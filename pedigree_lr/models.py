from __future__ import annotations
import sys
from copy import deepcopy
from dataclasses import dataclass, field
from datetime import timedelta, datetime
from decimal import Decimal
from io import StringIO
from itertools import combinations, permutations
from pathlib import Path
from random import Random
from typing import Mapping
import networkx as nx
import logging
import pandas as pd
import json

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stdout))

"""
This module contains classes and methods for representing and manipulating pedigrees.

The module provides data structures for markers, alleles, individuals, relationships, and pedigrees.
"""


@dataclass
class Marker:
    """
        Represents a genetic marker with a mutation rate and an optional number of copies.

        Attributes:
            name (str): The name of the marker.
            mutation_rate (float): The probability of mutation occurring at this marker.
            number_of_copies (int | None): The number of copies of this marker.
        """
    name: str
    mutation_rate: float
    number_of_copies: int | None = None


@dataclass
class Allele:
    """
        Represents an allele associated with a genetic marker.

        Attributes:
            marker (Marker): The associated genetic marker.
            value (int): The allele value.
            intermediate_value (int | None): An optional intermediate value.
            mutation_value (int | None): The number of mutation steps compared to parent.
            mutation_probability (float | None): The probability of mutation for this allele.
        """
    marker: Marker
    value: int
    intermediate_value: int | None = None
    mutation_value: int | None = None
    mutation_probability: float | None = None

    def __str__(self):
        return f"{self.value}" if self.intermediate_value is None else f"{self.value}.{self.intermediate_value}"


@dataclass
class Haplotype:
    """
        Represents a collection of alleles grouped by genetic markers.
        """

    def __init__(
            self
    ):
        """
            Initializes an empty haplotype.
            A haplotype is a collection of alleles grouped by genetic markers.
            Each marker can have multiple alleles (multi-copy), and the alleles are stored in a dictionary.
            """
        self.alleles: dict[str, list[Allele]] = {}

    def add_allele(
            self,
            marker: Marker,
            value: int,
            intermediate_value: int = None,
    ):
        """
                Adds an allele for a specific marker to the haplotype.

                Args:
                    marker (Marker): The genetic marker associated with the allele.
                    value (int): The allele value.
                    intermediate_value (int, optional): An intermediate allele value.
                """
        if marker.name not in self.alleles:
            self.alleles[marker.name] = []
        self.alleles[marker.name].append(Allele(marker, value, intermediate_value))

    def get_alleles_by_marker_name(
            self,
            marker_name: str,
    ) -> list[Allele]:
        """
                Retrieves a sorted list of alleles for a specific marker.

                Args:
                    marker_name (str): The name of the genetic marker.

                Returns:
                    list[Allele]: A sorted list of alleles associated with the marker.
                """
        return sorted(self.alleles.get(marker_name, []), key=lambda x: x.value)

    def __eq__(
            self,
            other: Haplotype
    ) -> bool:
        """
                Checks if two haplotypes are equal based on their alleles.

                Args:
                    other (Haplotype): The haplotype to compare with.

                Returns:
                    bool: True if the haplotypes are equivalent, False otherwise.
                """
        if self.alleles.keys() != other.alleles.keys():
            return False
        for marker, alleles in self.alleles.items():
            other_alleles = other.alleles.get(marker, [])
            if len(alleles) != len(other_alleles):
                return False
            # Compare alleles one-to-one
            for allele, other_allele in zip(sorted(alleles, key=lambda x: x.value),
                                            sorted(other_alleles, key=lambda x: x.value)):
                if (allele.value != other_allele.value or
                        allele.intermediate_value != other_allele.intermediate_value):
                    return False
        return True

    def allelic_difference(
            self,
            other: Haplotype
    ) -> int:
        """
                Calculates the allelic difference between two haplotypes.

                Args:
                    other (Haplotype): The haplotype to compare with.

                Returns:
                    int: The number of needed "mutations" to get from one haplotype to another.
                """
        if self.alleles.keys() != other.alleles.keys():
            return -1
        difference = 0
        for marker, alleles in self.alleles.items():
            other_alleles = other.alleles.get(marker, [])
            for allele, other_allele in zip(sorted(alleles, key=lambda x: x.value),
                                            sorted(other_alleles, key=lambda x: x.value)):
                difference += abs(allele.value - other_allele.value)
        return difference

    def __repr__(self):
        return f"Haplotype({self.alleles})"


@dataclass(frozen=False)
class Individual:
    """
        Represents an individual with a haplotype and additional metadata.
        """
    id: int
    name: str
    haplotype: Haplotype = field(default_factory=lambda: Haplotype())
    haplotype_class: str = "unknown"
    exclude: bool = False
    picking_probability: Decimal | None = None

    def add_allele(
            self,
            marker: Marker,
            value: int,
            intermediate_value: int = None,
    ):
        """
                Adds an allele to the individual's haplotype.
                """
        self.haplotype.add_allele(marker, value, intermediate_value)

    def get_alleles_by_marker_name(
            self,
            marker_name: str,
    ) -> list[Allele] | None:
        """
                Retrieves alleles by marker name from the individual's haplotype.
                """
        return self.haplotype.get_alleles_by_marker_name(marker_name)


@dataclass
class Relationship:
    """
        Represents a parent-child relationship between individuals.

        Attributes:
            parent_id (int): The ID of the parent individual.
            child_id (int): The ID of the child individual.
            edge_class (str): The type of relationship (default: "unknown").
        """
    parent_id: int
    child_id: int
    edge_class: str = "unknown"


class MarkerSet:
    """
        A collection of genetic markers with utility functions for management.
        """

    def __init__(self):
        """
                Initializes an empty marker set.
                """
        self.markers: list[Marker] = []

    def add_marker(
            self,
            marker: Marker,
    ):
        """
                Adds a marker to the marker set.
                """
        self.markers.append(marker)

    def get_marker_by_name(
            self,
            marker_name: str,
    ) -> Marker | None:
        """
                Retrieves a marker by its name.
                """
        for marker in self.markers:
            if marker.name == marker_name:
                return marker
        return None

    def read_marker_set_from_file(
            self,
            file,
    ):
        """
                Reads marker set data from a file.
                """
        header = next(file)
        if header.strip() != "marker,mutation_rate":
            logger.error(f"Invalid header in marker set file: {header}")

        for line in file:
            if "," not in line:
                logger.error(f"Invalid separator in marker set file: {line}")
                continue

            try:
                marker_name, mutation_rate = line.split(",")
            except ValueError:
                logger.error(f"Invalid line in marker set file: {line}")
                continue

            marker_name = marker_name.strip()
            mutation_rate = mutation_rate.strip()

            if not marker_name:
                logger.error(f"Empty marker name in line: {line}")
                continue

            try:
                mutation_rate = float(mutation_rate)
            except ValueError:
                logger.error(f"Invalid mutation rate value: {mutation_rate} in line: {line}")
                continue

            if mutation_rate < 0 or mutation_rate > 1:
                logger.error(f"Mutation rate out of bounds: {mutation_rate} in line: {line}")
                continue

            self.add_marker(Marker(marker_name, mutation_rate))

    def load_markers_from_database(
            self,
            markers: list[str],
    ):
        mutation_rates_path = Path(__file__).resolve().parent.parent / "data" / "mutation_rates.csv"
        mutation_rates_df = pd.read_csv(mutation_rates_path, index_col=0, header=0, names=["Marker", "Mutation rate"])

        for marker_name in markers:
            if marker_name not in mutation_rates_df.index:
                logger.error(f"Marker {marker_name} not found in mutation rates database.")
                continue
            mutation_rate = mutation_rates_df.loc[marker_name, "Mutation rate"]
            self.add_marker(Marker(marker_name, mutation_rate))


@dataclass
class Pedigree:
    """
        Represents a pedigree structure containing individuals and their relationships.

        Attributes:
            individuals (list[Individual]): List of individuals in the pedigree.
            relationships (list[Relationship]): List of parent-child relationships.
        """
    individuals: list[Individual] = field(default_factory=lambda: [])
    relationships: list[Relationship] = field(default_factory=lambda: [])
    picking_probabilities: dict[int, dict[tuple[int], Decimal]] = field(default_factory=lambda: {})

    def add_individual(
            self,
            individual_id: int | str,
            name: str,
    ):
        """
                Adds a new individual to the pedigree.

                Args:
                    individual_id (int or str): Unique identifier for the individual.
                    name (str): Name of the individual.
                """
        if any(individual.name == name for individual in self.individuals):
            logger.warning(f"Individual with name {name} already exists.")
        individual = Individual(individual_id, name)
        self.individuals.append(individual)

    def remove_individual(
            self,
            individual_id: int,
    ):
        """
                Removes an individual and their relationships from the pedigree.

                Args:
                    individual_id (int): Unique identifier of the individual to remove.
                """
        individual = self.get_individual_by_id(individual_id)
        if individual:
            self.individuals.remove(individual)
            self.relationships = [
                relationship
                for relationship in self.relationships
                if relationship.parent_id != individual_id
                   and relationship.child_id != individual_id
            ]

    def add_relationship(
            self,
            parent_id: int | str,
            child_id: int | str,
    ):
        """
                Adds a parent-child relationship to the pedigree.

                Args:
                    parent_id (int): ID of the parent individual.
                    child_id (int): ID of the child individual.
                """
        relationship = Relationship(parent_id, child_id)
        self.relationships.append(relationship)

    def read_tgf(
            self,
            file,
    ):
        """
                Reads a pedigree structure from a TGF (Trivial Graph Format) file.

                Args:
                    file: File-like object containing the pedigree data.
                """
        current_section = "node"
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line == "#":  # Switch to edge section
                current_section = "edge"
                continue
            if current_section == "node":
                individual_values = line.split()
                if len(individual_values) == 2:
                    individual_id, individual_name = line.split()
                elif len(individual_values) == 1:
                    individual_id = individual_values[0]
                    individual_name = individual_values[0]
                else:
                    logger.error(f"Invalid line in TGF file: {line}")
                    continue

                self.add_individual(str(individual_id), str(individual_name))
            elif current_section == "edge":
                parent_id, child_id = line.split()
                self.add_relationship(str(parent_id), str(child_id))

    def read_ped(
            self,
            file,
    ):
        """
                Reads a pedigree structure from a PED file (used in genetic studies).

                Args:
                    file: File-like object containing the pedigree data.
                """
        relationships = []
        for line in file:
            line = line.strip()
            if not line:
                continue
            family_id, individual_id, paternal_id, maternal_id, sex, phenotype = line.split()
            if sex == "2":
                continue
            elif sex == "1":
                self.add_individual(int(individual_id), str(individual_id))
                relationships.append((int(paternal_id), int(individual_id)))
        for paternal_id, child_id in relationships:
            if paternal_id == 0 or child_id == 0:
                continue
            self.add_relationship(paternal_id, child_id)

    def read_pedigree_from_file(
            self,
            file: StringIO,
            file_extension: str,
    ):
        """
                Reads a pedigree from a file based on its extension.

                Args:
                    file (StringIO): File-like object containing pedigree data.
                    file_extension (str): File extension indicating format (".tgf" or ".ped").
                """
        if file_extension == ".tgf":
            self.read_tgf(file)
        elif file_extension == ".ped":
            self.read_ped(file)

    def write_to_tgf(
            self
    ) -> bytes:
        """
                Serializes the pedigree into a TGF-formatted byte string.

                Returns:
                    bytes: TGF-formatted pedigree representation.
                """
        lines = []
        for individual in self.individuals:
            lines.append(f"{individual.id} {individual.name}")
        lines.append("#")
        for relationship in self.relationships:
            lines.append(f"{relationship.parent_id} {relationship.child_id}")
        return "\n".join(lines).encode("utf-8")

    def read_known_haplotypes_from_file(
            self,
            file: StringIO,
            marker_set: MarkerSet,
    ):
        json_data = json.load(file)

        for individual_name in json_data.keys():
            try:
                individual = self.get_individual_by_name(individual_name)
            except ValueError:
                logger.error(f"Individual {individual_name} not found in pedigree")
                return
            individual.haplotype_class = "known"
            individual_alleles = json_data[individual_name]

            for marker_name, values in individual_alleles.items():
                try:
                    marker = marker_set.get_marker_by_name(marker_name)
                except ValueError:
                    logger.error(f"Marker {marker_name} not found in marker set")
                    continue

                if not marker:
                    logger.error(f"Marker {marker_name} not found in marker set")
                    continue

                alleles = values.split(";")  # Use ";" as delimiter for multiple alleles
                number_of_copies = len(alleles)
                if not marker.number_of_copies:
                    marker.number_of_copies = number_of_copies
                elif marker.number_of_copies != number_of_copies:
                    logger.error(f"Number of copies mismatch for marker {marker_name}")
                    continue

                for allele in alleles:
                    if "." in allele:  # Intermediate allele
                        allele, intermediate_value = allele.split(".")
                        try:
                            allele = int(allele)
                            intermediate_value = int(intermediate_value)
                        except ValueError:
                            logger.error(f"Invalid allele or intermediate value: {allele}.{intermediate_value}")
                        individual.add_allele(marker, allele, intermediate_value)
                    else:
                        individual.add_allele(marker, int(allele))

    def get_individual_by_name(
            self,
            individual_name: str,
    ) -> Individual | None:
        for individual in self.individuals:
            if individual.name == individual_name:
                return individual
        return None

    def get_individual_by_id(
            self,
            individual_id: int,
    ) -> Individual | None:
        for individual in self.individuals:
            if individual.id == individual_id:
                return individual
        return None

    def get_unknown_individuals(
            self
    ) -> list[Individual]:
        return [
            individual
            for individual in self.individuals
            if individual.haplotype_class == "unknown"
        ]

    def get_known_individuals(
            self
    ) -> list[Individual]:
        return [
            individual
            for individual in self.individuals
            if individual.haplotype_class == "known"
        ]

    def set_suspect(
            self,
            suspect_name: str,
    ):
        previous_suspect = self.get_suspect()
        if previous_suspect:
            previous_suspect.haplotype_class = "known"

        suspect_individual = self.get_individual_by_name(suspect_name)
        if suspect_individual:
            suspect_individual.haplotype_class = "suspect"

    def reroot_pedigree(
            self,
            new_root_name: str,
    ):
        previous_root = self.get_suspect()
        if previous_root:
            previous_root.haplotype_class = "known"
        try:
            new_root = self.get_individual_by_name(new_root_name)
        except ValueError:
            logger.error(f"Individual {new_root_name} not found in pedigree")
            return
        new_root.haplotype_class = "suspect"
        current_graph = create_nx_graph(self).to_undirected()
        rerooted_graph = nx.DiGraph(nx.dfs_tree(current_graph, source=new_root.id))
        self.relationships = [
            Relationship(parent_id, child_id)
            for parent_id, child_id in rerooted_graph.edges()
        ]

    def exclude_individuals(
            self,
            excluded_individuals: list[str],
    ):
        for individual in self.individuals:
            individual.exclude = False

        for individual_name in excluded_individuals:
            individual = self.get_individual_by_name(individual_name)
            if individual:
                individual.exclude = True

    def get_level_order_traversal(
            self,
            source_name: str,
    ) -> list[Individual]:
        source = self.get_individual_by_name(source_name)
        ordered_individuals = []
        for level in nx.bfs_layers(create_nx_graph(self), sources=source.id):
            for individual_id in level:
                individual = self.get_individual_by_id(individual_id)
                if individual:
                    ordered_individuals.append(individual)
        return ordered_individuals

    def calculate_allele_probabilities(
            self,
            marker_set: MarkerSet,
            two_step_mutation_factor: float,
    ):
        for relationship in self.relationships:
            parent = self.get_individual_by_id(relationship.parent_id)
            child = self.get_individual_by_id(relationship.child_id)
            for marker in marker_set.markers:
                parent_alleles = parent.get_alleles_by_marker_name(marker.name)
                child_alleles = child.get_alleles_by_marker_name(marker.name)

                mutation_probability = calculate_mutation_probability(parent_alleles, child_alleles, marker,
                                                                      two_step_mutation_factor)

                for child_allele in child_alleles:
                    child_allele.mutation_probability += mutation_probability

    def to_string(
            self
    ):
        lines = ["Pedigree"]

        nx_graph = create_nx_graph(self)
        root = self.get_individual_by_id(1)
        node_depths = {
            node: depth
            for depth, nodes in enumerate(nx.bfs_layers(nx_graph, sources=root.id))
            for node in nodes
        }

        for node_id in nx.dfs_preorder_nodes(nx_graph):
            prefix = "\t" * node_depths[node_id]
            individual = self.get_individual_by_id(node_id)
            lines.append(f"{prefix}{individual.name}, {individual.haplotype_class}")
            for allele in individual.haplotype.alleles.values():
                lines.append(f"{prefix}.{allele}")

        return "\n".join(lines)

    def get_suspect(
            self
    ):
        for individual in self.individuals:
            if individual.haplotype_class == "suspect":
                return individual
        return None

    def check_known_haplotypes(
            self
    ):
        known_haplotypes = [
            individual
            for individual in self.individuals
            if individual.haplotype_class == "known"
        ]

        if not known_haplotypes:
            return

        marker_names = known_haplotypes[0].haplotype.alleles.keys()
        for individual in known_haplotypes[1:]:
            if individual.haplotype.alleles.keys() != marker_names:
                logger.error("Known haplotypes have different markers")
                return

    def check_pedigree_structure(
            self
    ):
        graph = create_nx_graph(self)
        if not nx.is_directed_acyclic_graph(graph):
            logger.error("Pedigree contains cycles")
        if not nx.is_connected(graph.to_undirected()):
            logger.error("Pedigree contains disconnected components")
        if not nx.is_tree(graph):
            logger.error("Pedigree is not a tree")
        pass

    def calculate_picking_probabilities(
            self,
            marker_set: MarkerSet,
    ):
        pedigree_deep_copy = deepcopy(self)
        suspect_haplotype = pedigree_deep_copy.get_suspect().haplotype
        p = nx.shortest_path(create_nx_graph(pedigree_deep_copy), source=pedigree_deep_copy.get_suspect().id)
        unknown_individuals = [i.id for i in pedigree_deep_copy.get_unknown_individuals()]  # if not i.exclude
        l_combinations_dict = {}

        # set all unknown haplotypes to the closest (ancestor) haplotype
        for individual in unknown_individuals:
            shortest_path = p[individual]
            for i, node_id in enumerate(shortest_path):
                individual = pedigree_deep_copy.get_individual_by_id(node_id)
                if individual.haplotype_class == "unknown":
                    individual.haplotype = deepcopy(
                        pedigree_deep_copy.get_individual_by_id(shortest_path[i - 1]).haplotype)
                    individual.haplotype_class = "estimated"

        for l in range(1, 2):
            l_combinations = list(combinations(unknown_individuals, l))
            l_combinations_dict[l] = {}

            for l_comb in l_combinations:
                l_pedigree_deepcopy = deepcopy(pedigree_deep_copy)
                total_comb_needed_mutations = 0
                l_comb_tuple = tuple(sorted([individual for individual in l_comb]))

                for comb_individual_id in l_comb:
                    comb_individual = l_pedigree_deepcopy.get_individual_by_id(comb_individual_id)
                    comb_individual.haplotype_class = "fixed"
                    comb_individual.haplotype = deepcopy(suspect_haplotype)

                for unknown_ind_id in unknown_individuals:
                    parent_ind_id = p[unknown_ind_id][-2]
                    unknown_ind = l_pedigree_deepcopy.get_individual_by_id(unknown_ind_id)
                    parent_ind = l_pedigree_deepcopy.get_individual_by_id(parent_ind_id)
                    number_of_mutations = unknown_ind.haplotype.allelic_difference(parent_ind.haplotype)
                    if number_of_mutations == -1:
                        print("error")
                    else:
                        total_comb_needed_mutations += number_of_mutations

                # print(f"l_comb: {l_comb_tuple}\t total_comb_needed_mutations: {total_comb_needed_mutations}")
                total_comb_needed_mutations += 1
                l_combinations_dict[l][l_comb_tuple] = Decimal(total_comb_needed_mutations)

            values = l_combinations_dict[l].values()
            max_value = max(values)
            for l_comb_tuple, total_comb_needed_mutations in l_combinations_dict[l].items():
                l_combinations_dict[l][l_comb_tuple] = Decimal(
                    ((max_value - l_combinations_dict[l][l_comb_tuple]) + 1) / (max_value + 1))
                # print(f"l_comb: {l_comb_tuple}\t picking probability: {l_combinations_dict[l][l_comb_tuple]}")
        self.picking_probabilities = l_combinations_dict

    def extend_pedigree(
            self,
    ):
        root = list(nx.topological_sort(create_nx_graph(self)))[0]
        generations = list(nx.bfs_layers(create_nx_graph(self), root))

        new_root_id = 0
        while new_root_id in [individual.id for individual in self.individuals]:
            new_root_id += 1

        self.add_individual(new_root_id, "new_root")
        self.add_relationship(new_root_id, root)

        previous_parent = new_root_id
        for i in range(1, len(generations) + 1):
            new_child = 0
            while new_child in [individual.id for individual in self.individuals]:
                new_child += 1

            self.add_individual(new_child, f"new_child_{i}")
            self.get_individual_by_id(new_child).haplotype_class = "unknown"
            self.add_relationship(previous_parent, new_child)
            previous_parent = new_child

        last_child_name = f"new_child_{len(generations)}"
        return last_child_name


@dataclass(frozen=True)
class IterationResult:
    probability: Decimal
    edge_probabilities: Mapping[tuple[int, int], Decimal]
    mutated_haplotypes: Mapping[int, Haplotype]
    fixed_individuals_ids: set[int] | None


@dataclass(frozen=True)
class SimulationParameters:
    number_of_iterations: int
    two_step_mutation_factor: float
    stability_window: int
    stability_min_iterations: int
    stability_threshold: float
    model_validity_threshold: float
    simulation_name: str
    number_of_threads: int
    results_path: Path
    random_seed: int | None = None


@dataclass(frozen=True)
class SimulationResult:
    pedigree: Pedigree
    marker_set: MarkerSet
    suspect_name: str
    simulation_parameters: SimulationParameters
    random: Random
    average_pedigree_probability: Decimal
    proposal_distribution: Mapping[int, Decimal]
    l_needed_iterations: Mapping[int, list[int]]
    l_model_probabilities: Mapping[int, list[Decimal]]
    l_model_validities: Mapping[int, bool]
    outside_match_probability: Decimal
    outside_needed_iterations: list[int]
    outside_model_probabilities: list[Decimal]
    outside_model_is_valid: bool
    run_time_pedigree_probability: timedelta
    run_time_proposal_distribution: timedelta
    total_run_time: timedelta

    # def download_results(
    #         self,
    #         simulation_parameters: SimulationParameters,
    # ) -> bytes:
    #     from pedigree_lr.reporting import create_report_bytes
    #     return create_report_bytes(self, simulation_parameters=simulation_parameters)


def create_nx_graph(
        pedigree: Pedigree
) -> nx.DiGraph:
    graph = nx.DiGraph()
    for individual in pedigree.individuals:
        graph.add_node(individual.id)
    for relationship in pedigree.relationships:
        graph.add_edge(relationship.parent_id, relationship.child_id)
    return graph


def get_mutation_probability(
        mutation_rate: float,
        mutation_value: float,
        two_step_mutation_factor: float
) -> float:
    no_mutation_rate = 1 - mutation_rate
    single_step_mutation_rate = (mutation_rate * (1 - two_step_mutation_factor)) / 2
    two_step_mutation_rate = (mutation_rate * two_step_mutation_factor) / 2

    if mutation_value == 0:
        return no_mutation_rate
    elif mutation_value == 1 or mutation_value == -1:
        return single_step_mutation_rate
    elif mutation_value == 2 or mutation_value == -2:
        return two_step_mutation_rate
    else:
        return 0.0


def get_single_copy_mutation_rate(
        mutation_rate: float,
        number_of_copies: int
) -> float:
    """
    P(at least one mutation) = mu_all - (1 - mu_1)^n
    mu_1 = 1 - (1 - mu_all)^(1/n)
    """
    return 1 - (1 - mutation_rate) ** (1 / number_of_copies)


def calculate_mutation_probability(
        parent_alleles: list[Allele],
        child_alleles: list[Allele],
        marker: Marker,
        two_step_mutation_factor: float
) -> Decimal:
    mutation_probability = Decimal(0)

    combs = [list(zip(parent_alleles, perm)) for perm in permutations(child_alleles)]

    # TODO: make more efficient by only calculating mutation probability for unique combinations
    for combination in combs:
        combination_probability = Decimal(1)
        for parent_allele, child_allele in combination:
            if child_allele.intermediate_value != parent_allele.intermediate_value:
                combination_probability = Decimal(0)  # No mutation for intermediate values mismatch
                break
            else:
                mutation_value = child_allele.value - parent_allele.value
                combination_probability *= Decimal(
                    get_mutation_probability(
                        get_single_copy_mutation_rate(marker.mutation_rate,
                                                      marker.number_of_copies),
                        mutation_value,
                        two_step_mutation_factor
                    )
                )
        mutation_probability += combination_probability

    return mutation_probability
