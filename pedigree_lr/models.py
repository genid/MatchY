from __future__ import annotations
import sys
from dataclasses import dataclass, field
from datetime import timedelta, datetime
from decimal import Decimal
from io import StringIO
from itertools import permutations
from pathlib import Path
from random import Random
from typing import Mapping
import networkx as nx
import logging

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

    def __eq__(self,
               other: "Haplotype"
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
        self.markers = []

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
        if header.strip() != "locus,mutation_rate":
            logger.error(f"Invalid header in marker set file: {header}")
            return

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

    def add_individual(
            self,
            individual_id: int,
            name: str,
    ):
        """
                Adds a new individual to the pedigree.

                Args:
                    individual_id (int): Unique identifier for the individual.
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
            parent_id: int,
            child_id: int,
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

                self.add_individual(int(individual_id), str(individual_name))
            elif current_section == "edge":
                parent_id, child_id = line.split()
                self.add_relationship(int(parent_id), int(child_id))

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

    def read_known_haplotype_from_file(
            self,
            individual_name: str,
            file,
            marker_set: MarkerSet,
    ):
        try:
            individual = self.get_individual_by_name(individual_name)
        except ValueError:
            logger.error(f"Individual {individual_name} not found in pedigree")
            return
        individual.haplotype_class = "known"
        header = next(file)  # Skip header
        if header.strip() != "marker,allele":
            logger.error(f"Invalid header in known haplotype file: {header}")

        for line in file:
            marker_name, values = line.split(",")
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
    ):
        for relationship in self.relationships:
            parent = self.get_individual_by_id(relationship.parent_id)
            child = self.get_individual_by_id(relationship.child_id)
            for marker in marker_set.markers:
                parent_alleles = parent.get_alleles_by_marker_name(marker.name)
                child_alleles = child.get_alleles_by_marker_name(marker.name)

                mutation_probability = calculate_mutation_probability(parent_alleles, child_alleles, marker)

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
        unknown_individuals = self.get_unknown_individuals()
        p = nx.shortest_path(create_nx_graph(self), source=self.get_suspect().id)

        for individual in unknown_individuals:
            shortest_path = p[individual.id][::-1]  # list includes source [0] and target [-1]

            intermediate_mutated_nodes = []
            last_haplotype = self.get_suspect().haplotype  # unknown individual's haplotype is suspect's haplotype
            steps = 1

            for i in range(1, len(shortest_path)):
                node = shortest_path[i]

                if self.get_individual_by_id(node).haplotype_class != "unknown":
                    haplotype = self.get_individual_by_id(node).haplotype
                    intermediate_mutated_nodes.append({
                        "source-target-haplotype": (last_haplotype, haplotype),
                        "steps": steps,
                    })
                    last_haplotype = haplotype
                    steps = 1
                else:
                    steps += 1

            total_mutation_probability = Decimal(1)
            for jump in intermediate_mutated_nodes:
                source_haplotype, target_haplotype = jump["source-target-haplotype"]
                steps = jump["steps"]

                mutation_probability = Decimal(1)
                for marker_name in source_haplotype.alleles.keys():
                    source_alleles = source_haplotype.get_alleles_by_marker_name(marker_name)
                    target_alleles = target_haplotype.get_alleles_by_marker_name(marker_name)
                    marker = marker_set.get_marker_by_name(marker_name)
                    mutation_probability *= Decimal(
                        calculate_mutation_probability(source_alleles, target_alleles, marker)
                    )

                mutation_probability **= Decimal(steps)
                total_mutation_probability *= mutation_probability

            individual.picking_probability = total_mutation_probability

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


@dataclass(frozen=True)
class SimulationResult:
    pedigree: Pedigree
    marker_set: MarkerSet
    suspect_name: str
    number_of_iterations: int
    random: Random
    average_pedigree_probability: Decimal
    proposal_distribution: Mapping[int, Decimal]
    outside_match_probability: Decimal
    run_time_pedigree_probability: timedelta
    run_time_proposal_distribution: timedelta
    total_run_time: timedelta

    def download_results(
            self,
            random_seed: int,
    ) -> bytes:
        bytes_data = StringIO()

        bytes_data.write("Simulation results\n")
        bytes_data.write("match-Y version 1.0.0\n\n")  # TODO: Remove hard coded version

        bytes_data.write(f"Date and time of report: \t{datetime.now()}\n")
        bytes_data.write(f"Number of iterations: \t{self.number_of_iterations}\n")
        bytes_data.write(f"Random seed: \t{random_seed}\n\n")

        bytes_data.write("Marker set with mutation rate\n")
        for marker in self.marker_set.markers:
            bytes_data.write(f"{marker.name}: {marker.mutation_rate}\n")

        bytes_data.write(f"\nNumber of nodes in pedigree: \t{len(self.pedigree.individuals)}\n")
        bytes_data.write(f"Number of edges in pedigree: \t{len(self.pedigree.relationships)}\n")

        bytes_data.write("\nIndividuals:\n")
        bytes_data.write("ID, Name, Haplotype class\n")
        for individual in self.pedigree.individuals:
            bytes_data.write(f"{individual.id}, {individual.name}, {individual.haplotype_class}\n")

        bytes_data.write("\nRelationships\n")
        bytes_data.write("Parent ID -> Child ID\n")
        for relationship in self.pedigree.relationships:
            bytes_data.write(f"{relationship.parent_id} -> {relationship.child_id}\n")

        bytes_data.write(f"\nSuspect: \t{self.suspect_name}\n\n")

        bytes_data.write(f"Average pedigree probability: \t{self.average_pedigree_probability}\n\n")
        bytes_data.write(f"Run time average pedigree probability: \t{self.run_time_pedigree_probability}\n")
        bytes_data.write(f"Run time proposal distribution: \t{self.run_time_proposal_distribution}\n")
        bytes_data.write(f"Total run time: \t{self.total_run_time}\n\n")

        bytes_data.write("\nMatch probabilities\n")

        for key in sorted(self.proposal_distribution.keys()):
            bytes_data.write(f"{key}: {self.proposal_distribution[key]:.4f}\n")

        bytes_data.write(f"\nOutside match probability: \t{self.outside_match_probability:.4f}\n")

        return bytes_data.getvalue().encode("utf-8")


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
) -> float:
    if mutation_value == 0:
        return 1 - mutation_rate
    elif mutation_value == 1 or mutation_value == -1:
        return (mutation_rate * 0.97) / 2
    elif mutation_value == 2 or mutation_value == -2:  # TODO: Remove hard coded value for 2-step mutation
        return (mutation_rate * 0.03) / 2
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
        marker: Marker
) -> Decimal:
    mutation_probability = Decimal(0)

    combinations = [list(zip(parent_alleles, perm)) for perm in permutations(child_alleles)]

    # TODO: make more efficient by only calculating mutation probability for unique combinations
    for combination in combinations:
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
                        mutation_value
                    )
                )
        mutation_probability += combination_probability

    return mutation_probability
