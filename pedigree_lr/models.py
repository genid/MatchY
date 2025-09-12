from __future__ import annotations
import sys
from copy import deepcopy
from dataclasses import dataclass, field
from datetime import timedelta
from decimal import Decimal
from io import StringIO
from itertools import permutations, product
from pathlib import Path
from typing import Mapping, Generator
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

class InvalidAveragePedigreeProbability(Exception):
    pass


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

    def __repr__(self):
        return f"{self.name} (mutation rate: {self.mutation_rate}, copies: {self.number_of_copies})"

    def __eq__(self, other):
        if not isinstance(other, Marker):
            return NotImplemented
        return (
                self.name == other.name and
                self.mutation_rate == other.mutation_rate and
                self.number_of_copies == other.number_of_copies
        )

    def __hash__(self):
        return hash((self.name, self.mutation_rate, self.number_of_copies))


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

    # Allele string representation, e.g. "10" or "10.2"
    def __repr__(self):
        return f"{self.value}" if self.intermediate_value is None else f"{self.value}.{self.intermediate_value}"

    def __eq__(self, other):
        if not isinstance(other, Allele):
            return NotImplemented
        return (self.marker, self.value, self.intermediate_value) == (other.marker, other.value,
                                                                      other.intermediate_value)

    def __hash__(self):
        return hash((self.marker, self.value, self.intermediate_value))

    def __lt__(self, other):
        if not isinstance(other, Allele):
            return NotImplemented
        return (self.value, self.intermediate_value or 0) < (other.value, other.intermediate_value or 0)


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
    id: str
    name: str
    haplotype: Haplotype = field(default_factory=lambda: Haplotype())
    haplotype_class: str = "unknown"
    exclude: bool = False
    picking_probability: Decimal | None = None
    closest_known_individuals: list[Individual] = field(default_factory=lambda: [])
    closest_known_distance: int | None = None

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

    def __str__(self):
        """
                Returns a string representation of the individual.
                """
        return f"Individual(id={self.id}, name={self.name}, haplotype_class={self.haplotype_class})"

    def __repr__(self):
        """
                Returns a detailed string representation of the individual.
                """
        return f"Individual(id={self.id}, name={self.name}, haplotype_class={self.haplotype_class})"

@dataclass
class Relationship:
    """
        Represents a parent-child relationship between individuals.

        Attributes:
            parent_id (int): The ID of the parent individual.
            child_id (int): The ID of the child individual.
            edge_class (str): The type of relationship (default: "unknown").
        """
    parent_id: str
    child_id: str
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
    picking_probabilities: dict[str, Decimal] = field(default_factory=lambda: {})

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
            individual_id: str,
    ):
        """
        Removes an individual and all their descendants (children, grandchildren, ...)
        and the corresponding relationships from the pedigree.

        Args:
            individual_id (str): Unique identifier of the individual to remove.
        """
        # Find the individual object
        individual = self.get_individual_by_id(individual_id)

        if not individual:
            return  # Nothing to remove if individual doesn't exist

        # Find all direct children of this individual
        children = [
            relationship.child_id
            for relationship in self.relationships
            if str(relationship.parent_id) == str(individual_id)
        ]

        # Recursively remove all children
        for child_id in children:
            self.remove_individual(child_id)

        # Finally, remove the individual itself
        if individual in self.individuals:
            self.individuals.remove(individual)

        # Remove all relationships involving this individual
        self.relationships = [
            relationship
            for relationship in self.relationships
            if (
                    str(relationship.parent_id) != str(individual_id)
                    and str(relationship.child_id) != str(individual_id)
            )
        ]

    def add_relationship(
            self,
            parent_id: str,
            child_id: str,
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
            individual = self.get_individual_by_name(individual_name)
            if not individual:
                logger.error(f"Individual {individual_name} not found in pedigree. Skipping known haplotype assignment.")
                return
            individual.haplotype_class = "known"
            individual_alleles = json_data[individual_name]

            for marker_name, values in individual_alleles.items():
                marker = marker_set.get_marker_by_name(marker_name)
                if not marker:
                    logger.error(f"Marker {marker_name} not found in marker set. Skipping allele assignment.")
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
            if str(individual.name) == str(individual_name):
                return individual
        return None

    def get_individual_by_id(
            self,
            individual_id: int | str,
    ) -> Individual | None:
        for individual in self.individuals:
            if str(individual.id) == str(individual_id):
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
        new_root = self.get_individual_by_name(new_root_name)
        if not new_root:
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

    def get_closest_known_individuals(
            self
    ):
        unknown_individuals = self.get_unknown_individuals()
        known_individuals = self.get_known_individuals() + [self.get_suspect()]

        undirected_graph = create_nx_graph(self).to_undirected()

        for unknown_individual in unknown_individuals:
            unknown_known_distance_dict = {}
            for known_individual in known_individuals:
                shortest_path = nx.shortest_path(undirected_graph, source=unknown_individual.id, target=known_individual.id)
                unknown_known_distance_dict[(unknown_individual.id, shortest_path[1], known_individual.id)] = len(shortest_path) - 1
            lowest_distance = min(unknown_known_distance_dict.values())
            closest_known_individual_ids = [(unknown_id, via_id)
                for (unknown_id, via_id, known_id), distance in unknown_known_distance_dict.items()
                if distance == lowest_distance
            ]

            for (unknown_id, via_id, known_id), distance in unknown_known_distance_dict.items():
                unknown_ind = self.get_individual_by_id(unknown_id)
                via_ind = self.get_individual_by_id(via_id)

                if (unknown_id, via_id) in closest_known_individual_ids:
                    if via_ind not in unknown_ind.closest_known_individuals:
                        unknown_ind.closest_known_individuals.append(via_ind)
                        unknown_ind.closest_known_distance = distance

    def calculate_allele_probabilities(
            self,
            marker_set: MarkerSet,
            two_step_mutation_factor: float,
            is_average_pedigree: bool = False,
    ):
        for relationship in self.relationships:
            parent = self.get_individual_by_id(relationship.parent_id)
            child = self.get_individual_by_id(relationship.child_id)
            for marker in marker_set.markers:
                parent_alleles = parent.get_alleles_by_marker_name(marker.name)
                child_alleles = child.get_alleles_by_marker_name(marker.name)

                mutation_probability = calculate_mutation_probability(
                    parent_alleles=parent_alleles,
                    child_alleles=child_alleles,
                    marker=marker,
                    two_step_mutation_factor=two_step_mutation_factor,
                    is_average_pedigree=is_average_pedigree
                )

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
    ):
        """
            Calculates a priori picking probabilities for all individuals with unknown haplotypes in the pedigree.

            This method estimates the likelihood that each unknown individual could plausibly be assigned
            the suspect's haplotype. It does so by:

            1. Estimating unknown haplotypes based on the closest known ancestor along the shortest path
               from the suspect.
            2. For each unknown individual, simulating the scenario where that individual is fixed to the
               suspect's haplotype and all others remain estimated.
            3. Counting the total number of mutations required across the pedigree for that configuration.
            4. Transforming the mutation counts into normalized picking probabilities, where fewer required
               mutations implies a higher probability.

            The resulting probabilities are stored in `self.picking_probabilities`, with each key being the
            individual's ID and the value being a Decimal between 0 and 1, representing the relative likelihood
            of being assigned the suspect haplotype.

            Note:
                - This method modifies internal state (`self.picking_probabilities`) but does not affect
                  the actual haplotypes in `self`.
                - Assumes that all individuals form a connected graph and that a suspect has been set.

            Raises:
                - Prints a warning if allelic difference computation fails (returns -1).
            """
        pedigree_deep_copy = deepcopy(self)
        suspect_haplotype = pedigree_deep_copy.get_suspect().haplotype
        if not suspect_haplotype:
            logger.error("Suspect does not exist or does not have a known haplotype.")
            return
        p = nx.shortest_path(create_nx_graph(pedigree_deep_copy), source=pedigree_deep_copy.get_suspect().id)
        unknown_individual_ids = [i.id for i in pedigree_deep_copy.get_unknown_individuals()]
        picking_probabilities = {}

        # Set all unknown haplotypes to the closest (ancestor) haplotype
        for individual in unknown_individual_ids:
            shortest_path = p[individual]
            for i, node_id in enumerate(shortest_path):
                individual = pedigree_deep_copy.get_individual_by_id(node_id)
                if individual.haplotype_class == "unknown":
                    individual.haplotype = deepcopy(
                        pedigree_deep_copy.get_individual_by_id(shortest_path[i - 1]).haplotype)
                    individual.haplotype_class = "estimated"

        for unknown_individual_id in unknown_individual_ids:
            pedigree_deepcopy = deepcopy(pedigree_deep_copy)
            total_needed_mutations = 0

            unknown_individual = pedigree_deepcopy.get_individual_by_id(unknown_individual_id)
            unknown_individual.haplotype_class = "fixed"
            unknown_individual.haplotype = deepcopy(suspect_haplotype)

            for unknown_ind_id in unknown_individual_ids:
                parent_ind_id = p[unknown_ind_id][-2]
                unknown_ind = pedigree_deepcopy.get_individual_by_id(unknown_ind_id)
                parent_ind = pedigree_deepcopy.get_individual_by_id(parent_ind_id)
                number_of_mutations = unknown_ind.haplotype.allelic_difference(parent_ind.haplotype)
                if number_of_mutations == -1:
                    logger.error("Allelic difference computation failed.")
                else:
                    total_needed_mutations += number_of_mutations

            total_needed_mutations += 1
            picking_probabilities[unknown_individual_id] = Decimal(total_needed_mutations)

        values = picking_probabilities.values()
        max_value = max(values)
        for unknown_individual_id, total_needed_mutations in picking_probabilities.items():
            picking_probabilities[unknown_individual_id] = Decimal(
                ((max_value - picking_probabilities[unknown_individual_id]) + 1) / (max_value + 1))

        self.picking_probabilities = picking_probabilities

    def extend_pedigree(
            self,
    ):
        G = create_nx_graph(self)
        root = [n for n,d in G.in_degree() if d==0][0]  # Find the root node (individual with no parents)
        generations = dict(enumerate(nx.bfs_layers(G, root)))

        highest_level_with_known_individual = None
        for level, individual_ids in generations.items():
            for individual_id in individual_ids:
                individual = self.get_individual_by_id(individual_id)
                if individual.haplotype_class == "known":
                    highest_level_with_known_individual = level + 1
                    break
            if highest_level_with_known_individual:
                break

        if not highest_level_with_known_individual:
            highest_level_with_known_individual = len(generations.keys())

        new_root_id = 0
        while str(new_root_id) in [str(individual.id) for individual in self.individuals]:
            new_root_id += 1  # Find a new unique ID for the root

        self.add_individual(new_root_id, "new_root")
        self.add_relationship(new_root_id, root)

        previous_parent = new_root_id
        for i in range(1, highest_level_with_known_individual + 1):
            new_child = 0
            while str(new_child) in [str(individual.id) for individual in self.individuals]:
                new_child += 1  # Find a new unique ID for the child

            self.add_individual(new_child, f"new_child_{i}")
            self.get_individual_by_id(new_child).haplotype_class = "unknown"
            self.add_relationship(previous_parent, new_child)
            previous_parent = new_child

        last_child_name = f"new_child_{highest_level_with_known_individual}"
        return last_child_name

    def remove_irrelevant_individuals(
            self,
            inside: bool = True,
            last_child_name: str | None = None,
    ) -> str | None:
        """
        Removes irrelevant individuals from the pedigree based on their haplotype class and relationships.
        If all descendants of an excluded individual are also excluded, the individual and its descendants are removed.
        If a suspect individual is removed, the closest known ancestor is returned.

        Args:
            inside (bool): If True, it is assumed to be the known pedigree.
                           If False, it is assumed to be the extended pedigree
            last_child_name (str | None): The name of the last child individual, if applicable.
        Returns:
            str | None: The name of the closest known ancestor if a suspect is removed, otherwise the name of the suspect.
        """
        # remove excluded individuals that have no non-excluded children
        individuals_to_remove = []
        unknown_individuals = self.get_unknown_individuals()
        G = create_nx_graph(self)

        if inside:
            for individual in sorted(self.individuals, key=lambda x: str(x.id)):
                descendants = nx.descendants(G, individual.id)
                if individual.exclude:
                    # if all descendants are excluded remove all individuals in the subtree
                    if all(self.get_individual_by_id(descendant_id).exclude for descendant_id in descendants):
                        individuals_to_remove.append(str(individual.id))
                        for descendant_id in descendants:
                            individuals_to_remove.append(str(descendant_id))
                elif individual.haplotype_class == "known" or individual.haplotype_class == "suspect":
                    irrelevant_in_path = True
                    for unknown_ind in unknown_individuals:
                        shortest_path = nx.shortest_path(G.to_undirected(), source=individual.id, target=unknown_ind.id)
                        # check if any of the individuals in the path has haplotype_class "known"
                        if not any(self.get_individual_by_id(node_id).haplotype_class == "known" for node_id in shortest_path[1:-1]):
                            irrelevant_in_path = False
                    if irrelevant_in_path:
                        individuals_to_remove.append(str(individual.id))

        else:
            nodes_to_keep = []
            known_individuals = self.get_known_individuals()
            suspect = self.get_suspect()
            if suspect:
                known_individuals.append(suspect)
            last_child_id = self.get_individual_by_name(last_child_name).id if last_child_name else None
            for known_individual in known_individuals:
                shortest_path = nx.shortest_path(G.to_undirected(), source=last_child_id, target=known_individual.id)
                if all(self.get_individual_by_id(node_id).haplotype_class == "unknown" for node_id in shortest_path[1:-1]):
                    # if the path contains no other known individuals, keep the path
                    for node_id in shortest_path:
                        individual = self.get_individual_by_id(node_id)
                        if str(individual.id) not in nodes_to_keep:
                            nodes_to_keep.append(str(individual.id))

            # remove all individuals that are not in the nodes_to_keep list
            for individual in sorted(self.individuals, key=lambda x: str(x.id)):
                if str(individual.id) not in nodes_to_keep:
                    individuals_to_remove.append(str(individual.id))

        suspect = self.get_suspect()
        closest_known_ancestor = None
        if suspect and suspect.id in individuals_to_remove:
            for ancestor_id in sorted(nx.ancestors(G, suspect.id), key=lambda x: str(x)):
                ancestor = self.get_individual_by_id(ancestor_id)
                if ancestor.haplotype_class == "known":
                    closest_known_ancestor = ancestor
                    break

        individuals_to_remove = list(set(individuals_to_remove))

        for individual_id in individuals_to_remove:
            self.remove_individual(individual_id)

        if closest_known_ancestor:
            return closest_known_ancestor.name
        else:
            return suspect.name


@dataclass(frozen=True)
class IterationResult:
    probability: Decimal
    edge_probabilities: Mapping[tuple[str, str], Decimal]
    mutated_haplotypes: Mapping[str, Haplotype]
    fixed_individual_id: str | int | None


@dataclass(frozen=False)
class SimulationParameters:
    max_number_of_iterations: int
    two_step_mutation_factor: float
    stability_window: int
    stability_min_iterations: int
    stability_threshold: float
    model_validity_threshold: float
    simulation_name: str
    number_of_threads: int
    results_path: Path
    random_seed: int | None = None
    user_name: str | None = None


@dataclass(frozen=True)
class SimulationResult:
    pedigree: Pedigree
    marker_set: MarkerSet
    root_name: str
    simulation_parameters: SimulationParameters
    random_seed: int

    average_pedigree_probability: Decimal
    extended_average_pedigree_probability: Decimal
    inside_match_probability: Mapping[int, Decimal]
    outside_match_probability: Decimal

    average_pedigree_needed_iterations: list[int]
    extended_needed_iterations: list[int]
    inside_needed_iterations: list[int]
    outside_needed_iterations: list[int]

    average_pedigree_model_pedigree_probabilities: list[Decimal]
    extended_model_pedigree_probabilities: list[Decimal]
    inside_model_probabilities: list[Decimal]
    outside_model_probabilities: list[Decimal]

    average_pedigree_model_is_valid: bool
    extended_model_is_valid: bool
    inside_model_is_valid: bool
    outside_model_is_valid: bool

    average_used_probabilities: list[Decimal]
    extended_used_probabilities: list[Decimal]
    inside_used_probabilities: list[Decimal]
    outside_used_probabilities: list[Decimal]

    run_time_pedigree_probability: timedelta
    run_time_proposal_distribution: timedelta
    run_time_extended_average_pedigree_probability: timedelta
    run_time_outside_match_probability: timedelta
    total_run_time: timedelta


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
    return 1 - ((1 - mutation_rate) ** (1 / number_of_copies))


def generate_unique_matchings(
        parent_alleles: list[Allele],
        child_alleles: list[Allele],
        is_average_pedigree: bool,
        )\
        -> Generator[list[tuple[Allele, Allele]]]:
    """
    Generate unique parent→child allele matchings
    without overcounting duplicates.
    """
    n = len(parent_alleles)
    seen = set()

    # for child_perm in permutations(child_alleles, n):
    #     yield list(zip(parent_alleles, child_perm))

    for child_perm in permutations(child_alleles, n):
        if child_perm in seen:
            continue
        seen.add(child_perm)
        yield list(zip(parent_alleles, child_perm))
    #
    # #
    # # if is_average_pedigree:
    # # for child_perm in permutations(child_alleles, n):
    # #     yield list(zip(parent_alleles, child_perm))
    # # #
    # # # else:
    # # #     for child_perm in permutations(child_alleles, n):
    # # #         if child_perm in seen:
    # # #             continue
    # # #         seen.add(child_perm)
    # # #         yield list(zip(parent_alleles, child_perm))
    #
    # n = len(parent_alleles)
    #
    # # if is_average_pedigree:
    # #     # keep all permutations
    # #     for child_perm in permutations(child_alleles, n):
    # #         yield list(zip(parent_alleles, child_perm))
    # # else:
    # seen = set()
    # for child_perm in permutations(child_alleles, n):
    #     for parent_pern in permutations(parent_alleles, n):
    #         yield list(zip(parent_pern, child_perm))
    #     combo = list(zip(parent_alleles, child_perm))
    #
    #     normalized = []
    #     for parent_value in set(parent_alleles):
    #         assigned_children = [c for p, c in combo if p == parent_value]
    #
    #         if len(set(assigned_children)) == 1:
    #             # all children same → collapse to one representative
    #             normalized.append((parent_value, (assigned_children[0],) * len(assigned_children)))
    #         else:
    #             # children differ → keep order (no collapsing)
    #             normalized.append((parent_value, tuple(assigned_children)))
    #
    #     normalized = tuple(sorted(normalized, key=lambda x: hash(x[0])))
    #
    #     if normalized not in seen:
    #         seen.add(normalized)
    #         yield combo


def calculate_mutation_probability(
        parent_alleles: list[Allele],
        child_alleles: list[Allele],
        marker: Marker,
        two_step_mutation_factor: float,
        is_average_pedigree: bool,
) -> Decimal:
    mutation_probability = Decimal(0)

    unique_matchings = list(generate_unique_matchings(parent_alleles, child_alleles, is_average_pedigree))

    for i, combination in enumerate(unique_matchings):
        combination_probability = Decimal(1)
        for parent_allele, child_allele in combination:
            if child_allele.intermediate_value != parent_allele.intermediate_value:
                combination_probability = Decimal(0)  # No mutation for intermediate values mismatch
                break
            else:
                mutation_value = child_allele.value - parent_allele.value
                parent_child_probability = Decimal(
                    get_mutation_probability(
                        get_single_copy_mutation_rate(marker.mutation_rate,
                                                      marker.number_of_copies),
                        mutation_value,
                        two_step_mutation_factor
                    ))
                combination_probability *= parent_child_probability
        mutation_probability += combination_probability

    return mutation_probability
