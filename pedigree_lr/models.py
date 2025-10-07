from __future__ import annotations
from copy import deepcopy
from dataclasses import dataclass, field
from datetime import timedelta
from decimal import Decimal
from functools import lru_cache
from io import StringIO
from itertools import permutations, combinations
from pathlib import Path
from typing import Mapping, Generator, List
import networkx as nx
import logging
import pandas as pd
import json
import numpy as np
from scipy.optimize import linear_sum_assignment

logger = logging.getLogger(__name__)

"""
This module contains classes and methods for representing and manipulating pedigrees.
The module provides data structures for markers, alleles, individuals, relationships, and pedigrees.
"""


class InvalidAveragePedigreeProbability(Exception):
    pass


@dataclass
class Bias:
    marker: Marker
    copy_nr: int
    direction: str  # "up" or "down"
    target_mass: float


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
        return (self.marker, self.value, self.intermediate_value) == (other.marker, other.value, other.intermediate_value)

    def __hash__(self):
        return hash((self.marker, self.value, self.intermediate_value))

    def __lt__(self, other):
        if not isinstance(other, Allele):
            return NotImplemented
        return (self.value, self.intermediate_value or 0) < (other.value, other.intermediate_value or 0)


@dataclass
class Haplotype:
    """Represents a collection of alleles grouped by genetic markers."""

    def __init__(self):
        self.alleles: dict[str, list[Allele]] = {}

    def __hash__(self):
        """
        Deterministic hash independent of insertion order and robust to None intermediate values.
        We sort markers by name, sort alleles using Allele.__lt__ (which treats None as 0),
        and normalise intermediate_value=None to 0 in the hashed representation.
        """
        normalised = tuple(
            (
                marker_name,
                tuple(
                    (a.value, (0 if a.intermediate_value is None else a.intermediate_value))
                    for a in sorted(alleles)
                ),
            )
            for marker_name, alleles in sorted(self.alleles.items(), key=lambda kv: kv[0])
        )
        return hash(normalised)

    def add_allele(self, marker: Marker, value: int, intermediate_value: int | None = None):
        if marker.name not in self.alleles:
            self.alleles[marker.name] = []
        self.alleles[marker.name].append(Allele(marker, value, intermediate_value))

    def get_alleles_by_marker_name(self, marker_name: str) -> list[Allele]:
        return sorted(self.alleles.get(marker_name, []), key=lambda x: x.value)

    def __eq__(self, other: 'Haplotype') -> bool:
        if self.alleles.keys() != other.alleles.keys():
            return False
        for marker, alleles in self.alleles.items():
            other_alleles = other.alleles.get(marker, [])
            if len(alleles) != len(other_alleles):
                return False
            # Compare multi-copy alleles one-to-one
            for allele, other_allele in zip(sorted(alleles, key=lambda x: x.value),
                                            sorted(other_alleles, key=lambda x: x.value)):
                if (allele.value != other_allele.value or
                        allele.intermediate_value != other_allele.intermediate_value):
                    return False
        return True

    def allelic_difference(self, other: 'Haplotype') -> int:
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
    id: str
    name: str
    haplotype: Haplotype = field(default_factory=lambda: Haplotype())
    haplotype_class: str = "unknown"
    exclude: bool = False
    picking_probability: Decimal | None = None
    closest_known_individuals: list[Individual] = field(default_factory=lambda: [])
    closest_known_distance: int | None = None

    def add_allele(self, marker: Marker, value: int, intermediate_value: int | None = None):
        self.haplotype.add_allele(marker, value, intermediate_value)

    def get_alleles_by_marker_name(self, marker_name: str) -> list[Allele] | None:
        return self.haplotype.get_alleles_by_marker_name(marker_name)

    def __str__(self):
        return f"Individual(id={self.id}, name={self.name}, haplotype_class={self.haplotype_class})"

    def __repr__(self):
        return f"Individual(id={self.id}, name={self.name}, haplotype_class={self.haplotype_class})"


@dataclass
class Relationship:
    parent_id: str
    child_id: str
    edge_class: str = "unknown"


class MarkerSet:
    def __init__(self):
        self.markers: list[Marker] = []

    def __hash__(self):
        return hash(tuple((m.name, m.mutation_rate, m.number_of_copies) for m in sorted(self.markers, key=lambda m: m.name)))

    def add_marker(self, marker: Marker):
        self.markers.append(marker)

    def get_marker_by_name(self, marker_name: str) -> Marker | None:
        for marker in self.markers:
            if marker.name == marker_name:
                return marker
        return None

    def read_marker_set_from_file(self, file):
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

    def load_markers_from_database(self, markers: list[str]):
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
    individuals: list[Individual] = field(default_factory=lambda: [])
    relationships: list[Relationship] = field(default_factory=lambda: [])
    picking_probabilities: dict[str, Decimal] = field(default_factory=lambda: {})

    def __hash__(self):
        return hash((
            tuple(sorted((str(ind.id), ind.name, ind.haplotype_class) for ind in self.individuals)),
            tuple(sorted((str(rel.parent_id), str(rel.child_id)) for rel in self.relationships))
        ))

    def add_individual(self, individual_id: int | str, name: str):
        if any(individual.name == name for individual in self.individuals):
            logger.warning(f"Individual with name {name} already exists.")
        individual = Individual(individual_id, name)
        self.individuals.append(individual)

    def remove_individual(self, individual_id: str):
        individual = self.get_individual_by_id(individual_id)
        if not individual:
            return
        children = [
            relationship.child_id
            for relationship in self.relationships
            if str(relationship.parent_id) == str(individual_id)
        ]
        for child_id in children:
            self.remove_individual(child_id)
        if individual in self.individuals:
            self.individuals.remove(individual)
        self.relationships = [
            relationship
            for relationship in self.relationships
            if (
                str(relationship.parent_id) != str(individual_id)
                and str(relationship.child_id) != str(individual_id)
            )
        ]

    def add_relationship(self, parent_id: str, child_id: str):
        relationship = Relationship(parent_id, child_id)
        self.relationships.append(relationship)

    def read_tgf(self, file):
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

    def read_ped(self, file):
        relationships: list[tuple[int, int]] = []
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

    def read_pedigree_from_file(self, file: StringIO, file_extension: str):
        if file_extension == ".tgf":
            self.read_tgf(file)
        elif file_extension == ".ped":
            self.read_ped(file)

    def write_to_tgf(self) -> bytes:
        lines = []
        for individual in self.individuals:
            lines.append(f"{individual.id} {individual.name}")
        lines.append("#")
        for relationship in self.relationships:
            lines.append(f"{relationship.parent_id} {relationship.child_id}")
        return "\n".join(lines).encode("utf-8")

    def read_known_haplotypes_from_file(self, file: StringIO, marker_set: 'MarkerSet'):
        json_data = json.load(file)
        for individual_name in json_data.keys():
            individual = self.get_individual_by_name(individual_name)
            if not individual:
                logger.warning(f"Individual {individual_name} not found in pedigree. Skipping known haplotype assignment.")
                return
            individual.haplotype_class = "known"
            individual_alleles = json_data[individual_name]
            for marker_name, values in individual_alleles.items():
                marker = marker_set.get_marker_by_name(marker_name)
                if not marker:
                    logger.warning(f"Marker {marker_name} not found in marker set. Skipping allele assignment.")
                    continue
                alleles = values.split(";")
                number_of_copies = len(alleles)
                if not marker.number_of_copies:
                    marker.number_of_copies = number_of_copies
                elif marker.number_of_copies != number_of_copies:
                    logger.error(f"Number of copies mismatch for marker {marker_name}")
                    continue
                for allele in alleles:
                    if "." in allele:
                        allele_val, intermediate_value = allele.split(".")
                        try:
                            allele_int = int(allele_val)
                            intermediate_int = int(intermediate_value)
                            individual.add_allele(marker, allele_int, intermediate_int)
                        except ValueError:
                            logger.error(f"Invalid allele or intermediate value: {allele_val}.{intermediate_value}")
                    else:
                        try:
                            individual.add_allele(marker, int(allele))
                        except ValueError:
                            logger.error(f"Invalid allele value: {allele}")

    def get_individual_by_name(self, individual_name: str) -> Individual | None:
        for individual in self.individuals:
            if str(individual.name) == str(individual_name):
                return individual
        return None

    def get_individual_by_id(self, individual_id: int | str) -> Individual | None:
        for individual in self.individuals:
            if str(individual.id) == str(individual_id):
                return individual
        return None

    def get_parent_by_child_id(self, child_id: str) -> Individual | None:
        for relationship in self.relationships:
            if str(relationship.child_id) == str(child_id):
                return self.get_individual_by_id(relationship.parent_id)
        return None

    def get_unknown_individuals(self) -> list[Individual]:
        return [
            individual for individual in self.individuals if individual.haplotype_class == "unknown"
        ]

    def get_known_individuals(self) -> list[Individual]:
        return [
            individual for individual in self.individuals if individual.haplotype_class == "known"
        ]

    def set_suspect(self, suspect_name: str):
        previous_suspect = self.get_suspect()
        if previous_suspect:
            previous_suspect.haplotype_class = "known"
        suspect_individual = self.get_individual_by_name(suspect_name)
        if suspect_individual:
            suspect_individual.haplotype_class = "suspect"

    def reroot_pedigree(self, new_root_name: str):
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

    def exclude_individuals(self, excluded_individuals: list[str]):
        for individual in self.individuals:
            individual.exclude = False
        for individual_name in excluded_individuals:
            individual = self.get_individual_by_name(individual_name)
            if individual:
                individual.exclude = True

    def get_level_order_traversal(self, source_name: str) -> list[Individual]:
        source = self.get_individual_by_name(source_name)
        ordered_individuals: list[Individual] = []
        for level in nx.bfs_layers(create_nx_graph(self), sources=source.id):
            for individual_id in level:
                individual = self.get_individual_by_id(individual_id)
                if individual:
                    ordered_individuals.append(individual)
        return ordered_individuals

    def get_closest_known_individuals(self):
        unknown_individuals = self.get_unknown_individuals()
        known_individuals = self.get_known_individuals() + [self.get_suspect()]
        undirected_graph = create_nx_graph(self).to_undirected()
        for unknown_individual in unknown_individuals:
            unknown_known_distance_dict: dict[tuple[str, str, str], int] = {}
            for known_individual in known_individuals:
                shortest_path = nx.shortest_path(undirected_graph, source=unknown_individual.id, target=known_individual.id)
                unknown_known_distance_dict[(unknown_individual.id, shortest_path[1], known_individual.id)] = len(shortest_path) - 1
            lowest_distance = min(unknown_known_distance_dict.values())
            closest_known_individual_ids = [
                (unknown_id, via_id)
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

    def calculate_allele_probabilities(self, marker_set: MarkerSet, two_step_mutation_factor: float, is_average_pedigree: bool = False):
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

    def to_string(self):
        lines = ["Pedigree"]
        nx_graph = create_nx_graph(self)
        root = self.get_individual_by_id(1)
        node_depths = {
            node: depth
            for depth, nodes in enumerate(nx.bfs_layers(nx_graph, sources=root.id))
            for node in nodes
        }
        for node_id in nx.dfs_preorder_nodes(nx_graph):
            prefix = "	" * node_depths[node_id]
            individual = self.get_individual_by_id(node_id)
            lines.append(f"{prefix}{individual.name}, {individual.haplotype_class}")
            # List of allele lists (per marker) for readability
            for alleles in individual.haplotype.alleles.values():
                lines.append(f"{prefix}." + ",".join(str(a) for a in alleles))
        return "\n".join(lines)

    def get_suspect(self):
        for individual in self.individuals:
            if individual.haplotype_class == "suspect":
                return individual
        return None

    def check_known_haplotypes(self):
        known_haplotypes = [
            individual for individual in self.individuals if individual.haplotype_class == "known"
        ]
        if not known_haplotypes:
            return
        marker_names = known_haplotypes[0].haplotype.alleles.keys()
        for individual in known_haplotypes[1:]:
            if individual.haplotype.alleles.keys() != marker_names:
                logger.error("Known haplotypes have different markers")
                return

    def check_pedigree_structure(self) -> bool:
        graph = create_nx_graph(self)
        if not nx.is_directed_acyclic_graph(graph):
            logger.error("Pedigree contains cycles")
            return False
        if not nx.is_connected(graph.to_undirected()):
            logger.error("Pedigree contains disconnected components")
            return False
        if not nx.is_tree(graph):
            logger.error("Pedigree is not a tree")
            return False
        for node in graph.nodes:
            if graph.in_degree(node) > 1:
                logger.error("Pedigree contains son(s) with multiple fathers")
                return False
        return True

    def calculate_picking_probabilities(self):
        pedigree_deep_copy = deepcopy(self)
        suspect_haplotype = pedigree_deep_copy.get_suspect().haplotype
        if not suspect_haplotype:
            logger.error("Suspect does not exist or does not have a known haplotype.")
            return
        p = nx.shortest_path(create_nx_graph(pedigree_deep_copy), source=pedigree_deep_copy.get_suspect().id)
        unknown_individual_ids = [i.id for i in pedigree_deep_copy.get_unknown_individuals()]
        picking_probabilities: dict[str, Decimal] = {}
        for individual in unknown_individual_ids:
            shortest_path = p[individual]
            for i, node_id in enumerate(shortest_path):
                individual_obj = pedigree_deep_copy.get_individual_by_id(node_id)
                if individual_obj.haplotype_class == "unknown":
                    individual_obj.haplotype = deepcopy(
                        pedigree_deep_copy.get_individual_by_id(shortest_path[i - 1]).haplotype
                    )
                    individual_obj.haplotype_class = "estimated"
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
        for unknown_individual_id in picking_probabilities.keys():
            picking_probabilities[unknown_individual_id] = Decimal(((max_value - picking_probabilities[unknown_individual_id]) + 1) / (max_value + 1))
        self.picking_probabilities = picking_probabilities

    def extend_pedigree(self):
        G = create_nx_graph(self)
        root = [n for n, d in G.in_degree() if d == 0][0]
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
            new_root_id += 1
        self.add_individual(new_root_id, "new_root")
        self.add_relationship(new_root_id, root)
        previous_parent = new_root_id
        for i in range(1, highest_level_with_known_individual + 1):
            new_child = 0
            while str(new_child) in [str(individual.id) for individual in self.individuals]:
                new_child += 1
            self.add_individual(new_child, f"new_child_{i}")
            self.get_individual_by_id(new_child).haplotype_class = "unknown"
            self.add_relationship(previous_parent, new_child)
            previous_parent = new_child
        last_child_name = f"new_child_{highest_level_with_known_individual}"
        return last_child_name

    def remove_irrelevant_individuals(self, inside: bool = True, last_child_name: str | None = None) -> str | None:
        individuals_to_remove: list[str] = []
        unknown_individuals = self.get_unknown_individuals()
        G = create_nx_graph(self)
        if inside:
            for individual in sorted(self.individuals, key=lambda x: str(x.id)):
                descendants = nx.descendants(G, individual.id)
                if individual.exclude:
                    if all(self.get_individual_by_id(descendant_id).exclude for descendant_id in descendants):
                        individuals_to_remove.append(str(individual.id))
                        for descendant_id in descendants:
                            individuals_to_remove.append(str(descendant_id))
                elif individual.haplotype_class in ("known", "suspect"):
                    irrelevant_in_path = True
                    for unknown_ind in unknown_individuals:
                        shortest_path = nx.shortest_path(G.to_undirected(), source=individual.id, target=unknown_ind.id)
                        # check if any of the individuals in the path has haplotype_class "known"
                        if not any(self.get_individual_by_id(node_id).haplotype_class == "known" for node_id in shortest_path[1:-1]):
                            irrelevant_in_path = False
                    if irrelevant_in_path:
                        individuals_to_remove.append(str(individual.id))
        else:
            nodes_to_keep: list[str] = []
            known_individuals = self.get_known_individuals()
            suspect = self.get_suspect()
            if suspect:
                known_individuals.append(suspect)
            last_child_id = self.get_individual_by_name(last_child_name).id if last_child_name else None
            for known_individual in known_individuals:
                shortest_path = nx.shortest_path(G.to_undirected(), source=last_child_id, target=known_individual.id)
                if all(self.get_individual_by_id(node_id).haplotype_class == "unknown" for node_id in shortest_path[1:-1]):
                    for node_id in shortest_path:
                        individual = self.get_individual_by_id(node_id)
                        if str(individual.id) not in nodes_to_keep:
                            nodes_to_keep.append(str(individual.id))
            for individual in sorted(self.individuals, key=lambda x: str(x.id)):
                if str(individual.id) not in nodes_to_keep:
                    individuals_to_remove.append(str(individual.id))
        suspect = self.get_suspect()
        closest_known_ancestor = None
        if suspect and str(suspect.id) in individuals_to_remove:
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
            return suspect.name if suspect else None

    def get_known_descendants(self, individual_id: str) -> list[Individual]:
        G = create_nx_graph(self)
        if individual_id not in G:
            raise ValueError(f"Individual ID {individual_id} not found in pedigree.")
        return [
            descendant
            for descendant_id in nx.descendants(G, individual_id)
            if (descendant := self.get_individual_by_id(descendant_id)).haplotype_class != "unknown"
        ]

    def get_mrca(self, individual_ids: list[str]) -> Individual | None:
        G = create_nx_graph(self)
        source = [n for n, d in G.in_degree() if d == 0][0]
        if not individual_ids:
            return None
        if any(ind_id not in G for ind_id in individual_ids):
            raise ValueError("One or more individual IDs not found in pedigree.")
        paths = [nx.ancestors(G, ind_id).union({ind_id}) for ind_id in individual_ids]
        common_ancestors = set.intersection(*paths)
        if not common_ancestors:
            return None

        mrca_id = max(common_ancestors, key=lambda n: len(nx.shortest_path(G, source=source, target=n)))

        return self.get_individual_by_id(mrca_id)

    class SimpleDifferenceMatrix:
        """Compute minimal mutation distance using Hungarian algorithm."""
        INTERMEDIATE_MISMATCH_PENALTY: int = 1000

        def __init__(self, allele1: List[Allele], allele2: List[Allele]):
            if not allele1 or not allele2:
                raise ValueError("Cannot compute distance against empty locus")
            if len(allele1) != len(allele2):
                raise ValueError("Alleles must have the same length")
            self.allele1 = sorted(allele1)
            self.allele2 = sorted(allele2)
            self.matrix = self._create_matrix()

        def _allele_distance(self, a1: Allele, a2: Allele) -> int:
            if a1 == a2:
                return 0
            diff = abs(a1.value - a2.value)
            if (a1.intermediate_value or 0) != (a2.intermediate_value or 0):
                diff += self.INTERMEDIATE_MISMATCH_PENALTY
            return diff

        def _create_matrix(self) -> List[List[int]]:
            return [[self._allele_distance(a1, a2) for a2 in self.allele2] for a1 in self.allele1]

        def calculate_mutations(self) -> List[int]:
            np_matrix = np.array(self.matrix)
            rows, cols = linear_sum_assignment(np_matrix)
            assignments = sorted(zip(rows, cols), key=lambda x: x[0])
            signed_mutations: list[int] = []
            for r, c in assignments:
                a1, a2 = self.allele1[r], self.allele2[c]
                if a1 == a2:
                    signed_mutations.append(0)
                    continue
                diff = a2.value - a1.value
                if (a1.intermediate_value or 0) != (a2.intermediate_value or 0):
                    diff = (1 if diff > 0 else -1) * (abs(diff) + 1)
                signed_mutations.append(diff)
            return signed_mutations

        def __str__(self) -> str:
            return "\n".join(map(str, self.matrix))

    @lru_cache(maxsize=128)
    def get_biases(
            self,
            individual_id: str,
            marker_set: MarkerSet,
            haplotypes_tuple: tuple[tuple[str, Haplotype], ...],
            bias_value: float | None = None
    ) -> list[Bias]:
        haplotypes = {hap[0]: hap[1] for hap in haplotypes_tuple}
        parent = self.get_parent_by_child_id(individual_id)
        if parent is None:
            return []
        parent_haplotype = haplotypes.get(parent.id)
        if parent_haplotype is None:
            return []
        known_descendants = self.get_known_descendants(individual_id)
        mrca_of_known_descendants = self.get_mrca([desc.id for desc in known_descendants])
        if mrca_of_known_descendants:
            distance_to_mrca = nx.shortest_path_length(
                create_nx_graph(self).to_undirected(),
                source=individual_id,
                target=mrca_of_known_descendants.id
            )
            if mrca_of_known_descendants.haplotype_class == "unknown":
                distance_to_mrca += 1
            distance_to_mrca = max(distance_to_mrca, 1)
        else:
            distance_to_mrca = 0

        if not known_descendants:
            return []
        biases: list[Bias] = []
        for marker in marker_set.markers:
            mutation_lists = [
                tuple(
                    "up" if mut > 0 else "down" if mut < 0 else "none"
                    for mut in self.SimpleDifferenceMatrix(
                        parent_haplotype.get_alleles_by_marker_name(marker.name),
                        haplotypes[desc.id].get_alleles_by_marker_name(marker.name),
                    ).calculate_mutations()
                )
                for desc in known_descendants
            ]
            mutation_set = set(mutation_lists)
            if len(mutation_set) != 1:
                continue
            mutation_tuple = next(iter(mutation_set))
            for copy_nr, direction in enumerate(mutation_tuple):
                if direction in {"up", "down"}:
                    if bias_value:
                        biases.append(Bias(marker, copy_nr, direction, bias_value))
                    else:
                        bias_value = min(max(0.1, (0.8 / (1 + distance_to_mrca))), 0.4)
                        biases.append(Bias(marker, copy_nr, direction, bias_value))

        # if not bias_value:
        #     new_bias_value = min(0.4, 0.3 * (len(biases) / (1 + distance_to_mrca)))
        #     for bias in biases:
        #         bias.target_mass = new_bias_value

        return biases


@dataclass(frozen=True)
class IterationResult:
    probability: Decimal
    edge_probabilities: Mapping[tuple[str, str], Decimal]
    edge_weight_factor: Decimal
    mutated_haplotypes: Mapping[str, Haplotype]
    fixed_individual_id: str | int | None


@dataclass(frozen=False)
class SimulationParameters:
    two_step_mutation_factor: float
    stability_window: int
    model_validity_threshold: float
    number_of_threads: int
    simulation_name: str
    results_path: Path
    bias: float | None = None
    user_name: str | None = None


@dataclass(frozen=True)
class SimulationResult:
    pedigree: Pedigree
    marker_set: MarkerSet
    root_name: str
    simulation_parameters: SimulationParameters
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
    per_individual_probabilities: Mapping[str, Decimal]
    run_time_pedigree_probability: timedelta
    run_time_proposal_distribution: timedelta
    run_time_extended_average_pedigree_probability: timedelta
    run_time_outside_match_probability: timedelta
    total_run_time: timedelta


def create_nx_graph(pedigree: Pedigree) -> nx.DiGraph:
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
    elif mutation_value in (1, -1):
        return single_step_mutation_rate
    elif mutation_value in (2, -2):
        return two_step_mutation_rate
    else:
        return 0.0


def get_single_copy_mutation_rate(mutation_rate: float, number_of_copies: int) -> float:
    """
    P(at least one mutation) = mu_all - (1 - mu_1)^n
    mu_1 = 1 - (1 - mu_all)^(1/n)
    """
    return 1 - ((1 - mutation_rate ** 1 / number_of_copies))


def generate_unique_matchings(
    parent_alleles: list[Allele],
    child_alleles: list[Allele],
    is_average_pedigree: bool,
) -> Generator[list[tuple[Allele, Allele]], None, None]:
    n = len(parent_alleles)
    seen = set()
    for child_perm in permutations(child_alleles, n):
        if child_perm in seen:
            continue
        seen.add(child_perm)
        yield list(zip(parent_alleles, child_perm))


def calculate_mutation_probability(
    parent_alleles: list[Allele],
    child_alleles: list[Allele],
    marker: Marker,
    two_step_mutation_factor: float,
    is_average_pedigree: bool,
) -> Decimal:
    mutation_probability = Decimal(0)
    unique_matchings = list(generate_unique_matchings(parent_alleles, child_alleles, is_average_pedigree))
    for combination in unique_matchings:
        combination_probability = Decimal(1)
        for parent_allele, child_allele in combination:
            if child_allele.intermediate_value != parent_allele.intermediate_value:
                combination_probability = Decimal(0)
                break
            else:
                mutation_value = child_allele.value - parent_allele.value
                parent_child_probability = Decimal(
                    get_mutation_probability(
                        get_single_copy_mutation_rate(marker.mutation_rate, marker.number_of_copies),
                        mutation_value,
                        two_step_mutation_factor
                    )
                )
                combination_probability *= parent_child_probability
        mutation_probability += combination_probability
    return mutation_probability