from __future__ import annotations
import sys
from dataclasses import dataclass, field
from datetime import timedelta
from decimal import Decimal
from io import StringIO
from itertools import permutations
from pathlib import Path
from typing import Mapping
import networkx as nx
import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(sys.stdout))

@dataclass
class Marker:
    name: str
    mutation_rate: float
    number_of_copies: int | None = None


@dataclass
class Allele:
    marker: Marker
    value: int
    intermediate_value: int | None = None
    mutation_value: int | None = None
    mutation_probability: float | None = None


@dataclass
class Haplotype:
    def __init__(self):
        self.alleles: dict[str, list[Allele]] = {}

    def add_allele(self, marker: Marker, value: int, intermediate_value: int = None):
        if marker.name not in self.alleles:
            self.alleles[marker.name] = []
        self.alleles[marker.name].append(Allele(marker, value, intermediate_value))

    def get_alleles_by_marker_name(self, marker_name: str) -> list[Allele]:
        return sorted(self.alleles.get(marker_name, []), key=lambda x: x.value)

    def __eq__(self, other: "Haplotype") -> bool:
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
    id: int
    name: str
    haplotype: Haplotype = field(default_factory=lambda: Haplotype())
    haplotype_class: str = "unknown"
    exclude: bool = False

    def add_allele(self, marker: Marker, value: int, intermediate_value: int = None):
        self.haplotype.add_allele(marker, value, intermediate_value)

    def get_alleles_by_marker_name(self, marker_name: str) -> list[Allele] | None:
        return self.haplotype.get_alleles_by_marker_name(marker_name)


@dataclass
class Relationship:
    parent_id: int
    child_id: int
    edge_class: str = "unknown"


class MarkerSet:
    def __init__(self):
        self.markers = []

    def add_marker(self, marker: Marker):
        self.markers.append(marker)

    def get_marker_by_name(self, marker_name: str) -> Marker | None:
        for marker in self.markers:
            if marker.name == marker_name:
                return marker
        return None

    def read_marker_set_from_file(self, file):
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
    individuals: list[Individual] = field(default_factory=lambda: [])
    relationships: list[Relationship] = field(default_factory=lambda: [])

    def add_individual(self, individual_id: int, name: str):
        if any(individual.name == name for individual in self.individuals):
            logger.warning(f"Individual with name {name} already exists.")
        individual = Individual(individual_id, name)
        self.individuals.append(individual)

    def remove_individual(self, individual_id: int):
        individual = self.get_individual_by_id(individual_id)
        if individual:
            self.individuals.remove(individual)
            self.relationships = [
                relationship
                for relationship in self.relationships
                if relationship.parent_id != individual_id
                and relationship.child_id != individual_id
            ]

    def add_relationship(self, parent_id: int, child_id: int):
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

                self.add_individual(int(individual_id), str(individual_name))
            elif current_section == "edge":
                parent_id, child_id = line.split()
                self.add_relationship(int(parent_id), int(child_id))

    def read_ped(self, file):
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

    def read_pedigree_from_file(self, file: StringIO,
                                file_extension: str):
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

    def read_known_haplotype_from_file(
        self, individual_name: str, file, marker_set: MarkerSet
    ):
        try:
            individual = self.get_individual_by_name(individual_name)
        except ValueError:
            logger.error(f"Individual {individual_name} not found in pedigree")
            return
        individual.haplotype_class = "known"
        header = next(file) # Skip header
        if header.strip() != "marker,allele":
            logger.error(f"Invalid header in known haplotype file: {header}")

        for line in file:
            marker_name, values = line.split(",")
            try:
                marker = marker_set.get_marker_by_name(marker_name)
            except ValueError:
                logger.error(f"Marker {marker_name} not found in marker set")
                continue

            alleles = values.split(";") # Use ";" as delimiter for multiple alleles
            number_of_copies = len(alleles)
            if not marker.number_of_copies:
                marker.number_of_copies = number_of_copies
            elif marker.number_of_copies != number_of_copies:
                logger.error(f"Number of copies mismatch for marker {marker_name}")
                continue

            for allele in alleles:
                if "." in allele: # Intermediate allele
                    allele, intermediate_value = allele.split(".")
                    try:
                        allele = int(allele)
                        intermediate_value = int(intermediate_value)
                    except ValueError:
                        logger.error(f"Invalid allele or intermediate value: {allele}.{intermediate_value}")
                    individual.add_allele(marker, allele, intermediate_value)
                else:
                    individual.add_allele(marker, int(allele))

    def get_individual_by_name(self, individual_name: str) -> Individual | None:
        for individual in self.individuals:
            if individual.name == individual_name:
                return individual
        return None

    def get_individual_by_id(self, individual_id: int) -> Individual | None:
        for individual in self.individuals:
            if individual.id == individual_id:
                return individual
        return None

    def get_unknown_individuals(self) -> list[Individual]:
        return [
            individual
            for individual in self.individuals
            if individual.haplotype_class == "unknown"
        ]

    def reroot_pedigree(self, new_root_name: str):
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

    def exclude_individuals(self, excluded_individuals: list[str]):
        for individual_name in excluded_individuals:
            individual = self.get_individual_by_name(individual_name)
            if individual:
                individual.exclude = True

    def get_level_order_traversal(self, source_name: str) -> list[Individual]:
        source = self.get_individual_by_name(source_name)
        ordered_individuals = []
        for level in nx.bfs_layers(create_nx_graph(self), sources=source.id):
            for individual_id in level:
                individual = self.get_individual_by_id(individual_id)
                if individual:
                    ordered_individuals.append(individual)
        return ordered_individuals

    def calculate_allele_probabilities(self, marker_set: MarkerSet):
        for relationship in self.relationships:
            parent = self.get_individual_by_id(relationship.parent_id)
            child = self.get_individual_by_id(relationship.child_id)
            for marker in marker_set.markers:
                parent_alleles = parent.get_alleles_by_marker_name(marker.name)
                child_alleles = child.get_alleles_by_marker_name(marker.name)

                mutation_probability = calculate_mutation_probability(parent_alleles, child_alleles, marker)

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
            prefix = "\t" * node_depths[node_id]
            individual = self.get_individual_by_id(node_id)
            lines.append(f"{prefix}{individual.name}, {individual.haplotype_class}")
            for allele in individual.haplotype.alleles.values():
                lines.append(f"{prefix}.{allele}")

        return "\n".join(lines)

    def get_suspect(self):
        for individual in self.individuals:
            if individual.haplotype_class == "suspect":
                return individual
        return None

    def check_known_haplotypes(self):
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

    def check_pedigree_structure(self):
        graph = create_nx_graph(self)
        if not nx.is_directed_acyclic_graph(graph):
            logger.error("Pedigree contains cycles")
        if not nx.is_connected(graph.to_undirected()):
            logger.error("Pedigree contains disconnected components")
        if not nx.is_tree(graph):
            logger.error("Pedigree is not a tree")
        pass


@dataclass(frozen=True)
class IterationResult:
    pedigree: Pedigree | None
    probability: Decimal
    edge_probabilities: Mapping[tuple[int, int], Decimal]


@dataclass(frozen=True)
class SimulationResult:
    average_pedigree_probability: Decimal
    proposal_distribution: Mapping[int, Decimal]
    run_time_pedigree_probability: timedelta
    run_time_proposal_distribution: timedelta


def create_nx_graph(pedigree: Pedigree) -> nx.DiGraph:
    graph = nx.DiGraph()
    for individual in pedigree.individuals:
        graph.add_node(individual.id)
    for relationship in pedigree.relationships:
        graph.add_edge(relationship.parent_id, relationship.child_id)
    return graph


def get_mutation_probability(mutation_rate: float, mutation_value: float) -> float:
    if mutation_value == 0:
        return 1 - mutation_rate
    elif mutation_value == 1 or mutation_value == -1:
        return (mutation_rate * 0.97) / 2
    elif mutation_value == 2 or mutation_value == -2: # TODO: Remove hard coded value for 2-step mutation
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
                combination_probability = Decimal(0) # No mutation for intermediate values mismatch
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