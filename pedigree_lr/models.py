from __future__ import annotations

from dataclasses import dataclass, field
from datetime import timedelta
from decimal import Decimal
from io import StringIO
from pathlib import Path
from typing import Mapping

import networkx as nx


@dataclass(frozen=True)
class Marker:
    name: str
    mutation_rate: float


@dataclass
class Allele:
    marker: Marker
    value: int
    intermediate_value: int | None = None
    parent_value: int | None = None
    parent_intermediate_value: int | None = None
    mutation_value: int | None = None
    mutation_probability: float | None = None


@dataclass
class Haplotype:
    def __init__(self):
        self.alleles: dict[str, Allele] = {}

    def __eq__(self, other: "Haplotype") -> bool:
        if len(self.alleles) != len(other.alleles):
            return False
        for allele in self.alleles.values():
            other_allele = other.alleles.get(allele.marker.name)
            if (other_allele and allele.value != other_allele.value or
                    allele.intermediate_value != other_allele.intermediate_value):
                return False
        return True


@dataclass(frozen=False)
class Individual:
    id: int
    name: str
    haplotype: Haplotype = field(default_factory=lambda: Haplotype())
    haplotype_class: str = "unknown"

    def add_allele(self, marker: Marker, value: int, intermediate_value: int = None):
        self.haplotype.alleles[marker.name] = Allele(marker, value, intermediate_value)

    def get_allele_by_marker_name(self, marker_name: str) -> Allele | None:
        return self.haplotype.alleles.get(marker_name)


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
        # Skip header
        next(file)
        for line in file:
            marker_name, mutation_rate = line.split(",")
            self.add_marker(Marker(str(marker_name), float(mutation_rate)))


@dataclass
class Pedigree:
    individuals: list[Individual] = field(default_factory=lambda: [])
    relationships: list[Relationship] = field(default_factory=lambda: [])

    def add_individual(self, individual_id: int, name: str):
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
                individual_id, individual_name = line.split()
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
        individual = self.get_individual_by_name(individual_name)
        individual.haplotype_class = "known"
        next(file) # Skip header
        for line in file:
            marker_name, value = line.split(",")
            marker = marker_set.get_marker_by_name(marker_name)
            if "." in value: # Intermediate allele
                value, intermediate_value = value.split(".")
                individual.add_allele(marker, int(value), int(intermediate_value))
            else:
                individual.add_allele(marker, int(value))

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
        new_root = self.get_individual_by_name(new_root_name)
        new_root.haplotype_class = "suspect"
        current_graph = create_nx_graph(self).to_undirected()
        rerooted_graph = nx.DiGraph(nx.dfs_tree(current_graph, source=new_root.id))
        self.relationships = [
            Relationship(parent_id, child_id)
            for parent_id, child_id in rerooted_graph.edges()
        ]

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
                parent_allele = parent.get_allele_by_marker_name(marker.name)
                child_allele = child.get_allele_by_marker_name(marker.name)

                child_allele.parent_value = parent_allele.value
                child_allele.parent_intermediate_value = parent_allele.intermediate_value

                if child_allele.intermediate_value != parent_allele.intermediate_value:
                    child_allele.mutation_probability = 0 # TODO: Implement intermediate allele mutation probability
                else:
                    child_allele.mutation_value = child_allele.value - parent_allele.value
                    child_allele.mutation_probability = get_mutation_probability(
                        marker.mutation_rate, child_allele.mutation_value
                    )

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
        return 0
