from dataclasses import dataclass, field
from datetime import timedelta
from decimal import Decimal
from pathlib import Path
from typing import Collection, Iterator, Mapping

import networkx as nx


@dataclass(frozen=True)
class Marker:
    name: str
    mutation_rate: float


@dataclass
class Allele:
    marker: Marker
    value: int
    parent_value: int | None = None
    mutation_value: int | None = None
    mutation_probability: float | None = None


@dataclass
class Haplotype:
    def __init__(self):
        self.alleles: dict[str, Allele] = {}


@dataclass(frozen=False)
class Individual:
    id: int
    name: str
    haplotype: Haplotype = field(default_factory=lambda: Haplotype())
    haplotype_class: str = "unknown"

    def add_allele(self, marker: Marker, value: int):
        self.haplotype.alleles[marker.name] = Allele(marker, value)

    def get_allele_by_marker_name(self, marker_name: str) -> Allele | None:
        return self.haplotype.alleles.get(marker_name)

    def has_same_haplotype_as(self, other_individual):
        if len(self.haplotype.alleles) != len(other_individual.haplotype.alleles):
            return False
        for allele in self.haplotype.alleles.values():
            other_allele = other_individual.get_allele_by_marker_name(
                allele.marker.name
            )
            if other_allele and allele.value != other_allele.value:
                return False
        return True


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

    def read_marker_set_from_file(self, path: Path):
        with path.open() as file:
            # Skip header
            next(file)
            for line in file:
                marker_name, mutation_rate = line.split(",")
                self.add_marker(Marker(str(marker_name), float(mutation_rate)))


class Pedigree:
    def __init__(self):
        self.graph = nx.DiGraph()
        self.individuals = []
        self.relationships = []

    def print(self):
        for individual in self.individuals:
            print(f"{individual.name},{individual.haplotype_class}")
            for allele in individual.haplotype.alleles:
                print(f"\t{allele}")

    def add_individual(self, individual_id: int, name: str):
        individual = Individual(individual_id, name)
        self.individuals.append(individual)
        self.graph.add_node(individual_id)

    def add_relationship(self, parent_id: int, child_id: int):
        relationship = Relationship(parent_id, child_id)
        self.relationships.append(relationship)
        self.graph.add_edge(parent_id, child_id)

    def read_pedigree_from_file(self, path: Path):
        with path.open() as file:
            current_section = "node"

            for line in file:
                # Remove leading and trailing whitespaces
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                # Switch to edge section when encountering '#'
                if line == "#":
                    current_section = "edge"
                    continue

                # Process lines based on the current section
                if current_section == "node":
                    # Split the line into id and name
                    individual_id, individual_name = line.split()
                    self.add_individual(int(individual_id), str(individual_name))
                elif current_section == "edge":
                    # Split the line into id1 and id2
                    parent_id, child_id = line.split()
                    self.add_relationship(int(parent_id), int(child_id))

    def get_edges(self) -> Iterator[tuple[Individual, Individual]]:
        for parent_id, child_id in self.graph.edges():
            parent = self.get_individual_by_id(parent_id)
            child = self.get_individual_by_id(child_id)
            yield parent, child

    def get_known_individuals_names(self) -> list[str]:
        return [
            individual.name
            for individual in self.individuals
            if individual.haplotype_class == "known"
        ]

    def get_surrounding_known_individuals(
        self, individual: Individual
    ) -> list[Individual]:
        surrounding_known_individuals = []
        for parent_id, child_id in self.graph.edges():
            if parent_id == individual.id:
                child = self.get_individual_by_id(child_id)
                if child.haplotype_class != "unknown":
                    surrounding_known_individuals.append(child)
            elif child_id == individual.id:
                parent = self.get_individual_by_id(parent_id)
                if parent.haplotype_class != "unknown":
                    surrounding_known_individuals.append(parent)
        return surrounding_known_individuals

    def read_known_haplotype_from_file(
        self, individual_name: str, path: Path, marker_set: MarkerSet
    ):
        individual = self.get_individual_by_name(individual_name)
        individual.haplotype_class = "known"
        with path.open() as file:
            # Skip header
            next(file)
            for line in file:
                marker_name, value = line.split(",")
                marker = marker_set.get_marker_by_name(marker_name)
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
        new_root = self.get_individual_by_name(new_root_name)
        new_root.haplotype_class = "suspect"
        self.graph = nx.DiGraph(
            nx.dfs_tree(self.graph.to_undirected(), source=new_root.id)
        )
        self.relationships = [
            Relationship(parent_id, child_id)
            for parent_id, child_id in self.graph.edges()
        ]

    def get_level_order_traversal(self, source: str) -> list[Individual]:
        source = self.get_individual_by_name(source)
        ordered_node_ids = []
        for level in nx.bfs_layers(self.graph, sources=source.id):
            for individual_id in level:
                individual = self.get_individual_by_id(individual_id)
                ordered_node_ids.append(individual)
        return ordered_node_ids

    def get_edges_by_node_id(self, node_id) -> list[tuple[int, int]]:
        edges = []
        for parent_id, child_id in self.graph.edges():
            if parent_id == node_id or child_id == node_id:
                edges.append((parent_id, child_id))
        return edges

    def get_unused_edges(self) -> list[Relationship]:
        return [
            relationship
            for relationship in self.relationships
            if relationship.edge_class == "unknown"
        ]

    def get_simulated_edges(self) -> list[Relationship]:
        return [
            relationship
            for relationship in self.relationships
            if relationship.edge_class == "simulated"
        ]

    def get_all_edges(self) -> list[Relationship]:
        return [relationship for relationship in self.relationships]

    def get_edges_with_one_unknown_and_one_known_individual(self) -> list[Relationship]:
        for relationship in self.relationships:
            parent = self.get_individual_by_id(relationship.parent_id)
            child = self.get_individual_by_id(relationship.child_id)
            if (
                parent.haplotype_class != "unknown"
                and child.haplotype_class == "unknown"
            ):
                yield child, parent
            elif (
                parent.haplotype_class == "unknown"
                and child.haplotype_class != "unknown"
            ):
                yield parent, child

    def set_relationship_class(
        self, parent: Individual, child: Individual, relationship_class: str
    ):
        for relationship in self.relationships:
            if (
                relationship.parent_id == parent.id
                and relationship.child_id == child.id
            ):
                relationship.edge_class = relationship_class

    def get_parent(self, individual):
        for parent_id, child_id in self.graph.edges():
            if child_id == individual.id:
                return self.get_individual_by_id(parent_id)

    def get_pedigree_probability(self, marker_set: MarkerSet, edge_type: str) -> float:
        pedigree_probability = 1

        if edge_type == "all":
            edges = self.get_all_edges()
        elif edge_type == "unused":
            edges = self.get_unused_edges()
        elif edge_type == "simulated":
            edges = self.get_simulated_edges()
        else:
            raise ValueError(f"Invalid edge type: {edge_type}")

        for relationship in edges:
            child = self.get_individual_by_id(relationship.child_id)
            parent = self.get_individual_by_id(relationship.parent_id)
            edge_probability = get_edge_probability(parent, child, marker_set)
            pedigree_probability *= edge_probability

        return pedigree_probability

    def calculate_allele_probabilities(self, marker_set: MarkerSet):
        for parent_id, child_id in self.graph.edges():
            parent = self.get_individual_by_id(parent_id)
            child = self.get_individual_by_id(child_id)
            for marker in marker_set.markers:
                parent_allele = parent.get_allele_by_marker_name(marker.name)
                child_allele = child.get_allele_by_marker_name(marker.name)

                child_allele.parent_value = parent_allele.value
                child_allele.mutation_value = child_allele.value - parent_allele.value
                child_allele.mutation_probability = get_mutation_probability(
                    marker.mutation_rate, child_allele.mutation_value
                )


@dataclass(frozen=True)
class IterationResult:
    pedigree_probability: float
    matching_haplotype_count: int


@dataclass(frozen=True)
class SimulationResult:
    average_pedigree_probability: Decimal
    proposal_distribution: Mapping[int, Decimal]
    run_time_pedigree_probability: timedelta
    run_time_proposal_distribution: timedelta


def get_mutation_probability(mutation_rate: float, mutation_value: float) -> float:
    if mutation_value == 0:
        return 1 - mutation_rate
    elif mutation_value == 1 or mutation_value == -1:
        return (mutation_rate * 0.97) / 2
    elif mutation_value == 2 or mutation_value == -2:
        return (mutation_rate * 0.03) / 2
    else:
        return 0


def get_edge_probability(
    known: Individual, unknown: Individual, marker_set: MarkerSet
) -> float:
    edge_probability = 1

    for marker in marker_set.markers:
        known_allele = known.get_allele_by_marker_name(marker.name).value
        unknown_allele = unknown.get_allele_by_marker_name(marker.name).value

        mutation_value = unknown_allele - known_allele
        mutation_probability = get_mutation_probability(
            marker.mutation_rate, mutation_value
        )

        edge_probability *= mutation_probability

    return edge_probability
