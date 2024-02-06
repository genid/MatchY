import random

import networkx as nx
from streamlit_agraph import agraph, Node, Edge, Config
import streamlit as st

# Constants for color codes
NODE_COLOR_KNOWN_HAPLOTYPE = "#b2d3c2"
NODE_COLOR_UNKNOWN_HAPLOTYPE = "#eeeeee"
NODE_COLOR_SUSPECT = "#ff0000"

EDGE_COLOR = "#aaaaaa"


class Marker:
    def __init__(self, marker_name: str, mutation_rate: float):
        self.name = marker_name
        self.mutation_rate = mutation_rate


class Allele:
    def __init__(self, marker: Marker,
                 value: int):
        self.marker = marker
        self.value = value
        self.parent_value: int = None
        self.mutation_value: int = None
        self.mutation_probability: float = None


class Haplotype:
    def __init__(self):
        self.alleles = []


class Individual:
    def __init__(self, individual_id: int, name: str):
        self.id = individual_id
        self.name = name
        self.haplotype: Haplotype = Haplotype()
        self.haplotype_class: str = "unknown"

    def add_allele(self, marker: Marker,
                   value: int):
        self.haplotype.alleles.append(Allele(marker, value))

    def get_allele_by_marker_name(self, marker_name: str) -> Allele:
        for allele in self.haplotype.alleles:
            if allele.marker.name == marker_name:
                return allele
        return None

    def mutate_allele(self, marker: Marker, source_allele: Allele):
        mutation_rate = marker.mutation_rate
        mutation_step = random.choices([0, 1], weights=[1 - mutation_rate, mutation_rate])[0]
        mutation_direction = random.choice([-1, 1])
        target_allele_value = source_allele.value + (mutation_step * mutation_direction)
        mutation_probability = (mutation_rate / 2) if mutation_step == 1 else 1 - mutation_rate
        self.add_allele(marker, target_allele_value)

    @property
    def node_color(self):
        if self.haplotype_class == "known":
            return NODE_COLOR_KNOWN_HAPLOTYPE
        elif self.haplotype_class == "unknown":
            return NODE_COLOR_UNKNOWN_HAPLOTYPE
        elif self.haplotype_class == "suspect":
            return NODE_COLOR_SUSPECT


class MarkerSet:
    def __init__(self):
        self.markers = []

    def add_marker(self, marker: Marker):
        self.markers.append(marker)

    def get_marker_by_name(self, marker_name: str) -> Marker:
        for marker in self.markers:
            if marker.name == marker_name:
                return marker
        return None

    def read_marker_set_from_file(self, filename: str):
        with open(filename, 'r') as file:
            # Skip header
            next(file)
            for line in file:
                marker_name, mutation_rate = line.split(",")
                self.add_marker(Marker(str(marker_name), float(mutation_rate)))


def get_mutation_probability(mutation_rate, mutation_value):
    if mutation_value == 0:
        return 1 - mutation_rate
    elif mutation_value == 1 or mutation_value == -1:
        return mutation_rate / 2
    else:
        return 0


class Pedigree:
    def __init__(self):
        self.graph = nx.DiGraph()
        self.individuals = []

    def add_individual(self, individual_id: int, name: str):
        individual = Individual(individual_id, name)
        self.individuals.append(individual)
        self.graph.add_node(individual_id)

    def add_relationship(self, parent_id: int, child_id: int):
        self.graph.add_edge(parent_id, child_id)

    def read_pedigree_from_file(self, filename: str):
        with open(filename, 'r') as file:
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

    def get_edges(self) -> tuple[Individual, Individual]:
        for parent_id, child_id in self.graph.edges():
            parent = self.get_individual_by_id(parent_id)
            child = self.get_individual_by_id(child_id)
            yield parent, child

    def get_known_individuals_names(self) -> list[str]:
        return [individual.name for individual in self.individuals if individual.haplotype_class == "known"]

    def get_surrounding_known_individuals(self, individual: Individual) -> list[Individual]:
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

    def read_known_haplotype_from_file(self, individual_name: str,
                                       filename: str,
                                       marker_set: MarkerSet):
        individual = self.get_individual_by_name(individual_name)
        individual.haplotype_class = "known"
        with open(filename, 'r') as file:
            # Skip header
            next(file)
            for line in file:
                marker_name, value = line.split(",")
                marker = marker_set.get_marker_by_name(marker_name)
                individual.add_allele(marker, int(value))

    def get_individual_by_name(self, individual_name: str) -> Individual:
        for individual in self.individuals:
            if individual.name == individual_name:
                return individual
        return None

    def get_individual_by_id(self, individual_id: int) -> Individual:
        for individual in self.individuals:
            if individual.id == individual_id:
                return individual
        return None

    def reroot_pedigree(self, new_root_name: str):
        new_root = self.get_individual_by_name(new_root_name)
        new_root.haplotype_class = "suspect"
        self.graph = nx.DiGraph(nx.dfs_tree(self.graph.to_undirected(), source=new_root.id))

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

    def calculate_allele_probabilities(self, marker_set: MarkerSet):
        for parent_id, child_id in self.graph.edges():
            parent = self.get_individual_by_id(parent_id)
            child = self.get_individual_by_id(child_id)
            for marker in marker_set.markers:
                parent_allele = parent.get_allele_by_marker_name(marker.name)
                child_allele = child.get_allele_by_marker_name(marker.name)

                child_allele.parent_value = parent_allele.value
                child_allele.mutation_value = child_allele.value - parent_allele.value
                child_allele.mutation_probability = get_mutation_probability(marker.mutation_rate,
                                                                             child_allele.mutation_value)

    def get_pedigree_probability(self, marker_set: MarkerSet) -> float:
        pedigree_probability = 1
        for parent_id, child_id in self.graph.edges():
            parent = self.get_individual_by_id(parent_id)
            child = self.get_individual_by_id(child_id)
            for marker in marker_set.markers:
                child_allele = child.get_allele_by_marker_name(marker.name)
                pedigree_probability *= child_allele.mutation_probability
        return pedigree_probability

    def visualize_pedigree(self) -> int:
        config = Config(
            width=1200,
            height=700,
            directed=True,
            hierarchical=True,
            direction="UD",
            sortMethod="directed",
            physics=False,
            nodeSpacing=150,
        )

        nodes = [Node(id=individual.id,
                      label=individual.name,
                      color=individual.node_color)
                 for individual in self.individuals]

        edges = [Edge(source=parent_id,
                      target=child_id,
                      color=EDGE_COLOR)
                 for (parent_id, child_id) in self.graph.edges()]

        selected_node_id = agraph(nodes=nodes, edges=edges, config=config)
        return selected_node_id


class Simulation:
    def __init__(self):
        self.simulation_results = []
