from models import Simulation, Pedigree, MarkerSet, Individual
import random
import streamlit as st
from copy import deepcopy


def mutate_haplotype(source, target, marker_set):
    target.haplotype_class = "simulated"
    for marker in marker_set.markers:
        target.mutate_allele(marker, source.get_allele_by_marker_name(marker.name))


def calculate_pedigree_probability(pedigree, marker_set, suspect, number_of_iterations, simulation):
    average_pedigree_probability = 0

    for i in range(number_of_iterations):
        pedigree_deep_copy = deepcopy(pedigree)

        level_order_traversal = pedigree_deep_copy.get_level_order_traversal(suspect)
        for individual in level_order_traversal:
            if individual.haplotype_class == "unknown":
                surrounding_known_individuals = pedigree_deep_copy.get_surrounding_known_individuals(individual)
                # randomly select one of the nodes
                random_known_node = random.choice(surrounding_known_individuals)
                mutate_haplotype(source=random_known_node, target=individual, marker_set=marker_set)

        pedigree_deep_copy.calculate_allele_probabilities(marker_set)
        pedigree_probability = pedigree_deep_copy.get_pedigree_probability(marker_set)
        average_pedigree_probability = ((average_pedigree_probability * i) + pedigree_probability) / (i + 1)
    st.write(f"Average pedigree probability after {i} iterations: {average_pedigree_probability}")


def run_simulation(pedigree: Pedigree, marker_set: MarkerSet, suspect: str, number_of_iterations: int):
    simulation = Simulation()

    # Step 1: calculate pedigree probability
    calculate_pedigree_probability(pedigree, marker_set, suspect, number_of_iterations, simulation)
