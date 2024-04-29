from models import Simulation, Pedigree, MarkerSet, Individual, get_edge_probability
import random
import streamlit as st
from copy import deepcopy
from stqdm import stqdm


def mutate_haplotype(source, target, marker_set):
    target.haplotype_class = "simulated"
    for marker in marker_set.markers:
        target.mutate_allele(marker, source.get_allele_by_marker_name(marker.name))


def calculate_pedigree_probability(pedigree, marker_set, suspect, number_of_iterations, simulation):
    """
    This function calculates the pedigree probability after a number of iterations. This corresponds to Pr(Hkn=hkn).
    This probability is only calculated once for a given pedigree. The result is stored in the simulation object.
    The calculation is performed on the rearranged pedigree, where the suspect is the root node.
    During each iteration, the haplotypes of the unknown individuals are predicted based on the haplotypes
    of the surrounding known individuals, and the given mutation rates.

    The order in which the unknown individuals are processed is determined by the level order traversal
    of the rearranged pedigree. This means that the (unknown) children of the suspect are processed first,
    then the (unknown) children of the children, and so on. For each unknown individual, the haplotypes of
    the surrounding known individuals are used to predict the haplotypes of the unknown individual.
    If an unknown individual is surrounded by multiple known individuals, one of them is randomly selected.

    After all unknown individuals have been processed, the pedigree probability is calculated. This is done by
    multiplying the probabilities of the individual edges in the rearranged pedigree that were not used to predict
    the haplotypes of the unknown individuals. An edge is the connection between a parent and a child in the pedigree.
    The probability of an edge is the probability of the child's haplotype given the parent's haplotype.
    The probability of an edge is calculated by multiplying the probabilities of the child's alleles given the parent's
    alleles. The probability of an allele is calculated by the mutation rate of the allele's marker.
    """

    average_pedigree_probability = 0

    for i in stqdm(range(number_of_iterations)):
        pedigree_deep_copy = deepcopy(pedigree)

        level_order_traversal = pedigree_deep_copy.get_level_order_traversal(suspect)
        for individual in level_order_traversal:
            if individual.haplotype_class == "unknown":
                surrounding_known_individuals = pedigree_deep_copy.get_surrounding_known_individuals(individual)
                random_known_node = random.choice(surrounding_known_individuals)
                mutate_haplotype(source=random_known_node, target=individual, marker_set=marker_set)
                pedigree_deep_copy.set_relationship_class(random_known_node, individual, "simulated")

        pedigree_deep_copy.calculate_allele_probabilities(marker_set)
        pedigree_probability = pedigree_deep_copy.get_pedigree_probability(marker_set, "unused")
        # st.write(pedigree_deep_copy.print_pedigree())
        average_pedigree_probability = ((average_pedigree_probability * i) + pedigree_probability) / (i + 1)

        del pedigree_deep_copy

    simulation.average_pedigree_probability = average_pedigree_probability
    return average_pedigree_probability


def calculate_l_matching_haplotypes(pedigree: Pedigree,
                                    marker_set: MarkerSet,
                                    suspect_name: str,
                                    number_of_iterations: int,
                                    simulation: Simulation,
                                    l: int,
                                    ):
    """
    This function simulates the pedigree and counts the number of times the simulation results in l matching haplotypes.
    First l number of unknown individuals are selected, and their haplotypes are set to the same haplotype
    as the suspect. Then the haplotypes of the remaining unknown individuals are predicted based on the
    haplotypes of the surrounding known individuals, and the given mutation rates. Each time a random edge with
    one unknown individual and one known individual is selected. The haplotype of the unknown individual is predicted
    based on the haplotype of the known individual.

    This calculation keeps track of two things:
    - The simulation probability. This corresponds to (Pr(H̃uk = huk|Hkn = hkn))
    - The pedigree probability.

    The simulation probability is the probability of this specific iteration of the simulation. It is calculated by
    the following steps:
    - Set the haplotypes of l unknown individuals to the same haplotype as the suspect.
    This results in a probability of l/total number of unknown individuals.
    """

    suspect = pedigree.get_individual_by_name(suspect_name)

    st.write(f"Calculating l={l} matching haplotypes")

    total_l_probability = 0
    total_l_probability_unweighted = 0

    for i in stqdm(range(number_of_iterations)):
        simulation_probability = 1

        pedigree_deep_copy = deepcopy(pedigree)
        unknown_individuals = pedigree_deep_copy.get_unknown_individuals()

        random_unknown_individuals = random.sample(unknown_individuals, l)

        simulation_probability *= l / len(unknown_individuals)

        for individual in random_unknown_individuals:
            individual.haplotype_class = "fixed"
            for marker in marker_set.markers:
                individual.add_allele(marker, suspect.get_allele_by_marker_name(marker.name).value)

        while len(pedigree_deep_copy.get_unknown_individuals()) > 0:
            edges_with_one_unknown_and_one_known_individual = list(pedigree_deep_copy.get_edges_with_one_unknown_and_one_known_individual())
            # simulation_probability *= (1 / len(edges_with_one_unknown_and_one_known_individual))
            random_unknown_individual, random_known_individual = random.choice(edges_with_one_unknown_and_one_known_individual)
            mutate_haplotype(source=random_known_individual, target=random_unknown_individual, marker_set=marker_set)
            pedigree_deep_copy.set_relationship_class(random_known_individual, random_unknown_individual, "simulated")
            edge_probability = get_edge_probability(marker_set, random_known_individual, random_unknown_individual)
            simulation_probability *= edge_probability

        # Calculate the probability of the entire pedigree (Pr(H = h))
        pedigree_deep_copy.calculate_allele_probabilities(marker_set)
        pedigree_probability = pedigree_deep_copy.get_pedigree_probability(marker_set, "all")

        # Calculate (Pr(Huk = huk|Hkn = hkn))
        conditional_probability = pedigree_probability / simulation.average_pedigree_probability

        # Check if the total number of matching haplotypes is equal to l (excluding the suspect)
        number_of_matching_haplotypes = 0
        for individual in pedigree_deep_copy.individuals:
            if individual.name != suspect_name:
                if individual.has_same_haplotype_as(suspect):
                    number_of_matching_haplotypes += 1

        if number_of_matching_haplotypes == l:
            total_l_probability += (conditional_probability / simulation_probability)
            total_l_probability_unweighted += 1

        del pedigree_deep_copy

    simulation.l_matching_haplotypes_probability[l] = (total_l_probability / number_of_iterations)
    st.write(f"Probability of {l} matching haplotypes: {simulation.l_matching_haplotypes_probability[l]}")
    st.write(f"Unweighted probability of {l} matching haplotypes: {total_l_probability_unweighted / number_of_iterations}")


def calculate_proposal_distribution(pedigree: Pedigree,
                                    marker_set: MarkerSet,
                                    suspect: str,
                                    number_of_iterations: int,
                                    simulation: Simulation,
                                    ):
    unknown_individuals = pedigree.get_unknown_individuals()
    number_of_unknown_individuals = len(unknown_individuals)

    for l in range(1, number_of_unknown_individuals + 1):  # l is the number of matching haplotypes
        calculate_l_matching_haplotypes(pedigree, marker_set, suspect, number_of_iterations, simulation, l)


def run_simulation(pedigree: Pedigree, marker_set: MarkerSet, suspect: str, number_of_iterations: int):
    """
    The total simulation is based on a Monte-Carlo simulation with Importance Sampling.
    This means that the simulation uses a new distribution where the probability of observing the suspect's
    haplotype is high. This simulation uses a proposal distribution where the number of matching haplotypes
    is set to a fixed number.
    """

    simulation = Simulation()

    """Step 1: calculate pedigree probability. This needs to be done only once for the original pedigree
    This corresponds to Pr(Hkn=hkn)"""
    average_pedigree_probability = calculate_pedigree_probability(pedigree, marker_set, suspect, number_of_iterations, simulation)
    st.write(f"Average pedigree probability Pr(Hkn=hkn) after {number_of_iterations} iterations: {average_pedigree_probability}")

    """Step 2: calculate the number of matching haplotypes with the suspect.
    This is done for the total number of individuals in the pedigree.
    This corresponds to Pr(Huk=hi, Hkn=hkn)"""
    proposal_distribution = calculate_proposal_distribution(pedigree, marker_set, suspect, number_of_iterations, simulation)
