from copy import deepcopy
from datetime import datetime, timedelta
from decimal import Decimal
from math import comb
from random import Random
from typing import Mapping, Collection

from pedigree_lr.models import (
    Allele,
    Haplotype,
    Individual,
    Marker,
    MarkerSet,
    Pedigree,
    SimulationResult,
    get_edge_probability,
)
from pedigree_lr.reporting import Reporter


def mutate_allele(marker: Marker, source_allele: Allele, random: Random) -> Allele:
    mutation_rate = marker.mutation_rate
    mutation_step = random.choices([0, 1], weights=[1 - mutation_rate, mutation_rate])[0]
    mutation_direction = random.choice([-1, 1])
    mutated_allele_value = source_allele.value + (mutation_step * mutation_direction)
    return Allele(marker, mutated_allele_value)


def mutate_haplotype(source: Individual, marker_set: MarkerSet, random: Random) -> Haplotype:
    haplotype = Haplotype()

    for marker in marker_set.markers:
        source_allele = source.get_allele_by_marker_name(marker.name)
        haplotype.alleles[marker.name] = mutate_allele(marker, source_allele, random)

    return haplotype


def simulate_pedigree_probability(
    pedigree: Pedigree,
    ordered_unknown_ids: Collection[int],
    marker_set: MarkerSet,
    random: Random
) -> Decimal:
    for individual_id in ordered_unknown_ids:
        individual = pedigree.get_individual_by_id(individual_id)
        parent = pedigree.get_parent(individual)

        individual.haplotype = mutate_haplotype(
            source=parent, marker_set=marker_set, random=random
        )
        individual.haplotype_class = "simulated"

        pedigree.set_relationship_class(parent, individual, "simulated")

    pedigree.calculate_allele_probabilities(marker_set)

    pedigree_probability = pedigree.get_pedigree_probability(marker_set, "unused")

    return Decimal(pedigree_probability)


def calculate_pedigree_probability(
    pedigree: Pedigree,
    suspect_name: str,
    marker_set: MarkerSet,
    number_of_iterations: int,
    random: Random,
    reporter: Reporter,
) -> Decimal:
    """
    This function calculates the pedigree probability after a number of iterations. This
    corresponds to Pr(Hkn=hkn). This probability is only calculated once for a given pedigree.
    The result is stored in the simulation object. The calculation is performed on the rearranged
    pedigree, where the suspect is the root node. During each iteration, the haplotypes of the
    unknown individuals are predicted based on the haplotypes of the surrounding known
    individuals, and the given mutation rates.

    The order in which the unknown individuals are processed is determined by the level order
    traversal of the rearranged pedigree. This means that the (unknown) children of the suspect
    are processed first, then the (unknown) children of the children, and so on. For each unknown
    individual, the haplotype of its father (parent node) is used to predict its own haplotype.

    After all unknown individuals have been processed, the pedigree probability is calculated.
    This is done by multiplying the probabilities of the individual edges in the rearranged
    pedigree that were not used to predict the haplotypes of the unknown individuals. An edge is
    the connection between a parent and a child in the pedigree. The probability of an edge is the
    probability of the child's haplotype given the parent's haplotype. The probability of an edge
    is calculated by multiplying the probabilities of the child's alleles given the parent's
    alleles. The probability of an allele is calculated by the mutation rate of the allele's
    marker.
    """
    average_pedigree_probability = Decimal(0)

    progress_bar = reporter.progress_bar(
        total=number_of_iterations, desc="Calculating pedigree probability"
    )

    ordered_unknown_ids = [
        individual.id
        for individual in pedigree.get_level_order_traversal(suspect_name)
        if individual.haplotype_class == "unknown"
    ]

    with progress_bar:
        for i in range(number_of_iterations):
            pedigree_deep_copy = deepcopy(pedigree)

            simulation_probability = simulate_pedigree_probability(
                pedigree=pedigree_deep_copy,
                ordered_unknown_ids=ordered_unknown_ids,
                marker_set=marker_set,
                random=random,
            )

            average_pedigree_probability = (
                (average_pedigree_probability * i) + simulation_probability
            ) / (i + 1)

            progress_bar.update(i)

    reporter.log(
        f"Average pedigree probability Pr(Hkn=hkn) after {number_of_iterations} "
        f"iterations: {average_pedigree_probability}"
    )

    return average_pedigree_probability


def simulate_l_matching_haplotypes(
    pedigree: Pedigree,
    suspect_id: int,
    ordered_unknown_ids: Collection[int],
    marker_set: MarkerSet,
    l: int,
    average_pedigree_probability: Decimal,
    random: Random,
) -> float:
    simulation_probability = Decimal(1)

    suspect = pedigree.get_individual_by_id(suspect_id)
    unknown_individuals = pedigree.get_unknown_individuals()

    if l > 0:
        random_unknown_individuals = random.sample(unknown_individuals, l)
        simulation_probability /= comb(len(unknown_individuals), l)

        for individual in random_unknown_individuals:
            individual.haplotype_class = "fixed"
            for marker in marker_set.markers:
                individual.add_allele(
                    marker, suspect.get_allele_by_marker_name(marker.name).value
                )

    for individual_id in ordered_unknown_ids:
        individual = pedigree.get_individual_by_id(individual_id)

        if individual.haplotype_class == "fixed":
            continue

        parent = pedigree.get_parent(individual)

        individual.haplotype = mutate_haplotype(
            source=parent, marker_set=marker_set, random=random
        )
        individual.haplotype_class = "simulated"

        pedigree.set_relationship_class(parent, individual, "simulated")

        edge_probability = get_edge_probability(parent, individual, marker_set)
        simulation_probability *= Decimal(edge_probability)

    # Calculate the probability of the entire pedigree (Pr(H = h))
    pedigree.calculate_allele_probabilities(marker_set)
    pedigree_probability = pedigree.get_pedigree_probability(marker_set, "all")

    # Calculate (Pr(Huk = huk|Hkn = hkn))
    conditional_probability = (
        Decimal(pedigree_probability) / average_pedigree_probability
    )

    # Check if the total number of matching haplotypes is equal to l (excluding the suspect)
    number_of_matching_haplotypes = 0
    for individual in pedigree.individuals:
        if individual.name != suspect.name and individual.haplotype_class != "known":
            if individual.has_same_haplotype_as(suspect):
                number_of_matching_haplotypes += 1

    return (
        conditional_probability / simulation_probability
        if number_of_matching_haplotypes == l and simulation_probability
        else 0
    )


def calculate_l_matching_haplotypes(
    pedigree: Pedigree,
    marker_set: MarkerSet,
    suspect_name: str,
    number_of_iterations: int,
    l: int,
    average_pedigree_probability: Decimal,
    random: Random,
    reporter: Reporter,
) -> Decimal:
    """
    This function simulates the pedigree and counts the number of times the simulation results in
    l matching haplotypes. First l number of unknown individuals are selected, and their
    haplotypes are set to the same haplotype as the suspect. Then the haplotypes of the remaining
    unknown individuals are predicted based on the haplotypes of the surrounding known
    individuals, and the given mutation rates. Each time a random edge with one unknown individual
    and one known individual is selected. The haplotype of the unknown individual is predicted
    based on the haplotype of the known individual.

    This calculation keeps track of two things:
    - The simulation probability. This corresponds to (Pr(H̃uk = huk|Hkn = hkn))
    - The pedigree probability.

    The simulation probability is the probability of this specific iteration of the simulation.
    It is calculated by the following steps:
    - Set the haplotypes of l unknown individuals to the same haplotype as the suspect.

    This results in a probability of l/total number of unknown individuals.
    """
    l_probability = Decimal()

    suspect = pedigree.get_individual_by_name(suspect_name)

    progress_bar = reporter.progress_bar(
        total=number_of_iterations, desc=f"Calculating {l} matching haplotypes"
    )

    ordered_unknown_ids = [
        individual.id
        for individual in pedigree.get_level_order_traversal(suspect.name)
        if individual.haplotype_class == "unknown"
    ]

    with progress_bar:
        for i in range(number_of_iterations):
            pedigree_deep_copy = deepcopy(pedigree)

            simulation_probability = simulate_l_matching_haplotypes(
                pedigree=pedigree_deep_copy,
                suspect_id=suspect.id,
                ordered_unknown_ids=ordered_unknown_ids,
                marker_set=marker_set,
                l=l,
                average_pedigree_probability=average_pedigree_probability,
                random=random,
            )

            l_probability = ((l_probability * i) + simulation_probability) / (i + 1)

            progress_bar.update(i)

    reporter.log(f"Probability of {l} matching haplotypes: {l_probability}")

    return l_probability


def calculate_proposal_distribution(
    pedigree: Pedigree,
    marker_set: MarkerSet,
    suspect_name: str,
    average_pedigree_probability: Decimal,
    number_of_iterations: int,
    random: Random,
    reporter: Reporter,
) -> Mapping[int, Decimal]:
    unknown_individuals = pedigree.get_unknown_individuals()
    number_of_unknowns = len(unknown_individuals)
    proposal_distribution: dict[int, Decimal] = {}

    for l in range(0, number_of_unknowns + 1):  # l is the number of matching haplotypes
        proposal_distribution[l] = calculate_l_matching_haplotypes(
            pedigree=pedigree,
            marker_set=marker_set,
            suspect_name=suspect_name,
            l=l,
            average_pedigree_probability=average_pedigree_probability,
            number_of_iterations=number_of_iterations,
            random=random,
            reporter=reporter,
        )

    return proposal_distribution


def run_simulation(
    pedigree: Pedigree,
    suspect_name: str,
    marker_set: MarkerSet,
    number_of_iterations: int,
    random: Random,
    reporter: Reporter,
) -> SimulationResult:
    """
    The total simulation is based on a Monte-Carlo simulation with Importance Sampling. This
    means that the simulation uses a new distribution where the probability of observing the
    suspect's haplotype is high. This simulation uses a proposal distribution where the number
    of matching haplotypes is set to a fixed number.
    """

    """
    Step 1: calculate pedigree probability. 
    This needs to be done only once for the original pedigree
    This corresponds to Pr(Hkn=hkn)
    """
    start_time = datetime.now()

    average_pedigree_probability = calculate_pedigree_probability(
        pedigree, suspect_name, marker_set, number_of_iterations, random, reporter
    )

    run_time_pedigree_probability = start_time.now() - start_time

    if average_pedigree_probability == 0:
        return SimulationResult(
            average_pedigree_probability=average_pedigree_probability,
            proposal_distribution={},
            run_time_pedigree_probability=run_time_pedigree_probability,
            run_time_proposal_distribution=timedelta(0),
        )

    """
    Step 2: calculate the number of matching haplotypes with the suspect.
    This is done for the total number of individuals in the pedigree.
    This corresponds to Pr(Huk=hi, Hkn=hkn)
    """
    start_time = datetime.now()

    l_matching_haplotypes_probability = calculate_proposal_distribution(
        pedigree=pedigree,
        marker_set=marker_set,
        suspect_name=suspect_name,
        average_pedigree_probability=average_pedigree_probability,
        number_of_iterations=number_of_iterations,
        random=random,
        reporter=reporter,
    )

    run_time_proposal_distribution = start_time.now() - start_time

    return SimulationResult(
        average_pedigree_probability=average_pedigree_probability,
        proposal_distribution=l_matching_haplotypes_probability,
        run_time_pedigree_probability=run_time_pedigree_probability,
        run_time_proposal_distribution=run_time_proposal_distribution,
    )
