import functools
from copy import deepcopy
from multiprocessing import Pool
from random import Random
from typing import Callable

from pedigree.models import MarkerSet, Pedigree, SimulationResult, Individual, IterationResult, \
    Marker, Allele


def mutate_allele(marker: Marker, source_allele: Allele, random: Random) -> int:
    mutation_rate = marker.mutation_rate
    mutation_step = random.choices(
        [0, 1, 2],
        weights=[1 - mutation_rate, mutation_rate * 0.97, mutation_rate * 0.03],
    )[0]
    mutation_direction = random.choice([-1, 1])
    return source_allele.value + (mutation_step * mutation_direction)


def mutate_haplotype(
    source: Individual, target: Individual, marker_set: MarkerSet, random: Random
) -> None:
    target.haplotype_class = "simulated"
    for marker in marker_set.markers:
        source_allele = source.get_allele_by_marker_name(marker.name)
        mutated_allele_value = mutate_allele(marker, source_allele, random)
        target.add_allele(marker, mutated_allele_value)


def calculate_iteration(
    pedigree: Pedigree, marker_set: MarkerSet, suspect: Individual, random: Random, _: int,
) -> IterationResult:
    pedigree_deep_copy = deepcopy(pedigree)

    while len(pedigree_deep_copy.get_unknown_individuals()) > 0:
        edges_with_one_unknown_and_one_known_individual = list(
            pedigree_deep_copy.get_edges_with_one_unknown_and_one_known_individual()
        )
        random_unknown_individual, random_known_individual = random.choice(
            edges_with_one_unknown_and_one_known_individual
        )
        mutate_haplotype(
            source=random_known_individual,
            target=random_unknown_individual,
            marker_set=marker_set,
            random=random
        )
        pedigree_deep_copy.set_relationship_class(
            random_known_individual, random_unknown_individual, "simulated"
        )

    pedigree_deep_copy.calculate_allele_probabilities(marker_set)
    pedigree_probability = pedigree_deep_copy.get_pedigree_probability(
        marker_set, "unused"
    )

    matching_haplotypes = 0
    for individual in pedigree_deep_copy.individuals:
        if individual.name != suspect.name and individual.name not in []:
            if individual.has_same_haplotype_as(suspect):
                matching_haplotypes += 1

    del pedigree_deep_copy

    return IterationResult(pedigree_probability, matching_haplotypes)


def calculate_l_matching_haplotypes(
    pedigree: Pedigree,
    marker_set: MarkerSet,
    suspect_name: str,
    number_of_iterations: int,
    random: Random,
    show_progress: Callable[[int, int], None] | None = None
) -> SimulationResult:
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
    suspect = pedigree.get_individual_by_name(suspect_name)

    with Pool(30) as pool:
        partial = functools.partial(
            calculate_iteration, pedigree, marker_set, suspect, random
        )

        results = []
        for i, result in enumerate(pool.imap_unordered(partial, range(number_of_iterations))):
            results.append(result)

            if show_progress:
                show_progress(i + 1, number_of_iterations)

    total_correct = 0
    total_l_matches = {}

    for result in results:
        if result.pedigree_probability > 0:
            total_correct += 1
            if result.matching_haplotype_count in total_l_matches:
                total_l_matches[result.matching_haplotype_count] += result.pedigree_probability
            else:
                total_l_matches[result.matching_haplotype_count] = result.pedigree_probability

    total_l_matches_sum = sum([v for k, v in total_l_matches.items()])
    total_l_matches_normalized = {
        k: v / total_l_matches_sum for k, v in total_l_matches.items()
    }
    total_0_count = total_l_matches.get(0, 0)
    total_not_0_count = sum([v for k, v in total_l_matches.items() if k != 0])

    return SimulationResult(
        iteration_results=results,
        total_correct=total_correct,
        total_l_matches_normalized=total_l_matches_normalized,
        total_0_count=total_0_count,
        total_not_0_count=total_not_0_count,
        average_pedigree_probability=0,  # TODO: what should this be?
        l_matching_haplotypes_probability={},  # TODO: what should this be?
    )


def run_simulation(
    pedigree: Pedigree,
    marker_set: MarkerSet,
    suspect: str,
    number_of_iterations: int,
    random: Random,
    show_progress: Callable[[int, int], None] | None = None
) -> SimulationResult:
    """
    The total simulation is based on a Monte-Carlo simulation with Importance Sampling.
    This means that the simulation uses a new distribution where the probability of observing the
    suspect's haplotype is high. This simulation uses a proposal distribution where the number of
    matching haplotypes is set to a fixed number.

    Calculate the number of matching haplotypes with the suspect. This is done for the total
    number of individuals in the pedigree. This corresponds to Pr(Huk=hi, Hkn=hkn)
    """
    return calculate_l_matching_haplotypes(
        pedigree, marker_set, suspect, number_of_iterations, random, show_progress
    )
