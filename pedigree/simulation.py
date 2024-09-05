from copy import deepcopy
from datetime import datetime
from decimal import Decimal
from itertools import combinations
from random import Random
from typing import Mapping

from pedigree.models import MarkerSet, Pedigree, SimulationResult, Individual, Marker, Allele, get_edge_probability
from pedigree.reporting import Reporter


def mutate_allele(marker: Marker, source_allele: Allele, random: Random) -> int:
    mutation_rate = marker.mutation_rate
    mutation_step = random.choices([0, 1], weights=[1 - mutation_rate, mutation_rate])[0]
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


def calculate_pedigree_probability(
    pedigree: Pedigree,
    marker_set: MarkerSet,
    suspect_name: str,
    number_of_iterations: int,
    random: Random,
    reporter: Reporter
) -> float:
    """
    This function calculates the pedigree probability after a number of iterations. This corresponds to Pr(Hkn=hkn).
    This probability is only calculated once for a given pedigree. The result is stored in the simulation object.
    The calculation is performed on the rearranged pedigree, where the suspect is the root node.
    During each iteration, the haplotypes of the unknown individuals are predicted based on the haplotypes
    of the surrounding known individuals, and the given mutation rates.

    The order in which the unknown individuals are processed is determined by the level order traversal
    of the rearranged pedigree. This means that the (unknown) children of the suspect are processed first,
    then the (unknown) children of the children, and so on. For each unknown individual, the haplotype of its
    father (parent node) is used to predict its own haplotype.

    After all unknown individuals have been processed, the pedigree probability is calculated. This is done by
    multiplying the probabilities of the individual edges in the rearranged pedigree that were not used to predict
    the haplotypes of the unknown individuals. An edge is the connection between a parent and a child in the pedigree.
    The probability of an edge is the probability of the child's haplotype given the parent's haplotype.
    The probability of an edge is calculated by multiplying the probabilities of the child's alleles given the parent's
    alleles. The probability of an allele is calculated by the mutation rate of the allele's marker.
    """

    average_pedigree_probability = 0

    progress_bar = reporter.progress_bar(
        total=number_of_iterations, desc="calculating pedigree probability"
    )

    with progress_bar:
        for i in range(number_of_iterations):
            pedigree_deep_copy = deepcopy(pedigree)

            level_order_traversal = pedigree_deep_copy.get_level_order_traversal(suspect_name)
            for individual in level_order_traversal:
                if individual.haplotype_class == "unknown":
                    parent = pedigree_deep_copy.get_parent(individual)
                    mutate_haplotype(
                        source=parent, target=individual, marker_set=marker_set, random=random
                    )

                    # surrounding_known_individuals = pedigree_deep_copy.get_surrounding_known_individuals(individual)
                    # random_known_node = random.choice(surrounding_known_individuals)
                    # mutate_haplotype(source=random_known_node, target=individual, marker_set=marker_set)
                    pedigree_deep_copy.set_relationship_class(parent, individual, "simulated")

            pedigree_deep_copy.calculate_allele_probabilities(marker_set)
            pedigree_probability = pedigree_deep_copy.get_pedigree_probability(marker_set, "unused")
            # alleles = pedigree_deep_copy.get_individual_by_name("Unknown").haplotype.get_alleles()
            # total_alleles_counter[alleles[0].value] += 1
            # total_probability_counter[alleles[0].value] += pedigree_probability

            average_pedigree_probability = ((average_pedigree_probability * i) + pedigree_probability) / (i + 1)

            progress_bar.update(i)

        # pedigree_deep_copy.print_pedigree()
        del pedigree_deep_copy

    reporter.log(
        f"Average pedigree probability Pr(Hkn=hkn) after {number_of_iterations} "
        f"iterations: {average_pedigree_probability}"
    )

    # for key, value in total_probability_counter.items():
    #     reporter.log((f"Allele {key}: {total_alleles_counter[key]/100000}")
    #     reporter.log((f"Allele {key}: {value / total_alleles_counter[key]}")

    return average_pedigree_probability


def calculate_l_matching_haplotypes(
    pedigree: Pedigree,
    marker_set: MarkerSet,
    suspect_name: str,
    number_of_iterations: int,
    l: int,
    average_pedigree_probability: float,
    random: Random,
    reporter: Reporter
) -> float:
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

    total_l_probability = Decimal()

    progress_bar = reporter.progress_bar(
        total=number_of_iterations, desc=f"calculating {l} matching haplotypes"
    )

    with progress_bar:
        for i in range(number_of_iterations):
            simulation_probability = 1

            pedigree_deep_copy = deepcopy(pedigree)
            unknown_individuals = pedigree_deep_copy.get_unknown_individuals()

            if l > 0:
                random_unknown_individuals = random.sample(unknown_individuals, l)
                simulation_probability = 1 / len(list(combinations(unknown_individuals, l)))

                for individual in random_unknown_individuals:
                    individual.haplotype_class = "fixed"
                    for marker in marker_set.markers:
                        individual.add_allele(marker, suspect.get_allele_by_marker_name(marker.name).value)

            # while len(pedigree_deep_copy.get_unknown_individuals()) > 0:
            #     edges_with_one_unknown_and_one_known_individual = list(pedigree_deep_copy.get_edges_with_one_unknown_and_one_known_individual())
            #     simulation_probability *= (1 / len(edges_with_one_unknown_and_one_known_individual))
            #     random_unknown_individual, random_known_individual = random.choice(edges_with_one_unknown_and_one_known_individual)
            #     mutate_haplotype(source=random_known_individual, target=random_unknown_individual, marker_set=marker_set)
            #     pedigree_deep_copy.set_relationship_class(random_known_individual, random_unknown_individual, "simulated")
            #     edge_probability = get_edge_probability(marker_set, random_known_individual, random_unknown_individual)
            #     simulation_probability *= edge_probability

            level_order_traversal = pedigree_deep_copy.get_level_order_traversal(suspect_name)
            for individual in level_order_traversal:
                if individual.haplotype_class == "unknown":
                    parent = pedigree_deep_copy.get_parent(individual)
                    mutate_haplotype(
                        source=parent, target=individual, marker_set=marker_set, random=random
                    )
                    pedigree_deep_copy.set_relationship_class(parent, individual, "simulated")
                    edge_probability = get_edge_probability(marker_set, parent, individual)
                    simulation_probability *= edge_probability

            # Calculate the probability of the entire pedigree (Pr(H = h))
            pedigree_deep_copy.calculate_allele_probabilities(marker_set)
            pedigree_probability = pedigree_deep_copy.get_pedigree_probability(marker_set, "all")

            # Calculate (Pr(Huk = huk|Hkn = hkn))
            try:
                conditional_probability = Decimal(pedigree_probability) / Decimal(average_pedigree_probability)
            except:
                conditional_probability = 0

            # Check if the total number of matching haplotypes is equal to l (excluding the suspect)
            number_of_matching_haplotypes = 0
            for individual in pedigree_deep_copy.individuals:
                if individual.name != suspect_name and individual.haplotype_class != "known":
                    if individual.has_same_haplotype_as(suspect):
                        number_of_matching_haplotypes += 1

            if number_of_matching_haplotypes == l:
                try:
                    total_l_probability += (Decimal(conditional_probability) / Decimal(simulation_probability))
                except:
                    total_l_probability += 0

            del pedigree_deep_copy

            progress_bar.update(i)

    l_probability = (total_l_probability / number_of_iterations) if number_of_iterations else 0

    reporter.log(
        f"Probability of {l} matching haplotypes: {l_probability}"
    )

    return l_probability


def calculate_proposal_distribution(
    pedigree: Pedigree,
    marker_set: MarkerSet,
    suspect_name: str,
    number_of_iterations: int,
    average_pedigree_probability: float,
    random: Random,
    reporter: Reporter,
) -> Mapping[int, float]:
    unknown_individuals = pedigree.get_unknown_individuals()
    number_of_unknown_individuals = len(unknown_individuals)
    proposal_distribution: dict[int, float] = {}

    for l in range(0, number_of_unknown_individuals + 1):  # l is the number of matching haplotypes
        proposal_distribution[l] = calculate_l_matching_haplotypes(
            pedigree=pedigree,
            marker_set=marker_set,
            suspect_name=suspect_name,
            number_of_iterations=number_of_iterations,
            l=l,
            average_pedigree_probability=average_pedigree_probability,
            random=random,
            reporter=reporter
        )

    return proposal_distribution


def run_simulation(
    pedigree: Pedigree,
    marker_set: MarkerSet,
    suspect_name: str,
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
    start_time = datetime.now()

    """
    Step 1: calculate pedigree probability. 
    This needs to be done only once for the original pedigree
    This corresponds to Pr(Hkn=hkn)
    """
    average_pedigree_probability = calculate_pedigree_probability(
        pedigree, marker_set, suspect_name, number_of_iterations, random, reporter
    )

    run_time_pedigree_probability = start_time.now() - start_time
    start_time = datetime.now()

    """
    Step 2: calculate the number of matching haplotypes with the suspect.
    This is done for the total number of individuals in the pedigree.
    This corresponds to Pr(Huk=hi, Hkn=hkn)
    """
    l_matching_haplotypes_probability = calculate_proposal_distribution(
        pedigree=pedigree,
        marker_set=marker_set,
        suspect_name=suspect_name,
        number_of_iterations=number_of_iterations,
        average_pedigree_probability=average_pedigree_probability,
        random=random,
        reporter=reporter
    )

    run_time_proposal_distribution = start_time.now() - start_time

    return SimulationResult(
        average_pedigree_probability=average_pedigree_probability,
        proposal_distribution=l_matching_haplotypes_probability,
        run_time_pedigree_probability=run_time_pedigree_probability,
        run_time_proposal_distribution=run_time_proposal_distribution,
    )
