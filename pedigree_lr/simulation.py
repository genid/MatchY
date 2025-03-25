import math
import operator
from collections.abc import Set
from copy import deepcopy
from datetime import datetime, timedelta
from decimal import Decimal
from functools import reduce
from idlelib.replace import replace
from math import comb, factorial
import numpy as np
from random import Random
from typing import Collection, Mapping, Sequence
from pedigree_lr.reporting import Reporter
from pedigree_lr.models import (
    Allele,
    Haplotype,
    Individual,
    IterationResult,
    Marker,
    MarkerSet,
    Pedigree,
    Relationship,
    SimulationResult,
    calculate_mutation_probability,
    get_single_copy_mutation_rate
)


def mutate_allele(
        marker: Marker,
        source_alleles: list[Allele],
        random: Random,
        two_step_mutation_factor: float = 0.03,  # TODO: make this a parameter
) -> list[Allele]:
    """
        Applies mutation to a list of alleles based on the given marker's mutation rate.

        This function determines the mutation rate per allele, including a small chance
        for two-step mutations. Each allele undergoes mutation with a probability
        derived from the mutation rate. The mutation step is chosen randomly, and its
        direction (increase or decrease) is also randomly assigned.

        Args:
            marker (Marker): The genetic marker associated with the alleles, which
                determines the mutation rate.
            source_alleles (list[Allele]): The list of alleles to mutate.
            random (Random): A random number generator used for mutation decisions.
            two_step_mutation_factor (float): The factor that determines the probability
                of a two-step mutation. Defaults to 0.03.

        Returns:
            list[Allele]: A new list of mutated alleles with potentially changed values.

        Notes:
            - The two-step mutation rate is currently fixed at 3% of the base mutation rate.
            - Alleles cannot mutate below a value of 1 to avoid negative values.
            - Intermediate allele values are preserved without mutation.
        """

    mutation_rate = get_single_copy_mutation_rate(marker.mutation_rate, marker.number_of_copies)
    two_step_mutation_rate = mutation_rate * two_step_mutation_factor
    mutation_rate = mutation_rate * (1 - two_step_mutation_factor)

    mutated_alleles = []
    for source_allele in source_alleles:
        mutation_step = random.choices([0, 1, 2],
                                       weights=[1 - mutation_rate, mutation_rate, two_step_mutation_rate])[0]
        mutation_direction = random.choice([-1, 1])  # Assumption that direction is symmetric
        mutated_allele_value = max(1, source_allele.value + (mutation_step * mutation_direction))  # Account for negative values

        # Intermediate allele values are preserved without mutation for now
        mutated_intermediate_allele_value = source_allele.intermediate_value

        mutated_alleles.append(Allele(marker, mutated_allele_value, mutated_intermediate_allele_value))
    return mutated_alleles


def mutate_haplotype(
        source: Haplotype,
        marker_set: MarkerSet,
        random: Random
) -> Haplotype:
    haplotype = Haplotype()

    for marker in marker_set.markers:
        source_alleles = source.alleles[marker.name]
        haplotype.alleles[marker.name] = mutate_allele(marker, source_alleles, random)

    return haplotype


def get_edge_probability(
        known: Haplotype,
        unknown: Haplotype,
        marker_set: MarkerSet
) -> Decimal:
    edge_probability = Decimal(1)

    for marker in marker_set.markers:
        known_alleles = sorted(known.alleles[marker.name], key=lambda allele: allele.value)
        unknown_alleles = sorted(unknown.alleles[marker.name], key=lambda allele: allele.value)

        mutation_probability = calculate_mutation_probability(known_alleles, unknown_alleles, marker)
        edge_probability *= mutation_probability

    return edge_probability


def get_edge_probabilities(
        haplotypes: Mapping[int, Haplotype],
        relationships: Collection[Relationship],
        marker_set: MarkerSet,
) -> Mapping[tuple[int, int], Decimal]:
    return {
        (relationship.parent_id, relationship.child_id): get_edge_probability(
            haplotypes[relationship.parent_id],
            haplotypes[relationship.child_id],
            marker_set=marker_set,
        )
        for relationship in relationships
    }


def update_average(
        new_probability: Decimal,
        old_probability: Decimal,
        i: int
) -> Decimal:
    return ((new_probability * i) + old_probability) / (
            i + 1
    )


def stable(
        probabilities: list[Decimal],
        i: int,
        window: int = 500,
        min_iterations: int = 2000,
        threshold: float = 0.0001
) -> bool:
    if i > window and i > min_iterations:
        if all(
                abs((probabilities[i] - probabilities[i - 1]) / probabilities[i - 1]) < threshold
                if probabilities[i - 1] != 0 else probabilities[i] == 0
                for i in range(i - window, i)
        ):
            print(f"Stable after {i} iterations")
            return True
    return False


def model_valid(
        probabilities: list[Decimal],
        threshold: float = 0.005
) -> bool:
    return all(
        abs((probabilities[i] - probabilities[i - 1]) / probabilities[i - 1]) < threshold
        if probabilities[i - 1] != 0 else probabilities[i] == 0
        for i in range(1, len(probabilities))
    )


def create_simulated_pedigree(
        individuals: Mapping[int, Individual],
        relationships: Collection[Relationship],
        haplotypes: Mapping[int, Haplotype],
        simulated_individual_ids: Set[int],
        simulated_relationship_ids: Set[tuple[int, int]],
        fixed_individual_ids: Set[int],
        marker_set: MarkerSet,
) -> Pedigree:
    pedigree = Pedigree(
        individuals=[
            Individual(
                id=individual.id,
                name=individual.name,
                haplotype=haplotypes[individual.id],
                haplotype_class="simulated"
                if individual.id in simulated_individual_ids
                else (
                    "fixed"
                    if individual.id in fixed_individual_ids
                    else individual.haplotype_class
                ),
            )
            for individual in individuals.values()
        ],
        relationships=[
            Relationship(
                parent_id=relationship.parent_id,
                child_id=relationship.child_id,
                edge_class="simulated"
                if (relationship.parent_id, relationship.child_id)
                   in simulated_relationship_ids
                else relationship.edge_class,
            )
            for relationship in relationships
        ],
    )
    pedigree.calculate_allele_probabilities(marker_set)
    return pedigree


def simulate_pedigree_probability(
        individuals: Mapping[int, Individual],
        relationships: Collection[Relationship],
        parents: Mapping[int, int],
        ordered_unknown_ids: Collection[int],
        marker_set: MarkerSet,
        random: Random,
        show_simulated_pedigrees: bool = False,
) -> IterationResult:
    simulated_individual_ids: set[int] = set()
    simulated_relationship_ids: set[tuple[int, int]] = set()

    haplotypes = {
        individual.id: individual.haplotype for individual in individuals.values()
    }

    for individual_id in ordered_unknown_ids:
        parent_id = parents[individual_id]
        parent_haplotype = haplotypes[parent_id]

        haplotypes[individual_id] = mutate_haplotype(
            source=parent_haplotype, marker_set=marker_set, random=random
        )

        simulated_individual_ids.add(individual_id)
        simulated_relationship_ids.add((parent_id, individual_id))

    unused_relationships = [
        relationship
        for relationship in relationships
        if (relationship.parent_id, relationship.child_id)
           not in simulated_relationship_ids
    ]

    edge_probabilities = get_edge_probabilities(
        haplotypes, unused_relationships, marker_set
    )

    # pedigree_probability = reduce(operator.mul, edge_probabilities.values(), 1) # alternative way to calculate pedigree probability

    if any(probability == 0 for probability in edge_probabilities.values()):
        pedigree_probability = Decimal(0)
    else:
        log_sum = sum(math.log(p) for p in edge_probabilities.values())  # if len edge_probabilities > 0 else 1
        pedigree_probability = math.exp(log_sum)  # TODO: Solve underflow

    return IterationResult(
        pedigree=create_simulated_pedigree(
            individuals=individuals,
            relationships=relationships,
            haplotypes=haplotypes,
            simulated_individual_ids=simulated_individual_ids,
            simulated_relationship_ids=simulated_relationship_ids,
            fixed_individual_ids=set(),
            marker_set=marker_set,
        )
        if show_simulated_pedigrees
        else None,
        probability=Decimal(pedigree_probability),
        edge_probabilities=edge_probabilities,
    )


def calculate_pedigree_probability(
        pedigree: Pedigree,
        suspect_name: str,
        marker_set: MarkerSet,
        number_of_iterations: int,
        random: Random,
        reporter: Reporter,
        show_simulated_pedigrees: bool,
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

    progress_bar = reporter.progress_bar(
        total=number_of_iterations, desc="Calculating pedigree probability"
    )

    individuals = {individual.id: individual for individual in pedigree.individuals}

    parents = {
        relationship.child_id: relationship.parent_id
        for relationship in pedigree.relationships
    }

    ordered_unknown_ids = [
        individual.id
        for individual in pedigree.get_level_order_traversal(suspect_name)
        if individual.haplotype_class == "unknown"
    ]

    average_pedigree_probability = Decimal()
    pedigree_probabilities: list[Decimal] = []
    edge_probabilities: dict[tuple[int, int], list[Decimal]] = {
        (relationship.parent_id, relationship.child_id): []
        for relationship in pedigree.relationships
    }

    with (progress_bar):
        for i in range(number_of_iterations):  # TODO: dynamics number of iterations
            iteration_result = simulate_pedigree_probability(
                individuals=individuals,
                relationships=pedigree.relationships,
                parents=parents,
                ordered_unknown_ids=ordered_unknown_ids,
                marker_set=marker_set,
                random=random,
                show_simulated_pedigrees=show_simulated_pedigrees,
            )

            if iteration_result.pedigree:
                reporter.log(iteration_result.pedigree.to_string())

            average_pedigree_probability = update_average(
                average_pedigree_probability, iteration_result.probability, i
            )

            pedigree_probabilities.append(average_pedigree_probability)

            # if stable(pedigree_probabilities, i):
            #     break

            edge_probabilities = {
                edge: probabilities + [iteration_result.edge_probabilities.get(edge, 0)]
                for edge, probabilities in edge_probabilities.items()
            }
            progress_bar.update(i)

    average_edge_probabilities = {
        edge: sum(edge_probabilities) / number_of_iterations
        for edge, edge_probabilities in edge_probabilities.items()
    }

    reporter.log(
        f"Average pedigree probability Pr(Hkn=hkn) after {number_of_iterations} iterations: "
        f"{average_pedigree_probability} (log {10 * math.log10(average_pedigree_probability)})"
    )

    # reporter.log(f"Average edge probabilities after {number_of_iterations} iterations:")
    # for (child_id, parent_id), edge_probability in average_edge_probabilities.items():
    #     reporter.log(
    #         f"{individuals[child_id].name} -> {individuals[parent_id].name} = {edge_probability}"
    #     )

    return average_pedigree_probability


def simulate_l_matching_haplotypes(
        individuals: Mapping[int, Individual],
        relationships: Collection[Relationship],
        parents: Mapping[int, int],
        suspect: Individual,
        ordered_unknown_ids: Sequence[int],
        marker_set: MarkerSet,
        l: int,
        average_pedigree_probability: Decimal,
        random: Random,
        show_simulated_pedigrees: bool,
) -> IterationResult:
    simulated_individual_ids: set[int] = set()
    simulated_relationship_ids: set[tuple[int, int]] = set()
    fixed_individual_ids: set[int] = set()

    simulation_probability = Decimal(1)

    if l > 0:
        total_picking_probability = Decimal(1)

        for _ in range(l):
            # get a list of all remaining unknown individuals that are not excluded and not fixed
            available_individuals_ids = [id for id in ordered_unknown_ids if
                                         not individuals[id].exclude and id not in fixed_individual_ids]
            # get the picking probabilities of the available individuals
            picking_probabilities = {id:
                                         individuals[id].picking_probability for id in available_individuals_ids
                                     }
            # normalize the picking probabilities
            normalized_picking_probabilities = {
                id: probability / sum(picking_probabilities.values()) for id, probability in
                picking_probabilities.items()
            }
            # pick an individual based on the normalized picking probabilities
            picked_individual_id = random.choices(
                available_individuals_ids,
                weights=[float(normalized_picking_probabilities[id]) for id in available_individuals_ids]
            )[0]
            # add the picked individual to the fixed individuals
            fixed_individual_ids.add(picked_individual_id)
            # get the picking probability of the picked individual
            picking_probability = normalized_picking_probabilities[picked_individual_id]
            # multiply the total picking probability with the picking probability of the picked individual
            total_picking_probability *= picking_probability

        unordered_picking_score = (total_picking_probability * factorial(l))
        simulation_probability *= unordered_picking_score


    haplotypes = {
        individual.id: (
            suspect.haplotype
            if individual.id in fixed_individual_ids
            else individual.haplotype
        )
        for individual in individuals.values()
    }

    for individual_id in ordered_unknown_ids:
        if individual_id in fixed_individual_ids:
            continue

        parent_id = parents[individual_id]
        parent_haplotype = haplotypes[parent_id]

        haplotypes[individual_id] = mutate_haplotype(
            source=parent_haplotype, marker_set=marker_set, random=random
        )

        simulated_individual_ids.add(individual_id)
        simulated_relationship_ids.add((parent_id, individual_id))

    # Calculate the probability of the entire pedigree (Pr(H = h))
    all_edge_probabilities = get_edge_probabilities(
        haplotypes, relationships, marker_set
    )
    pedigree_probability = reduce(operator.mul, all_edge_probabilities.values(), 1)

    # Calculate the probability of the simulated edges
    simulated_edge_probabilities = (
        all_edge_probabilities[edge] for edge in simulated_relationship_ids
    )
    simulation_probability = reduce(
        operator.mul, simulated_edge_probabilities, simulation_probability
    )

    # Calculate (Pr(Huk = huk|Hkn = hkn))
    if average_pedigree_probability == 0:
        conditional_probability = Decimal(0)
    else:
        conditional_probability = pedigree_probability / average_pedigree_probability

    # Check if the total number of matching haplotypes is equal to l (excluding the suspect)
    number_of_matching_haplotypes = sum(
        haplotypes[individual_id] == suspect.haplotype and not individuals[individual_id].exclude
        # TODO: check if this is correct
        for individual_id in simulated_individual_ids | fixed_individual_ids
    )

    probability = (
        conditional_probability / simulation_probability
        if number_of_matching_haplotypes == l and simulation_probability and simulation_probability != 0
        else Decimal(0)
    )  # TODO: use average of three models

    return IterationResult(
        pedigree=create_simulated_pedigree(
            individuals=individuals,
            relationships=relationships,
            haplotypes=haplotypes,
            simulated_individual_ids=simulated_individual_ids,
            simulated_relationship_ids=simulated_relationship_ids,
            fixed_individual_ids=fixed_individual_ids,
            marker_set=marker_set,
        )
        if show_simulated_pedigrees
        else None,
        probability=probability,
        edge_probabilities=all_edge_probabilities,
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
        show_simulated_pedigrees: bool,
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

    average_edge_probabilities: dict[tuple[int, int], Decimal] = {
        (relationship.parent_id, relationship.child_id): Decimal(0)
        for relationship in pedigree.relationships
    }

    progress_bar = reporter.progress_bar(
        total=number_of_iterations * 3, desc=f"Calculating {l} matching haplotypes"
    )

    individuals = {individual.id: individual for individual in pedigree.individuals}

    parents = {
        relationship.child_id: relationship.parent_id
        for relationship in pedigree.relationships
    }

    ordered_unknown_ids = [
        individual.id
        for individual in pedigree.get_level_order_traversal(suspect_name)
        if individual.haplotype_class == "unknown"
    ]

    suspect = pedigree.get_individual_by_name(suspect_name)

    with (progress_bar):
        model_probabilities = []
        for m in range(3):  # Model validation
            l_probabilities = []
            for i in range(number_of_iterations):
                iteration_result = simulate_l_matching_haplotypes(
                    individuals=individuals,
                    relationships=pedigree.relationships,
                    parents=parents,
                    suspect=suspect,
                    ordered_unknown_ids=ordered_unknown_ids,
                    marker_set=marker_set,
                    l=l,
                    average_pedigree_probability=average_pedigree_probability,
                    random=random,
                    show_simulated_pedigrees=show_simulated_pedigrees,
                )

                l_probability = update_average(l_probability, iteration_result.probability, i)
                l_probabilities.append(l_probability)

                if stable(l_probabilities, i):
                    model_probabilities.append(l_probability)
                    break

                average_edge_probabilities = {
                    edge: update_average(average_edge_probability, iteration_result.edge_probabilities.get(edge, 0), i)
                    for edge, average_edge_probability in average_edge_probabilities.items()
                }

                progress_bar.update(i)

            with open(f"probability_test_results/l_{l}_model_{m}_probabilities.txt", "w") as f:
                for prob in l_probabilities:
                    f.write(f"{prob}\n")

        # check if all three model probabilities are stable (<0.5% change)
        if model_valid(model_probabilities):
            print(f"Model is valid")
        else:
            print(f"Model is not valid!")
            print(model_probabilities)

    reporter.log(f"Probability of {l} matching haplotypes: {l_probability}")

    reporter.log(
        f"Average edge probabilities of {l} matching haplotypes after "
        f"{number_of_iterations} iterations:"
    )
    for (child_id, parent_id), edge_probability in average_edge_probabilities.items():
        reporter.log(
            f"{individuals[child_id].name} -> {individuals[parent_id].name} = {edge_probability}"
        )

    return l_probability


def calculate_proposal_distribution(
        pedigree: Pedigree,
        marker_set: MarkerSet,
        suspect_name: str,
        average_pedigree_probability: Decimal,
        number_of_iterations: int,
        random: Random,
        reporter: Reporter,
        show_simulated_pedigrees: bool,
) -> Mapping[int, Decimal]:
    unknown_individuals = pedigree.get_unknown_individuals()
    number_of_unknowns = len(unknown_individuals)
    proposal_distribution: dict[int, Decimal] = {}
    number_of_excluded_individuals = len([individual for individual in pedigree.individuals if individual.exclude])

    # TODO: check if this is correct
    # TODO: l=0 (no matching haplotypes) should be calculated by 1 - sum(proposal_distribution)
    for l in range(1,
                   number_of_unknowns + 1 - number_of_excluded_individuals):  # l is the number of matching haplotypes
        proposal_distribution[l] = calculate_l_matching_haplotypes(
            pedigree=pedigree,
            marker_set=marker_set,
            suspect_name=suspect_name,
            l=l,
            average_pedigree_probability=average_pedigree_probability,
            number_of_iterations=number_of_iterations,
            random=random,
            reporter=reporter,
            show_simulated_pedigrees=show_simulated_pedigrees,
        )

    proposal_distribution[0] = Decimal(1) - sum(proposal_distribution.values())
    # TODO: check if proposal_distribution adds up to unity
    return proposal_distribution


def run_simulation(
        pedigree: Pedigree,
        suspect_name: str,
        marker_set: MarkerSet,
        number_of_iterations: int,
        random: Random,
        reporter: Reporter,
        show_simulated_pedigrees: bool,
) -> SimulationResult:
    """
    The total simulation is based on a Monte-Carlo simulation with Importance Sampling. This
    means that the simulation uses a new distribution where the probability of observing the
    suspect's haplotype is high. This simulation uses a proposal distribution where the number
    of matching haplotypes is set to a fixed number.
    """

    # Reroot the pedigree to have the suspect as the most recent common ancestor
    # This is important for the level order traversal in which the unknown individuals are processed
    pedigree.reroot_pedigree(suspect_name)

    # Calculate the a priori match probabilities for all unknown individuals in the pedigree.
    # These (normalized) probabilities are used to "randomly" pick the l number of individuals
    # which will match the suspect haplotype
    pedigree.calculate_picking_probabilities(marker_set)

    """
    Step 1: calculate pedigree probability. 
    This needs to be done only once for the original pedigree
    This corresponds to Pr(Hkn=hkn)
    """
    start_time = datetime.now()

    average_pedigree_probability = calculate_pedigree_probability(
        pedigree=pedigree,
        suspect_name=suspect_name,
        marker_set=marker_set,
        number_of_iterations=number_of_iterations,
        random=random,
        reporter=reporter,
        show_simulated_pedigrees=show_simulated_pedigrees,
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
        show_simulated_pedigrees=show_simulated_pedigrees,
    )

    run_time_proposal_distribution = start_time.now() - start_time

    extended_pedigree = deepcopy(pedigree)
    extended_pedigree.extend_pedigree()

    # TODO: run calculate_pedigree_probability for extended pedigree as well for outside match probability



    return SimulationResult(
        average_pedigree_probability=average_pedigree_probability,
        proposal_distribution=l_matching_haplotypes_probability,
        run_time_pedigree_probability=run_time_pedigree_probability,
        run_time_proposal_distribution=run_time_proposal_distribution,
    )
