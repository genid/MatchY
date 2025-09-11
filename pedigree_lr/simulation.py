from configparser import ConfigParser
from decimal import Decimal
import math
import multiprocessing
import operator
import pickle
from functools import partial
from itertools import combinations
from copy import deepcopy
from datetime import datetime, timedelta
from functools import reduce
from pathlib import Path
from random import Random, SystemRandom
from typing import Collection, Mapping, Sequence
import pandas as pd
from pedigree_lr.reporting import Reporter, create_html_pdf_report, ProgressBar
from pedigree_lr.visualization import plot_probabilities, save_pedigree_to_png
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
    get_single_copy_mutation_rate,
    SimulationParameters, InvalidAveragePedigreeProbability
)


def mutate_allele(
        marker: Marker,
        source_alleles: list[Allele],
        random_seed: int,
        two_step_mutation_factor: float,
) -> list[Allele]:
    """
        Applies mutation to a list of alleles based on the given marker's mutation rate.

        This function determines the mutation rate per allele, including a small chance
        for two-step mutations. Each allele undergoes mutation with a probability
        derived from the mutation rate. The mutation step is chosen randomly, and its
        direction (increase or decrease) is also randomly assigned.

        Args:
            marker (Marker): A marker from the marker set, with corresponding mutation rate.
            source_alleles (list[Allele]): The list of alleles to mutate.
            random_seed (Random): A random number generator used for mutation decisions.
            two_step_mutation_factor (float): The factor that determines the probability
                of a two-step mutation. Defaults to 0.03.

        Returns:
            list[Allele]: A new list of mutated alleles with potentially changed values.

        Notes:
            - The two-step mutation rate is currently fixed at 3% of the base mutation rate.
            - Alleles cannot mutate below a value of 1 to avoid negative values.
            - Intermediate allele values are preserved without mutation.
        """
    # TODO: add correct random seed handling
    random = SystemRandom()  # Ensure each mutation is reproducible with the same seed

    mutation_rate = get_single_copy_mutation_rate(marker.mutation_rate, marker.number_of_copies)
    two_step_mutation_rate = mutation_rate * two_step_mutation_factor
    single_step_mutation_rate = mutation_rate * (1 - two_step_mutation_factor)
    no_mutation_rate = 1 - mutation_rate

    mutated_alleles = []
    for source_allele in source_alleles:
        mutation_step = random.choices([0, 1, 2],
                                       weights=[no_mutation_rate, single_step_mutation_rate, two_step_mutation_rate])[0]
        mutation_direction = random.choice([-1, 1])  # Assumption that direction is symmetric
        mutated_allele_value = max(1, source_allele.value + (
                    mutation_step * mutation_direction))  # Account for negative values

        # Intermediate allele values are preserved without mutation for now
        mutated_intermediate_allele_value = source_allele.intermediate_value

        mutated_alleles.append(Allele(marker, mutated_allele_value, mutated_intermediate_allele_value))
    return sorted(mutated_alleles)


def mutate_haplotype(
        source: Haplotype,
        marker_set: MarkerSet,
        random_seed: int,
        two_step_mutation_factor: float,
) -> Haplotype:
    """
        Generates a mutated version of a given haplotype by applying mutations to its alleles.

        Args:
            source (Haplotype): The original haplotype to be mutated.
            marker_set (MarkerSet): A set of markers with corresponding mutation rates.
            random_seed (Random): A random number generator instance to introduce stochasticity in mutation.
            two_step_mutation_factor (float): The factor that determines the probability
                of a two-step mutation. Defaults to 0.03.

        Returns:
            Haplotype: A new haplotype instance with mutated alleles.

        Note:
            - Each marker in `marker_set` is iterated over, and its corresponding alleles in `source`
              undergo mutation through the `mutate_allele` function.
            - The mutation process is controlled by the `random` instance.
        """

    target_haplotype = Haplotype()  # Initialize a new haplotype instance

    for marker in marker_set.markers:
        source_alleles = source.alleles[marker.name]
        target_haplotype.alleles[marker.name] = mutate_allele(
            marker=marker,
            source_alleles=source_alleles,
            random_seed=random_seed,
            two_step_mutation_factor=two_step_mutation_factor,
        )


    return target_haplotype


def serialize_haplotype(
        haplotype: Haplotype
) -> tuple:
    """
        Serializes a haplotype object into a tuple format for efficient caching

        Args:
            haplotype (Haplotype): The haplotype to be serialized.

        Returns:
            tuple: A tuple representation of the haplotype, where each marker's alleles are sorted
                   by their values and intermediate values.
        """
    return tuple(
        (marker_name, tuple(sorted((allele.value, allele.intermediate_value) for allele in alleles)))
        for marker_name, alleles in sorted(haplotype.alleles.items())
    )


def get_edge_probability(
        source: Haplotype,
        target: Haplotype,
        marker_set: MarkerSet,
        two_step_mutation_factor: float,
        is_average_pedigree: bool
) -> Decimal:
    """
        Computes the probability of mutating from a source haplotype to a target haplotype
        based on a given set of markers.

        Args:
            source (Haplotype): The original haplotype before mutation.
            target (Haplotype): The haplotype after potential mutation.
            marker_set (MarkerSet): A set of markers with corresponding mutation rates.
            two_step_mutation_factor (float): A factor that determines the probability of a two-step mutation. Defaults to 0.03.

        Returns:
            Decimal: The cumulative probability of mutation across all markers.

        Note:
            - The alleles for each marker are sorted based on their values before comparison.
            - The mutation probability for each marker is determined using `calculate_mutation_probability`.
            - The final edge probability is the product of all individual mutation probabilities.
        """

    edge_probability = Decimal(1)

    for marker in marker_set.markers:
        source_alleles = sorted(source.alleles[marker.name], key=lambda allele: allele.value)
        target_alleles = sorted(target.alleles[marker.name], key=lambda allele: allele.value)

        mutation_probability = calculate_mutation_probability(
            parent_alleles=source_alleles,
            child_alleles=target_alleles,
            marker=marker,
            two_step_mutation_factor=two_step_mutation_factor,
            is_average_pedigree=is_average_pedigree,
        )

        edge_probability *= mutation_probability

    return edge_probability


def get_edge_probabilities(
        haplotypes: Mapping[str, Haplotype],
        relationships: Collection[Relationship],
        marker_set: MarkerSet,
        two_step_mutation_factor: float,
        is_average_pedigree: bool,
) -> Mapping[tuple[str, str], Decimal]:
    """
        Computes mutation probabilities for all parent-child haplotype pairs in a given set of relationships.

        Args:
            haplotypes (Mapping[int, Haplotype]): A mapping of individual IDs to their respective haplotypes.
            relationships (Collection[Relationship]): A collection of relationships defining parent-child pairs.
            marker_set (MarkerSet): A set of markers with corresponding mutation rates.
            two_step_mutation_factor (float): A factor that determines the probability of a two-step mutation. Defaults to 0.03.

        Returns:
            Mapping[tuple[int, int], Decimal]: A dictionary where keys are (parent_id, child_id) tuples
            and values are the corresponding mutation probabilities.

        Note:
            - Each relationship defines a parent-child connection.
            - The function calculates the probability of mutating from the parent's haplotype to the child's haplotype
              using `get_edge_probability`.
            - The final result is a dictionary mapping each parent-child pair to their computed mutation probability.
        """

    return {
        (relationship.parent_id, relationship.child_id): get_edge_probability(
            source=haplotypes[relationship.parent_id],
            target=haplotypes[relationship.child_id],
            marker_set=marker_set,
            two_step_mutation_factor=two_step_mutation_factor,
            is_average_pedigree=is_average_pedigree,
        )
        for relationship in relationships
    }


def update_average(
        new_probability: Decimal,
        old_probability: Decimal | None,
        i: int
) -> Decimal:
    """
        Updates the running average of a probability value using an incremental approach.

        Args:
            new_probability (Decimal): The newly computed probability to incorporate.
            old_probability (Decimal): The current running average probability.
            i (int): The number of values previously averaged.

        Returns:
            Decimal: The updated average probability.

        Note:
            - This function implements an incremental mean formula:
              updated_avg = ((old_avg * count) + new_value) / (count + 1).
            - It ensures numerical stability when updating an average iteratively.
        """
    if old_probability is None:
        return new_probability
    return ((old_probability * i) + new_probability) / (
            i + 1
    )


def is_stable(
        probabilities: list[Decimal],
        threshold: float,
        reporter: Reporter,
) -> bool:
    """
        Determines whether a sequence of probabilities has stabilized based on recent changes.

        Args:
            probabilities (list[Decimal]): A list of probability values recorded over iterations.
            threshold (float, optional): The maximum allowed relative change between consecutive probabilities
                                         for stability. Default is 0.0001.
            reporter (Reporter): A reporting tool for logging and tracking progress.

        Returns:
            bool: True if the probability values have stabilized, False otherwise.

        Note:
            - Stability is checked only if `i` exceeds both `window` and `min_iterations`.
            - Stability is determined by ensuring that the relative change in probability
              over the last `window` iterations remains below `threshold`.
            - If the previous probability is zero, stability is achieved only if the current probability is also zero.
        """

    min_value = float(min(probabilities))
    max_value = float(max(probabilities))
    if min_value > 0:
        if (max_value - min_value) / min_value < threshold:
            return True
    if min_value == 0 and max_value == 0:
        reporter.log("\nProbabilities are all zero. Pedigree is impossible.")
        return True
    return False


def is_model_valid(
        probabilities: list[Decimal],
        threshold: float = 0.005
) -> tuple[bool, Decimal, list[Decimal], float]:
    """
        Checks whether three independent simulation runs produce consistent results.

        Args:
            probabilities (list[Decimal]): A list of final probability values from three separate simulation runs.
            threshold (float, optional): The maximum allowed relative difference between results
                                         for the model to be considered valid. Default is 0.005.

        Returns:
            bool: True if all three runs produce similar probabilities within the threshold, False otherwise.

        Note:
            - The function assumes `probabilities` contains exactly three values, each corresponding to a separate run.
            - The model is considered valid if the relative difference between consecutive results
              is within the specified `threshold`.
            - If a previous probability is zero, validity is ensured only if the current probability is also zero.
        """
    # check if there is a combination of three probabilities that are within the threshold
    if len(probabilities) < 3:
        raise ValueError("At least three probabilities are required to validate the model.")

    combs = combinations(probabilities, 3)

    lowest_comb = None
    lowest_variance_sum = None
    lowest_mean = None

    for comb in combs:
        mean_probability = Decimal(sum(comb) / len(comb))
        abs_deviation_sum = sum([abs(comb[i] - mean_probability) for i in range(len(comb))])
        model_is_valid = all(
        abs((comb[i] - mean_probability) / mean_probability) < threshold
        if mean_probability != 0 else comb[i] == 0
        for i in range(len(comb))
        )

        if model_is_valid:
            if lowest_variance_sum is None or abs_deviation_sum < lowest_variance_sum:
                lowest_variance_sum = abs_deviation_sum
                lowest_comb = comb
                lowest_mean = mean_probability

    if lowest_comb is not None and lowest_mean is not None:
        if lowest_mean > 1.0:
            return False, Decimal(sum(probabilities) / len(probabilities)), [], threshold / 2
        return True, Decimal(lowest_mean), list(lowest_comb), threshold

    else:
        return False, Decimal(sum(probabilities) / len(probabilities)), [], threshold


def simulate_pedigree_probability(
        individuals: Mapping[str, Individual],
        relationships: Collection[Relationship],
        parents: Mapping[str, str],
        ordered_unknown_ids: Collection[str],
        marker_set: MarkerSet,
        random_seed: int,
        two_step_mutation_factor: float
) -> IterationResult:
    """
        Simulates the probabilities of a pedigree based on the provided relationships and markers.

        Args:
            individuals (Mapping[int, Individual]): A dictionary mapping individual IDs to their respective `Individual` objects.
            relationships (Collection[Relationship]): A collection of `Relationship` objects representing parent-child relationships.
            parents (Mapping[int, int]): A mapping of individual IDs to their parent IDs.
            ordered_unknown_ids (Collection[int]): A collection of individual IDs whose haplotypes are unknown and need to be simulated.
            marker_set (MarkerSet): A set of markers used for allele probability calculations.
            random_seed (Random): A random number generator used to simulate the mutation of haplotypes.
            two_step_mutation_factor (float): A factor that determines the probability of a two-step mutation. Defaults to 0.03.

        Returns:
            IterationResult: An object containing the calculated pedigree probability, and the edge probabilities for the relationships.

        Note:
            - The function simulates the haplotypes of individuals with unknown haplotypes by mutating their parent's haplotype.
            - It calculates the edge probabilities for parent-child relationships and uses them to compute the overall pedigree probability.
            - If any edge probability is zero, the pedigree probability is set to zero to avoid underflow.
        """

    random = SystemRandom()

    simulated_individual_ids: set[str] = set()
    simulated_relationship_ids: set[tuple[str, str]] = set()

    haplotypes = {
        individual.id: individual.haplotype for individual in individuals.values()
    }

    for individual_id in ordered_unknown_ids:
        parent_id = parents[individual_id]
        parent_haplotype = haplotypes[parent_id]

        # Use the parents haplotype to predict the unknown individuals haplotype
        haplotypes[individual_id] = mutate_haplotype(
            source=parent_haplotype,
            marker_set=marker_set,
            random_seed=random_seed,
            two_step_mutation_factor=two_step_mutation_factor,
        )

        simulated_individual_ids.add(individual_id)
        simulated_relationship_ids.add((parent_id, individual_id))

    # distances = set(sorted([individuals[ind].closest_known_distance for ind in ordered_unknown_ids]))
    #
    # for distance in distances:
    #     unknown_ids_with_distance = [ind for ind in ordered_unknown_ids if individuals.get(ind).closest_known_distance == distance]
    #     for individual_id in unknown_ids_with_distance:
    #         individual = individuals.get(individual_id)
    #         closest_known_individuals = individual.closest_known_individuals
    #
    #         selected_closest_known_individual = random.choice(closest_known_individuals)
    #         selected_closest_known_haplotype = haplotypes[selected_closest_known_individual.id]
    #
    #         # Use the closest known individual's haplotype to predict the unknown individual's haplotype
    #         haplotypes[individual_id] = mutate_haplotype(
    #             source=selected_closest_known_haplotype,
    #             marker_set=marker_set,
    #             random_seed=random_seed,
    #             two_step_mutation_factor=two_step_mutation_factor,
    #         )
    #
    #         simulated_individual_ids.add(individual_id)
    #         simulated_relationship_ids.add((selected_closest_known_individual.id, individual_id))

    unused_relationships = [
        relationship
        for relationship in relationships
        if (relationship.parent_id, relationship.child_id)
           not in simulated_relationship_ids and (relationship.child_id, relationship.parent_id)
           not in simulated_relationship_ids
    ]

    edge_probabilities = get_edge_probabilities(
        haplotypes, unused_relationships, marker_set, two_step_mutation_factor, True
    )

    if any(probability == 0 for probability in edge_probabilities.values()):
        pedigree_probability = Decimal(0)
    else:
        pedigree_probability = reduce(operator.mul, edge_probabilities.values(), 1)

    return IterationResult(
        probability=Decimal(pedigree_probability),
        edge_probabilities=edge_probabilities,
        mutated_haplotypes=haplotypes,
        fixed_individual_id=None,
    )


def simulate_pedigree_iteration(
        i: int,
        individuals: Mapping[str, Individual],
        relationships: Collection[Relationship],
        parents: Mapping[str, str],
        ordered_unknown_ids: Sequence[str],
        marker_set: MarkerSet,
        random_seed: int,
        two_step_mutation_factor: float
) -> IterationResult:
    return simulate_pedigree_probability(
        individuals=individuals,
        relationships=relationships,
        parents=parents,
        ordered_unknown_ids=ordered_unknown_ids,
        marker_set=marker_set,
        random_seed=random_seed + i,  # Ensure each process has a unique random seed
        two_step_mutation_factor=two_step_mutation_factor,
    )


def process_iteration_results(
        simulate_func: callable,
        simulation_parameters: SimulationParameters,
        progress_bar: ProgressBar,
        reporter: Reporter,
        out_file_name: str,
        number_of_threads: int = 1,
        is_outside: bool = False
):
    with (progress_bar):
        valid = False
        trial = 1
        tightening = 0

        model_iterations = {m: 0 for m in range(3)}
        model_probabilities = {m: [] for m in range(3)}

        while not valid:
            for m in range(3):  # Model validation
                with open(
                        f"{simulation_parameters.results_path}/{out_file_name}_m_{m}_outside_{is_outside}.txt",
                        "a", 1) as f:
                    with multiprocessing.Pool(number_of_threads) as pool:
                        print(f"\nStarting trial {trial}, model {m + 1} with {number_of_threads} threads...")
                        #chunksize = max(1, simulation_parameters.stability_window // (number_of_threads * 40))
                        chunk_size = 25
                        iteration_results = pool.imap_unordered(simulate_func,
                                                                range(0, simulation_parameters.stability_window),
                                                                chunksize=chunk_size)
                        for i, iteration_result in enumerate(iteration_results):
                            model_iterations[m] += 1
                            if len(model_probabilities[m]) == 0:
                                old_probability = None
                            else:
                                old_probability = model_probabilities[m][-1]

                            model_new_mean = update_average(
                                old_probability=old_probability,
                                new_probability=iteration_result.probability,
                                i=model_iterations[m]
                            )
                            model_probabilities[m].append(model_new_mean)

                            if i % 100 == 0:
                                f.write(f"{model_new_mean}\n")

                            progress_bar.update(1)

            total_model_probability_mean = Decimal(
                sum(model_probabilities[m][-1] for m in range(3)) / 3
            )

            # all values in the model_probabilities should be within simulation_parameters.model_validity_threshold of the total_model_probability_mean
            number_of_probabilities_outside_threshold = sum(
                sum(
                    abs(prob - total_model_probability_mean) / total_model_probability_mean > simulation_parameters.model_validity_threshold
                    for prob in model_probabilities[m])
                for m in range(3)
            )

            if all(
                    all(
                        abs(prob - total_model_probability_mean) / total_model_probability_mean < simulation_parameters.model_validity_threshold
                        for prob in model_probabilities[m])
                    for m in range(3)
            ):
                if total_model_probability_mean > 1.0:
                    reporter.log(
                        f"\nModel is stable after {trial} trials. However, current mean: {total_model_probability_mean} (log {10 * math.log10(total_model_probability_mean)}) is greater than unity. Starting new trial with tightened threshold: {simulation_parameters.model_validity_threshold / 2}"
                    )

                    if out_file_name == "match_probabilities" and tightening == 2:
                        reporter.log(
                            "\nMaximum number of tightening reached. Model is not valid. Restarting simulation."
                        )
                        raise InvalidAveragePedigreeProbability("Average pedigree probability likely not valid.")

                    simulation_parameters.model_validity_threshold /= 1.5  # If the mean probability is greater than 1, we need to tighten the threshold
                    tightening += 1

                else:
                    valid = True
                    reporter.log(
                        f"\nModel is valid after {trial} trials! Current mean: {total_model_probability_mean} (log {10 * math.log10(total_model_probability_mean)})")
            else:
                reporter.log(
                    f"\nModel is not valid after {trial} trials. Current mean: {total_model_probability_mean} (log {10 * math.log10(total_model_probability_mean)}). {number_of_probabilities_outside_threshold} data points fall outside range. Starting new trial...")
                trial += 1
                model_probabilities = {m: model_probabilities[m][-1:] for m in
                                       range(3)}  # Keep only the last probability of each model

    return total_model_probability_mean, [int(trial * simulation_parameters.stability_window * 3)], [
        Decimal(model_probabilities[m][-1]) for m in range(3)], True, [Decimal(model_probabilities[m][-1]) for m in
                                                                       range(3)]


def calculate_average_pedigree_probability(
        pedigree: Pedigree,
        root_name: str,
        suspect_haplotype: Haplotype,
        marker_set: MarkerSet,
        simulation_parameters: SimulationParameters,
        random_seed: int,
        reporter: Reporter,
        is_outside: bool,
        number_of_threads: int = 1,
) -> tuple[Decimal, list[int], list[Decimal], bool, list[Decimal]]:
    """
    Calculates the average probability of a given pedigree using a Monte Carlo simulation.

    This function estimates the probability of the pedigree, denoted as P(Hv), by simulating the inheritance
    of haplotypes through the pedigree over a given number of iterations. The suspect is set as the root of
    the pedigree to facilitate a level order traversal, ensuring that unknown individuals are processed in
    a structured manner.

    ### Steps:
    1. The pedigree is traversed in level order, ensuring unknown individuals are assigned haplotypes
       based on the haplotypes of surrounding known individuals and mutation rates.
    2. Each unknown individual's haplotype is predicted using the haplotype of its parent.
    3. After all unknown individuals have been processed, the pedigree probability is calculated by
       multiplying the probabilities of the edges that were not used to determine haplotypes.
    4. The edge probabilities are derived from allele probabilities, which depend on mutation rates.
    5. The process is repeated for a set number of iterations, with an average pedigree probability
       computed over all runs.
    6. The function stops early if the probability estimate stabilizes.

    Args:
        pedigree (Pedigree): The pedigree structure containing individuals and relationships.
        root_name (str): The name of the suspect, used as the root of the pedigree.
        suspect_haplotype (Haplotype): The haplotype of the suspect individual, used as a reference for predictions.
        marker_set (MarkerSet): A set of genetic markers used in haplotype prediction.
        simulation_parameters (Mapping[str, any]): A dictionary containing simulation parameters such as
        random_seed (Random): A random number generator for controlled stochastic processes.
        reporter (Reporter): A reporting tool for logging and tracking progress.
        is_outside (bool): A flag indicating whether the simulation is performed on the extended pedigree.
        number_of_threads (int): The number of threads to use for parallel processing. Default is 1.

    Returns:
        Decimal: The estimated average pedigree probability P(Hv).
        list[int]: A list of the number of iterations needed for stability across three models.
        list[Decimal]: A list of model probabilities from three independent runs.
        bool: A boolean indicating whether the model is valid based on the stability of the probabilities.

    Notes:
        - The function maintains a running average of pedigree probabilities to detect stabilization.
        - If the probability estimate stabilizes before reaching `number_of_iterations`, the function
          terminates early.
        - Edge probabilities are logged, showing how likely specific parent-child haplotype transmissions are.
        - The final probability is reported both as a raw value and a log-transformed value.
    """

    progress_bar = reporter.progress_bar(
        total=simulation_parameters.max_number_of_iterations * 3, desc="Calculating average pedigree probability"
    )

    individuals = {individual.id: individual for individual in pedigree.individuals}

    child_parent_dict = {
        relationship.child_id: relationship.parent_id
        for relationship in pedigree.relationships
    }

    ordered_unknown_ids = [
        individual.id
        for individual in pedigree.get_level_order_traversal(root_name)
        if individual.haplotype_class == "unknown"
    ]

    known_individuals = [
        individual for individual in pedigree.get_level_order_traversal(root_name)
        if individual.haplotype_class != "unknown"
    ]

    if len(known_individuals) == 1:
        # If only the suspect is known, return a probability of 1
        reporter.log("\nRoot is only known individual in (reduced) pedigree. Average pedigree probability is 1.")
        return Decimal(1), [0,0,0], [Decimal(1),Decimal(1),Decimal(1)], True, [Decimal(1)]  # If only the suspect is known, return a probability of 1

    simulate_func = partial(
        simulate_pedigree_iteration,
        individuals=individuals,
        relationships=pedigree.relationships,
        parents=child_parent_dict,
        ordered_unknown_ids=ordered_unknown_ids,
        marker_set=marker_set,
        random_seed=random_seed,
        two_step_mutation_factor=simulation_parameters.two_step_mutation_factor,
    )

    return process_iteration_results(
        simulate_func=simulate_func,
        simulation_parameters=simulation_parameters,
        progress_bar=progress_bar,
        reporter=reporter,
        out_file_name="average_pedigree_probabilities",
        number_of_threads=number_of_threads,
        is_outside=is_outside,
    )


def simulate_matching_haplotypes(
        individuals: Mapping[str, Individual],
        relationships: Collection[Relationship],
        parents: Mapping[str, str],
        suspect_haplotype: Haplotype,
        ordered_unknown_ids: Sequence[str],
        marker_set: MarkerSet,
        average_pedigree_probability: Decimal,
        random_seed: int,
        two_step_mutation_factor: float,
        picking_probabilities: dict[str, Decimal],
        is_outside: bool,
) -> IterationResult:
    """
    Simulates the inheritance of haplotypes in a pedigree and computes the probability of obtaining
    exactly `l` matching haplotypes.

    This function performs a stochastic simulation to determine the likelihood of `l` unknown individuals
    inheriting the same haplotype as the suspect. It follows these steps:

    ### Steps:
    1. **Selection of `l` Individuals:**
       - Selects `l` unknown individuals at random, with a probability based on predefined weights.
       - These individuals are assigned the same haplotype as the suspect.

    2. **Haplotype Prediction:**
       - Remaining unknown individuals have their haplotypes determined by simulating inheritance
         from their parents, incorporating mutation probabilities.

    3. **Probability Computation:**
       - Computes the probability of the specific haplotype configuration (`Pr(H = h)`).
       - Derives the conditional probability (`Pr(Huk = huk | Hkn = hkn)`) by normalizing against the
         average pedigree probability.
       - Ensures the final configuration contains exactly `l` matching haplotypes (excluding the suspect).

    4. **Result Compilation:**
       - If requested, generates a simulated pedigree showing assigned haplotypes.
       - Returns the probability of achieving `l` matching haplotypes along with edge probabilities.

    ### Args:
        individuals (Mapping[int, Individual]): A dictionary mapping individual IDs to `Individual` objects.
        relationships (Collection[Relationship]): A collection of parent-child relationships in the pedigree.
        parents (Mapping[int, int]): A dictionary mapping child IDs to their respective parent IDs.
        suspect (Individual): The suspect whose haplotype is being matched.
        ordered_unknown_ids (Sequence[int]): A sequence of unknown individual IDs in level-order traversal.
        marker_set (MarkerSet): A set of genetic markers used in haplotype inheritance simulation.
        average_pedigree_probability (Decimal): The average probability of the pedigree, used for normalization.
        random_seed (Random): A random number generator to ensure controlled stochastic behavior.

    ### Returns:
        IterationResult: A result object containing:
        - The computed probability of obtaining exactly `l` matching haplotypes.
        - Edge probabilities indicating haplotype transmission likelihoods.

    ### Notes:
        - If `average_pedigree_probability` is zero, the conditional probability is set to zero.
        - The selection of `l` individuals follows weighted probabilities to reflect real-world inheritance patterns.
        - Mutation rates influence haplotype assignments for non-fixed individuals.
        - The simulation probability accounts for both individual selection and haplotype propagation.
        - The function ensures that the final result adheres to `l` matches exactly before returning a probability.
    """
    # TODO: add correct random seed handling
    random = SystemRandom()

    simulated_individual_ids: set[str] = set()
    simulated_relationship_ids: set[tuple[str, str]] = set()
    fixed_individual_id: str

    simulation_probability = Decimal(1)

    normalized_picking_probabilities = {
        individual: probability / sum(picking_probabilities.values()) for individual, probability in
        picking_probabilities.items()
    }
    # pick a combination of individuals based on the normalized picking probabilities
    picked_individual = random.choices(
        population=list(normalized_picking_probabilities.keys()),
        weights=[float(normalized_picking_probabilities[ind]) for ind in picking_probabilities.keys()],
    )[0]
    # add the picked individual to the fixed individuals
    fixed_individual_id = picked_individual
    # get the picking probability of the picked individual
    picking_probability = normalized_picking_probabilities[picked_individual]
    simulation_probability *= picking_probability

    haplotypes = {
        individual.id: (
            deepcopy(suspect_haplotype)
            if individual.id == fixed_individual_id
            else deepcopy(individual.haplotype)
        )
        for individual in individuals.values()
    }

    for i, individual_id in enumerate(ordered_unknown_ids):
        if individual_id == fixed_individual_id:
            continue

        parent_id = parents[individual_id]
        parent_haplotype = haplotypes[parent_id]

        haplotypes[individual_id] = deepcopy(mutate_haplotype(
            source=parent_haplotype,
            marker_set=marker_set,
            random_seed=random_seed+i,
            two_step_mutation_factor=two_step_mutation_factor,
        ))

        simulated_individual_ids.add(individual_id)
        simulated_relationship_ids.add((parent_id, individual_id))

    # distances = set(sorted([individuals[ind].closest_known_distance for ind in ordered_unknown_ids]))
    #
    # for distance in distances:
    #     unknown_ids_with_distance = [ind for ind in ordered_unknown_ids if individuals.get(ind).closest_known_distance == distance]
    #
    #     for individual_id in unknown_ids_with_distance:
    #         if individual_id == fixed_individual_id:
    #             continue
    #
    #         individual = individuals.get(individual_id)
    #         closest_known_individuals = individual.closest_known_individuals
    #
    #         selected_closest_known_individual = random.choice(closest_known_individuals)
    #         selected_closest_known_haplotype = haplotypes[selected_closest_known_individual.id]
    #
    #         # Use the closest known individual's haplotype to predict the unknown individual's haplotype
    #         haplotypes[individual_id] = mutate_haplotype(
    #             source=selected_closest_known_haplotype,
    #             marker_set=marker_set,
    #             random_seed=random_seed,
    #             two_step_mutation_factor=two_step_mutation_factor,
    #         )
    #
    #         simulated_individual_ids.add(individual_id)
    #         simulated_relationship_ids.add((selected_closest_known_individual.id, individual_id))

    # Calculate the probability of the entire pedigree (Pr(H = h))
    all_edge_probabilities = get_edge_probabilities(
        haplotypes, relationships, marker_set, two_step_mutation_factor, False
    )

    pedigree_probability = reduce(operator.mul, all_edge_probabilities.values(), 1)

    # Calculate the probability of the simulated edges
    simulated_edge_probabilities = []
    for (source, target) in simulated_relationship_ids:
        if (source, target) in all_edge_probabilities:
            simulated_edge_probabilities.append(all_edge_probabilities[(source, target)])
        elif (target, source) in all_edge_probabilities:
            simulated_edge_probabilities.append(all_edge_probabilities[(target, source)])

    simulation_probability = reduce(
        operator.mul, simulated_edge_probabilities, simulation_probability
    )

    # Calculate (Pr(Huk = huk|Hkn = hkn))
    if average_pedigree_probability == 0:
        conditional_probability = Decimal(0)
    else:
        conditional_probability = pedigree_probability / average_pedigree_probability

    # Check if the total number of matching haplotypes is equal to l (excluding the suspect and excluded individuals)
    non_excluded_matching_ids = [
        individual_id for individual_id in ordered_unknown_ids
        if haplotypes[individual_id] == suspect_haplotype and not individuals[individual_id].exclude
    ]

    all_matching_ids = [
        individual_id for individual_id in ordered_unknown_ids
        if haplotypes[individual_id] == suspect_haplotype
    ]

    total_matching_ids = {fixed_individual_id}.union(set(all_matching_ids))
    total_number_of_matching_haplotypes = len(total_matching_ids)

    number_of_non_excluded_matching_haplotypes = len(non_excluded_matching_ids)

    # number of matching has to be at least 1
    # if outside: S* is always matching
    if not is_outside:
        probability = (
           conditional_probability / (simulation_probability * total_number_of_matching_haplotypes)
            if number_of_non_excluded_matching_haplotypes >= 1 and simulation_probability != 0
            else Decimal(0)
        )
    else:
        probability = (
            conditional_probability / simulation_probability
            if simulation_probability != 0
            else Decimal(0)
        )

    return IterationResult(
        probability=probability,
        edge_probabilities=all_edge_probabilities,
        mutated_haplotypes=haplotypes,
        fixed_individual_id=fixed_individual_id,
    )


def simulate_iteration(
        i: int,
        individuals: Mapping[str, Individual],
        relationships: Collection[Relationship],
        parents: Mapping[str, str],
        suspect_haplotype: Haplotype,
        ordered_unknown_ids: Sequence[str],
        marker_set: MarkerSet,
        average_pedigree_probability: Decimal,
        random_seed: int,
        two_step_mutation_factor: float,
        picking_probabilities: dict[str, Decimal],
        is_outside: bool,
) -> IterationResult:
    """Wrapper to allow multiprocessing without lambda."""
    return simulate_matching_haplotypes(
        individuals=individuals,
        relationships=relationships,
        parents=parents,
        suspect_haplotype=suspect_haplotype,
        ordered_unknown_ids=ordered_unknown_ids,
        marker_set=marker_set,
        average_pedigree_probability=average_pedigree_probability,
        random_seed=random_seed + i,  # Ensure each process has a unique random seed
        two_step_mutation_factor=two_step_mutation_factor,
        picking_probabilities=picking_probabilities,
        is_outside=is_outside,
    )


def calculate_matching_haplotypes(
        pedigree: Pedigree,
        marker_set: MarkerSet,
        root_name: str,
        suspect_haplotype: Haplotype,
        simulation_parameters: SimulationParameters,
        average_pedigree_probability: Decimal,
        reporter: Reporter,
        is_outside: bool,
        random_seed: int = 42,
        number_of_threads: int = 1,
) -> tuple[Decimal, list[int], list[Decimal], bool, list[Decimal]]:
    """
    Simulates the pedigree to estimate the probability of obtaining exactly `l` matching haplotypes
    between unknown individuals and the suspect.

    This function follows a Monte-Carlo simulation approach with importance sampling. The process involves:

    1. Selecting `l` unknown individuals and setting their haplotypes to match the suspect's haplotype.
    2. Predicting the haplotypes of the remaining unknown individuals based on known haplotypes and mutation rates.
    3. Iteratively updating the simulation probability and pedigree probability over multiple iterations.
    4. Validating the stability of the simulation results across three independent runs.

    ### Simulation Steps:
    - The simulation tracks the probability of the specific iteration (`Pr(H̃uk = huk | Hkn = hkn)`).
    - The haplotypes of `l` unknown individuals are forcibly set to match the suspect's haplotype.
    - The haplotypes of the remaining unknown individuals are inferred based on surrounding known individuals.
    - The simulation probability is averaged across iterations to estimate the final probability.

    ### Stability and Model Validation:
    - The simulation runs three independent times to ensure the results are stable.
    - If all three model probabilities are within 0.5% variation, the model is considered valid.
    - If stability is not reached, increasing `number_of_iterations` is recommended.

    Args:
        pedigree (Pedigree): The pedigree structure containing individuals and relationships.
        marker_set (MarkerSet): A set of genetic markers used for allele probability calculations.
        root_name (str): The name of the suspect whose haplotype is used as a reference.
        suspect_haplotype (Haplotype): The haplotype of the suspect individual, used for matching.
        simulation_parameters (Mapping[str, any]): A dictionary containing simulation parameters such as
        average_pedigree_probability (Decimal): The precomputed average probability of the pedigree.
        random_seed (int): A seed for the random number generator to ensure reproducibility.
        reporter (Reporter): A reporting tool for tracking progress and logging results.
        is_outside (bool): A flag indicating whether the simulation is performed on the extended pedigree.
        number_of_threads (int): The number of threads to use for parallel processing. Default is 1.

    Returns:
        Decimal: The estimated probability of obtaining exactly `l` matching haplotypes.

    Notes:
        - A progress bar is displayed to track iterations.
        - Edge probabilities are computed and logged to analyze relationship strength.
        - If the model is deemed invalid, additional iterations should be performed.
    """

    progress_bar = reporter.progress_bar(
        total=simulation_parameters.max_number_of_iterations * 3, desc=f"Calculating matching haplotypes"
    )

    individuals = {individual.id: individual for individual in pedigree.individuals}

    parents = {
        relationship.child_id: relationship.parent_id
        for relationship in pedigree.relationships
    }

    ordered_unknown_ids = [
        individual.id
        for individual in pedigree.get_level_order_traversal(root_name)
        if individual.haplotype_class == "unknown"
    ]

    simulate_func = partial(
        simulate_iteration,
        individuals=individuals,
        relationships=pedigree.relationships,
        parents=parents,
        suspect_haplotype=suspect_haplotype,
        ordered_unknown_ids=ordered_unknown_ids,
        marker_set=marker_set,
        average_pedigree_probability=average_pedigree_probability,
        random_seed=random_seed,
        two_step_mutation_factor=simulation_parameters.two_step_mutation_factor,
        picking_probabilities=pedigree.picking_probabilities,
        is_outside=is_outside,
    )

    return process_iteration_results(
        simulate_func=simulate_func,
        simulation_parameters=simulation_parameters,
        progress_bar=progress_bar,
        reporter=reporter,
        out_file_name="match_probabilities",
        number_of_threads=number_of_threads,
        is_outside=is_outside
    )


def calculate_proposal_distribution(
        pedigree: Pedigree,
        marker_set: MarkerSet,
        root_name: str,
        suspect_haplotype: Haplotype,
        average_pedigree_probability: Decimal,
        simulation_parameters: SimulationParameters,
        random_seed: int,
        reporter: Reporter,
) -> tuple[Mapping[int, Decimal], list[int], list[Decimal], bool, list[Decimal]]:
    """
        Calculates the proposal distribution for the number of matching haplotypes between unknown individuals
        in a pedigree and a given suspect's haplotype. This distribution is used in importance sampling for
        Monte-Carlo simulations.

        The function iterates through the possible values of `l`, representing the number of matching haplotypes,
        and for each value, it calculates the corresponding probability using a Monte-Carlo simulation. It then
        computes the proposal distribution, which represents the likelihood of obtaining `l` matching haplotypes
        from the total number of unknown individuals in the pedigree.

        ### Steps:
        1. For each possible value of `l` (number of matching haplotypes), the probability is calculated using
           `calculate_l_matching_haplotypes()`.
        2. The probability for `l = 0` (no matching haplotypes) is calculated as the remaining probability
           (`1 - sum(proposal_distribution.values())`).
        3. If any individuals are marked as excluded, they are not considered in the calculation.

        Args:
            pedigree (Pedigree): The pedigree structure containing individuals and relationships.
            marker_set (MarkerSet): A set of genetic markers used for allele probability calculations.
            root_name (str): The name of the suspect whose haplotype is used as a reference.
            suspect_haplotype (Haplotype): The haplotype of the suspect individual, used for matching.
            average_pedigree_probability (Decimal): The precomputed average probability of the pedigree.
            simulation_parameters (Mapping[str, any]): A dictionary containing simulation parameters such as
            random_seed (int): A seed for the random number generator to ensure reproducibility.
            reporter (Reporter): A reporting tool for tracking progress and logging results.

        Returns:
            Mapping[int, Decimal]: A dictionary where the keys are integers representing the number of matching
            haplotypes (`l`), and the values are the corresponding probabilities.

        Notes:
            - The function assumes that the total number of unknown individuals is greater than or equal to `l`.
            - If the computed probability for `l = 0` is negative, a warning is logged indicating that the proposal
              distribution does not add up to unity.
            - Individuals marked with `exclude=True` are ignored in the calculation.
        """

    unknown_individuals = pedigree.get_unknown_individuals()
    proposal_distribution: dict[int, Decimal] = {}
    needed_iterations: list[int] = []
    model_probabilities: list[Decimal] = []
    model_validity: bool = False
    used_probabilities: list[Decimal] = []
    number_of_non_excluded_unknowns = len([ind for ind in unknown_individuals if not ind.exclude])

    # There has to be at least one non-excluded unknown individual in the pedigree
    if number_of_non_excluded_unknowns == 0:
        reporter.log("\nWarning! No (non-excluded) unknown individuals left in the pedigree. Only calculating outside match probability.")
        proposal_distribution[0] = Decimal(1)
        proposal_distribution[1] = Decimal(0)
        return proposal_distribution, needed_iterations, model_probabilities, model_validity, used_probabilities

    max_number_of_threads = multiprocessing.cpu_count()
    number_of_threads = min(simulation_parameters.number_of_threads, max_number_of_threads)

    proposal_distribution[1], needed_iterations, model_probabilities, model_validity, used_probabilities = calculate_matching_haplotypes(
        pedigree=pedigree,
        marker_set=marker_set,
        root_name=root_name,
        suspect_haplotype=suspect_haplotype,
        average_pedigree_probability=average_pedigree_probability,
        simulation_parameters=simulation_parameters,
        reporter=reporter,
        is_outside=False,
        random_seed=random_seed,
        number_of_threads=number_of_threads,
    )

    # The non-match probability (i.e. the probability of 0 matching individuals) is 1 - the probability of at least 1 matching individual
    proposal_distribution[0] = Decimal(1) - proposal_distribution[1]
    if proposal_distribution[0] < 0:
        reporter.log("\nWarning! Proposal distribution does not add up to unity.")
    return proposal_distribution, needed_iterations, model_probabilities, model_validity, used_probabilities


def calculate_outside_match_probability(
        pedigree: Pedigree,
        marker_set: MarkerSet,
        root_name: str,
        suspect_haplotype: Haplotype,
        average_pedigree_probability: Decimal,
        simulation_parameters: SimulationParameters,
        random_seed: int,
        reporter: Reporter,
        number_of_threads: int = 1,
)   -> tuple[Decimal, list[int], list[Decimal], bool, list[Decimal]]:
    outside_match_probability, needed_iterations, model_probabilities, model_is_valid, used_probabilities = calculate_matching_haplotypes(
        pedigree=pedigree,
        marker_set=marker_set,
        root_name=root_name,
        suspect_haplotype=suspect_haplotype,
        average_pedigree_probability=average_pedigree_probability,
        simulation_parameters=simulation_parameters,
        random_seed=random_seed,
        reporter=reporter,
        is_outside=True,
        number_of_threads=number_of_threads,
    )
    return outside_match_probability, needed_iterations, model_probabilities, model_is_valid, used_probabilities


def run_simulation(
        input_pedigree: Pedigree,
        suspect_name: str,
        marker_set: MarkerSet,
        simulation_parameters: SimulationParameters,
        random_seed: int,
        reporter: Reporter,
        skip_inside: bool = False,
        skip_outside: bool = False,
) -> SimulationResult:
    """
        Runs a Monte-Carlo simulation with Importance Sampling to calculate match probabilities
        between a suspect and individuals in a pedigree.

        The goal is to compute the probability of having 0, 1, 2, ... matching haplotypes with the suspect.
        The simulation is performed using a proposal distribution where the probability of observing the
        suspect's haplotype is high, and the number of matching haplotypes is set to a fixed number.

        The simulation is divided into three main steps:
        1. **Calculate the average pedigree probability (P(hv))**:
           This step computes the average probability of the pedigree based on the number of iterations.
        2. **Calculate the proposal distribution**:
           The proposal distribution calculates the number of matching haplotypes in the pedigree relative to the suspect.
        3. **Calculate the outside match probability**:
           The pedigree is extended, and an outside match probability is calculated.

        Args:
            input_pedigree (Pedigree): The pedigree object containing the family structure and genetic information.
            suspect_name (str): The name of the suspect whose haplotype is being compared to others.
            marker_set (MarkerSet): A set of genetic markers used for calculating allele probabilities.
            simulation_parameters (Mapping[str, float]): A dictionary containing simulation parameters such as
                - `number_of_iterations`: The number of iterations for the Monte-Carlo simulation.
            random_seed (Random): A random number generator used for random selection and mutation.
            reporter (Reporter): A reporter object used to track the progress and output of the simulation.
            skip_inside (bool): If True, skips the inside match probability calculation.
            skip_outside (bool): If True, skips the outside match probability calculation.

        Returns:
            SimulationResult: An object containing the following:
                - `average_pedigree_probability`: The average probability of the original pedigree.
                - `proposal_distribution`: The calculated proposal distribution for the number of matching haplotypes.
                - `run_time_pedigree_probability`: The time taken to calculate the average pedigree probability.
                - `run_time_proposal_distribution`: The time taken to calculate the proposal distribution.

        Note:
            - The simulation assumes that the pedigree is re-rooted to have the suspect as the most recent common ancestor.
            - The simulation involves several time-consuming steps, such as calculating pedigree probabilities and proposal distributions.
            - If the average pedigree probability is zero, the simulation is considered impossible, and the function returns early.
        """
    config_path = Path(__file__).resolve().parent.parent / "data" / "config.ini"
    global_config = ConfigParser()
    global_config.optionxform = str  # type: ignore
    global_config.read(config_path)

    # Create deep copies of the input pedigree to preserve the original
    pedigree = deepcopy(input_pedigree)

    # write pedigree to tgf
    pedigree_bytes_data = pedigree.write_to_tgf()
    with open(f"{simulation_parameters.results_path}/original_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf", "wb") as f:
        f.write(pedigree_bytes_data)

    save_pedigree_to_png(pedigree=pedigree,
                         global_config=global_config,
                         results_path=simulation_parameters.results_path,
                         pedigree_name="original_pedigree")

    extended_pedigree = deepcopy(pedigree)

    # Extend the pedigree and store the name of the last added individual (used for edge cases)
    last_child_name = extended_pedigree.extend_pedigree()
    extended_pedigree_bytes_data = extended_pedigree.write_to_tgf()
    with open(f"{simulation_parameters.results_path}/extended_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf", "wb") as f:
        f.write(extended_pedigree_bytes_data)

    save_pedigree_to_png(pedigree=extended_pedigree,
                         global_config=global_config,
                         results_path=simulation_parameters.results_path,
                         pedigree_name="extended_pedigree")

    # Extract and store the suspect's haplotype before modifying the pedigree
    suspect_haplotype = deepcopy(pedigree.get_individual_by_name(suspect_name).haplotype)

    # Remove irrelevant individuals and re-root the pedigree with the most recent informative ancestor
    root_name = pedigree.remove_irrelevant_individuals(
        inside=True
    )
    relevant_pedigree_bytes_data = pedigree.write_to_tgf()
    with open(
            f"{simulation_parameters.results_path}/relevant_original_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf",
            "wb") as f:
        f.write(relevant_pedigree_bytes_data)

    pedigree.reroot_pedigree(root_name)
    pedigree.get_closest_known_individuals()
    rerooted_pedigree_bytes_data = pedigree.write_to_tgf()
    with open(f"{simulation_parameters.results_path}/rerooted_original_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf",
              "wb") as f:
        f.write(rerooted_pedigree_bytes_data)

    # Perform the same process on the extended pedigree (excluding 'inside' nodes)
    extended_root_name = extended_pedigree.remove_irrelevant_individuals(
        inside=False,
        last_child_name=last_child_name)
    extended_relevant_pedigree_bytes_data = extended_pedigree.write_to_tgf()
    with open(
            f"{simulation_parameters.results_path}/relevant_extended_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf",
            "wb") as f:
        f.write(extended_relevant_pedigree_bytes_data)

    extended_pedigree.reroot_pedigree(extended_root_name)
    extended_pedigree.get_closest_known_individuals()
    extended_rerooted_pedigree_bytes_data = extended_pedigree.write_to_tgf()
    with open(
            f"{simulation_parameters.results_path}/rerooted_extended_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf",
            "wb") as f:
        f.write(extended_rerooted_pedigree_bytes_data)

    # Compute a priori match probabilities for unknown individuals
    # These normalized probabilities are used to probabilistically assign the suspect haplotype
    if len(pedigree.individuals) > 0:
        pedigree.calculate_picking_probabilities()

    # Save the processed pedigree to disk for reproducibility/debugging
    timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
    output_file = f"{simulation_parameters.results_path}/pedigree_{timestamp}.pkl"
    with open(output_file, "wb") as f:
        # noinspection PyTypeChecker
        pickle.dump(pedigree, f)

    # In the extended pedigree, fix the last added child to have 100% picking probability
    last_child = extended_pedigree.get_individual_by_name(last_child_name)
    last_child.picking_probability = Decimal(1)

    # Set picking probability of all other individuals to 0
    for individual in extended_pedigree.individuals:
        if individual.id != last_child.id:
            individual.picking_probability = Decimal(0)

    extended_pedigree.picking_probabilities = {
        individual.id: individual.picking_probability
        for individual in extended_pedigree.individuals}

    """
    Step 1: calculate average pedigree probability. 
    This needs to be done only once for the original pedigree.
    This corresponds to P(hv)
    """
    max_number_of_threads = multiprocessing.cpu_count()
    number_of_threads = min(simulation_parameters.number_of_threads, max_number_of_threads)

    start_time_average_pedigree_probability = datetime.now()

    while True:  # loop to allow retrying in case of BadInputError
        try:
            if skip_inside:
                reporter.log("\nSkipping inside match probability calculation.")
                average_pedigree_probability = Decimal(1)
                average_pedigree_needed_iterations = [0, 0, 0]
                average_pedigree_model_pedigree_probabilities = [Decimal(1), Decimal(1), Decimal(1)]
                average_pedigree_model_is_valid = True
                average_pedigree_used_probabilities = [Decimal(1)]
            else:
                average_pedigree_probability, average_pedigree_needed_iterations, average_pedigree_model_pedigree_probabilities, average_pedigree_model_is_valid, average_pedigree_used_probabilities = calculate_average_pedigree_probability(
                    pedigree=pedigree,
                    root_name=root_name,
                    suspect_haplotype=suspect_haplotype,
                    marker_set=marker_set,
                    simulation_parameters=simulation_parameters,
                    random_seed=random_seed,
                    reporter=reporter,
                    is_outside=False,
                    number_of_threads=number_of_threads,
                )
            run_time_pedigree_probability = datetime.now() - start_time_average_pedigree_probability

            # If the average pedigree probability is zero, the simulation ends here,
            # because this means that the current pedigree is impossible (e.g. because of a 3-step mutation)
            if average_pedigree_probability == Decimal(0):
                simulation_result = SimulationResult(
                    pedigree=input_pedigree,
                    marker_set=marker_set,
                    root_name=root_name,
                    simulation_parameters=simulation_parameters,
                    random_seed=random_seed,

                    average_pedigree_probability=average_pedigree_probability,
                    extended_average_pedigree_probability=Decimal(0),
                    inside_match_probability={1: Decimal(0)},
                    outside_match_probability=Decimal(0),

                    average_pedigree_needed_iterations=average_pedigree_needed_iterations,
                    extended_needed_iterations=[],
                    inside_needed_iterations=[],
                    outside_needed_iterations=[],

                    average_pedigree_model_pedigree_probabilities=average_pedigree_model_pedigree_probabilities,
                    extended_model_pedigree_probabilities=[],
                    inside_model_probabilities=[],
                    outside_model_probabilities=[],

                    average_pedigree_model_is_valid=average_pedigree_model_is_valid,
                    extended_model_is_valid=False,
                    inside_model_is_valid=False,
                    outside_model_is_valid=False,

                    average_used_probabilities=average_pedigree_used_probabilities,
                    extended_used_probabilities=[],
                    inside_used_probabilities=[],
                    outside_used_probabilities=[],

                    run_time_pedigree_probability=run_time_pedigree_probability,
                    run_time_proposal_distribution=timedelta(0),
                    run_time_extended_average_pedigree_probability=timedelta(0),
                    run_time_outside_match_probability=timedelta(0),
                    total_run_time=run_time_pedigree_probability,
                )

                with open(f"{simulation_parameters.results_path}/simulation_result_{datetime.now().strftime('%Y%m%d%H%M%S')}.pkl", "wb") as f:
                    # noinspection PyTypeChecker
                    pickle.dump(simulation_result, f)
                return simulation_result

            """
            Step 2: calculate the number of matching haplotypes with the suspect.
            This is done for the total number of unknown individuals (x) in the pedigree.
            This corresponds to px=P(m(Hu)=x|hv).
            """
            start_time_proposal_distribution = datetime.now()
            if skip_inside:
                reporter.log("\nSkipping proposal distribution calculation.")
                inside_match_probability = {1: Decimal(0)}
                inside_needed_iterations = [0, 0, 0]
                inside_model_probabilities = [Decimal(0), Decimal(0), Decimal(0)]
                inside_model_validity = False
                inside_used_probabilities = [Decimal(0)]
            else:
                inside_match_probability, inside_needed_iterations, inside_model_probabilities, inside_model_validity, inside_used_probabilities = calculate_proposal_distribution(
                    pedigree=pedigree,
                    marker_set=marker_set,
                    root_name=root_name,
                    suspect_haplotype=suspect_haplotype,
                    average_pedigree_probability=average_pedigree_probability,
                    simulation_parameters=simulation_parameters,
                    random_seed=random_seed,
                    reporter=reporter,
                )
            run_time_proposal_distribution = datetime.now() - start_time_proposal_distribution
            break

        except InvalidAveragePedigreeProbability:
            continue

    """
    Step 3: calculate the outside match probability.
    This is done by extending the pedigree with an additional branch. One extra generation is added, above the most
    recent common ancestor (the current root of the tree), and then a branch is added all the way back to the last 
    generation. This last generation is the last child in the extended pedigree. 
    The outside match probability is calculated as the probability that this last child has the same 
    haplotype as the suspect.
    """

    reporter.log(f"\nCalculating outside match probability...")

    start_time_extended_average_pedigree_probability = datetime.now()

    while True:  # loop to allow retrying in case of BadInputError
        try:
            if skip_outside:
                reporter.log("\nSkipping extended average pedigree probability calculation.")
                extended_average_pedigree_probability = Decimal(1)
                extended_needed_iterations = [0, 0, 0]
                extended_model_pedigree_probabilities = [Decimal(1), Decimal(1), Decimal(1)]
                extended_model_is_valid = True
                extended_used_probabilities = [Decimal(1)]
            else:
                extended_average_pedigree_probability, extended_needed_iterations, extended_model_pedigree_probabilities, extended_model_is_valid, extended_used_probabilities = calculate_average_pedigree_probability(
                    pedigree=extended_pedigree,
                    root_name=extended_root_name,
                    suspect_haplotype=suspect_haplotype,
                    marker_set=marker_set,
                    simulation_parameters=simulation_parameters,
                    random_seed=random_seed,
                    reporter=reporter,
                    is_outside=True,
                    number_of_threads=number_of_threads,
                )
            run_time_extended_average_pedigree_probability = datetime.now() - start_time_extended_average_pedigree_probability

            start_time_outside_match_probability = datetime.now()
            if skip_outside:
                reporter.log("\nSkipping outside match probability calculation.")
                outside_match_probability = Decimal(0)
                outside_needed_iterations = [0, 0, 0]
                outside_model_probabilities = [Decimal(0), Decimal(0), Decimal(0)]
                outside_model_is_valid = False
                outside_used_probabilities = [Decimal(0)]
            else:
                outside_match_probability, outside_needed_iterations, outside_model_probabilities, outside_model_is_valid, outside_used_probabilities = calculate_outside_match_probability(
                    pedigree=extended_pedigree,
                    marker_set=marker_set,
                    root_name=extended_root_name,
                    suspect_haplotype=suspect_haplotype,
                    average_pedigree_probability=extended_average_pedigree_probability,
                    simulation_parameters=simulation_parameters,
                    random_seed=random_seed,
                    reporter=reporter,
                    number_of_threads=number_of_threads,
                )
            run_time_outside_match_probability = datetime.now() - start_time_outside_match_probability
            break
        except InvalidAveragePedigreeProbability:
            simulation_parameters.model_validity_threshold /= 3.375
            continue

    total_run_time = datetime.now() - start_time_average_pedigree_probability

    simulation_result = SimulationResult(
        pedigree=input_pedigree,
        root_name=root_name,
        marker_set=marker_set,
        simulation_parameters=simulation_parameters,
        random_seed=random_seed,

        average_pedigree_probability=average_pedigree_probability,
        extended_average_pedigree_probability=extended_average_pedigree_probability,
        inside_match_probability=inside_match_probability,
        outside_match_probability=outside_match_probability,

        average_pedigree_needed_iterations=average_pedigree_needed_iterations,
        extended_needed_iterations=extended_needed_iterations,
        inside_needed_iterations=inside_needed_iterations,
        outside_needed_iterations=outside_needed_iterations,

        average_pedigree_model_pedigree_probabilities=average_pedigree_model_pedigree_probabilities,
        extended_model_pedigree_probabilities=extended_model_pedigree_probabilities,
        inside_model_probabilities=inside_model_probabilities,
        outside_model_probabilities=outside_model_probabilities,

        average_pedigree_model_is_valid=average_pedigree_model_is_valid,
        extended_model_is_valid=extended_model_is_valid,
        inside_model_is_valid=inside_model_validity,
        outside_model_is_valid=outside_model_is_valid,

        average_used_probabilities=average_pedigree_used_probabilities,
        extended_used_probabilities=extended_used_probabilities,
        inside_used_probabilities=inside_used_probabilities,
        outside_used_probabilities=outside_used_probabilities,

        run_time_pedigree_probability=run_time_pedigree_probability,
        run_time_proposal_distribution=run_time_proposal_distribution,
        run_time_extended_average_pedigree_probability=run_time_extended_average_pedigree_probability,
        run_time_outside_match_probability=run_time_outside_match_probability,
        total_run_time=total_run_time,
    )

    plot_probabilities(
        simulation_result=simulation_result,
        results_path=simulation_parameters.results_path
    )

    with open(f"{simulation_parameters.results_path}/simulation_result_{datetime.now().strftime('%Y%m%d%H%M%S')}.pkl",
              "wb") as f:
    # noinspection PyTypeChecker
        pickle.dump(simulation_result, f)

    pdf_data = create_html_pdf_report(
        result=simulation_result
    )

    with open(f"{simulation_parameters.results_path}/pdf_report_{datetime.now().strftime('%Y%m%d%H%M%S')}.pdf",
              "wb") as f:
        f.write(pdf_data)

    return simulation_result
