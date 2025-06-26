from functools import lru_cache
import math
import multiprocessing
import operator
import pickle
from collections import defaultdict
from copy import deepcopy
from datetime import datetime, timedelta
from functools import reduce
from random import Random
from typing import Collection, Mapping, Sequence
import pandas as pd
from _decimal import Decimal
from pedigree_lr.reporting import Reporter, create_html_pdf_report
from pedigree_lr.visualization import plot_probabilities
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
    get_single_copy_mutation_rate, SimulationParameters
)


def mutate_allele(
        marker: Marker,
        source_alleles: list[Allele],
        random: Random,
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
    return mutated_alleles


def mutate_haplotype(
        source: Haplotype,
        marker_set: MarkerSet,
        random: Random,
        two_step_mutation_factor: float,
) -> Haplotype:
    """
        Generates a mutated version of a given haplotype by applying mutations to its alleles.

        Args:
            source (Haplotype): The original haplotype to be mutated.
            marker_set (MarkerSet): A set of markers with corresponding mutation rates.
            random (Random): A random number generator instance to introduce stochasticity in mutation.
            two_step_mutation_factor (float): The factor that determines the probability
                of a two-step mutation. Defaults to 0.03.

        Returns:
            Haplotype: A new haplotype instance with mutated alleles.

        Note:
            - Each marker in `marker_set` is iterated over, and its corresponding alleles in `source`
              undergo mutation through the `mutate_allele` function.
            - The mutation process is controlled by the `random` instance.
        """

    target_haplotype = Haplotype()

    for marker in marker_set.markers:
        source_alleles = source.alleles[marker.name]
        target_haplotype.alleles[marker.name] = mutate_allele(
            marker=marker,
            source_alleles=source_alleles,
            random=random,
            two_step_mutation_factor=two_step_mutation_factor,
        )

    return target_haplotype


def serialize_haplotype(haplotype: Haplotype) -> tuple:
    return tuple(
        (marker_name, tuple(sorted((allele.value, allele.intermediate_value) for allele in alleles)))
        for marker_name, alleles in sorted(haplotype.alleles.items())
    )


@lru_cache(maxsize=100_000)
def _cached_edge_probability(
    serialized_source: tuple,
    serialized_target: tuple,
    marker_names: tuple,
    marker_rates: tuple,
    marker_copies: tuple,
    two_step_mutation_factor: float
) -> Decimal:
    edge_probability = Decimal(1)

    for i, marker_name in enumerate(marker_names):
        rate = marker_rates[i]
        copies = marker_copies[i]
        marker = Marker(name=marker_name, mutation_rate=rate, number_of_copies=copies)

        source_alleles = [Allele(marker, v, iv) for v, iv in serialized_source[i][1]]
        target_alleles = [Allele(marker, v, iv) for v, iv in serialized_target[i][1]]

        mutation_probability = calculate_mutation_probability(
            sorted(source_alleles, key=lambda a: a.value),
            sorted(target_alleles, key=lambda a: a.value),
            marker,
            two_step_mutation_factor
        )
        edge_probability *= mutation_probability

    return edge_probability


def get_edge_probability(
    source: Haplotype,
    target: Haplotype,
    marker_set: MarkerSet,
    two_step_mutation_factor: float
) -> Decimal:
    serialized_source = serialize_haplotype(source)
    serialized_target = serialize_haplotype(target)
    marker_names = tuple(marker.name for marker in marker_set.markers)
    marker_rates = tuple(marker.mutation_rate for marker in marker_set.markers)
    marker_copies = tuple(marker.number_of_copies for marker in marker_set.markers)

    return _cached_edge_probability(
        serialized_source,
        serialized_target,
        marker_names,
        marker_rates,
        marker_copies,
        two_step_mutation_factor
    )


def get_edge_probabilities(
        haplotypes: Mapping[int, Haplotype],
        relationships: Collection[Relationship],
        marker_set: MarkerSet,
        two_step_mutation_factor: float
) -> Mapping[tuple[int, int], Decimal]:
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
        )
        for relationship in relationships
    }


def update_average(
        new_probability: Decimal,
        old_probability: Decimal,
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
        reporter.log("Probabilities are all zero. Pedigree is impossible.")
        return True
    return False


def is_model_valid(
        probabilities: list[Decimal],
        threshold: float = 0.005
) -> bool:
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

    mean_probability = Decimal(sum(probabilities) / len(probabilities))
    return all(
        abs((probabilities[i] - mean_probability) / mean_probability) < threshold
        if mean_probability != 0 else probabilities[i] == 0
        for i in range(len(probabilities))
    )


def simulate_pedigree_probability(
        individuals: Mapping[int, Individual],
        relationships: Collection[Relationship],
        parents: Mapping[int, int],
        ordered_unknown_ids: Collection[int],
        marker_set: MarkerSet,
        random: Random,
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
            random (Random): A random number generator used to simulate the mutation of haplotypes.
            two_step_mutation_factor (float): A factor that determines the probability of a two-step mutation. Defaults to 0.03.

        Returns:
            IterationResult: An object containing the calculated pedigree probability, and the edge probabilities for the relationships.

        Note:
            - The function simulates the haplotypes of individuals with unknown haplotypes by mutating their parent's haplotype.
            - It calculates the edge probabilities for parent-child relationships and uses them to compute the overall pedigree probability.
            - If any edge probability is zero, the pedigree probability is set to zero to avoid underflow.
        """

    simulated_individual_ids: set[int] = set()
    simulated_relationship_ids: set[tuple[int, int]] = set()

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
            random=random,
            two_step_mutation_factor=two_step_mutation_factor,
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
        haplotypes, unused_relationships, marker_set, two_step_mutation_factor
    )

    if any(probability == 0 for probability in edge_probabilities.values()):
        pedigree_probability = Decimal(0)
    else:
        log_sum = sum(math.log(p) for p in edge_probabilities.values())
        pedigree_probability = math.exp(log_sum)  # TODO: Solve underflow

    return IterationResult(
        probability=Decimal(pedigree_probability),
        edge_probabilities=edge_probabilities,
        mutated_haplotypes=haplotypes,
        fixed_individuals_ids=None,
    )


def calculate_average_pedigree_probability(
        pedigree: Pedigree,
        suspect_name: str,
        marker_set: MarkerSet,
        simulation_parameters: SimulationParameters,
        random: Random,
        reporter: Reporter,
        is_outside: bool,
) -> Decimal:
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
        suspect_name (str): The name of the suspect, used as the root of the pedigree.
        marker_set (MarkerSet): A set of genetic markers used in haplotype prediction.
        simulation_parameters (Mapping[str, any]): A dictionary containing simulation parameters such as
        random (Random): A random number generator for controlled stochastic processes.
        reporter (Reporter): A reporting tool for logging and tracking progress.
        is_outside (bool): A flag indicating whether the simulation is performed on the extended pedigree.

    Returns:
        Decimal: The estimated average pedigree probability after running the simulation.

    Notes:
        - The function maintains a running average of pedigree probabilities to detect stabilization.
        - If the probability estimate stabilizes before reaching `number_of_iterations`, the function
          terminates early.
        - Edge probabilities are logged, showing how likely specific parent-child haplotype transmissions are.
        - The final probability is reported both as a raw value and a log-transformed value.
    """

    progress_bar = reporter.progress_bar(
        total=simulation_parameters.number_of_iterations, desc="Calculating average pedigree probability"
    )

    individuals = {individual.id: individual for individual in pedigree.individuals}

    child_parent_dict = {
        relationship.child_id: relationship.parent_id
        for relationship in pedigree.relationships
    }

    ordered_unknown_ids = [
        individual.id
        for individual in pedigree.get_level_order_traversal(suspect_name)
        if individual.haplotype_class == "unknown"
    ]

    running_average_pedigree_probability = Decimal()
    average_pedigree_probabilities_list: list[Decimal] = []  # This list is used to check for stability
    edge_probabilities: dict[tuple[int, int], list[Decimal]] = {
        (relationship.parent_id, relationship.child_id): []
        for relationship in pedigree.relationships
    }

    with (progress_bar):
        with open(f"{simulation_parameters.results_path}/average_pedigree_probabilities_outside_{is_outside}_{datetime.now().strftime('%Y%m%d%H%M%S')}.txt", "w") as f:
            # TODO: implement model validation
            for iteration_id in range(simulation_parameters.number_of_iterations):
                iteration_result = simulate_pedigree_probability(
                    individuals=individuals,
                    relationships=pedigree.relationships,
                    parents=child_parent_dict,
                    ordered_unknown_ids=ordered_unknown_ids,
                    marker_set=marker_set,
                    random=random,
                    two_step_mutation_factor=simulation_parameters.two_step_mutation_factor,
                )

                running_average_pedigree_probability = update_average(
                    old_probability=running_average_pedigree_probability,
                    new_probability=iteration_result.probability,
                    i=iteration_id
                )

                average_pedigree_probabilities_list.append(running_average_pedigree_probability)

                # TODO: is_stable is very inefficient
                if iteration_id % 100 == 0:
                    f.write(f"{running_average_pedigree_probability}\n")

                if iteration_id % simulation_parameters.stability_window == 0 and iteration_id > simulation_parameters.stability_min_iterations and running_average_pedigree_probability < Decimal(1.0):
                    if is_stable(
                            probabilities=average_pedigree_probabilities_list,
                            threshold=simulation_parameters.stability_threshold,
                            reporter=reporter
                    ):
                        reporter.log(f"Average pedigree probability is stable after {iteration_id} iterations.")
                        break
                    average_pedigree_probabilities_list = []

                # edge_probabilities = {
                #     edge: probabilities + [iteration_result.edge_probabilities.get(edge, 0)]
                #     for edge, probabilities in edge_probabilities.items()
                # }

                progress_bar.update(1)

                if iteration_id == simulation_parameters.number_of_iterations - 1:
                    reporter.log(f"Average pedigree probability is not stable after {simulation_parameters.number_of_iterations} iterations. "
                                 f"Increase the number of iterations.")

    # average_edge_probabilities = {
    #     edge: sum(edge_probabilities) / simulation_parameters['number_of_iterations']
    #     for edge, edge_probabilities in edge_probabilities.items()
    # }

    reporter.log(
        f"Average pedigree probability P(Hv): "
        f"{running_average_pedigree_probability} (log {10 * math.log10(running_average_pedigree_probability)})"
    )

    return running_average_pedigree_probability


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
        two_step_mutation_factor: float,
        picking_probabilities: dict[int, dict[tuple[int], Decimal]],
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
        l (int): The number of unknown individuals that should match the suspect’s haplotype.
        average_pedigree_probability (Decimal): The average probability of the pedigree, used for normalization.
        random (Random): A random number generator to ensure controlled stochastic behavior.

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

    simulated_individual_ids: set[int] = set()
    simulated_relationship_ids: set[tuple[int, int]] = set()
    fixed_individual_ids: set[int] = set()

    simulation_probability = Decimal(1)

    if l > 0:
        l_picking_probabilities = picking_probabilities[l]
        normalized_picking_probabilities = {
            combination: probability / sum(l_picking_probabilities.values()) for combination, probability in
            l_picking_probabilities.items()
        }
        # print(normalized_picking_probabilities)
        # pick a combination of individuals based on the normalized picking probabilities
        picked_combination = random.choices(
            population=list(normalized_picking_probabilities.keys()),
            weights=[float(normalized_picking_probabilities[comb]) for comb in l_picking_probabilities.keys()],
        )[0]
        # add the picked individual to the fixed individuals
        for individual_id in picked_combination:
            fixed_individual_ids.add(individual_id)
        # get the picking probability of the picked individual
        picking_probability = normalized_picking_probabilities[picked_combination]
        simulation_probability *= picking_probability

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
            source=parent_haplotype,
            marker_set=marker_set,
            random=random,
            two_step_mutation_factor=two_step_mutation_factor,
        )

        simulated_individual_ids.add(individual_id)
        simulated_relationship_ids.add((parent_id, individual_id))

    # Calculate the probability of the entire pedigree (Pr(H = h))
    all_edge_probabilities = get_edge_probabilities(
        haplotypes, relationships, marker_set, two_step_mutation_factor
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

    # Check if the total number of matching haplotypes is equal to l (excluding the suspect and excluded individuals)
    non_excluded_matching_ids = [
        individual_id for individual_id in ordered_unknown_ids
        if haplotypes[individual_id] == suspect.haplotype and not individuals[individual_id].exclude
    ]

    number_of_non_excluded_matching_haplotypes = len(non_excluded_matching_ids)

    total_matching_ids = fixed_individual_ids.union(set(non_excluded_matching_ids))
    total_number_of_matching_haplotypes = len(total_matching_ids)

    if is_outside:
        # number of matching has to be exactly 1
        probability = (
            conditional_probability / (simulation_probability * total_number_of_matching_haplotypes)
            if number_of_non_excluded_matching_haplotypes == l and simulation_probability != 0
            else Decimal(0)
        )
    else:
        # number of matching has to be at least 1
        probability = (
            conditional_probability / (simulation_probability * total_number_of_matching_haplotypes)
            if number_of_non_excluded_matching_haplotypes >= l and simulation_probability != 0
            else Decimal(0)
        )

    return IterationResult(
        probability=probability,
        edge_probabilities=all_edge_probabilities,
        mutated_haplotypes=haplotypes,
        fixed_individuals_ids=fixed_individual_ids,
    )


def calculate_l_matching_haplotypes(
        pedigree: Pedigree,
        marker_set: MarkerSet,
        suspect_name: str,
        simulation_parameters: SimulationParameters,
        l: int,
        average_pedigree_probability: Decimal,
        random: Random,
        reporter: Reporter,
        is_outside: bool,
        number_of_threads: int = 1,
) -> tuple[Decimal, list[int], list[Decimal], bool]:
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
        suspect_name (str): The name of the suspect whose haplotype is used as a reference.
        simulation_parameters (Mapping[str, any]): A dictionary containing simulation parameters such as
        l (int): The fixed number of unknown individuals whose haplotypes are set to match the suspect.
        average_pedigree_probability (Decimal): The precomputed average probability of the pedigree.
        random (Random): A random number generator for controlled stochastic processes.
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

    l_probability = Decimal()

    average_edge_probabilities: dict[tuple[int, int], Decimal] = {
        (relationship.parent_id, relationship.child_id): Decimal(0)
        for relationship in pedigree.relationships
    }

    progress_bar = reporter.progress_bar(
        total=simulation_parameters.number_of_iterations * 3, desc=f"Calculating {l} matching haplotypes"
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
    picking_probabilities = pedigree.picking_probabilities

    with (progress_bar):
        model_probabilities = []
        needed_iterations = []

        for m in range(3):  # Model validation
            l_probabilities = []

            with open(f"{simulation_parameters.results_path}/{l}_pedigree_probabilities_model_{m}_outside_{is_outside}_{datetime.now().strftime('%Y%m%d%H%M%S')}.txt", "w", 1) as f:
                for i in range(simulation_parameters.number_of_iterations):
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
                        two_step_mutation_factor=simulation_parameters.two_step_mutation_factor,
                        picking_probabilities=picking_probabilities,
                        is_outside=is_outside,
                    )

                    l_probability = update_average(old_probability=l_probability,
                                                   new_probability=iteration_result.probability,
                                                   i=i)
                    if i > simulation_parameters.stability_window:
                        if len(l_probabilities) > 1:
                            if abs((l_probability - l_probabilities[-1]) / l_probabilities[-1]) > 0.2:
                                print(f"Iteration {i} - Model {m + 1} - l={l}")
                                print(iteration_result.probability)
                                print(iteration_result.edge_probabilities)
                                print(iteration_result.fixed_individuals_ids)

                                alleles_dict = defaultdict(lambda: defaultdict(list))
                                for id, haplotype in iteration_result.mutated_haplotypes.items():
                                    for marker, alleles in haplotype.alleles.items():
                                        for allele in alleles:
                                            alleles_dict[marker][id].append(str(allele))

                                df = pd.DataFrame.from_dict(alleles_dict, orient="index")
                                print(df)
                                df.index.name = "Marker"
                                df.to_csv(f"{simulation_parameters.results_path}/alleles_model_{m + 1}_l_{l}_iteration_{i}.csv")

                                with open(f"{simulation_parameters.results_path}/alleles_model_{m + 1}_l_{l}_iteration_{i}.txt", "w") as logf:
                                    logf.write(f"Model {m + 1} - l={l} - Iteration {i}\n")
                                    logf.write(f"Probability: {iteration_result.probability}\n")
                                    logf.write(f"Edge probabilities: {iteration_result.edge_probabilities}\n")
                                    logf.write(f"Fixed individuals: {iteration_result.fixed_individuals_ids}\n")
                                    logf.write(f"Suspect id: {suspect.id}\n")

                    l_probabilities.append(l_probability)

                    if i % 100 == 0:
                        f.write(f"{l_probability}\n")

                    if i % simulation_parameters.stability_window == 0 and i > simulation_parameters.stability_min_iterations and l_probability < Decimal(
                            1.0):
                        if is_stable(
                                probabilities=l_probabilities,
                                threshold=simulation_parameters.stability_threshold,
                                reporter=reporter
                        ):
                            model_probabilities.append(l_probability)
                            needed_iterations.append(i)
                            reporter.log(f"Model {m + 1} for l={l} is stable after {i} iterations.")
                            progress_bar.update_total(
                                simulation_parameters.number_of_iterations - i
                            )
                            break
                        l_probabilities = []

                    average_edge_probabilities = {
                        edge: update_average(
                            old_probability=average_edge_probability,
                            new_probability=iteration_result.edge_probabilities.get(edge, 0),
                            i=i
                        )
                        for edge, average_edge_probability in average_edge_probabilities.items()
                    }

                    progress_bar.update(1)

                    if i == simulation_parameters.number_of_iterations - 1:
                        model_probabilities.append(l_probability)
                        needed_iterations.append(i)
                        reporter.log(f"Model {m + 1} for l={l} is not stable after {simulation_parameters.number_of_iterations} iterations. "
                                     f"Increase the number of iterations.")

        # check if all three model probabilities are stable (<0.5% change)
        if len(model_probabilities) == 3:
            mean_probability = Decimal(sum(model_probabilities) / len(model_probabilities))

        model_is_valid = is_model_valid(
                probabilities=model_probabilities,
                threshold=simulation_parameters.model_validity_threshold,
        )

        if model_is_valid:
            reporter.log(f"Model for l={l} is valid. Probability of {l} matching haplotypes: {mean_probability}")
        else:
            reporter.log(f"Model for l={l} is invalid! This may be the effect of individual models reaching a local optimum. "
                         f"Try to increase the stability window or decrease the stability threshold.")
            reporter.log(f"Model probabilities: {model_probabilities}")


    return mean_probability, needed_iterations, model_probabilities, model_is_valid


def calculate_proposal_distribution(
        pedigree: Pedigree,
        marker_set: MarkerSet,
        suspect_name: str,
        average_pedigree_probability: Decimal,
        simulation_parameters: SimulationParameters,
        random: Random,
        reporter: Reporter,
) -> tuple[Mapping[int, Decimal], Mapping[int, list[int]], Mapping[int, list[Decimal]], Mapping[int, bool]]:
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
            suspect_name (str): The name of the suspect whose haplotype is used as a reference.
            average_pedigree_probability (Decimal): The precomputed average probability of the pedigree.
            simulation_parameters (Mapping[str, any]): A dictionary containing simulation parameters such as
            random (Random): A random number generator for controlled stochastic processes.
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
    number_of_unknowns = len(unknown_individuals)
    proposal_distribution: dict[int, Decimal] = {}
    needed_iterations: dict[int, list[int]] = {}
    model_probabilities: dict[int, list[Decimal]] = {}
    model_validities: dict[int, bool] = {}
    number_of_excluded_individuals = len([individual for individual in pedigree.individuals if individual.exclude])

    # There has to be at least one non-excluded unknown individual in the pedigree
    if number_of_unknowns <= number_of_excluded_individuals:
        reporter.log("Warning! No unknown individuals left in the pedigree. Only calculating outside match probability.")
        proposal_distribution[0] = Decimal(1)
        return proposal_distribution, needed_iterations, model_probabilities, model_validities

    max_number_of_threads = multiprocessing.cpu_count()
    number_of_threads = min(simulation_parameters.number_of_threads, max_number_of_threads)

    proposal_distribution[1], needed_iterations[1], model_probabilities[1], model_validities[1] = calculate_l_matching_haplotypes(
        pedigree=pedigree,
        marker_set=marker_set,
        suspect_name=suspect_name,
        l=1,
        average_pedigree_probability=average_pedigree_probability,
        simulation_parameters=simulation_parameters,
        random=random,
        reporter=reporter,
        is_outside=False,
        number_of_threads=number_of_threads,
    )

    proposal_distribution[0] = Decimal(1) - sum(proposal_distribution.values())
    if proposal_distribution[0] < 0:
        reporter.log("Warning! Proposal distribution does not add up to unity.")
    return proposal_distribution, needed_iterations, model_probabilities, model_validities


def calculate_outside_match_probability(
        pedigree: Pedigree,
        marker_set: MarkerSet,
        suspect_name: str,
        average_pedigree_probability: Decimal,
        simulation_parameters: SimulationParameters,
        random: Random,
        reporter: Reporter,
)   -> tuple[Decimal, list[int], list[Decimal], bool]:
    outside_match_probability, needed_iterations, model_probabilities, model_is_valid = calculate_l_matching_haplotypes(
        pedigree=pedigree,
        marker_set=marker_set,
        suspect_name=suspect_name,
        l=1,
        average_pedigree_probability=average_pedigree_probability,
        simulation_parameters=simulation_parameters,
        random=random,
        reporter=reporter,
        is_outside=True,
    )
    return outside_match_probability, needed_iterations, model_probabilities, model_is_valid


def run_simulation(
        pedigree: Pedigree,
        suspect_name: str,
        marker_set: MarkerSet,
        simulation_parameters: SimulationParameters,
        random: Random,
        reporter: Reporter,
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
            pedigree (Pedigree): The pedigree object containing the family structure and genetic information.
            suspect_name (str): The name of the suspect whose haplotype is being compared to others.
            marker_set (MarkerSet): A set of genetic markers used for calculating allele probabilities.
            simulation_parameters (Mapping[str, float]): A dictionary containing simulation parameters such as
                - `number_of_iterations`: The number of iterations for the Monte-Carlo simulation.
            random (Random): A random number generator used for random selection and mutation.
            reporter (Reporter): A reporter object used to track the progress and output of the simulation.

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

    extended_pedigree = deepcopy(pedigree)
    last_child_name = extended_pedigree.extend_pedigree()

    # The pedigree is re-rooted to have the suspect as the most recent common ancestor.
    # This is important for the level order traversal in which the unknown individuals are processed
    pedigree.reroot_pedigree(suspect_name)
    extended_pedigree.reroot_pedigree(suspect_name)

    # The a priori match probabilities is calculated for all unknown individuals in the pedigree.
    # These (normalized) probabilities are used to "randomly" pick the l number of individuals
    # which will be fixed to the suspect haplotype
    pedigree.calculate_picking_probabilities(marker_set)

    with open(f"{simulation_parameters.results_path}/pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.pkl", "wb") as f:
        pickle.dump(pedigree, f)

    last_child_individual = extended_pedigree.get_individual_by_name(last_child_name)
    last_child_individual.picking_probability = Decimal(1)
    for individual in extended_pedigree.individuals:
        if individual.id != last_child_individual.id:
            individual.picking_probability = Decimal(0)

    extended_pedigree.picking_probabilities = {1: {
        tuple([individual.id]): individual.picking_probability
        for individual in extended_pedigree.individuals}
    }

    """
    Step 1: calculate average pedigree probability. 
    This needs to be done only once for the original pedigree.
    This corresponds to P(hv)
    """
    start_time_average_pedigree_probability = datetime.now()
    average_pedigree_probability = calculate_average_pedigree_probability(
        pedigree=pedigree,
        suspect_name=suspect_name,
        marker_set=marker_set,
        simulation_parameters=simulation_parameters,
        random=random,
        reporter=reporter,
        is_outside=False,
    )
    run_time_pedigree_probability = datetime.now() - start_time_average_pedigree_probability

    # If the average pedigree probability is zero, the simulation ends here,
    # because this means that the current pedigree is impossible (i.e. because of a 3-step mutation)
    if average_pedigree_probability == Decimal(0):
        simulation_result = SimulationResult(
            pedigree=pedigree,
            suspect_name=suspect_name,
            marker_set=marker_set,
            simulation_parameters=simulation_parameters,
            random=random,
            average_pedigree_probability=average_pedigree_probability,
            proposal_distribution={},
            outside_match_probability=Decimal(0),
            run_time_pedigree_probability=run_time_pedigree_probability,
            run_time_proposal_distribution=timedelta(0),
            total_run_time=timedelta(0),
            l_needed_iterations={},
            l_model_probabilities={},
            l_model_validities={},
            outside_needed_iterations=[],
            outside_model_probabilities=[],
            outside_model_is_valid=False,
        )

        with open(f"{simulation_parameters.results_path}/simulation_result_{datetime.now().strftime('%Y%m%d%H%M%S')}.pkl", "wb") as f:
            pickle.dump(simulation_result, f)
        return simulation_result

    """
    Step 2: calculate the number of matching haplotypes with the suspect.
    This is done for the total number of unknown individuals (x) in the pedigree.
    This corresponds to px=P(m(Hu)=x|hv).
    """
    start_time_proposal_distribution = datetime.now()

    l_matching_haplotypes_probability, l_needed_iterations, l_model_probabilities, l_model_validities = calculate_proposal_distribution(
        pedigree=pedigree,
        marker_set=marker_set,
        suspect_name=suspect_name,
        average_pedigree_probability=average_pedigree_probability,
        simulation_parameters=simulation_parameters,
        random=random,
        reporter=reporter,
    )

    run_time_proposal_distribution = datetime.now() - start_time_proposal_distribution

    """
    Step 3: calculate the outside match probability.
    This is done by extending the pedigree with an additional branch. One extra generation is added, above the most
    recent common ancestor (the current root of the tree), and then a branch is added all the way back to the last 
    generation. This last generation is the last child in the extended pedigree. 
    The outside match probability is calculated as the probability that this last child has the same 
    haplotype as the suspect.
    """

    # reporter.log(f"Calculating outside match probability...")

    # extended_pedigree_probability = calculate_average_pedigree_probability(
    #     pedigree=extended_pedigree,
    #     suspect_name=suspect_name,
    #     marker_set=marker_set,
    #     simulation_parameters=simulation_parameters,
    #     random=random,
    #     reporter=reporter,
    #     is_outside=True,
    # )

    # outside_match_probability, outside_needed_iterations, outside_model_probabilities, outside_model_is_valid = calculate_outside_match_probability(
    #     pedigree=extended_pedigree,
    #     marker_set=marker_set,
    #     suspect_name=suspect_name,
    #     average_pedigree_probability=extended_pedigree_probability,
    #     simulation_parameters=simulation_parameters,
    #     random=random,
    #     reporter=reporter,
    # )

    outside_match_probability, outside_needed_iterations, outside_model_probabilities, outside_model_is_valid = Decimal(0), [], [], False

    total_run_time = datetime.now() - start_time_average_pedigree_probability

    plot_probabilities(simulation_parameters.results_path, list(l_matching_haplotypes_probability.keys()))

    simulation_result = SimulationResult(
        pedigree=pedigree,
        suspect_name=suspect_name,
        marker_set=marker_set,
        simulation_parameters=simulation_parameters,
        random=random,
        average_pedigree_probability=average_pedigree_probability,
        proposal_distribution=l_matching_haplotypes_probability,
        l_needed_iterations=l_needed_iterations,
        l_model_probabilities=l_model_probabilities,
        l_model_validities=l_model_validities,
        outside_match_probability=outside_match_probability,
        outside_needed_iterations=outside_needed_iterations,
        outside_model_probabilities=outside_model_probabilities,
        outside_model_is_valid=outside_model_is_valid,
        run_time_pedigree_probability=run_time_pedigree_probability,
        run_time_proposal_distribution=run_time_proposal_distribution,
        total_run_time=total_run_time,
    )

    with open(f"{simulation_parameters.results_path}/simulation_result_{datetime.now().strftime('%Y%m%d%H%M%S')}.pkl",
              "wb") as f:
        pickle.dump(simulation_result, f)

    # with open(f"{simulation_parameters.results_path}/report_{datetime.now().strftime('%Y%m%d%H%M%S')}.txt",
    #           "wb") as f:
    #     bytes_data = simulation_result.download_results(simulation_parameters=simulation_parameters)
    #     f.write(bytes_data)

    pdf_data = create_html_pdf_report(
        result=simulation_result
    )

    with open(f"{simulation_parameters.results_path}/pdf_report_{datetime.now().strftime('%Y%m%d%H%M%S')}.pdf",
              "wb") as f:
        f.write(pdf_data)

    return simulation_result
