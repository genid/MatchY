import csv
from collections import defaultdict
from configparser import ConfigParser
from decimal import Decimal
import math
import multiprocessing
import operator
import pickle
from functools import partial
from copy import deepcopy
from datetime import datetime, timedelta
from functools import reduce
from pathlib import Path
from random import Random
import os
from typing import Collection, Mapping, Sequence, List, Dict
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
    SimulationParameters,
    InvalidAveragePedigreeProbability,
    Bias,
)


def _rng() -> Random:
    global _RND
    if _RND is None:
        _RND = Random(os.getpid())
    return _RND


def _neutral_step_probabilities(marker: Marker, two_step_mutation_factor: float) -> Dict[int, float]:
    mu = float(get_single_copy_mutation_rate(marker.mutation_rate, marker.number_of_copies))
    two = max(mu * float(two_step_mutation_factor), 0.0)
    one = max(mu - two, 0.0)
    p0 = max(1.0 - mu, 0.0)
    return {-2: two / 2.0, -1: one / 2.0, 0: p0, +1: one / 2.0, +2: two / 2.0}


def mutate_alleles(marker: Marker, source_alleles: List[Allele], two_step_mutation_factor: float, biases: List[Bias]):
    rnd = _rng()
    mutated_alleles: List[Allele] = []
    weighted_prob_prod = 1.0
    unweighted_prob_prod = 1.0
    for copy_nr, source_allele in enumerate(source_alleles):
        unweighted = _neutral_step_probabilities(marker, two_step_mutation_factor)
        weighted = unweighted.copy()
        for bias in biases:
            if bias.marker == marker and bias.copy_nr == copy_nr:
                if bias.direction == "up":
                    weighted[1] += bias.target_mass
                elif bias.direction == "down":
                    weighted[-1] += bias.target_mass
                weighted[0] -= bias.target_mass
        step = rnd.choices(list(weighted.keys()), weights=list(weighted.values()))[0]
        mutated_value = max(1, source_allele.value + step)
        mutated_alleles.append(Allele(marker, mutated_value, source_allele.intermediate_value))
        weighted_prob_prod *= weighted[step]
        unweighted_prob_prod *= unweighted[step]
    return sorted(mutated_alleles), weighted_prob_prod, unweighted_prob_prod


def mutate_haplotype(source: Haplotype, marker_set: MarkerSet, two_step_mutation_factor: float, biases: List[Bias]):
    target = Haplotype()
    w_edge = 1.0
    u_edge = 1.0
    for marker in marker_set.markers:
        source_alleles = source.alleles[marker.name]
        alleles, w, u = mutate_alleles(marker, source_alleles, two_step_mutation_factor, biases)
        target.alleles[marker.name] = alleles
        w_edge *= w
        u_edge *= u
    return target, w_edge, u_edge


def get_edge_probability(source: Haplotype, target: Haplotype, marker_set: MarkerSet, two_step_mutation_factor: float,
                         is_average_pedigree: bool) -> Decimal:
    edge = Decimal(1)
    for marker in marker_set.markers:
        src = sorted(source.alleles[marker.name], key=lambda a: a.value)
        dst = sorted(target.alleles[marker.name], key=lambda a: a.value)
        p = calculate_mutation_probability(parent_alleles=src, child_alleles=dst, marker=marker,
                                           two_step_mutation_factor=two_step_mutation_factor,
                                           is_average_pedigree=is_average_pedigree)
        edge *= p
    return edge


def get_edge_probabilities(haplotypes: Mapping[str, Haplotype], relationships: Collection[Relationship],
                           marker_set: MarkerSet, two_step_mutation_factor: float, is_average_pedigree: bool):
    return {
        (rel.parent_id, rel.child_id): get_edge_probability(haplotypes[rel.parent_id], haplotypes[rel.child_id],
                                                            marker_set, two_step_mutation_factor, is_average_pedigree)
        for rel in relationships
    }


def simulate_pedigree_probability(
        pedigree: Pedigree,
        root_name: str,
        marker_set: MarkerSet,
        two_step_mutation_factor: float,
        bias_value: float = 0.1
) -> IterationResult:
    child_of = {rel.child_id: rel.parent_id for rel in pedigree.relationships}
    ordered_unknown = [ind.id for ind in pedigree.get_level_order_traversal(root_name) if
                       ind.haplotype_class == "unknown"]
    haplotypes = {ind.id: ind.haplotype for ind in pedigree.individuals}
    simulated = set()
    w_total = 1.0
    u_total = 1.0

    for uid in ordered_unknown:
        pid = child_of[uid]
        parent = haplotypes[pid]
        hap_tuple = tuple(sorted(haplotypes.items(), key=lambda kv: str(kv[0])))
        biases = pedigree.get_biases(uid, marker_set, hap_tuple, bias_value)
        haplotypes[uid], w, u = mutate_haplotype(parent, marker_set, two_step_mutation_factor, biases)
        w_total *= w
        u_total *= u
        simulated.add((pid, uid))

    unused = [rel for rel in pedigree.relationships if
              (rel.parent_id, rel.child_id) not in simulated and (rel.child_id, rel.parent_id) not in simulated]
    edges = get_edge_probabilities(haplotypes, unused, marker_set, two_step_mutation_factor, True)
    w_factor = Decimal(0) if w_total == 0 else Decimal(u_total) / Decimal(w_total)
    if any(p == 0 for p in edges.values()):
        ped_p = Decimal(0)
    else:
        ped_p = reduce(operator.mul, edges.values(), 1)
    return IterationResult(probability=Decimal(ped_p), edge_probabilities=edges, edge_weight_factor=w_factor,
                           mutated_haplotypes=haplotypes, fixed_individual_id=None)


def simulate_pedigree_iteration(i: int, root_name: str, pedigree: Pedigree, marker_set: MarkerSet,
                                two_step_mutation_factor: float, bias_value: float = 0.1) -> IterationResult:
    return simulate_pedigree_probability(pedigree, root_name, marker_set, two_step_mutation_factor, bias_value)


def _a1_pool_init(sim_fn, ctx: dict):
    global _A1_SIM_FUNC, _A1_SIM_CTX, _RND
    _A1_SIM_FUNC = sim_fn
    _A1_SIM_CTX = dict(ctx or {})
    seed_base = int(_A1_SIM_CTX.get('random_seed', _A1_SIM_CTX.get('seed_base', 0)))
    _RND = Random(seed_base + os.getpid())


def _a1_run_batch(args):
    start, count = args
    out_pairs = []
    per_ind_weighted = {}
    per_pair_weighted = {}
    sim = _A1_SIM_FUNC
    kw = _A1_SIM_CTX or {}
    suspect = kw.get('suspect_haplotype', None)
    ordered_unknown_ids = kw.get('ordered_unknown_ids', [])
    individuals = kw.get('individuals', {})
    is_inside = (suspect is not None) and (kw.get('is_outside', False) is False)
    for j in range(start, start + count):
        it = sim(i=j, **kw)
        prob = it.probability
        w = Decimal(it.edge_weight_factor)
        out_pairs.append((prob, w))
        if is_inside and prob != 0:
            haplotypes = it.mutated_haplotypes
            try:
                matches = [uid for uid in ordered_unknown_ids if
                           (haplotypes.get(uid) == suspect) and (not individuals[uid].exclude)]
            except KeyError:
                matches = []
            total_matching_ids = {it.fixed_individual_id}
            for uid in ordered_unknown_ids:
                if haplotypes.get(uid) == suspect:
                    total_matching_ids.add(uid)
            K = len(total_matching_ids)
            if K > 0:
                base = prob
                for uid in matches:
                    per_ind_weighted[uid] = per_ind_weighted.get(uid, Decimal(0)) + (base * w)
                mlen = len(matches)
                if mlen >= 2:
                    for a in range(mlen):
                        ua = matches[a]
                        for b in range(a + 1, mlen):
                            ub = matches[b]
                            key = (ua, ub) if str(ua) < str(ub) else (ub, ua)
                            per_pair_weighted[key] = per_pair_weighted.get(key, Decimal(0)) + (base * w)
    return out_pairs, per_ind_weighted, per_pair_weighted


def _safe_log10dB(x: Decimal) -> float:
    try:
        xf = float(x)
        if xf <= 0.0:
            return float('-inf')
        return 10.0 * math.log10(xf)
    except Exception:
        return float('-inf')


def process_iteration_results(
        simulate_func: callable,
        simulation_parameters: SimulationParameters,
        progress_bar: ProgressBar,
        reporter: Reporter,
        out_file_name: str,
        number_of_threads: int = 1,
        is_outside: bool = False
) -> (Decimal, List[int], List[Decimal], Dict[str, Decimal]):
    with (progress_bar):
        valid = False
        trial = 1
        tightening = 0
        model_iterations = {m: 0 for m in range(3)}
        model_probabilities = {m: [] for m in range(3)}
        weight_sums = {m: Decimal(0) for m in range(3)}
        weighted_sums = {m: Decimal(0) for m in range(3)}
        since_write = {m: 0 for m in range(3)}
        per_individual_weighted_sums = {m: defaultdict(Decimal) for m in range(3)}
        per_pair_weighted_sums = {m: defaultdict(Decimal) for m in range(3)}
        threads = max(1, int(number_of_threads))
        window = max(1, int(simulation_parameters.stability_window))
        batch_size = min(max(256, window // max(1, threads * 4)), 512)

        while not valid:
            for m in range(3):  # Model validation
                with open(
                        f"{simulation_parameters.results_path}/{out_file_name}_m_{m}_outside_{is_outside}.txt",
                        "a", 1
                ) as f:
                    tasks = []
                    remaining = window
                    start = 0
                    while remaining > 0:
                        c = min(batch_size, remaining)
                        tasks.append((start, c))
                        start += c
                        remaining -= c
                    with multiprocessing.Pool(threads, initializer=_a1_pool_init,
                              initargs=(simulate_func.func, simulate_func.keywords or {})) as pool:
                        print("\nStarting trial {} , model {} with {} threads (batch_size={})...".format(trial, m + 1,
                                                                                                         threads,
                                                                                                         batch_size))
                        for batch_out, per_ind_dict, per_pair_dict in pool.imap_unordered(_a1_run_batch, tasks,
                                                                                          chunksize=1):
                            n_new = len(batch_out)
                            for prob, w in batch_out:
                                model_iterations[m] += 1
                                weight_sums[m] += Decimal(w)
                                weighted_sums[m] += (prob * Decimal(w))
                                current_mean = (weighted_sums[m] / weight_sums[m]) if weight_sums[m] != 0 else Decimal(
                                    0)
                                model_probabilities[m].append(current_mean)
                            for uid, w_sum in per_ind_dict.items():
                                per_individual_weighted_sums[m][uid] += Decimal(w_sum)
                            for key, w_sum in per_pair_dict.items():
                                per_pair_weighted_sums[m][key] += Decimal(w_sum)
                            progress_bar.update(n_new)
                            since_write[m] += n_new
                            while since_write[m] >= 100:
                                f.write(str(model_probabilities[m][-1]) + "\n")
                                since_write[m] -= 100

            total_mean = Decimal(sum(model_probabilities[m][-1] for m in range(3)) / 3)
            number_outside = sum(sum(
                abs(prob - total_mean) / total_mean > simulation_parameters.model_validity_threshold for prob in
                model_probabilities[m]) for m in range(3))
            if all(all(abs(prob - total_mean) / total_mean < simulation_parameters.model_validity_threshold for prob in
                       model_probabilities[m]) for m in range(3)):
                if total_mean > 1.0:
                    reporter.log(
                        ("\nModel is stable after {} trials. However, current mean: {} (log {}) is greater than unity. "
                         "Starting new trial with tightened threshold: {}").format(trial, total_mean,
                                                                                   _safe_log10dB(total_mean), (
                                                                                           simulation_parameters.model_validity_threshold / 2)))
                    if out_file_name == "match_probabilities" and tightening == 2:
                        reporter.log(
                            "\nMaximum number of tightening reached. Model is not valid. Restarting simulation.")
                        raise InvalidAveragePedigreeProbability("Average pedigree probability likely not valid.")
                    simulation_parameters.model_validity_threshold /= 1.5
                    tightening += 1
                else:
                    valid = True
                    reporter.log(
                        f"\nModel is valid after {trial} trials! Current mean: {total_mean} (log {_safe_log10dB(total_mean)})")
            else:
                reporter.log((
                                 "\nModel is not valid after {} trials. Current mean: {} (log {}). {} data points fall outside range. "
                                 "Starting new trial...").format(trial, total_mean, _safe_log10dB(total_mean),
                                                                 number_outside))

            if not valid:
                trial += 1
                model_probabilities = {m: model_probabilities[m][-1:] for m in range(3)}
        try:
            if out_file_name == "match_probabilities" and (not is_outside):
                ctx = simulate_func.keywords or {}
                individuals = ctx.get("individuals", {})
                ordered_unknown_ids = ctx.get("ordered_unknown_ids", [])
                per_model_means = []
                for m in range(3):
                    denom = weight_sums[m] if weight_sums[m] != 0 else Decimal(0)
                    model_means = {}
                    if denom != 0:
                        for uid, w_sum in per_individual_weighted_sums[m].items():
                            model_means[uid] = (w_sum / denom)
                    else:
                        model_means = {uid: Decimal(0) for uid in per_individual_weighted_sums[m].keys()}
                    per_model_means.append(model_means)
                final_per_individual = defaultdict(Decimal)
                all_ids = set().union(*[set(d.keys()) for d in per_model_means]) if per_model_means else set()
                for uid in all_ids:
                    vals = [d.get(uid, Decimal(0)) for d in per_model_means]
                    final_per_individual[uid] = sum(vals) / Decimal(3)
                ts = datetime.now().strftime('%Y%m%d%H%M%S')
                csv_path = "{}/per_individual_marginal_probabilities_{}.csv".format(simulation_parameters.results_path,
                                                                                    ts)
                with open(csv_path, "w", newline="") as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(["individual_id", "individual_name", "probability"])
                    for uid in ordered_unknown_ids:
                        ind = individuals.get(uid)
                        if ind is None or getattr(ind, "exclude", False):
                            continue
                        name = getattr(ind, "name", str(uid))
                        writer.writerow([uid, name, str(final_per_individual.get(uid, Decimal(0)))])
                reporter.log("Saved per-individual marginal probabilities to: {}".format(csv_path))
        except Exception as _csv_err:
            reporter.log("[warn] Could not write per-individual/pairwise CSVs: {}".format(_csv_err))
        return (
            total_mean,
            [int(trial * simulation_parameters.stability_window * 3)],
            [Decimal(model_probabilities[m][-1]) for m in range(3)],
            final_per_individual if (out_file_name == "match_probabilities" and (not is_outside)) else {}
        )


def calculate_average_pedigree_probability(pedigree: Pedigree, root_name: str, marker_set: MarkerSet,
                                           simulation_parameters: SimulationParameters, reporter: Reporter,
                                           is_outside: bool, number_of_threads: int = 1) -> (Decimal, List[int], List[Decimal], Dict[str, Decimal]):
    progress_bar = reporter.progress_bar(desc="Calculating average pedigree probability")
    known = [ind for ind in pedigree.get_level_order_traversal(root_name) if ind.haplotype_class != "unknown"]
    if len(known) == 1:
        reporter.log("\nRoot is only known individual in (reduced) pedigree. Average pedigree probability is 1.")
        return Decimal(1), [0], [Decimal(1), Decimal(1), Decimal(1)], {}
    simulate_func = partial(simulate_pedigree_iteration, root_name=root_name, pedigree=pedigree, marker_set=marker_set,
                            two_step_mutation_factor=simulation_parameters.two_step_mutation_factor, bias_value=simulation_parameters.bias)
    return process_iteration_results(simulate_func, simulation_parameters, progress_bar, reporter,
                                     "average_pedigree_probabilities", number_of_threads, is_outside)


def simulate_matching_haplotypes(pedigree: Pedigree, individuals: Mapping[str, Individual],
                                 relationships: Collection[Relationship], parents: Mapping[str, str],
                                 suspect_haplotype: Haplotype, ordered_unknown_ids: Sequence[str],
                                 marker_set: MarkerSet, average_pedigree_probability: Decimal,
                                 two_step_mutation_factor: float, picking_probabilities: Dict[str, Decimal],
                                 is_outside: bool,
                                 bias_value: float = 0.1
                                 ) -> IterationResult:
    rnd = _rng()
    simulation_probability = Decimal(1)
    total_pick = sum(picking_probabilities.values())
    norm_pick = {k: (v / total_pick) if total_pick != 0 else Decimal(0) for k, v in picking_probabilities.items()}
    pop = list(norm_pick.keys())
    wts = [float(norm_pick[k]) for k in pop]
    fixed_individual_id = rnd.choices(population=pop, weights=wts, k=1)[0]
    simulation_probability *= norm_pick[fixed_individual_id]

    haplotypes = {
        ind.id: (deepcopy(suspect_haplotype) if ind.id == fixed_individual_id else individuals[ind.id].haplotype) for
        ind in individuals.values()}

    _orig = None
    try:
        if fixed_individual_id in individuals:
            _orig = individuals[fixed_individual_id].haplotype_class
            individuals[fixed_individual_id].haplotype_class = "known"
    except Exception:
        _orig = None

    try:
        fh = haplotypes[fixed_individual_id]
        missing = []
        for m in marker_set.markers:
            if not fh.alleles.get(m.name):
                if suspect_haplotype.alleles.get(m.name):
                    fh.alleles[m.name] = deepcopy(suspect_haplotype.alleles[m.name])
                else:
                    missing.append(m.name)
        if missing:
            print("[warn] fixed individual {} missing loci for bias: {}".format(fixed_individual_id, missing))
    except Exception:
        pass

    simulated = set()  # type: ignore
    w_total = 1.0
    u_total = 1.0
    for uid in ordered_unknown_ids:
        if uid == fixed_individual_id:
            continue
        pid = parents[uid]
        parent = haplotypes[pid]
        hap_tuple = tuple(sorted(haplotypes.items(), key=lambda kv: str(kv[0])))
        try:
            biases = pedigree.get_biases(uid, marker_set, hap_tuple, bias_value)
        except ValueError as _e:
            if 'empty locus' in str(_e).lower():
                biases = []
                print("[warn] get_biases empty locus for id={}; continuing without bias".format(uid))
            else:
                raise
        haplotypes[uid], w, u = mutate_haplotype(parent, marker_set, two_step_mutation_factor, biases)
        w_total *= w
        u_total *= u
        simulated.add((pid, uid))

    all_edges = get_edge_probabilities(haplotypes, relationships, marker_set, two_step_mutation_factor, False)
    w_factor = Decimal(0) if w_total == 0 else Decimal(u_total) / Decimal(w_total)
    ped_p = reduce(operator.mul, all_edges.values(), 1)

    sim_edges = []
    for (s, t) in simulated:
        if (s, t) in all_edges:
            sim_edges.append(all_edges[(s, t)])
        elif (t, s) in all_edges:
            sim_edges.append(all_edges[(t, s)])
    simulation_probability = reduce(operator.mul, sim_edges, simulation_probability)

    if average_pedigree_probability == 0:
        cond = Decimal(0)
    else:
        cond = ped_p / average_pedigree_probability

    non_excl_match = [uid for uid in ordered_unknown_ids if
                      (haplotypes[uid] == suspect_haplotype) and (not individuals[uid].exclude)]
    total_match = [uid for uid in ordered_unknown_ids if haplotypes[uid] == suspect_haplotype]
    total_number = len({fixed_individual_id}.union(set(total_match)))
    number_non_excl = len(non_excl_match)

    if not is_outside:
        probability = (cond / (simulation_probability * total_number)) if (
                number_non_excl >= 1 and simulation_probability != 0) else Decimal(0)
    else:
        probability = (cond / simulation_probability) if simulation_probability != 0 else Decimal(0)

    if _orig is not None:
        try:
            individuals[fixed_individual_id].haplotype_class = _orig
        except Exception:
            pass

    return IterationResult(
        probability=probability,
        edge_probabilities=all_edges,
        edge_weight_factor=w_factor,
        mutated_haplotypes=haplotypes,
        fixed_individual_id=fixed_individual_id,
    )


def simulate_iteration(
        i: int,
        pedigree: Pedigree,
        individuals: Mapping[str, Individual],
        relationships: Collection[Relationship],
        parents: Mapping[str, str],
        suspect_haplotype: Haplotype,
        ordered_unknown_ids: Sequence[str],
        marker_set: MarkerSet,
        average_pedigree_probability: Decimal,
        two_step_mutation_factor: float,
        picking_probabilities: dict[str, Decimal],
        is_outside: bool,
        bias_value: float = 0.1
) -> IterationResult:
    """Wrapper to allow multiprocessing without lambda."""
    return simulate_matching_haplotypes(
        pedigree=pedigree,
        individuals=individuals,
        relationships=relationships,
        parents=parents,
        suspect_haplotype=suspect_haplotype,
        ordered_unknown_ids=ordered_unknown_ids,
        marker_set=marker_set,
        average_pedigree_probability=average_pedigree_probability,
        two_step_mutation_factor=two_step_mutation_factor,
        picking_probabilities=picking_probabilities,
        is_outside=is_outside,
        bias_value=bias_value
    )


def calculate_matching_haplotypes(pedigree: Pedigree, marker_set: MarkerSet, root_name: str,
                                  suspect_haplotype: Haplotype, simulation_parameters: SimulationParameters,
                                  average_pedigree_probability: Decimal, reporter: Reporter, is_outside: bool,
                                  number_of_threads: int = 1) -> (Decimal, List[int], List[Decimal], Dict[str, Decimal]):
    progress_bar = reporter.progress_bar(desc="Calculating matching haplotypes")
    individuals = {ind.id: ind for ind in pedigree.individuals}
    parents = {rel.child_id: rel.parent_id for rel in pedigree.relationships}
    ordered_unknown_ids = [ind.id for ind in pedigree.get_level_order_traversal(root_name) if
                           ind.haplotype_class == "unknown"]
    simulate_func = partial(simulate_iteration, pedigree=pedigree, individuals=individuals,
                            relationships=pedigree.relationships, parents=parents, suspect_haplotype=suspect_haplotype,
                            ordered_unknown_ids=ordered_unknown_ids, marker_set=marker_set,
                            average_pedigree_probability=average_pedigree_probability,
                            two_step_mutation_factor=simulation_parameters.two_step_mutation_factor,
                            picking_probabilities=pedigree.picking_probabilities, is_outside=is_outside,
                            bias_value=simulation_parameters.bias)
    return process_iteration_results(simulate_func, simulation_parameters, progress_bar, reporter,
                                     "match_probabilities", number_of_threads, is_outside)


def calculate_proposal_distribution(pedigree: Pedigree, marker_set: MarkerSet, root_name: str,
                                    suspect_haplotype: Haplotype, average_pedigree_probability: Decimal,
                                    simulation_parameters: SimulationParameters, reporter: Reporter) -> (
        Dict[int, Decimal], List[int], List[Decimal], Dict[str, Decimal]):
    unknown_individuals = pedigree.get_unknown_individuals()
    proposal: Dict[int, Decimal] = {}
    needed: List[int] = []
    model_probs: List[Decimal] = []
    model_valid = False
    used_probs: List[Decimal] = []
    non_excl = len([ind for ind in unknown_individuals if not ind.exclude])
    if non_excl == 0:
        reporter.log(
            "\nWarning! No (non-excluded) unknown individuals left in the pedigree. Only calculating outside match probability.")
        proposal[0] = Decimal(1)
        proposal[1] = Decimal(0)
        return proposal, needed, model_probs, model_valid, used_probs

    max_threads = multiprocessing.cpu_count()
    n_threads = min(simulation_parameters.number_of_threads, max_threads)
    sim_params_inside = simulation_parameters
    proposal[1], needed, model_probs, per_ind_probs = calculate_matching_haplotypes(pedigree, marker_set,
                                                                                              root_name,
                                                                                              suspect_haplotype,
                                                                                              sim_params_inside,
                                                                                              average_pedigree_probability,
                                                                                              reporter,
                                                                                              is_outside=False,
                                                                                              number_of_threads=n_threads)
    proposal[0] = Decimal(1) - proposal[1]
    if proposal[0] < 0:
        reporter.log("\nWarning! Proposal distribution does not add up to unity.")
    return proposal, needed, model_probs, per_ind_probs


def calculate_outside_match_probability(pedigree: Pedigree, marker_set: MarkerSet, root_name: str,
                                        suspect_haplotype: Haplotype, average_pedigree_probability: Decimal,
                                        simulation_parameters: SimulationParameters, reporter: Reporter,
                                        number_of_threads: int = 1):
    return calculate_matching_haplotypes(pedigree, marker_set, root_name, suspect_haplotype, simulation_parameters,
                                         average_pedigree_probability, reporter, is_outside=True,
                                         number_of_threads=number_of_threads)


def run_simulation(
        input_pedigree: Pedigree,
        suspect_name: str,
        marker_set: MarkerSet,
        simulation_parameters: SimulationParameters,
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
    with open(f"{simulation_parameters.results_path}/original_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf",
              "wb") as f:
        f.write(pedigree_bytes_data)

    save_pedigree_to_png(pedigree=pedigree,
                         global_config=global_config,
                         results_path=simulation_parameters.results_path,
                         pedigree_name="original_pedigree")

    extended_pedigree = deepcopy(pedigree)

    # Extend the pedigree and store the name of the last added individual (used for edge cases)
    last_child_name = extended_pedigree.extend_pedigree()
    extended_pedigree_bytes_data = extended_pedigree.write_to_tgf()
    with open(f"{simulation_parameters.results_path}/extended_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf",
              "wb") as f:
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
    with open(
            f"{simulation_parameters.results_path}/rerooted_original_pedigree_{datetime.now().strftime('%Y%m%d%H%M%S')}.tgf",
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
            else:
                average_pedigree_probability, average_pedigree_needed_iterations, average_pedigree_model_pedigree_probabilities, _ = calculate_average_pedigree_probability(
                    pedigree=pedigree,
                    root_name=root_name,
                    marker_set=marker_set,
                    simulation_parameters=simulation_parameters,
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

                    per_individual_probabilities={},

                    run_time_pedigree_probability=run_time_pedigree_probability,
                    run_time_proposal_distribution=timedelta(0),
                    run_time_extended_average_pedigree_probability=timedelta(0),
                    run_time_outside_match_probability=timedelta(0),
                    total_run_time=run_time_pedigree_probability,
                )

                with open(
                        f"{simulation_parameters.results_path}/simulation_result_{datetime.now().strftime('%Y%m%d%H%M%S')}.pkl",
                        "wb") as f:
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
                inside_needed_iterations = [0]
                inside_model_probabilities = [Decimal(0), Decimal(0), Decimal(0)]
            else:
                inside_match_probability, inside_needed_iterations, inside_model_probabilities, inside_model_per_individual = calculate_proposal_distribution(
                    pedigree=pedigree,
                    marker_set=marker_set,
                    root_name=root_name,
                    suspect_haplotype=suspect_haplotype,
                    average_pedigree_probability=average_pedigree_probability,
                    simulation_parameters=simulation_parameters,
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

    if skip_outside:
        reporter.log("\nSkipping extended average pedigree probability calculation.")
        extended_average_pedigree_probability = Decimal(1)
        outside_match_probability = Decimal(0)
        extended_needed_iterations = [0]
        outside_needed_iterations = [0]
        extended_model_pedigree_probabilities = [Decimal(1), Decimal(1), Decimal(1)]
        outside_model_probabilities = [Decimal(0), Decimal(0), Decimal(0)]
        run_time_extended_average_pedigree_probability = timedelta(0)
        run_time_outside_match_probability = timedelta(0)
    else:
        while True:  # loop to allow retrying in case of BadInputError
            try:
                extended_average_pedigree_probability, extended_needed_iterations, extended_model_pedigree_probabilities, _  = calculate_average_pedigree_probability(
                    pedigree=extended_pedigree,
                    root_name=extended_root_name,
                    marker_set=marker_set,
                    simulation_parameters=simulation_parameters,
                    reporter=reporter,
                    is_outside=True,
                    number_of_threads=number_of_threads,
                )
                run_time_extended_average_pedigree_probability = datetime.now() - start_time_extended_average_pedigree_probability
                start_time_outside_match_probability = datetime.now()

                outside_match_probability, outside_needed_iterations, outside_model_probabilities, _ = calculate_outside_match_probability(
                    pedigree=extended_pedigree,
                    marker_set=marker_set,
                    root_name=extended_root_name,
                    suspect_haplotype=suspect_haplotype,
                    average_pedigree_probability=extended_average_pedigree_probability,
                    simulation_parameters=simulation_parameters,
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

        per_individual_probabilities=inside_model_per_individual if not skip_inside else {},

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
