from argparse import ArgumentParser, Namespace
from decimal import Decimal
from pathlib import Path
from pedigree_lr.config import load_config
from pedigree_lr.data import load_marker_set_from_config
import json
from main import simulate
import pickle
from datetime import timedelta
from pedigree_lr.models import SimulationResult


def split_configs(args: Namespace):
    config_path = args.config_path
    folder = Path(config_path).resolve().parent
    per_marker_folder = folder / f"{Path(config_path).stem}_per_marker"
    per_marker_folder.mkdir(exist_ok=True)

    config = load_config(Path(config_path))
    marker_set = load_marker_set_from_config(config)

    unmutated_markers = {}

    for marker in marker_set.markers:
        new_marker_set_file = per_marker_folder / f"{marker.name}.csv"
        with open(new_marker_set_file, "w") as file:
            file.write("marker,mutation_rate\n")
            file.write(f"{marker.name},{marker.mutation_rate}\n")

        new_known_haplotypes_file = per_marker_folder / f"{marker.name}_known_haplotypes.json"
        with open(config.known_haplotypes) as infile, open(new_known_haplotypes_file, "w") as outfile:
            known_haplotypes = json.load(infile)
            filtered_haplotypes = {ind: {marker.name: hap for m, hap in markers.items() if m == marker.name}
                                    for ind, markers in known_haplotypes.items()}
            seen_haps = set()
            for ind in list(filtered_haplotypes.keys()):
                hap = filtered_haplotypes[ind][marker.name]
                seen_haps.add(hap)
            if len(seen_haps) == 1:
                if args.combine_unmutated:
                    for ind in filtered_haplotypes.keys():
                        if ind not in unmutated_markers:
                            unmutated_markers[ind] = filtered_haplotypes[ind]
                        else:
                            unmutated_markers[ind][marker.name] = filtered_haplotypes[ind][marker.name]

                if args.only_mutated or args.combine_unmutated:
                    print(f"Skipping marker {marker.name} as it has no variation in known haplotypes.")
                    continue

            json.dump(filtered_haplotypes, outfile, indent=4)

        with open(per_marker_folder / f"{marker.name}.ini", "w") as file:
            file.write(f"[pedigree]\n")
            file.write(f"path = {config.pedigree}\n")
            file.write(f"suspect = {config.suspect}\n")
            file.write(f"marker_set = {new_marker_set_file}\n")
            file.write(f"number_of_iterations = {config.simulation_parameters.max_number_of_iterations}\n")
            file.write(f"two_step_mutation_factor = {config.simulation_parameters.two_step_mutation_factor}\n")
            file.write(f"stability_window = {config.simulation_parameters.stability_window}\n")
            file.write(f"stability_min_iterations = {config.simulation_parameters.stability_min_iterations}\n")
            file.write(f"stability_threshold = {config.simulation_parameters.stability_threshold}\n")
            file.write(f"model_validity_threshold = {config.simulation_parameters.model_validity_threshold}\n")
            file.write(f"random_seed = {config.random_seed if config.random_seed is not None else ''}\n")
            file.write(f"known_haplotypes = {new_known_haplotypes_file}\n")
            exclude_individuals = ",".join(config.exclude_individuals) if config.exclude_individuals else ""
            file.write(f"exclude_individuals = {exclude_individuals}\n")
            file.write(f"simulation_name = {config.simulation_parameters.simulation_name}_{marker.name}\n")
            file.write(f"number_of_threads = {config.simulation_parameters.number_of_threads}\n")
            file.write(f"results_path = {config.simulation_parameters.results_path}\n")

        simulate(str(per_marker_folder / f"{marker.name}.ini"),
                 skip_inside=args.skip_inside,
                 skip_outside=args.skip_outside
                 )

    if args.combine_unmutated and unmutated_markers:
        new_marker_set_file = per_marker_folder / f"combined_unmutated.csv"
        with open(new_marker_set_file, "w") as file:
            file.write("marker,mutation_rate\n")
            for marker in marker_set.markers:
                if all(ind in unmutated_markers and marker.name in unmutated_markers[ind] for ind in unmutated_markers):
                    file.write(f"{marker.name},{marker.mutation_rate}\n")

        new_known_haplotypes_file = per_marker_folder / f"combined_unmutated_known_haplotypes.json"
        with open(new_known_haplotypes_file, "w") as outfile:
            json.dump(unmutated_markers, outfile, indent=4)

        with open(per_marker_folder / f"combined_unmutated.ini", "w") as file:
            file.write(f"[pedigree]\n")
            file.write(f"path = {config.pedigree}\n")
            file.write(f"suspect = {config.suspect}\n")
            file.write(f"marker_set = {new_marker_set_file}\n")
            file.write(f"number_of_iterations = {config.simulation_parameters.max_number_of_iterations}\n")
            file.write(f"two_step_mutation_factor = {config.simulation_parameters.two_step_mutation_factor}\n")
            file.write(f"stability_window = {config.simulation_parameters.stability_window}\n")
            file.write(f"stability_min_iterations = {config.simulation_parameters.stability_min_iterations}\n")
            file.write(f"stability_threshold = {config.simulation_parameters.stability_threshold}\n")
            file.write(f"model_validity_threshold = {config.simulation_parameters.model_validity_threshold}\n")
            file.write(f"random_seed = {config.random_seed if config.random_seed is not None else ''}\n")
            file.write(f"known_haplotypes = {new_known_haplotypes_file}\n")
            exclude_individuals = ",".join(config.exclude_individuals) if config.exclude_individuals else ""
            file.write(f"exclude_individuals = {exclude_individuals}\n")
            file.write(f"simulation_name = {config.simulation_parameters.simulation_name}_combined_unmutated\n")
            file.write(f"number_of_threads = {config.simulation_parameters.number_of_threads}\n")
            file.write(f"results_path = {config.simulation_parameters.results_path}\n")

        simulate(str(per_marker_folder / f"combined_unmutated.ini"),
                 skip_inside=args.skip_inside,
                 skip_outside=args.skip_outside
                 )

    combined_inside_match_probability = Decimal(1)
    combined_outside_match_probability = Decimal(1)
    combined_needed_iterations = 0
    combined_run_time = timedelta()

    results_glob = Path(config.simulation_parameters.results_path).rglob("simulation_result_*.pkl")

    for result_file in results_glob:
        if result_file.exists():
            with open(result_file, 'rb') as f:
                result: SimulationResult = pickle.load(f)
                combined_inside_match_probability *= Decimal(result.inside_match_probability[1])
                combined_outside_match_probability *= Decimal(result.outside_match_probability)
                combined_needed_iterations += sum(result.average_pedigree_needed_iterations) + sum(
                    result.extended_needed_iterations) + sum(result.inside_needed_iterations) + sum(
                    result.outside_needed_iterations)
                combined_run_time += result.total_run_time

    with open(config.simulation_parameters.results_path / "combined_results.txt", "w") as f:
        f.write(f"{combined_inside_match_probability}\n")
        f.write(f"{combined_outside_match_probability}\n")
        f.write(f"{combined_needed_iterations}\n")
        f.write(f"{combined_run_time}\n")

if __name__ == '__main__':
    """
    Entry point for the script. Parses command-line arguments and runs the simulation.

    Command-line arguments:
    -c, --config-path: Path to the configuration INI file. Defaults to "config.ini".
    """
    parser = ArgumentParser(description='Run the pedigreeLR application')
    parser.add_argument(
        "-c",
        '--config-path',
        default="config.ini",
        type=str,
        required=False,
        help='The path to the config ini file.'
    )
    parser.add_argument(
        "-i",
        "--skip-inside",
        action="store_true",
        help='Skip inside pedigree probabilities.'
    )
    parser.add_argument(
        "-o",
        "--skip-outside",
        action="store_true",
        help='Skip outside pedigree probabilities.'
    )

    parser.add_argument(
        "-m",
        "--only-mutated",
        action="store_true",
        help='Only compute probabilities for mutated markers.'
    )

    parser.add_argument(
        "-x",
        "--combine-unmutated",
        action="store_true",
        help='Combine unmutated markers into a single analysis.'
    )

    args = parser.parse_args()

    split_configs(args)
