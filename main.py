from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path
from pprint import pprint
from random import Random

from pedigree.config import load_config
from pedigree.data import load_marker_set, load_pedigree
from pedigree.simulation import run_simulation

from tqdm import tqdm


def simulate(
    config_path: str = "config.ini"
):
    config = load_config(Path(config_path))
    marker_set = load_marker_set(config)
    pedigree = load_pedigree(config, marker_set)

    # pedigree.print()

    start = datetime.now()

    with tqdm(total=config.number_of_iterations, desc="Simulating") as bar:
        result = run_simulation(
            pedigree=pedigree,
            marker_set=marker_set,
            suspect=config.suspect,
            number_of_iterations=config.number_of_iterations,
            random=Random(config.random_seed),
            show_progress=lambda count, total: bar.update(1)
        )

    print(f"total_l_matches_normalized: {result.total_l_matches_normalized}")
    print(
        f"Total correct: "
        f"{result.total_correct}/{config.number_of_iterations}="
        f"{result.total_correct / config.number_of_iterations}"
    )
    lr = result.total_0_count / result.total_not_0_count if result.total_not_0_count else "NaN"
    print(
        f"Likelihood ratio: {result.total_0_count}/{result.total_not_0_count}={lr}"
    )

    print(datetime.now() - start)


if __name__ == '__main__':
    parser = ArgumentParser(description='Run the pedigreeLR application')
    parser.add_argument(
        "-c",
        '--config-path',
        default="config.ini",
        type=str,
        required=False,
        help='The path to the config ini file.'
    )
    args = parser.parse_args()

    simulate(args.config_path)
