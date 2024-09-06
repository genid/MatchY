import json
from argparse import ArgumentParser
from pathlib import Path
from pprint import pprint
from random import Random

from pedigree_lr.config import load_config
from pedigree_lr.data import load_marker_set, load_pedigree
from pedigree_lr.reporting import ConsoleReporter
from pedigree_lr.simulation import run_simulation


def simulate(
    config_path: str = "config.ini"
):
    config = load_config(Path(config_path))
    marker_set = load_marker_set(config)
    pedigree = load_pedigree(config, marker_set)

    # pedigree.print()

    reporter = ConsoleReporter()

    result = run_simulation(
        pedigree=pedigree,
        suspect_name=config.suspect,
        marker_set=marker_set,
        number_of_iterations=config.number_of_iterations,
        random=Random(config.random_seed),
        reporter=reporter
    )

    pprint(result)


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
