from argparse import ArgumentParser
from pathlib import Path
from random import Random
from pedigree_lr.config import load_config
from pedigree_lr.data import load_marker_set_from_config, load_pedigree_from_config
from pedigree_lr.reporting import ConsoleReporter
from pedigree_lr.simulation import run_simulation

"""
This module provides functionality to simulate a pedigree-based likelihood ratio (LR) analysis
using configurations specified in an INI file. It includes functions for loading configurations,
setting up pedigrees, and running simulations.

Dependencies:
- argparse: For parsing command-line arguments.
- pathlib: For handling file paths.
- random: For generating random numbers.
- pedigree_lr.config: For loading configuration files.
- pedigree_lr.data: For loading marker sets and pedigrees.
- pedigree_lr.reporting: For reporting simulation results.
- pedigree_lr.simulation: For running the simulation.

Usage:
Run this script directly to execute the simulation with a specified configuration file.
"""


def simulate(config_path: str = "config.ini"):
    """
    Simulates a pedigree-based likelihood ratio analysis.

    Args:
        config_path (str): Path to the configuration INI file. Defaults to "config.ini".

    The function performs the following steps:
    1. Loads the configuration file.
    2. Sets up the marker set and pedigree based on the configuration.
    3. Excludes individuals from the pedigree as specified in the configuration.
    4. Initializes a console reporter for simulation results.
    5. Runs the simulation using the provided parameters and random seed.

    Returns:
        None
    """
    config = load_config(Path(config_path))
    marker_set = load_marker_set_from_config(config)
    pedigree = load_pedigree_from_config(config, marker_set)
    pedigree.exclude_individuals(config.exclude_individuals)

    reporter = ConsoleReporter()

    simulation_result = run_simulation(
        input_pedigree=pedigree,
        suspect_name=config.suspect,
        marker_set=marker_set,
        simulation_parameters=config.simulation_parameters,
        reporter=reporter,
        random_seed=config.random_seed
    )


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
    args = parser.parse_args()

    simulate(args.config_path)
