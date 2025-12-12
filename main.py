from argparse import ArgumentParser
from pathlib import Path
import logging
from pedigree_lr.config import load_config
from pedigree_lr.data import load_marker_set_from_config, load_pedigree_from_config, load_trace_from_file
from pedigree_lr.reporting import ConsoleReporter
from pedigree_lr.simulation import run_simulation

logger = logging.getLogger(__name__)

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


def simulate(
        config_path: str = "config.ini",
        skip_inside: bool = False,
        skip_outside: bool = False,
        trace_mode: bool = False,
        adaptive_bias: bool = False
) -> None:
    """
    Simulates a pedigree-based likelihood ratio analysis.

    Args:
        config_path (str): Path to the configuration INI file. Defaults to "config.ini".
        skip_inside (bool): If True, skips inside pedigree probabilities. Defaults to False.
        skip_outside (bool): If True, skips outside pedigree probabilities. Defaults to False.
        trace_mode (bool): If True, uses uniform picking probabilities. Defaults to False.

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

    # If adaptive bias is enabled, set bias to None to trigger adaptive mode
    if adaptive_bias:
        config.simulation_parameters.bias = None

    marker_set = load_marker_set_from_config(config)
    pedigree = load_pedigree_from_config(config, marker_set)
    pedigree.exclude_individuals(config.exclude_individuals)

    # Priority: JSON trace > CSV trace > config suspect (auto-detection)
    trace = None
    suspect_name = config.suspect

    if hasattr(pedigree, '_trace_from_json') and pedigree._trace_from_json is not None:
        # TRACE found in JSON file
        trace = pedigree._trace_from_json
        trace_mode = True  # Auto-enable trace mode
        suspect_name = None  # No suspect in trace mode
        logger.info("TRACE profile detected in JSON file. Enabling trace mode.")
        if config.trace:
            logger.warning("Both JSON TRACE and CSV trace file specified. Using JSON TRACE (CSV ignored).")
    elif config.trace:
        # CSV trace file specified
        trace = load_trace_from_file(config.trace, marker_set)
        trace_mode = True  # Auto-enable trace mode
        suspect_name = None  # No suspect in trace mode
        logger.info("Loading trace from CSV file (legacy mode). Enabling trace mode.")
    elif suspect_name is None:
        # No trace and no suspect specified - error
        logger.error("No suspect or trace profile specified. Please provide either a suspect name or a TRACE profile.")
        raise ValueError("No suspect or trace profile specified in configuration.")

    reporter = ConsoleReporter()

    simulation_result = run_simulation(
        input_pedigree=pedigree,
        suspect_name=suspect_name,
        trace=trace,
        marker_set=marker_set,
        simulation_parameters=config.simulation_parameters,
        reporter=reporter,
        skip_inside=skip_inside,
        skip_outside=skip_outside,
        trace_mode=trace_mode
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
    parser.add_argument(
        "-i",
        "--skip-inside",
        action="store_true",
        default=False,
        required=False,
        help='Skip inside pedigree probabilities.'
    )
    parser.add_argument(
        "-o",
        "--skip-outside",
        action="store_true",
        default=False,
        required=False,
        help='Skip outside pedigree probabilities.'
    )

    parser.add_argument(
        "-t",
        "--trace-mode",
        action="store_true",
        default=False,
        required=False,
        help='Use uniform picking probabilities.'
    )

    parser.add_argument(
        "-a",
        "--adaptive-bias",
        action="store_true",
        default=False,
        required=False,
        help='Enable adaptive bias mode (dynamically adjusts bias per model based on performance).'
    )

    args = parser.parse_args()

    simulate(args.config_path,
             skip_inside=args.skip_inside,
             skip_outside=args.skip_outside,
             trace_mode=args.trace_mode,
             adaptive_bias=args.adaptive_bias
             )
