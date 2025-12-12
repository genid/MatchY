"""
Generate config.ini files from Streamlit session state for CLI execution.
"""
from configparser import ConfigParser
from pathlib import Path
import subprocess
import sys
import threading
from typing import Optional, Callable
from pedigree_lr.models import SimulationParameters, MarkerSet


def generate_config_from_streamlit(
    simulation_name: str,
    user_name: str,
    pedigree_file_path: str,
    haplotypes_file_path: str,
    marker_set: MarkerSet,
    suspect_name: Optional[str],
    excluded_individuals: list[str],
    simulation_parameters: SimulationParameters,
    output_dir: Path,
    trace_mode: bool = False,
) -> Path:
    """
    Generate a config.ini file for CLI execution based on Streamlit settings.

    Args:
        simulation_name: Name for the simulation
        user_name: User's name for report
        pedigree_file_path: Path to saved pedigree file
        haplotypes_file_path: Path to saved haplotypes JSON
        marker_set: MarkerSet object to save as CSV
        suspect_name: Name of suspect individual (None for trace mode)
        excluded_individuals: List of excluded individual names
        simulation_parameters: SimulationParameters object
        output_dir: Directory to save results (parent, not simulation-specific)
        trace_mode: Whether trace mode is enabled

    Returns:
        Path to generated config.ini file
    """
    config = ConfigParser()
    config.optionxform = str  # Preserve case

    # Save marker set to temporary CSV file
    temp_dir = Path(output_dir) / "temp_files"
    temp_dir.mkdir(exist_ok=True, parents=True)
    marker_set_path = temp_dir / "marker_set.csv"

    with open(marker_set_path, "w") as f:
        f.write("marker,mutation_rate\n")
        for marker in marker_set.markers:
            f.write(f"{marker.name},{marker.mutation_rate}\n")

    # Get parent directory for results_path (load_config will create simulation-specific folder)
    results_parent = Path(output_dir).parent

    # Pedigree section
    config["pedigree"] = {
        "path": pedigree_file_path,
        "known_haplotypes": haplotypes_file_path,
        "marker_set": str(marker_set_path),
        "simulation_name": simulation_name,
        "user_name": user_name,
        "exclude_individuals": ",".join(excluded_individuals) if excluded_individuals else "",
        "two_step_mutation_factor": str(simulation_parameters.two_step_mutation_factor),
        "stability_window": str(simulation_parameters.stability_window),
        "model_validity_threshold": str(simulation_parameters.model_validity_threshold),
        "number_of_threads": str(simulation_parameters.number_of_threads),
        "results_path": str(results_parent),
    }

    # Add suspect or trace based on mode
    if trace_mode:
        # In trace mode, suspect is not needed (TRACE from JSON)
        pass
    elif suspect_name:
        config["pedigree"]["suspect"] = suspect_name

    # Add bias if specified (not adaptive)
    if simulation_parameters.bias is not None:
        config["pedigree"]["bias"] = str(simulation_parameters.bias)

    # Save config file
    config_path = output_dir / "streamlit_generated_config.ini"
    with open(config_path, "w") as f:
        config.write(f)

    return config_path


def run_cli_subprocess(
    config_path: Path,
    skip_inside: bool = False,
    skip_outside: bool = False,
    trace_mode: bool = False,
    adaptive_bias: bool = False,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> tuple[int, str, str]:
    """
    Run main.py CLI in a subprocess with the generated config.

    Args:
        config_path: Path to config.ini file
        skip_inside: Skip inside pedigree probabilities
        skip_outside: Skip outside pedigree probabilities
        trace_mode: Enable trace mode
        adaptive_bias: Enable adaptive bias mode
        progress_callback: Optional callback function for progress updates

    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    # Build command
    cmd = [
        sys.executable,  # Use same Python interpreter
        "main.py",
        "--config-path", str(config_path),
    ]

    if skip_inside:
        cmd.append("--skip-inside")
    if skip_outside:
        cmd.append("--skip-outside")
    if trace_mode:
        cmd.append("--trace-mode")
    if adaptive_bias:
        cmd.append("--adaptive-bias")

    # Run subprocess with real-time output capture
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,  # Line buffered
        universal_newlines=True
    )

    # Capture output in real-time
    stdout_lines = []
    stderr_lines = []

    def read_stdout():
        for line in iter(process.stdout.readline, ''):
            if line:
                stdout_lines.append(line)
                if progress_callback:
                    progress_callback(line.strip())

    def read_stderr():
        for line in iter(process.stderr.readline, ''):
            if line:
                stderr_lines.append(line)

    # Start threads to read output
    stdout_thread = threading.Thread(target=read_stdout)
    stderr_thread = threading.Thread(target=read_stderr)
    stdout_thread.start()
    stderr_thread.start()

    # Wait for completion
    return_code = process.wait()

    # Wait for output threads to finish
    stdout_thread.join()
    stderr_thread.join()

    return return_code, ''.join(stdout_lines), ''.join(stderr_lines)
