from __future__ import annotations
from configparser import ConfigParser
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from pedigree_lr.models import SimulationParameters

_CONFIG_LOCATION = Path("config.yaml")


@dataclass(frozen=True)
class Config:
    pedigree: Path
    suspect: str | None
    trace: Path | None
    marker_set: Path
    known_haplotypes: Path
    simulation_parameters: SimulationParameters
    exclude_individuals: list[str] = None


def load_config(path: Path) -> Config:
    config = ConfigParser()
    config.optionxform = str  # type: ignore
    config.read(path)

    # Check if results_path already contains a datetime folder (from GUI)
    # Pattern: ends with _YYYYMMDDHHMMSS
    import re
    results_path_config = Path(config["pedigree"]["results_path"]).resolve()

    # Check if path ends with datetime pattern (YYYYMMDDHHMMSS)
    if re.search(r'_\d{14}$', results_path_config.name):
        # Path already includes datetime folder (created by GUI)
        results_path = results_path_config
    else:
        # Path is parent directory, create subfolder with datetime (CLI mode)
        folder_name = f"{config['pedigree']['simulation_name'].replace(' ', '_').lower()}_{datetime.now().strftime('%Y%m%d%H%M%S')}"
        results_path = results_path_config / folder_name

    results_path.mkdir(parents=True, exist_ok=True)

    if "bias" in config["pedigree"]:
        bias = float(config["pedigree"]["bias"])
        if not (0.0 <= bias <= 0.5):
            raise ValueError("Bias must be between 0 and 0.5.")
    else:
        bias = None

    # Both suspect and trace are now optional
    trace = config["pedigree"].get("trace", None)
    suspect = config["pedigree"].get("suspect", None)

    return Config(
        pedigree=Path(config["pedigree"]["path"]),
        suspect=suspect,
        marker_set=Path(config["pedigree"]["marker_set"]),
        known_haplotypes=Path(config["pedigree"]["known_haplotypes"]),
        trace=Path(trace) if trace else None,
        simulation_parameters=SimulationParameters(
            two_step_mutation_factor=float(config["pedigree"]["two_step_mutation_fraction"]),
            stability_window=int(config["pedigree"]["batch_length"]),
            model_validity_threshold=float(config["pedigree"]["convergence_criterion"]),
            bias=bias,
            simulation_name=str(config["pedigree"]["simulation_name"]),
            number_of_threads=int(config["pedigree"]["number_of_threads"]),
            results_path=results_path,
            user_name=str(config["pedigree"]["user_name"]),
        ),
        exclude_individuals=config["pedigree"]["exclude_individuals"].split(","),
    )