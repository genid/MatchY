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

    folder_name = f"{config['pedigree']['simulation_name'].replace(' ', '_').lower()}_{datetime.now().strftime('%Y%m%d%H%M%S')}"
    results_path = Path(config["pedigree"]["results_path"]).resolve() / folder_name
    results_path.mkdir(parents=True, exist_ok=True)

    if "bias" in config["pedigree"]:
        bias = float(config["pedigree"]["bias"])
        if not (0.0 <= bias <= 0.5):
            raise ValueError("Bias must be between 0 and 0.5.")
    else:
        bias = None

    if "trace" in config["pedigree"]:
        trace = config["pedigree"]["trace"]
    else:
        trace = None

    return Config(
        pedigree=Path(config["pedigree"]["path"]),
        suspect=config["pedigree"]["suspect"],
        marker_set=Path(config["pedigree"]["marker_set"]),
        known_haplotypes=Path(config["pedigree"]["known_haplotypes"]),
        trace=Path(trace) if trace else None,
        simulation_parameters=SimulationParameters(
            two_step_mutation_factor=float(config["pedigree"]["two_step_mutation_factor"]),
            stability_window=int(config["pedigree"]["stability_window"]),
            model_validity_threshold=float(config["pedigree"]["model_validity_threshold"]),
            bias=bias,
            simulation_name=str(config["pedigree"]["simulation_name"]).replace(' ', '_').lower(),
            number_of_threads=int(config["pedigree"]["number_of_threads"]),
            results_path=results_path,
            user_name=str(config["pedigree"]["user_name"]).replace(' ', '_').lower(),
        ),
        exclude_individuals=config["pedigree"]["exclude_individuals"].split(","),
    )