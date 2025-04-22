from __future__ import annotations

from configparser import ConfigParser
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

_CONFIG_LOCATION = Path("config.yaml")


@dataclass(frozen=True)
class Config:
    pedigree: Path
    suspect: str
    marker_set: Path
    known_haplotypes: Path
    simulation_parameters: Mapping[str, any]
    random_seed: int | None
    exclude_individuals: list[str] = None


def load_config(path: Path) -> Config:
    config = ConfigParser()
    config.optionxform = str  # type: ignore
    config.read(path)

    return Config(
        pedigree=Path(config["pedigree"]["path"]),
        suspect=config["pedigree"]["suspect"],
        marker_set=Path(config["pedigree"]["marker_set"]),
        known_haplotypes=Path(config["pedigree"]["known_haplotypes"]),
        simulation_parameters={
            "number_of_iterations": int(config["pedigree"]["number_of_iterations"]),
            "two_step_mutation_factor": float(config["pedigree"]["two_step_mutation_factor"]),
            "stability_window": int(config["pedigree"]["stability_window"]),
            "stability_min_iterations": int(config["pedigree"]["stability_min_iterations"]),
            "stability_threshold": float(config["pedigree"]["stability_threshold"]),
            "model_validity_threshold": float(config["pedigree"]["model_validity_threshold"]),
            "simulation_name": str(config["pedigree"]["simulation_name"]),
        },
        random_seed=int(config["pedigree"]["random_seed"]) if config["pedigree"]["random_seed"] else None,
        exclude_individuals=config["pedigree"]["exclude_individuals"].split(","),
    )
