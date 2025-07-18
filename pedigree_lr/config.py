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
    suspect: str
    marker_set: Path
    known_haplotypes: Path
    simulation_parameters: SimulationParameters
    random_seed: int | None
    exclude_individuals: list[str] = None


def load_config(path: Path) -> Config:
    config = ConfigParser()
    config.optionxform = str  # type: ignore
    config.read(path)

    folder_name = f"{config['pedigree']['simulation_name'].replace(' ', '_').lower()}_{datetime.now().strftime('%Y%m%d%H%M%S')}"
    results_path = Path(config["pedigree"]["results_path"]).resolve() / folder_name
    results_path.mkdir(parents=True, exist_ok=True)

    return Config(
        pedigree=Path(config["pedigree"]["path"]),
        suspect=config["pedigree"]["suspect"],
        marker_set=Path(config["pedigree"]["marker_set"]),
        known_haplotypes=Path(config["pedigree"]["known_haplotypes"]),
        simulation_parameters=SimulationParameters(
            max_number_of_iterations=int(config["pedigree"]["number_of_iterations"]),
            two_step_mutation_factor=float(config["pedigree"]["two_step_mutation_factor"]),
            stability_window=int(config["pedigree"]["stability_window"]),
            stability_min_iterations=int(config["pedigree"]["stability_min_iterations"]),
            stability_threshold=float(config["pedigree"]["stability_threshold"]),
            model_validity_threshold=float(config["pedigree"]["model_validity_threshold"]),
            simulation_name=str(config["pedigree"]["simulation_name"]).replace(' ', '_').lower(),
            number_of_threads=int(config["pedigree"]["number_of_threads"]),
            results_path=results_path,
            random_seed=int(config["pedigree"]["random_seed"]) if config["pedigree"]["random_seed"] else None,
            user_name="admin",
        ),
        random_seed=int(config["pedigree"]["random_seed"]) if config["pedigree"]["random_seed"] else None,
        exclude_individuals=config["pedigree"]["exclude_individuals"].split(","),
    )