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
    known_haplotypes: Mapping[str, Path]
    number_of_iterations: int
    random_seed: int | None
    show_simulated_pedigrees: bool


def load_config(path: Path) -> Config:
    config = ConfigParser()
    config.optionxform = str  # type: ignore
    config.read(path)

    return Config(
        pedigree=Path(config["pedigree"]["path"]),
        suspect=config["pedigree"]["suspect"],
        marker_set=Path(config["pedigree"]["marker_set"]),
        known_haplotypes={
            name: Path(path)
            for name, path in config["pedigree.known_haplotypes"].items()
        },
        number_of_iterations=int(config["pedigree"]["number_of_iterations"]),
        random_seed=int(config["pedigree"]["random_seed"])
        if config["pedigree"]["random_seed"]
        else None,
        show_simulated_pedigrees=config["pedigree"]["show_simulated_pedigrees"].lower()
        == "true",
    )
