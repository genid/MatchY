from io import StringIO
from pedigree_lr.config import Config
from pedigree_lr.models import MarkerSet, Pedigree
import json
from pathlib import Path


def load_marker_set_from_config(config: Config) -> MarkerSet:
    marker_set = MarkerSet()
    with open(config.marker_set) as file:
        marker_set.read_marker_set_from_file(file)
    return marker_set


def get_marker_set_names() -> list[str]:
    kits_path = Path(__file__).resolve().parent.parent / "data" / "kits.json"
    with open(kits_path, "r") as file:
        kits = json.load(file)
    return list(kits.keys())


def load_marker_set_from_database(marker_set_name: str) -> MarkerSet:
    kits_path = Path(__file__).resolve().parent.parent / "data" / "kits.json"
    with open(kits_path, "r") as file:
        kits = json.load(file)

    if marker_set_name not in kits:
        raise ValueError(f"Marker set {marker_set_name} not found in database.")

    markers = kits[marker_set_name]

    marker_set = MarkerSet()
    marker_set.load_markers_from_database(markers)
    return marker_set

def load_marker_set_from_upload(marker_set_file) -> MarkerSet:
    marker_set = MarkerSet()
    marker_set.read_marker_set_from_file(marker_set_file)
    return marker_set


def load_pedigree_from_config(config: Config, marker_set: MarkerSet) -> Pedigree:
    pedigree = Pedigree()
    file_extension = config.pedigree.suffix
    with open(config.pedigree, "r") as file:
        pedigree.read_pedigree_from_file(file=file,
                                         file_extension=file_extension) # TODO: resolve pedigree load from config

    with open(config.known_haplotypes) as file:
        pedigree.read_known_haplotypes_from_file(file=file,
                                                 marker_set=marker_set)

    pedigree.check_known_haplotypes()
    pedigree.reroot_pedigree(config.suspect)
    return pedigree


def load_pedigree_from_upload(pedigree_file: StringIO,
                              file_extension: str) -> Pedigree:
    pedigree = Pedigree()
    pedigree.read_pedigree_from_file(pedigree_file, file_extension)
    pedigree.check_pedigree_structure()
    return pedigree
