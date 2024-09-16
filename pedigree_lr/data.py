from pedigree_lr.config import Config
from pedigree_lr.models import MarkerSet, Pedigree


def load_marker_set_from_config(config: Config) -> MarkerSet:
    marker_set = MarkerSet()
    with open(config.marker_set) as file:
        marker_set.read_marker_set_from_file(file)
    return marker_set


def load_marker_set_from_upload(marker_set_file) -> MarkerSet:
    marker_set = MarkerSet()
    marker_set.read_marker_set_from_file(marker_set_file)
    return marker_set


def load_pedigree_from_config(config: Config, marker_set: MarkerSet) -> Pedigree:
    pedigree = Pedigree()
    with open(config.pedigree) as file:
        pedigree.read_pedigree_from_file(file)

    for name, path in config.known_haplotypes.items():
        with open(path) as file:
            pedigree.read_known_haplotype_from_file(name, file, marker_set)

    pedigree.reroot_pedigree(config.suspect)

    return pedigree


def load_pedigree_from_upload(pedigree_file, marker_set: MarkerSet) -> Pedigree:
    pedigree = Pedigree()
    pedigree.read_pedigree_from_file(pedigree_file)
    return pedigree
