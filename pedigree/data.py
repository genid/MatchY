from pedigree.config import Config
from pedigree.models import MarkerSet, Pedigree


def load_marker_set(config: Config) -> MarkerSet:
    marker_set = MarkerSet()
    marker_set.read_marker_set_from_file(config.marker_set)
    return marker_set


def load_pedigree(config: Config, marker_set: MarkerSet) -> Pedigree:
    pedigree = Pedigree()
    pedigree.read_pedigree_from_file(config.pedigree)

    for name, path in config.known_haplotypes.items():
        pedigree.read_known_haplotype_from_file(name, path, marker_set)

    pedigree.reroot_pedigree(config.suspect)

    return pedigree
