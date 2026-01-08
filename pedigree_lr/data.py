import json
from io import StringIO
from pedigree_lr.config import Config
from pedigree_lr.models import MarkerSet, Pedigree, Haplotype
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
    with open(config.pedigree, "r") as file: # type: StringIO
        pedigree.read_pedigree_from_file(file=file,
                                         file_extension=file_extension)

    with open(config.known_haplotypes) as file: # type: StringIO
        trace_from_json = pedigree.read_known_haplotypes_from_file(file=file,
                                                                    marker_set=marker_set)
        # Store trace temporarily on pedigree object for CLI retrieval
        if trace_from_json is not None:
            pedigree._trace_from_json = trace_from_json

    pedigree.check_known_haplotypes()
    if not pedigree.check_pedigree_structure():
        raise ValueError("Pedigree structure is invalid.")

    pedigree.set_suspect(config.suspect)

    return pedigree


def load_pedigree_from_upload(pedigree_file: StringIO,
                              file_extension: str) -> Pedigree | None:
    pedigree = Pedigree()
    pedigree.read_pedigree_from_file(pedigree_file, file_extension)
    if not pedigree.check_pedigree_structure():
        return None
    else:
        return pedigree


def load_trace_from_file(trace_path: Path,
                         marker_set: MarkerSet) -> Haplotype:
    trace_profile = Haplotype()
    with open(trace_path, "r") as file:
        for line in file.readlines():
            marker_name, alleles_values = line.strip().split(",")
            marker = marker_set.get_marker_by_name(marker_name)
            if not marker:
                raise ValueError(f"Marker {marker_name} not found in marker set.")

            alleles = alleles_values.split(";")
            number_of_copies = len(alleles)
            if not marker.number_of_copies:
                marker.number_of_copies = number_of_copies
            elif marker.number_of_copies != number_of_copies:
                raise ValueError(f"Marker {marker_name} has inconsistent number of copies.")
            for allele in alleles:
                if "." in allele:
                    allele_val, intermediate_value = allele.split(".")
                    try:
                        allele_int = int(allele_val)
                        intermediate_int = int(intermediate_value)
                        trace_profile.add_allele(marker, allele_int, intermediate_int)
                    except ValueError:
                        raise ValueError(f"Invalid allele value: {allele}")
                else:
                    try:
                        allele_int = int(allele)
                        trace_profile.add_allele(marker, allele_int)
                    except ValueError:
                        raise ValueError(f"Invalid allele value: {allele}")
    return trace_profile