import sys
from abc import ABC

from stqdm import stqdm
from streamlit.delta_generator import DeltaGenerator
from tqdm import tqdm
from io import StringIO
from datetime import datetime


def create_report_bytes(
        self,
        random_seed: int,
) -> bytes:
    bytes_data = StringIO()

    bytes_data.write("Simulation results\n")
    bytes_data.write("match-Y version 1.0.0\n\n")  # TODO: Remove hard coded version

    bytes_data.write(f"Date and time of report: \t{datetime.now()}\n")
    bytes_data.write(f"Number of iterations: \t{self.simulation_parameters['number_of_iterations']}\n")
    bytes_data.write(f"Random seed: \t{random_seed}\n\n")

    bytes_data.write("Marker set with mutation rate\n")
    for marker in self.marker_set.markers:
        bytes_data.write(f"{marker.name}: {marker.mutation_rate}\n")

    bytes_data.write(f"\nNumber of nodes in pedigree: \t{len(self.pedigree.individuals)}\n")
    bytes_data.write(f"Number of edges in pedigree: \t{len(self.pedigree.relationships)}\n")

    bytes_data.write("\nIndividuals:\n")
    bytes_data.write("ID\tName\tClass\tExcluded\n")
    for individual in self.pedigree.individuals:
        bytes_data.write(f"{individual.id}\t{individual.name}\t{individual.haplotype_class}\t{individual.exclude}\n")

    bytes_data.write("\nKnown haplotypes:\n")
    for individual in self.pedigree.individuals:
        if individual.haplotype_class != "unknown":
            bytes_data.write(f"{individual.name}:\n")
            for marker, alleles in individual.haplotype.alleles.items():
                bytes_data.write(f"{marker}: {','.join(str(allele) for allele in alleles)}\n")

    bytes_data.write("\nRelationships\n")
    bytes_data.write("Parent ID -> Child ID\n")
    for relationship in self.pedigree.relationships:
        bytes_data.write(f"{relationship.parent_id} -> {relationship.child_id}\n")

    bytes_data.write(f"\nSuspect: \t{self.suspect_name}\n\n")

    bytes_data.write(f"Average pedigree probability: \t{self.average_pedigree_probability}\n\n")
    bytes_data.write(f"Run time average pedigree probability: \t{self.run_time_pedigree_probability}\n")
    bytes_data.write(f"Run time proposal distribution: \t{self.run_time_proposal_distribution}\n")
    bytes_data.write(f"Total run time: \t{self.total_run_time}\n\n")

    bytes_data.write("\nMatch probabilities\n")
    bytes_data.write("i\tProbability\tNeeded iterations\tModel probabilities\n")

    for key in sorted(self.proposal_distribution.keys()):
        if key in self.l_needed_iterations:
            bytes_data.write(
                f"{key}\t{self.proposal_distribution[key]:.4E}\t{self.l_needed_iterations[key]}\t{self.l_model_probabilities[key]}\n")
        else:
            bytes_data.write(
                f"{key}\t{self.proposal_distribution[key]:.4E}\n")

    bytes_data.write(f"\nOutside match probability: \t{self.outside_match_probability:.4E}\n")
    bytes_data.write(f"Outside match probability needed iterations: {self.outside_needed_iterations}\n")
    bytes_data.write(f"Outside match probability model probabilities: {self.outside_model_probabilities}\n")

    return bytes_data.getvalue().encode("utf-8")


class ProgressBar(ABC):
    def update(self, count: int):
        raise NotImplementedError

    def update_total(self, total: int):
        raise NotImplementedError

    def __enter__(self) -> "ProgressBar":
        raise NotImplementedError

    def __exit__(self, *args) -> None:
        raise NotImplementedError


class ConsoleProgressBar(ProgressBar):
    def __init__(self, total: int, desc: str) -> None:
        self.tqdm = tqdm(total=total, desc=desc, file=sys.__stdout__)

    def update(self, count: int):
        self.tqdm.update(count)

    def update_total(self, total: int):
        self.tqdm.total -= total
        self.tqdm.refresh()

    def __enter__(self) -> ProgressBar:
        return self

    def __exit__(self, *args) -> None:
        self.tqdm.__exit__(*args)


class StreamlitProgressBar(ProgressBar):
    def __init__(self, total: int, desc: str, st_container: DeltaGenerator) -> None:
        self.stqdm = stqdm(total=total, desc=desc, st_container=st_container)

    def update(self, count: int):
        self.stqdm.update(count)

    def update_total(self, total: int):
        self.stqdm.total -= total
        self.stqdm.refresh()

    def __enter__(self) -> ProgressBar:
        return self

    def __exit__(self, *args) -> None:
        self.stqdm.__exit__(*args)


class Reporter(ABC):
    def log(self, msg: str) -> None:
        raise NotImplementedError

    def progress_bar(self, total: int, desc: str) -> ProgressBar:
        raise NotImplementedError


class ConsoleReporter(Reporter):
    def log(self, msg: str) -> None:
        print(msg)

    def progress_bar(self, total: int, desc: str) -> ProgressBar:
        return ConsoleProgressBar(total, desc)


class StreamlitReporter(Reporter):
    def __init__(
        self, log_container: DeltaGenerator, progress_container: DeltaGenerator
    ) -> None:
        self.log_container = log_container
        self.progress_container = progress_container

    def log(self, msg: str) -> None:
        self.log_container.text(msg)

    def progress_bar(self, total: int, desc: str) -> ProgressBar:
        return StreamlitProgressBar(total, desc, self.progress_container)
