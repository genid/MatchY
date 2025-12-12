import sys
from abc import ABC
from stqdm import stqdm
from streamlit.delta_generator import DeltaGenerator
from tqdm import tqdm
from datetime import datetime
from pedigree_lr.models import SimulationResult
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from pathlib import Path
import logging


def setup_logger_cli(name=None, level=logging.INFO):
    """
    Setup a logger for CLI use.
    Logs go to stdout.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        cli_handler = logging.StreamHandler(sys.stdout)
        cli_handler.setLevel(level)
        formatter = logging.Formatter("[%(levelname)s] %(name)s: %(message)s")
        cli_handler.setFormatter(formatter)
        logger.addHandler(cli_handler)

    return logger


def setup_logger_streamlit(name=None, level=logging.INFO):
    """
    Setup a logger for Streamlit pages.
    Logs go to stdout AND Streamlit colored boxes.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # CLI handler
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        cli_handler = logging.StreamHandler(sys.stdout)
        cli_handler.setLevel(level)
        formatter = logging.Formatter("[%(levelname)s]: %(message)s")
        cli_handler.setFormatter(formatter)
        logger.addHandler(cli_handler)

    # Streamlit handler
    try:
        import streamlit as st

        if not any(getattr(h, "_is_streamlit_handler", False) for h in logger.handlers):

            class StreamlitHandler(logging.Handler):
                _is_streamlit_handler = True  # mark handler to prevent duplicates

                def emit(self, record):
                    msg = self.format(record)
                    if record.levelno >= logging.ERROR:
                        st.error(msg)
                    elif record.levelno >= logging.WARNING:
                        st.warning(msg)
                    else:
                        st.info(msg)

            sh = StreamlitHandler()
            sh.setLevel(level)
            sh.setFormatter(logging.Formatter("[%(levelname)s]: %(message)s"))
            logger.addHandler(sh)

    except ModuleNotFoundError:
        # Not running in Streamlit
        pass

    return logger


def create_html_pdf_report(
        result: SimulationResult
) -> bytes:
    report_template_folder = Path(__file__).resolve().parent.parent / "data"
    env = Environment(loader=FileSystemLoader(searchpath=report_template_folder))  # assumes template in cwd
    template = env.get_template("report_template.html")

    results_path = Path(result.simulation_parameters.results_path)
    images = list(results_path.glob("*.png"))
    images = [str(image) for image in images if image.suffix == ".png"]

    html_out = template.render(
        title="Simulation Report",
        subtitle="Pedigree based match probability results",
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        simulation_parameters=result.simulation_parameters,
        result=result,
        images=images,
        logo_path=str(Path(__file__).resolve().parent.parent / "logo.png"),
    )

    # Generate PDF
    html = HTML(string=html_out, base_url=".")
    pdf_bytes = html.write_pdf()

    return pdf_bytes


def normalize_probabilities(per_individual_probabilities: dict, outside_match_probability=None) -> list:
    """
    Normalize per-individual match probabilities to sum to 100%, including outside match probability.

    Args:
        per_individual_probabilities: Dictionary mapping individual ID to match probability (as Decimal)
        outside_match_probability: Outside pedigree match probability (as Decimal), if available

    Returns:
        List of tuples (individual_id, normalized_percentage) sorted by probability (highest first)
        Includes an "Outside Pedigree" entry if outside_match_probability is provided
    """
    if not per_individual_probabilities:
        return []

    # Calculate total probability including outside match probability
    total_prob = sum(float(prob) for prob in per_individual_probabilities.values())

    if outside_match_probability is not None:
        total_prob += float(outside_match_probability)

    # Avoid division by zero
    if total_prob == 0:
        result = [(ind_id, 0.0) for ind_id in per_individual_probabilities.keys()]
        if outside_match_probability is not None:
            result.append(("Outside Pedigree", 0.0))
        return result

    # Normalize each probability to percentage (sum = 100%)
    normalized = [
        (ind_id, (float(prob) / total_prob) * 100.0)
        for ind_id, prob in per_individual_probabilities.items()
    ]

    # Add outside match probability if provided
    if outside_match_probability is not None:
        outside_normalized = (float(outside_match_probability) / total_prob) * 100.0
        normalized.append(("Outside Pedigree", outside_normalized))

    # Sort by probability (highest first)
    normalized.sort(key=lambda x: x[1], reverse=True)

    return normalized


def create_trace_mode_report(
        result: SimulationResult
) -> bytes:
    """
    Generate a simplified PDF report for trace mode analysis.

    Focuses on identifying the most likely trace donor by showing normalized
    per-individual match probabilities. Omits sections not relevant to trace analysis.

    Args:
        result: SimulationResult object containing simulation results

    Returns:
        PDF report as bytes
    """
    report_template_folder = Path(__file__).resolve().parent.parent / "data"
    env = Environment(loader=FileSystemLoader(searchpath=report_template_folder))
    template = env.get_template("trace_report_template.html")

    results_path = Path(result.simulation_parameters.results_path)
    images = list(results_path.glob("*.png"))
    images = [str(image) for image in images if image.suffix == ".png"]

    # Normalize probabilities to sum to 100% (including outside match probability)
    normalized_probs = normalize_probabilities(
        result.per_individual_probabilities,
        result.outside_match_probability
    )

    # Extract summary statistics
    most_likely_donor = normalized_probs[0][0] if normalized_probs else "N/A"
    highest_probability = normalized_probs[0][1] if normalized_probs else 0.0

    html_out = template.render(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        simulation_parameters=result.simulation_parameters,
        result=result,
        normalized_probabilities=normalized_probs,
        most_likely_donor=most_likely_donor,
        highest_probability=highest_probability,
        total_run_time=result.total_run_time,
        images=images,
        logo_path=str(Path(__file__).resolve().parent.parent / "logo.png"),
    )

    # Generate PDF
    html = HTML(string=html_out, base_url=".")
    pdf_bytes = html.write_pdf()

    return pdf_bytes


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
    def __init__(self, desc: str) -> None:
        # self.tqdm = tqdm(total=total, desc=desc, file=sys.__stdout__)
        self.tqdm = tqdm(desc=desc, file=sys.__stdout__)

    def update(self, count: int):
        self.tqdm.update(count)

    def update_total(self, total: int):
        # self.tqdm.total -= total
        self.tqdm.refresh()

    def __enter__(self) -> ProgressBar:
        return self

    def __exit__(self, *args) -> None:
        self.tqdm.__exit__(*args)


class StreamlitProgressBar(ProgressBar):
    def __init__(self, desc: str, st_container: DeltaGenerator) -> None:
        self.stqdm = stqdm(desc=desc, st_container=st_container)

    def update(self, count: int):
        self.stqdm.update(count)

    def update_total(self, total: int):
        self.stqdm.refresh()

    def __enter__(self) -> ProgressBar:
        return self

    def __exit__(self, *args) -> None:
        self.stqdm.__exit__(*args)


class Reporter(ABC):
    def log(self, msg: str) -> None:
        raise NotImplementedError

    def progress_bar(self, desc: str) -> ProgressBar:
        raise NotImplementedError


class ConsoleReporter(Reporter):
    def log(self, msg: str) -> None:
        print(msg)

    def progress_bar(self, desc: str) -> ProgressBar:
        return ConsoleProgressBar(desc)


class StreamlitReporter(Reporter):
    def __init__(
        self, log_container: DeltaGenerator, progress_container: DeltaGenerator
    ) -> None:
        self.log_container = log_container
        self.progress_container = progress_container

    def log(self, msg: str) -> None:
        self.log_container.text(msg)

    def progress_bar(self, desc: str) -> ProgressBar:
        return StreamlitProgressBar(desc=desc, st_container=self.progress_container)
