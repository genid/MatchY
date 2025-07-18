import sys
from abc import ABC
from stqdm import stqdm
from streamlit.delta_generator import DeltaGenerator
from tqdm import tqdm
from io import StringIO
from datetime import datetime
from pedigree_lr.models import SimulationParameters, SimulationResult
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from pathlib import Path


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
        title="Match-Y Simulation Report",
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
