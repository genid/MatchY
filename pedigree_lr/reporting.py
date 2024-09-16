import sys
from abc import ABC

from stqdm import stqdm
from streamlit.delta_generator import DeltaGenerator
from tqdm import tqdm


class ProgressBar(ABC):
    def update(self, count: int):
        raise NotImplementedError

    def __enter__(self) -> "ProgressBar":
        raise NotImplementedError

    def __exit__(self, *args) -> None:
        raise NotImplementedError


class ConsoleProgressBar(ProgressBar):
    def __init__(self, total: int, desc: str) -> None:
        self.tqdm = tqdm(total=total, desc=desc, file=sys.__stdout__)

    def update(self, count: int):
        self.tqdm.update(1)

    def __enter__(self) -> ProgressBar:
        return self

    def __exit__(self, *args) -> None:
        self.tqdm.__exit__(*args)


class StreamlitProgressBar(ProgressBar):
    def __init__(self, total: int, desc: str, st_container: DeltaGenerator) -> None:
        self.stqdm = stqdm(total=total, desc=desc, st_container=st_container)

    def update(self, count: int):
        self.stqdm.update(1)

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
