"""
Compare Python vs Rust simulation outputs.

Simple pedigree: Father → Suspect + Brother
5 markers, Father has one allelic difference (DYF393S1: 14 vs 15) from Suspect.
Brother's haplotype is unknown.

Usage (from repository root):
    python compare_run.py
"""
import sys
import subprocess
import json
from pathlib import Path
from copy import deepcopy
from unittest.mock import MagicMock

repo_root = Path(__file__).resolve().parent
examples = repo_root / "examples"
sys.path.insert(0, str(repo_root))

# Stub UI/rendering dependencies before any pedigree_lr import (Windows multiprocessing
# re-runs top-level imports in worker processes, so stubs must be at module level).
_STUBS = [
    "weasyprint", "stqdm", "streamlit", "streamlit.delta_generator",
    "streamlit.components", "streamlit.components.v1",
    "streamlit_agraph", "jinja2",
]
for _mod in _STUBS:
    sys.modules[_mod] = MagicMock()

from pedigree_lr.models import SimulationParameters, Pedigree, MarkerSet
import pedigree_lr.visualization as _viz_mod
_viz_mod.save_pedigree_to_png = lambda **kw: None
_viz_mod.plot_probabilities = lambda *a, **kw: None
from pedigree_lr.simulation import run_simulation


class SilentReporter:
    def log(self, msg): pass
    def progress_bar(self, desc):
        class _PB:
            def __enter__(self): return self
            def __exit__(self, *a): pass
            def update(self, *a, **kw): pass
        return _PB()


def load_data():
    marker_set = MarkerSet()
    with open(examples / "simple_markers.csv") as f:
        marker_set.read_marker_set_from_file(f)

    pedigree = Pedigree()
    with open(examples / "simple.tgf") as f:
        pedigree.read_pedigree_from_file(file=f, file_extension=".tgf")
    with open(examples / "simple.json") as f:
        pedigree.read_known_haplotypes_from_file(file=f, marker_set=marker_set)

    pedigree.check_known_haplotypes()
    pedigree.set_suspect("Suspect")
    return marker_set, pedigree


def run_python(marker_set, pedigree):
    results_path = repo_root / "results" / "simple_compare"
    results_path.mkdir(parents=True, exist_ok=True)
    params = SimulationParameters(
        two_step_mutation_factor=0.03,
        stability_window=10000,
        model_validity_threshold=0.05,
        number_of_threads=1,
        simulation_name="simple_compare",
        results_path=results_path,
        bias=None,
    )
    return run_simulation(
        input_pedigree=deepcopy(pedigree),
        suspect_name="Suspect",
        trace=None,
        marker_set=marker_set,
        simulation_parameters=params,
        reporter=SilentReporter(),
    )


def run_rust():
    rust_exe = Path("C:/cargo-target/matchy/release/matchy.exe")
    config_toml = examples / "simple_suspect.toml"
    proc = subprocess.run(
        [str(rust_exe), "--config-path", str(config_toml)],
        capture_output=True, text=True, cwd=str(examples),
    )
    if proc.returncode != 0:
        print("[stderr]", proc.stderr[:3000], file=sys.stderr)
        return None

    # Parse the human-readable stdout summary
    result = {}
    for line in proc.stdout.splitlines():
        line = line.strip()
        if line.startswith("Average pedigree probability:"):
            result["average_pedigree_probability"] = line.split(":", 1)[1].strip()
        elif line.startswith("P(matches = "):
            k = int(line.split("P(matches = ")[1].split(")")[0])
            val = line.split("=", 2)[-1].strip()
            result.setdefault("inside_match_probability", {})[str(k)] = val
        elif line.startswith("P(match outside) ="):
            result["outside_match_probability"] = line.split("=", 1)[1].strip()
    return result


if __name__ == "__main__":
    marker_set, pedigree = load_data()

    print("=" * 60)
    print("Running Python simulation...")
    print("  Pedigree: Father -> Suspect + Brother")
    print("  Suspect known, Father known (1 mutation), Brother unknown")
    print("=" * 60)
    py = run_python(marker_set, pedigree)

    print("\n--- Python Results ---")
    print(f"average_pedigree_probability : {py.average_pedigree_probability}")
    print(f"outside_match_probability    : {py.outside_match_probability}")
    print("inside_match_probability (per k):")
    for k, v in sorted(py.inside_match_probability.items()):
        print(f"  k={k}: {v}")
    print("per_individual_probabilities:")
    for name, v in sorted(py.per_individual_probabilities.items()):
        print(f"  {name}: {v}")

    print("\n" + "=" * 60)
    print("Running Rust simulation...")
    print("=" * 60)
    rust = run_rust()

    if rust:
        print("\n--- Rust Results ---")
        print(f"average_pedigree_probability : {rust.get('average_pedigree_probability', 'n/a')}")
        print(f"outside_match_probability    : {rust.get('outside_match_probability', 'n/a')}")
        inside = rust.get("inside_match_probability", {})
        if inside:
            print("inside_match_probability (per k):")
            for k, v in sorted(inside.items(), key=lambda x: int(x[0])):
                print(f"  k={k}: {v}")
    else:
        print("[Rust simulation failed — check output above]")
