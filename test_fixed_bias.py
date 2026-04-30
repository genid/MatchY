"""
Test fixed bias vs auto bias in the Python simulation.
Runs same pedigree with bias=None (auto), bias=0.2 (fixed), bias=0.05 (fixed small).
Results should be numerically close since all are valid IS proposals.
"""
import sys
from pathlib import Path
from copy import deepcopy
from unittest.mock import MagicMock
from decimal import Decimal

repo_root = Path(__file__).resolve().parent
examples = repo_root / "examples"
sys.path.insert(0, str(repo_root))

_STUBS = ["weasyprint", "stqdm", "streamlit", "streamlit.delta_generator",
          "streamlit.components", "streamlit.components.v1",
          "streamlit_agraph", "jinja2"]
for _mod in _STUBS:
    sys.modules[_mod] = MagicMock()

from pedigree_lr.models import SimulationParameters, Pedigree, MarkerSet
import pedigree_lr.visualization as _viz_mod
_viz_mod.save_pedigree_to_png = lambda **kw: None
_viz_mod.plot_probabilities = lambda *a, **kw: None
from pedigree_lr.simulation import run_simulation


class SilentReporter:
    def log(self, msg): print(msg, flush=True)
    def progress_bar(self, desc):
        class _PB:
            def __enter__(self): return self
            def __exit__(self, *a): pass
            def update(self, *a, **kw): pass
        return _PB()


def load_simple():
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
    return marker_set, pedigree, "Suspect"


def load_mockfo():
    from pedigree_lr.data import load_marker_set_from_database
    marker_set = load_marker_set_from_database("RMplex")
    pedigree = Pedigree()
    with open(examples / "mockfo.tgf") as f:
        pedigree.read_pedigree_from_file(file=f, file_extension=".tgf")
    with open(examples / "mockfo.json") as f:
        pedigree.read_known_haplotypes_from_file(file=f, marker_set=marker_set)
    pedigree.check_known_haplotypes()
    pedigree.set_suspect("Verdachte")
    return marker_set, pedigree, "Verdachte"


def run_with_bias(marker_set, pedigree, suspect_name, bias_value, label):
    results_path = repo_root / "results" / f"bias_test_{label}"
    results_path.mkdir(parents=True, exist_ok=True)
    params = SimulationParameters(
        two_step_mutation_factor=0.03,
        stability_window=10000,
        model_validity_threshold=0.05,
        number_of_threads=3,
        simulation_name=f"bias_test_{label}",
        results_path=results_path,
        bias=bias_value,
    )
    result = run_simulation(
        input_pedigree=deepcopy(pedigree),
        suspect_name=suspect_name,
        trace=None,
        marker_set=marker_set,
        simulation_parameters=params,
        reporter=SilentReporter(),
    )
    return result


if __name__ == "__main__":
    import sys
    use_mockfo = "--mockfo" in sys.argv
    print("Loading data...", flush=True)
    if use_mockfo:
        marker_set, pedigree, suspect = load_mockfo()
        print("Using mockfo pedigree", flush=True)
    else:
        marker_set, pedigree, suspect = load_simple()
        print("Using simple pedigree", flush=True)

    results = {}
    for bias_val, label in [(None, "auto"), (0.2, "fixed_0.2"), (0.05, "fixed_0.05"), (0.4, "fixed_0.4")]:
        print(f"\n{'='*60}", flush=True)
        print(f"Running with bias={bias_val} ({label})", flush=True)
        r = run_with_bias(marker_set, pedigree, bias_val, label)
        r._suspect = suspect
        results[label] = r
        print(f"  avg_ped_p   = {r.average_pedigree_probability:.6E}")
        print(f"  outside_p   = {r.outside_match_probability:.6E}")
        print(f"  inside (k=1)= {r.inside_match_probability.get(1, 'N/A')}")
        print(f"  per_indiv   = {dict(r.per_individual_probabilities)}")

    print("\n" + "="*60)
    print("COMPARISON:")
    auto = results["auto"]
    for label, r in results.items():
        if label == "auto":
            continue
        ratio_ped = float(r.average_pedigree_probability) / float(auto.average_pedigree_probability) if auto.average_pedigree_probability else float('inf')
        ratio_out = float(r.outside_match_probability) / float(auto.outside_match_probability) if auto.outside_match_probability else float('inf')
        k1_auto = auto.inside_match_probability.get(1, Decimal(0))
        k1_this = r.inside_match_probability.get(1, Decimal(0))
        ratio_in = float(k1_this) / float(k1_auto) if k1_auto else float('inf')
        print(f"  {label} vs auto: avg_ped_p ratio={ratio_ped:.3f}, outside ratio={ratio_out:.3f}, inside k=1 ratio={ratio_in:.3f}")
