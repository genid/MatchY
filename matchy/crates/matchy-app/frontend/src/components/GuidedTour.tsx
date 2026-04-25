import { useEffect, useState } from "react";

const TOUR_KEY = "matchy_tour_done";

const STEPS = [
  {
    title: "Welcome to MatchY",
    content:
      "MatchY computes Y-STR haplotype match probabilities using Monte Carlo simulation on a pedigree. This short tour walks you through the five-step workflow.",
    icon: "🧬",
  },
  {
    title: "Step 1 — Build or load a pedigree",
    content:
      "Go to the Pedigree tab. Click '+ New pedigree' to start from scratch, or import a TGF / PED file. Draw relationships by dragging the handles on each node. Ctrl+Z undoes the last action. The pedigree must be a directed acyclic graph (DAG).",
    icon: "🌳",
  },
  {
    title: "Step 2 — Select a marker set",
    content:
      "Go to Marker Sets to pick a built-in kit (e.g. RMplex, YFiler) or use the Custom Set Builder to hand-pick markers. You can add custom markers to the pool, override mutation rates (shown in orange when changed), and save sets for reuse. Custom markers and rate changes persist across sessions.",
    icon: "📊",
  },
  {
    title: "Step 3 — Load haplotype data",
    content:
      "Go to the Haplotypes tab. Import a JSON file, or use 'Paste from Excel…' to paste tab-separated data directly from a spreadsheet. Column headers should be individual names; row labels should be marker names. Use the ⚠ icon to spot individuals missing from the pedigree, and the 'exclude' toggle per column to exclude them from the calculation.",
    icon: "🔬",
  },
  {
    title: "Step 4 — Configure and run",
    content:
      "On the Run tab: the data summary bar shows what's loaded. Select a suspect (or enable Trace mode), adjust parameters (or reset to defaults), then press Run Simulation or Ctrl+Enter. Each convergence chart shows the inter-model variance % and dashed lines for the current average and the criterion band.",
    icon: "⚙️",
  },
  {
    title: "Step 5 — Interpret results",
    content:
      "Headline probability cards show pedigree probability, inside-match, and outside-match at a glance. The per-individual match probability table is sorted by probability. Click 'Copy results' to copy all numbers to the clipboard.",
    icon: "📈",
  },
  {
    title: "Step 6 — Report and sessions",
    content:
      "Click Generate Report to produce an interactive HTML report with convergence charts, haplotype tables, and pedigree diagrams. Enable Auto-Save in Settings to automatically store every run. Use Save session / Load session to restore your full workspace, or New session to start fresh.",
    icon: "💾",
  },
];

export function GuidedTour({ onClose }: { onClose: () => void }) {
  const [step, setStep] = useState(0);
  const current = STEPS[step];
  const isLast = step === STEPS.length - 1;

  const finish = () => {
    localStorage.setItem(TOUR_KEY, "1");
    onClose();
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40">
      <div className="bg-white rounded-xl shadow-2xl w-full max-w-md mx-4 overflow-hidden">
        {/* Progress bar */}
        <div className="h-1 bg-gray-100">
          <div
            className="h-1 bg-blue-500 transition-all duration-300"
            style={{ width: `${((step + 1) / STEPS.length) * 100}%` }}
          />
        </div>

        <div className="p-6">
          <div className="text-4xl mb-3">{current.icon}</div>
          <h2 className="text-lg font-bold text-gray-900 mb-2">{current.title}</h2>
          <p className="text-sm text-gray-600 leading-relaxed">{current.content}</p>
        </div>

        <div className="px-6 pb-5 flex items-center justify-between">
          <button
            onClick={finish}
            className="text-xs text-gray-400 hover:text-gray-600 transition-colors"
          >
            Skip tour
          </button>
          <div className="flex items-center gap-2">
            {/* Step dots */}
            <div className="flex gap-1 mr-3">
              {STEPS.map((_, i) => (
                <button
                  key={i}
                  onClick={() => setStep(i)}
                  className={`w-1.5 h-1.5 rounded-full transition-colors ${
                    i === step ? "bg-blue-500" : "bg-gray-300"
                  }`}
                />
              ))}
            </div>
            {step > 0 && (
              <button
                onClick={() => setStep((s) => s - 1)}
                className="text-sm border border-gray-300 text-gray-600 hover:bg-gray-50 px-3 py-1.5 rounded transition-colors"
              >
                Back
              </button>
            )}
            <button
              onClick={() => (isLast ? finish() : setStep((s) => s + 1))}
              className="text-sm bg-blue-600 hover:bg-blue-700 text-white font-medium px-4 py-1.5 rounded transition-colors"
            >
              {isLast ? "Get started" : "Next"}
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}

export function useTour() {
  const [show, setShow] = useState(false);

  useEffect(() => {
    if (!localStorage.getItem(TOUR_KEY)) {
      setShow(true);
    }
  }, []);

  return {
    show,
    open: () => setShow(true),
    close: () => setShow(false),
  };
}
