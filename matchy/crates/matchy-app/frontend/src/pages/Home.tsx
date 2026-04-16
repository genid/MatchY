import { useState } from "react";
import { useAppStore } from "../store/appStore";
import { useSimulation } from "../hooks/useSimulation";
import type { SimulationRequest } from "../types/matchy";

// Default simulation parameters
const DEFAULTS = {
  twoStepMutationFraction: 0.03,
  batchLength: 10000,
  convergenceCriterion: 0.02,
  bias: null as number | null,
  numberOfThreads: 4,
  skipInside: false,
  skipOutside: false,
  traceMode: false,
  adaptiveBias: false,
  simulationName: "simulation",
  userName: "",
};

export default function Home() {
  const { pedigree, haplotypes, markers, suspect, setSuspect, exclude, setExclude, simulation } =
    useAppStore();
  const { startSimulation, cancelSimulation } = useSimulation();
  const [params, setParams] = useState(DEFAULTS);

  const unknownIndividuals =
    pedigree?.individuals.filter((i) => i.haplotypeClass === "unknown" && !i.exclude) ?? [];

  const knownIndividuals =
    pedigree?.individuals.filter((i) => !i.exclude) ?? [];

  const canRun = !!pedigree && haplotypes !== null && markers.length > 0;

  const handleRun = () => {
    if (!canRun) return;
    startSimulation({
      ...params,
    });
  };

  return (
    <div className="p-6 max-w-5xl mx-auto space-y-6">
      <h1 className="text-2xl font-bold text-gray-900">Run Simulation</h1>

      {/* Status banners */}
      {!pedigree && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          No pedigree loaded. Go to <strong>Pedigree</strong> to load one.
        </div>
      )}
      {pedigree && haplotypes === null && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          No haplotypes loaded. Go to <strong>Haplotypes</strong> to load them.
        </div>
      )}
      {pedigree && haplotypes !== null && markers.length === 0 && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          No marker set loaded. Go to <strong>Marker Sets</strong> to select one.
        </div>
      )}

      <div className="grid grid-cols-2 gap-6">
        {/* Left: configuration */}
        <div className="space-y-4">
          <section className="bg-white rounded-lg border p-4 space-y-3">
            <h2 className="font-semibold text-gray-700">Mode</h2>
            <label className="flex items-center gap-2 text-sm">
              <input
                type="checkbox"
                checked={params.traceMode}
                onChange={(e) => setParams((p) => ({ ...p, traceMode: e.target.checked }))}
              />
              Trace mode (identify most likely donor)
            </label>

            {!params.traceMode && (
              <div>
                <label className="block text-sm font-medium text-gray-600 mb-1">Suspect</label>
                <select
                  className="w-full border rounded px-2 py-1.5 text-sm"
                  value={suspect ?? ""}
                  onChange={(e) => setSuspect(e.target.value || null)}
                >
                  <option value="">— select —</option>
                  {knownIndividuals.map((i) => (
                    <option key={i.id} value={i.name}>
                      {i.name}
                    </option>
                  ))}
                </select>
              </div>
            )}
          </section>

          <section className="bg-white rounded-lg border p-4 space-y-3">
            <h2 className="font-semibold text-gray-700">Parameters</h2>
            <div className="grid grid-cols-2 gap-3 text-sm">
              <div>
                <label className="block text-gray-600 mb-1">Two-step fraction</label>
                <input
                  type="number"
                  step="0.001"
                  min="0"
                  max="1"
                  className="w-full border rounded px-2 py-1"
                  value={params.twoStepMutationFraction}
                  onChange={(e) =>
                    setParams((p) => ({ ...p, twoStepMutationFraction: parseFloat(e.target.value) }))
                  }
                />
              </div>
              <div>
                <label className="block text-gray-600 mb-1">Batch length</label>
                <input
                  type="number"
                  min="100"
                  className="w-full border rounded px-2 py-1"
                  value={params.batchLength}
                  onChange={(e) =>
                    setParams((p) => ({ ...p, batchLength: parseInt(e.target.value, 10) }))
                  }
                />
              </div>
              <div>
                <label className="block text-gray-600 mb-1">Convergence criterion</label>
                <input
                  type="number"
                  step="0.001"
                  min="0"
                  className="w-full border rounded px-2 py-1"
                  value={params.convergenceCriterion}
                  onChange={(e) =>
                    setParams((p) => ({ ...p, convergenceCriterion: parseFloat(e.target.value) }))
                  }
                />
              </div>
              <div>
                <label className="block text-gray-600 mb-1">Threads</label>
                <input
                  type="number"
                  min="1"
                  max="32"
                  className="w-full border rounded px-2 py-1"
                  value={params.numberOfThreads}
                  onChange={(e) =>
                    setParams((p) => ({ ...p, numberOfThreads: parseInt(e.target.value, 10) }))
                  }
                />
              </div>
            </div>
            <div className="flex gap-4 text-sm">
              <label className="flex items-center gap-2">
                <input
                  type="checkbox"
                  checked={params.skipInside}
                  onChange={(e) => setParams((p) => ({ ...p, skipInside: e.target.checked }))}
                />
                Skip inside
              </label>
              <label className="flex items-center gap-2">
                <input
                  type="checkbox"
                  checked={params.skipOutside}
                  onChange={(e) => setParams((p) => ({ ...p, skipOutside: e.target.checked }))}
                />
                Skip outside
              </label>
              <label className="flex items-center gap-2">
                <input
                  type="checkbox"
                  checked={params.adaptiveBias}
                  onChange={(e) => setParams((p) => ({ ...p, adaptiveBias: e.target.checked }))}
                />
                Adaptive bias
              </label>
            </div>
          </section>

          <div className="flex gap-3">
            <button
              disabled={!canRun || simulation.running}
              onClick={handleRun}
              className="flex-1 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-300 text-white font-medium py-2 px-4 rounded transition-colors"
            >
              {simulation.running ? "Running…" : "Run Simulation"}
            </button>
            {simulation.running && (
              <button
                onClick={cancelSimulation}
                className="bg-red-100 hover:bg-red-200 text-red-700 font-medium py-2 px-4 rounded"
              >
                Cancel
              </button>
            )}
          </div>
        </div>

        {/* Right: progress / results */}
        <div className="space-y-4">
          {simulation.error && (
            <div className="rounded bg-red-50 border border-red-200 p-3 text-sm text-red-800">
              {simulation.error}
            </div>
          )}
          {simulation.progress.length > 0 && (
            <section className="bg-white rounded-lg border p-4">
              <h2 className="font-semibold text-gray-700 mb-3">Progress</h2>
              <div className="text-xs font-mono space-y-1 max-h-64 overflow-y-auto">
                {simulation.progress.slice(-20).map((e, i) => (
                  <div key={i} className="text-gray-600">
                    Trial {e.trial} · Model {e.model} · Iter {e.iteration} · Mean {e.currentMean}
                    {e.converged && " ✓"}
                  </div>
                ))}
              </div>
            </section>
          )}
          {simulation.result && (
            <section className="bg-white rounded-lg border p-4">
              <h2 className="font-semibold text-gray-700 mb-2">Results</h2>
              <pre className="text-xs text-gray-600 overflow-auto max-h-64">
                {JSON.stringify(simulation.result, null, 2)}
              </pre>
            </section>
          )}
        </div>
      </div>
    </div>
  );
}
