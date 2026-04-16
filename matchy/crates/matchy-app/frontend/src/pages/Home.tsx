import { useState } from "react";
import { useAppStore } from "../store/appStore";
import { useSimulation } from "../hooks/useSimulation";
import { useSettings } from "./Settings";

export default function Home() {
  const { pedigree, haplotypes, markers, suspect, setSuspect, simulation } =
    useAppStore();
  const { startSimulation, cancelSimulation } = useSimulation();
  const savedSettings = useSettings();
  const [params, setParams] = useState({
    twoStepMutationFraction: savedSettings.defaultTwoStepFraction,
    batchLength: savedSettings.defaultBatchLength,
    convergenceCriterion: savedSettings.defaultConvergenceCriterion,
    bias: null as number | null,
    numberOfThreads: savedSettings.defaultThreads,
    skipInside: false,
    skipOutside: false,
    traceMode: false,
    adaptiveBias: false,
    simulationName: savedSettings.simulationName,
    userName: savedSettings.userName,
  });

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
            <section className="bg-white rounded-lg border p-4 space-y-2">
              <h2 className="font-semibold text-gray-700">Results</h2>
              <p className="text-sm text-gray-600">
                Converged: <strong>{simulation.result.converged ? "Yes" : "No"}</strong>
                {" · "}Trials: <strong>{simulation.result.trials}</strong>
              </p>
              {simulation.result.inside_match_probabilities && (
                <div className="text-sm">
                  <p className="font-medium text-gray-700">Inside-pedigree match probabilities</p>
                  <p className="text-gray-500 text-xs">
                    Avg. pedigree probability: {simulation.result.inside_match_probabilities.average_pedigree_probability}
                  </p>
                  <table className="text-xs mt-1 w-full border-collapse">
                    <thead>
                      <tr className="bg-gray-50">
                        <th className="border px-2 py-1 text-left">Matches</th>
                        <th className="border px-2 py-1 text-left">Probability</th>
                      </tr>
                    </thead>
                    <tbody>
                      {Object.entries(simulation.result.inside_match_probabilities.probabilities)
                        .sort(([a], [b]) => parseInt(a) - parseInt(b))
                        .map(([n, p]) => (
                          <tr key={n}>
                            <td className="border px-2 py-1">{n}</td>
                            <td className="border px-2 py-1">{p}</td>
                          </tr>
                        ))}
                    </tbody>
                  </table>
                </div>
              )}
              {simulation.result.outside_match_probability && (
                <p className="text-sm">
                  Outside match probability:{" "}
                  <strong>{simulation.result.outside_match_probability}</strong>
                </p>
              )}
              {!simulation.result.inside_match_probabilities && !simulation.result.outside_match_probability && (
                <p className="text-sm text-gray-500">
                  No suspect specified — pedigree probability only.
                </p>
              )}
            </section>
          )}
        </div>
      </div>
    </div>
  );
}
