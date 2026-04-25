import { useRef, useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import { useSimulation } from "../hooks/useSimulation";
import { useSettings } from "./Settings";
import {
  ConvergenceChart,
  type ConvergenceChartRef,
} from "../components/ConvergenceChart";
import { renderPedigreeSvgDataUrl } from "../utils/pedigreeSvg";
import { saveSession, loadSession } from "../utils/session";
import type { PedigreeData } from "../types/matchy";

// Format a number with at least 2 decimal places; 2 sig figs for small numbers;
// scientific notation (2 dp) outside [1e-3, 1e3).
function fmt2sig(n: number): string {
  if (n === 0) return "0.00";
  const abs = Math.abs(n);
  if (abs >= 1e-3 && abs < 1e3) {
    const mag = Math.floor(Math.log10(abs));
    const decimals = Math.max(2, 1 - mag); // at least 2 decimal places
    return n.toFixed(decimals);
  }
  return n.toExponential(2);
}

function fmtPct(v: number | string): string {
  const n = typeof v === "string" ? parseFloat(v) : v;
  if (!isFinite(n)) return "—";
  return fmt2sig(n * 100) + " %";
}

function fmtLr(lr: number | null): string {
  if (lr === null) return "∞";
  return fmt2sig(lr);
}


export default function Home() {
  const { pedigree, pedigreeTgf, haplotypes, haplotypesJson, markers, markerSetCsv, selectedKitName, suspect, setSuspect, exclude, setExclude, simulationName, setSimulationName, userName, setUserName, simulation, setPedigree, setHaplotypes, setMarkerSet, setSimulationResult, setSimulationProgress, resetAll } =
    useAppStore();
  const { startSimulation, cancelSimulation } = useSimulation();
  const savedSettings = useSettings();
  const [params, setParams] = useState({
    twoStepMutationFraction: savedSettings.defaultTwoStepFraction,
    batchLength: savedSettings.defaultBatchLength,
    convergenceCriterion: savedSettings.defaultConvergenceCriterion,
    bias: null as number | null,
    seed: null as number | null,
    numberOfThreads: savedSettings.defaultThreads,
    skipInside: false,
    skipOutside: false,
    traceMode: false,
    adaptiveBias: false,
  });
  const [cpuCount, setCpuCount] = useState<number>(256);
  useEffect(() => {
    invoke<number>("get_cpu_count").then(setCpuCount).catch(() => {});
  }, []);

  const [reportError, setReportError] = useState<string | null>(null);
  const [reportGenerating, setReportGenerating] = useState(false);
  const [sessionError, setSessionError] = useState<string | null>(null);
  const [autoSaveMsg, setAutoSaveMsg] = useState<string | null>(null);
  const [copyMsg, setCopyMsg] = useState(false);
  const [sessionBusy, setSessionBusy] = useState(false);

  // Persist sim name and analyst name across tab navigation
const pedigreeChartRef = useRef<ConvergenceChartRef>(null);
  const extendedPedigreeChartRef = useRef<ConvergenceChartRef>(null);
  const insideChartRef = useRef<ConvergenceChartRef>(null);
  const outsideChartRef = useRef<ConvergenceChartRef>(null);

  const navigate = useNavigate();

  // Ctrl+Enter to run
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.ctrlKey || e.metaKey) && e.key === "Enter" && canRun && !simulation.running) {
        e.preventDefault();
        handleRun();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  });

  const allIndividuals = pedigree?.individuals ?? [];
  const haplotypeNames = new Set(Object.keys(haplotypes?.haplotypeTable ?? {}));
  const knownIndividuals = allIndividuals.filter(
    (i) => !i.exclude && haplotypeNames.has(i.name)
  );

  const canRun = !!pedigree && haplotypes !== null && markers.length > 0;

  const handleRun = () => {
    if (!canRun) return;
    startSimulation({ ...params, simulationName, userName });
  };

  const handleSaveSession = async () => {
    if (!pedigreeTgf || !haplotypesJson) return;
    setSessionError(null);
    setSessionBusy(true);
    try {
      await saveSession(
        {
          version: 1,
          pedigreeTgf,
          haplotypesJson,
          selectedKitName,
          markerSetCsv,
          suspect,
          exclude,
          params: { ...params, simulationName, userName },
          simulationResult: simulation.result,
          simulationProgress: simulation.progress,
        },
        simulationName || "session",
      );
    } catch (e) {
      setSessionError(String(e));
    } finally {
      setSessionBusy(false);
    }
  };

  const handleLoadSession = async () => {
    setSessionError(null);
    setSessionBusy(true);
    try {
      const sess = await loadSession();
      if (!sess) return;

      setPedigree(sess.pedigreeData as PedigreeData, sess.pedigreeTgf);
      setHaplotypes(sess.haplotypesData, sess.haplotypesJson);
      setMarkerSet(sess.selectedKitName, sess.markers, sess.markerSetCsv);
      setSuspect(sess.suspect);
      setExclude(sess.exclude);
      const { simulationName: sn, userName: un, ...restParams } = sess.params;
      setParams(restParams);
      setSimulationName(sn ?? "");
      setUserName(un ?? "");
      if (sess.simulationProgress.length > 0) {
        setSimulationProgress(sess.simulationProgress);
      }
      if (sess.simulationResult) {
        setSimulationResult({ success: true, error: null, result: sess.simulationResult });
      }
    } catch (e) {
      setSessionError(String(e));
    } finally {
      setSessionBusy(false);
    }
  };

  const handleGenerateReport = async () => {
    if (!simulation.result) return;
    setReportError(null);
    setReportGenerating(true);
    try {
      // Collect chart images
      const chartImages: Record<string, string> = {};
      const pedigreeImg = pedigreeChartRef.current?.toBase64Image();
      const extendedPedigreeImg = extendedPedigreeChartRef.current?.toBase64Image();
      const insideImg = insideChartRef.current?.toBase64Image();
      const outsideImg = outsideChartRef.current?.toBase64Image();
      if (pedigreeImg) chartImages["pedigree_probability"] = pedigreeImg;
      if (extendedPedigreeImg) chartImages["extended_pedigree_probability"] = extendedPedigreeImg;
      if (insideImg) chartImages["inside_match_probability"] = insideImg;
      if (outsideImg) chartImages["outside_match_probability"] = outsideImg;

      const knownNames = new Set<string>(
        haplotypes ? Object.keys(haplotypes.haplotypeTable) : []
      );
      const pedigreeImage = pedigree
        ? renderPedigreeSvgDataUrl(pedigree, suspect, exclude, knownNames)
        : null;

      // Build extended pedigree SVG (outside-match pedigree)
      let extendedPedigreeImage: string | null = null;
      if (pedigreeTgf && haplotypesJson) {
        try {
          const extPedigree = await invoke<{ individuals: { id: string; name: string; haplotypeClass: string; exclude: boolean }[]; relationships: { parentId: string; childId: string }[] }>(
            "build_extended_pedigree",
            {
              pedigreeTgf,
              haplotypesJson,
              markerSetName: selectedKitName,
              markerSetCsv,
              suspect,
              traceMode: params.traceMode,
            }
          );
          // Extended pedigree: new_child nodes are unknown; all original nodes retain their class
          const extKnownNames = new Set<string>(
            extPedigree.individuals
              .filter((i) => i.haplotypeClass === "known" || i.haplotypeClass === "suspect")
              .map((i) => i.name)
          );
          extendedPedigreeImage = renderPedigreeSvgDataUrl(
            extPedigree as Parameters<typeof renderPedigreeSvgDataUrl>[0],
            suspect,
            exclude,
            extKnownNames,
          );
        } catch {
          // Non-fatal: report without extended pedigree image
        }
      }
      if (extendedPedigreeImage) chartImages["extended_pedigree"] = extendedPedigreeImage;

      // Enrich pedigree with effective haplotypeClass / exclude flag before
      // passing to the report — the store keeps raw "unknown" for all nodes.
      const enrichedPedigree = pedigree ? {
        ...pedigree,
        individuals: pedigree.individuals.map((ind) => {
          let haplotypeClass: string;
          if (exclude.includes(ind.name)) {
            haplotypeClass = "excluded";
          } else if (suspect === ind.name) {
            haplotypeClass = "suspect";
          } else if (haplotypes?.haplotypeTable[ind.name]) {
            haplotypeClass = "known";
          } else if (haplotypes) {
            haplotypeClass = "unknown";
          } else {
            haplotypeClass = ind.haplotypeClass;
          }
          return { ...ind, haplotypeClass, exclude: exclude.includes(ind.name) };
        }),
      } : null;

      const html = await invoke<string>("generate_report", {
        resultJson: JSON.stringify(simulation.result),
        pedigreeImage,
        chartImages,
        pedigreeJson: enrichedPedigree ? JSON.stringify(enrichedPedigree) : null,
        haplotypesJson: haplotypes ? JSON.stringify(haplotypes) : null,
        markersJson: markers.length > 0 ? JSON.stringify(markers) : null,
        reportDate: new Date().toLocaleDateString("nl-NL"),
        progressEventsJson: JSON.stringify(simulation.progress),
      });

      await invoke<string>("save_and_open_report", {
        html,
        simulationName: simulationName,
      });

      // Auto-save to runs folder if configured
      if (savedSettings.runsFolder) {
        try {
          await invoke("save_run", {
            folder: savedSettings.runsFolder,
            html,
            resultJson: JSON.stringify(simulation.result),
            simulationName: simulationName,
          });
          setAutoSaveMsg("Saved to runs folder ✓");
          setTimeout(() => setAutoSaveMsg(null), 3000);
        } catch (e) {
          console.warn("Auto-save run failed:", e);
          setAutoSaveMsg("Auto-save failed");
          setTimeout(() => setAutoSaveMsg(null), 3000);
        }
      }
    } catch (e) {
      setReportError(String(e));
    } finally {
      setReportGenerating(false);
    }
  };

  const hasProgress = simulation.progress.length > 0;
  const hasResult = !!simulation.result;

  return (
    <div className="p-6 max-w-5xl mx-auto space-y-6">
      <div className="flex items-center justify-between">
        <h1 className="text-2xl font-bold text-gray-900">Run Simulation</h1>
        <div className="flex items-center gap-2">
          <button
            onClick={() => {
              if (pedigree || haplotypes || simulation.result) {
                if (confirm("Start a new session? All unsaved data will be lost.")) resetAll();
              } else {
                resetAll();
              }
            }}
            className="text-sm border border-gray-300 bg-white hover:bg-gray-50 text-gray-700 font-medium py-1.5 px-3 rounded transition-colors"
          >
            New session
          </button>
          <button
            onClick={handleLoadSession}
            disabled={sessionBusy}
            className="text-sm border border-gray-300 bg-white hover:bg-gray-50 disabled:opacity-50 text-gray-700 font-medium py-1.5 px-3 rounded transition-colors"
          >
            Load session…
          </button>
          <button
            onClick={handleSaveSession}
            disabled={sessionBusy || (!pedigreeTgf && !haplotypesJson)}
            className="text-sm border border-gray-300 bg-white hover:bg-gray-50 disabled:opacity-50 text-gray-700 font-medium py-1.5 px-3 rounded transition-colors"
          >
            Save session…
          </button>
        </div>
      </div>

      {sessionError && (
        <div className="rounded bg-red-50 border border-red-200 p-3 text-sm text-red-800">
          Session error: {sessionError}
        </div>
      )}

      {/* Status banners */}
      {!pedigree && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          No pedigree loaded.{" "}
          <button onClick={() => navigate("/pedigree")} className="font-semibold underline hover:text-yellow-900">
            Go to Pedigree →
          </button>
        </div>
      )}
      {pedigree && markers.length === 0 && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          No marker set loaded.{" "}
          <button onClick={() => navigate("/markers")} className="font-semibold underline hover:text-yellow-900">
            Go to Marker Sets →
          </button>
        </div>
      )}
      {pedigree && markers.length > 0 && haplotypes === null && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          No haplotypes loaded.{" "}
          <button onClick={() => navigate("/haplotypes")} className="font-semibold underline hover:text-yellow-900">
            Go to Haplotypes →
          </button>
        </div>
      )}

      {/* Data summary */}
      {(pedigree || haplotypes || markers.length > 0) && (
        <div className="rounded bg-white border border-gray-200 px-4 py-2.5 flex flex-wrap gap-x-5 gap-y-1 text-sm text-gray-600">
          <span>
            <span className="text-gray-400 mr-1">Pedigree:</span>
            {pedigree
              ? <strong className="text-gray-800">{pedigree.individuals.length} individuals</strong>
              : <span className="text-gray-400 italic">none</span>}
          </span>
          <span>
            <span className="text-gray-400 mr-1">Markers:</span>
            {markers.length > 0
              ? <strong className="text-gray-800">{markers.length}{selectedKitName ? ` (${selectedKitName})` : " (custom)"}</strong>
              : <span className="text-gray-400 italic">none</span>}
          </span>
          <span>
            <span className="text-gray-400 mr-1">Haplotypes:</span>
            {haplotypes
              ? <strong className="text-gray-800">{Object.keys(haplotypes.haplotypeTable).length} known</strong>
              : <span className="text-gray-400 italic">none</span>}
          </span>
        </div>
      )}

      {/* Marker overlap warnings */}
      {markers.length > 0 && haplotypes && (() => {
        const markerNames = new Set(markers.map((m) => m.name));
        const haploMarkers = new Set(haplotypes.markerNames);
        const inSetNotInHaplo = markers.filter((m) => !haploMarkers.has(m.name)).map((m) => m.name);
        const inHaploNotInSet = haplotypes.markerNames.filter((m) => !markerNames.has(m));
        return (
          <>
            {inSetNotInHaplo.length > 0 && (
              <div className="rounded bg-orange-50 border border-orange-200 px-4 py-2 text-xs text-orange-800">
                ⚠ {inSetNotInHaplo.length} marker{inSetNotInHaplo.length > 1 ? "s" : ""} in the active set not found in haplotypes (will use zero-frequency fallback):{" "}
                <span className="font-mono">{inSetNotInHaplo.slice(0, 5).join(", ")}{inSetNotInHaplo.length > 5 ? `…+${inSetNotInHaplo.length - 5}` : ""}</span>
              </div>
            )}
            {inHaploNotInSet.length > 0 && (
              <div className="rounded bg-blue-50 border border-blue-200 px-4 py-2 text-xs text-blue-800">
                ℹ {inHaploNotInSet.length} marker{inHaploNotInSet.length > 1 ? "s" : ""} in haplotypes not in active marker set (will be ignored):{" "}
                <span className="font-mono">{inHaploNotInSet.slice(0, 5).join(", ")}{inHaploNotInSet.length > 5 ? `…+${inHaploNotInSet.length - 5}` : ""}</span>
              </div>
            )}
          </>
        );
      })()}

      <div className="grid grid-cols-2 gap-6">
        {/* Left: configuration */}
        <div className="space-y-4">
          <section className="bg-white rounded-lg border p-4 space-y-3">
            <h2 className="font-semibold text-gray-700">Identification</h2>
            <div className="grid grid-cols-2 gap-3 text-sm">
              <div>
                <label className="block text-gray-600 mb-1">Simulation name</label>
                <input
                  type="text"
                  className="w-full border rounded px-2 py-1"
                  value={simulationName}
                  onChange={(e) => setSimulationName(e.target.value)}
                />
              </div>
              <div>
                <label className="block text-gray-600 mb-1">Analyst name</label>
                <input
                  type="text"
                  className="w-full border rounded px-2 py-1"
                  value={userName}
                  onChange={(e) => setUserName(e.target.value)}
                />
              </div>
            </div>
          </section>

          <section className="bg-white rounded-lg border p-4 space-y-3">
            <h2 className="font-semibold text-gray-700">Mode</h2>
            <label
              className={`flex items-center gap-2 text-sm ${!haplotypes?.traceHaplotype ? "opacity-40 cursor-not-allowed" : ""}`}
              title={!haplotypes?.traceHaplotype ? "Load a TRACE profile in the Haplotype Editor to enable trace mode" : "Find the most likely donor for the TRACE profile"}
            >
              <input
                type="checkbox"
                checked={params.traceMode}
                disabled={!haplotypes?.traceHaplotype}
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
            <div className="flex items-center justify-between">
              <h2 className="font-semibold text-gray-700">Parameters</h2>
              <button
                onClick={() => setParams((p) => ({
                  ...p,
                  twoStepMutationFraction: savedSettings.defaultTwoStepFraction,
                  batchLength: savedSettings.defaultBatchLength,
                  convergenceCriterion: savedSettings.defaultConvergenceCriterion,
                  numberOfThreads: savedSettings.defaultThreads,
                  bias: null,
                  seed: null,
                  skipInside: false,
                  skipOutside: false,
                  adaptiveBias: false,
                }))}
                className="text-xs text-gray-400 hover:text-gray-600 border border-gray-200 hover:border-gray-400 rounded px-2 py-0.5 transition-colors"
                title="Reset all parameters to defaults from Settings"
              >
                Reset to defaults
              </button>
            </div>
            <div className="grid grid-cols-2 gap-3 text-sm">
              <div>
                <label className="block text-gray-600 mb-1" title="Fraction of mutations that are two-step (±2 repeat units) rather than one-step (±1). Typical value: 0.03.">Two-step fraction</label>
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
                <label className="block text-gray-600 mb-1" title="Number of Monte Carlo trials per batch. Each batch produces one data point on the convergence chart. Larger values give smoother convergence but slower updates.">Batch length</label>
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
                <label className="block text-gray-600 mb-1" title="Maximum allowed relative spread between the three independent model estimates before the result is accepted. Lower = more precise but slower. Typical: 0.02 (2%).">Convergence criterion</label>
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
                <label className="block text-gray-600 mb-1" title={`Number of parallel threads used for the three-model ensemble. Higher = faster on multi-core machines. Max on this machine: ${cpuCount}.`}>
                  Threads <span className="text-gray-400 font-normal">(max {cpuCount})</span>
                </label>
                <input
                  type="number"
                  min="1"
                  max={cpuCount}
                  className="w-full border rounded px-2 py-1"
                  value={params.numberOfThreads}
                  onChange={(e) =>
                    setParams((p) => ({ ...p, numberOfThreads: parseInt(e.target.value, 10) }))
                  }
                />
              </div>
            </div>
            <div className="flex gap-4 text-sm">
              <label className="flex items-center gap-2" title="Skip the pedigree-probability simulation stage (P inside pedigree). Use when only outside-pedigree matching is needed.">
                <input
                  type="checkbox"
                  checked={params.skipInside}
                  onChange={(e) => setParams((p) => ({ ...p, skipInside: e.target.checked }))}
                />
                Skip inside
              </label>
              <label className="flex items-center gap-2" title="Skip the outside-pedigree simulation stage (P outside pedigree / random match probability). Use when only within-pedigree analysis is needed.">
                <input
                  type="checkbox"
                  checked={params.skipOutside}
                  onChange={(e) => setParams((p) => ({ ...p, skipOutside: e.target.checked }))}
                />
                Skip outside
              </label>
              <label className="flex items-center gap-2" title="Automatically tune the importance-sampling bias parameter during the run to improve simulation efficiency. Recommended for most cases.">
                <input
                  type="checkbox"
                  checked={params.adaptiveBias}
                  onChange={(e) => setParams((p) => ({ ...p, adaptiveBias: e.target.checked }))}
                />
                Adaptive bias
              </label>
            </div>
            {!params.adaptiveBias && (
              <div className="text-sm">
                <label className="block text-gray-600 mb-1" title="Fixed importance-sampling bias factor. Leave blank to let the engine choose automatically. Only relevant when Adaptive bias is off.">
                  Bias (leave blank for auto)
                </label>
                <input
                  type="number"
                  step="0.01"
                  min="0"
                  max="1"
                  placeholder="auto"
                  className="w-32 border rounded px-2 py-1"
                  value={params.bias ?? ""}
                  onChange={(e) =>
                    setParams((p) => ({
                      ...p,
                      bias: e.target.value === "" ? null : parseFloat(e.target.value),
                    }))
                  }
                />
              </div>
            )}
            <div className="text-sm">
              <label className="block text-gray-600 mb-1" title="Random seed for reproducible results. Leave blank to use a different random seed each run.">
                Seed (leave blank for default)
              </label>
              <input
                type="number"
                min="0"
                step="1"
                placeholder="default"
                className="w-32 border rounded px-2 py-1"
                value={params.seed ?? ""}
                onChange={(e) =>
                  setParams((p) => ({
                    ...p,
                    seed: e.target.value === "" ? null : parseInt(e.target.value, 10),
                  }))
                }
              />
            </div>
          </section>

          {allIndividuals.length > 0 && (
            <section className="bg-white rounded-lg border p-4 space-y-2">
              <h2 className="font-semibold text-gray-700">Excluded Individuals</h2>
              <p className="text-xs text-gray-400">
                Individuals excluded from the match probability calculation.
              </p>
              <div className="space-y-1 max-h-36 overflow-y-auto">
                {[...allIndividuals].sort((a, b) => a.name.localeCompare(b.name)).map((ind) => (
                  <label key={ind.id} className="flex items-center gap-2 text-sm">
                    <input
                      type="checkbox"
                      checked={exclude.includes(ind.name)}
                      onChange={(e) => {
                        if (e.target.checked) {
                          setExclude([...exclude, ind.name]);
                        } else {
                          setExclude(exclude.filter((n) => n !== ind.name));
                        }
                      }}
                    />
                    {ind.name}
                  </label>
                ))}
              </div>
            </section>
          )}

          {autoSaveMsg && (
            <div className={`rounded px-3 py-2 text-xs font-medium ${autoSaveMsg.includes("failed") ? "bg-red-50 text-red-700 border border-red-200" : "bg-emerald-50 text-emerald-700 border border-emerald-200"}`}>
              {autoSaveMsg}
            </div>
          )}

          <div className="flex gap-3">
            <button
              disabled={!canRun || simulation.running}
              onClick={handleRun}
              title={
                !pedigree ? "Missing: pedigree" :
                haplotypes === null ? "Missing: haplotypes" :
                markers.length === 0 ? "Missing: marker set" :
                simulation.running ? "Simulation in progress" :
                "Run simulation (Ctrl+Enter)"
              }
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

          {simulation.running && !hasProgress && (
            <section className="bg-white rounded-lg border p-5 flex items-center gap-3 text-sm text-gray-500">
              <div className="w-5 h-5 rounded-full border-2 border-blue-500 border-t-transparent animate-spin flex-shrink-0" />
              Starting simulation — waiting for first batch…
            </section>
          )}

          {hasProgress && (
            <section className="bg-white rounded-lg border p-4 space-y-3">
              <h2 className="font-semibold text-gray-700">Convergence</h2>
              <ConvergenceChart
                ref={pedigreeChartRef}
                events={simulation.progress}
                stage="pedigree_probability"
                title="Pedigree probability"
                convergenceCriterion={params.convergenceCriterion}
              />
              {!params.skipInside && (
                <ConvergenceChart
                  ref={insideChartRef}
                  events={simulation.progress}
                  stage="inside_match_probability"
                  title="Inside-pedigree match"
                  convergenceCriterion={params.convergenceCriterion}
                />
              )}
              {!params.skipOutside && (
                <ConvergenceChart
                  ref={extendedPedigreeChartRef}
                  events={simulation.progress}
                  stage="extended_pedigree_probability"
                  title="Extended pedigree probability"
                  convergenceCriterion={params.convergenceCriterion}
                />
              )}
              {!params.skipOutside && (
                <ConvergenceChart
                  ref={outsideChartRef}
                  events={simulation.progress}
                  stage="outside_match_probability"
                  title="Outside-pedigree match"
                  convergenceCriterion={params.convergenceCriterion}
                />
              )}
            </section>
          )}

          {hasResult && simulation.result && (
            <section className="bg-white rounded-lg border p-4 space-y-3">
              <div className="flex items-center justify-between">
                <h2 className="font-semibold text-gray-700">Results</h2>
                <span
                  title={`${simulation.result.trials} convergence batches × ${params.batchLength} = ${(simulation.result.trials * params.batchLength).toLocaleString()} MC samples`}
                  className={`text-xs font-medium px-2 py-0.5 rounded-full cursor-default ${simulation.result.converged ? "bg-green-100 text-green-700" : "bg-orange-100 text-orange-700"}`}
                >
                  {simulation.result.converged ? "Converged" : "Not converged"} · {simulation.result.trials} batches
                </span>
              </div>

              {/* Headline probability cards */}
              {(() => {
                const perInd = simulation.result.per_individual_probabilities;
                const lrValues = perInd
                  ? Object.values(perInd)
                      .map((p) => parseFloat(p as string))
                      .filter((n) => isFinite(n) && n > 0)
                      .map((n) => 1 / n)
                  : [];
                const avgLr = lrValues.length > 0
                  ? lrValues.reduce((a, b) => a + b, 0) / lrValues.length
                  : null;
                return (
                  <div className="grid grid-cols-2 gap-2">
                    <div className="rounded-lg bg-blue-50 border border-blue-100 px-3 py-2">
                      <p className="text-xs text-blue-500 mb-0.5">Pedigree probability</p>
                      <p className="text-lg font-bold text-blue-800 font-mono">
                        {simulation.result.inside_match_probabilities
                          ? fmtPct(simulation.result.inside_match_probabilities.average_pedigree_probability)
                          : "—"}
                      </p>
                    </div>
                    {simulation.result.inside_match_probabilities && (
                      <div className="rounded-lg bg-emerald-50 border border-emerald-100 px-3 py-2">
                        <p className="text-xs text-emerald-600 mb-0.5">Inside-pedigree match</p>
                        <p className="text-lg font-bold text-emerald-800 font-mono">
                          {fmtPct(
                            Object.values(simulation.result.inside_match_probabilities.probabilities ?? {})
                              .reduce((s, v) => s + parseFloat(v as string), 0)
                          )}
                        </p>
                      </div>
                    )}
                    {simulation.result.outside_match_probability && (
                      <div className="rounded-lg bg-purple-50 border border-purple-100 px-3 py-2">
                        <p className="text-xs text-purple-600 mb-0.5">Outside-pedigree match</p>
                        <p className="text-lg font-bold text-purple-800 font-mono">
                          {fmtPct(simulation.result.outside_match_probability)}
                        </p>
                      </div>
                    )}
                    {avgLr !== null && (
                      <div className="rounded-lg bg-amber-50 border border-amber-100 px-3 py-2">
                        <p className="text-xs text-amber-600 mb-0.5">Average LR</p>
                        <p className="text-lg font-bold text-amber-800 font-mono">
                          {fmtLr(avgLr)}
                        </p>
                      </div>
                    )}
                  </div>
                );
              })()}

              {/* Per-individual match probabilities — the main table of interest */}
              {simulation.result.per_individual_probabilities &&
                Object.keys(simulation.result.per_individual_probabilities).length > 0 && (() => {
                  const sorted = Object.entries(simulation.result.per_individual_probabilities!)
                    .sort(([, a], [, b]) => parseFloat(b as string) - parseFloat(a as string));
                  const lrList = sorted.map(([, p]) => {
                    const n = parseFloat(p as string);
                    return isFinite(n) && n > 0 ? 1 / n : null;
                  });
                  const validLrs = lrList.filter((v): v is number => v !== null);
                  const avgLr = validLrs.length > 0
                    ? validLrs.reduce((a, b) => a + b, 0) / validLrs.length
                    : null;
                  return (
                    <div className="text-sm">
                      <p className="font-medium text-gray-700 mb-1">Match probability per individual</p>
                      <table className="text-xs w-full border-collapse">
                        <thead>
                          <tr className="bg-gray-50">
                            <th className="border px-2 py-1 text-left">Individual</th>
                            <th className="border px-2 py-1 text-right">P(match)</th>
                            <th className="border px-2 py-1 text-right">LR (1/P)</th>
                          </tr>
                        </thead>
                        <tbody>
                          {sorted.map(([name, prob], i) => (
                            <tr key={name} className="hover:bg-gray-50">
                              <td className="border px-2 py-1 font-medium">{name}</td>
                              <td className="border px-2 py-1 text-right font-mono">
                                {fmt2sig(parseFloat(prob as string))}
                              </td>
                              <td className="border px-2 py-1 text-right font-mono">
                                {fmtLr(lrList[i])}
                              </td>
                            </tr>
                          ))}
                          {avgLr !== null && (
                            <tr className="bg-amber-50 font-semibold">
                              <td className="border px-2 py-1 text-gray-600">Average LR</td>
                              <td className="border px-2 py-1" />
                              <td className="border px-2 py-1 text-right font-mono text-amber-800">
                                {fmtLr(avgLr)}
                              </td>
                            </tr>
                          )}
                        </tbody>
                      </table>
                    </div>
                  );
                })()}

              {!simulation.result.inside_match_probabilities &&
                !simulation.result.outside_match_probability && (
                  <p className="text-sm text-gray-500">
                    No suspect specified — pedigree probability only.
                  </p>
                )}

              <button
                onClick={() => {
                  const r = simulation.result!;
                  const lines: string[] = [`Simulation: ${simulationName || "—"} | Analyst: ${userName || "—"}`];
                  lines.push(`Converged: ${r.converged ? "Yes" : "No"} | Batches: ${r.trials}`);
                  if (r.inside_match_probabilities) {
                    lines.push(`Pedigree probability: ${fmtPct(r.inside_match_probabilities.average_pedigree_probability)}`);
                    const insideSum = Object.values(r.inside_match_probabilities.probabilities ?? {})
                      .reduce((s, v) => s + parseFloat(v as string), 0);
                    lines.push(`Inside-pedigree match: ${fmtPct(insideSum)}`);
                  }
                  if (r.outside_match_probability) lines.push(`Outside-pedigree match: ${fmtPct(r.outside_match_probability)}`);
                  if (r.per_individual_probabilities) {
                    const sorted = Object.entries(r.per_individual_probabilities)
                      .sort(([, a], [, b]) => parseFloat(b as string) - parseFloat(a as string));
                    const lrs = sorted.map(([, p]) => { const n = parseFloat(p as string); return isFinite(n) && n > 0 ? 1 / n : null; });
                    const validLrs = lrs.filter((v): v is number => v !== null);
                    const avgLr = validLrs.length > 0 ? validLrs.reduce((a, b) => a + b, 0) / validLrs.length : null;
                    lines.push("Per-individual match probabilities:");
                    sorted.forEach(([name, prob], i) => {
                      const lr = lrs[i];
                      lines.push(`  ${name}: P=${fmt2sig(parseFloat(prob as string))}  LR=${fmtLr(lr)}`);
                    });
                    if (avgLr !== null) lines.push(`  Average LR: ${fmtLr(avgLr)}`);
                  }
                  navigator.clipboard.writeText(lines.join("\n"));
                  setCopyMsg(true);
                  setTimeout(() => setCopyMsg(false), 2000);
                }}
                className="text-xs border border-gray-200 rounded px-3 py-1.5 hover:bg-gray-50 text-gray-600 transition-colors"
              >
                {copyMsg ? "Copied ✓" : "Copy results"}
              </button>

              {reportError && (
                <div className="rounded bg-red-50 border border-red-200 p-2 text-xs text-red-700">
                  {reportError}
                </div>
              )}
              <div className="flex gap-2 mt-2">
                <button
                  onClick={handleGenerateReport}
                  disabled={reportGenerating}
                  className="flex-1 bg-emerald-600 hover:bg-emerald-700 disabled:bg-gray-300 text-white font-medium py-2 px-4 rounded text-sm transition-colors"
                >
                  {reportGenerating ? "Generating…" : "Generate Report"}
                </button>
                {simulation.result?.per_individual_probabilities && (
                  <button
                    onClick={() => navigate("/pedigree", { state: { showProbs: true } })}
                    title="Show per-individual match probabilities on the pedigree nodes"
                    className="flex-1 bg-indigo-600 hover:bg-indigo-700 text-white font-medium py-2 px-3 rounded text-sm transition-colors"
                  >
                    Show in pedigree
                  </button>
                )}
              </div>
            </section>
          )}
        </div>
      </div>
    </div>
  );
}
