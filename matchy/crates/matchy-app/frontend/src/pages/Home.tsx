import { useRef, useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import { useSimulation } from "../hooks/useSimulation";
import { useSettings } from "./Settings";
import { useT, getTranslations } from "../i18n";
import {
  ConvergenceChart,
  type ConvergenceChartRef,
} from "../components/ConvergenceChart";
import { renderPedigreeSvgDataUrl } from "../utils/pedigreeSvg";
import { saveSession, loadSession } from "../utils/session";
import type { PedigreeData } from "../types/matchy";

// 3 significant figures: fixed notation for [1e-3, 1e3), scientific otherwise.
function fmt2sig(n: number): string {
  if (n === 0) return "0";
  const abs = Math.abs(n);
  if (abs >= 1e-3 && abs < 1e3) {
    return parseFloat(n.toPrecision(3)).toString();
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
  const { pedigree, pedigreeTgf, haplotypes, haplotypesJson, markers, markerSetCsv, selectedKitName, suspect, setSuspect, exclude, setExclude, simulationName, setSimulationName, userName, setUserName, simulation, setPedigree, setHaplotypes, setMarkerSet, setSimulationResult, setSimulationProgress, resetAll, locale, simParams: params, setSimParams: setParams } =
    useAppStore();
  const { startSimulation, cancelSimulation } = useSimulation();
  const savedSettings = useSettings();
  const t = useT();
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

  const canRun = !!pedigree && haplotypes !== null && markers.length > 0 && (params.traceMode || !!suspect);

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
      const svgClassLabels = {
        unknown: t("ped_class_unknown"),
        known: t("ped_class_known"),
        suspect: t("ped_class_suspect"),
        excluded: t("ped_class_excluded"),
      };
      const pedigreeImage = pedigree
        ? renderPedigreeSvgDataUrl(pedigree, suspect, exclude, knownNames, svgClassLabels)
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
            svgClassLabels,
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
        langJson: JSON.stringify(getTranslations(locale)),
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
          setAutoSaveMsg(t("run_auto_saved"));
          setTimeout(() => setAutoSaveMsg(null), 3000);
        } catch (e) {
          console.warn("Auto-save run failed:", e);
          setAutoSaveMsg(t("run_auto_save_failed"));
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
        <h1 className="text-2xl font-bold text-gray-900">{t("run_title")}</h1>
        <div className="flex items-center gap-2">
          <button
            onClick={() => {
              if (pedigree || haplotypes || simulation.result) {
                if (confirm(t("run_new_session") + "?")) resetAll();
              } else {
                resetAll();
              }
            }}
            className="text-sm border border-gray-300 bg-white hover:bg-gray-50 text-gray-700 font-medium py-1.5 px-3 rounded transition-colors"
          >
            {t("run_new_session")}
          </button>
          <button
            onClick={handleLoadSession}
            disabled={sessionBusy}
            className="text-sm border border-gray-300 bg-white hover:bg-gray-50 disabled:opacity-50 text-gray-700 font-medium py-1.5 px-3 rounded transition-colors"
          >
            {t("run_load_session")}
          </button>
          <button
            onClick={handleSaveSession}
            disabled={sessionBusy || (!pedigreeTgf && !haplotypesJson)}
            className="text-sm border border-gray-300 bg-white hover:bg-gray-50 disabled:opacity-50 text-gray-700 font-medium py-1.5 px-3 rounded transition-colors"
          >
            {t("run_save_session")}
          </button>
        </div>
      </div>

      {sessionError && (
        <div className="rounded bg-red-50 border border-red-200 p-3 text-sm text-red-800">
          {t("run_session_error")} {sessionError}
        </div>
      )}

      {/* Status banners */}
      {!pedigree && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          {t("run_no_pedigree")}{" "}
          <button onClick={() => navigate("/pedigree")} className="font-semibold underline hover:text-yellow-900">
            {t("run_go_pedigree")}
          </button>
        </div>
      )}
      {pedigree && markers.length === 0 && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          {t("run_no_markers")}{" "}
          <button onClick={() => navigate("/markers")} className="font-semibold underline hover:text-yellow-900">
            {t("run_go_markers")}
          </button>
        </div>
      )}
      {pedigree && markers.length > 0 && haplotypes === null && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          {t("run_no_haplotypes")}{" "}
          <button onClick={() => navigate("/haplotypes")} className="font-semibold underline hover:text-yellow-900">
            {t("run_go_haplotypes")}
          </button>
        </div>
      )}
      {pedigree && markers.length > 0 && haplotypes !== null && !params.traceMode && !suspect && (
        <div className="rounded bg-yellow-50 border border-yellow-200 p-3 text-sm text-yellow-800">
          {t("run_no_suspect_warning")}
        </div>
      )}

      {/* Data summary */}
      {(pedigree || haplotypes || markers.length > 0) && (
        <div className="rounded bg-white border border-gray-200 px-4 py-2.5 flex flex-wrap gap-x-5 gap-y-1 text-sm text-gray-600">
          <span>
            <span className="text-gray-400 mr-1">{t("run_pedigree")}:</span>
            {pedigree
              ? <strong className="text-gray-800">{pedigree.individuals.length} {t("run_individuals")}</strong>
              : <span className="text-gray-400 italic">{t("run_none")}</span>}
          </span>
          <span>
            <span className="text-gray-400 mr-1">{t("run_markers")}:</span>
            {markers.length > 0
              ? <strong className="text-gray-800">{markers.length}{selectedKitName ? ` (${selectedKitName})` : ` (${t("run_custom")})`}</strong>
              : <span className="text-gray-400 italic">{t("run_none")}</span>}
          </span>
          <span>
            <span className="text-gray-400 mr-1">{t("run_haplotypes")}:</span>
            {haplotypes
              ? <strong className="text-gray-800">{Object.keys(haplotypes.haplotypeTable).length} {t("run_known")}</strong>
              : <span className="text-gray-400 italic">{t("run_none")}</span>}
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
            <h2 className="font-semibold text-gray-700">{t("run_section_identification")}</h2>
            <div className="grid grid-cols-2 gap-3 text-sm">
              <div>
                <label className="block text-gray-600 mb-1">{t("run_simulation_name")}</label>
                <input
                  type="text"
                  className="w-full border rounded px-2 py-1"
                  value={simulationName}
                  onChange={(e) => setSimulationName(e.target.value)}
                />
              </div>
              <div>
                <label className="block text-gray-600 mb-1">{t("run_analyst_name")}</label>
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
            <h2 className="font-semibold text-gray-700">{t("run_section_mode")}</h2>
            <label
              className={`flex items-center gap-2 text-sm ${!haplotypes?.traceHaplotype ? "opacity-40 cursor-not-allowed" : ""}`}
              title={!haplotypes?.traceHaplotype ? t("run_trace_mode_tooltip") : t("run_trace_mode")}
            >
              <input
                type="checkbox"
                checked={params.traceMode}
                disabled={!haplotypes?.traceHaplotype}
                onChange={(e) => setParams({ ...params, traceMode: e.target.checked })}
              />
              {t("run_trace_mode")}
            </label>

            {!params.traceMode && (
              <div>
                <label className="block text-sm font-medium text-gray-600 mb-1">{t("run_suspect")}</label>
                <select
                  className="w-full border rounded px-2 py-1.5 text-sm"
                  value={suspect ?? ""}
                  onChange={(e) => setSuspect(e.target.value || null)}
                >
                  <option value="">{t("run_select")}</option>
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
              <h2 className="font-semibold text-gray-700">{t("run_section_parameters")}</h2>
              <button
                onClick={() => setParams({
                  ...params,
                  twoStepMutationFraction: savedSettings.defaultTwoStepFraction,
                  batchLength: savedSettings.defaultBatchLength,
                  convergenceCriterion: savedSettings.defaultConvergenceCriterion,
                  numberOfThreads: savedSettings.defaultThreads,
                  bias: null,
                  seed: null,
                  skipInside: false,
                  skipOutside: false,
                  adaptiveBias: false,
                })}
                className="text-xs text-gray-400 hover:text-gray-600 border border-gray-200 hover:border-gray-400 rounded px-2 py-0.5 transition-colors"
                title={t("run_reset_defaults")}
              >
                {t("run_reset_defaults")}
              </button>
            </div>
            <div className="grid grid-cols-2 gap-3 text-sm">
              <div>
                <label className="block text-gray-600 mb-1" title={t("run_tooltip_two_step")}>{t("run_two_step_fraction")}</label>
                <input
                  type="number"
                  step="0.001"
                  min="0"
                  max="1"
                  className="w-full border rounded px-2 py-1"
                  value={params.twoStepMutationFraction}
                  onChange={(e) =>
                    setParams({ ...params, twoStepMutationFraction: parseFloat(e.target.value) })
                  }
                />
              </div>
              <div>
                <label className="block text-gray-600 mb-1" title={t("run_tooltip_batch_length")}>{t("run_batch_length")}</label>
                <input
                  type="number"
                  min="100"
                  className="w-full border rounded px-2 py-1"
                  value={params.batchLength}
                  onChange={(e) =>
                    setParams({ ...params, batchLength: parseInt(e.target.value, 10) })
                  }
                />
              </div>
              <div>
                <label className="block text-gray-600 mb-1" title={t("run_tooltip_convergence")}>{t("run_convergence_criterion")}</label>
                <input
                  type="number"
                  step="0.001"
                  min="0"
                  className="w-full border rounded px-2 py-1"
                  value={params.convergenceCriterion}
                  onChange={(e) =>
                    setParams({ ...params, convergenceCriterion: parseFloat(e.target.value) })
                  }
                />
              </div>
              <div>
                <label className="block text-gray-600 mb-1" title={t("run_tooltip_threads").replace("{n}", String(cpuCount))}>
                  {t("run_threads")} <span className="text-gray-400 font-normal">(max {cpuCount})</span>
                </label>
                <input
                  type="number"
                  min="1"
                  max={cpuCount}
                  className="w-full border rounded px-2 py-1"
                  value={params.numberOfThreads}
                  onChange={(e) =>
                    setParams({ ...params, numberOfThreads: parseInt(e.target.value, 10) })
                  }
                />
              </div>
            </div>
            <div className="flex gap-4 text-sm">
              <label className="flex items-center gap-2" title={t("run_tooltip_skip_inside")}>
                <input
                  type="checkbox"
                  checked={params.skipInside}
                  onChange={(e) => setParams({ ...params, skipInside: e.target.checked })}
                />
                {t("run_skip_inside")}
              </label>
              <label className="flex items-center gap-2" title={t("run_tooltip_skip_outside")}>
                <input
                  type="checkbox"
                  checked={params.skipOutside}
                  onChange={(e) => setParams({ ...params, skipOutside: e.target.checked })}
                />
                {t("run_skip_outside")}
              </label>
              <label className="flex items-center gap-2" title={t("run_tooltip_adaptive_bias")}>
                <input
                  type="checkbox"
                  checked={params.adaptiveBias}
                  onChange={(e) => setParams({ ...params, adaptiveBias: e.target.checked })}
                />
                {t("run_adaptive_bias")}
              </label>
            </div>
            {!params.adaptiveBias && (
              <div className="text-sm">
                <label className="block text-gray-600 mb-1" title={t("run_tooltip_bias")}>
                  {t("run_bias")}
                </label>
                <input
                  type="number"
                  step="0.01"
                  min="0"
                  max="1"
                  placeholder={t("run_bias_auto")}
                  className="w-32 border rounded px-2 py-1"
                  value={params.bias ?? ""}
                  onChange={(e) =>
                    setParams({ ...params, bias: e.target.value === "" ? null : parseFloat(e.target.value) })
                  }
                />
              </div>
            )}
            <div className="text-sm">
              <label className="block text-gray-600 mb-1" title={t("run_tooltip_seed")}>
                {t("run_seed")}
              </label>
              <input
                type="number"
                min="0"
                step="1"
                placeholder={t("run_seed_default")}
                className="w-32 border rounded px-2 py-1"
                value={params.seed ?? ""}
                onChange={(e) =>
                  setParams({ ...params, seed: e.target.value === "" ? null : parseInt(e.target.value, 10) })
                }
              />
            </div>
          </section>

          {allIndividuals.length > 0 && (
            <section className="bg-white rounded-lg border p-4 space-y-2">
              <h2 className="font-semibold text-gray-700">{t("run_exclude_section")}</h2>
              <p className="text-xs text-gray-400">{t("run_exclude_desc")}</p>
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
            <div className={`rounded px-3 py-2 text-xs font-medium ${autoSaveMsg.includes("failed") || autoSaveMsg.includes("mislukt") || autoSaveMsg.includes("fehlgeschlagen") ? "bg-red-50 text-red-700 border border-red-200" : "bg-emerald-50 text-emerald-700 border border-emerald-200"}`}>
              {autoSaveMsg}
            </div>
          )}

          <div className="flex gap-3">
            <button
              disabled={!canRun || simulation.running}
              onClick={handleRun}
              title={
                !pedigree ? t("run_no_pedigree") :
                haplotypes === null ? t("run_no_haplotypes") :
                markers.length === 0 ? t("run_no_markers") :
                !params.traceMode && !suspect ? t("run_no_suspect_warning") :
                simulation.running ? t("run_running") :
                "Ctrl+Enter"
              }
              className="flex-1 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-300 text-white font-medium py-2 px-4 rounded transition-colors"
            >
              {simulation.running ? t("run_running") : t("run_start")}
            </button>
            {simulation.running && (
              <button
                onClick={cancelSimulation}
                className="bg-red-100 hover:bg-red-200 text-red-700 font-medium py-2 px-4 rounded"
              >
                {t("run_cancel")}
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
              {t("run_starting")}
            </section>
          )}

          {hasProgress && (
            <section className="bg-white rounded-lg border p-4 space-y-3">
              <h2 className="font-semibold text-gray-700">{t("run_convergence_title")}</h2>
              <ConvergenceChart
                ref={pedigreeChartRef}
                events={simulation.progress}
                stage="pedigree_probability"
                title={t("run_pedigree_prob_card")}
                convergenceCriterion={params.convergenceCriterion}
                batchLength={params.batchLength}
              />
              {!params.skipInside && (
                <ConvergenceChart
                  ref={insideChartRef}
                  events={simulation.progress}
                  stage="inside_match_probability"
                  title={t("run_inside_match_card")}
                  convergenceCriterion={params.convergenceCriterion}
                  batchLength={params.batchLength}
                />
              )}
              {!params.skipOutside && (
                <ConvergenceChart
                  ref={extendedPedigreeChartRef}
                  events={simulation.progress}
                  stage="extended_pedigree_probability"
                  title={t("run_ext_pedigree_card")}
                  convergenceCriterion={params.convergenceCriterion}
                  batchLength={params.batchLength}
                />
              )}
              {!params.skipOutside && (
                <ConvergenceChart
                  ref={outsideChartRef}
                  events={simulation.progress}
                  stage="outside_match_probability"
                  title={t("run_outside_match_card")}
                  convergenceCriterion={params.convergenceCriterion}
                  batchLength={params.batchLength}
                />
              )}
            </section>
          )}

          {hasResult && simulation.result && (
            <section className="bg-white rounded-lg border p-4 space-y-3">
              <div className="flex items-center justify-between">
                <h2 className="font-semibold text-gray-700">{t("run_results_title")}</h2>
                <span
                  title={`${simulation.result.trials} ${t("run_convergence_title")} ${t("run_batches")} × ${params.batchLength} = ${(simulation.result.trials * params.batchLength).toLocaleString()} MC samples`}
                  className={`text-xs font-medium px-2 py-0.5 rounded-full cursor-default ${simulation.result.converged ? "bg-green-100 text-green-700" : "bg-orange-100 text-orange-700"}`}
                >
                  {simulation.result.converged ? t("run_converged") : t("run_not_converged")} · {simulation.result.trials} {t("run_batches")}
                </span>
              </div>

              {/* Headline probability cards */}
              {(() => {
                const perInd = simulation.result.per_individual_probabilities;
                const probs = perInd
                  ? Object.values(perInd)
                      .map((p) => parseFloat(p as string))
                      .filter((n) => isFinite(n) && n > 0)
                  : [];
                const avgLr = probs.length > 0
                  ? probs.length / probs.reduce((a, b) => a + b, 0)
                  : null;
                return (
                  <div className="grid grid-cols-2 gap-2">
                    <div className="rounded-lg bg-blue-50 border border-blue-100 px-3 py-2">
                      <p className="text-xs text-blue-500 mb-0.5">{t("run_pedigree_prob_card")}</p>
                      <p className="text-lg font-bold text-blue-800 font-mono">
                        {simulation.result.inside_match_probabilities
                          ? fmtPct(simulation.result.inside_match_probabilities.average_pedigree_probability)
                          : "—"}
                      </p>
                    </div>
                    {simulation.result.inside_match_probabilities && (
                      <div className="rounded-lg bg-emerald-50 border border-emerald-100 px-3 py-2">
                        <p className="text-xs text-emerald-600 mb-0.5">{t("run_inside_match_card")}</p>
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
                        <p className="text-xs text-purple-600 mb-0.5">{t("run_outside_match_card")}</p>
                        <p className="text-lg font-bold text-purple-800 font-mono">
                          {fmtPct(simulation.result.outside_match_probability)}
                        </p>
                      </div>
                    )}
                    {avgLr !== null && (
                      <div className="rounded-lg bg-amber-50 border border-amber-100 px-3 py-2">
                        <p className="text-xs text-amber-600 mb-0.5">{t("run_avg_lr")}</p>
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
                  const validProbs = sorted
                    .map(([, p]) => parseFloat(p as string))
                    .filter((n) => isFinite(n) && n > 0);
                  const avgLr = validProbs.length > 0
                    ? validProbs.length / validProbs.reduce((a, b) => a + b, 0)
                    : null;
                  return (
                    <div className="text-sm">
                      <p className="font-medium text-gray-700 mb-1">{t("run_per_individual_title")}</p>
                      <table className="text-xs w-full border-collapse">
                        <thead>
                          <tr className="bg-gray-50">
                            <th className="border px-2 py-1 text-left">{t("run_individual")}</th>
                            <th className="border px-2 py-1 text-right">{t("run_match_probability")}</th>
                            <th className="border px-2 py-1 text-right">{t("run_lr")}</th>
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
                              <td className="border px-2 py-1 text-gray-600">{t("run_avg_lr")}</td>
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
                    {t("run_no_suspect_msg")}
                  </p>
                )}

              <button
                onClick={() => {
                  const r = simulation.result!;
                  const lines: string[] = [`${t("run_simulation_name")}: ${simulationName || "—"} | ${t("run_analyst_name")}: ${userName || "—"}`];
                  lines.push(`${t("run_converged")}: ${r.converged ? "Yes" : "No"} | ${t("run_batches")}: ${r.trials}`);
                  if (r.inside_match_probabilities) {
                    lines.push(`${t("run_pedigree_prob_card")}: ${fmtPct(r.inside_match_probabilities.average_pedigree_probability)}`);
                    const insideSum = Object.values(r.inside_match_probabilities.probabilities ?? {})
                      .reduce((s, v) => s + parseFloat(v as string), 0);
                    lines.push(`${t("run_inside_match_card")}: ${fmtPct(insideSum)}`);
                  }
                  if (r.outside_match_probability) lines.push(`${t("run_outside_match_card")}: ${fmtPct(r.outside_match_probability)}`);
                  if (r.per_individual_probabilities) {
                    const sorted = Object.entries(r.per_individual_probabilities)
                      .sort(([, a], [, b]) => parseFloat(b as string) - parseFloat(a as string));
                    const lrs = sorted.map(([, p]) => { const n = parseFloat(p as string); return isFinite(n) && n > 0 ? 1 / n : null; });
                    const clipProbs = sorted.map(([, p]) => parseFloat(p as string)).filter((n) => isFinite(n) && n > 0);
                    const avgLr = clipProbs.length > 0 ? clipProbs.length / clipProbs.reduce((a, b) => a + b, 0) : null;
                    lines.push(`${t("run_per_individual_title")}:`);
                    sorted.forEach(([name, prob], i) => {
                      const lr = lrs[i];
                      lines.push(`  ${name}: P=${fmt2sig(parseFloat(prob as string))}  LR=${fmtLr(lr)}`);
                    });
                    if (avgLr !== null) lines.push(`  ${t("run_avg_lr")}: ${fmtLr(avgLr)}`);
                  }
                  navigator.clipboard.writeText(lines.join("\n"));
                  setCopyMsg(true);
                  setTimeout(() => setCopyMsg(false), 2000);
                }}
                className="text-xs border border-gray-200 rounded px-3 py-1.5 hover:bg-gray-50 text-gray-600 transition-colors"
              >
                {copyMsg ? t("run_copied") : t("run_copy")}
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
                  {reportGenerating ? t("run_generating") : t("run_generate_report")}
                </button>
                {simulation.result?.per_individual_probabilities && (
                  <button
                    onClick={() => navigate("/pedigree", { state: { showProbs: true } })}
                    title={t("run_show_in_pedigree")}
                    className="flex-1 bg-indigo-600 hover:bg-indigo-700 text-white font-medium py-2 px-3 rounded text-sm transition-colors"
                  >
                    {t("run_show_in_pedigree")}
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
