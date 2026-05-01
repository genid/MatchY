import { useEffect, useMemo, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import { useT } from "../i18n";
import type { MarkerInfo } from "../types/matchy";

const SAVED_SETS_KEY = "matchy_saved_marker_sets";
const MARKER_POOL_KEY = "matchy_marker_pool_extras";

interface SavedMarkerSet {
  name: string;
  markers: MarkerInfo[];
}

interface MarkerPoolExtras {
  customMarkers: MarkerInfo[];
  rateOverrides: Record<string, number>;
  copiesOverrides: Record<string, number>;
}

function loadSavedSets(): SavedMarkerSet[] {
  try {
    const stored = localStorage.getItem(SAVED_SETS_KEY);
    if (stored) return JSON.parse(stored) as SavedMarkerSet[];
  } catch {}
  return [];
}

function persistSavedSets(sets: SavedMarkerSet[]) {
  localStorage.setItem(SAVED_SETS_KEY, JSON.stringify(sets));
}

function loadMarkerPoolExtras(): MarkerPoolExtras {
  try {
    const stored = localStorage.getItem(MARKER_POOL_KEY);
    if (stored) return JSON.parse(stored) as MarkerPoolExtras;
  } catch {}
  return { customMarkers: [], rateOverrides: {}, copiesOverrides: {} };
}

function persistMarkerPoolExtras(extras: MarkerPoolExtras) {
  localStorage.setItem(MARKER_POOL_KEY, JSON.stringify(extras));
}

export default function MarkerSets() {
  const { selectedKitName, markers, setMarkerSet, haplotypes } = useAppStore();
  const t = useT();
  const [kitNames, setKitNames] = useState<string[]>([]);
  const [markerPool, setMarkerPool] = useState<MarkerInfo[]>([]);
  const [loading, setLoading] = useState(false);
  const [poolLoading, setPoolLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Custom builder state
  const [customMode, setCustomMode] = useState(false);
  const [selected, setSelected] = useState<Set<string>>(new Set());
  const [rateOverrides, setRateOverrides] = useState<Map<string, number>>(() => {
    const extras = loadMarkerPoolExtras();
    return new Map(Object.entries(extras.rateOverrides));
  });
  const [copiesOverrides, setCopiesOverrides] = useState<Map<string, number>>(() => {
    const extras = loadMarkerPoolExtras();
    return new Map(Object.entries(extras.copiesOverrides));
  });
  const [customMarkers, setCustomMarkers] = useState<MarkerInfo[]>(() => loadMarkerPoolExtras().customMarkers);
  const [search, setSearch] = useState("");

  // Add-to-pool form
  const [newName, setNewName] = useState("");
  const [newRate, setNewRate] = useState("0.004");
  const [newCopies, setNewCopies] = useState("1");

  // Saved sets
  const [savedSets, setSavedSets] = useState<SavedMarkerSet[]>(loadSavedSets);
  const [saveSetName, setSaveSetName] = useState("");
  const [saveSetDialogOpen, setSaveSetDialogOpen] = useState(false);

  useEffect(() => {
    invoke<string[]>("list_kits")
      .then((names) => setKitNames(names.sort()))
      .catch((e) => setError(String(e)));

    setPoolLoading(true);
    invoke<MarkerInfo[]>("list_all_markers")
      .then((pool) => {
        const extras = loadMarkerPoolExtras();
        const baseNames = new Set(pool.map((m) => m.name));
        const merged = [
          ...pool,
          ...extras.customMarkers.filter((m) => !baseNames.has(m.name)),
        ].sort((a, b) => a.name.localeCompare(b.name));
        setMarkerPool(merged);
        setPoolLoading(false);
      })
      .catch((e) => { setError(String(e)); setPoolLoading(false); });
  }, []);

  const applySelection = (
    sel: Set<string>,
    rates: Map<string, number>,
    copies: Map<string, number>,
    pool: MarkerInfo[]
  ) => {
    const active: MarkerInfo[] = [...sel].map((name) => {
      const base = pool.find((m) => m.name === name);
      const detected = detectedCopies.get(name);
      const effectiveCopies = typeof detected === "number"
        ? detected
        : copies.get(name) ?? base?.numberOfCopies ?? 1;
      return {
        name,
        mutationRate: rates.get(name) ?? base?.mutationRate ?? 0.004,
        numberOfCopies: effectiveCopies,
      };
    });
    const csv =
      active.length > 0
        ? "marker,mutation_rate,number_of_copies\n" +
          active.map((m) => `${m.name},${m.mutationRate},${m.numberOfCopies ?? 1}`).join("\n")
        : null;
    setMarkerSet(null, active, csv);
  };

  const handleSelectKit = async (name: string) => {
    if (selectedKitName === name) return;
    setCustomMode(false);
    setLoading(true);
    setError(null);
    try {
      const markerInfos = await invoke<MarkerInfo[]>("load_kit", { name });
      setMarkerSet(name, markerInfos, null);
    } catch (e) {
      setError(String(e));
    } finally {
      setLoading(false);
    }
  };

  const handleEnterCustomMode = () => {
    setCustomMode(true);
    setError(null);
    // Pre-populate from current custom set (if not a built-in kit)
    if (!selectedKitName && markers.length > 0) {
      const sel = new Set(markers.map((m) => m.name));
      const rates = new Map(markers.map((m) => [m.name, m.mutationRate]));
      const cops = new Map(markers.map((m) => [m.name, m.numberOfCopies ?? 1]));
      setSelected(sel);
      setRateOverrides(rates);
      setCopiesOverrides(cops);
    } else {
      setSelected(new Set());
      // Keep persisted rate/copies overrides so custom edits survive re-entering custom mode
    }
  };

  const toggleMarker = (name: string) => {
    const newSel = new Set(selected);
    if (newSel.has(name)) newSel.delete(name); else newSel.add(name);
    setSelected(newSel);
    applySelection(newSel, rateOverrides, copiesOverrides, markerPool);
  };

  const toggleAll = (check: boolean) => {
    const visible = filteredPool.map((m) => m.name);
    const newSel = new Set(selected);
    visible.forEach((n) => (check ? newSel.add(n) : newSel.delete(n)));
    setSelected(newSel);
    applySelection(newSel, rateOverrides, copiesOverrides, markerPool);
  };

  const persistExtras = (
    rates: Map<string, number>,
    cops: Map<string, number>,
    customs: MarkerInfo[]
  ) => {
    persistMarkerPoolExtras({
      customMarkers: customs,
      rateOverrides: Object.fromEntries(rates),
      copiesOverrides: Object.fromEntries(cops),
    });
  };

  const handleRateChange = (name: string, value: string) => {
    const rate = parseFloat(value);
    if (isNaN(rate)) return;
    const newRates = new Map(rateOverrides);
    newRates.set(name, rate);
    setRateOverrides(newRates);
    persistExtras(newRates, copiesOverrides, customMarkers);
    if (selected.has(name)) applySelection(selected, newRates, copiesOverrides, markerPool);
  };

  const handleCopiesChange = (name: string, value: string) => {
    const copies = parseInt(value, 10);
    if (isNaN(copies) || copies < 1) return;
    const newCops = new Map(copiesOverrides);
    newCops.set(name, copies);
    setCopiesOverrides(newCops);
    persistExtras(rateOverrides, newCops, customMarkers);
    if (selected.has(name)) applySelection(selected, rateOverrides, newCops, markerPool);
  };

  const handleResetToFactory = () => {
    if (!window.confirm(t("markers_confirm_reset"))) return;
    setPoolLoading(true);
    invoke<MarkerInfo[]>("list_all_markers")
      .then((pool) => {
        const emptyRates = new Map<string, number>();
        const emptyCops = new Map<string, number>();
        setMarkerPool(pool);
        setCustomMarkers([]);
        setRateOverrides(emptyRates);
        setCopiesOverrides(emptyCops);
        persistMarkerPoolExtras({ customMarkers: [], rateOverrides: {}, copiesOverrides: {} });
        // Remove any selected custom markers that no longer exist in pool
        const poolNames = new Set(pool.map((m) => m.name));
        const newSel = new Set([...selected].filter((n) => poolNames.has(n)));
        setSelected(newSel);
        if (newSel.size > 0) applySelection(newSel, emptyRates, emptyCops, pool);
        setPoolLoading(false);
      })
      .catch((e) => { setError(String(e)); setPoolLoading(false); });
  };

  const handleAddToPool = () => {
    const name = newName.trim();
    const rate = parseFloat(newRate);
    const copies = parseInt(newCopies, 10);
    if (!name || isNaN(rate) || isNaN(copies)) return;
    if (markerPool.some((m) => m.name === name)) {
      setError(t("markers_error_already_in_pool").replace("{name}", name));
      return;
    }
    const newMarker: MarkerInfo = { name, mutationRate: rate, numberOfCopies: copies };
    const newPool = [...markerPool, newMarker].sort((a, b) => a.name.localeCompare(b.name));
    const newCustoms = [...customMarkers, newMarker];
    setMarkerPool(newPool);
    setCustomMarkers(newCustoms);
    // Auto-select the new marker
    const newSel = new Set(selected);
    newSel.add(name);
    const newRates = new Map(rateOverrides);
    newRates.set(name, rate);
    const newCops = new Map(copiesOverrides);
    newCops.set(name, copies);
    setSelected(newSel);
    setRateOverrides(newRates);
    setCopiesOverrides(newCops);
    persistExtras(newRates, newCops, newCustoms);
    applySelection(newSel, newRates, newCops, newPool);
    setNewName("");
    setError(null);
  };

  const handleOpenSaveDialog = () => {
    setSaveSetName("");
    setSaveSetDialogOpen(true);
  };

  const handleSaveSet = () => {
    const name = saveSetName.trim();
    if (!name || selected.size === 0) return;
    const activeMarkers: MarkerInfo[] = [...selected].map((mName) => {
      const base = markerPool.find((m) => m.name === mName);
      const detected = detectedCopies.get(mName);
      const effectiveCopies = typeof detected === "number"
        ? detected
        : copiesOverrides.get(mName) ?? base?.numberOfCopies ?? 1;
      return {
        name: mName,
        mutationRate: rateOverrides.get(mName) ?? base?.mutationRate ?? 0.004,
        numberOfCopies: effectiveCopies,
      };
    });
    const updated = savedSets.filter((s) => s.name !== name);
    updated.push({ name, markers: activeMarkers });
    updated.sort((a, b) => a.name.localeCompare(b.name));
    setSavedSets(updated);
    persistSavedSets(updated);
    setSaveSetDialogOpen(false);
    setSaveSetName("");
  };

  const handleLoadSavedSet = (set: SavedMarkerSet) => {
    setCustomMode(false);
    const csv =
      set.markers.length > 0
        ? "marker,mutation_rate,number_of_copies\n" +
          set.markers.map((m) => `${m.name},${m.mutationRate},${m.numberOfCopies ?? 1}`).join("\n")
        : null;
    setMarkerSet(null, set.markers, csv);
  };

  const handleDeleteSavedSet = (name: string) => {
    const updated = savedSets.filter((s) => s.name !== name);
    setSavedSets(updated);
    persistSavedSets(updated);
  };

  // Derive copy counts from loaded haplotypes: marker → detected count or "conflict"
  const detectedCopies = useMemo<Map<string, number | "conflict">>(() => {
    const result = new Map<string, number | "conflict">();
    if (!haplotypes) return result;
    for (const markerName of haplotypes.markerNames) {
      let detected: number | null = null;
      let conflict = false;
      for (const indName of Object.keys(haplotypes.haplotypeTable)) {
        const allele = haplotypes.haplotypeTable[indName]?.[markerName];
        if (!allele || allele.trim() === "") continue;
        const count = allele.trim().split(";").length;
        if (detected === null) {
          detected = count;
        } else if (detected !== count) {
          conflict = true;
          break;
        }
      }
      if (detected !== null) result.set(markerName, conflict ? "conflict" : detected);
    }
    return result;
  }, [haplotypes]);

  const filteredPool = markerPool.filter((m) =>
    m.name.toLowerCase().includes(search.toLowerCase())
  );
  const allVisibleSelected =
    filteredPool.length > 0 && filteredPool.every((m) => selected.has(m.name));

  const overallRate =
    markers.length > 0
      ? 1 - markers.reduce((acc, m) => acc * (1 - m.mutationRate), 1)
      : null;

  return (
    <div className="p-6 max-w-4xl mx-auto space-y-5">
      <h1 className="text-xl font-bold text-gray-900">{t("markers_title")}</h1>

      {error && (
        <div className="rounded bg-red-50 border border-red-200 p-3 text-sm text-red-800">
          {error}
        </div>
      )}

      {markers.length > 0 && (
        <div className="rounded bg-blue-50 border border-blue-200 p-3 text-sm text-blue-800">
          <strong>{t("markers_active_prefix")}</strong>{" "}
          {selectedKitName ? (
            <>
              {t("markers_builtin_kit")} <em>{selectedKitName}</em>
            </>
          ) : (
            <>{t("markers_custom_set_active")}</>
          )}{" "}
          ({markers.length} markers
          {overallRate !== null && (
            <>
              , {t("markers_overall_rate")}{" "}
              <strong>{overallRate.toExponential(4)}</strong>
            </>
          )}
          )
        </div>
      )}

      {/* Built-in kits */}
      <section className="bg-white rounded-lg border p-4">
        <h2 className="font-semibold text-gray-700 mb-3">{t("markers_embedded")}</h2>
        {kitNames.length === 0 ? (
          <p className="text-sm text-gray-400">{t("markers_loading")}</p>
        ) : (
          <div className="grid grid-cols-2 gap-2">
            {kitNames.map((name) => (
              <button
                key={name}
                onClick={() => handleSelectKit(name)}
                disabled={loading}
                className={`text-left px-3 py-2 rounded border text-sm transition-colors ${
                  selectedKitName === name
                    ? "border-blue-400 bg-blue-50 text-blue-800 font-medium"
                    : "border-gray-200 hover:bg-gray-50 text-gray-700"
                }`}
              >
                {name}
                {selectedKitName === name && (
                  <span className="ml-2 text-xs text-blue-500">{t("markers_active_badge")}</span>
                )}
              </button>
            ))}
          </div>
        )}
      </section>

      {/* Saved custom sets */}
      {savedSets.length > 0 && (
        <section className="bg-white rounded-lg border p-4">
          <h2 className="font-semibold text-gray-700 mb-3">{t("markers_saved_title")}</h2>
          <div className="grid grid-cols-2 gap-2">
            {savedSets.map((set) => {
              const isActive = !selectedKitName && markers.length > 0 &&
                markers.length === set.markers.length &&
                set.markers.every((m) => markers.some((am) => am.name === m.name));
              return (
                <div key={set.name} className={`flex items-center gap-1 px-2 py-1.5 rounded border text-sm ${isActive ? "border-blue-400 bg-blue-50" : "border-gray-200"}`}>
                  <button
                    onClick={() => handleLoadSavedSet(set)}
                    className={`flex-1 text-left truncate ${isActive ? "text-blue-800 font-medium" : "text-gray-700 hover:text-gray-900"}`}
                    title={`${set.markers.length} markers`}
                  >
                    {set.name}
                    {isActive && <span className="ml-2 text-xs text-blue-500">✓ active</span>}
                  </button>
                  <span className="text-xs text-gray-400 whitespace-nowrap">{set.markers.length}m</span>
                  <button
                    onClick={() => handleDeleteSavedSet(set.name)}
                    className="text-gray-300 hover:text-red-500 ml-1 text-xs leading-none"
                    title={t("markers_delete_set_tooltip")}
                  >
                    ✕
                  </button>
                </div>
              );
            })}
          </div>
        </section>
      )}

      {/* Custom set builder */}
      <section className="bg-white rounded-lg border p-4">
        <div className="flex items-center justify-between mb-3">
          <h2 className="font-semibold text-gray-700">
            {t("markers_custom_builder")}
            {customMode && (
              <span className="ml-2 text-xs font-normal text-blue-600">
                — {selected.size} {t("markers_selected_suffix")}
                {!selectedKitName && markers.length > 0 && selected.size > 0 && (
                  <span className="ml-1 text-green-600">{t("markers_active_badge")}</span>
                )}
              </span>
            )}
          </h2>
          {!customMode ? (
            <div className="flex items-center gap-2">
              <button
                onClick={handleEnterCustomMode}
                className="text-sm bg-blue-600 hover:bg-blue-700 text-white px-3 py-1.5 rounded"
              >
                {t("markers_build_custom")}
              </button>
              {(customMarkers.length > 0 || rateOverrides.size > 0 || copiesOverrides.size > 0) && (
                <button
                  onClick={handleResetToFactory}
                  className="text-xs text-red-500 hover:text-red-700 border border-red-200 hover:border-red-400 px-2 py-1.5 rounded"
                  title={t("markers_reset_factory_tooltip")}
                >
                  {t("markers_reset_factory")}
                </button>
              )}
            </div>
          ) : (
            <div className="flex items-center gap-2">
              {selected.size > 0 && (
                <button
                  onClick={handleOpenSaveDialog}
                  className="text-xs bg-emerald-600 hover:bg-emerald-700 text-white px-2 py-1 rounded"
                >
                  {t("markers_save_set_open_btn")}
                </button>
              )}
              <button
                onClick={handleResetToFactory}
                className="text-xs text-red-500 hover:text-red-700 border border-red-200 hover:border-red-400 px-2 py-1 rounded"
                title={t("markers_reset_factory_tooltip")}
              >
                {t("markers_reset_factory")}
              </button>
              <button
                onClick={() => setCustomMode(false)}
                className="text-xs text-gray-500 hover:text-gray-700"
              >
                {t("markers_close_builder")}
              </button>
            </div>
          )}
        </div>

        {/* Save-set name dialog */}
        {saveSetDialogOpen && (
          <div className="mb-3 p-3 bg-gray-50 rounded border border-gray-200 flex gap-2 items-end">
            <div className="flex-1">
              <label className="block text-xs text-gray-600 mb-1">{t("markers_save_set_name_label")}</label>
              <input
                autoFocus
                type="text"
                className="w-full border rounded px-2 py-1 text-sm"
                placeholder={t("markers_save_set_example")}
                value={saveSetName}
                onChange={(e) => setSaveSetName(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") handleSaveSet();
                  if (e.key === "Escape") setSaveSetDialogOpen(false);
                }}
              />
            </div>
            <button
              onClick={handleSaveSet}
              disabled={!saveSetName.trim()}
              className="bg-emerald-600 hover:bg-emerald-700 disabled:opacity-40 text-white text-sm px-3 py-1.5 rounded"
            >
              {t("markers_save_btn")}
            </button>
            <button
              onClick={() => setSaveSetDialogOpen(false)}
              className="bg-white border text-gray-700 text-sm px-3 py-1.5 rounded hover:bg-gray-50"
            >
              {t("ped_cancel")}
            </button>
          </div>
        )}

        {!customMode && (
          <p className="text-sm text-gray-400">{t("markers_custom_desc")}</p>
        )}

        {customMode && (
          <>
            {/* Search + select-all */}
            <div className="flex gap-2 mb-2 items-center">
              <input
                type="text"
                placeholder={t("markers_search_placeholder")}
                className="flex-1 border rounded px-3 py-1.5 text-sm"
                value={search}
                onChange={(e) => setSearch(e.target.value)}
              />
              <button
                onClick={() => toggleAll(!allVisibleSelected)}
                className="text-xs border rounded px-3 py-1.5 bg-gray-50 hover:bg-gray-100 text-gray-700 whitespace-nowrap"
              >
                {allVisibleSelected ? t("markers_deselect_all") : t("markers_select_all")}
              </button>
            </div>

            {/* Copy-count conflict warning */}
            {[...detectedCopies.entries()].some(([, v]) => v === "conflict") && (
              <div className="mb-2 rounded bg-orange-50 border border-orange-200 p-2 text-xs text-orange-800">
                {t("markers_copy_conflict")}
              </div>
            )}

            {/* Pool table */}
            {poolLoading ? (
              <p className="text-sm text-gray-400 py-4">{t("markers_pool_loading")}</p>
            ) : (
              <div className="overflow-auto max-h-72 border rounded">
                <table className="text-xs w-full border-collapse">
                  <thead className="sticky top-0 bg-gray-50">
                    <tr>
                      <th className="border px-2 py-2 w-8" />
                      <th className="border px-3 py-2 text-left">{t("markers_marker")}</th>
                      <th className="border px-3 py-2 text-right">{t("markers_mutation_rate")}</th>
                      <th className="border px-3 py-2 text-right">{t("markers_copies")}</th>
                    </tr>
                  </thead>
                  <tbody>
                    {filteredPool.map((m, i) => {
                      const isChecked = selected.has(m.name);
                      const rate = rateOverrides.get(m.name) ?? m.mutationRate;
                      const detected = detectedCopies.get(m.name);
                      const copies = typeof detected === "number"
                        ? detected
                        : copiesOverrides.get(m.name) ?? m.numberOfCopies ?? 1;
                      const hasConflict = detected === "conflict";
                      return (
                        <tr
                          key={m.name}
                          className={
                            isChecked
                              ? "bg-blue-50"
                              : i % 2 === 0
                              ? "bg-white"
                              : "bg-gray-50"
                          }
                        >
                          <td className="border px-2 py-1 text-center">
                            <input
                              type="checkbox"
                              checked={isChecked}
                              onChange={() => toggleMarker(m.name)}
                            />
                          </td>
                          <td className="border px-3 py-1 font-mono">{m.name}</td>
                          <td className="border px-3 py-1 text-right">
                            <div className="flex items-center justify-end gap-1">
                              {rateOverrides.has(m.name) && rateOverrides.get(m.name) !== m.mutationRate && (
                                <span className="text-orange-400 text-xs" title={`Modified from factory default (${m.mutationRate})`}>●</span>
                              )}
                              <input
                                type="number"
                                step="0.0001"
                                min="0"
                                max="1"
                                className={`w-28 border rounded px-1 py-0.5 text-right font-mono text-xs ${rateOverrides.has(m.name) && rateOverrides.get(m.name) !== m.mutationRate ? "border-orange-300 text-orange-700" : ""}`}
                                value={rate}
                                onChange={(e) => handleRateChange(m.name, e.target.value)}
                              />
                            </div>
                          </td>
                          <td className="border px-3 py-1 text-right">
                            {typeof detected === "number" ? (
                              <span className="inline-flex items-center gap-1">
                                <span className="font-mono text-xs text-gray-700">{detected}</span>
                                <span className="text-xs text-gray-400">{t("markers_auto_copies")}</span>
                              </span>
                            ) : hasConflict ? (
                              <span className="inline-flex items-center gap-1">
                                <input
                                  type="number"
                                  min="1"
                                  max="4"
                                  className="w-14 border border-orange-400 rounded px-1 py-0.5 text-right text-xs"
                                  value={copiesOverrides.get(m.name) ?? m.numberOfCopies ?? 1}
                                  onChange={(e) => handleCopiesChange(m.name, e.target.value)}
                                />
                                <span title="Individuals have different copy counts for this marker" className="text-orange-500 cursor-help">⚠</span>
                              </span>
                            ) : (
                              <input
                                type="number"
                                min="1"
                                max="4"
                                className="w-14 border rounded px-1 py-0.5 text-right text-xs"
                                value={copies}
                                onChange={(e) => handleCopiesChange(m.name, e.target.value)}
                              />
                            )}
                          </td>
                        </tr>
                      );
                    })}
                    {filteredPool.length === 0 && (
                      <tr>
                        <td colSpan={4} className="text-center text-gray-400 py-4">
                          {t("markers_no_match").replace("{search}", search)}
                        </td>
                      </tr>
                    )}
                  </tbody>
                </table>
              </div>
            )}

            {/* Add new marker to pool */}
            <div className="mt-3 pt-3 border-t">
              <p className="text-xs font-medium text-gray-600 mb-2">{t("markers_add_marker")}:</p>
              <div className="flex gap-2 items-end flex-wrap">
                <div>
                  <label className="block text-xs text-gray-500 mb-0.5">{t("markers_marker_name")}</label>
                  <input
                    type="text"
                    className="border rounded px-2 py-1 text-xs w-32 font-mono"
                    placeholder="DYS390"
                    value={newName}
                    onChange={(e) => setNewName(e.target.value)}
                    onKeyDown={(e) => e.key === "Enter" && handleAddToPool()}
                  />
                </div>
                <div>
                  <label className="block text-xs text-gray-500 mb-0.5">{t("markers_rate")}</label>
                  <input
                    type="number"
                    step="0.0001"
                    min="0"
                    max="1"
                    className="border rounded px-2 py-1 text-xs w-24"
                    value={newRate}
                    onChange={(e) => setNewRate(e.target.value)}
                  />
                </div>
                <div>
                  <label className="block text-xs text-gray-500 mb-0.5">{t("markers_copies")}</label>
                  <input
                    type="number"
                    min="1"
                    max="4"
                    className="border rounded px-2 py-1 text-xs w-16"
                    value={newCopies}
                    onChange={(e) => setNewCopies(e.target.value)}
                  />
                </div>
                <button
                  onClick={handleAddToPool}
                  className="bg-blue-600 hover:bg-blue-700 text-white text-xs font-medium px-3 py-1.5 rounded"
                >
                  {t("markers_add_marker")}
                </button>
              </div>
            </div>
          </>
        )}
      </section>
    </div>
  );
}
