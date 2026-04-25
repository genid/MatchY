import { useEffect, useRef, useState } from "react";
import { open, save } from "@tauri-apps/plugin-dialog";
import { readTextFile, writeTextFile } from "@tauri-apps/plugin-fs";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

function isValidAllele(val: string): boolean {
  if (!val.trim()) return false;
  return val.split(";").every((part) =>
    part.split(".").every((n) => /^\d+$/.test(n.trim()))
  );
}

interface ValidationError {
  message: string;
  marker?: string;
  individual?: string;
}

function validate(
  table: Record<string, Record<string, string>>,
  trace: Record<string, string> | null,
  markerNames: string[],
  columns: string[]
): ValidationError[] {
  const errors: ValidationError[] = [];
  for (const marker of markerNames) {
    const copyNums = new Set<number>();
    for (const ind of columns) {
      const val = table[ind]?.[marker] ?? "";
      if (!val) {
        errors.push({ message: `Empty cell at ${ind} / ${marker}`, marker, individual: ind });
        continue;
      }
      if (!isValidAllele(val)) {
        errors.push({ message: `Invalid allele '${val}' at ${ind} / ${marker}`, marker, individual: ind });
      }
      copyNums.add(val.split(";").length);
    }
    if (trace) {
      const val = trace[marker] ?? "";
      if (val && isValidAllele(val)) copyNums.add(val.split(";").length);
    }
    if (copyNums.size > 1) {
      errors.push({ message: `Inconsistent copy number at marker ${marker}`, marker });
    }
  }
  return errors;
}

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------

export default function HaplotypeEditor() {
  const { haplotypes, pedigree, markers, selectedKitName, setHaplotypes, suspect, setSuspect, exclude, setExclude } =
    useAppStore();

  const [table, setTable] = useState<Record<string, Record<string, string>>>({});
  const [trace, setTrace] = useState<Record<string, string> | null>(null);
  const [columns, setColumns] = useState<string[]>([]);

  const [localMarkerNames, setLocalMarkerNames] = useState<string[]>([]);
  const [error, setError] = useState<string | null>(null);
  const [feedback, setFeedback] = useState<string | null>(null);
  const [importing, setImporting] = useState(false);

  // Add-individual panel state
  const [addName, setAddName] = useState("");
  const [addFromPedigree, setAddFromPedigree] = useState("");

  // Paste from spreadsheet
  const [pasteOpen, setPasteOpen] = useState(false);
  const [pasteText, setPasteText] = useState("");
  const [pasteError, setPasteError] = useState<string | null>(null);

  // Marker names: prefer loaded kit, then local (post-import), then store
  const kitMarkerNames = markers.map((m) => m.name);
  const markerNames =
    kitMarkerNames.length > 0
      ? kitMarkerNames
      : localMarkerNames.length > 0
      ? localMarkerNames
      : haplotypes?.markerNames ?? [];

  // Initialise local state once on mount from store (empty deps avoids circular loop
  // when auto-save writes back to the store and would re-trigger initialisation).
  // Re-initialises on unmount/remount (e.g. navigation away and back).
  // eslint-disable-next-line react-hooks/exhaustive-deps
  useEffect(() => {
    if (haplotypes) {
      setTable(structuredClone(haplotypes.haplotypeTable));
      setTrace(haplotypes.traceHaplotype ? structuredClone(haplotypes.traceHaplotype) : null);
      setColumns(Object.keys(haplotypes.haplotypeTable));
      setLocalMarkerNames(haplotypes.markerNames);
    }
  }, []);

  // Auto-save to simulation (debounced) whenever local state changes.
  useEffect(() => {
    if (columns.length === 0) return;
    const errs = validate(table, trace, markerNames, columns);
    if (errs.length > 0) return;
    const timer = setTimeout(async () => {
      try {
        const json = await invoke<string>("export_haplotypes_json", {
          haplotypeTable: table,
          traceHaplotype: trace,
        });
        setHaplotypes({ haplotypeTable: table, traceHaplotype: trace, markerNames }, json);
      } catch {}
    }, 300);
    return () => clearTimeout(timer);
  }, [table, trace, columns]); // markerNames derives from state already tracked

  const emptyRow = (): Record<string, string> =>
    Object.fromEntries(markerNames.map((m) => [m, ""]));

  const copyLastRow = (): Record<string, string> =>
    columns.length > 0
      ? structuredClone(table[columns[columns.length - 1]] ?? emptyRow())
      : emptyRow();

  const showFeedback = (msg: string) => {
    setFeedback(msg);
    setTimeout(() => setFeedback(null), 2500);
  };

  // ------------------------------------------------------------------
  // Cell editing
  // ------------------------------------------------------------------

  const handleCellChange = (ind: string, marker: string, val: string) =>
    setTable((prev) => ({ ...prev, [ind]: { ...prev[ind], [marker]: val } }));

  const handleTraceChange = (marker: string, val: string) =>
    setTrace((prev) => ({ ...prev!, [marker]: val }));

  // Keyboard: Tab → next cell, Enter → cell below
  const cellRef = useRef<Map<string, HTMLInputElement>>(new Map());
  const cellKey = (col: string, row: number) => `${col}::${row}`;

  const handleKeyDown = (
    e: React.KeyboardEvent,
    colIdx: number,
    rowIdx: number,
    allCols: string[]
  ) => {
    let nextCol = colIdx;
    let nextRow = rowIdx;
    if (e.key === "Tab") {
      e.preventDefault();
      nextCol = e.shiftKey ? colIdx - 1 : colIdx + 1;
      if (nextCol < 0) { nextCol = allCols.length - 1; nextRow--; }
      if (nextCol >= allCols.length) { nextCol = 0; nextRow++; }
    } else if (e.key === "Enter") {
      e.preventDefault();
      nextRow = e.shiftKey ? rowIdx - 1 : rowIdx + 1;
    } else {
      return;
    }
    nextRow = Math.max(0, Math.min(markerNames.length - 1, nextRow));
    nextCol = Math.max(0, Math.min(allCols.length - 1, nextCol));
    cellRef.current.get(cellKey(allCols[nextCol], nextRow))?.focus();
  };

  // ------------------------------------------------------------------
  // Add / Remove
  // ------------------------------------------------------------------

  const doAddIndividual = (name: string, copyFrom?: string) => {
    name = name.trim();
    if (!name) { setError("Enter a name"); return; }
    if (columns.includes(name)) { setError(`'${name}' already in table`); return; }
    if (name.toUpperCase() === "TRACE") { setError("Use 'Add TRACE' for TRACE profile"); return; }
    setError(null);
    const values = copyFrom && table[copyFrom] ? structuredClone(table[copyFrom]) : copyLastRow();
    setColumns((prev) => [...prev, name]);
    setTable((prev) => ({ ...prev, [name]: values }));
    setAddName("");
    setAddFromPedigree("");
    showFeedback(`Added: ${name}`);
  };

  const handleRemoveIndividual = (name: string) => {
    setColumns((prev) => prev.filter((c) => c !== name));
    setTable((prev) => {
      const next = { ...prev };
      delete next[name];
      return next;
    });
    showFeedback(`Removed: ${name}`);
  };

  const handleImportCsvForIndividual = async (indName: string) => {
    setError(null);
    try {
      const filePath = await open({
        multiple: false,
        filters: [{ name: "CSV", extensions: ["csv", "txt"] }],
      });
      if (!filePath || typeof filePath !== "string") return;
      const content = await readTextFile(filePath);

      const lines = content.split(/\r?\n/).map((l) => l.trim()).filter(Boolean);
      if (lines.length < 2) { setError("CSV has no data rows."); return; }

      // Detect delimiter (comma or semicolon between first two fields)
      const header = lines[0];
      const delim = header.includes(";") && !header.includes(",") ? ";" : ",";
      const headerCols = header.split(delim).map((h) => h.trim().toLowerCase());
      const markerIdx = headerCols.findIndex((h) => h === "marker_name" || h === "marker");
      const allelesIdx = headerCols.findIndex((h) => h === "alleles" || h === "allele");
      if (markerIdx === -1 || allelesIdx === -1) {
        setError("CSV must have columns 'marker_name' and 'alleles'.");
        return;
      }

      const updates: Record<string, string> = {};
      const unknown: string[] = [];
      for (const line of lines.slice(1)) {
        const parts = line.split(delim);
        const marker = parts[markerIdx]?.trim();
        const alleles = parts[allelesIdx]?.trim() ?? "";
        if (!marker) continue;
        if (markerNames.includes(marker)) {
          updates[marker] = alleles;
        } else {
          unknown.push(marker);
        }
      }

      setTable((prev) => ({
        ...prev,
        [indName]: { ...prev[indName], ...updates },
      }));

      const msg = `Imported ${Object.keys(updates).length} markers for ${indName}`;
      showFeedback(unknown.length > 0 ? `${msg} (${unknown.length} unrecognised markers skipped)` : msg);
    } catch (e) {
      setError(String(e));
    }
  };

  const handleAddTrace = () => {
    setTrace(copyLastRow());
    showFeedback("TRACE profile added");
  };

  const handleRemoveTrace = () => {
    setTrace(null);
    showFeedback("TRACE profile removed");
  };

  // ------------------------------------------------------------------
  // Import / Export
  // ------------------------------------------------------------------

  const handleImport = async () => {
    setError(null);
    setImporting(true);
    try {
      const filePath = await open({
        multiple: false,
        filters: [{ name: "Haplotypes JSON", extensions: ["json"] }],
      });
      if (!filePath || typeof filePath !== "string") return;
      const content = await readTextFile(filePath);
      const parsed = JSON.parse(content) as Record<string, Record<string, string>>;

      const newTrace = parsed["TRACE"] ?? null;
      const newTable = Object.fromEntries(
        Object.entries(parsed).filter(([k]) => k !== "TRACE")
      );
      const newCols = Object.keys(newTable);
      const names =
        kitMarkerNames.length > 0
          ? kitMarkerNames
          : Object.keys(Object.values(newTable)[0] ?? {});

      // Ensure each individual has all markers
      for (const ind of newCols) {
        for (const m of names) {
          if (!newTable[ind][m]) newTable[ind][m] = "";
        }
      }

      setTable(newTable);
      setTrace(newTrace);
      setColumns(newCols);
      setLocalMarkerNames(names);
      showFeedback(`Imported ${newCols.length} individuals, ${names.length} markers`);
    } catch (e) {
      setError(String(e));
    } finally {
      setImporting(false);
    }
  };

  const handleExport = async () => {
    setError(null);
    try {
      const json = await invoke<string>("export_haplotypes_json", {
        haplotypeTable: table,
        traceHaplotype: trace,
      });
      const filePath = await save({
        filters: [{ name: "Haplotypes JSON", extensions: ["json"] }],
        defaultPath: "haplotypes.json",
      });
      if (!filePath) return;
      await writeTextFile(filePath, json);
      showFeedback("Exported");
    } catch (e) {
      setError(String(e));
    }
  };

  // ------------------------------------------------------------------
  // Paste from spreadsheet (Excel/Sheets tab-separated format)
  // Expected: first row = header (blank cell + individual names)
  //           first column = marker names, rest = allele values
  // ------------------------------------------------------------------

  const handlePasteImport = () => {
    setPasteError(null);
    const lines = pasteText.split(/\r?\n/).map((l) => l.trimEnd()).filter((l) => l.trim());
    if (lines.length < 2) { setPasteError("Need at least a header row and one data row."); return; }

    const header = lines[0].split("\t");
    // First cell may be "Marker", "marker_name", blank, etc. — skip it
    const indNames = header.slice(1).map((s) => s.trim()).filter(Boolean);
    if (indNames.length === 0) { setPasteError("No individual names found in the header row."); return; }

    const newTable: Record<string, Record<string, string>> = {};
    const pastedMarkerNames: string[] = [];

    for (const line of lines.slice(1)) {
      const cells = line.split("\t");
      const markerName = cells[0]?.trim();
      if (!markerName) continue;
      pastedMarkerNames.push(markerName);
      for (let i = 0; i < indNames.length; i++) {
        const ind = indNames[i];
        if (!newTable[ind]) newTable[ind] = {};
        newTable[ind][markerName] = cells[i + 1]?.trim() ?? "";
      }
    }

    if (pastedMarkerNames.length === 0) { setPasteError("No marker rows found."); return; }

    // Fill in missing markers for each individual
    const effectiveMarkers = kitMarkerNames.length > 0 ? kitMarkerNames : pastedMarkerNames;
    for (const ind of indNames) {
      for (const m of effectiveMarkers) {
        if (!newTable[ind][m]) newTable[ind][m] = "";
      }
    }

    // Merge with existing table: new individuals are added, existing ones updated
    const mergedTable = { ...table };
    const mergedCols = [...columns];
    for (const ind of indNames) {
      if (!mergedCols.includes(ind)) mergedCols.push(ind);
      mergedTable[ind] = { ...(mergedTable[ind] ?? {}), ...newTable[ind] };
    }

    setTable(mergedTable);
    setColumns(mergedCols);
    if (kitMarkerNames.length === 0) setLocalMarkerNames(pastedMarkerNames);
    setPasteOpen(false);
    setPasteText("");
    showFeedback(`Pasted ${indNames.length} individuals, ${pastedMarkerNames.length} markers`);
  };

  // ------------------------------------------------------------------
  // Derived state
  // ------------------------------------------------------------------

  const validationErrors = validate(table, trace, markerNames, columns);
  const hasErrors = validationErrors.length > 0;
  const errorSet = new Set(validationErrors.map((e) => `${e.individual}::${e.marker}`));
  const markerErrors = new Set(validationErrors.map((e) => e.marker).filter(Boolean));

  // Individuals in pedigree not yet in table
  const pedigreeNames = pedigree?.individuals.map((i) => i.name) ?? [];
  const availableFromPedigree = pedigreeNames.filter(
    (n) => !columns.includes(n) && n.toUpperCase() !== "TRACE"
  );

  // All display columns: TRACE first (if present), then individuals
  const allCols = trace ? ["__TRACE__", ...columns] : columns;

  // No marker names at all → ask to select a kit
  if (markerNames.length === 0 && columns.length === 0) {
    return (
      <div className="h-full flex flex-col items-center justify-center gap-4 text-gray-400 p-8">
        <p className="text-sm font-medium">No marker set selected</p>
        <p className="text-xs text-center max-w-xs">
          Go to <strong className="text-gray-600">Marker Sets</strong> to select a kit or upload
          a custom CSV, or import a haplotypes JSON directly.
        </p>
        <button
          onClick={handleImport}
          disabled={importing}
          className="text-sm bg-blue-600 text-white rounded px-4 py-2 hover:bg-blue-700 disabled:opacity-50"
        >
          Import haplotypes JSON
        </button>
      </div>
    );
  }

  return (
    <div className="h-full flex flex-col overflow-hidden">
      {/* ── Toolbar ── */}
      <div className="bg-white border-b px-4 py-2 flex items-center gap-2 flex-wrap flex-shrink-0">
        <button
          onClick={handleImport}
          disabled={importing}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50 disabled:opacity-50"
        >
          {importing ? "Importing…" : "Import JSON"}
        </button>
        <button
          onClick={handleExport}
          disabled={columns.length === 0}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50 disabled:opacity-50"
        >
          Export JSON
        </button>
        <button
          onClick={() => { setPasteOpen(true); setPasteText(""); setPasteError(null); }}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50"
          title="Paste tab-separated data copied from Excel or Google Sheets"
        >
          Paste from Excel…
        </button>
        <span className="text-xs text-gray-400 ml-1">
          {markerNames.length} markers · {columns.length} individual{columns.length !== 1 ? "s" : ""}
          {trace ? " · TRACE" : ""}
          {selectedKitName ? ` · kit: ${selectedKitName}` : ""}
          {columns.length > 0 && !hasErrors && (
            <span className="ml-2 text-emerald-600">● auto-saved</span>
          )}
        </span>

        {feedback && (
          <span className="ml-auto text-xs text-green-600 font-medium">{feedback}</span>
        )}
      </div>

      {pasteOpen && (
        <div className="fixed inset-0 bg-black/30 z-50 flex items-center justify-center p-6">
          <div className="bg-white rounded-xl shadow-2xl border w-[560px] max-h-[80vh] flex flex-col">
            <div className="px-5 pt-4 pb-3 border-b">
              <p className="font-semibold text-gray-800">Paste from Excel / Google Sheets</p>
              <p className="text-xs text-gray-500 mt-1">
                Copy a range from your spreadsheet and paste it below. Expected format: first row = individual names (skip first cell), first column = marker names.
              </p>
            </div>
            <textarea
              autoFocus
              className="flex-1 font-mono text-xs p-3 resize-none focus:outline-none border-b min-h-[200px]"
              placeholder={"(blank)\tFather\tSon\tUncle\nDYS19\t15\t15\t16\nDYS389I\t13\t13\t14"}
              value={pasteText}
              onChange={(e) => setPasteText(e.target.value)}
            />
            {pasteError && (
              <p className="px-4 py-2 text-xs text-red-700 bg-red-50 border-b">{pasteError}</p>
            )}
            <div className="px-4 py-3 flex gap-2 justify-end">
              <button
                onClick={() => { setPasteOpen(false); setPasteText(""); setPasteError(null); }}
                className="text-sm border rounded px-4 py-1.5 hover:bg-gray-50 text-gray-700"
              >
                Cancel
              </button>
              <button
                onClick={handlePasteImport}
                disabled={!pasteText.trim()}
                className="text-sm bg-blue-600 hover:bg-blue-700 disabled:opacity-40 text-white rounded px-4 py-1.5"
              >
                Import
              </button>
            </div>
          </div>
        </div>
      )}

      {error && (
        <div className="bg-red-50 border-b border-red-200 px-4 py-2 text-xs text-red-700 flex items-center gap-2 flex-shrink-0">
          {error}
          <button className="underline ml-1" onClick={() => setError(null)}>dismiss</button>
        </div>
      )}

      {hasErrors && (
        <div className="bg-amber-50 border-b border-amber-200 px-4 py-1.5 text-xs text-amber-700 flex-shrink-0">
          {validationErrors.slice(0, 3).map((e, i) => (
            <div key={i}>⚠ {e.message}</div>
          ))}
          {validationErrors.length > 3 && (
            <div>…and {validationErrors.length - 3} more error(s)</div>
          )}
        </div>
      )}

      {/* ── Body: left panel + table ── */}
      <div className="flex-1 flex flex-row overflow-hidden">
        {/* Left panel */}
        <div className="w-56 flex-shrink-0 bg-white border-r flex flex-col overflow-y-auto">
          {/* Add Individual */}
          <div className="p-3 border-b">
            <p className="text-xs font-semibold text-gray-600 mb-2 uppercase tracking-wide">
              👤 Add Individual
            </p>

            {availableFromPedigree.length > 0 ? (
              <>
                <label className="text-xs text-gray-500 block mb-1">From pedigree</label>
                <select
                  value={addFromPedigree}
                  onChange={(e) => setAddFromPedigree(e.target.value)}
                  className="w-full text-xs border rounded px-2 py-1.5 mb-2 focus:outline-none focus:ring-1 focus:ring-blue-400"
                >
                  <option value="">— select —</option>
                  {availableFromPedigree.map((n) => (
                    <option key={n} value={n}>{n}</option>
                  ))}
                </select>
                <button
                  onClick={() => doAddIndividual(addFromPedigree)}
                  disabled={!addFromPedigree}
                  className="w-full text-xs bg-blue-600 text-white rounded px-2 py-1.5 hover:bg-blue-700 disabled:opacity-40 mb-2"
                >
                  Add from pedigree
                </button>
                <p className="text-xs text-gray-400 mb-1">Or add manually:</p>
              </>
            ) : (
              <p className="text-xs text-gray-400 mb-1">Enter name:</p>
            )}

            <input
              type="text"
              value={addName}
              onChange={(e) => setAddName(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && doAddIndividual(addName)}
              placeholder="e.g. Father"
              className="w-full text-xs border rounded px-2 py-1.5 mb-2 focus:outline-none focus:ring-1 focus:ring-blue-400"
            />
            <button
              onClick={() => doAddIndividual(addName)}
              disabled={!addName.trim()}
              className="w-full text-xs bg-blue-600 text-white rounded px-2 py-1.5 hover:bg-blue-700 disabled:opacity-40"
            >
              ➕ Add
            </button>
          </div>

          {/* TRACE */}
          <div className="p-3 border-b">
            <p className="text-xs font-semibold text-gray-600 mb-2 uppercase tracking-wide">
              🔬 TRACE Profile
            </p>
            {trace ? (
              <button
                onClick={handleRemoveTrace}
                className="w-full text-xs bg-red-50 text-red-700 border border-red-200 rounded px-2 py-1.5 hover:bg-red-100"
              >
                🗑️ Remove TRACE
              </button>
            ) : (
              <button
                onClick={handleAddTrace}
                disabled={columns.length === 0}
                className="w-full text-xs bg-purple-600 text-white rounded px-2 py-1.5 hover:bg-purple-700 disabled:opacity-40"
              >
                Add TRACE Profile
              </button>
            )}
            <p className="text-xs text-gray-400 mt-1.5">
              {trace
                ? "TRACE present — shown as first column"
                : "TRACE is used for donor identification (trace mode)"}
            </p>
          </div>

          {/* Individuals list: import CSV + remove */}
          {columns.length > 0 && (
            <div className="p-3 border-b">
              <p className="text-xs font-semibold text-gray-600 mb-2 uppercase tracking-wide">
                Individuals
              </p>
              <div className="flex flex-col gap-1">
                {columns.map((name) => (
                  <div key={name} className="flex items-center justify-between gap-1">
                    <span className="text-xs text-gray-700 truncate flex-1">{name}</span>
                    <button
                      onClick={() => handleImportCsvForIndividual(name)}
                      className="text-xs text-blue-500 hover:text-blue-700 px-1 flex-shrink-0"
                      title={`Import CSV for ${name}`}
                    >
                      📥
                    </button>
                    <button
                      onClick={() => handleRemoveIndividual(name)}
                      className="text-xs text-red-500 hover:text-red-700 px-1 flex-shrink-0"
                      title={`Remove ${name}`}
                    >
                      ✕
                    </button>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Format guide */}
          <div className="p-3 mt-auto">
            <p className="text-xs font-semibold text-gray-500 mb-1">Allele formats</p>
            <p className="text-xs text-gray-400 leading-relaxed">
              Single: <code>15</code><br />
              Decimal: <code>15.2</code><br />
              Multi-copy: <code>15;16</code><br />
              Decimal multi: <code>15.2;16</code>
            </p>
          </div>
        </div>

        {/* ── Table ── */}
        {columns.length === 0 && !trace ? (
          <div className="flex-1 flex flex-col items-center justify-center gap-3 text-gray-400">
            <p className="text-sm">No individuals in table</p>
            <p className="text-xs text-center max-w-xs">Add individuals from the left panel, import a JSON file, or use <strong className="text-gray-500">Paste from Excel…</strong> in the toolbar.</p>
          </div>
        ) : (
          <div className="flex-1 overflow-auto p-3">
            <table className="text-xs border-collapse">
              <thead>
                <tr className="sticky top-0 z-10 bg-gray-50">
                  <th className="sticky left-0 z-20 bg-gray-50 border px-3 py-2 text-left font-semibold text-gray-600 whitespace-nowrap min-w-[120px]">
                    Marker
                  </th>
                  {trace && (
                    <th className="border px-2 py-2 text-center font-semibold text-purple-700 bg-purple-50 whitespace-nowrap min-w-[100px]">
                      TRACE
                    </th>
                  )}
                  {columns.map((name) => {
                    const isExcluded = exclude.includes(name);
                    const isSuspect = suspect === name;
                    const inPedigree = !pedigree || pedigree.individuals.some((i) => i.name === name);
                    const cls = isExcluded ? "excluded" : isSuspect ? "suspect" : "known";
                    const colClass = isExcluded
                      ? "text-gray-400 bg-gray-50"
                      : isSuspect
                      ? "text-blue-800 bg-blue-50"
                      : "text-green-800 bg-green-50";
                    return (
                      <th
                        key={name}
                        className={`border px-2 py-1.5 text-center font-medium whitespace-nowrap min-w-[100px] ${colClass}`}
                      >
                        <div className="flex items-center justify-center gap-1">
                          {!inPedigree && (
                            <span title="Individual not found in loaded pedigree" className="text-orange-400 cursor-help">⚠</span>
                          )}
                          {name}
                        </div>
                        <div className="font-normal text-gray-400 text-xs">{cls}</div>
                        <button
                          onClick={() => {
                            if (isExcluded) {
                              setExclude(exclude.filter((n) => n !== name));
                            } else {
                              setExclude([...exclude, name]);
                            }
                          }}
                          className={`mt-0.5 text-xs rounded px-1.5 py-0 border transition-colors ${
                            isExcluded
                              ? "bg-gray-200 text-gray-500 border-gray-300 hover:bg-red-50 hover:text-red-600"
                              : "bg-white text-gray-400 border-gray-200 hover:bg-red-50 hover:text-red-600 hover:border-red-200"
                          }`}
                          title={isExcluded ? "Click to include in calculation" : "Click to exclude from calculation"}
                        >
                          {isExcluded ? "excluded" : "exclude"}
                        </button>
                        <button
                          onClick={() => setSuspect(isSuspect ? null : name)}
                          className={`mt-0.5 text-xs rounded px-1.5 py-0 border transition-colors ${
                            isSuspect
                              ? "bg-red-100 text-red-700 border-red-300 hover:bg-red-50"
                              : "bg-white text-gray-400 border-gray-200 hover:bg-red-50 hover:text-red-700 hover:border-red-200"
                          }`}
                          title={isSuspect ? "Clear as suspect" : "Set as suspect"}
                        >
                          {isSuspect ? "⚑ suspect" : "suspect"}
                        </button>
                      </th>
                    );
                  })}
                </tr>
              </thead>
              <tbody>
                {markerNames.map((marker, rowIdx) => {
                  const knownAlleles = columns
                    .map((n) => table[n]?.[marker])
                    .filter(Boolean);
                  const unique = new Set(knownAlleles);
                  const hasDiversity = unique.size > 1;
                  const rowBg = hasDiversity
                    ? "bg-red-50"
                    : rowIdx % 2 === 0
                    ? "bg-white"
                    : "bg-gray-50";

                  return (
                    <tr key={marker} className={rowBg}>
                      <td className="sticky left-0 z-10 bg-inherit border px-3 py-1 font-mono text-gray-800 whitespace-nowrap">
                        {markerErrors.has(marker) && (
                          <span className="mr-1 text-red-500" title="Error on this marker">●</span>
                        )}
                        {hasDiversity && !markerErrors.has(marker) && (
                          <span className="mr-1 text-amber-400" title="Allelic diversity">◆</span>
                        )}
                        {marker}
                      </td>

                      {trace && (
                        <td className="border p-0 bg-purple-50">
                          <input
                            ref={(el) => {
                              if (el) cellRef.current.set(cellKey("__TRACE__", rowIdx), el);
                            }}
                            type="text"
                            value={trace[marker] ?? ""}
                            onChange={(e) => handleTraceChange(marker, e.target.value)}
                            onKeyDown={(e) => handleKeyDown(e, 0, rowIdx, allCols)}
                            className={`w-full px-2 py-1 font-mono bg-transparent text-center text-purple-800 focus:outline-none focus:bg-purple-100 ${
                              trace[marker] && !isValidAllele(trace[marker])
                                ? "text-red-600 bg-red-50"
                                : ""
                            }`}
                          />
                        </td>
                      )}

                      {columns.map((ind, colIdx) => {
                        const val = table[ind]?.[marker] ?? "";
                        const isErr = errorSet.has(`${ind}::${marker}`);
                        const traceOffset = trace ? 1 : 0;
                        return (
                          <td key={ind} className="border p-0">
                            <input
                              ref={(el) => {
                                if (el) cellRef.current.set(cellKey(ind, rowIdx), el);
                              }}
                              type="text"
                              value={val}
                              onChange={(e) => handleCellChange(ind, marker, e.target.value)}
                              onKeyDown={(e) =>
                                handleKeyDown(e, colIdx + traceOffset, rowIdx, allCols)
                              }
                              className={`w-full px-2 py-1 font-mono bg-transparent text-center focus:outline-none focus:bg-yellow-50 ${
                                isErr ? "text-red-600 bg-red-50" : "text-gray-800"
                              }`}
                            />
                          </td>
                        );
                      })}
                    </tr>
                  );
                })}
              </tbody>
            </table>
            <p className="text-xs text-gray-400 mt-2">
              Tab / Shift+Tab to move between cells · Enter / Shift+Enter to move up/down ·
              ◆ = allelic diversity · ● = validation error
            </p>
          </div>
        )}
      </div>
    </div>
  );
}
