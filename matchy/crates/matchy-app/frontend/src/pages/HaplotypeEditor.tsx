import { useState } from "react";
import { open, save } from "@tauri-apps/plugin-dialog";
import { readTextFile, writeTextFile } from "@tauri-apps/plugin-fs";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import type { HaplotypesParseResult } from "../types/matchy";

export default function HaplotypeEditor() {
  const {
    haplotypes, haplotypesJson,
    pedigree, pedigreeTgf,
    selectedKitName, markerSetCsv,
    setHaplotypes,
  } = useAppStore();
  const [error, setError] = useState<string | null>(null);
  const [importing, setImporting] = useState(false);

  // ------------------------------------------------------------------
  // Import haplotypes from JSON file
  // ------------------------------------------------------------------
  const handleImport = async () => {
    if (!pedigree) {
      setError("Load a pedigree first before importing haplotypes.");
      return;
    }
    if (!selectedKitName && !markerSetCsv) {
      setError("Select a marker set (Marker Sets page) before importing haplotypes.");
      return;
    }

    setError(null);
    setImporting(true);
    try {
      const filePath = await open({
        multiple: false,
        filters: [{ name: "Haplotypes JSON", extensions: ["json"] }],
      });
      if (!filePath || typeof filePath !== "string") return;

      const content = await readTextFile(filePath);
      const result = await invoke<HaplotypesParseResult>("parse_haplotypes_json", {
        jsonContent: content,
        pedigreeTgf,
        markerSetName: selectedKitName ?? undefined,
        markerSetCsv: markerSetCsv ?? undefined,
      });
      setHaplotypes(result, content);
    } catch (e) {
      setError(String(e));
    } finally {
      setImporting(false);
    }
  };

  // ------------------------------------------------------------------
  // Export haplotypes to JSON file
  // ------------------------------------------------------------------
  const handleExport = async () => {
    if (!haplotypes) return;
    setError(null);
    try {
      const json = await invoke<string>("export_haplotypes_json", {
        haplotypeTable: haplotypes.haplotypeTable,
        traceHaplotype: haplotypes.traceHaplotype,
      });
      const filePath = await save({
        filters: [{ name: "Haplotypes JSON", extensions: ["json"] }],
        defaultPath: "haplotypes.json",
      });
      if (!filePath) return;
      await writeTextFile(filePath, json);
    } catch (e) {
      setError(String(e));
    }
  };

  // ------------------------------------------------------------------
  // Render
  // ------------------------------------------------------------------

  if (!pedigree) {
    return (
      <div className="p-8 flex flex-col items-center justify-center gap-3 text-gray-400 h-full">
        <p className="text-sm">No pedigree loaded.</p>
        <p className="text-xs text-gray-400">Go to <strong className="text-gray-500">Pedigree</strong> to load one first.</p>
      </div>
    );
  }

  const markerNames = haplotypes?.markerNames ?? [];
  const individualNames = pedigree.individuals.map((i) => i.name);
  const knownNames = pedigree.individuals
    .filter((i) => i.haplotypeClass !== "unknown" && i.haplotypeClass !== "excluded")
    .map((i) => i.name);

  return (
    <div className="h-full flex flex-col">
      {/* Toolbar */}
      <div className="bg-white border-b px-4 py-2 flex items-center gap-2">
        <button
          onClick={handleImport}
          disabled={importing}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50 disabled:opacity-50"
        >
          {importing ? "Importing…" : "Import JSON"}
        </button>
        <button
          onClick={handleExport}
          disabled={!haplotypes}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50 disabled:opacity-50"
        >
          Export JSON
        </button>

        {haplotypes && (
          <span className="text-xs text-gray-500 ml-2">
            {markerNames.length} markers · {Object.keys(haplotypes.haplotypeTable).length} individuals typed
            {haplotypes.traceHaplotype && " · TRACE profile present"}
          </span>
        )}
      </div>

      {error && (
        <div className="bg-red-50 border-b border-red-200 px-4 py-2 text-sm text-red-700">
          {error}
        </div>
      )}

      {!haplotypes ? (
        <div className="flex-1 flex flex-col items-center justify-center gap-4 text-gray-400">
          <p className="text-sm">No haplotypes loaded</p>
          <button
            onClick={handleImport}
            disabled={importing}
            className="text-sm bg-blue-600 text-white rounded px-4 py-2 hover:bg-blue-700 disabled:opacity-50"
          >
            Import haplotypes JSON
          </button>
          {(!selectedKitName && !markerSetCsv) && (
            <p className="text-xs text-amber-600">
              Select a marker set on the <strong>Marker Sets</strong> page first.
            </p>
          )}
        </div>
      ) : (
        <div className="flex-1 overflow-auto p-4">
          {/* Spreadsheet-style table: rows = markers, columns = individuals */}
          <div className="overflow-auto border rounded shadow-sm">
            <table className="text-xs border-collapse min-w-full">
              <thead>
                <tr className="bg-gray-50 sticky top-0 z-10">
                  <th className="sticky left-0 z-20 bg-gray-50 border px-3 py-2 text-left font-semibold text-gray-600 whitespace-nowrap">
                    Marker
                  </th>
                  {haplotypes.traceHaplotype && (
                    <th className="border px-3 py-2 text-center font-semibold text-purple-700 bg-purple-50 whitespace-nowrap">
                      TRACE
                    </th>
                  )}
                  {individualNames.map((name) => {
                    const ind = pedigree.individuals.find((i) => i.name === name);
                    const cls = ind?.haplotypeClass ?? "unknown";
                    const isKnown = cls !== "unknown" && cls !== "excluded";
                    return (
                      <th
                        key={name}
                        className={`border px-3 py-2 text-center font-medium whitespace-nowrap ${
                          isKnown ? "text-green-800 bg-green-50" : "text-gray-500"
                        }`}
                      >
                        {name}
                        <div className="font-normal text-gray-400 text-xs">{cls}</div>
                      </th>
                    );
                  })}
                </tr>
              </thead>
              <tbody>
                {markerNames.map((marker, rowIdx) => {
                  // Check for allelic diversity across known individuals
                  const knownAlleles = knownNames
                    .map((n) => haplotypes.haplotypeTable[n]?.[marker])
                    .filter(Boolean);
                  const unique = new Set(knownAlleles);
                  const hasDiversity = unique.size > 1;

                  return (
                    <tr
                      key={marker}
                      className={`${hasDiversity ? "bg-red-50" : rowIdx % 2 === 0 ? "bg-white" : "bg-gray-50"} hover:bg-yellow-50`}
                    >
                      <td className="sticky left-0 z-10 bg-inherit border px-3 py-1.5 font-mono text-gray-800 whitespace-nowrap">
                        {marker}
                        {hasDiversity && (
                          <span className="ml-1 text-red-500 text-xs" title="Allelic diversity detected">●</span>
                        )}
                      </td>
                      {haplotypes.traceHaplotype && (
                        <td className="border px-2 py-1 text-center bg-purple-50 font-mono text-purple-800">
                          {haplotypes.traceHaplotype[marker] ?? "—"}
                        </td>
                      )}
                      {individualNames.map((name) => {
                        const value = haplotypes.haplotypeTable[name]?.[marker];
                        return (
                          <td
                            key={name}
                            className={`border px-2 py-1 text-center font-mono ${
                              value ? "text-gray-800" : "text-gray-300"
                            }`}
                          >
                            {value ?? "—"}
                          </td>
                        );
                      })}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
          <p className="text-xs text-gray-400 mt-2">
            Rows highlighted in red have allelic diversity across typed individuals.
          </p>
        </div>
      )}
    </div>
  );
}
