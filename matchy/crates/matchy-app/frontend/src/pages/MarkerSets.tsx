import { useEffect, useState } from "react";
import { open } from "@tauri-apps/plugin-dialog";
import { readTextFile } from "@tauri-apps/plugin-fs";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import type { MarkerInfo } from "../types/matchy";

export default function MarkerSets() {
  const { selectedKitName, markers, setMarkerSet } = useAppStore();
  const [kitNames, setKitNames] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [customFileName, setCustomFileName] = useState<string | null>(null);

  // Load kit list on mount
  useEffect(() => {
    invoke<string[]>("list_kits")
      .then((names) => setKitNames(names.sort()))
      .catch((e) => setError(String(e)));
  }, []);

  const handleSelectKit = async (name: string) => {
    if (selectedKitName === name) return; // already active
    setLoading(true);
    setError(null);
    setCustomFileName(null);
    try {
      const markerInfos = await invoke<MarkerInfo[]>("load_kit", { name });
      setMarkerSet(name, markerInfos, null);
    } catch (e) {
      setError(String(e));
    } finally {
      setLoading(false);
    }
  };

  const handleCustomCsv = async () => {
    setError(null);
    try {
      const filePath = await open({
        multiple: false,
        filters: [{ name: "Marker set CSV", extensions: ["csv"] }],
      });
      if (!filePath || typeof filePath !== "string") return;

      const content = await readTextFile(filePath);
      const markerInfos = await invoke<MarkerInfo[]>("load_custom_csv", {
        csvContent: content,
      });

      const fileName = filePath.split(/[\\/]/).pop() ?? filePath;
      setCustomFileName(fileName);
      setMarkerSet(null, markerInfos, content);
    } catch (e) {
      setError(String(e));
    }
  };

  const isCustomActive = !selectedKitName && markers.length > 0;

  return (
    <div className="p-6 max-w-3xl mx-auto space-y-5">
      <h1 className="text-xl font-bold text-gray-900">Marker Sets</h1>

      {error && (
        <div className="rounded bg-red-50 border border-red-200 p-3 text-sm text-red-800">
          {error}
        </div>
      )}

      {/* Active marker set summary */}
      {markers.length > 0 && (
        <div className="rounded bg-blue-50 border border-blue-200 p-3 text-sm text-blue-800">
          <strong>Active:</strong>{" "}
          {selectedKitName ? (
            <>Built-in kit <em>{selectedKitName}</em></>
          ) : (
            <>Custom CSV{customFileName ? ` — ${customFileName}` : ""}</>
          )}{" "}
          ({markers.length} markers)
        </div>
      )}

      {/* Built-in kits */}
      <section className="bg-white rounded-lg border p-4">
        <h2 className="font-semibold text-gray-700 mb-3">Built-in Kits</h2>
        {kitNames.length === 0 ? (
          <p className="text-sm text-gray-400">Loading…</p>
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
                  <span className="ml-2 text-xs text-blue-500">✓ active</span>
                )}
              </button>
            ))}
          </div>
        )}
      </section>

      {/* Custom CSV */}
      <section className="bg-white rounded-lg border p-4">
        <h2 className="font-semibold text-gray-700 mb-1">Custom CSV</h2>
        <p className="text-sm text-gray-500 mb-3">
          CSV with columns: <code className="bg-gray-100 px-1 rounded">marker,mutation_rate</code>
        </p>
        <button
          onClick={handleCustomCsv}
          className={`text-sm border rounded px-3 py-1.5 transition-colors ${
            isCustomActive
              ? "border-blue-400 bg-blue-50 text-blue-700"
              : "bg-white hover:bg-gray-50 text-gray-700"
          }`}
        >
          {isCustomActive
            ? `✓ ${customFileName ?? "Custom CSV"} active`
            : "Upload custom CSV…"}
        </button>
      </section>

      {/* Marker table */}
      {markers.length > 0 && (
        <section className="bg-white rounded-lg border p-4">
          <h2 className="font-semibold text-gray-700 mb-3">
            Markers ({markers.length})
          </h2>
          <div className="overflow-auto max-h-80 border rounded">
            <table className="text-xs w-full border-collapse">
              <thead className="sticky top-0">
                <tr className="bg-gray-50">
                  <th className="border px-3 py-2 text-left">Marker</th>
                  <th className="border px-3 py-2 text-right">Mutation Rate</th>
                  <th className="border px-3 py-2 text-right">Copies</th>
                </tr>
              </thead>
              <tbody>
                {markers.map((m, i) => (
                  <tr key={m.name} className={i % 2 === 0 ? "bg-white" : "bg-gray-50"}>
                    <td className="border px-3 py-1 font-mono">{m.name}</td>
                    <td className="border px-3 py-1 text-right font-mono">
                      {m.mutationRate.toExponential(3)}
                    </td>
                    <td className="border px-3 py-1 text-right">
                      {m.numberOfCopies ?? 1}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      )}
    </div>
  );
}
