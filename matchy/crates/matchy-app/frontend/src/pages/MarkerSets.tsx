import { useEffect, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import type { MarkerInfo } from "../types/matchy";

export default function MarkerSets() {
  const { selectedKitName, markers, setMarkerSet } = useAppStore();
  const [kitNames, setKitNames] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    invoke<string[]>("list_kits")
      .then(setKitNames)
      .catch((e) => setError(String(e)));
  }, []);

  const handleSelectKit = async (name: string) => {
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

  return (
    <div className="p-6 max-w-3xl mx-auto space-y-5">
      <h1 className="text-xl font-bold text-gray-900">Marker Sets</h1>

      {error && (
        <div className="rounded bg-red-50 border border-red-200 p-3 text-sm text-red-800">
          {error}
        </div>
      )}

      <section className="bg-white rounded-lg border p-4">
        <h2 className="font-semibold text-gray-700 mb-3">Built-in Kits</h2>
        <div className="space-y-2">
          {kitNames.map((name) => (
            <div
              key={name}
              className={`flex items-center justify-between px-3 py-2 rounded border cursor-pointer transition-colors ${
                selectedKitName === name
                  ? "border-blue-400 bg-blue-50"
                  : "border-gray-200 hover:bg-gray-50"
              }`}
              onClick={() => handleSelectKit(name)}
            >
              <span className="font-medium text-sm">{name}</span>
              {selectedKitName === name && (
                <span className="text-xs text-blue-600 font-medium">Active ({markers.length} markers)</span>
              )}
            </div>
          ))}
          {kitNames.length === 0 && !loading && (
            <p className="text-sm text-gray-500">No kits available.</p>
          )}
          {loading && <p className="text-sm text-gray-500">Loading…</p>}
        </div>
      </section>

      {selectedKitName && markers.length > 0 && (
        <section className="bg-white rounded-lg border p-4">
          <h2 className="font-semibold text-gray-700 mb-3">
            {selectedKitName} ({markers.length} markers)
          </h2>
          <div className="overflow-auto max-h-64">
            <table className="text-xs w-full border-collapse">
              <thead>
                <tr className="bg-gray-50">
                  <th className="border px-3 py-2 text-left">Marker</th>
                  <th className="border px-3 py-2 text-right">Mutation Rate</th>
                  <th className="border px-3 py-2 text-right">Copies</th>
                </tr>
              </thead>
              <tbody>
                {markers.map((m) => (
                  <tr key={m.name} className="hover:bg-gray-50">
                    <td className="border px-3 py-1 font-mono">{m.name}</td>
                    <td className="border px-3 py-1 text-right">{m.mutationRate.toExponential(3)}</td>
                    <td className="border px-3 py-1 text-right">{m.numberOfCopies ?? 1}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      )}

      <section className="bg-white rounded-lg border p-4">
        <h2 className="font-semibold text-gray-700 mb-2">Custom CSV</h2>
        <p className="text-sm text-gray-500 mb-3">
          Upload a CSV file with columns: <code>marker,mutation_rate</code>
        </p>
        {/* TODO: file upload button wired to Tauri dialog + load_custom_csv */}
        <button className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50">
          Upload custom CSV
        </button>
      </section>
    </div>
  );
}
