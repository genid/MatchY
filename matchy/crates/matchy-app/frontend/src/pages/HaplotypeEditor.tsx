import { useAppStore } from "../store/appStore";

export default function HaplotypeEditor() {
  const { haplotypes, markers, pedigree } = useAppStore();

  if (!pedigree) {
    return (
      <div className="p-6 text-gray-500 text-sm">Load a pedigree first.</div>
    );
  }

  if (!haplotypes) {
    return (
      <div className="p-6">
        <p className="text-gray-500 text-sm mb-4">No haplotypes loaded.</p>
        {/* TODO: file import button wired to Tauri dialog + parse_haplotypes_json */}
        <button className="text-sm bg-blue-600 text-white rounded px-3 py-1.5">
          Import JSON
        </button>
      </div>
    );
  }

  const individualNames = pedigree.individuals.map((i) => i.name);
  const markerNames = haplotypes.markerNames;

  return (
    <div className="p-4 overflow-auto">
      <div className="flex items-center justify-between mb-3">
        <h1 className="text-xl font-bold text-gray-900">Haplotype Editor</h1>
        <div className="flex gap-2">
          <button className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50">
            Import JSON
          </button>
          <button className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50">
            Export JSON
          </button>
        </div>
      </div>

      {/* Spreadsheet-style table: rows = markers, columns = individuals + TRACE */}
      <div className="overflow-auto border rounded">
        <table className="text-xs border-collapse min-w-full">
          <thead>
            <tr className="bg-gray-50">
              <th className="sticky left-0 bg-gray-50 border px-3 py-2 text-left font-semibold text-gray-600">
                Marker
              </th>
              {haplotypes.traceHaplotype && (
                <th className="border px-3 py-2 text-center font-semibold text-purple-700 bg-purple-50">
                  TRACE
                </th>
              )}
              {individualNames.map((name) => (
                <th key={name} className="border px-3 py-2 text-center font-medium text-gray-700">
                  {name}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {markerNames.map((marker) => (
              <tr key={marker} className="hover:bg-gray-50">
                <td className="sticky left-0 bg-white border px-3 py-1.5 font-mono text-gray-800">
                  {marker}
                </td>
                {haplotypes.traceHaplotype && (
                  <td className="border px-2 py-1 text-center bg-purple-50 font-mono">
                    {haplotypes.traceHaplotype[marker] ?? "—"}
                  </td>
                )}
                {individualNames.map((name) => (
                  <td key={name} className="border px-2 py-1 text-center font-mono">
                    {/* TODO: make cells editable (TanStack Table) */}
                    {haplotypes.haplotypeTable[name]?.[marker] ?? "—"}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
