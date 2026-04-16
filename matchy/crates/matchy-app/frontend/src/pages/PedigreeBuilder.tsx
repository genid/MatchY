import { useCallback, useEffect } from "react";
import {
  ReactFlow,
  addEdge,
  useNodesState,
  useEdgesState,
  Controls,
  Background,
  type Node,
  type Edge,
  type Connection,
} from "@xyflow/react";
import "@xyflow/react/dist/style.css";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import type { PedigreeData } from "../types/matchy";

// Node color by haplotype class (matches Python visualization.py color scheme)
const CLASS_COLORS: Record<string, string> = {
  unknown: "#e2e8f0",
  known: "#c6f6d5",
  suspect: "#fed7d7",
  estimated: "#fefcbf",
  fixed: "#bee3f8",
  excluded: "#718096",
};

function pedigreeToFlow(pedigree: PedigreeData): { nodes: Node[]; edges: Edge[] } {
  const nodes: Node[] = pedigree.individuals.map((ind) => ({
    id: ind.id,
    data: { label: ind.name, haplotypeClass: ind.haplotypeClass },
    position: { x: 0, y: 0 }, // will be laid out by dagre
    style: {
      background: CLASS_COLORS[ind.haplotypeClass] ?? "#e2e8f0",
      border: "1px solid #cbd5e0",
      borderRadius: "6px",
      padding: "6px 12px",
      fontSize: "12px",
    },
  }));
  const edges: Edge[] = pedigree.relationships.map((rel) => ({
    id: `${rel.parentId}-${rel.childId}`,
    source: rel.parentId,
    target: rel.childId,
    type: "smoothstep",
  }));
  return { nodes, edges };
}

export default function PedigreeBuilder() {
  const { pedigree, setPedigree } = useAppStore();
  const [nodes, setNodes, onNodesChange] = useNodesState([]);
  const [edges, setEdges, onEdgesChange] = useEdgesState([]);

  useEffect(() => {
    if (pedigree) {
      const { nodes: n, edges: e } = pedigreeToFlow(pedigree);
      setNodes(n);
      setEdges(e);
    }
  }, [pedigree]);

  const onConnect = useCallback(
    (connection: Connection) => setEdges((eds) => addEdge(connection, eds)),
    [setEdges]
  );

  const handleFileImport = async () => {
    // TODO: use Tauri dialog to open file, then invoke parse_tgf or parse_ped
  };

  const handleExport = async () => {
    if (!pedigree) return;
    try {
      const tgf = await invoke<string>("export_tgf", { data: pedigree });
      // TODO: trigger file download via Tauri
      console.log("TGF output:", tgf);
    } catch (err) {
      console.error(err);
    }
  };

  return (
    <div className="h-full flex flex-col">
      <div className="bg-white border-b px-4 py-2 flex gap-2">
        <button
          onClick={handleFileImport}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50"
        >
          Import TGF / PED
        </button>
        <button
          onClick={handleExport}
          disabled={!pedigree}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50 disabled:opacity-50"
        >
          Export TGF
        </button>
        <div className="ml-auto flex gap-2 text-xs items-center">
          {Object.entries(CLASS_COLORS).map(([cls, color]) => (
            <span key={cls} className="flex items-center gap-1">
              <span style={{ background: color }} className="inline-block w-3 h-3 rounded border border-gray-300" />
              {cls}
            </span>
          ))}
        </div>
      </div>
      <div className="flex-1">
        <ReactFlow
          nodes={nodes}
          edges={edges}
          onNodesChange={onNodesChange}
          onEdgesChange={onEdgesChange}
          onConnect={onConnect}
          fitView
        >
          <Controls />
          <Background />
        </ReactFlow>
      </div>
    </div>
  );
}
