import { useCallback, useEffect, useState } from "react";
import {
  ReactFlow,
  addEdge,
  useNodesState,
  useEdgesState,
  Controls,
  Background,
  MiniMap,
  type Node,
  type Edge,
  type Connection,
} from "@xyflow/react";
import "@xyflow/react/dist/style.css";
import dagre from "dagre";
import { open, save } from "@tauri-apps/plugin-dialog";
import { readTextFile, writeTextFile } from "@tauri-apps/plugin-fs";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import type { PedigreeData } from "../types/matchy";

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const NODE_WIDTH = 120;
const NODE_HEIGHT = 40;

const CLASS_COLORS: Record<string, string> = {
  unknown:   "#e2e8f0",
  known:     "#c6f6d5",
  suspect:   "#fed7d7",
  estimated: "#fefcbf",
  fixed:     "#bee3f8",
  excluded:  "#e2e8f0",
};

// ---------------------------------------------------------------------------
// Dagre layout
// ---------------------------------------------------------------------------

function applyDagreLayout(
  nodes: Node[],
  edges: Edge[]
): Node[] {
  const g = new dagre.graphlib.Graph();
  g.setDefaultEdgeLabel(() => ({}));
  g.setGraph({ rankdir: "TB", nodesep: 40, ranksep: 60 });

  nodes.forEach((n) => g.setNode(n.id, { width: NODE_WIDTH, height: NODE_HEIGHT }));
  edges.forEach((e) => g.setEdge(e.source, e.target));

  dagre.layout(g);

  return nodes.map((n) => {
    const pos = g.node(n.id);
    return {
      ...n,
      position: {
        x: pos.x - NODE_WIDTH / 2,
        y: pos.y - NODE_HEIGHT / 2,
      },
    };
  });
}

// ---------------------------------------------------------------------------
// Conversion helpers
// ---------------------------------------------------------------------------

function pedigreeToFlow(pedigree: PedigreeData): { nodes: Node[]; edges: Edge[] } {
  const rawNodes: Node[] = pedigree.individuals.map((ind) => ({
    id: ind.id,
    data: { label: ind.name },
    position: { x: 0, y: 0 },
    style: {
      background: CLASS_COLORS[ind.haplotypeClass] ?? "#e2e8f0",
      border: ind.exclude ? "2px dashed #a0aec0" : "1px solid #cbd5e0",
      borderRadius: "6px",
      padding: "4px 10px",
      fontSize: "11px",
      width: `${NODE_WIDTH}px`,
      textAlign: "center" as const,
      opacity: ind.exclude ? 0.6 : 1,
    },
  }));

  const edges: Edge[] = pedigree.relationships.map((rel) => ({
    id: `${rel.parentId}-${rel.childId}`,
    source: rel.parentId,
    target: rel.childId,
    type: "smoothstep",
    style: { stroke: "#94a3b8" },
  }));

  return { nodes: applyDagreLayout(rawNodes, edges), edges };
}

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------

export default function PedigreeBuilder() {
  const { pedigree, pedigreeTgf, setPedigree } = useAppStore();
  const [nodes, setNodes, onNodesChange] = useNodesState([]);
  const [edges, setEdges, onEdgesChange] = useEdgesState([]);
  const [error, setError] = useState<string | null>(null);
  const [importing, setImporting] = useState(false);

  // Re-layout whenever the pedigree changes
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

  // Import TGF or PED via Tauri file dialog
  const handleImport = async () => {
    setError(null);
    setImporting(true);
    try {
      const filePath = await open({
        multiple: false,
        filters: [
          { name: "Pedigree files", extensions: ["tgf", "ped"] },
          { name: "All files", extensions: ["*"] },
        ],
      });
      if (!filePath || typeof filePath !== "string") return;

      const content = await readTextFile(filePath);
      const ext = filePath.split(".").pop()?.toLowerCase();
      const command = ext === "ped" ? "parse_ped" : "parse_tgf";
      const data = await invoke<PedigreeData>(command, {
        [ext === "ped" ? "pedContent" : "tgfContent"]: content,
      });
      setPedigree(data, content);
    } catch (e) {
      setError(String(e));
    } finally {
      setImporting(false);
    }
  };

  // Export current pedigree as TGF
  const handleExport = async () => {
    if (!pedigree) return;
    setError(null);
    try {
      const tgf = await invoke<string>("export_tgf", { data: pedigree });
      const filePath = await save({
        filters: [{ name: "TGF Pedigree", extensions: ["tgf"] }],
        defaultPath: "pedigree.tgf",
      });
      if (!filePath) return;
      await writeTextFile(filePath, tgf);
    } catch (e) {
      setError(String(e));
    }
  };

  return (
    <div className="h-full flex flex-col">
      {/* Toolbar */}
      <div className="bg-white border-b px-4 py-2 flex items-center gap-2 flex-wrap">
        <button
          onClick={handleImport}
          disabled={importing}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50 disabled:opacity-50"
        >
          {importing ? "Importing…" : "Import TGF / PED"}
        </button>
        <button
          onClick={handleExport}
          disabled={!pedigree}
          className="text-sm bg-white border rounded px-3 py-1.5 hover:bg-gray-50 disabled:opacity-50"
        >
          Export TGF
        </button>

        {pedigree && (
          <span className="text-xs text-gray-500 ml-2">
            {pedigree.individuals.length} individuals ·{" "}
            {pedigree.relationships.length} relationships
          </span>
        )}

        {/* Legend */}
        <div className="ml-auto flex gap-3 text-xs items-center">
          {Object.entries(CLASS_COLORS).map(([cls, color]) => (
            <span key={cls} className="flex items-center gap-1">
              <span
                style={{ background: color }}
                className="inline-block w-3 h-3 rounded border border-gray-300"
              />
              <span className="text-gray-600">{cls}</span>
            </span>
          ))}
        </div>
      </div>

      {error && (
        <div className="bg-red-50 border-b border-red-200 px-4 py-2 text-sm text-red-700">
          {error}
        </div>
      )}

      {/* React Flow canvas */}
      <div className="flex-1">
        {pedigree ? (
          <ReactFlow
            nodes={nodes}
            edges={edges}
            onNodesChange={onNodesChange}
            onEdgesChange={onEdgesChange}
            onConnect={onConnect}
            fitView
            fitViewOptions={{ padding: 0.2 }}
            minZoom={0.2}
          >
            <Controls />
            <MiniMap
              nodeStrokeColor="#cbd5e0"
              nodeColor={(n) => (n.style?.background as string) ?? "#e2e8f0"}
            />
            <Background gap={16} color="#e2e8f0" />
          </ReactFlow>
        ) : (
          <div className="h-full flex flex-col items-center justify-center gap-4 text-gray-400">
            <svg className="w-16 h-16 opacity-30" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2" />
            </svg>
            <p className="text-sm">No pedigree loaded</p>
            <button
              onClick={handleImport}
              className="text-sm bg-blue-600 text-white rounded px-4 py-2 hover:bg-blue-700"
            >
              Import TGF or PED file
            </button>
          </div>
        )}
      </div>
    </div>
  );
}
