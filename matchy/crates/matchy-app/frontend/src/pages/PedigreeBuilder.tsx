import { useCallback, useEffect, useRef, useState } from "react";
import { useLocation } from "react-router-dom";
import {
  ReactFlow,
  useNodesState,
  useEdgesState,
  Controls,
  Background,
  MiniMap,
  NodeToolbar,
  Handle,
  Position,
  type Node,
  type Edge,
  type NodeProps,
  type OnConnect,
  type OnConnectStart,
  type ReactFlowInstance,
} from "@xyflow/react";
import "@xyflow/react/dist/style.css";
import dagre from "dagre";
import { open, save } from "@tauri-apps/plugin-dialog";
import { readTextFile, writeTextFile } from "@tauri-apps/plugin-fs";
import { invoke } from "@tauri-apps/api/core";
import { useAppStore } from "../store/appStore";
import { useT } from "../i18n";
import type { PedigreeData, RelationshipData } from "../types/matchy";

// ---------------------------------------------------------------------------
// Constants & helpers
// ---------------------------------------------------------------------------

const NODE_WIDTH = 148;
const NODE_HEIGHT = 44;
const NODE_HEIGHT_PROB = 66; // taller when probability badge is shown

function fmtProb(p: number): string {
  if (p === 0) return "0.00";
  const abs = Math.abs(p);
  if (abs >= 1e-3 && abs < 1e3) {
    const mag = Math.floor(Math.log10(abs));
    return p.toFixed(Math.max(2, 1 - mag));
  }
  return p.toExponential(2);
}

function probHeatColor(t: number): string {
  // t: 0 = lowest (green), 1 = highest (red)
  const h = Math.round(120 - 120 * t);
  return `hsl(${h}, 65%, 75%)`;
}

const CLASS_COLORS: Record<string, string> = {
  unknown:   "#e2e8f0",
  known:     "#c6f6d5",
  suspect:   "#fed7d7",
  estimated: "#fefcbf",
  fixed:     "#bee3f8",
  excluded:  "#e9ecf0",
};

// Muted versions used when the probability heat overlay is active so the footer strip dominates
const CLASS_COLORS_MUTED: Record<string, string> = {
  unknown:   "#e2e8f0",
  known:     "#e8f5ec",
  suspect:   "#fce8e8",
  estimated: "#fefdf0",
  fixed:     "#eef6fc",
  excluded:  "#e9ecf0",
};

const CLASS_COLORS_DARK: Record<string, string> = {
  unknown:   "#334155",
  known:     "#14532d",
  suspect:   "#7f1d1d",
  estimated: "#713f12",
  fixed:     "#1e3a5f",
  excluded:  "#1e293b",
};

const CLASS_COLORS_MUTED_DARK: Record<string, string> = {
  unknown:   "#2d3748",
  known:     "#1a3d26",
  suspect:   "#4a1a1a",
  estimated: "#3d2810",
  fixed:     "#1a2d42",
  excluded:  "#1e293b",
};

function applyDagreLayout(nodes: Node[], edges: Edge[], nodeHeight = NODE_HEIGHT): Node[] {
  const g = new dagre.graphlib.Graph();
  g.setDefaultEdgeLabel(() => ({}));
  g.setGraph({ rankdir: "TB", nodesep: 60, ranksep: 80 });
  nodes.forEach((n) => g.setNode(n.id, { width: NODE_WIDTH, height: nodeHeight }));
  edges.forEach((e) => g.setEdge(e.source, e.target));
  dagre.layout(g);
  return nodes.map((n) => {
    const pos = g.node(n.id);
    return { ...n, position: { x: pos.x - NODE_WIDTH / 2, y: pos.y - nodeHeight / 2 } };
  });
}

function collectDescendants(rootId: string, rels: RelationshipData[]): Set<string> {
  const result = new Set([rootId]);
  const queue = [rootId];
  while (queue.length) {
    const cur = queue.shift()!;
    for (const r of rels) {
      if (r.parentId === cur && !result.has(r.childId)) {
        result.add(r.childId);
        queue.push(r.childId);
      }
    }
  }
  return result;
}

function wouldCreateCycle(parentId: string, childId: string, rels: RelationshipData[]): boolean {
  const visited = new Set<string>();
  const queue = [parentId];
  while (queue.length) {
    const cur = queue.shift()!;
    if (cur === childId) return true;
    for (const r of rels) {
      if (r.childId === cur && !visited.has(r.parentId)) {
        visited.add(r.parentId);
        queue.push(r.parentId);
      }
    }
  }
  return false;
}

const EMPTY_PEDIGREE: PedigreeData = { individuals: [], relationships: [] };

// ---------------------------------------------------------------------------
// Custom node
// ---------------------------------------------------------------------------

type PedigreeNodeCallbacks = {
  onRemove: (id: string) => void;
  onRename: (id: string, newName: string) => void;
  onToggleExclude: (id: string) => void;
  onToggleSuspect: (id: string) => void;
};

type PedigreeNodeData = {
  label: string;
  haplotypeClass: string;
  exclude: boolean;
  isSuspect: boolean;
  hasHaplotype: boolean;
  probOverlayActive: boolean;
  cb: React.MutableRefObject<PedigreeNodeCallbacks>;
  matchProb?: number | null;
  probColor?: string;
};

type PedigreeNodeType = Node<PedigreeNodeData, "pedigree">;

function PedigreeNode({ id, data, selected }: NodeProps<PedigreeNodeType>) {
  const [hovered, setHovered] = useState(false);
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(data.label);
  const t = useT();
  const darkMode = useAppStore((s) => s.darkMode);
  const hideTimer = useRef<ReturnType<typeof setTimeout> | null>(null);

  const showToolbar = selected || hovered;

  const cancelHide = () => {
    if (hideTimer.current !== null) {
      clearTimeout(hideTimer.current);
      hideTimer.current = null;
    }
  };

  const scheduleHide = () => {
    cancelHide();
    hideTimer.current = setTimeout(() => setHovered(false), 900);
  };

  const commitRename = () => {
    const name = draft.trim();
    if (name && name !== data.label) {
      data.cb.current.onRename(id, name);
    } else {
      setDraft(data.label);
    }
    setEditing(false);
  };

  return (
    <>
      <NodeToolbar isVisible={showToolbar} position={Position.Top}>
        <div
          className="flex items-center gap-1 bg-white border border-gray-200 rounded-lg shadow-lg px-2 py-1.5"
          onMouseEnter={() => { cancelHide(); setHovered(true); }}
          onMouseLeave={scheduleHide}
        >
          <button
            onMouseDown={(e) => e.stopPropagation()}
            onClick={() => data.cb.current.onToggleExclude(id)}
            className={`text-xs rounded px-1.5 py-0.5 border whitespace-nowrap transition-colors ${
              data.exclude
                ? "bg-gray-200 text-gray-600 border-gray-300 hover:bg-green-50 hover:text-green-700 hover:border-green-300"
                : "bg-white text-gray-500 border-gray-200 hover:bg-orange-50 hover:text-orange-700 hover:border-orange-200"
            }`}
            title={data.exclude ? t("ped_include") : t("ped_exclude")}
          >
            {data.exclude ? t("ped_include") : t("ped_exclude")}
          </button>
          {data.hasHaplotype && (
            <button
              onMouseDown={(e) => e.stopPropagation()}
              onClick={() => data.cb.current.onToggleSuspect(id)}
              className={`text-xs rounded px-1.5 py-0.5 border whitespace-nowrap transition-colors ${
                data.isSuspect
                  ? "bg-red-100 text-red-700 border-red-300 hover:bg-red-50"
                  : "bg-white text-gray-500 border-gray-200 hover:bg-red-50 hover:text-red-700 hover:border-red-200"
              }`}
              title={data.isSuspect ? t("ped_clear_suspect") : t("ped_suspect")}
            >
              {data.isSuspect ? t("ped_suspect_active") : t("ped_suspect")}
            </button>
          )}
          <button
            onMouseDown={(e) => e.stopPropagation()}
            onClick={() => data.cb.current.onRemove(id)}
            className="text-xs bg-red-50 border border-red-200 text-red-600 rounded px-1.5 py-0.5 hover:bg-red-100"
            title={t("ped_remove")}
          >
            {t("ped_remove")}
          </button>
        </div>
      </NodeToolbar>

      {/* Top handle — drag from here to add/connect a parent */}
      <Handle
        type="target"
        position={Position.Top}
        style={{ background: "#94a3b8", width: 10, height: 10, top: -5 }}
        title={t("ped_drag_parent")}
      />

      <div
        onMouseEnter={() => { cancelHide(); setHovered(true); }}
        onMouseLeave={scheduleHide}
        style={{
          background: (data.probOverlayActive
            ? (darkMode ? CLASS_COLORS_MUTED_DARK : CLASS_COLORS_MUTED)[data.haplotypeClass]
            : (darkMode ? CLASS_COLORS_DARK : CLASS_COLORS)[data.haplotypeClass]) ?? (darkMode ? "#334155" : "#e2e8f0"),
          border: selected
            ? "2px solid #3b82f6"
            : data.exclude
            ? "2px dashed #94a3b8"
            : "1.5px solid #cbd5e0",
          borderRadius: 8,
          padding: data.matchProb !== undefined ? "7px 14px 0" : "8px 14px",
          minWidth: NODE_WIDTH,
          textAlign: "center",
          opacity: data.exclude ? 0.65 : 1,
          cursor: "default",
          userSelect: "none",
          overflow: "hidden",
        }}
      >
        {editing ? (
          <input
            autoFocus
            value={draft}
            onChange={(e) => setDraft(e.target.value)}
            onBlur={commitRename}
            onKeyDown={(e) => {
              if (e.key === "Enter") commitRename();
              if (e.key === "Escape") { setDraft(data.label); setEditing(false); }
              e.stopPropagation();
            }}
            className="text-xs text-center bg-transparent border-b border-gray-500 outline-none w-full"
          />
        ) : (
          <span
            className="text-xs font-medium block"
            style={{ color: darkMode ? "#e2e8f0" : "#1a202c" }}
            title={t("ped_double_click_rename")}
            onDoubleClick={() => { setDraft(data.label); setEditing(true); }}
          >
            {data.label}
          </span>
        )}
        {data.matchProb !== undefined && (
          <div
            style={{
              margin: "5px -14px 0",
              padding: "2px 6px 4px",
              backgroundColor: data.probColor ?? "#e2e8f0",
            }}
            title={`P(match) = ${data.matchProb ?? "n/a"}`}
          >
            <span className="text-[10px] font-mono font-semibold tracking-tight" style={{ color: "#1a202c" }}>
              {data.matchProb !== null ? fmtProb(data.matchProb) : "—"}
            </span>
          </div>
        )}
      </div>

      {/* Bottom handle — drag from here to add/connect a child */}
      <Handle
        type="source"
        position={Position.Bottom}
        style={{ background: "#94a3b8", width: 10, height: 10, bottom: -5 }}
        title={t("ped_drag_child")}
      />
    </>
  );
}

const NODE_TYPES = { pedigree: PedigreeNode };

// ---------------------------------------------------------------------------
// Main component
// ---------------------------------------------------------------------------

type DropDialog = {
  /** Screen-relative position for the dialog */
  x: number;
  y: number;
  fromNodeId: string;
  fromNodeName: string;
  /** 'source' = bottom handle → new node is a child
   *  'target' = top handle → new node is a parent */
  handleType: "source" | "target";
};

type ConfirmDialog = {
  message: string;
  onConfirm: () => void | Promise<void>;
};

export default function PedigreeBuilder() {
  const { pedigree, setPedigree, clearPedigree, haplotypes, suspect, setSuspect, exclude, setExclude, simulation, darkMode } = useAppStore();
  const location = useLocation();
  const t = useT();
  const perIndProbs = simulation.result?.per_individual_probabilities ?? null;
  const [showProbOverlay, setShowProbOverlay] = useState<boolean>(
    !!(location.state as { showProbs?: boolean } | null)?.showProbs && !!perIndProbs
  );
  const [nodes, setNodes, onNodesChange] = useNodesState<Node>([]);
  const [edges, setEdges, onEdgesChange] = useEdgesState<Edge>([]);
  const [error, setError] = useState<string | null>(null);
  const [feedback, setFeedback] = useState<string | null>(null);
  const [importing, setImporting] = useState(false);
  const [history, setHistory] = useState<PedigreeData[]>([]);
  const rfRef = useRef<ReactFlowInstance | null>(null);

  // Founder prompt dialog (shown when canvas is empty and user clicks "Start building")
  const [founderDialogOpen, setFounderDialogOpen] = useState(false);
  const [founderName, setFounderName] = useState("");
  const founderInputRef = useRef<HTMLInputElement>(null);

  // "Drop-to-add" dialog (from dragging a handle to empty canvas)
  const [dropDialog, setDropDialog] = useState<DropDialog | null>(null);
  const [dropName, setDropName] = useState("");
  const dropInputRef = useRef<HTMLInputElement>(null);
  const canvasRef = useRef<HTMLDivElement>(null);

  const [confirmDialog, setConfirmDialog] = useState<ConfirmDialog | null>(null);

  const showConfirmDialog = (message: string, onConfirm: () => void | Promise<void>) => {
    setConfirmDialog({ message, onConfirm });
  };

  // Track which handle the user started dragging from
  const connectStartRef = useRef<{ nodeId: string; handleType: "source" | "target" } | null>(null);

  // Ref always pointing to current pedigree for stable callbacks
  const pedRef = useRef<PedigreeData>(pedigree ?? EMPTY_PEDIGREE);
  useEffect(() => { pedRef.current = pedigree ?? EMPTY_PEDIGREE; }, [pedigree]);

  const showFeedback = (msg: string) => {
    setFeedback(msg);
    setTimeout(() => setFeedback(null), 2500);
  };

  const commitPedigree = useCallback(
    async (updated: PedigreeData, skipHistory = false) => {
      try {
        if (!skipHistory) setHistory((prev) => [...prev.slice(-30), pedRef.current]);
        if (updated.individuals.length === 0) {
          clearPedigree();
        } else {
          const tgf = await invoke<string>("export_tgf", { data: updated });
          setPedigree(updated, tgf);
        }
      } catch (e) {
        setError(String(e));
      }
    },
    [setPedigree, clearPedigree]
  );

  const handleUndo = useCallback(async () => {
    setHistory((prev) => {
      if (prev.length === 0) return prev;
      const restored = prev[prev.length - 1];
      invoke<string>("export_tgf", { data: restored })
        .then((tgf) => { setPedigree(restored, tgf); setFeedback(t("ped_feedback_undone")); setTimeout(() => setFeedback(null), 2500); })
        .catch((e) => setError(String(e)));
      return prev.slice(0, -1);
    });
  }, [setPedigree, t]);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.ctrlKey || e.metaKey) && e.key === "z" && !e.shiftKey) {
        e.preventDefault();
        handleUndo();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [handleUndo]);

  // Stable callbacks ref
  const cbRef = useRef<PedigreeNodeCallbacks>({
    onRemove: async () => {},
    onRename: async () => {},
    onToggleExclude: () => {},
    onToggleSuspect: () => {},
  });

  cbRef.current = {
    onRemove: async (id) => {
      const cur = pedRef.current;
      const directChildren = cur.relationships.filter((r) => r.parentId === id);
      const hasParent = cur.relationships.some((r) => r.childId === id);
      const name = cur.individuals.find((i) => i.id === id)?.name ?? id;

      // Founder with exactly one child: just remove this node, child becomes new root
      if (!hasParent && directChildren.length === 1) {
        const newInds = cur.individuals.filter((i) => i.id !== id);
        const newRels = cur.relationships.filter((r) => r.parentId !== id);
        await commitPedigree({ individuals: newInds, relationships: newRels });
        showFeedback(t("ped_feedback_removed").replace("{name}", name));
        return;
      }

      const toRemove = collectDescendants(id, cur.relationships);

      const doRemove = async () => {
        const newInds = cur.individuals.filter((i) => !toRemove.has(i.id));
        const newRels = cur.relationships.filter(
          (r) => !toRemove.has(r.parentId) && !toRemove.has(r.childId)
        );
        await commitPedigree({ individuals: newInds, relationships: newRels });
        const extra = toRemove.size > 1 ? " " + t("ped_feedback_descendants").replace("{n}", String(toRemove.size - 1)) : "";
        showFeedback(t("ped_feedback_removed").replace("{name}", name) + extra);
      };

      if (toRemove.size > 1) {
        showConfirmDialog(t("ped_confirm_remove_descendants"), doRemove);
        return;
      }
      await doRemove();
    },

    onRename: async (id, newNameVal) => {
      const cur = pedRef.current;
      if (cur.individuals.some((i) => i.id !== id && (i.name === newNameVal || i.id === newNameVal))) {
        setError(t("ped_error_already_exists").replace("{name}", newNameVal));
        return;
      }
      const newInds = cur.individuals.map((i) =>
        i.id === id ? { ...i, id: newNameVal, name: newNameVal } : i
      );
      const newRels = cur.relationships.map((r) => ({
        parentId: r.parentId === id ? newNameVal : r.parentId,
        childId: r.childId === id ? newNameVal : r.childId,
      }));
      await commitPedigree({ individuals: newInds, relationships: newRels });
      showFeedback(t("ped_feedback_renamed").replace("{name}", newNameVal));
    },

    onToggleExclude: (id) => {
      const cur = pedRef.current;
      const name = cur.individuals.find((i) => i.id === id)?.name ?? id;
      setExclude(
        exclude.includes(name) ? exclude.filter((n) => n !== name) : [...exclude, name]
      );
    },

    onToggleSuspect: (id) => {
      const cur = pedRef.current;
      const name = cur.individuals.find((i) => i.id === id)?.name ?? id;
      setSuspect(suspect === name ? null : name);
    },
  };

  // Build ReactFlow nodes/edges from pedigree data, overlaying live store state for colors
  const toFlow = useCallback(
    (ped: PedigreeData): { nodes: Node[]; edges: Edge[] } => {
      // Compute probability scale for overlay
      const activeProbs = showProbOverlay && perIndProbs ? perIndProbs : null;
      const probNums = activeProbs
        ? Object.values(activeProbs).map((v) => parseFloat(v as string)).filter((n) => isFinite(n))
        : [];
      const minP = probNums.length > 0 ? Math.min(...probNums) : 0;
      const maxP = probNums.length > 0 ? Math.max(...probNums) : 1;
      const probRange = maxP - minP;

      const nodeH = activeProbs ? NODE_HEIGHT_PROB : NODE_HEIGHT;

      const rawNodes: Node[] = ped.individuals.map((ind) => {
        // Derive effective class from live store values (highest-priority first)
        let cls: string;
        if (exclude.includes(ind.name)) {
          cls = "excluded";
        } else if (suspect === ind.name) {
          cls = "suspect";
        } else if (haplotypes?.haplotypeTable[ind.name]) {
          cls = "known";
        } else if (haplotypes) {
          cls = "unknown"; // haplotypes loaded but this individual has none
        } else {
          cls = ind.haplotypeClass; // no haplotypes loaded yet — keep imported class
        }

        let matchProb: number | null | undefined = undefined;
        let probColor: string | undefined = undefined;
        if (activeProbs) {
          const raw = activeProbs[ind.name];
          if (raw !== undefined) {
            matchProb = parseFloat(raw as string);
            const t = probRange > 0 ? (matchProb - minP) / probRange : 0.5;
            probColor = probHeatColor(t);
          } else {
            matchProb = null;
            probColor = "#e2e8f0";
          }
        }

        return {
          id: ind.id,
          type: "pedigree",
          data: { label: ind.name, haplotypeClass: cls, exclude: exclude.includes(ind.name), isSuspect: suspect === ind.name, hasHaplotype: !!(haplotypes?.haplotypeTable[ind.name]), probOverlayActive: !!activeProbs, cb: cbRef, matchProb, probColor },
          position: { x: 0, y: 0 },
          style: { width: NODE_WIDTH },
        };
      });
      const edges: Edge[] = ped.relationships.map((r) => ({
        id: `${r.parentId}-${r.childId}`,
        source: r.parentId,
        target: r.childId,
        type: "smoothstep",
        style: { stroke: "#94a3b8" },
      }));
      return { nodes: applyDagreLayout(rawNodes, edges, nodeH), edges };
    },
    [haplotypes, suspect, exclude, showProbOverlay, perIndProbs]
  );

  // Sync ReactFlow whenever the pedigree store changes
  useEffect(() => {
    const { nodes: n, edges: e } = toFlow(pedigree ?? EMPTY_PEDIGREE);
    setNodes(n);
    setEdges(e);
    if (n.length > 0) setTimeout(() => rfRef.current?.fitView({ padding: 0.25, duration: 300 }), 60);
  }, [pedigree, toFlow, setNodes, setEdges]);

  // Record which node/handle the drag started from
  const onConnectStart = useCallback<OnConnectStart>((_event, { nodeId, handleType }) => {
    if (nodeId && handleType) {
      connectStartRef.current = {
        nodeId,
        handleType: handleType as "source" | "target",
      };
    }
  }, []);

  // Drag ended on empty canvas → prompt to add a new node
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const onConnectEnd = useCallback((event: MouseEvent | TouchEvent, connectionState: any) => {
    // If the connection was completed to an existing node, onConnect handles it
    if (connectionState?.isValid) {
      connectStartRef.current = null;
      return;
    }
    const info = connectStartRef.current;
    connectStartRef.current = null;
    if (!info) return;

    const cur = pedRef.current;
    const fromNode = cur.individuals.find((i) => i.id === info.nodeId);
    if (!fromNode) return;

    const canvas = canvasRef.current;
    const rect = canvas?.getBoundingClientRect() ?? { left: 0, top: 0 };
    const clientEvent = event instanceof MouseEvent ? event : (event as TouchEvent).changedTouches[0];
    const x = clientEvent.clientX - rect.left;
    const y = clientEvent.clientY - rect.top;

    setDropDialog({
      x,
      y,
      fromNodeId: info.nodeId,
      fromNodeName: fromNode.name,
      handleType: info.handleType,
    });
    setDropName("");
    setTimeout(() => dropInputRef.current?.focus(), 40);
  }, []);

  // Drag handle → handle: add relationship between existing nodes
  const onConnect = useCallback<OnConnect>(
    async ({ source, target }) => {
      connectStartRef.current = null;
      if (!source || !target || source === target) return;
      const cur = pedRef.current;
      if (cur.relationships.some((r) => r.parentId === source && r.childId === target)) return;
      if (cur.relationships.some((r) => r.childId === target)) {
        setError(t("ped_error_already_parent").replace("{name}", target));
        return;
      }
      if (wouldCreateCycle(source, target, cur.relationships)) {
        setError(t("ped_error_cycle"));
        return;
      }
      await commitPedigree({ ...cur, relationships: [...cur.relationships, { parentId: source, childId: target }] });
      showFeedback(t("ped_feedback_connected").replace("{src}", source).replace("{tgt}", target));
    },
    [commitPedigree]
  );

  // Edge deleted via Delete key → remove relationship
  const onEdgesDelete = useCallback(
    async (deleted: Edge[]) => {
      const cur = pedRef.current;
      const deletedSet = new Set(deleted.map((e) => `${e.source}||${e.target}`));
      const newRels = cur.relationships.filter(
        (r) => !deletedSet.has(`${r.parentId}||${r.childId}`)
      );
      await commitPedigree({ ...cur, relationships: newRels });
      showFeedback(t("ped_feedback_connections").replace("{n}", String(deleted.length)));
    },
    [commitPedigree, t]
  );

  // Node deleted via Delete key → smart cascade with warning
  const onNodesDelete = useCallback(
    async (deleted: Node[]) => {
      const cur = pedRef.current;

      // Check whether any deletion would cascade (excluding the founder-with-one-child case)
      let wouldCascade = false;
      for (const n of deleted) {
        const directChildren = cur.relationships.filter((r) => r.parentId === n.id);
        const hasParent = cur.relationships.some((r) => r.childId === n.id);
        if (!(!hasParent && directChildren.length === 1)) {
          if (collectDescendants(n.id, cur.relationships).size > 1) { wouldCascade = true; break; }
        }
      }
      if (wouldCascade) {
        showConfirmDialog(t("ped_confirm_remove_descendants"), async () => {
          const toRemoveK = new Set<string>();
          for (const n of deleted) {
            const directChildren = cur.relationships.filter((r) => r.parentId === n.id);
            const hasParent = cur.relationships.some((r) => r.childId === n.id);
            if (!hasParent && directChildren.length === 1) {
              toRemoveK.add(n.id);
            } else {
              collectDescendants(n.id, cur.relationships).forEach((id) => toRemoveK.add(id));
            }
          }
          const newInds = cur.individuals.filter((i) => !toRemoveK.has(i.id));
          const newRels = cur.relationships.filter(
            (r) => !toRemoveK.has(r.parentId) && !toRemoveK.has(r.childId)
          );
          await commitPedigree({ individuals: newInds, relationships: newRels });
          showFeedback(t("ped_feedback_individuals").replace("{n}", String(toRemoveK.size)));
        });
        return;
      }

      const toRemove = new Set<string>();
      for (const n of deleted) {
        const directChildren = cur.relationships.filter((r) => r.parentId === n.id);
        const hasParent = cur.relationships.some((r) => r.childId === n.id);
        if (!hasParent && directChildren.length === 1) {
          toRemove.add(n.id); // founder with one child: only remove the founder
        } else {
          collectDescendants(n.id, cur.relationships).forEach((id) => toRemove.add(id));
        }
      }

      const newInds = cur.individuals.filter((i) => !toRemove.has(i.id));
      const newRels = cur.relationships.filter(
        (r) => !toRemove.has(r.parentId) && !toRemove.has(r.childId)
      );
      await commitPedigree({ individuals: newInds, relationships: newRels });
      showFeedback(t("ped_feedback_individuals").replace("{n}", String(toRemove.size)));
    },
    [commitPedigree, t]
  );

  // ------------------------------------------------------------------
  // Panel / dialog actions
  // ------------------------------------------------------------------

  const openFounderDialog = () => {
    setFounderName("");
    setFounderDialogOpen(true);
    setTimeout(() => founderInputRef.current?.focus(), 40);
  };

  const handleAddFounder = async () => {
    const name = founderName.trim();
    if (!name) return;
    const cur = pedRef.current;
    if (cur.individuals.some((i) => i.name === name || i.id === name)) {
      setError(t("ped_error_already_exists").replace("{name}", name));
      return;
    }
    setError(null);
    const newInd = { id: name, name, haplotypeClass: "unknown" as const, exclude: false };
    await commitPedigree({ ...cur, individuals: [...cur.individuals, newInd] });
    setFounderDialogOpen(false);
    setFounderName("");
    showFeedback(t("ped_feedback_added_founder").replace("{name}", name));
  };

  const addIndividualWithRelationship = async (
    name: string,
    relationship: { parentId: string; childId: string } | null
  ) => {
    const cur = pedRef.current;
    if (cur.individuals.some((i) => i.name === name || i.id === name)) {
      setError(t("ped_error_already_exists").replace("{name}", name));
      return false;
    }
    setError(null);
    const newInd = { id: name, name, haplotypeClass: "unknown" as const, exclude: false };
    const newInds = [...cur.individuals, newInd];
    const newRels = relationship ? [...cur.relationships, relationship] : [...cur.relationships];
    await commitPedigree({ individuals: newInds, relationships: newRels });
    showFeedback(t("ped_feedback_added").replace("{name}", name));
    return true;
  };

  const handleConfirmDrop = async () => {
    if (!dropDialog) return;
    const name = dropName.trim();
    if (!name) return;
    const relationship =
      dropDialog.handleType === "source"
        ? { parentId: dropDialog.fromNodeId, childId: name }   // bottom handle → child
        : { parentId: name, childId: dropDialog.fromNodeId };  // top handle → parent

    // Validate for parent case (dragging top handle to empty canvas → new node becomes parent)
    if (dropDialog.handleType === "target") {
      const cur = pedRef.current;
      if (cur.relationships.some((r) => r.childId === dropDialog.fromNodeId)) {
        setError(t("ped_error_already_parent").replace("{name}", dropDialog.fromNodeName));
        return;
      }
      if (wouldCreateCycle(name, dropDialog.fromNodeId, cur.relationships)) {
        setError(t("ped_error_add_cycle"));
        return;
      }
    }

    const ok = await addIndividualWithRelationship(name, relationship);
    if (ok) {
      setDropDialog(null);
      setDropName("");
    }
  };

  const handleImport = async () => {
    setError(null);
    setImporting(true);
    try {
      const filePath = await open({
        multiple: false,
        filters: [{ name: "Pedigree files", extensions: ["tgf", "ped"] }],
      });
      if (!filePath || typeof filePath !== "string") return;
      const content = await readTextFile(filePath);
      const ext = filePath.split(".").pop()?.toLowerCase();
      const command = ext === "ped" ? "parse_ped" : "parse_tgf";
      const data = await invoke<PedigreeData>(command, {
        [ext === "ped" ? "pedContent" : "tgfContent"]: content,
      });
      setPedigree(data, content);
      showFeedback(t("ped_feedback_imported").replace("{n}", String(data.individuals.length)));
    } catch (e) {
      setError(String(e));
    } finally {
      setImporting(false);
    }
  };

  const handleExport = async () => {
    const cur = pedRef.current;
    if (!cur.individuals.length) return;
    setError(null);
    try {
      const tgf = await invoke<string>("export_tgf", { data: cur });
      const filePath = await save({
        filters: [{ name: "TGF Pedigree", extensions: ["tgf"] }],
        defaultPath: "pedigree.tgf",
      });
      if (!filePath) return;
      await writeTextFile(filePath, tgf);
      showFeedback(t("ped_feedback_saved"));
    } catch (e) {
      setError(String(e));
    }
  };

  const current = pedigree ?? EMPTY_PEDIGREE;
  const individuals = current.individuals;

  // ------------------------------------------------------------------
  // Render
  // ------------------------------------------------------------------

  return (
    <div className="h-full flex flex-row overflow-hidden">
      {/* ── Left panel ── */}
      <div className="w-48 flex-shrink-0 bg-white border-r flex flex-col overflow-y-auto text-xs">

        <div className="p-3 border-b">
          <p className="font-semibold text-gray-700 text-sm mb-2">{t("ped_pedigree")}</p>
          <div className="flex gap-1 mb-1">
            <button onClick={handleImport} disabled={importing}
              className="flex-1 bg-white border rounded px-2 py-1.5 hover:bg-gray-50 disabled:opacity-50">
              {importing ? "…" : t("ped_import")}
            </button>
            <button onClick={handleExport} disabled={!individuals.length}
              className="flex-1 bg-white border rounded px-2 py-1.5 hover:bg-gray-50 disabled:opacity-50">
              {t("ped_export")}
            </button>
          </div>
          <div className="flex gap-1 mt-1">
            <button
              onClick={() => {
                const doNew = () => {
                  commitPedigree(EMPTY_PEDIGREE, true);
                  setHistory([]);
                  openFounderDialog();
                };
                if (individuals.length === 0) {
                  doNew();
                } else {
                  showConfirmDialog(t("ped_confirm_clear"), doNew);
                }
              }}
              className="flex-1 text-xs bg-blue-600 text-white rounded px-2 py-1.5 hover:bg-blue-700"
            >
              {t("ped_new")}
            </button>
            <button
              onClick={handleUndo}
              disabled={history.length === 0}
              title={t("ped_undo")}
              className="text-xs bg-white border rounded px-2 py-1.5 hover:bg-gray-50 disabled:opacity-40"
            >
              {t("ped_undo")}
            </button>
          </div>
        </div>

        {feedback && (
          <div className="mx-2 mt-2 px-2 py-1 bg-green-50 border border-green-200 rounded text-green-700">
            {feedback}
          </div>
        )}
        {error && (
          <div className="mx-2 mt-2 px-2 py-1 bg-red-50 border border-red-200 rounded text-red-700">
            {error}
            <button className="ml-1 underline" onClick={() => setError(null)}>✕</button>
          </div>
        )}

        {/* Tips */}
        <div className="p-3 border-b text-gray-400 space-y-1 leading-relaxed">
          <p className="font-semibold text-gray-500">{t("ped_how_to_use")}</p>
          <p>{t("ped_tip_hover")}</p>
          <p>{t("ped_tip_remove")}</p>
          <p>{t("ped_tip_rename")}</p>
          <p>{t("ped_tip_drag_bottom")}</p>
          <p>{t("ped_tip_drag_top")}</p>
          <p>{t("ped_tip_drag_connect")}</p>
          <p>{t("ped_tip_delete")}</p>
          <p>{t("ped_tip_ctrl_z")}</p>
        </div>

        {/* Probability overlay toggle */}
        {perIndProbs && (
          <div className="p-3 border-t">
            <button
              onClick={() => setShowProbOverlay((v) => !v)}
              className={`w-full text-xs rounded px-2 py-1.5 border transition-colors ${
                showProbOverlay
                  ? "bg-indigo-600 text-white border-indigo-600 hover:bg-indigo-700"
                  : "bg-white text-gray-700 border-gray-300 hover:bg-gray-50"
              }`}
            >
              {showProbOverlay ? t("ped_hide_probs") : t("ped_show_probs")}
            </button>
            {showProbOverlay && (
              <div className="mt-1.5 flex justify-between text-xs text-gray-400">
                <span style={{ color: probHeatColor(0) }}>{t("ped_prob_low")}</span>
                <span style={{ color: probHeatColor(1) }}>{t("ped_prob_high")}</span>
              </div>
            )}
          </div>
        )}

        {/* Legend */}
        <div className="p-3 mt-auto border-t">
          <p className="font-semibold text-gray-500 mb-1.5 uppercase tracking-wide">{t("ped_legend")}</p>
          <div className="space-y-1">
            {Object.entries(CLASS_COLORS)
              .filter(([cls]) => cls !== "estimated" && cls !== "fixed")
              .map(([cls, color]) => (
                <span key={cls} className="flex items-center gap-1.5 text-gray-600">
                  <span style={{ background: color }}
                    className="inline-block w-3 h-3 rounded border border-gray-300 flex-shrink-0" />
                  {t(`ped_class_${cls}` as Parameters<typeof t>[0])}
                </span>
              ))}
          </div>
          {individuals.length > 0 && (
            <p className="text-gray-400 mt-2">
              {t("ped_count_summary")
                .replace("{n}", String(individuals.length))
                .replace("{r}", String(current.relationships.length))}
            </p>
          )}
        </div>
      </div>

      {/* ── Canvas ── */}
      <div className="flex-1 relative" ref={canvasRef}>

        {/* Confirm dialog */}
        {confirmDialog && (
          <div className="absolute inset-0 flex items-center justify-center bg-black/30 z-50">
            <div className="bg-white rounded-xl shadow-2xl p-5 w-72 border">
              <p className="text-sm text-gray-800 mb-4">{confirmDialog.message}</p>
              <div className="flex gap-2">
                <button
                  onClick={() => {
                    const fn = confirmDialog.onConfirm;
                    setConfirmDialog(null);
                    fn();
                  }}
                  className="flex-1 bg-red-600 text-white text-sm rounded px-3 py-1.5 hover:bg-red-700"
                >
                  {t("ped_ok")}
                </button>
                <button
                  onClick={() => setConfirmDialog(null)}
                  className="flex-1 bg-white border text-gray-700 text-sm rounded px-3 py-1.5 hover:bg-gray-50"
                >
                  {t("ped_cancel")}
                </button>
              </div>
            </div>
          </div>
        )}

        {/* Founder prompt dialog */}
        {founderDialogOpen && (
          <div
            className="absolute inset-0 flex items-center justify-center bg-black/20 z-50"
            onClick={(e) => { if (e.target === e.currentTarget) setFounderDialogOpen(false); }}
          >
            <div className="bg-white rounded-xl shadow-2xl p-5 w-64 border">
              <p className="text-sm font-semibold text-gray-800 mb-1">{t("ped_founder_title")}</p>
              <p className="text-xs text-gray-500 mb-3">{t("ped_founder_desc")}</p>
              <input
                ref={founderInputRef}
                autoFocus
                type="text"
                value={founderName}
                onChange={(e) => setFounderName(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") handleAddFounder();
                  if (e.key === "Escape") setFounderDialogOpen(false);
                }}
                placeholder={t("ped_founder_placeholder")}
                className="w-full text-sm border rounded px-3 py-2 mb-3 focus:outline-none focus:ring-2 focus:ring-blue-400"
              />
              <div className="flex gap-2">
                <button onClick={handleAddFounder} disabled={!founderName.trim()}
                  className="flex-1 bg-blue-600 text-white text-sm rounded px-3 py-1.5 hover:bg-blue-700 disabled:opacity-40">
                  {t("ped_start")}
                </button>
                <button onClick={() => setFounderDialogOpen(false)}
                  className="flex-1 bg-white border text-gray-700 text-sm rounded px-3 py-1.5 hover:bg-gray-50">
                  {t("ped_cancel")}
                </button>
              </div>
            </div>
          </div>
        )}

        {/* Drop-to-add dialog (appears near where user released the drag) */}
        {dropDialog && (
          <div
            className="absolute inset-0 z-50"
            onClick={(e) => { if (e.target === e.currentTarget) { setDropDialog(null); setDropName(""); } }}
          >
            <div
              className="absolute bg-white rounded-xl shadow-2xl p-4 w-52 border"
              style={{
                left: Math.min(dropDialog.x, (canvasRef.current?.clientWidth ?? 600) - 220),
                top: Math.min(dropDialog.y, (canvasRef.current?.clientHeight ?? 400) - 160),
              }}
            >
              <p className="text-sm font-semibold text-gray-800 mb-0.5">
                {dropDialog.handleType === "source" ? t("ped_add_child_drop") : t("ped_add_parent_drop")}
              </p>
              <p className="text-xs text-gray-500 mb-2">
                {dropDialog.handleType === "source"
                  ? <>{t("ped_child_of")} <strong>{dropDialog.fromNodeName}</strong></>
                  : <>{t("ped_parent_of")} <strong>{dropDialog.fromNodeName}</strong></>
                }
              </p>
              <input
                ref={dropInputRef}
                autoFocus
                type="text"
                value={dropName}
                onChange={(e) => setDropName(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") handleConfirmDrop();
                  if (e.key === "Escape") { setDropDialog(null); setDropName(""); }
                }}
                placeholder={t("ped_name_placeholder")}
                className="w-full text-sm border rounded px-2 py-1.5 mb-2 focus:outline-none focus:ring-2 focus:ring-blue-400"
              />
              <div className="flex gap-1.5">
                <button onClick={handleConfirmDrop} disabled={!dropName.trim()}
                  className="flex-1 bg-blue-600 text-white text-xs rounded px-2 py-1.5 hover:bg-blue-700 disabled:opacity-40">
                  {t("ped_add")}
                </button>
                <button onClick={() => { setDropDialog(null); setDropName(""); }}
                  className="flex-1 bg-white border text-gray-700 text-xs rounded px-2 py-1.5 hover:bg-gray-50">
                  {t("ped_cancel")}
                </button>
              </div>
            </div>
          </div>
        )}

        {individuals.length > 0 ? (
          <ReactFlow
            nodes={nodes}
            edges={edges}
            onNodesChange={onNodesChange}
            onEdgesChange={onEdgesChange}
            onConnect={onConnect}
            onConnectStart={onConnectStart}
            onConnectEnd={onConnectEnd}
            onNodesDelete={onNodesDelete}
            onEdgesDelete={onEdgesDelete}
            nodeTypes={NODE_TYPES}
            onInit={(instance) => { rfRef.current = instance; }}
            fitView
            fitViewOptions={{ padding: 0.25 }}
            minZoom={0.15}
            nodesDraggable
            nodesConnectable
            deleteKeyCode={["Delete", "Backspace"]}
          >
            <Controls />
            <MiniMap
              nodeStrokeColor={darkMode ? "#475569" : "#cbd5e0"}
              nodeColor={(n) => {
                const d = n.data as PedigreeNodeData;
                const colorMap = darkMode ? CLASS_COLORS_DARK : CLASS_COLORS;
                return (showProbOverlay && d.probColor) ? d.probColor : (colorMap[d.haplotypeClass] ?? (darkMode ? "#334155" : "#e2e8f0"));
              }}
            />
            <Background gap={16} color={darkMode ? "#1e293b" : "#e8ecf0"} />
          </ReactFlow>
        ) : (
          <div className="h-full flex flex-col items-center justify-center gap-4 text-gray-400 p-8">
            <svg className="w-16 h-16 opacity-30" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5}
                d="M17 20h5v-2a3 3 0 00-5.356-1.857M17 20H7m10 0v-2c0-.656-.126-1.283-.356-1.857M7 20H2v-2a3 3 0 015.356-1.857M7 20v-2c0-.656.126-1.283.356-1.857m0 0a5.002 5.002 0 019.288 0" />
            </svg>
            <p className="text-sm font-medium">{t("ped_empty_title")}</p>
            <p className="text-xs text-center max-w-xs">{t("ped_empty_hint")}</p>
            <div className="flex gap-2">
              <button onClick={openFounderDialog}
                className="text-sm bg-blue-600 text-white rounded px-4 py-2 hover:bg-blue-700">
                {t("ped_start_building")}
              </button>
              <button onClick={handleImport}
                className="text-sm bg-white border rounded px-4 py-2 hover:bg-gray-50">
                {t("ped_import_file")}
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
