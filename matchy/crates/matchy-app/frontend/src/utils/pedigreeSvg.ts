import dagre from "dagre";
import type { PedigreeData } from "../types/matchy";

const NODE_W = 148;  // matches PedigreeBuilder NODE_WIDTH
const NODE_H = 44;   // matches PedigreeBuilder NODE_HEIGHT
const PADDING = 24;

// Mirrors CLASS_COLORS from PedigreeBuilder
const CLASS_FILL: Record<string, string> = {
  unknown:   "#e2e8f0",
  known:     "#c6f6d5",
  suspect:   "#fed7d7",
  estimated: "#fefcbf",
  fixed:     "#bee3f8",
  excluded:  "#e9ecf0",
};

export function renderPedigreeSvgDataUrl(
  pedigree: PedigreeData,
  suspect: string | null,
  exclude: string[],
  knownNames: Set<string>,
  classLabels?: { unknown: string; known: string; suspect: string; excluded: string },
): string {
  const labels = classLabels ?? { unknown: "unknown", known: "known", suspect: "suspect", excluded: "excluded" };
  const { individuals, relationships } = pedigree;

  // Use dagre for layout — same algorithm as PedigreeBuilder
  const g = new dagre.graphlib.Graph();
  g.setDefaultEdgeLabel(() => ({}));
  g.setGraph({ rankdir: "TB", nodesep: 60, ranksep: 80 });

  individuals.forEach((ind) => g.setNode(ind.id, { width: NODE_W, height: NODE_H }));
  relationships.forEach((rel) => g.setEdge(rel.parentId, rel.childId));

  dagre.layout(g);

  // Extract top-left positions from dagre centre coordinates
  const pos = new Map<string, { x: number; y: number }>();
  individuals.forEach((ind) => {
    const n = g.node(ind.id);
    if (n) pos.set(ind.id, { x: n.x - NODE_W / 2, y: n.y - NODE_H / 2 });
  });

  // Compute canvas bounds and apply padding offset
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  pos.forEach(({ x, y }) => {
    minX = Math.min(minX, x);
    minY = Math.min(minY, y);
    maxX = Math.max(maxX, x + NODE_W);
    maxY = Math.max(maxY, y + NODE_H);
  });
  const ox = PADDING - minX;
  const oy = PADDING - minY;
  const totalWidth  = maxX - minX + PADDING * 2;
  const totalHeight = maxY - minY + PADDING * 2;

  // Draw edges first (so nodes sit on top) — L-shaped like PedigreeBuilder smoothstep
  let paths = "";
  for (const rel of relationships) {
    const p = pos.get(rel.parentId);
    const c = pos.get(rel.childId);
    if (!p || !c) continue;
    const px = p.x + ox + NODE_W / 2;
    const py = p.y + oy + NODE_H;
    const cx = c.x + ox + NODE_W / 2;
    const cy = c.y + oy;
    const mid = (py + cy) / 2;
    paths += `<path d="M${px},${py} L${px},${mid} L${cx},${mid} L${cx},${cy}" stroke="#94a3b8" stroke-width="1.5" fill="none"/>`;
  }

  // Draw nodes
  let svgNodes = "";
  for (const ind of individuals) {
    const p = pos.get(ind.id);
    if (!p) continue;
    const x = p.x + ox;
    const y = p.y + oy;

    const isExcluded = exclude.includes(ind.name);
    const isSuspect  = suspect === ind.name;
    const isKnown    = knownNames.has(ind.name);

    const cls = isExcluded ? "excluded" : isSuspect ? "suspect" : isKnown ? "known" : "unknown";
    const fill    = CLASS_FILL[cls] ?? "#e2e8f0";
    // Border: dashed gray for excluded, solid blue for suspect, solid green for known, solid gray for unknown
    const stroke      = isExcluded ? "#94a3b8" : isSuspect ? "#3b82f6" : isKnown ? "#4ade80" : "#cbd5e0";
    const strokeW     = isExcluded ? "1.5" : isSuspect ? "2" : "1.5";
    const dash        = isExcluded ? 'stroke-dasharray="5 3"' : "";
    const opacity     = isExcluded ? ' opacity="0.65"' : "";
    const clsLabel    = isExcluded ? labels.excluded : isSuspect ? labels.suspect : isKnown ? labels.known : labels.unknown;

    const name = ind.name.length > 16 ? ind.name.slice(0, 15) + "…" : ind.name;
    // Text color: dark on all backgrounds (matches PedigreeBuilder light mode)
    const txtCol = "#1a202c";

    svgNodes += `<g${opacity}>`;
    svgNodes += `<rect x="${x}" y="${y}" width="${NODE_W}" height="${NODE_H}" rx="8" fill="${fill}" stroke="${stroke}" stroke-width="${strokeW}" ${dash}/>`;
    svgNodes += `<text x="${x + NODE_W / 2}" y="${y + 18}" text-anchor="middle" fill="${txtCol}" font-weight="500" font-size="11">${escSvg(name)}</text>`;
    svgNodes += `<text x="${x + NODE_W / 2}" y="${y + 33}" text-anchor="middle" fill="#64748b" font-size="9">${escSvg(clsLabel)}</text>`;
    svgNodes += `</g>`;
  }

  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${totalWidth}" height="${totalHeight}" viewBox="0 0 ${totalWidth} ${totalHeight}"><style>text{font-family:ui-sans-serif,system-ui,sans-serif}</style>${paths}${svgNodes}</svg>`;
  return `data:image/svg+xml;base64,${btoa(unescape(encodeURIComponent(svg)))}`;
}

function escSvg(s: string): string {
  return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}
