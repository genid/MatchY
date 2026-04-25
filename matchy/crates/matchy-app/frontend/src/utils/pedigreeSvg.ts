import type { PedigreeData } from "../types/matchy";

const NODE_W = 90;
const NODE_H = 32;
const H_GAP = 24;
const V_GAP = 48;
const PADDING = 20;

/**
 * @param knownNames - set of individual names that have a typed haplotype
 *   (derived from haplotypeTable keys). Any individual NOT in this set is
 *   treated as unknown and rendered in gray.
 */
export function renderPedigreeSvgDataUrl(
  pedigree: PedigreeData,
  suspect: string | null,
  exclude: string[],
  knownNames: Set<string>,
): string {
  const { individuals, relationships } = pedigree;

  const childrenOf = new Map<string, string[]>();
  const parentsOf = new Map<string, string[]>();
  individuals.forEach((i) => {
    childrenOf.set(i.id, []);
    parentsOf.set(i.id, []);
  });
  relationships.forEach((r) => {
    childrenOf.get(r.parentId)?.push(r.childId);
    parentsOf.get(r.childId)?.push(r.parentId);
  });

  // BFS from roots
  const roots = individuals
    .filter((i) => (parentsOf.get(i.id) ?? []).length === 0)
    .map((i) => i.id);
  const layer = new Map<string, number>();
  const queue: { id: string; depth: number }[] = roots.map((id) => ({ id, depth: 0 }));
  while (queue.length > 0) {
    const item = queue.shift()!;
    if (layer.has(item.id)) continue;
    layer.set(item.id, item.depth);
    for (const childId of childrenOf.get(item.id) ?? []) {
      queue.push({ id: childId, depth: item.depth + 1 });
    }
  }
  individuals.forEach((i) => {
    if (!layer.has(i.id)) layer.set(i.id, 0);
  });

  const maxLayer = Math.max(...layer.values(), 0);
  const layerNodes: string[][] = Array.from({ length: maxLayer + 1 }, () => []);
  layer.forEach((l, id) => layerNodes[l].push(id));

  const totalWidth =
    Math.max(...layerNodes.map((ns) => ns.length * (NODE_W + H_GAP) - H_GAP), NODE_W) +
    PADDING * 2;
  const totalHeight = (maxLayer + 1) * (NODE_H + V_GAP) - V_GAP + PADDING * 2;

  const pos = new Map<string, { x: number; y: number }>();
  layerNodes.forEach((nodes, l) => {
    const rowWidth = nodes.length * (NODE_W + H_GAP) - H_GAP;
    const startX = (totalWidth - rowWidth) / 2;
    nodes.forEach((id, i) => {
      pos.set(id, {
        x: startX + i * (NODE_W + H_GAP),
        y: PADDING + l * (NODE_H + V_GAP),
      });
    });
  });

  let paths = "";
  for (const rel of relationships) {
    const p = pos.get(rel.parentId);
    const c = pos.get(rel.childId);
    if (!p || !c) continue;
    const px = p.x + NODE_W / 2;
    const py = p.y + NODE_H;
    const cx = c.x + NODE_W / 2;
    const cy = c.y;
    const mid = (py + cy) / 2;
    paths += `<path d="M${px},${py} L${px},${mid} L${cx},${mid} L${cx},${cy}" stroke="#94a3b8" stroke-width="1.5" fill="none"/>`;
  }

  let nodes = "";
  for (const ind of individuals) {
    const p = pos.get(ind.id);
    if (!p) continue;
    const isExcluded = exclude.includes(ind.name);
    const isSuspect = suspect === ind.name;
    const isKnown = knownNames.has(ind.name);

    const fill = isExcluded
      ? "#f3f4f6"
      : isSuspect
      ? "#dbeafe"
      : isKnown
      ? "#dcfce7"
      : "#ffffff";
    const stroke = isExcluded
      ? "#9ca3af"
      : isSuspect
      ? "#3b82f6"
      : isKnown
      ? "#16a34a"
      : "#9ca3af";
    const textColor = isExcluded
      ? "#6b7280"
      : isSuspect
      ? "#1e40af"
      : isKnown
      ? "#166534"
      : "#6b7280";
    const label = isExcluded
      ? "excluded"
      : isSuspect
      ? "suspect"
      : isKnown
      ? "known"
      : "unknown";
    const strokeDash = !isKnown && !isSuspect && !isExcluded ? 'stroke-dasharray="4 2"' : "";

    const name = ind.name.length > 11 ? ind.name.slice(0, 10) + "…" : ind.name;
    nodes += `<rect x="${p.x}" y="${p.y}" width="${NODE_W}" height="${NODE_H}" rx="4" fill="${fill}" stroke="${stroke}" stroke-width="1.5" ${strokeDash}/>`;
    nodes += `<text x="${p.x + NODE_W / 2}" y="${p.y + 14}" text-anchor="middle" fill="${textColor}" font-weight="600" font-size="10">${escSvg(name)}</text>`;
    nodes += `<text x="${p.x + NODE_W / 2}" y="${p.y + 26}" text-anchor="middle" fill="${textColor}" font-size="8" opacity="0.8">${label}</text>`;
  }

  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${totalWidth}" height="${totalHeight}" viewBox="0 0 ${totalWidth} ${totalHeight}"><style>text{font-family:ui-sans-serif,system-ui,sans-serif}</style>${paths}${nodes}</svg>`;
  return `data:image/svg+xml;base64,${btoa(unescape(encodeURIComponent(svg)))}`;
}

function escSvg(s: string): string {
  return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}
