/**
 * Session save/load for MatchY.
 *
 * A session file (.matchy.json) bundles all inputs, parameters and an optional
 * previous result so any simulation can be resumed exactly where it was left off.
 */

import { open, save } from "@tauri-apps/plugin-dialog";
import { readTextFile, writeTextFile } from "@tauri-apps/plugin-fs";
import { invoke } from "@tauri-apps/api/core";
import type { HaplotypesParseResult, MarkerInfo, ProgressEvent, SimulationResult } from "../types/matchy";

export interface SessionParams {
  twoStepMutationFraction: number;
  batchLength: number;
  convergenceCriterion: number;
  bias: number | null;
  numberOfThreads: number;
  skipInside: boolean;
  skipOutside: boolean;
  traceMode: boolean;
  adaptiveBias: boolean;
  simulationName: string;
  userName: string;
  seed: number | null;
}

export interface SessionData {
  version: number;
  pedigreeTgf: string;
  haplotypesJson: string;
  selectedKitName: string | null;
  markerSetCsv: string | null;
  suspect: string | null;
  exclude: string[];
  params: SessionParams;
  simulationResult: SimulationResult | null;
  simulationProgress: ProgressEvent[];
}

export interface LoadedSession {
  pedigreeData: Awaited<ReturnType<typeof invoke<{ individuals: unknown[]; relationships: unknown[] }>>>;
  haplotypesData: HaplotypesParseResult;
  haplotypesJson: string;
  pedigreeTgf: string;
  markers: MarkerInfo[];
  selectedKitName: string | null;
  markerSetCsv: string | null;
  suspect: string | null;
  exclude: string[];
  params: SessionParams;
  simulationResult: SimulationResult | null;
  simulationProgress: ProgressEvent[];
}

export async function saveSession(
  session: SessionData,
  defaultName: string,
): Promise<void> {
  const filePath = await save({
    defaultPath: `${defaultName}.matchy.json`,
    filters: [{ name: "MatchY Session", extensions: ["matchy.json", "json"] }],
  });
  if (!filePath) return;
  await writeTextFile(filePath, JSON.stringify(session, null, 2));
}

export async function loadSession(): Promise<LoadedSession | null> {
  const filePath = await open({
    multiple: false,
    filters: [{ name: "MatchY Session", extensions: ["matchy.json", "json"] }],
  });
  if (!filePath || typeof filePath !== "string") return null;

  const raw = await readTextFile(filePath);
  const session: SessionData = JSON.parse(raw);

  if (!session.pedigreeTgf || !session.haplotypesJson) {
    throw new Error("Invalid session file: missing pedigree or haplotypes.");
  }

  // Parse pedigree
  const pedigreeData = await invoke<{ individuals: unknown[]; relationships: unknown[] }>(
    "parse_tgf",
    { tgfContent: session.pedigreeTgf }
  );

  // Parse markers first (needed for parse_haplotypes_json)
  let markers: MarkerInfo[] = [];
  let selectedKitName = session.selectedKitName ?? null;
  let markerSetCsv = session.markerSetCsv ?? null;

  if (selectedKitName) {
    try {
      markers = await invoke<MarkerInfo[]>("load_kit", { name: selectedKitName });
    } catch {
      selectedKitName = null;
    }
  }
  if (markers.length === 0 && markerSetCsv) {
    try {
      markers = await invoke<MarkerInfo[]>("load_custom_csv", { csvContent: markerSetCsv });
    } catch {
      markerSetCsv = null;
    }
  }

  // Parse haplotypes (requires pedigree TGF + marker set)
  const haplotypesData = await invoke<HaplotypesParseResult>("parse_haplotypes_json", {
    jsonContent: session.haplotypesJson,
    pedigreeTgf: session.pedigreeTgf,
    markerSetName: selectedKitName,
    markerSetCsv: selectedKitName ? null : markerSetCsv,
  });

  return {
    pedigreeData,
    haplotypesData,
    haplotypesJson: session.haplotypesJson,
    pedigreeTgf: session.pedigreeTgf,
    markers,
    selectedKitName,
    markerSetCsv,
    suspect: session.suspect ?? null,
    exclude: session.exclude ?? [],
    params: session.params,
    simulationResult: session.simulationResult ?? null,
    simulationProgress: session.simulationProgress ?? [],
  };
}
