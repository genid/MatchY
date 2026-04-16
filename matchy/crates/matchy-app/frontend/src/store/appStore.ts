import { create } from "zustand";
import type {
  PedigreeData,
  HaplotypesParseResult,
  MarkerInfo,
  ProgressEvent,
  SimulationResponse,
  SimulationResult,
} from "../types/matchy";

interface SimulationStatus {
  running: boolean;
  progress: ProgressEvent[];
  response: SimulationResponse | null;
  result: SimulationResult | null;
  error: string | null;
}

interface AppStore {
  // Pedigree
  pedigree: PedigreeData | null;
  pedigreeTgf: string;
  setPedigree: (data: PedigreeData, tgf: string) => void;

  // Haplotypes
  haplotypes: HaplotypesParseResult | null;
  haplotypesJson: string;
  setHaplotypes: (data: HaplotypesParseResult, json: string) => void;

  // Marker set
  selectedKitName: string | null;
  markers: MarkerInfo[];
  markerSetCsv: string | null;
  setMarkerSet: (kitName: string | null, markers: MarkerInfo[], csv: string | null) => void;

  // Simulation config
  suspect: string | null;
  setSuspect: (name: string | null) => void;
  exclude: string[];
  setExclude: (names: string[]) => void;

  // Simulation status
  simulation: SimulationStatus;
  setSimulationRunning: (running: boolean) => void;
  addProgressEvent: (event: ProgressEvent) => void;
  setSimulationResult: (response: SimulationResponse | null) => void;
  setSimulationError: (error: string | null) => void;
  resetSimulation: () => void;
}

export const useAppStore = create<AppStore>((set) => ({
  pedigree: null,
  pedigreeTgf: "",
  setPedigree: (data, tgf) => set({ pedigree: data, pedigreeTgf: tgf }),

  haplotypes: null,
  haplotypesJson: "",
  setHaplotypes: (data, json) => set({ haplotypes: data, haplotypesJson: json }),

  selectedKitName: null,
  markers: [],
  markerSetCsv: null,
  setMarkerSet: (kitName, markers, csv) =>
    set({ selectedKitName: kitName, markers, markerSetCsv: csv }),

  suspect: null,
  setSuspect: (name) => set({ suspect: name }),
  exclude: [],
  setExclude: (names) => set({ exclude: names }),

  simulation: {
    running: false,
    progress: [],
    response: null,
    result: null,
    error: null,
  },
  setSimulationRunning: (running) =>
    set((s) => ({ simulation: { ...s.simulation, running } })),
  addProgressEvent: (event) =>
    set((s) => ({
      simulation: {
        ...s.simulation,
        progress: [...s.simulation.progress, event],
      },
    })),
  setSimulationResult: (response) =>
    set((s) => ({
      simulation: {
        ...s.simulation,
        response,
        result: response?.result ?? null,
      },
    })),
  setSimulationError: (error) =>
    set((s) => ({ simulation: { ...s.simulation, error } })),
  resetSimulation: () =>
    set((s) => ({
      simulation: { running: false, progress: [], response: null, result: null, error: null },
    })),
}));
