import { create } from "zustand";
import type {
  PedigreeData,
  HaplotypesParseResult,
  MarkerInfo,
  ProgressEvent,
  SimulationResponse,
  SimulationResult,
} from "../types/matchy";
import type { Locale } from "../i18n/types";
import { detectSystemLocale } from "../i18n";

const PREFS_KEY = "matchy_preferences";
const SETTINGS_KEY = "matchy_settings";

export interface SimParams {
  twoStepMutationFraction: number;
  batchLength: number;
  convergenceCriterion: number;
  bias: number | null;
  seed: number | null;
  numberOfThreads: number;
  skipInside: boolean;
  skipOutside: boolean;
  traceMode: boolean;
  adaptiveBias: boolean;
}

const DEFAULT_SIM_PARAMS: SimParams = {
  twoStepMutationFraction: 0.03,
  batchLength: 10000,
  convergenceCriterion: 0.02,
  bias: null,
  seed: null,
  numberOfThreads: 3,
  skipInside: false,
  skipOutside: false,
  traceMode: false,
  adaptiveBias: false,
};

function loadSimParamDefaults(): SimParams {
  try {
    const stored = localStorage.getItem(SETTINGS_KEY);
    if (stored) {
      const s = JSON.parse(stored);
      return {
        ...DEFAULT_SIM_PARAMS,
        twoStepMutationFraction: s.defaultTwoStepFraction ?? DEFAULT_SIM_PARAMS.twoStepMutationFraction,
        batchLength: s.defaultBatchLength ?? DEFAULT_SIM_PARAMS.batchLength,
        convergenceCriterion: s.defaultConvergenceCriterion ?? DEFAULT_SIM_PARAMS.convergenceCriterion,
        numberOfThreads: s.defaultThreads ?? DEFAULT_SIM_PARAMS.numberOfThreads,
      };
    }
  } catch {}
  return { ...DEFAULT_SIM_PARAMS };
}

interface Preferences {
  darkMode: boolean;
  locale: Locale;
}

function loadPreferences(): Preferences {
  try {
    const stored = localStorage.getItem(PREFS_KEY);
    if (stored) {
      const p = JSON.parse(stored) as Partial<Preferences>;
      return {
        darkMode: p.darkMode ?? false,
        locale: (p.locale as Locale) ?? detectSystemLocale(),
      };
    }
  } catch {}
  return { darkMode: false, locale: detectSystemLocale() };
}

function savePreferences(prefs: Preferences) {
  try {
    localStorage.setItem(PREFS_KEY, JSON.stringify(prefs));
  } catch {}
}

interface SimulationStatus {
  running: boolean;
  progress: ProgressEvent[];
  response: SimulationResponse | null;
  result: SimulationResult | null;
  error: string | null;
}

interface AppStore {
  // Appearance / locale (persisted, not reset by resetAll)
  darkMode: boolean;
  locale: Locale;
  setDarkMode: (v: boolean) => void;
  setLocale: (v: Locale) => void;

  // Pedigree
  pedigree: PedigreeData | null;
  pedigreeTgf: string;
  setPedigree: (data: PedigreeData, tgf: string) => void;
  clearPedigree: () => void;

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
  simulationName: string;
  setSimulationName: (name: string) => void;
  userName: string;
  setUserName: (name: string) => void;
  simParams: SimParams;
  setSimParams: (params: SimParams) => void;

  // Simulation status
  simulation: SimulationStatus;
  setSimulationRunning: (running: boolean) => void;
  addProgressEvent: (event: ProgressEvent) => void;
  setSimulationResult: (response: SimulationResponse | null) => void;
  setSimulationError: (error: string | null) => void;
  setSimulationProgress: (events: ProgressEvent[]) => void;
  resetSimulation: () => void;
  resetAll: () => void;
}

const initialPrefs = loadPreferences();

export const useAppStore = create<AppStore>((set, get) => ({
  darkMode: initialPrefs.darkMode,
  locale: initialPrefs.locale,
  setDarkMode: (v) => {
    set({ darkMode: v });
    savePreferences({ darkMode: v, locale: get().locale });
  },
  setLocale: (v) => {
    set({ locale: v });
    savePreferences({ darkMode: get().darkMode, locale: v });
  },

  pedigree: null,
  pedigreeTgf: "",
  setPedigree: (data, tgf) => set({ pedigree: data, pedigreeTgf: tgf }),
  clearPedigree: () => set({ pedigree: null, pedigreeTgf: "" }),

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
  simulationName: "",
  setSimulationName: (name) => set({ simulationName: name }),
  userName: "",
  setUserName: (name) => set({ userName: name }),
  simParams: loadSimParamDefaults(),
  setSimParams: (params) => set({ simParams: params }),

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
  setSimulationProgress: (events) =>
    set((s) => ({ simulation: { ...s.simulation, progress: events } })),
  resetSimulation: () =>
    set({ simulation: { running: false, progress: [], response: null, result: null, error: null } }),
  resetAll: () =>
    set({
      pedigree: null, pedigreeTgf: "",
      haplotypes: null, haplotypesJson: "",
      selectedKitName: null, markers: [], markerSetCsv: null,
      suspect: null, exclude: [],
      simulationName: "", userName: "",
      simParams: loadSimParamDefaults(),
      simulation: { running: false, progress: [], response: null, result: null, error: null },
      // darkMode and locale intentionally not reset
    }),
}));
