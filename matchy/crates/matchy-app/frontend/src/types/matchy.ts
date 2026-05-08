// TypeScript type mirrors of the Rust structs exposed via Tauri commands.

export interface IndividualData {
  id: string;
  name: string;
  haplotypeClass: "unknown" | "known" | "suspect" | "estimated" | "fixed" | "excluded";
  exclude: boolean;
}

export interface RelationshipData {
  parentId: string;
  childId: string;
}

export interface PedigreeData {
  individuals: IndividualData[];
  relationships: RelationshipData[];
}

export interface ValidationResult {
  isValid: boolean;
  isDag: boolean;
  isConnected: boolean;
  error: string | null;
}

export interface MarkerInfo {
  name: string;
  mutationRate: number;
  numberOfCopies: number | null;
}

export interface HaplotypesParseResult {
  /** individual_name → marker_name → allele_string */
  haplotypeTable: Record<string, Record<string, string>>;
  /** marker_name → allele_string for TRACE profile */
  traceHaplotype: Record<string, string> | null;
  /** Ordered list of marker names */
  markerNames: string[];
}

export interface SimulationRequest {
  pedigreeTgf: string;
  haplotypesJson: string;
  markerSetName: string | null;
  markerSetCsv: string | null;
  suspect: string | null;
  exclude: string[];
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
  autoBiasStrength: number | null;
  autoBiasMin: number | null;
  autoBiasMax: number | null;
  debugZeroProbSamples: number | null;
  debugZeroProbPath?: string | null;
}

/** Mirrors Rust's MatchProbabilities (snake_case keys — no rename_all) */
export interface MatchProbabilities {
  probabilities: Record<string, string>;       // match_count (as string key) → Decimal string
  average_pedigree_probability: string;
}

/** Mirrors Rust's SimulationParameters (snake_case keys, as serialised into SimulationResult). */
export interface SimulationParametersSnapshot {
  two_step_mutation_fraction: number;
  batch_length: number;
  convergence_criterion: number;
  bias: number | null;
  number_of_threads: number;
  suspect: string | null;
  exclude: string[];
  skip_inside: boolean;
  skip_outside: boolean;
  trace_mode: boolean;
  adaptive_bias: boolean;
  simulation_name: string;
  user_name: string;
  seed: number | null;
  auto_bias_strength: number;
  auto_bias_min: number;
  auto_bias_max: number;
  debug_zero_prob_samples: number | null;
  debug_zero_prob_path: string | null;
}

/** Mirrors Rust's SimulationResult (snake_case keys) */
export interface SimulationResult {
  /** MatchY version that produced this result (frozen at simulation time). */
  app_version?: string;
  /** Parameters frozen at simulation time — the authoritative source for reports. */
  parameters?: SimulationParametersSnapshot;
  inside_match_probabilities: MatchProbabilities | null;
  outside_match_probability: string | null;    // Decimal as string
  per_individual_probabilities: Record<string, string> | null;
  trials: number;
  converged: boolean;
}

export interface SimulationResponse {
  success: boolean;
  error: string | null;
  result: SimulationResult | null;
}

export interface ProgressEvent {
  trial: number;
  model: number;
  iteration: number;
  currentMean: string;
  stage: "pedigree_probability" | "extended_pedigree_probability" | "inside_match_probability" | "outside_match_probability";
  converged: boolean;
}
