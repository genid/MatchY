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
}

export interface SimulationResponse {
  success: boolean;
  error: string | null;
}

export interface ProgressEvent {
  trial: number;
  model: number;
  iteration: number;
  currentMean: string;
  stage: "pedigree_probability" | "inside_match_probability" | "outside_match_probability";
  converged: boolean;
}
