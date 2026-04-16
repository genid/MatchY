import { useCallback } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import { useAppStore } from "../store/appStore";
import type { ProgressEvent, SimulationRequest, SimulationResponse } from "../types/matchy";

export function useSimulation() {
  const {
    pedigreeTgf,
    haplotypesJson,
    selectedKitName,
    markerSetCsv,
    suspect,
    exclude,
    setSimulationRunning,
    addProgressEvent,
    setSimulationResult,
    setSimulationError,
    resetSimulation,
  } = useAppStore();

  const startSimulation = useCallback(
    async (params: Omit<SimulationRequest, "pedigreeTgf" | "haplotypesJson" | "markerSetName" | "markerSetCsv" | "suspect" | "exclude">) => {
      resetSimulation();
      setSimulationRunning(true);

      // Subscribe to progress events
      const unlisten = await listen<ProgressEvent>("simulation-progress", (event) => {
        addProgressEvent(event.payload);
      });

      try {
        const request: SimulationRequest = {
          pedigreeTgf,
          haplotypesJson,
          markerSetName: selectedKitName,
          markerSetCsv,
          suspect,
          exclude,
          ...params,
        };

        const response = await invoke<SimulationResponse>("run_simulation", {
          request,
        });
        setSimulationResult(response);

        if (!response.success && response.error) {
          setSimulationError(response.error);
        }
      } catch (err) {
        setSimulationError(String(err));
      } finally {
        setSimulationRunning(false);
        unlisten();
      }
    },
    [pedigreeTgf, haplotypesJson, selectedKitName, markerSetCsv, suspect, exclude]
  );

  const cancelSimulation = useCallback(async () => {
    try {
      await invoke("cancel_simulation");
    } catch (err) {
      console.error("Failed to cancel simulation:", err);
    }
  }, []);

  return { startSimulation, cancelSimulation };
}
