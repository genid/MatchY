import { useEffect, useState } from "react";
import { save } from "@tauri-apps/plugin-dialog";
import { writeTextFile } from "@tauri-apps/plugin-fs";

const STORAGE_KEY = "matchy_settings";

interface AppSettings {
  defaultThreads: number;
  defaultBatchLength: number;
  defaultConvergenceCriterion: number;
  defaultTwoStepFraction: number;
  simulationName: string;
  userName: string;
}

const DEFAULTS: AppSettings = {
  defaultThreads: 3,
  defaultBatchLength: 10000,
  defaultConvergenceCriterion: 0.02,
  defaultTwoStepFraction: 0.03,
  simulationName: "simulation",
  userName: "",
};

function loadSettings(): AppSettings {
  try {
    const stored = localStorage.getItem(STORAGE_KEY);
    if (stored) return { ...DEFAULTS, ...JSON.parse(stored) };
  } catch {}
  return { ...DEFAULTS };
}

export function useSettings(): AppSettings {
  return loadSettings();
}

export default function Settings() {
  const [settings, setSettings] = useState<AppSettings>(DEFAULTS);
  const [saved, setSaved] = useState(false);
  const [exportError, setExportError] = useState<string | null>(null);

  // Load on mount
  useEffect(() => {
    setSettings(loadSettings());
  }, []);

  const handleSave = () => {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(settings));
    setSaved(true);
    setTimeout(() => setSaved(false), 2000);
  };

  // Export as a TOML config file the CLI can use
  const handleExportToml = async () => {
    setExportError(null);
    const toml = `[simulation]
name = "${settings.simulationName}"
user = "${settings.userName}"
results_path = "./results"

[files]
pedigree = "./pedigree.tgf"
known_haplotypes = "./haplotypes.json"
marker_set = "RMplex"

[parameters]
two_step_mutation_fraction = ${settings.defaultTwoStepFraction}
batch_length = ${settings.defaultBatchLength}
convergence_criterion = ${settings.defaultConvergenceCriterion}
number_of_threads = ${settings.defaultThreads}

[mode]
# suspect = "IndividualName"
`;
    try {
      const filePath = await save({
        filters: [{ name: "TOML Config", extensions: ["toml"] }],
        defaultPath: `${settings.simulationName}.toml`,
      });
      if (!filePath) return;
      await writeTextFile(filePath, toml);
    } catch (e) {
      setExportError(String(e));
    }
  };

  const field = (
    label: string,
    key: keyof AppSettings,
    type: "text" | "number" = "number",
    extra?: React.InputHTMLAttributes<HTMLInputElement>
  ) => (
    <div>
      <label className="block text-sm text-gray-600 mb-1">{label}</label>
      <input
        type={type}
        className="w-full border rounded px-2 py-1.5 text-sm"
        value={settings[key] as string | number}
        onChange={(e) =>
          setSettings((s) => ({
            ...s,
            [key]: type === "number"
              ? (key === "defaultBatchLength" || key === "defaultThreads"
                  ? parseInt(e.target.value, 10)
                  : parseFloat(e.target.value))
              : e.target.value,
          }))
        }
        {...extra}
      />
    </div>
  );

  return (
    <div className="p-6 max-w-xl mx-auto space-y-6">
      <h1 className="text-xl font-bold text-gray-900">Settings</h1>

      <section className="bg-white rounded-lg border p-4 space-y-4">
        <h2 className="font-semibold text-gray-700">Identification</h2>
        <div className="grid grid-cols-2 gap-4">
          {field("Simulation name", "simulationName", "text")}
          {field("Analyst name", "userName", "text")}
        </div>
      </section>

      <section className="bg-white rounded-lg border p-4 space-y-4">
        <h2 className="font-semibold text-gray-700">Default Simulation Parameters</h2>
        <div className="grid grid-cols-2 gap-4">
          {field("Threads (1–3)", "defaultThreads", "number", { min: 1, max: 3 })}
          {field("Batch length", "defaultBatchLength", "number", { min: 100 })}
          {field("Convergence criterion", "defaultConvergenceCriterion", "number", { step: 0.001, min: 0.001, max: 0.2 })}
          {field("Two-step fraction", "defaultTwoStepFraction", "number", { step: 0.001, min: 0, max: 1 })}
        </div>
        <p className="text-xs text-gray-400">
          Threads controls how many of the 3 Monte Carlo models run in parallel.
          The convergence criterion is the maximum relative deviation between models
          before the result is accepted (lower = more precise but slower).
        </p>
      </section>

      {exportError && (
        <div className="rounded bg-red-50 border border-red-200 p-3 text-sm text-red-700">
          {exportError}
        </div>
      )}

      <div className="flex items-center gap-3 flex-wrap">
        <button
          onClick={handleSave}
          className="bg-blue-600 hover:bg-blue-700 text-white font-medium py-2 px-4 rounded text-sm"
        >
          Save Settings
        </button>
        <button
          onClick={handleExportToml}
          className="bg-white border hover:bg-gray-50 text-gray-700 font-medium py-2 px-4 rounded text-sm"
        >
          Export TOML config…
        </button>
        {saved && <span className="text-sm text-green-600 font-medium">Saved ✓</span>}
      </div>

      <section className="bg-gray-50 rounded-lg border p-4 text-xs text-gray-500 space-y-1">
        <p><strong>CLI usage:</strong></p>
        <p>
          After exporting a TOML config, run the CLI from the project directory:
        </p>
        <code className="block bg-gray-100 rounded px-2 py-1 font-mono text-gray-700 mt-1">
          matchy --config-path case.toml
        </code>
        <p className="mt-2">
          Settings are stored locally in your browser and persist between sessions.
        </p>
      </section>
    </div>
  );
}
