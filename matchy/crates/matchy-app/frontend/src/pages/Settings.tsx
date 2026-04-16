import { useState } from "react";

const DEFAULT_SETTINGS = {
  defaultThreads: 4,
  defaultBatchLength: 10000,
  defaultConvergenceCriterion: 0.02,
  defaultTwoStepFraction: 0.03,
};

export default function Settings() {
  const [settings, setSettings] = useState(DEFAULT_SETTINGS);
  const [saved, setSaved] = useState(false);

  const handleSave = () => {
    // TODO: persist to app_config_dir()/settings.toml via Tauri fs plugin
    setSaved(true);
    setTimeout(() => setSaved(false), 2000);
  };

  return (
    <div className="p-6 max-w-xl mx-auto space-y-6">
      <h1 className="text-xl font-bold text-gray-900">Settings</h1>

      <section className="bg-white rounded-lg border p-4 space-y-4">
        <h2 className="font-semibold text-gray-700">Default Simulation Parameters</h2>
        <div className="grid grid-cols-2 gap-4 text-sm">
          <div>
            <label className="block text-gray-600 mb-1">Threads</label>
            <input
              type="number"
              min="1"
              max="32"
              className="w-full border rounded px-2 py-1.5"
              value={settings.defaultThreads}
              onChange={(e) =>
                setSettings((s) => ({ ...s, defaultThreads: parseInt(e.target.value, 10) }))
              }
            />
          </div>
          <div>
            <label className="block text-gray-600 mb-1">Batch length</label>
            <input
              type="number"
              min="100"
              className="w-full border rounded px-2 py-1.5"
              value={settings.defaultBatchLength}
              onChange={(e) =>
                setSettings((s) => ({ ...s, defaultBatchLength: parseInt(e.target.value, 10) }))
              }
            />
          </div>
          <div>
            <label className="block text-gray-600 mb-1">Convergence criterion</label>
            <input
              type="number"
              step="0.001"
              min="0"
              className="w-full border rounded px-2 py-1.5"
              value={settings.defaultConvergenceCriterion}
              onChange={(e) =>
                setSettings((s) => ({ ...s, defaultConvergenceCriterion: parseFloat(e.target.value) }))
              }
            />
          </div>
          <div>
            <label className="block text-gray-600 mb-1">Two-step fraction</label>
            <input
              type="number"
              step="0.001"
              min="0"
              max="1"
              className="w-full border rounded px-2 py-1.5"
              value={settings.defaultTwoStepFraction}
              onChange={(e) =>
                setSettings((s) => ({ ...s, defaultTwoStepFraction: parseFloat(e.target.value) }))
              }
            />
          </div>
        </div>
      </section>

      <div className="flex items-center gap-3">
        <button
          onClick={handleSave}
          className="bg-blue-600 hover:bg-blue-700 text-white font-medium py-2 px-4 rounded"
        >
          Save Settings
        </button>
        {saved && <span className="text-sm text-green-600">Saved!</span>}
      </div>
    </div>
  );
}
