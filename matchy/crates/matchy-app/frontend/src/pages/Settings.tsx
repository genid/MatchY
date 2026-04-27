import { useEffect, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { save, open } from "@tauri-apps/plugin-dialog";
import { writeTextFile } from "@tauri-apps/plugin-fs";
import { useAppStore } from "../store/appStore";
import { useT } from "../i18n";
import type { Locale } from "../i18n/types";

const STORAGE_KEY = "matchy_settings";

interface AppSettings {
  defaultThreads: number;
  defaultBatchLength: number;
  defaultConvergenceCriterion: number;
  defaultTwoStepFraction: number;
  runsFolder: string;
}

const DEFAULTS: AppSettings = {
  defaultThreads: 3,
  defaultBatchLength: 10000,
  defaultConvergenceCriterion: 0.02,
  defaultTwoStepFraction: 0.03,
  runsFolder: "",
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

const LOCALES: { value: Locale; labelKey: "settings_lang_en" | "settings_lang_nl" | "settings_lang_de" | "settings_lang_es" | "settings_lang_fr" | "settings_lang_pt" }[] = [
  { value: "en", labelKey: "settings_lang_en" },
  { value: "nl", labelKey: "settings_lang_nl" },
  { value: "de", labelKey: "settings_lang_de" },
  { value: "es", labelKey: "settings_lang_es" },
  { value: "fr", labelKey: "settings_lang_fr" },
  { value: "pt", labelKey: "settings_lang_pt" },
];

export default function Settings() {
  const [settings, setSettings] = useState<AppSettings>(DEFAULTS);
  const [saved, setSaved] = useState(false);
  const [exportError, setExportError] = useState<string | null>(null);
  const [cpuCount, setCpuCount] = useState<number>(256);

  const { darkMode, setDarkMode, locale, setLocale } = useAppStore();
  const t = useT();

  useEffect(() => {
    setSettings(loadSettings());
    invoke<number>("get_cpu_count").then(setCpuCount).catch(() => {});
  }, []);

  const handleSave = () => {
    localStorage.setItem(STORAGE_KEY, JSON.stringify(settings));
    setSaved(true);
    setTimeout(() => setSaved(false), 2000);
  };

  const handlePickRunsFolder = async () => {
    const result = await open({ directory: true, multiple: false });
    if (result && typeof result === "string") {
      setSettings((s) => ({ ...s, runsFolder: result }));
    }
  };

  const handleExportToml = async () => {
    setExportError(null);
    const toml = `[simulation]
name = "simulation"
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
        defaultPath: "matchy.toml",
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
      <h1 className="text-xl font-bold text-gray-900">{t("settings_title")}</h1>

      {/* ── Appearance ── */}
      <section className="bg-white rounded-lg border p-4 space-y-4">
        <h2 className="font-semibold text-gray-700">{t("settings_appearance")}</h2>
        <div className="flex items-center justify-between">
          <span className="text-sm text-gray-700">{t("settings_dark_mode")}</span>
          <button
            role="switch"
            aria-checked={darkMode}
            onClick={() => setDarkMode(!darkMode)}
            className={`relative inline-flex h-6 w-11 items-center rounded-full transition-colors focus:outline-none ${
              darkMode ? "bg-blue-600" : "bg-gray-300"
            }`}
          >
            <span
              className={`inline-block h-4 w-4 transform rounded-full bg-white shadow transition-transform ${
                darkMode ? "translate-x-6" : "translate-x-1"
              }`}
            />
          </button>
        </div>
        <div className="flex items-center justify-between">
          <span className="text-sm text-gray-700">{t("settings_language")}</span>
          <select
            value={locale}
            onChange={(e) => setLocale(e.target.value as Locale)}
            className="border rounded px-2 py-1 text-sm"
          >
            {LOCALES.map(({ value, labelKey }) => (
              <option key={value} value={value}>{t(labelKey)}</option>
            ))}
          </select>
        </div>
      </section>

      {/* ── Simulation defaults ── */}
      <section className="bg-white rounded-lg border p-4 space-y-4">
        <h2 className="font-semibold text-gray-700">{t("settings_simulation_defaults")}</h2>
        <div className="grid grid-cols-2 gap-4">
          {field(`${t("settings_threads")} (${t("settings_threads_max")} ${cpuCount})`, "defaultThreads", "number", { min: 1, max: cpuCount })}
          {field(t("settings_batch_length"), "defaultBatchLength", "number", { min: 100 })}
          {field(t("settings_convergence_criterion"), "defaultConvergenceCriterion", "number", { step: 0.001, min: 0.001, max: 0.2 })}
          {field(t("settings_two_step_fraction"), "defaultTwoStepFraction", "number", { step: 0.001, min: 0, max: 1 })}
        </div>
        <p className="text-xs text-gray-400">{t("settings_convergence_desc")}</p>
      </section>

      {/* ── Auto-save ── */}
      <section className="bg-white rounded-lg border p-4 space-y-4">
        <h2 className="font-semibold text-gray-700">{t("settings_auto_save")}</h2>
        <p className="text-xs text-gray-400">{t("settings_auto_save_desc")}</p>
        <div className="flex items-center gap-2">
          <input
            type="text"
            readOnly
            placeholder={t("settings_no_folder")}
            className="flex-1 border rounded px-2 py-1.5 text-sm bg-gray-50 text-gray-700"
            value={settings.runsFolder}
          />
          <button
            onClick={handlePickRunsFolder}
            className="border rounded px-3 py-1.5 text-sm bg-white hover:bg-gray-50 text-gray-700 whitespace-nowrap"
          >
            {t("settings_choose_folder")}
          </button>
          {settings.runsFolder && (
            <button
              onClick={() => setSettings((s) => ({ ...s, runsFolder: "" }))}
              className="text-gray-400 hover:text-red-500 text-sm px-1"
              title="Clear"
            >
              ✕
            </button>
          )}
        </div>
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
          {t("settings_save")}
        </button>
        <button
          onClick={handleExportToml}
          className="bg-white border hover:bg-gray-50 text-gray-700 font-medium py-2 px-4 rounded text-sm"
        >
          {t("settings_export_toml")}
        </button>
        {saved && <span className="text-sm text-green-600 font-medium">{t("settings_saved")}</span>}
      </div>

      <section className="bg-gray-50 rounded-lg border p-4 text-xs text-gray-500 space-y-1">
        <p><strong>{t("settings_cli_usage")}</strong></p>
        <p>{t("settings_cli_desc")}</p>
        <code className="block bg-gray-100 rounded px-2 py-1 font-mono text-gray-700 mt-1">
          matchy --config-path case.toml
        </code>
        <p className="mt-2">{t("settings_local_storage_note")}</p>
      </section>
    </div>
  );
}
