import { useEffect, useRef, useState } from "react";
import { Routes, Route, NavLink } from "react-router-dom";
import { getCurrentWindow } from "@tauri-apps/api/window";
import { getVersion } from "@tauri-apps/api/app";
import { invoke } from "@tauri-apps/api/core";
import { save } from "@tauri-apps/plugin-dialog";
import { writeTextFile } from "@tauri-apps/plugin-fs";
import { saveSession } from "./utils/session";
import Home from "./pages/Home";
import PedigreeBuilder from "./pages/PedigreeBuilder";
import HaplotypeEditor from "./pages/HaplotypeEditor";
import MarkerSets from "./pages/MarkerSets";
import Settings from "./pages/Settings";
import { GuidedTour, useTour } from "./components/GuidedTour";
import { useAppStore } from "./store/appStore";
import { useT } from "./i18n";

function App() {
  const tour = useTour();
  const darkMode = useAppStore((s) => s.darkMode);
  const t = useT();

  const [closeDialogOpen, setCloseDialogOpen] = useState(false);
  const [saving, setSaving] = useState(false);
  const [appVersion, setAppVersion] = useState<string>("");
  const destroyRef = useRef<(() => Promise<void>) | null>(null);

  useEffect(() => {
    getVersion().then(setAppVersion).catch(() => {});
  }, []);

  // Apply/remove `dark` class on <html> whenever darkMode changes
  useEffect(() => {
    document.documentElement.classList.toggle("dark", darkMode);
  }, [darkMode]);

  useEffect(() => {
    let cleanup: (() => void) | undefined;
    getCurrentWindow().onCloseRequested((event) => {
      const { pedigree, haplotypes, simulation } = useAppStore.getState();
      const hasData = !!(pedigree || haplotypes || simulation.result);
      if (hasData) {
        event.preventDefault();
        destroyRef.current = () => getCurrentWindow().destroy();
        setCloseDialogOpen(true);
      }
    }).then((unlisten) => { cleanup = unlisten; });
    return () => { cleanup?.(); };
  }, []);

  const handleQuit = async () => {
    setCloseDialogOpen(false);
    await destroyRef.current?.();
  };

  const handleSaveAndQuit = async () => {
    setSaving(true);
    try {
      const {
        pedigree, haplotypes, haplotypesJson, pedigreeTgf,
        selectedKitName, markerSetCsv, suspect, exclude,
        simulationName, userName, simParams,
        simulation,
      } = useAppStore.getState();

      if (pedigree) {
        const tgf = await invoke<string>("export_tgf", { data: pedigree });
        const filePath = await save({
          filters: [{ name: "TGF Pedigree", extensions: ["tgf"] }],
          defaultPath: "pedigree.tgf",
        });
        if (filePath) await writeTextFile(filePath, tgf);
      }

      if (haplotypes && haplotypesJson) {
        const filePath = await save({
          filters: [{ name: "JSON Haplotypes", extensions: ["json"] }],
          defaultPath: "haplotypes.json",
        });
        if (filePath) await writeTextFile(filePath, haplotypesJson);
      }

      if (pedigreeTgf && haplotypesJson) {
        await saveSession(
          {
            version: 1,
            pedigreeTgf,
            haplotypesJson,
            selectedKitName,
            markerSetCsv,
            suspect,
            exclude,
            params: { ...simParams, simulationName, userName },
            simulationResult: simulation.result,
            simulationProgress: simulation.progress,
          },
          simulationName || "session",
        );
      }
    } catch {
      // If save fails or user cancels, don't close
      setSaving(false);
      return;
    }
    setSaving(false);
    setCloseDialogOpen(false);
    await destroyRef.current?.();
  };

  const navItems = [
    { path: "/",           label: t("nav_run") },
    { path: "/pedigree",   label: t("nav_pedigree") },
    { path: "/markers",    label: t("nav_marker_sets") },
    { path: "/haplotypes", label: t("nav_haplotypes") },
    { path: "/settings",   label: t("nav_settings") },
  ];

  return (
    <div className="flex flex-col h-screen">
      {tour.show && <GuidedTour onClose={tour.close} />}

      {/* Close confirmation dialog */}
      {closeDialogOpen && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/50">
          <div className="bg-white rounded-lg shadow-xl border border-gray-200 p-6 max-w-sm w-full mx-4 space-y-4">
            <h2 className="text-base font-semibold text-gray-900">{t("app_close_title")}</h2>
            <p className="text-sm text-gray-600">{t("app_close_msg")}</p>
            <div className="flex flex-col gap-2 pt-1">
              <button
                onClick={handleSaveAndQuit}
                disabled={saving}
                className="w-full bg-blue-600 hover:bg-blue-700 disabled:opacity-60 text-white font-medium py-2 px-4 rounded text-sm"
              >
                {saving ? "…" : t("app_close_save_quit")}
              </button>
              <button
                onClick={handleQuit}
                disabled={saving}
                className="w-full bg-white hover:bg-gray-50 disabled:opacity-60 border border-gray-300 text-gray-700 font-medium py-2 px-4 rounded text-sm"
              >
                {t("app_close_btn")}
              </button>
              <button
                onClick={() => setCloseDialogOpen(false)}
                disabled={saving}
                className="w-full text-gray-400 hover:text-gray-600 disabled:opacity-60 py-1.5 px-4 rounded text-sm"
              >
                {t("ped_cancel")}
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Top navigation bar */}
      <header className="bg-white border-b border-gray-200 px-4 py-1 flex items-center gap-6 shadow-sm">
        <img src={darkMode ? "/logo_minimal_white.png" : "/logo.png"} alt="MatchY" className="h-20 w-auto" />
        <nav className="flex gap-1">
          {navItems.map(({ path, label }) => (
            <NavLink
              key={path}
              to={path}
              end={path === "/"}
              className={({ isActive }) =>
                `px-3 py-1.5 rounded text-sm font-medium transition-colors ${
                  isActive
                    ? "bg-blue-100 text-blue-700"
                    : "text-gray-600 hover:bg-gray-100 hover:text-gray-900"
                }`
              }
            >
              {label}
            </NavLink>
          ))}
        </nav>
        <div className="ml-auto flex items-center gap-3">
          {appVersion && (
            <span className="text-xs text-gray-400 select-none">v{appVersion}</span>
          )}
          <button
            onClick={tour.open}
            title={t("app_tour_tooltip")}
            className="text-xs text-gray-400 hover:text-gray-700 border border-gray-200 hover:border-gray-400 rounded-full w-6 h-6 flex items-center justify-center font-bold transition-colors"
          >
            ?
          </button>
        </div>
      </header>

      {/* Page content */}
      <main className="flex-1 overflow-auto bg-gray-50">
        <Routes>
          <Route path="/" element={<Home />} />
          <Route path="/pedigree" element={<PedigreeBuilder />} />
          <Route path="/haplotypes" element={<HaplotypeEditor />} />
          <Route path="/markers" element={<MarkerSets />} />
          <Route path="/settings" element={<Settings />} />
        </Routes>
      </main>
    </div>
  );
}

export default App;
