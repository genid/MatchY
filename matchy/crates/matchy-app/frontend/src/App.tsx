import { useEffect } from "react";
import { Routes, Route, NavLink } from "react-router-dom";
import { getCurrentWindow } from "@tauri-apps/api/window";
import { ask } from "@tauri-apps/plugin-dialog";
import Home from "./pages/Home";
import PedigreeBuilder from "./pages/PedigreeBuilder";
import HaplotypeEditor from "./pages/HaplotypeEditor";
import MarkerSets from "./pages/MarkerSets";
import Settings from "./pages/Settings";
import { GuidedTour, useTour } from "./components/GuidedTour";
import { useAppStore } from "./store/appStore";

const NAV_ITEMS = [
  { path: "/", label: "Run" },
  { path: "/pedigree", label: "Pedigree" },
  { path: "/markers", label: "Marker Sets" },
  { path: "/haplotypes", label: "Haplotypes" },
  { path: "/settings", label: "Settings" },
];

function App() {
  const tour = useTour();

  useEffect(() => {
    let cleanup: (() => void) | undefined;
    getCurrentWindow().onCloseRequested((event) => {
      const { pedigree, haplotypes, simulation } = useAppStore.getState();
      const hasData = !!(pedigree || haplotypes || simulation.result);
      if (hasData) {
        event.preventDefault();
        ask("Unsaved session data will be lost. Close MatchY?", {
          title: "Close MatchY",
          kind: "warning",
          okLabel: "Close",
          cancelLabel: "Cancel",
        }).then(async (confirmed) => {
          if (confirmed) await getCurrentWindow().destroy();
        });
      }
    }).then((unlisten) => { cleanup = unlisten; });
    return () => { cleanup?.(); };
  }, []);

  return (
    <div className="flex flex-col h-screen">
      {tour.show && <GuidedTour onClose={tour.close} />}

      {/* Top navigation bar */}
      <header className="bg-white border-b border-gray-200 px-4 py-1 flex items-center gap-6 shadow-sm">
        <img src="/logo.png" alt="MatchY" className="h-20 w-auto" />
        <nav className="flex gap-1">
          {NAV_ITEMS.map(({ path, label }) => (
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
        <div className="ml-auto">
          <button
            onClick={tour.open}
            title="Show guided tour"
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
