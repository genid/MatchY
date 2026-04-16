import { Routes, Route, NavLink } from "react-router-dom";
import Home from "./pages/Home";
import PedigreeBuilder from "./pages/PedigreeBuilder";
import HaplotypeEditor from "./pages/HaplotypeEditor";
import MarkerSets from "./pages/MarkerSets";
import Settings from "./pages/Settings";

const NAV_ITEMS = [
  { path: "/", label: "Run" },
  { path: "/pedigree", label: "Pedigree" },
  { path: "/haplotypes", label: "Haplotypes" },
  { path: "/markers", label: "Marker Sets" },
  { path: "/settings", label: "Settings" },
];

function App() {
  return (
    <div className="flex flex-col h-screen">
      {/* Top navigation bar */}
      <header className="bg-white border-b border-gray-200 px-4 py-2 flex items-center gap-6 shadow-sm">
        <div className="flex items-center gap-2">
          <span className="text-lg font-bold text-blue-700">🧬 MatchY</span>
        </div>
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
