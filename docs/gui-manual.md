# MatchY — Desktop Application User Manual

## Overview

MatchY is a forensic genetics tool for estimating Y-STR haplotype match probabilities using Monte Carlo simulation. The desktop application provides an interactive interface for building pedigrees, loading haplotype data, running simulations, and generating professional PDF-ready reports.

The application follows a five-step workflow: **Pedigree → Marker Sets → Haplotypes → Run → Report**.

---

## Installation

Download the appropriate installer from the [GitHub Releases](https://github.com/genid/MatchY/releases/latest) page:

| Platform | File | Notes |
|----------|------|-------|
| Windows | `MatchY-{version}-x64-windows-setup.exe` | Recommended — NSIS installer |
| Windows | `MatchY-{version}-x64-windows.msi` | MSI installer |
| Windows | `MatchY-{version}-x64-windows-portable.zip` | No installation required — unzip and run `matchy-app.exe` |
| Linux | `MatchY-{version}-x86_64-linux.AppImage` | Make executable, then run |
| macOS Intel | `MatchY-{version}-x86_64-macos.dmg` | See macOS note below |
| macOS ARM | `MatchY-{version}-aarch64-macos.dmg` | See macOS note below |

> **macOS:** Binaries are unsigned. Right-click → Open the first time, or run `xattr -dr com.apple.quarantine MatchY.app` in Terminal.

---

## Interface Overview

The application has five tabs accessible from the top navigation bar:

| Tab | Purpose |
|-----|---------|
| **Run** | Configure and execute simulations, view live results |
| **Pedigree** | Build or import the Y-chromosome lineage |
| **Marker Sets** | Select or build the marker panel |
| **Haplotypes** | Enter known Y-STR profiles |
| **Settings** | Application preferences and defaults |

The top-right corner shows the application version (`v1.0.0`) and a **?** button to launch the guided tour.

---

## Step 1 — Pedigree

### Building a pedigree

1. Go to the **Pedigree** tab.
2. Click **+ New pedigree** and type the name of the oldest ancestor (e.g. "Grandfather") in the dialog that appears.
3. Hover over or click a node to reveal the toolbar above it. The toolbar contains buttons to **Exclude/Include** the individual, mark them as **PoI** (only shown when a haplotype is loaded for that individual), and **✕** to remove them.
4. To add a **child**: drag the **bottom handle (●)** to an empty area of the canvas — a dialog appears asking for the new individual's name.
5. To add a **parent**: drag the **top handle (●)** to an empty area of the canvas — a dialog appears asking for the parent's name.
6. To **connect two existing nodes**: drag a handle (●) from one node and release it on the opposite handle of another node.
7. Double-click any label to rename an individual.
8. Select a node and press **Delete** (or click **✕** in the toolbar) to remove it. If removing the individual would also delete descendants, a confirmation dialog appears first.
9. Press **Ctrl+Z** to undo the last action.

Constraints:
- The pedigree must be a **directed acyclic graph (DAG)** — no cycles.
- Each individual may have at most one parent in the pedigree (Y-chromosome inheritance is paternal).

### Importing and exporting

Click **Import** to load a pedigree from:
- `.tgf` — Tab-separated graph format (two-column: parent ID → child ID)
- `.ped` — Standard pedigree format

Click **Export** to save the current pedigree as a `.tgf` file.

### Legend

| Colour | Meaning |
|--------|---------|
| Light gray | Unknown individual (no haplotype loaded) |
| Green | Known individual (haplotype loaded) |
| Pink | Person of interest (PoI) — set via the **PoI** toolbar button |
| Gray + dashed border | Excluded from calculation |

After running a simulation, click **Show in pedigree** on the Run tab to navigate to the Pedigree tab with match probabilities overlaid on each node. The overlay can also be toggled directly on the Pedigree tab via the **Show P(match)** button in the left panel.

---

## Step 2 — Marker Sets

Go to the **Marker Sets** tab to select which Y-STR markers are used in the simulation.

### Embedded kits

Click **Load kit** next to any built-in kit:

| Kit | Description |
|-----|-------------|
| RMplex | 30-marker Rapidly Mutating (RM) Y-STR panel |
| Yfiler plus | Thermo Fisher Yfiler Plus — 27 Y-STR markers |
| PowerPlex Y23 | Promega PowerPlex Y23 — 23 Y-STR markers |
| Combined | Union of RMplex, Yfiler Plus, and PowerPlex Y23 |

The loaded kit becomes the active marker set, shown in the **Active marker set** panel.

### Custom marker set

Two ways to create a custom set:

1. **Load CSV** — Upload a CSV file with columns `marker_name`, `mutation_rate`, and optionally `copies`.
2. **Custom Set Builder** — Click **Build custom set…** to open the marker pool. Search, select individual markers, override mutation rates (shown in orange when changed), and save the set for later reuse.

Saved sets persist between sessions and can be loaded, renamed, or deleted from the **Saved Sets** panel.

---

## Step 3 — Haplotypes

Go to the **Haplotypes** tab to enter the known Y-STR profiles.

### Adding individuals

- Click **Add from pedigree** to populate the table with individuals from the current pedigree.
- Click **➕ Add** and type a name to add an individual not in the pedigree.
- Click **Add TRACE Profile** to add a crime-scene profile for trace mode.

### Entering allele data

- Click any cell and type the allele value.
- Use **Tab** / **Shift+Tab** to move between cells horizontally.
- Use **Enter** / **Shift+Enter** to move between cells vertically.
- For multi-copy markers, separate allele values with a semicolon (e.g. `35;37`).

### Importing data

- **Import JSON** — Load a previously exported JSON file.
- **Paste from Excel…** — Copy a range from Excel or Google Sheets and paste it. Expected format: first row = individual names, first column = marker names.
- **Import CSV for {name}** — Import per-individual allele data from a CSV with columns `marker_name` and `alleles`.

### Excluding individuals

Toggle **exclude** on any individual's column header to exclude them from the match probability calculation. Excluded individuals are still shown in pedigree diagrams but are not considered as potential matching individuals.

---

## Step 4 — Run Simulation

Go to the **Run** tab.

### Data summary bar

Before running, confirm the data summary bar shows all three items loaded:
- **Pedigree** — number of individuals and relationships
- **Markers** — active marker set name and count
- **Haplotypes** — number of known profiles

### Mode

| Option | Description |
|--------|-------------|
| **Person of interest (PoI)** | Select the suspect from the dropdown. The simulation calculates the probability that another pedigree member matches the PoI's haplotype. |
| **Trace mode** | Automatically identifies the most likely donor from all pedigree members. Requires a TRACE profile loaded in the Haplotypes tab. |

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Threads** | max | Number of parallel CPU threads. Higher = faster on multi-core machines. |
| **Batch length** | 10,000 | Monte Carlo samples per batch. Each batch adds one data point to the convergence chart. |
| **Convergence criterion** | 0.02 (2%) | Maximum relative spread between the three independent model estimates before the result is accepted. Lower = more precise but slower. |
| **Two-step fraction** | 0.03 | Fraction of mutations that jump ±2 repeat units instead of ±1. |
| **Adaptive bias** | off | Differentiates the importance-sampling bias across the three ensemble models based on their past performance. Off by default; leave off unless you have a specific reason to enable it. |
| **Bias** | auto | Fixed importance-sampling bias factor. Leave blank to let the engine choose automatically. |
| **Seed** | 0 | Random seed. Leave blank to use the default seed (0), which gives reproducible results. Enter a different integer to obtain a different reproducible sequence. |

### Skip options

- **Skip inside** — Omit the within-pedigree match probability stage.
- **Skip outside** — Omit the outside-pedigree random match probability stage.

### Running

Click **Start simulation** or press **Ctrl+Enter**. The simulation runs continuously until convergence. Click **Cancel** to stop early.

### Live results

While running, the convergence charts show each model's running mean probability. The **inter-model variance** badge shows the current spread between models; when it drops below the criterion, the simulation converges and stops automatically.

### Results cards

After convergence, headline cards display:

| Card | Meaning |
|------|---------|
| **Pedigree probability** | P(observed pedigree state given all typed haplotypes) |
| **Inside-pedigree match** | P(at least one other non-excluded pedigree member has the same haplotype as the PoI) |
| **Outside-pedigree match** | P(a random male outside this pedigree has the same haplotype as the PoI) |
| **Pedigree odds** | P(PoI is the only match in pedigree) / P(at least one other also matches) — this is an odds ratio, not a likelihood ratio |
| **Average LR** | 1 / mean P(match) across all non-excluded unknown pedigree members (trace mode) |

> **Note:** Some cards are hidden in trace mode (they are not meaningful when no PoI is specified).

### Per-individual results

The **per-individual match probability table** lists each non-excluded unknown pedigree member with:

| Column | Meaning |
|--------|---------|
| **Match probability** | P(this individual has the same haplotype as the PoI), estimated by Monte Carlo simulation |
| **%** | Match probability expressed as a percentage |
| **LR** | Likelihood ratio = 1 / P(match). Quantifies how many times more likely an observed match is under the hypothesis that this individual is the true donor, compared to a random pedigree member coincidentally sharing the haplotype |

Individuals are sorted by match probability, highest first. The LR is the primary forensic metric: a higher LR indicates stronger evidence that the individual is related to the PoI.

Click **Copy results** to copy all numbers to the clipboard.

---

## Step 5 — Reports and Sessions

### Generating a report

Click **Generate Report** after a simulation completes. An interactive HTML report opens in your browser containing:
- Match probability summary table
- Narrative interpretation paragraphs
- Convergence charts
- Performance metrics
- Simulation parameters
- Pedigree diagrams
- Known haplotype tables

Click **⬇ Download PDF** in the report to save a PDF copy.

### Auto-save

Enable **Auto-Save Runs** in Settings and select a folder. Each run is automatically saved as an HTML report in a dated subfolder.

### Sessions

| Action | Description |
|--------|-------------|
| **Save session…** | Export the full workspace (pedigree, haplotypes, marker set, parameters) to a JSON file. |
| **Load session…** | Restore a previously saved workspace. |
| **New session** | Clear all data and start fresh. |

Sessions allow you to reproduce a simulation exactly or share a case with a colleague.

### Exporting TOML config

In Settings, click **Export TOML config…** to save the current parameters as a `.toml` file for use with the CLI.

---

## Settings

| Setting | Description |
|---------|-------------|
| **Dark mode** | Toggle dark/light theme. |
| **Language** | UI language: English, Nederlands, Deutsch, Español, Français, Português, 中文. |
| **Default simulation parameters** | Set default threads, batch length, convergence criterion, two-step fraction. |
| **Auto-Save Runs** | Folder for automatic report saving. |

Settings are stored locally and persist between sessions.

---

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| **Ctrl+Enter** | Start simulation (on Run tab) |
| **Ctrl+Z** | Undo (on Pedigree tab) |
| **Tab / Shift+Tab** | Move between haplotype cells |
| **Enter / Shift+Enter** | Move up/down in haplotype table |
| **Delete** | Remove selected pedigree node |
