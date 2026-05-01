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
| Windows | `MatchY-{version}-x64-windows-portable.zip` | No installation required |
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
2. Click **+ New pedigree** and name the oldest ancestor (e.g. "Grandfather").
3. Hover over any node to reveal the toolbar. Use **+ Child** to add a descendant, or drag the **bottom handle (●)** onto another node to connect them.
4. To add a parent, drag the **top handle (●)** upward.
5. Double-click any label to rename an individual.
6. Select a node and press **Delete** (or click **✕**) to remove it and all its descendants.
7. Press **Ctrl+Z** to undo the last action.

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
| White / default | Unknown individual (no haplotype loaded) |
| Blue | Known individual (haplotype loaded) |
| Orange flag **⚑** | Person of interest (PoI) |
| Red | Excluded from calculation |
| Dashed border | Estimated haplotype |

After running a simulation, click **Show P(match)** on the Run tab to overlay match probabilities on the pedigree nodes.

---

## Step 2 — Marker Sets

Go to the **Marker Sets** tab to select which Y-STR markers are used in the simulation.

### Embedded kits

Click **Load kit** next to any built-in kit:

| Kit | Markers |
|-----|---------|
| RMplex | Rearranged markers panel |
| Yfiler plus | Promega Yfiler Plus |
| PowerPlex Y23 | PowerPlex Y23 |
| Combined | Merged superset |

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
- For multi-copy markers, separate allele values with a comma or space.

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
| **Adaptive bias** | on | Automatically tunes the importance-sampling bias for efficiency. Recommended. |
| **Bias** | auto | Fixed bias factor. Only relevant when Adaptive bias is off. |
| **Seed** | random | Fixed random seed for reproducible results. Leave blank for a different seed each run. |

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
| **Pedigree probability** | P(observed pedigree state) |
| **Inside-pedigree match** | P(at least one other pedigree member matches PoI) |
| **Outside-pedigree match** | P(random outside male matches PoI) |
| **Pedigree odds** | LR: PoI is unique in pedigree vs. at least one other matches |
| **Average LR** | Average likelihood ratio across pedigree members (trace mode) |

The **per-individual match probability table** lists each unknown individual with their match probability, percentage, and LR, sorted highest first.

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
