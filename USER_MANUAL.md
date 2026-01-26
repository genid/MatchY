# MatchY User Manual

**Version 1.0**
**For Forensic Practitioners and Scientists**

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Getting Started](#2-getting-started)
3. [Understanding the Science](#3-understanding-the-science)
4. [Using the GUI (Graphical Interface)](#4-using-the-gui-graphical-interface)
5. [Using the CLI (Command Line)](#5-using-the-cli-command-line)
6. [Helper Tools](#6-helper-tools)
7. [Interpreting Results](#7-interpreting-results)
8. [Best Practices](#8-best-practices)
9. [Frequently Asked Questions](#9-frequently-asked-questions)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Introduction

### What is MatchY?

MatchY is a specialized forensic software tool designed to calculate **match probabilities** for Y-chromosome STR (Y-STR) profiles within pedigrees. It uses Monte Carlo simulation methods to assess the likelihood that a DNA profile belongs to a specific individual within a family tree.

### Who Should Use MatchY?

- **Forensic Practitioners**: Case officers, forensic analysts, and investigators working with Y-STR evidence
- **Forensic Scientists**: Researchers and experts requiring statistical evaluation of Y-STR matchescases

### Key Features

✅ **Two Operating Modes**:
- **Standard Mode**: Calculate match probabilities for a known suspect/person of interest (PoI) within a pedigree
- **Trace Mode**: Identify the most likely donor of an unknown Y-STR profile (trace) from pedigree members

✅ **Flexible Execution**:
- Browser-based GUI (simple, visual)
- Command-line interface (faster, multi-threaded)

✅ **Built-in Tools**:
- Pedigree Builder (visual family tree creator)
- Haplotype Editor (Y-STR profile manager)

✅ **PDF Reports**:
- PDF reports with visualizations
- Match probability tables
- Quality metrics and diagnostics

---

## 2. Getting Started

### System Requirements

- **Operating System**: Windows, macOS, or Linux
- **Python**: Version 3.8 or higher
- **RAM**: Minimum 4GB (8GB+ recommended for large pedigrees)
- **Storage**: 500MB free space

### Installation

#### Docker image

The easiest way to install MatchY is via the Docker Image.

1. **Pull Docker image**
    ```bash
    docker pull dionzand/matchy:latest
    ```

2. **Run Docker with mounted results folder (recommended)**
    ```bash
    # Windows (PowerShell)
    docker run -p 8501:8501 -v C:\Users\YourName\MatchY-Results:/results dionzand/matchy:latest

    # macOS/Linux
    docker run -p 8501:8501 -v ~/MatchY-Results:/results dionzand/matchy:latest
    ```

    Or without mounted folder:
    ```bash
    docker run -p 8501:8501 dionzand/matchy:latest
    ```

    > **Note**: Mounting the results folder with `-v` allows you to access generated reports and output files directly on your local machine. Replace the path before the `:` with your desired local directory.

3. **Go to GUI**
    - Open `http://localhost:8501` (or specified port) in your favorite browser

#### Docker CLI image

For command-line usage, use the CLI Docker image:

1. **Pull Docker CLI image**
    ```bash
    docker pull dionzand/match-cli:latest
    ```

2. **Run with mounted folders**
    ```bash
    # Windows (PowerShell)
    docker run --rm `
      -v C:\Users\YourName\MatchY-Data:/app/data `
      -v C:\Users\YourName\MatchY-Results:/app/results `
      dionzand/match-cli:latest --config data/config.ini

    # macOS/Linux
    docker run --rm \
      -v ~/MatchY-Data:/app/data \
      -v ~/MatchY-Results:/app/results \
      dionzand/match-cli:latest --config data/config.ini
    ```

    This mounts:
    - Your local data folder (containing config files) to `/app/data` in the container
    - Your local results folder to `/app/results` in the container for easy access to outputs


#### Build from source

If you prefer to run MatchY locally, follow these steps:

1. **Install python**
2. **Clone repository**:
   ```bash
   git clone https://github.com/genid/MatchY.git
   ```
3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
4. **Run the application**:
  - CLI interface: ```python main.py```
  - Streamlit dashboard: ```streamlit run 1_🧬_Home.py```

Your default web browser will open automatically. If not, navigate to `http://localhost:8501`

> **Note**: The GUI runs locally on your computer. No data is sent to external servers.

---

## 3. Understanding the Science

### 3.1 What is Y-STR Analysis?

Y-chromosome Short Tandem Repeats (Y-STRs) are genetic markers inherited exclusively through the paternal line. Unlike autosomal DNA, Y-STRs:
- Pass from father to son virtually unchanged
- Undergo mutations at known rates (rapidly mutating (RM) Y-STRs have a mutation rate of >1% per generation)
- Allow tracing of paternal lineages over multiple generations

### 3.2 The Match Probability Problem

When a Y-STR profile from evidence matches (or partially matches) a reference profile, we need to answer:

> **"What is the probability that this profile belongs to a specific person within a family, given what we know about their relatives?"**

MatchY solves this by simulating thousands of possible scenarios through the family tree, accounting for:
- Mutation rates at each Y-STR marker
- The pedigree structure (who is related to whom)
- Known Y-STR profiles of family members
- The possibility of mutations occurring over generations

### 3.3 Weight of Evidence (Standard Mode)

#### What It Measures

The **weight of evidence** compares two hypotheses:

- **Hp (Prosecution Hypothesis)**: The suspect is the source of the evidence profile
- **Hd (Defense Hypothesis)**: Someone else (inside or outside the pedigree) is the source

#### Key Output Values

**1. Inside Match Probability Distribution**

This tells you: *"What is the probability that at least one other non-excluded, untyped individual within the pedigree shares the suspect's haplotype, given the pedigree structure and known haplotypes of typed relatives?"*

**Interpretation**:
The higher the Inside Match Probability, the higher the likelihood that another individual matches the trace profile, and thus the lower the evidential weight against the suspect.

> **📌 Forensic Note**: The match probability is highly dependent on the pedigree size and structure, as well as the number and position of typed individuals within the pedigree. Genotyping the father of the suspect (and indicating at least one allelic difference) is the best way to decrease the match probability and increase the evidential weight against the suspect.

**2. Outside Match Probability**

This represents: *"What is the probability that an unrelated male from (just) outside the current pedigree would match by coincidence?"*

> **📌 Forensic Note**: In some cases it might be that the outside match probability is actually higher than the inside match probability. Especially in the trace mode this means that the true trace donor is more likely to be outside the current pedigree.

### 3.4 Trace Mode (Donor Identification)

#### What It Does

**Trace mode** answers a different question:

> **"Given an unknown Y-STR profile (trace), which person in this pedigree is the most likely source?"**

#### How It Works

1. The trace profile is compared against all possible positions in the pedigree
2. For each individual, MatchY calculates: *"What's the probability the trace came from this person?"*
3. Results are **normalized** so all probabilities sum to 100%

#### Key Output

**Normalized Match Probabilities**

| Rank | Individual      | Normalized Probability |
|------|-----------------|------------------------|
| 1    | John Smith      | 45.2%                  |
| 2    | Robert Smith    | 28.7%                  |
| 3    | Outside Pedigree| 15.3%                  |
| 4    | Michael Smith   | 10.8%                  |

**Interpretation**:
- **Rank 1 (John Smith, 45.2%)**: Most likely source, but not certain
- **Outside Pedigree (15.3%)**: ~1 in 7 chance the donor isn't in the pedigree at all
- **Low top probability**: Suggests markers can't discriminate well between relatives


### 3.5 Mutation Modeling

#### Single-Step vs Two-Step Mutations

Y-STR mutations typically involve:
- **Single-step**: Gain or loss of one repeat unit (most common)
- **Two-step**: Gain or loss of two repeat units (less common)

**Two-Step Mutation Fraction** (default: 0.03)
- Controls the proportion of mutations that are two-step
- Based on published Y-STR mutation studies
- Affects calculations when evidence shows 2+ mutations from reference

> **📌 Technical Note**: Unless you have population-specific mutation data, use the default value (0.03). This is based on large-scale studies of Y-STR mutation patterns.

#### Mutation Rates by Marker

Each Y-STR marker has its own mutation rate. MatchY uses:
- **Built-in databases**: Pre-loaded marker sets with published rates
- **Custom rates**: Upload your own mutation rate files if needed

Common marker sets:
- **Yfiler™**: 17 markers
- **Yfiler™ Plus**: 25 markers
- **PowerPlex® Y23**: 22 markers
- **RMplex**: 30 markers

> **📌 Forensic Note**: The use of a kit with rapidly mutating (RM) Y-STR markers, such as RMplex, is adviced. Since the match probability is highly dependent on observed allelic differences within the pedigree, a marker set with elevated mutation rates increase the probability of observing such allelic differences. This does not, however, mean that the match probability will be higher when using a marker set with higher mutation rates - on the contrary actually: an allelic difference in a slowly mutating marker will result in a lower match probability (i.e. higher evidential value). The likelihood of observing a mutation in such a marker, however, is of course lower.

---

## 4. Using the GUI (Graphical Interface)

### 4.1 Overview

The GUI is designed for users who prefer a visual, point-and-click interface. It's ideal for:
- Single-case analysis
- Learning how MatchY works
- Cases not requiring high-performance computing

### 4.2 Main Workflow

```
Upload Files → Set Parameters → Run Simulation → Review Results → Download Report
```

### 4.3 Step-by-Step Guide

#### Step 1: Select Marker Set

**Location**: Left sidebar

1. Click the **"Select marker set"** dropdown
2. Choose your Y-STR panel (e.g., "Yfiler", "PowerPlex_Y23", "RMplex")

> **📌 What to Choose**: Select the marker set that matches your laboratory's Y-STR typing kit. If your kit isn't listed, you can create a custom marker set file (see Section 6.2).

#### Step 2: Upload Pedigree File

**Location**: Left sidebar

1. Click **"Upload pedigree file"**
2. Select a `.tgf` or `.ped` file

> **💡 Tip**: Don't have a pedigree file? Use the **Pedigree Builder** tool (see Section 6.1) to create one visually.

**File Formats**:

**TGF (Trivial Graph Format)** - Recommended for most users
```
1 GrandFather
2 Father
3 Son
4 Uncle
#
1 2
1 4
2 3
```
- Lines before `#`: Node definitions (ID, Name)
- Lines after `#`: Edges (Parent Child relationships)

**PED (Pedigree Format)** - Standard genetic format
```
Family1 1 0 0 1 0
Family1 2 1 0 1 0
Family1 3 2 0 1 0
```
- Columns: FamilyID, IndividualID, FatherID, MotherID, Sex, Phenotype
- Sex: 1=Male, 2=Female (MatchY only uses males)

> **⚠️ Important**:
> - Only **male** individuals are considered when opening .ped files
> - Use consistent names between pedigree and haplotype files
> - Ensure the pedigree forms a valid tree (no cycles, no disconnected branches)

#### Step 3: Upload Haplotypes File

**Location**: Left sidebar

1. Click **"Upload haplotypes file"**
2. Select a `.json` file containing Y-STR profiles

**JSON Format**:
```json
{
  "GrandFather": {
    "DYS391": "10",
    "DYS389I": "13",
    "DYS439": "12",
    "DYS389II": "29"
  },
  "Father": {
    "DYS391": "10",
    "DYS389I": "13",
    "DYS439": "11",
    "DYS389II": "29"
  }
}
```

**For Trace Mode**, include a special "TRACE" key:
```json
{
  "TRACE": {
    "DYS391": "10",
    "DYS389I": "13",
    "DYS439": "12",
    "DYS389II": "29"
  },
  "Father": {
    "DYS391": "10",
    ...
  }
}
```

**Multi-Copy Markers** (e.g., DYS385):
```json
{
  "Individual1": {
    "DYS385": "11;14"
  }
}
```
- Use semicolon `;` to separate multiple copies
- Order doesn't matter (MatchY will sort them)

**Intermediate Alleles** (e.g., 13.2):
```json
{
  "Individual1": {
    "DYS458": "17.2"
  }
}
```

> **💡 Tip**: Use the **Haplotype Editor** (see Section 6.2) to create and validate your JSON files.

#### Step 4: Upload Files

**Location**: Left sidebar

1. After selecting both files, click **"Upload files"**
2. Wait for confirmation messages:
   - ✅ "Pedigree file uploaded successfully"
   - ✅ "Haplotypes file uploaded successfully"

The pedigree visualization will appear in the main area.

#### Step 5: Configure Simulation

**Location**: Main area

##### A. Select Mode

**Standard Mode** (Default):
- Uncheck "Use trace mode"
- Select the **suspect** from the dropdown (the individual whose profile matches the evidence)

> **📌 Forensic Note**: Only a typed individual can be selected as suspect. This means that at least one typed individual should be present in the pedigree.

**Trace Mode**:
- Check "Use trace mode"
- Suspect selection becomes disabled (TRACE profile will be used instead)

> **📌 When to Use Each Mode**:
> - **Standard Mode**: You have a suspect in mind and want to assess match probability
> - **Trace Mode**: You have unknown DNA and want to identify which pedigree member it came from

##### B. Exclude Individuals (Optional)

**"Choose individuals to exclude from the simulation"**

Use this when certain individuals are:
- Deceased before the crime
- Proven to be elsewhere (alibi)
- Eliminated by other evidence

> **⚠️ Important**: Excluded individuals are not considered in the match probability calculation, but their (simulated) haplotypes will still be used for the simulation process.

> **📌 Forensic Note**: The match probability will usually decrease when excluding individuals (and therefore the weight of evidence against the suspect will increase). It is therefore adviced to use exclusion of individuals wisely.

##### C. Set Simulation Parameters

Click **"Set simulation parameters"** to expand advanced options:

**Two-step mutation fraction** (Default: 0.03)
- Range: 0.00 - 1.00
- Controls proportion of two-step mutations
- **When to change**: Can be set to 0 when no two-step mutations are allowed

**Batch length** (Default: 10,000)
- Number of iterations where results must remain stable
- Higher = more confident results, longer runtime
- **Recommended**: 10,000 for most cases, increase for more reliability

**Convergence criterion** (Default: 0.0500)
- Maximum allowed variation between convergence checks
- Lower = stricter convergence, more iterations
- **Recommended**: 0.0500 for most cases, increase for more reliability

**Skip inside pedigree probabilities**
- Check to disable inside pedigree calculations
- Only calculate outside match probability
- **Use when**: Only interested in outside match probability

**Skip outside pedigree probabilities**
- Check to disable outside match calculations
- Only calculate inside pedigree distribution
- **Use when**: Confident the source is within the pedigree.

> **💡 Tip**: Skipping the inside/outside pedigree probabilities can significantly reduce the runtime of the analysis. Only use any of the two when necessary.

**Set bias mode**
- Choose from any of the following: fixed bias, dynamic bias, adaptive bias mode.
- **Recommended**: Adaptive bias mode will converge the fastest, but might lead to a local optimum.

> **📌 What is bias?**: Bias is introduced in the simulation to steer the mutation process in a certain direction. Without bias, the mutation process is completely random, which might result in the simulation to follow irrelevant paths. Setting bias usually results in faster convergence, but might results in certain edge-cases to not be considered.

**Understanding the bias modes**

Each bias mode nudges the simulation in a given direction. When simulating the haplotype of an unknown individual, it takes into account the haplotypes of other relevant individuals. Take for example:

Suspect [10] -> unknown 1 [?]-> unknown 2 [?] -> known A [11]

Since we know that the suspect has allele [10] and known A has allele [11], it makes more sense that unknown 1 and 2 also are either allele [10] or [11], and not [9] or [8]. Setting a bias value nudges the simulation in the most likely direction.

- **Fixed bias**: When using a fixed bias value, the same bias value is applied to all edges in the pedigree.
- **Dynamic bias**: The bias value is determined by the number of allelic differences and the number of meioses between the unknown individual that is currently being simulated, and the closest known individual. The more allelic differences and/or the less meioses still left, the higher the bias value will be. For instance, when an allelic difference of [2] needs to be resolved in only 2 meioses, the bias value will be higher than when only [1] allelic difference needs to be resolve in 3 meioses.
- **Adaptive bias**: Automatically adjusts mutation bias based on convergence of the model.


> **📌 For Most Cases**: Use the default values unless you have specific reasons to adjust them. The defaults are based on extensive testing of the software.

##### D. Choose Execution Mode

**Execution Mode** options:

**1. Run in Browser (Single-threaded)**
- Runs in your web browser
- Slower (uses one CPU core)
- Real-time progress updates
- **Best for**: Small pedigrees, learning, quick checks

**2. Send CLI job from Browser (Multi-threaded)**
- Intermediate between GUI and CLI
- Set up your simulation via browser
- Send the simulation job to CLI
- **Best for**: Setting up large simulations visually

**3. Submit to CLI (Multi-threaded)**
- Runs as separate process
- Much faster (uses multiple CPU cores)
- Shows output after completion
- **Best for**: Large pedigrees, production casework

If choosing CLI mode:
- Adjust **"Number of threads"** slider (default: 4)
- More threads = faster runtime (up to your CPU core count)

> **💡 Performance Tip**: Runtime is dependent on a lot of factors: pedigree size and structure, number of typed and untyped individuals in the pedigree, number of typed markers and their mutation rates.

##### E. Name Your Simulation

**"Give this simulation a name"**: Enter a case-relevant identifier
- Examples: "Case2024-0123", "Smith_Family_Analysis", "Trace_Donor_Case456"
- This name appears in the PDF report filename

**"Your name"**: Enter the name of the user
- Appears in the PDF report header
- Important for chain of custody and quality assurance

#### Step 6: Run Simulation

1. Click **"Start simulation"** (blue button)
2. Progress indicators will appear:
   - Browser mode: Real-time progress bars and log messages
   - CLI mode: Output log after completion

> **⚠️ Do Not**: Close your browser or navigate away during simulation. You'll see a warning message while it runs.

#### Step 7: Review Results

##### Standard Mode Results

**Match Probabilities**:
- Table showing the inside and outside match probability

##### Trace Mode Results

**Most Likely Donor**:
- Name and probability displayed prominently
- Color-coded severity indicator

**Ranked Table**:
- All individuals sorted by match probability
- Includes "Outside Pedigree" option
- Normalized to sum to 100%

#### Step 8: Download Report

1. Click **"Download simulation report"** (or "Download trace donor identification report")
2. PDF file saves to your downloads folder
3. Filename format: `[SimulationName]_report.pdf` or `[SimulationName]_trace_report.pdf`

**Report Contents**:
- Summary statistics
- Match probability tables/charts
- Known haplotypes table
- Pedigree diagrams
- Simulation parameters used
- User name and date
- Quality metrics (convergence plots)

> **📌 For Case Files**: This PDF is designed to be court-ready. It includes all necessary information for peer review and legal proceedings.

### 4.4 GUI Tips and Tricks

**Tip 1: Use Example Files**
- Click "Download example pedigree file in TGF format" in the sidebar
- Great for learning the file formats

**Tip 2: Test with Simple Cases First**
- Start with 3-5 individuals
- Add complexity once you understand the workflow

**Tip 3: Check the Pedigree Visualization**
- After upload, verify the tree structure looks correct
- Known individuals are highlighted
- Suspect has a distinct color

**Tip 4: Save Your Configuration**
- Keep your JSON haplotype files in a secure location
- Document your marker set choice
- Maintain version control for complex cases

---

## 5. Using the CLI (Command Line)

### 5.1 Why Use the CLI?

The command-line interface is preferred for:
- **Batch processing**: Run multiple analyses automatically
- **High performance**: Full multi-threading support
- **Integration**: Embed MatchY in automated workflows
- **Reproducibility**: Script your exact analysis steps

### 5.2 Basic Usage

#### Local Installation

```bash
python main.py --config-path path/to/config.ini
```

#### Docker CLI

When using the Docker CLI image, mount your data and results folders:

```bash
# Windows (PowerShell)
docker run --rm `
  -v C:\Users\YourName\MatchY-Data:/app/data `
  -v C:\Users\YourName\MatchY-Results:/app/results `
  dionzand/match-cli:latest --config data/config.ini

# macOS/Linux
docker run --rm \
  -v ~/MatchY-Data:/app/data \
  -v ~/MatchY-Results:/app/results \
  dionzand/match-cli:latest --config data/config.ini
```

**Benefits of volume mounting:**
- 📁 Access your configuration files from your local machine
- 📊 Results are saved directly to your local results folder
- 🔄 No need to copy files in/out of the container

### 5.3 Creating a Configuration File

#### Method 1: Generate from GUI

The easiest way to create a config file:
1. Use the GUI to set up your analysis
2. Choose "Submit to CLI (Multi-threaded)"
3. MatchY creates `streamlit_generated_config.ini` in the results folder
4. Use this file as a template for future CLI runs

#### Method 2: Manual Creation

Create a file named `config.ini`:

```ini
[pedigree]
path = data/pedigrees/smith_family.tgf
known_haplotypes = data/haplotypes/smith_profiles.json
marker_set = data/marker_sets/yfiler_plus.csv
suspect = Father
simulation_name = Smith_Case_2024
user_name = Dr. Jane Doe
exclude_individuals = Uncle,Cousin
two_step_mutation_fraction = 0.03
batch_length = 1000
convergence_criterion = 0.01
number_of_threads = 4
results_path = results/
bias = 0.1
```

**Required Fields**:
- `path`: Path to pedigree file (.tgf or .ped)
- `known_haplotypes`: Path to JSON haplotypes file
- `marker_set`: Path to CSV marker set file (or use GUI to generate)
- `simulation_name`: Name for this analysis
- `user_name`: User name for report
- `results_path`: Directory to save results (created automatically)

**Optional Fields**:
- `suspect`: Individual name (omit for trace mode with TRACE in JSON)
- `exclude_individuals`: Comma-separated list of names to exclude
- `bias`: Mutation bias value (omit to disable)

**Marker Set CSV Format**:
```csv
marker,mutation_rate
DYS391,0.00251
DYS389I,0.00251
DYS439,0.00251
DYS389II,0.00251
DYS438,0.00074
```

### 5.4 Command-Line Options

**Full syntax**:
```bash
python main.py [OPTIONS]
```

**Available Options**:

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--config-path` | `-c` | Path to configuration file | `config.ini` |
| `--skip-inside` | `-i` | Skip inside pedigree calculations | Disabled |
| `--skip-outside` | `-o` | Skip outside pedigree calculations | Disabled |
| `--trace-mode` | `-t` | Enable trace mode | Disabled |
| `--adaptive-bias` | `-a` | Enable adaptive bias mode | Disabled |

**Examples**:

**Standard Analysis**:
```bash
# Local installation
python main.py -c cases/case123.ini

# Docker CLI
docker run --rm -v ~/data:/app/data -v ~/results:/app/results \
  dionzand/match-cli:latest --config data/case123.ini
```

**Trace Mode**:
```bash
# Local installation
python main.py -c cases/trace_analysis.ini --trace-mode

# Docker CLI
docker run --rm -v ~/data:/app/data -v ~/results:/app/results \
  dionzand/match-cli:latest --config data/trace_analysis.ini --trace-mode
```

**Skip Outside Probability** (faster, inside-only):
```bash
# Local installation
python main.py -c cases/case123.ini --skip-outside

# Docker CLI
docker run --rm -v ~/data:/app/data -v ~/results:/app/results \
  dionzand/match-cli:latest --config data/case123.ini --skip-outside
```

**High-Performance Mode** (adaptive bias, 8 threads):
```bash
# Local installation
python main.py -c cases/complex_pedigree.ini --adaptive-bias

# Docker CLI
docker run --rm -v ~/data:/app/data -v ~/results:/app/results \
  dionzand/match-cli:latest --config data/complex_pedigree.ini --adaptive-bias
```
(Set `number_of_threads = 8` in config.ini)

**Skip Inside Probability**:
```bash
# Local installation
python main.py -c cases/case123.ini --skip-inside

# Docker CLI
docker run --rm -v ~/data:/app/data -v ~/results:/app/results \
  dionzand/match-cli:latest --config data/case123.ini --skip-inside
```

---

## 6. Helper Tools

### 6.1 Pedigree Builder

**Location**: Navigate to "📊 Pedigree Builder" in the sidebar

The Pedigree Builder is a visual tool for creating family trees without writing TGF or PED files manually.

#### Step-by-Step Guide

**Step 1: Add Root Individual**
1. Enter name in "Individual name" field
2. Click "Add individual"
3. This becomes the top of your pedigree (typically the oldest known ancestor)

**Step 2: Add Children**
1. Select parent from "Select parent" dropdown
2. Enter child's name
3. Click "Add child to selected parent"
4. Repeat for all children

**Step 3: Build the Tree**
- Continue adding generations downward
- The tool automatically creates parent-child relationships
- You can have multiple children per parent

**Step 4: Visualize**
- The pedigree diagram updates automatically
- Review the structure for accuracy

**Step 5: Export**
1. Click "Download pedigree as TGF"
2. Save the `.tgf` file
3. Use this file in the main MatchY interface

#### Example Walkthrough

**Goal**: Create a 3-generation pedigree

```
Grandfather (root)
├── Father
│   ├── Son1
│   └── Son2
└── Uncle
    └── Cousin
```

**Actions**:
1. Add "Grandfather" (root)
2. Select "Grandfather", add child "Father"
3. Select "Grandfather", add child "Uncle"
4. Select "Father", add child "Son1"
5. Select "Father", add child "Son2"
6. Select "Uncle", add child "Cousin"
7. Download TGF

> **📌 Remember**: Only include males in the pedigree.

#### Common Issues

**Issue**: "Individual already exists"
- **Solution**: Each name must be unique. Use last names or suffixes (e.g., "John Sr.", "John Jr.")

**Issue**: Can't select a parent
- **Solution**: You must add at least one individual first. Start with the root.

### 6.2 Haplotype Editor

**Location**: Navigate to "🧬 Haplotype Editor" in the sidebar

The Haplotype Editor helps create and validate Y-STR profile JSON files.

#### Step-by-Step Guide

**Step 1: Select Marker Set**
- Choose the same marker set you'll use in the main analysis
- This loads the correct marker names

**Step 2: Add Individuals**
- Enter individual name (must match pedigree exactly)
- Click "Add individual"
- Repeat for all known profiles

**Step 3: Enter Y-STR Data**

For each individual:
1. Locate each marker row
2. Enter allele values in the input boxes

**Single-Copy Markers** (e.g., DYS391):
- Enter one value: `14`

**Multi-Copy Markers** (e.g., DYS385):
- Enter multiple values separated by semicolons: `11;14`
- Order doesn't matter

**Intermediate Alleles** (e.g., DYS458):
- Use decimal notation: `17.2`

**Step 4: Add TRACE Profile (Optional)**

For trace mode:
1. Click "Add TRACE profile"
2. A special "TRACE" individual appears
3. Enter alleles for the unknown profile

**Step 5: Validate**
- Click "Validate haplotypes"
- Checks for:
  - Missing markers
  - Invalid allele formats
  - Inconsistent multi-copy markers
  - Duplicate individual names

**Step 6: Export**
1. Click "Download haplotypes as JSON"
2. Save the `.json` file
3. Use in the main MatchY interface

#### Validation Rules

✅ **Valid Entries**:
```
14              (single copy)
11;14           (multi-copy, two copies)
11;13;15        (multi-copy, three copies)
17.2            (intermediate allele)
11.1;14.3       (multi-copy with intermediate alleles)
```

❌ **Invalid Entries**:
```
14.5.2          (multiple decimal points)
abc             (non-numeric)
11, 14          (comma instead of semicolon)
11; 14          (space after semicolon)
```

#### Best Practices

**1. Copy-Paste from Lab Reports**
- Use spreadsheet software to organize data first
- Copy values column by column into the editor
- Always validate afterward

**2. Use Consistent Naming**
- Individual names in Haplotype Editor must match Pedigree Builder exactly
- Case-sensitive: "John" ≠ "john"

**3. Handle Missing Markers**
- MatchY doesn't allow missing values
- If a marker could not be analyzed in one individual, that marker should be removed from the marker set

**4. Document Allele Calls**
- Keep original lab reports alongside JSON files
- Note any ambiguous calls or stutter peaks

**5. Trace Profile Quality**
- Only include unambiguous loci in TRACE
- Exclude markers with mixtures or artifacts

---

## 7. Interpreting Results

### 7.1 Standard Mode Interpretation

#### Average pedigree probability
```
The average pedigree probability gives the user an indication on how likely the current pedigree is that was provided. This value does not directly have an influence on the weight of evidence. It can however give an indication on how likely it is that the pedigree was correctly reconstructed. When the pedigree contains impossible haplotype combinations (e.g. allelic differences that would require three or more-step mutations) the average pedigree probability will be 0.
```

#### Extended average pedigree probability
```
Same as average pedigree probability, but for the extended pedigree (i.e. the outside match probability pedigree).
```

#### Inside Match probability

Example: `0.0002` (or 0.02%)

**What This Tells You**:

This is the probability that at least one other non-excluded, untyped individual in the pedigree matches the suspect's haplotype.

**Example Statement**:
```
"The inside match probability of 0.0002 indicates that there is a probability of approximately 0.02% that at least one other non-excluded, untyped individual in the pedigree matches the suspect's haplotype."
```

#### Outside Match probability

Example: `0.0005` (or 0.05%)

**What This Tells You**:

This is the probability that a random individual (just) outside the provided pedigree matches the suspect's haplotype.

**Example Statement**:
```
"The outside match probability of 0.0005 indicates that there is a probability of approximately 0.05% that a random individual (just) outside the provided pedigree matches the suspect's haplotype."
```

> **📌 Forensic Note**: Match probabilities depends on:
> - Pedigree structure
> - Number and location of typed individuals
> - Markers used (and their mutation rates)



#### The Convergence Plots

Your PDF report includes several convergence plots that show how the probabilities stabilize over iterations
- **Good convergence**: Line flattens out
- **Poor convergence**: Line continues oscillating

**What to Look For**:
✅ **Flat lines in the right half of the plot**: Good convergence, reliable results
❌ **Oscillating lines throughout**: Poor convergence, increase batch length
❌ **Sudden jumps**: When your pedigree contains large allelic differences, spikes might occur. Use adaptive bias

> **📌 Quality Assurance**: Always review convergence plots before finalizing reports. Poor convergence means unreliable probabilities.

### 7.2 Trace Mode Interpretation

#### The Normalized Probability Table

Example output:

| Rank | Individual       | Original Probability | Normalized Probability |
|------|------------------|----------------------|------------------------|
| 1    | John Smith       | 4.52E-15             | 45.2%                  |
| 2    | Robert Smith     | 2.87E-15             | 28.7%                  |
| 3    | Outside Pedigree | 1.53E-15             | 15.3%                  |
| 4    | Michael Smith    | 1.08E-15             | 10.8%                  |

**Understanding the Columns**:

**Original Probability**:
- Raw probability from simulation
- Very small numbers (scientific notation)
- Not directly interpretable alone

**Normalized Probability**:
- Scaled so all probabilities sum to 100%
- **This is what you report**
- Represents relative likelihood among options

**What This Tells You**:

**Scenario 1: High Inside Pedigree Probability**
```
Interpretation: Strong indication that an individual inside the pedigree is the most
likely source, though not conclusive proof.

Example Statement: "John Smith shows a normalized match probability
of 65.4%, indicating he is the most likely source of the trace DNA
among pedigree members. The next closest individual (Robert Smith, 22.1%)
is approximately 3 times less likely."
```


**Scenario 2: High Outside Pedigree Probability**
```
Interpretation: There's a substantial chance the true source is not
in the pedigree at all.

Example Statement: "While John Smith shows the highest match probability
among pedigree members (38.7%), the 'Outside Pedigree' option accounts
for 42.5% probability. This suggests the trace donor may be an unrelated
male or a paternal relative not included in the pedigree."
```


## 8. Best Practices

### 8.1 Quality Assurance

#### Pre-Analysis Checklist

Before running a simulation, verify:

- [ ] Pedigree visualizes correctly (no missing/extra individuals)
- [ ] All individual names match between pedigree and haplotype files
- [ ] Known individuals have complete Y-STR profiles (no missing markers)
- [ ] Marker set matches laboratory typing kit
- [ ] Suspect or TRACE is correctly identified
- [ ] Excluded individuals are appropriate for case context
- [ ] Simulation name is case-relevant and unique
- [ ] User name is correct for report

#### Post-Analysis Checklist

After simulation completes, check:

- [ ] All convergence plots show stable lines
- [ ] No warning messages in console output
- [ ] Results make sense given pedigree structure
- [ ] PDF report contains all expected sections
- [ ] Known haplotypes table matches input data
- [ ] Simulation parameters in report match intent

#### Validation Run

For critical cases, consider running the simulation twice:

   - The second run will have a different seed, thus resulting in slightly different results
   - Results should be nearly identical
   - Large differences indicate poor convergence

### 8.3 Performance Optimization

#### When to Adjust Parameters

**Increase number of cores** (e.g., 4 → 20):
- Increasing the number of cores is the easiest way to decrease runtime
- Runtime decreases roughly linear with increasing number of cores
- Consider moving to high performance computing (HPC) clusters for (very) large cases

**Decrease pedigree size**:
- Start with a smaller pedigree (only close relatives)
- See how the results change when more distant relatives are added to the pedigree

**Increase Batch Length** (e.g., 1000 → 5000):
- Convergence plots show oscillation
- Results vary between runs
- Complex pedigree (many unknown individuals)

**Decrease Model Validity** (e.g., 0.01 → 0.005):
- Need higher precision for publication
- Close probabilities require better resolution
- Willing to wait for longer runtime

**Enable Adaptive Bias**:
- When the standard analysis shows poor convergence
- When the standard analysis run for very long

**Skip Inside or Outside**:
- Testing/debugging (faster turnaround)
- Only need one type of probability
- Preliminary analysis before full run


## 10. Troubleshooting

### File Upload Issues

#### Problem: "Invalid pedigree file format"

**Cause**: File doesn't match TGF or PED format

**Solution**:
1. Open file in text editor or graph editor (such as yEd)
2. Check format matches examples in Section 4.3
3. Common issues:
   - Missing `#` separator in TGF files
   - Extra spaces or tabs
   - Windows vs Linux line endings (use `dos2unix` if needed)

**Fix for TGF**:
```
# Correct:
1 Father
2 Son
#
1 2

# Incorrect (no separator):
1 Father
2 Son
1 2
```

---

#### Problem: "Individual [Name] not found in pedigree"

**Cause**: Name mismatch between pedigree and haplotype files

**Solution**:
1. Compare names in both files (case-sensitive!)
2. Check for:
   - Extra spaces: "John " vs "John"
   - Different capitalization: "John" vs "john"
   - Special characters: "José" vs "Jose"
3. Edit files to make names exactly match

---

#### Problem: "Marker [Name] not found in marker set"

**Cause**: Haplotype file contains markers not in selected marker set

**Solution**:
1. Verify marker set matches your typing kit
2. Check marker name spelling (case-sensitive)
3. Common issues:
   - "DYS389-I" vs "DYS389I" (hyphen differences)
   - "DYSS391" vs "DYS391" (typo)

---

#### Problem: "Number of copies mismatch for marker [Name]"

**Cause**: Different individuals have different numbers of alleles for a multi-copy marker

**Solution**:
1. Check all entries for that marker in JSON file
2. Example problem:

    ***incorrect:***
   ```json
   "Individual1": {"DYS385": "11;14"},      // 2 copies
   "Individual2": {"DYS385": "11;14;15"}    // 3 copies ❌
   ```
    ***correct:***
   ```json
   "Individual1": {"DYS385": "11;14;14"},   // 3 copies
   "Individual2": {"DYS385": "11;14;15"}    // 3 copies
   ```

3. Verify laboratory data—this may indicate a data entry error


> **📌 What is the correct number of copies?**: The number of copies is not fixed for a given marker. However, it is necessary that within a pedigree the exact number of copies is shared between all pedigree members. Overlapping/stacked alleles should be entered twice (e.g., [11;14;14] in the example above).


## Appendix A: Glossary

**Allele**: A specific value (repeat count) at a Y-STR marker (e.g., "14" at DYS391).

**Haplotype**: The full set of alleles across all Y-STR markers for an individual.

**Intermediate Allele**: An allele with fractional repeats (e.g., 17.2), typically due to incomplete repeat units.

**Inside Match Probability**: Probability that at least one other non-excluded, untyped individual within the provided pedigree matches the suspect's haplotype.

**Mutation**: Change in repeat count passed from father to son.

**Outside Match Probability**: Probability that a random individual (just) outside the provided pedigree matches the suspect's haplotype.

**Pedigree**: Family tree showing paternal relationships.

**Batch Length**: Number of simulation iterations required to have stable results before stopping.

**Suspect**: In standard mode, the individual whose profile is being evaluated.

**TRACE**: In trace mode, the unknown Y-STR profile being compared against pedigree members.

**Two-Step Mutation**: A mutation involving gain or loss of 2 repeat units (less common than single-step).

**Y-STR**: Y-chromosome Short Tandem Repeat marker; inherited paternally.

---

## Appendix B: Quick Reference Cards

### GUI Quick Start

1. **Select marker set** (sidebar)
2. **Upload pedigree** (.tgf or .ped file)
3. **Upload haplotypes** (.json file)
4. Click **"Upload files"**
5. **Configure simulation**:
   - Select suspect OR enable trace mode
   - Set execution mode (browser/CLI)
   - Name simulation
6. Click **"Start simulation"**
7. **Download report** when complete

### CLI Quick Start

1. Create `config.ini` with required fields
2. Run: `python main.py -c config.ini`
3. Wait for completion
4. Find PDF report in `results/` folder

### Trace Mode Checklist

- [ ] JSON file contains `"TRACE"` key with haplotype
- [ ] "Use trace mode" checkbox enabled (GUI)
- [ ] OR `--trace-mode` flag used (CLI)
- [ ] Report shows normalized probabilities table
- [ ] TRACE column appears in known haplotypes table

### Troubleshooting Flowchart

```
Problem?
├─ Upload fails
│  ├─ Check file format (TGF/PED for pedigree, JSON for haplotypes)
│  └─ Verify names match between files
│
├─ Simulation won't start
│  ├─ Check suspect is selected (or trace mode enabled)
│  └─ Verify all individuals have haplotype data
│
├─ Simulation runs forever
│  ├─ Reduce batch length
│  └─ Use CLI with more threads
│
├─ Poor convergence
│  ├─ Increase batch length
│  └─ Enable adaptive bias
│
└─ Report issues
   ├─ Update WeasyPrint
   └─ Check results/ folder for plot images
```

---

## Appendix C: Example Files

### Example Pedigree (TGF)

**File: `example_pedigree.tgf`**
```
1 Grandfather
2 Father
3 Uncle
4 Son1
5 Son2
6 Cousin
#
1 2
1 3
2 4
2 5
3 6
```

### Example Haplotypes (JSON)

**File: `example_haplotypes.json`**
```json
{
  "Grandfather": {
    "DYS391": "10",
    "DYS389I": "13",
    "DYS439": "12",
    "DYS389II": "29",
    "DYS438": "12",
    "DYS385": "11;14"
  },
  "Father": {
    "DYS391": "10",
    "DYS389I": "13",
    "DYS439": "11",
    "DYS389II": "29",
    "DYS438": "12",
    "DYS385": "11;14"
  },
  "Son": {
    "DYS391": "10",
    "DYS389I": "13",
    "DYS439": "11",
    "DYS389II": "29",
    "DYS438": "12",
    "DYS385": "11;14"
  }
}
```

### Example Config (CLI)

**File: `example_config.ini`**
```ini
[pedigree]
path = examples/example_pedigree.tgf
known_haplotypes = examples/example_haplotypes.json
marker_set = data/marker_sets/yfiler.csv
suspect = Son
simulation_name = Example_Analysis
user_name = User Name
exclude_individuals =
two_step_mutation_fraction = 0.10
batch_length = 1000
convergence_criterion = 0.01
number_of_threads = 4
results_path = results/
bias = 0.1
```

---

## Feedback and Support

For questions, bug reports, or suggestions:
- **GitHub**: [Repository link]
- **Email**: [Support email]
- **Documentation**: [Online docs link]

---

**End of User Manual**

*MatchY - Forensic Y-STR Match Probability Analysis*
*Empowering forensic science with rigorous statistical methods*
