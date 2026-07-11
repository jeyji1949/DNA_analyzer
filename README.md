# 🧬 DNA Analyzer — BioSeq Lab

> Python application for bioinformatic analysis of bacterial DNA sequences.
> Full graphical interface, double-strand analysis, 7 types of analyses, interactive visualizations, and multi-format export.

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)
![Tkinter](https://img.shields.io/badge/GUI-Tkinter-1B7A4A?style=flat-square)
![Matplotlib](https://img.shields.io/badge/Plots-Matplotlib-orange?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Features](#-features)
- [What's New & Improvements](#-whats-new--improvements)
- [Project Structure](#-project-structure)
- [Installation](#-installation)
- [Running the App](#-running-the-app)
- [User Guide](#-user-guide)
- [Double-Strand Analysis](#-double-strand-analysis-new-feature)
- [Analysis Descriptions](#-analysis-descriptions)
- [Exporting Results](#-exporting-results)
- [Technical Architecture](#-technical-architecture)
- [Dependencies](#-dependencies)

---

## 🔭 Overview

DNA Analyzer is a desktop application built in **Python / Tkinter** for analyzing bacterial DNA sequences from either direct input or a FASTA file.

It covers the fundamental bioinformatics analyses with a modular architecture:

- **ORF detection** across all 6 reading frames starting from the 1st nucleotide (sense and antisense strands)
- **Identification of the most likely coding ORF** with a score out of 100 and protein translation
- **Simultaneous double-strand analysis** with a colored visual diagram and optimal result
- **Regulatory motif search**: promoters (-10 / -35), GC-rich Rho-independent terminators, Shine-Dalgarno sequences
- **Physicochemical analyses**: GC%, Tm, molecular mass
- **Interactive graphical visualizations** via Matplotlib
- **Export** to CSV, Excel, JSON, TXT, FASTA

![Main interface — ORF detection](screenshots/screenshot_orf_main.png)
*Main interface: ORF detection across the 6 reading frames, best candidate identified (Frame +1, ★ score)*

---

## ✨ Features

### Biological analyses
| Analysis | Description |
|---|---|
| 📋 Reading frames | Codon-by-codon segmentation across the 6 frames starting from the 1st nucleotide |
| 🏆 Coding ORF | Score /100: length, ATG, protein size, strand, stop — best candidate identified |
| ⚙️ Protein translation | Translation of the main ORF, colored by biochemical property |
| 📍 Promoters | Box-35 / Box-10 pairs with mismatch tolerance |
| 🎯 Shine-Dalgarno | AGGAGG consensus located 7–9 bp upstream of ATG |
| 🔚 Terminators | **GC-rich** stem-loop structures (stem ≥ 50% GC) + polyT |
| ✂️ Restriction sites | 15 common enzymes (EcoRI, BamHI, HindIII…) |
| 📊 Statistics | GC%, Wallace Tm, Nearest-Neighbor Tm, molecular mass |

### Interface & usability
- Graphical interface with a **biology-themed** design (emerald green palette)
- Real-time DNA double-helix animation in the title bar
- **Real-time input validation** (empty, too short, non-ATCG characters)
- **Threaded analysis** — non-blocking UI on long FASTA sequences
- FASTA file import (`.fasta`, `.fa`, `.txt`, `.seq`)
- One-click **reverse complement strand** calculation
- Configurable minimum ORF length threshold (default: 30 bp)

---

## 🆕 What's New & Improvements

### Critical fixes

**Fix A1 — Analysis button**
The `RUN ANALYSIS` button was missing its `command=self._run` argument and did not respond to clicks. Fixed in `gui/input_frame.py`.

**Fix A2 — Minimum ORF length threshold**
The default value (90 bp) was too restrictive for short sequences. Reduced to **30 bp** in `orf_finder.py` and `input_frame.py`.

### New features

**Improvement A3 — User input validation**
Three checks before any analysis:
- Empty sequence → inline red error message
- Sequence < 15 nucleotides → warning showing the length
- Non-ATCG characters → list of invalid characters

Validation also happens **in real time** while typing (red border if invalid characters are present). Messages disappear automatically after 5 seconds.

**Improvement A4 — ORF tabs redesign**
Two new distinct tabs replace the single ORF tab:
- **📋 Reading Frames**: all segments between stop codons across the 6 frames from the 1st nucleotide, with an indication of ATG presence
- **🏆 Coding ORF**: best candidate with a score /100, detailed justification, colored DNA sequence, and translated protein

**Improvement A5 — Biological filter for terminators**
Active Rho-independent terminators have a G-C-rich stem. Only terminators with a **stem GC ≥ 50%** are kept in the main analyses (others remain available as a fallback).

**Improvement A6 — Double-Strand Analysis Window** ⭐
See the dedicated section below.

---

## 📁 Project Structure

```
dna_analyzer/
│
├── main.py                    # Entry point — run with: python main.py
│
├── analysis/
│   ├── __init__.py
│   ├── orf_finder.py          # 6-frame ORFs + coding score [IMPROVED]
│   ├── motif_finder.py        # Promoters, SD, GC-rich terminators, restriction
│   └── statistics.py         # Composition, Tm, molecular mass, sliding GC
│
├── data/
│   ├── __init__.py
│   └── codon_table.py        # Codon table + restriction sites
│
├── gui/
│   ├── __init__.py
│   ├── app.py                # Main window + analysis orchestration [IMPROVED]
│   ├── input_frame.py        # Input, validation, action buttons [IMPROVED]
│   ├── results_frame.py      # 7 result tabs [IMPROVED]
│   ├── visualize_frame.py    # Matplotlib charts
│   └── dual_strand_window.py # Double-strand analysis window [NEW]
│
├── export/
│   ├── __init__.py
│   └── exporter.py           # CSV, Excel, JSON, TXT, FASTA export
│
└── requirements.txt          # Python dependencies
```

---

## 🛠 Installation

### Prerequisites

- Python **3.10 or higher**
- Up-to-date `pip`

### 1. Clone the repository

```bash
git clone https://github.com/your-username/dna_analyzer.git
cd dna_analyzer
```

### 2. Create a virtual environment (recommended)

```bash
python -m venv bioenv
source bioenv/bin/activate        # Linux / macOS
bioenv\Scripts\activate           # Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

#### Main dependencies

```
matplotlib>=3.7.0
openpyxl>=3.1.0
```

> **Note**: `tkinter` ships with standard Python.
> On Ubuntu/Debian, if tkinter is missing:
> ```bash
> sudo apt install python3-tk
> ```

---

## 🚀 Running the App

```bash
python main.py
```

The interface opens with a demo sequence pre-loaded, ready to be analyzed.

---

## 📖 User Guide

### Left panel — Input

```
┌─────────────────────────────────┐
│  DNA Sequence (5'→3')           │
│  ┌───────────────────────────┐  │
│  │ ATGCCCGGGTAA...           │  │
│  └───────────────────────────┘  │
│  300 bp  •  GC: 48.3%           │
│                                 │
│  [📁 Import FASTA] [↔ RC] [✕]   │
│                                 │
│  Analyses to run:               │
│  ☑ ORF Detection (6 frames)     │
│  ☑ Promoters -10 / -35          │
│  ☑ Shine-Dalgarno               │
│  ☑ Rho-independent terminators  │
│  ☑ Restriction sites            │
│                                 │
│  Min. ORF length (bp): [30]     │
│                                 │
│  [▶  RUN ANALYSIS]              │
│  [🧬  DOUBLE-STRAND DIAGRAM]    │
└─────────────────────────────────┘
```

| Element | Role |
|---|---|
| Text box | Paste the DNA sequence directly (ATCG only) — real-time validation |
| `📁 Import FASTA` | Load a FASTA file — `>` header lines are ignored |
| `↔ RC` | Replaces the sequence with its reverse complement strand |
| `✕` | Clears the input box |
| Checkboxes | Enable/disable each analysis individually |
| Min. ORF length | Threshold in bp — `30` for all ORFs, `90+` for biological candidates |
| `▶ RUN ANALYSIS` | Starts the analysis with input validation |
| `🧬 DOUBLE-STRAND DIAGRAM` | **[NEW]** Opens the window for visual analysis of both strands |

> **Automatic validation:** if the sequence is empty, too short, or contains non-ATCG characters, a red error message is shown and the analysis does not start.

---

### Right panel — Result tabs

#### 📋 Reading Frames *(new)*
Shows all segments between stop codons across the **6 reading frames**, starting from the 1st nucleotide, on both strands. Each row shows:
- Frame (+1 to +3 sense, -1 to -3 antisense) and strand
- Start and end position (1-based)
- Length in bp, number of amino acids
- Presence of an internal ATG (★)
- Terminal stop codon

The **best segment** (longest with ATG) is highlighted in green.

#### 🏆 Coding ORF *(new)*
Identifies the best coding candidate with:
- **Score /100** calculated from 7 biological criteria
- Detailed point-by-point justification
- Summary cards: frame, strand, position, length, ATG, stop, molecular mass
- ORF DNA sequence (by triplets)
- Translated protein colored (hydrophobic/polar/charged/other)
- Ranking of the top 10 candidates

#### 📍 Promoters
Table of detected Box-35 / Box-10 pairs. Quality levels: `Perfect consensus` > `Strong` > `Moderate` > `Weak`.

#### 🎯 Shine-Dalgarno
SD sites located with their associated ATG, spacing, and number of mismatches.

#### 🔚 Terminators
Stem-loop structures (GC-rich) + polyT. Columns: position, Arm1, Loop, Arm2, PolyT, stem length.

#### ✂️ Restriction
Enzymes present in the sequence, number of cuts, and 1-based positions.

#### 📊 Statistics
Length, GC%, AT%, G/C ratio, Wallace Tm, Nearest-Neighbor Tm, ss/ds molecular mass, colored visual composition bar.

---

## 🧬 Double-Strand Analysis *(new feature)*

Accessible via the **"🧬 DOUBLE-STRAND DIAGRAM"** button in the left panel.

### Biological principle
A double-stranded DNA molecule has two reading orientations:
- **Sense strand 5'→3' (+)** — the sequence as entered
- **Antisense strand 3'→5' (-)** — the reverse complement strand

A gene can be encoded on either strand. This window **simultaneously** analyzes both orientations and identifies the best candidate.

![Double-strand map — simultaneous analysis of both strands] ("2d137923-fcf8-4cee-995b-a141a202f439.png")
*Linear map of both strands: ORFs (green/orange), Box-35 promoters (blue), Box-10 (orange), GC-rich terminators. The sense strand (+) is identified as the best candidate (score 90/100).*

### Window tabs

#### 📊 Double-Strand Map
Linear representation of both strands on two stacked Matplotlib axes, in the same style as the main sequence map:
- 🟩 **Green/teal rectangles** — ORFs (gold for the best one)
- 🔵 **Blue lines** — Box-35 promoters
- 🔴 **Red lines** — Box-10 promoters
- 🟠 **Orange zones** — GC-rich terminators
- 🔺 **Purple triangles** — Shine-Dalgarno sites
- Matplotlib toolbar (zoom, PNG save)

#### 🔤 Strand + and Strand − Sequences
Each nucleotide is colored according to its biological function:

| Color | Element |
|---|---|
| Dark green `#A9DFBF` | Best-candidate ORF |
| Pale green `#D5F5E3` | Standard ORF |
| Pale blue `#AED6F1` | Box -35 promoter |
| Pale red `#F5B7B1` | Box -10 promoter |
| Pale orange `#FDEBD0` | Rho-independent terminator (GC-rich) |
| Pale purple `#E8DAEF` | Shine-Dalgarno site |
| **Bold green text** | ATG codon |
| **Bold red text** | Stop codon |

Each 60-base line also shows position annotations on the right (`─35@pos`, `TERM@pos(GC:62%)`).

#### ★ Optimal Result

![Optimal result — double-strand analysis](screenshots/screenshot_optimal_result.png)
*Optimal result: sense strand (+) selected with a score of 90/100. Comparison table of strand (+) vs strand (−), translated protein of 47 amino acids colored by biochemical property.*

- **Gold banner** with the 12 key metrics of the winning strand
- **Mini-map** of the optimal strand with legend
- **Translated protein** colored by biochemical property
- **Comparison table** of strand (+) vs strand (-) with winning cells highlighted in green

---

### Visualization Tab (main panel)

![Sequence map — linear view](screenshots/screenshot_sequence_map.png)
*Linear view of the sequence (300 bp): ORF detected in green (positions 202–234), Box-35 promoters (blue lines) and Box-10 (orange lines) distributed across the sequence.*

Use the dropdown menu to choose among 4 charts:

| Chart | Description |
|---|---|
| Sequence map | Linear view with ORFs, promoters, SD, and terminators |
| GC% (sliding window) | GC-content curve over a 100 bp window with an average GC% line |
| Nucleotide composition | Histogram + pie chart of A/T/G/C |
| ORF lengths | Horizontal bars for the first 20 ORFs, colored by strand |

---

## 💾 Exporting Results

Accessible via **File** in the menu bar, after running an analysis.

| Format | Content |
|---|---|
| 📊 Excel (`.xlsx`) | 6 tabs: ORFs, Promoters, SD, Terminators, Restriction, Statistics |
| 📄 CSV (`.csv`) | All results in a tab-delimited text file |
| 🗄️ JSON (`.json`) | Full structured export, ideal for programmatic processing |
| 📝 TXT Report (`.txt`) | Readable, formatted report ready to be shared |
| 🧬 FASTA (`.fasta`) | Analyzed sequence in standard FASTA format |

---

## 🏗 Technical Architecture

### Analysis flow

```
Sequence input
      │
      ▼
input_frame.py ──► Validation (empty / too short / non-ATCG)
      │
      ▼
app.py (_run_analysis)
      │
      ├──► find_reading_frames()  →  find_best_coding_orf()   [6 frames + score]
      ├──► find_orfs()                                         [classic ATG→Stop]
      ├──► find_promoters()                                    [Box-35 / Box-10]
      ├──► find_shine_dalgarno()                               [AGGAGG]
      ├──► find_terminators() + GC ≥ 50% filter                [stem-loop]
      ├──► find_restriction_sites()                            [15 enzymes]
      └──► nucleotide_composition() + Tm + MW                 [statistics]
               │
               ▼
      results_frame.display()    visualize_frame.draw()
```

### Coding potential score

```python
score = 0
score += min(40, orf_length / seq_length * 200)       # relative length
score += 20 if has_atg else 0                         # start codon
score += 20 if num_aa >= 50 else 10 if num_aa >= 30 else 0  # size
score += 5 if strand == '+' else 0                    # sense strand
score += 10 if stop_codon present else 0              # complete ORF
score += 5 if '*' not in protein else 0               # no internal stop
```

### Threading

```python
import threading

def _worker():
    results = run_all_analyses(seq)
    self.root.after(0, lambda: self._finish_analysis(results))

threading.Thread(target=_worker, daemon=True).start()
```

---

## 📦 Dependencies

```
matplotlib>=3.7.0       # Charts embedded in Tkinter
openpyxl>=3.1.0         # Excel export (.xlsx)
```

> The `tkinter`, `threading`, `csv`, `json`, `math`, and `collections` modules are included in Python's standard library.

Installation:
```bash
pip install matplotlib openpyxl
```

---

## 🧪 Test Sequence

The following sequence (300 bp) is used as a demo sequence. It contains a main ORF on the **antisense strand** (score 90/100), promoters, and GC-rich terminators:

```
CATTTCTCTTAAGATTTATTCTATCTTAACACAACAACTTTTAATAAAAGATATGTAGAT
TACAATTTAAATAGATTGTAATATTTGTAACACTAACATTAATATAGTTGTTATTTTTGT
TACATAAACCACTAATAACTCATAATCTTTTAAAACTTATATTTGAGATAACATCAACTT
TACATTACAAGTTATAAAACAAAAGAAGTGGGACACAGAATTCGTCTTGAACACTGTGTC
CCACCTCGTCCCCAAAACTTGCTCTGTCCGTAGAAAAATAAAAAGGGGCCCCCTTTGTTG
```

**Expected result:**

| Criterion | Strand + | Strand − |
|---|---|---|
| Coding score | 76/100 | **90/100 ★** |
| Best ORF | 54 bp / 18 aa | **141 bp / 47 aa** |
| ATG codon | ✓ | ✓ |
| GC-rich terminators | 2 | 2 |

---

## 👥 Authors

me xoxo

---

*DNA Analyzer — BioSeq Lab | Python / Tkinter / Matplotlib*
