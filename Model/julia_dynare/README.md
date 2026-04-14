# NK-IOSOE Julia/Dynare — Free MATLAB Replacement

This folder contains a complete Julia translation of `main_SOE_gap.m` for the
**12-sector New Keynesian Input-Output Small Open Economy (NK-IOSOE)** model
of the Chilean economy.

Run the full simulation pipeline without MATLAB by installing two free tools:
**Julia** and **Octave + Dynare**.

---

## Quick Start (3 steps)

### Step 1 — Install Julia

Download from <https://julialang.org/downloads/> (≥ 1.9).

### Step 2 — Install Octave and Dynare

**Octave** is a free, drop-in MATLAB replacement:

| OS | Command |
|----|---------|
| Linux (Ubuntu/Debian) | `sudo apt install octave` |
| macOS (Homebrew) | `brew install octave` |
| Windows | Download from <https://octave.org/download> |

**Dynare** — install for Octave use:

| OS | Instructions |
|----|-------------|
| Linux (Ubuntu) | `sudo apt install dynare` (or <https://dynare.org/download>) |
| macOS | Download `.pkg` from <https://dynare.org/download>; choose "Octave" version |
| Windows | Download installer from <https://dynare.org/download>; choose "Octave" version |

After installing Dynare, note the path to its `matlab/` subfolder:
- Linux: `/usr/share/dynare/matlab`
- macOS: `/Applications/Dynare/6.x/matlab`
- Windows: `C:\Program Files\Dynare\6.x\matlab`

### Step 3 — Install Julia packages and run

```bash
# From this folder:
cd Model/julia_dynare

# Install all Julia dependencies (first time only, takes ~2 min)
julia --project=. -e "import Pkg; Pkg.instantiate()"

# Edit DYNARE_MATLAB_PATH in main_SOE_gap.jl to match your install,
# or set the environment variable:
export DYNARE_HOME="/usr/share/dynare/matlab"   # Linux
# export DYNARE_HOME="/Applications/Dynare/6.4/matlab"  # macOS
# set DYNARE_HOME=C:\Program Files\Dynare\6.4\matlab     # Windows (cmd)

# Run the model
julia --project=. main_SOE_gap.jl
```

---

## Data Files Required

The script searches for these files in common locations (sibling `modelo_chile/`
folder or via `NKIOSOE_DATA_DIR` environment variable):

| File | Notes |
|------|-------|
| `Stata_to_excel_few_industries_chile.xlsx` | **Save the original `.xls` as `.xlsx`** (Excel → Save As → `.xlsx`) |
| `IO_2021_chile.csv` | 12×12 input-output matrix |
| `fpa_vector_few_industries_chile.csv` | Price adjustment frequencies |

Set `NKIOSOE_DATA_DIR` if files are elsewhere:

```bash
export NKIOSOE_DATA_DIR="/path/to/your/data"
```

---

## Selecting an Exercise

Edit the `EXERCISE` variable near the top of `main_SOE_gap.jl`:

```julia
EXERCISE = 0   # 0 = Baseline (all shocks, matches SMM calibration)
               # 1 = Preference shock only
               # 2 = Manufacturing TFP shock only   ← default
               # 3 = Monetary policy shock only
```

---

## File Structure

```
julia_dynare/
├── Project.toml                  # Julia package dependencies
├── README.md                     # This file
├── main_SOE_gap.jl              # Main script — run this
├── steady_ntwsoe.jl             # Outer steady-state residuals (outer NLsolve)
├── steady_ntwsoe_system.jl      # Inner production-block system (inner NLsolve)
├── utils.jl                     # Helpers: Spearman, Lyapunov, run_dynare_octave
└── mod/                         # Dynare model files
    ├── NK_SOE_lev_gap2.mod      # Main model (levels + output-gap version)
    ├── definition_block_nsec.mod
    ├── definition_block_io.mod
    ├── definition_block_lab.mod
    └── solution_block.mod
```

---

## How It Works

| Step | MATLAB original | Julia equivalent |
|------|-----------------|-----------------|
| Data I/O | `readtable` (Excel) | `XLSX.readxlsx` |
| CSV I/O | `readmatrix` | `CSV.read` |
| `.mat` files | `load` / `save` | `MAT.matread` / `MAT.matwrite` |
| Nonlinear solve | `fsolve` | `NLsolve.nlsolve` (trust-region) |
| Dynare run | `dynare model.mod` | System call to Octave + Dynare |
| Spearman corr | `corr(..., 'Spearman')` | `StatsBase.corspearman` |
| Lyapunov eq | `dlyap` | `local_dlyap` (doubling algorithm) |
| Parameter save | `save params_val_ul.mat` | `MAT.matwrite` (v5 format) |

The Dynare `.mod` files are **unchanged** — Julia saves `params_val_ul.mat` that
Dynare's `load params_val_ul.mat` reads, exactly as in the MATLAB workflow.

---

## SMM Estimates

If `smm_estimates.mat` (or `smm_best_so_far.mat`) exists in the `modelo_chile/`
sibling folder, Julia loads it automatically and overrides the default parameter
values — identical behavior to the MATLAB script.

---

## Environment Variables Reference

| Variable | Default | Meaning |
|----------|---------|---------|
| `DYNARE_HOME` | `(edit in script)` | Path to Dynare's `matlab/` subfolder |
| `OCTAVE_EXE` | `octave` | Octave executable name/path |
| `NKIOSOE_DATA_DIR` | `(auto-search)` | Folder containing the three data files |
| `SMM_EXERCISE` | `(not set)` | Forces exercise index; used by SMM estimation |

---

## Troubleshooting

**"Required data file not found"**
- Convert `.xls` → `.xlsx` and place next to this script or set `NKIOSOE_DATA_DIR`.

**"Dynare results file not found"**
- Check that Octave and Dynare are on your PATH.
- Verify `DYNARE_HOME` points to the correct `matlab/` subfolder.
- Run `octave --eval "addpath('$DYNARE_HOME'); dynare('mod/NK_SOE_lev_gap2.mod')"` manually to see error output.

**"Outer steady-state solver did not converge"**
- The model is at the calibrated parameter values from `params_val.mat`.
  If parameters were changed, the initial guesses may need adjustment.

**Julia package errors**
- Run `julia --project=. -e "import Pkg; Pkg.status()"` to verify all packages are installed.
- Run `julia --project=. -e "import Pkg; Pkg.instantiate()"` to install missing ones.
