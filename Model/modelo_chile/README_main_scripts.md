# Comparison: main_SOE.m vs main_SOE_ul.m

## Overview

Both `main_SOE.m` and `main_SOE_ul.m` are MATLAB driver scripts for the 12-sector New Keynesian Small Open Economy (SOE) DSGE model for Chile based on Ferrante, Graves & Iacoviello (2023). The key difference lies in **model specification and flexibility**:

| Aspect | main_SOE.m | main_SOE_ul.m |
|--------|-----------|-------------|
| **Model Type** | Sticky price (baseline) | Flexible price with flexible labor |
| **Primary Use** | Sticky price equilibrium analysis | Output gap analysis (flexible vs sticky) |
| **Key Plotting Script** | `figs_SOE_ss.m` | `figs_SOE_gap.m` |
| **Output Metrics** | Deviations from steady state | Output gaps (Y vs Y_f) |
| **File Size** | 752 lines | 872 lines |

---

## Detailed Differences

### 1. **Model Specification**

**main_SOE.m (Sticky Price Model)**
- Implements the standard NK-IOSOE model with sticky prices
- Firms cannot freely adjust prices (governed by Calvo or Rotemberg pricing)
- Labor market has frictions and adjustment costs
- Represents the actual economy with nominal rigidities
- Comment: `% NK_IOSOE_lev: 12-sector SOE DSGE model for Chile`

**main_SOE_ul.m (Flexible Price / "Unconstrained Labor" Model)**
- `ul` suffix likely indicates "unconstrained labor" or flexible labor reference equilibrium
- Includes flexible-price allocations (Y_f) for comparison
- Allows firms to adjust prices freely and immediately
- Serves as the counterfactual benchmark equilibrium
- Used to calculate output gaps as deviations from this flexible-price equilibrium

### 2. **Plotting and Output Analysis**

**main_SOE.m: Steady State Deviations**
```matlab
% Line 716 - Common plots for all exercises (STEADY STATE DEVIATIONS)
figs_SOE_ss  % Plots deviations from steady state
```
- Plots show how variables deviate from **steady state** (Y_ss)
- Focuses on the magnitude of business cycle fluctuations
- No reference to flexible-price allocation
- Useful for: Understanding cyclical behavior, impulse responses relative to trend

**main_SOE_ul.m: Output Gaps**
```matlab
% Line 832 - Common plots for all exercises (OUTPUT GAPS)
figs_SOE_gap  % Plots output gaps: deviation from flexible price equilibrium
```
- Plots show **output gaps** = deviation from flexible-price equilibrium (Y - Y_f)
- Measures the distortion caused by nominal rigidities and labor frictions
- Separates trend changes from cyclical deviations with pricing/labor constraints
- Useful for: Analyzing monetary policy effectiveness, welfare implications of sticky prices

### 3. **Data Path Handling**

**main_SOE.m (Hardcoded Paths)**
```matlab
% Lines 72-74 - Fixed paths
data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';
filename_industries = fullfile(data_path, 'Stata_to_excel_few_industries_chile.xls');
```
- Requires exact path hardcoding
- May fail if user has different OneDrive structure or username
- Less portable across different machines

**main_SOE_ul.m (Adaptive Path Search)**
```matlab
% Lines 71-122 - Multi-level fallback search
candidates = {
    pwd,                                % current working directory
    fullfile(repo_root, 'Process Data Codes'),
    fullfile(repo_root, 'Model', 'modelo_chile'),
    'C:\Users\mgiarda\OneDrive - Banco Central de Chile\...',
    'C:\Users\adiaz\OneDrive - Banco Central de Chile\...'
};
```
- Dynamically searches multiple folders:
  1. Current working directory (pwd)
  2. Repository root relative paths
  3. Multiple user-specific paths
  4. Script directory
  5. Plain filenames in current folder
- **Much more robust** for different user setups
- Provides informative error messages listing all attempted paths
- **Recommended for collaboration and shared projects**

### 4. **Dynare Integration**

**main_SOE.m (Basic)**
```matlab
% Lines 70-71 - Fixed path
dynare_path = 'C:\Program Files\dynare\6.2\matlab';
addpath(dynare_path);
```

**main_SOE_ul.m (Robust with Fallbacks)**
```matlab
% Lines 582-622 - Comprehensive Dynare detection
if isempty(which('dynare'))
    fprintf('Dynare not found on MATLAB path. Searching common install locations...\n');
    % Check environment variable
    dyn_home = getenv('DYNARE_HOME');
    % Search common installation folders
    roots = {'C:\Program Files\Dynare', 'C:\Program Files (x86)\Dynare', 'C:\Dynare'};
    % Multiple version subdirectory search...
```
- Checks if Dynare is already on path
- Reads `DYNARE_HOME` environment variable
- Searches common installation directories
- Tests multiple Dynare version subdirectories
- More compatible with different Dynare installations

---

## When to Use Which Script

### Use **main_SOE.m** when:
- ✓ Analyzing sticky-price model responses only
- ✓ You want steady-state deviation plots (deviations from trend)
- ✓ You're not interested in flexible-price counterfactuals
- ✓ You have a standard setup with fixed OneDrive paths
- ✓ Simpler analysis focused on nominal rigidity effects

### Use **main_SOE_ul.m** when:
- ✓ You need output gaps (Y vs Y_f comparison)
- ✓ Analyzing welfare costs of prices/labor frictions
- ✓ Monetary policy effectiveness analysis
- ✓ Contributing to a collaborative project
- ✓ You have a non-standard folder setup or different user
- ✓ You want more robust Dynare detection
- ✓ Your Dynare installation is in a non-standard location

---

## Exercise Structure

Both scripts support the same 4 exercise modes (controlled by `EXERCISE` variable):

```matlab
EXERCISE = 2;  % <<<< CHANGE THIS TO SELECT EXERCISE
% 0 = Baseline (all shocks)
% 1 = Preference shock only
% 2 = Manufacturing TFP shock only
% 3 = Monetary policy shock only
```

The difference is only in **how results are visualized** (steady-state deviations vs. output gaps).

---

## File Organization

```
Model/modelo_chile/
├── main_SOE.m                          ← Sticky price model (basic)
├── main_SOE_ul.m                       ← Flexible price comparison (robust)
│
├── figs_SOE_ss.m                       ← Used by main_SOE
│                                          (steady state deviations)
├── figs_SOE_gap.m                      ← Used by main_SOE_ul
│                                          (output gaps)
│
├── NK_IOSOE_lev.mod                    ← Dynare model file (both use this)
│
├── main_SOE_ul/                        ← Exercise-specific plots
│   ├── plot_shock_effects.m
│   ├── plot_manufacturing_shock.m
│   └── plot_figure7_manufacturing_shock.m
│
└── utils/                              ← Shared utilities
```

---

## Recent Updates

As of March 2026:
- **figs_SOE_gap.m** has been updated to correctly handle models **without** flexible-price allocations
- Figures 4-7 now display **deviations from steady state (Y_ss)** rather than output gaps (Y_f)
- This makes the plotting consistent with actual model capabilities

---

## Troubleshooting

### Path Error: "Cannot find data files"
**Solution:** Use `main_SOE_ul.m` instead - it has better path detection

### Dynare Not Found Error  
**Solution:** Use `main_SOE_ul.m` - it includes Dynare path searching

### Need to modify paths
**In main_SOE_ul.m**, edit the `candidates` cell array (lines 75-80) to add your custom paths

---

## References

- **Base Model:** Ferrante, Graves & Iacoviello (2023), *Journal of Monetary Economics*, 140, pp. S64-S81
- **Chile Model:** 12-sector NK-IOSOE model with input-output structure
- **Extensions:** Open economy features, sectoral heterogeneity, network effects
