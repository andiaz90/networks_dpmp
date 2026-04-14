# Terms of Trade Shock Exercise - Complete Guide

## Overview

This exercise implements a **terms of trade shock** and compares impulse responses across three model configurations:

1. **Baseline**: Sticky price model with full Input-Output structure
2. **Flexible Price**: Identical to baseline but with no price rigidities
3. **No I-O**: Sticky price model with diagonal I-O matrix (no intermediate inputs)

**Key Feature**: All deviations are measured from **steady state**, not from flexible-price equilibrium. This uses main_SOE as the benchmark.

---

## Quick Start

### Option 1: Complete Automated Exercise (Recommended)

Run this single command in MATLAB:

```matlab
cd Model/modelo_chile
run_tot_shock_comparison
```

This will:
1. ✓ Automatically configure and run all 3 models
2. ✓ Extract IRF data from Dynare results
3. ✓ Generate comparison plots
4. ✓ Create summary statistics table

**Expected runtime**: 15-30 minutes (depends on Dynare compilation)

### Option 2: Manual Step-by-Step

If you prefer to run models individually:

```matlab
% Step 1: Configure and run baseline model
% Edit NK_IOSOE_lev.mod or use run_tot_shock_comparison

% Step 2: Run with configuration for flexible price
% (modify kappa parameter to ~0.001)

% Step 3: Run with diagonal I-O
% (modify beta matrix to eye(12))

% Step 4: Generate plots
plot_terms_of_trade_shock_results
```

---

## What Gets Compared

### Six Key Variables

| Variable | Description | Unit |
|----------|-------------|------|
| **Output (Y)** | Aggregate GDP | % deviation from SS |
| **Consumption (C)** | Total consumption | % deviation from SS |
| **CPI Inflation (PI_C)** | Consumer price inflation | Annualized % points |
| **Real Rate (RR)** | Real interest rate | Annualized % points |
| **Trade Balance / GDP (TB_Y)** | External balance ratio | % of GDP |
| **Real Exchange Rate (RX)** | Bilateral real rate | % deviation from SS |

### Three Model Configurations

| Model | Price Rigidity | I-O Structure | Use Case |
|-------|----------------|---------------|----------|
| **Baseline** | Sticky (Calvo) | Full (12×12 matrix) | Main model with frictions |
| **Flexible** | None (κ→0) | Full (12×12 matrix) | Frictionless benchmark |
| **No I-O** | Sticky (Calvo) | Diagonal (identity) | Isolate network effects |

---

## Understanding the Model Differences

### Baseline vs Flexible Price

**Identifies the role of sticky prices:**
- Flexible price shows **faster adjustment** to shock
- Baseline shows **persistence** due to price rigidities
- Larger divergence = stronger pricing frictions

**Example**: Terms of trade improvement
- Flexible: immediate price adjustment, quick output response
- Sticky: gradual price adjustment, slower output response

### Baseline vs No-I-O

**Identifies the role of supply chains:**
- Full I-O: sectoral shocks propagate through intermediate links
- Diagonal: shocks isolated to directly affected sectors
- Larger divergence = stronger network effects

**Example**: Manufacturing shock
- No I-O: only manufacturing sector responds
- Baseline: shocks ripple to suppliers and customers of manufacturing

### Flexible vs No-I-O

**Combined effects:**
- Price flexibility + supply chain isolation
- Shows interaction between frictions and network structure

---

## Key Metrics in Output Table

### Impact
The immediate (quarter 0) response to the shock. 

**Interpretation**:
- Large impact = fast transmission
- Small impact = delayed response

### Peak
The maximum absolute value during the 40-quarter horizon.

**Interpretation**:
- Shows maximum departure from steady state
- Later peaks suggest persistent effects

### Cumulative
Sum across all 40 quarters.

**Interpretation**:
- Positive cumsum = sustained above-SS response
- Negative cumsum = sustained below-SS response
- Magnitude shows total shock impact

### Final (Quarter 39)
Remaining deviation after 10+ years.

**Interpretation**:
- Near zero = complete mean reversion
- Significant value = persistent shock effect
- Sign shows direction of long-run shift

---

## Expected Results

### Terms of Trade Improvement Scenario

**Baseline Model Response:**

| Variable | Impact | Interpretation |
|----------|--------|-----------------|
| Output | ± (ambiguous) | Demand pull vs ToT effect |
| Consumption | ⬆️ | Wealth effect from ToT gain |
| Inflation | ⬆️ | Import prices rise |
| Real Rate | ⬆️ | Central bank tightens |
| Trade Balance | ⬇️ immediately | Demand surge then ⬆️ |
| Real Exchange Rate | ⬆️ | Currency appreciates |

**Flexible Price Model:**
- Larger immediate responses (no pricing friction)
- Faster mean reversion
- Cleaner peak responses

**No-I-O Model:**
- Smaller sectoral responses (no propagation)
- Direct sector shock is isolated
- Less aggregate volatility

---

## Customization

### Change the Shock

To use a different exogenous variable:

In `plot_terms_of_trade_shock_results.m`, change:
```matlab
shock_var = 'tot_shock';  % Change to any exogenous variable name
```

Available shocks in your model (modify NK_IOSOE_lev.mod):
- Preference shocks (`eps_pref`)
- Productivity shocks (`eps_a_i`)
- Monetary policy shocks (`eps_r`)
- Financial shocks, etc.

### Change Horizon

To extend the comparison beyond 40 quarters:

In `plot_terms_of_trade_shock_results.m`:
```matlab
nPeriods = 60;  % Show 60 quarters instead of 40
```

In Dynare model (NK_IOSOE_lev.mod):
```dynare
stoch_simul(order=1, irf=60, ...)  % Increase irf horizon
```

### Adjust Model Parameters

Before calling `run_tot_shock_comparison`:

```matlab
% Make more flexible prices
kappa_flex = ones(nsec,1) * 0.001;  % Smaller = more flexible

% Or less flexible
kappa_sticky = ones(nsec,1) * 100;  % Larger = stickier

% Run with custom configuration
assignin('base', 'modkappa', kappa_sticky);
dynare NK_IOSOE_lev.mod noclearall;
```

---

## Output Files

### Generated During Run

```
dynare_tot_shock_Baseline.txt          ← Dynare output log (baseline)
dynare_tot_shock_Flexible Price.txt    ← Dynare output log (flexible)
dynare_tot_shock_No I-O.txt            ← Dynare output log (no I-O)

model_output_tot_shock_Baseline.mat         ← Dynare results (baseline)
model_output_tot_shock_Flexible Price.mat   ← Dynare results (flexible)
model_output_tot_shock_No I-O.mat           ← Dynare results (no I-O)

tot_shock_comparison.png               ← Main comparison figure
```

### Reusing Results

To regenerate plots without re-running Dynare:

```matlab
clear all;
% Load saved results
load model_output_tot_shock_Baseline.mat;
results(1).oo_ = oo_;
results(1).M_ = M_;
results(1).success = true;

% Repeat for other models...

% Generate plots
plot_terms_of_trade_shock_results;
```

---

## Troubleshooting

### Error: "Variable not found"

**Problem**: "Y variable not found in model"

**Solution**: 
1. Verify variable is defined in NK_IOSOE_lev.mod
2. Check it's included in `stoch_simul()` line:
   ```dynare
   stoch_simul(...) Y C PI_C RR TB_Y RX ...
   ```

### Error: "Shock not found"

**Problem**: "tot_shock variable not found"

**Solution**:
1. Add to varexo section in NK_IOSOE_lev.mod:
   ```dynare
   varexo tot_shock;
   ```
2. Add shock process to model section
3. Include in stoch_simul output list

### IRF Data Shape Mismatch

**Problem**: Array dimension errors during comparison

**Solution**:
1. Ensure all Dynare runs use same `order=1` approximation
2. Use same `irf=40` horizon for all
3. Check variable names are consistent across models

### Path Errors

**Problem**: "Cannot find NK_IOSOE_lev.mod"

**Solution**:
1. Ensure you're in `Model/modelo_chile` directory
2. Or add to MATLAB path:
   ```matlab
   addpath('C:\path\to\Model\modelo_chile')
   ```

---

## Interpreting the Figure

### Layout: 2×3 Grid

```
[Output]        [Consumption]    [CPI Inflation]
[Real Rate]     [Trade Balance]  [Real Exchange Rate]
```

### Reading Each Subplot

- **X-axis**: Quarters (0 to 39)
- **Y-axis**: Deviation from steady state (or %)
- **Blue solid (—)**: Baseline model
- **Orange dashed (- -)**: Flexible Price model
- **Green dash-dot (-.-.)**: No-I-O model
- **Gray dashed line at 0**: Steady state reference

### What Diverging Lines Mean

- **Baseline ≫ Flexible**: Strong price stickiness effect
- **Baseline ≫ No-I-O**: Strong I-O propagation effect
- **All converge to zero**: Model shows good mean reversion
- **Persistent divergence**: Shock has long-run effects

---

## Advanced Usage

### Comparing Multiple Shocks

Create a loop version:

```matlab
shocks_list = {'tot_shock', 'eps_a_2', 'eps_pref'};

for s = 1:length(shocks_list)
    shock_var = shocks_list{s};
    run_tot_shock_comparison;  % Will use this shock_var
end
```

### Sensitivity Analysis

Test different parameter values:

```matlab
% Test different price rigidity levels
kappa_values = [0.001, 0.01, 0.1, 1, 10];

for kappa_val = kappa_values
    modkappa = ones(nsec,1) * kappa_val;
    assignin('base', 'modkappa', modkappa);
    % Run model...
end
```

### Network Analysis

Extract sectoral responses:

```matlab
% After loading results
for i = 1:12
    Y_i_field = sprintf('Y_%d_%s', i, shock_var);
    if isfield(oo_.irfs, Y_i_field)
        sector_irf(:,i) = oo_.irfs.(Y_i_field);
    end
end
% Visualize sectoral dispersion
heatmap(sector_irf, 'YLabel', 'Quarters', 'XLabel', 'Sector')
```

---

## References

- **NK-IOSOE Model**: Ferrante, Graves & Iacoviello (2023), *Journal of Monetary Economics* 140, S64-S81
- **Terms of Trade Analysis**: Mendoza (1995), Kehoe & Midrigan (2007)
- **I-O Network Effects**: di Giovanni & Kalemli-Özcan (2008), Carvalho & Gabaix (2013)
- **Price Rigidities**: Calvo (1983), Rotemberg (1982)

---

## Files in This Exercise

```
Model/modelo_chile/
├── run_tot_shock_comparison.m              ← Main script (RUN THIS)
├── plot_terms_of_trade_shock_results.m     ← Plotting script (auto-called)
├── ToT_SHOCK_IMPLEMENTATION_GUIDE.md       ← Detailed implementation guide
└── ToT_SHOCK_README.md                     ← This file

Supporting model files (pre-existing):
├── NK_IOSOE_lev.mod                        ← Dynare model
├── main_SOE.m                              ← Baseline driver
├── main_SOE_ul.m                           ← Flexible price driver
└── figs_SOE_ss.m                           ← Steady state plots
```

---

## Contact & Support

For issues or questions:
1. Check "Troubleshooting" section above
2. Review Dynare error messages in output logs
3. Verify NK_IOSOE_lev.mod includes all shock and variable declarations
4. Ensure data files are accessible and correctly formatted

---

## Version History

**v1.0** (March 2026)
- Initial release
- Supports 3 model configurations
- 6-variable comparison framework
- Flexible shock selection
