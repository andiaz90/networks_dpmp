# Terms of Trade Shock Exercise - Implementation Guide

## Overview

This guide explains how to add a **terms of trade shock** to your existing NK_IOSOE_lev.mod model and run comparative analysis across 3 models:

1. **Baseline**: Sticky price with Input-Output structure
2. **Flexible Price**: Same I-O structure but with no price rigidities  
3. **No I-O**: Sticky price with diagonal I-O matrix (no intermediate inputs)

---

## What is a Terms of Trade Shock?

The terms of trade (ToT) is the ratio of export prices to import prices:
$$\text{ToT} = \frac{P_X}{P_M}$$

A **terms of trade shock** is an exogenous change in this ratio, often driven by:
- Global commodity price changes
- External demand shocks
- Currency movements

---

## Implementation Steps

### Step 1: Add Terms of Trade Variables to Dynare Model

Add these lines to **NK_IOSOE_lev.mod** in the variable declarations section:

```dynare
// Terms of Trade variables
var Tot_shock;      % Exogenous terms of trade index
var Tot;            % Log real terms of trade
var Pf_ss;          % Foreign prices (steady state reference)
```

### Step 2: Add to Exogenous Variables

In the `varexo` section of your Dynare model, add:

```dynare
varexo 
    tot_shock;      % Terms of trade shock (exogenous)
    % ... other shocks ...
end;
```

### Step 3: Add Shock Process Definition

In the `model()` section, define the ToT shock process:

```dynare
% Terms of Trade shock process (AR(1))
Tot_shock = rho_tot * Tot_shock(-1) + e_tot;

% Log terms of trade (with shock influence)
Tot = Tot_ss + Tot_shock;
```

where `rho_tot` is the persistence parameter (0.95 is typical).

### Step 4: Add Impulse Response Configuration

In the `shocks;` section, specify the shock size:

```dynare
shocks;
    var tot_shock = 0.01;   % 1% trade shock
end;
```

### Step 5: Run Dynare with Estimation Settings

Add stoch_simul for IRF generation:

```dynare
stoch_simul(order=1, irf=40, periods=0) Y C PI_C RR TB_Y RX Tot_shock;
```

---

## Variable Mapping

The comparison plots expect these variables in your Dynare output:

| Variable | Description | Expected Unit |
|----------|-------------|----------------|
| `Y` | Aggregate output | Log deviation from SS |
| `C` | Consumption | Log deviation from SS |
| `PI_C` | CPI inflation | Annualized percentage points |
| `RR` | Real interest rate | Annualized percentage points |
| `TB_Y` | Trade balance / GDP | Percentage of GDP |
| `RX` | Real exchange rate | Log deviation from SS |

**Important**: All variables should be **log deviations from steady state**, NOT deviations from flexible price equilibrium.

---

## Running the Exercise

### Option A: Use the Provided Script

```matlab
% In MATLAB, from the modelo_chile folder:
terms_of_trade_shock_comparison;
```

This script will:
1. Run the baseline model (main_SOE parameters)
2. Run the flexible price version
3. Run the diagonal I-O version
4. Generate comparison plots automatically

### Option B: Manual Multi-Model Run

If you prefer to run each model separately:

```matlab
% 1. Run baseline
main_SOE;  % Modify EXERCISE = 0 (baseline)

% 2. Run flexible price version
main_SOE_ul;  % Already has flexible price setup

% 3. Run diagonal I-O version
% Modify NK_IOSOE_lev.mod: set beta = diag(ones(12,1))
% Then run main_SOE again

% 4. Generate comparison
plot_tot_shock_comparison;
```

---

## Understanding the Results

### Key Comparisons

**Baseline vs Flexible Price:**
- Shows effect of sticky prices on adjustment
- Flexible prices should show faster return to SS

**Baseline vs Diagonal I-O:**
- Shows importance of supply chain linkages
- Diagonal means no intermediate input sharing
- Should show smaller sectoral reallocation

### Expected Response Patterns

**Terms of Trade Improvement** (ToT increases):

1. **Output**: 
   - May increase (export sector boom) or decrease (relative price effect)
   - Flexible price responds faster

2. **Consumption**:
   - Usually increases (terms of trade wealth effect)
   - Sticky price: delayed response

3. **Inflation**:
   - Import-weighted CPI rises in ToT improvement
   - Flexible price: larger immediate jump

4. **Real Interest Rate**:
   - Central bank may raise rates to moderate demand
   - Flexible price: no nominal rigidity buffer

5. **Trade Balance**:
   - Initially worsens (demand for imports rises)
   - Then improves as export boom materializes

6. **Real Exchange Rate**:
   - Currency appreciates with ToT improvement
   - Response slower under sticky prices

---

## Customization Options

### Adjust Shock Persistence

In Dynare, modify `rho_tot`:
```dynare
rho_tot = 0.80;  % Less persistent shock (0-1)
```

Lower values = faster mean reversion  
Higher values = prolonged shock effect

### Change Shock Size

In Dynare `shocks;` block:
```dynare
var tot_shock = 0.05;  % 5% shock (larger impact)
```

### Change IRF Horizon

In MATLAB script or Dynare:
```matlab
nPeriods = 60;  % Show 60 quarters instead of 40
```

---

## Model-Specific Notes

### Baseline Model Setup
- Uses main_SOE.m parameters
- Sticky prices active (Calvo/Rotemberg)
- Full I-O matrix structure
- Shows deviations from steady state

### Flexible Price Setup
- Zero or near-zero price adjustment costs (κ → 0)
- Otherwise identical to baseline
- Represents frictionless reference equilibrium
- **Deviations still measured from steady state**

### Diagonal I-O Setup
- Replaces β matrix with identity matrix
- Each sector only produces for final demand
- No intermediate input relationships
- Isolates network effects of I-O structure

---

## Troubleshooting

### "Variable not found" Error

**Problem**: "Y variable not found in model"

**Solution**: Ensure the variable is defined in NK_IOSOE_lev.mod and included in stoch_simul list

### IRF Data Shape Mismatch

**Problem**: Dimension errors when comparing IRFs

**Solution**: Ensure all Dynare runs use same:
- Order of approximation (order=1)
- IRF horizon (irf=40)
- Same variable list in stoch_simul

### Path Errors

**Problem**: "model_output_tot_shock_baseline.mat not found"

**Solution**: 
1. Ensure script is run from modelo_chile folder
2. Check that Dynare runs completed successfully
3. Verify save paths match between scripts

---

## References

- **NK_IOSOE Model**: Ferrante, Graves & Iacoviello (2023)
- **Terms of Trade Shocks**: Mendoza (1995), Kehoe & Midrigan (2007)
- **I-O Network Effects**: di Giovanni & Kalemli-Özcan (2008)

---

## Files Generated

```
Output Files:
├── terms_of_trade_shock_comparison.m      ← Main script (run this)
├── plot_tot_shock_comparison.m            ← Plotting script (called automatically)
├── model_output_tot_shock_baseline.mat    ← Baseline model results
├── model_output_tot_shock_flexible.mat    ← Flexible price results
├── model_output_tot_shock_no_io.mat       ← Diagonal I-O results
├── terms_of_trade_shock_comparison.png    ← Comparison figure
└── dynare_tot_shock_*.txt                 ← Dynare output logs
```

---

## Next Steps

1. **Add ToT variables** to your NK_IOSOE_lev.mod
2. **Test single model run** with terms of trade shock
3. **Verify IRF output** matches expected variable list
4. **Run full comparison** with terms_of_trade_shock_comparison.m
5. **Interpret results** using the "Understanding the Results" section above
