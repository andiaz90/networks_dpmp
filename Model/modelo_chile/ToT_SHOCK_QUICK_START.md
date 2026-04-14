# Terms of Trade Shock Exercise - QUICK REFERENCE

## 🚀 Get Started in 2 Minutes

### Step 1: Navigate to correct folder
```matlab
cd Model/modelo_chile
```

### Step 2: Validate setup (optional but recommended)
```matlab
validate_tot_shock_setup
```

### Step 3: Run the exercise
```matlab
run_tot_shock_comparison
```

**That's it!** The script will automatically:
- ✓ Configure 3 model versions
- ✓ Run Dynare for each configuration
- ✓ Extract and compare IRFs
- ✓ Generate 6-variable comparison figure
- ✓ Produce summary statistics table

**⏱️ Expected time: 15-30 minutes**

---

## 📁 Files Created by This Exercise

After running, you'll have:

```
Dynare Outputs:
├── dynare_tot_shock_Baseline.txt
├── dynare_tot_shock_Flexible Price.txt
├── dynare_tot_shock_No I-O.txt

MATLAB Outputs:
├── model_output_tot_shock_Baseline.mat
├── model_output_tot_shock_Flexible Price.mat
├── model_output_tot_shock_No I-O.mat

Figures:
└── tot_shock_comparison.png               ← THE MAIN RESULT
```

---

## 📊 What You'll Get

### Main Output: tot_shock_comparison.png

A 2×3 figure with IRF comparisons for:

```
[Output]              [Consumption]         [Inflation]
[Real Interest Rate]  [Trade Balance/GDP]   [Real Exchange Rate]
```

Each subplot shows 3 lines:
- **Blue solid** = Baseline (sticky price + I-O)
- **Orange dashed** = Flexible Price (+ I-O)
- **Green dash-dot** = No I-O (sticky price only)

### Console Output: Summary Table

```
Output
===================================================================
Model                               Impact    Peak   Peak Qtr
-------------------------------------------------------------------
Baseline                           -0.5234  -1.2341     12
Flexible Price                     -0.8234  -1.0234      8
No I-O                             -0.3421  -0.7234     15
```

---

## 🔎 Key Interpretations

### When Baseline > Flexible Price
👉 **Price stickiness matters** - sticky prices slow adjustment

### When Baseline > No-I-O  
👉 **Supply chains matter** - I-O linkages amplify shock

### When Flexible ≈ No-I-O
👉 **Frictions dominate** - pricing matters more than networks

### When All lines converge to zero
👉 **Good mean reversion** - shock is temporary

### When lines stay diverged
👉 **Long-run effects** - persistent shock impact

---

## ⚙️ Customization

### Use a Different Shock

In `plot_terms_of_trade_shock_results.m`, line ~20:
```matlab
shock_var = 'tot_shock';    % Change this
```

Available shocks (from NK_IOSOE_lev.mod):
- `eps_a_i` - Productivity shocks
- `eps_pref` - Preference shocks  
- `eps_r` - Monetary policy shocks
- Any varexo in your model

### Extend the Horizon

Change in `plot_terms_of_trade_shock_results.m`, line ~18:
```matlab
nPeriods = 60;   % Default is 40, increase for longer view
```

And in NK_IOSOE_lev.mod:
```dynare
stoch_simul(order=1, irf=60, ...)   % Match to nPeriods
```

### Compare Only 2 Models

Comment out unwanted model in `run_tot_shock_comparison.m`:
```matlab
% Skip flexible price run:
% run_terms_of_trade_shock_model('flexible', ...);
```

---

## 🐛 Troubleshooting

| Problem | Solution |
|---------|----------|
| "Y variable not found" | Add to NK_IOSOE_lev.mod stoch_simul list |
| "tot_shock not found" | Add `varexo tot_shock;` to .mod file |
| "File not found" | Check data_path in run_tot_shock_comparison |
| Dynare error | Review dynare_tot_shock_*.txt files for details |
| MATLAB crashes | Reduce nPeriods or check available RAM |

**First time troubleshooting?**
→ See `validate_tot_shock_setup` for automated diagnosis

---

## 📝 Understanding the Variables

| Variable | What it measures | Interpretation |
|----------|------------------|-----------------|
| **Y** | Total output/GDP | Economy size relative to SS |
| **C** | Total consumption | Living standard change |
| **PI_C** | Consumer price inflation | Cost of living increase/decrease |
| **RR** | Real interest rate | Return on savings (inflation-adjusted) |
| **TB_Y** | Trade balance as % GDP | External solvency indicator |
| **RX** | Real exchange rate | International competitiveness |

---

## 🧮 The Three Models Compared

### Model 1: BASELINE
- 🔒 **Sticky prices** (Calvo parameter θ)
- 🕸️ **Full I-O** (12×12 matrix from actual Chilean data)
- 📈 **Most realistic** for current economy
- Use for: Actual policy effects

### Model 2: FLEXIBLE PRICE
- 🔓 **No price frictions** (κ → 0, θ → 0)
- 🕸️ **Same I-O structure** (unchanged)
- 📉 **Perfectly competitive** benchmark
- Use for: Show importance of rigidities

### Model 3: NO I-O
- 🔒 **Sticky prices** (same as baseline)
- ◼️ **Diagonal I-O** (identity matrix)
- 🏭 **Self-sufficient sectors** (no intermediate trade)
- Use for: Show importance of supply chains

---

## 📚 Documentation

| File | Purpose |
|------|---------|
| `ToT_SHOCK_README.md` | Complete user guide (15 min read) |
| `ToT_SHOCK_IMPLEMENTATION_GUIDE.md` | How to add ToT shock to model |
| `validate_tot_shock_setup.m` | Auto-check system before running |
| This page | Quick reference for experienced users |

---

## 🎯 Advanced Usage

### Re-run plots without re-running Dynare
```matlab
% Load previous results
load model_output_tot_shock_Baseline.mat
results(1).oo_ = oo_;
results(1).M_ = M_;
results(1).success = true;

% Repeat for other models...

% Just plot
plot_terms_of_trade_shock_results
```

### Extract specific sector response
```matlab
% Inside plot_terms_of_trade_shock_results.m, add:
Y_2_field = 'Y_2_tot_shock';  % Sector 2 (Mining)
if isfield(oo_.irfs, Y_2_field)
    sector_2_irf = oo_.irfs.(Y_2_field);
    plot(sector_2_irf, 'k-', 'LineWidth', 2);
end
```

### Batch comparison of multiple shocks
```matlab
shocks = {'tot_shock', 'eps_a_2', 'eps_pref'};
for s = 1:length(shocks)
    shock_var = shocks{s};
    fprintf('Analyzing %s...\n', shock_var);
    plot_terms_of_trade_shock_results;
    saveas(gcf, sprintf('comparison_%s.png', shock_var));
end
```

---

## 💡 Common Questions

**Q: Why compare to steady state, not flexible price?**  
A: Your model (main_SOE) is sticky-price without flexible-price counterpart. Using SS as benchmark allows consistent measurement across all variants.

**Q: Which model should I trust most?**  
A: Baseline - it includes your actual Chilean I-O data and price rigidities. Flexible/No-I-O are useful **comparisons**.

**Q: Can I add a new shock type?**  
A: Yes! Add to NK_IOSOE_lev.mod `varexo` block, define process in equation block, and list in `stoch_simul()`.

**Q: How do I interpret trade balance response?**  
A: ToT improvement usually causes initial trade deficit (demand for imports rises), then improves (exports boom). Sticky prices delay this.

---

## 📈 Expected Output Patterns

### Terms of Trade Improvement (1% shock)

| Model | Y effect | C effect | Effect size |
|-------|----------|----------|-------------|
| Baseline | ±0.2-0.5% | +0.5-1.0% | Medium |
| Flexible | +0.3-0.7% | +0.8-1.5% | Larger |
| No-I-O | ±0.1-0.3% | +0.4-0.8% | Smaller |

**Why different sizes?**
- Flexible price: faster, stronger response
- No-I-O: isolated to directly affected sectors
- Baseline: intermediate due to frictions + networks

---

## ✅ Verification Checklist

Before asking for help, verify:

- [ ] Running from `Model/modelo_chile` folder
- [ ] All data files accessible (validate_tot_shock_setup confirms)
- [ ] Dynare installed and on path (validate script checks)
- [ ] NK_IOSOE_lev.mod includes Y, C, PI_C, RR, TB_Y, RX variables
- [ ] varexo section exists in .mod file
- [ ] Recent MATLAB version (R2019b or later recommended)

---

## 🚨 When Something Goes Wrong

### Step 1: Run validation
```matlab
validate_tot_shock_setup
```

### Step 2: Check Dynare output
```matlab
% Read the error log
type dynare_tot_shock_Baseline.txt
```

### Step 3: Review requirements in .mod file
```dynare
% Should contain:
var Y C PI_C RR TB_Y RX ...;
varexo tot_shock ...;
% ... model equations ...
stoch_simul(order=1, irf=40, periods=0) Y C PI_C RR TB_Y RX ...;
```

### Step 4: Check data files
```matlab
% Verify files exist:
ls C:\Users\adiaz\OneDrive*\data_github\*.xls
ls C:\Users\adiaz\OneDrive*\data_github\*.csv
```

---

## 📞 Support Resources

- **Complete guide**: `ToT_SHOCK_README.md`
- **Implementation details**: `ToT_SHOCK_IMPLEMENTATION_GUIDE.md`  
- **System check**: Run `validate_tot_shock_setup`
- **Dynare help**: See dynare_tot_shock_*.txt output files
- **MATLAB help**: `help dynare` or check Dynare documentation

---

## 🎓 Learning Path

1. **Start here** (this page) - 5 min
2. Run `validate_tot_shock_setup` - 2 min
3. Run `run_tot_shock_comparison` - 20 min
4. Inspect `tot_shock_comparison.png` - 10 min
5. Read console table output - 5 min
6. Optional: Read `ToT_SHOCK_README.md` for theory - 15 min

**Total: ~1 hour to complete exercise**

---

Last updated: March 2026  
Compatible with: MATLAB R2019b+, Dynare 6.2+, Chile 12-sector model

