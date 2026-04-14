# CES Calibration Guide: Mapping Parameters to Expenditure Shares

## Overview

This guide explains how CES (Constant Elasticity of Substitution) function parameters in the model map to expenditure shares observed in the data. Understanding this mapping is essential for calibrating the model to match input-output tables and consumption data.

---

## Key CES Functions in the Model

The model uses several nested CES aggregators:

1. **Production Function** (by sector)
2. **Intermediate Input Aggregator** (by sector)  
3. **Consumption: Goods vs. Services** (at household level)
4. **Within-category Consumption** (goods basket, services basket)
5. **Home vs. Foreign Varieties** (by good/sector)

---

## 1. Production Function (CES over Labor, Intermediates, Imports)

### Functional Form

For sector $i$, production is:

$$Y_i = A_i \left[ \alpha_i^{1/\varepsilon_i^Y} M_i^{(\varepsilon_i^Y-1)/\varepsilon_i^Y} + \alpha_i^{V,1/\varepsilon_i^Y} V_i^{(\varepsilon_i^Y-1)/\varepsilon_i^Y} + (1-\alpha_i-\alpha_i^V)^{1/\varepsilon_i^Y} L_i^{(\varepsilon_i^Y-1)/\varepsilon_i^Y} \right]^{\varepsilon_i^Y/(\varepsilon_i^Y-1)}$$

Where:
- $M_i$ = composite intermediate inputs  
- $L_i$ = labor
- $V_i$ = imported intermediates
- $\alpha_i$ = intermediate input share parameter
- $\alpha_i^V$ = import share parameter
- $\varepsilon_i^Y$ = elasticity of substitution in production

### First-Order Conditions (Cost Minimization)

The FOCs give us input demand functions:

$$M_i = \left(\frac{MC_i}{PM_i}\right)^{\varepsilon_i^Y} \alpha_i \cdot Y_i$$

$$L_i = \left(\frac{MC_i}{w}\right)^{\varepsilon_i^Y} (1-\alpha_i-\alpha_i^V) \cdot Y_i$$

$$V_i = \left(\frac{MC_i}{P^V}\right)^{\varepsilon_i^Y} \alpha_i^V \cdot Y_i$$

### Mapping to Expenditure Shares

**In steady state with all prices equal to 1**, the expenditure shares are:

$$\text{Share of intermediates in total costs} = \frac{PM_i \cdot M_i}{MC_i \cdot Y_i} = \alpha_i$$

$$\text{Share of labor in total costs} = \frac{w \cdot L_i}{MC_i \cdot Y_i} = 1 - \alpha_i - \alpha_i^V$$

$$\text{Share of imports in total costs} = \frac{P^V \cdot V_i}{MC_i \cdot Y_i} = \alpha_i^V$$

**Key Insight:** When $\varepsilon_i^Y = 1$ (Cobb-Douglas case), the CES share parameters **equal** the expenditure shares **exactly**, regardless of prices.

**For general CES ($\varepsilon_i^Y \neq 1$):** The relationship is:

$$\text{Expenditure Share} = \frac{\alpha_i \cdot (PX_i/PM_i)^{\varepsilon_i^Y-1}}{\alpha_i \cdot (PX_i/PM_i)^{\varepsilon_i^Y-1} + \alpha_i^V \cdot (PX_i/P^V)^{\varepsilon_i^Y-1} + (1-\alpha_i-\alpha_i^V) \cdot (PX_i/w)^{\varepsilon_i^Y-1}}$$

Where $PX_i$ is the common normalization (often the output price).

### Calibration Strategy

Since $\varepsilon_i^Y$ is typically small (e.g., 0.8 in your model), the model is "close" to Cobb-Douglas, so:

```matlab
% From data:
share_intermediate = total_intermediate_costs / total_costs
share_labor = total_labor_costs / total_costs  
share_imports = total_import_costs / total_costs

% Set parameters:
alpha(i) = share_intermediate(i)
alphaV(i) = share_imports(i)
% Automatically: (1-alpha(i)-alphaV(i)) = share_labor(i)
```

This is what happens in [main_SOE.m](main_SOE.m#L87):
```matlab
alpha      = table2array(alldata(:,12));  % Intermediate share
alpha_V    = table2array(alldata(:,13));  % Import share
```

---

## 2. Intermediate Input Aggregator (CES over input varieties)

### Functional Form

Sector $i$ combines inputs from all sectors into a composite intermediate:

$$M_i = \left[ \sum_{j=1}^{n} \beta_{i,j} \cdot (M_{ji})^{(\varepsilon_i^M-1)/\varepsilon_i^M} \right]^{\varepsilon_i^M/(\varepsilon_i^M-1)}$$

Where:
- $M_{ji}$ = amount of good $j$ used as intermediate by sector $i$
- $\beta_{i,j}$ = share parameter for input $j$ in sector $i$'s intermediate bundle
- $\varepsilon_i^M$ = elasticity of substitution between intermediate inputs

### Price Index

The corresponding price index is:

$$PM_i = \left[ \sum_{j=1}^{n} \beta_{i,j} \cdot (PH_j)^{1-\varepsilon_i^M} \right]^{1/(1-\varepsilon_i^M)}$$

### Demand for Each Input

Cost minimization gives:

$$M_{ji} = \beta_{i,j} \left(\frac{PM_i}{PH_j}\right)^{\varepsilon_i^M} M_i$$

### Mapping to Expenditure Shares

The expenditure on input $j$ by sector $i$ is:

$$\text{Expenditure}_{ji} = PH_j \cdot M_{ji} = \beta_{i,j} \left(\frac{PH_j}{PM_i}\right)^{1-\varepsilon_i^M} \cdot PM_i \cdot M_i$$

The **expenditure share** of input $j$ in sector $i$'s total intermediate spending is:

$$s_{ji} = \frac{PH_j \cdot M_{ji}}{PM_i \cdot M_i} = \beta_{i,j} \left(\frac{PH_j}{PM_i}\right)^{1-\varepsilon_i^M}$$

**Key Insights:**

1. **When $\varepsilon_i^M = 1$ (Cobb-Douglas):** 
   $$s_{ji} = \beta_{i,j}$$
   The share parameters **equal** expenditure shares exactly.

2. **When prices are equal ($PH_j = PM_i$ for all $j$):**
   $$s_{ji} = \beta_{i,j}$$
   Again, parameters equal shares.

3. **General case:** The relationship depends on relative prices, but for **calibration purposes**, we typically assume steady-state prices are close to normalization, so:

$$\beta_{i,j} \approx s_{ji}^{\text{data}}$$

### Calibration Strategy

From the Input-Output table, we observe expenditure matrix $IO_{ij}$ where row $i$ = supplier, column $j$ = purchaser.

```matlab
% Read IO table: IO(supplier, purchaser) = expenditure
IO_matrix = readmatrix('IO_2021_chile.csv');

% Normalize by columns: each column sums to 1
beta_matrix = IO_matrix ./ sum(IO_matrix, 1);

% beta_matrix(i,j) = share of input i in sector j's intermediate bundle
% This matches the definition of beta_{j,i} in the model
% So we transpose to get beta_{i,j}:
modbeta = beta_matrix';  % Now modbeta(i,j) = beta_{i,j}
```

This is exactly what [main_SOE.m](main_SOE.m#L100-L103) does:
```matlab
betax  = betaio ./ sum(betaio,1);  % Column normalization
modbeta  = betax';  % Transpose
```

**Important:** The constraint is $\sum_{j=1}^{n} \beta_{i,j} = 1$ for each sector $i$.

---

## 3. Consumption Aggregator: Goods vs. Services (CES)

### Functional Form

Household preferences aggregate goods and services:

$$C = \left[ \omega^{1/\sigma} C_g^{(\sigma-1)/\sigma} + (1-\omega)^{1/\sigma} C_s^{(\sigma-1)/\sigma} \right]^{\sigma/(\sigma-1)}$$

Where:
- $C_g$ = goods consumption aggregate
- $C_s$ = services consumption aggregate  
- $\omega$ = weight parameter on goods
- $\sigma$ = elasticity of substitution (typically not estimated, often set to 1 for Cobb-Douglas)

**In your model:** The preferences are actually specified in **log form** (Cobb-Douglas):

$$\log C = \omega_g \log C_g + \omega_s \log C_s$$

Where $\omega_g + \omega_s = 1$.

### Expenditure Shares

For Cobb-Douglas preferences:

$$\frac{p_g \cdot C_g}{p_g \cdot C_g + p_s \cdot C_s} = \omega_g$$

So the **expenditure shares equal the weight parameters** exactly.

### Calibration Strategy

From consumption expenditure data:

```matlab
% From data:
spend_good = consumption_expenditure_on_goods_sectors;
spend_serv = consumption_expenditure_on_services_sectors;

% Calculate shares:
omega_g = sum(spend_good) / (sum(spend_good) + sum(spend_serv));
omega_s = 1 - omega_g;
```

This is in [main_SOE.m](main_SOE.m#L177):
```matlab
ombar_val = 0.57;  % Steady-state share of goods
```

---

## 4. Within-Category Consumption (Cobb-Douglas)

### Functional Form

Within goods and services, consumption follows Cobb-Douglas:

$$C_g = \prod_{i \in \text{goods}} C_{gi}^{\gamma_i^g}$$

$$C_s = \prod_{i \in \text{services}} C_{si}^{\gamma_i^s}$$

Where $\sum_{i \in \text{goods}} \gamma_i^g = 1$ and $\sum_{i \in \text{services}} \gamma_i^s = 1$.

### Expenditure Shares

By Cobb-Douglas properties:

$$\gamma_i^g = \frac{P_i \cdot C_{gi}}{p_g \cdot C_g}$$

The share parameters **equal** expenditure shares.

### Calibration Strategy

```matlab
% Identify sectors:
goods_sectors = (spend_good > spend_serv);
services_sectors = (spend_serv > spend_good);

% Calculate Cobb-Douglas weights:
gammag = spend_good ./ sum(spend_good);
gammas = spend_serv ./ sum(spend_serv);
```

This is in [main_SOE.m](main_SOE.m#L147-L152):
```matlab
goods    = spend_good > spend_serv;
services = spend_serv > spend_good;
gammag   = spend_good / sum(spend_good);
gammas   = spend_serv / sum(spend_serv);
```

---

## 5. Home vs. Foreign Varieties (CES)

### Functional Form

For each sector/good $i$, consumers choose between home and foreign varieties:

$$C_i = \left[ \varrho_i^{1/\sigma^H} (C^H_i)^{(\sigma^H-1)/\sigma^H} + (1-\varrho_i)^{1/\sigma^H} (C^F_i)^{(\sigma^H-1)/\sigma^H} \right]^{\sigma^H/(\sigma^H-1)}$$

Where:
- $C^H_i$ = home variety consumption
- $C^F_i$ = foreign variety consumption
- $\varrho_i$ = home bias parameter
- $\sigma^H$ = Armington elasticity (substitution between home/foreign)

### Expenditure Shares

The share of expenditure on home variety is:

$$s^H_i = \frac{p^H_i \cdot C^H_i}{P_i \cdot C_i} = \varrho_i \left(\frac{p^H_i}{P_i}\right)^{1-\sigma^H}$$

Where $P_i = \left[ \varrho_i (p^H_i)^{1-\sigma^H} + (1-\varrho_i) (p^F_i)^{1-\sigma^H} \right]^{1/(1-\sigma^H)}$ is the price index.

### Calibration Strategy

**This is endogenous!** The model jointly solves for:
- $\varrho_i$ (home bias parameters)
- Prices  
- Quantities

Subject to matching observed:
- Trade balance / GDP
- Home consumption shares by sector

This happens in the steady-state calibration routine `steady_ntwsoe_calib.m`.

---

## Summary Table: CES Parameters and Data Mapping

| CES Function | Parameter | Data Input | Mapping Formula | Notes |
|--------------|-----------|------------|-----------------|-------|
| **Production** | $\alpha_i$ | Intermediate cost share | $\alpha_i \approx \frac{\text{Intermediate costs}}{\text{Total costs}}$ | Direct from IO table |
| **Production** | $\alpha_i^V$ | Import cost share | $\alpha_i^V \approx \frac{\text{Import costs}}{\text{Total costs}}$ | Direct from trade data |
| **Intermediates** | $\beta_{i,j}$ | IO expenditure matrix | $\beta_{i,j} = \frac{IO_{ji}}{\sum_k IO_{kj}}$ | Column-normalized IO table |
| **Consumption (G vs S)** | $\omega_g$ | Goods expenditure share | $\omega_g = \frac{\text{Goods spending}}{\text{Total spending}}$ | Direct from consumption data |
| **Within Goods** | $\gamma_i^g$ | Sector spending (goods) | $\gamma_i^g = \frac{\text{Spending on good }i}{\text{Total goods spending}}$ | Cobb-Douglas weights |
| **Within Services** | $\gamma_i^s$ | Sector spending (services) | $\gamma_i^s = \frac{\text{Spending on service }i}{\text{Total services spending}}$ | Cobb-Douglas weights |
| **Home/Foreign** | $\varrho_i$ | Trade flows | Endogenously calibrated | Solved jointly with prices |

---

## Important Notes on Elasticities

### When does the mapping become exact?

1. **Cobb-Douglas case ($\varepsilon = 1$):** Parameters = Expenditure shares **exactly**, regardless of prices.

2. **Normalized prices:** If steady-state prices are normalized to 1 (or close), then for small elasticities ($\varepsilon$ close to 1), the approximation is very good:
   $$\text{parameter} \approx \text{expenditure share}$$

3. **General CES:** The relationship involves relative prices:
   $$\text{expenditure share} = \frac{\text{parameter} \cdot (\text{relative price})^{1-\varepsilon}}{\sum_j \text{parameter}_j \cdot (\text{relative price}_j)^{1-\varepsilon}}$$

### Practical Calibration

For your model with $\varepsilon^Y = 0.8$ and $\varepsilon^M = 0.1$:

- Both are **close to Cobb-Douglas** (1.0)
- In steady state with normalized prices, the approximation is excellent
- Therefore: **Set CES parameters equal to observed expenditure shares**

This is the standard approach in multi-sector DSGE models.

---

## Verification

To verify the calibration is correct, check in steady state:

```matlab
% After solving steady state:

% 1. Check production expenditure shares
intermediate_share_model = (PMi_ss .* M_ss) ./ (MCi_ss .* Yi_ss);
labor_share_model = (w_ss .* L_ss) ./ (MCi_ss .* Yi_ss);
import_share_model = (PV_ss .* Vi_ss) ./ (MCi_ss .* Yi_ss);

% Compare to parameters:
fprintf('Intermediate share: model=%.3f, param=%.3f\n', mean(intermediate_share_model), mean(alpha));

% 2. Check IO shares
for i = 1:nsec
    for j = 1:nsec
        io_share_model(i,j) = (pH_ss(j) * Mji_ss(j,i)) / (PMi_ss(i) * M_ss(i));
    end
end
% Compare to modbeta

% 3. Check consumption shares
goods_share_model = (p_g_ss * C_g_ss) / (p_g_ss * C_g_ss + p_s_ss * C_s_ss);
% Compare to ombar_val
```

---

## References

- See [main_SOE.m](main_SOE.m) for implementation
- See [steady_ntwsoe_calib.m](steady_ntwsoe_calib.m) for steady-state system
- See [NK_IOSOE_lev.mod](NK_IOSOE_lev.mod) for model equations
- See [IO_NOTATION_EXPLANATION.md](IO_NOTATION_EXPLANATION.md) for index conventions

---

## Common Pitfalls

1. **Index ordering:** Be careful with $(i,j)$ notation in IO tables vs. model code
2. **Normalization:** Make sure column sums equal 1 for $\beta$ matrix
3. **Elasticity values:** For $\varepsilon$ far from 1, the parameter ≠ expenditure share mapping breaks down
4. **Price normalization:** Verify steady-state prices are reasonably close to 1
5. **Completeness:** Ensure $\alpha_i + \alpha_i^V < 1$ (otherwise no room for labor!)

