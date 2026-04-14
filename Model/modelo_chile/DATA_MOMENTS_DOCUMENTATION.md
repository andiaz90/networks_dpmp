# Data Moments for SMM Estimation — NK-IOSOE Chile Model

**Script:** `compute_data_moments.m`  
**Output:** `data_moments_chile.mat` (struct `dm_chile`)  
**Estimation sample:** 2006Q1 – 2023Q4 (72 quarters)  
**Filtering:** All series HP-filtered at quarterly frequency with λ = 1600 before computing moments. This is required because the model is linearised around its steady state, so its solution describes cyclical (HP-filtered) deviations only.

---

## Moment Vector (15 elements)

| # | Moment | Variable | Description |
|---|--------|----------|-------------|
| 1 | `d_std_Yg`     | $\bar{\sigma}^Y_G$       | Output-weighted avg std dev of sectoral output, **goods** sectors (1–5) |
| 2 | `d_std_PHg`    | $\bar{\sigma}^{P_H}_G$   | Output-weighted avg std dev of sectoral price deflator, goods |
| 3 | `d_std_Lg`     | $\bar{\sigma}^L_G$       | Output-weighted avg std dev of sectoral employment, goods |
| 4 | `d_std_Ys`     | $\bar{\sigma}^Y_S$       | Output-weighted avg std dev of sectoral output, **services** sectors (6–12) |
| 5 | `d_std_PHs`    | $\bar{\sigma}^{P_H}_S$   | Output-weighted avg std dev of sectoral price deflator, services |
| 6 | `d_std_Ls`     | $\bar{\sigma}^L_S$       | Output-weighted avg std dev of sectoral employment, services |
| 7 | `0`            | $1 - \varrho(y^d, y^m)$  | Rank-correlation of **output** std devs (data vs model) — target = 0 (perfect rank order) |
| 8 | `0`            | $1 - \varrho(p^d, p^m)$  | Rank-correlation of **price** std devs (data vs model) — target = 0 |
| 9 | `0`            | $1 - \varrho(l^d, l^m)$  | Rank-correlation of **employment** std devs (data vs model) — target = 0 |
| 10 | `d_std_GDP`   | $\sigma^{GDP}$           | Std dev of aggregate real GDP cycle |
| 11 | `d_std_pi`    | $\sigma^\pi$             | Std dev of aggregate inflation (QoQ GDP deflator growth) |
| 12 | `d_corr_GDPpi`| $\text{corr}(GDP,\pi)$   | Correlation of aggregate output and inflation |
| 13 | `d_omG`       | $\bar\omega_G$           | Goods expenditure share (steady-state calibration target) |
| 14 | `d_std_Q`     | $\sigma^Q$               | Std dev of real exchange rate cycle |
| 15 | `d_TBGDP`     | $\overline{TB/GDP}$      | Average trade balance as share of nominal GDP |

Moments 7–9 are not directly loaded from data. The model's moment function (`smm_model_moments.m`) computes `1 − ρ(data_std, model_std)` for each of the three variable types, and the data target is always **0** (meaning data and model rank the sectors in the same volatility order).

---

## Data Sources and Construction

### 1. Sectoral Output — `y_d` [12×1]

| Item | Detail |
|------|--------|
| **File** | `pib_sectorial_bc.xlsx`, sheet **"Cuadro"** |
| **Producer** | Banco Central de Chile (BCCh) |
| **Frequency** | Quarterly |
| **Units** | Thousands of millions of chained 2018 pesos |
| **Coverage** | 2009Q1 – 2026Q2 (70 obs; sample window used: 2009Q1–2023Q4) |
| **Transformation** | Log → HP-filter (λ=1600) → `std()` |

**Column mapping** (xlsx column index, col 1 = date):

| Model sector | Name | xlsx cols used |
|---|---|---|
| 1 | Agriculture & Fishing | col 2 (Agropecuario-silvícola) **+** col 3 (Pesca) — summed |
| 2 | Mining | col 4 (Minería total aggregate) |
| 3 | Manufacturing | col 7 (Industria Manufacturera aggregate) |
| 4 | Utilities | col 17 (Electricidad, gas y agua) |
| 5 | Construction | col 18 (Construcción) |
| 6 | Trade/Hotels | col 19 (Comercio, restaurantes y hoteles aggregate) |
| 7 | Transport/Comm | col 22 (Transporte) **+** col 23 (Comunicaciones e info) — summed |
| 8 | Finance | col 25 (Servicios financieros) |
| 9 | Real Estate | col 27 (Servicios de vivienda e inmobiliarios) |
| 10 | Business Services | col 26 (Servicios empresariales) |
| 11 | Personal Services | col 28 (Servicios personales) |
| 12 | Public Admin | col 29 (Administración pública) |

Sub-sector columns (e.g., copper mining, manufacturing sub-industries) are skipped; aggregate totals are used to avoid double-counting.

---

### 2. Sectoral Price Deflators — `p_d` [12×1]

| Item | Detail |
|------|--------|
| **File** | `deflactor_pib.csv` |
| **Producer** | Banco Central de Chile (BCCh) |
| **Frequency** | Quarterly |
| **Units** | Price index (2018 = 100) |
| **Coverage** | 1996Q1 – 2023Q4 (112 obs) |
| **Date format** | `MMMYYYY` — MAR=Q1, JUN=Q2, SEP=Q3, DIC=Q4 |
| **Transformation** | Log → HP-filter (λ=1600) → `std()`; isolated NaNs filled by linear interpolation before filtering |

The CSV has 28 series in columns 2–29. Mapping to 12 model sectors:

| Model sector | deflactor_pib cols used | Rule |
|---|---|---|
| 1 | 1 (Agropecuario) + 2 (Pesca) | simple mean |
| 2 | 3 (Minería total) | direct |
| 3 | 6 (Manufactura total) | direct |
| 4 | 15 (Electricidad, gas, agua) | direct |
| 5 | 16 (Construcción) | direct |
| 6 | 17 (Comercio) + 18 (Restaurantes y hoteles) | simple mean |
| 7 | 19 (Transporte) + 20 (Comunicaciones) | simple mean |
| 8 | 21 (Servicios financieros) | direct |
| 9 | 23 (Servicios de vivienda) | direct |
| 10 | 22 (Servicios empresariales) | direct |
| 11 | 24 (Servicios personales) | direct |
| 12 | 25 (Administración pública) | direct |

Column 28 (Deflactor del PIB total) is used for aggregate inflation (moments 11–12).

---

### 3. Sectoral Employment — `l_d` [12×1]

| Item | Detail |
|------|--------|
| **File** | `count_workers_by_sector.csv` |
| **Producer** | INE / BCCh (derived from ENE/NENE survey) |
| **Frequency** | Monthly |
| **Units** | Worker headcount |
| **Coverage** | 2006M1 onwards |
| **Transformation** | Monthly → quarterly (mean of 3 months per quarter) → log → HP-filter (λ=1600) → `std()` |

The CSV has 12 sector columns in **alphabetical order**. They are reordered to model sector order using:

| CSV column | Sector name (English) | Model sector |
|---|---|---|
| 2 | Agriculture, forestry and fishing | 1 |
| 3 | Business services | 10 |
| 4 | Construction | 5 |
| 5 | Electricity, gas, water and waste management | 4 |
| 6 | Financial intermediation | 8 |
| 7 | Manufacturing | 3 |
| 8 | Mining | 2 |
| 9 | Personal services | 11 |
| 10 | Public administration | 12 |
| 11 | Real estate and housing services | 9 |
| 12 | Transport, communications and information services | 7 |
| 13 | Wholesale and retail trade; accommodation | 6 |

---

### 4. Cross-Sectional Averages — Moments 1–6

The sectoral `y_d`, `p_d`, `l_d` vectors are aggregated to goods/services scalar moments using **steady-state output weights** loaded from `params_val.mat` (field `Yi_ss`):

$$\bar\sigma^Y_G = \sum_{i \in G} w^G_i \, \sigma^Y_i, \qquad w^G_i = \frac{Y^{ss}_i}{\sum_{j\in G} Y^{ss}_j}$$

and analogously for services and for price/employment. Equal weights are used as a fallback if `params_val.mat` is not found.

---

### 5. Aggregate Output — `d_std_GDP` (Moment 10)

| Item | Detail |
|------|--------|
| **File** | `pib_sectorial_bc.xlsx`, sheet "Cuadro", **column 32** (PIB total) |
| **Transformation** | Log → HP-filter (λ=1600) → `std()` |
| **Note** | This is the BCCh official aggregate GDP series, not the sum of sectoral rows, to avoid possible rounding discrepancies |

---

### 6. Aggregate Inflation — `d_std_pi` (Moment 11) and `d_corr_GDPpi` (Moment 12)

| Item | Detail |
|------|--------|
| **File** | `deflactor_pib.csv`, column 28 (Deflactor del PIB total) |
| **Construction** | $\pi_t = P_t / P_{t-1}$ (quarter-on-quarter gross inflation rate) |
| **Transformation** | log → HP-filter (λ=1600) → `std()` for moment 11; `corr(GDP_{hp}, \pi_{hp})` for moment 12 |
| **Alignment** | GDP and π series both trimmed to the common estimation sample end before computing correlation |

---

### 7. Real Exchange Rate — `d_std_Q` (Moment 14)

| Item | Detail |
|------|--------|
| **File** | `reer_chile_bis.xlsx`, sheet **"Real"** |
| **Producer** | Bank for International Settlements (BIS) |
| **Series** | **RBCL** — Real Broad REER for Chile (CPI-based, 2020 = 100) |
| **URL** | https://www.bis.org/statistics/eer/broad.xlsx |
| **Frequency** | Monthly (column K in the "Real" sheet, after date column) |
| **Coverage** | ~1994M1 onwards |
| **Transformation** | Monthly → quarterly average → log → HP-filter (λ=1600) → `std()` |
| **Fallback** | 0.052 (literature estimate for Chile) if file read fails |

---

### 8. Trade Balance-to-GDP Ratio — `d_TBGDP` (Moment 15)

| Item | Detail |
|------|--------|
| **File** | `datos_CCNN_mayo2025.xlsx`, sheet **"Nom trimestral"** |
| **Producer** | Banco Central de Chile (BCCh), National Accounts |
| **Columns used** | Col I: Exports of goods and services (`6. Exportaciones de Bienes y Servicios`); Col J: Imports (`7. Importaciones`); Col K: Nominal GDP (`8. PIB a precios corrientes`) |
| **Units** | Thousands of millions of current pesos |
| **Construction** | $\overline{TB/GDP} = \text{mean}\left(\frac{X_t - M_t}{GDP_t}\right)$ over the estimation sample — a **level** moment (mean, not std dev) |
| **Fallback** | −0.02 (historical average for Chile) if file read fails |

---

### 9. Goods Expenditure Share — `d_omG` (Moment 13)

| Item | Detail |
|------|--------|
| **Value** | 0.57 (hardcoded) |
| **Source** | SS calibration target (`ombar_val` in `params_val.mat`) derived from BCCh National Accounts consumption breakdown (goods = agricultural + manufacturing + utilities + construction share of household expenditure) |
| **Note** | This is a calibrated steady-state value, not an HP-filtered second moment. It is included in the moment vector to pin the goods/services expenditure split. |

---

## Sector Classification

| Model # | Name | Good / Service |
|---|---|---|
| 1 | Agriculture & Fishing | Good |
| 2 | Mining | Good |
| 3 | Manufacturing | Good |
| 4 | Utilities (Electricity, Gas, Water) | Good |
| 5 | Construction | Good |
| 6 | Trade, Restaurants & Hotels | Service |
| 7 | Transport & Communications | Service |
| 8 | Financial Services | Service |
| 9 | Real Estate & Housing | Service |
| 10 | Business Services | Service |
| 11 | Personal Services | Service |
| 12 | Public Administration | Service |

---

## Data Files Summary

| File | Location | Used for |
|---|---|---|
| `count_workers_by_sector.csv` | `data_github/` | Employment moments (3, 6, 9) |
| `deflactor_pib.csv` | `data_github/` | Price moments (2, 5, 8) and inflation moments (11, 12) |
| `pib_sectorial_bc.xlsx` | `data_github/` | Output moments (1, 4, 7) and aggregate GDP (10) |
| `reer_chile_bis.xlsx` | `data_github/` | Real exchange rate moment (14) |
| `datos_CCNN_mayo2025.xlsx` | `data_github/` | Trade balance moment (15) |
| `params_val.mat` | `Model/modelo_chile/` | Steady-state weights for cross-sectional averages and goods share (13) |
