# NK-IOSOE Calibration Audit — 19 Aug 2026

Audit of the baseline run `julia Model/julia_dynare/main_SOE_gap.jl`
(θ from `estimation_results/smm_checkpoint.csv`, written 2026-07-24 10:21).

Every number below is traced to a file. Sources: `Data/2021_Cuadros_12x12.xlsx`
(Cuadro 23, valor agregado), `Data/sector_calibration.csv`,
`Model/julia_dynare/utils.jl`, `smm_estimation.jl`, `estimation_results/*`.

---

## Summary

Three things are blocking, in the sense that no result in the current run is
quotable until they are fixed:

1. The reported θ is the CMA-ES **starting point**, not an estimate.
2. The objective **does not reproduce** across the two code paths (21.47 vs 15.09).
3. There is **no inference** — all standard errors are exactly zero.

Beyond that, one calibration choice has first-order consequences for the paper's
central claim: the intermediate-input shares are built on a denominator that
excludes capital income, which inflates the network multiplier by roughly 70%.

---

## A. Blocking

### A1. The parameter vector was never estimated

`smm_estimation.jl:801` restricts the search to seven dimensions:

```julia
_FREE = [1, 4, 5, 6, 33, 34, 36]   # ilabcosts, log κ_V, ρ_om, ρ_A, ρ_ξ, σ_ξ, κ_w
```

All seven sit **exactly** at their initial values in the checkpoint:

| Param | Checkpoint | Reads as |
|---|---|---|
| `ilabcosts` | 1.0 | round default |
| `log(kappaV)` | 13.815510557964274 | = log(10⁶) exactly |
| `rho_om` | 0.5 | round default |
| `rho_A` | 0.5 | round default |
| `rho_xi` | 0.7 | round default |
| `sigma_xi` | 0.01 | round default |
| `kappaw` | 100.0 | round default |

`estimation_results/smm_progress_log.csv` has two rows: eval 50 at 9.5 s and
eval 100 at 19.7 s, with `best_obj` frozen at 21.471848 across both, and 6
failed evaluations in the second batch. The run lasted 20 seconds out of a
300,000-evaluation budget.

The remaining 29 parameters are pinned by design (elasticities at Atalay, 24
sectoral shock sizes measured from data, external ρ/σ from
`external_shock_calibration.csv`). That is the intended redesign. But it means
**the current model is calibrated at a guess in every free dimension**, and the
`main_SOE_gap.jl` header "Parameters (smm_checkpoint.csv)" reads as if these
were estimated. Relabel the printout or re-run.

### A2. The objective is not reproducible — **RESOLVED 2026-08-19**

> **Diagnosis: there is no live disagreement between the two code paths. The
> 21.47 is a stale number computed on a model that no longer exists.**
>
> File timeline for 2026-07-24:
>
> | Time | Event |
> |---|---|
> | **10:21** | `smm_checkpoint.csv` written, **obj = 21.4718** |
> | 10:48 | `NK_SOE_lev_gap2.mod` → **`eps_pc` added** (world copper price) |
> | 11:16 | `.mod` → `choice3` |
> | 13:23 | `.mod` → `ownership` |
> | 13:40 | `.mod` → **`GDP_vol` added** (volume GDP) |
> | 14:11–14:12 | `smm_model_moments.jl`, `smm_estimation.jl` updated |
>
> Verified by grepping the dated `.mod` backups:
>
> | file | `eps_pc` | `GDP_vol` |
> |---|---:|---:|
> | `NK_SOE_lev_gap2.BACKUP_precopper.mod` (model as of 10:21) | 0 | 0 |
> | `NK_SOE_lev_gap2.mod` (current) | 4 | 6 |
>
> So the run that produced 21.47 had **27 active shocks, not 28** — no copper
> shock — and computed its aggregate and labour moments on nominal `GDP` rather
> than `GDP_vol`. `GDP_vol` enters 5 of the 60 moments: `std(GDP)`,
> `corr(GDP,π)`, `corr(GDP,Q)`, `corr(N,GDP)`, `corr(N,GDP/N)`. That is exactly
> the set of blocks that moved.
>
> Two things ruled out:
> - **Data moments are not the cause.** `sectoral_moments.csv` and
>   `aggregate_moments.csv` are both dated 2026-07-23 14:24 — they predate the
>   checkpoint and have not changed since.
> - **The weighting matrix is not the cause.**
>   `diff` of `build_weighting_matrix` between `utils.BACKUP_precopper.jl`
>   (10:50) and the current `utils.jl` is empty.
>
> The estimator and `main_SOE_gap.jl` compute the 60 moments with line-for-line
> identical formulas (`smm_estimation.jl:456-633` vs `main_SOE_gap.jl:1197-1345`),
> activate the same 28 shocks by name, and read the same steady state. The only
> remaining difference is Klein vs Dynare QZ, which should be numerical noise.
>
> **A latent trap was also found and closed.** `main_SOE_gap.jl` can substitute
> the estimator's own stored moments from `smm_results.csv` for the fit table.
> That file currently has a correct 60 rows and a `model` column, so the
> existing row-count guard passes — but its contents were computed at 10:21 on
> the pre-copper model (it reports `std(Y_1)` = 0.0301 against today's 0.0181).
> Only the `smm_param_source == "smm_estimates.csv"` condition prevented those
> stale moments from being printed as this run's fit. Row count is necessary,
> not sufficient.
>
> **Fix applied.** A provenance stamp, dependency-free (no new `using`):
> - `utils.jl` — new `objective_dep_files`, `objective_fingerprint`,
>   `write_objective_provenance`, `check_objective_provenance`. The fingerprint
>   is `basename:bytes:mtime` over the `.mod`, `utils.jl`, `smm_estimation.jl`,
>   `smm_model_moments.jl` and the two data-moment CSVs.
> - `smm_estimation.jl` — `save_checkpoint` now writes
>   `estimation_results/objective_provenance.txt` on every checkpoint.
> - `main_SOE_gap.jl` — on load, compares the stamp (falling back to mtimes for
>   pre-existing checkpoints) and prints either
>   `stored obj = X — comparable` or
>   `stored obj = X ** NOT COMPARABLE **` plus the list of changed files. The
>   same flag now also gates the `smm_results.csv` substitution.
>
> Against the current repo the guard correctly returns
> `comparable = false, changed = [NK_SOE_lev_gap2.mod, utils.jl,
> smm_estimation.jl, smm_model_moments.jl]` and correctly does *not* flag the
> data moments.
>
> **Still to verify:** that the two paths agree *today*. Run a short estimation
> on the current model (it writes `smm_estimates.csv`, `smm_results.csv` and the
> new provenance stamp), then run `main_SOE_gap.jl`. The provenance check will
> pass, `main` will load the estimator's moments, and the existing
> "Dynare (QZ) diagnostic comparison" block will print both objectives and their
> difference. Anything beyond ~1e-6 is a real Klein/QZ discrepancy worth chasing;
> agreement closes A2 completely.

---

#### Original finding (superseded by the diagnosis above)

Same θ, two code paths, different answers — in every block:

| Block | Estimator (`smm_progress_log.csv`) | `main_SOE_gap.jl` | Ratio |
|---|---:|---:|---:|
| Y | 0.280 | 0.488 | 1.74 |
| PH | 4.180 | 1.742 | 0.42 |
| L | 0.812 | 0.866 | 1.07 |
| Agg | 1.313 | 0.837 | 0.64 |
| Rank | 8.639 | 4.777 | 0.55 |
| CorrYP | 4.285 | 3.926 | 0.92 |
| NLab | 1.963 | 2.456 | 1.25 |
| **Total** | **21.472** | **15.091** | **0.70** |

`smm_estimation.jl:165-169` explicitly claims the two share one objective
("moved 2026-07-10 … so the printed fit is the SAME objective the estimator
minimizes"). It is not. Suspects, in order:

- **Solution method.** The estimator solves via Klein from the g1 Jacobian
  (`klein_hit_pct` in the log); `main` reports "Moment fit (Dynare (QZ))".
- **Active shock set.** `main` prints "Active shocks (by name): 28 of 30
  exogenous". Confirm the estimator activates the same 28.
- **HP-filter path** on the theoretical spectrum.

Diagnose by evaluating both paths on the same θ with a single shock active and
comparing `Γ` element by element.

### A3. No standard errors

`estimation_results/smm_inference.csv`: every `std_err` is `0.0` and every
`t_stat` is `Inf`, for all 36 parameters. The numerical Jacobian is degenerate.
Expected for the 29 pinned dimensions, but it also fails for the 7 free ones,
which means the step size or the finite-difference loop is broken. No
overidentification (J) statistic is reported anywhere.

---

## B. Calibration — first order

### B1. Intermediate shares drop capital income — **RESOLVED 2026-08-19**

> **Resolution: adopt the FGI convention. `sector_calibration.csv` regenerated.**
>
> First, a citation correction. Two different papers were being conflated:
>
> | file | actual paper |
> |---|---|
> | `CF_KeynesianSupply_latest.pdf` | Cesà-Bianchi & Ferrero (2021), *The Transmission of Keynesian Supply Shocks* |
> | `FGI_jme_2023.pdf` | **Ferrante, Graves & Iacoviello (2023, JME)**, *The inflationary effects of sectoral reallocation* |
>
> The bibliography maps `ferrante2023inflationary` to the Cesà-Bianchi–Ferrero
> PDF. **Fix this before submission** — CF has no Ferrante on it.
>
> FGI is the actual template: same CES nest (ε_Y between labour and
> intermediates, ε_M across intermediates), same Γ_ij normalised to unit row
> sums, same κ_i converted from Pastén et al., same γᵢᵍ/γᵢˢ split, and Table 1
> carries γ = 2, χ̄ = 1, ψ = 1, **ε = 10** — precisely the NK-IOSOE "from
> literature" row. They also *estimate* ε_M = 0.13 and ε_Y = 0.82, which is
> essentially what this project pins (0.20 / 0.80, attributed to Atalay).
>
> **FGI Table 1:** `Intermediate Input Share (Range)  α_i  0.11 to 0.83  BEA`,
> with the note *"Industries with lowest and highest values of α_i are 'Housing'
> and 'Funds, Trusts, and Other Financial Vehicles.'"*
>
> Housing lowest is a clean falsification test, because housing is where the two
> conventions diverge most — imputed rents are almost all operating surplus with
> a near-zero wage bill:
>
> | Servicios inmobiliarios | value | rank of 12 |
> |---|---:|---|
> | revenue-based `CI/GO` | 0.215 | **1st — lowest**, matches FGI |
> | cost-based `CI/(CI+REM)` | 0.887 | **12th — highest** |
>
> So FGI use α_i = intermediate inputs / **gross output**, with the residual
> (1 − α_i) = VA/GO going entirely to labour. They have no capital — *"as our
> model does not include capital, these hiring costs capture a variety of
> frictions affecting a firm's ability to expand its productive capacity"* —
> accept the resulting high labour share, and keep ε = 10.
>
> **The inconsistency in one sentence:** the NK-IOSOE imported FGI's model
> structure, FGI's ε = 10 and FGI's elasticities, but Rubbo's α convention.
> Rubbo's convention is right for Phillips-curve sufficient statistics with no
> steady state reported; it is the wrong pairing with ε = 10 and a calibrated,
> reported steady state.
>
> **Applied.** `Data/build_sector_calibration.py` (new) regenerates α and α_V
> from `2021_Cuadros_12x12.xlsx` (Cuadro 23) and `2021_MIP_12x12.xlsx` (sheet 1,
> domestic activity × activity), with the convention on one visible line
> (`DENOMINATOR = "GO"`; set `"CI+REM"` to reproduce the old file). It rewrites
> only α and α_V — the other five columns are carried through unchanged because
> their provenance has not been reconstructed. Old file backed up to
> `sector_calibration.csv.BACKUP_precostbase_20260819`.
>
> Result:
>
> | | intermediate share (GO-wtd) | implied GO/VA |
> |---|---:|---:|
> | was | 0.700 | 3.34 |
> | now | **0.480** | **1.92** |
> | Chilean accounts | 0.480 | 1.92 |
>
> FGI check passes: lowest α is Servicios inmobiliarios at 0.215; range 0.22–0.72
> against FGI's 0.11–0.83 across 66 US sectors.
>
> Cost of the convention, stated plainly for the paper: labour now absorbs 100%
> of value added against 41.6% in the Chilean accounts. That is inherent to a
> no-capital model and is what FGI accept. A fixed factor with weight EBE/GO
> would fix it and is worth one robustness paragraph — Rubbo's caveat about
> "rental payments to fixed factors" bites harder in Chile than the US, since
> mining is 24.8% of EBE and real estate a further 13.4%.
>
> **Everything downstream moves**: the steady state, `gdp_va_shares_baseline.csv`,
> `vae_baseline.csv`, every amplification number, and θ itself (the estimates
> were obtained under the old shares). Re-run `main_SOE_gap.jl`, then re-estimate.
>
> **B1b — the network matrix was the total, not the domestic one. FIXED.**
> `modbeta` was built from `IO_2021_chile.csv`, whose column sums equal Cuadro 23
> *total* intermediate consumption, while `main_SOE_gap.jl:876` asserted "modbeta
> is the domestic IO matrix" and the model treats imports as a separate bundle
> V_i. So imported inputs were booked to domestic suppliers.
>
> I first called this second-order because Γ is row-normalised. **That was wrong** —
> it is a 40% change:
>
> ```
> ||Γ_dom − Γ_tot||_F / ||Γ_tot||_F = 0.397     max cell change 0.301
> ```
>
> | purchaser ← supplier | total | domestic |
> |---|---:|---:|
> | Utilities ← Mining | 0.314 | **0.013** |
> | Serv. personales ← Manufactura | 0.433 | 0.205 |
> | Agro ← Agro | 0.139 | 0.362 |
> | Construcción ← Manufactura | 0.562 | 0.365 |
>
> Utilities←Mining is the clearest case: Chile *imports* the fuel it burns, and
> the total matrix credited that to domestic mining. That single cell would
> distort every oil and mining shock exercise in the paper. And now that α_i is
> explicitly the *domestic* intermediate share, the composition of that bundle
> has to be domestic for internal consistency.
>
> A domestic matrix already existed as `mip_12x12.csv` but was written with
> European thousands separators ("3.541" = 3541), so the loader could never read
> it. `build_sector_calibration.py` now emits `IO_2021_chile_domestic.csv` in the
> loader's format (`;`-separated, no header, rows = supplier, cols = purchaser),
> with a row-total integrity check against the workbook.
>
> `main_SOE_gap.jl` and all five `*_shock_analysis.jl` scripts now prefer the
> domestic file and fall back to the total one with a loud warning.
> `smoke_test.jl` checks for it.

---

### B1c. Other flags cleared before the re-run (2026-08-19)

| item | was | now |
|---|---|---|
| `ygdp_target` (`main_SOE_gap.jl:743`) | `2.0`, a round "IMF target" placeholder | `1.924`, from Cuadro 23 — and only now a meaningful check, since under the old cost base the model ratio was 3.29 and could never match |
| `SECTOR_VA_SHARE` (`utils.jl:478`) | hardcoded from an old model run under the *old* α convention — stale twice over | Cuadro 23 value-added shares |
| weight-matrix comment | justified the design with "Public Admin (0.17% of VA)" — a model artifact | corrected; Cuadro 23 puts public admin at 5.08%. The 1/d² argument stands on its own |
| provenance dependency list | model + moment code + data moments only | also `sector_calibration.csv`, both IO matrices, the fpa vector and both shock-calibration files — all of which move the steady state and hence every model moment |

The last one matters: without it, today's α change would have silently
invalidated θ with no warning, which is the exact failure mode A2 was about.

---

#### B1 background (framing superseded twice; arithmetic throughout is correct)

#### B1-mid. Reframing against Rubbo / Afrouzi–Bhattarai

> **Correction.** The original text below called this convention "the worst of
> the options" and "indefensible in a replication package." That was wrong.
> Netting non-labour value added out of costs is the *standard* convention in
> exactly the production-network literature this paper builds on:
>
> - **Afrouzi & Bhattarai (2024)**, Data Appendix E.1, explicitly:
>   `Total Industry Cost_j = Total Intermediate_j + Compensation of Employees_j`
>   — gross operating surplus and production taxes are excluded from cost and
>   carried as a revenue/cost wedge `(1+ω_j)`.
> - **Rubbo (2023, Econometrica)**, §6.1: "we net out the non-labor component of
>   value added from total costs in the computation of labor and input shares."
>   She flags the caveat directly: this is innocuous if non-labour VA is profit,
>   but "if instead nonlabor value added also includes rental payments to fixed
>   factors (such as capital and land), our calibration is an upper bound."
> - **Ferrante et al. (2023)**, the paper this project follows most closely, uses
>   `Y_kt(j) = e^{a_kt} Z_kt(j)^{α_k} N_kt(j)^{1−α_k}` — intermediates and labour,
>   **no capital** — with α_k taken from Pastén, Schoenle & Weber (2020).
>
> So `α + α_V = CI/(CI+REM)` is the cost-based input–output matrix that the
> theory actually calls for. It is not an error.
>
> **The real inconsistency is the markup, and it is sharper.** Those papers carry
> the wedge explicitly and let the data set it. The NK-IOSOE imposes ε = 10
> everywhere, i.e. an 11.1% markup and a wedge of 1.111. Under the same
> convention the Chilean data imply:
>
> | | 1+ω = GO/(CI+REM) | markup | implied ε |
> |---|---:|---:|---:|
> | **Aggregate** | **1.435** | **43.5%** | **3.30** |
> | Servicios inmobiliarios | 4.124 | 312% | 1.32 |
> | Minería | 2.379 | 138% | 1.72 |
> | Utilities | 1.574 | 57% | 2.74 |
> | Manufactura | 1.214 | 21% | 5.68 |
> | Intermediación financiera | 1.158 | 16% | 7.34 |
>
> The model's ε = 10 markup accounts for 11.1% of gross output where the
> calibration implies 43.5% — leaving **22.6% of gross output** as revenue with
> no counterpart in the model. That, not a data-construction mistake, is what
> produces gross output/GDP = 3.29 against 1.92 in the accounts.
>
> **Rubbo's caveat binds hard in Chile.** Mining is 24.8% of Chilean EBE and
> real estate a further 13.4% — 38% of operating surplus between them. In both
> the surplus is a rent to a natural resource or a housing stock, i.e. exactly
> the "rental payments to fixed factors" case in which she says the calibration
> becomes an upper bound. The US papers can more comfortably read EBE as profit;
> Chile cannot.
>
> **Two objects, not one.** The cost-based ratio is the right input to the
> sectoral Phillips curves and to pass-through. It is *not* the national-accounts
> gross-output-to-GDP ratio. The `IMF target: 2.0` check printed by
> `main_SOE_gap.jl` compares the two and will always fail; either drop it or
> relabel it. The same caveat propagates to `tables/vae_baseline.csv` and to any
> amplification claim quoted from the steady state.
>
> Revised options are in the "Suggested order of work" section.

---

#### Original finding (framing superseded; the arithmetic below is still correct)

#### B1-old. Intermediate shares drop capital income; the network multiplier is ~70% too large

`sector_calibration.csv` sets αᵢ + α_Vᵢ = CIᵢ / (CIᵢ + Remuneracionesᵢ), i.e.
**excedente bruto de explotación is deleted from the denominator** and its mass
redistributed pro rata to materials and labor. Verified against Cuadro 23:

- corr(α + α_V, CI/(CI+REM)) = **0.9998**, mean absolute error 0.008
- corr(α + α_V, CI/(CI+VA)) = 0.496 (the data concept)
- EBE is **56.8%** of Chilean value added in 2021

Consequence:

| | Intermediate share (GO-weighted) | Implied gross output / VA |
|---|---:|---:|
| Chilean IO 2021 | 0.480 | **1.92** |
| `sector_calibration.csv` | 0.700 | **3.34** |
| Model prints | — | 3.29 (nominal) / 4.49 (unweighted) |

1/(1 − 0.700) = 3.34 accounts for the printed 3.29 almost exactly. The "IMF
target 2.0" flag in the run output is not a tuning miss — it is this identity.

Per-sector, worst first:

| Sector | CI/GO data | α+α_V calib | gap |
|---|---:|---:|---:|
| Inmobiliario | 0.215 | 0.882 | +0.667 |
| Minería | 0.355 | 0.847 | +0.491 |
| Utilities | 0.542 | 0.853 | +0.311 |
| Agro | 0.576 | 0.825 | +0.250 |
| Transporte | 0.522 | 0.717 | +0.195 |

**Why this matters more than the other items.** The paper's contribution is
network amplification. The Leontief multiplier is a direct function of the
intermediate share. At 0.70 instead of 0.48 the model's Domar weights, its
sectoral-to-aggregate pass-through, and every VAE number in
`tables/vae_baseline.csv` are inflated by a factor that has nothing to do with
Chilean production structure. A referee at any of the target journals will
compute gross output / GDP from the calibration table in ten seconds.

**This is a modeling choice, not a coding bug** — the model has no capital, so
EBE must go somewhere. But pro-rata is the worst of the options:

- *(a)* Reassign all EBE to labor: α = CI/GO, labor share = 1 − CI/GO. Preserves
  the data intermediate share and GO/VA = 1.92 exactly. Cheapest fix; costs an
  implausibly high labor share (0.52 aggregate vs 0.42 in the data).
- *(b)* Add a fixed factor in the CES with share EBE/GO. Preserves both the
  intermediate share and the labor share; costs one more term in the production
  function and marginal cost.
- *(c)* Absorb part of EBE into pure profit. With ε = 10 the markup delivers only
  ~10% of revenue, so this covers at most a fifth of EBE. Not sufficient alone.

(b) is the right answer for a top-journal submission; (a) is the right answer if
you need the number fixed this week.

### B2. Sectoral VA shares are far from Chilean data

Model vs Cuadro 23, % of total VA:

| Sector | Data | Model | Gap |
|---|---:|---:|---:|
| Manufactura | 9.54 | 20.74 | **+11.20** |
| Serv. empresariales | 9.75 | 21.72 | **+11.97** |
| Comercio, hoteles, rest. | 12.66 | 3.01 | **−9.65** |
| Inmobiliario | 8.55 | 2.90 | −5.65 |
| Adm. pública | 5.08 | 0.14 | −4.94 |
| Serv. personales | 13.27 | 8.53 | −4.74 |
| Construcción | 6.45 | 3.01 | −3.44 |
| Financiera | 3.91 | 6.07 | +2.16 |
| Transporte | 8.31 | 10.30 | +1.99 |
| Agro | 3.92 | 5.64 | +1.72 |
| Minería | 15.78 | 14.56 | −1.22 |
| Utilities | 2.79 | 3.37 | +0.58 |

Public administration is effectively absent from the model: Y_ss = 0.0022, or
0.06% of gross output against 3.63% in the data.

Partly downstream of B1, partly the consumption vector. `sector_calibration.csv`
allocates **49.6% of all household consumption to manufacturing**, and ~0 to
mining, construction and public administration (`spend_good` = 1.086e-8 for
sectors 2 and 5; `spend_serv` = 0.0011 for sector 12). Comercio, which is the
largest consumption category in the Chilean accounts, gets 5.2%. Check whether
the consumption vector was built on an activity-by-activity table without
reallocating trade and transport margins.

### B3. The weight matrix uses model VA shares, not data

`utils.jl:478`:

```julia
const SECTOR_VA_SHARE = [5.07, 12.63, 15.87, 2.90, 3.14, 3.34,
                         10.93, 6.71, 2.13, 26.03, 11.08, 0.17]
```

These are an earlier *model* run, not BCCh. Ratio of correct to current weight:

| Sector | current | data | ×  |
|---|---:|---:|---:|
| Adm. pública | 0.17% | 5.08% | **×29.9** |
| Inmobiliario | 2.13% | 8.55% | ×4.0 |
| Comercio | 3.34% | 12.66% | ×3.8 |
| Construcción | 3.14% | 6.45% | ×2.1 |
| Serv. empresariales | 26.03% | 9.75% | ×0.37 |
| Manufactura | 15.87% | 9.54% | ×0.60 |

The estimator is systematically downweighting the sectors the model gets most
wrong — the comment at `utils.jl:488` even justifies the choice by citing
"Public Admin (0.17% of VA)", treating a model artifact as the fact.

Honest scoring: correcting this moves the total objective only 15.09 → 15.13, so
it is a correctness fix, not a fit fix. But it changes *which* sectors the
optimizer chases, and it is indefensible in a replication package.

### B4. No investment, no government

C/GDP = 97.7% in steady state (GDP = C + TB). Chile: C ≈ 62%, I ≈ 24%, G ≈ 15%.

Structural, not a calibration error — but it is why construction and public
administration have no demand base, and it removes the investment channel
through which sectoral and monetary shocks normally propagate in a production-
network model. This needs an explicit paragraph in the paper, not silence.

### B5. Foreign debt position

Q·B*/GDP = 161% of quarterly GDP = **40% of annual GDP**. Chile's gross external
debt ≈ 70% of GDP; net IIP ≈ −20%. In a one-asset model the net concept is the
right target, so the calibration is roughly 2× too levered. Pick the target
explicitly and state it in the paper.

---

## C. Fit diagnostics

Conditional on A2 being resolved — the numbers below are from the `main` path.

### C1. Labor is the model's largest failure (≈37% of the objective)

| Moment | Data | Model |
|---|---:|---:|
| corr(N, GDP) | 0.860 | **0.110** |
| corr(N, GDP/N) | 0.286 | **−0.530** |
| rank corr, sectoral labor vol | 1.000 | **−0.252** |

The rank correlation is *negative*: the model orders sectors by employment
volatility backwards. That single moment is 20.8% of the objective; the two
comovement moments add 16.3%.

The proximate cause is visible in the sectoral table: the model puts
std(L₁₁) = 0.093 in Servicios personales against 0.020 in the data (the largest
miss anywhere in the fit), and only 0.014 in Minería against 0.045. σ_tfp,₁₁ =
0.0737 is the largest measured shock in the table, loaded onto the least
volatile large sector. Check `compute_sectoral_shocks.jl` for sector 11 — this
looks like a measurement artifact, not a fit failure.

Under GHH with κ_w = 100 (never estimated — see A1) labor supply responds to the
wage alone. Getting corr(N, GDP) from 0.11 to 0.86 is the single highest-return
modeling change available.

### C2. Sectoral prices are far too smooth

| | data | model |
|---|---:|---:|
| std(PH₁) Agro | 0.128 | 0.018 |
| std(PH₄) Utilities | 0.127 | 0.021 |
| std(PH₃) Manuf | 0.053 | 0.014 |
| std(PH₂) Minería | 0.167 | 0.109 |

Only mining tracks, and only because it receives `eps_pvstar` directly. Every
other sector's price volatility collapses toward the aggregate — the model has
no sector-specific price driver of the right magnitude.

Note κ_V = 10⁶ (import price adjustment). At that value import prices are
effectively fixed and exchange-rate pass-through is switched off, which is hard
to square with citing Romero (2025 JIE) as the pass-through anchor. The skill
notes the paper reports κ_V = 10¹³ — reconcile the two before the calibration
table goes in.

### C3. The preference-shock block does nothing

From the IRF table, impact on GDP: `eps_om_2` = −0.0000, `eps_om_5` = +0.0000,
largest is `eps_om_7` at −0.065%. Twelve σ_om parameters are being measured from
data and fed into a model where they have no traction.

This is also the likely reason the corr(Yᵢ, PHᵢ) block (26% of the objective)
fails: model correlations span −0.26 to +0.05, data span −0.88 to +0.48. Without
a working demand shock there is no cross-sector variation in the supply/demand
mix for the model to match. Either the preference block is mis-wired in the
`.mod`, or drop the 12 parameters and say so.

### C4. Real exchange rate persistence

autocorr(Q) = 0.19 vs 0.72 in the data (3.6% of the objective). With ρ_pvstar
pinned at 0.695 externally and no other persistent real driver, this is now
unfixable by estimation.

### C5. Things to check before any figure ships

- **TB impulse responses.** `eps_pc` moves TB by **+108.06**, `eps_pvst` by
  **−74.55**. If those are percent deviations, a single shock moving the trade
  balance by 100+ percent is a units/scale problem — most likely TB_ss = 0.018
  as a denominator. Verify before the IRF figures go in the paper.
- **Gross output / GDP prints twice with different values** in the same run:
  3.291 in the steady-state accounting block, 4.4885 in the summary. Nominal vs
  unweighted real sum. Pick one; the paper will quote whichever lands in the
  table.
- **"Simulation: no data (0 periods)."** Moments come from the theoretical
  spectrum, so this is GMM on second moments, not SMM. HP-filtering a theoretical
  spectrum and HP-filtering a 72-quarter sample are not the same object; the
  small-sample bias is uncorrected. Either simulate at T = 72 with the same
  filter, or state the approximation.

---

## Roadmap after the 2026-08-19 session

Objective composition at the last full run (C: cl=19.1, reallocation shock on,
purchaser-price consumption weights — superseded by the basic-price rebuild):

| block | value | % |
|---|---:|---:|
| corr(Y_i,PH_i) | 5.434 | **36.4** |
| rank correlations | 4.138 | **27.7** |
| std(PH_i) | 1.475 | 9.9 |
| labour comovement | 1.276 | 8.6 |
| std(Y_i) | 1.167 | 7.8 |
| aggregates | 0.850 | 5.7 |
| std(L_i) | 0.569 | 3.8 |
| goods expenditure share | 0.008 | 0.1 |

Single largest terms: `corr(Y₁₁,PH₁₁)` 22.4%, `rank corr labour` 14.7%,
`corr(Y₈,PH₈)` 8.8%, `rank corr output` 6.7%, `rank corr prices` 6.3%,
`corr(N,GDP)` 6.0%, `std(Y₅)` construction 5.5%.

**0. Re-run and re-baseline.** Everything below is conditioned on a composition
that may have shifted: today changed the α convention, the IO matrix, the
consumption weights, ombar, cl (was inert), added a shock and added a moment.

**1. Re-estimate — nothing has ever been estimated on this model.** All seven
free parameters still sit at their CMA-ES starting values (A1); the July run
managed 100 evaluations in 20 seconds. `cl` was inert until today, so it has
never been searched at all. This is mechanical and almost certainly the largest
single reduction available. Two prerequisites: decide `cl` first (it is not
interior-identified — the objective was still falling at 40 against a bound of
50, so a free search returns the bound), and fix the Jacobian so the run
produces standard errors (A3 — every one is currently 0.0).

**2. Sectoral variance decomposition.** The aggregate one (N, GDP) was decisive:
it showed `epsA_11` is 26% of employment variance and 3.7% of GDP, and `epsA_2`
the mirror image. The same calculation per sector would say *which shock* drives
`corr(Y₁₁,PH₁₁)` — 22% of the objective and currently undiagnosed. Cheap:
the state-space matrices are already on disk.

**3. `isigma_tfp_11 = 0.0737`.** The largest measured sectoral shock, in the
sector with the *lowest* data employment volatility (0.0197). Data say personal
services has volatile output and stable employment; the model cannot do both, so
a shock sized off output volatility blows up employment (model 0.095 vs 0.020).
Implicated in the labour rank correlation and `corr(N,GDP)` — together ~21% of
the objective. Check `compute_sectoral_shocks.jl` for sector 11.

**4. Investment.** Construction has *no final demand*: no household consumption
(its output is investment) and no investment in the model. VA 0.84% vs 6.45% in
data, gross output share 0.9% vs 7.8%, `VAE/VA` = 8.5x (a visible outlier), and
`std(Y₅)` = 0.168 vs 0.055 under the reallocation shock — 5.5% of the objective
from one sector the model structurally cannot represent. Also fixes C/GDP =
97.7% against 62% in data. Major `.mod` work.

**5. Government consumption.** Public administration 0.29% vs 5.08%. G is ~15%
of Chilean GDP and entirely absent. Cheaper than investment.

**6. Fixed factor for capital/resource rents.** Mining is 21.5% of model
employment against 1.6% of Chilean headcount, because the FGI convention books
copper rent (58% of mining gross output) as labour income. Also restores
labour/VA to 0.416 from 1.0. More a credibility and validation fix than an
objective reduction, but mining is 15.8% of value added and the number is
indefensible in a seminar. See B1.

**7. The weighting matrix is deciding what "fit" means.** 64% of the objective
sits in two cross-sectional pattern-matching blocks (`corr(Y_i,PH_i)` 36% + rank
correlations 28%), and the rank targets are 1.0 — unreachable by construction.
That is what pushes `cl` to its bound and what made flexible wages fit better
than the Chilean literature's sticky ones. Before optimising harder against this
criterion, it is worth asking whether it scores the right thing. FGI's precedent
is to target the model-data correlation of the cross-section with an identity
weight matrix, which is the same idea but a different relative weight.

**8. Smaller items.** Drop or re-specify the 12 `om_i` taste shifters (0.5% of N
variance, 0.1% of GDP, two exactly zero by construction). `autocorr(Q)` 0.18 vs
0.72 is 3.8% and unfixable by estimation while ρ_pvstar is pinned. Delete the
dead `kappa` column from `sector_calibration.csv` before it misleads someone —
it nearly misled me.

---

## Original order of work (superseded by the roadmap above)

1. ~~**A2** — close the 21.47 / 15.09 gap.~~ **Done.** Diagnosed as staleness, not
   a code-path bug; provenance guard added. Verification run still pending.
2. ~~**B1** — capital income.~~ **Done.** FGI convention adopted;
   `sector_calibration.csv` regenerated via `Data/build_sector_calibration.py`.
   Requires a `main_SOE_gap.jl` re-run and then re-estimation.
3. ~~**B3** — swap in data VA shares.~~ **Done**, along with the domestic IO
   matrix (B1b), the gross-output target, and the provenance dependency list.
   See B1c.
4. **C3** — determine whether the preference block is broken or should be dropped.
   Twelve parameters is a lot of the θ vector doing nothing.
5. **A1** — re-run CMA-ES properly with a real evaluation budget.
6. **A3** — fix the Jacobian so there are standard errors.
7. **C1** — the labor comovement problem is the substantive modeling question,
   and worth attacking only once 1–6 are settled.
