# Unified model + SMM fixes + inference — change log & validation

**Date:** 2026-06-18 · **Branch:** `agustin_development`

This implements three things you asked for: (1) a single canonical model file,
(2) the C1/C2 correctness fixes from the audit, and (3) GMM inference (standard
errors + overidentification J-test).

---

## 1. One model file

`NK_SOE_lev_gap2.mod` is now the **single** model. `NK_SOE_lev_gap2_smm.mod` and
the stale `nk_iosoe_smm_context.jls` were deleted. The unified model =
the oil/lev model **plus** the 12 sectoral demand (taste) shocks:

- **Added** endogenous: `om_1..om_12`, `norm_g`, `norm_s` (14 vars).
- **Added** equations: 12 sectoral AR(1) taste processes `om_i = rho_om1·om_i(-1) + sigma_om_i·eps_om_i`, plus the two normalizers `norm_g`, `norm_s` (14 eqs).
- **Added** shocks: `eps_om_1..12`; **removed** the aggregate `eps_om`.
- **Kept** the full oil block (`POstar`, `VOil`/`VNon`, `PIV_i`, `eps_postar`) — toggleable via `shock_eps_postar`, off during estimation.
- Sectoral/total consumption now carry the taste factor `gammag_i·exp(om_i)/norm_g` (and the `_s`, `_f` variants), with `norm_g/norm_s` preserving the goods/services budgets.
- **Removed** parameters: `sigma_om` (scalar), `rho_om2`, `shock_eps_om`. `utils.jl::write_params_mod` updated to match (no longer writes them).

**Key safety property:** at the steady state `om_i = 0` and `norm_g = norm_s = 1`, so every modified equation reduces to its original form — **the validated steady state is unchanged by construction.**

Static validation performed (no Dynare needed):

- Preprocessor macros balanced (`@#for`/`@#endfor` = 82/82).
- No dangling references to removed symbols (`eps_om`, `sigma_om`, `rho_om2`, `shock_eps_om`).
- All 29 declared shocks present in the `shocks;` block.
- **Square system:** 556→570 raw endogenous, +14 vars / +14 eqs exactly (the old raw count 556 matches Dynare's `orig_endo_nbr`, confirming the counter).
- **Parameter consistency:** every parameter the `.mod` declares is written by `write_params_mod`, and vice-versa (zero mismatches).

## 2. C1/C2 correctness fixes (SMM now matches the model)

- **C1 (shock layout):** `Σe` is now built **by shock name** (`active_shock_indices` in `smm_model_moments.jl`, resolved once in `build_baseline`), not by hardcoded indices. Active set: `eps_i`, `eps_pvstar`, `eps_xi`, `epsA_1..12`, `eps_om_1..12` (27 shocks); oil and labour-supply off. A startup check aborts if the loaded context doesn't yield exactly 27 active shocks (i.e., the old model is still loaded).
- **C2 (dead demand channel):** the 12 estimated `sigma_om_i` now exist in the model, so they actually drive sectoral demand; previously `set_param!` was a silent no-op. Also fixed: the common TFP persistence `rho_A` now sets `rho_tfp1_i` (the model uses sector-specific `rho_tfp1_i`), not the ignored scalar `rho_tfp1`.
- `run_smm_estimation.jl` loads only the single unified context.

## 3. Inference (`smm_inference.jl`)

GMM sandwich standard errors and the overidentification J-test, run automatically
at the end of estimation; writes `Data/smm_inference.csv` (param, estimate, se, t).

- Moment Jacobian `G = ∂model_moments/∂θ` by central finite differences at θ̂.
- `V = (G'WG)^{-1} G'W S W G (G'WG)^{-1}`, `SE = sqrt(diag(V)/T)`, `T=72`.
- `J = T·g'S^{-1}g ~ χ²(K−P)`, df = 58−35 = 23. (χ² p-value via a self-contained incomplete-gamma; validated against known critical values.)
- **`S` (data-moment covariance):** uses `Data/moment_cov.csv` if present (rigorous route — produce it via a moving-block bootstrap of the data series). If absent, falls back to principled delta-method diagonal variances (σ̂: σ²/2; ρ̂: (1−ρ²)²), clearly flagged as **indicative**. For the paper, supply the bootstrap `moment_cov.csv`.

---

## What you must do to use this (I cannot run Dynare here)

1. **Recompile the unified model** (this rebuilds the context with 29 shocks / 570 vars and the compiled Jacobian):
   ```
   julia --project=. main_SOE_gap.jl        # EXERCISE=0
   ```
2. **Run the smoke test** — it parses all sources and checks the setup:
   ```
   julia --project=. smoke_test.jl
   ```
3. **Confirm the solve** (from the audit recipe): steady-state residual < 1e-10 and all state-transition eigenvalues < 1.
4. **Estimate:**
   ```
   julia --threads=auto --project=. run_smm_estimation.jl
   ```
   It will print `Active shocks (by name): 27 of 29` — if not, it aborts telling you to recompile.
5. **(For publication SEs)** generate `Data/moment_cov.csv` via a block bootstrap of the data moments; otherwise SEs are reported as indicative.

## 4. Shock-analysis scripts (name-based shock columns)

All five IRF scripts (`oil_`, `agr_`, `min_`, `mfg_`, `agrmin_shock_analysis.jl`)
previously hardcoded the shock column from the old varexo order (`eps_postar`=5,
`epsA_k`=5+k). With the unified model `eps_postar`=4 and `epsA_k`=4+k, so these
were all off by one. They now resolve the `ghu` column **by name** via a shared
helper `exo_col(name, mod_dir)` in `shock_plots_common.jl`, which reads the
model's `modfile.json` exogenous list. This is robust to any future shock-order
change and self-adjusts to whichever model is compiled (it errors clearly if the
named shock is absent). All scripts parse cleanly.

## Things to watch / open items

- `run_dynare_smm_subprocess.jl` referenced the retired `_smm` model; it's now obsolete (the unified model is the only one). Remove or repoint if you use it.
- The hand-rolled Klein **fallback** still has stale hardcoded dims (audit C3/C4); it's dormant when Dynare.jl's native solver runs. De-hardcode it if you ever rely on the fallback.
