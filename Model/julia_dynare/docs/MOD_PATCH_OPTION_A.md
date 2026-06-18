# MOD File Patch — Option A: 12 Sectoral Demand Shocks

**File to edit:** `mod/NK_SOE_lev_gap2_smm.mod`  
**After editing, recompile with:** `julia --project=. main_SOE_gap.jl` (EXERCISE=0)

---

## 1. VAREXO block

**Remove** `eps_om` (aggregate preference shock).  
**Add** `eps_om_1` through `eps_om_12`.  
**New order must match the Julia Σe indices exactly:**

```
varexo
    eps_i       % 1  monetary policy
    epschi      % 2  labor supply (inactive — Σe[2,2]=0)
    eps_pvstar  % 3  import price
    epsA_1  epsA_2  epsA_3  epsA_4  epsA_5  epsA_6
    epsA_7  epsA_8  epsA_9  epsA_10 epsA_11 epsA_12  % 4:15  sectoral TFP
    eps_xi      % 16 aggregate demand
    eps_om_1  eps_om_2  eps_om_3  eps_om_4  eps_om_5  eps_om_6
    eps_om_7  eps_om_8  eps_om_9  eps_om_10 eps_om_11 eps_om_12  % 17:28 sectoral demand
;
```

---

## 2. PARAMETERS block

**Remove:** `sigma_om`  
**Add:** `sigma_om_1 sigma_om_2 ... sigma_om_12`  
**Keep:** `rho_om1` (common AR(1) persistence, shared across all sectoral demand shocks)

```dynare
parameters ... sigma_om_1 sigma_om_2 sigma_om_3 sigma_om_4 sigma_om_5 sigma_om_6
               sigma_om_7 sigma_om_8 sigma_om_9 sigma_om_10 sigma_om_11 sigma_om_12 ...;
```

---

## 3. PARAMETER CALIBRATION block

**Remove:** `sigma_om = ...;`  
**Add initial values** (will be overwritten by SMM):

```dynare
sigma_om_1 = 0.03;  sigma_om_2 = 0.03;  sigma_om_3 = 0.03;
sigma_om_4 = 0.03;  sigma_om_5 = 0.03;  sigma_om_6 = 0.03;
sigma_om_7 = 0.03;  sigma_om_8 = 0.03;  sigma_om_9 = 0.03;
sigma_om_10 = 0.03; sigma_om_11 = 0.03; sigma_om_12 = 0.03;
```

---

## 4. TASTE SHOCK AR(1) EQUATIONS

**Remove** the single aggregate AR(1):
```dynare
% REMOVE THIS:
eps_om = rho_om1 * eps_om(-1) + sigma_om * e_om;
```

**Add** 12 sector-specific AR(1) processes:
```dynare
% ADD THESE (one per sector):
eps_om_1  = rho_om1 * eps_om_1(-1)  + sigma_om_1  * eps_om_1;
eps_om_2  = rho_om1 * eps_om_2(-1)  + sigma_om_2  * eps_om_2;
eps_om_3  = rho_om1 * eps_om_3(-1)  + sigma_om_3  * eps_om_3;
eps_om_4  = rho_om1 * eps_om_4(-1)  + sigma_om_4  * eps_om_4;
eps_om_5  = rho_om1 * eps_om_5(-1)  + sigma_om_5  * eps_om_5;
eps_om_6  = rho_om1 * eps_om_6(-1)  + sigma_om_6  * eps_om_6;
eps_om_7  = rho_om1 * eps_om_7(-1)  + sigma_om_7  * eps_om_7;
eps_om_8  = rho_om1 * eps_om_8(-1)  + sigma_om_8  * eps_om_8;
eps_om_9  = rho_om1 * eps_om_9(-1)  + sigma_om_9  * eps_om_9;
eps_om_10 = rho_om1 * eps_om_10(-1) + sigma_om_10 * eps_om_10;
eps_om_11 = rho_om1 * eps_om_11(-1) + sigma_om_11 * eps_om_11;
eps_om_12 = rho_om1 * eps_om_12(-1) + sigma_om_12 * eps_om_12;
```

> **Note on naming:** In Dynare, the `varexo` name and the shock variable in the equation are the same symbol. If your current `.mod` uses `e_om` as the innovation and `eps_om` as the state, replicate that pattern: each `eps_om_i` is the state variable driven by the `varexo` `eps_om_i` directly (unit-variance innovation, amplitude scaled by `sigma_om_i`).

---

## 5. SECTORAL DEMAND EQUATIONS

Each sector's demand (taste/preference) shifter currently uses the single `eps_om`. Replace with sector-specific `eps_om_i`.

Find all occurrences of `exp(eps_om)` (or `eps_om` in linearized form) inside sectoral demand or utility equations and replace:

```dynare
% BEFORE (sector i demand):
... * exp(eps_om) ...

% AFTER (sector i demand, where i=1..12):
... * exp(eps_om_i) ...
```

The exact location depends on how taste shocks enter the CES demand system. In the NK-IOSOE model, `eps_om` typically appears in the sectoral expenditure share or the marginal utility of consumption for each good. Replace the shared `eps_om` with the sector-indexed `eps_om_$(i)`.

---

## 6. STEADY STATE / INITVAL

Remove any `eps_om = 0;` initialization.  
Add `eps_om_1 = 0; ... eps_om_12 = 0;` (or they will default to 0).

---

## 7. POST-EDIT CHECKLIST

After editing the `.mod` file:

1. Run `julia --project=. main_SOE_gap.jl` (EXERCISE=0) to recompile and check Blanchard-Kahn conditions.
2. Verify `n_exo = 28` in the Dynare output (was 17 before).
3. Verify the varexo ordering matches the Σe indices in `smm_estimation.jl` and `smm_model_moments.jl`:
   - Index 1 = `eps_i`, Index 3 = `eps_pvstar`, Index 4:15 = `epsA_1:12`, Index 16 = `eps_xi`, Index 17:28 = `eps_om_1:12`
4. Run `julia --project=. compute_data_moments.jl` to regenerate `sectoral_moments.csv` with the new `corr_YPH` column.
5. Run `julia --threads=auto --project=. run_smm_estimation.jl` to start estimation (35 params, 58 moments).
