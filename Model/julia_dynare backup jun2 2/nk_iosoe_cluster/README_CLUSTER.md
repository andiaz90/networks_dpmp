# NK-IOSOE — Cluster bundle (estimation + shock plots)

Lightweight, self-contained package to run the SMM estimation of the NK-IOSOE
Chile model on an HPC cluster and then generate the shock-plot figures **using
the estimated parametrization**. All heavy artifacts (the compiled Dynare
context `nk_iosoe_context.jls` ≈ 20 MB and the sparse Jacobians) are **built on
the cluster**, not shipped — so this bundle is only ~1.5 MB.

## What's inside

```
nk_iosoe_cluster/
├── Project.toml                 # full dependency set (estimation + Plots/StatsPlots)
├── run_pipeline.sh              # the run sequence (single source of truth)
├── cluster/
│   ├── setup_cluster.jl         # one-time package install + precompile
│   ├── run_full_pipeline.sh     # SLURM: moments→compile→estimate→shocks
│   ├── run_estimation.sh        # SLURM: moments→compile→estimate only
│   └── run_shocks.sh            # SLURM: shocks only (after estimation)
├── Model/julia_dynare/          # all Julia source + mod/ source (.mod only)
│   ├── compute_data_moments.jl  # builds data moments from raw data
│   ├── main_SOE_gap.jl          # compiles the Dynare model → context.jls
│   ├── run_smm_estimation.jl    # CMA-ES SMM → Data/smm_estimates.csv
│   ├── run_all_shocks.jl        # oil/agr/min/mfg shock analyses
│   ├── *_shock_analysis.jl, shock_plots_common.jl, smm_*.jl, utils.jl, ...
│   └── mod/                     # NK_SOE_lev_gap2.mod + definition/solution blocks
└── Data/                        # raw inputs needed to run (NOT generated)
    ├── fpa_vector_few_industries_chile.csv   # price-adjustment freqs (κ)
    ├── IO_2021_chile.csv, sector_calibration.csv
    ├── count_workers_by_sector.csv, deflactor_pib.csv, pib_sectorial_bc.csv
    ├── reer_chile_bis.xlsx, datos_CCNN_mayo2025.xlsx
    └── computed/export_shares_chile.csv
```

Generated **on the cluster** (not shipped): `Data/{sectoral,aggregate}_moments.csv`,
`mod/nk_iosoe_context.jls`, `mod/NK_SOE_lev_gap2/...` Jacobians, `Data/smm_estimates.csv`,
all figures and LaTeX tables.

## Quick start

```bash
# 1. upload + extract
scp nk_iosoe_cluster_lite_*.tar.gz user@cluster:~/
ssh user@cluster
tar -xzf nk_iosoe_cluster_lite_*.tar.gz
cd nk_iosoe_cluster

# 2. one-time package install (login node)
module load julia            # or your cluster's exact module name
julia --project=. cluster/setup_cluster.jl

# 3a. everything in one job (estimation + shocks):
sbatch cluster/run_full_pipeline.sh

# 3b. OR split it:
sbatch cluster/run_estimation.sh      # writes Data/smm_estimates.csv
sbatch cluster/run_shocks.sh          # after the estimation job finishes
```

Monitor: `squeue -u $USER` and `tail -f nkiosoe_full_<JOBID>.out`.

## How the estimated parametrization reaches the shocks

`run_smm_estimation.jl` writes the estimated parameter vector to
`Data/smm_estimates.csv`. Every shock script (`oil/agr/min/mfg_shock_analysis.jl`)
checks for that file on startup and, if present, loads it as the structural
parametrization before solving — so the IRFs are produced at the estimated θ
automatically. No manual step needed; just run estimation before shocks (the
full-pipeline script enforces this order).

## Outputs

| File | Contents |
|------|----------|
| `Data/smm_estimates.csv` | estimated parameter vector θ̂ |
| `Data/smm_results.csv` | full moment-fit table (data vs model) |
| `Data/smm_checkpoint.csv` | warm-start checkpoint (resume an interrupted run) |
| `results/Figures/<shock>/*.pdf` | shock IRF + decomposition figures |
| `results/tables/*.tex` | shock LaTeX tables |
| `Model/julia_dynare/figures/`, `.../tables/` | local copies of the same |

## Resuming an interrupted estimation

The best θ found so far is saved to `Data/smm_checkpoint.csv` and reloaded as a
warm start on the next run. Just resubmit the same job.

## Configuration

Edit the `#SBATCH` headers in `cluster/*.sh` for your site: `--partition`,
`--cpus-per-task` (= Julia threads for the estimation), `--mem`, `--time`, and
the `module load` Julia version. The estimation scales well to ~48 threads.

## Notes

- Threads, not Distributed: the estimation uses `Threads.@threads`, so SLURM uses
  `--cpus-per-task=N`, not `-p N`. BLAS threads are pinned to 1.
- A local Julia depot (`.julia_depot/`) is created inside the package to keep the
  install self-contained and reusable across jobs.
- Shock figures are collected under `results/` (the scripts' Dropbox/Overleaf
  default path is overridden via `OVERLEAF_ROOT` in `run_pipeline.sh`).
