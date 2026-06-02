# NK-IOSOE SMM Estimation — HPC Cluster Guide

## Overview

This folder contains everything needed to run the SMM estimation on an HPC cluster with SLURM. The estimation uses Julia's native threading (`Threads.@threads`) — **not** distributed workers — so the SLURM script uses `--threads=N`, not `-p N`.

---

## Quick Start (4 steps)

### Step 1 — Prepare locally (run on your Mac/laptop)

Make sure these two scripts have been run first:
```bash
# In julia_dynare/
julia --project=. main_SOE_gap.jl          # generates nk_iosoe_context.jls
julia --project=. compute_data_moments.jl  # generates sectoral/aggregate_moments.csv
```

Then bundle everything:
```bash
bash cluster/prepare_cluster_package.sh
# Creates: nk_iosoe_cluster_YYYYMMDD.tar.gz
```

### Step 2 — Upload to cluster

```bash
scp nk_iosoe_cluster_YYYYMMDD.tar.gz user@cluster.bcentral.cl:~/
```

### Step 3 — Install packages on cluster (once)

```bash
ssh user@cluster.bcentral.cl
tar -xzf nk_iosoe_cluster_YYYYMMDD.tar.gz
cd nk_iosoe_cluster
julia --project=cluster cluster/setup_cluster.jl
```

This downloads and precompiles all Julia packages (~5–15 minutes first time, instant afterwards).

### Step 4 — Submit job

```bash
sbatch cluster/run_smm_hpc.sh
```

Monitor with:
```bash
squeue -u $USER
tail -f smm_JOBID.out
```

---

## Configuration

Edit `run_smm_hpc.sh` before submitting:

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--cpus-per-task` | 32 | Number of Julia threads. More = faster up to ~48. |
| `--mem` | 128G | Memory. 64G is sufficient; 128G gives headroom. |
| `--time` | 7 days | SMM with 30,000 evals takes 20–40 hours on 32 cores. |
| `--partition` | all | Change to your cluster's partition name. |
| `module load` | Julia/1.11.6 | Match your cluster's available Julia version. |

---

## File Structure

```
nk_iosoe_cluster/
├── julia_code/
│   ├── run_smm_estimation.jl      ← entry point
│   ├── smm_estimation.jl          ← CMA-ES + moment function
│   ├── smm_model_moments.jl       ← Klein solver + HP filter
│   ├── utils.jl                   ← hp_filtered_cov_fast, Lyapunov
│   ├── steady_ntwsoe.jl           ← steady-state outer system
│   ├── steady_ntwsoe_system.jl    ← steady-state inner system
│   └── mod/
│       ├── nk_iosoe_context.jls   ← compiled Dynare model (15 MB)
│       ├── dynare_endo_names.csv  ← endogenous variable names
│       ├── params_jl.mod          ← calibrated parameter values
│       └── NK_SOE_lev_gap2/model/julia/  ← compiled Jacobians
│           ├── SparseDynamicG1!.jl
│           ├── SparseDynamicResid!.jl
│           └── ...
├── Data/
│   ├── sectoral_moments.csv       ← HP-filtered data moments by sector
│   └── aggregate_moments.csv      ← Aggregate data moments
├── cluster/
│   ├── run_smm_hpc.sh             ← SLURM submission script
│   ├── setup_cluster.jl           ← one-time package install
│   ├── Project.toml               ← Julia package manifest
│   └── prepare_cluster_package.sh ← local bundling script
```

---

## Output Files

All outputs are written to `Data/`:

| File | Contents |
|------|----------|
| `smm_results.csv` | Full moment-fit table (data vs model, diff, rel%) |
| `smm_estimates.csv` | 23 estimated parameter values |
| `smm_checkpoint.csv` | Warm-start checkpoint (re-use in next run) |
| `smm_results_JOBID.csv` | Job-stamped copy of results |
| `smm_JOBID.out` | Full SLURM stdout log |

---

## Resuming an Interrupted Run

If the job times out or is cancelled, the best solution found so far is saved to `Data/smm_checkpoint.csv`. To resume:

1. Copy `smm_checkpoint.csv` back to your local `Data/` folder
2. Re-run `prepare_cluster_package.sh` (it includes the checkpoint automatically)
3. Re-upload and resubmit

The warm-start loads the checkpoint at the beginning of estimation.

---

## Performance Notes

- **32 threads** on a 3 GHz Xeon = ~55–70 ms/eval with HP filter
- **30,000 evaluations** (current `max_evals`) takes approximately **30–45 hours**
- Klein cache is 0% during active search (expected — every generation has new structural params)
- If the loss stagnates early, the fix is already applied (`ftol=1e-10`, parameter rescaling)
- `fail` count should be < 5% of evaluations; > 20% means bounds are too wide

## Dependencies

| Package | Purpose |
|---------|---------|
| `Dynare` | Klein (2000) first-order perturbation, model compilation |
| `NLsolve` | Steady-state solver (trust-region method) |
| `CMAEvolutionStrategy` | CMA-ES optimizer |
| `CSV`, `DataFrames` | Data I/O |
| `StatsBase` | Spearman rank correlation |
| `XLSX` | (optional) reading raw Excel data if re-computing moments |
