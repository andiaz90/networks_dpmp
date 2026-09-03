#!/usr/bin/env julia
"""
compare_kappa_micro.jl — SMM-estimated sectoral price rigidity against the
micro-calibrated (frequency-of-price-adjustment) vector.

THIS IS THE EXHIBIT THE EXERCISE EXISTS FOR. The model's kappa_i are currently
calibrated from Pasten, Schoenle & Weber (2020) US PPI frequencies rescaled by
a single constant to one observed Chilean anchor (Albagli, Grigoli, Luttini,
Quevedo & Rojas 2026, regular PPI manufacturing median 0.252/month). Eleven of
the twelve sectors therefore carry an AMERICAN cross-sectional pattern. This
script asks whether Chilean aggregate data agree with it.

Reads (whichever exist):
  estimation_results/smm_estimates_kappa_affine.csv       level + dispersion
  estimation_results/smm_estimates_kappa_free_sh*.csv     eleven free kappa_i
  estimation_results/smm_estimates_lk*.csv                the lambda_kappa profile
  mod/params_jl.mod                                       kappa^cal (baseline)
  ../../Data/fpa_vector_few_industries_chile.csv          theta_m^cal

Writes:
  ../../Data/computed/kappa_comparison.csv
  ../../Data/computed/kappa_comparison.tex     (booktabs, paper-ready)

Reported in DURATIONS as well as kappa. A referee reads "prices reset every 3.4
months", not "kappa = 11.0", and the duration is the object the micro evidence
is measured in — which is the whole point of the comparison.

USAGE:  julia --project=. compare_kappa_micro.jl
"""

using CSV, DataFrames, Printf, Statistics

const JD = @__DIR__

# Honour SMM_ESTIMATION_DIR, the same redirect the estimation uses, so a scratch
# laptop run can be inspected without touching estimation_results/. Both folders
# are searched: the redirected one first, then the canonical one, so a scratch
# run that has only produced a checkpoint can still be compared against
# whatever completed estimates already exist.
const EST_CANON = joinpath(JD, "estimation_results")
const EST_DIRS  = let d = get(ENV, "SMM_ESTIMATION_DIR", "")
    dirs = String[]
    if !isempty(d)
        p = isabspath(d) ? d : joinpath(JD, d)
        isdir(p) && push!(dirs, p)
    end
    isdir(EST_CANON) && push!(dirs, EST_CANON)
    unique(dirs)
end
const DATA = abspath(joinpath(JD, "..", "..", "Data"))
const OUT  = joinpath(DATA, "computed")
mkpath(OUT)

const NSEC   = 12
const NKPC   = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]   # mining has no Phillips curve
const NAMES  = ["Agriculture", "Mining", "Manufacturing", "Utilities", "Construction",
                "Trade", "Transport", "Financial", "Real estate", "Business serv.",
                "Personal serv.", "Public admin."]
const EPS    = 10.0
const BETA   = 0.986

"""
    kappa_to_duration(kappa) -> quarters

Invert the Calvo->Rotemberg map kappa = s(eps-1)/[(1-s)(1-s*beta)], s = 1 - theta_q,
by bisection, and return the implied mean price duration 1/theta_q in quarters.
Monotone in s, so bisection is safe and exact to machine precision in 200 steps.
"""
function kappa_to_duration(kappa::Real)
    (!isfinite(kappa) || kappa <= 0) && return NaN
    lo, hi = 1e-9, 1 - 1e-12
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        km  = mid * (EPS - 1) / ((1 - mid) * (1 - mid * BETA))
        km < kappa ? (lo = mid) : (hi = mid)
    end
    s = 0.5 * (lo + hi)
    return 1.0 / max(1 - s, 1e-12)
end

"Read kappa_i out of params_jl.mod — the calibrated vector actually in the model."
function kappa_calibrated()
    p = joinpath(JD, "mod", "params_jl.mod")
    isfile(p) || error("""
        $(p) not found. It is written by main_SOE_gap.jl; run that first
        (EXERCISE=0) so the calibrated kappa vector exists to compare against.
        """)
    k = fill(NaN, NSEC)
    for ln in eachline(p)
        m = match(r"^\s*kappa_(\d+)\s*=\s*([-0-9.eE+]+)\s*;", ln)
        m === nothing && continue
        i = parse(Int, m.captures[1])
        1 <= i <= NSEC && (k[i] = parse(Float64, m.captures[2]))
    end
    any(isnan, k[NKPC]) && error("params_jl.mod is missing kappa_i for some NKPC sector.")
    return k
end

"Apply the utils.jl map to a named estimate dictionary."
function kappa_from_estimates(est::Dict{String,Float64}, kcal::Vector{Float64})
    haskey(est, "log_kappa_bar") && haskey(est, "lambda_kappa") || return nothing
    lk  = [log(kcal[i]) for i in NKPC]
    mu  = mean(lk)
    out = copy(kcal)
    for (k, i) in enumerate(NKPC)
        out[i] = exp(est["log_kappa_bar"] + est["lambda_kappa"] * (lk[k] - mu) +
                     get(est, "dlog_kappa_$(i)", 0.0))
    end
    return out
end

function read_est(path)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    (hasproperty(df, :param) && hasproperty(df, :value)) || return nothing
    return Dict(String(r.param) => Float64(r.value) for r in eachrow(df)
                if !ismissing(r.value) && isfinite(Float64(r.value)))
end

# --------------------------------------------------------------------------- #
kcal = kappa_calibrated()
fpa  = let f = joinpath(DATA, "fpa_vector_few_industries_chile.csv")
    isfile(f) ? [parse(Float64, l) for l in eachline(f) if !isempty(strip(l))] : fill(NaN, NSEC)
end

@printf "\n%s\n  SECTORAL PRICE RIGIDITY: SMM ESTIMATE vs MICRO CALIBRATION\n%s\n\n" repeat("=",74) repeat("=",74)
@printf "  Calibrated vector: Pasten-Schoenle-Weber (2020) US PPI monthly frequencies,\n"
@printf "  all durations rescaled by one constant so that MANUFACTURING matches Albagli\n"
@printf "  et al. (2026) Chile, 0.252/month. Sector 2 (Mining) has no Phillips curve.\n\n"

# Collect every available estimate across the searched folders. A scratch run
# usually has only smm_checkpoint.csv (written on every improvement) and no
# smm_estimates.csv until it finishes — so the CHECKPOINT is included, and
# labelled, rather than the script reporting "nothing found" mid-run.
cands = Tuple{String,String}[]
for dir in EST_DIRS
    tag = dir == EST_CANON ? "" : " [" * basename(dir) * "]"
    push!(cands, ("affine (level+dispersion)" * tag, joinpath(dir, "smm_estimates_kappa_affine.csv")))
    for f in sort(filter(x -> occursin(r"^smm_estimates_kappa_free_sh", x), readdir(dir)))
        push!(cands, ("free, ridge=" * replace(f, r"^.*_sh" => "", ".csv" => "") * tag, joinpath(dir, f)))
    end
    for f in sort(filter(x -> occursin(r"^smm_estimates_lk", x), readdir(dir)))
        push!(cands, ("profile " * replace(f, r"^.*_lk" => "lambda_kappa=", ".csv" => "") * tag, joinpath(dir, f)))
    end
    push!(cands, ("smm_estimates.csv" * tag, joinpath(dir, "smm_estimates.csv")))
    push!(cands, ("smm_checkpoint.csv (IN PROGRESS — not a final estimate)" * tag,
                  joinpath(dir, "smm_checkpoint.csv")))
end
@printf "  searching: %s\n\n" join(map(basename, EST_DIRS), ", ")

cols   = DataFrame(sector = 1:NSEC, name = NAMES,
                   theta_m_cal = fpa, kappa_cal = kcal,
                   dur_q_cal = kappa_to_duration.(kcal))
labels = String[]

for (lab, path) in cands
    est = read_est(path)
    est === nothing && continue
    kv = kappa_from_estimates(est, kcal)
    kv === nothing && continue
    push!(labels, lab)
    cols[!, "kappa_" * lab] = kv
    cols[!, "dur_q_" * lab] = kappa_to_duration.(kv)
    @printf "  loaded: %-28s  log kappa_bar = %7.4f   lambda_kappa = %6.4f\n" (
        lab) est["log_kappa_bar"] est["lambda_kappa"]
end

if isempty(labels)
    @printf "\n  No estimate carrying theta[39:51] was found in: %s\n" join(EST_DIRS, ", ")
    @printf "\n  Nothing has been ESTIMATED yet — the ladder's rungs 0-2 only verify the\n"
    @printf "  wiring and score the CALIBRATED kappa; they write no theta[39:51].\n"
    @printf "  Produce one first:\n"
    @printf "      bash run_kappa_ladder.sh 7                  # local, bounded, scratch output\n"
    @printf "      sbatch cluster/run_kappa_smm.sh affine      # the headline run\n\n"
    exit(0)
end

# --------------------------------------------------------------------------- #
lab = labels[1]
kest = cols[!, "kappa_" * lab]
@printf "\n  %-16s %8s %9s %9s %9s %9s %8s\n" "sector" "theta_m" "kappa cal" "kappa est" "dur cal" "dur est" "ratio"
@printf "  %s\n" repeat("-", 76)
for i in 1:NSEC
    if i ∉ NKPC
        @printf "  %-16s %8.3f %9.1f %9s %9.2f %9s %8s   (no NKPC)\n" (
            NAMES[i]) fpa[i] kcal[i] "--" cols.dur_q_cal[i] "--" "--"
        continue
    end
    @printf "  %-16s %8.3f %9.1f %9.1f %9.2f %9.2f %8.2f\n" NAMES[i] fpa[i] kcal[i] kest[i] (
        cols.dur_q_cal[i]) cols[!, "dur_q_"*lab][i] (kest[i] / kcal[i])
end

lc = [log(kcal[i]) for i in NKPC]
le = [log(kest[i]) for i in NKPC]
sp = let ra = sortperm(sortperm(lc)), rb = sortperm(sortperm(le))
    cor(Float64.(ra), Float64.(rb))
end
@printf "\n  %s\n" repeat("-", 76)
@printf "  corr( log kappa^cal , log kappa^est )  = %+.3f   (level)\n" cor(lc, le)
@printf "  rank corr                              = %+.3f   (ordering — the object the\n" sp
@printf "                                                    paper's heterogeneity claim rests on)\n"
@printf "  mean log kappa:  calibrated %.3f  ->  estimated %.3f   (x%.2f on the level)\n" (
    mean(lc)) mean(le) exp(mean(le) - mean(lc))
@printf "  sd   log kappa:  calibrated %.3f  ->  estimated %.3f   (x%.2f on the dispersion)\n" (
    std(lc)) std(le) (std(le) / max(std(lc), 1e-12))
@printf "  mean duration :  calibrated %.2f q ->  estimated %.2f q\n" (
    mean(cols.dur_q_cal[NKPC])) mean(cols[!, "dur_q_"*lab][NKPC])

@printf "\n  HOW TO READ THIS.\n"
@printf "    lambda_kappa near 1 and a rank corr near 1: Chilean macro data endorse the\n"
@printf "      US cross-section, and the calibration stands as it is.\n"
@printf "    lambda_kappa near 0: the data want HOMOGENEOUS stickiness — a rejection of\n"
@printf "      the paper's central mechanism, and it has to be reported as one.\n"
@printf "    lambda_kappa > 1: Chilean rigidity is MORE dispersed than the US pattern.\n"
@printf "    A large level shift with the ordering intact says the Albagli anchor is off,\n"
@printf "      not the Pasten pattern — those are separate claims and this separates them.\n"
@printf "\n  NONE of this substitutes for measuring the frequencies in the Chilean IPP\n"
@printf "  microdata. It tests whether aggregate data are CONSISTENT with the borrowed\n"
@printf "  cross-section; it cannot tell you the true one.\n\n"

# --------------------------------------------------------------------------- #
CSV.write(joinpath(OUT, "kappa_comparison.csv"), cols)

open(joinpath(OUT, "kappa_comparison.tex"), "w") do io
    println(io, "% Generated by compare_kappa_micro.jl — do not edit by hand.")
    println(io, "\\begin{tabular}{lrrrr}")
    println(io, "\\toprule")
    println(io, " & \\multicolumn{2}{c}{Calibrated (micro)} & \\multicolumn{2}{c}{Estimated (SMM)} \\\\")
    println(io, "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}")
    println(io, "Sector & \$\\kappa_i\$ & Duration (q) & \$\\kappa_i\$ & Duration (q) \\\\")
    println(io, "\\midrule")
    for i in 1:NSEC
        if i ∉ NKPC
            @printf(io, "%s & %.1f & %.2f & \\multicolumn{2}{c}{---} \\\\\n",
                    NAMES[i], kcal[i], cols.dur_q_cal[i])
        else
            @printf(io, "%s & %.1f & %.2f & %.1f & %.2f \\\\\n",
                    NAMES[i], kcal[i], cols.dur_q_cal[i], kest[i], cols[!, "dur_q_"*lab][i])
        end
    end
    println(io, "\\midrule")
    @printf(io, "Mean & %.1f & %.2f & %.1f & %.2f \\\\\n",
            mean(kcal[NKPC]), mean(cols.dur_q_cal[NKPC]),
            mean(kest[NKPC]), mean(cols[!, "dur_q_"*lab][NKPC]))
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

@printf "  wrote %s\n" joinpath(OUT, "kappa_comparison.csv")
@printf "  wrote %s\n\n" joinpath(OUT, "kappa_comparison.tex")
@printf "  For the scatter (log kappa^cal on x, log kappa^est on y, 45-degree line),\n"
@printf "  build it with the ipom-policy-plots house style rather than a default theme.\n\n"
