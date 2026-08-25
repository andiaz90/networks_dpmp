"""
compute_data_moments.jl
=======================
Compute all empirical moments needed for SMM estimation of the NK-IOSOE
Chile model. Saves results to modelo_chile/data_moments_chile.mat and also
returns the dm_chile NamedTuple so it can be used directly in Julia.

Translated from compute_data_moments.m (MATLAB).

DATA SOURCES  (all expected in Data/ folder at repo root)
  Employment (L) : count_workers_by_sector.csv   (monthly, 2006M1–)
  Prices     (P) : deflactor_pib.csv             (quarterly, 1996Q1–2023Q4)
  Output     (Y) : pib_sectorial_bc.xlsx         (quarterly real GDP, BCCh)
  REER       (Q) : reer_chile_bis.xlsx           (BIS Real Broad REER, monthly)
  TB/GDP         : datos_CCNN_mayo2025.xlsx       (BCCh CCNN, sheet "Nom trimestral")

METHODOLOGY
  - Monthly employment → quarterly average (mean of 3 months per quarter)
  - All series HP-filtered at quarterly frequency (λ = 1600)
  - Sample aligned to common window (default: 2006Q1 – 2023Q4)
  - std dev computed on HP-filtered log deviations

Usage:
  julia --project=. compute_data_moments.jl
or from REPL:
  include("compute_data_moments.jl")
"""

using CSV, DataFrames, XLSX
using Statistics, LinearAlgebra, SparseArrays
using Dates, Printf

# =========================================================================== #
#  PATHS                                                                       #
# =========================================================================== #

SCRIPT_DIR  = @__DIR__
REPO_ROOT   = abspath(joinpath(SCRIPT_DIR, "..", ".."))
DATA_DIR    = joinpath(REPO_ROOT, "Data")

# Output: plain CSV files in Data/ folder (no .mat, no MATLAB dependency)
OUT_SECTORAL  = joinpath(DATA_DIR, "sectoral_moments.csv")
OUT_AGGREGATE = joinpath(DATA_DIR, "aggregate_moments.csv")

@printf "\n%s\n" repeat("=", 61)
@printf "  Computing data moments for SMM estimation (Chile)\n"
@printf "  Output: %s\n          %s\n" OUT_SECTORAL OUT_AGGREGATE
@printf "%s\n\n" repeat("=", 61)

function _main()   # wrapped in function so we can use `return` for early exit

# =========================================================================== #
#  FAST PATH: if CSV outputs already exist, load them and skip all xlsx work  #
# =========================================================================== #
# FORCE_RECOMPUTE=1 skips the fast path (2026-08-24). Without it the only way to
# recompute was to delete the CSVs — and deleting the WRONG one destroys moments
# 61-63, which §2b then cannot carry forward. That trap fired on 2026-08-21 and
# again today: an edit to §9c produced no effect at all because the fast path
# silently loaded the previous run's file. Delete nothing; set the flag.
# (plain local, NOT const — this whole file lives inside _main())
FORCE_RECOMPUTE = get(ENV, "FORCE_RECOMPUTE", "0") != "0"
if FORCE_RECOMPUTE
    @printf "\n%s\n  FORCE_RECOMPUTE=1 — ignoring existing CSVs, recomputing from raw sources.\n%s\n\n" repeat("=",61) repeat("=",61)
end

if !FORCE_RECOMPUTE && isfile(OUT_SECTORAL) && isfile(OUT_AGGREGATE)
    sec = CSV.read(OUT_SECTORAL,  DataFrame)
    agg = CSV.read(OUT_AGGREGATE, DataFrame)

    # STALENESS CHECK (2026-07-08): files written before the corr(N,GDP) /
    # corr(N,GDP/N) additions also predate the sector-8/10 valid-window fix
    # (their std_Y was biased −13% by constant backfill). Never load them —
    # fall through and recompute from the raw sources.
    _agg_keys = Set(String.(agg.moment))
    if !("corr_NGDP" in _agg_keys && "corr_NAPL" in _agg_keys)
        @printf "\n%s\n  Existing CSV moment files are STALE (pre-2026-07-08:\n" repeat("=",61)
        @printf "  no corr_NGDP/corr_NAPL; sector-8/10 stds biased by backfill).\n"
        @printf "  Recomputing from raw sources and overwriting...\n%s\n" repeat("=",61)
        @goto recompute
    end

    # STALENESS CHECK (2026-08-21): files written before section 9c have no
    # idiosyncratic columns, so compute_sectoral_shocks.jl would silently fall
    # back to TOTAL volatility — the very double-counting 9c exists to remove.
    if !(hasproperty(sec, :std_Y_idio) && hasproperty(sec, :corr_YPH_idio))
        @printf "\n%s\n  Existing CSV moment files are STALE (pre-2026-08-21:\n" repeat("=",61)
        @printf "  no std_Y_idio/corr_YPH_idio common-factor columns).\n"
        @printf "  Recomputing from raw sources and overwriting...\n%s\n" repeat("=",61)
        @goto recompute
    end
    # Same date: files written before the COVID exclusion and the common
    # sectoral window carry moments computed on a different sample. They are
    # not comparable with anything produced now, so never load them.
    if !("covid_excluded" in _agg_keys)
        @printf "\n%s\n  Existing CSV moment files are STALE (pre-2026-08-21: no\n" repeat("=",61)
        @printf "  covid_excluded flag, so they predate the COVID exclusion and the\n"
        @printf "  common sectoral window fix). Recomputing and overwriting...\n%s\n" repeat("=",61)
        @goto recompute
    end

    @printf "\n%s\n  CSV moment files already exist (current format) — loading directly.\n" repeat("=",61)
    @printf "  (Delete and re-run to recompute from raw Excel/CSV sources.)\n"
    @printf "%s\n\n" repeat("=",61)

    # 'name' column may not exist if CSV was written by bootstrap_csv.jl
    _names = hasproperty(sec, :name) ? sec.name :
             ["Agriculture","Mining","Manufacturing","Utilities","Construction",
              "Trade/Hotels","Transport/Comm","Finance","Real Estate",
              "Business Serv.","Personal Serv.","Public Admin."]

    @printf "  %-20s  %8s  %8s  %8s\n" "Sector" "std(Y)" "std(PH)" "std(L)"
    @printf "  %s\n" repeat("-", 50)
    for (i, r) in enumerate(eachrow(sec))
        @printf "  %-20s  %8.4f  %8.4f  %8.4f\n" _names[i] r.std_Y r.std_PH r.std_L
    end
    agg_d = Dict(String(r.moment) => Float64(r.value) for r in eachrow(agg))
    @printf "\n  std(GDP)=%.4f  std(pi)=%.4f  corr(GDP,pi)=%.4f\n" agg_d["std_GDP"] agg_d["std_pi"] agg_d["corr_GDPpi"]
    @printf "  std(Q)=%.4f   autocorr(Q)=%.4f  TB/GDP=%.4f\n\n" agg_d["std_Q"] agg_d["autocorr_Q"] agg_d["TBGDP"]
    @printf "  Loaded from: %s\n\n" OUT_SECTORAL
    return   # done — no raw file processing needed
end

@label recompute
@printf "\n%s\n  Computing data moments from raw files (Chile)\n%s\n\n" repeat("=",61) repeat("=",61)
@printf "  (Tip: if pib_sectorial_bc.xlsx fails to load, open it in Excel\n"
@printf "   and File → Save As → .xlsx to fix XLSX.jl compatibility.)\n\n"

# =========================================================================== #
#  SETTINGS                                                                    #
# =========================================================================== #

NSEC   = 12
LAMBDA = 1600.0     # HP filter smoothing parameter (quarterly)

# Estimation sample
SAMPLE_START = (year=2006, q=1)
SAMPLE_END   = (year=2023, q=4)

# Sector classification (1-based indices)
GOODS    = [1,2,3,4,5]      # sectors 1-5
SERVICES = [6,7,8,9,10,11,12]  # sectors 6-12

SECTOR_NAMES = [
    "Agriculture", "Mining", "Manufacturing", "Utilities", "Construction",
    "Trade/Hotels", "Transport/Comm", "Finance", "Real Estate",
    "Business Serv.", "Personal Serv.", "Public Admin.",
]


# =========================================================================== #
#  HELPER FUNCTIONS                                                            #
# =========================================================================== #

"""
    hp_cycle(y, λ) -> Vector{Float64}

Standard Hodrick-Prescott filter. Solves (I + λ D'D) trend = y where D is
the (T-2)×T second-difference matrix. Returns the cycle component y - trend.
"""
function hp_cycle(y::AbstractVector{<:Real}, λ::Real)
    T  = length(y)
    T < 4 && return zeros(T)

    # Build the (T-2)×T second-difference matrix D in sparse form
    m   = T - 2
    ri  = vcat(1:m, 1:m, 1:m)
    ci  = vcat(1:m, 2:m+1, 3:m+2)
    val = vcat(fill(1.0, m), fill(-2.0, m), fill(1.0, m))
    D   = sparse(ri, ci, val, m, T)

    # Solve the HP system: (I + λ D'D) trend = y
    H     = sparse(1.0I, T, T) + λ * (D'D)
    trend = H \ collect(Float64, y)
    return collect(Float64, y) .- trend
end

"""
    fillmissing_linear(x) -> Vector{Float64}

Replace NaN entries by piecewise-linear interpolation (forward/backward fill
at boundaries). Equivalent to MATLAB's fillmissing(x, 'linear').
"""
function fillmissing_linear(x::AbstractVector{<:Real})
    result = collect(Float64, x)
    good   = findall(!isnan, result)
    isempty(good) && return result

    # Fill left edge
    result[1:good[1]-1]    .= result[good[1]]
    # Fill right edge
    result[good[end]+1:end] .= result[good[end]]
    # Linear interpolation between valid points
    for k in 1:length(good)-1
        i1, i2 = good[k], good[k+1]
        if i2 - i1 > 1
            slope = (result[i2] - result[i1]) / (i2 - i1)
            for j in i1+1:i2-1
                result[j] = result[i1] + slope * (j - i1)
            end
        end
    end
    return result
end

"""
    nanmean(x) -> Float64

Mean of a vector, ignoring NaN entries.
"""
nanmean(x) = mean(filter(!isnan, x))

"""
    monthly_to_quarterly(Y_m, yr_m, mth_m)
      -> (Y_q::Matrix, yr_out::Vector, qt_out::Vector)

Convert a monthly matrix to quarterly by averaging the 3 months within each
quarter. Partial quarters (1–2 months) are averaged over available months.
"""
function monthly_to_quarterly(Y_m::Matrix{<:Real},
                               yr_m::Vector{<:Integer},
                               mth_m::Vector{<:Integer})
    qt_m   = ceil.(Int, mth_m ./ 3)
    yq_idx = collect(zip(yr_m, qt_m))
    uq     = unique(yq_idx)
    sort!(uq)

    nQ  = length(uq)
    nS  = size(Y_m, 2)
    Y_q = fill(NaN, nQ, nS)

    for (t, (y, q)) in enumerate(uq)
        idx = findall(==(( y, q)), yq_idx)
        if length(idx) >= 1
            Y_q[t, :] = vec(mean(Y_m[idx, :], dims=1))
        end
    end

    yr_out = [yq[1] for yq in uq]
    qt_out = [yq[2] for yq in uq]
    return Y_q, yr_out, qt_out
end

# =========================================================================== #
#  COVID EXCLUSION  (2026-08-21)                                              #
# =========================================================================== #
#
# WHY. Under the Ley de Protección al Empleo (21.227, April 2020) workers whose
# contracts were SUSPENDED stayed on the AFC (Seguro de Cesantía) register and
# drew unemployment-insurance funds. They are legally employed and counted, but
# they were not working. So AFC headcount barely fell in exactly the sectors
# whose output collapsed:
#
#     sector             2020Q2 vs 2019Q4, log pp
#                          output      AFC employment
#     Personal services    -20.6            -4.5
#     Trade / hotels       -24.6           -10.5
#     Transport / comms    -18.1            -5.7
#
# Feeding that into an SMM objective asks the model to reproduce a comovement
# that is an artefact of a furlough scheme, not of labour demand. No setting of
# c or kappa_w can deliver it, and trying corrupts both. The pandemic quarters
# also dominate the aggregate second moments through one enormous common shock.
#
# WHAT WE DO. Filter the FULL series — the HP trend is estimated on every
# observation, so the cycle is well defined through the pandemic — then drop
# 2020Q1-2021Q2 when forming the moments. This is the standard treatment and is
# what Data/diagnose_sectoral_employment.py already did for its diagnostics.
#
# Set COVID_EXCLUDE=0 to keep the pandemic in, for the robustness table.

COVID_EXCLUDE = get(ENV, "COVID_EXCLUDE", "1") != "0"
COVID_FIRST   = (2020, 1)     # inclusive
COVID_LAST    = (2021, 2)     # inclusive

_qidx(y, q) = 4y + q

"""
    covid_keep(n) -> BitVector

Mask of length `n` for a quarterly series that ENDS at SAMPLE_END, false on the
excluded pandemic quarters. Every series here is clipped at SAMPLE_END by
`align_sample`, which masks rather than pads, so the END is what they share —
count backwards from it rather than assuming a common start.
"""
function covid_keep(n::Integer)
    keep = trues(n)
    COVID_EXCLUDE || return keep
    qend = _qidx(SAMPLE_END.year, SAMPLE_END.q)
    lo, hi = _qidx(COVID_FIRST...), _qidx(COVID_LAST...)
    for t in 1:n
        qt = qend - (n - t)
        (lo <= qt <= hi) && (keep[t] = false)
    end
    return keep
end

"Std of a cyclical series with the pandemic quarters dropped."
mstd(x::AbstractVector{<:Real}) = std(x[covid_keep(length(x))])

"""
    mac1(x) -> Float64

First-order autocorrelation of a cyclical series, pandemic quarters dropped.

A (t, t+1) PAIR is kept only if BOTH quarters survive the mask. Deleting the
excluded quarters and then correlating x[1:end-1] with x[2:end] would splice
2019Q4 onto 2021Q3 and count that as a one-quarter transition, which biases the
estimate toward zero exactly where the series has its largest movements.
"""
function mac1(x::AbstractVector{<:Real})
    k    = covid_keep(length(x)) .& isfinite.(x)
    pair = k[1:end-1] .& k[2:end]
    sum(pair) < 3 && return NaN
    return cor(x[1:end-1][pair], x[2:end][pair])
end

"""
Correlation of two cyclical series with the pandemic quarters dropped. The two
must already be the same length and calendar-aligned (the callers truncate to a
common tail first).
"""
function mcor(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    length(a) == length(b) || error("mcor: lengths $(length(a)) and $(length(b)) differ — align the tails before calling.")
    k = covid_keep(length(a))
    return complete_cor(a[k], b[k])
end

"""
    align_sample(values, yrs, qts, y0, q0, y1, q1)

Return the subset of `values` (rows) whose (year, quarter) lies in [y0Q q0, y1Q q1].
"""
function align_sample(values, yrs, qts, y0, q0, y1, q1)
    mask = [((y, q) >= (y0, q0)) && ((y, q) <= (y1, q1))
            for (y, q) in zip(yrs, qts)]
    return values[mask, :], yrs[mask], qts[mask]
end

"""
    complete_cor(a, b) -> Float64

Pearson correlation using only rows where both a[i] and b[i] are finite.
Equivalent to MATLAB's corr(a, b, 'rows', 'complete').
"""
function complete_cor(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    ok = isfinite.(a) .& isfinite.(b)
    sum(ok) < 3 && return NaN
    return cor(a[ok], b[ok])
end


# =========================================================================== #
#  1. SECTORAL EMPLOYMENT                                                      #
#     count_workers_by_sector.csv  (monthly, columns in alphabetical order)   #
# =========================================================================== #

@printf "--- 1. Loading employment data ---\n"

# CSV columns 2–13 are sectors in ALPHABETICAL order.
# Mapping: CSV-alpha position → model sector number
ALPHA_TO_MODEL = [1, 10, 5, 4, 8, 3, 2, 11, 12, 9, 7, 6]

fname_emp = joinpath(DATA_DIR, "count_workers_by_sector_sa.csv")  # SA (raw: count_workers_by_sector.csv)
emp_df    = CSV.read(fname_emp, DataFrame)

# First column is date (YYYY-MM-DD); columns 2–13 are sector employment
date_col   = emp_df[:, 1]
yr_m_emp   = year.(date_col)
mth_m_emp  = month.(date_col)

L_alpha = Matrix{Float64}(emp_df[:, 2:13])
L_monthly = fill(NaN, size(L_alpha, 1), NSEC)
for k in 1:NSEC
    L_monthly[:, ALPHA_TO_MODEL[k]] = L_alpha[:, k]
end

L_qrt, yr_q_emp, qt_q_emp = monthly_to_quarterly(L_monthly, yr_m_emp, mth_m_emp)
@printf "  Employment: %d monthly → %d quarterly obs  (%dQ%d – %dQ%d)\n" size(L_monthly,1) size(L_qrt,1) yr_q_emp[1] qt_q_emp[1] yr_q_emp[end] qt_q_emp[end]


# =========================================================================== #
#  2. SECTORAL PRICE DEFLATORS                                                 #
#     deflactor_pib.csv  (quarterly, date format "MMMYYYY", 28 series)        #
# =========================================================================== #

@printf "\n--- 2. Loading price deflators ---\n"

fname_defl = joinpath(DATA_DIR, "deflactor_pib_sa.csv")  # SA (raw: deflactor_pib.csv)

# File has 3 metadata rows (SERIES, DESCRIPCION, UNIDAD), then data rows.
# We read everything as strings, skip the header rows.
raw_defl = CSV.read(fname_defl, DataFrame, header=false, types=String,
                    silencewarnings=true)
# Rows 1–3 are metadata; data starts at row 4
data_rows_defl = eachrow(Matrix{Union{String,Missing}}(raw_defl[4:end, :]))
n_defl = nrow(raw_defl) - 3

yr_defl = zeros(Int, n_defl)
q_defl  = zeros(Int, n_defl)

MONTH_TO_Q = Dict("MAR"=>1, "JUN"=>2, "SEP"=>3, "DIC"=>4)

for (t, row) in enumerate(data_rows_defl)
    s   = strip(replace(string(row[1]), '"' => ""))
    mon = uppercase(s[1:3])
    yr  = parse(Int, strip(s[4:end]))
    yr_defl[t] = yr
    q_defl[t]  = get(MONTH_TO_Q, mon, 0)
end

P_raw = fill(NaN, n_defl, 28)
for (t, row) in enumerate(data_rows_defl)
    for c in 1:28
        v = string(get(row, c+1, ""))
        v_stripped = strip(replace(v, '"' => ""))
        if v_stripped != "_" && v_stripped != "" && v_stripped != "missing"
            parsed = tryparse(Float64, v_stripped)
            parsed !== nothing && (P_raw[t, c] = parsed)
        end
    end
end

@printf "  Deflators: %d quarterly obs  (%dQ%d – %dQ%d)\n" n_defl yr_defl[1] q_defl[1] yr_defl[end] q_defl[end]

# Map 28 deflator columns → 12 model sectors (some sectors average sub-components)
DEFL_COLS = [
    [1, 2],   # Sector  1: Agro + Pesca (simple mean)
    [3],      # Sector  2: Mining (aggregate)
    [6],      # Sector  3: Manufacturing (aggregate)
    [15],     # Sector  4: Utilities
    [16],     # Sector  5: Construction
    [17, 18], # Sector  6: Retail + Hotels
    [19, 20], # Sector  7: Transport + Comm
    [21],     # Sector  8: Finance
    [23],     # Sector  9: Real Estate
    [22],     # Sector 10: Business services
    [24],     # Sector 11: Personal services
    [25],     # Sector 12: Public admin
]

P_sec = fill(NaN, n_defl, NSEC)
for i in 1:NSEC
    cols = DEFL_COLS[i]
    P_sec[:, i] = [nanmean(P_raw[t, cols]) for t in 1:n_defl]
end
P_agg = P_raw[:, 28]   # total GDP deflator (column 28)


# =========================================================================== #
#  3. SECTORAL REAL OUTPUT                                                     #
#     pib_sectorial_bc.xlsx  Sheet "Cuadro"                                   #
# =========================================================================== #

@printf "\n--- 3. Loading sectoral real output (pib_sectorial_bc.csv) ---\n"

# pib_sectorial_bc.csv was generated from the BCCh xlsx by Python/openpyxl
# (XLSX.jl cannot handle the non-standard XML in that file).
# CSV layout: row 1 = title, row 2 = empty, row 3 = column headers, rows 4+ = data.
# Col 1 = date (YYYY-MM-DD HH:MM:SS), cols 2–33 = activity values.
# Sector mapping (pib_num col k = CSV data col k, i.e. overall CSV col k+1):
#   pib_num 1  = Agropecuario-silvícola      → sector 1 (+ pib_num 2 Pesca)
#   pib_num 2  = Pesca
#   pib_num 3  = Minería                      → sector 2
#   pib_num 6  = Ind. Manufacturera           → sector 3
#   pib_num 16 = Electricidad/gas/agua        → sector 4
#   pib_num 17 = Construcción                 → sector 5
#   pib_num 18 = Comercio, rest. y hoteles    → sector 6
#   pib_num 21 = Transporte                   → sector 7 (+ pib_num 22 Comunicaciones)
#   pib_num 22 = Comunicaciones
#   pib_num 24 = Serv. financieros            → sector 8
#   pib_num 25 = Serv. empresariales          → sector 10
#   pib_num 26 = Vivienda e inmobiliarios     → sector 9
#   pib_num 27 = Serv. personales             → sector 11
#   pib_num 28 = Admón. pública               → sector 12
#   pib_num 31 = PIB total (col 32)           → GDP

fname_pib_csv = joinpath(DATA_DIR, "pib_sectorial_bc_sa.csv")  # SA (raw: pib_sectorial_bc.csv)
Y_sec_raw     = Matrix{Float64}(undef, 0, NSEC)
GDP_data      = Float64[]
yr_y          = Int[]
qt_y          = Int[]
pib_error_msg = ""

if isfile(fname_pib_csv)
    try
        # Read raw CSV — all as String to handle mixed date/numeric columns
        raw_df   = CSV.read(fname_pib_csv, DataFrame; header=false, types=String,
                            silencewarnings=true)
        # Skip 3 header rows (title, blank, column labels); data starts at row 4
        n_pib    = nrow(raw_df) - 3
        n_cols   = ncol(raw_df)
        @printf "  CSV rows (data): %d,  columns: %d\n" n_pib n_cols

        # --- Parse dates from col 1 (format "YYYY-MM-DD HH:MM:SS") ---------- #
        yr_y_pib = zeros(Int, n_pib)
        qt_y_pib = zeros(Int, n_pib)
        for t in 1:n_pib
            s = strip(string(raw_df[t+3, 1]))
            # Extract year and month from leading "YYYY-MM-" substring
            if length(s) >= 7
                yr_  = tryparse(Int, s[1:4])
                mon_ = tryparse(Int, s[6:7])
                if yr_ !== nothing && mon_ !== nothing
                    yr_y_pib[t] = yr_
                    qt_y_pib[t] = ceil(Int, mon_ / 3)
                end
            end
        end
        @printf "  Date range: %dQ%d – %dQ%d\n" yr_y_pib[1] qt_y_pib[1] yr_y_pib[end] qt_y_pib[end]

        # --- Parse numeric data (cols 2–33 of CSV = pib_num cols 1–32) ------ #
        n_data_cols = min(32, n_cols - 1)
        pib_num = fill(NaN, n_pib, n_data_cols)
        for t in 1:n_pib
            for c in 1:n_data_cols
                v = strip(string(raw_df[t+3, c+1]))
                isempty(v) && continue
                p = tryparse(Float64, v)
                p !== nothing && isfinite(p) && (pib_num[t, c] = p)
            end
        end
        n_ok = sum(!isnan, pib_num)
        @printf "  Numeric parse: %d / %d cells valid (%.0f%%)\n" n_ok length(pib_num) (100*n_ok/length(pib_num))
        n_ok == 0 && error("All numeric cells parsed as NaN — check CSV structure.")

        # --- Aggregate to 12 model sectors ---------------------------------- #
        Y_sec_raw = fill(NaN, n_pib, NSEC)
        nc = n_data_cols
        nc >= 2  && (Y_sec_raw[:, 1]  = pib_num[:, 1] .+ pib_num[:, 2])   # Agro+Pesca
        nc >= 3  && (Y_sec_raw[:, 2]  = pib_num[:, 3])                      # Mining
        nc >= 6  && (Y_sec_raw[:, 3]  = pib_num[:, 6])                      # Manufacturing
        nc >= 16 && (Y_sec_raw[:, 4]  = pib_num[:, 16])                     # Utilities
        nc >= 17 && (Y_sec_raw[:, 5]  = pib_num[:, 17])                     # Construction
        nc >= 18 && (Y_sec_raw[:, 6]  = pib_num[:, 18])                     # Trade+Hotels
        nc >= 22 && (Y_sec_raw[:, 7]  = pib_num[:, 21] .+ pib_num[:, 22])  # Transport+Comm
        nc >= 24 && (Y_sec_raw[:, 8]  = pib_num[:, 24])                     # Finance
        nc >= 26 && (Y_sec_raw[:, 9]  = pib_num[:, 26])                     # Real Estate
        nc >= 25 && (Y_sec_raw[:, 10] = pib_num[:, 25])                     # Business Serv
        nc >= 27 && (Y_sec_raw[:, 11] = pib_num[:, 27])                     # Personal Serv
        nc >= 28 && (Y_sec_raw[:, 12] = pib_num[:, 28])                     # Public Admin
        nc >= 31 && (GDP_data         = pib_num[:, 31])                      # PIB total

        yr_y = yr_y_pib
        qt_y = qt_y_pib

    catch e
        err_str = sprint(showerror, e)
        pib_error_msg = "  EXCEPTION: $err_str\n"
        @printf "\n  %s\n  EXCEPTION reading pib_sectorial_bc.csv:\n  %s\n  %s\n\n" repeat("!",58) err_str repeat("!",58)
    end
else
    pib_error_msg = "  File not found: $fname_pib_csv\n  → Regenerate with: python3 scripts/xlsx_to_csv.py\n"
    @printf "  File not found: %s\n" fname_pib_csv
    @printf "  Run: python3 scripts/xlsx_to_csv.py  to regenerate from pib_sectorial_bc.xlsx\n"
end


# =========================================================================== #
#  4. ALIGN ALL SERIES TO COMMON SAMPLE                                        #
# =========================================================================== #

@printf "\n--- 4. Aligning to common sample (%dQ%d – %dQ%d) ---\n" SAMPLE_START.year SAMPLE_START.q SAMPLE_END.year SAMPLE_END.q

L_sample, _, _ = align_sample(L_qrt,  yr_q_emp, qt_q_emp,
                                SAMPLE_START.year, SAMPLE_START.q,
                                SAMPLE_END.year,   SAMPLE_END.q)

P_sample, _, _ = align_sample(P_sec,  yr_defl, q_defl,
                                SAMPLE_START.year, SAMPLE_START.q,
                                SAMPLE_END.year,   SAMPLE_END.q)

Pa_sample, _, _ = align_sample(reshape(P_agg, :, 1), yr_defl, q_defl,
                                 SAMPLE_START.year, SAMPLE_START.q,
                                 SAMPLE_END.year,   SAMPLE_END.q)
Pa_sample = vec(Pa_sample)

Y_qrt = Matrix{Float64}(undef, 0, NSEC)
GDP_sample = Float64[]

if !isempty(yr_y)
    Y_tmp, yr_y_s, _ = align_sample(Y_sec_raw, yr_y, qt_y,
                                      SAMPLE_START.year, SAMPLE_START.q,
                                      SAMPLE_END.year,   SAMPLE_END.q)
    if !isempty(Y_tmp)
        Y_qrt      = Y_tmp
        GDP_sample = align_sample(reshape(GDP_data, :, 1), yr_y, qt_y,
                                   SAMPLE_START.year, SAMPLE_START.q,
                                   SAMPLE_END.year,   SAMPLE_END.q)[1] |> vec
        @printf "  Y_qrt aligned: %d quarters\n" size(Y_qrt, 1)
    end
end

nT = size(L_sample, 1)
@printf "  Common sample: %d quarters (L), %d quarters (P)\n" nT size(P_sample,1)


# --------------------------------------------------------------------------- #
#  4b. SECTORAL MOMENTS ON ONE COMMON WINDOW  (2026-08-21 — BUG FIX)           #
# --------------------------------------------------------------------------- #
#
# THE BUG. align_sample clips each series to SAMPLE_START..SAMPLE_END, but the
# series do not all START there. The seasonally-adjusted output vintage
# (pib_sectorial_bc_sa.csv) begins 2013Q1, while employment and the deflators
# begin 2006Q1. The result, in every sectoral_moments.csv written before today:
#
#     std_Y     2013Q1-2023Q4   (44 quarters)
#     corr_YPH  2013Q1-2023Q4   (44 quarters)
#     std_PH    2006Q1-2023Q4   (72 quarters)   <- different sample
#     std_L     2006Q1-2023Q4   (72 quarters)   <- different sample
#
# That is fatal for moments 44-46, the cross-sectional RANK correlations, which
# compare the ordering of std_L (or std_PH) against the ordering of std_Y. Two
# orderings measured over different decades are not comparable, and the model
# produces both from a single stationary distribution. It also makes the level
# blocks internally inconsistent.
#
# HOW MUCH IT MATTERED. The rank correlation between data std_L and data std_Y:
#     mismatched (72q vs 44q):  -0.231
#     aligned    (44q vs 44q):  +0.000
# The -0.231 was an artefact. It comes almost entirely from the sectors with the
# strongest AFC coverage phase-in over 2006-2012 — Mining (std_L 0.045 -> 0.017
# once the early years are dropped), Finance (0.046 -> 0.015), Public Admin
# (0.043 -> 0.028). The Seguro de Cesantia covers only contracts signed after
# October 2002, so its sectoral counts grow mechanically through the early
# sample at very different rates by sector (Utilities +142% log, Mining +116%,
# vs Agriculture +21%). HP(1600) removes a smooth trend but not that.
#
# FIX. Restrict every sectoral series to the window where ALL THREE are
# observed, before any filtering. Set SECTORAL_COMMON_WINDOW=0 to reproduce the
# old (mismatched) behaviour for comparison.

SECTORAL_COMMON_WINDOW = get(ENV, "SECTORAL_COMMON_WINDOW", "1") != "0"

if SECTORAL_COMMON_WINDOW && !isempty(Y_qrt)
    # CAREFUL: align_sample MASKS rows rather than padding them, so a series
    # that starts late simply comes back SHORTER. Row t is therefore NOT the
    # same calendar quarter across L_sample (72 rows), P_sample (72) and
    # Y_qrt (44). What they share is the END: every one is clipped at
    # SAMPLE_END. So take the last nc rows of each — the same convention the
    # aggregate block already uses via `end-nmin+1:end`.
    nc = minimum((size(L_sample,1), size(P_sample,1), size(Y_qrt,1)))

    # Within that tail window, drop any leading rows where some sector is still
    # unobserved (a few deflator series start later than the output vintage).
    good = trues(nc)
    for i in 1:NSEC, M in (L_sample, P_sample, Y_qrt)
        col = M[end-nc+1:end, i]
        good .&= .!(isnan.(col) .| (col .<= 0))
    end
    w1 = findfirst(good); w2 = findlast(good)

    if w1 === nothing || w2 === nothing || (w2 - w1 + 1) < 20
        @printf "  WARNING: common sectoral window has < 20 quarters — leaving series on their own samples.\n"
        @printf "  Rank moments 44–46 are then NOT comparable; treat them as untargeted.\n"
    else
        L_sample = L_sample[end-nc+1:end, :][w1:w2, :]
        P_sample = P_sample[end-nc+1:end, :][w1:w2, :]
        Y_qrt    = Y_qrt[   end-nc+1:end, :][w1:w2, :]
        # GDP_sample comes from the same align_sample call as Y_qrt, so it is
        # already the same length; keep it in step.
        if length(GDP_sample) >= nc
            GDP_sample = GDP_sample[end-nc+1:end][w1:w2]
        end
        nT = size(L_sample, 1)
        @printf "  Sectoral common window: %d quarters (tail-aligned at %dQ%d), all of Y/PH/L observed for all %d sectors.\n" nT SAMPLE_END.year SAMPLE_END.q NSEC
        @printf "  std_Y, std_PH, std_L and corr_YPH are now ALL computed on this window.\n"
    end
else
    @printf "  SECTORAL_COMMON_WINDOW=0 — sectoral series left on their own samples (rank moments 44–46 not comparable).\n"
end


# =========================================================================== #
#  5. HP-FILTER AND COMPUTE SECTORAL STD DEVS                                 #
# =========================================================================== #

@printf "\n--- 5. HP filtering and computing std devs ---\n"

if COVID_EXCLUDE
    _k = covid_keep(nT)
    @printf "  COVID EXCLUSION: %dQ%d–%dQ%d dropped from every moment (%d of %d sectoral quarters kept).\n" (
        COVID_FIRST[1]) (COVID_FIRST[2]) (COVID_LAST[1]) (COVID_LAST[2]) sum(_k) nT
    @printf "  Series are HP-filtered on the FULL sample; only the moment sums skip those quarters.\n"
else
    @printf "  COVID_EXCLUDE=0 — pandemic quarters retained.\n"
end

# Employment std devs
l_d = fill(NaN, NSEC)
for i in 1:NSEC
    x = log.(L_sample[:, i])
    any(isnan.(x) .| isinf.(x)) && (@printf "  WARNING: sector %d employment has missing values, skipping.\n" i; continue)
    l_d[i] = mstd(hp_cycle(x, LAMBDA))
end

# Price deflator std devs
p_d = fill(NaN, NSEC)
for i in 1:NSEC
    x   = log.(P_sample[:, i])
    bad = isnan.(x) .| isinf.(x)
    if sum(.!bad) < 20
        @printf "  WARNING: sector %d price has too few valid obs, skipping.\n" i
        continue
    end
    p_d[i] = mstd(hp_cycle(fillmissing_linear(x), LAMBDA))
end

@printf "  Employment std devs (%%): %s\n" join([@sprintf("%.3f", v*100) for v in l_d], "  ")
@printf "  Price      std devs (%%): %s\n" join([@sprintf("%.3f", v*100) for v in p_d], "  ")

# Sectoral output std devs, and their first-order autocorrelations.
#
# PERSISTENCE MOMENTS (2026-08-21). Before today the 63-moment vector contained
# exactly ONE autocorrelation, autocorr(Q), while the estimation had three free
# persistence parameters (rho_A, rho_om, rho_zeta). They were identified only
# indirectly, through the way the HP filter reweights variance across
# frequencies, which is why rho_A ran to its upper bound of 0.99: nothing in the
# objective pushed back against maximum persistence. autocorr(GDP), autocorr(pi)
# and these twelve sectoral autocorrelations discipline it directly.
#
# The sectoral spread is also informative in its own right. In the data it runs
# from ~0.11 (Agriculture) to ~0.89 (Finance), which a single COMMON rho_A
# cannot generate; targeting these twelve is what makes that restriction
# testable rather than assumed.
y_d  = fill(NaN, NSEC)
ac_y = fill(NaN, NSEC)
if !isempty(Y_qrt)
    for i in 1:NSEC
        x   = Y_qrt[:, i]
        bad = (x .<= 0) .| isnan.(x)
        if sum(.!bad) < 20
            @printf "  WARNING: sector %d output has too few valid obs (%d), skipping.\n" i sum(.!bad)
            continue
        end
        # Trim to the sector's valid window before filtering. Sectors 8 (Finance)
        # and 10 (Business services) only start in 2013Q1; backfilling their
        # leading NaNs with a constant (old behaviour) shrank their HP std by
        # ~13% and corrupted the cross-sector volatility ranking used by the
        # rank-correlation moments. fillmissing_linear now handles interior
        # gaps only.
        i1 = findfirst(!, bad); i2 = findlast(!, bad)
        x_log = log.(fillmissing_linear(x[i1:i2]))
        y_d[i]  = mstd(hp_cycle(x_log, LAMBDA))
        ac_y[i] = mac1(hp_cycle(x_log, LAMBDA))
        (i1 > 1 || i2 < length(x)) &&
            @printf "  NOTE: sector %d output std computed on valid window (obs %d–%d of %d).\n" i i1 i2 length(x)
    end
    @printf "  Output     std devs (%%): %s\n" join([@sprintf("%.3f", v*100) for v in y_d], "  ")
else
    @printf "  Output data not loaded; y_d set to NaN.\n"
end

# Aggregate inflation: QoQ growth of total GDP deflator
pi_series = Pa_sample[2:end] ./ Pa_sample[1:end-1]
pi_series[pi_series .<= 0] .= NaN
pi_hp    = hp_cycle(fillmissing_linear(log.(pi_series)), LAMBDA)
d_std_pi = mstd(pi_hp)
d_ac_pi  = mac1(pi_hp)
@printf "  std(pi) = %.5f\n" d_std_pi


# =========================================================================== #
#  6. AGGREGATE OUTPUT MOMENTS                                                 #
# =========================================================================== #

@printf "\n--- 6. Aggregate moments ---\n"

d_std_GDP    = NaN
d_ac_GDP     = NaN
d_corr_GDPpi = NaN
GDP_hp       = Float64[]

if !isempty(GDP_sample) && sum(GDP_sample .> 0) >= 20
    GDP_hp       = hp_cycle(fillmissing_linear(log.(GDP_sample)), LAMBDA)
    d_std_GDP    = mstd(GDP_hp)
    d_ac_GDP     = mac1(GDP_hp)
    nmin         = min(length(GDP_hp)-1, length(pi_hp))
    d_corr_GDPpi = mcor(GDP_hp[end-nmin+1:end], pi_hp[end-nmin+1:end])
    @printf "  std(GDP)      = %.5f  (BCCh pib_sectorial_bc)\n" d_std_GDP
    @printf "  corr(GDP, pi) = %.4f\n" d_corr_GDPpi
elseif !isempty(Y_qrt)
    GDP_real     = vec(sum(Y_qrt, dims=2))
    GDP_hp       = hp_cycle(log.(GDP_real), LAMBDA)
    d_std_GDP    = mstd(GDP_hp)
    d_ac_GDP     = mac1(GDP_hp)
    nmin         = min(length(GDP_hp)-1, length(pi_hp))
    d_corr_GDPpi = mcor(GDP_hp[end-nmin+1:end], pi_hp[end-nmin+1:end])
    @printf "  std(GDP)      = %.5f  (sum of sectoral VA)\n" d_std_GDP
    @printf "  corr(GDP, pi) = %.4f\n" d_corr_GDPpi
else
    d_std_GDP    = 0.021
    d_ac_GDP     = 0.75
    d_corr_GDPpi = -0.15
    @printf "  std(GDP) / corr: using defaults (no output data).\n"
end

# --- Employment comovement moments -------------------------------------- #
# corr(N, GDP): aggregate employment–output comovement (strongly positive in
# Chilean data, ≈ +0.70). corr(N, GDP/N): employment–labor-productivity
# correlation (≈ +0.07, the Galí-1999 fact). Together these discipline the
# demand vs supply shock mix and the wage stickiness kappaw.
d_corr_NGDP = NaN
d_corr_NAPL = NaN
if !isempty(GDP_sample) && size(L_sample, 1) >= 20 && sum(GDP_sample .> 0) >= 20
    N_series = vec(sum(L_sample, dims=2))
    nmin_n   = min(length(GDP_sample), length(N_series))   # both end at SAMPLE_END
    g_log    = log.(fillmissing_linear(GDP_sample[end-nmin_n+1:end]))
    n_log    = log.(N_series[end-nmin_n+1:end])
    n_hp     = hp_cycle(n_log, LAMBDA)
    g_hp     = hp_cycle(g_log, LAMBDA)
    apl_hp   = hp_cycle(g_log .- n_log, LAMBDA)
    d_corr_NGDP = mcor(n_hp, g_hp)
    d_corr_NAPL = mcor(n_hp, apl_hp)
    @printf "  corr(N, GDP)   = %.4f  (%d quarters)\n" d_corr_NGDP nmin_n
    @printf "  corr(N, GDP/N) = %.4f\n" d_corr_NAPL
else
    @printf "  corr(N,GDP) / corr(N,GDP/N): insufficient data, set to NaN.\n"
end


# =========================================================================== #
#  7a. REAL EXCHANGE RATE  (BIS Real Broad REER for Chile, RBCL)              #
#      reer_chile_bis.xlsx  Sheet "Real"                                      #
# =========================================================================== #

@printf "\n--- 7a. Real exchange rate (BIS REER) ---\n"

d_std_Q      = NaN
d_autocorr_Q = NaN
reer_hp_aligned = Float64[]

fname_reer = joinpath(DATA_DIR, "reer_chile_bis.xlsx")
if isfile(fname_reer)
    try
        wb_reer   = XLSX.readxlsx(fname_reer)
        ws_reer   = wb_reer["Real"]

        # Row 5 = country codes (RBDZ, RBAR, …, RBCL), row 6 onward = data
        # Find column with header "RBCL" (Chile)
        n_header_cols = XLSX.get_dimension(ws_reer).stop.column_number
        cl_col = nothing
        for c in 1:n_header_cols
            v = ws_reer[5, c]
            if !ismissing(v) && string(v) == "RBCL"
                cl_col = c
                break
            end
        end
        cl_col === nothing && (cl_col = 11)   # fallback: column K

        # Read monthly data from row 6 onward
        last_data_row = 6
        while !ismissing(ws_reer[last_data_row + 1, 1]) && last_data_row < 1000
            last_data_row += 1
        end
        n_reer = last_data_row - 5

        reer_dates = Date[]
        reer_vals  = Float64[]
        for t in 1:n_reer
            d_cell = ws_reer[5 + t, 1]
            v_cell = ws_reer[5 + t, cl_col]
            d_cell === missing && continue
            v_cell === missing && continue

            d_parsed = d_cell isa Date ? d_cell :
                       tryparse(Date, string(d_cell), dateformat"mm/dd/yyyy")
            d_parsed === nothing && continue

            v_parsed = v_cell isa Number ? Float64(v_cell) :
                       tryparse(Float64, string(v_cell))
            v_parsed === nothing && continue

            push!(reer_dates, d_parsed)
            push!(reer_vals,  v_parsed)
        end

        # Monthly → quarterly (average within each quarter)
        yr_reer = year.(reer_dates)
        qt_reer = ceil.(Int, month.(reer_dates) ./ 3)
        uq_reer = unique(collect(zip(yr_reer, qt_reer)))
        sort!(uq_reer)

        reer_qrt = [nanmean(reer_vals[findall(i -> (yr_reer[i], qt_reer[i]) == yq, eachindex(reer_vals))])
                    for yq in uq_reer]
        yr_reer_q  = [yq[1] for yq in uq_reer]
        qt_reer_q  = [yq[2] for yq in uq_reer]

        # Align to sample
        reer_s_mat, _, _ = align_sample(reshape(reer_qrt, :, 1),
                                          yr_reer_q, qt_reer_q,
                                          SAMPLE_START.year, SAMPLE_START.q,
                                          SAMPLE_END.year,   SAMPLE_END.q)
        reer_s = vec(reer_s_mat)

        if sum(!isnan, reer_s) >= 20
            reer_hp      = hp_cycle(log.(reer_s), LAMBDA)   # log-transform: std in fraction units
            d_std_Q      = mstd(reer_hp)
            # autocorr with the pandemic out: keep a (t,t+1) PAIR only if BOTH
            # quarters survive the mask, rather than deleting the excluded
            # quarters and splicing 2019Q4 onto 2021Q3 as a spurious transition.
            let k = covid_keep(length(reer_hp)),
                ok = isfinite.(reer_hp),
                pair = (k .& ok)[1:end-1] .& (k .& ok)[2:end]
                d_autocorr_Q = sum(pair) > 2 ?
                    cor(reer_hp[1:end-1][pair], reer_hp[2:end][pair]) : NaN
            end
            reer_hp_aligned = reer_hp
            @printf "  std(Q)      = %.5f  (BIS REER, %d quarterly obs)\n" d_std_Q sum(!isnan, reer_s)
            @printf "  autocorr(Q) = %.4f  (AR(1) of HP-filtered log REER)\n" d_autocorr_Q
        else
            @printf "  WARNING: insufficient REER obs in sample.\n"
        end
    catch e
        @printf "  WARNING: could not read reer_chile_bis.xlsx — %s\n" string(e)
    end
else
    @printf "  File not found: %s\n" fname_reer
end

if isnan(d_std_Q)
    d_std_Q      = 0.052
    d_autocorr_Q = 0.75
    @printf "  std(Q) / autocorr(Q): using fallback values (%.4f / %.4f)\n" d_std_Q d_autocorr_Q
end

d_corr_GDPQ = NaN


# =========================================================================== #
#  7b. TRADE BALANCE  (datos_CCNN_mayo2025.xlsx, sheet "Nom trimestral")      #
# =========================================================================== #

@printf "\n--- 7b. Trade balance (TB/GDP) ---\n"

d_TBGDP     = NaN
d_std_TBGDP = NaN
fname_ccnn = joinpath(DATA_DIR, "datos_CCNN_mayo2025.xlsx")

if isfile(fname_ccnn)
    try
        wb_tb  = XLSX.readxlsx(fname_ccnn)
        ws_tb  = wb_tb["Nom trimestral"]

        # Row 3 = headers, data from row 4.
        # Identify Export (col I≈9), Import (col J≈10), GDP (col K≈11)
        # Try to find by header name first
        n_cols_tb = 15
        i_x = i_m = i_g = nothing
        for c in 1:n_cols_tb
            cell_val = ws_tb[3, c]
            hdr = cell_val === missing ? "" : string(cell_val)
            occursin(r"[Ee]xport|x6_"i, hdr) && (i_x = c)
            occursin(r"[Ii]mport|x7_"i, hdr) && (i_m = c)
            occursin(r"PIB|x8_"i,        hdr) && (i_g = c)
        end
        i_x === nothing && (i_x = 9)
        i_m === nothing && (i_m = 10)
        i_g === nothing && (i_g = 11)

        # Parse date column (col A) and data columns
        last_tb = 4
        while !ismissing(ws_tb[last_tb + 1, 1]) && last_tb < 500
            last_tb += 1
        end
        n_tb  = last_tb - 3

        yr_tb = zeros(Int, n_tb)
        qt_tb = zeros(Int, n_tb)
        X_tb  = fill(NaN, n_tb)
        M_tb  = fill(NaN, n_tb)
        G_tb  = fill(NaN, n_tb)

        for t in 1:n_tb
            d_cell = ws_tb[3 + t, 1]
            if d_cell isa Date
                yr_tb[t] = year(d_cell)
                qt_tb[t] = ceil(Int, month(d_cell) / 3)
            elseif d_cell isa Number
                # Excel serial
                d_conv   = Date(1899, 12, 30) + Day(round(Int, d_cell))
                yr_tb[t] = year(d_conv)
                qt_tb[t] = ceil(Int, month(d_conv) / 3)
            end
            for (vec_, col) in [(X_tb, i_x), (M_tb, i_m), (G_tb, i_g)]
                v = ws_tb[3 + t, col]
                if v isa Number && !ismissing(v)
                    vec_[t] = Float64(v)
                end
            end
        end

        # Align to sample and compute TB/GDP
        mask_tb = [((y, q) >= (SAMPLE_START.year, SAMPLE_START.q)) &&
                   ((y, q) <= (SAMPLE_END.year,   SAMPLE_END.q))
                   for (y, q) in zip(yr_tb, qt_tb)]
        ok_tb = mask_tb .& (!isnan).(X_tb) .& (!isnan).(M_tb) .& (G_tb .> 0)

        if sum(ok_tb) >= 20
            tb_gdp_series = (X_tb[ok_tb] .- M_tb[ok_tb]) ./ G_tb[ok_tb]
            d_TBGDP     = mean(tb_gdp_series)
            # HP-filter the TB/GDP ratio (level, not log — can be negative) and take std dev
            d_std_TBGDP = mstd(hp_cycle(tb_gdp_series, LAMBDA))
            @printf "  TB/GDP      = %.4f  (BCCh CCNN, %d quarterly obs)\n" d_TBGDP sum(ok_tb)
            @printf "  std(TB/GDP) = %.5f  (HP-filtered)\n" d_std_TBGDP
        else
            @printf "  WARNING: insufficient CCNN obs for TB/GDP.\n"
        end
    catch e
        @printf "  WARNING: could not compute TB/GDP — %s\n" string(e)
    end
else
    @printf "  File not found: %s\n" fname_ccnn
end

isnan(d_TBGDP)     && (d_TBGDP = -0.02;  @printf "  TB/GDP: using fallback %.4f\n" d_TBGDP)
isnan(d_std_TBGDP) && (d_std_TBGDP = 0.025; @printf "  std(TB/GDP): using fallback %.4f\n" d_std_TBGDP)


# =========================================================================== #
#  7c. corr(GDP, Q)                                                           #
# =========================================================================== #

if !isempty(reer_hp_aligned) && !isempty(GDP_hp)
    nr = min(length(reer_hp_aligned), length(GDP_hp))
    nr >= 20 && (d_corr_GDPQ = mcor(GDP_hp[end-nr+1:end], reer_hp_aligned[end-nr+1:end]))
    @printf "  corr(GDP,Q)  = %.4f\n" d_corr_GDPQ
end
isnan(d_corr_GDPQ) && (d_corr_GDPQ = -0.15; @printf "  corr(GDP,Q): using fallback %.4f\n" d_corr_GDPQ)


# =========================================================================== #
#  8. GOODS EXPENDITURE SHARE (from calibration)                              #
# =========================================================================== #

d_omG = 0.57   # from BCCh CCNN calibration (ombar_val)
@printf "\n  omG = %.2f  (from SS calibration)\n" d_omG


# =========================================================================== #
#  7c. SOE GREAT RATIOS: consumption and investment  (moments 81-84)          #
# =========================================================================== #
#
# WHY (2026-08-21). Until now nothing in the 77-moment objective disciplined
# the COMPOSITION of demand — there was no consumption or investment moment at
# all, only GDP. That is a gap relative to every small-open-economy paper we
# compare to: Garcia-Cicco, Pancrazi & Uribe (2010 AER) Table 4 and the Central
# Bank's own XMAS Table 5 both report, for each of C and I, the triple
# {standard deviation, correlation with output, first-order autocorrelation}.
# With no such moment the model is free to get the level of GDP volatility right
# through any mix of C and I, and the aggregate block misbehaves accordingly.
#
# We target the two volatility RATIOS rather than the levels. std(C)/std(GDP)
# and std(I)/std(GDP) are the scale-free, conventionally reported objects, and
# they cannot be traded off against the overall shock scale the way levels can.
#
# SOURCE: datos_CCNN_mayo2025.xlsx, sheet "Real trimestral" — chained-volume
# quarterly national accounts. Column 3 is household consumption (2.1 Consumo de
# hogares e IPSFL), column 8 gross fixed capital formation (4. Formacion bruta
# de capital fijo), column 11 real GDP (7. Producto Interno Bruto). Note this is
# a DIFFERENT sheet from the nominal one used for the trade balance above.

@printf "\n--- 7c. Great ratios (real C and I, CCNN) ---\n"

d_ratio_stdC = NaN; d_ratio_stdI = NaN
d_corr_CGDP  = NaN; d_corr_IGDP  = NaN

if isfile(fname_ccnn)
    try
        ws_r = XLSX.readxlsx(fname_ccnn)["Real trimestral"]
        last_r = 4
        while !ismissing(ws_r[last_r + 1, 1]) && last_r < 500
            last_r += 1
        end
        n_r = last_r - 3
        yr_r = zeros(Int, n_r); qt_r = zeros(Int, n_r)
        C_r  = fill(NaN, n_r);  I_r = fill(NaN, n_r);  G_r = fill(NaN, n_r)
        for t in 1:n_r
            dv = ws_r[t + 3, 1]
            if dv isa Dates.Date || dv isa Dates.DateTime
                yr_r[t] = Dates.year(dv); qt_r[t] = ceil(Int, Dates.month(dv) / 3)
            else
                ss = strip(string(dv))
                length(ss) >= 7 && (yr_r[t] = something(tryparse(Int, ss[1:4]), 0);
                                    qt_r[t] = ceil(Int, something(tryparse(Int, ss[6:7]), 0) / 3))
            end
            _num(c) = let v = ws_r[t + 3, c]
                v === missing ? NaN : (v isa Number ? Float64(v) : something(tryparse(Float64, strip(string(v))), NaN))
            end
            C_r[t] = _num(3); I_r[t] = _num(8); G_r[t] = _num(11)
        end

        CIG, _, _ = align_sample(hcat(C_r, I_r, G_r), yr_r, qt_r,
                                 SAMPLE_START.year, SAMPLE_START.q,
                                 SAMPLE_END.year,   SAMPLE_END.q)
        ok_r = vec(all(x -> isfinite(x) && x > 0, CIG, dims=2))
        if sum(ok_r) >= 20
            r1 = findfirst(ok_r); r2 = findlast(ok_r)
            c_hp = hp_cycle(log.(CIG[r1:r2, 1]), LAMBDA)
            i_hp = hp_cycle(log.(CIG[r1:r2, 2]), LAMBDA)
            g_hp = hp_cycle(log.(CIG[r1:r2, 3]), LAMBDA)
            sg = mstd(g_hp)
            d_ratio_stdC = sg > 1e-10 ? mstd(c_hp) / sg : NaN
            d_ratio_stdI = sg > 1e-10 ? mstd(i_hp) / sg : NaN
            d_corr_CGDP  = mcor(c_hp, g_hp)
            d_corr_IGDP  = mcor(i_hp, g_hp)
            @printf "  Real CCNN sample: %d quarters (%dQ%d-%dQ%d)\n" (r2 - r1 + 1) yr_r[1] qt_r[1] SAMPLE_END.year SAMPLE_END.q
            @printf "  std(C)/std(GDP) = %.4f    corr(C,GDP) = %.4f\n" d_ratio_stdC d_corr_CGDP
            @printf "  std(I)/std(GDP) = %.4f    corr(I,GDP) = %.4f\n" d_ratio_stdI d_corr_IGDP
            @printf "  (this GDP is the CCNN expenditure-side volume measure, used only to\n"
            @printf "   normalise these four moments; std(GDP) itself stays on the BCCh series)\n"
        else
            @printf "  WARNING: fewer than 20 valid quarters in 'Real trimestral' — great ratios NaN.\n"
        end
    catch e
        @printf "  WARNING: could not read 'Real trimestral' (%s) — great ratios NaN.\n" sprint(showerror, e)
    end
else
    @printf "  %s not found — great ratios NaN.\n" basename(fname_ccnn)
end


# =========================================================================== #
#  9. OUTPUT-WEIGHTED CROSS-SECTIONAL AVERAGES                                #
# =========================================================================== #

@printf "\n--- 9. Output-weighted cross-sectional moments ---\n"

# SS sector-size weights, in order of preference (2026-07-08):
#   1. sector_calibration.csv, column Yi_ss (and PH_ss if present → value weights)
#   2. mod/params_jl.mod — parse Y_ss<i> and PH_ss<i>; weights = PH_ss·Y_ss
#      (constant-price VALUE shares, consistent with the model's constant-price
#      Y/VA aggregates; PH_ss_i ≠ 1, so raw quantities misweight sectors)
#   3. equal weights (loud warning — averages then NOT citable)
Y_ss_vec = ones(NSEC)
_wsrc = "equal (FALLBACK — do not cite the G/S averages)"
calib_file = joinpath(DATA_DIR, "sector_calibration.csv")
if isfile(calib_file)
    try
        calib_df = CSV.read(calib_file, DataFrame)
        _yi = hasproperty(calib_df, :Yi_ss) ? Float64.(calib_df.Yi_ss) :
              hasproperty(calib_df, :yi_ss) ? Float64.(calib_df.yi_ss) : nothing
        if _yi !== nothing
            _ph = hasproperty(calib_df, :PH_ss) ? Float64.(calib_df.PH_ss) : ones(NSEC)
            Y_ss_vec = _ph .* _yi
            _wsrc = "sector_calibration.csv (PH_ss·Yi_ss value shares)"
        end
    catch e
        @printf "  Could not read sector_calibration.csv (%s).\n" string(e)
    end
end
if startswith(_wsrc, "equal")
    params_mod = joinpath(SCRIPT_DIR, "mod", "params_jl.mod")
    if isfile(params_mod)
        _ys  = fill(NaN, NSEC); _phs = fill(NaN, NSEC)
        for ln in eachline(params_mod)
            m = match(r"^Y_ss(\d+)\s*=\s*([-0-9.eE+]+)\s*;", ln)
            m !== nothing && (k = parse(Int, m[1]); k <= NSEC && (_ys[k]  = parse(Float64, m[2])))
            m = match(r"^PH_ss(\d+)\s*=\s*([-0-9.eE+]+)\s*;", ln)
            m !== nothing && (k = parse(Int, m[1]); k <= NSEC && (_phs[k] = parse(Float64, m[2])))
        end
        if !any(isnan, _ys) && !any(isnan, _phs)
            Y_ss_vec = _phs .* _ys
            _wsrc = "params_jl.mod (PH_ss·Y_ss constant-price value shares)"
        else
            @printf "  params_jl.mod parsed but Y_ss/PH_ss incomplete (%d/%d, %d/%d).\n" sum(!isnan, _ys) NSEC sum(!isnan, _phs) NSEC
        end
    end
end
@printf "  Sector weights: %s\n" _wsrc
startswith(_wsrc, "equal") &&
    @printf "  WARNING: equal weights in use — regenerate params_jl.mod (run main_SOE_gap.jl)\n           or add Yi_ss/PH_ss columns to sector_calibration.csv.\n"

goods_idx    = GOODS
services_idx = SERVICES
w_g = Y_ss_vec[goods_idx]    ./ sum(Y_ss_vec[goods_idx])
w_s = Y_ss_vec[services_idx] ./ sum(Y_ss_vec[services_idx])

d_std_Yg  = dot(w_g, y_d[goods_idx])
d_std_PHg = dot(w_g, p_d[goods_idx])
d_std_Lg  = dot(w_g, l_d[goods_idx])
d_std_Ys  = dot(w_s, y_d[services_idx])
d_std_PHs = dot(w_s, p_d[services_idx])
d_std_Ls  = dot(w_s, l_d[services_idx])


# =========================================================================== #
#  9b. SECTORAL PRICE–OUTPUT CORRELATIONS  corr(Y_i, PH_i)                   #
# =========================================================================== #
# Sign interpretation:
#   corr < 0 → TFP shock dominates (supply ↑ → quantity ↑, price ↓)
#   corr > 0 → demand shock dominates (demand ↑ → quantity ↑, price ↑)
# This is the key moment for identifying the supply vs demand decomposition
# per sector in the NK-IOSOE Option-A estimation (Option A adds 12 sectoral
# demand shocks sigma_om_1:12).

@printf "\n--- 9b. corr(Y_i, PH_i): supply vs demand identification ---\n"

corr_YPH_d = fill(NaN, NSEC)
if !isempty(Y_qrt) && size(Y_qrt, 1) >= 20 && size(P_sample, 1) >= 20
    nT_Y = size(Y_qrt, 1)
    nT_P = size(P_sample, 1)
    nT_c = min(nT_Y, nT_P)   # both aligned to same sample period → should be equal
    for i in 1:NSEC
        xY = Y_qrt[end-nT_c+1:end, i]
        xP = P_sample[end-nT_c+1:end, i]
        bad = (xY .<= 0) .| isnan.(xY) .| (xP .<= 0) .| isnan.(xP)
        if sum(.!bad) < 20
            @printf "  WARNING: sector %d has too few valid obs for corr(Y,PH).\n" i
            continue
        end
        # Trim to the jointly valid window (sectors 8/10 start 2013Q1) so the
        # HP cycle is not distorted by constant-backfilled leading segments.
        i1 = findfirst(!, bad); i2 = findlast(!, bad)
        y_hp = hp_cycle(fillmissing_linear(log.(xY[i1:i2])), LAMBDA)
        p_hp = hp_cycle(fillmissing_linear(log.(xP[i1:i2])), LAMBDA)
        corr_YPH_d[i] = mcor(y_hp, p_hp)
    end
    @printf "  %-4s  %-20s  %10s\n" "Sec" "Name" "corr(Y,PH)"
    @printf "  %s\n" repeat("-", 38)
    for i in 1:NSEC
        tag = i in GOODS ? "[G]" : "[S]"
        @printf "  %-2d %-17s%s  %10.4f\n" i SECTOR_NAMES[i] tag corr_YPH_d[i]
    end
else
    @printf "  WARNING: Y_qrt unavailable — corr_YPH_d set to NaN (will use 0 in estimation).\n"
end


# =========================================================================== #
#  9c. IDIOSYNCRATIC SECTORAL MOMENTS — COMMON-FACTOR REMOVAL                 #
#      (2026-08-21; refinement 1 of compute_sectoral_shocks.jl's footer)      #
# =========================================================================== #
#
# WHY. std_Y / corr_YPH above are TOTAL sectoral volatility: the idiosyncratic
# sector-specific component PLUS the common factor (aggregate demand, external
# / terms-of-trade, monetary). compute_sectoral_shocks.jl inverts the SIZE of
# the sector-specific TFP and taste shocks from these numbers. Using the TOTAL
# double-counts, because the model ALSO has aggregate and external shocks that
# generate the common factor on top of the sectoral ones. That is what
# lambda_A = 0.59 in the 2026-08-21 estimation run was silently undoing.
#
# METHOD (Foerster, Sarte & Watson 2011 JPE, sec. III). Stack the HP-filtered
# sectoral output and price cycles into one T x 2N panel, standardise each
# column, take the first NFAC principal components, and project every series
# on them. The residual is the idiosyncratic component.
#
#   * The factors are estimated JOINTLY from [Y_panel P_panel] rather than
#     separately from each. The SAME factor is then removed from Y_i and P_i,
#     so corr(resid_Y_i, resid_P_i) is still a well-defined supply/demand
#     diagnostic. Removing a Y-factor from Y and a different P-factor from P
#     would contaminate exactly the correlation the inversion depends on.
#   * NFAC = 1 by default (set FACTOR_NFAC to override). One factor absorbs
#     ~33% of the panel variance; a second takes it to ~51% but starts to
#     soak up genuinely sectoral variation, which is what we are trying to
#     keep. Treat NFAC = 2 as the robustness check, not the baseline.
#
# OUTPUT. Three NEW columns in sectoral_moments.csv — std_Y_idio, std_PH_idio,
# corr_YPH_idio. The existing std_Y / std_PH / corr_YPH columns are UNCHANGED
# and remain the SMM targets: the model must reproduce TOTAL volatility, since
# it has both shock types. Only the shock INVERSION uses the _idio columns.

# NOTE: no `const` — this whole file body sits inside _main(), and `const` on a
# local is a syntax error in Julia. Same reason LAMBDA and SAMPLE_START above
# are plain assignments.
FACTOR_NFAC = parse(Int, get(ENV, "FACTOR_NFAC", "1"))

@printf "\n--- 9c. Idiosyncratic moments (common-factor removal, NFAC=%d) ---\n" FACTOR_NFAC

std_Y_idio    = fill(NaN, NSEC)
std_PH_idio   = fill(NaN, NSEC)
corr_YPH_idio = fill(NaN, NSEC)
# First-order autocorrelation of the IDIOSYNCRATIC output cycle (2026-08-24).
# Feeds the sector-specific TFP persistence rho_{A,i}; see the note where it is
# filled below.
ac_Y_idio     = fill(NaN, NSEC)
d_rbar_YY     = NaN   # moment 78: average pairwise cross-sectoral output correlation

"""
    common_factors(panel, k) -> (F, varshare)

First `k` principal components of `panel` (T x N), computed on the
column-standardised data. Returns the T x k factor scores and the share of
total panel variance carried by each component.
"""
function common_factors(panel::AbstractMatrix{<:Real}, k::Int)
    T, N = size(panel)
    Z = similar(Matrix{Float64}(panel))
    for j in 1:N
        col = @view panel[:, j]
        s   = std(col)
        Z[:, j] = s > 1e-12 ? (col .- mean(col)) ./ s : zeros(T)
    end
    U, S, _ = svd(Z)
    F = U[:, 1:k] .* S[1:k]'
    return F, (S .^ 2) ./ sum(S .^ 2)
end

"""
    project_out(x, F) -> residual

Residual from an OLS regression of `x` on a constant and the factors `F`.
"""
function project_out(x::AbstractVector{<:Real}, F::AbstractMatrix{<:Real})
    X = hcat(ones(length(x)), F)
    return x .- X * (X \ collect(float.(x)))
end

"""
    rbar_pairwise(panel, keep) -> (rbar, npairs)

Average pairwise correlation across the columns of a T x N panel, using only the
rows flagged in `keep`. This is the cross-sectional comovement statistic of the
production-network literature — FSW (2011 JPE) Table 7, Atalay (2017) eq. 18.

Factored out of the inline loop that used to compute moment 78, because that loop
called `complete_cor` directly and so was the ONLY sectoral moment in the file
that did not drop the pandemic quarters (fixed 2026-08-24).
"""
function rbar_pairwise(panel::AbstractMatrix{<:Real}, keep::AbstractVector{Bool})
    N = size(panel, 2)
    acc = 0.0; np = 0
    for i in 1:N, j in (i+1):N
        c = complete_cor(panel[keep, i], panel[keep, j])
        isnan(c) && continue
        acc += c; np += 1
    end
    return (np > 0 ? acc / np : NaN), np
end

"Share of total panel variance carried by the first principal component."
function pc1_share(panel::AbstractMatrix{<:Real})
    size(panel, 1) < 5 && return NaN
    _, vs = common_factors(panel, 1)
    return vs[1]
end

if !isempty(Y_qrt) && size(Y_qrt, 1) >= 20 && size(P_sample, 1) >= 20
    nT_Y = size(Y_qrt, 1); nT_P = size(P_sample, 1)
    nT_c = min(nT_Y, nT_P)
    Yc = Y_qrt[end-nT_c+1:end, :]
    Pc = P_sample[end-nT_c+1:end, :]

    # A factor model needs a BALANCED panel, so restrict to the window where
    # every sector has valid Y and P. Sectors 8 (Finance) and 10 (Business
    # services) are the binding constraint; the SA output vintage starts
    # 2013Q1, so in practice this is the whole output sample.
    ok = trues(nT_c)
    for i in 1:NSEC
        ok .&= .!((Yc[:, i] .<= 0) .| isnan.(Yc[:, i]) .|
                  (Pc[:, i] .<= 0) .| isnan.(Pc[:, i]))
    end
    t1 = findfirst(ok); t2 = findlast(ok)
    nbal = (t1 === nothing || t2 === nothing) ? 0 : sum(ok[t1:t2])

    if nbal < 20 || nbal < t2 - t1 + 1
        @printf "  WARNING: no contiguous balanced window (%d usable obs) — _idio columns left NaN.\n" nbal
        @printf "  → compute_sectoral_shocks.jl will fall back to TOTAL moments and warn.\n"
    else
        y_pan = hcat([hp_cycle(log.(Yc[t1:t2, i]), LAMBDA) for i in 1:NSEC]...)
        p_pan = hcat([hp_cycle(log.(Pc[t1:t2, i]), LAMBDA) for i in 1:NSEC]...)
        @printf "  Balanced panel: %d quarters x %d sectors (x2 blocks).\n" (t2 - t1 + 1) NSEC

        # `covid_keep(n)` dates its mask by counting BACKWARDS from SAMPLE_END,
        # which is only correct for a series that ends AT SAMPLE_END. The
        # balanced panel ends at t2, and t2 < nT_c whenever the last quarter is
        # unbalanced for any sector. In that case every mask applied to y_pan or
        # p_pan — mstd, mcor, and rbar — would drop the wrong quarters, silently.
        # Build the panel's own mask by padding out to SAMPLE_END and clipping.
        # (Guard added 2026-08-24; pre-existing exposure, not a new one.)
        pan_tail_lag = nT_c - t2
        pan_keep     = covid_keep(size(y_pan, 1) + pan_tail_lag)[1:size(y_pan, 1)]
        if pan_tail_lag > 0
            @printf "  NOTE: panel ends %d quarter(s) before SAMPLE_END; pandemic mask shifted accordingly.\n" pan_tail_lag
        end

        F, vshare = common_factors(hcat(y_pan, p_pan), FACTOR_NFAC)
        pc_txt = join([@sprintf("%.0f%%", 100 * vshare[k]) for k in 1:FACTOR_NFAC], ", ")
        @printf "  Variance share of retained PCs: %s  (cumulative %.0f%%)\n" pc_txt (100 * sum(vshare[1:FACTOR_NFAC]))

        std_Y_tot_pan  = [mstd(y_pan[:, i])                       for i in 1:NSEC]
        corr_tot_pan   = [mcor(y_pan[:, i], p_pan[:, i])  for i in 1:NSEC]
        for i in 1:NSEC
            ry = project_out(y_pan[:, i], F)
            rp = project_out(p_pan[:, i], F)
            std_Y_idio[i]    = mstd(ry)
            std_PH_idio[i]   = mstd(rp)
            corr_YPH_idio[i] = mcor(ry, rp)
            # PERSISTENCE OF THE IDIOSYNCRATIC CYCLE (2026-08-24).
            #
            # This measures rho_{A,i}, the sector-specific persistence of the TFP
            # shock, under exactly the same unit-loading assumption that
            # compute_sectoral_shocks.jl already makes for the shock SIZES: with
            # the common factor projected out and a unit shock-to-output loading,
            # the idiosyncratic output cycle IS the shock, so its autocorrelation
            # is the shock's. Sizes and persistence are then measured off the same
            # object under the same assumption, which is the point.
            #
            # It is deliberately NOT the total-cycle autocorr_Y already in this
            # file (moments 66-77). That one contains the common factor, whose
            # persistence the model generates from its own aggregate and external
            # shocks; feeding it back into rho_A would double-count exactly as
            # using total volatility for the sigmas would.
            ac_Y_idio[i]     = mac1(ry)
        end

        @printf "\n  %-4s %-18s %9s %9s %7s %10s %11s\n" "Sec" "Name" "std_Y" "std_Y_id" "share" "corr_YPH" "corr_id"
        @printf "  %s\n" repeat("-", 74)
        for i in 1:NSEC
            shr = std_Y_tot_pan[i] > 0 ? std_Y_idio[i] / std_Y_tot_pan[i] : NaN
            @printf "  %-2d   %-18s %9.5f %9.5f %6.2f %10.4f %11.4f\n" i SECTOR_NAMES[i] std_Y_tot_pan[i] std_Y_idio[i] shr corr_tot_pan[i] corr_YPH_idio[i]
        end
        @printf "  %s\n" repeat("-", 74)
        # AVERAGE PAIRWISE CROSS-SECTORAL OUTPUT CORRELATION (moment 78, 2026-08-21)
        #
        # rbar = mean over i<j of corr(Y_i, Y_j). This is THE cross-sectional
        # summary statistic of the production-network literature — Foerster,
        # Sarte & Watson (2011 JPE) Table 7 and Atalay (2017) eq. 18 both use it —
        # and it is what the input substitution elasticities govern: with
        # complementary inputs a sector-specific shock propagates and rbar is
        # high; with substitutable inputs buyers reallocate and rbar collapses.
        #
        # It replaces the three Spearman rank correlations in the objective. Those
        # have no precedent (the word "Spearman" does not appear in FSW, Atalay,
        # Rubbo, FGI, Baqaee-Farhi, Luttini-Pasten-Rubbo, GCPU, SGU or XMAS) and
        # they are step functions of theta, so locally flat and globally cliffed —
        # a bad object for CMA-ES. rbar is smooth, structural and standard.
        #
        # Computed on the TOTAL cycles, not the idiosyncratic residuals: the model
        # counterpart is corr(Y_i,Y_j) with every shock on, and removing a common
        # factor would drive it mechanically toward zero on the data side only.
        #
        # TWO CHANGES, 2026-08-24.
        #
        # (a) THE PANDEMIC MASK WAS MISSING. This loop called `complete_cor`
        #     directly rather than `mcor`, so rbar was the one sectoral moment
        #     still computed over 2020Q1-2021Q2. COVID is a textbook common
        #     factor — it moves all twelve sectors at once — so leaving it in
        #     inflates rbar, and the model was being asked to match a target
        #     contaminated by exactly the comovement it is accused of
        #     over-producing. Both versions are printed so the size of the
        #     artefact is on the record; the TARGETED value now excludes COVID.
        #
        # (b) A TRANSFORMATION DIAGNOSTIC. The entire "excess comovement"
        #     diagnosis rests on this one number, and it had never been checked
        #     against an alternative filter. Foerster, Sarte & Watson (2011 JPE)
        #     report that ~90% of US IP variation is common-factor driven, on
        #     QUARTERLY GROWTH RATES of 117 sectors; here it is HP-filtered log
        #     LEVELS of 12 broad sectors. If rbar and the PC1 share move a lot
        #     between the two transformations, the gap to the model is partly a
        #     measurement choice and not a structural failure. Growth rates are
        #     DIAGNOSTIC ONLY — the model counterpart is HP-filtered, so the
        #     targeted moment must stay HP-filtered.
        let T = size(y_pan, 1)
            keep_all = trues(T)
            keep_cov = pan_keep          # panel-dated mask, see the NOTE above

            rb_in,  np_in  = rbar_pairwise(y_pan, keep_all)
            rb_out, np_out = rbar_pairwise(y_pan, keep_cov)

            # Log growth rates. A (t-1, t) difference survives only if BOTH
            # quarters do — otherwise 2019Q4 gets spliced onto 2021Q3 and
            # counted as one quarter's growth, the same trap `mac1` avoids.
            ly     = log.(Yc[t1:t2, :])
            g_pan  = ly[2:end, :] .- ly[1:end-1, :]
            keep_g = keep_cov[1:end-1] .& keep_cov[2:end]
            rb_g, np_g = rbar_pairwise(g_pan, keep_g)

            d_rbar_YY = rb_out          # <-- moment 78, pandemic excluded

            @printf "\n  Average pairwise cross-sectoral output correlation (rbar)\n"
            @printf "  %s\n" repeat("-", 68)
            @printf "  %-34s %8s %8s %8s\n" "transformation" "rbar" "PC1" "pairs"
            @printf "  %-34s %8.4f %8.2f %8d\n" "HP log levels, COVID IN"  rb_in  pc1_share(y_pan[keep_all, :]) np_in
            @printf "  %-34s %8.4f %8.2f %8d  <- TARGETED\n" "HP log levels, COVID OUT" rb_out pc1_share(y_pan[keep_cov, :]) np_out
            @printf "  %-34s %8.4f %8.2f %8d  (diagnostic)\n" "log growth rates, COVID OUT" rb_g pc1_share(g_pan[keep_g, :]) np_g
            @printf "  %s\n" repeat("-", 68)
            @printf "  COVID contribution to rbar: %+.4f\n" (rb_in - rb_out)
            @printf "  FSW (2011 JPE) report ~0.90 PC1 share for US IP growth, 117 sectors.\n"
            @printf "  A large level-vs-growth gap here means the model/data comparison is\n"
            @printf "  partly a filter choice, not only a structural miss.\n"
        end

        nflip = count(i -> !isnan(corr_tot_pan[i]) && !isnan(corr_YPH_idio[i]) &&
                           corr_tot_pan[i] * corr_YPH_idio[i] < 0, 1:NSEC)
        @printf "  %d of %d sectors change the SIGN of corr(Y,PH) once the common factor is out.\n" nflip NSEC
        @printf "  (A positive total correlation driven by the common demand factor is NOT\n"
        @printf "   evidence of a sector-specific demand shock — that is the whole point.)\n"
    end
else
    @printf "  Y_qrt unavailable — _idio columns left NaN.\n"
end


# =========================================================================== #
#  10. SUMMARY AND VALIDATION                                                  #
# =========================================================================== #

@printf "\n%s\n  DATA MOMENTS SUMMARY\n%s\n\n" repeat("=",61) repeat("=",61)
@printf "%-20s  %8s  %8s  %8s  %10s\n" "Sector" "std(Y)" "std(PH)" "std(L)" "corr(Y,PH)"
@printf "%s\n" repeat("-", 66)
for i in 1:NSEC
    tag = i in GOODS ? "[G]" : "[S]"
    @printf "%-2d %-17s%s  %8.4f  %8.4f  %8.4f  %10.4f\n" i SECTOR_NAMES[i] tag y_d[i] p_d[i] l_d[i] corr_YPH_d[i]
end
@printf "%s\n" repeat("-", 66)
@printf "  Goods    avg: std(Y)=%6.4f  std(P)=%6.4f  std(L)=%6.4f\n" d_std_Yg  d_std_PHg  d_std_Lg
@printf "  Services avg: std(Y)=%6.4f  std(P)=%6.4f  std(L)=%6.4f\n\n" d_std_Ys  d_std_PHs  d_std_Ls
@printf "  std(GDP)      = %.5f\n" d_std_GDP
@printf "  std(pi)       = %.5f\n" d_std_pi
@printf "  corr(GDP,pi)  = %.4f\n" d_corr_GDPpi
@printf "  omG           = %.2f\n"  d_omG
@printf "  std(Q)        = %.5f\n" d_std_Q
@printf "  autocorr(Q)   = %.4f\n" d_autocorr_Q
@printf "  corr(GDP,Q)   = %.4f\n" d_corr_GDPQ
@printf "  TB/GDP        = %.4f\n" d_TBGDP
@printf "  std(TB/GDP)   = %.5f\n" d_std_TBGDP
@printf "  corr(N,GDP)   = %.4f\n" d_corr_NGDP
@printf "  corr(N,GDP/N) = %.4f\n" d_corr_NAPL

# Validate: all 47 moments must be non-NaN
d_vec_check = [y_d; p_d; l_d; d_std_GDP; d_std_pi; d_corr_GDPpi;
               d_omG; d_std_Q; d_TBGDP; d_std_TBGDP; d_autocorr_Q; d_corr_GDPQ;
               d_corr_NGDP; d_corr_NAPL]
nan_idx = findall(isnan, d_vec_check)
if !isempty(nan_idx)
    moment_labels = [["std(Y_$i)" for i in 1:NSEC];
                     ["std(PH_$i)" for i in 1:NSEC];
                     ["std(L_$i)"  for i in 1:NSEC];
                     ["std(GDP)","std(pi)","corr(GDP,pi)","omG",
                      "std(Q)","TB/GDP","std(TB/GDP)","autocorr(Q)","corr(GDP,Q)",
                      "corr(N,GDP)","corr(N,GDP/N)"]]

    # Build a readable list — included directly in the error() so it appears
    # in the exception message even when stdout has scrolled past.
    nan_lines = join(["  [$(k)] $(moment_labels[k])" for k in nan_idx], "\n")

    # Likely cause diagnostics
    causes = String[]
    any(isnan, y_d)  && push!(causes, "→ std(Y_i) NaN: pib_sectorial_bc.xlsx loading diagnostic:\n$(pib_error_msg)")
    any(isnan, p_d)  && push!(causes, "→ std(PH_i) NaN: deflactor_pib.csv columns didn't parse. Check file path and column count.")
    any(isnan, l_d)  && push!(causes, "→ std(L_i) NaN: count_workers_by_sector.csv missing or sector alignment wrong.")
    isnan(d_std_Q)   && push!(causes, "→ std(Q) NaN: reer_chile_bis.xlsx not found or RBCL column missing.")
    isnan(d_TBGDP)     && push!(causes, "→ TB/GDP NaN: datos_CCNN_mayo2025.xlsx not found or column mapping failed.")
    isnan(d_std_TBGDP) && push!(causes, "→ std(TB/GDP) NaN: TB/GDP series has fewer than 20 observations.")
    cause_lines = isempty(causes) ? "" : "\n\nLikely causes:\n" * join(causes, "\n")

    @printf "\n%s\n" repeat("!", 62)
    @printf "  NaN MOMENTS (%d of 47):\n" length(nan_idx)
    println(nan_lines)
    @printf "%s\n" repeat("!", 62)

    error("""
    $(length(nan_idx)) of 47 moments are NaN.

    NaN moments:
    $(nan_lines)
    $(cause_lines)

    DATA_DIR searched: $(DATA_DIR)
    Files present: $(filter(f -> endswith(f, r"\.xlsx|\.csv"), readdir(DATA_DIR, join=true) .|> basename))
    """)
end
@assert length(d_vec_check) == 47 "BUG: expected 47 moments, got $(length(d_vec_check))"
@printf "\nValidation passed: all 47 aggregate/sectoral moments are non-NaN.\n"

# Soft check for corr_YPH_d — NaN is acceptable if Y data is unavailable
n_nan_corr = sum(isnan.(corr_YPH_d))
if n_nan_corr > 0
    @printf "  NOTE: %d / %d corr(Y_i,PH_i) values are NaN (Y data unavailable).\n" n_nan_corr NSEC
    @printf "  → Estimation will treat these as 0 (via hasproperty fallback in smm_estimation.jl).\n"
else
    @printf "  corr(Y_i,PH_i) validation passed: all %d values non-NaN.\n" NSEC
end


# =========================================================================== #
#  SAVE — plain CSV (no .mat, no MATLAB dependency)                           #
# =========================================================================== #

# --- 1. Sectoral moments (one row per sector) ----------------------------
df_sec = DataFrame(
    sector   = 1:NSEC,
    name     = SECTOR_NAMES,
    std_Y    = y_d,
    std_PH   = p_d,
    std_L    = l_d,
    corr_YPH = corr_YPH_d,   # key identifier for supply vs demand decomposition (Option-A)
    # --- Idiosyncratic (common factor projected out) — section 9c ----------
    # These are the INVERSION inputs, not SMM targets. See 9c and the header of
    # compute_sectoral_shocks.jl.
    std_Y_idio    = std_Y_idio,
    std_PH_idio   = std_PH_idio,
    corr_YPH_idio = corr_YPH_idio,
    autocorr_Y    = ac_y,     # persistence moments (66-77) — see section 5
    # Sector-specific TFP persistence rho_{A,i}, measured on the idiosyncratic
    # cycle. Read by compute_sectoral_shocks.jl; NOT itself a targeted moment.
    autocorr_Y_idio = ac_Y_idio,
    std_Yg   = [i in GOODS    ? d_std_Yg  : NaN for i in 1:NSEC],
    std_PHg  = [i in GOODS    ? d_std_PHg : NaN for i in 1:NSEC],
    std_Lg   = [i in GOODS    ? d_std_Lg  : NaN for i in 1:NSEC],
    std_Ys   = [i in SERVICES ? d_std_Ys  : NaN for i in 1:NSEC],
    std_PHs  = [i in SERVICES ? d_std_PHs : NaN for i in 1:NSEC],
    std_Ls   = [i in SERVICES ? d_std_Ls  : NaN for i in 1:NSEC],
)
CSV.write(OUT_SECTORAL, df_sec)

# --- 2. Aggregate moments (key-value, easy to read in any language) ------
df_agg = DataFrame(
    moment = ["std_GDP",   "std_pi",    "corr_GDPpi",
              "omG",       "std_Q",     "autocorr_Q",
              "corr_GDPQ", "TBGDP",     "std_TBGDP",
              "corr_NGDP", "corr_NAPL",
              "autocorr_GDP", "autocorr_pi",
              "rbar_YY", "ratio_stdC", "ratio_stdI", "corr_CGDP", "corr_IGDP",
              "sample_start_year", "sample_start_q",
              "sample_end_year",   "sample_end_q",
              # Sample provenance (2026-08-21). Written so the estimation log,
              # and the paper's data section, can state the effective sample
              # without anyone having to re-derive it from the raw files.
              "covid_excluded", "covid_first_year", "covid_first_q",
              "covid_last_year", "covid_last_q",
              "sectoral_common_window", "n_sectoral_quarters"],
    value  = [d_std_GDP,   d_std_pi,    d_corr_GDPpi,
              d_omG,       d_std_Q,     d_autocorr_Q,
              d_corr_GDPQ, d_TBGDP,    d_std_TBGDP,
              d_corr_NGDP, d_corr_NAPL,
              d_ac_GDP,    d_ac_pi,
              d_rbar_YY, d_ratio_stdC, d_ratio_stdI, d_corr_CGDP, d_corr_IGDP,
              Float64(SAMPLE_START.year), Float64(SAMPLE_START.q),
              Float64(SAMPLE_END.year),   Float64(SAMPLE_END.q),
              COVID_EXCLUDE ? 1.0 : 0.0,
              Float64(COVID_FIRST[1]), Float64(COVID_FIRST[2]),
              Float64(COVID_LAST[1]),  Float64(COVID_LAST[2]),
              SECTORAL_COMMON_WINDOW ? 1.0 : 0.0,
              Float64(sum(covid_keep(nT)))],
)

# --- 2b. CARRY FORWARD the moments this script does not compute ----------- #
#
# Moments 61-63 — std_omG, std_pigap, corr_pigap_om — and the LEVEL omG are
# produced by Data/build_reallocation_calibration.py, not here. Until today this
# script simply overwrote aggregate_moments.csv without them, so deleting the
# file and regenerating it silently dropped three targeted moments and replaced
# omG with a hardcoded 0.57.
#
# That is not a hypothetical: it happened on 2026-08-21. The estimator caught
# std_omG and refused to run, but main_SOE_gap.jl ran anyway and printed a fit
# table in which std(pi_g - pi_s) was 26% of the objective — entirely because
# the data value was missing (0.0) and that moment carries inverse-squared-data
# weighting, so w = 1/0.01^2 = 10,000 against a model value of 0.024.
#
# So: preserve whatever the reallocation script last wrote, and say loudly when
# there is nothing to preserve.
let prior = Dict{String,Float64}()
    if isfile(OUT_AGGREGATE)
        try
            _old = CSV.read(OUT_AGGREGATE, DataFrame)
            for r in eachrow(_old)
                prior[String(r.moment)] = Float64(r.value)
            end
        catch e
            @printf "  WARNING: could not read the existing %s (%s) — reallocation moments not carried forward.\n" basename(OUT_AGGREGATE) sprint(showerror, e)
        end
    end

    carried  = String[]
    missing_ = String[]
    for k in ("omG", "std_omG", "std_pigap", "corr_pigap_om")
        if haskey(prior, k) && isfinite(prior[k])
            if k in df_agg.moment
                df_agg[findfirst(==(k), df_agg.moment), :value] = prior[k]
            else
                push!(df_agg, (k, prior[k]))
            end
            push!(carried, k)
        else
            push!(missing_, k)
        end
    end

    isempty(carried) || @printf "\n  Carried forward from the previous aggregate_moments.csv: %s\n" join(carried, ", ")
    if !isempty(missing_)
        @printf "\n  %s\n" repeat("!", 72)
        @printf "  MISSING moments this script does not compute: %s\n" join(missing_, ", ")
        @printf "  These are moments 61-63 (and the omG level). They come from a separate\n"
        @printf "  script. The SMM estimator will REFUSE to run without them, and any fit\n"
        @printf "  table built meanwhile will be badly distorted. Run now:\n\n"
        @printf "      python3 Data/build_reallocation_calibration.py\n\n"
        @printf "  It preserves every other row, so it is safe to run straight after this.\n"
        @printf "  %s\n" repeat("!", 72)
    end
end
CSV.write(OUT_AGGREGATE, df_agg)

@printf "\nMoments saved to:\n  %s\n  %s\n" OUT_SECTORAL OUT_AGGREGATE
@printf "\nAll moments computed from data. Ready to run smm_estimation.\n\n"

end  # function _main()

Base.invokelatest(_main)  # Julia 1.12: world-age fix
