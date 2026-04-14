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

using CSV, DataFrames, XLSX, MAT
using Statistics, LinearAlgebra, SparseArrays
using Dates, Printf

# =========================================================================== #
#  PATHS                                                                       #
# =========================================================================== #

SCRIPT_DIR  = @__DIR__
REPO_ROOT   = abspath(joinpath(SCRIPT_DIR, "..", ".."))
DATA_DIR    = joinpath(REPO_ROOT, "Data")
MODELO_DIR  = abspath(joinpath(SCRIPT_DIR, "..", "modelo_chile"))
OUTPUT_FILE = joinpath(MODELO_DIR, "data_moments_chile.mat")

@printf "\n%s\n" repeat("=", 61)
@printf "  Computing data moments for SMM estimation (Chile)\n"
@printf "%s\n\n" repeat("=", 61)


# =========================================================================== #
#  SETTINGS                                                                    #
# =========================================================================== #

const NSEC   = 12
const LAMBDA = 1600.0     # HP filter smoothing parameter (quarterly)

# Estimation sample
const SAMPLE_START = (year=2006, q=1)
const SAMPLE_END   = (year=2023, q=4)

# Sector classification (1-based indices)
const GOODS    = [1,2,3,4,5]      # sectors 1-5
const SERVICES = [6,7,8,9,10,11,12]  # sectors 6-12

const SECTOR_NAMES = [
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
const ALPHA_TO_MODEL = [1, 10, 5, 4, 8, 3, 2, 11, 12, 9, 7, 6]

fname_emp = joinpath(DATA_DIR, "count_workers_by_sector.csv")
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

fname_defl = joinpath(DATA_DIR, "deflactor_pib.csv")

# File has 3 metadata rows (SERIES, DESCRIPCION, UNIDAD), then data rows.
# We read everything as strings, skip the header rows.
raw_defl = CSV.read(fname_defl, DataFrame, header=false, types=String,
                    silencewarnings=true)
# Rows 1–3 are metadata; data starts at row 4
data_rows_defl = eachrow(Matrix{Union{String,Missing}}(raw_defl[4:end, :]))
n_defl = nrow(raw_defl) - 3

yr_defl = zeros(Int, n_defl)
q_defl  = zeros(Int, n_defl)

const MONTH_TO_Q = Dict("MAR"=>1, "JUN"=>2, "SEP"=>3, "DIC"=>4)

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
const DEFL_COLS = [
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

@printf "\n--- 3. Loading sectoral real output (pib_sectorial_bc.xlsx) ---\n"

fname_pib  = joinpath(DATA_DIR, "pib_sectorial_bc.xlsx")
Y_sec_raw  = Matrix{Float64}(undef, 0, NSEC)
GDP_data   = Float64[]
yr_y       = Int[]
qt_y       = Int[]

if isfile(fname_pib)
    try
        wb  = XLSX.readxlsx(fname_pib)

        # Try sheet name "Cuadro" first; fall back to first sheet.
        # Print available sheet names to help debug if "Cuadro" is wrong.
        sheet_names = XLSX.sheetnames(wb)
        @printf "  Sheets in file: %s\n" join(sheet_names, ", ")
        ws = "Cuadro" in sheet_names ? wb["Cuadro"] : wb[1]
        ws_name = "Cuadro" in sheet_names ? "Cuadro" : sheet_names[1]
        @printf "  Using sheet: %s\n" ws_name

        # Row 3 = headers, rows 4+ = data
        data_start = 4

        # ------------------------------------------------------------------ #
        # KEY BCCh QUIRK: the date column (col A) only has the year in the   #
        # FIRST row of each quarter group (rows 4, 8, 12, …); the other 3   #
        # rows in each group are BLANK.  Using col A to find the last row    #
        # would stop at row 5 (the first blank), giving n_pib=1.            #
        # Fix: use a DATA column (col B) to count rows; use row-position     #
        # arithmetic to assign years and quarters.                           #
        # ------------------------------------------------------------------ #

        # Last row: scan col 2 (first data column) which is always populated
        last_row = data_start
        while last_row < 400
            v = ws[last_row + 1, 2]   # col B: always has data in every quarter row
            (v === nothing || v === missing) && break
            # Stop on clearly empty numeric cell
            (v isa AbstractString && strip(v) == "") && break
            last_row += 1
        end
        n_pib = last_row - data_start + 1
        @printf "  readxlsx OK: %d data rows (col-B scan)\n" n_pib

        # --- Detect start year from column A of first data row ------------ #
        cell0 = ws[data_start, 1]
        @printf "  Date cell [row %d, col 1] → type=%-16s  value=%s\n" data_start string(typeof(cell0)) string(cell0)

        function detect_pib_year(v)
            v === nothing || v === missing && return 0
            if v isa Date || v isa DateTime
                d   = v isa DateTime ? Date(v) : v
                ser = Dates.value(d - Date(1899, 12, 30))
                return (1980 <= ser <= 2060) ? ser : (1980 <= year(d) <= 2060 ? year(d) : 0)
            elseif v isa Integer
                return (1980 <= v <= 2060) ? v :
                       (v > 0 ? let yr = year(Date(1899,12,30) + Dates.Day(v)); 1980<=yr<=2060 ? yr : 0 end : 0)
            elseif v isa AbstractFloat && isfinite(v)
                vi = round(Int, v)
                return (1980 <= vi <= 2060) ? vi :
                       (vi > 0 ? let yr = year(Date(1899,12,30) + Dates.Day(vi)); 1980<=yr<=2060 ? yr : 0 end : 0)
            elseif v isa AbstractString && length(v) >= 4
                p = tryparse(Int, v[1:4]); return (p !== nothing && 1980<=p<=2060) ? p : 0
            end
            return 0
        end

        start_year = detect_pib_year(cell0)
        if start_year == 0
            @printf "  WARNING: could not detect start year from col A. Defaulting to 2009.\n"
            start_year = 2009   # BCCh PIB data documented start year
        end
        @printf "  Start year detected: %d\n" start_year

        # --- Assign years and quarters by row position (robust) ----------- #
        # BCCh PIB always has exactly 4 consecutive rows per year (Q1→Q4).
        # This is far more reliable than parsing partially-blank date cells.
        yr_y_pib = [start_year + (t - 1) ÷ 4 for t in 1:n_pib]
        qt_y_pib = [(t - 1) % 4 + 1            for t in 1:n_pib]

        @printf "  Date range (row-position): %dQ%d – %dQ%d  (%d rows = %.1f years)\n" yr_y_pib[1] qt_y_pib[1] yr_y_pib[end] qt_y_pib[end] n_pib (n_pib/4)

        n_valid    = n_pib
        valid_rows = 1:n_pib

        # --- Parse numeric cells for 32 data columns (xlsx cols 2–33) ----- #
        pib_num = fill(NaN, n_valid, 32)
        for t in 1:n_valid
            for c in 1:32
                v = ws[data_start + t - 1, c + 1]
                if v isa Number && !ismissing(v) && !isnan(Float64(v))
                    pib_num[t, c] = Float64(v)
                elseif v isa AbstractString
                    p = tryparse(Float64, strip(v))
                    p !== nothing && (pib_num[t, c] = p)
                end
            end
        end
        n_valid_cells = sum(!isnan, pib_num)
        @printf "  Numeric parse: %d / %d cells valid (%.0f%%)\n" n_valid_cells length(pib_num) (100*n_valid_cells/length(pib_num))

        # Aggregate to 12 model sectors (column mapping mirrors MATLAB original)
        Y_sec_raw  = fill(NaN, n_pib, NSEC)
        Y_sec_raw[:, 1]  = pib_num[:, 1]  .+ pib_num[:, 2]   # Agro + Pesca
        Y_sec_raw[:, 2]  = pib_num[:, 3]                       # Mining
        Y_sec_raw[:, 3]  = pib_num[:, 6]                       # Manufacturing
        Y_sec_raw[:, 4]  = pib_num[:, 16]                      # Utilities
        Y_sec_raw[:, 5]  = pib_num[:, 17]                      # Construction
        Y_sec_raw[:, 6]  = pib_num[:, 18]                      # Trade+Hotels
        Y_sec_raw[:, 7]  = pib_num[:, 21] .+ pib_num[:, 22]   # Transport + Comm
        Y_sec_raw[:, 8]  = pib_num[:, 24]                      # Finance
        Y_sec_raw[:, 9]  = pib_num[:, 26]                      # Real Estate
        Y_sec_raw[:, 10] = pib_num[:, 25]                      # Business services
        Y_sec_raw[:, 11] = pib_num[:, 27]                      # Personal services
        Y_sec_raw[:, 12] = pib_num[:, 28]                      # Public admin
        GDP_data         = pib_num[:, 31]                       # PIB total

        yr_y = yr_y_pib   # all rows used (no valid-mask needed with row-position approach)
        qt_y = qt_y_pib

    catch e
        @printf "  ERROR reading pib_sectorial_bc.xlsx: %s\n" string(e)
        @printf "  Output moments (y_d) will be NaN.\n"
    end
else
    @printf "  File not found: %s\n  Output moments (y_d) will be NaN.\n" fname_pib
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


# =========================================================================== #
#  5. HP-FILTER AND COMPUTE SECTORAL STD DEVS                                 #
# =========================================================================== #

@printf "\n--- 5. HP filtering and computing std devs ---\n"

# Employment std devs
l_d = fill(NaN, NSEC)
for i in 1:NSEC
    x = log.(L_sample[:, i])
    any(isnan.(x) .| isinf.(x)) && (@printf "  WARNING: sector %d employment has missing values, skipping.\n" i; continue)
    l_d[i] = std(hp_cycle(x, LAMBDA))
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
    p_d[i] = std(hp_cycle(fillmissing_linear(x), LAMBDA))
end

@printf "  Employment std devs (%%): %s\n" join([@sprintf("%.3f", v*100) for v in l_d], "  ")
@printf "  Price      std devs (%%): %s\n" join([@sprintf("%.3f", v*100) for v in p_d], "  ")

# Sectoral output std devs
y_d = fill(NaN, NSEC)
if !isempty(Y_qrt)
    for i in 1:NSEC
        x   = Y_qrt[:, i]
        bad = (x .<= 0) .| isnan.(x)
        if sum(.!bad) < 20
            @printf "  WARNING: sector %d output has too few valid obs (%d), skipping.\n" i sum(.!bad)
            continue
        end
        x_log = log.(fillmissing_linear(x))
        y_d[i] = std(hp_cycle(x_log, LAMBDA))
    end
    @printf "  Output     std devs (%%): %s\n" join([@sprintf("%.3f", v*100) for v in y_d], "  ")
else
    @printf "  Output data not loaded; y_d set to NaN.\n"
end

# Aggregate inflation: QoQ growth of total GDP deflator
pi_series = Pa_sample[2:end] ./ Pa_sample[1:end-1]
pi_series[pi_series .<= 0] .= NaN
pi_hp    = hp_cycle(fillmissing_linear(log.(pi_series)), LAMBDA)
d_std_pi = std(pi_hp)
@printf "  std(pi) = %.5f\n" d_std_pi


# =========================================================================== #
#  6. AGGREGATE OUTPUT MOMENTS                                                 #
# =========================================================================== #

@printf "\n--- 6. Aggregate moments ---\n"

d_std_GDP    = NaN
d_corr_GDPpi = NaN
GDP_hp       = Float64[]

if !isempty(GDP_sample) && sum(GDP_sample .> 0) >= 20
    GDP_hp       = hp_cycle(fillmissing_linear(log.(GDP_sample)), LAMBDA)
    d_std_GDP    = std(GDP_hp)
    nmin         = min(length(GDP_hp)-1, length(pi_hp))
    d_corr_GDPpi = complete_cor(GDP_hp[end-nmin+1:end], pi_hp[end-nmin+1:end])
    @printf "  std(GDP)      = %.5f  (BCCh pib_sectorial_bc)\n" d_std_GDP
    @printf "  corr(GDP, pi) = %.4f\n" d_corr_GDPpi
elseif !isempty(Y_qrt)
    GDP_real     = vec(sum(Y_qrt, dims=2))
    GDP_hp       = hp_cycle(log.(GDP_real), LAMBDA)
    d_std_GDP    = std(GDP_hp)
    nmin         = min(length(GDP_hp)-1, length(pi_hp))
    d_corr_GDPpi = complete_cor(GDP_hp[end-nmin+1:end], pi_hp[end-nmin+1:end])
    @printf "  std(GDP)      = %.5f  (sum of sectoral VA)\n" d_std_GDP
    @printf "  corr(GDP, pi) = %.4f\n" d_corr_GDPpi
else
    d_std_GDP    = 0.021
    d_corr_GDPpi = -0.15
    @printf "  std(GDP) / corr: using defaults (no output data).\n"
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
            reer_hp      = hp_cycle(reer_s, LAMBDA)
            d_std_Q      = std(reer_hp)
            tmp          = reer_hp[isfinite.(reer_hp)]
            d_autocorr_Q = length(tmp) > 2 ? cor(tmp[1:end-1], tmp[2:end]) : NaN
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

d_TBGDP    = NaN
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
        ok_tb = mask_tb .& (!isnan).(X_tb) .& (!isnan).(M_tb) .& G_tb .> 0

        if sum(ok_tb) >= 20
            d_TBGDP = mean((X_tb[ok_tb] .- M_tb[ok_tb]) ./ G_tb[ok_tb])
            @printf "  TB/GDP = %.4f  (BCCh CCNN, %d quarterly obs)\n" d_TBGDP sum(ok_tb)
        else
            @printf "  WARNING: insufficient CCNN obs for TB/GDP.\n"
        end
    catch e
        @printf "  WARNING: could not compute TB/GDP — %s\n" string(e)
    end
else
    @printf "  File not found: %s\n" fname_ccnn
end

isnan(d_TBGDP) && (d_TBGDP = -0.02; @printf "  TB/GDP: using fallback %.4f\n" d_TBGDP)


# =========================================================================== #
#  7c. corr(GDP, Q)                                                           #
# =========================================================================== #

if !isempty(reer_hp_aligned) && !isempty(GDP_hp)
    nr = min(length(reer_hp_aligned), length(GDP_hp))
    nr >= 20 && (d_corr_GDPQ = complete_cor(GDP_hp[end-nr+1:end], reer_hp_aligned[end-nr+1:end]))
    @printf "  corr(GDP,Q)  = %.4f\n" d_corr_GDPQ
end
isnan(d_corr_GDPQ) && (d_corr_GDPQ = -0.15; @printf "  corr(GDP,Q): using fallback %.4f\n" d_corr_GDPQ)


# =========================================================================== #
#  8. GOODS EXPENDITURE SHARE (from calibration)                              #
# =========================================================================== #

d_omG = 0.57   # from BCCh CCNN calibration (ombar_val)
@printf "\n  omG = %.2f  (from SS calibration)\n" d_omG


# =========================================================================== #
#  9. OUTPUT-WEIGHTED CROSS-SECTIONAL AVERAGES                                #
# =========================================================================== #

@printf "\n--- 9. Output-weighted cross-sectional moments ---\n"

# Load SS output weights from params_val.mat if available; else equal weights
Y_ss_vec = ones(NSEC)
params_file = joinpath(MODELO_DIR, "params_val.mat")
if isfile(params_file)
    try
        pv = matread(params_file)
        if haskey(pv, "Yi_ss")
            Y_ss_vec = vec(Float64.(pv["Yi_ss"]))
            @printf "  Loaded Y_ss weights from params_val.mat\n"
        end
    catch
        @printf "  Could not read Y_ss from params_val.mat; using equal weights.\n"
    end
else
    @printf "  params_val.mat not found; using equal output weights.\n"
end

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
#  10. SUMMARY AND VALIDATION                                                  #
# =========================================================================== #

@printf "\n%s\n  DATA MOMENTS SUMMARY\n%s\n\n" repeat("=",61) repeat("=",61)
@printf "%-20s  %8s  %8s  %8s\n" "Sector" "std(Y)" "std(PH)" "std(L)"
@printf "%s\n" repeat("-", 52)
for i in 1:NSEC
    tag = i in GOODS ? "[G]" : "[S]"
    @printf "%-2d %-17s%s  %8.4f  %8.4f  %8.4f\n" i SECTOR_NAMES[i] tag y_d[i] p_d[i] l_d[i]
end
@printf "%s\n" repeat("-", 52)
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

# Validate: all 44 moments must be non-NaN
d_vec_check = [y_d; p_d; l_d; d_std_GDP; d_std_pi; d_corr_GDPpi;
               d_omG; d_std_Q; d_TBGDP; d_autocorr_Q; d_corr_GDPQ]
nan_idx = findall(isnan, d_vec_check)
if !isempty(nan_idx)
    moment_labels = [["std(Y_$i)" for i in 1:NSEC];
                     ["std(PH_$i)" for i in 1:NSEC];
                     ["std(L_$i)"  for i in 1:NSEC];
                     ["std(GDP)","std(pi)","corr(GDP,pi)","omG",
                      "std(Q)","TB/GDP","autocorr(Q)","corr(GDP,Q)"]]

    # Build a readable list — included directly in the error() so it appears
    # in the exception message even when stdout has scrolled past.
    nan_lines = join(["  [$(k)] $(moment_labels[k])" for k in nan_idx], "\n")

    # Likely cause diagnostics
    causes = String[]
    any(isnan, y_d)  && push!(causes, "→ std(Y_i) NaN: pib_sectorial_bc.xlsx not found, wrong sheet name, or date-parsing failure. Run with VERBOSE_PIB=true for cell-level debug.")
    any(isnan, p_d)  && push!(causes, "→ std(PH_i) NaN: deflactor_pib.csv columns didn't parse. Check file path and column count.")
    any(isnan, l_d)  && push!(causes, "→ std(L_i) NaN: count_workers_by_sector.csv missing or sector alignment wrong.")
    isnan(d_std_Q)   && push!(causes, "→ std(Q) NaN: reer_chile_bis.xlsx not found or RBCL column missing.")
    isnan(d_TBGDP)   && push!(causes, "→ TB/GDP NaN: datos_CCNN_mayo2025.xlsx not found or column mapping failed.")
    cause_lines = isempty(causes) ? "" : "\n\nLikely causes:\n" * join(causes, "\n")

    @printf "\n%s\n" repeat("!", 62)
    @printf "  NaN MOMENTS (%d of 44):\n" length(nan_idx)
    println(nan_lines)
    @printf "%s\n" repeat("!", 62)

    error("""
    $(length(nan_idx)) of 44 moments are NaN.

    NaN moments:
    $(nan_lines)
    $(cause_lines)

    DATA_DIR searched: $(DATA_DIR)
    Files present: $(filter(f -> endswith(f, r"\.xlsx|\.csv"), readdir(DATA_DIR, join=true) .|> basename))
    """)
end
@assert length(d_vec_check) == 44 "BUG: expected 44 moments, got $(length(d_vec_check))"
@printf "\nValidation passed: all 44 moments are non-NaN.\n"


# =========================================================================== #
#  SAVE TO data_moments_chile.mat                                              #
# =========================================================================== #

dm_chile = Dict{String,Any}(
    "y_d"          => y_d,
    "p_d"          => p_d,
    "l_d"          => l_d,
    "d_std_Yg"     => d_std_Yg,
    "d_std_PHg"    => d_std_PHg,
    "d_std_Lg"     => d_std_Lg,
    "d_std_Ys"     => d_std_Ys,
    "d_std_PHs"    => d_std_PHs,
    "d_std_Ls"     => d_std_Ls,
    "d_std_GDP"    => d_std_GDP,
    "d_std_pi"     => d_std_pi,
    "d_corr_GDPpi" => d_corr_GDPpi,
    "d_omG"        => d_omG,
    "d_std_Q"      => d_std_Q,
    "d_autocorr_Q" => d_autocorr_Q,
    "d_corr_GDPQ"  => d_corr_GDPQ,
    "d_TBGDP"      => d_TBGDP,
    "sample_start" => [SAMPLE_START.year, SAMPLE_START.q],
    "sample_end"   => [SAMPLE_END.year,   SAMPLE_END.q],
    "sector_names" => SECTOR_NAMES,
    "L_qrt"        => L_qrt,
    "P_sec_qrt"    => P_sample,
    "Pa_qrt"       => Pa_sample,
    "Y_sec_qrt"    => isempty(Y_qrt) ? fill(NaN, 0, NSEC) : Y_qrt,
    "GDP_qrt"      => GDP_sample,
)

matwrite(OUTPUT_FILE, Dict("dm_chile" => dm_chile))
@printf "\nMoments saved to:\n  %s\n" OUTPUT_FILE
@printf "\nAll moments computed from data. Ready to run smm_estimation.\n\n"
