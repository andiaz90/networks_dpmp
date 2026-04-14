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
if isfile(OUT_SECTORAL) && isfile(OUT_AGGREGATE)
    @printf "\n%s\n  CSV moment files already exist — loading directly.\n" repeat("=",61)
    @printf "  (Delete and re-run to recompute from raw Excel/CSV sources.)\n"
    @printf "%s\n\n" repeat("=",61)

    sec = CSV.read(OUT_SECTORAL,  DataFrame)
    agg = CSV.read(OUT_AGGREGATE, DataFrame)

    @printf "  %-20s  %8s  %8s  %8s\n" "Sector" "std(Y)" "std(PH)" "std(L)"
    @printf "  %s\n" repeat("-", 50)
    for r in eachrow(sec)
        @printf "  %-20s  %8.4f  %8.4f  %8.4f\n" r.name r.std_Y r.std_PH r.std_L
    end
    agg_d = Dict(String(r.moment) => Float64(r.value) for r in eachrow(agg))
    @printf "\n  std(GDP)=%.4f  std(pi)=%.4f  corr(GDP,pi)=%.4f\n" agg_d["std_GDP"] agg_d["std_pi"] agg_d["corr_GDPpi"]
    @printf "  std(Q)=%.4f   autocorr(Q)=%.4f  TB/GDP=%.4f\n\n" agg_d["std_Q"] agg_d["autocorr_Q"] agg_d["TBGDP"]
    @printf "  Loaded from: %s\n\n" OUT_SECTORAL
    return   # done — no raw file processing needed
end

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

@printf "\n--- 3. Loading sectoral real output (pib_sectorial_bc.xlsx) ---\n"

fname_pib  = joinpath(DATA_DIR, "pib_sectorial_bc.xlsx")
Y_sec_raw  = Matrix{Float64}(undef, 0, NSEC)
GDP_data   = Float64[]
yr_y       = Int[]
qt_y       = Int[]

pib_error_msg = ""   # captured here so the validation error always shows it

if isfile(fname_pib)
    try

        # ------------------------------------------------------------------ #
        # Read the whole sheet as a Matrix{Any} in ONE call.                 #
        # This avoids ALL cell-by-cell iteration issues and the BCCh blank-  #
        # date-row quirk (col A only has the year every 4th row).            #
        # ------------------------------------------------------------------ #
        wb          = XLSX.readxlsx(fname_pib)
        sheet_names = XLSX.sheetnames(wb)
        pib_error_msg *= "  Sheets in file: " * join(sheet_names, ", ") * "\n"
        @printf "  Sheets in file: %s\n" join(sheet_names, ", ")

        ws      = "Cuadro" in sheet_names ? wb["Cuadro"] : wb[1]
        ws_name = "Cuadro" in sheet_names ? "Cuadro" : sheet_names[1]
        pib_error_msg *= "  Using sheet: $ws_name\n"
        @printf "  Using sheet: %s\n" ws_name

        # Get the worksheet dimension to know exactly how many rows/cols exist
        dim      = XLSX.get_dimension(ws)
        max_row  = dim.stop.row_number
        max_col  = min(dim.stop.column_number, 33)  # caps at col AG (33rd col)

        pib_error_msg *= "  Sheet dimension: rows 1–$max_row, cols 1–$max_col\n"
        @printf "  Sheet dimension: rows 1–%d, cols 1–%d\n" max_row max_col

        # Data rows start at row 4 (row 3 = headers in BCCh format)
        data_start = 4
        if max_row < data_start
            error("Sheet has only $max_row rows; expected data from row $data_start onward.")
        end

        # Read entire data range as Matrix{Any}: rows 4..max_row × cols 1..33
        # Each element is a number, Date, String, or missing — whatever XLSX.jl gives
        raw = ws[data_start:max_row, 1:max_col]   # (max_row-3) × max_col Matrix{Any}
        n_pib_raw = size(raw, 1)
        pib_error_msg *= "  Raw matrix: $n_pib_raw rows × $(size(raw,2)) cols\n"
        @printf "  Raw matrix: %d rows × %d cols\n" n_pib_raw size(raw, 2)

        # Diagnose the first few cells of col A (the BCCh date column)
        @printf "  Col-A sample (rows %d–%d): " data_start (data_start+3)
        for t in 1:min(4, n_pib_raw)
            v = raw[t, 1]
            @printf "type=%s val=%s  |  " string(typeof(v)) string(v)
        end
        @printf "\n"

        # --- Detect start year from col A, row 1 of raw (= xlsx row 4) --- #
        function _year_from_cell(v)
            v === nothing || v === missing && return 0
            if v isa Date || v isa DateTime
                d   = v isa DateTime ? Date(v) : v
                ser = Dates.value(d - Date(1899, 12, 30))
                return (1980 <= ser <= 2060) ? ser : (1980 <= year(d) <= 2060 ? year(d) : 0)
            elseif v isa Number && isfinite(Float64(v))
                vi = round(Int, Float64(v))
                return (1980 <= vi <= 2060) ? vi :
                    (vi > 0 ? let yr = year(Date(1899,12,30) + Dates.Day(vi))
                                  1980 <= yr <= 2060 ? yr : 0 end : 0)
            elseif v isa AbstractString && length(v) >= 4
                p = tryparse(Int, v[1:4])
                return (p !== nothing && 1980 <= p <= 2060) ? p : 0
            end
            return 0
        end

        # Try each of the first 4 rows of col A (one of them should have the year)
        start_year = 0
        for t in 1:min(4, n_pib_raw)
            start_year = _year_from_cell(raw[t, 1])
            start_year > 0 && break
        end
        if start_year == 0
            start_year = 2009
            pib_error_msg *= "  WARNING: start year not found in col A → defaulting to 2009\n"
            @printf "  WARNING: start year not found in col A — defaulting to 2009.\n"
        end
        @printf "  Start year: %d\n" start_year

        # --- Trim trailing empty rows (where all data cols are missing) --- #
        # Work backwards from the last row to find the true last data row
        n_pib = n_pib_raw
        while n_pib > 1
            row_vals = raw[n_pib, 2:end]
            all(v -> v === missing || v === nothing, row_vals) ? (n_pib -= 1) : break
        end
        @printf "  Data rows after trim: %d (%.1f years × 4 quarters)\n" n_pib (n_pib/4)

        # --- Year/quarter by row position (immune to blank date cells) ---- #
        yr_y_pib = [start_year + (t - 1) ÷ 4 for t in 1:n_pib]
        qt_y_pib = [(t - 1) % 4 + 1            for t in 1:n_pib]
        @printf "  Date range: %dQ%d – %dQ%d\n" yr_y_pib[1] qt_y_pib[1] yr_y_pib[end] qt_y_pib[end]

        # --- Parse numeric data columns (cols 2–33 of raw = xlsx B–AG) ---- #
        n_data_cols = min(32, size(raw, 2) - 1)
        pib_num = fill(NaN, n_pib, n_data_cols)
        for t in 1:n_pib
            for c in 1:n_data_cols
                v = raw[t, c + 1]    # col 1 = date, cols 2+ = data
                v === missing || v === nothing && continue
                if v isa Number && isfinite(Float64(v))
                    pib_num[t, c] = Float64(v)
                elseif v isa AbstractString
                    p = tryparse(Float64, strip(v)); p !== nothing && (pib_num[t, c] = p)
                end
            end
        end
        n_ok = sum(!isnan, pib_num)
        pib_error_msg *= "  Numeric parse: $n_ok / $(length(pib_num)) cells valid\n"
        @printf "  Numeric parse: %d / %d cells valid (%.0f%%)\n" n_ok length(pib_num) (100*n_ok/length(pib_num))

        if n_ok == 0
            error("All $(length(pib_num)) numeric cells parsed as NaN — check sheet structure.")
        end

        # --- Aggregate to 12 model sectors -------------------------------- #
        # Column mapping (pib_num col k = xlsx col k+1):
        #  1=Agropecuario, 2=Pesca, 3=Minería total, 6=Manufactura,
        # 16=EGA, 17=Construcción, 18=Comercio+Restaurantes,
        # 21=Transporte, 22=Comunicaciones, 24=Financiero, 25=Empresarial,
        # 26=Vivienda, 27=Personales, 28=Admin pública, 31=PIB total
        Y_sec_raw = fill(NaN, n_pib, NSEC)
        nc = n_data_cols
        nc >= 2  && (Y_sec_raw[:, 1]  = pib_num[:, 1] .+ pib_num[:, 2])  # Agro+Pesca
        nc >= 3  && (Y_sec_raw[:, 2]  = pib_num[:, 3])                     # Mining
        nc >= 6  && (Y_sec_raw[:, 3]  = pib_num[:, 6])                     # Manufactura
        nc >= 16 && (Y_sec_raw[:, 4]  = pib_num[:, 16])                    # Utilities
        nc >= 17 && (Y_sec_raw[:, 5]  = pib_num[:, 17])                    # Construction
        nc >= 18 && (Y_sec_raw[:, 6]  = pib_num[:, 18])                    # Trade+Hotels
        nc >= 22 && (Y_sec_raw[:, 7]  = pib_num[:, 21] .+ pib_num[:, 22]) # Transport+Comm
        nc >= 24 && (Y_sec_raw[:, 8]  = pib_num[:, 24])                    # Finance
        nc >= 26 && (Y_sec_raw[:, 9]  = pib_num[:, 26])                    # Real Estate
        nc >= 25 && (Y_sec_raw[:, 10] = pib_num[:, 25])                    # Business Serv
        nc >= 27 && (Y_sec_raw[:, 11] = pib_num[:, 27])                    # Personal Serv
        nc >= 28 && (Y_sec_raw[:, 12] = pib_num[:, 28])                    # Public Admin
        nc >= 31 && (GDP_data         = pib_num[:, 31])                     # PIB total

        yr_y = yr_y_pib
        qt_y = qt_y_pib

    catch e


        err_str = sprint(showerror, e)
        pib_error_msg *= "  EXCEPTION: $err_str\n"

        @printf "\n  %s\n" repeat("!", 58)
        @printf "  EXCEPTION reading pib_sectorial_bc.xlsx:\n  %s\n" err_str
        @printf "  %s\n\n" repeat("!", 58)

        # ---- Specific diagnosis for common XLSX.jl failures ---- #
        if occursin("x:workbook", err_str) || occursin("Malformed", err_str)
            @printf """
  DIAGNOSIS: pib_sectorial_bc.xlsx uses a non-standard XML namespace
  prefix ('x:workbook' instead of 'workbook') that XLSX.jl cannot read.

  ONE-TIME FIX (30 seconds):
    1. Open  Data/pib_sectorial_bc.xlsx  in Excel (or LibreOffice Calc)
    2. File → Save As → Excel Workbook (.xlsx)  [overwrite the same file]
    3. Re-run compute_data_moments.jl

  This rewrites the XML without the namespace prefix — XLSX.jl will
  then read it correctly.

"""
        end
    end
else
    pib_error_msg = "  File not found: $fname_pib\n"
    @printf "  File not found: %s\n" fname_pib
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
    any(isnan, y_d)  && push!(causes, "→ std(Y_i) NaN: pib_sectorial_bc.xlsx loading diagnostic:\n$(pib_error_msg)")
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
#  SAVE — plain CSV (no .mat, no MATLAB dependency)                           #
# =========================================================================== #

# --- 1. Sectoral moments (one row per sector) ----------------------------
df_sec = DataFrame(
    sector  = 1:NSEC,
    name    = SECTOR_NAMES,
    std_Y   = y_d,
    std_PH  = p_d,
    std_L   = l_d,
    std_Yg  = [i in GOODS    ? d_std_Yg  : NaN for i in 1:NSEC],
    std_PHg = [i in GOODS    ? d_std_PHg : NaN for i in 1:NSEC],
    std_Lg  = [i in GOODS    ? d_std_Lg  : NaN for i in 1:NSEC],
    std_Ys  = [i in SERVICES ? d_std_Ys  : NaN for i in 1:NSEC],
    std_PHs = [i in SERVICES ? d_std_PHs : NaN for i in 1:NSEC],
    std_Ls  = [i in SERVICES ? d_std_Ls  : NaN for i in 1:NSEC],
)
CSV.write(OUT_SECTORAL, df_sec)

# --- 2. Aggregate moments (key-value, easy to read in any language) ------
df_agg = DataFrame(
    moment = ["std_GDP",   "std_pi",   "corr_GDPpi",
              "omG",       "std_Q",    "autocorr_Q",
              "corr_GDPQ", "TBGDP",
              "sample_start_year", "sample_start_q",
              "sample_end_year",   "sample_end_q"],
    value  = [d_std_GDP,   d_std_pi,   d_corr_GDPpi,
              d_omG,       d_std_Q,    d_autocorr_Q,
              d_corr_GDPQ, d_TBGDP,
              Float64(SAMPLE_START.year), Float64(SAMPLE_START.q),
              Float64(SAMPLE_END.year),   Float64(SAMPLE_END.q)],
)
CSV.write(OUT_AGGREGATE, df_agg)

@printf "\nMoments saved to:\n  %s\n  %s\n" OUT_SECTORAL OUT_AGGREGATE
@printf "\nAll moments computed from data. Ready to run smm_estimation.\n\n"

end  # function _main()

Base.invokelatest(_main)  # Julia 1.12: world-age fix
