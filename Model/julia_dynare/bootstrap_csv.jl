"""
bootstrap_csv.jl
================
One-time script: reads the existing MATLAB-generated .mat files and writes
clean CSV files that replace all .mat dependencies going forward.

Run ONCE:
  julia --project=. bootstrap_csv.jl

After this, delete or archive the .mat files — they are no longer needed.
"""

using MAT, CSV, DataFrames, Printf

SCRIPT_DIR = @__DIR__
MODELO_DIR = abspath(joinpath(SCRIPT_DIR, "..", "modelo_chile"))
DATA_DIR   = abspath(joinpath(SCRIPT_DIR, "..", "..", "Data"))

@printf "\n%s\nBootstrapping CSV files from MATLAB .mat outputs\n%s\n\n" repeat("=",52) repeat("=",52)

# Helper: extract a named vector from the MAT dict, always returning Vector{Float64}
function gv(d::Dict, k::String)
    v = d[k]
    v isa AbstractMatrix ? vec(Float64.(v)) :
    v isa AbstractVector ? Float64.(v)      : [Float64(v)]
end

# =========================================================================== #
#  1.  Data/sector_calibration.csv                                            #
#      Source: modelo_chile/params_val_ul.mat                                 #
#      Contains: alpha, alpha_V, var_rho, spend_good, spend_serv, kappa,     #
#                is_goods flag — everything read from Stata_to_excel*.xls     #
# =========================================================================== #

params_file = joinpath(MODELO_DIR, "params_val_ul.mat")
if !isfile(params_file)
    @printf "  SKIP sector_calibration.csv — %s not found\n" params_file
else
    p    = matread(params_file)
    nsec = 12

    sector_names = [
        "Agropecuario-silvicola y Pesca",
        "Mineria",
        "Industria manufacturera",
        "Electricidad gas agua y gestion de desechos",
        "Construccion",
        "Comercio hoteles y restaurantes",
        "Transporte comunicaciones y servicios de informacion",
        "Intermediacion financiera",
        "Servicios inmobiliarios y de vivienda",
        "Servicios empresariales",
        "Servicios personales",
        "Administracion publica",
    ]

    modalpha   = gv(p, "modalpha")
    modalphaV  = gv(p, "modalphaV")
    modvarrho  = gv(p, "modvarrho")
    modgammag  = gv(p, "modgammag")
    modgammas  = gv(p, "modgammas")
    modkappa   = gv(p, "modkappa")
    goods_v    = gv(p, "goods")
    ombar      = Float64(p["ombar_val"])

    # spend_good / spend_serv: reconstruct proportional basket weights
    # modgammag[i] = γ^g_i = share of sector i within the goods basket
    # spend_good[i] ∝ modgammag[i] × ombar  (goods share of total)
    spend_good = modgammag .* ombar
    spend_serv = modgammas .* (1.0 - ombar)

    df = DataFrame(
        sector     = 1:nsec,
        name       = sector_names,
        alpha      = modalpha,
        alpha_V    = modalphaV,
        var_rho    = modvarrho,
        spend_good = spend_good,
        spend_serv = spend_serv,
        kappa      = modkappa,
        is_goods   = Int.(round.(goods_v)),
    )

    out = joinpath(DATA_DIR, "sector_calibration.csv")
    CSV.write(out, df)
    @printf "  Written: %s\n" out
    @printf "  Columns: %s\n\n" join(names(df), ", ")
end

# =========================================================================== #
#  2.  Data/sectoral_moments.csv                                              #
#      Source: modelo_chile/data_moments_chile.mat                            #
#      Contains: std(Y_i), std(PH_i), std(L_i) for i = 1..12                #
# =========================================================================== #

dm_file = joinpath(MODELO_DIR, "data_moments_chile.mat")
if !isfile(dm_file)
    @printf "  SKIP sectoral_moments.csv — %s not found\n" dm_file
else
    dm_raw = matread(dm_file)
    dm     = dm_raw["dm_chile"]

    y_d = gv(dm, "y_d")
    p_d = gv(dm, "p_d")
    l_d = gv(dm, "l_d")
    nsec = length(y_d)

    df_sec = DataFrame(sector=1:nsec, std_Y=y_d, std_PH=p_d, std_L=l_d)
    out = joinpath(DATA_DIR, "sectoral_moments.csv")
    CSV.write(out, df_sec)
    @printf "  Written: %s\n" out
    show(df_sec); println("\n")

    # ---- aggregate moments ------------------------------------------------
    agg_keys = ["std_GDP","std_pi","corr_GDPpi","omG",
                "std_Q","autocorr_Q","corr_GDPQ","TBGDP"]
    agg_vals = Float64[
        dm["d_std_GDP"], dm["d_std_pi"], dm["d_corr_GDPpi"],
        dm["d_omG"],     dm["d_std_Q"],  dm["d_autocorr_Q"],
        dm["d_corr_GDPQ"], dm["d_TBGDP"],
    ]
    df_agg = DataFrame(moment=agg_keys, value=agg_vals)
    out2 = joinpath(DATA_DIR, "aggregate_moments.csv")
    CSV.write(out2, df_agg)
    @printf "  Written: %s\n" out2
    show(df_agg); println("\n")
end

@printf "\n%s\nBootstrap complete. MAT files are no longer needed.\n%s\n\n" repeat("=",52) repeat("=",52)
