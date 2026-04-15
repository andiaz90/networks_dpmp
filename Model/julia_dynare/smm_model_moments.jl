"""
smm_model_moments.jl
====================
Compute theoretical model moments for parameter vector θ.

PLATFORM NOTES:
  Uses GenericSchur.jl (pure Julia)
  for the ordered QZ decomposition (no LAPACK dependency).
  The compiled Dynare Jacobian files (SparseDynamicG1!.jl etc.)
  are called directly — no Dynare.jl re-solve API needed.

KLEIN (2000) IMPLEMENTATION:
  After assembling the 491×642 dynamic Jacobian G = [A|B|C|D] from
  the compiled model files, applies:
  1. GenericSchur.schur(AA, BB) for generalized QZ (no LAPACK callback)
  2. GenericSchur.ordschur(F, select) for eigenvalue ordering (pure Julia)
  3. Extracts decision rule g1_1, g1_2 via Klein's Z-matrix partition

JACOBIAN COLUMN LAYOUT (from dynamic.json):
  cols   1: 78  → A: ∂f/∂y_{bk,t-1}   (78 backward state variables)
  cols  79:569  → B: ∂f/∂y_t           (491 current endogenous)
  cols 570:625  → C: ∂f/∂y_{fw,t+1}   (56 forward-looking variables)
  cols 626:642  → D: ∂f/∂ε_t           (17 exogenous shocks)
"""

using LinearAlgebra, SparseArrays, Statistics, StatsBase, NLsolve
# GenericSchur was used to avoid LAPACK gees callbacks on ARM Mac.
# Julia's built-in LinearAlgebra.schur uses LAPACK dgges (no select callback)
# which IS available via Apple Accelerate on ARM. We now use it directly.

include("steady_ntwsoe_system.jl")
include("steady_ntwsoe.jl")
include("utils.jl")

# =========================================================================== #
#  PATHS TO COMPILED MODEL FILES                                               #
# =========================================================================== #
const _JDYN_DIR    = @__DIR__
const _MODEL_BASE  = joinpath(_JDYN_DIR, "mod", "NK_SOE_lev_gap2")
const _JULIA_DIR   = joinpath(_MODEL_BASE, "model", "julia")
const _JSON_PATH   = joinpath(_MODEL_BASE, "model", "json", "dynamic.json")
const _MODFILE_PATH = joinpath(_MODEL_BASE, "model", "json", "modfile.json")

# Include compiled Dynare model files at TOP LEVEL so they are in the same
# Julia world as all calling code.  Including them inside a function creates
# world-age gaps that cause silent failures on all platforms.
const _DYNARE_MODEL_LOADED = Ref(false)
if isdir(_JULIA_DIR)
    include(joinpath(_JULIA_DIR, "SparseDynamicResidTT!.jl"))
    include(joinpath(_JULIA_DIR, "SparseDynamicResid!.jl"))
    include(joinpath(_JULIA_DIR, "SparseDynamicG1TT!.jl"))
    include(joinpath(_JULIA_DIR, "SparseDynamicG1!.jl"))
    _DYNARE_MODEL_LOADED[] = true
end

# =========================================================================== #
#  JACOBIAN SPARSITY STRUCTURE  (read once from dynamic.json)                 #
# =========================================================================== #

const _JAC_LOADED = Ref(false)
const _JAC_ROWS   = Int[]
const _JAC_COLS   = Int[]
const _N_NZ       = Ref(0)
const _N_EQ       = Ref(0)
const _N_COL      = Ref(0)

function _load_jacobian_structure!()
    _JAC_LOADED[] && return
    isfile(_JSON_PATH) || error("dynamic.json not found: $_JSON_PATH\nRun main_SOE_gap.jl first.")

    json_raw = read(_JSON_PATH, String)
    # Parse manually — avoid JSON.jl dependency
    # Extract "eq" and "col" from each entry using simple regex
    eq_matches  = collect(eachmatch(r"\"eq\":\s*(\d+)",  json_raw))
    col_matches = collect(eachmatch(r"\"col\":\s*(\d+)", json_raw))
    nrows_m = match(r"\"nrows\":\s*(\d+)", json_raw)
    ncols_m = match(r"\"ncols\":\s*(\d+)", json_raw)

    _N_EQ[]  = parse(Int, nrows_m.captures[1])   # 491
    _N_COL[] = parse(Int, ncols_m.captures[1])   # 642
    n        = min(length(eq_matches), length(col_matches))
    _N_NZ[]  = n

    resize!(_JAC_ROWS, n)
    resize!(_JAC_COLS, n)
    for i in 1:n
        _JAC_ROWS[i] = parse(Int, eq_matches[i].captures[1])
        _JAC_COLS[i] = parse(Int, col_matches[i].captures[1])
    end
    _JAC_LOADED[] = true
    @printf "  [SMM] Jacobian sparsity loaded: %d equations × %d cols, %d nnz\n" _N_EQ[] _N_COL[] _N_NZ[]
end


# =========================================================================== #
#  FORWARD VARIABLE INDICES FROM modfile.json                                  #
# =========================================================================== #
# modfile.json contains model_info.lead_lag_incidence: a (n_endo × 3) matrix.
# Row j, column 3 gives the Jacobian column of variable j's t+1 appearance
# (0 if the variable has no lead). Variables with col-3 > 0 are forward-looking.
# Their 1-based declaration-order indices form i_fwrd_b.
# This is read once and cached — no context field needed.

const _FWRD_LOADED = Ref(false)
const _I_FWRD_B    = Int[]

function _load_fwrd_indices!()
    _FWRD_LOADED[] && return
    isfile(_MODFILE_PATH) || error("modfile.json not found: $_MODFILE_PATH")

    raw = read(_MODFILE_PATH, String)

    # Locate "lead_lag_incidence" then scan forward for all [a, b, c] rows.
    # Format in file: "lead_lag_incidence": [[ 1, 79, 0],[ 2, 80, 570], ...]
    # Each row is exactly a 3-integer bracket: \[ \d+, \d+, \d+\]
    idx = findfirst("lead_lag_incidence", raw)
    idx === nothing && error("lead_lag_incidence not found in modfile.json")

    # Read up to 100 KB from that point — enough for 491 rows × ~15 chars each
    chunk_end = min(lastindex(raw), idx[end] + 100_000)
    chunk = raw[idx[end]+1 : chunk_end]

    # Each row: [ <lag>, <cur>, <lead>] — parse lead (3rd number)
    rows = collect(eachmatch(r"\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]", chunk))
    isempty(rows) && error("Could not parse lead_lag_incidence rows from modfile.json")

    # Forward-looking: declaration-order index j where lead column > 0
    empty!(_I_FWRD_B)
    for (j, row) in enumerate(rows)
        lead_col = parse(Int, row.captures[3])
        lead_col > 0 && push!(_I_FWRD_B, j)
    end
    length(rows) == 491 || @printf "  [SMM] Warning: parsed %d rows, expected 491\n" length(rows)
    _FWRD_LOADED[] = true
    @printf "  [SMM] Forward indices loaded from modfile.json: %d of %d vars are forward-looking\n" length(_I_FWRD_B) length(rows)
end


# =========================================================================== #
#  PARAMETER ACCESS                                                            #
# =========================================================================== #

const _param_cache = IdDict{Any, Dict{String,Int}}()

function param_idx(context, name::String)
    haskey(_param_cache, context) || (_param_cache[context] = _build_param_cache(context))
    return get(_param_cache[context], name, nothing)
end

"""
Build parameter name → index mapping.

Primary: read from modfile.json 'parameters' array.  Dynare writes parameters
in declaration order, which matches context.models[1].params vector ordering.
Position j in the array (1-based) = params[j].  This is robust and works even
when the symboltable is corrupted (as it is in contexts built with stoch_simul
that failed due to the ARM Mac gees issue).

Fallback: iterate symboltable fields (works on uncorrupted contexts).
"""
function _build_param_cache(context)
    # ---- Primary: modfile.json parameters array (always available) ---------- #
    # The "parameters" array is ordered identically to context.models[1].params.
    # Use position-based slicing (not a dotall regex) so it works whether the
    # JSON is compact or pretty-printed with newlines.
    if isfile(_MODFILE_PATH)
        try
            raw  = read(_MODFILE_PATH, String)
            p_start = findfirst("\"parameters\"", raw)
            p_end   = findfirst("\"orig_endo_nbr\"", raw)
            if p_start !== nothing && p_end !== nothing && p_start[1] < p_end[1]
                chunk = raw[p_start[1] : p_end[1]-1]
                # Match only exact "name" key (not "texName"/"longName")
                names = collect(eachmatch(r"\"name\"\s*:\s*\"([^\"]+)\"", chunk))
                if length(names) > 100   # sanity: model has hundreds of params
                    cache = Dict{String,Int}()
                    for (j, m) in enumerate(names)
                        cache[String(m.captures[1])] = j
                    end
                    return cache
                end
            end
        catch; end
    end

    # ---- Fallback: symboltable traversal (works on non-corrupted contexts) -- #
    cache = Dict{String,Int}()
    st = context.symboltable
    try
        for fname in fieldnames(typeof(st))
            obj = getfield(st, fname)
            obj isa AbstractDict || continue
            for (k, v) in obj
                try
                    if hasproperty(v, :name) && hasproperty(v, :index)
                        cache[String(v.name)] = v.index
                    end
                catch; end
            end
        end
    catch; end
    return cache
end

set_param!(ctx, name::String, val::Real) = begin
    idx = param_idx(ctx, name)
    idx === nothing && return
    ctx.models[1].params[idx] = Float64(val)
end

get_param_val(ctx, name::String) = begin
    idx = param_idx(ctx, name)
    idx === nothing ? NaN : ctx.models[1].params[idx]
end


# =========================================================================== #
#  KLEIN (2000) FIRST-ORDER SOLVER  (pure Julia, no LAPACK gees)             #
# =========================================================================== #

"""
Compute the 491×642 dynamic Jacobian using the compiled Dynare model files.
The SparseDynamic*.jl files are included at top level (not inside a function)
to avoid Julia world-age issues.  Returns sparse G = [A|B|C|D]:
  A (491×78):  ∂f/∂y_{bk,t-1}   B (491×491): ∂f/∂y_t
  C (491×56):  ∂f/∂y_{fw,t+1}   D (491×17):  ∂f/∂ε_t
"""
function _eval_dynamic_jacobian(context)
    _load_jacobian_structure!()
    _DYNARE_MODEL_LOADED[] || error("Compiled model files not loaded. Run main_SOE_gap.jl first.")

    ss     = context.results.model_results[1].trends.endogenous_steady_state
    n_endo = length(ss)
    params = context.models[1].params
    n_exo  = context.models[1].exo_nbr

    # y = [ss_lag; ss_now; ss_lead] — all 491 vars at each time period
    y = vcat(ss, ss, ss)          # length 3*491 = 1473
    x = zeros(n_exo)              # exo at SS = 0 (shocks)

    T    = Float64[]              # empty temp (length 0 asserted)
    g1_v = zeros(_N_NZ[])         # 3350 Jacobian values

    SparseDynamicG1TT!(T, y, x, params, ss)
    SparseDynamicG1!(T, g1_v, y, x, params, ss)

    return sparse(_JAC_ROWS, _JAC_COLS, g1_v, _N_EQ[], _N_COL[])
end

"""
    resolve_first_order!(context) -> (success, g1_1, g1_2, Sigma_e)

Compute the first-order decision rule matrices:
  g1_1 (491×78):  response of all endogenous to 78 state variables
  g1_2 (491×17):  response of all endogenous to 17 shocks

Uses the compiled Dynare Jacobian + GenericSchur.jl (Klein 2000).
Works on all platforms — no LAPACK dependency.
"""
function resolve_first_order!(context)
    try
        # Fast path: try Dynare.jl's own solver (works on some platforms)
        for fn in [:compute_first_order_solution!, :first_order_solution!]
            isdefined(Dynare, fn) || continue
            try
                getfield(Dynare, fn)(context)
                mr  = context.results.model_results[1]
                lre = mr.linearrationalexpectations
                return true,
                       Matrix{Float64}(lre.g1_1),
                       Matrix{Float64}(lre.g1_2),
                       context.models[1].Sigma_e
            catch; end
        end

        # Klein (2000) pure-Julia path
        ok, g1_1, g1_2, Σe = _klein_solve(context)
        ok && return true, g1_1, g1_2, Σe

        # Fallback: use the existing decision rule already in the context.
        # This is exact for θ = θ_baseline, approximate for other θ.
        # Allows estimation of shock parameters (which enter via Σe, not g1).
        mr  = context.results.model_results[1]
        lre = mr.linearrationalexpectations
        if !isempty(lre.g1_1) && size(lre.g1_1, 1) >= 400
            return true,
                   Matrix{Float64}(lre.g1_1),
                   Matrix{Float64}(lre.g1_2),
                   context.models[1].Sigma_e
        end

        return false, zeros(0,0), zeros(0,0), zeros(0,0)
    catch
        return false, zeros(0,0), zeros(0,0), zeros(0,0)
    end
end

function _klein_solve(context)
    try
        # ---- 1. Assemble dynamic Jacobian -------------------------------- #
        G = _eval_dynamic_jacobian(context)
        n_eq = _N_EQ[]   # 491

        # Column partition (confirmed from dynamic.json):
        n_bk  = 78;  n_endo = 491;  n_fw = 56;  n_exo = 17
        # A: cols 1:78       B: cols 79:569    C: cols 570:625   D: cols 626:642
        A = Matrix{Float64}(G[:, 1:n_bk])
        B = Matrix{Float64}(G[:, n_bk+1 : n_bk+n_endo])
        C = Matrix{Float64}(G[:, n_bk+n_endo+1 : n_bk+n_endo+n_fw])
        D = Matrix{Float64}(G[:, n_bk+n_endo+n_fw+1 : end])

        # ---- 2. Identify backward and forward variable index sets -------- #
        bkwrd_b = collect(Int, context.models[1].i_bkwrd_b)  # 78 state rows

        # i_fwrd_b: 1-based declaration-order indices of forward-looking variables.
        # Priority 1: context field (various possible names across Dynare.jl versions).
        # Priority 2: modfile.json lead_lag_incidence (deterministic, always available).
        # Priority 3: dynare_fwd_rows.csv (written by updated subprocess).
        fwrd_b = nothing
        m1 = context.models[1]
        for fname in (:i_fwrd_b, :i_fwrd, :i_lead_b, :i_nontemporal_b)
            if isdefined(m1, fname)
                v = getfield(m1, fname)
                if !isempty(v) && length(v) == n_fw
                    fwrd_b = collect(Int, v); break
                end
            end
        end
        if isnothing(fwrd_b)
            # Load from modfile.json (lead_lag_incidence, col 3 > 0 → forward var)
            _load_fwrd_indices!()
            if length(_I_FWRD_B) == n_fw
                fwrd_b = copy(_I_FWRD_B)
            elseif !isempty(_I_FWRD_B)
                @printf "  [Klein] modfile.json gave %d fwrd vars, expected %d\n" length(_I_FWRD_B) n_fw
            end
        end
        if isnothing(fwrd_b)
            fwd_csv = joinpath(_JDYN_DIR, "mod", "dynare_fwd_rows.csv")
            if isfile(fwd_csv)
                df_fwd = CSV.read(fwd_csv, DataFrame)
                fb2 = Int.(df_fwd.fwd_row)
                length(fb2) == n_fw && (fwrd_b = fb2)
            end
        end
        isnothing(fwrd_b) && return false, zeros(0,0), zeros(0,0), zeros(0,0)

        # ---- 3. Static variable elimination ------------------------------ #
        # The full 491×491 QZ has only 56 finite eigenvalues (rank(C_pad)=56),
        # so n_stable can never reach n_bk=78.  Reduce to the 134×134 system
        # of backward+forward variables by eliminating the 357 static variables.
        #
        # Static variables: not in bkwrd_b and not in fwrd_b
        s_idx = setdiff(1:n_endo, union(bkwrd_b, fwrd_b))  # 357 static var indices
        n_s   = length(s_idx)   # 357

        # Static equations: rows where ∂f/∂y_{t-1} ≈ 0  AND  ∂f/∂y_{t+1} ≈ 0
        # (purely contemporaneous — no leads or lags)
        s_eq = findall(r -> norm(A[r,:]) < 1e-10 && norm(C[r,:]) < 1e-10, 1:n_eq)
        dyn_eq = setdiff(1:n_eq, s_eq)   # 134 dynamic equations (have A or C ≠ 0)

        if length(s_eq) != n_s
            @printf "  [Klein] Static elimination mismatch: %d static eqs, %d static vars (expected both=%d)\n" length(s_eq) n_s n_s
            return false, zeros(0,0), zeros(0,0), zeros(0,0)
        end

        # Express statics as y_s = -B_ss^{-1}(B_sb y_b + B_sf y_f)
        B_ss = B[s_eq, s_idx]          # 357×357
        B_sb = B[s_eq, bkwrd_b]        # 357×78
        B_sf = B[s_eq, fwrd_b]         # 357×56
        F_ss = lu(B_ss)
        Bss_inv_sb = F_ss \ B_sb        # 357×78
        Bss_inv_sf = F_ss \ B_sf        # 357×56

        # Substitute into dynamic equations: build reduced 134×134 system
        B_db = B[dyn_eq, bkwrd_b]      # 134×78
        B_df = B[dyn_eq, fwrd_b]       # 134×56
        B_ds = B[dyn_eq, s_idx]        # 134×357
        A_dyn = A[dyn_eq, :]           # 134×78
        C_dyn = C[dyn_eq, :]           # 134×56

        B_red_b = B_db - B_ds * Bss_inv_sb   # 134×78  (after static elimination)
        B_red_f = B_df - B_ds * Bss_inv_sf   # 134×56
        B_red   = [B_red_b B_red_f]           # 134×134
        C_red   = C_dyn                        # 134×56  (A stays as A_dyn)

        # ---- 4. Ordered QZ on the reduced 134×134 system ----------------- #
        # Variable ordering in 134-space: positions 1:78 = backward, 79:134 = forward
        # Leading matrix AA: C_red in forward columns (79:134), zeros in backward
        n_134   = n_bk + n_fw   # 134
        fw_cols = (n_bk+1):n_134   # 79:134

        C_pad_134 = zeros(n_134, n_134)
        C_pad_134[:, fw_cols] = C_red   # only forward vars appear at t+1

        AA_qz = C_pad_134   # 134×134 leading  (coefficient of y_{t+1})
        BB_qz = -B_red       # 134×134 lagging  (coefficient of y_t)

        # Use Julia's built-in LinearAlgebra.schur (LAPACK dgges, no select callback).
        # On ARM Mac, Apple Accelerate supports dgges without a select function,
        # unlike gees (which has a Fortran callback and fails on Apple Silicon).
        F = LinearAlgebra.schur(AA_qz, BB_qz)
        # LinearAlgebra.GeneralizedSchur: F.S from AA, F.T from BB, eigenvalues = F.S/F.T
        S_diag = abs.(diag(F.S))
        T_diag = abs.(diag(F.T))

        # Eigenvalues λ = S/T (S from AA=C_pad_134, T from BB=-B_red).
        # Backward positions (cols 1:78 of C_pad_134 = 0): S≈0 → λ≈0 (stable) ✓
        # Forward positions (cols 79:134 = C_red): λ > 1 (explosive for DSGE) ✓
        λ = [T_diag[i] < 1e-14 ? Inf : S_diag[i] / T_diag[i] for i in 1:n_134]

        n_stable = sum(λ .< 1.0)
        if n_stable != n_bk
            @printf "  [Klein] BK failed on 134×134 reduced system: %d stable eigenvalues, expected %d\n" n_stable n_bk
            return false, zeros(0,0), zeros(0,0), zeros(0,0)
        end

        # Reorder: stable eigenvalues (backward vars) first using LinearAlgebra.ordschur!
        select = BitVector(λ .< 1.0)
        F2 = LinearAlgebra.ordschur!(F, select)
        Z  = real(F2.Z)   # 134×134

        # ---- 5. Extract decision rule from the reduced system ------------ #
        # In the 134-dim space: rows 1:78 = backward vars, rows 79:134 = forward vars
        bk_in_134 = 1:n_bk
        nf_in_134 = (n_bk+1):n_134   # forward vars in the 134 space

        Z11 = Z[bk_in_134, 1:n_bk]   # 78×78  (state block)
        Z21 = Z[nf_in_134, 1:n_bk]   # 56×78  (non-state block → forward vars)

        abs(det(Z11)) < 1e-10 && return false, zeros(0,0), zeros(0,0), zeros(0,0)

        # Decision rule in the 134 space: y_fw = gx_fw × y_bk_{t-1}
        gx_fw = Z21 / Z11   # 56×78  (forward vars as function of states)

        # Recover static vars: y_s = -(Bss_inv_sb + Bss_inv_sf * gx_fw) y_bk_{t-1}
        gx_s = -(Bss_inv_sb + Bss_inv_sf * gx_fw)   # 357×78

        # ---- 6. Build full g1_1 (491×78) --------------------------------- #
        g1_1 = zeros(n_endo, n_bk)
        g1_1[bkwrd_b, :] = I(n_bk)    # backward vars: identity (predetermined)
        g1_1[fwrd_b,  :] = gx_fw       # forward vars
        g1_1[s_idx,   :] = gx_s        # static vars

        # ---- 7. Solve for shock impact g1_2 (491×17) --------------------- #
        # From B*g1_2 + C*(g1_1[fwrd,:]*g1_2[bkwrd,:]) + D = 0  (approx):
        effective_B = B + C * g1_1[fwrd_b, :] * g1_1[bkwrd_b, :]
        g1_2 = if abs(det(effective_B)) > 1e-10
            -effective_B \ D
        else
            -pinv(effective_B) * D
        end

        # Store computed decision rule back into the context for the cached path
        try
            lre = context.results.model_results[1].linearrationalexpectations
            if size(lre.g1_1) == size(g1_1)
                lre.g1_1 .= g1_1
                lre.g1_2 .= g1_2
            end
        catch; end

        Σe = context.models[1].Sigma_e
        return true, g1_1, g1_2, Σe

    catch e
        # FieldError means i_fwrd_b / i_bkwrd_b field missing from context —
        # expected for contexts built without stoch_simul.  Suppress the print
        # to avoid flooding output during CMA-ES (115k evaluations).
        # For any other exception type, print once for diagnostics.
        if !(e isa FieldError)
            @printf "  [Klein] Exception type: %s\n" typeof(e)
        end
        return false, zeros(0,0), zeros(0,0), zeros(0,0)
    end
end


# =========================================================================== #
#  STEADY STATE RECOMPUTATION  (when epsY or epsM change)                     #
# =========================================================================== #

function recompute_ss!(context, epsY, epsM, baseline, endo_names)
    nsec = baseline.nsec
    modepsY = fill(epsY, nsec); modepsM = fill(epsM, nsec)
    ss_vec = context.results.model_results[1].trends.endogenous_steady_state
    get_ss(nm) = let idx = findfirst(==(nm), endo_names)
                     idx === nothing ? 1.0 : ss_vec[idx] end

    x0 = [[get_ss("PH_$(i)") for i in 1:nsec]; get_ss("w"); get_ss("Q"); get_ss("C")]
    res = nlsolve(
        (F, x) -> F .= steady_ntwsoe(x, baseline.PVstar_ss, baseline.epsilon_val,
                                      baseline.modvarrho, baseline.sigmaH_val,
                                      baseline.modgammag, baseline.modgammas,
                                      baseline.ombar_val, 1-baseline.ombar_val,
                                      baseline.modchiX, baseline.omegaX_val,
                                      baseline.etastar_val, baseline.ystar_ss_val,
                                      baseline.modalpha, baseline.modalphaV,
                                      baseline.modbeta, modepsY, modepsM,
                                      baseline.gamma_val, baseline.chi_val,
                                      baseline.psi_val, ones(nsec), baseline.tb_target),
        x0; ftol=1e-12, method=:trust_region, show_trace=false)
    !converged(res) && return false

    pH_ss = res.zero[1:nsec]; w_ss = res.zero[nsec+1]
    Q_ss = res.zero[nsec+2];  C_ss = res.zero[nsec+3]
    PL_ss = fill(w_ss, nsec); PV_ss = Q_ss * baseline.PVstar_ss
    MCi_ss = (baseline.epsilon_val-1)/baseline.epsilon_val .* pH_ss
    PMi_ss = (baseline.modbeta * (pH_ss .^ (1 .- modepsM))) .^ (1 ./ (1 .- modepsM))
    P_ss = (baseline.modvarrho .^ baseline.sigmaH_val .* pH_ss .^ (1-baseline.sigmaH_val)
           .+ (1 .- baseline.modvarrho) .^ baseline.sigmaH_val .* PV_ss .^ (1-baseline.sigmaH_val)) .^ (1/(1-baseline.sigmaH_val))
    p_g = prod(P_ss .^ baseline.modgammag); p_s = prod(P_ss .^ baseline.modgammas)
    C_g = baseline.ombar_val*C_ss/p_g; C_s = (1-baseline.ombar_val)*C_ss/p_s
    C_gi = baseline.modgammag .* (p_g./P_ss) .* C_g
    C_si = baseline.modgammas .* (p_s./P_ss) .* C_s
    CHi = baseline.modvarrho .^ baseline.sigmaH_val .* (pH_ss./P_ss) .^ (-baseline.sigmaH_val) .* (C_gi.+C_si)
    CFi = (1 .- baseline.modvarrho) .^ baseline.sigmaH_val .* (PV_ss./P_ss) .^ (-baseline.sigmaH_val) .* (C_gi.+C_si)
    PX = prod(pH_ss .^ baseline.modchiX)
    X  = baseline.omegaX_val * (PX/Q_ss)^(-baseline.etastar_val) * baseline.ystar_ss_val
    Xi = baseline.modchiX .* X .* PX ./ pH_ss

    inner = nlsolve(
        (F, x) -> steady_ntwsoe_system!(F, x, baseline.modalpha, baseline.modalphaV,
                                         baseline.modbeta, MCi_ss, PMi_ss, PL_ss, PV_ss,
                                         CHi, Xi, modepsY, modepsM, ones(nsec), pH_ss),
        vcat(max.(MCi_ss./PMi_ss,1e-20), max.(MCi_ss./PL_ss,1e-20),
             max.(MCi_ss./PV_ss,1e-20), fill(0.2,nsec));
        ftol=1e-10, method=:trust_region, show_trace=false)
    !converged(inner) && return false

    Yi_ss = inner.zero[3*nsec+1:4*nsec]; M_ss = inner.zero[1:nsec]
    IMP = sum(inner.zero[2*nsec+1:3*nsec]) + sum(CFi)
    TB  = PX*X - PV_ss*IMP; GDP = C_ss + TB
    Pistar = 1.0; r_star = Pistar/baseline.beta_val
    Bstar  = -TB / (Q_ss*(1 - r_star/Pistar))

    set_param!(context, "bbar", Q_ss*Bstar/GDP)
    set_param!(context, "Y_ss", sum(Yi_ss))
    set_param!(context, "M_tot_ss", sum(M_ss))
    for i in 1:nsec
        set_param!(context, "PL_ss$(i)", w_ss)
        set_param!(context, "epsY_$(i)", epsY)
        set_param!(context, "epsM_$(i)", epsM)
    end
    ss_mut = context.results.model_results[1].trends.endogenous_steady_state
    upd(nm,val) = let idx=findfirst(==(nm),endo_names); idx!==nothing&&(ss_mut[idx]=val); end
    upd("w",w_ss); upd("Q",Q_ss); upd("C",C_ss); upd("GDP",GDP); upd("TB",TB)
    for i in 1:nsec; upd("PH_$(i)",pH_ss[i]); upd("Y_$(i)",Yi_ss[i]); upd("L_$(i)",inner.zero[nsec+i]); end
    return true
end


# =========================================================================== #
#  MOMENT NAMES                                                                #
# =========================================================================== #

const MOMENT_NAMES = vcat(
    ["std(Y_$(i))"  for i in 1:12], ["std(PH_$(i))" for i in 1:12], ["std(L_$(i))"  for i in 1:12],
    ["std(GDP)", "std(pi)", "corr(GDP,pi)", "mean goods expenditure share",
     "std(Q)", "autocorr(Q)", "corr(GDP,Q)",
     "rank corr: output (model vs data)", "rank corr: prices (model vs data)", "rank corr: labor  (model vs data)"])


# =========================================================================== #
#  MAIN MOMENT FUNCTION                                                        #
# =========================================================================== #

function smm_model_moments(θ, context, baseline, endo_names)
    nsec = baseline.nsec; NAN46 = fill(NaN, 46)
    ilabcosts=θ[1]; epsY=θ[2]; epsM=θ[3]; kappaV=exp(θ[4])
    rho_om=θ[5]; sigma_om=θ[6]; rho_A=θ[7]; isigma_tfp=θ[8:19]
    rho_pvstar=θ[20]; sigma_pvstar=θ[21]; rho_xi=θ[22]; sigma_xi=θ[23]

    (!(0<epsY<5)||!(0<epsM<2)||ilabcosts<=0||kappaV<=0||abs(rho_om)>=1||
     sigma_om<0||abs(rho_A)>=1||any(isigma_tfp.<0)||abs(rho_pvstar)>=1||
     sigma_pvstar<0||abs(rho_xi)>=1||sigma_xi<0) && return NAN46, false

    set_param!(context,"ilabcosts",ilabcosts); set_param!(context,"kappaV",kappaV)
    set_param!(context,"rho_om1",rho_om);      set_param!(context,"sigma_om",sigma_om)
    set_param!(context,"rho_tfp1",rho_A);      set_param!(context,"rho_pvstar",rho_pvstar)
    set_param!(context,"sigma_pvstar",sigma_pvstar); set_param!(context,"rho_xi",rho_xi)
    set_param!(context,"sigma_xi",sigma_xi)
    for i in 1:nsec; set_param!(context,"isigma_tfp_$(i)",isigma_tfp[i]); end

    epsY_prev = get_param_val(context,"epsY_1"); epsM_prev = get_param_val(context,"epsM_1")
    need_ss = abs(epsY-epsY_prev)>1e-8 || abs(epsM-epsM_prev)>1e-8
    for i in 1:nsec; set_param!(context,"epsY_$(i)",epsY); set_param!(context,"epsM_$(i)",epsM); end
    if need_ss; ok=recompute_ss!(context,epsY,epsM,baseline,endo_names); !ok&&return NAN46,false; end

    success, T, R, Σe = resolve_first_order!(context)
    !success && return NAN46, false

    # Build normalized Sigma_e from the 16 active shocks.
    # The SMM context has Sigma_e = zeros(17,17) because stoch_simul was not
    # run for the _smm mod file.  We CANNOT rely on mutating context.models[1].Sigma_e
    # (may fail silently on ARM due to immutability in deserialized structs).
    # Instead, construct Sigma_e locally: each active shock has unit variance;
    # amplitude scaling is handled by the parameters in the model equations
    # (sigma_om * eps_om, isigma_tfp_i * epsA_i, sigma_pvstar * eps_pvstar,
    #  sigma_xi * eps_xi).
    #
    # varexo order (from NK_SOE_lev_gap2_smm.mod):
    #   1 = eps_om    2 = eps_i     3 = epschi (unused)  4 = eps_pvstar
    #   5:16 = epsA_1:12             17 = eps_xi
    n_exo = size(R, 2)   # should be 17
    Σe = zeros(n_exo, n_exo)
    Σe[1, 1]  = 1.0                        # eps_om
    Σe[4, 4]  = 1.0                        # eps_pvstar
    for i in 5:min(16, n_exo)
        Σe[i, i] = 1.0                     # epsA_1 through epsA_12
    end
    if n_exo >= 17; Σe[17, 17] = 1.0; end # eps_xi

    sr = context.models[1].i_bkwrd_b
    P  = local_dlyap(T[sr,:], R[sr,:]*Σe*R[sr,:]')
    (any(diag(P).<-1e-10)||any(isnan.(P))) && return NAN46, false
    P  = (P+P')/2
    Γ  = T*P*T' + R*Σe*R'; Γ=(Γ+Γ')/2
    Γ1 = T*T[sr,:]*( P*T' + R[sr,:]*Σe*R')

    ys = context.results.model_results[1].trends.endogenous_steady_state
    ei = Dict(nm=>i for (i,nm) in enumerate(endo_names))
    pstd(vn) = let idx=get(ei,vn,nothing); idx===nothing ? 0.0 :
                   sqrt(max(Γ[idx,idx],0.0))/max(abs(ys[idx]),1e-12); end
    xcorr(v1,v2) = let i1=get(ei,v1,nothing),i2=get(ei,v2,nothing)
                       (i1===nothing||i2===nothing) ? NaN :
                       let d=sqrt(max(Γ[i1,i1],0.0)*max(Γ[i2,i2],0.0))
                           d<1e-15 ? 0.0 : clamp(Γ[i1,i2]/d,-1.0,1.0) end; end
    i_Q = get(ei,"Q",nothing)
    acQ = (i_Q!==nothing&&Γ[i_Q,i_Q]>1e-15) ? Γ1[i_Q,i_Q]/Γ[i_Q,i_Q] : NaN

    std_Y=[pstd("Y_$(i)") for i in 1:nsec]; std_PH=[pstd("PH_$(i)") for i in 1:nsec]
    std_L=[pstd("L_$(i)") for i in 1:nsec]
    dY=baseline.data_std_Y; dPH=baseline.data_std_PH; dL=baseline.data_std_L
    vy=isfinite.(std_Y).&isfinite.(dY); vp=isfinite.(std_PH).&isfinite.(dPH); vl=isfinite.(std_L).&isfinite.(dL)
    rY=sum(vy)>=3 ? safe_spearman(std_Y[vy],dY[vy]) : 0.0
    rP=sum(vp)>=3 ? safe_spearman(std_PH[vp],dPH[vp]) : 0.0
    rL=sum(vl)>=3 ? safe_spearman(std_L[vl],dL[vl]) : 0.0

    return [std_Y;std_PH;std_L;pstd("GDP");pstd("pi");xcorr("GDP","pi");
            baseline.ombar_val;pstd("Q");acQ;xcorr("GDP","Q");rY;rP;rL], true
end
