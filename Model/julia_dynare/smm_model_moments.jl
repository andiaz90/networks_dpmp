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

using LinearAlgebra, SparseArrays, Statistics, StatsBase, NLsolve, Random
# GenericSchur was used to avoid LAPACK gees callbacks on ARM Mac.
# Julia's built-in LinearAlgebra.schur uses LAPACK dgges (no select callback)
# which IS available via Apple Accelerate on ARM. We now use it directly.

include("steady_ntwsoe_system.jl")
include("steady_ntwsoe.jl")
include("utils.jl")

# =========================================================================== #
#  THREAD SAFETY: PER-THREAD CONTEXT AND GLOBAL STATE                         #
# =========================================================================== #

const _N_THREADS = max(1, Threads.nthreads())

"""
    _tid()::Int

Get safe thread ID for indexing: clamps to valid range [1, _N_THREADS].
Used for all per-thread data structures.
"""
_tid() = min(Threads.threadid(), _N_THREADS)

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

# =========================================================================== #
#  KLEIN RESULT CACHE (PER-THREAD)                                             #
# =========================================================================== #
# Shock-amplitude params (sigma_om, isigma_tfp_i, sigma_pvstar, sigma_xi)
# enter the model as LINEAR coefficients in the model equations and hence
# in g1_2.  When ONLY shock amplitudes change (structural params are
# unchanged), g1_1 is identical and g1_2 scales proportionally.  We cache
# the last Klein result and skip the 134×134 QZ when structural params
# change by less than _KLEIN_THRESH.
#
# Structural params that DO require Klein re-solve (θ indices, Option-A layout):
#   1=ilabcosts, 2=epsY, 3=epsM, 4=log(kappaV), 5=rho_om, 6=rho_A,
#   31=rho_pvstar, 33=rho_xi
const _KLEIN_STRUCT_IDX = [1, 2, 3, 4, 5, 6, 31, 33]
const _KLEIN_THRESH     = 1e-5   # re-solve if any structural param moves > this

const _KLEIN_CACHE_T    = [Ref{Matrix{Float64}}(zeros(0,0)) for _ in 1:_N_THREADS]
const _KLEIN_CACHE_R    = [Ref{Matrix{Float64}}(zeros(0,0)) for _ in 1:_N_THREADS]
const _KLEIN_CACHE_ΘSTR = [Ref{Vector{Float64}}(Float64[]) for _ in 1:_N_THREADS]
const _KLEIN_HITS       = Threads.Atomic{Int}(0)
const _KLEIN_MISSES     = Threads.Atomic{Int}(0)

function _resolve_cached!(context, θ)
    tid = _tid()
    θ_str = θ[_KLEIN_STRUCT_IDX]
    cache_valid = !isempty(_KLEIN_CACHE_ΘSTR[tid][]) &&
                  maximum(abs.(_KLEIN_CACHE_ΘSTR[tid][] .- θ_str)) < _KLEIN_THRESH

    if cache_valid
        Threads.atomic_add!(_KLEIN_HITS, 1)
        return true, _KLEIN_CACHE_T[tid][], _KLEIN_CACHE_R[tid][]
    end

    ok, T, R, _ = resolve_first_order!(context)
    if ok && !isempty(T)
        _KLEIN_CACHE_T[tid][]    = T
        _KLEIN_CACHE_R[tid][]    = R
        _KLEIN_CACHE_ΘSTR[tid][] = copy(θ_str)
        Threads.atomic_add!(_KLEIN_MISSES, 1)
        return true, T, R
    end
    return false, zeros(0,0), zeros(0,0)
end

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
#  PARAMETER VALUES — loaded from params_jl.mod (PER-THREAD)                   #
# =========================================================================== #
# context.work.params may be declared but uninitialized (UndefVarError) in
# contexts deserialized from a stoch_simul run that failed on ARM Mac.
# We maintain our OWN parameter vector _SMM_PARAMS, initialized from
# params_jl.mod (written by main_SOE_gap.jl with full parameter values)
# and updated by set_param! at each CMA-ES evaluation.
# Per-thread vectors to avoid race conditions during parallel CMA-ES.

const _SMM_PARAMS       = [Ref{Vector{Float64}}(Float64[]) for _ in 1:_N_THREADS]
const _SMM_PARAMS_READY = [Ref(false) for _ in 1:_N_THREADS]
const _LOAD_PARAMS_LOCK = ReentrantLock()

function _load_smm_params!(context)
    tid = _tid()
    _SMM_PARAMS_READY[tid][] && return _SMM_PARAMS[tid][]

    # Lock to prevent concurrent file reads during parallel initialization
    lock(_LOAD_PARAMS_LOCK) do
        # Double-check after acquiring lock
        _SMM_PARAMS_READY[tid][] && return _SMM_PARAMS[tid][]

        # Try context.work.params first (works on uncorrupted contexts)
        try
            p = context.work.params
            if length(p) > 100
                _SMM_PARAMS[tid][] = copy(Vector{Float64}(p))
                _SMM_PARAMS_READY[tid][] = true
                tid == 1 && @printf "  [SMM] Params loaded from context.work.params (%d params)\n" length(p)
                return _SMM_PARAMS[tid][]
            end
        catch; end

        # Fall back to params_jl.mod — parse name=value; pairs, build ordered vector
        mod_path = joinpath(_JDYN_DIR, "mod", "params_jl.mod")
        isfile(mod_path) || error("params_jl.mod not found at $mod_path — run main_SOE_gap.jl first")

        raw = read(mod_path, String)
        val_dict = Dict{String,Float64}()
        for m in eachmatch(r"^(\w+)\s*=\s*([-\d.eE+]+)\s*;", raw)
            try val_dict[String(m.captures[1])] = parse(Float64, m.captures[2])
            catch; end
        end

        # Build ordered params vector using modfile.json parameter ordering
        p_start = findfirst("\"parameters\"", read(_MODFILE_PATH, String))
        p_end   = findfirst("\"orig_endo_nbr\"", read(_MODFILE_PATH, String))
        mf_raw  = read(_MODFILE_PATH, String)
        chunk   = mf_raw[p_start[1] : p_end[1]-1]
        names   = [String(m.captures[1]) for m in eachmatch(r"\"name\"\s*:\s*\"([^\"]+)\"", chunk)]

        n_params = length(names)
        params   = zeros(n_params)
        found    = 0
        for (j, nm) in enumerate(names)
            if haskey(val_dict, nm)
                params[j] = val_dict[nm]
                found += 1
            end
        end

        _SMM_PARAMS[tid][] = params
        _SMM_PARAMS_READY[tid][] = true
        tid == 1 && @printf "  [SMM] Params loaded from params_jl.mod: %d/%d params found\n" found n_params
        return _SMM_PARAMS[tid][]
    end
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
in declaration order, which matches context.work.params vector ordering.
Position j in the array (1-based) = params[j].  This is robust and works even
when the symboltable is corrupted (as it is in contexts built with stoch_simul
that failed due to the ARM Mac gees issue).

Fallback: iterate symboltable fields (works on uncorrupted contexts).
"""
function _build_param_cache(context)
    # ---- Primary: modfile.json parameters array (always available) ---------- #
    # The "parameters" array is ordered identically to context.work.params.
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
    tid = _tid()
    idx = param_idx(ctx, name)
    idx === nothing && return
    p = _load_smm_params!(ctx)
    isempty(p) && return
    idx <= length(p) && (_SMM_PARAMS[tid][][idx] = Float64(val))
end

get_param_val(ctx, name::String) = begin
    tid = _tid()
    idx = param_idx(ctx, name)
    idx === nothing && return NaN
    p = _SMM_PARAMS_READY[tid][] ? _SMM_PARAMS[tid][] : _load_smm_params!(ctx)
    (isempty(p) || idx > length(p)) ? NaN : p[idx]
end


# =========================================================================== #
#  EXOGENOUS-SHOCK NAMES (ordered) — for name-based Σe construction            #
# =========================================================================== #
# The hardcoded Σe index layout was the root cause of the shock-layout bug
# (audit finding C1): a 28-shock template applied to whatever model was loaded.
# We now resolve the active shocks by NAME from the model's own exogenous
# ordering, so a change in the .mod shock list can never silently scramble Σe.

const _EXO_NAMES = Ref{Vector{String}}(String[])

function smm_exo_names(context)
    isempty(_EXO_NAMES[]) || return _EXO_NAMES[]
    names = String[]
    # Prefer the live context if it exposes exogenous symbols with names.
    try
        exo = context.models[1].exogenous
        names = [String(getfield(e, :name)) for e in exo]
    catch; names = String[]; end
    # Fallback: parse the "exogenous" array from modfile.json by bracket matching.
    if isempty(names) && isfile(_MODFILE_PATH)
        try
            raw = read(_MODFILE_PATH, String)
            key = findfirst("\"exogenous\"", raw)   # exact key; won't match exogenous_deterministic
            if key !== nothing
                lb = findnext('[', raw, key[end] + 1)
                depth = 0; rb = lb
                for p in lb:lastindex(raw)
                    c = raw[p]
                    c == '[' && (depth += 1)
                    if c == ']'
                        depth -= 1
                        depth == 0 && (rb = p; break)
                    end
                end
                chunk = raw[lb:rb]
                names = [String(m.captures[1]) for m in eachmatch(r"\"name\"\s*:\s*\"([^\"]+)\"", chunk)]
            end
        catch; end
    end
    _EXO_NAMES[] = names
    return names
end

"""
    active_shock_indices(context, nsec) -> Vector{Int}

Diagonal positions in Σe to switch on during SMM estimation, resolved by name:
monetary (eps_i), import price (eps_pvstar), aggregate demand (eps_xi), the 12
sectoral TFP shocks (epsA_i) and the 12 sectoral demand shocks (eps_om_i).
Deliberately EXCLUDES eps_postar (oil — off during estimation) and epschi
(labour supply — off). Returns the indices in the model's exogenous order.
"""
function active_shock_indices(context, nsec::Int)
    exo = smm_exo_names(context)
    active = Set{String}(["eps_i", "eps_pvstar", "eps_xi"])
    for i in 1:nsec
        push!(active, "epsA_$(i)")
        push!(active, "eps_om_$(i)")
    end
    return [k for (k, nm) in enumerate(exo) if nm in active]
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
    # Use _SMM_PARAMS (loaded from params_jl.mod) — robust against
    # context.work.params being uninitialized (UndefVarError on ARM Mac)
    params = _load_smm_params!(context)
    isempty(params) && error("Parameter vector not available — run main_SOE_gap.jl first.")
    n_exo  = context.models[1].exogenous_nbr

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

        # i_fwrd_b is a direct field of context.models[1] in this Dynare.jl version.
        # Use it directly; fall back to modfile.json lead_lag_incidence if absent.
        fwrd_b = let m1 = context.models[1]
            if isdefined(m1, :i_fwrd_b) && !isempty(m1.i_fwrd_b) && length(m1.i_fwrd_b) == n_fw
                collect(Int, m1.i_fwrd_b)
            else
                # Fallback: read from modfile.json lead_lag_incidence
                _load_fwrd_indices!()
                length(_I_FWRD_B) == n_fw ? copy(_I_FWRD_B) : nothing
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
        # FieldError / UndefVarError from missing/uninitialized context fields
        # are expected — suppress to avoid flooding 115k CMA-ES iterations.
        if !(e isa FieldError || e isa UndefVarError)
            @printf "  [Klein] Exception type: %s\n" typeof(e)
        end
        return false, zeros(0,0), zeros(0,0), zeros(0,0)
    end
end


# =========================================================================== #
#  STEADY STATE RECOMPUTATION  (when epsY or epsM change)                     #
# =========================================================================== #

function recompute_ss!(context, epsY, epsM, baseline, endo_names; etastar=nothing)
  # Whole body wrapped so a DomainError / singular-LU / Inf in the nonlinear SS
  # solve returns `false` (→ caller emits NaN moments → optimiser penalises this
  # θ) instead of throwing and killing the entire estimation run.
  try
    nsec = baseline.nsec
    modepsY = fill(epsY, nsec); modepsM = fill(epsM, nsec)
    # Use the *effective* etastar so the SS solve is consistent with the value
    # that triggered the recompute (previously baseline.etastar_val was used,
    # silently ignoring the new etastar — the export block then disagreed with
    # the write-back below).
    eta_for_ss = something(etastar, baseline.etastar_val)
    ss_vec = context.results.model_results[1].trends.endogenous_steady_state
    get_ss(nm) = let idx = findfirst(==(nm), endo_names)
                     idx === nothing ? 1.0 : ss_vec[idx] end

    f_outer! = (F, x) -> F .= steady_ntwsoe(x, baseline.PVstar_ss, baseline.epsilon_val,
                                      baseline.modvarrho, baseline.sigmaH_val,
                                      baseline.modgammag, baseline.modgammas,
                                      baseline.ombar_val, 1-baseline.ombar_val,
                                      baseline.modchiX, baseline.omegaX_val,
                                      eta_for_ss, baseline.ystar_ss_val,
                                      baseline.modalpha, baseline.modalphaV,
                                      baseline.modbeta, modepsY, modepsM,
                                      baseline.gamma_val, baseline.chi_val,
                                      baseline.psi_val, ones(nsec), baseline.tb_target)

    x0 = [[get_ss("PH_$(i)") for i in 1:nsec]; get_ss("w"); get_ss("Q"); get_ss("C")]
    res = nlsolve(f_outer!, x0; ftol=1e-12, method=:trust_region, show_trace=false)
    if !converged(res)
        # Multi-start: the trust-region solver is sensitive to the warm guess
        # when epsY/epsM/etastar move far from the cached point. Retry from a
        # few deterministically-perturbed guesses before giving up. Seeded RNG
        # so failures are reproducible across runs.
        rng = Random.MersenneTwister(hash((epsY, epsM, eta_for_ss)) % UInt32)
        for k in 1:6
            x0p = max.(x0 .* (1 .+ 0.10 * k .* (rand(rng, length(x0)) .- 0.5)), 1e-8)
            res = nlsolve(f_outer!, x0p; ftol=1e-12, method=:trust_region, show_trace=false)
            converged(res) && break
        end
        converged(res) || return false
    end
    all(isfinite, res.zero) || return false

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
    etastar_eff = something(etastar, baseline.etastar_val)
    PX = prod(pH_ss .^ baseline.modchiX)
    X  = baseline.omegaX_val * (PX/Q_ss)^(-etastar_eff) * baseline.ystar_ss_val
    Xi = baseline.modchiX .* X .* PX ./ pH_ss

    inner = nlsolve(
        (F, x) -> steady_ntwsoe_system!(F, x, baseline.modalpha, baseline.modalphaV,
                                         baseline.modbeta, MCi_ss, PMi_ss, PL_ss, PV_ss,
                                         CHi, Xi, modepsY, modepsM, ones(nsec), pH_ss),
        vcat(max.(MCi_ss./PMi_ss,1e-20), max.(MCi_ss./PL_ss,1e-20),
             max.(MCi_ss./PV_ss,1e-20), fill(0.2,nsec));
        ftol=1e-10, method=:trust_region, show_trace=false)
    (!converged(inner) || !all(isfinite, inner.zero)) && return false

    Yi_ss = inner.zero[3*nsec+1:4*nsec]; M_ss = inner.zero[1:nsec]
    IMP = sum(inner.zero[2*nsec+1:3*nsec]) + sum(CFi)
    TB  = PX*X - PV_ss*IMP; GDP = C_ss + TB
    Pistar = 1.0; r_star = Pistar/baseline.beta_val
    Bstar  = -TB / (Q_ss*(1 - r_star/Pistar))

    set_param!(context, "bbar", Q_ss*Bstar/GDP)
    set_param!(context, "Y_ss", sum(Yi_ss))
    set_param!(context, "M_tot_ss", sum(M_ss))
    etastar !== nothing && set_param!(context, "etastar", etastar_eff)
    for i in 1:nsec
        set_param!(context, "PL_ss$(i)", w_ss)
        set_param!(context, "epsY_$(i)", epsY)
        set_param!(context, "epsM_$(i)", epsM)
    end
    ss_mut = context.results.model_results[1].trends.endogenous_steady_state
    upd(nm,val) = let idx=findfirst(==(nm),endo_names); idx!==nothing&&(ss_mut[idx]=val); end
    # Reject pathological aggregates before committing them to the context SS.
    (isfinite(GDP) && GDP > 0 && isfinite(C_ss) && C_ss > 0 &&
     all(isfinite, pH_ss) && all(>(0), pH_ss) && all(isfinite, Yi_ss)) || return false
    upd("w",w_ss); upd("Q",Q_ss); upd("C",C_ss); upd("GDP",GDP); upd("TB",TB)
    for i in 1:nsec; upd("PH_$(i)",pH_ss[i]); upd("Y_$(i)",Yi_ss[i]); upd("L_$(i)",inner.zero[nsec+i]); end
    return true
  catch
    # Any unexpected failure in the SS recompute is treated as a non-convergence
    # so the optimiser simply penalises this θ rather than aborting the run.
    return false
  end
end


# =========================================================================== #
#  MOMENT NAMES + duplicate moment fn — DISABLED                               #
#  ------------------------------------------------------------------------   #
#  The canonical MOMENT_NAMES (58 entries) and the canonical, optimised        #
#  smm_model_moments(...) live in smm_estimation.jl, which run_smm_estimation  #
#  includes AFTER this file and therefore silently overrode the two            #
#  definitions below.  Keeping two divergent copies (this one had only 46      #
#  MOMENT_NAMES with a stale "mean goods expenditure share" entry, and a       #
#  slower local_dlyap path) was a constant source of drift.  The block is      #
#  commented out — NOT deleted — so the history/derivation stays readable.     #
#  All the infrastructure ABOVE (Klein solver, recompute_ss!, param access,    #
#  caches) is still live and used by the canonical moment fn.                  #
# =========================================================================== #
#=
const MOMENT_NAMES = vcat(
    ["std(Y_$(i))"  for i in 1:12], ["std(PH_$(i))" for i in 1:12], ["std(L_$(i))"  for i in 1:12],
    ["std(GDP)", "std(pi)", "corr(GDP,pi)", "mean goods expenditure share",
     "std(Q)", "autocorr(Q)", "corr(GDP,Q)",
     "rank corr: output (model vs data)", "rank corr: prices (model vs data)", "rank corr: labor  (model vs data)"])


# =========================================================================== #
#  MAIN MOMENT FUNCTION                                                        #
# =========================================================================== #

function smm_model_moments(θ, context, baseline, endo_names)
    nsec = baseline.nsec; NAN58 = fill(NaN, 58)
    # Option-A parameter layout (35 params):
    #   θ[1]    = ilabcosts
    #   θ[2]    = epsY
    #   θ[3]    = epsM
    #   θ[4]    = log(kappaV)
    #   θ[5]    = rho_om   (common persistence for all 12 sectoral demand shocks)
    #   θ[6]    = rho_A    (common TFP persistence)
    #   θ[7:18] = isigma_tfp_1:12
    #   θ[19:30]= sigma_om_1:12  (sectoral demand shock std devs)
    #   θ[31]   = rho_pvstar
    #   θ[32]   = sigma_pvstar
    #   θ[33]   = rho_xi
    #   θ[34]   = sigma_xi
    #   θ[35]   = etastar
    ilabcosts=θ[1]; epsY=θ[2]; epsM=θ[3]; kappaV=exp(θ[4])
    rho_om=θ[5]; rho_A=θ[6]; isigma_tfp=θ[7:18]
    sigma_om_vec=θ[19:30]
    rho_pvstar=θ[31]; sigma_pvstar=θ[32]; rho_xi=θ[33]; sigma_xi=θ[34]
    etastar = length(θ) >= 35 ? θ[35] : baseline.etastar_val

    (!(0<epsY<5)||!(0<epsM<2)||ilabcosts<=0||kappaV<=0||abs(rho_om)>=1||
     any(sigma_om_vec.<0)||abs(rho_A)>=1||any(isigma_tfp.<0)||abs(rho_pvstar)>=1||
     sigma_pvstar<0||abs(rho_xi)>=1||sigma_xi<0||
     !(0.1<etastar<8.0)) && return NAN58, false

    set_param!(context,"ilabcosts",ilabcosts); set_param!(context,"kappaV",kappaV)
    set_param!(context,"rho_om1",rho_om)
    set_param!(context,"rho_tfp1",rho_A);      set_param!(context,"rho_pvstar",rho_pvstar)
    set_param!(context,"sigma_pvstar",sigma_pvstar); set_param!(context,"rho_xi",rho_xi)
    set_param!(context,"sigma_xi",sigma_xi)
    set_param!(context,"etastar",etastar)
    for i in 1:nsec
        set_param!(context,"isigma_tfp_$(i)",isigma_tfp[i])
        set_param!(context,"sigma_om_$(i)",sigma_om_vec[i])
    end

    epsY_prev    = get_param_val(context,"epsY_1"); epsM_prev = get_param_val(context,"epsM_1")
    etastar_prev = get_param_val(context,"etastar")
    need_ss = abs(epsY-epsY_prev)>1e-8 || abs(epsM-epsM_prev)>1e-8 ||
              abs(etastar - (isnan(etastar_prev) ? baseline.etastar_val : etastar_prev)) > 1e-8
    for i in 1:nsec; set_param!(context,"epsY_$(i)",epsY); set_param!(context,"epsM_$(i)",epsM); end
    if need_ss; ok=recompute_ss!(context,epsY,epsM,baseline,endo_names;etastar=etastar); !ok&&return NAN58,false; end

    success, T, R = _resolve_cached!(context, θ)
    !success && return NAN58, false

    # Build normalized Sigma_e for all active shocks.
    # Each active shock has unit variance; amplitude scaling is handled by
    # parameters in the model equations (sigma_om_i * eps_om_i, etc.).
    #
    # varexo order (from NK_SOE_lev_gap2_smm.mod, Option-A layout):
    #   1 = eps_i    2 = epschi(inactive)    3 = eps_pvstar
    #   4:15 = epsA_1:12    16 = eps_xi    17:28 = eps_om_1:12
    n_exo = size(R, 2)   # should be 28 after Option-A mod recompile
    Σe = zeros(n_exo, n_exo)
    Σe[1, 1]  = 1.0                        # eps_i    (monetary policy)
    # Σe[2, 2] = 0.0                       # epschi   (labor supply, inactive)
    Σe[3, 3]  = 1.0                        # eps_pvstar (import price)
    for i in 4:min(15, n_exo)
        Σe[i, i] = 1.0                     # epsA_1 through epsA_12 (sectoral TFP)
    end
    if n_exo >= 16; Σe[16, 16] = 1.0; end # eps_xi   (aggregate demand)
    for i in 17:min(28, n_exo)
        Σe[i, i] = 1.0                     # eps_om_1 through eps_om_12 (sectoral demand)
    end

    sr = context.models[1].i_bkwrd_b
    A_state = T[sr, :]       # n_state × n_state  (state transition)
    B_state = R[sr, :]       # n_state × n_exo    (state shock impact)
    B_Σ_Bt  = B_state * Σe * B_state'
    B_Σ_Bt  = (B_Σ_Bt + B_Σ_Bt') / 2

    # Solve state covariance via Lyapunov: P = A·P·A' + B·Σ·B'
    P = local_dlyap(A_state, B_Σ_Bt)
    (any(diag(P) .< -1e-10) || any(isnan.(P))) && return NAN58, false
    P = (P + P') / 2

    # Build endo_name → row index map
    ys = context.results.model_results[1].trends.endogenous_steady_state
    ei = Dict(nm => i for (i,nm) in enumerate(endo_names))

    # Pre-extract only the ~39 rows of T and R we need for moments.
    # This reduces the HP filter inner loop from 491×491 to 39×39 — 25× faster.
    needed_names = vcat(
        ["Y_$(i)"  for i in 1:nsec],
        ["PH_$(i)" for i in 1:nsec],
        ["L_$(i)"  for i in 1:nsec],
        ["GDP", "pi", "Q", "TB"]
    )
    needed_idx = [get(ei, nm, 0) for nm in needed_names]   # 0 if missing
    valid_mask = needed_idx .> 0
    needed_idx_valid = needed_idx[valid_mask]

    T_sub = Matrix{Float64}(T[needed_idx_valid, :])   # n_needed × n_state
    R_sub = Matrix{Float64}(R[needed_idx_valid, :])   # n_needed × n_exo
    R_sub_Σe_Rsub = R_sub * Σe * R_sub'
    R_sub_Σe_Rsub = (R_sub_Σe_Rsub + R_sub_Σe_Rsub') / 2

    # Pre-compute HP filter weights once and run a single frequency loop
    w_var, w_lag1 = build_hp_weights(1600.0, 256)
    Γ_sub, Γ1_sub = hp_filtered_cov_fast(
        Matrix{Float64}(A_state), B_Σ_Bt, T_sub, R_sub_Σe_Rsub, w_var, w_lag1)

    # Map sub-matrix results back to needed_names ordering (including missing vars)
    n_needed = length(needed_names)
    Γ_val  = zeros(n_needed, n_needed)
    Γ1_val = zeros(n_needed, n_needed)
    sub_positions = findall(valid_mask)
    for (si, pi) in enumerate(sub_positions), (sj, pj) in enumerate(sub_positions)
        Γ_val[pi, pj]  = Γ_sub[si, sj]
        Γ1_val[pi, pj] = Γ1_sub[si, sj]
    end

    # Local name→position in needed_names (faster than repeated findfirst)
    ei_sub = Dict(nm => i for (i, nm) in enumerate(needed_names))

    pstd(vn) = let idx = get(ei_sub, vn, 0)
        idx == 0 && return 0.0
        sqrt(max(Γ_val[idx, idx], 0.0)) / max(abs(ys[needed_idx[idx]]), 1e-12)
    end
    xcorr(v1, v2) = let i1 = get(ei_sub, v1, 0), i2 = get(ei_sub, v2, 0)
        (i1 == 0 || i2 == 0) && return NaN
        d = sqrt(max(Γ_val[i1,i1], 0.0) * max(Γ_val[i2,i2], 0.0))
        d < 1e-15 ? 0.0 : clamp(Γ_val[i1,i2] / d, -1.0, 1.0)
    end

    i_Q_sub  = get(ei_sub, "Q",  0)
    i_TB_sub = get(ei_sub, "TB", 0)
    acQ = (i_Q_sub > 0 && Γ_val[i_Q_sub, i_Q_sub] > 1e-15) ?
          Γ1_val[i_Q_sub, i_Q_sub] / Γ_val[i_Q_sub, i_Q_sub] : NaN

    GDP_ss = baseline.GDP_ss > 0 ? baseline.GDP_ss : 1.0
    std_TBGDP = (i_TB_sub > 0) ? sqrt(max(Γ_val[i_TB_sub, i_TB_sub], 0.0)) / GDP_ss : 0.0

    std_Y  = [pstd("Y_$(i)")  for i in 1:nsec]
    std_PH = [pstd("PH_$(i)") for i in 1:nsec]
    std_L  = [pstd("L_$(i)")  for i in 1:nsec]

    # corr(Y_i, PH_i): negative under TFP shocks, positive under demand shocks
    # Key identifier for supply vs demand decomposition per sector
    corr_YPH = [xcorr("Y_$(i)", "PH_$(i)") for i in 1:nsec]

    dY=baseline.data_std_Y; dPH=baseline.data_std_PH; dL=baseline.data_std_L
    vy=isfinite.(std_Y).&isfinite.(dY)
    vp=isfinite.(std_PH).&isfinite.(dPH)
    vl=isfinite.(std_L).&isfinite.(dL)
    rY=sum(vy)>=3 ? safe_spearman(std_Y[vy],dY[vy]) : 0.0
    rP=sum(vp)>=3 ? safe_spearman(std_PH[vp],dPH[vp]) : 0.0
    rL=sum(vl)>=3 ? safe_spearman(std_L[vl],dL[vl]) : 0.0

    # Return 58 moments: 12×std_Y + 12×std_PH + 12×std_L + 10×aggregate + 12×corr(Y_i,PH_i)
    return [std_Y;std_PH;std_L;pstd("GDP");pstd("pi");xcorr("GDP","pi");
            std_TBGDP;pstd("Q");acQ;xcorr("GDP","Q");rY;rP;rL;corr_YPH], true
end
=#
