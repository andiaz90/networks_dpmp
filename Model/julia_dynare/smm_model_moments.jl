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

"""
    get_power_deriv(x, p, k) -> d^k(x^p)/dx^k

Dynare's helper for derivatives of a power. The generated SparseDynamicG1!.jl
CALLS it (e.g. `get_power_deriv(y[594], -params[1], 1)`) but no generated file
DEFINES it, and it is not exported by the installed Dynare.jl either — Dynare
normally supplies it inside its own module, whereas these files are included
into Main. So every call raised

    UndefVarError: `get_power_deriv` not defined in `Main`

ROOT CAUSE OF THE FLAT OBJECTIVE (found 2026-08-19). _eval_dynamic_jacobian
therefore ALWAYS threw, the Klein path could never run, and
resolve_first_order! fell through to its last resort: the decision rule already
in the context. That fallback returns success while returning a θ-INDEPENDENT
solution, so the SMM objective was exactly constant in all 7 free parameters.
CMA-ES "converged" in ~126 evaluations with θ unchanged, every Jacobian column
was zero (std err 0.0, t = Inf), and both the July and August runs estimated
nothing while reporting success. main_SOE_gap.jl was never affected: it rewrites
params_jl.mod and re-runs Dynare in a subprocess.

Semantics match Dynare's own implementation exactly, including the guard that
returns 0 for a zero base with a positive integer exponent below the derivative
order (where x^(p-k) would otherwise be a division by zero).
"""
function get_power_deriv(x::Real, p::Real, k::Integer)
    if abs(x) < 1e-12 && p > 0 && k > p && abs(p - round(p)) < 1e-12
        return 0.0
    end
    dxp = float(x)^(float(p) - k)
    pp  = float(p)
    for _ in 1:k
        dxp *= pp
        pp  -= 1.0
    end
    return dxp
end

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
# Shock-amplitude params (sigma_om, isigma_tfp_i, sigma_pvstar, sigma_zeta)
# enter the model as LINEAR coefficients in the model equations and hence
# in g1_2.  When ONLY shock amplitudes change (structural params are
# unchanged), g1_1 is identical and g1_2 scales proportionally.  We cache
# the last Klein result and skip the 134×134 QZ when structural params
# change by less than _KLEIN_THRESH.
#
# Structural params that DO require Klein re-solve (θ indices, Option-A layout):
#   1=ilabcosts, 2=epsY, 3=epsM, 4=log(kappaV), 5=rho_om, 6=rho_A,
#   31=rho_pvstar, 33=rho_zeta, 36=kappaw (wage Rotemberg cost enters the wage-PC
#   Jacobian directly — a stale cache here would silently freeze kappaw)
#
# 2026-08-21 — THE ABOVE REASONING WAS HALF RIGHT, AND THE HALF THAT WAS WRONG
# SILENTLY FROZE THE SHOCK SIZES.
#
# It is true that g1_1 (the decision rule T) is invariant to shock standard
# deviations. It is NOT true that this makes the cache safe, because
# _resolve_cached! returns T *and R* from the same cache, and R is where the
# sigmas live. In NK_SOE_lev_gap2.mod the shock sizes are parameters INSIDE the
# equations:
#
#     om_@{i}    = rho_om1*om_@{i}(-1) + sigma_om_@{i}*eps_om_@{i};        (:820)
#     exp(A_@{i}) = ... + isigma_tfp_@{i}*epsA_@{i};                       (:830)
#
# and the moment code sets Sigma_e to the IDENTITY on active shocks. So every
# sigma reaches the moments through R, and R alone.
#
# With sigmas excluded from the cache key, any evaluation that changes ONLY
# shock sizes hit the cache and returned the previous R — i.e. the previous
# shock sizes. Symptom observed 2026-08-21: calibrate_sectoral_shocks.jl, whose
# whole job is to vary theta[7:30] and nothing else, reported a byte-identical
# rescale ratio across all four iterations while the shock sizes moved by a
# factor of 0.6. It could never converge, because the model never saw the update.
#
# CMA-ES mostly escaped this because its candidates perturb several structural
# parameters at once, so the key almost always missed — and the pre-flight
# sensitivity scan escaped it by accident, because each perturbation leaves the
# cache keyed on the PREVIOUS parameter's perturbed value, forcing a miss on the
# next one. Neither is a defence; both are luck.
#
# Every element of theta affects T or R, so the only correct key is the whole
# vector. The cache now serves its remaining honest purpose — skipping a
# re-solve when the identical theta is evaluated twice, which happens on the
# decomposition re-evaluation and on repeated candidates. Measured cost is nil:
# the Klein cache hit rate was already 0-2% during estimation.
const _KLEIN_STRUCT_IDX = collect(1:38)

# Steady-state solve tolerance. Was hardcoded 1e-12 in three places, which is
# far tighter than anything downstream needs: the SS is the point a FIRST-ORDER
# approximation is taken around, and the moments are accurate to O(1e-6) at
# best. 1e-10 is still two orders tighter than the model's own accuracy and
# cuts trust-region iterations materially — which matters now that freeing epsY
# makes this run on nearly every evaluation. SMM_SS_FTOL overrides.
const _SS_FTOL = parse(Float64, get(ENV, "SMM_SS_FTOL", "1e-10"))
const _KLEIN_THRESH     = 1e-5   # re-solve if any structural param moves > this

const _KLEIN_CACHE_T    = [Ref{Matrix{Float64}}(zeros(0,0)) for _ in 1:_N_THREADS]
const _KLEIN_CACHE_R    = [Ref{Matrix{Float64}}(zeros(0,0)) for _ in 1:_N_THREADS]
const _KLEIN_CACHE_ΘSTR = [Ref{Vector{Float64}}(Float64[]) for _ in 1:_N_THREADS]
const _KLEIN_HITS       = Threads.Atomic{Int}(0)
const _KLEIN_MISSES     = Threads.Atomic{Int}(0)
# Counts evaluations that fell back to the decision rule already in the context,
# i.e. the one at θ_baseline rather than at the requested θ. Any non-zero value
# means some evaluations did not reflect θ (2026-08-19).
const _STALE_DR_USES    = Threads.Atomic{Int}(0)
# Klein aborts. The Klein path is dead code whenever Dynare's solver works, so
# its diagnostics are printed ONCE and then counted — printing per evaluation
# buried the log. Totals are reported at the end of the run: silent failure is
# what made the flat-objective bug survive two estimation runs, so the
# information is suppressed, never discarded.
const _KLEIN_ABORTS     = Threads.Atomic{Int}(0)
# One-shot flag so fast-path diagnostics print once, not once per evaluation.
const _FASTPATH_REPORTED = Ref(false)

# Per-thread Dynare solver workspaces (2026-08-19). Constructing these is not
# free, so they are built once per thread and reused. They depend only on model
# DIMENSIONS, not on parameter values, so caching them across θ is safe.
const _DYN_WS     = [Ref{Any}(nothing) for _ in 1:_N_THREADS]
const _CSS_WS     = [Ref{Any}(nothing) for _ in 1:_N_THREADS]
const _SS_OPTS    = Ref{Any}(nothing)

"""
    _dynare_solve!(context, params) -> (ok, g1_1, g1_2, Sigma_e)

Call Dynare.jl's own first-order solver.

The previous code invoked `compute_first_order_solution!(context)`, which has no
such method — it takes nine arguments in this version (Dynare TkH8b,
perturbations.jl:637). The MethodError was swallowed by a bare `catch`, so the
fast path appeared to "not work on this platform" when it had simply never been
called correctly.

Why this rather than the hand-rolled Klein: `_klein_solve` assumes the backward
and forward variable sets are DISJOINT (it forms `[y_b y_f]` and treats it as
square). That held for the old 491-variable model (78 + 56 = 134) but not for
the current one, where |bkwrd| = 92, |fwrd| = 58 and the union is 112 — i.e. 38
variables are MIXED. Handling mixed variables correctly is the fiddly heart of
any Klein implementation, and Dynare already does it.
"""
function _dynare_solve!(context, params)
    tid = _tid()
    m   = context.models[1]
    res = context.results.model_results[1]

    if _DYN_WS[tid][] === nothing
        _DYN_WS[tid][] = Dynare.DynamicWs(context; order = 1)
        _CSS_WS[tid][] = Dynare.ComputeStochSimulWs(context)
    end
    if _SS_OPTS[] === nothing
        _SS_OPTS[] = Dynare.StochSimulOptions(Dict{String,Any}("order" => 1, "irf" => 0))
    end

    # Argument convention copied verbatim from Dynare's own caller,
    # compute_stoch_simul! (perturbations.jl:547-551):
    #
    #     endogenous  = results.trends.endogenous_steady_state
    #     set_endogenous3!(css_ws.endogenous3, endogenous)
    #     exogenous   = results.trends.exogenous_steady_state
    #     compute_first_order_solution!(context, endogenous3, exogenous,
    #                                   endogenous, params, model, ws, css_ws, options)
    #
    # i.e. the first argument is the STACKED [y_{t-1}; y_t; y_{t+1}] buffer
    # css_ws.endogenous3 (3*591 = 1773 — hence "AssertionError: length(y) ==
    # 1773" when the flat steady state was passed), the third is the flat one.
    # Dynare also zeroes the exogenous steady state first, so we do too.
    css     = _CSS_WS[tid][]
    endo_ss = res.trends.endogenous_steady_state
    exo_ss  = res.trends.exogenous_steady_state
    fill!(exo_ss, 0.0)
    Dynare.set_endogenous3!(css.endogenous3, endo_ss)

    Dynare.compute_first_order_solution!(
        context, css.endogenous3, exo_ss, endo_ss, params, m,
        _DYN_WS[tid][], css, _SS_OPTS[];
        variance_decomposition = false)

    lre = res.linearrationalexpectations
    (isempty(lre.g1_1) || isempty(lre.g1_2)) && return false, zeros(0,0), zeros(0,0), zeros(0,0)
    return true, Matrix{Float64}(lre.g1_1), Matrix{Float64}(lre.g1_2), m.Sigma_e
end

function _resolve_cached!(context, θ)
    tid = _tid()
    # Tolerate legacy θ vectors shorter than the newest struct index (e.g. a
    # 35-length θ without kappaw): drop out-of-range entries instead of erroring.
    θ_str = θ[filter(<=(length(θ)), _KLEIN_STRUCT_IDX)]
    cache_valid = !isempty(_KLEIN_CACHE_ΘSTR[tid][]) &&
                  length(_KLEIN_CACHE_ΘSTR[tid][]) == length(θ_str) &&
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
    # Cross-check against the Jacobian's own equation count rather than a
    # hardcoded 491 (stale since the model grew to 591 equations, 2026-08-19).
    _N_EQ[] > 0 && length(rows) != _N_EQ[] &&
        @printf "  [SMM] Warning: parsed %d lead_lag rows, Jacobian has %d equations\n" length(rows) _N_EQ[]
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
monetary (eps_i), import price (eps_pvstar), aggregate demand (eps_zeta), world copper
price (eps_pc), the 12
sectoral TFP shocks (epsA_i) and the 12 sectoral demand shocks (eps_om_i).
Deliberately EXCLUDES eps_postar (oil — off during estimation) and epschi
(labour supply — off). Returns the indices in the model's exogenous order.
"""
function active_shock_indices(context, nsec::Int)
    exo = smm_exo_names(context)
    # eps_pc = world copper price; eps_omg = goods/services demand reallocation
    # (FGI 2023 omega_t, added 2026-08-19). MUST match the active_set in
    # main_SOE_gap.jl — the two are separate literals and have drifted before.
    active = Set{String}(["eps_i", "eps_pvstar", "eps_zeta", "eps_pc", "eps_omg"])
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
    # ---------------------------------------------------------------------- #
    # CRITICAL (2026-08-19): push our parameter vector INTO the context before
    # any solver runs.
    #
    # set_param! writes to _SMM_PARAMS (see the note at ~line 204: context.work
    # .params was unreliable on ARM Macs, so this file keeps its own vector).
    # But the fast path below calls Dynare's OWN solver, which reads
    # context.work.params. When that field is readable — as it is on this
    # machine, cf. "[SMM] Params loaded from context.work.params (655 params)" —
    # the fast path succeeds, returns immediately, and Klein never runs. Every
    # evaluation then returns the decision rule at the ORIGINAL parameters,
    # whatever θ says.
    #
    # Consequence: the SMM objective was EXACTLY CONSTANT in all 7 free
    # parameters. A full-range probe confirmed zero spread for every one of
    # them, including sigma_zeta = 0 (which removes a shock worth 10.4% of
    # employment variance) and kappaw 0->400. That is why CMA-ES terminated
    # after ~126 evaluations with best_obj never improving, why every Jacobian
    # column was zero (std err 0.0, t = Inf), and why the July run behaved
    # identically. main_SOE_gap.jl was unaffected because it rewrites
    # params_jl.mod and re-runs Dynare in a subprocess.
    # ---------------------------------------------------------------------- #
    try
        sp = _SMM_PARAMS[_tid()][]
        if !isempty(sp)
            wp = context.work.params
            length(wp) == length(sp) && copyto!(wp, sp)
        end
    catch
        # context.work.params unreadable (the ARM case the design anticipated):
        # the fast path will fail too and we fall through to Klein, which reads
        # _SMM_PARAMS directly. Nothing to do.
    end

    try
        # Fast path: try Dynare.jl's own solver (works on some platforms).
        # Failures were silently swallowed; they are now reported once each, so
        # we can tell whether this path is unavailable or merely untried
        # (2026-08-19). Dynare's solver demonstrably works on this machine in
        # main_SOE_gap.jl's subprocess, so knowing WHY it fails in-process is
        # worth more than reimplementing Klein.
        try
            ok, g1_1, g1_2, Se = _dynare_solve!(context, _load_smm_params!(context))
            if ok
                _FASTPATH_REPORTED[] || (@printf "  [SMM] using Dynare's own first-order solver\n"; _FASTPATH_REPORTED[] = true)
                return true, g1_1, g1_2, Se
            end
            _FASTPATH_REPORTED[] || (@printf "  [SMM] Dynare solver returned an empty decision rule; falling back\n"; _FASTPATH_REPORTED[] = true)
        catch err
            if !_FASTPATH_REPORTED[]
                @printf "  [SMM] Dynare solver FAILED -> %s\n" first(split(sprint(showerror, err), '\n'))
                _FASTPATH_REPORTED[] = true
            end
        end

        # Klein (2000) pure-Julia path
        ok, g1_1, g1_2, Σe = _klein_solve(context)
        ok && return true, g1_1, g1_2, Σe

        # Fallback: the decision rule ALREADY in the context — i.e. the one at
        # θ_baseline, NOT at the requested θ.
        #
        # This is the same hazard as the fast path above (2026-08-19): it
        # reports success while returning a solution that does not depend on θ,
        # so the objective goes flat and the estimator silently converges at its
        # starting point. The original note claimed it "allows estimation of
        # shock parameters (which enter via Σe, not g1)", but in THIS model the
        # shock sizes are in-equation parameters (isigma_tfp_i * epsA_i), so
        # they enter g1_2 — and the Σe returned here is stale too.
        #
        # Kept as a last resort so a single bad evaluation cannot kill a long
        # run, but it now announces itself. If this fires often, the estimates
        # are not trustworthy.
        mr  = context.results.model_results[1]
        lre = mr.linearrationalexpectations
        if !isempty(lre.g1_1) && size(lre.g1_1, 1) >= 400
            Threads.atomic_add!(_STALE_DR_USES, 1)
            if _STALE_DR_USES[] <= 3 || _STALE_DR_USES[] % 500 == 0
                @printf "  [WARN] resolve_first_order! fell back to the STALE decision rule (use #%d) — this evaluation does NOT reflect θ.\n" _STALE_DR_USES[]
            end
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
        n_eq = _N_EQ[]

        # ------------------------------------------------------------------ #
        # Column partition, DERIVED from the loaded model (2026-08-19).
        #
        # WAS hardcoded: n_bk=78; n_endo=491; n_fw=56; n_exo=17 — the dimensions
        # of a much older model. The current one is 92/591/58/31. The immediate
        # consequence was that the fwrd_b length check below (== n_fw) failed,
        # _klein_solve returned false, and resolve_first_order! fell through to
        # the stale decision rule — which reports success while ignoring θ. That
        # is what made the SMM objective exactly constant in all 7 free
        # parameters and caused two estimation runs to "converge" instantly
        # having estimated nothing.
        #
        # Everything is now read from the model, and cross-checked against the
        # Jacobian's own column count so a future model change fails loudly
        # instead of silently degrading to a θ-independent solution.
        # ------------------------------------------------------------------ #
        n_endo = length(context.results.model_results[1].trends.endogenous_steady_state)
        n_exo  = context.models[1].exogenous_nbr
        n_bk   = length(context.models[1].i_bkwrd_b)
        n_fw   = _N_COL[] - n_bk - n_endo - n_exo
        if n_fw <= 0 || n_bk + n_endo + n_fw + n_exo != _N_COL[]
            Threads.atomic_add!(_KLEIN_ABORTS, 1)
            _KLEIN_ABORTS[] == 1 && @printf "  [SMM] Klein: inconsistent column partition (bk=%d endo=%d fw=%d exo=%d, Jacobian has %d cols)\n" n_bk n_endo n_fw n_exo _N_COL[]
            return false, zeros(0,0), zeros(0,0), zeros(0,0)
        end
        # A: 1:n_bk | B: next n_endo | C: next n_fw | D: last n_exo
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
        if isnothing(fwrd_b)
            Threads.atomic_add!(_KLEIN_ABORTS, 1)
            _KLEIN_ABORTS[] == 1 && @printf "  [SMM] Klein: could not resolve %d forward-looking indices (i_fwrd_b=%d, modfile=%d) — ABORTING rather than returning a stale rule.\n" n_fw (isdefined(context.models[1], :i_fwrd_b) ? length(context.models[1].i_fwrd_b) : -1) length(_I_FWRD_B)
            return false, zeros(0,0), zeros(0,0), zeros(0,0)
        end

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
            Threads.atomic_add!(_KLEIN_ABORTS, 1)
            _KLEIN_ABORTS[] == 1 && @printf "  [Klein] Static elimination mismatch: %d static eqs, %d static vars (expected both=%d) — Klein path unusable; further occurrences counted silently\n" length(s_eq) n_s n_s
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
            Threads.atomic_add!(_KLEIN_ABORTS, 1)
            _KLEIN_ABORTS[] == 1 && @printf "  [Klein] BK failed on the reduced system: %d stable eigenvalues, expected %d\n" n_stable n_bk
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
            Threads.atomic_add!(_KLEIN_ABORTS, 1)
            _KLEIN_ABORTS[] == 1 && @printf "  [Klein] Exception type: %s\n" typeof(e)
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

    # USE-BEFORE-DEFINITION BUG, fixed 2026-08-21.
    #
    # _sub was assigned ~40 lines below, in the write-back block, but it is
    # captured by the f_outer! closure just below and therefore read the moment
    # nlsolve first calls it. Julia raises UndefVarError, the try/catch wrapping
    # this function turns that into `return false`, and the caller emits NaN
    # moments. Result: recompute_ss! failed 100% of the time, silently.
    #
    # It went unnoticed because this function only runs when epsY, epsM or
    # etastar MOVE, and all three have been pinned since the theta-reduction —
    # so the whole path was dead code. Freeing epsY on 2026-08-21 executed it for
    # the first time and every single evaluation failed (110/110, ~30 ms each,
    # far too fast for a real nonlinear solve — that timing was the tell).
    _sub = hasproperty(baseline, :subsMC_val) && baseline.subsMC_val > 0 ?
           baseline.subsMC_val : 1.0

    f_outer! = (F, x) -> F .= steady_ntwsoe(x, baseline.PVstar_ss,
                                      (baseline.epsilon_val-1)/(baseline.epsilon_val*_sub),
                                      baseline.modvarrho, baseline.sigmaH_val,
                                      baseline.modgammag, baseline.modgammas,
                                      baseline.ombar_val, 1-baseline.ombar_val,
                                      baseline.modchiX, baseline.omegaX_val,
                                      eta_for_ss, baseline.ystar_ss_val,
                                      baseline.modalpha, baseline.modalphaV,
                                      baseline.modbeta, modepsY, modepsM,
                                      baseline.modgammaG, baseline.gshare_target,
                                      baseline.modalphaK, baseline.modchiI,
                                      baseline.nuK_val,
                                      baseline.gamma_val, baseline.chi_val,
                                      baseline.psi_val, ones(nsec), baseline.tb_target)

    # +1 unknown since 2026-08-20: the government-demand scale. Seeded from the
    # solved Gi_ss in params_jl.mod so the warm start stays warm.
    _g0 = sum(baseline.modGi) > 0 ? sum(baseline.modGi) : 0.20
    x0 = [[get_ss("PH_$(i)") for i in 1:nsec]; get_ss("w"); get_ss("Q"); get_ss("C"); _g0]
    res = nlsolve(f_outer!, x0; ftol=_SS_FTOL, method=:trust_region, show_trace=false)
    if !converged(res)
        # Multi-start: the trust-region solver is sensitive to the warm guess
        # when epsY/epsM/etastar move far from the cached point. Retry from a
        # few deterministically-perturbed guesses before giving up. Seeded RNG
        # so failures are reproducible across runs.
        rng = Random.MersenneTwister(hash((epsY, epsM, eta_for_ss)) % UInt32)
        for k in 1:6
            x0p = max.(x0 .* (1 .+ 0.10 * k .* (rand(rng, length(x0)) .- 0.5)), 1e-8)
            res = nlsolve(f_outer!, x0p; ftol=_SS_FTOL, method=:trust_region, show_trace=false)
            converged(res) && break
        end
        converged(res) || return false
    end
    all(isfinite, res.zero) || return false

    pH_ss = res.zero[1:nsec]; w_ss = res.zero[nsec+1]
    Q_ss = res.zero[nsec+2];  C_ss = res.zero[nsec+3]
    PL_ss = fill(w_ss, nsec); PV_ss = Q_ss * baseline.PVstar_ss
    # LPR markup subsidy: MC/PH = (eps-1)/(eps*subsMC). subsMC = 1 is the
    # no-subsidy case, so a params file predating the subsidy still works.
    # _sub is defined above, before the f_outer! closure that captures it.
    MCi_ss = (baseline.epsilon_val-1)/(baseline.epsilon_val*_sub) .* pH_ss
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

    # ESTIMATION mode (see steady_ntwsoe_system.jl): the endowment Kbar is held
    # FIXED at its calibrated value, so the system stays 4*nsec and the capital
    # cost share is free to drift with the estimated parameters — as in
    # Baqaee-Farhi and LPR, where Kbar_f is a datum.
    inner = nlsolve(
        (F, x) -> steady_ntwsoe_system!(F, x, baseline.modalpha, baseline.modalphaV,
                                         baseline.modbeta, MCi_ss, PMi_ss, PL_ss, PV_ss,
                                         CHi, Xi, modepsY, modepsM, ones(nsec), pH_ss,
                                         baseline.modalphaK, baseline.modchiI,
                                         baseline.nuK_val, baseline.modKbar),
        vcat(max.(MCi_ss./PMi_ss,1e-20), max.(MCi_ss./PL_ss,1e-20),
             max.(MCi_ss./PV_ss,1e-20), fill(0.2,nsec));
        ftol=_SS_FTOL, method=:trust_region, show_trace=false)
    (!converged(inner) || !all(isfinite, inner.zero)) && return false

    Yi_ss = inner.zero[3*nsec+1:4*nsec]; M_ss = inner.zero[1:nsec]
    IMP = sum(inner.zero[2*nsec+1:3*nsec]) + sum(CFi)
    TB  = PX*X - PV_ss*IMP; GDP = C_ss + TB
    Pistar = 1.0; r_star = Pistar/baseline.beta_val
    Bstar  = -TB / (Q_ss*(1 - r_star/Pistar))

    set_param!(context, "bbar", Q_ss*Bstar/GDP)
    set_param!(context, "Y_ss", sum(pH_ss .* Yi_ss))  # constant-SS-price gross output (matches Y in .mod)
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
    # GHH preferences: rescale the disutility so the recomputed SS still satisfies
    # MRS = chi0*N^psi = w, then refresh the composite Zc = C - chi0*N^(1+psi)/(1+psi)
    # and aggregate labor N so the linearization point solves the new preference
    # block exactly (chi = 1 at SS). chi0 tracks C_ss as epsY/epsM/etastar move.
    N_new    = sum(@view inner.zero[nsec+1:2*nsec])
    chi0_new = baseline.chi_val * C_ss^baseline.gamma_val
    Zc_new   = C_ss - chi0_new * N_new^(1 + baseline.psi_val) / (1 + baseline.psi_val)
    (isfinite(chi0_new) && isfinite(Zc_new) && Zc_new > 0) || return false
    set_param!(context, "chi0", chi0_new)
    upd("N", N_new); upd("Zc", Zc_new); upd("Zc_f", Zc_new)
    return true
  catch
    # Any unexpected failure in the SS recompute is treated as a non-convergence
    # so the optimiser simply penalises this θ rather than aborting the run.
    return false
  end
end


# =========================================================================== #
#  NOTE: this file holds INFRASTRUCTURE ONLY.                                  #
#  ------------------------------------------------------------------------   #
#  The canonical MOMENT_NAMES and smm_model_moments(...) live in               #
#  smm_estimation.jl. A second, divergent copy of both used to sit here,       #
#  block-commented since the 2026-07 cleanup; deleted 2026-08-21 (recover      #
#  from git if the derivation is ever needed). Everything ABOVE — the Klein    #
#  solver, recompute_ss!, param access, the caches — is live and used by the   #
#  canonical moment function.                                                  #
# =========================================================================== #
