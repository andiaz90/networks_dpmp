"""
utils.jl
========
Helper functions for NK-IOSOE Julia translation.

Contents:
  - safe_spearman   : Spearman rank correlation that returns 0 when model vector is uniform
  - local_dlyap     : Discrete Lyapunov equation solver (for variance decomposition)
  - ternary_str     : Conditional string selection (replaces MATLAB ternary_str)
  - run_dynare_octave : Run Dynare via Octave system call
  - load_dynare_results : Read oo_ / M_ structures back from Dynare's output .mat
"""

using LinearAlgebra, StatsBase, MAT

# --------------------------------------------------------------------------- #
#  Rank correlation                                                             #
# --------------------------------------------------------------------------- #

"""
    safe_spearman(x, y) -> Float64

Spearman rank correlation between `x` and `y`.
Returns 0.0 when `x` is nearly uniform (model std devs identical across sectors),
mirroring the behaviour of smm_model_moments.m.
"""
function safe_spearman(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    # If all model values are (nearly) identical the rank order is undefined
    if std(x) < 1e-12 * (abs(mean(x)) + 1e-12)
        return 0.0
    end
    return corspearman(x, y)
end


# --------------------------------------------------------------------------- #
#  Discrete Lyapunov equation                                                  #
# --------------------------------------------------------------------------- #

"""
    local_dlyap(A, Q) -> P

Solve the discrete Lyapunov equation  P = A P Aᵀ + Q  by doubling iterations.
Equivalent to MATLAB's `dlyap(A', Q)` (note: dlyap solves Aᵀ P A − P + Q = 0).
Returns the symmetric positive-semidefinite solution.
"""
function local_dlyap(A::AbstractMatrix{<:Real}, Q::AbstractMatrix{<:Real}; maxiter::Int=500, tol::Real=1e-14)
    P = copy(Q)
    Ai = copy(A)
    for _ in 1:maxiter
        P_new = Ai * P * Ai' + P
        err   = norm(P_new - P, Inf) / (norm(P, Inf) + 1e-14)
        P     = (P_new + P_new') / 2   # enforce symmetry
        Ai    = Ai * Ai
        if err < tol
            break
        end
    end
    return (P + P') / 2
end


# --------------------------------------------------------------------------- #
#  Misc                                                                        #
# --------------------------------------------------------------------------- #

"""
    ternary_str(cond, s_true, s_false) -> String

Return `s_true` if `cond` is true, else `s_false`.
Replaces MATLAB's inline ternary_str function.
"""
ternary_str(cond::Bool, s_true::String, s_false::String) = cond ? s_true : s_false


# --------------------------------------------------------------------------- #
#  Dynare via Octave                                                           #
# --------------------------------------------------------------------------- #

"""
    run_dynare_octave(mod_dir, mod_name, dynare_matlab_path; octave_exe="octave")

Run `dynare('mod_name.mod', 'noclearall')` via Octave and wait for completion.
After Dynare finishes it saves `mod_name_results.mat` in `mod_dir`.

Arguments:
  mod_dir            : absolute path to the directory containing .mod files
  mod_name           : name WITHOUT extension (e.g. "NK_SOE_lev_gap2")
  dynare_matlab_path : path to Dynare's matlab/ subfolder
                       (e.g. "/usr/share/dynare/matlab" on Linux,
                             "C:/Program Files/Dynare/6.4/matlab" on Windows)
  octave_exe         : name or full path of the Octave executable (default "octave")
"""
function run_dynare_octave(
    mod_dir            :: String,
    mod_name           :: String,
    dynare_matlab_path :: String;
    octave_exe         :: String = "octave",
)
    # Build Octave eval string
    # We cd into mod_dir so that all @#include paths resolve correctly
    octave_script = """
    addpath('$(replace(dynare_matlab_path, "\\" => "/"))');
    cd('$(replace(mod_dir, "\\" => "/"))');
    dynare('$(mod_name).mod', 'noclearall');
    """

    cmd = Cmd([octave_exe, "--no-gui", "--eval", octave_script])
    println("  Running: $octave_exe --no-gui --eval \"...\"")
    flush(stdout)
    run(cmd)   # blocks until Dynare finishes
end


"""
    load_dynare_results(mod_dir, mod_name) -> (oo, M_model, options)

Load the `oo_`, `M_` and `options_` structs that Dynare writes to
`<mod_name>_results.mat` after a successful run.

Returns a named tuple with fields:
  .oo       : Dynare output structure (Dict)
  .M_model  : Dynare model structure  (Dict)
  .options  : Dynare options structure (Dict)
"""
function load_dynare_results(mod_dir::String, mod_name::String)
    results_file = joinpath(mod_dir, "$(mod_name)_results.mat")
    if !isfile(results_file)
        # Dynare also writes inside an output subfolder in some versions
        alt = joinpath(mod_dir, mod_name, "Output", "$(mod_name)_results.mat")
        if isfile(alt)
            results_file = alt
        else
            error("Dynare results file not found: $results_file\n" *
                  "Make sure Dynare ran successfully.")
        end
    end

    data = matread(results_file)
    oo      = get(data, "oo_",      Dict())
    M_model = get(data, "M_",       Dict())
    options = get(data, "options_", Dict())

    return (oo = oo, M_model = M_model, options = options)
end
