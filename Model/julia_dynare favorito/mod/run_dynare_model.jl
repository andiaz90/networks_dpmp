# run_dynare_model.jl
# =====================================================================
# This file is include()d by main_SOE_gap.jl at RUNTIME (inside _main())
# AFTER params_jl.mod has been written to this directory.
#
# Key: include() called inside a Julia function re-compiles and executes
# this file at that call point, so @dynare "NK_SOE_lev_gap2" runs the
# Dynare preprocessor with CWD = this directory (set by the caller).
# =====================================================================

using Dynare

context = @dynare "NK_SOE_lev_gap2"
