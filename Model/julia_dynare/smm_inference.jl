"""
smm_inference.jl
================
Asymptotic inference for the NK-IOSOE analytical-moment GMM estimator:
standard errors (GMM sandwich) and the overidentifying-restrictions J-test.

THEORY (analytical-moment GMM)
  θ̂ = argmin g(θ)' W g(θ),   g(θ) = data_moments − model_moments(θ),  K moments, P params.
  Because the model moments are an exact function of θ, all sampling uncertainty
  is in the DATA moments. Let S = asymptotic Cov(√T · data_moments) and
  G = ∂ model_moments/∂θ' (K×P, sign irrelevant in the sandwich). Then

     √T (θ̂ − θ0) →d N(0, V),
     V = (G'WG)^{-1} G'W S W G (G'WG)^{-1},     SE(θ̂_p) = sqrt(V_pp / T).

  Overidentification test (efficient weighting W = S^{-1}):
     J = T · g(θ̂)' S^{-1} g(θ̂)  →d  χ²(K − P).

INPUTS
  S (the K×K data-moment covariance) is what makes these SEs publication-grade.
  Rigorous route: provide Data/moment_cov.csv (K×K), e.g. from a moving-block
  bootstrap of the underlying series (resample the data, recompute the K moments
  B≥K times, take the sample covariance × T).
  Fallback (no file): a diagonal S from textbook delta-method variances per
  moment type — INDICATIVE only; reported with a clear warning.

OUTPUT
  estimation_results/smm_inference.csv  (param, estimate, se, t_stat) and a printed J-test.
"""

using LinearAlgebra, Printf, CSV, DataFrames

const SAMPLE_T = 72   # 2006Q1–2023Q4 quarters (matches compute_data_moments.jl)

# ---- χ² upper-tail p-value via regularized incomplete gamma Q(a,x) --------- #
# Numerical Recipes gammq (series + continued fraction). No Distributions dep.
function _gammq(a::Float64, x::Float64)
    (x < 0 || a <= 0) && return NaN
    x == 0 && return 1.0
    if x < a + 1.0                       # series for P(a,x), Q = 1 − P
        ap = a; sum = 1.0/a; del = sum
        for _ in 1:1000
            ap += 1.0; del *= x/ap; sum += del
            abs(del) < abs(sum)*1e-12 && break
        end
        return 1.0 - sum*exp(-x + a*log(x) - lgamma_(a))
    else                                  # continued fraction for Q(a,x)
        b = x + 1.0 - a; c = 1e300; d = 1.0/b; h = d
        for i in 1:1000
            an = -i*(i - a); b += 2.0
            d = an*d + b; abs(d) < 1e-300 && (d = 1e-300)
            c = b + an/c; abs(c) < 1e-300 && (c = 1e-300)
            d = 1.0/d; del = d*c; h *= del
            abs(del - 1.0) < 1e-12 && break
        end
        return exp(-x + a*log(x) - lgamma_(a))*h
    end
end
# log Γ via Lanczos (avoid SpecialFunctions dependency)
function lgamma_(x::Float64)
    g = 7.0
    c = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
         771.32342877765313, -176.61502916214059, 12.507343278686905,
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
    x -= 1.0
    a = c[1]; t = x + g + 0.5
    for i in 2:9; a += c[i]/(x + (i-1)); end
    return 0.5*log(2π) + (x + 0.5)*log(t) - t + log(a)
end
chisq_pvalue(stat::Real, dof::Integer) = _gammq(dof/2, stat/2)


# ---- diagonal fallback S: delta-method asymptotic variances of √T·moments -- #
# std-dev moment σ̂: √T(σ̂−σ)→N(0, σ²/2). correlation ρ̂: √T(ρ̂−ρ)→N(0,(1−ρ²)²).
# Rank correlations (cross-sectional over NSEC sectors): variance ≈ 1/(NSEC−1).
function _fallback_S(dm::Vector{<:Real})
    K = length(dm); s = ones(K)
    stdpos = vcat(1:36, [37, 38, 40, 41])            # std-dev type moments
    corrpos = vcat([39, 42, 43], 47:58, [59, 60])    # correlation type moments (59–60 = corr(N,GDP), corr(N,GDP/N))
    rankpos = [44, 45, 46]                           # cross-sectional rank corr
    for k in stdpos;  s[k] = max(dm[k]^2/2, 1e-10);              end
    for k in corrpos; s[k] = max((1 - dm[k]^2)^2, 1e-4);        end
    for k in rankpos; s[k] = 1.0/max(NSEC - 1, 1);              end
    return Matrix(Diagonal(s))
end

function _load_or_build_S(dm::Vector{<:Real})
    f = joinpath(DATA_DIR, "moment_cov.csv")
    if isfile(f)
        M = Matrix{Float64}(CSV.read(f, DataFrame))
        if size(M) == (length(dm), length(dm))
            @printf "  Inference: using bootstrap moment covariance %s\n" basename(f)
            return (M + M')/2, true
        end
        @printf "  [warn] %s has wrong size %s; using fallback S.\n" basename(f) string(size(M))
    end
    @printf "  Inference: NO moment_cov.csv — using delta-method diagonal S (INDICATIVE SEs).\n"
    @printf "             For publication SEs, supply Data/moment_cov.csv (block bootstrap).\n"
    return _fallback_S(dm), false
end


"""
    compute_smm_inference(θ_hat, m_hat, context, baseline, endo_names; T=SAMPLE_T)

Finite-differences the moment Jacobian at θ̂, forms the GMM sandwich covariance,
and computes SEs and the overidentifying J-test. Saves estimation_results/smm_inference.csv.
"""
function compute_smm_inference(θ_hat, m_hat, context, baseline, endo_names; T::Int=SAMPLE_T)
    @printf "\n%s\n  SMM INFERENCE (standard errors + J-test)\n%s\n" repeat("=",60) repeat("=",60)
    W = build_weighting_matrix(data_moments)
    g = data_moments .- m_hat
    S, have_boot = _load_or_build_S(data_moments)

    # Moment Jacobian G = ∂ model_moments / ∂θ  (K×P), central differences.
    K = N_MOMENTS; P = N_THETA
    G = zeros(K, P)
    for p in 1:P
        h  = 1e-4 * max(1.0, abs(θ_hat[p]))
        θp = copy(θ_hat); θm = copy(θ_hat)
        θp[p] = min(θ_hat[p] + h, UB[p]); θm[p] = max(θ_hat[p] - h, LB[p])
        denom = θp[p] - θm[p]
        mp, okp = smm_model_moments(θp, context, baseline, endo_names)
        mm, okm = smm_model_moments(θm, context, baseline, endo_names)
        if okp && okm && denom > 0 && !any(isnan, mp) && !any(isnan, mm)
            G[:, p] = (mp .- mm) ./ denom
        else
            @printf "  [warn] Jacobian column %d (%s) unreliable (model failed at θ±h).\n" p PARAM_LABELS[p]
        end
    end

    # Sandwich covariance of √T(θ̂−θ0).
    bread = Symmetric(G' * W * G)
    bread_inv = pinv(Matrix(bread))
    meat = G' * W * S * W * G
    Vθ = bread_inv * meat * bread_inv
    se = sqrt.(max.(diag(Vθ) ./ T, 0.0))
    tstat = θ_hat ./ se

    # Overidentification J-test (efficient weighting S^{-1}).
    Sinv = pinv(S)
    Jstat = T * (g' * Sinv * g)
    dof = max(K - P, 1)
    pval = chisq_pvalue(Jstat, dof)

    # Report
    @printf "\n  %-3s  %-16s  %10s  %10s  %8s\n" "#" "Parameter" "Estimate" "Std.Err" "t"
    @printf "  %s\n" repeat("-", 54)
    for p in 1:P
        @printf "  %-3d  %-16s  %10.4f  %10.4f  %8.2f\n" p PARAM_LABELS[p] θ_hat[p] se[p] tstat[p]
    end
    @printf "\n  J-test of overidentifying restrictions:\n"
    @printf "    J = %.3f,  df = %d (= %d moments − %d params),  p-value = %.4f%s\n" Jstat dof K P pval (have_boot ? "" : "  [INDICATIVE — diagonal S]")
    @printf "    %s\n" (pval < 0.05 ? "Overidentifying restrictions REJECTED at 5% (model misspecification)." :
                                      "Fail to reject overid restrictions at 5%.")

    df = DataFrame(param=PARAM_LABELS, estimate=collect(Float64, θ_hat),
                   std_err=se, t_stat=tstat)
    try
        atomic_write_csv(joinpath(ESTIMATION_DIR, "smm_inference.csv"), df)
        @printf "\n  Saved: %s\n\n" joinpath(ESTIMATION_DIR, "smm_inference.csv")
    catch err
        @printf "  [warn] could not save smm_inference.csv: %s\n" sprint(showerror, err)
    end
    return se, Jstat, dof, pval
end
