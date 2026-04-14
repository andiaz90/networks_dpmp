function [moments, info_flag, var_names_out] = smm_model_moments(theta, M_in, options_in, oo_in, baseline)
% SMM_MODEL_MOMENTS  Compute theoretical model moments for a given parameter vector.
%
% Uses Dynare's linearised state-space to obtain exact unconditional second
% moments without simulation, then maps them to the 42-element target vector
% described in the estimation section of the paper.
%
% INPUTS
%   theta    - [10x1] parameter vector (see order below)
%   M_in     - Dynare M_ structure (after initial compilation)
%   options_in - Dynare options_ structure
%   oo_in    - Dynare oo_ structure (after initial compilation)
%   baseline - struct with fixed calibrated objects:
%              .nsec, .goods, .services, .Y_ss, .GDP_ss
%              .modalpha, .modalphaV, .modbeta, .modgammag, .modgammas
%              .modvarrho, .modchiX, .modkappa
%              .gamma_val, .psi_val, .chi_val, .epsilon_val, .beta_val
%              .ombar_val, .tb_target, .PVstar_ss, .sigmaH_val
%              .etastar_val, .omegaX_val, .ystar_ss_val
%              .steady_fn  - handle to steady_ntwsoe (or [])
%
% OUTPUTS
%   moments      - [42x1] model moment vector (same ordering as data moments)
%   info_flag    - 0 = success, >0 = failure code from Dynare resol
%   var_names_out - cell array describing the 15 moments (for diagnostics)
%
% THETA ORDERING
%   theta(1)     ilabcosts      inverse aggregate labour adjustment cost
%   theta(2)     epsY           elasticity of substitution in production (common)
%   theta(3)     epsM           elasticity of substitution between materials (common)
%   theta(4)     log_kappaV     log of import price adjustment cost kappa_v
%   theta(5)     rho_om         AR(1) persistence of goods-services reallocation shock
%   theta(6)     sigma_om       std dev of goods-services reallocation shock
%   theta(7)     rho_A          AR(1) persistence of sectoral TFP shocks (common)
%   theta(8:19)  isigma_tfp_i   std dev of each sector's TFP shock (i=1,...,12)
%   theta(20)    rho_pvstar     AR(1) persistence of import price (PVstar) shock
%   theta(21)    sigma_pvstar   std dev of import price shock

    nsec     = baseline.nsec;
    goods    = baseline.goods;      % logical, 1 for goods sectors (1-5)
    services = baseline.services;   % logical, 1 for services sectors (6-12)

    % ------------------------------------------------------------------ %
    %  1. Unpack theta and enforce bounds                                 %
    % ------------------------------------------------------------------ %
    ilabcosts  = theta(1);
    epsY       = theta(2);
    epsM       = theta(3);
    kappaV     = exp(theta(4));   % stored as log for numerical stability
    rho_om     = theta(5);
    sigma_om   = theta(6);
    rho_A      = theta(7);
    isigma_tfp = theta(8:19);   % [nsec×1] sector-specific TFP shock std devs
    rho_pvstar   = theta(20);   % AR(1) persistence of import price shock
    sigma_pvstar = theta(21);   % std dev of import price shock
    rho_xi       = theta(22);   % AR(1) persistence of preference/demand shock
    sigma_xi     = theta(23);   % std dev of preference/demand shock

    % Hard feasibility checks (return large penalty if violated)
    if epsY <= 0 || epsY >= 5 || epsM <= 0 || epsM >= 2 || ...
       ilabcosts <= 0 || kappaV <= 0 || abs(rho_om) >= 1 || ...
       sigma_om < 0 || abs(rho_A) >= 1 || any(isigma_tfp < 0) || ...
       abs(rho_pvstar) >= 1 || sigma_pvstar < 0 || ...
       abs(rho_xi) >= 1 || sigma_xi < 0
        moments  = NaN(46,1);
        info_flag = 99;
        var_names_out = get_moment_names();
        return
    end

    % ------------------------------------------------------------------ %
    %  2. Copy Dynare structures (avoid contaminating globals)            %
    % ------------------------------------------------------------------ %
    M_loc      = M_in;
    oo_loc     = oo_in;
    options_loc = options_in;

    % ------------------------------------------------------------------ %
    %  3. Update M_.params with new theta                                 %
    % ------------------------------------------------------------------ %
    % IMPORTANT: do NOT use an anonymous-function helper like
    %   set_param = @(name,val) set_dynare_param(M_loc, name, val)
    % MATLAB anonymous functions capture variables BY VALUE at creation
    % time, so the closure would permanently hold the original M_loc = M_in.
    % Every call would return M_in + one change, silently discarding all
    % previous updates.  Call set_dynare_param(M_loc, ...) directly so
    % that each call sees the already-updated M_loc.

    % -- parameters that do NOT require SS recomputation --
    M_loc = set_dynare_param(M_loc, 'ilabcosts', ilabcosts);
    M_loc = set_dynare_param(M_loc, 'kappaV',    kappaV);
    M_loc = set_dynare_param(M_loc, 'rho_om1',   rho_om);
    M_loc = set_dynare_param(M_loc, 'sigma_om',  sigma_om);
    M_loc = set_dynare_param(M_loc, 'rho_tfp1',  rho_A);
    M_loc = set_dynare_param(M_loc, 'rho_pvstar',   rho_pvstar);
    M_loc = set_dynare_param(M_loc, 'sigma_pvstar', sigma_pvstar);
    M_loc = set_dynare_param(M_loc, 'rho_xi',       rho_xi);
    M_loc = set_dynare_param(M_loc, 'sigma_xi',     sigma_xi);
    % Activate the PVstar shock (unit variance innovation; size controlled by sigma_pvstar)
    M_loc.Sigma_e(4, 4) = 1.0;
    % Activate the preference/demand shock (eps_xi is the 17th exogenous variable
    % in the varexo ordering: eps_om, eps_i, epschi, eps_pvstar, epsA_1..12, eps_xi)
    M_loc.Sigma_e(17, 17) = 1.0;

    % sectoral labour adjustment cost cl_i = 1/ilabcosts  (cl_i>0 → friction)
    % In the mod cl_@{i} is the *level* parameter; ilabcosts scales the cost in the
    % economy-wide FOC for labour.  We keep the cl_i profiles fixed and only
    % move the aggregate scaling ilabcosts.
    for i = 1:nsec
        M_loc = set_dynare_param(M_loc, sprintf('isigma_tfp_%d',i), isigma_tfp(i));
    end

    % -- parameters that require SS recomputation (epsY, epsM) --
    need_ss_update = (abs(epsY - get_param(M_in,'epsY_1')) > 1e-8) || ...
                     (abs(epsM - get_param(M_in,'epsM_1')) > 1e-8);

    for i = 1:nsec
        M_loc = set_dynare_param(M_loc, sprintf('epsY_%d',i), epsY);
        M_loc = set_dynare_param(M_loc, sprintf('epsM_%d',i), epsM);
    end

    % ------------------------------------------------------------------ %
    %  4. Recompute steady state if epsY or epsM changed                 %
    % ------------------------------------------------------------------ %
    if need_ss_update && ~isempty(baseline.steady_fn)
        [oo_loc, M_loc] = recompute_ss(epsY, epsM, M_loc, oo_loc, baseline);
        if isempty(oo_loc)
            moments  = NaN(46,1);
            info_flag = 98;
            var_names_out = get_moment_names();
            return
        end
    end

    % ------------------------------------------------------------------ %
    %  5. Re-solve linearised model (no recompilation)                   %
    %  Dynare 6.2 signature:                                             %
    %    resol(check_flag, M_, options_, dr_in,                         %
    %          endo_steady_state, exo_steady_state, exo_det_ss)         %
    %  Returns: [dr, info, params]                                       %
    % ------------------------------------------------------------------ %
    try
        endo_ss = oo_loc.dr.ys;
        exo_ss  = zeros(M_loc.exo_nbr, 1);
        if isfield(M_loc, 'exo_det_nbr') && M_loc.exo_det_nbr > 0
            exo_det_ss = zeros(M_loc.exo_det_nbr, 1);
        else
            exo_det_ss = zeros(0, 1);
        end
        [dr_new, info_flag, ~] = resol(0, M_loc, options_loc, oo_loc.dr, endo_ss, exo_ss, exo_det_ss);
    catch ME
        persistent resol_warned;
        if isempty(resol_warned)
            fprintf('[smm_model_moments] resol() threw: %s\n', ME.message);
            resol_warned = 1;
        end
        moments  = NaN(46,1);
        info_flag = 97;
        var_names_out = get_moment_names();
        return
    end

    if info_flag(1) ~= 0
        persistent bk_warned;
        if isempty(bk_warned)
            fprintf('[smm_model_moments] resol() info=%d (BK/numerical failure).\n', info_flag(1));
            bk_warned = 1;
        end
        moments  = NaN(46,1);
        var_names_out = get_moment_names();
        return
    end

    % ------------------------------------------------------------------ %
    %  6. Compute theoretical variance-covariance via Lyapunov equation  %
    % ------------------------------------------------------------------ %
    T = dr_new.ghx;          % ny × n_state  (impact of lagged states)
    R = dr_new.ghu;          % ny × n_exog   (impact of innovations)
    % M_.Sigma_e holds shock variances declared in shocks block (all = 1.0)
    % The sigma scaling is already embedded in ghu via the linearisation.
    Qe = M_loc.Sigma_e;      % n_exog × n_exog (diagonal, unit variances)

    % State-variable rows in decision-rule ordering
    % dr.state_var(k) = declaration-order index of k-th predetermined variable
    % dr.order_var(i) = declaration-order index of variable at DR row i
    % So DR row for state k = find(order_var == state_var(k))
    n_state = size(T, 2);
    state_dr_rows = zeros(n_state, 1);
    for k = 1:n_state
        pos = find(dr_new.order_var == dr_new.state_var(k), 1);
        if isempty(pos)
            moments  = NaN(46,1);
            info_flag = 96;
            var_names_out = get_moment_names();
            return
        end
        state_dr_rows(k) = pos;
    end

    A = T(state_dr_rows, :);      % n_state × n_state
    B = R(state_dr_rows, :);      % n_state × n_exog

    % Solve discrete Lyapunov: P = A*P*A' + B*Qe*B'
    % dlyap_smm falls back to Smith doubling if CST not available.
    Q_lyap = B * Qe * B';
    Q_lyap = (Q_lyap + Q_lyap') / 2;   % enforce symmetry
    try
        P_state = dlyap_smm(A, Q_lyap);
    catch ME_lyap
        persistent lyap_warned;
        if isempty(lyap_warned)
            fprintf('[smm_model_moments] Lyapunov solve failed: %s\n', ME_lyap.message);
            fprintf('  A %dx%d, max|eig(A)|=%.6f\n', size(A,1), size(A,2), max(abs(eig(A))));
            lyap_warned = 1;
        end
        moments  = NaN(46,1);
        info_flag = 96;
        var_names_out = get_moment_names();
        return
    end

    % Check for numerical issues
    if any(diag(P_state) < -1e-10) || any(isnan(P_state(:))) || any(isinf(P_state(:)))
        moments  = NaN(46,1);
        info_flag = 95;
        var_names_out = get_moment_names();
        return
    end
    P_state = (P_state + P_state') / 2;  % enforce symmetry before propagating

    % Contemporaneous variance-covariance of ALL endogenous variables (in DR order)
    Gamma = T * P_state * T' + R * Qe * R';

    % Lag-1 cross-covariance: Gamma_1(i,j) = Cov(y_{t}(i), y_{t-1}(j))
    % Derivation: E[y_t y_{t-1}'] = T*A*P_state*T' + T*A*B*Qe*R'
    Gamma_1 = T * A * (P_state * T' + B * Qe * R');

    % ------------------------------------------------------------------ %
    %  7. Extract model moments                                           %
    % ------------------------------------------------------------------ %
    % Use endo_ss (the SS we passed into resol) for all normalisations,
    % NOT dr_new.ys.  Dynare's bytecode steady-state solver is called
    % inside every resol() invocation and converges to its tolerance
    % (~1e-6), introducing numerical noise in dr_new.ys.  With the
    % default FD step (~sqrt(eps) ≈ 1.5e-8) that noise is ~67× larger
    % than the signal, making FD gradients unreliable.  endo_ss is the
    % algebraically-exact SS we prepared before calling resol, so it
    % carries no solver noise and keeps the objective smooth.
    ys = endo_ss;                % steady-state vector (declaration order, length ny)
    ov = dr_new.order_var;       % DR ordering: ov(k) = index in M_.endo_names

    % -- Helper: percentage std dev of variable varname --
    % (level std / SS level; for level vars like Y, GDP, Q, etc.)
    pct_std = @(varname) get_pct_std(varname, Gamma, ov, M_loc, ys);

    % Sectoral std devs (percentage)
    std_Y  = arrayfun(@(i) pct_std(sprintf('Y_%d',i)),  1:nsec)';  % [nsec×1]
    std_PH = arrayfun(@(i) pct_std(sprintf('PH_%d',i)), 1:nsec)';
    std_L  = arrayfun(@(i) pct_std(sprintf('L_%d',i)),  1:nsec)';

    % Aggregate moments
    m_std_GDP    = pct_std('GDP');
    m_std_pi     = pct_std('pi');
    m_corr_GDPpi = get_corr('GDP', 'pi', Gamma, ov, M_loc);
    m_om_g       = baseline.ombar_val;
    m_std_Q      = pct_std('Q');

    % AR(1) autocorrelation of Q (identifies rho_pvstar)
    i_Q   = find(dr_new.order_var == find(strcmp(cellstr(M_loc.endo_names),'Q'),1), 1);
    i_GDP = find(dr_new.order_var == find(strcmp(cellstr(M_loc.endo_names),'GDP'),1), 1);
    if ~isempty(i_Q) && Gamma(i_Q, i_Q) > 0
        m_autocorr_Q = Gamma_1(i_Q, i_Q) / Gamma(i_Q, i_Q);
    else
        m_autocorr_Q = NaN;
    end
    % corr(GDP, Q) contemporaneous
    if ~isempty(i_GDP) && ~isempty(i_Q) && Gamma(i_GDP,i_GDP) > 0 && Gamma(i_Q,i_Q) > 0
        m_corr_GDPQ = Gamma(i_GDP, i_Q) / sqrt(Gamma(i_GDP,i_GDP) * Gamma(i_Q,i_Q));
    else
        m_corr_GDPQ = NaN;
    end

    % ------------------------------------------------------------------ %
    %  8. Rank correlations (cross-sectional, model vs data)             %
    %   [44] Spearman rho: rank(std_Y_model) vs rank(std_Y_data)        %
    %   [45] Spearman rho: rank(std_PH_model) vs rank(std_PH_data)      %
    %   [46] Spearman rho: rank(std_L_model)  vs rank(std_L_data)       %
    % ------------------------------------------------------------------ %
    d_std_Y  = baseline.data_std_Y;   % [nsec×1] from data
    d_std_PH = baseline.data_std_PH;
    d_std_L  = baseline.data_std_L;

    valid_y = isfinite(std_Y)  & isfinite(d_std_Y);
    valid_p = isfinite(std_PH) & isfinite(d_std_PH);
    valid_l = isfinite(std_L)  & isfinite(d_std_L);

    % safe_spearman returns 0 (not NaN) when either vector is constant.
    % This happens at theta0 when only one sector has active TFP shocks
    % (e.g. params_val_ul.mat was saved during a single-shock exercise):
    % I-O network propagation makes all model std devs nearly identical
    % → Spearman rank-corr is undefined. Returning 0 means the optimizer
    % sees a penalty of (1 - 0)^2 = 1, encouraging sector discrimination.
    if sum(valid_y) >= 3
        m_rho_Y = safe_spearman(std_Y(valid_y),  d_std_Y(valid_y));
    else
        m_rho_Y = 0;
    end
    if sum(valid_p) >= 3
        m_rho_PH = safe_spearman(std_PH(valid_p), d_std_PH(valid_p));
    else
        m_rho_PH = 0;
    end
    if sum(valid_l) >= 3
        m_rho_L = safe_spearman(std_L(valid_l),  d_std_L(valid_l));
    else
        m_rho_L = 0;
    end

    % ------------------------------------------------------------------ %
    %  9. Stack into moment vector (46 elements)                         %
    %   [1-12]  std(Y_i)    [13-24] std(PH_i)   [25-36] std(L_i)        %
    %   [37]    std(GDP)    [38]    std(pi)      [39]    corr(GDP,pi)   %
    %   [40]    mean omG    [41]    std(Q)       [42]    autocorr(Q)    %
    %   [43]    corr(GDP,Q) [44]    rho_Y        [45]    rho_PH         %
    %   [46]    rho_L                                                   %
    % ------------------------------------------------------------------ %
    moments = [
        std_Y;
        std_PH;
        std_L;
        m_std_GDP;
        m_std_pi;
        m_corr_GDPpi;
        m_om_g;
        m_std_Q;
        m_autocorr_Q;
        m_corr_GDPQ;
        m_rho_Y;
        m_rho_PH;
        m_rho_L
    ];

    var_names_out = get_moment_names();
end

% ======================================================================= %
%  Local helper functions                                                   %
% ======================================================================= %

function M_out = set_dynare_param(M_in, pname, val)
% Works with both char-array (Dynare <=5) and cell-array (Dynare 6) param_names.
    M_out = M_in;
    pn = M_in.param_names;
    if ischar(pn)
        % char array: each row is a name padded with spaces
        idx = find(strcmp(cellstr(pn), pname), 1);
    else
        idx = find(strcmp(pn, pname), 1);
    end
    if ~isempty(idx)
        M_out.params(idx) = val;
    end
end

function val = get_param(M_in, pname)
    pn = M_in.param_names;
    if ischar(pn)
        idx = find(strcmp(cellstr(pn), pname), 1);
    else
        idx = find(strcmp(pn, pname), 1);
    end
    if isempty(idx)
        val = NaN;
    else
        val = M_in.params(idx);
    end
end

function std_pct = get_pct_std(varname, Gamma, order_var, M_loc, ys)
    pos = find(strcmp(cellstr(M_loc.endo_names), varname), 1);
    if isempty(pos)
        std_pct = 0;
        return
    end
    dr_pos = find(order_var == pos, 1);
    if isempty(dr_pos)
        std_pct = 0;
        return
    end
    ss_val = abs(ys(pos));
    var_v  = Gamma(dr_pos, dr_pos);
    std_pct = sqrt(max(0, var_v)) / max(ss_val, 1e-12);
end

function c = get_corr(vname1, vname2, Gamma, order_var, M_loc)
    pos1 = find(strcmp(cellstr(M_loc.endo_names), vname1), 1);
    pos2 = find(strcmp(cellstr(M_loc.endo_names), vname2), 1);
    dr1  = find(order_var == pos1, 1);
    dr2  = find(order_var == pos2, 1);
    cov12 = Gamma(dr1, dr2);
    var1  = Gamma(dr1, dr1);
    var2  = Gamma(dr2, dr2);
    denom = sqrt(max(var1, 0) * max(var2, 0));
    if denom < 1e-15
        c = 0;
    else
        c = cov12 / denom;
        c = max(-1, min(1, c));
    end
end

function names = get_moment_names()
    n = cell(46,1);
    for ii = 1:12, n{ii}    = sprintf('std(Y_%d)',  ii); end
    for ii = 1:12, n{12+ii} = sprintf('std(PH_%d)', ii); end
    for ii = 1:12, n{24+ii} = sprintf('std(L_%d)',  ii); end
    n{37} = 'std(GDP)';
    n{38} = 'std(pi)';
    n{39} = 'corr(GDP, pi)';
    n{40} = 'mean goods expenditure share';
    n{41} = 'std(Q)';
    n{42} = 'autocorr(Q)';
    n{43} = 'corr(GDP, Q)';
    n{44} = 'rank corr: output (model vs data)';
    n{45} = 'rank corr: prices (model vs data)';
    n{46} = 'rank corr: labor  (model vs data)';
    names = n;
end

function [oo_out, M_out] = recompute_ss(epsY, epsM, M_loc, oo_loc, baseline)
% Recompute the FULL steady state when epsY or epsM change, and update BOTH
% oo_loc.steady_state (endogenous SS) and M_loc.params (Dynare parameters
% that enter the model equations and affect the linearisation).
%
% The .mod file loads params_val_ul.mat at compile time and embeds several
% SS scalars as Dynare parameters:
%   Global : Ctot_ss Ctotg_ss Ctots_ss VA_ss M_tot_ss Y_ss IMP_ss bbar
%   Sector : PL_ss_i  (= w_ss, one per sector, 12 total)
% These MUST be updated in M_.params before calling resol(); otherwise the
% decision rule is linearised around the wrong steady state.
%
% Returns [] (empty) on failure so the caller can issue info=98.
    oo_out = [];
    M_out  = [];
    nsec   = baseline.nsec;

    % ---- unpack calibration objects (fixed across SMM iterations) --------
    modepsY    = ones(nsec,1) * epsY;
    modepsM    = ones(nsec,1) * epsM;
    alpha_vec  = baseline.modalpha;
    alphaV_vec = baseline.modalphaV;
    beta_mat   = baseline.modbeta;
    varrho_val = baseline.modvarrho;
    gammag_vec = baseline.modgammag;
    gammas_vec = baseline.modgammas;
    chiX_vec   = baseline.modchiX;
    A_vec      = ones(nsec,1);
    om_g       = baseline.ombar_val;
    om_s       = 1 - om_g;
    PVstar_ss  = baseline.PVstar_ss;
    Pistar_ss  = 1.0;           % always 1 in this model
    epsilon    = baseline.epsilon_val;
    sigmaH     = baseline.sigmaH_val;
    etastar    = baseline.etastar_val;
    omegaX     = baseline.omegaX_val;
    Ystar      = baseline.ystar_ss_val;
    gamma_val  = baseline.gamma_val;
    chi_val    = baseline.chi_val;
    psi_val    = baseline.psi_val;
    beta_val   = baseline.beta_val;
    tb_target  = baseline.tb_target;
    chii_b     = get_param(M_loc, 'chii_b');

    % ---- Step 1: solve for pH_ss, w_ss, Q_ss, C_ss ----------------------
    en_cell   = cellstr(M_loc.endo_names);
    C_ss_prev = oo_loc.steady_state(find(strcmp(en_cell,'C'),1));
    w_ss_prev = oo_loc.steady_state(find(strcmp(en_cell,'w'),1));
    Q_ss_prev = oo_loc.steady_state(find(strcmp(en_cell,'Q'),1));
    pH_guess  = arrayfun(@(i) oo_loc.steady_state(find(strcmp(en_cell,...
                    sprintf('PH_%d',i)),1)), 1:nsec)';

    x0   = [pH_guess; w_ss_prev; Q_ss_prev; C_ss_prev];
    opts = optimoptions('fsolve','TolFun',1e-12,'Display','off','MaxIterations',500);

    try
        [ss_sol, fval, exitflag] = fsolve(@(x) baseline.steady_fn(x, PVstar_ss, epsilon,...
            varrho_val, sigmaH, gammag_vec, gammas_vec, om_g, om_s, chiX_vec, omegaX,...
            etastar, Ystar, alpha_vec, alphaV_vec, beta_mat, modepsY, modepsM,...
            gamma_val, chi_val, psi_val, A_vec, tb_target), x0, opts);
    catch
        return
    end
    if exitflag <= 0 || norm(fval) > 1e-6
        return
    end

    pH_ss = ss_sol(1:nsec);
    w_ss  = ss_sol(nsec+1);
    Q_ss  = ss_sol(nsec+2);
    C_ss  = ss_sol(nsec+3);

    % ---- Step 2: derive all remaining SS quantities ----------------------
    PL_ss  = ones(nsec,1) * w_ss;
    PV_ss  = Q_ss * PVstar_ss;
    MCi_ss = (epsilon-1)/epsilon * pH_ss;

    PMi_ss = zeros(nsec,1);
    for ii = 1:nsec
        PMi_ss(ii) = (sum(beta_mat(ii,:) .* pH_ss'.^(1-modepsM(ii))))^(1/(1-modepsM(ii)));
    end

    P_ss   = (varrho_val.^sigmaH .* pH_ss.^(1-sigmaH) + ...
              (1-varrho_val).^sigmaH .* PV_ss.^(1-sigmaH)).^(1/(1-sigmaH));
    p_g_ss = prod(P_ss.^gammag_vec);
    p_s_ss = prod(P_ss.^gammas_vec);

    C_g_ss  = om_g * C_ss / p_g_ss;
    C_s_ss  = om_s * C_ss / p_s_ss;
    C_gi_ss = gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss;
    C_si_ss = gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss;

    CHg_ss = varrho_val.^sigmaH .* (pH_ss./P_ss).^(-sigmaH) .* C_gi_ss;
    CHs_ss = varrho_val.^sigmaH .* (pH_ss./P_ss).^(-sigmaH) .* C_si_ss;
    CFg_ss = (1-varrho_val).^sigmaH .* (PV_ss./P_ss).^(-sigmaH) .* C_gi_ss;
    CFs_ss = (1-varrho_val).^sigmaH .* (PV_ss./P_ss).^(-sigmaH) .* C_si_ss;
    CHi_ss = CHg_ss + CHs_ss;
    CFi_ss = CFg_ss + CFs_ss;

    PX_ss = prod(pH_ss.^chiX_vec);
    X_ss  = omegaX * (PX_ss/Q_ss)^(-etastar) * Ystar;
    Xi_ss = chiX_vec .* X_ss .* PX_ss ./ pH_ss;

    % ---- Step 3: solve production system for M, L, Vi, Yi ---------------
    % Initial guess from current SS
    M_prev  = arrayfun(@(i) oo_loc.steady_state(find(strcmp(en_cell,sprintf('M_%d',i)),1)),  1:nsec)';
    L_prev  = arrayfun(@(i) oo_loc.steady_state(find(strcmp(en_cell,sprintf('L_%d',i)),1)),  1:nsec)';
    Vi_prev = arrayfun(@(i) oo_loc.steady_state(find(strcmp(en_cell,sprintf('V_%d',i)),1)),  1:nsec)';
    Yi_prev = arrayfun(@(i) oo_loc.steady_state(find(strcmp(en_cell,sprintf('Y_%d',i)),1)),  1:nsec)';
    % Replace any zeros/NaN in guesses with analytical approximation
    for ii = 1:nsec
        if ~isfinite(M_prev(ii))  || M_prev(ii)  <= 0, M_prev(ii)  = 0.1; end
        if ~isfinite(L_prev(ii))  || L_prev(ii)  <= 0, L_prev(ii)  = 0.1; end
        if ~isfinite(Vi_prev(ii)) || Vi_prev(ii) <= 0, Vi_prev(ii) = 0.01; end
        if ~isfinite(Yi_prev(ii)) || Yi_prev(ii) <= 0, Yi_prev(ii) = 0.1; end
    end
    x0_sys = [M_prev; L_prev; Vi_prev; Yi_prev];

    has_sys = isfield(baseline,'steady_sys_fn') && ~isempty(baseline.steady_sys_fn);
    if has_sys
        opts_sys = optimoptions('fsolve','TolFun',1e-10,'Display','off','MaxIterations',500);
        try
            [x_sol, fval2, ef2] = fsolve(@(x) baseline.steady_sys_fn(x, alpha_vec, alphaV_vec,...
                beta_mat, MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss,...
                modepsY, modepsM, A_vec, pH_ss), x0_sys, opts_sys);
        catch
            x_sol = x0_sys;  ef2 = -1;
        end
        if ef2 > 0 && norm(fval2) < 1e-6
            M_ss  = x_sol(1:nsec);
            L_ss  = x_sol(nsec+1:2*nsec);
            Vi_ss = x_sol(2*nsec+1:3*nsec);
            Yi_ss = x_sol(3*nsec+1:4*nsec);
        else
            % Fallback: analytical first-order approximation (no intermediates)
            M_ss  = (MCi_ss./PMi_ss).^modepsY .* alpha_vec    .* (CHi_ss + Xi_ss);
            L_ss  = (MCi_ss./PL_ss ).^modepsY .* (1-alpha_vec-alphaV_vec) .* (CHi_ss + Xi_ss);
            Vi_ss = (MCi_ss./PV_ss  ).^modepsY .* alphaV_vec  .* (CHi_ss + Xi_ss);
            Yi_ss = A_vec .* (alpha_vec.^(1./modepsY).*M_ss.^((modepsY-1)./modepsY) + ...
                              alphaV_vec.^(1./modepsY).*Vi_ss.^((modepsY-1)./modepsY) + ...
                              (1-alphaV_vec-alpha_vec).^(1./modepsY).*L_ss.^((modepsY-1)./modepsY)).^(modepsY./(modepsY-1));
        end
    else
        % No production-system solver available — use analytical approximation
        M_ss  = (MCi_ss./PMi_ss).^modepsY .* alpha_vec    .* (CHi_ss + Xi_ss);
        L_ss  = (MCi_ss./PL_ss ).^modepsY .* (1-alpha_vec-alphaV_vec) .* (CHi_ss + Xi_ss);
        Vi_ss = (MCi_ss./PV_ss  ).^modepsY .* alphaV_vec  .* (CHi_ss + Xi_ss);
        Yi_ss = A_vec .* (alpha_vec.^(1./modepsY).*M_ss.^((modepsY-1)./modepsY) + ...
                          alphaV_vec.^(1./modepsY).*Vi_ss.^((modepsY-1)./modepsY) + ...
                          (1-alphaV_vec-alpha_vec).^(1./modepsY).*L_ss.^((modepsY-1)./modepsY)).^(modepsY./(modepsY-1));
    end

    % ---- Step 4: aggregate quantities ------------------------------------
    IMP_tot_ss = sum(Vi_ss) + sum(CFi_ss);
    TB_ss      = PX_ss*X_ss - PV_ss*IMP_tot_ss;
    GDP_ss     = C_ss + TB_ss;
    N_ss       = sum(L_ss);
    Y_tot_ss   = sum(Yi_ss);
    M_tot_ss   = sum(M_ss);
    VA_ss      = sum(Yi_ss - M_ss);
    Ctotg_ss   = sum(gammag_vec .* (p_g_ss./P_ss) .* C_g_ss);
    Ctots_ss   = sum(gammas_vec .* (p_s_ss./P_ss) .* C_s_ss);
    Ctot_ss    = Ctotg_ss + Ctots_ss;

    r_star_ss  = Pistar_ss / beta_val;
    Bstar_ss   = -TB_ss / (Q_ss * (1 - r_star_ss/Pistar_ss));
    bbar_new   = Q_ss * Bstar_ss / GDP_ss;

    % ---- Step 5: update M_loc.params (Dynare parameters used in model) --
    M_out = M_loc;
    % PL_ss_i = w_ss for every sector
    for i = 1:nsec
        M_out = set_dynare_param(M_out, sprintf('PL_ss%d',i), w_ss);
    end
    M_out = set_dynare_param(M_out, 'Ctot_ss',  Ctot_ss);
    M_out = set_dynare_param(M_out, 'Ctotg_ss', Ctotg_ss);
    M_out = set_dynare_param(M_out, 'Ctots_ss', Ctots_ss);
    M_out = set_dynare_param(M_out, 'VA_ss',    VA_ss);
    M_out = set_dynare_param(M_out, 'M_tot_ss', M_tot_ss);
    M_out = set_dynare_param(M_out, 'Y_ss',     Y_tot_ss);
    M_out = set_dynare_param(M_out, 'IMP_ss',   IMP_tot_ss);
    M_out = set_dynare_param(M_out, 'bbar',     bbar_new);

    % ---- Step 6: update oo_loc.steady_state for endogenous variables -----
    oo_out = oo_loc;
    upd = @(nm, v) update_ss_var(oo_out, M_out, nm, v);

    oo_out = upd('w',   w_ss);
    oo_out = upd('Q',   Q_ss);
    oo_out = upd('C',   C_ss);
    oo_out = upd('N',   N_ss);
    oo_out = upd('GDP', GDP_ss);
    oo_out = upd('TB',  TB_ss);
    for i = 1:nsec
        oo_out = upd(sprintf('PH_%d',i), pH_ss(i));
        oo_out = upd(sprintf('MC_%d',i), MCi_ss(i));
        oo_out = upd(sprintf('PM_%d',i), PMi_ss(i));
        oo_out = upd(sprintf('PL_%d',i), PL_ss(i));
        oo_out = upd(sprintf('P_%d', i), P_ss(i));
        oo_out = upd(sprintf('Y_%d', i), Yi_ss(i));
        oo_out = upd(sprintf('L_%d', i), L_ss(i));
        oo_out = upd(sprintf('M_%d', i), M_ss(i));
        oo_out = upd(sprintf('V_%d', i), Vi_ss(i));
    end
end

function oo_out = update_ss_var(oo_in, M_loc, varname, val)
    oo_out = oo_in;
    idx = find(strcmp(cellstr(M_loc.endo_names), varname), 1);
    if ~isempty(idx)
        oo_out.steady_state(idx) = val;
        % Also update dr.ys if it exists
        if isfield(oo_out, 'dr') && isfield(oo_out.dr, 'ys')
            oo_out.dr.ys(idx) = val;
        end
    end
end

function P = dlyap_smm(A, Q)
% Solve discrete Lyapunov equation  P = A*P*A' + Q  for stable A.
% Uses MATLAB's dlyap (Control System Toolbox) if available; otherwise
% falls back to the Smith doubling algorithm which converges in O(log(1/eps))
% matrix multiplications and requires only max|eig(A)| < 1.
    try
        P = dlyap(A, Q);
        return
    catch
    end
    % Smith doubling: accumulate S = sum_{k>=0} A^k * Q * (A^k)'
    % which equals the solution P = A*P*A' + Q.
    P  = Q;
    Ak = A;
    for k = 1:200
        P_new = P + Ak * P * Ak';
        Ak    = Ak * Ak;
        if norm(P_new - P, 'fro') < 1e-13 * (norm(P_new, 'fro') + 1e-30)
            P = P_new;
            break
        end
        P = P_new;
        if max(abs(Ak(:))) < 1e-15
            break
        end
    end
    P = (P + P') / 2;  % enforce symmetry
end

function rc = safe_spearman(a, b)
% Spearman rank-correlation of a and b.
% Returns 0 (=worst fit, target is 1) when either vector is constant
% so the optimizer penalises degenerate solutions without crashing.
    if std(double(a)) < 1e-12 || std(double(b)) < 1e-12
        rc = 0;
    else
        rc = corr(double(a), double(b), 'Type', 'Spearman');
    end
end

function rc = safe_rank_corr(data_vec, model_vec)
% 1 - Spearman rank-corr(data_vec, model_vec).
% Returns 0 (no penalty) when either vector is constant, which happens when
% data_moments_chile.mat is missing and the XLS fallback loaded all-zero data.
    if std(double(data_vec)) < 1e-12 || std(double(model_vec)) < 1e-12
        rc = 0;
    else
        rc = 1 - corr(data_vec, model_vec, 'type', 'Spearman');
    end
end
