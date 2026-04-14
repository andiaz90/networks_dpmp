% smm_estimation.m
%
% Simulated Method of Moments (SMM) estimation for the NK-IOSOE Chile model.
%
% STRATEGY
%   The model is solved at first order by Dynare (NK_SOE_lev_gap2.mod).
%   For each candidate parameter vector theta, Dynare's internal `resol`
%   function re-solves the linearised system WITHOUT recompiling the model.
%   Theoretical unconditional second moments are then derived analytically
%   from the state-space representation via the discrete Lyapunov equation,
%   matching the 42 target moments described in the paper.
%
% PREREQUISITES
%   1. Run main_SOE_gap.m at least ONCE.  This compiles NK_SOE_lev_gap2.mod and
%      leaves M_, oo_, options_ in the base workspace plus params_val_ul.mat.
%      After this one-time compilation, SMM re-solves with resol() only.
%   2. Populate the data-moment placeholders in Section 2 below.
%
% ESTIMATED PARAMETERS  theta (21-element vector)
%   [1]     ilabcosts                aggregate labour adjustment cost (inverse)
%   [2]     epsY                     elast. of subst. in production
%   [3]     epsM                     elast. of subst. between materials
%   [4]     log(kappaV)              log of import price adj. cost
%   [5]     rho_om                   AR persistence, goods-services shock
%   [6]     sigma_om                 std dev, goods-services shock
%   [7]     rho_A                    AR persistence, sectoral TFP shocks (common)
%   [8-19]  isigma_tfp_1,...,12      std dev of each sector's TFP shock
%   [20]    rho_pvstar               AR persistence, import price (PVstar) shock
%   [21]    sigma_pvstar             std dev of import price shock

clear all 
close all
clc

fprintf('\n============================================================\n');
fprintf('  SMM ESTIMATION: NK-IOSOE Chile Model\n');
fprintf('============================================================\n\n');

%% ================================================================== %%
%%  0. CHECK PREREQUISITES  (auto-run main_SOE if needed)              %%
%% ================================================================== %%
if ~exist('M_','var') || ~exist('oo_','var') || ~exist('options_','var') || ...
        ~isfield(oo_,'dr') || ~isfield(oo_.dr,'ghx')

    fprintf('Dynare workspace not found. Running main_SOE_gap.m (EXERCISE=0: Baseline)...\n');
    fprintf('This compiles NK_SOE_lev_gap2.mod and may take ~1 minute.\n\n');

    % Use env var to pass EXERCISE=0 across clear all in main_SOE_gap.m
    setenv('SMM_EXERCISE', '0');
    main_soe_path = fullfile(fileparts(mfilename('fullpath')), 'main_SOE_gap.m');
    try
        run(main_soe_path);
    catch ME_setup
        setenv('SMM_EXERCISE', '');
        error('main_SOE_gap.m failed: %s\nFix the model setup before running SMM.', ME_setup.message);
    end
    setenv('SMM_EXERCISE', '');   % clear the override flag

    % Verify Dynare variables are now present
    if ~exist('M_','var') || ~exist('oo_','var') || ~exist('options_','var')
        error('main_SOE_gap.m ran but M_/oo_/options_ still not in workspace.');
    end
    if ~isfield(oo_,'dr') || ~isfield(oo_.dr,'ghx')
        error('oo_.dr not initialised after main_SOE_gap.m. Check stoch_simul output.');
    end
    fprintf('\nModel initialised. Proceeding with SMM estimation...\n\n');
end

nsec     = 12;
goods    = logical([1;1;1;1;1;0;0;0;0;0;0;0]);   % sectors 1-5
services = logical([0;0;0;0;0;1;1;1;1;1;1;1]);   % sectors 6-12

%% ================================================================== %%
%%  1. LOAD SECTORAL DATA MOMENTS                                       %%
%%     Priority: data_moments_chile.mat (from compute_data_moments.m)  %%
%%     Fallback: Stata_to_excel_few_industries_chile.xls               %%
%% ================================================================== %%
fprintf('Loading sectoral data moments...\n');

data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';

dm_file = fullfile(fileparts(mfilename('fullpath')), 'data_moments_chile.mat');
if exist(dm_file, 'file')
    fprintf('  Loading from data_moments_chile.mat\n');
    tmp   = load(dm_file);
    dm    = tmp.dm_chile;
    y_d   = dm.y_d;       % [12x1] sectoral output std devs
    p_d   = dm.p_d;       % [12x1] sectoral price  std devs
    l_d   = dm.l_d;       % [12x1] sectoral employ std devs (from count_workers_by_sector.csv)
    d_std_GDP    = dm.d_std_GDP;
    d_std_pi     = dm.d_std_pi;
    d_corr_GDPpi = dm.d_corr_GDPpi;
    d_omG        = dm.d_omG;
    d_std_Q      = dm.d_std_Q;
    d_autocorr_Q = dm.d_autocorr_Q;
    d_corr_GDPQ  = dm.d_corr_GDPQ;
    d_TBGDP      = dm.d_TBGDP;
    fprintf('  Sample: %dQ%d - %dQ%d\n', dm.sample_start(1), dm.sample_start(2), ...
            dm.sample_end(1), dm.sample_end(2));
else
    fprintf('  data_moments_chile.mat not found. Falling back to Stata XLS.\n');
    fprintf('  Run compute_data_moments.m first for fully data-based moments.\n');
    filename_industries = fullfile(data_path, 'Stata_to_excel_few_industries_chile.xls');
    varNames = {'kk','BEAName','Nameshort','dp','dy','dl','dri','share', ...
                'industrytype','spend_good','spend_serv','alpha','alpha_V','var_rho'};
    varTypes = {'double','char','char','double','double','double','double', ...
                'double','categorical','double','double','double','double','double'};
    opts_xls = spreadsheetImportOptions('NumVariables',14,'VariableNames',varNames, ...
               'VariableTypes',varTypes,'DataRange','A2');
    alldata  = readtable(filename_industries, opts_xls);
    y_d = table2array(alldata(1:nsec, 5));
    p_d = table2array(alldata(1:nsec, 4));
    l_d = table2array(alldata(1:nsec, 6));

    %% Aggregate moment placeholders (fill from BCCh/IFS data)
    d_std_GDP    = 0.0215;
    d_std_pi     = 0.0041;
    d_corr_GDPpi = -0.15;
    d_omG        = 0.57;
    d_std_Q      = 0.0520;
    d_autocorr_Q = 0.75;   % AR(1) of HP-filtered log REER (lit. estimate for Chile)
    d_corr_GDPQ  = -0.15;  % corr(GDP, Q): typically negative in SOE models
    d_TBGDP      = -0.02;
end

%% ================================================================== %%
%%  3. BUILD DATA MOMENT VECTOR AND BASELINE PARAMETER STRUCTURE       %%
%% ================================================================== %%

% Steady-state output weights for cross-sectional averages
get_ss   = @(nm) oo_.dr.ys(find(strcmp(cellstr(M_.endo_names), nm), 1));
Y_ss_vec = arrayfun(@(i) get_ss(sprintf('Y_%d',i)), 1:nsec)';
w_g = Y_ss_vec(goods)   / sum(Y_ss_vec(goods));
w_s = Y_ss_vec(services)/ sum(Y_ss_vec(services));

% Full 46-element data moment vector.
%   [1-12]  std(Y_i)     sectoral output volatility (identifies isigma_tfp_i)
%   [13-24] std(PH_i)    sectoral price  volatility (over-id: tests IO + stickiness)
%   [25-36] std(L_i)     sectoral labor  volatility (over-id: tests IO + ilabcosts)
%   [37]    std(GDP)     [38] std(pi)    [39] corr(GDP,pi)
%   [40]    mean omG     [41] std(Q)
%   [42]    autocorr(Q)  [43] corr(GDP,Q)
%   [44]    rank corr output (model vs data)   target = 1
%   [45]    rank corr prices (model vs data)   target = 1
%   [46]    rank corr labor  (model vs data)   target = 1
data_moments = [
    y_d;                        % 1-12   std(Y_i)
    p_d;                        % 13-24  std(PH_i)
    l_d;                        % 25-36  std(L_i)
    d_std_GDP;                  % 37
    d_std_pi;                   % 38
    d_corr_GDPpi;               % 39
    d_omG;                      % 40     SS target (passive)
    d_std_Q;                    % 41
    d_autocorr_Q;               % 42     identifies rho_pvstar
    d_corr_GDPQ;                % 43     over-identifies pvstar shock
    1.0;                        % 44     rank corr output: target = perfect
    1.0;                        % 45     rank corr prices: target = perfect
    1.0                         % 46     rank corr labor:  target = perfect
];
assert(numel(data_moments) == 46, ...
    'BUG: expected 46 data moments, got %d. Check that d_autocorr_Q and d_corr_GDPQ are loaded.', ...
    numel(data_moments));

% -------- baseline: fixed objects passed into the moment function --------
load params_val_ul.mat;

baseline = struct();
baseline.nsec        = nsec;
baseline.goods       = goods;
baseline.services    = services;
baseline.Y_ss        = Y_ss_vec;
baseline.GDP_ss      = get_ss('GDP');
baseline.ombar_val   = ombar_val;
baseline.tb_target   = d_TBGDP;   % from BCCh CCNN data (was hardcoded -0.02)
baseline.modalpha    = modalpha;
baseline.modalphaV   = modalphaV;
baseline.modbeta     = modbeta;
baseline.modgammag   = modgammag;
baseline.modgammas   = modgammas;
baseline.modvarrho   = modvarrho;
baseline.modchiX     = modchiX;
baseline.modkappa    = modkappa;
baseline.gamma_val   = gamma_val;
baseline.psi_val     = 1;
baseline.chi_val     = 1;
baseline.epsilon_val = 10;
baseline.beta_val    = beta_val;
baseline.PVstar_ss   = PVstar_ss;
baseline.sigmaH_val  = sigmaH_val;
baseline.etastar_val = etastar_val;
baseline.omegaX_val  = omegaX_val;
baseline.ystar_ss_val= ystar_ss_val;

if exist('steady_ntwsoe','file')
    baseline.steady_fn = @steady_ntwsoe;
else
    warning('steady_ntwsoe not on path -- SS will not update when epsY/epsM change.');
    baseline.steady_fn = [];
end

if exist('steady_ntwsoe_system','file')
    baseline.steady_sys_fn = @steady_ntwsoe_system;
else
    warning('steady_ntwsoe_system not on path -- production system SS will use analytical approximation.');
    baseline.steady_sys_fn = [];
end

baseline.data_std_Y  = y_d;
baseline.data_std_PH = p_d;
baseline.data_std_L  = l_d;

fprintf('Data moments:\n');
fprintf('  std(Y):  '); fprintf('%.4f ', y_d'); fprintf('\n');
fprintf('  std(PH): '); fprintf('%.4f ', p_d'); fprintf('\n');
fprintf('  std(L):  '); fprintf('%.4f ', l_d'); fprintf('\n');
fprintf('  std(GDP)=%.4f  std(pi)=%.4f  corr(GDP,pi)=%.3f\n', ...
        d_std_GDP, d_std_pi, d_corr_GDPpi);
fprintf('  std(Q)=%.4f  autocorr(Q)=%.3f  corr(GDP,Q)=%.3f\n', ...
        d_std_Q, d_autocorr_Q, d_corr_GDPQ);
fprintf('  TB/GDP SS target=%.3f  omG=%.2f\n\n', d_TBGDP, d_omG);

%% ================================================================== %%
%%  4. PARAMETER BOUNDS AND INITIAL VALUES                             %%
%% ================================================================== %%
theta0 = [
    ilabcosts_val;          % 1   ilabcosts
    modepsY(1);             % 2   epsY
    modepsM(1);             % 3   epsM
    log(kappaV_val);        % 4   log(kappaV)
    rho_om1_val;            % 5   rho_om
    sigma_om_val;           % 6   sigma_om
    rho_tfp1_val;           % 7   rho_A
    isigma_tfp_val(:);      % 8-19  isigma_tfp_1,...,12
    rho_pvstar_val;         % 20  rho_pvstar
    sigma_pvstar_val;       % 21  sigma_pvstar
    rho_xi_val;             % 22  rho_xi   (preference/demand shock persistence)
    sigma_xi_val            % 23  sigma_xi (preference/demand shock std dev)
];

% Fix 2: enforce rho_pvstar >= 0.50 (foreign price shocks are persistent) and
%         sigma_pvstar >= 0.005 (keep import price channel economically active).
%         Previously both collapsed to near-zero, shutting off cost-push inflation.
lb = [1e-3; 0.10; 0.01; log(1e3); -0.99; 1e-5; -0.99; ...
      1e-5*ones(nsec,1); 0.50; 0.005; ...   % rho_pvstar, sigma_pvstar
      0.00; 0];                              % rho_xi, sigma_xi
ub = [100;  3.00; 1.50; log(1e16); 0.99; 0.50;  0.99; ...
      0.50*ones(nsec,1);  0.99; 0.50; ...   % rho_pvstar, sigma_pvstar
      0.99; 0.20];                          % rho_xi, sigma_xi

param_labels = [{'ilabcosts','epsY','epsM','log(kappaV)', ...
                 'rho_om','sigma_om','rho_A'}, ...
                 arrayfun(@(i) sprintf('isigma_tfp_%d',i), 1:nsec, 'UniformOutput',false), ...
                 {'rho_pvstar','sigma_pvstar','rho_xi','sigma_xi'}];
n_theta = numel(theta0);

%% ================================================================== %%
%%  5. WEIGHTING MATRIX                                                %%
%%                                                                     %%
%%  Diagonal relative-scale weights: W_ii = 1 / max(|d_i|, floor)^2  %%
%%  This makes every moment contribute equally in percentage terms,   %%
%%  so large-magnitude moments (e.g. std(Y_1)=0.40) do not swamp     %%
%%  small ones (e.g. std(L_12)=0.05). The floor (0.01) prevents      %%
%%  near-zero moments from getting extreme weights.                   %%
%%                                                                     %%
%  Moments that are pinned by the SS calibration (omG) and the       %%
%  corr(GDP,pi) structural mismatch are down-weighted so they do not %%
%  dominate the objective with unfixable residuals. TB/GDP is not    %%
%  included in the SMM moment vector because it is imposed           %%
%  exogenously as the steady-state closure target.                   %%
%%                                                                     %%
%%  Fix 1: corr(GDP,pi) reduced from 0.10x to 0.02x.  The data value %%
%%  is -0.02 (near-zero), giving a baseline weight 1/0.02^2 = 2490.  %%
%%  Even at 0.10x this moment contributed 88% of total loss (79/90), %%
%%  pinning the optimizer to a local minimum and preventing it from   %%
%%  improving inflation and exchange-rate moments.  At 0.02x its      %%
%%  contribution falls to ~16, proportionate with other moments.      %%
%%                                                                     %%
%%  Fix 3: small-magnitude moments (|d_i| < 0.02) are capped at 0.25x%%
%%  their natural weight.  The floor_w=0.01 is insufficient for price %%
%%  moments like std(PH_10)=0.0166 which get weight 1/0.0166^2=3628, %%
%%  amplifying tiny absolute errors into dominant loss terms.         %%
%% ================================================================== %%
floor_w = 0.01;
w_diag  = 1 ./ max(abs(data_moments), floor_w).^2;
W = diag(w_diag);

% Fix 3: cap weights for small-magnitude moments to prevent mechanical dominance
small_moments_idx = abs(data_moments) < 0.02;
W(small_moments_idx, small_moments_idx) = W(small_moments_idx, small_moments_idx) * 0.25;

% Fix 1: further reduce corr(GDP,pi) — was 0.10x, now 0.02x
W(39,39) = W(39,39) * 0.02;   % corr(GDP,pi): supply-only model structurally cannot match
W(40,40) = W(40,40) * 0.10;   % mean omG: SS target, passive
W(43,43) = W(43,43) * 0.50;   % corr(GDP,Q): over-id, mild downweight
% Rank correlation moments: downweight slightly (bounded [-1,1]; typical residuals ~0.2-0.5)
W(44,44) = W(44,44) * 0.50;   % rank corr output
W(45,45) = W(45,45) * 0.50;   % rank corr prices
W(46,46) = W(46,46) * 0.50;   % rank corr labor

fprintf('Weighting matrix notes: omG downweighted; TB/GDP removed from SMM moments because it is fixed by the SS closure.\n');

%% ================================================================== %%
%%  6. OPTIMISATION                                                     %%
%%  smm_obj and smm_obj_full are local functions at end of this file.  %%
%% ================================================================== %%
obj_fn = @(t) smm_obj(t, data_moments, W, baseline, M_, options_, oo_);

%% ================================================================== %%
%%  5b. PRE-FLIGHT DIAGNOSTIC — errors out if anything is broken      %%
%%      so the optimizer never runs with a broken configuration.       %%
%% ================================================================== %%
fprintf('=== PRE-FLIGHT CHECK ===\n');

% Reset persistent variables inside smm_model_moments so all messages show
clear smm_model_moments

% ---- (A) Parameter names ---------------------------------------------------
pn = M_.param_names;
if ischar(pn), pn_cell = cellstr(pn); else, pn_cell = pn(:); end
fprintf('  M_.param_names: %s (%d entries)\n', ...
    class(M_.param_names), numel(pn_cell));

smm_pcheck = {'ilabcosts','kappaV','rho_om1','sigma_om','rho_tfp1', ...
              'isigma_tfp_1','epsY_1','epsM_1','rho_pvstar','sigma_pvstar'};
pf_missing_params = {};
for kp = 1:numel(smm_pcheck)
    idx_p = find(strcmp(pn_cell, smm_pcheck{kp}), 1);
    if isempty(idx_p)
        fprintf('  [MISSING PARAM] %s\n', smm_pcheck{kp});
        pf_missing_params{end+1} = smm_pcheck{kp}; %#ok<AGROW>
    else
        fprintf('  [param OK] %-20s = %.6g\n', smm_pcheck{kp}, M_.params(idx_p));
    end
end

% ---- (B) Endogenous variable names ----------------------------------------
en_cell = cellstr(M_.endo_names);
fprintf('  M_.endo_names: %d variables. Checking required names...\n', numel(en_cell));
nsec_pf = baseline.nsec;
pf_missing_vars = {};
for kv = 1:nsec_pf
    for vp = {'Y','PH','L'}
        vn = sprintf('%s_%d', vp{1}, kv);
        if ~any(strcmp(en_cell, vn))
            fprintf('  [MISSING VAR] %s\n', vn);
            pf_missing_vars{end+1} = vn; %#ok<AGROW>
        end
    end
end
for vn = {'GDP','pi','Q','TB'}
    if ~any(strcmp(en_cell, vn{1}))
        fprintf('  [MISSING VAR] %s\n', vn{1});
        pf_missing_vars{end+1} = vn{1}; %#ok<AGROW>
    end
end
if isempty(pf_missing_vars)
    fprintf('  All required endo vars found.\n');
end

% ---- (C) resol at original M_ (as compiled by dynare) --------------------
pf_exo_ss = zeros(M_.exo_nbr, 1);
pf_ned = 0; if isfield(M_,'exo_det_nbr'), pf_ned = M_.exo_det_nbr; end
pf_exo_det_ss = zeros(pf_ned, 1);
fprintf('  which(resol): %s\n', which('resol'));
fprintf('  resol() at original M_ ... ');
try
    [~, info_orig, ~] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, pf_exo_ss, pf_exo_det_ss);
    if info_orig(1) == 0
        fprintf('OK (info=0)\n');
    else
        fprintf('FAILED (info=%d)\n', info_orig(1));
    end
catch ME_orig
    fprintf('THREW: %s\n', ME_orig.message);
    info_orig = [97];
end

% ---- (D) resol at theta0-updated M_ (manual, no smm_model_moments) -------
fprintf('  resol() at theta0-updated M_ ... ');
M_pf = M_;
pf_scalar_params = {'ilabcosts','kappaV','rho_om1','sigma_om','rho_tfp1'};
pf_scalar_vals   = [theta0(1), exp(theta0(4)), theta0(5), theta0(6), theta0(7)];
for ks = 1:numel(pf_scalar_params)
    idx_s = find(strcmp(pn_cell, pf_scalar_params{ks}), 1);
    if ~isempty(idx_s), M_pf.params(idx_s) = pf_scalar_vals(ks); end
end
for ki = 1:nsec_pf
    idx_s = find(strcmp(pn_cell, sprintf('isigma_tfp_%d',ki)), 1);
    if ~isempty(idx_s), M_pf.params(idx_s) = theta0(7+ki); end
    idx_s = find(strcmp(pn_cell, sprintf('epsY_%d',ki)), 1);
    if ~isempty(idx_s), M_pf.params(idx_s) = theta0(2); end
    idx_s = find(strcmp(pn_cell, sprintf('epsM_%d',ki)), 1);
    if ~isempty(idx_s), M_pf.params(idx_s) = theta0(3); end
end
try
    [~, info_t0, ~] = resol(0, M_pf, options_, oo_.dr, oo_.dr.ys, pf_exo_ss, pf_exo_det_ss);
    if info_t0(1) == 0
        fprintf('OK (info=0)\n');
    else
        fprintf('FAILED (info=%d) — BK/numerical failure at theta0\n', info_t0(1));
    end
catch ME_t0
    fprintf('THREW: %s\n', ME_t0.message);
    info_t0 = [97];
end

% ---- (E) Full smm_model_moments test at theta0 ---------------------------
fprintf('  smm_model_moments(theta0) ... ');
[m_test, info_smm, ~] = smm_model_moments(theta0, M_, options_, oo_, baseline);
pf_ok = (info_smm(1) == 0) && ~any(isnan(m_test));

if pf_ok
    obj_test = (data_moments - m_test)' * W * (data_moments - m_test);
    fprintf('OK  obj(theta0) = %.6f\n', obj_test);
    mnames_pf = smm_moment_names();
    fprintf('\n  %-30s  %8s  %8s\n','Moment','Data','Model');
    for mm = 1:46
        fprintf('  %-30s  %8.5f  %8.5f\n', mnames_pf{mm}, data_moments(mm), m_test(mm));
    end
    fprintf('\n=== PRE-FLIGHT PASSED \u2014 launching optimizer ===\n\n');
else
    fprintf('FAILED  info=%d  NaN_idx=%s\n', info_smm(1), mat2str(find(isnan(m_test))'));
    fprintf('\n  --- DIAGNOSIS ---\n');
    if ~isempty(pf_missing_params)
        fprintf('  Missing params: %s\n', strjoin(pf_missing_params,', '));
        fprintf('  --> Recompile NK_SOE_lev_gap2.mod (run main_SOE_gap.m manually with EXERCISE=0).\n');
    end
    if ~isempty(pf_missing_vars)
        fprintf('  Missing endo vars: %s\n', strjoin(pf_missing_vars,', '));
        fprintf('  --> Variable name mismatch. Check Y_@{i}, PH_@{i}, L_@{i} in .mod file.\n');
    end
    if info_orig(1) ~= 0
        fprintf('  resol FAILS even at original M_ (info=%d).\n', info_orig(1));
        fprintf('  --> The compiled +NK_SOE_lev_gap2/ files are inconsistent with M_.\n');
        fprintf('  --> Delete +NK_SOE_lev_gap2/ folder, then run main_SOE_gap.m manually.\n');
    elseif info_t0(1) ~= 0
        fprintf('  resol FAILS specifically at theta0 (info=%d).\n', info_t0(1));
        fprintf('  --> These theta0 parameter values violate B-K conditions.\n');
        fprintf('  --> Reset theta0 to compiled calibration values:\n');
        fprintf('      kappaV=M_.params(find(strcmp(pn_cell,''kappaV''),1))\n');
    end
    if info_smm(1) == 97
        fprintf('  resol() threw inside smm_model_moments (info=97).\n');
        fprintf('  --> Most likely model not recompiled after structural change.\n');
        fprintf('  --> Close MATLAB, run main_SOE_gap.m freshly (not via SMM auto-run).\n');
    elseif info_smm(1) == 0 && any(isnan(m_test))
        nan_idx_diag = find(isnan(m_test));
        if all(ismember(nan_idx_diag, [44 45 46]))
            fprintf('  Rank-corr moments (44-46) are NaN.\n');
            fprintf('  Likely cause: params_val_ul.mat was saved from a single-shock exercise\n');
            fprintf('  (e.g. EXERCISE=2), so most isigma_tfp_i = 0. With only one active TFP\n');
            fprintf('  shock, I-O propagation makes all model std devs nearly identical and\n');
            fprintf('  Spearman rank-corr is undefined.\n');
            fprintf('  FIX: run main_SOE_gap.m with EXERCISE=0 (Baseline) to regenerate\n');
            fprintf('       params_val_ul.mat with all TFP shocks active, then re-run\n');
            fprintf('       smm_estimation.m. smm_model_moments has now been patched to\n');
            fprintf('       return 0 (not NaN) for undefined rank-corrs — if the issue\n');
            fprintf('       persists, clear the workspace and restart MATLAB.\n');
        else
            fprintf('  resol succeeded but moments are NaN.\n');
            fprintf('  NaN moments: %s\n', mat2str(nan_idx_diag'));
            fprintf('  --> Check endo_names above; missing variables return 0 std dev.\n');
        end
    end
    error('SMM:preflight', ...
          'Pre-flight FAILED (info=%d). Fix issues above, then re-run.', info_smm(1));
end

% --- Warm start: prefer best-so-far checkpoint, then final results ---
best_ckpt_file = fullfile(fileparts(mfilename('fullpath')), 'smm_best_so_far.mat');
prev_results_file = fullfile(fileparts(mfilename('fullpath')), 'smm_results.mat');
warm_started = false;

if exist(best_ckpt_file, 'file')
    try
        tmp_ckpt = load(best_ckpt_file, 'smm_best_so_far');
        theta_prev = tmp_ckpt.smm_best_so_far.theta_best;
        if numel(theta_prev) == n_theta && all(theta_prev >= lb) && all(theta_prev <= ub)
            theta0 = theta_prev;
            warm_started = true;
            fprintf('Warm start: loaded theta_best from smm_best_so_far.mat (obj=%.6f, saved=%s).\n\n', ...
                tmp_ckpt.smm_best_so_far.obj_best, tmp_ckpt.smm_best_so_far.saved_at);
        end
    catch
        fprintf('Could not load smm_best_so_far.mat checkpoint.\n');
    end
end

if ~warm_started && exist(prev_results_file, 'file')
    try
        tmp_prev = load(prev_results_file, 'smm_results');
        theta_prev = tmp_prev.smm_results.theta_hat;
        if numel(theta_prev) == n_theta && all(theta_prev >= lb) && all(theta_prev <= ub)
            theta0 = theta_prev;
            warm_started = true;
            fprintf('Warm start: loaded theta_hat from smm_results.mat (obj=%.6f).\n\n', ...
                tmp_prev.smm_results.obj_hat);
        end
    catch
        fprintf('Could not load smm_results.mat for warm start.\n');
    end
end

if ~warm_started
    fprintf('No valid warm-start file found - using default theta0.\n\n');
end

fprintf('Starting optimisation (%d parameters, %d moments)...\n\n', n_theta, 46);

% CMA-ES: Covariance Matrix Adaptation Evolution Strategy (Hansen, 2012).
% Derivative-free, robust to the ~1e-6 Dynare noise, adapts search
% covariance to the problem geometry — well suited for 21 heterogeneous
% parameters on different scales.
% insigma = (ub-lb)/6 so that the initial ±3-sigma interval covers the
% full feasible range; bounds are enforced natively by cmaes.m.
smm_fail_count = 0;
x0_col   = max(lb, min(ub, theta0(:)));   % column vector required by cmaes
insigma  = (ub(:) - lb(:)) / 6;

cma_opts = struct();
cma_opts.LBounds     = lb(:);
cma_opts.UBounds     = ub(:);
cma_opts.MaxFunEvals = 5000 * n_theta;
cma_opts.TolFun      = 1e-6;
cma_opts.TolX        = 1e-6;
cma_opts.Seed        = 42;
cma_opts.DispFinals  = 'on';
cma_opts.DispModulo  = 50;   % print every 50 iterations

fprintf('--- CMA-ES (blackbox, derivative-free) ---\n');
try
    [~, ~, ~, stopflag, ~, bestever] = cmaes('smm_obj_for_cmaes', x0_col, insigma, ...
        cma_opts, data_moments, W, baseline, M_, options_, oo_);
    theta_hat = max(lb(:), min(ub(:), bestever.x));
    fprintf('CMA-ES done.  obj = %.6f  stopflag: %s\n', bestever.f, strjoin(stopflag, ', '));
    fprintf('  (resol/NaN failures: %d)\n\n', smm_fail_count);
catch ME
    warning('SMM:cmaes', 'CMA-ES failed: %s. Falling back to theta0.', ME.message);
    theta_hat = max(lb, min(ub, theta0));
end


[obj_hat, psi_hat, moments_hat] = smm_obj_full(theta_hat, data_moments, W, baseline, M_, options_, oo_);

%% ================================================================== %%
%%  7. RESULTS TABLE                                                    %%
%% ================================================================== %%
fprintf('\n============================================================\n');
fprintf('  SMM RESULTS\n');
fprintf('============================================================\n\n');
fprintf('Objective value at theta_hat: %.6f\n\n', obj_hat);

moment_names = smm_moment_names();

fprintf('%-6s  %-15s  %10s  %10s\n','Idx','Parameter','Initial','Estimate');
fprintf('%s\n', repmat('-',48,1));
for k = 1:n_theta
    if k == 4
        fprintf('  %2d  %-15s  %10.4f  %10.4f  [log]\n', ...
            k, param_labels{k}, theta0(k), theta_hat(k));
        fprintf('  --  %-15s  %10.2e  %10.2e  [level]\n', ...
            'kappaV', exp(theta0(k)), exp(theta_hat(k)));
    else
        fprintf('  %2d  %-15s  %10.4f  %10.4f\n', k, param_labels{k}, theta0(k), theta_hat(k));
    end
end

fprintf('\nMoment fit:\n');
fprintf('%-34s  %9s  %9s  %9s\n','Moment','Data','Model','Diff');
fprintf('%s\n', repmat('-',66,1));
for m = 1:numel(moment_names)
    fprintf('%-34s  %9.5f  %9.5f  %+9.5f\n', ...
        moment_names{m}, data_moments(m), moments_hat(m), psi_hat(m));
end

%% ================================================================== %%
%%  8. FINAL SUMMARY                                                    %%
%% ================================================================== %%

fprintf('\n============================================================\n');
fprintf('  ESTIMATES\n');
fprintf('============================================================\n\n');
fprintf('%-6s  %-15s  %12s\n','Idx','Parameter','Estimate');
fprintf('%s\n', repmat('-',38,1));
for k = 1:n_theta
    if k == 4
        fprintf('  %2d  %-15s  %12.4f  [log]\n', ...
            k, param_labels{k}, theta_hat(k));
    else
        fprintf('  %2d  %-15s  %12.4f\n', ...
            k, param_labels{k}, theta_hat(k));
    end
end

%% ================================================================== %%
%%  9. SAVE AND APPLY ESTIMATES                                         %%
%% ================================================================== %%
% Save M_ so main_SOE_gap.m can use the SAME linearization for its moment
% table.  The .mod file embeds steady-state parameters (Ctot_ss, w_ss,
% PH_ss_i, ...) via load params_val_ul.mat at compile time.  When
% main_SOE_gap re-runs Dynare with estimated epsY/epsM, those SS parameters
% change → different linearization → different moments even at the same
% theta_hat.  Using M_smm ensures both scripts evaluate the identical model.
M_smm       = M_;
options_smm = options_;
oo_smm      = oo_;

smm_results = struct( ...
    'theta_hat',    theta_hat, ...
    'obj_hat',      obj_hat, ...
    'psi_hat',      psi_hat, ...
    'moments_hat',  moments_hat, ...
    'data_moments', data_moments, ...
    'param_labels', {param_labels}, ...
    'moment_names', {moment_names}, ...
    'W',            W);
save('smm_results.mat', 'smm_results', 'M_smm', 'options_smm', 'oo_smm');
fprintf('\nResults saved to smm_results.mat\n');

% Write estimates to smm_estimates.mat so main_SOE.m picks them up automatically.
% Variable names match those used in main_SOE.m; loading the file there
% (EXERCISE == 0) directly overrides the hard-coded defaults.
ilabcosts_val  = theta_hat(1);
modepsY        = ones(nsec,1) * theta_hat(2);
modepsM        = ones(nsec,1) * theta_hat(3);
kappaV_val     = exp(theta_hat(4));
rho_om1_val    = theta_hat(5);
sigma_om_val   = theta_hat(6);
rho_tfp1_val   = theta_hat(7);
isigma_tfp_val = theta_hat(8:19);
rho_pvstar_val   = theta_hat(20);
sigma_pvstar_val = theta_hat(21);
rho_xi_val       = theta_hat(22);
sigma_xi_val     = theta_hat(23);
save('smm_estimates.mat', 'ilabcosts_val', 'modepsY', 'modepsM', ...
    'kappaV_val', 'rho_om1_val', 'sigma_om_val', 'rho_tfp1_val', ...
    'isigma_tfp_val', 'rho_pvstar_val', 'sigma_pvstar_val', ...
    'rho_xi_val', 'sigma_xi_val');
fprintf('smm_estimates.mat saved.  Re-run main_SOE_gap.m (EXERCISE=0) to apply estimates.\n\n');

% ==================================================================== %
%  LOCAL FUNCTIONS  (MATLAB R2016b+: allowed at end of script files)   %
% ==================================================================== %

function obj = smm_obj(theta, d_mom, W_mat, bas, M_in, opt_in, oo_in)
% Scalar SMM objective for use with fminsearch / fmincon.
% smm_fail_count in the caller's workspace is incremented on solver failure
% so the user can see how often the model rejects candidate parameters.
    [m_model, info] = smm_model_moments(theta, M_in, opt_in, oo_in, bas);
    if info(1) ~= 0 || any(isnan(m_model))
        % Increment diagnostic counter (lives in smm_estimation.m's workspace)
        try
            assignin('caller', 'smm_fail_count', ...
                evalin('caller','smm_fail_count') + 1);
        catch
        end
        obj = 1e8;
        return
    end
    psi = d_mom - m_model;
    obj = psi' * W_mat * psi;
end

function [obj, psi, m_model] = smm_obj_full(theta, d_mom, W_mat, bas, M_in, opt_in, oo_in)
% Full-output version for reporting.
    [m_model, info] = smm_model_moments(theta, M_in, opt_in, oo_in, bas);
    if info(1) ~= 0 || any(isnan(m_model))
        obj = 1e8;  psi = NaN(numel(d_mom),1);  m_model = NaN(numel(d_mom),1);
        return
    end
    psi = d_mom - m_model;
    obj = psi' * W_mat * psi;
end

function names = smm_moment_names()
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
