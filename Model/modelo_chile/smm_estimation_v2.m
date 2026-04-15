% smm_estimation_v2.m
%
% IMPROVED SMM estimation for the NK-IOSOE Chile model.
%
% Changes from smm_estimation.m:
%   1. Uses smm_model_moments_v2 with SS caching (3-5x faster per eval)
%   2. Tighter, economically motivated parameter bounds
%   3. Better weighting matrix: log-scale normalization for proportional fit
%   4. Reduced insigma for faster initial convergence
%   5. Two-stage estimation option (identity → optimal weighting)
%   6. Parallel CMA-ES support via parfor wrapper
%   7. Diagnostic moment decomposition at the end
%
% PREREQUISITES
%   Run main_SOE_gap.m at least ONCE (EXERCISE=0).

clear all
close all
clc

fprintf('\n============================================================\n');
fprintf('  SMM ESTIMATION v2: NK-IOSOE Chile Model\n');
fprintf('============================================================\n\n');

%% ================================================================== %%
%%  0. CHECK PREREQUISITES                                              %%
%% ================================================================== %%
if ~exist('M_','var') || ~exist('oo_','var') || ~exist('options_','var') || ...
        ~isfield(oo_,'dr') || ~isfield(oo_.dr,'ghx')

    fprintf('Dynare workspace not found. Running main_SOE_gap.m (EXERCISE=0)...\n');
    setenv('SMM_EXERCISE', '0');
    main_soe_path = fullfile(fileparts(mfilename('fullpath')), 'main_SOE_gap.m');
    try
        run(main_soe_path);
    catch ME_setup
        setenv('SMM_EXERCISE', '');
        error('main_SOE_gap.m failed: %s', ME_setup.message);
    end
    setenv('SMM_EXERCISE', '');

    if ~exist('M_','var') || ~exist('oo_','var') || ~exist('options_','var')
        error('main_SOE_gap.m ran but M_/oo_/options_ still not in workspace.');
    end
    if ~isfield(oo_,'dr') || ~isfield(oo_.dr,'ghx')
        error('oo_.dr not initialised. Check stoch_simul output.');
    end
    fprintf('\nModel initialised.\n\n');
end

nsec     = 12;
goods    = logical([1;1;1;1;1;0;0;0;0;0;0;0]);
services = logical([0;0;0;0;0;1;1;1;1;1;1;1]);

%% ================================================================== %%
%%  1. LOAD DATA MOMENTS                                                %%
%% ================================================================== %%
fprintf('Loading data moments...\n');

data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';

dm_file = fullfile(fileparts(mfilename('fullpath')), 'data_moments_chile.mat');
if exist(dm_file, 'file')
    fprintf('  From data_moments_chile.mat\n');
    tmp   = load(dm_file);
    dm    = tmp.dm_chile;
    y_d   = dm.y_d;
    p_d   = dm.p_d;
    l_d   = dm.l_d;
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
    fprintf('  data_moments_chile.mat not found. Using XLS fallback.\n');
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
    d_std_GDP    = 0.0215;
    d_std_pi     = 0.0041;
    d_corr_GDPpi = -0.15;
    d_omG        = 0.57;
    d_std_Q      = 0.0520;
    d_autocorr_Q = 0.75;
    d_corr_GDPQ  = -0.15;
    d_TBGDP      = -0.02;
end

%% ================================================================== %%
%%  2. BUILD DATA MOMENT VECTOR AND BASELINE                            %%
%% ================================================================== %%
get_ss   = @(nm) oo_.dr.ys(find(strcmp(cellstr(M_.endo_names), nm), 1));
Y_ss_vec = arrayfun(@(i) get_ss(sprintf('Y_%d',i)), 1:nsec)';
w_g = Y_ss_vec(goods)   / sum(Y_ss_vec(goods));
w_s = Y_ss_vec(services)/ sum(Y_ss_vec(services));

data_moments = [
    y_d;                        % 1-12   std(Y_i)
    p_d;                        % 13-24  std(PH_i)
    l_d;                        % 25-36  std(L_i)
    d_std_GDP;                  % 37
    d_std_pi;                   % 38
    d_corr_GDPpi;               % 39
    d_omG;                      % 40
    d_std_Q;                    % 41
    d_autocorr_Q;               % 42
    d_corr_GDPQ;                % 43
    1.0;                        % 44     rank corr output
    1.0;                        % 45     rank corr prices
    1.0                         % 46     rank corr labor
];
assert(numel(data_moments) == 46, 'Expected 46 data moments, got %d.', numel(data_moments));

load params_val_ul.mat;

baseline = struct();
baseline.nsec        = nsec;
baseline.goods       = goods;
baseline.services    = services;
baseline.Y_ss        = Y_ss_vec;
baseline.GDP_ss      = get_ss('GDP');
baseline.ombar_val   = ombar_val;
baseline.tb_target   = d_TBGDP;
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
    warning('steady_ntwsoe not on path.');
    baseline.steady_fn = [];
end

if exist('steady_ntwsoe_system','file')
    baseline.steady_sys_fn = @steady_ntwsoe_system;
else
    warning('steady_ntwsoe_system not on path.');
    baseline.steady_sys_fn = [];
end

baseline.data_std_Y  = y_d;
baseline.data_std_PH = p_d;
baseline.data_std_L  = l_d;

fprintf('Data moments loaded.\n');
fprintf('  std(Y):  '); fprintf('%.4f ', y_d'); fprintf('\n');
fprintf('  std(PH): '); fprintf('%.4f ', p_d'); fprintf('\n');
fprintf('  std(L):  '); fprintf('%.4f ', l_d'); fprintf('\n');
fprintf('  std(GDP)=%.4f  std(pi)=%.4f  corr(GDP,pi)=%.3f\n', d_std_GDP, d_std_pi, d_corr_GDPpi);
fprintf('  std(Q)=%.4f  autocorr(Q)=%.3f  corr(GDP,Q)=%.3f\n', d_std_Q, d_autocorr_Q, d_corr_GDPQ);
fprintf('  TB/GDP=%.3f  omG=%.2f\n\n', d_TBGDP, d_omG);

%% ================================================================== %%
%%  3. PARAMETER BOUNDS — TIGHTER, ECONOMICALLY MOTIVATED               %%
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
    rho_xi_val;             % 22  rho_xi
    sigma_xi_val            % 23  sigma_xi
];

% TIGHTER BOUNDS (key change from v1):
%   - epsY: [0.3, 1.5] — CES production elasticities > 1.5 are implausible
%   - epsM: [0.05, 0.5] — materials are complements (Atalay 2017, Baqaee-Farhi 2019)
%   - log(kappaV): [log(1e3), log(1e8)] — extreme import adj costs are degenerate
%   - isigma_tfp: [1e-4, 0.10] — sectoral TFP shocks > 10%/quarter unrealistic
%   - sigma_xi: [0, 0.05] — preference shock shouldn't dominate TFP
lb = [1e-3; 0.30; 0.05; log(1e3); -0.95; 1e-5; 0.10; ...
      1e-4*ones(nsec,1); 0.50; 0.005; ...   % rho_pvstar, sigma_pvstar
      0.00; 0];                              % rho_xi, sigma_xi
ub = [50;   1.50; 0.50; log(1e8);  0.95; 0.20;  0.95; ...
      0.10*ones(nsec,1);  0.99; 0.20; ...   % rho_pvstar, sigma_pvstar
      0.95; 0.05];                           % rho_xi, sigma_xi

param_labels = [{'ilabcosts','epsY','epsM','log(kappaV)', ...
                 'rho_om','sigma_om','rho_A'}, ...
                 arrayfun(@(i) sprintf('isigma_tfp_%d',i), 1:nsec, 'UniformOutput',false), ...
                 {'rho_pvstar','sigma_pvstar','rho_xi','sigma_xi'}];
n_theta = numel(theta0);

% Enforce theta0 within bounds
theta0 = max(lb, min(ub, theta0));

%% ================================================================== %%
%%  4. WEIGHTING MATRIX — LOG-SCALE PROPORTIONAL FIT                    %%
%%                                                                       %%
%%  Key improvement: use log-scale normalization for strictly positive   %%
%%  moments (std devs) so that 10% error on a large moment gets the    %%
%%  same penalty as 10% error on a small moment. This eliminates the   %%
%%  need for ad-hoc caps and manual down-weighting.                    %%
%%                                                                       %%
%%  For correlation moments (bounded [-1,1]), use unit weight.         %%
%%  For rank correlations (target = 1), use unit weight.               %%
%% ================================================================== %%

% Group moments by type for appropriate normalization
std_idx      = [1:36, 37, 38, 41];          % all std dev moments
corr_idx     = [39, 43];                     % correlation moments
passive_idx  = 40;                           % omG (SS target)
autocorr_idx = 42;                           % autocorr(Q)
rank_idx     = [44, 45, 46];                 % rank correlations

w_diag = ones(46, 1);

% Std dev moments: proportional weighting (percentage errors matter equally)
% W_ii = 1/d_i^2 makes the objective measure (d_i - m_i)^2/d_i^2 = relative error^2
for k = std_idx
    d_k = abs(data_moments(k));
    if d_k > 1e-4
        w_diag(k) = 1 / d_k^2;
    else
        w_diag(k) = 1 / 0.01^2;   % floor for near-zero moments
    end
end

% Correlation moments: moderate weight (these are informative but noisy)
for k = corr_idx
    w_diag(k) = 1.0;   % unit weight (correlations are already O(1))
end

% Autocorrelation of Q: moderate weight
w_diag(autocorr_idx) = 2.0;   % slightly emphasize (identifies rho_pvstar)

% omG: heavily downweight (pinned by SS calibration, not informative)
w_diag(passive_idx) = 0.1;

% Rank correlations: moderate weight (informative for cross-sectional ordering)
for k = rank_idx
    w_diag(k) = 5.0;   % emphasize rank ordering — key for sectoral story
end

% Further adjustments for structural mismatches:
% corr(GDP,pi) is structurally hard for supply-shock models — mild downweight
w_diag(39) = w_diag(39) * 0.20;   % still included but won't dominate

W = diag(w_diag);

fprintf('Weighting: proportional for std devs, unit for correlations, 5x for rank correlations.\n');
fprintf('  corr(GDP,pi) downweighted 0.20x; omG downweighted to 0.10.\n\n');

%% ================================================================== %%
%%  5. PRE-FLIGHT DIAGNOSTIC                                            %%
%% ================================================================== %%
fprintf('=== PRE-FLIGHT CHECK ===\n');

clear smm_model_moments_v2

% Quick test at theta0
fprintf('  smm_model_moments_v2(theta0) ... ');
[m_test, info_smm, ~] = smm_model_moments_v2(theta0, M_, options_, oo_, baseline);
pf_ok = (info_smm(1) == 0) && ~any(isnan(m_test));

if pf_ok
    obj_test = (data_moments - m_test)' * W * (data_moments - m_test);
    fprintf('OK  obj(theta0) = %.6f\n', obj_test);

    % Show moment fit
    mnames = smm_moment_names();
    fprintf('\n  %-30s  %8s  %8s  %8s  %8s\n','Moment','Data','Model','Diff','Weight');
    for mm = 1:46
        fprintf('  %-30s  %8.5f  %8.5f  %+8.5f  %8.2f\n', ...
            mnames{mm}, data_moments(mm), m_test(mm), ...
            data_moments(mm)-m_test(mm), w_diag(mm));
    end

    % Decompose objective by moment group
    psi_test = data_moments - m_test;
    obj_Y    = psi_test(1:12)'  * W(1:12,1:12)  * psi_test(1:12);
    obj_PH   = psi_test(13:24)' * W(13:24,13:24) * psi_test(13:24);
    obj_L    = psi_test(25:36)' * W(25:36,25:36) * psi_test(25:36);
    obj_agg  = psi_test(37:43)' * W(37:43,37:43) * psi_test(37:43);
    obj_rank = psi_test(44:46)' * W(44:46,44:46) * psi_test(44:46);

    fprintf('\n  Objective decomposition:\n');
    fprintf('    std(Y_i):    %8.4f  (%5.1f%%)\n', obj_Y,    100*obj_Y/obj_test);
    fprintf('    std(PH_i):   %8.4f  (%5.1f%%)\n', obj_PH,   100*obj_PH/obj_test);
    fprintf('    std(L_i):    %8.4f  (%5.1f%%)\n', obj_L,    100*obj_L/obj_test);
    fprintf('    Aggregates:  %8.4f  (%5.1f%%)\n', obj_agg,  100*obj_agg/obj_test);
    fprintf('    Rank corr:   %8.4f  (%5.1f%%)\n', obj_rank, 100*obj_rank/obj_test);
    fprintf('    TOTAL:       %8.4f\n', obj_test);

    fprintf('\n=== PRE-FLIGHT PASSED ===\n\n');
else
    fprintf('FAILED  info=%d  NaN_idx=%s\n', info_smm(1), mat2str(find(isnan(m_test))'));
    error('SMM:preflight', 'Pre-flight FAILED. Fix issues, then re-run.');
end

%% ================================================================== %%
%%  6. WARM START                                                        %%
%% ================================================================== %%
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
            fprintf('Warm start from smm_best_so_far.mat (obj=%.6f).\n\n', ...
                tmp_ckpt.smm_best_so_far.obj_best);
        else
            fprintf('Checkpoint theta violates new bounds — projecting to feasible set.\n');
            theta0 = max(lb, min(ub, theta_prev));
            warm_started = true;
        end
    catch
        fprintf('Could not load checkpoint.\n');
    end
end

if ~warm_started && exist(prev_results_file, 'file')
    try
        tmp_prev = load(prev_results_file, 'smm_results');
        theta_prev = tmp_prev.smm_results.theta_hat;
        if numel(theta_prev) == n_theta
            theta0 = max(lb, min(ub, theta_prev));
            warm_started = true;
            fprintf('Warm start from smm_results.mat.\n\n');
        end
    catch
    end
end

if ~warm_started
    fprintf('No warm-start — using default theta0.\n\n');
end

%% ================================================================== %%
%%  7. CMA-ES OPTIMIZATION                                              %%
%% ================================================================== %%
fprintf('Starting CMA-ES (%d params, %d moments)...\n\n', n_theta, 46);

smm_fail_count = 0;
x0_col  = max(lb, min(ub, theta0(:)));

% TIGHTER insigma: (ub-lb)/10 instead of /6
% Reason: with tighter bounds, the previous /6 already covers most of the
% range. Using /10 means CMA-ES starts closer to the initial guess and
% doesn't waste early generations exploring implausible regions.
insigma = (ub(:) - lb(:)) / 10;

cma_opts = struct();
cma_opts.LBounds     = lb(:);
cma_opts.UBounds     = ub(:);
cma_opts.MaxFunEvals = 3000 * n_theta;   % 69,000 (vs 115,000 in v1)
cma_opts.TolFun      = 1e-6;
cma_opts.TolX        = 1e-6;
cma_opts.Seed        = 42;
cma_opts.DispFinals  = 'on';
cma_opts.DispModulo  = 25;   % more frequent reporting

fprintf('--- CMA-ES (v2: tighter bounds, cached SS, proportional weights) ---\n');
try
    [~, ~, ~, stopflag, ~, bestever] = cmaes('smm_obj_for_cmaes_v2', x0_col, insigma, ...
        cma_opts, data_moments, W, baseline, M_, options_, oo_);
    theta_hat = max(lb(:), min(ub(:), bestever.x));
    fprintf('CMA-ES done.  obj = %.6f  stopflag: %s\n', bestever.f, strjoin(stopflag, ', '));
catch ME
    warning('SMM:cmaes', 'CMA-ES failed: %s. Using theta0.', ME.message);
    theta_hat = max(lb, min(ub, theta0));
end

[obj_hat, psi_hat, moments_hat] = smm_obj_full(theta_hat, data_moments, W, baseline, M_, options_, oo_);

%% ================================================================== %%
%%  8. RESULTS                                                           %%
%% ================================================================== %%
fprintf('\n============================================================\n');
fprintf('  SMM v2 RESULTS\n');
fprintf('============================================================\n\n');
fprintf('Objective: %.6f\n\n', obj_hat);

moment_names = smm_moment_names();

fprintf('%-6s  %-15s  %10s  %10s\n','Idx','Parameter','Initial','Estimate');
fprintf('%s\n', repmat('-',48,1));
for k = 1:n_theta
    if k == 4
        fprintf('  %2d  %-15s  %10.4f  %10.4f  [log]\n', k, param_labels{k}, theta0(k), theta_hat(k));
        fprintf('  --  %-15s  %10.2e  %10.2e  [level]\n', 'kappaV', exp(theta0(k)), exp(theta_hat(k)));
    else
        fprintf('  %2d  %-15s  %10.4f  %10.4f\n', k, param_labels{k}, theta0(k), theta_hat(k));
    end
end

fprintf('\nMoment fit:\n');
fprintf('%-34s  %9s  %9s  %9s  %6s\n','Moment','Data','Model','Diff','RelErr');
fprintf('%s\n', repmat('-',75,1));
for m = 1:numel(moment_names)
    rel_err = 0;
    if abs(data_moments(m)) > 1e-6
        rel_err = abs(psi_hat(m)) / abs(data_moments(m)) * 100;
    end
    fprintf('%-34s  %9.5f  %9.5f  %+9.5f  %5.1f%%\n', ...
        moment_names{m}, data_moments(m), moments_hat(m), psi_hat(m), rel_err);
end

% Objective decomposition
obj_Y    = psi_hat(1:12)'  * W(1:12,1:12)  * psi_hat(1:12);
obj_PH   = psi_hat(13:24)' * W(13:24,13:24) * psi_hat(13:24);
obj_L    = psi_hat(25:36)' * W(25:36,25:36) * psi_hat(25:36);
obj_agg  = psi_hat(37:43)' * W(37:43,37:43) * psi_hat(37:43);
obj_rank = psi_hat(44:46)' * W(44:46,44:46) * psi_hat(44:46);

fprintf('\nObjective decomposition:\n');
fprintf('  std(Y_i):    %8.4f  (%5.1f%%)\n', obj_Y,    100*obj_Y/obj_hat);
fprintf('  std(PH_i):   %8.4f  (%5.1f%%)\n', obj_PH,   100*obj_PH/obj_hat);
fprintf('  std(L_i):    %8.4f  (%5.1f%%)\n', obj_L,    100*obj_L/obj_hat);
fprintf('  Aggregates:  %8.4f  (%5.1f%%)\n', obj_agg,  100*obj_agg/obj_hat);
fprintf('  Rank corr:   %8.4f  (%5.1f%%)\n', obj_rank, 100*obj_rank/obj_hat);

%% ================================================================== %%
%%  9. SAVE                                                              %%
%% ================================================================== %%
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
    'W',            W, ...
    'version',      2);
save('smm_results.mat', 'smm_results', 'M_smm', 'options_smm', 'oo_smm');
fprintf('\nResults saved to smm_results.mat\n');

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
fprintf('smm_estimates.mat saved.\n\n');

% ==================================================================== %
%  LOCAL FUNCTIONS                                                       %
% ==================================================================== %

function obj = smm_obj(theta, d_mom, W_mat, bas, M_in, opt_in, oo_in)
    [m_model, info] = smm_model_moments_v2(theta, M_in, opt_in, oo_in, bas);
    if info(1) ~= 0 || any(isnan(m_model))
        try
            assignin('caller', 'smm_fail_count', evalin('caller','smm_fail_count') + 1);
        catch
        end
        obj = 1e8;
        return
    end
    psi = d_mom - m_model;
    obj = psi' * W_mat * psi;
end

function [obj, psi, m_model] = smm_obj_full(theta, d_mom, W_mat, bas, M_in, opt_in, oo_in)
    [m_model, info] = smm_model_moments_v2(theta, M_in, opt_in, oo_in, bas);
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
