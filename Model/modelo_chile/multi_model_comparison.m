%% ===========================================================================
% TERMS OF TRADE SHOCK - COMPARE MULTIPLE MODEL VARIANTS
% ===========================================================================
% This script runs 4 model configurations and compares IRFs (output gaps):
% 1. Baseline:   Sticky price + Full I-O + Labor adjustment costs
% 2. Flexible:   Flexible price + Full I-O + No labor costs
% 3. NoLabor:    Sticky price + Full I-O + No labor costs
% 4. NearClosed: Sticky price + Full I-O + Minimal trade
% Uses NK_SOE_lev_gap2.mod: IRFs reported as output gaps (deviation from flex-price equilibrium).
% ===========================================================================

clear all;
close all;
clc;

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║  TERMS OF TRADE SHOCK: MULTI-MODEL COMPARISON                     ║\n');
fprintf('║  Comparing: Baseline vs Flexible Price vs No I-O                  ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

restoredefaultpath;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultFigureWindowStyle','normal');

%% Add paths
this_dir = fileparts(mfilename('fullpath'));
pths = { ...
    this_dir ...
    fullfile(this_dir, '..', 'utils') ...
    'C:\Program Files\Dynare\6.4\matlab' ...
    'C:\Program Files\Dynare\6.2\matlab' ...
    'C:\Users\mgiarda\OneDrive - Banco Central de Chile\Escritorio\Github_BC\DPMP-DME-Modelos\Model\utils' ...
    'C:\Users\mgiarda\OneDrive - Banco Central de Chile\Escritorio\Proyectos\ntw_data_github' ...
    'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github' ...
};
for ip = 1:numel(pths)
    if isfolder(pths{ip})
        addpath(pths{ip});
    end
end

% Ensure working directory is modelo_chile so save/load params_val.mat are consistent
cd(this_dir);
fprintf('Working directory: %s\n', this_dir);
fprintf('Run logs: dynare_comparison_*.txt | Figures: *_comparison_irf.png\n\n');

% Verify required helper functions are on the path
assert(~isempty(which('steady_ntwsoe')), ...
    'steady_ntwsoe not found on path. Check this_dir = %s', this_dir);
assert(~isempty(which('steady_ntwsoe_system')), ...
    'steady_ntwsoe_system not found on path.');

nsec = 12;

%% Load data (same as main_SOE.m)
data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';
filename_industries = fullfile(data_path, 'Stata_to_excel_few_industries_chile.xls');
filename_io         = fullfile(data_path, 'IO_2021_chile.csv');
filename_fpa        = fullfile(datade tsi _path, 'fpa_vector_few_industries_chile.csv');

fprintf('Reading data files...\n');
varNames = {'kk','BEAName','Nameshort','dp','dy','dl','dri','share','industrytype','spend_good','spend_serv','alpha','alpha_V', 'var_rho'};
varTypes = {'double','char','char','double','double','double','double','double','categorical','double','double','double','double','double'};
opts = spreadsheetImportOptions('NumVariables',14,'VariableNames',varNames,'VariableTypes',varTypes,'DataRange','A2');

alldata    = readtable(filename_industries, opts);
names      = table2array(alldata(:,3));
alpha      = table2array(alldata(:,12));
alpha_V    = table2array(alldata(:,13));
var_rho    = table2array(alldata(:,14));
spend_good = table2array(alldata(:,10));
spend_serv = table2array(alldata(:,11));

% Data guards
alpha(~isfinite(alpha))     = 0;
alpha_V(~isfinite(alpha_V)) = 0;
var_rho(~isfinite(var_rho)) = 0.5;
alpha   = max(alpha, 0);
alpha_V = max(alpha_V, 0);
for ii = 1:numel(alpha)
    if alpha(ii) + alpha_V(ii) >= 0.98
        scale = 0.98 / max(alpha(ii) + alpha_V(ii), 1e-12);
        alpha(ii)   = alpha(ii) * scale;
        alpha_V(ii) = alpha_V(ii) * scale;
    end
end
var_rho = min(max(var_rho, 1e-4), 1-1e-4);

% I-O matrix (full)
betaio = readmatrix(filename_io);
betaio = betaio(1:nsec,1:nsec);
betaio(~isfinite(betaio) | betaio < 0) = 0;
betax = betaio ./ sum(betaio,1);
modbeta_full = betax';

% Diagonal I-O (preserve diagonal share, renormalize)
modbeta_diag = diag(diag(modbeta_full));
for jj = 1:nsec
    cs = sum(modbeta_diag(:,jj));
    if cs > 0
        modbeta_diag(:,jj) = modbeta_diag(:,jj)/cs;
    else
        modbeta_diag(jj,jj) = 1.0;
    end
end

% Price flexibility
fpa = readmatrix(filename_fpa);
fpa = fpa(1:nsec);
fpa = ones(nsec,1) * mean(fpa);
theta = (1-fpa).^3;

beta_val     = 0.986;
epsilon_val  = 10;
kappa_calib  = theta*epsilon_val./((1-theta).*(1-theta*beta_val));
kappa_sticky = kappa_calib;
kappa_flex   = ones(nsec,1) * 0.001;

% Spending shares
gammag = spend_good / sum(spend_good);
gammas = spend_serv / sum(spend_serv);

%% Scalar parameters (mirroring main_tot_shock.m)
phi_val        = 1.2;
rhoi_val       = 0.6;
rhoirule_val   = 0.74;
ilabcosts_val  = 0.1;
ombar_val      = 0.57;

% Shock scenario to run: 'tot' or 'manufacturing'
shock_scenario = 'manufacturing';

% Baseline non-ToT shocks
sigma_i_val     = 0.0;
rho_om1_val     = 0.1;  rho_om2_val  = 0;    sigma_om_val  = 0.0;
rho_tfp1_val    = 0.95; rho_tfp2_val = 0.0;
rho_val         = 0.1;
sigma_L_agg_val = 0.0;
isigma_tfp_val  = zeros(nsec,1);

% Terms of Trade shock calibration (used when shock_scenario='tot')
rho_pvstar_val   = 0.9;
sigma_pvstar_val = -0.10;

% Exogenous shock activation controls used in NK_IOSOE_tot.mod shocks block
shock_eps_pvstar_val = 1.0;
shock_eps_om_val     = 0.0;
shock_eps_i_val      = 0.0;
shock_epschi_val     = 0.0;
shock_epsA_val       = zeros(nsec,1);

% Manufacturing TFP shock setup
manufacturing_keywords = {'manufact', 'manuf', 'manufactura', 'industria', 'industry'};
manufacturing_idx = [];
name_strings = strings(nsec,1);
for ii = 1:nsec
    name_strings(ii) = lower(string(names{ii}));
end
for kk = 1:numel(manufacturing_keywords)
    idx_try = find(contains(name_strings, manufacturing_keywords{kk}), 1, 'first');
    if ~isempty(idx_try)
        manufacturing_idx = idx_try;
        break;
    end
end
if isempty(manufacturing_idx)
    manufacturing_idx = 3;  % fallback if label matching fails
end

manufacturing_tfp_size = -0.10;   % 10% negative TFP innovation in manufacturing sector

switch lower(shock_scenario)
    case 'tot'
        % ToT shock active, sectoral TFP shocks inactive
        shock_eps_pvstar_val = 1.0;
        shock_epsA_val       = zeros(nsec,1);
        isigma_tfp_val       = zeros(nsec,1);
        fprintf('Shock scenario: Terms of Trade (eps_pvstar)\n');

    case 'manufacturing'
        % Manufacturing-sector TFP shock active, ToT shock inactive
        sigma_pvstar_val = 0.0;
        shock_eps_pvstar_val = 0.0;
        shock_epsA_val       = zeros(nsec,1);
        shock_epsA_val(manufacturing_idx) = 1.0;
        isigma_tfp_val       = zeros(nsec,1);
        isigma_tfp_val(manufacturing_idx) = manufacturing_tfp_size;
        fprintf('Shock scenario: Manufacturing TFP (epsA_%d), sector="%s"\n', ...
            manufacturing_idx, string(names{manufacturing_idx}));

    otherwise
        error('Unknown shock_scenario="%s". Use ''tot'' or ''manufacturing''.', shock_scenario);
end

% Production
modepsM    = ones(nsec,1) * 0.1;
modepsY    = ones(nsec,1) * 0.8;
modchiX    = ones(nsec,1) * 1/nsec;
modvarrho  = var_rho;
modA       = ones(nsec,1);
modalphaV  = alpha_V;

% SOE parameters
Pistar_ss_val  = 1.00;
Rworld_ss_val  = Pistar_ss_val / beta_val;
PVstar_ss      = 1;
kappaV_val     = 1e13;
epsilonV_val   = 1e13;
epsilonX_val   = 1;
omegaX_val     = 1;
chii_b_val     = 0.001;
etastar_val    = 1;
ystar_ss_val   = 1;
sigmaH_val     = 0.999;
gamma_val      = 2;
tb_target      = -0.02;

goods    = spend_good > spend_serv;
services = spend_serv > spend_good;
modgammag = gammag;
modgammas = gammas;

%% Define model configurations
fprintf('\nDefining 4 model configurations (stable models only)...\n');
fprintf('Note: Diagonal I-O configs dropped due to Dynare solver issues with flexible+diagonal combo\n\n');

configs = struct();

% Config 1: NEW BASELINE - STICKY + FULL I-O + LABOR ADJUSTMENT COSTS
configs(1).name = 'Baseline';
configs(1).description = 'Sticky Price + Full I-O + Labor Adjustment Costs';
configs(1).modbeta = modbeta_full;
configs(1).modkappa = kappa_sticky;
configs(1).modalpha = alpha;
configs(1).modcl = ones(nsec,1) * 10;  % Labor costs enabled
configs(1).modclneg = ones(nsec,1) * 10;
configs(1).modalphaV = modalphaV;
configs(1).modvarrho = modvarrho;
configs(1).omegaX = omegaX_val;

% Config 2: FLEXIBLE PRICE + FULL I-O + NO LABOR COSTS
configs(2).name = 'Flexible';
configs(2).description = 'Flexible Price + Full I-O + No Labor Costs';
configs(2).modbeta = modbeta_full;
configs(2).modkappa = kappa_flex;
configs(2).modalpha = alpha;
configs(2).modcl = zeros(nsec,1);
configs(2).modclneg = zeros(nsec,1);
configs(2).modalphaV = modalphaV;
configs(2).modvarrho = modvarrho;
configs(2).omegaX = omegaX_val;

% Config 3: LABOR COSTS REMOVED (Baseline without labor adjustment costs)
configs(3).name = 'NoLabor';
configs(3).description = 'Sticky Price + Full I-O + No Labor Costs';
configs(3).modbeta = modbeta_full;
configs(3).modkappa = kappa_sticky;
configs(3).modalpha = alpha;
configs(3).modcl = zeros(nsec,1);
configs(3).modclneg = zeros(nsec,1);
configs(3).modalphaV = modalphaV;
configs(3).modvarrho = modvarrho;
configs(3).omegaX = omegaX_val;

% Config 4: NEAR-CLOSED ECONOMY
configs(4).name = 'NearClosed';
configs(4).description = 'Near-Closed Economy (Sticky + Full I-O + Minimal Trade)';
configs(4).modbeta = modbeta_full;
configs(4).modkappa = kappa_sticky;
configs(4).modalpha = alpha;
configs(4).modcl = zeros(nsec,1);
configs(4).modclneg = zeros(nsec,1);
configs(4).modalphaV = modalphaV * 0.05;        % 5% of calibrated import share
configs(4).modvarrho = 1 - (1 - modvarrho)*0.05; % 5% of calibrated home bias gap
configs(4).omegaX = omegaX_val * 0.05;           % 5% of calibrated export scale
configs(4).tb_target_cfg = 0.0;                   % Balanced trade (near-zero trade → TB≈0)

%% Run each model configuration
results = struct('success', {});

for cfg = 1:length(configs)
    fprintf('════════════════════════════════════════════════\n');
    fprintf('Configuration %d: %s\n', cfg, configs(cfg).description);
    fprintf('════════════════════════════════════════════════\n\n');
    
    % Set per-config parameters in base workspace
    modbeta  = configs(cfg).modbeta;
    modkappa = configs(cfg).modkappa;
    modalpha = configs(cfg).modalpha;
    modcl    = configs(cfg).modcl;
    modclneg = configs(cfg).modclneg;
    modcm    = zeros(nsec,1);
    
    % Override SOE parameters if specified in config (for closed economy)
    if isfield(configs(cfg), 'modalphaV')
        modalphaV_cfg = configs(cfg).modalphaV;
    else
        modalphaV_cfg = modalphaV;
    end
    
    if isfield(configs(cfg), 'modvarrho')
        modvarrho_cfg = configs(cfg).modvarrho;
    else
        modvarrho_cfg = modvarrho;
    end
    
    if isfield(configs(cfg), 'omegaX')
        omegaX_cfg = configs(cfg).omegaX;
    else
        omegaX_cfg = omegaX_val;
    end

    assignin('base', 'nsec',      nsec);
    assignin('base', 'modbeta',   modbeta);
    assignin('base', 'modkappa',  modkappa);
    assignin('base', 'modalpha',  modalpha);
    assignin('base', 'modalphaV', modalphaV_cfg);
    assignin('base', 'modgammag', modgammag);
    assignin('base', 'modgammas', modgammas);
    assignin('base', 'names',     names);
    assignin('base', 'modcl',     modcl);
    assignin('base', 'modclneg',  modclneg);
    assignin('base', 'modcm',     modcm);
    assignin('base', 'modepsM',   modepsM);
    assignin('base', 'modepsY',   modepsY);
    assignin('base', 'modchiX',   modchiX);
    assignin('base', 'modvarrho', modvarrho_cfg);
    % Only override omegaX_val when it actually differs (e.g. near-closed economy)
    if abs(omegaX_cfg - omegaX_val) > 1e-12
        assignin('base', 'omegaX_val', omegaX_cfg);
    end

    try
        % Capture Dynare output
        diary_file = sprintf('dynare_comparison_%s.txt', configs(cfg).name);
        diary(diary_file);

        % Compute steady state for this configuration (inline, same as main_SOE.m)
        fprintf('Computing steady state for %s...\n', configs(cfg).description);

        % Set aliases expected by the SS routines (matching main_SOE.m conventions)
        epsilon    = epsilon_val;
        gamma      = gamma_val;
        psi        = 1;
        chi        = 1;
        Pistar_ss  = Pistar_ss_val;
        Rworld_ss  = Rworld_ss_val;
        sigmaH     = sigmaH_val;
        varrho_val = modvarrho_cfg;
        
        % Allow config-specific tb_target (e.g. for near-closed economy)
        if isfield(configs(cfg), 'tb_target_cfg') && ~isempty(configs(cfg).tb_target_cfg)
            tb_target_use = configs(cfg).tb_target_cfg;
        else
            tb_target_use = tb_target;
        end
        gammag_vec = modgammag;
        gammas_vec = modgammas;
        om_g       = ombar_val;
        om_s       = 1 - ombar_val;
        chiX_vec   = modchiX;
        omegaX     = omegaX_cfg;
        etastar    = etastar_val;
        Ystar      = ystar_ss_val;
        alpha_vec  = modalpha;
        alphaV_vec = modalphaV_cfg;
        beta_mat   = modbeta;
        epsY_vec   = modepsY;
        epsM_vec   = modepsM;
        A_vec      = modA;

        % Initial guesses
        pHvec_guess = 1.0*ones(nsec,1);
        C_guess     = 1.0;
        w_guess     = 1.0;
        Q_guess     = 1.0;
        x_guess     = [pHvec_guess; w_guess; Q_guess; C_guess];

        fprintf('  Testing residual with initial guess...\n');
        resid = steady_ntwsoe(x_guess, PVstar_ss, epsilon, varrho_val, sigmaH, ...
            gammag_vec, gammas_vec, om_g, om_s, chiX_vec, omegaX, etastar, Ystar, ...
            alpha_vec, alphaV_vec, beta_mat, epsY_vec, epsM_vec, gamma, chi, psi, A_vec, tb_target_use);
        fprintf('  Initial residual norm: %.6f\n', norm(resid));

        fsolve_opts = optimoptions('fsolve', 'TolFun', 1e-14, 'MaxIter', 2000, 'Display', 'off');
        [ss_sol, fval, exitflag] = fsolve(@(x_vec) steady_ntwsoe(x_vec, PVstar_ss, epsilon, varrho_val, sigmaH, ...
            gammag_vec, gammas_vec, om_g, om_s, chiX_vec, omegaX, etastar, Ystar, ...
            alpha_vec, alphaV_vec, beta_mat, epsY_vec, epsM_vec, gamma, chi, psi, A_vec, tb_target_use), ...
            x_guess, fsolve_opts);
        fprintf('  fsolve exit flag: %d. Final residual norm: %.6e\n', exitflag, norm(fval));
        
        % If first attempt failed, try randomized restarts
        if exitflag <= 0
            fprintf('  Trying randomized restarts for %s...\n', configs(cfg).description);
            rng(42);
            best_norm = norm(fval);
            best_sol  = ss_sol;
            best_flag = exitflag;
            for attempt = 1:20
                x_try = [abs(x_guess(1:nsec) .* (0.7 + 0.6*rand(nsec,1))); ...
                         0.8 + 0.4*rand(); 0.8 + 0.4*rand(); 0.8 + 0.4*rand()];
                [sol_try, fval_try, flag_try] = fsolve(@(x_vec) steady_ntwsoe(x_vec, PVstar_ss, epsilon, varrho_val, sigmaH, ...
                    gammag_vec, gammas_vec, om_g, om_s, chiX_vec, omegaX, etastar, Ystar, ...
                    alpha_vec, alphaV_vec, beta_mat, epsY_vec, epsM_vec, gamma, chi, psi, A_vec, tb_target_use), ...
                    x_try, fsolve_opts);
                if norm(fval_try) < best_norm
                    best_norm = norm(fval_try); best_sol = sol_try; best_flag = flag_try;
                end
                if flag_try > 0 && norm(fval_try) < 1e-8; break; end
            end
            ss_sol   = best_sol;
            fval     = fsolve(@(x_vec) steady_ntwsoe(x_vec, PVstar_ss, epsilon, varrho_val, sigmaH, ...
                gammag_vec, gammas_vec, om_g, om_s, chiX_vec, omegaX, etastar, Ystar, ...
                alpha_vec, alphaV_vec, beta_mat, epsY_vec, epsM_vec, gamma, chi, psi, A_vec, tb_target_use), ss_sol, fsolve_opts);
            exitflag = best_flag;
            fprintf('  Best restart residual norm: %.6e\n', norm(fval));
        end
        
        if norm(fval) > 1e-6
            error('fsolve did not converge for %s (residual %.2e).', configs(cfg).description, norm(fval));
        end

        pH_ss = ss_sol(1:nsec,1);
        w_ss  = ss_sol(nsec+1,1);
        Q_ss  = ss_sol(nsec+2,1);
        C_ss  = ss_sol(nsec+3,1);

        %% Evaluate the steady state (identical to main_SOE.m)
        PL_ss  = ones(nsec,1)*w_ss;
        PV_ss  = Q_ss*PVstar_ss;

        MCi_ss = (epsilon-1)/epsilon*pH_ss;
        PMi_ss = zeros(nsec,1);
        for ii = 1:nsec
            PMi_ss(ii) = (sum(beta_mat(ii,:).*(pH_ss'.^(1-epsM_vec(ii)))))^(1/(1-epsM_vec(ii)));
        end

        P_ss   = (varrho_val.^(sigmaH).*pH_ss.^(1-sigmaH) + (1-varrho_val).^(sigmaH).*PV_ss.^(1-sigmaH)).^(1/(1-sigmaH));
        p_g_ss = prod(P_ss.^gammag_vec);
        p_s_ss = prod(P_ss.^gammas_vec);

        C_g_ss  = om_g*C_ss/p_g_ss;
        C_s_ss  = om_s*C_ss/p_s_ss;
        C_gi_ss = gammag_vec.*(p_g_ss./P_ss).*C_g_ss;
        C_si_ss = gammas_vec.*(p_s_ss./P_ss).*C_s_ss;

        CHg_ss = varrho_val.^(sigmaH).*(pH_ss./P_ss).^(-sigmaH).*C_gi_ss;
        CHs_ss = varrho_val.^(sigmaH).*(pH_ss./P_ss).^(-sigmaH).*C_si_ss;
        CFg_ss = (1-varrho_val).^(sigmaH).*(PV_ss./P_ss).^(-sigmaH).*C_gi_ss;
        CFs_ss = (1-varrho_val).^(sigmaH).*(PV_ss./P_ss).^(-sigmaH).*C_si_ss;

        CHi_ss = CHg_ss + CHs_ss;
        CFi_ss = CFg_ss + CFs_ss;

        PX_ss  = prod(pH_ss.^chiX_vec);
        X_ss   = omegaX*(PX_ss/Q_ss)^(-etastar)*Ystar;
        Xi_ss  = chiX_vec.*X_ss.*PX_ss./pH_ss;

        % Intermediate use initial guess
        intermediate_use_guess = zeros(nsec,1);
        for i = 1:nsec
            for j = 1:nsec
                intermediate_use_guess(i) = intermediate_use_guess(i) + beta_mat(j,i);
            end
        end
        intermediate_use_guess = intermediate_use_guess .* mean(CHi_ss + Xi_ss);

        M_ss = (MCi_ss./PMi_ss).^epsY_vec.*alpha_vec.*(CHi_ss + Xi_ss + intermediate_use_guess);

        intermediate_use_ss = zeros(nsec,1);
        for i = 1:nsec
            for j = 1:nsec
                intermediate_use_ss(i) = intermediate_use_ss(i) + beta_mat(j,i) * (PMi_ss(j)/pH_ss(i))^epsM_vec(j) * M_ss(j);
            end
        end

        L_ss  = (MCi_ss./PL_ss).^epsY_vec.*(1-alpha_vec-alphaV_vec).*(CHi_ss + Xi_ss + intermediate_use_ss);
        Vi_ss = (MCi_ss./PV_ss).^epsY_vec.*alphaV_vec.*(CHi_ss + Xi_ss + intermediate_use_ss);
        Yi_ss = A_vec.*(alpha_vec.^(1./epsY_vec).*M_ss.^((epsY_vec-1)./epsY_vec) ...
            + alphaV_vec.^(1./epsY_vec).*Vi_ss.^((epsY_vec-1)./epsY_vec) ...
            + (1-alphaV_vec-alpha_vec).^(1./epsY_vec).*L_ss.^((epsY_vec-1)./epsY_vec)).^(epsY_vec./(epsY_vec-1));

        x_guess1 = [M_ss; L_ss; Vi_ss; Yi_ss];
        opts1    = optimoptions('fsolve','TolFun',1e-10,'Display','off');
        x_sol    = fsolve(@(x_vec) steady_ntwsoe_system(x_vec, alpha_vec, alphaV_vec, beta_mat, ...
            MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss, epsY_vec, epsM_vec, A_vec, pH_ss), ...
            x_guess1, opts1);
        fprintf('  Production system residual norm: %.6e\n', norm(steady_ntwsoe_system(x_sol, alpha_vec, alphaV_vec, beta_mat, MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss, epsY_vec, epsM_vec, A_vec, pH_ss)));

        M_ss  = x_sol(1:nsec);
        L_ss  = x_sol(nsec+1:2*nsec);
        Vi_ss = x_sol(2*nsec+1:3*nsec);
        Yi_ss = x_sol(3*nsec+1:4*nsec);

        V_ss   = sum(Vi_ss);
        CF_ss  = sum(CFi_ss);
        mkupV  = 1;
        IMP_ss = mkupV*(V_ss + CF_ss);

        TB_ss  = PX_ss*X_ss - PV_ss*IMP_ss;
        GDP_ss = C_ss + TB_ss;
        N_ss   = sum(L_ss);
        r_star_ss = Rworld_ss;
        pi_ss  = 1;
        r_ss   = 1/beta_val;
        tbgdp  = TB_ss/GDP_ss;

        Bstar_ss = -TB_ss / (Q_ss*(1 - r_star_ss/Pistar_ss));
        bbar_val = Q_ss*Bstar_ss/GDP_ss;

        Pistar_ss_val  = Pistar_ss;
        Rworld_ss_val  = Rworld_ss;
        PVstar_ss_val  = PVstar_ss;

        % Validate SS
        fprintf('  SS: GDP=%.4f  C=%.4f  TB/GDP=%.3f  Q=%.4f  IMP=%.4f\n', ...
            GDP_ss, C_ss, tbgdp, Q_ss, IMP_ss);
        if ~isfinite(GDP_ss) || ~isfinite(IMP_ss) || GDP_ss <= 0
            error('SS values are non-finite or non-positive. Steady-state solver failed to converge.');
        end

        % Create indexed steady state variables for Dynare initval block
        for i = 1:nsec
            eval(sprintf('Cgi_ss%d = C_gi_ss(%d);', i, i));
            eval(sprintf('Csi_ss%d = C_si_ss(%d);', i, i));
            eval(sprintf('CFg_ss%d = CFg_ss(%d);', i, i));
            eval(sprintf('CFs_ss%d = CFs_ss(%d);', i, i));
            eval(sprintf('CHg_ss%d = CHg_ss(%d);', i, i));
            eval(sprintf('CHs_ss%d = CHs_ss(%d);', i, i));
            eval(sprintf('Vi_ss%d = Vi_ss(%d);', i, i));
            eval(sprintf('PH_ss%d = pH_ss(%d);', i, i));
            eval(sprintf('MC_ss%d = MCi_ss(%d);', i, i));
            eval(sprintf('Y_ss%d = Yi_ss(%d);', i, i));
            eval(sprintf('L_ss%d = L_ss(%d);', i, i));
            eval(sprintf('P_ss%d = P_ss(%d);', i, i));
            eval(sprintf('PMi_ss%d = PMi_ss(%d);', i, i));
            eval(sprintf('Mi_ss%d = M_ss(%d);', i, i));
            eval(sprintf('PL_ss%d = PL_ss(%d);', i, i));
            eval(sprintf('Xi_ss%d = Xi_ss(%d);', i, i));
        end

        CF_ss         = sum(CFi_ss);
        CFg_total_ss  = sum(CFg_ss);
        CFs_total_ss  = sum(CFs_ss);
        M_tot_ss      = sum(M_ss);
        Y_tot_ss      = sum(Yi_ss);
        VA_ss         = sum(Yi_ss - M_ss);
        Ctotg_ss      = sum(gammag_vec.*(p_g_ss./P_ss).*C_g_ss);
        Ctots_ss      = sum(gammas_vec.*(p_s_ss./P_ss).*C_s_ss);
        Ctot_ss       = Ctotg_ss + Ctots_ss;
        M_tot_ss_val  = M_tot_ss;
        VA_ss_val     = VA_ss;
        Ctotg_ss_val  = Ctotg_ss;
        Ctots_ss_val  = Ctots_ss;
        Ctot_ss_val   = Ctot_ss;
        IMP_ss_val    = IMP_ss;

        save params_val_ul.mat;

        % Ensure ToT shock parameters are in params_val_ul.mat
        save('params_val_ul.mat', 'rho_pvstar_val', 'sigma_pvstar_val', 'shock_eps_pvstar_val', 'shock_epsA_val', '-append');

        fprintf('Running Dynare for %s...\n', configs(cfg).description);
        dynare NK_SOE_lev_gap2.mod noclearall;
        
        diary off;
        
        % Store results
        results(cfg).name = configs(cfg).name;
        results(cfg).description = configs(cfg).description;
        results(cfg).M_ = M_;
        results(cfg).oo_ = oo_;
        results(cfg).success = true;
        
        fprintf('✓ %s completed successfully\n', configs(cfg).description);
        
        % Save results
        output_filename = sprintf('model_comparison_%s.mat', configs(cfg).name);
        save(output_filename, 'oo_', 'M_', 'configs', '-v7.3');
        fprintf('✓ Saved: %s\n\n', output_filename);
        
    catch ME
        fprintf('✗ FAILED: %s\n', ME.message);
        fprintf('Stack:\n%s\n', ME.getReport('extended'));
        diary off;
        results(cfg).success = false;
        results(cfg).error = ME.message;
    end
end

%% Check if we have enough successful runs
n_success = sum([results.success]);
fprintf('\n════════════════════════════════════════════════\n');
fprintf('Results: %d/%d models completed successfully\n', n_success, length(configs));
fprintf('════════════════════════════════════════════════\n\n');

if n_success < 1
    fprintf('✗ No successful model runs. Cannot generate plots.\n');
    return;
end

%% Extract and compare IRFs
fprintf('Generating comparison plots...\n\n');

%% Extract IRF data and create plots
nPeriods = 40;

% Optional normalization: scale all IRFs so Q impact is a target percent change
% Keep FALSE to preserve the structural ToT shock size set by sigma_pvstar_val.
normalize_to_Q_impact = false;
target_Q_impact_pct = -10;   % e.g., -10 means a 10% fall in Q on impact

% Variables to compare (8 key variables) - using actual model variable names from NK_IOSOE_lev.mod
variables = {
    'GDP',    'log_dev',  'Output (GDP)';
    'C',      'log_dev',  'Consumption';
    'pi',     'ann_pct',  'CPI Inflation (annualized)';
    'r',      'ann_pct',  'Nominal Interest Rate (annualized)';
    'TB',     'level',    'Trade Balance (level dev.)';
    'Q',      'log_dev',  'Real Exchange Rate';
    'IMP',    'log_dev',  'Imports';
    'X',      'log_dev',  'Exports';
    'Ygap',   'log_dev',  'Output Gap (NK vs Flex-Price)';
    'GDPgap', 'log_dev',  'GDP Gap (NK vs Flex-Price)'
};

% Inflation comparison plot
inflation_vars = {
    'pi',   'ann_pct',  'Aggregate Inflation (annualized)';
    'pi_s', 'ann_pct',  'Services Inflation (annualized)';
    'pi_g', 'ann_pct',  'Goods Inflation (annualized)';
    'pi_e', 'ann_pct',  'Foreign Goods Inflation (annualized)'
};

% Model display info
model_names = {};
colors = {};
linestyles = {};

for m = 1:length(results)
    if results(m).success
        model_names{end+1} = results(m).name;
    end
end

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.9290 0.6940 0.1250]};
linestyles = {'-', '--', '-.', ':'};

%% Identify available shock - use 'GDP_' prefix (GDP has no sectoral variants)
fprintf('Searching for available shocks...\n');

selected_shock = '';

% First, prefer eps_pvstar (Terms of Trade shock)
for m = 1:length(results)
    if results(m).success && isfield(results(m).oo_, 'irfs')
        if isfield(results(m).oo_.irfs, 'GDP_eps_pvstar')
            selected_shock = 'eps_pvstar';
            break;
        end
    end
end

% Second fallback: scan for epsA shocks (manufacturing TFP)
if isempty(selected_shock)
    for m = 1:length(results)
        if results(m).success && isfield(results(m).oo_, 'irfs')
            irf_names = fieldnames(results(m).oo_.irfs);
            for k = 1:length(irf_names)
                fname = irf_names{k};
                if strncmp(fname, 'GDP_epsA_', 9)
                    selected_shock = fname(5:end);  % e.g., 'epsA_3'
                    break;
                end
            end
            if ~isempty(selected_shock), break; end
        end
    end
end

% Final fallback: any non-trivial GDP_ shock
if isempty(selected_shock)
    for m = 1:length(results)
        if results(m).success && isfield(results(m).oo_, 'irfs')
            irf_names = fieldnames(results(m).oo_.irfs);
            for k = 1:length(irf_names)
                fname = irf_names{k};
                if strncmp(fname, 'GDP_', 4)
                    selected_shock = fname(5:end);
                    break;
                end
            end
            if ~isempty(selected_shock), break; end
        end
    end
end

if isempty(selected_shock)
    fprintf('✗ No shocks with IRFs found in model outputs.\n');
    fprintf('Ensure stoch_simul is configured in NK_IOSOE_lev.mod\n');
    return;
end

fprintf('✓ Selected shock: %s\n\n', selected_shock);

%% Extract IRF data
irfs_data = struct();
sector_irfs = struct();

for m = 1:length(results)
    if ~results(m).success
        continue;
    end
    
    oo_ = results(m).oo_;
    M_ = results(m).M_;
    
    fprintf('Extracting IRFs for: %s\n', results(m).description);
    
    if ~isfield(oo_, 'irfs')
        fprintf('  ✗ No IRFs available\n');
        continue;
    end
    
    % Extract each variable
    for v = 1:size(variables, 1)
        var_name = variables{v, 1};
        var_type = variables{v, 2};
        display_name = variables{v, 3};
        
        % Build IRF field name
        irf_field = sprintf('%s_%s', var_name, selected_shock);
        
        if isfield(oo_.irfs, irf_field)
            irf_vals = oo_.irfs.(irf_field);
            
            % Convert to percentage
            if strcmp(var_type, 'log_dev')
                irf_vals = 100 * irf_vals;
            elseif strcmp(var_type, 'ann_pct')
                % Quarterly model: annualize by multiplying by 400
                irf_vals = 400 * irf_vals;
            elseif strcmp(var_type, 'level')
                % Keep raw level deviations (for variables near 0 at SS)
                irf_vals = irf_vals;
            end
            
            % Pad to nPeriods
            if length(irf_vals) < nPeriods
                irf_vals = [irf_vals; zeros(nPeriods - length(irf_vals), 1)];
            else
                irf_vals = irf_vals(1:nPeriods);
            end
            
            irfs_data.(sprintf('m%d_%s', m, var_name)) = irf_vals;
            fprintf('  ✓ %s\n', var_name);
        else
            fprintf('  ✗ %s not found\n', var_name);
        end
    end

    % Extract inflation variables (pi, pi_g, pi_s, pi_e)
    for v = 1:size(inflation_vars, 1)
        var_name = inflation_vars{v, 1};
        var_type = inflation_vars{v, 2};

        irf_field = sprintf('%s_%s', var_name, selected_shock);
        if isfield(oo_.irfs, irf_field)
            irf_vals = oo_.irfs.(irf_field);

            if strcmp(var_type, 'log_dev')
                irf_vals = 100 * irf_vals;
            elseif strcmp(var_type, 'ann_pct')
                irf_vals = 400 * irf_vals;
            elseif strcmp(var_type, 'level')
                irf_vals = irf_vals;
            end

            if length(irf_vals) < nPeriods
                irf_vals = [irf_vals; zeros(nPeriods - length(irf_vals), 1)];
            else
                irf_vals = irf_vals(1:nPeriods);
            end

            irfs_data.(sprintf('m%d_%s', m, var_name)) = irf_vals;
        end
    end

    % Extract sectoral IRFs (production, consumption of goods, labor demand)
    sector_sets = {
        'Y',  'log_dev',  'Production';
        'Cg', 'log_dev',  'Consumption of goods';
        'L',  'log_dev',  'Labor demand'
    };

    for s = 1:size(sector_sets, 1)
        var_prefix = sector_sets{s, 1};
        var_type = sector_sets{s, 2};
        data_mat = NaN(nPeriods, nsec);

        for i = 1:nsec
            irf_field = sprintf('%s_%d_%s', var_prefix, i, selected_shock);
            if isfield(oo_.irfs, irf_field)
                irf_vals = oo_.irfs.(irf_field);

                if strcmp(var_type, 'log_dev')
                    irf_vals = 100 * irf_vals;
                elseif strcmp(var_type, 'ann_pct')
                    irf_vals = 400 * irf_vals;
                end

                if length(irf_vals) < nPeriods
                    irf_vals = [irf_vals; zeros(nPeriods - length(irf_vals), 1)];
                else
                    irf_vals = irf_vals(1:nPeriods);
                end

                data_mat(:, i) = irf_vals;
            end
        end

        sector_irfs.(sprintf('m%d_%s', m, var_prefix)) = data_mat;
    end
    fprintf('\n');
end

% Normalize all IRFs by model so each one corresponds to target Q impact
if normalize_to_Q_impact
    fprintf('Normalizing IRFs to target Q impact of %.2f%%...\n', target_Q_impact_pct);
    for m = 1:length(results)
        if ~results(m).success
            continue;
        end

        q_key = sprintf('m%d_Q', m);
        if ~isfield(irfs_data, q_key)
            fprintf('  %s: Q IRF not found, skipping normalization\n', results(m).description);
            continue;
        end

        q_impact = irfs_data.(q_key)(1);
        if ~isfinite(q_impact) || abs(q_impact) < 1e-10
            fprintf('  %s: Q impact too small/invalid (%.3e), skipping normalization\n', ...
                results(m).description, q_impact);
            continue;
        end

        scale_factor = target_Q_impact_pct / q_impact;

        for v = 1:size(variables, 1)
            var_name = variables{v, 1};
            irf_key = sprintf('m%d_%s', m, var_name);
            if isfield(irfs_data, irf_key)
                irfs_data.(irf_key) = irfs_data.(irf_key) * scale_factor;
            end
        end

        for v = 1:size(inflation_vars, 1)
            var_name = inflation_vars{v, 1};
            irf_key = sprintf('m%d_%s', m, var_name);
            if isfield(irfs_data, irf_key)
                irfs_data.(irf_key) = irfs_data.(irf_key) * scale_factor;
            end
        end

        % Normalize sectoral IRFs too
        sector_prefixes = {'Y', 'Cg', 'L'};
        for s = 1:numel(sector_prefixes)
            s_key = sprintf('m%d_%s', m, sector_prefixes{s});
            if isfield(sector_irfs, s_key)
                sector_irfs.(s_key) = sector_irfs.(s_key) * scale_factor;
            end
        end

        fprintf('  %s: Q impact %.4f%% -> %.4f%% (scale x%.4f)\n', ...
            results(m).description, q_impact, target_Q_impact_pct, scale_factor);
    end
    fprintf('\n');
end

%% Create comparison figure
fprintf('Creating comparison figure (2x4 layout)...\n\n');

fig = figure('Name', 'Multi-Model Comparison', ...
    'NumberTitle', 'off', ...
    'Position', [50 50 1500 900], ...
    'Visible', 'on');
drawnow;

% Plot each variable
for v = 1:size(variables, 1)
    var_name = variables{v, 1};
    display_name = variables{v, 3};
    
    subplot(2, 4, v);
    hold on;
    grid on;
    grid minor;
    
    % Plot each model
    model_idx = 0;
    for m = 1:length(results)
        if results(m).success
            model_idx = model_idx + 1;
            
            irf_key = sprintf('m%d_%s', m, var_name);
            if isfield(irfs_data, irf_key)
                irf_vals = irfs_data.(irf_key);
                plot(0:nPeriods-1, irf_vals, ...
                    'Color', colors{model_idx}, ...
                    'LineStyle', linestyles{model_idx}, ...
                    'LineWidth', 2.5, ...
                    'DisplayName', model_names{model_idx});
            end
        end
    end
    
    % Add zero line
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    
    var_type_v = variables{v, 2};
    if strcmp(var_type_v, 'level')
        ylabel_str = sprintf('%s (level dev.)', display_name);
    elseif strcmp(var_type_v, 'ann_pct')
        ylabel_str = sprintf('%s (ann. %%, p.p.)', display_name);
    else
        ylabel_str = sprintf('%s (%% dev.)', display_name);
    end
    xlabel('Quarters', 'FontSize', 10);
    ylabel(ylabel_str, 'FontSize', 10);
    title(display_name, 'FontSize', 11, 'FontWeight', 'bold');
    
    if v == 4
        leg = legend('Location', 'best', 'FontSize', 8);
        leg.Box = 'on';
    end
    
    xlim([0, nPeriods-1]);
    enforce_min_yaxis(gca);
    hold off;
end

% Title
n_success = sum([results.success]);
sgtitle(sprintf('%d-Model Comparison (%d successful): %s Shock', length(configs), n_success, selected_shock),...
    'FontSize', 13, 'FontWeight', 'bold');

% Save figure
try
    print(fig, '-dpng', '-r300', 'model_comparison_irf.png');
    fprintf('✓ Figure saved: model_comparison_irf.png\n\n');
catch
    fprintf('⚠ Could not save figure\n\n');
end

%% Print summary table
fprintf('\n════════════════════════════════════════════════════════════════════\n');
fprintf('SUMMARY STATISTICS - 4 MODEL COMPARISON\n');
fprintf('Models: Baseline (Labor Costs) vs Flexible vs No-Labor vs Near-Closed Economy\n');
fprintf('════════════════════════════════════════════════════════════════════\n\n');

%% Inflation comparison figure
fprintf('Creating inflation comparison figure (2x2 layout)...\n\n');

fig = figure('Name', 'Inflation Comparison', ...
    'NumberTitle', 'off', ...
    'Position', [80 80 1200 800], ...
    'Visible', 'on');
drawnow;

for v = 1:size(inflation_vars, 1)
    var_name = inflation_vars{v, 1};
    display_name = inflation_vars{v, 3};

    subplot(2, 2, v);
    hold on;
    grid on;
    grid minor;

    model_idx = 0;
    for m = 1:length(results)
        if results(m).success
            model_idx = model_idx + 1;

            irf_key = sprintf('m%d_%s', m, var_name);
            if isfield(irfs_data, irf_key)
                irf_vals = irfs_data.(irf_key);
                plot(0:nPeriods-1, irf_vals, ...
                    'Color', colors{model_idx}, ...
                    'LineStyle', linestyles{model_idx}, ...
                    'LineWidth', 2.5, ...
                    'DisplayName', model_names{model_idx});
            end
        end
    end

    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
        'HandleVisibility', 'off');

    ylabel_str = sprintf('%s (ann. %%, p.p.)', display_name);
    xlabel('Quarters', 'FontSize', 10);
    ylabel(ylabel_str, 'FontSize', 10);
    title(display_name, 'FontSize', 11, 'FontWeight', 'bold');

    if v == 2
        leg = legend('Location', 'best', 'FontSize', 8);
        leg.Box = 'on';
    end

    xlim([0, nPeriods-1]);
    enforce_min_yaxis(gca);
    hold off;
end

sync_ylims(gcf);
sgtitle(sprintf('Inflation Responses: %s', selected_shock),...
    'FontSize', 12, 'FontWeight', 'bold');

try
    print(fig, '-dpng', '-r300', 'inflation_comparison_irf.png');
    fprintf('✓ Figure saved: inflation_comparison_irf.png\n\n');
catch
    fprintf('⚠ Could not save inflation figure\n\n');
end

%% Sectoral IRF figures
sector_figs = {
    'Y',  'Production by sector',           'Y_';
    'Cg', 'Consumption of goods by sector', 'Cg_';
    'L',  'Labor demand by sector',         'L_'
};

for sf = 1:size(sector_figs, 1)
    prefix = sector_figs{sf, 1};
    fig_title = sector_figs{sf, 2};

    fig = figure('Name', fig_title, ...
        'NumberTitle', 'off', ...
        'Position', [50 50 1600 1000], ...
        'Visible', 'on');
    drawnow;

    for i = 1:nsec
        subplot(3, 4, i);
        hold on;
        grid on;
        grid minor;

        model_idx = 0;
        for m = 1:length(results)
            if results(m).success
                model_idx = model_idx + 1;
                s_key = sprintf('m%d_%s', m, prefix);
                if isfield(sector_irfs, s_key)
                    series = sector_irfs.(s_key)(:, i);
                    if all(~isfinite(series))
                        continue;
                    end
                    plot(0:nPeriods-1, series, ...
                        'Color', colors{model_idx}, ...
                        'LineStyle', linestyles{model_idx}, ...
                        'LineWidth', 2.0, ...
                        'DisplayName', model_names{model_idx});
                end
            end
        end

        yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
            'HandleVisibility', 'off');

        if ~isempty(names)
            title_str = names{i};
        else
            title_str = sprintf('Sector %d', i);
        end
        title(title_str, 'FontSize', 9, 'FontWeight', 'bold');
        xlim([0, nPeriods-1]);

        if i == 1
            ylabel('Percent dev.', 'FontSize', 9);
        end

        if i == 4
            leg = legend('Location', 'best', 'FontSize', 7);
            leg.Box = 'on';
        end

        enforce_min_yaxis(gca);
        hold off;
    end

    sync_ylims(gcf);
    sgtitle(sprintf('%s: %s', fig_title, selected_shock),...
        'FontSize', 12, 'FontWeight', 'bold');

    out_name = sprintf('sector_irf_%s.png', prefix);
    try
        print(fig, '-dpng', '-r300', out_name);
        fprintf('✓ Figure saved: %s\n', out_name);
    catch
        fprintf('⚠ Could not save figure: %s\n', out_name);
    end
    fprintf('\n');
end

for v = 1:size(variables, 1)
    var_name = variables{v, 1};
    display_name = variables{v, 3};
    
    fprintf('%s\n', display_name);
    fprintf(repmat('=', 1, 110));
    fprintf('\n%-35s %12s %12s %12s\n', 'Model', 'Impact', 'Peak', 'Peak Qtr');
    fprintf(repmat('-', 1, 110));
    fprintf('\n');
    
    model_idx = 0;
    for m = 1:length(results)
        if results(m).success
            model_idx = model_idx + 1;
            irf_key = sprintf('m%d_%s', m, var_name);
            
            if isfield(irfs_data, irf_key)
                irf_vals = irfs_data.(irf_key);
                
                impact = irf_vals(1);
                [peak_abs, peak_idx] = max(abs(irf_vals));
                peak = irf_vals(peak_idx);
                
                fprintf('%-35s %12.4f %12.4f %12d\n', ...
                    model_names{model_idx}, impact, peak, peak_idx-1);
            end
        end
    end
    fprintf('\n');
end

fprintf('════════════════════════════════════════════════════════════════════\n\n');

%% =========================================================================
%%  LOCAL FUNCTIONS  (copied from main_tot_shock.m)
%% =========================================================================

function ss = tot_compute_ss(modbeta, beta_val, Pistar_ss, ...
        modgammag, modgammas, ombar_val, modchiX, omegaX_val, etastar_val, ...
        ystar_ss_val, modalpha, modalphaV, modepsY, modepsM, ...
        gamma_val, chi, psi, modA, tb_target, PVstar_ss, epsilon, ...
        modvarrho, sigmaH_val, nsec)
    % Mirrors main_SOE.m: calls steady_ntwsoe DIRECTLY (no penalty wrapper)
    % so fsolve can compute meaningful Jacobians.

    sigmaH     = sigmaH_val;
    om_g       = ombar_val;
    om_s       = 1 - ombar_val;
    varrho_val = modvarrho;
    beta_mat   = modbeta;
    epsY_vec   = modepsY;
    epsM_vec   = modepsM;
    alpha_vec  = modalpha;
    alphaV_vec = modalphaV;
    chiX_vec   = modchiX;
    gammag_vec = modgammag;
    gammas_vec = modgammas;
    omegaX     = omegaX_val;
    etastar    = etastar_val;
    Ystar      = ystar_ss_val;
    A_vec      = modA;

    % Multiple starting points (same as main_SOE.m uses all-ones)
    x0_list = [ ...
        [ones(nsec,1);       1.00; 1.00; 1.00], ...
        [ones(nsec,1)*0.90;  0.95; 1.05; 0.90], ...
        [ones(nsec,1)*1.10;  1.05; 0.95; 1.10], ...
        [ones(nsec,1)*0.75;  0.90; 1.10; 0.80]  ...
    ];
    % Call steady_ntwsoe directly — identical to main_SOE.m
    ss_fn = @(x) steady_ntwsoe(x, PVstar_ss, epsilon, ...
        varrho_val, sigmaH, gammag_vec, gammas_vec, om_g, om_s, ...
        chiX_vec, omegaX, etastar, Ystar, alpha_vec, alphaV_vec, ...
        beta_mat, epsY_vec, epsM_vec, gamma_val, chi, psi, A_vec, tb_target);

    opt = optimoptions('fsolve', 'TolFun', 1e-14, 'Display', 'off', ...
        'MaxFunctionEvaluations', 1e5, 'MaxIterations', 4e3);

    ss_sol    = [];
    best_fnorm = inf;
    best_ef    = -inf;
    for k = 1:size(x0_list, 2)
        x0 = x0_list(:, k);
        [ss_try, f_try, ef] = fsolve(ss_fn, x0, opt);
        f_norm = norm(f_try);
        fprintf('    SS attempt %d: exitflag=%d  residual=%.2e\n', k, ef, f_norm);
        if isfinite(f_norm) && (ef > best_ef || f_norm < best_fnorm)
            best_fnorm = f_norm;
            best_ef    = ef;
            ss_sol     = ss_try;
        end
        if ef > 0 && f_norm < 1e-8
            break;   % converged well enough
        end
    end

    if isempty(ss_sol) || best_fnorm > 1e-4
        error('tot_compute_ss:no_solution', ...
            'Steady-state solver failed (best residual=%.2e). Check parameters.', best_fnorm);
    end

    pH_ss  = ss_sol(1:nsec);
    w_ss   = ss_sol(nsec+1);
    Q_ss   = ss_sol(nsec+2);
    C_ss   = ss_sol(nsec+3);

    PL_ss  = ones(nsec,1)*w_ss;
    PV_ss  = Q_ss*PVstar_ss;
    MCi_ss = (epsilon-1)/epsilon*pH_ss;

    PMi_ss = zeros(nsec,1);
    for ii = 1:nsec
        PMi_ss(ii) = (sum(beta_mat(ii,:).*(pH_ss'.^(1-epsM_vec(ii)))))^(1/(1-epsM_vec(ii)));
    end
    PMi_ss = max(real(PMi_ss), 1e-10);

    P_ss    = (varrho_val.^sigmaH.*pH_ss.^(1-sigmaH) + ...
               (1-varrho_val).^sigmaH.*PV_ss.^(1-sigmaH)).^(1/(1-sigmaH));
    p_g_ss  = prod(P_ss.^gammag_vec);
    p_s_ss  = prod(P_ss.^gammas_vec);
    C_g_ss  = om_g*C_ss/p_g_ss;
    C_s_ss  = om_s*C_ss/p_s_ss;
    C_gi_ss = gammag_vec.*(p_g_ss./P_ss).*C_g_ss;
    C_si_ss = gammas_vec.*(p_s_ss./P_ss).*C_s_ss;
    CHg_ss  = varrho_val.^sigmaH.*(pH_ss./P_ss).^(-sigmaH).*C_gi_ss;
    CHs_ss  = varrho_val.^sigmaH.*(pH_ss./P_ss).^(-sigmaH).*C_si_ss;
    CFg_ss  = (1-varrho_val).^sigmaH.*(PV_ss./P_ss).^(-sigmaH).*C_gi_ss;
    CFs_ss  = (1-varrho_val).^sigmaH.*(PV_ss./P_ss).^(-sigmaH).*C_si_ss;
    PX_ss   = prod(pH_ss.^chiX_vec);
    X_ss    = omegaX*(PX_ss/Q_ss)^(-etastar)*Ystar;
    Xi_ss   = chiX_vec.*X_ss.*PX_ss./pH_ss;

    int_g = zeros(nsec,1);
    for i = 1:nsec
        for j = 1:nsec, int_g(i) = int_g(i) + beta_mat(j,i); end
    end
    int_g = int_g * max(mean(CHg_ss+CHs_ss+Xi_ss), 1e-10);
    M0 = (MCi_ss./PMi_ss).^epsY_vec .* alpha_vec .* max(CHg_ss+CHs_ss+Xi_ss+int_g, 1e-10);
    M0 = max(real(M0), 1e-10);
    int2 = zeros(nsec,1);
    for i = 1:nsec
        for j = 1:nsec
            int2(i) = int2(i) + beta_mat(j,i)*(PMi_ss(j)/pH_ss(i))^epsM_vec(j)*M0(j);
        end
    end
    scale_rhs = max(CHg_ss+CHs_ss+Xi_ss+int2, 1e-10);
    L0 = (MCi_ss./PL_ss).^epsY_vec .* (1-alpha_vec-alphaV_vec) .* scale_rhs;
    V0 = (MCi_ss./PV_ss).^epsY_vec .* alphaV_vec .* scale_rhs;
    L0 = max(real(L0), 1e-10);
    V0 = max(real(V0), 1e-10);
    Y0 = A_vec.*(alpha_vec.^(1./epsY_vec).*M0.^((epsY_vec-1)./epsY_vec) ...
         + alphaV_vec.^(1./epsY_vec).*V0.^((epsY_vec-1)./epsY_vec) ...
         + (1-alphaV_vec-alpha_vec).^(1./epsY_vec).*L0.^((epsY_vec-1)./epsY_vec)) ...
         .^(epsY_vec./(epsY_vec-1));
    Y0 = max(real(Y0), 1e-10);

    opt2 = optimoptions('fsolve','TolFun',1e-10,'Display','off', ...
        'MaxFunctionEvaluations', 1e5, 'MaxIterations', 4e3);
    x0_prod = [M0;L0;V0;Y0];
    x0_prod(~isfinite(x0_prod)|~isreal(x0_prod)) = 1e-4;
    x0_prod = max(real(x0_prod), 1e-10);

    % Call steady_ntwsoe_system directly — same as main_SOE.m
    sys_fn = @(x) steady_ntwsoe_system(x, alpha_vec, alphaV_vec, beta_mat, ...
        MCi_ss, PMi_ss, PL_ss, PV_ss, CHg_ss+CHs_ss, Xi_ss, ...
        epsY_vec, epsM_vec, A_vec, pH_ss);
    [x_sol, f_sol, ef_sys] = fsolve(sys_fn, x0_prod, opt2);
    fprintf('    Production system: exitflag=%d  residual=%.2e\n', ef_sys, norm(f_sol));
    if isempty(x_sol), x_sol = x0_prod; end

    M_ss  = max(real(x_sol(1:nsec)), 1e-10);
    L_ss  = max(real(x_sol(nsec+1:2*nsec)), 1e-10);
    Vi_ss = max(real(x_sol(2*nsec+1:3*nsec)), 1e-10);
    Yi_ss = max(real(x_sol(3*nsec+1:4*nsec)), 1e-10);

    V_ss   = sum(Vi_ss);
    CF_ss  = sum(CFg_ss + CFs_ss);
    IMP_ss = V_ss + CF_ss;
    TB_ss  = PX_ss*X_ss - PV_ss*IMP_ss;
    GDP_ss = C_ss + TB_ss;
    N_ss   = sum(L_ss);
    r_ss   = 1/beta_val;
    pi_ss  = 1;
    Bstar_ss  = -TB_ss / (Q_ss*(1 - (Pistar_ss/beta_val)/Pistar_ss));
    bbar_val  = Q_ss*Bstar_ss/GDP_ss;

    ss.pH_ss=pH_ss; ss.w_ss=w_ss; ss.Q_ss=Q_ss; ss.C_ss=C_ss;
    ss.P_ss=P_ss; ss.p_g_ss=p_g_ss; ss.p_s_ss=p_s_ss;
    ss.C_g_ss=C_g_ss; ss.C_s_ss=C_s_ss;
    ss.C_gi_ss=C_gi_ss; ss.C_si_ss=C_si_ss;
    ss.CHg_ss=CHg_ss; ss.CHs_ss=CHs_ss;
    ss.CFg_ss=CFg_ss; ss.CFs_ss=CFs_ss;
    ss.Vi_ss=Vi_ss; ss.MCi_ss=MCi_ss; ss.PMi_ss=PMi_ss;
    ss.Yi_ss=Yi_ss; ss.L_ss=L_ss; ss.M_ss=M_ss; ss.PL_ss=PL_ss;
    ss.Xi_ss=Xi_ss; ss.PX_ss=PX_ss; ss.X_ss=X_ss;
    ss.V_ss=V_ss; ss.CF_ss=CF_ss; ss.IMP_ss=IMP_ss;
    ss.TB_ss=TB_ss; ss.GDP_ss=GDP_ss; ss.N_ss=N_ss;
    ss.r_ss=r_ss; ss.pi_ss=pi_ss;
    ss.Bstar_ss=Bstar_ss; ss.bbar_val=bbar_val;
    ss.Pistar_ss_val=Pistar_ss; ss.Rworld_ss_val=Pistar_ss/beta_val;
    ss.PVstar_ss_val=PVstar_ss;
end

function tot_load_ss(ss, nsec)
    assignin('caller','pH_ss',         ss.pH_ss);
    assignin('caller','w_ss',          ss.w_ss);
    assignin('caller','Q_ss',          ss.Q_ss);
    assignin('caller','C_ss',          ss.C_ss);
    assignin('caller','P_ss',          ss.P_ss);
    assignin('caller','p_g_ss',        ss.p_g_ss);
    assignin('caller','p_s_ss',        ss.p_s_ss);
    assignin('caller','C_g_ss',        ss.C_g_ss);
    assignin('caller','C_s_ss',        ss.C_s_ss);
    assignin('caller','MCi_ss',        ss.MCi_ss);
    assignin('caller','PMi_ss',        ss.PMi_ss);
    assignin('caller','PL_ss',         ss.PL_ss);
    assignin('caller','Yi_ss',         ss.Yi_ss);
    assignin('caller','L_ss',          ss.L_ss);
    assignin('caller','M_ss',          ss.M_ss);
    assignin('caller','Vi_ss',         ss.Vi_ss);
    assignin('caller','C_gi_ss',       ss.C_gi_ss);
    assignin('caller','C_si_ss',       ss.C_si_ss);
    assignin('caller','CFg_ss',        ss.CFg_ss);
    assignin('caller','CFs_ss',        ss.CFs_ss);
    assignin('caller','CHg_ss',        ss.CHg_ss);
    assignin('caller','CHs_ss',        ss.CHs_ss);
    assignin('caller','PX_ss',         ss.PX_ss);
    assignin('caller','X_ss',          ss.X_ss);
    assignin('caller','V_ss',          ss.V_ss);
    assignin('caller','CF_ss',         ss.CF_ss);
    assignin('caller','IMP_ss',        ss.IMP_ss);
    assignin('caller','GDP_ss',        ss.GDP_ss);
    assignin('caller','N_ss',          ss.N_ss);
    assignin('caller','r_ss',          ss.r_ss);
    assignin('caller','pi_ss',         ss.pi_ss);
    assignin('caller','Bstar_ss',      ss.Bstar_ss);
    assignin('caller','bbar_val',      ss.bbar_val);
    assignin('caller','Pistar_ss_val', ss.Pistar_ss_val);
    assignin('caller','Rworld_ss_val', ss.Rworld_ss_val);
    assignin('caller','PVstar_ss_val', ss.PVstar_ss_val);
    for i = 1:nsec
        assignin('caller', sprintf('Cgi_ss%d',  i), ss.C_gi_ss(i));
        assignin('caller', sprintf('Csi_ss%d',  i), ss.C_si_ss(i));
        assignin('caller', sprintf('CFg_ss%d',  i), ss.CFg_ss(i));
        assignin('caller', sprintf('CFs_ss%d',  i), ss.CFs_ss(i));
        assignin('caller', sprintf('CHg_ss%d',  i), ss.CHg_ss(i));
        assignin('caller', sprintf('CHs_ss%d',  i), ss.CHs_ss(i));
        assignin('caller', sprintf('Vi_ss%d',   i), ss.Vi_ss(i));
        assignin('caller', sprintf('PH_ss%d',   i), ss.pH_ss(i));
        assignin('caller', sprintf('MC_ss%d',   i), ss.MCi_ss(i));
        assignin('caller', sprintf('Y_ss%d',    i), ss.Yi_ss(i));
        assignin('caller', sprintf('L_ss%d',    i), ss.L_ss(i));
        assignin('caller', sprintf('P_ss%d',    i), ss.P_ss(i));
        assignin('caller', sprintf('PMi_ss%d',  i), ss.PMi_ss(i));
        assignin('caller', sprintf('Mi_ss%d',   i), ss.M_ss(i));
        assignin('caller', sprintf('PL_ss%d',   i), ss.PL_ss(i));
        assignin('caller', sprintf('Xi_ss%d',   i), ss.Xi_ss(i));
    end
end

function residual = tot_safe_steady_residual(x, ...
    PVstar_ss, epsilon, varrho_val, sigmaH, gammag_vec, gammas_vec, ...
    om_g, om_s, chiX_vec, omegaX, etastar, Ystar, alpha_vec, alphaV_vec, ...
    beta_mat, epsY_vec, epsM_vec, gamma_val, chi, psi, A_vec, tb_target, nsec)
    penalty = 1e6 * ones(nsec+3,1);
    if any(~isfinite(x)) || ~isreal(x), residual = penalty; return; end
    x = max(real(x), 1e-10);
    residual = steady_ntwsoe(x, PVstar_ss, epsilon, ...
        varrho_val, sigmaH, gammag_vec, gammas_vec, om_g, om_s, ...
        chiX_vec, omegaX, etastar, Ystar, alpha_vec, alphaV_vec, ...
        beta_mat, epsY_vec, epsM_vec, gamma_val, chi, psi, A_vec, tb_target);
    if any(~isfinite(residual)) || ~isreal(residual)
        residual = penalty;
    else
        residual = real(residual);
    end
end

function residual = tot_safe_system_residual(x, ...
    alpha_vec, alphaV_vec, beta_mat, MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss, ...
    epsY_vec, epsM_vec, A_vec, pH_ss, nsec)
    penalty = 1e6 * ones(4*nsec,1);
    if any(~isfinite(x)) || ~isreal(x), residual = penalty; return; end
    x = max(real(x), 1e-10);
    residual = steady_ntwsoe_system(x, alpha_vec, alphaV_vec, beta_mat, ...
        max(MCi_ss,1e-10), max(PMi_ss,1e-10), max(PL_ss,1e-10), max(PV_ss,1e-10), ...
        max(CHi_ss,1e-10), max(Xi_ss,1e-10), epsY_vec, epsM_vec, A_vec, max(pH_ss,1e-10));
    if any(~isfinite(residual)) || ~isreal(residual)
        residual = penalty;
    else
        residual = real(residual);
    end
end
