% main_tot_shock.m
% =========================================================================
% Terms of Trade Shock Exercise
% Compares three model variants:
%   1. Baseline NK-IOSOE model (calibrated I-O matrix, sticky prices)
%   2. Flexible price model  (kappa = 0, costless price adjustment)
%   3. No Input-Output model (diagonal I-O matrix, sticky prices)
%
% For each variant Dynare computes IRFs to a 1-standard-deviation shock to
% PVstar (the foreign price of the imported bundle = terms of trade shock).
% A positive shock raises import prices, worsening the SOE's terms of trade.
%
% Figure shows: Output (GDP), Consumption, Inflation, Real Interest Rate,
%               Trade Balance / GDP, Real Exchange Rate.
% =========================================================================

clear all; close all; clc;

fprintf('Running script file: %s\n', mfilename('fullpath'));
which main_tot_shock -all
which steady_ntwsoe_system -all

restoredefaultpath;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultFigureWindowStyle','docked');

%% ========== DATA FILENAMES (defined early for path resolution) ==========
filename_industries = 'Stata_to_excel_few_industries_chile.xls';
filename_io         = 'IO_2021_chile.csv';
filename_fpa        = 'fpa_vector_few_industries_chile.csv';

%% ========== ADD PATHS ==========
pths = { ...
    'C:\Program Files\Dynare\6.4\matlab' ...
    'C:\Users\mgiarda\OneDrive - Banco Central de Chile\Escritorio\Github_BC\DPMP-DME-Modelos\Model\utils' ...
    'C:\Users\mgiarda\OneDrive - Banco Central de Chile\Escritorio\Proyectos\ntw_data_github' ...
    'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github' ...
    'C:\github\27-9-25\old files\network\DPMP-DME-Modelos\carpeta_Jose' ...
};
for ip = 1:numel(pths)
    if isfolder(pths{ip}), addpath(pths{ip}); end
end


%% ========== SECTOR / DATA SETUP ==========
nsec = 12;

varNames = {'kk','BEAName','Nameshort','dp','dy','dl','dri','share','industrytype',...
            'spend_good','spend_serv','alpha','alpha_V','var_rho'};
varTypes = {'double','char','char','double','double','double','double','double',...
            'categorical','double','double','double','double','double'};
opts = spreadsheetImportOptions('NumVariables',14,'VariableNames',varNames,...
                                'VariableTypes',varTypes,'DataRange','A2');

alldata    = readtable(filename_industries, opts);
alpha      = table2array(alldata(:,12));
alpha_V    = table2array(alldata(:,13));
var_rho    = table2array(alldata(:,14));
spend_good = table2array(alldata(:,10));
spend_serv = table2array(alldata(:,11));

% Data guards for numerical robustness
alpha(~isfinite(alpha))       = 0;
alpha_V(~isfinite(alpha_V))   = 0;
var_rho(~isfinite(var_rho))   = 0.5;
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

%% -------- I-O MATRIX (FULL) --------
betaio        = readmatrix(filename_io);
betaio        = betaio(1:nsec, 1:nsec);
betaio(~isfinite(betaio) | betaio < 0) = 0;
betaio_colsum = sum(betaio, 1);
betaio_colsum(betaio_colsum <= 0) = 1;
betax         = betaio ./ betaio_colsum;   % column-normalize
modbeta_full  = betax';                     % modbeta(i,j) = share of input i in sector j

%% ========== COMMON PARAMETERS (shared across variants) ==========
beta_val       = 0.986;


% Policy parameters
phi_val       = 2.5;
rhoi_val      = 0.6;
rhoirule_val  = 0.74;
ilabcosts_val = 0.1;
ombar_val     = 0.57;

% Set all non-ToT shocks to zero for isolated exercise
sigma_i_val     = 0.0;
rho_om1_val     = 0.1;   rho_om2_val   = 0;   sigma_om_val  = 0.0;
rho_tfp1_val    = 0.5;   rho_tfp2_val  = 0.0;
rho_val         = 0.1;
sigma_L_agg_val = 0.0;
isigma_tfp_val  = zeros(nsec,1);

% ---- TERMS OF TRADE SHOCK parameters ----
rho_pvstar_val   = 0.9;    % AR(1) persistence of PVstar
sigma_pvstar_val = 0.01;   % 1 std-dev = 1% shock to PVstar

% Production cost parameters
modcl    = zeros(nsec,1);
modclneg = zeros(nsec,1);
modcm    = zeros(nsec,1);

% Elasticities
modepsM   = ones(nsec,1) * 0.1;
modepsY   = ones(nsec,1) * 0.8;
modchiX   = ones(nsec,1) * 1/nsec;
modvarrho = var_rho;
modA      = ones(nsec,1);
modalpha  = alpha;
modalphaV = alpha_V;

% SOE parameters
Pistar_ss      = 1.00;
Rworld_ss      = Pistar_ss / beta_val;
kappaV_val     = 1e13;
epsilonV_val   = 1e13;
epsilonX_val   = 1;
omegaX_val     = 1;
chii_b_val     = 0.001;
etastar_val    = 1;
ystar_ss_val   = 1;
PVstar_ss      = 1;
sigmaH_val     = 0.999;
gamma_val      = 2;
gamma          = gamma_val;
psi            = 1;
chi            = 1;
epsilon        = 10;
tb_target      = -0.02;

% Price-stickiness (Calvo -> Rotemberg)
fpa            = readmatrix(filename_fpa);
fpa            = fpa(1:nsec);
fpa            = ones(nsec,1) * mean(fpa);
theta          = (1-fpa).^3;
kappa_calib    = theta*(epsilon) ./ ((1-theta).*(1-theta*beta_val));

% Consumption basket
goods      = spend_good > spend_serv;
services   = spend_serv > spend_good;
modgammag  = spend_good / sum(spend_good);
modgammas  = spend_serv / sum(spend_serv);


%% =========================================================================
%%  VARIANT 1 - BASELINE  (full IO, calibrated kappa)
%% =========================================================================
fprintf('\n=== VARIANT 1: BASELINE ===\n');

modbeta  = modbeta_full;
modkappa = kappa_calib;

ss_base = tot_compute_ss(modbeta, beta_val, Pistar_ss, modgammag, modgammas, ...
    ombar_val, modchiX, omegaX_val, etastar_val, ystar_ss_val, ...
    modalpha, modalphaV, modepsY, modepsM, gamma_val, chi, psi, ...
    modA, tb_target, PVstar_ss, epsilon, modvarrho, sigmaH_val, nsec);

tot_load_ss(ss_base, nsec);
save params_val.mat;

fprintf('Running Dynare (Baseline)...\n');
dynare NK_IOSOE_tot.mod noclearall nolog

irf_base = tot_extract_irfs(oo_, ss_base.GDP_ss);
oo_base  = oo_;
fprintf('  Baseline GDP_ss = %.4f,  TB/GDP = %.3f\n', ...
         ss_base.GDP_ss, ss_base.TB_ss/ss_base.GDP_ss);

%% =========================================================================
%%  VARIANT 2 - FLEXIBLE PRICES  (full IO, kappa = 0)
%% =========================================================================
fprintf('\n=== VARIANT 2: FLEXIBLE PRICES ===\n');

modkappa = zeros(nsec, 1);   % zero price-adjustment cost

% Steady state unchanged by kappa -- reuse ss_base
tot_load_ss(ss_base, nsec);
save params_val.mat;

fprintf('Running Dynare (Flexible Prices)...\n');
dynare NK_IOSOE_tot.mod noclearall nolog

irf_flex = tot_extract_irfs(oo_, ss_base.GDP_ss);
oo_flex  = oo_;

%% =========================================================================
%%  VARIANT 3 - NO IO  (diagonal IO matrix, calibrated kappa)
%% =========================================================================
fprintf('\n=== VARIANT 3: NO INPUT-OUTPUT (diagonal IO) ===\n');

% Build diagonal IO: off-diagonal elements -> 0, then renormalize columns
modbeta_diag = diag(diag(modbeta_full));
for jj = 1:nsec
    cs = sum(modbeta_diag(:, jj));
    if cs > 0
        modbeta_diag(:, jj) = modbeta_diag(:, jj) / cs;
    else
        modbeta_diag(jj, jj) = 1.0;
    end
end
modbeta  = modbeta_diag;
modkappa = kappa_calib;   % restore sticky prices

ss_noio = tot_compute_ss(modbeta, beta_val, Pistar_ss, modgammag, modgammas, ...
    ombar_val, modchiX, omegaX_val, etastar_val, ystar_ss_val, ...
    modalpha, modalphaV, modepsY, modepsM, gamma_val, chi, psi, ...
    modA, tb_target, PVstar_ss, epsilon, modvarrho, sigmaH_val, nsec);

tot_load_ss(ss_noio, nsec);
save params_val.mat;

fprintf('Running Dynare (No IO)...\n');
dynare NK_IOSOE_tot.mod noclearall nolog

irf_noio = tot_extract_irfs(oo_, ss_noio.GDP_ss);
oo_noio  = oo_;
fprintf('  No-IO GDP_ss = %.4f,  TB/GDP = %.3f\n', ...
         ss_noio.GDP_ss, ss_noio.TB_ss/ss_noio.GDP_ss);

%% =========================================================================
%%  FIGURE - IRF comparison
%% =========================================================================
fprintf('\nBuilding IRF comparison figure...\n');

T_plot  = 40;
periods = 1:T_plot;

col_base = [0.12 0.35 0.70];   % deep blue
col_flex = [0.85 0.18 0.10];   % red
col_noio = [0.10 0.60 0.25];   % green
lw = 2.5;

% Panel definitions: {irf_field, panel_title, y-label, display_scale}
panels = { ...
    'GDP',    'Output (GDP)',          '% dev. from SS',   100; ...
    'C',      'Consumption',           '% dev. from SS',   100; ...
    'pi',     'Inflation',             '% dev. from SS',   100; ...
    'r_real', 'Real Interest Rate',    'pp deviation',     100; ...
    'TB_GDP', 'Trade Balance / GDP',   'pp change',        100; ...
    'Q',      'Real Exchange Rate',    '% dev. from SS',   100; ...
};

fig = figure('Name','ToT Shock - IRF Comparison', ...
             'Position',[40 40 1440 840], 'Visible','on');

for pp = 1:6
    fld = panels{pp,1};   ttl = panels{pp,2};
    ylb = panels{pp,3};   scl = panels{pp,4};

    yb = scl * tot_pad(irf_base.(fld), T_plot);
    yf = scl * tot_pad(irf_flex.(fld), T_plot);
    yn = scl * tot_pad(irf_noio.(fld), T_plot);

    subplot(2, 3, pp);
    h1 = plot(periods, yb, '-',  'Color', col_base, 'LineWidth', lw);  hold on;
    h2 = plot(periods, yf, '--', 'Color', col_flex, 'LineWidth', lw);
    h3 = plot(periods, yn, ':',  'Color', col_noio, 'LineWidth', lw);
    plot(periods, zeros(1, T_plot), 'k-', 'LineWidth', 0.8);

    title(ttl, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ylb, 'FontSize', 10);
    xlabel('Quarters', 'FontSize', 10);
    xlim([1 T_plot]);   grid on;   box on;
    set(gca, 'FontSize', 10);
end

% Shared legend at the bottom
lgd = legend([h1 h2 h3], ...
    {'Baseline', 'Flexible Prices', 'No Input-Output'}, ...
    'Orientation', 'horizontal', 'FontSize', 12, 'Box', 'on');
lgd.Position = [0.23 0.005 0.54 0.04];

sgtitle({ ...
    'Impulse Responses to a Terms of Trade Shock', ...
    sprintf('1%% increase in import price P^{V*}  (\\rho_{ToT} = %.1f)', rho_pvstar_val)}, ...
    'FontSize', 13, 'FontWeight', 'bold');

saveas(fig, 'tot_shock_irf_comparison.png');
saveas(fig, 'tot_shock_irf_comparison.fig');
fprintf('Figure saved: tot_shock_irf_comparison.png / .fig\n');
fprintf('\n=== main_tot_shock.m finished ===\n');


%% =========================================================================
%%  LOCAL FUNCTIONS  (must appear after all script statements in MATLAB)
%% =========================================================================

% -------------------------------------------------------------------------
% tot_compute_ss: solve for steady state given an IO matrix and kappa config
% -------------------------------------------------------------------------
function ss = tot_compute_ss(modbeta, beta_val, Pistar_ss, ...
        modgammag, modgammas, ombar_val, modchiX, omegaX_val, etastar_val, ...
        ystar_ss_val, modalpha, modalphaV, modepsY, modepsM, ...
        gamma_val, chi, psi, modA, tb_target, PVstar_ss, epsilon, ...
        modvarrho, sigmaH_val, nsec)

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

    x0_list = [ ...
        [ones(nsec,1);      1.00; 1.00; 1.00], ...
        [ones(nsec,1)*0.90; 0.95; 1.05; 0.90], ...
        [ones(nsec,1)*1.10; 1.05; 0.95; 1.10], ...
        [ones(nsec,1)*0.75; 0.90; 1.10; 0.80]  ...
    ];
    opt = optimoptions('fsolve', ...
        'TolFun',1e-14, 'Display','off', ...
        'MaxFunctionEvaluations', 1e5, 'MaxIterations', 4e3);

    ss_sol = [];
    best_fnorm = inf;
    for k = 1:size(x0_list,2)
        x0 = x0_list(:,k);
        try
            [ss_try, f_try, exitflag] = fsolve(@(x) tot_safe_steady_residual(x, ...
                PVstar_ss, epsilon, varrho_val, sigmaH, gammag_vec, gammas_vec, ...
                om_g, om_s, chiX_vec, omegaX, etastar, Ystar, alpha_vec, alphaV_vec, ...
                beta_mat, epsY_vec, epsM_vec, gamma_val, chi, psi, A_vec, tb_target, nsec), ...
                x0, opt);

            f_norm = norm(f_try);
            if isfinite(f_norm) && f_norm < best_fnorm
                best_fnorm = f_norm;
                ss_sol = ss_try;
            end
            if exitflag > 0 && isfinite(f_norm)
                ss_sol = ss_try;
                break;
            end
        catch
        end
    end

    if isempty(ss_sol)
        warning('tot_compute_ss:outer_solve_failed', ...
            'No feasible outer steady-state solution found; using safe fallback initial state.');
        ss_sol = x0_list(:,1);
    end

    pH_ss = max(ss_sol(1:nsec), 1e-8);
    w_ss  = max(ss_sol(nsec+1), 1e-8);
    Q_ss  = max(ss_sol(nsec+2), 1e-8);
    C_ss  = max(ss_sol(nsec+3), 1e-8);

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
    CHi_ss  = CHg_ss + CHs_ss;
    CFi_ss  = CFg_ss + CFs_ss;
    PX_ss   = prod(pH_ss.^chiX_vec);
    X_ss    = omegaX*(PX_ss/Q_ss)^(-etastar)*Ystar;
    Xi_ss   = chiX_vec.*X_ss.*PX_ss./pH_ss;

    % Initial guess for production inputs
    int_g = zeros(nsec,1);
    for i = 1:nsec
        for j = 1:nsec, int_g(i) = int_g(i) + beta_mat(j,i); end
    end
    int_g = int_g * max(mean(CHi_ss + Xi_ss), 1e-10);
    M0 = (MCi_ss./PMi_ss).^epsY_vec .* alpha_vec .* max(CHi_ss + Xi_ss + int_g, 1e-10);
    M0 = max(real(M0), 1e-10);

    int2 = zeros(nsec,1);
    for i = 1:nsec
        for j = 1:nsec
            int2(i) = int2(i) + beta_mat(j,i)*(PMi_ss(j)/pH_ss(i))^epsM_vec(j)*M0(j);
        end
    end
    scale_rhs = max(CHi_ss+Xi_ss+int2, 1e-10);
    L0 = (MCi_ss./PL_ss).^epsY_vec .* (1-alpha_vec-alphaV_vec) .* scale_rhs;
    V0 = (MCi_ss./PV_ss).^epsY_vec .* alphaV_vec .* scale_rhs;
    L0 = max(real(L0), 1e-10);
    V0 = max(real(V0), 1e-10);
    Y0 = A_vec.*(alpha_vec.^(1./epsY_vec).*M0.^((epsY_vec-1)./epsY_vec) ...
         + alphaV_vec.^(1./epsY_vec).*V0.^((epsY_vec-1)./epsY_vec) ...
         + (1-alphaV_vec-alpha_vec).^(1./epsY_vec).*L0.^((epsY_vec-1)./epsY_vec)) ...
         .^(epsY_vec./(epsY_vec-1));
    Y0 = max(real(Y0), 1e-10);

    opt2  = optimoptions('fsolve','TolFun',1e-10,'Display','off', ...
        'MaxFunctionEvaluations', 1e5, 'MaxIterations', 4e3);
    x0_prod = [M0;L0;V0;Y0];
    x0_prod(~isfinite(x0_prod) | ~isreal(x0_prod)) = 1e-4;
    x0_prod = max(real(x0_prod), 1e-10);
    x0_prod_alt = [ones(nsec,1)*1e-2; ones(nsec,1)*1e-2; ones(nsec,1)*1e-2; ones(nsec,1)*1e-2];

    x_sol = [];
    for kk = 1:2
        try
            if kk == 1
                x_init = x0_prod;
            else
                x_init = x0_prod_alt;
            end
            [x_try, f_try, ef] = fsolve(@(x) tot_safe_system_residual(x, ...
                alpha_vec, alphaV_vec, beta_mat, MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss, ...
                epsY_vec, epsM_vec, A_vec, pH_ss, nsec), x_init, opt2);
            if ef > 0 && all(isfinite(f_try))
                x_sol = x_try;
                break;
            elseif isempty(x_sol) && all(isfinite(f_try))
                x_sol = x_try;
            end
        catch
        end
    end

    if isempty(x_sol)
        warning('tot_compute_ss:inner_solve_failed', ...
            'Production-system solve failed; using safe fallback quantities.');
        x_sol = x0_prod;
    end

    M_ss  = max(real(x_sol(1:nsec)), 1e-10);
    L_ss  = max(real(x_sol(nsec+1:2*nsec)), 1e-10);
    Vi_ss = max(real(x_sol(2*nsec+1:3*nsec)), 1e-10);
    Yi_ss = max(real(x_sol(3*nsec+1:4*nsec)), 1e-10);

    V_ss   = sum(Vi_ss);
    CF_ss  = sum(CFi_ss);
    IMP_ss = V_ss + CF_ss;
    TB_ss  = PX_ss*X_ss - PV_ss*IMP_ss;
    GDP_ss = C_ss + TB_ss;
    N_ss   = sum(L_ss);
    r_ss   = 1/beta_val;
    pi_ss  = 1;
    r_star_ss = Pistar_ss/beta_val;
    Bstar_ss  = -TB_ss / (Q_ss*(1 - r_star_ss/Pistar_ss));
    bbar_val  = Q_ss*Bstar_ss/GDP_ss;

    ss.pH_ss      = pH_ss;   ss.w_ss       = w_ss;
    ss.Q_ss       = Q_ss;    ss.C_ss       = C_ss;
    ss.P_ss       = P_ss;    ss.p_g_ss     = p_g_ss;
    ss.p_s_ss     = p_s_ss;  ss.C_g_ss     = C_g_ss;
    ss.C_s_ss     = C_s_ss;  ss.C_gi_ss    = C_gi_ss;
    ss.C_si_ss    = C_si_ss; ss.CHg_ss     = CHg_ss;
    ss.CHs_ss     = CHs_ss;  ss.CFg_ss     = CFg_ss;
    ss.CFs_ss     = CFs_ss;  ss.Vi_ss      = Vi_ss;
    ss.MCi_ss     = MCi_ss;  ss.PMi_ss     = PMi_ss;
    ss.Yi_ss      = Yi_ss;   ss.L_ss       = L_ss;
    ss.M_ss       = M_ss;    ss.PL_ss      = PL_ss;
    ss.Xi_ss      = Xi_ss;   ss.PX_ss      = PX_ss;
    ss.X_ss       = X_ss;    ss.V_ss       = V_ss;
    ss.CF_ss      = CF_ss;   ss.IMP_ss     = IMP_ss;
    ss.TB_ss      = TB_ss;   ss.GDP_ss     = GDP_ss;
    ss.N_ss       = N_ss;    ss.r_ss       = r_ss;
    ss.pi_ss      = pi_ss;   ss.Bstar_ss   = Bstar_ss;
    ss.bbar_val   = bbar_val;
    ss.Pistar_ss_val = Pistar_ss;
    ss.Rworld_ss_val = Pistar_ss/beta_val;
    ss.PVstar_ss_val = PVstar_ss;
end

% -------------------------------------------------------------------------
% tot_load_ss: push steady-state variables into the calling workspace
% -------------------------------------------------------------------------
function tot_load_ss(ss, nsec)
    assignin('caller','pH_ss',          ss.pH_ss);
    assignin('caller','w_ss',           ss.w_ss);
    assignin('caller','Q_ss',           ss.Q_ss);
    assignin('caller','C_ss',           ss.C_ss);
    assignin('caller','P_ss',           ss.P_ss);
    assignin('caller','p_g_ss',         ss.p_g_ss);
    assignin('caller','p_s_ss',         ss.p_s_ss);
    assignin('caller','C_g_ss',         ss.C_g_ss);
    assignin('caller','C_s_ss',         ss.C_s_ss);
    assignin('caller','MCi_ss',         ss.MCi_ss);
    assignin('caller','PMi_ss',         ss.PMi_ss);
    assignin('caller','PL_ss',          ss.PL_ss);
    assignin('caller','Yi_ss',          ss.Yi_ss);
    assignin('caller','L_ss',           ss.L_ss);
    assignin('caller','M_ss',           ss.M_ss);
    assignin('caller','Vi_ss',          ss.Vi_ss);
    assignin('caller','C_gi_ss',        ss.C_gi_ss);
    assignin('caller','C_si_ss',        ss.C_si_ss);
    assignin('caller','CFg_ss',         ss.CFg_ss);
    assignin('caller','CFs_ss',         ss.CFs_ss);
    assignin('caller','CHg_ss',         ss.CHg_ss);
    assignin('caller','CHs_ss',         ss.CHs_ss);
    assignin('caller','PX_ss',          ss.PX_ss);
    assignin('caller','X_ss',           ss.X_ss);
    assignin('caller','V_ss',           ss.V_ss);
    assignin('caller','CF_ss',          ss.CF_ss);
    assignin('caller','IMP_ss',         ss.IMP_ss);
    assignin('caller','GDP_ss',         ss.GDP_ss);
    assignin('caller','N_ss',           ss.N_ss);
    assignin('caller','r_ss',           ss.r_ss);
    assignin('caller','pi_ss',          ss.pi_ss);
    assignin('caller','Bstar_ss',       ss.Bstar_ss);
    assignin('caller','bbar_val',       ss.bbar_val);
    assignin('caller','Pistar_ss_val',  ss.Pistar_ss_val);
    assignin('caller','Rworld_ss_val',  ss.Rworld_ss_val);
    assignin('caller','PVstar_ss_val',  ss.PVstar_ss_val);

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

% -------------------------------------------------------------------------
% tot_extract_irfs: extract the eps_pvstar IRFs from oo_.irfs
% -------------------------------------------------------------------------
function irf = tot_extract_irfs(oo_, GDP_ss)
    shock = 'eps_pvstar';
    T     = 40;
    vars  = {'GDP','C','pi','r','TB','Q'};
    for v = 1:numel(vars)
        fname = [vars{v} '_' shock];
        if isfield(oo_.irfs, fname)
            raw           = oo_.irfs.(fname);
            irf.(vars{v}) = raw(1:min(T, numel(raw)));
        else
            irf.(vars{v}) = zeros(1, T);
            warning('tot_extract_irfs: field not found: %s', fname);
        end
    end
    % Real interest rate: r_dev - E[pi_{t+1}]  (approximate Fisher decomp.)
    r_v  = irf.r;
    pi_v = irf.pi;
    irf.r_real = r_v - [pi_v(2:end), 0];
    % Change in Trade Balance / GDP ratio: IRF(TB) / GDP_ss
    irf.TB_GDP = irf.TB / GDP_ss;
end

% -------------------------------------------------------------------------
% tot_safe_steady_residual: guards fsolve from NaN/Inf/complex evaluations
% -------------------------------------------------------------------------
function residual = tot_safe_steady_residual(x, ...
    PVstar_ss, epsilon, varrho_val, sigmaH, gammag_vec, gammas_vec, ...
    om_g, om_s, chiX_vec, omegaX, etastar, Ystar, alpha_vec, alphaV_vec, ...
    beta_mat, epsY_vec, epsM_vec, gamma_val, chi, psi, A_vec, tb_target, nsec)

    penalty = 1e6 * ones(nsec+3,1);

    if any(~isfinite(x)) || ~isreal(x)
        residual = penalty;
        return;
    end

    x = max(real(x), 1e-10);

    try
        residual = steady_ntwsoe(x, PVstar_ss, epsilon, ...
            varrho_val, sigmaH, gammag_vec, gammas_vec, om_g, om_s, ...
            chiX_vec, omegaX, etastar, Ystar, alpha_vec, alphaV_vec, ...
            beta_mat, epsY_vec, epsM_vec, gamma_val, chi, psi, A_vec, tb_target);

        if any(~isfinite(residual)) || ~isreal(residual)
            residual = penalty;
        else
            residual = real(residual);
        end
    catch
        residual = penalty;
    end
end

% -------------------------------------------------------------------------
% tot_safe_system_residual: guards inner production-system solver
% -------------------------------------------------------------------------
function residual = tot_safe_system_residual(x, ...
    alpha_vec, alphaV_vec, beta_mat, MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss, ...
    epsY_vec, epsM_vec, A_vec, pH_ss, nsec)

    penalty = 1e6 * ones(4*nsec,1);

    if any(~isfinite(x)) || ~isreal(x)
        residual = penalty;
        return;
    end

    x = max(real(x), 1e-10);

    try
        residual = steady_ntwsoe_system(x, alpha_vec, alphaV_vec, beta_mat, ...
            max(MCi_ss,1e-10), max(PMi_ss,1e-10), max(PL_ss,1e-10), max(PV_ss,1e-10), ...
            max(CHi_ss,1e-10), max(Xi_ss,1e-10), epsY_vec, epsM_vec, A_vec, max(pH_ss,1e-10));

        if any(~isfinite(residual)) || ~isreal(residual)
            residual = penalty;
        else
            residual = real(residual);
        end
    catch
        residual = penalty;
    end
end

% -------------------------------------------------------------------------
% tot_pad: pad or trim vector to length T
% -------------------------------------------------------------------------
function y = tot_pad(x, T)
    n = numel(x);
    if n >= T
        y = x(1:T);
    else
        y = [x, zeros(1, T-n)];
    end
end
