
% Based on: Ferrante, Graves & Iacoviello (2023) JME 140 S64-S81
% NK_SOE_lev_gap2: 12-sector SOE DSGE model for Chile (flexible-price + output-gap version)
% Chile sectors (12):
% 1 Agropecuario-silvícola y Pesca
% 2 Minería
% 3 Industria manufacturera
% 4 Electricidad, gas, agua y gestión de desechos
% 5 Construcción
% 6 Comercio, hoteles y restaurantes
% 7 Transporte, comunicaciones y servicios de información
% 8 Intermediación financiera
% 9 Servicios inmobiliarios y de vivienda
% 10 Servicios empresariales
% 11 Servicios personales
% 12 Administración pública

% beta_mat(i,j) = share of input j used by sector j (from I-O table)


clear all;
close all;
clc;

%% ========== EXERCISE SELECTOR ==========
% Choose which exercise to run:
% 0 = Baseline (all shocks active with default parameters)
% 1 = Exercise 1: Consumption preference shock only
% 2 = Exercise 2: Manufacturing TFP shock only
% 3 = Exercise 3: Monetary policy shock only

% When called from smm_estimation.m, EXERCISE is forced to 0 (Baseline)
% via the environment variable SMM_EXERCISE (survives clear all).
smm_ex_env = getenv('SMM_EXERCISE');
smm_called  = ~isempty(smm_ex_env);   % true when invoked from smm_estimation.m
if smm_called
    EXERCISE = str2double(smm_ex_env);
    set(0, 'DefaultFigureVisible', 'off');  % suppress all figures during SMM
else
    EXERCISE = 0;  % <<<< CHANGE THIS TO SELECT EXERCISE (0=Baseline matches SMM calibration)
end

exercise_labels = {'Baseline (all shocks)', 'Exercise 1: Preference shock', ...
    'Exercise 2: TFP shock to Manufacturing', 'Exercise 3: Monetary policy shock'};
fprintf('\n============================================================\n');
fprintf('  NK-SOE 12-sector model for Chile\n');
if smm_called
    fprintf('  [SMM mode] Exercise %d: %s\n', EXERCISE, exercise_labels{EXERCISE+1});
else
    fprintf('  %s\n', exercise_labels{EXERCISE+1});
end
fprintf('============================================================\n\n');
if ~ismember(EXERCISE, 0:3)
    error('Invalid EXERCISE value. Must be 0, 1, 2, or 3.');
end

restoredefaultpath
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultFigureWindowStyle','docked');

%% Add paths
pths = { ...
    'C:\Program Files\Dynare\6.4\matlab' ...
    'C:\Users\mgiarda\OneDrive - Banco Central de Chile\Escritorio\Github_BC\DPMP-DME-Modelos\Model\utils' ...
    'C:\Users\mgiarda\OneDrive - Banco Central de Chile\Escritorio\Proyectos\ntw_data_github' ...
};
for ip = 1:numel(pths)
    if isfolder(pths{ip})
        addpath(pths{ip});
    else
        warning('Path not found (skipped): %s', pths{ip});
    end
end

%% Basic setup
nsec = 12;

% Filenames for data (try to locate on user's PC)
default_files = {'Stata_to_excel_few_industries_chile.xls', 'IO_2021_chile.csv', 'fpa_vector_few_industries_chile.csv'};

% Candidate folders to search (add your personal paths here if needed)
repo_root = fileparts(fileparts(mfilename('fullpath')));
candidates = { ...
    pwd, ...                                % current working directory
    fullfile(repo_root, 'Process Data Codes'), ...
    fullfile(repo_root, 'Model', 'modelo_chile'), ...
    'C:\Users\mgiarda\OneDrive - Banco Central de Chile\Escritorio\Proyectos\ntw_data_github', ...
    'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github' ...
};

fullpaths = cell(size(default_files));
found = false(size(default_files));
for i = 1:numel(candidates)
    for j = 1:numel(default_files)
        if ~found(j)
            try
                cand = fullfile(candidates{i}, default_files{j});
            catch
                cand = '';
            end
            if ~isempty(cand) && isfile(cand)
                fullpaths{j} = cand;
                found(j) = true;
            end
        end
    end
end

% Also check paths relative to this script location
script_dir = fileparts(which('main_SOE_gap'));
for j = 1:numel(default_files)
    if ~found(j)
        cand = fullfile(script_dir, default_files{j});
        if isfile(cand)
            fullpaths{j} = cand;
            found(j) = true;
        end
    end
end

% Final fallback: check plain filenames in current folder
for j = 1:numel(default_files)
    if ~found(j) && isfile(default_files{j})
        fullpaths{j} = fullfile(pwd, default_files{j});
        found(j) = true;
    end
end

if any(~found)
    missing = default_files(~found);
    msg = sprintf('Missing required data files: %s\nTried these folders:\n', strjoin(missing, ', '));
    for i = 1:numel(candidates)
        msg = [msg sprintf('  %s\n', candidates{i})];
    end
    error(msg);
end

filename_industries = fullpaths{1};
filename_io         = fullpaths{2};
filename_fpa        = fullpaths{3};

%% READ DATA (12 sectors)
varNames = {'kk','BEAName','Nameshort','dp','dy','dl','dri','share','industrytype','spend_good','spend_serv','alpha','alpha_V', 'var_rho'};
varTypes = {'double','char','char','double','double','double','double','double','categorical','double','double','double','double','double'};
opts = spreadsheetImportOptions('NumVariables',14,'VariableNames',varNames,'VariableTypes',varTypes,'DataRange','A2');

alldata    = readtable(filename_industries, opts);
names      = table2array(alldata(:,3));
%share      = table2array(alldata(:,8));
p_d        = table2array(alldata(:,4));
y_d        = table2array(alldata(:,5));
l_d        = table2array(alldata(:,6));
ri_d       = table2array(alldata(:,7));
spend_good = table2array(alldata(:,10));
spend_serv = table2array(alldata(:,11));
alpha      = table2array(alldata(:,12));
alpha_V    = table2array(alldata(:,13));
var_rho    = table2array(alldata(:,14));

%% SET PARAMETERS FROM DATA
betaio = readmatrix(filename_io);
betaio = betaio(1:nsec,1:nsec);
% Normalize columns: each column j shows input shares for sector j
% Column sum = 1 means all intermediate inputs sum to 100% of sector j's intermediate bundle
betax  = betaio ./ sum(betaio,1);
% Transpose so beta_mat(i,j) = share of input i used by sector j
modbeta  = betax';      % I-O matrix: rows=supplying product, cols=purchasing sector
modalpha = alpha;  % intermediates share (vector) - will be modified
modalphaV = alpha_V; %import share 




% Override discount factor for SOE
beta_val = 0.986;  % Quarterly discount factor (annualized: ~5.7% discount rate)

% Price rigidities
theta = readmatrix(filename_fpa);  % Frequency of price adjustment (quarterly)
kappa = theta*(10-1)./((1-theta).*(1-theta*beta_val)); % convert to Rotemberg adjustment cost parameter
modkappa = kappa;  % Price adjustment costs by sector (κ_i)

% Spending categories
goods    = spend_good > spend_serv;  % Indicator: 1 if sector produces goods
services = spend_serv > spend_good;  % Indicator: 1 if sector produces services
gammag   = spend_good / sum(spend_good);  % Share of each good in goods consumption basket (γ^g_i)
modgammag = gammag;
gammas    = spend_serv / sum(spend_serv);  % Share of each service in services consumption basket (γ^s_i)
modgammas = gammas;



% Cost parameters
modcl     = zeros(nsec,1);  % Labor adjustment costs (set to zero)
modclneg  = zeros(nsec,1);  % Negative labor adjustment costs (set to zero)
modcm     = zeros(nsec,1);  % Material adjustment costs (set to zero)

% Elasticities
modepsM = ones(nsec,1)*0.1;  % Substitution elasticity for intermediate inputs (ε^m_i)
modepsY = ones(nsec,1)*0.8;  % Substitution elasticity in production (ε^Y_i)

% Labor market parameters
modpsil = ones(nsec,1)*(-1000);  % Inverse of sectoral labor adjustment costs (ψ^l_i, large negative = no adjustment)
modpsim = ones(nsec,1)*(-1000);  % Inverse of material adjustment costs (ψ^m_i, large negative = no adjustment)

% Policy parameters (Monetary Policy Rule)
phi_val        = 2.5;            % Taylor rule coefficient on inflation (φ_π) - response to inflation deviations
rhoi_val       = 0.6;            % Interest rate persistence in Taylor rule (ρ_i) - autoregressive parameter
rhoirule_val   = 0.74;           % Interest rate smoothing parameter in Taylor rule (alternative specification)

% Structural parameters
gammaind_val   = 0;              % Sectoral price indexation parameter (0 = no indexation to past inflation)
ilabcosts_val  = 0.1;            % Inverse of aggregate labor adjustment cost parameter (1/ψ^L)
ombar_val      = 0.57;           % Steady-state share of goods in total consumption basket (ω̄)

% ========== SHOCK PARAMETERS ==========
% Set shocks based on selected exercise
if EXERCISE == 0
    % BASELINE: All shocks active with default parameters
    sigma_i_val    = 0.001;          % Monetary policy shock
    rho_om1_val    = 0.1;            % Preference shock persistence
    rho_om2_val    = 0;
    sigma_om_val   = 0.0001;         % Preference shock std dev
    rho_tfp1_val   = 0.5;            % TFP shock persistence
    rho_tfp2_val   = 0.0;
    tfpshock       = zeros(nsec,1);
    isigma_tfp_val = 0.05*ones(nsec,1);  % TFP shocks for all sectors
    
elseif EXERCISE == 1
    % EXERCISE 1: Preference shock only
    sigma_i_val    = 0.0;            % NO monetary shock
    rho_om1_val    = 0.95;           % High persistence for preference shock
    rho_om2_val    = 0;
    sigma_om_val   = 0.01;           % PREFERENCE SHOCK ACTIVATED
    rho_tfp1_val   = 0.5;
    rho_tfp2_val   = 0.0;
    tfpshock       = zeros(nsec,1);
    isigma_tfp_val = zeros(nsec,1);  % NO TFP shocks
    
elseif EXERCISE == 2
    % EXERCISE 2: Manufacturing TFP shock only
    manufacturing_sector = 3;  % Sector 3 = Manufacturing
    sigma_i_val    = 0.0;            % NO monetary shock
    rho_om1_val    = 0.1;
    rho_om2_val    = 0;
    sigma_om_val   = 0.0;            % NO preference shock
    rho_tfp1_val   = 0.95;           % High persistence for TFP shock
    rho_tfp2_val   = 0.0;
    tfpshock       = zeros(nsec,1);
    tfpshock(manufacturing_sector) = 0.01;  % TFP shock to manufacturing
    isigma_tfp_val = zeros(nsec,1);
    isigma_tfp_val(manufacturing_sector) = 0.01;  % ONLY MANUFACTURING GETS TFP SHOCK
    
elseif EXERCISE == 3
    % EXERCISE 3: Monetary policy shock only
    sigma_i_val    = 0.01;           % MONETARY SHOCK ACTIVATED
    rho_om1_val    = 0.1;
    rho_om2_val    = 0;
    sigma_om_val   = 0.0;            % NO preference shock
    rho_tfp1_val   = 0.5;
    rho_tfp2_val   = 0.0;
    tfpshock       = zeros(nsec,1);
    isigma_tfp_val = zeros(nsec,1);  % NO TFP shocks
end

% Common shock parameters (apply to all exercises)
% Monetary policy shock (when active)
% (sigma_i_val set above based on exercise)

% Consumption preference shock (goods vs. services) (when active)
% (rho_om1_val, rho_om2_val, sigma_om_val set above based on exercise)

% TFP shocks (sectoral productivity) (when active)
% (rho_tfp1_val, rho_tfp2_val, tfpshock, isigma_tfp_val set above based on exercise)


% Labor shocks (aggregate and sectoral)
rho_val          = 0.1;              % Persistence parameter for labor supply shocks (ρ_L)
sigma_L_agg_val  = 0.0;           % Standard deviation of aggregate labor shock (σ_L) - common shock to all sectors
sigma_L_het      = zeros(nsec,1);    % Standard deviations of sector-specific labor shocks (σ_{L,i})
sigma_L_het(services) = 0.0;         % Services sectors have no idiosyncratic labor shocks (only aggregate)

% Export and home bias parameters
modchiX    = ones(nsec,1)*1/nsec;  % Export shares by sector (χ^X_i) - uniform distribution across sectors
modvarrho  = var_rho;               % Home bias parameter in consumption (ϱ_i) - Armington weight on domestic goods
modA       = ones(nsec,1);          % Steady-state TFP levels by sector (Ā_i) - normalized to 1

%% SOE-specific parameters
% World interest rate must satisfy Euler equation in steady state
% From Euler: r_star = Pistar_ss/beta in SS
% And from risk premium: r_star = Rworld (when debt is at target)
% Therefore: Rworld = Pistar_ss/beta
Pistar_ss = 1.00;            % Steady-state foreign inflation (Π*) - zero inflation
Rworld_ss = Pistar_ss/beta_val;  % Steady-state world interest rate derived from Euler equation

% Import elasticities (extreme values = nearly fixed quantities)
kappaV_val     = 10000000000000;  % Import adjustment cost (κ^V) - very large = imports nearly fixed
epsilonV_val   = 10000000000000;  % Import demand elasticity (ε^V) - very large = fixed import demand
epsilonX_val   = 1;
omegaX_val     = 1;               % Export demand scale parameter (Ω^X) - shifts export demand

% Debt and risk premium parameters
chii_b_val     = .001;            % Debt-elastic risk premium parameter (χ_b) - sensitivity of interest rate to debt
etastar_val    = 3.5;               % Foreign demand elasticity (η*) - substitution between domestic exports

% Foreign economy parameters
xi_rstar_val   = 0.2;             % Persistence of world interest rate shocks (ξ_{R*})
ystar_ss_val   = 1;               % Steady-state foreign output (Y*) - normalized to 1
PVstar_ss      = 1;               % Steady-state foreign price level (P^V*) - normalized to 1
sigmaH_val     = 0.999;           % Armington aggregation elasticity (σ^H) - substitution between home/foreign goods



% Household preference parameters
gamma_val      = 2;               % Risk aversion / inverse of intertemporal elasticity of substitution (γ)
gamma          = gamma_val;       % (duplicate for compatibility)
psi            = 1;               % Inverse of Frisch elasticity of labor supply (ψ) - labor supply rigidity
chi            = 1;               % Labor disutility weight (χ) - scaling parameter for labor disutility
epsilon        = 10;              % Elasticity of substitution between goods varieties (ε) - determines markups

% Initial guesses for steady state (in levels, so must be positive)
pHvec_guess = 1.0*ones(nsec,1);  % Prices around 1
C_guess  = 1.0;                   % Consumption around 1
w_guess  = 1.0;                   % Wage around 1
Q_guess  = 1.0;                   % Exchange rate around 1

dm_file = fullfile(fileparts(mfilename('fullpath')), 'data_moments_chile.mat');
tb_target = 0.02;
if exist(dm_file, 'file')
    tmp_dm = load(dm_file, 'dm_chile');
    if isfield(tmp_dm, 'dm_chile') && isfield(tmp_dm.dm_chile, 'd_TBGDP') && isfinite(tmp_dm.dm_chile.d_TBGDP)
        tb_target = tmp_dm.dm_chile.d_TBGDP;
    end
end

calib_shares=0;

CHshare_target=ones(nsec,1)*0.98;
CHshare_target(2)=0.98;
X_share_target = 0.2;

options = optimoptions('fsolve','TolFun',1e-14);

% ---- Load SMM estimates (override defaults when available) ----
% smm_estimates.mat is updated by smm_estimation.m during CMA-ES whenever
% a better calibration is found, and overwritten at the end with theta_hat.
% Applied for all exercises.
smm_est_file = fullfile(fileparts(mfilename('fullpath')), 'smm_estimates.mat');
smm_ckpt_file = fullfile(fileparts(mfilename('fullpath')), 'smm_best_so_far.mat');
param_names = {'ilabcosts', 'epsY', 'epsM', 'kappaV', 'rho_om1', 'sigma_om', 'rho_tfp1', 'isigma_tfp(1)', 'rho_pvstar', 'sigma_pvstar'};
if exist(smm_est_file, 'file')
    load(smm_est_file);   % loads: ilabcosts_val, modepsY, modepsM, kappaV_val,
                          %        rho_om1_val, sigma_om_val, rho_tfp1_val,
                          %        rho_pvstar_val, sigma_pvstar_val,
                          %        rho_psi_val, sigma_psi_val, isigma_tfp_val
    smm_param_source = 'smm_estimates.mat';
elseif exist(smm_ckpt_file, 'file')
    tmp_ckpt = load(smm_ckpt_file, 'smm_best_so_far');
    theta_best = tmp_ckpt.smm_best_so_far.theta_best;
    ilabcosts_val  = theta_best(1);
    modepsY        = ones(nsec,1) * theta_best(2);
    modepsM        = ones(nsec,1) * theta_best(3);
    kappaV_val     = exp(theta_best(4));
    rho_om1_val    = theta_best(5);
    sigma_om_val   = theta_best(6);
    rho_tfp1_val   = theta_best(7);
    isigma_tfp_val = theta_best(8:19);
    rho_pvstar_val   = theta_best(20);
    sigma_pvstar_val = theta_best(21);
    rho_xi_val       = theta_best(22);
    sigma_xi_val     = theta_best(23);
    smm_param_source = sprintf('smm_best_so_far.mat (obj=%.6f)', tmp_ckpt.smm_best_so_far.obj_best);
else
    smm_param_source = 'hard-coded defaults';
end

fprintf('--- Parameters (%s) ---\n', smm_param_source);
param_vals = [ilabcosts_val, modepsY(1), modepsM(1), kappaV_val, rho_om1_val, sigma_om_val, rho_tfp1_val, isigma_tfp_val(1), rho_pvstar_val, sigma_pvstar_val];
for i = 1:length(param_names)
    fprintf('  %-16s %g\n', param_names{i}, param_vals(i));
end
fprintf('\n');

%% Compute Steady State
fprintf('--- Steady State ---\n');
if calib_shares==0
    x_guess=[pHvec_guess;w_guess;Q_guess;C_guess];
    sigmaH = sigmaH_val; varrho_val = modvarrho;
    gammag_vec = modgammag; gammas_vec = modgammas;
    om_g=ombar_val; om_s=1-ombar_val; chiX_vec = modchiX; omegaX = omegaX_val;
    etastar = etastar_val; Ystar = ystar_ss_val; alpha_vec = modalpha;
    alphaV_vec = modalphaV; beta_mat = modbeta; epsY_vec = modepsY;
    epsM_vec = modepsM; A_vec = modA;

    [ss_sol, fval, exitflag] = fsolve(@(x_vec) steady_ntwsoe(x_vec,PVstar_ss,epsilon,varrho_val,sigmaH,...
        gammag_vec,gammas_vec,om_g,om_s,chiX_vec,omegaX,etastar,Ystar,...
        alpha_vec,alphaV_vec,beta_mat,epsY_vec,epsM_vec,gamma,chi,psi,A_vec,tb_target),x_guess,options );

    if exitflag <= 0
        warning('fsolve did not converge (exit flag %d). Model may fail.', exitflag);
    end

    pH_ss      = ss_sol(1:nsec,1);
    w_ss       = ss_sol(nsec+1,1);
    Q_ss       = ss_sol(nsec+2,1);
    C_ss       = ss_sol(nsec+3,1);
end 

if calib_shares==1
    x_guess=[pHvec_guess;w_guess;Q_guess;C_guess;modvarrho;omegaX_val];

    sigmaH = sigmaH_val; 
    gammag_vec = modgammag; gammas_vec = modgammas;
    om_g=ombar_val; om_s=1-ombar_val; chiX_vec = modchiX; omegaX = omegaX_val;
    etastar = etastar_val; Ystar = ystar_ss_val; alpha_vec = modalpha;
    alphaV_vec = modalphaV; beta_mat = modbeta; epsY_vec = modepsY;
    epsM_vec = modepsM; A_vec = modA;

    resid=steady_ntwsoe_calib(x_guess,PVstar_ss,epsilon,sigmaH,...
        gammag_vec,gammas_vec,om_g,om_s,chiX_vec,etastar,Ystar,...
        alpha_vec,alphaV_vec,beta_mat,epsY_vec,epsM_vec,...
        gamma,chi,psi,A_vec,tb_target,CHshare_target,X_share_target);

    ss_sol=fsolve(@(x_vec) steady_ntwsoe_calib(x_vec,PVstar_ss,epsilon,sigmaH,...
        gammag_vec,gammas_vec,om_g,om_s,chiX_vec,etastar,Ystar,...
        alpha_vec,alphaV_vec,beta_mat,epsY_vec,epsM_vec,...
        gamma,chi,psi,A_vec,tb_target,CHshare_target,X_share_target),x_guess );

    pH_ss      = ss_sol(1:nsec,1);
    w_ss       = ss_sol(nsec+1,1);
    Q_ss       = ss_sol(nsec+2,1);
    C_ss       = ss_sol(nsec+3,1);
    varrho_val = ss_sol(nsec+4:2*nsec+3,1);
    omegaX     = ss_sol(2*nsec+4,1);
    
    % Update parameter vectors with solved values
    modvarrho  = varrho_val;
end

%% Evaluate the steady state
om_g = ombar_val;
om_s = 1-ombar_val;

PL_ss   = ones(nsec,1)*w_ss;
PV_ss   = Q_ss*PVstar_ss;

MCi_ss  = (epsilon-1)/epsilon*pH_ss; % Marginal costs given prices
PMi_ss  = zeros(nsec,1);
for ii = 1:nsec
    PMi_ss(ii) = (sum(beta_mat(ii,:).* (pH_ss'.^(1-modepsM(ii)))))^(1/(1-modepsM(ii)));
end

P_ss    = (varrho_val.^(sigmaH).*pH_ss.^(1-sigmaH)+(1-varrho_val).^(sigmaH).*PV_ss.^(1-sigmaH)).^(1/(1-sigmaH)); % price of the good by the side of consumers

p_g_ss  = prod(P_ss.^gammag_vec); %multiply Goods price index
p_s_ss  = prod(P_ss.^gammas_vec); %multiply Services price index

C_g_ss  = om_g*C_ss/p_g_ss; %consumption of goods
C_s_ss  = om_s*C_ss/p_s_ss; % consumption of services

C_gi_ss = gammag_vec.*(p_g_ss./P_ss).*C_g_ss; % consumption of each good i if it is in goods' set
C_si_ss = gammas_vec.*(p_s_ss./P_ss).*C_s_ss; % consumption of each good i if it is in services' set

CHg_ss = varrho_val.^(sigmaH).*(pH_ss./P_ss).^(-sigmaH).*C_gi_ss; % consumption of good i in goods if home produced
CHs_ss = varrho_val.^(sigmaH).*(pH_ss./P_ss).^(-sigmaH).*C_si_ss; % consumption of good i in services if home produced

CFg_ss = (1-varrho_val).^(sigmaH).*(PV_ss./P_ss).^(-sigmaH).*C_gi_ss;  % consumption of good i in goods if foreign produced
CFs_ss = (1-varrho_val).^(sigmaH).*(PV_ss./P_ss).^(-sigmaH).*C_si_ss;  % consumption of good i in services if foreign produced

CHi_ss = CHg_ss+CHs_ss; % home demand for good i
CFi_ss = CFg_ss+CFs_ss; % foreign demand for good i

PX_ss  = prod(pH_ss.^chiX_vec); % Exports price index

X_ss   = omegaX*(PX_ss/Q_ss)^(-etastar)*Ystar; % Total exports

Xi_ss  = chiX_vec.*X_ss.*PX_ss./pH_ss; % sector i exports

% Compute intermediate use correctly: sum_j beta(j,i) * (PM(j)/PH(i))^epsM(j) * M(j)
% First compute initial guess
intermediate_use_guess = zeros(nsec,1);
for i = 1:nsec
    for j = 1:nsec
        % Use a simple approximation for initial guess
        intermediate_use_guess(i) = intermediate_use_guess(i) + beta_mat(j,i);
    end
end
intermediate_use_guess = intermediate_use_guess .* mean(CHi_ss + Xi_ss);

% Solve for production inputs
M_ss   = (MCi_ss./PMi_ss).^epsY_vec.*alpha_vec.*(CHi_ss + Xi_ss + intermediate_use_guess);

% Now compute proper intermediate use with M_ss
intermediate_use_ss = zeros(nsec,1);
for i = 1:nsec
    for j = 1:nsec
        intermediate_use_ss(i) = intermediate_use_ss(i) + beta_mat(j,i) * (PMi_ss(j)/pH_ss(i))^epsM_vec(j) * M_ss(j);
    end
end

L_ss   = (MCi_ss./PL_ss).^epsY_vec.*(1-alpha_vec-alphaV_vec).*(CHi_ss + Xi_ss + intermediate_use_ss);
Vi_ss  = (MCi_ss./PV_ss).^epsY_vec.*alphaV_vec.*(CHi_ss + Xi_ss + intermediate_use_ss);
Yi_ss  = A_vec.*(alpha_vec.^(1./epsY_vec).*M_ss.^((epsY_vec-1)./epsY_vec)... 
    + alphaV_vec.^(1./epsY_vec).*Vi_ss.^((epsY_vec-1)./epsY_vec)...
    + (1-alphaV_vec-alpha_vec).^(1./epsY_vec).*L_ss.^((epsY_vec-1)./epsY_vec)).^(epsY_vec./(epsY_vec-1));
x_guess1=[M_ss; L_ss; Vi_ss; Yi_ss];
options1 = optimoptions('fsolve','TolFun',1e-10,'Display','off');

x_sol=fsolve(@(x_vec) steady_ntwsoe_system(x_vec,alpha_vec, alphaV_vec, beta_mat, MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss,epsY_vec, epsM_vec, A_vec,pH_ss),x_guess1,options1);




M_ss   = x_sol(1:nsec);
L_ss   = x_sol(nsec+1:2*nsec);
Vi_ss  = x_sol(2*nsec+1:3*nsec);
Yi_ss  = x_sol(3*nsec+1:4*nsec);

V_ss   = sum(Vi_ss);
CF_ss  = sum(CFi_ss);
mkupV  = 1;
IMP_tot_ss = mkupV*(V_ss+CF_ss);

TB_ss  = PX_ss*X_ss-PV_ss*IMP_tot_ss;
GDP_ss = C_ss+TB_ss;
N_ss   = sum(L_ss);  % Total labor = sum of sectoral labor (adjustment costs zero in SS)
Y_ss   = sum(Yi_ss); % Aggregate gross output
r_star_ss = Rworld_ss;
pi_ss  = 1;                % Domestic steady-state inflation (zero)
r_ss   = 1/beta_val;       % Gross interest rate (from Euler equation)

tbgdp  = TB_ss/GDP_ss;

% Compute steady state debt from debt accumulation equation
% In SS: Q*Bstar = -TB + r_star*Q*Bstar/Pistar
% => Bstar = -TB / [Q*(1 - r_star/Pistar)]
Bstar_ss = -TB_ss / (Q_ss*(1 - r_star_ss/Pistar_ss));

% Compute debt-to-GDP target bbar from r_star equation
% r_star = Rworld*exp(-chii_b*(bbar - Q*Bstar/GDP))
% In SS with r_star=Rworld: bbar = Q*Bstar/GDP
bbar_val = Q_ss*Bstar_ss/GDP_ss;

% IMF calibration targets for Chile:
%   Gross output / GDP = 2.0   (IMF, gross-output-based national accounts)
%   Foreign debt / GDP = 0.70  (IMF, external debt statistics)
ygdp_target  = 2.0;   % gross output to GDP ratio (IMF)
fdebt_target = 0.70;  % foreign (external) debt to GDP ratio (IMF)
ygdp_ss      = Y_ss / GDP_ss;        % model-implied gross output / GDP
fdebt_ss     = Q_ss * abs(Bstar_ss) / GDP_ss;  % model-implied |foreign debt| / GDP

Pistar_ss_val=Pistar_ss;
Rworld_ss_val=Rworld_ss;
PVstar_ss_val=PVstar_ss;

%% Create indexed steady state variables for Dynare initval block
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

CF_ss = sum(CFi_ss);
CFg_total_ss = sum(CFg_ss);
CFs_total_ss = sum(CFs_ss);
M_tot_ss     = sum(M_ss);        M_tot_ss_val  = M_tot_ss;
Y_tot_ss     = sum(Yi_ss);
VA_ss        = sum(Yi_ss - M_ss); VA_ss_val     = VA_ss;
Ctotg_ss     = sum(gammag_vec.*(p_g_ss./P_ss).*C_g_ss); Ctotg_ss_val = Ctotg_ss;
Ctots_ss     = sum(gammas_vec.*(p_s_ss./P_ss).*C_s_ss); Ctots_ss_val = Ctots_ss;
Ctot_ss      = Ctotg_ss + Ctots_ss;                     Ctot_ss_val  = Ctot_ss;
IMP_ss_val   = IMP_tot_ss;

% Foreign demand shock parameters (not in gap model, kept for smm_estimates.mat compatibility)
rho_psi_val   = 0.5;     % AR(1) persistence of foreign demand shock
sigma_psi_val = 0.001;   % Std dev of foreign demand shock

% Import price (PVstar) shock parameters
rho_pvstar_val       = 0.9;   % AR(1) persistence of import price shock
sigma_pvstar_val     = 0.03;  % Std dev of import price shock (initial value for SMM estimation)

% Preference/demand shock parameters (Fix 4: intertemporal taste shock)
% xi enters the Euler equation: xi*C^(-gamma) = beta*xi(+1)*C(+1)^(-gamma)*r/pi(+1)
% A positive xi shock raises desired current consumption relative to the future,
% creating demand-driven positive GDP-inflation covariance to offset the
% supply-shock-induced negative correlation observed in the baseline estimation.
rho_xi_val   = 0.80;   % AR(1) persistence of preference shock (initial for SMM)
sigma_xi_val = 0.005;  % Std dev of preference shock (initial for SMM)

% Workspace variables read by NK_SOE_lev_gap2.mod shocks block
% (shock_eps*_val = variance for stoch_simul; 1 = active, 0 = inactive)
shock_eps_om_val     = double(sigma_om_val > 0);   % active when preference shock on
shock_eps_i_val      = double(sigma_i_val  > 0);   % active when monetary shock on
% PVstar shock: activate for EXERCISE=0 when sigma_pvstar was estimated (>0).
% smm_model_moments forces M_loc.Sigma_e(4,4)=1 internally, so the full
% Dynare run must mirror that — otherwise the compiled model misses exchange
% rate dynamics that SMM matched against std(Q) and autocorr(Q).
shock_eps_pvstar_val = double(EXERCISE == 0 && sigma_pvstar_val > 0);
% Preference shock: same convention — activate during baseline compilation.
shock_eps_xi_val     = double(EXERCISE == 0 && sigma_xi_val > 0);
shock_epsA_val       = ones(nsec, 1);  % all TFP shocks active (amplitude controlled by isigma_tfp_val)

fprintf('  GDP=%.4f  TB/GDP=%.4f  X=%.4f  IMP=%.4f  bbar=%.4f\n', ...
    GDP_ss, tbgdp, X_ss, IMP_tot_ss, bbar_val);
fprintf('  TB/GDP target: %.4f  |  world rate (ann.): %.2f%%\n', ...
    tb_target, 400*(Rworld_ss-1));
fprintf('  Bstar=%.4f  Q*Bstar=%.4f  Q*Bstar/GDP=%.3f (%.1f%%)\n', ...
    Bstar_ss, Q_ss*Bstar_ss, Q_ss*Bstar_ss/GDP_ss, 100*Q_ss*Bstar_ss/GDP_ss);
fprintf('  Gross output/GDP : %.3f  (IMF target: %.1f)%s\n', ygdp_ss,  ygdp_target,  ...
    ternary_str(abs(ygdp_ss  - ygdp_target)  > 0.20, '  << off target', ''));
fprintf('  Foreign debt/GDP : %.3f  (IMF target: %.2f)%s\n', fdebt_ss, fdebt_target, ...
    ternary_str(abs(fdebt_ss - fdebt_target) > 0.10, '  << off target', ''));
fprintf('\n');

save params_val_ul.mat

%% Run Dynare
fprintf('--- Running Dynare ---\n');

% Ensure Dynare is available on the MATLAB path. Try common locations and
% the DYNARE_HOME environment variable if the `dynare` function is not found.
if isempty(which('dynare'))
    fprintf('Dynare not found on MATLAB path. Searching common install locations...\n');
    dynare_found = false;

    % Check environment variable first
    dyn_home = getenv('DYNARE_HOME');
    if ~isempty(dyn_home)
        cand = fullfile(dyn_home, 'matlab');
        if isfile(fullfile(cand, 'dynare.m'))
            addpath(cand);
            dynare_found = true;
            found_path = cand;
        end
    end

    % Search common installation roots
    if ~dynare_found
        roots = {'C:\Program Files\Dynare', 'C:\Program Files (x86)\Dynare', 'C:\Dynare'};
        for ri = 1:numel(roots)
            r = roots{ri};
            if isfolder(r)
                d = dir(r);
                for k = 1:numel(d)
                    if d(k).isdir && ~startsWith(d(k).name, '.')
                        cand = fullfile(r, d(k).name, 'matlab');
                        if isfile(fullfile(cand, 'dynare.m'))
                            addpath(cand);
                            dynare_found = true;
                            found_path = cand;
                            break;
                        end
                        % also allow if dynare.m lives directly under the root
                        cand2 = fullfile(r, 'matlab');
                        if isfile(fullfile(cand2, 'dynare.m'))
                            addpath(cand2);
                            dynare_found = true;
                            found_path = cand2;
                            break;
                        end
                    end
                end
            end
            if dynare_found, break; end
        end
    end

    % As a last resort, try the path entries already declared above (pths)
    if ~dynare_found && exist('pths', 'var')
        for i = 1:numel(pths)
            try
                cand = pths{i};
                if isfile(fullfile(cand, 'dynare.m'))
                    addpath(cand);
                    dynare_found = true;
                    found_path = cand;
                    break;
                end
            catch
            end
        end
    end

    rehash;
    if ~isempty(which('dynare'))
        dynare_found = true;
        if ~exist('found_path', 'var')
            found_path = fileparts(which('dynare'));
        end
    end

    if dynare_found
        fprintf('Added Dynare to MATLAB path: %s\n', found_path);
    else
        error(['Dynare function not found. Please install Dynare or add its matlab folder to your MATLAB path.\n' ...
               'Tried environment variable DYNARE_HOME and common locations under C:\\Program Files\\Dynare and C:\\Dynare.\n' ...
               'If Dynare is already installed on this machine, provide its matlab folder path and I can add it to the script.']);
    end
end

% Turn on diary to capture all Dynare output
diary_file = 'dynare_output_log.txt';
diary(diary_file);

try
    evalc('dynare NK_SOE_lev_gap2.mod noclearall');
catch ME
    diary off;
    fprintf('ERROR: Dynare failed – %s\n', ME.message);
    return;
end

diary off;

%% Check Dynare results
simulation_successful = false;

if ~exist('oo_', 'var')
    fprintf('ERROR: Dynare did not complete (oo_ not found).\n');
    return;
end

ss_ok  = isfield(oo_, 'steady_state') && ~isempty(oo_.steady_state) && ...
         ~any(isnan(oo_.steady_state)) && ~any(isinf(oo_.steady_state));
bk_ok  = false;
if isfield(oo_, 'dr') && isfield(oo_.dr, 'eigval') && ~isempty(oo_.dr.eigval)
    bk_ok = sum(abs(oo_.dr.eigval) > 1.000001) == M_.nsfwrd;
end
sim_ok = isfield(oo_, 'endo_simul') && ~isempty(oo_.endo_simul);

fprintf('--- Dynare results ---\n');
if ss_ok,  fprintf('  Steady state  : OK\n');   else, fprintf('  Steady state  : FAILED\n');  end
if bk_ok,  fprintf('  Blanchard-Kahn: OK\n');   else, fprintf('  Blanchard-Kahn: FAILED\n'); end
if sim_ok, fprintf('  Simulation    : OK\n');   else, fprintf('  Simulation    : FAILED\n'); end
fprintf('\n');

if ~sim_ok
    fprintf('Simulation failed – aborting.\n');
    return;
end
simulation_successful = true;

%% Extract simulation results for each endogenous variable (robust)
if exist('oo_', 'var') && isfield(oo_, 'endo_simul') && ~isempty(oo_.endo_simul)
    [nRows, nCols] = size(oo_.endo_simul);
    nEndo = M_.endo_nbr;
    % Determine orientation: rows==nEndo (variables x periods) or cols==nEndo (periods x variables)
    if nRows == nEndo
        nToExtract = nEndo;
        for i = 1:nToExtract
            name = strtrim(char(M_.endo_names(i, :)));
            try
                eval([name ' = oo_.endo_simul(i, :);']);
            catch ME
                warning('Could not assign simulation series to variable "%s": %s', name, ME.message);
            end
        end
    elseif nCols == nEndo
        % Transposed orientation: periods x variables
        nToExtract = nEndo;
        for i = 1:nToExtract
            name = strtrim(char(M_.endo_names(i, :)));
            try
                % take column i and transpose to row vector for consistency
                eval([name ' = oo_.endo_simul(:, i)'';']);
            catch ME
                warning('Could not assign simulation series (transposed) to variable "%s": %s', name, ME.message);
            end
        end
    else
        warning('Simulation results matrix has unexpected shape [%d x %d]; expected one dimension to equal %d (number of endogenous variables). Skipping extraction.', nRows, nCols, nEndo);
    end
else
    warning(['No simulation results found in ''oo_.endo_simul''. This can happen if Dynare failed before producing results ' ...
             'or if you are running the model without calling stoch_simul. Skipping extraction of simulated series.']);
end

%% Initialize variables for storing results
p = {};
Y = {};
A = {};
pl = {};
L = {};
M = {};
MUP = {};
MC = {};
p_m = NaN(nsec, 1);
y_m = NaN(nsec, 1);
A_m = NaN(nsec, 1);
l_m = NaN(nsec, 1);
m_m = NaN(nsec, 1);

Ymat = [];
pmat = [];
name_vec = strings(nsec, 1);

%% Loop through each sector to extract and process results
for i = 1:nsec
    try
        pl{i} = eval(sprintf('PL_%d', i));
        p{i} = eval(sprintf('P_%d', i));
        pm{i} = eval(sprintf('PM_%d', i));
        L{i} = eval(sprintf('L_%d', i));
        Y{i} = eval(sprintf('Y_%d', i));
        M{i} = eval(sprintf('M_%d', i));
        MUP{i} = exp(eval(sprintf('P_%d', i))) ./ exp(eval(sprintf('MC_%d', i)));
        MC{i} = eval(sprintf('MC_%d', i));
        Ymat(:, i) = Y{i};
        pmat(:, i) = p{i};
        A{i} = eval(sprintf('A_%d', i));
        name_vec(i) = names{i}(1:min(50, numel(names{i})));
        p_m(i, 1) = 100 * (mean(p{i}(5)) - mean(p{i}(1)));
        A_m(i, 1) = 100 * (mean(A{i}(5)) - mean(A{i}(1)));
        y_m(i, 1) = 100 * (mean(Y{i}(5)) - mean(Y{i}(1)));
        l_m(i, 1) = 100 * (mean(L{i}(5)) - mean(L{i}(1)));
        m_m(i, 1) = 100 * (mean(M{i}(5)) - mean(M{i}(1)));
    catch ME
        warning('Could not extract sector %d variables: %s', i, ME.message);
    end
end

%% Rank correlations: model vs data (sectoral output, prices, labor)
fprintf('--- Rank correlations (model vs data, Spearman) ---\n');

% Compute analytical percentage std devs from the decision rule (Lyapunov equation).
% This mirrors smm_model_moments.m exactly, and is robust to single-shock exercises
% (EXERCISE=2, etc.) where oo_.endo_simul sample std devs are nearly identical
% across sectors, making Spearman rank-corr undefined from simulated data alone.
std_Y_m  = NaN(nsec, 1);
std_PH_m = NaN(nsec, 1);
std_L_m  = NaN(nsec, 1);

rc_lyap_ok = false;
if isfield(oo_, 'dr') && isfield(oo_.dr, 'ghx') && ~isempty(oo_.dr.ghx)
    try
        T_rc  = oo_.dr.ghx;          % [n_endo_dr × n_state]
        R_rc  = oo_.dr.ghu;          % [n_endo_dr × n_exog]
        Qe_rc = M_.Sigma_e;           % [n_exog × n_exog]
        n_st  = size(T_rc, 2);
        state_rows = zeros(n_st, 1);
        for kk = 1:n_st
            state_rows(kk) = find(oo_.dr.order_var == oo_.dr.state_var(kk), 1);
        end
        A_rc = T_rc(state_rows, :);
        B_rc = R_rc(state_rows, :);
        P_st = local_dlyap(A_rc, B_rc * Qe_rc * B_rc');
        Gamma_rc = T_rc * P_st * T_rc' + R_rc * Qe_rc * R_rc';
        Gamma_rc = (Gamma_rc + Gamma_rc') / 2;

        en_rc  = cellstr(M_.endo_names);
        ov_rc  = oo_.dr.order_var;
        ys_rc  = oo_.dr.ys;          % steady state in declaration order
        for i = 1:nsec
            for kv = 1:3
                switch kv
                    case 1, vn = sprintf('Y_%d',  i);
                    case 2, vn = sprintf('PH_%d', i);
                    case 3, vn = sprintf('L_%d',  i);
                end
                dec_idx = find(strcmp(en_rc, strtrim(vn)), 1);
                if isempty(dec_idx), continue; end
                dr_idx  = find(ov_rc == dec_idx, 1);
                if isempty(dr_idx),  continue; end
                ss_val  = abs(ys_rc(dec_idx));
                if ss_val < 1e-12, ss_val = 1; end
                pstd = sqrt(max(Gamma_rc(dr_idx, dr_idx), 0)) / ss_val;
                switch kv
                    case 1, std_Y_m(i)  = pstd;
                    case 2, std_PH_m(i) = pstd;
                    case 3, std_L_m(i)  = pstd;
                end
            end
        end
        rc_lyap_ok = true;
    catch ME_rc
        fprintf('  [rank corr] Lyapunov approach failed (%s); falling back to simulation std devs.\n', ME_rc.message);
    end
end

% Fallback: simulation-based std devs when Lyapunov is unavailable
if ~rc_lyap_ok && isfield(oo_, 'endo_simul') && ~isempty(oo_.endo_simul)
    en_sim   = cellstr(M_.endo_names);
    n_endo   = numel(en_sim);
    sim_data = oo_.endo_simul;
    if size(sim_data, 1) ~= n_endo && size(sim_data, 2) == n_endo
        sim_data = sim_data';
    end
    ys_sim = oo_.dr.ys;
    for i = 1:nsec
        for kv = 1:3
            switch kv
                case 1, vn = sprintf('Y_%d',  i);
                case 2, vn = sprintf('PH_%d', i);
                case 3, vn = sprintf('L_%d',  i);
            end
            pos    = find(strcmp(en_sim, strtrim(vn)), 1);
            if isempty(pos), continue; end
            ts     = sim_data(pos, :);
            ss_val = abs(ys_sim(pos));
            if ss_val < 1e-12, ss_val = 1; end
            pstd   = std(double(ts)) / ss_val;
            switch kv
                case 1, std_Y_m(i)  = pstd;
                case 2, std_PH_m(i) = pstd;
                case 3, std_L_m(i)  = pstd;
            end
        end
    end
end

% Only compare sectors where both model and data are finite
valid_y = isfinite(std_Y_m)  & isfinite(y_d);
valid_p = isfinite(std_PH_m) & isfinite(p_d);
valid_l = isfinite(std_L_m)  & isfinite(l_d);

% safe_spearman: returns 0 when model vector is constant (undefined rank order),
% matching the behaviour in smm_model_moments.m.
if sum(valid_y) >= 3
    rho_y = safe_spearman_local(std_Y_m(valid_y), y_d(valid_y));
    pval_y = NaN;  % p-value not reliable for safe_spearman; omit
else
    rho_y = NaN; pval_y = NaN;
end
if sum(valid_p) >= 3
    rho_p = safe_spearman_local(std_PH_m(valid_p), p_d(valid_p));
    pval_p = NaN;
else
    rho_p = NaN; pval_p = NaN;
end
if sum(valid_l) >= 3
    rho_l = safe_spearman_local(std_L_m(valid_l), l_d(valid_l));
    pval_l = NaN;
else
    rho_l = NaN; pval_l = NaN;
end

fprintf('  %-30s  %8s\n', 'Variable', 'Spearman r');
fprintf('  %s\n', repmat('-', 1, 42));
fprintf('  %-30s  %8.4f\n', 'Sectoral output (Y)', rho_y);
fprintf('  %-30s  %8.4f\n', 'Prices (PH)',         rho_p);
fprintf('  %-30s  %8.4f\n', 'Labor (L)',           rho_l);
if rho_y == 0 || rho_p == 0 || rho_l == 0
    fprintf('  (0 = model std devs are uniform across sectors for this exercise)\n');
end
fprintf('\n');

rank_corr = struct( ...
    'rho_output', rho_y, 'pval_output', pval_y, ...
    'rho_price',  rho_p, 'pval_price',  pval_p, ...
    'rho_labor',  rho_l, 'pval_labor',  pval_l);

%% Save the model results
output_files = {'model_output_IOSOE', 'model_output_IOSOE_ex1', 'model_output_IOSOE_ex2', 'model_output_IOSOE_ex3'};
save(output_files{EXERCISE+1});

%% Generate plots and perform model checks
if smm_called
    fprintf('\n[SMM] Skipping plots (called from smm_estimation.m).\n');
    set(0, 'DefaultFigureVisible', 'on');  % restore figure visibility
else
    fprintf('--- Generating plots ---\n');
    nsec_val = nsec;
    % Common plots for all exercises (OUTPUT GAPS)
    figs_SOE_gap  % Plots output gaps: deviation from flexible price equilibrium

    % Exercise-specific plots
    if EXERCISE == 1
        plot_shock_effects
    elseif EXERCISE == 2
        plot_manufacturing_shock
        try
            plot_figure7_manufacturing_shock
        catch ME
            fprintf('  WARNING: Figure 7 failed – %s\n', ME.message);
        end
    elseif EXERCISE == 3
        plot_shock_effects
    else
        plot_shock_effects
    end

    steady_state_table
end

fprintf('============================================================\n');
fprintf('  Done: %s\n', exercise_labels{EXERCISE+1});
fprintf('============================================================\n\n');

%% ================================================================ %%
%%  MOMENT FIT TABLE                                                   %%
%%  Loads the definitive moments from smm_results.mat — computed at   %%
%%  theta_hat during the SMM optimisation.  This guarantees the table %%
%%  is identical to the smm_estimation.m output regardless of which   %%
%%  exercise is currently running (exercise-specific Dynare runs have  %%
%%  different active shocks and cannot reproduce the SMM moments).     %%
%% ================================================================ %%
%%  MOMENT FIT TABLE                                                   %%
%%  Computes model moments live from the current Dynare run using the  %%
%%  same Lyapunov-based method as smm_model_moments.m / smm_estimation %%
%% ================================================================ %%
fprintf('--- Moment fit (SMM) ---\n');

% ---- Load data moments (mirrors smm_estimation.m §1 and §3 exactly) ----
% Priority: data_moments_chile.mat (HP-filtered std devs from compute_data_moments.m)
% Fallback: Excel file columns dp/dy/dl (same as smm_estimation.m fallback)
dm_file_fit = fullfile(fileparts(mfilename('fullpath')), 'data_moments_chile.mat');
if exist(dm_file_fit, 'file')
    tmp_dm_fit   = load(dm_file_fit, 'dm_chile');
    dm_fit       = tmp_dm_fit.dm_chile;
    fit_y_d      = dm_fit.y_d;       % [12x1] sectoral output std devs (HP-filtered)
    fit_p_d      = dm_fit.p_d;       % [12x1] sectoral price  std devs
    fit_l_d      = dm_fit.l_d;       % [12x1] sectoral employ std devs
    d_std_GDP    = dm_fit.d_std_GDP;
    d_std_pi     = dm_fit.d_std_pi;
    d_corr_GDPpi = dm_fit.d_corr_GDPpi;
    d_omG        = dm_fit.d_omG;
    d_std_Q      = dm_fit.d_std_Q;
    d_autocorr_Q = dm_fit.d_autocorr_Q;
    d_corr_GDPQ  = dm_fit.d_corr_GDPQ;
    d_TBGDP_fit  = dm_fit.d_TBGDP;
else
    warning('data_moments_chile.mat not found — using Excel file fallback (run compute_data_moments.m for correct values).');
    fit_y_d      = y_d;
    fit_p_d      = p_d;
    fit_l_d      = l_d;
    d_std_GDP    = 0.04107;
    d_std_pi     = 0.02179;
    d_corr_GDPpi = -0.02003;
    d_omG        = ombar_val;
    d_std_Q      = 0.04348;
    d_autocorr_Q = 0.71643;
    d_corr_GDPQ  = 0.13814;
    d_TBGDP_fit  = tb_target;
end

% Stack 46-element data moment vector (same ordering as smm_estimation.m)
data_moments_fit = [
    fit_y_d;       % 1-12   std(Y_i)  — from data_moments_chile.mat
    fit_p_d;       % 13-24  std(PH_i)
    fit_l_d;       % 25-36  std(L_i)
    d_std_GDP;     % 37
    d_std_pi;      % 38
    d_corr_GDPpi;  % 39
    d_omG;         % 40
    d_std_Q;       % 41
    d_autocorr_Q;  % 42
    d_corr_GDPQ;   % 43
    1.0;           % 44  rank corr output
    1.0;           % 45  rank corr prices
    1.0            % 46  rank corr labor
];

% ---- Build baseline struct from current workspace ----
get_ss_fit = @(nm) oo_.dr.ys(find(strcmp(cellstr(M_.endo_names), nm), 1));
Y_ss_fit   = arrayfun(@(i) get_ss_fit(sprintf('Y_%d',i)), 1:nsec)';

baseline_fit = struct();
baseline_fit.nsec         = nsec;
baseline_fit.goods        = goods;
baseline_fit.services     = services;
baseline_fit.Y_ss         = Y_ss_fit;
baseline_fit.GDP_ss       = get_ss_fit('GDP');
baseline_fit.ombar_val    = ombar_val;
baseline_fit.tb_target    = d_TBGDP_fit;
baseline_fit.modalpha     = modalpha;
baseline_fit.modalphaV    = modalphaV;
baseline_fit.modbeta      = modbeta;
baseline_fit.modgammag    = modgammag;
baseline_fit.modgammas    = modgammas;
baseline_fit.modvarrho    = modvarrho;
baseline_fit.modchiX      = modchiX;
baseline_fit.modkappa     = modkappa;
baseline_fit.gamma_val    = gamma_val;
baseline_fit.psi_val      = psi;
baseline_fit.chi_val      = chi;
baseline_fit.epsilon_val  = epsilon;
baseline_fit.beta_val     = beta_val;
baseline_fit.PVstar_ss    = PVstar_ss;
baseline_fit.sigmaH_val   = sigmaH_val;
baseline_fit.etastar_val  = etastar_val;
baseline_fit.omegaX_val   = omegaX_val;
baseline_fit.ystar_ss_val = ystar_ss_val;
baseline_fit.data_std_Y   = fit_y_d;   % HP-filtered sectoral output std devs
baseline_fit.data_std_PH  = fit_p_d;
baseline_fit.data_std_L   = fit_l_d;
if exist('steady_ntwsoe','file')
    baseline_fit.steady_fn = @steady_ntwsoe;
else
    baseline_fit.steady_fn = [];
end
if exist('steady_ntwsoe_system','file')
    baseline_fit.steady_sys_fn = @steady_ntwsoe_system;
else
    baseline_fit.steady_sys_fn = [];
end

% ---- Build theta vector: always load from smm_estimates.mat ----
% We read directly from the file so the moment table always reflects the
% full estimated calibration (all 12 TFP shocks active, baseline EXERCISE=0)
% regardless of which exercise is currently running in the workspace.
smm_est_file_fit = fullfile(fileparts(mfilename('fullpath')), 'smm_estimates.mat');
smm_ckpt_file_fit = fullfile(fileparts(mfilename('fullpath')), 'smm_best_so_far.mat');
if exist(smm_est_file_fit, 'file')
    tmp_est = load(smm_est_file_fit);
    tf_ilabcosts    = tmp_est.ilabcosts_val;
    tf_epsY         = tmp_est.modepsY(1);
    tf_epsM         = tmp_est.modepsM(1);
    tf_kappaV       = log(tmp_est.kappaV_val);
    tf_rho_om       = tmp_est.rho_om1_val;
    tf_sigma_om     = tmp_est.sigma_om_val;
    tf_rho_tfp      = tmp_est.rho_tfp1_val;
    tf_isigma_tfp   = tmp_est.isigma_tfp_val(:);
    tf_rho_pvstar   = tmp_est.rho_pvstar_val;
    tf_sigma_pvstar = tmp_est.sigma_pvstar_val;
    tf_rho_xi       = tmp_est.rho_xi_val;
    tf_sigma_xi     = tmp_est.sigma_xi_val;
    fprintf('  (theta from smm_estimates.mat — full baseline calibration)\n');
elseif exist(smm_ckpt_file_fit, 'file')
    tmp_ckpt_fit = load(smm_ckpt_file_fit, 'smm_best_so_far');
    th_best = tmp_ckpt_fit.smm_best_so_far.theta_best;
    tf_ilabcosts    = th_best(1);
    tf_epsY         = th_best(2);
    tf_epsM         = th_best(3);
    tf_kappaV       = th_best(4);   % already log
    tf_rho_om       = th_best(5);
    tf_sigma_om     = th_best(6);
    tf_rho_tfp      = th_best(7);
    tf_isigma_tfp   = th_best(8:19);
    tf_rho_pvstar   = th_best(20);
    tf_sigma_pvstar = th_best(21);
    tf_rho_xi       = th_best(22);
    tf_sigma_xi     = th_best(23);
    fprintf('  (theta from smm_best_so_far.mat — best checkpoint)\n');
else
    % Fall back to workspace (exercise-specific — rank corr may be 0 for single-shock exercises)
    tf_ilabcosts    = ilabcosts_val;
    tf_epsY         = modepsY(1);
    tf_epsM         = modepsM(1);
    tf_kappaV       = log(kappaV_val);
    tf_rho_om       = rho_om1_val;
    tf_sigma_om     = sigma_om_val;
    tf_rho_tfp      = rho_tfp1_val;
    tf_isigma_tfp   = isigma_tfp_val(:);
    tf_rho_pvstar   = rho_pvstar_val;
    tf_sigma_pvstar = sigma_pvstar_val;
    tf_rho_xi       = rho_xi_val;
    tf_sigma_xi     = sigma_xi_val;
    fprintf('  WARNING: smm_estimates.mat not found — using current exercise parameters.\n');
    fprintf('  Rank correlations may be 0 for single-shock exercises.\n');
end

theta_fit = [
    tf_ilabcosts;       % 1   ilabcosts
    tf_epsY;            % 2   epsY
    tf_epsM;            % 3   epsM
    tf_kappaV;          % 4   log(kappaV)
    tf_rho_om;          % 5   rho_om
    tf_sigma_om;        % 6   sigma_om
    tf_rho_tfp;         % 7   rho_A
    tf_isigma_tfp;      % 8-19 isigma_tfp_i (all 12 sectors)
    tf_rho_pvstar;      % 20  rho_pvstar
    tf_sigma_pvstar;    % 21  sigma_pvstar
    tf_rho_xi;          % 22  rho_xi
    tf_sigma_xi         % 23  sigma_xi
];

% ---- Build weighting matrix W (same as smm_estimation.m) ----
floor_w  = 0.01;
w_diag_fit = 1 ./ max(abs(data_moments_fit), floor_w).^2;
W_fit    = diag(w_diag_fit);
small_idx = abs(data_moments_fit) < 0.02;
W_fit(small_idx, small_idx) = W_fit(small_idx, small_idx) * 0.25;
W_fit(39,39) = W_fit(39,39) * 0.02;
W_fit(40,40) = W_fit(40,40) * 0.10;
W_fit(43,43) = W_fit(43,43) * 0.50;
W_fit(44,44) = W_fit(44,44) * 0.50;
W_fit(45,45) = W_fit(45,45) * 0.50;
W_fit(46,46) = W_fit(46,46) * 0.50;

% ---- Compute model moments from current Dynare run ----
fprintf('  Computing model moments from current Dynare run...\n');
[m_fit, info_fit, ~] = smm_model_moments(theta_fit, M_, options_, oo_, baseline_fit);

if info_fit(1) ~= 0 || any(isnan(m_fit))
    fprintf('  WARNING: smm_model_moments returned info=%d with %d NaN(s). Table may be incomplete.\n', ...
        info_fit(1), sum(isnan(m_fit)));
end

% ---- Compute residuals, weighted loss, and print table ----
resid_fit  = data_moments_fit - m_fit;
wsq_fit    = diag(W_fit) .* resid_fit.^2;
obj_val    = resid_fit' * W_fit * resid_fit;
mnames_fit = smm_moment_names();

sep = repmat('-', 1, 72);
fprintf('\n%-32s  %9s  %9s  %9s  %10s\n', 'Moment', 'Data', 'Model', 'Diff', 'W*Diff^2');
fprintf('%s\n', sep);
for mm = 1:46
    marker = '';
    if isfinite(wsq_fit(mm)) && wsq_fit(mm) > 0.1 * obj_val, marker = ' *'; end
    if isfinite(wsq_fit(mm))
        fprintf('%-32s  %9.5f  %9.5f  %+9.5f  %10.6f%s\n', ...
            mnames_fit{mm}, data_moments_fit(mm), m_fit(mm), resid_fit(mm), wsq_fit(mm), marker);
    else
        fprintf('%-32s  %9.5f  %9.5f  %+9.5f  %10s\n', ...
            mnames_fit{mm}, data_moments_fit(mm), m_fit(mm), resid_fit(mm), 'NaN');
    end
end
fprintf('%s\n', sep);
fprintf('%-32s  %9s  %9s  %9s  %10.6f\n', 'TOTAL LOSS', '', '', '', obj_val);
fprintf('  (* = contributes >10%% of total loss)\n\n');

moment_fit_table = struct('moment_names', {mnames_fit}, ...
    'data', data_moments_fit, 'model', m_fit, 'residuals', resid_fit, ...
    'weighted_sq', wsq_fit, 'total_loss', obj_val);

%% ------------------------------------------------------------------ %%
%%  Local helpers                                                        %%
%% ------------------------------------------------------------------ %%
function s = ternary_str(cond, a, b)
% Return string a when cond is true, b otherwise.
if cond, s = a; else, s = b; end
end

function rho = safe_spearman_local(x, y)
% Spearman rank correlation; returns 0 (not NaN) when x is constant.
% Mirrors safe_spearman() in smm_model_moments.m: a constant model vector
% means sector std devs are indistinguishable (e.g. single-shock exercise).
if std(double(x)) < 1e-12
    rho = 0;
    return
end
try
    rho = corr(x(:), y(:), 'Type', 'Spearman');
    if isnan(rho), rho = 0; end
catch
    rho = 0;
end
end

function P = local_dlyap(A, Q)
% Solve discrete Lyapunov equation  P = A*P*A' + Q  for stable A.
% Tries MATLAB dlyap (Control System Toolbox); falls back to Smith doubling,
% which only requires max|eig(A)| < 1 and no additional toolboxes.
    try
        P = dlyap(A, Q);
        return
    catch
    end
    % Smith doubling: accumulate S = sum_{k>=0} A^k * Q * (A^k)'
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
    P = (P + P') / 2;
end

function names = smm_moment_names()
% Returns the 46-element moment name cell array (mirrors smm_estimation.m).
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


