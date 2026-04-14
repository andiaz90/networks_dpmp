%% ===========================================================================
% TERMS OF TRADE SHOCK COMPARISON EXERCISE
% ===========================================================================
% This script compares IRFs to a terms of trade shock across 3 models:
% 1. Baseline (sticky price with IO structure) - main_SOE
% 2. Flexible Price (main_SOE_ul)
% 3. No IO (diagonal IO matrix - only direct production)
%
% Variables plotted: Output, Consumption, Inflation, Real Interest Rate, 
%                    Trade Balance/GDP, Real Exchange Rate
%
% All deviations measured from steady state (not from flexible price)
% ===========================================================================

clear all;
close all;
clc;

fprintf('\n========================================\n');
fprintf('TERMS OF TRADE SHOCK COMPARISON\n');
fprintf('========================================\n\n');

restoredefaultpath;
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
dynare_path = 'C:\Program Files\dynare\6.2\matlab';
addpath(dynare_path);

%% Data paths
data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';
filename_industries = fullfile(data_path, 'Stata_to_excel_few_industries_chile.xls');
filename_io         = fullfile(data_path, 'IO_2021_chile.csv');
filename_fpa        = fullfile(data_path, 'fpa_vector_few_industries_chile.csv');

%% READ DATA (12 sectors)
varNames = {'kk','BEAName','Nameshort','dp','dy','dl','dri','share','industrytype','spend_good','spend_serv','alpha','alpha_V', 'var_rho'};
varTypes = {'double','char','char','double','double','double','double','double','categorical','double','double','double','double','double'};
opts = spreadsheetImportOptions('NumVariables',14,'VariableNames',varNames,'VariableTypes',varTypes,'DataRange','A2');

fprintf('Reading industry data...\n');
alldata    = readtable(filename_industries, opts);
names      = table2array(alldata(:,3));
p_d        = table2array(alldata(:,4));
y_d        = table2array(alldata(:,5));
l_d        = table2array(alldata(:,6));
ri_d       = table2array(alldata(:,7));
spend_good = table2array(alldata(:,10));
spend_serv = table2array(alldata(:,11));
alpha      = table2array(alldata(:,12));
alpha_V    = table2array(alldata(:,13));
var_rho    = table2array(alldata(:,14));

fprintf('Reading I-O matrix...\n');
betaio = readmatrix(filename_io);
betaio = betaio(1:nsec,1:nsec);
betax  = betaio ./ sum(betaio,1);
modbeta  = betax';
modalpha = alpha;
modalphaV = alpha_V;

fprintf('Reading price flexibility...\n');
fpa = readmatrix(filename_fpa);
fpa = fpa(1:nsec);
fpa = ones(nsec,1) * mean(fpa);
theta = (1-fpa).^3;

beta_val = 0.986;
kappa = theta*(10-1)./((1-theta).*(1-theta*beta_val));
modkappa = kappa;

% Spending categories
goods    = spend_good > spend_serv;
services = spend_serv > spend_good;
gammag   = spend_good / sum(spend_good);
modgammag = gammag;
gammas   = spend_serv / sum(spend_serv);
modgammas = gammas;

modcl     = zeros(nsec,1);
modclneg  = zeros(nsec,1);
modcm     = zeros(nsec,1);

fprintf('\n========================================\n');
fprintf('MODEL 1: BASELINE (Sticky Price + IO)\n');
fprintf('========================================\n');
run_terms_of_trade_shock_model('baseline', modbeta, modalpha, modalphaV, modkappa, ...
    modgammag, modgammas, modcl, modclneg, modcm, nsec);

fprintf('\n========================================\n');
fprintf('MODEL 2: FLEXIBLE PRICE + IO\n');
fprintf('========================================\n');
% For flexible price, set kappa to very small value
modkappa_flex = zeros(nsec,1) + 0.001;
run_tot_shock_comparison;

fprintf('\n========================================\n');
fprintf('MODEL 3: BASELINE (Sticky Price) + DIAGONAL IO\n');
fprintf('========================================\n');
% Create diagonal IO matrix (no intermediate inputs)
modbeta_diag = eye(nsec);
modalpha_diag = zeros(nsec,1);  % No intermediate inputs
modalphaV_diag = modalphaV;      % Keep import structure

fprintf('\n========================================\n');
fprintf('GENERATING COMPARISON PLOTS\n');
fprintf('========================================\n');

% Generate comparison plots
plot_terms_of_trade_shock_results;

fprintf('\n========================================\n');
fprintf('TERMS OF TRADE SHOCK COMPARISON COMPLETED\n');
fprintf('========================================\n\n');

