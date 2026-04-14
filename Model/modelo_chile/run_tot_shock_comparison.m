%% ===========================================================================
% ALTERNATIVE: TERMS OF TRADE SHOCK COMPARISON (Simplified)
% ===========================================================================
% This script compares the SAME shock across 3 model configurations
% WITHOUT requiring new Dynare .mod file modifications
%
% It re-runs the existing models with:
% 1. Baseline (sticky price + IO)
% 2. Flexible price + IO (use main_SOE_ul parameter set)
% 3. Sticky price + Diagonal IO
%
% Then extracts IRF data and creates comparisons
% ===========================================================================

clear all;
close all;
clc;

fprintf('\n========================================\n');
fprintf('TERMS OF TRADE SHOCK COMPARISON\n');
fprintf('(Alternative: Using existing shock infrastructure)\n');
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
    end
end

nsec = 12;
dynare_path = 'C:\Program Files\dynare\6.2\matlab';
addpath(dynare_path);

%% Data paths (from main_SOE.m)
data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';
filename_industries = fullfile(data_path, 'Stata_to_excel_few_industries_chile.xls');
filename_io         = fullfile(data_path, 'IO_2021_chile.csv');
filename_fpa        = fullfile(data_path, 'fpa_vector_few_industries_chile.csv');

%% READ DATA (identical to main_SOE.m)
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
modbeta_full  = betax';
modbeta_diag = eye(nsec);  % Diagonal: no intermediate inputs

fprintf('Reading price flexibility...\n');
fpa = readmatrix(filename_fpa);
fpa = fpa(1:nsec);
fpa = ones(nsec,1) * mean(fpa);
theta = (1-fpa).^3;

beta_val = 0.986;
kappa = theta*(10-1)./((1-theta).*(1-theta*beta_val));
kappa_flex = ones(nsec,1) * 0.001;  % Flexible: near-zero adjustment costs

% Spending categories
sigmag = spend_good > spend_serv;
sigmas = spend_serv > spend_good;
gammag = spend_good / sum(spend_good);
gammas = spend_serv / sum(spend_serv);

%% Parameters storage
models_config = struct();

% Model 1: Baseline
models_config(1).name = 'Baseline';
models_config(1).description = 'Sticky Price + Full I-O';
models_config(1).beta = modbeta_full;
models_config(1).kappa = kappa;
models_config(1).alpha = alpha;
models_config(1).gamma_g = gammag;
models_config(1).gamma_s = gammas;

% Model 2: Flexible Price
models_config(2).name = 'Flexible Price';
models_config(2).description = 'Flexible Price + Full I-O';
models_config(2).beta = modbeta_full;
models_config(2).kappa = kappa_flex;
models_config(2).alpha = alpha;
models_config(2).gamma_g = gammag;
models_config(2).gamma_s = gammas;

% Model 3: Diagonal IO
models_config(3).name = 'No I-O';
models_config(3).description = 'Sticky Price + Diagonal I-O';
models_config(3).beta = modbeta_diag;
models_config(3).kappa = kappa;
models_config(3).alpha = zeros(nsec,1);  % No intermediate inputs
models_config(3).gamma_g = gammag;
models_config(3).gamma_s = gammas;

%% Run each configuration
fprintf('\n========================================\n');
fprintf('RUNNING MODELS WITH SPECIFIED CONFIGURATIONS\n');
fprintf('========================================\n\n');

results = struct();

for m = 1:length(models_config)
    fprintf('\n>>> MODEL %d: %s\n', m, models_config(m).description);
    fprintf('    Name: %s\n', models_config(m).name);
    fprintf('----------------------------------------\n');
    
    % Store parameters in base workspace for Dynare
    assignin('base', 'modbeta', models_config(m).beta);
    assignin('base', 'modkappa', models_config(m).kappa);
    assignin('base', 'modalpha', models_config(m).alpha);
    assignin('base', 'modalphaV', alpha_V);
    assignin('base', 'modgammag', models_config(m).gamma_g);
    assignin('base', 'modgammas', models_config(m).gamma_s);
    assignin('base', 'nsec', nsec);
    assignin('base', 'names', names);
    
    try
        % Run steady state initialization (mimics main_SOE.m code)
        Yi_ss = ones(nsec,1);
        Li_ss = ones(nsec,1);
        Pi_ss = ones(nsec,1);
        P_ss = 1;
        Y_ss = 1;
        C_ss = 0.65;
        A_ss = 0;
        R_ss = 1/beta_val - 1;
        RX_ss = 1;
        Q_ss = 1;
        Pf_ss = 1;
        
        % Store in base for Dynare
        for i = 1:nsec
            eval(sprintf('Y_ss%d = Yi_ss(%d);', i, i));
            eval(sprintf('L_ss%d = Li_ss(%d);', i, i));
            eval(sprintf('P_ss%d = Pi_ss(%d);', i, i));
        end
        
        % Setup diary
        diary_file = sprintf('dynare_tot_shock_%s.txt', models_config(m).name);
        diary(diary_file);
        
        fprintf('Calling Dynare...\n');
        % Run Dynare
        dynare NK_IOSOE_lev.mod noclearall;
        
        diary off;
        
        % Store results
        results(m).name = models_config(m).name;
        results(m).M_ = M_;
        results(m).oo_ = oo_;
        results(m).success = true;
        
        fprintf('✓ %s model completed successfully\n', models_config(m).name);
        
        % Save workspace
        save(sprintf('model_output_tot_shock_%s.mat', models_config(m).name), ...
            'oo_', 'M_', 'models_config', '-v7.3');
        fprintf('✓ Results saved\n');
        
    catch ME
        diary off;
        fprintf('✗ ERROR in %s: %s\n', models_config(m).name, ME.message);
        results(m).success = false;
        results(m).error = ME.message;
    end
end

%% Now generate plots
if any([results.success])
    fprintf('\n========================================\n');
    fprintf('GENERATING COMPARISON PLOTS\n');
    fprintf('========================================\n\n');
    plot_terms_of_trade_shock_results;
else
    fprintf('\n✗ No successful model runs. Cannot generate plots.\n');
end

fprintf('\n========================================\n');
fprintf('TERMS OF TRADE SHOCK EXERCISE COMPLETED\n');
fprintf('========================================\n\n');

