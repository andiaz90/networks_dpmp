% run_exercises.m
% Runs three shock exercises and saves plots and steady state moments to output folder
% Exercise 1: Shock to preferences for good consumption
% Exercise 2: TFP shock to manufacturing (sector 3)
% Exercise 3: Monetary policy shock

clear all;
close all;
clc;

%% Add paths first (before exercises)
pths = { ...
    'C:\Program Files\Dynare\4.5.6\matlab' ...
    'C:\Program Files\Dynare\6.2\matlab' ...
    'C:\github\networks_dpmp_dme\DPMP-DME-Modelos\Model\utils' ...
    'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github' ...
};
for ip = 1:numel(pths)
    if isfolder(pths{ip})
        addpath(pths{ip});
    end
end

fprintf('\n========================================\n');
fprintf('RUNNING SHOCK EXERCISES\n');
fprintf('========================================\n\n');

% Create output directory structure
output_dir = 'output';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

ex1_dir = fullfile(output_dir, 'exercise1_preference_shock');
ex2_dir = fullfile(output_dir, 'exercise2_tfp_manufacturing');
ex3_dir = fullfile(output_dir, 'exercise3_monetary_policy');

if ~exist(ex1_dir, 'dir'), mkdir(ex1_dir); end
if ~exist(ex2_dir, 'dir'), mkdir(ex2_dir); end
if ~exist(ex3_dir, 'dir'), mkdir(ex3_dir); end

fprintf('Output directories created:\n');
fprintf('  - %s\n', ex1_dir);
fprintf('  - %s\n', ex2_dir);
fprintf('  - %s\n', ex3_dir);

%% ========== EXERCISE 1: PREFERENCE SHOCK ==========
fprintf('\n\n========================================\n');
fprintf('EXERCISE 1: PREFERENCE SHOCK\n');
fprintf('========================================\n\n');

% Run main model with preference shock
main_SOE_exercise1

% Check if Dynare ran successfully
if exist('oo_', 'var') && exist('M_', 'var')
    % Extract and save steady state table
    ss_table_ex1 = extract_steady_state_table('Exercise 1: Preference Shock', oo_, M_);
    writetable(ss_table_ex1, fullfile(ex1_dir, 'steady_state_moments.csv'));
    fprintf('Saved steady state table to: %s\n', fullfile(ex1_dir, 'steady_state_moments.csv'));
    
    % Save plots
    save_exercise_plots(ex1_dir, 'preference_shock');
    
    % Save workspace
    save(fullfile(ex1_dir, 'exercise1_results.mat'));
    fprintf('Saved workspace to: %s\n', fullfile(ex1_dir, 'exercise1_results.mat'));
else
    warning('Exercise 1 did not complete successfully. Skipping output save.');
    fprintf('oo_ exists: %d, M_ exists: %d\n', exist('oo_', 'var'), exist('M_', 'var'));
end

%% ========== EXERCISE 2: TFP SHOCK TO MANUFACTURING ==========
fprintf('\n\n========================================\n');
fprintf('EXERCISE 2: TFP SHOCK TO MANUFACTURING\n');
fprintf('========================================\n\n');

% Clear previous results but preserve directory paths
ex2_dir_temp = ex2_dir;
clearvars -except ex2_dir_temp ex3_dir
close all;
ex2_dir = ex2_dir_temp;
clear ex2_dir_temp;

% Re-add paths after clearing
pths = { ...
    'C:\Program Files\Dynare\4.5.6\matlab' ...
    'C:\Program Files\Dynare\6.2\matlab' ...
    'C:\github\networks_dpmp_dme\DPMP-DME-Modelos\Model\utils' ...
    'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github' ...
};
for ip = 1:numel(pths)
    if isfolder(pths{ip})
        addpath(pths{ip});
    end
end

% Run main model with TFP shock to manufacturing
main_SOE_exercise2

% Check if Dynare ran successfully
if exist('oo_', 'var') && exist('M_', 'var')
    % Extract and save steady state table
    ss_table_ex2 = extract_steady_state_table('Exercise 2: TFP Shock to Manufacturing', oo_, M_);
    writetable(ss_table_ex2, fullfile(ex2_dir, 'steady_state_moments.csv'));
    fprintf('Saved steady state table to: %s\n', fullfile(ex2_dir, 'steady_state_moments.csv'));
    
    % Save plots
    save_exercise_plots(ex2_dir, 'tfp_manufacturing');
    
    % Save workspace
    save(fullfile(ex2_dir, 'exercise2_results.mat'));
    fprintf('Saved workspace to: %s\n', fullfile(ex2_dir, 'exercise2_results.mat'));
else
    warning('Exercise 2 did not complete successfully. Skipping output save.');
    fprintf('oo_ exists: %d, M_ exists: %d\n', exist('oo_', 'var'), exist('M_', 'var'));
end

%% ========== EXERCISE 3: MONETARY POLICY SHOCK ==========
fprintf('\n\n========================================\n');
fprintf('EXERCISE 3: MONETARY POLICY SHOCK\n');
fprintf('========================================\n\n');

% Clear previous results but preserve directory path
ex3_dir_temp = ex3_dir;
clearvars -except ex3_dir_temp
close all;
ex3_dir = ex3_dir_temp;
clear ex3_dir_temp;

% Re-add paths after clearing
pths = { ...
    'C:\Program Files\Dynare\4.5.6\matlab' ...
    'C:\Program Files\Dynare\6.2\matlab' ...
    'C:\github\networks_dpmp_dme\DPMP-DME-Modelos\Model\utils' ...
    'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github' ...
};
for ip = 1:numel(pths)
    if isfolder(pths{ip})
        addpath(pths{ip});
    end
end

% Run main model with monetary policy shock
main_SOE_exercise3

% Check if Dynare ran successfully
if exist('oo_', 'var') && exist('M_', 'var')
    fprintf('Debug: oo_ class = %s, M_ class = %s\n', class(oo_), class(M_));
    fprintf('Debug: oo_ empty = %d, M_ empty = %d\n', isempty(oo_), isempty(M_));
    
    % Extract and save steady state table
    ss_table_ex3 = extract_steady_state_table('Exercise 3: Monetary Policy Shock', oo_, M_);
    writetable(ss_table_ex3, fullfile(ex3_dir, 'steady_state_moments.csv'));
    fprintf('Saved steady state table to: %s\n', fullfile(ex3_dir, 'steady_state_moments.csv'));
    
    % Save plots
    save_exercise_plots(ex3_dir, 'monetary_policy');
    
    % Save workspace
    save(fullfile(ex3_dir, 'exercise3_results.mat'));
    fprintf('Saved workspace to: %s\n', fullfile(ex3_dir, 'exercise3_results.mat'));
else
    warning('Exercise 3 did not complete successfully. Skipping output save.');
    fprintf('oo_ exists: %d, M_ exists: %d\n', exist('oo_', 'var'), exist('M_', 'var'));
end

%% Summary
fprintf('\n\n========================================\n');
fprintf('ALL EXERCISES COMPLETED SUCCESSFULLY\n');
fprintf('========================================\n\n');
fprintf('Results saved in:\n');
fprintf('  1. %s\n', 'output/exercise1_preference_shock');
fprintf('  2. %s\n', 'output/exercise2_tfp_manufacturing');
fprintf('  3. %s\n', 'output/exercise3_monetary_policy');
fprintf('\nEach folder contains:\n');
fprintf('  - steady_state_moments.csv (steady state table)\n');
fprintf('  - Plot figures (*.png or *.fig)\n');
fprintf('  - Full workspace (*.mat)\n\n');
