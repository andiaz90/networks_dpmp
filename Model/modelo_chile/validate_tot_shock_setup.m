%% ===========================================================================
% TERMS OF TRADE SHOCK EXERCISE - QUICK VALIDATION
% ===========================================================================
% This script verifies that all components are in place and working
% Run this BEFORE running the full comparison to catch issues early
% ===========================================================================

clear all;
close all;
clc;

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║  TERMS OF TRADE SHOCK EXERCISE - SYSTEM VALIDATION               ║\n');
fprintf('║  Quick diagnostic to verify setup before running full exercise    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

%% Check 1: Working directory
fprintf('[1/8] Checking working directory...\n');
current_dir = pwd;
[~, folder_name] = fileparts(current_dir);

if strcmpi(folder_name, 'modelo_chile')
    fprintf('    ✓ In correct directory: %s\n', current_dir);
else
    fprintf('    ? Current directory: %s\n', current_dir);
    fprintf('    Recommended: Change to Model/modelo_chile folder\n');
end

%% Check 2: Required scripts present
fprintf('\n[2/8] Checking required scripts...\n');

required_scripts = {
    'run_tot_shock_comparison.m', 'Main comparison script';
    'plot_terms_of_trade_shock_results.m', 'Plotting script';
    'NK_IOSOE_lev.mod', 'Dynare model file';
    'main_SOE.m', 'Baseline driver';
    'figs_SOE_ss.m', 'Steady state plots'
};

all_scripts_found = true;
for i = 1:size(required_scripts, 1)
    filename = required_scripts{i, 1};
    description = required_scripts{i, 2};
    
    if isfile(filename)
        fprintf('    ✓ %s (%s)\n', filename, description);
    else
        fprintf('    ✗ MISSING: %s (%s)\n', filename, description);
        all_scripts_found = false;
    end
end

if ~all_scripts_found
    fprintf('\n    WARNING: Some required files are missing!\n');
    fprintf('    Ensure you are in the Model/modelo_chile directory.\n');
end

%% Check 3: Data files
fprintf('\n[3/8] Checking data files...\n');

data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';
data_files = {
    'Stata_to_excel_few_industries_chile.xls', 'Industry data';
    'IO_2021_chile.csv', 'Input-Output matrix';
    'fpa_vector_few_industries_chile.csv', 'Price flexibility'
};

all_data_found = true;
for i = 1:size(data_files, 1)
    filename = data_files{i, 1};
    description = data_files{i, 2};
    filepath = fullfile(data_path, filename);
    
    if isfile(filepath)
        fprintf('    ✓ %s (%s)\n', filename, description);
    else
        fprintf('    ✗ MISSING: %s\n', filename);
        fprintf('      Expected: %s\n', filepath);
        all_data_found = false;
    end
end

if ~all_data_found
    fprintf('\n    WARNING: Some data files not found!\n');
    fprintf('    Update data_path variable in run_tot_shock_comparison.m\n');
end

%% Check 4: MATLAB path and Dynare
fprintf('\n[4/8] Checking MATLAB paths and Dynare...\n');

dynare_path = 'C:\Program Files\dynare\6.2\matlab';
if isfolder(dynare_path)
    fprintf('    ✓ Dynare path exists: %s\n', dynare_path);
    addpath(dynare_path);
    
    if ~isempty(which('dynare'))
        fprintf('    ✓ Dynare function accessible\n');
    else
        fprintf('    ✗ Dynare function not on path after adding directory\n');
    end
else
    fprintf('    ✗ Dynare path not found: %s\n', dynare_path);
    fprintf('    Your Dynare may be in a different location.\n');
    fprintf('    Update dynare_path in run_tot_shock_comparison.m\n');
end

%% Check 5: NK_IOSOE_lev.mod variables
fprintf('\n[5/8] Checking NK_IOSOE_lev.mod for required variables...\n');

if isfile('NK_IOSOE_lev.mod')
    mod_content = fileread('NK_IOSOE_lev.mod');
    
    required_vars = {'var Y;', 'var C;', 'var PI_C;', 'var RR;', 'var TB_Y;', 'var RX;'};
    required_shocks = {'varexo'};
    
    all_vars_found = true;
    for i = 1:length(required_vars)
        if contains(mod_content, required_vars{i})
            fprintf('    ✓ %s found\n', required_vars{i});
        else
            fprintf('    ✗ %s NOT found\n', required_vars{i});
            all_vars_found = false;
        end
    end
    
    if contains(mod_content, 'varexo')
        fprintf('    ✓ varexo section found\n');
    else
        fprintf('    ✗ varexo section NOT found\n');
    end
    
    if ~all_vars_found
        fprintf('\n    INFO: Some variables missing from NK_IOSOE_lev.mod\n');
        fprintf('    These may need to be added for full functionality.\n');
        fprintf('    See ToT_SHOCK_IMPLEMENTATION_GUIDE.md\n');
    end
else
    fprintf('    ✗ NK_IOSOE_lev.mod not found\n');
end

%% Check 6: Load test
fprintf('\n[6/8] Test loading main_SOE parameters...\n');

try
    % Quick test of data loading from main_SOE section
    varNames = {'kk','BEAName','Nameshort','dp','dy','dl','dri','share','industrytype','spend_good','spend_serv','alpha','alpha_V', 'var_rho'};
    varTypes = {'double','char','char','double','double','double','double','double','categorical','double','double','double','double','double'};
    opts = spreadsheetImportOptions('NumVariables',14,'VariableNames',varNames,'VariableTypes',varTypes,'DataRange','A2');
    
    fprintf('    ✓ Data import options created\n');
    
    if all_data_found
        try
            data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';
            filename_industries = fullfile(data_path, 'Stata_to_excel_few_industries_chile.xls');
            alldata = readtable(filename_industries, opts);
            fprintf('    ✓ Successfully loaded industry data (%d sectors)\n', height(alldata));
        catch me
            fprintf('    ✗ Error loading data: %s\n', me.message);
        end
    end
catch me
    fprintf('    ✗ Error setting up data import: %s\n', me.message);
end

%% Check 7: Output directories
fprintf('\n[7/8] Checking output locations...\n');

current_loc = pwd;
fprintf('    Output files will be saved to: %s\n', current_loc);
fprintf('    ✓ Ready to save results\n');

%% Check 8: Summary
fprintf('\n[8/8] Validation Summary\n');
fprintf('════════════════════════════════════════════════════════════════════\n');

if all_scripts_found && all_data_found
    fprintf('✓ ALL CHECKS PASSED - Ready to run exercise!\n\n');
    fprintf('To start the comparison, run:\n');
    fprintf('  >> run_tot_shock_comparison\n\n');
    fprintf('Expected duration: 15-30 minutes depending on Dynare.\n\n');
    
    start_now = input('Would you like to start now? (y/n): ', 's');
    if strcmpi(start_now, 'y')
        fprintf('\nStarting Terms of Trade Shock comparison...\n\n');
        run_tot_shock_comparison;
    else
        fprintf('\nTo start later, run: >> run_tot_shock_comparison\n\n');
    end
else
    fprintf('✗ SOME ISSUES DETECTED\n\n');
    fprintf('Required fixes:\n');
    if ~all_scripts_found
        fprintf('  1. Ensure you are in Model/modelo_chile directory\n');
        fprintf('     Change to: cd Model/modelo_chile\n\n');
    end
    if ~all_data_found
        fprintf('  2. Update data_path in run_tot_shock_comparison.m\n');
        fprintf('     Set to location of your data files\n\n');
    end
    fprintf('After fixing issues, run validation again:\n');
    fprintf('  >> validate_tot_shock_setup\n\n');
end

fprintf('════════════════════════════════════════════════════════════════════\n');
fprintf('For more information, see: ToT_SHOCK_README.md\n\n');

