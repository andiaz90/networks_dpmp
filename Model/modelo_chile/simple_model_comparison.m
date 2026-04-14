%% ===========================================================================
% SIMPLE COMPARISON: Extract and plot results from existing model runs
% ===========================================================================
% This script loads previously computed model results and creates comparison plots
% It works WITH your existing Dynare model structure
% ===========================================================================

clear all;
close all;
clc;

fprintf('\n========================================\n');
fprintf('TERMS OF TRADE SHOCK COMPARISON\n');
fprintf('Loading existing model outputs...\n');
fprintf('========================================\n\n');

%% Try to load baseline model output
fprintf('Loading baseline model results...\n');

baseline_files = {'model_output_IOSOE.mat', 'model_output_IOSOE_ex2.mat'};
baseline_loaded = false;

for i = 1:length(baseline_files)
    if isfile(baseline_files{i})
        load(baseline_files{i}, 'oo_', 'M_');
        results(1).oo_ = oo_;
        results(1).M_ = M_;
        results(1).name = 'Baseline (Sticky Price + I-O)';
        results(1).success = true;
        fprintf('✓ Loaded: %s\n', baseline_files{i});
        baseline_loaded = true;
        break;
    end
end

if ~baseline_loaded
    fprintf('✗ No baseline model output found. Run main_SOE first.\n');
    return;
end

%% Try to load flexible price model output (from main_SOE_ul)
fprintf('\nLoading flexible price model results...\n');

flex_files = {'model_output_IOSOE_ul.mat'};
flex_loaded = false;

for i = 1:length(flex_files)
    if isfile(flex_files{i})
        load(flex_files{i}, 'oo_', 'M_');
        results(2).oo_ = oo_;
        results(2).M_ = M_;
        results(2).name = 'Flexible Price + I-O';
        results(2).success = true;
        fprintf('✓ Loaded: %s\n', flex_files{i});
        flex_loaded = true;
        break;
    end
end

if ~flex_loaded
    fprintf('⚠ Flexible price model output not found (optional)\n');
    results(2).success = false;
end

%% Configuration for available models
n_models_available = sum([results.success]);
fprintf('\n========================================\n');
fprintf('Models loaded: %d\n', n_models_available);
fprintf('========================================\n\n');

if n_models_available < 1
    fprintf('✗ No model results available. Cannot generate plots.\n');
    return;
end

%% Extract available shocks and variables
fprintf('Analyzing available data in models...\n\n');

% Get info from first successful model
for m = 1:length(results)
    if results(m).success
        M_ = results(m).M_;
        oo_ = results(m).oo_;
        
        fprintf('Model: %s\n', results(m).name);
        fprintf('  Endogenous variables: %d\n', M_.endo_nbr);
        fprintf('  Exogenous shocks: %d\n', M_.exo_nbr);
        
        % List available shocks
        fprintf('  Available shocks:\n');
        for s = 1:min(10, M_.exo_nbr)
            fprintf('    - %s\n', M_.exo_names{s});
        end
        if M_.exo_nbr > 10
            fprintf('    ... and %d more\n', M_.exo_nbr - 10);
        end
        
        % Check for IRFs
        if isfield(oo_, 'irfs')
            irf_names = fieldnames(oo_.irfs);
            fprintf('  IRFs available: %d\n', length(irf_names));
            if ~isempty(irf_names)
                fprintf('  First few IRFs:\n');
                for i = 1:min(5, length(irf_names))
                    fprintf('    - %s\n', irf_names{i});
                end
            end
        else
            fprintf('  No IRFs available in oo_ structure\n');
        end
        fprintf('\n');
        break;
    end
end

%% Determine which shock to use
fprintf('Determining which shock to analyze...\n');

% Look for common shock names
shock_candidates = {'eps_a_2', 'eps_a_', 'e_a', 'a_shock', 'eps_pref', 'eps_r'};
selected_shock = '';

for m = 1:length(results)
    if results(m).success
        M_ = results(m).M_;
        
        % Try to find a shock that has IRFs
        if isfield(results(m).oo_, 'irfs')
            for s = 1:length(shock_candidates)
                irf_list = fieldnames(results(m).oo_.irfs);
                matching = contains(irf_list, shock_candidates{s});
                if any(matching)
                    selected_shock = shock_candidates{s};
                    fprintf('✓ Selected shock: %s\n', selected_shock);
                    break;
                end
            end
        end
        if ~isempty(selected_shock)
            break;
        end
    end
end

if isempty(selected_shock)
    fprintf('✗ Could not identify a shock with available IRFs.\n');
    fprintf('Please run main_SOE with stoch_simul to generate IRFs.\n');
    return;
end

%% Extract variable IRFs
fprintf('\nExtracting IRF data...\n\n');

nPeriods = 40;

% Variables to extract
vars_to_extract = {
    'Y', 'Output';
    'C', 'Consumption';
    'PI_C', 'Inflation';
    'RR', 'Real Rate';
    'TB_Y', 'Trade Balance/GDP';
    'RX', 'Real Exchange Rate'
};

irfs_all = struct();

for m = 1:length(results)
    if ~results(m).success
        continue;
    end
    
    oo_ = results(m).oo_;
    M_ = results(m).M_;
    
    fprintf('Extracting from: %s\n', results(m).name);
    
    if ~isfield(oo_, 'irfs')
        fprintf('  ✗ No IRFs in this model\n');
        results(m).success = false;
        continue;
    end
    
    for v = 1:size(vars_to_extract, 1)
        var_name = vars_to_extract{v, 1};
        var_display = vars_to_extract{v, 2};
        
        % Try to find IRF for this variable and shock
        irf_field = sprintf('%s_%s', var_name, selected_shock);
        
        if isfield(oo_.irfs, irf_field)
            irf_data = oo_.irfs.(irf_field);
            
            % Convert to percentage
            irf_data = 100 * irf_data;
            
            % Pad/trim to nPeriods
            if length(irf_data) < nPeriods
                irf_data = [irf_data; zeros(nPeriods - length(irf_data), 1)];
            else
                irf_data = irf_data(1:nPeriods);
            end
            
            irfs_all.(sprintf('model%d_%s', m, var_name)) = irf_data;
            fprintf('  ✓ %s\n', var_name);
        else
            fprintf('  ✗ %s_%s not found\n', var_name, selected_shock);
        end
    end
    fprintf('\n');
end

%% Create comparison figure
fprintf('Creating comparison figure...\n\n');

fig = figure('Name', 'Model Comparison', ...
    'NumberTitle', 'off', ...
    'Position', [50 50 1400 900], ...
    'Visible', 'on');
drawnow;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]};
linestyles = {'-', '--', '-.'};

% Plot 2x3 layout
for v = 1:size(vars_to_extract, 1)
    var_name = vars_to_extract{v, 1};
    var_display = vars_to_extract{v, 2};
    
    subplot(2, 3, v);
    hold on;
    grid on;
    grid minor;
    
    has_data_plotted = false;
    
    % Plot each model
    for m = 1:length(results)
        if results(m).success
            irf_key = sprintf('model%d_%s', m, var_name);
            
            if isfield(irfs_all, irf_key)
                irf_data = irfs_all.(irf_key);
                plot(0:nPeriods-1, irf_data, ...
                    'Color', colors{m}, ...
                    'LineStyle', linestyles{m}, ...
                    'LineWidth', 2.5, ...
                    'DisplayName', results(m).name);
                has_data_plotted = true;
            end
        end
    end
    
    % Add zero line
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    
    xlabel('Quarters', 'FontSize', 10);
    ylabel(sprintf('%s (% dev from SS)', var_display), 'FontSize', 10);
    title(var_display, 'FontSize', 11, 'FontWeight', 'bold');
    
    if v == 6
        lg = legend('Location', 'best', 'FontSize', 9);
        lg.Box = 'on';
    end
    
    xlim([0, nPeriods-1]);
    hold off;
end

% Title
sgtitle(sprintf('Model Comparison: %s shock (40 quarters)', selected_shock), ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save
try
    print(fig, '-dpng', '-r300', 'model_comparison.png');
    fprintf('✓ Figure saved: model_comparison.png\n\n');
catch
    fprintf('✗ Could not save figure\n\n');
end

%% Display summary
fprintf('========================================\n');
fprintf('Summary Statistics\n');
fprintf('========================================\n\n');

for v = 1:size(vars_to_extract, 1)
    var_name = vars_to_extract{v, 1};
    var_display = vars_to_extract{v, 2};
    
    fprintf('%s\n', var_display);
    fprintf(repmat('-', 1, 80));
    fprintf('\n%-30s %12s %12s %12s\n', 'Model', 'Impact', 'Peak', 'Peak Qtr');
    fprintf(repmat('-', 1, 80));
    fprintf('\n');
    
    for m = 1:length(results)
        if results(m).success
            irf_key = sprintf('model%d_%s', m, var_name);
            
            if isfield(irfs_all, irf_key)
                irf_data = irfs_all.(irf_key);
                
                impact = irf_data(1);
                [peak_abs, peak_idx] = max(abs(irf_data));
                peak_val = irf_data(peak_idx);
                
                fprintf('%-30s %12.4f %12.4f %12d\n', ...
                    results(m).name, impact, peak_val, peak_idx-1);
            end
        end
    end
    fprintf('\n');
end

fprintf('========================================\n');
fprintf('✓ Comparison Complete\n');
fprintf('========================================\n\n');

