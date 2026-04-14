%% ===========================================================================
% PLOT RESULTS: TERMS OF TRADE SHOCK COMPARISON
% ===========================================================================
% Creates comparison plots of IRFs across 3 models
% Compares: Baseline, Flexible Price, and No-I-O configurations
% ===========================================================================

fprintf('Generating Terms of Trade shock comparison plots...\n');

%% Setup
model_names = {'Baseline', 'Flexible Price', 'No I-O'};
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]};
linestyles = {'-', '--', '-.'};

nPeriods = 40;
shock_var = 'tot_shock';  % Terms of trade shock

% Variables to plot
vars_to_plot = {
    'Y',       'log_dev',  'Output';
    'C',       'log_dev',  'Consumption';
    'PI_C',    'pct',      'CPI Inflation (annualized)';
    'RR',      'pct',      'Real Interest Rate (annualized)';
    'TB_Y',    'level',    'Trade Balance / GDP';
    'RX',      'log_dev',  'Real Exchange Rate'
};

%% Extract IRFs from loaded models
n_models = length(model_names);
irfs = struct();

fprintf('\nExtracting IRFs from Dynare results...\n');

for m = 1:n_models
    if ~isfield(results, 'success') || m > length(results)
        fprintf('  ✗ Model %d results not available\n', m);
        continue;
    end
    
    if ~results(m).success
        fprintf('  ✗ Model %d (%s) failed\n', m, model_names{m});
        continue;
    end
    
    oo_ = results(m).oo_;
    M_ = results(m).M_;
    
    fprintf('  Processing %s...\n', model_names{m});
    
    % Get shock index
    tot_shock_idx = find(strcmp(M_.exo_names, shock_var), 1);
    
    if isempty(tot_shock_idx)
        fprintf('    ✗ Shock variable not found: %s\n', shock_var);
        continue;
    end
    
    % Extract IRFs for each variable
    for v = 1:size(vars_to_plot, 1)
        var_name = vars_to_plot{v, 1};
        var_type = vars_to_plot{v, 2};
        
        var_idx = find(strcmp(M_.endo_names, var_name), 1);
        
        if isempty(var_idx)
            fprintf('    ✗ Variable %s not found\n', var_name);
            continue;
        end
        
        % IRF naming convention: VAR_SHOCK
        irf_field = sprintf('%s_%s', var_name, shock_var);
        
        if isfield(oo_.irfs, irf_field)
            irf_data = oo_.irfs.(irf_field);
            
            % Convert to appropriate units
            if strcmp(var_type, 'log_dev')
                % Log deviation: multiply by 100 for percentage
                irf_data = 100 * irf_data;
            elseif strcmp(var_type, 'pct')
                % Already in percentage form, may need to annualize
                if contains(var_name, {'PI', 'RR'})  % Inflation or interest rate
                    irf_data = 400 * irf_data;  % Quarterly to annualized
                else
                    irf_data = 100 * irf_data;
                end
            elseif strcmp(var_type, 'level')
                % Percentage level (TB/Y already in %)
                irf_data = 100 * irf_data;
            end
            
            % Pad to nPeriods if needed
            if length(irf_data) < nPeriods
                irf_data = [irf_data; zeros(nPeriods - length(irf_data), 1)];
            else
                irf_data = irf_data(1:nPeriods);
            end
            
            irfs.(sprintf('m%d_%s', m, var_name)) = irf_data;
            
        else
            fprintf('    ✗ IRF field not found: %s\n', irf_field);
        end
    end
    
    fprintf('    ✓ Complete\n');
end

%% Create comparison figure
fprintf('\nCreating comparison figure...\n');

fig = figure('Name', 'Terms of Trade Shock Comparison', ...
    'NumberTitle', 'off', ...
    'Position', [50 50 1500 950], ...
    'Visible', 'on');
drawnow;

% 2x3 layout
for v = 1:size(vars_to_plot, 1)
    var_name = vars_to_plot{v, 1};
    var_display = vars_to_plot{v, 3};
    
    subplot(2, 3, v);
    hold on;
    grid on;
    grid minor;
    
    has_data = false;
    
    % Plot each model
    for m = 1:n_models
        irf_field = sprintf('m%d_%s', m, var_name);
        
        if isfield(irfs, irf_field)
            irf_data = irfs.(irf_field);
            plot(0:nPeriods-1, irf_data, ...
                'Color', colors{m}, ...
                'LineStyle', linestyles{m}, ...
                'LineWidth', 2.5, ...
                'DisplayName', model_names{m});
            has_data = true;
        end
    end
    
    % Add zero line
    yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    
    xlabel('Quarters', 'FontSize', 10);
    ylabel(var_display, 'FontSize', 10);
    title(var_display, 'FontSize', 11, 'FontWeight', 'bold');
    
    % Legend on bottom-right subplot
    if v == 6 && has_data
        leg = legend('Location', 'best', 'FontSize', 9);
        leg.Box = 'on';
    end
    
    xlim([0, nPeriods-1]);
    enforce_min_yaxis(gca);
    hold off;
end

% Master title
sgtitle({'Terms of Trade Shock: Model Comparison (40 Quarters)',...
         'Baseline (Sticky+IO) vs Flexible Price vs No I-O', ...
         'All deviations measured from steady state'}, ...
    'FontSize', 13, 'FontWeight', 'bold');

% Save
try
    print(fig, '-dpng', '-r300', 'tot_shock_comparison.png');
    fprintf('✓ Figure saved: tot_shock_comparison.png\n');
catch
    fprintf('✗ Could not save figure\n');
end

%% Create summary table
fprintf('\n========================================\n');
fprintf('IMPULSE RESPONSE SUMMARY TABLE\n');
fprintf('========================================\n');

for v = 1:size(vars_to_plot, 1)
    var_name = vars_to_plot{v, 1};
    var_display = vars_to_plot{v, 3};
    
    fprintf('\n%s\n', var_display);
    fprintf(repmat('=', 1, 95));
    fprintf('\n%-20s %12s %12s %12s %12s %12s\n', ...
        'Model', 'Impact', 'Peak', 'Peak Qtr', 'Cumulative', 'Final');
    fprintf(repmat('-', 1, 95));
    
    for m = 1:n_models
        irf_field = sprintf('m%d_%s', m, var_name);
        
        if isfield(irfs, irf_field)
            irf_data = irfs.(irf_field);
            
            impact = irf_data(1);
            [peak_abs, peak_idx] = max(abs(irf_data));
            peak_val = irf_data(peak_idx);
            cumsum_val = sum(irf_data);
            final_val = irf_data(end);
            
            fprintf('%-20s %12.4f %12.4f %12d %12.4f %12.4f\n', ...
                model_names{m}, impact, peak_val, peak_idx-1, cumsum_val, final_val);
        end
    end
end

fprintf('\n========================================\n');
fprintf('✓ COMPARISON ANALYSIS COMPLETE\n');
fprintf('========================================\n\n');

fprintf('Notes:\n');
fprintf('- Impact = first period (quarter 0) response\n');
fprintf('- Peak = maximum absolute value during 40-quarter horizon\n');
fprintf('- Peak Qtr = quarter when peak is reached\n');
fprintf('- Cumulative = sum of all 40 periods\n');
fprintf('- Final = response at quarter 39\n\n');

