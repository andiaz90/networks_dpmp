%% ===========================================================================
% PLOTTING SCRIPT: TERMS OF TRADE SHOCK COMPARISON
% ===========================================================================
% Compares IRFs across 3 models for a Terms of Trade shock
% Variables: Output, Consumption, Inflation, Real Interest Rate, 
%            Trade Balance/GDP, Real Exchange Rate
%
% All deviations measured from steady state
% ===========================================================================

fprintf('Generating comparison plots for Terms of Trade shock...\n\n');

%% Load model outputs
models = {'baseline', 'flexible', 'no_io'};
model_names = {'Baseline (Sticky Price + IO)', 'Flexible Price + IO', 'Sticky Price + Diagonal IO'};
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]};
linestyles = {'-', '--', '-.'};

results = struct();
for i = 1:length(models)
    try
        filename = sprintf('model_output_tot_shock_%s.mat', models{i});
        if isfile(filename)
            load(filename, 'oo_', 'M_');
            results.(models{i}).oo_ = oo_;
            results.(models{i}).M_ = M_;
            fprintf('✓ Loaded results for %s\n', models{i});
        else
            fprintf('✗ File not found: %s\n', filename);
        end
    catch ME
        fprintf('✗ Error loading %s: %s\n', models{i}, ME.message);
    end
end

if ~any(isfield(results, models))
    warning('No model results loaded. Cannot generate plots.');
    return;
end

%% Extract variables for comparison
% Initialize storage
nPeriods = 40;  % IRF horizon
variables = struct();

variable_list = {
    'Y',       'log_diff',  'Output (%)';
    'C',       'log_diff',  'Consumption (%)';
    'PI_C',    'level',     'Inflation (annualized %)';
    'RR',      'level',     'Real Interest Rate (annualized %)';
    'TB_Y',    'level',     'Trade Balance / GDP (%)';
    'RX',      'log_diff',  'Real Exchange Rate (%)'
};

fprintf('\nExtracting IRFs...\n');

for j = 1:size(variable_list, 1)
    var_name = variable_list{j, 1};
    var_type = variable_list{j, 2};
    var_display = variable_list{j, 3};
    
    variables.(var_name).display = var_display;
    variables.(var_name).type = var_type;
    variables.(var_name).irf = struct();
    
    for i = 1:length(models)
        if isfield(results, models{i}) && isfield(results.(models{i}), 'oo_')
            oo_ = results.(models{i}).oo_;
            M_ = results.(models{i}).M_;
            
            try
                % Find variable index
                var_idx = strmatch(var_name, M_.endo_names, 'exact');
                
                if ~isempty(var_idx)
                    % Extract IRF to terms of trade shock (assumed to be 'tot_shock')
                    tot_shock_idx = strmatch('tot_shock', M_.exo_names, 'exact');
                    
                    if ~isempty(tot_shock_idx)
                        if isfield(oo_, 'irfs')
                            irf_name = sprintf('%s_tot_shock', var_name);
                            if isfield(oo_.irfs, irf_name)
                                irf_data = oo_.irfs.(irf_name);
                                
                                % Convert to percentage deviations
                                if strcmp(var_type, 'log_diff')
                                    irf_data = 100 * irf_data;
                                elseif strcmp(var_type, 'level')
                                    % Annualize if interest rate or inflation
                                    if strcmpi(var_name, 'PI_C') || strcmpi(var_name, 'RR')
                                        irf_data = 100 * irf_data;  % already in decimal form
                                    end
                                end
                                
                                % Pad to nPeriods
                                if length(irf_data) < nPeriods
                                    irf_data = [irf_data; zeros(nPeriods - length(irf_data), 1)];
                                else
                                    irf_data = irf_data(1:nPeriods);
                                end
                                
                                variables.(var_name).irf.(models{i}) = irf_data;
                                fprintf('  ✓ %s - %s\n', models{i}, var_name);
                            else
                                fprintf('  ✗ IRF not found: %s\n', irf_name);
                            end
                        end
                    else
                        fprintf('  ✗ tot_shock not found in exogenous variables\n');
                    end
                else
                    fprintf('  ✗ Variable %s not found\n', var_name);
                end
            catch ME
                fprintf('  ✗ Error processing %s - %s: %s\n', models{i}, var_name, ME.message);
            end
        end
    end
end

%% Create comparison figure
fprintf('\nCreating comparison figure...\n');

fig = figure('Name', 'Terms of Trade Shock Comparison', 'NumberTitle', 'off', ...
    'Position', [100 100 1400 900], 'Visible', 'on');
drawnow;

% 2x3 subplot layout
subplot_positions = {[1 1], [1 2], [1 3], [2 1], [2 2], [2 3]};

for j = 1:size(variable_list, 1)
    var_name = variable_list{j, 1};
    var_display = variable_list{j, 3};
    
    ax = subplot(2, 3, j);
    hold on;
    grid on;
    
    % Plot each model
    for i = 1:length(models)
        if isfield(variables.(var_name).irf, models{i})
            irf_data = variables.(var_name).irf.(models{i});
            plot(0:nPeriods-1, irf_data, ...
                'Color', colors{i}, ...
                'LineStyle', linestyles{i}, ...
                'LineWidth', 2.5, ...
                'DisplayName', model_names{i});
        end
    end
    
    % Add zero line
    yline(0, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    xlabel('Quarters');
    ylabel(var_display);
    title(var_display, 'FontWeight', 'bold');
    
    % Legend only on first subplot
    if j == 1
        legend('Location', 'best', 'FontSize', 9);
    end
    
    enforce_min_yaxis(gca);
    hold off;
end

% Overall title
sgtitle({'Terms of Trade Shock: Model Comparison',...
         'All deviations from steady state (not from flexible price)'}, ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
print(fig, '-dpng', '-r300', 'terms_of_trade_shock_comparison.png');
fprintf('✓ Figure saved: terms_of_trade_shock_comparison.png\n');

%% Create detailed comparison tables
fprintf('\nGenerating comparison table (impact and peak response)...\n');

for j = 1:size(variable_list, 1)
    var_name = variable_list{j, 1};
    var_display = variable_list{j, 3};
    
    fprintf('\n%s\n', var_display);
    fprintf(repmat('=', 1, 80));
    fprintf('\nModel                                  Impact    Peak      Peak Quarter\n');
    fprintf(repmat('-', 1, 80));
    
    for i = 1:length(models)
        if isfield(variables.(var_name).irf, models{i})
            irf_data = variables.(var_name).irf.(models{i});
            
            impact = irf_data(1);  % First period
            [peak_val, peak_idx] = max(abs(irf_data));
            peak_val = irf_data(peak_idx);
            
            fprintf('%-36s  %8.4f  %8.4f  %6d\n', ...
                model_names{i}, impact, peak_val, peak_idx-1);
        end
    end
end

fprintf('\n========================================\n');
fprintf('✓ Terms of Trade Shock Comparison Complete\n');
fprintf('========================================\n\n');

