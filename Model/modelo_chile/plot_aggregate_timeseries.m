
% Set time horizon - use existing settings if available, otherwise default
if exist('horizon_scatter_model', 'var') && exist('t_end', 'var')
    time_end = t_end;
    time = 0:time_end-1;
else
    time_end = min(40, size(oo_.endo_simul, 2)); % Default to 40 periods or available data
    time = 0:time_end-1;
end

% Define key variables to extract and plot
key_vars = {
    'C',        'Total Consumption';
    'C_g',      'Consumption of Goods';
    'C_s',      'Consumption of Services';
    'Y',        'Output';
    'pi',       'Inflation';
    'r',        'Interest Rate';
    'N',        'Labor';
    'w',        'Wage';
    'p_g',      'Goods Price';
    'p_s',      'Services Price';
    'TB',       'Trade Balance';
    'Q',        'Exchange Rate';
    'X',        'Exports';
    'Vimp',     'Imports'
};

% Extract time series data
data = struct();
available_vars = {};
available_labels = {};

for i = 1:size(key_vars, 1)
    var_name = key_vars{i, 1};
    var_label = key_vars{i, 2};
    
    % Find variable index in Dynare results
    var_idx = strmatch(var_name, M_.endo_names, 'exact');
    
    if ~isempty(var_idx)
        % Extract time series
        data.(var_name) = oo_.endo_simul(var_idx, 1:time_end);
        available_vars{end+1} = var_name;
        available_labels{end+1} = var_label;
    else
        fprintf('Warning: Variable %s not found in simulation results.\n', var_name);
    end
end

if isempty(available_vars)
    error('No variables found in simulation results. Check model solution.');
end

% Create figure
fig = figure('Name', 'Aggregate Macroeconomic Time Series', ...
             'Position', [100 50 1600 1000], ...
             'Visible', 'on', ...
             'WindowStyle', 'docked');
clf(fig);  % Clear figure to prevent multiple overlapping plots

% Calculate subplot layout
n_vars = length(available_vars);
n_cols = 4;  % 4 columns like plot_12_figs
n_rows = ceil(n_vars / n_cols);

% Plot each variable
for i = 1:n_vars
    subplot(n_rows, n_cols, i);
    cla reset;  % Clear axis to prevent multiple overlapping plots
    
    var_name = available_vars{i};
    var_label = available_labels{i};
    var_data = data.(var_name);
    
    % Convert to appropriate units and calculate deviations
    if strcmp(var_name, 'pi')
        % Inflation: get steady state and convert to deviation in annualized percentage points
        pi_idx = strmatch('pi', M_.endo_names, 'exact');
        pi_ss = oo_.steady_state(pi_idx);
        plot_data = 400 * (var_data - pi_ss);  % Annualized deviation in percentage points
        ylabel_text = 'pp dev. (ann.)';
    elseif strcmp(var_name, 'r')
        % Interest rate: get steady state and convert to deviation in annualized percentage points
        r_idx = strmatch('r', M_.endo_names, 'exact');
        r_ss = oo_.steady_state(r_idx);
        plot_data = 400 * (var_data - r_ss);  % Annualized deviation in percentage points
        ylabel_text = 'pp dev. (ann.)';
    elseif contains(var_name, {'TB'})
        % Trade balance: levels
        plot_data = var_data;
        ylabel_text = 'Level';
    else
        % Most variables: percentage deviation from steady state
        var_idx = strmatch(var_name, M_.endo_names, 'exact');
        if ~isempty(var_idx)
            var_ss = oo_.steady_state(var_idx);
            plot_data = 100 * (var_data - var_ss) / var_ss;
        else
            % Fallback: use first period
            plot_data = 100 * (var_data - var_data(1)) / var_data(1);
        end
        ylabel_text = '% dev. from SS';
    end
    
    % Plot time series
    plot(time, plot_data(1:length(time)), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    
    % Add zero line for deviation variables
    if ~contains(var_name, {'TB'}) && ~strcmp(var_name, 'pi') && ~strcmp(var_name, 'r')
        hold on;
        plot(time, zeros(size(time)), 'k--', 'LineWidth', 0.5);
    end
    
    % Formatting
    title(var_label, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Quarters');
    ylabel(ylabel_text);
    grid on;
    
    % Set appropriate y-axis limits
    if max(abs(plot_data)) > 0
        ylim_padding = max(abs(plot_data)) * 0.1;  % 10% padding
        ylim([min(plot_data) - ylim_padding, max(plot_data) + ylim_padding]);
    end
    
    % Improve appearance
    set(gca, 'FontSize', 9);
    axis tight;
end

% Add overall title
sgtitle('Macroeconomic Time Series - Key Aggregate Variables', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Add information box
info_text = sprintf(['Simulation horizon: %d quarters\n' ...
                     'Deviations from initial steady state'], ...
                    time_end, length(available_vars), size(key_vars, 1));

annotation('textbox', [0.02 0.02 0.25 0.08], ...
           'String', info_text, ...
           'FontSize', 9, 'EdgeColor', 'black', 'BackgroundColor', 'white');

% Make sure plot is displayed
drawnow;

% Save figure
try
    print(fig, '-dpng', '-r300', 'aggregate_timeseries.png');
    fprintf('Aggregate time series plot displayed and saved as aggregate_timeseries.png\n');
catch
    fprintf('Aggregate time series plot displayed (saving failed but plot is visible)\n');
end

% Display summary statistics
fprintf('\n=== TIME SERIES SUMMARY ===\n');
fprintf('Simulation horizon: %d quarters\n', time_end);
fprintf('Variables plotted: %d\n\n', length(available_vars));

for i = 1:length(available_vars)
    var_name = available_vars{i};
    var_label = available_labels{i};
    var_data = data.(var_name);
    
    if strcmp(var_name, 'pi')
        % Inflation: deviation from steady state in annualized pp
        pi_idx = strmatch('pi', M_.endo_names, 'exact');
        pi_ss = oo_.steady_state(pi_idx);
        dev_data = 400 * (var_data - pi_ss);
        mean_val = mean(dev_data);
        std_val = std(dev_data);
        fprintf('  %s: Mean dev = %.2f pp, Std = %.2f pp\n', var_label, mean_val, std_val);
    elseif strcmp(var_name, 'r')
        % Interest rate: deviation from steady state in annualized pp
        r_idx = strmatch('r', M_.endo_names, 'exact');
        r_ss = oo_.steady_state(r_idx);
        dev_data = 400 * (var_data - r_ss);
        mean_val = mean(dev_data);
        std_val = std(dev_data);
        fprintf('  %s: Mean dev = %.2f pp, Std = %.2f pp\n', var_label, mean_val, std_val);
    elseif contains(var_name, {'TB'})
        % Level variables
        mean_val = mean(var_data);
        std_val = std(var_data);
        fprintf('  %s: Mean = %.3f, Std = %.3f\n', var_label, mean_val, std_val);
    else
        % Percentage deviations from steady state
        var_idx = strmatch(var_name, M_.endo_names, 'exact');
        if ~isempty(var_idx)
            var_ss = oo_.steady_state(var_idx);
            dev_data = 100 * (var_data - var_ss) / var_ss;
        else
            dev_data = 100 * (var_data - var_data(1)) / var_data(1);
        end
        mean_val = mean(dev_data);
        std_val = std(dev_data);
        max_val = max(abs(dev_data));
        fprintf('  %s: Mean dev = %.2f%%, Std = %.2f%%, Max |dev| = %.2f%%\n', ...
                var_label, mean_val, std_val, max_val);
    end
end

fprintf('\nNote: "dev" = deviation from steady state; "pp" = percentage points\n');
