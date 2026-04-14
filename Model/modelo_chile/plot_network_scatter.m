% plot_network_scatter.m
% Creates scatter plots showing the relationship between upstream connectedness 
% and various economic responses to shocks
%
% Generates scatter plots for:
% 1. Consumption response
% 2. Sector price
% 3. Production
% 4. Home consumption
% 5. Foreign consumption
% 6. Sector imports
% 7. Sector exports

%% Setup
if ~exist('nsec', 'var')
    nsec = 12;
end

if ~exist('names', 'var')
    warning('Sector names not loaded. Using generic names.');
    names = arrayfun(@(i) sprintf('Sector %d', i), 1:nsec, 'UniformOutput', false);
end

if ~exist('modbeta', 'var')
    error('I-O matrix (modbeta) not found. Please run main_SOE.m first.');
end

%% Calculate Upstream Connectedness
% Upstream connectedness: how much sector i uses inputs from all other sectors
% This is based on the Leontief inverse: (I - Beta)^(-1)
% where Beta is the I-O matrix

fprintf('\n=== Calculating Upstream Connectedness ===\n');

% Method 1: Simple input shares (direct upstream linkages)
upstream_direct = sum(modbeta, 1)';  % Sum across rows for each column (sector j's input use)

% Method 2: Leontief inverse (total upstream linkages including indirect)
I_matrix = eye(nsec);
Leontief_inverse = inv(I_matrix - modbeta);
upstream_total = sum(Leontief_inverse, 1)' - 1;  % Subtract 1 to exclude own sector

fprintf('Upstream connectedness (direct): min=%.3f, max=%.3f, mean=%.3f\n', ...
    min(upstream_direct), max(upstream_direct), mean(upstream_direct));
fprintf('Upstream connectedness (total): min=%.3f, max=%.3f, mean=%.3f\n', ...
    min(upstream_total), max(upstream_total), mean(upstream_total));

% Use total upstream connectedness for analysis
upstream_conn = upstream_total;

%% Identify shock and time horizon for analysis
shock_name = 'epsA_3';  % Default: Manufacturing TFP shock
time_point = 1;  % Impact response (quarter 0)

if exist('oo_', 'var') && isfield(oo_, 'irfs')
    irf_fields = fieldnames(oo_.irfs);
    
    % Try to identify available shock
    if any(contains(irf_fields, 'epsA_3'))
        shock_name = 'epsA_3';
    elseif any(contains(irf_fields, 'epsA'))
        for i = 1:length(irf_fields)
            tokens = regexp(irf_fields{i}, '_(epsA_\d+)$', 'tokens');
            if ~isempty(tokens)
                shock_name = tokens{1}{1};
                break;
            end
        end
    elseif any(contains(irf_fields, 'eps_i'))
        shock_name = 'eps_i';
    elseif any(contains(irf_fields, 'eps_om'))
        shock_name = 'eps_om';
    end
    fprintf('Using shock: %s\n', shock_name);
    fprintf('Analyzing impact response (quarter %d)\n', time_point-1);
else
    error('IRF results (oo_) not found. Please run the model first.');
end

%% Extract responses for each variable
fprintf('\n=== Extracting Sectoral Responses ===\n');

% Initialize response vectors
C_response = nan(nsec, 1);      % Total consumption
P_response = nan(nsec, 1);      % Prices
Y_response = nan(nsec, 1);      % Production
CH_response = nan(nsec, 1);     % Home consumption
CF_response = nan(nsec, 1);     % Foreign consumption
IMP_response = nan(nsec, 1);    % Imports (= CF)
% Note: Model has aggregate exports (X) only, not sectoral exports (X_i)

for i = 1:nsec
    % Price (PH_i - corrected variable name)
    var_name = sprintf('PH_%d_%s', i, shock_name);
    if isfield(oo_.irfs, var_name)
        P_response(i) = oo_.irfs.(var_name)(time_point);
    end
    
    % Production (Y_i)
    var_name = sprintf('Y_%d_%s', i, shock_name);
    if isfield(oo_.irfs, var_name)
        Y_response(i) = oo_.irfs.(var_name)(time_point);
    end
    
    % Home consumption (CHg_i or CHs_i depending on sector type)
    if goods(i)
        var_name = sprintf('CHg_%d_%s', i, shock_name);
    else
        var_name = sprintf('CHs_%d_%s', i, shock_name);
    end
    if isfield(oo_.irfs, var_name)
        CH_response(i) = oo_.irfs.(var_name)(time_point);
    end
    
    % Foreign consumption (CFg_i or CFs_i depending on sector type)
    % This equals imports for each sector
    if goods(i)
        var_name = sprintf('CFg_%d_%s', i, shock_name);
    else
        var_name = sprintf('CFs_%d_%s', i, shock_name);
    end
    if isfield(oo_.irfs, var_name)
        CF_response(i) = oo_.irfs.(var_name)(time_point);
        IMP_response(i) = CF_response(i);  % Imports = foreign consumption
    end
    
    % Total consumption per sector = Home + Foreign
    if ~isnan(CH_response(i)) && ~isnan(CF_response(i))
        C_response(i) = CH_response(i) + CF_response(i);
    elseif ~isnan(CH_response(i))
        C_response(i) = CH_response(i);
    elseif ~isnan(CF_response(i))
        C_response(i) = CF_response(i);
    end
end

fprintf('Valid responses: C=%d, P=%d, Y=%d, CH=%d, CF=%d, IMP=%d\n', ...
    sum(~isnan(C_response)), sum(~isnan(P_response)), sum(~isnan(Y_response)), ...
    sum(~isnan(CH_response)), sum(~isnan(CF_response)), ...
    sum(~isnan(IMP_response)));

%% Create scatter plots with linear fits
fprintf('\n=== Creating Scatter Plots ===\n');

fig = figure('Name', 'Network Connectedness and Economic Responses', ...
             'Position', [50, 50, 1600, 800], 'Color', 'white');

% Define variables to plot (6 plots total - no sectoral exports in model)
var_data = {C_response, P_response, Y_response, CH_response, CF_response, IMP_response};
var_names = {'Consumption (C_i)', 'Price (p^H_i)', 'Production (Y_i)', ...
             'Home Consumption (C^H_i)', 'Foreign Consumption (C^F_i)', ...
             'Imports (IMP_i)'};
var_labels = {'Consumption Response (%)', 'Price Response (%)', 'Production Response (%)', ...
              'Home Cons. Response (%)', 'Foreign Cons. Response (%)', ...
              'Import Response (%)'};

% Create subplots (2 rows x 3 columns)
for idx = 1:6
    subplot(2, 3, idx);
    
    response = var_data{idx};
    valid_idx = ~isnan(response) & ~isnan(upstream_conn);
    
    if sum(valid_idx) < 2
        text(0.5, 0.5, 'Insufficient data', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'FontSize', 10);
        title(var_names{idx}, 'FontSize', 11, 'FontWeight', 'bold');
        continue;
    end
    
    % Extract valid data
    x_data = upstream_conn(valid_idx);
    y_data = 100 * response(valid_idx);  % Convert to percentage
    sector_labels = names(valid_idx);  % Get sector names for valid data
    
    % Create scatter plot
    scatter(x_data, y_data, 80, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    
    % Add labels for each point
    for i = 1:length(x_data)
        text(x_data(i), y_data(i), ['  ' sector_labels{i}], ...
             'FontSize', 8, 'VerticalAlignment', 'middle', ...
             'HorizontalAlignment', 'left', 'Interpreter', 'none');
    end
    
    % Fit linear regression
    p = polyfit(x_data, y_data, 1);  % Linear fit
    y_fit = polyval(p, x_data);
    
    % Plot regression line
    [x_sorted, sort_idx] = sort(x_data);
    plot(x_sorted, y_fit(sort_idx), 'r-', 'LineWidth', 2);
    
    % Calculate R-squared
    y_mean = mean(y_data);
    SS_tot = sum((y_data - y_mean).^2);
    SS_res = sum((y_data - y_fit).^2);
    R_squared = 1 - SS_res/SS_tot;
    
    % Calculate correlation
    corr_coef = corr(x_data, y_data);
    
    % Formatting
    grid on;
    xlabel('Upstream Connectedness', 'FontSize', 10);
    ylabel(var_labels{idx}, 'FontSize', 10);
    title(var_names{idx}, 'FontSize', 11, 'FontWeight', 'bold');
    hold off;
end

% Add overall title
sgtitle(sprintf('Network Connectedness vs. Economic Responses to %s Shock (Impact)', ...
        strrep(shock_name, '_', '\_')), 'FontSize', 14, 'FontWeight', 'bold');

%% Save figure
saveas(fig, 'network_scatter_plots.png');
fprintf('\n✓ Scatter plots saved: network_scatter_plots.png\n');

%% Print summary statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('%-25s | Corr. | R² | Slope\n', 'Variable');
fprintf('-----------------------------------------------------------\n');

for idx = 1:6
    response = var_data{idx};
    valid_idx = ~isnan(response) & ~isnan(upstream_conn);
    
    if sum(valid_idx) >= 2
        x_data = upstream_conn(valid_idx);
        y_data = 100 * response(valid_idx);
        
        corr_coef = corr(x_data, y_data);
        p = polyfit(x_data, y_data, 1);
        
        y_mean = mean(y_data);
        SS_tot = sum((y_data - y_mean).^2);
        y_fit = polyval(p, x_data);
        SS_res = sum((y_data - y_fit).^2);
        R_squared = 1 - SS_res/SS_tot;
        
        fprintf('%-25s | %5.3f | %4.3f | %7.3f\n', ...
            var_names{idx}, corr_coef, R_squared, p(1));
    else
        fprintf('%-25s | N/A   | N/A  | N/A\n', var_names{idx});
    end
end

fprintf('\n');
