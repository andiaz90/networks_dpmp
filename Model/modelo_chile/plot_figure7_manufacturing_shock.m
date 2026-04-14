% plot_figure7_manufacturing_shock.m
% Figure 7: Impact of Manufacturing TFP Shock on Sectoral and Macro Variables
% Shows: Sector-level employment, consumption, imports, exports
%        Open economy variables: trade balance, exchange rate, interest rate

fprintf('\n=== GENERATING FIGURE 7: MANUFACTURING SHOCK EFFECTS ===\n');

% Manufacturing sector index
manufacturing_idx = 3;
shock_name = 'epsA_3';  % Manufacturing TFP shock

% Check if model has been solved
if ~exist('oo_', 'var') || ~isfield(oo_, 'irfs')
    error('Model not solved. Run main_SOE_exercise2.m first to generate IRFs.');
end

% Time horizon (quarters)
irf_fields = fieldnames(oo_.irfs);
sample_irf = oo_.irfs.(irf_fields{1});
T_full = length(sample_irf);
T = min(40, T_full);  % Show 40 quarters (10 years)
quarters = 0:T-1;

%% === FIGURE 1: MANUFACTURING SECTOR + OPEN ECONOMY VARIABLES ===
fprintf('Plotting Figure 7a: Manufacturing sector and open economy variables...\n');

% Define variables to plot for manufacturing sector + OE variables
manuf_oe_vars = {
    % Manufacturing sector variables (sector 3)
    'L_3',      'Manuf. Employment';
    'Y_3',      'Manuf. Output';
    'P_3',      'Manuf. Price';
    'Xi_3',     'Manuf. Exports';
    % Open Economy variables
    'TB',       'Trade Balance';
    'Q',        'Real Exchange Rate';
    'r',        'Interest Rate';
    'X',        'Total Exports'
};

fig1 = figure('Name', 'Figure 7a: Manufacturing Sector & Open Economy Variables', ...
              'Position', [50 50 1400 900], ...
              'Color', 'white', ...
              'Visible', 'on');
figure(fig1);  % Bring to front

n_vars = size(manuf_oe_vars, 1);
n_rows = 2;
n_cols = 4;

for i = 1:n_vars
    subplot(n_rows, n_cols, i);
    
    var_name = manuf_oe_vars{i, 1};
    var_label = manuf_oe_vars{i, 2};
    irf_name = [var_name '_' shock_name];
    
    % Special handling for computed variables
    if strcmp(var_name, 'Xi_3')
        % Compute manufacturing exports: Xi_3 = chiX_3 * X * PX / PH_3
        % Need IRFs for X, PX, and PH_3
        X_irf_name = ['X_' shock_name];
        PX_irf_name = ['PX_' shock_name];
        PH3_irf_name = ['PH_3_' shock_name];
        
        if ismember(X_irf_name, irf_fields) && ismember(PX_irf_name, irf_fields) && ismember(PH3_irf_name, irf_fields)
            % Get IRFs (in log deviations)
            X_irf = oo_.irfs.(X_irf_name)(1:T);
            PX_irf = oo_.irfs.(PX_irf_name)(1:T);
            PH3_irf = oo_.irfs.(PH3_irf_name)(1:T);
            
            % Xi_3 = chiX_3 * X * PX / PH_3
            % In log: log(Xi_3) ≈ log(X) + log(PX) - log(PH_3) + constant
            % For deviations: Xi_3_dev ≈ X_dev + PX_dev - PH3_dev
            irf_data = X_irf + PX_irf - PH3_irf;
        else
            irf_data = [];
        end
    elseif ismember(irf_name, irf_fields)
        irf_data = oo_.irfs.(irf_name)(1:T);
    else
        irf_data = [];
    end
    
    % Plot if data is available
    if ~isempty(irf_data)
        % Color code: red for manufacturing, green for OE variables
        if contains(var_name, '_3') || strcmp(var_name, 'Xi_3')
            line_color = [0.8 0.1 0.1];  % Red for manufacturing
            line_width = 2.8;
        else
            line_color = [0.1 0.5 0.2];  % Green for OE variables
            line_width = 2.5;
        end
        
        % Plot
        plot(quarters, 100*irf_data, 'LineWidth', line_width, 'Color', line_color);
        hold on;
        plot(quarters, zeros(size(quarters)), 'k--', 'LineWidth', 0.7);
        
        % Formatting
        title(var_label, 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Quarters', 'FontSize', 10);
        
        if ismember(var_name, {'r'})
            ylabel('Percentage Points', 'FontSize', 10);
        elseif strcmp(var_name, 'TB')
            ylabel('% of GDP', 'FontSize', 10);
        else
            ylabel('% Deviation', 'FontSize', 10);
        end
        
        grid on;
        xlim([0 T-1]);
        enforce_min_yaxis(gca);
        
        % Add peak value annotation
        [peak_val, peak_idx] = max(abs(irf_data));
        peak_impact = 100*irf_data(peak_idx);
        peak_quarter = peak_idx - 1;
        text(0.98, 0.95, sprintf('Peak: %.2f%% (Q%d)', peak_impact, peak_quarter), ...
             'Units', 'normalized', 'FontSize', 8, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'white', 'EdgeColor', 'none');
        
    else
        % Variable not available
        text(0.5, 0.5, sprintf('%s\nNot Available', var_label), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
        axis off;
    end
end

% Main title
sgtitle({'Figure 7a: Manufacturing TFP Shock',...
         'Effects on Manufacturing Sector and Open Economy Variables'}, ...
        'FontSize', 15, 'FontWeight', 'bold');

% Force display
drawnow;

% Save figure
print(fig1, '-dpng', '-r300', 'figure7a_manufacturing_and_OE.png');
fprintf('✓ Figure 7a saved: figure7a_manufacturing_and_OE.png\n');

%% === FIGURE 2: AGGREGATE MACROECONOMIC VARIABLES ===
fprintf('Plotting Figure 7b: Aggregate macroeconomic variables...\n');

% Aggregate variables to plot
agg_vars = {
    'N',        'Aggregate Employment';
    'C',        'Aggregate Consumption';
    'Y',        'Aggregate Output (GDP)';
    'pi',       'Aggregate Inflation';
    'pi_g',     'Goods Inflation';
    'pi_s',     'Services Inflation';
    'r',        'Interest Rate';
    'w',        'Real Wage';
    'Q',        'Real Exchange Rate';
    'TB',       'Trade Balance';
    'X',        'Total Exports';
    'IMP',      'Total Imports'
};

% Create figure for macro variables
fig2 = figure('Name', 'Figure 7b: Aggregate Macroeconomic Effects', ...
              'Position', [100 100 1600 900], ...
              'Color', 'white', ...
              'Visible', 'on');
figure(fig2);  % Bring to front

n_agg = size(agg_vars, 1);
n_rows = 3;
n_cols = 4;

for i = 1:n_agg
    subplot(n_rows, n_cols, i);
    
    var_name = agg_vars{i, 1};
    var_label = agg_vars{i, 2};
    irf_name = [var_name '_' shock_name];
    
    % Check if IRF exists
    if ismember(irf_name, irf_fields)
        irf_data = oo_.irfs.(irf_name)(1:T);
        
        % Plot IRF with blue color for aggregate variables
        plot(quarters, 100*irf_data, 'LineWidth', 2.5, 'Color', [0.1 0.3 0.7]);
        hold on;
        
        % Zero line
        plot(quarters, zeros(size(quarters)), 'k--', 'LineWidth', 0.8);
        
        % Formatting
        title(var_label, 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Quarters', 'FontSize', 10);
        
        % Y-axis label depends on variable type
        if ismember(var_name, {'pi', 'r'})
            ylabel('Percentage Points', 'FontSize', 10);
        elseif ismember(var_name, {'TB', 'NX'})
            ylabel('% of GDP', 'FontSize', 10);
        else
            ylabel('% Deviation', 'FontSize', 10);
        end
        
        grid on;
        xlim([0 T-1]);
        enforce_min_yaxis(gca);
        
        % Add impact statistics
        [peak_val, peak_idx] = max(abs(irf_data));
        peak_impact = 100 * irf_data(peak_idx);
        peak_quarter = peak_idx - 1;
        impact_q1 = 100 * irf_data(min(2, length(irf_data)));
        
        text(0.98, 0.95, sprintf('Q1: %.2f%%\nPeak: %.2f%% (Q%d)', ...
             impact_q1, peak_impact, peak_quarter), ...
             'Units', 'normalized', 'FontSize', 7.5, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'white', 'EdgeColor', 'none');
        
    else
        % Variable not available
        text(0.5, 0.5, sprintf('%s\nNot Available', var_label), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
        axis off;
    end
end

% Main title
sgtitle({'Figure 7b: Manufacturing TFP Shock',...
         'Effects on Aggregate Macroeconomic Variables'}, ...
        'FontSize', 15, 'FontWeight', 'bold');

% Force display
drawnow;

% Save figure
print(fig2, '-dpng', '-r300', 'figure7b_aggregate_effects.png');
fprintf('✓ Figure 7b saved: figure7b_aggregate_effects.png\n');

%% === FIGURE 3: SECTORAL IMPORTS ===
fprintf('Plotting Figure 7c: Sectoral imports (foreign consumption)...\n');

% Create figure for sectoral imports
fig3 = figure('Name', 'Figure 7c: Sectoral Imports', ...
              'Position', [150 150 1600 900], ...
              'Color', 'white', ...
              'Visible', 'on');
figure(fig3);  % Bring to front

n_rows = 3;
n_cols = 4;

% Get sector names from workspace
if exist('names', 'var')
    sector_names = names;
else
    sector_names = cell(nsec,1);
    for i = 1:nsec
        sector_names{i} = sprintf('Sector %d', i);
    end
end

for i = 1:nsec
    subplot(n_rows, n_cols, i);
    
    % Sectoral imports = CFg_i + CFs_i (foreign goods + foreign services consumption)
    CFg_irf_name = sprintf('CFg_%d_%s', i, shock_name);
    CFs_irf_name = sprintf('CFs_%d_%s', i, shock_name);
    
    % Check if both components exist
    if ismember(CFg_irf_name, irf_fields) && ismember(CFs_irf_name, irf_fields)
        CFg_data = oo_.irfs.(CFg_irf_name)(1:T);
        CFs_data = oo_.irfs.(CFs_irf_name)(1:T);
        
        % For small deviations, we can approximate total imports deviation as sum
        % More accurate would be to compute weighted average, but this is close enough
        import_data = CFg_data + CFs_data;
        
        % Highlight manufacturing sector
        if i == manufacturing_idx
            line_color = [0.8 0.1 0.1];  % Red for manufacturing
            line_width = 2.8;
        else
            line_color = [0.2 0.4 0.7];  % Blue for others
            line_width = 2.0;
        end
        
        % Plot
        plot(quarters, 100*import_data, 'LineWidth', line_width, 'Color', line_color);
        hold on;
        plot(quarters, zeros(size(quarters)), 'k--', 'LineWidth', 0.7);
        
        % Formatting
        sector_label = sector_names{i};
        if length(sector_label) > 20
            sector_label = sector_label(1:20);
        end
        title(sprintf('%s', sector_label), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel('Quarters', 'FontSize', 9);
        ylabel('% Deviation', 'FontSize', 9);
        
        grid on;
        xlim([0 T-1]);
        enforce_min_yaxis(gca);
        
        % Add peak value annotation
        [peak_val, peak_idx] = max(abs(import_data));
        peak_impact = 100*import_data(peak_idx);
        peak_quarter = peak_idx - 1;
        text(0.98, 0.95, sprintf('Peak: %.2f%%', peak_impact), ...
             'Units', 'normalized', 'FontSize', 7, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'white', 'EdgeColor', 'none');
        
    else
        % Variable not available
        text(0.5, 0.5, sprintf('Sector %d\nNot Available', i), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.5 0.5 0.5]);
        axis off;
    end
end

% Main title
sgtitle({'Figure 7c: Manufacturing TFP Shock',...
         'Effects on Sectoral Imports (Foreign Consumption)'}, ...
        'FontSize', 15, 'FontWeight', 'bold');

% Force display
drawnow;

% Save figure
print(fig3, '-dpng', '-r300', 'figure7c_sectoral_imports.png');
fprintf('✓ Figure 7c saved: figure7c_sectoral_imports.png\n');

%% === SUMMARY STATISTICS ===
fprintf('\n=== MANUFACTURING TFP SHOCK: IMPACT SUMMARY ===\n');
fprintf('Horizon: %d quarters (%.1f years)\n', T, T/4);
fprintf('Shock: 1%% positive TFP shock to Manufacturing (Sector 3)\n');
fprintf('Persistence: rho = %.2f\n\n', rho_tfp1_val);

% Combine variables from both figures for summary
all_summary_vars = [manuf_oe_vars; agg_vars];

fprintf('%-35s %12s %12s %12s\n', 'Variable', 'Impact (Q1)', 'Peak', 'At Quarter');
fprintf('%s\n', repmat('-', 1, 75));

for i = 1:size(all_summary_vars, 1)
    var_name = all_summary_vars{i, 1};
    var_label = all_summary_vars{i, 2};
    irf_name = [var_name '_' shock_name];
    
    if ismember(irf_name, irf_fields)
        irf_data = oo_.irfs.(irf_name);
        
        impact_q1 = 100 * irf_data(min(2, length(irf_data)));  % Quarter 1
        [peak_val, peak_idx] = max(abs(irf_data));
        peak_impact = 100 * irf_data(peak_idx);
        peak_quarter = peak_idx - 1;
        
        fprintf('%-35s %11.3f%% %11.3f%% %12d\n', ...
                var_label, impact_q1, peak_impact, peak_quarter);
    end
end

fprintf('\n✓ Figure 7 generation complete!\n');
fprintf('  - Figure 7a: figure7a_manufacturing_and_OE.png (Manufacturing + Open Economy)\n');
fprintf('  - Figure 7b: figure7b_aggregate_effects.png (Aggregate Variables)\n\n');
