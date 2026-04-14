
    % Define key variables to plot
    key_vars = {
        'C',        'Consumption';
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
        'IMP',     'Imports'
    };
    
    % Define shocks to compare
    shocks = {
        'eps_om', 'Consumption Preference Shock', [0.2 0.4 0.8];  % Blue
        'eps_i', 'Monetary Policy Shock', [0.8 0.2 0.2]   % Red
    };
    
    % Check which IRFs are available for both shocks
    irf_fields = fieldnames(oo_.irfs);
    available_vars = {};
    available_labels = {};
    
    % Check which variables have IRFs for both shocks
    for i = 1:size(key_vars, 1)
        var_name = key_vars{i, 1};
        has_both_shocks = true;
        
        % Check if both shocks have IRFs for this variable
        for j = 1:size(shocks, 1)
            shock_name = shocks{j, 1};
            irf_name = [var_name '_' shock_name];
            if ~ismember(irf_name, irf_fields)
                has_both_shocks = false;
                break;
            end
        end
        
        if has_both_shocks
            available_vars{end+1} = var_name;
            available_labels{end+1} = key_vars{i, 2};
        else
            fprintf('Warning: IRF for %s not available for both shocks.\n', var_name);
        end
    end
    
    if isempty(available_vars)
        error('No IRFs found for both consumption and monetary policy shocks. Check model solution.');
    end
    
    % Create figure - make it visible and bring to front
    n_vars = length(available_vars);
    n_rows = ceil(n_vars / 3);
    n_cols = min(3, n_vars);
    
    fig = figure('Name', 'Comparison: Consumption vs Monetary Policy Shocks', ...
                 'Position', [100 100 1400 900], ...
                 'Visible', 'on', ...
                 'WindowStyle', 'docked');  % Use docked for better visibility
    
    % Bring figure to front
    figure(fig);
    drawnow;
    
    % Time horizon (quarters) - use first available shock
    first_shock = shocks{1, 1};
    T_full = length(oo_.irfs.([available_vars{1} '_' first_shock]));
    T = min(20, T_full);  % Limit to 20 periods
    quarters = 0:T-1;
    
    % Plot each variable
    for i = 1:n_vars
        subplot(n_rows, n_cols, i);
        
        var_name = available_vars{i};
        max_abs_val = 0;
        
        % Plot both shocks for this variable
        for j = 1:size(shocks, 1)
            shock_name = shocks{j, 1};
            shock_label = shocks{j, 2};
            shock_color = shocks{j, 3};
            
            irf_name = [var_name '_' shock_name];
            irf_data = oo_.irfs.(irf_name)(1:T);  % Use only first T periods
            
            % Plot IRF
            plot(quarters, 100*irf_data, 'LineWidth', 2, 'Color', shock_color, ...
                 'DisplayName', shock_label);
            hold on;
            
            % Track maximum value for y-axis scaling
            max_abs_val = max(max_abs_val, max(abs(irf_data)));
        end
        
        % Add zero line
        plot(quarters, zeros(size(quarters)), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        
        % Formatting
        title(available_labels{i}, 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Quarters');
        if contains(var_name, {'pi', 'r'})
            ylabel('Percentage Points');
        else
            ylabel('Percent Deviation');
        end
        
        grid on;
        xlim([0 T-1]);
        
        % Add legend for first subplot only
        if i == 1
            legend('Location', 'best', 'FontSize', 10);
        end
        
        % Add some visual improvements
        if max_abs_val > 0
            ylim_range = max_abs_val * 110;  % 110% of max for some padding
            ylim([-ylim_range ylim_range]);
        end
        enforce_min_yaxis(gca);
    end

    % Add main title
    sgtitle('Impulse Response Functions: Consumption vs Monetary Policy Shocks',...
            'FontSize', 16, 'FontWeight', 'bold');
    % Add shock information (smaller textbox, bottom-right)
    annotation('textbox', [0.70 0.02 0.28 0.07], ... % [x y w h] in normalized units
               'String', sprintf(['Blue: Consumption Shock (eps_C, \\sigma = %.2f)\n' ...
                                  'Red: Monetary Policy Shock (eps_i, \\sigma = %.2f)'], ...
                                  0.01, 0.01), ...
               'Units', 'normalized', ...
               'FontSize', 8, ...
               'EdgeColor', 'black', ...
               'BackgroundColor', 'white', ...
               'HorizontalAlignment', 'left', ...
               'VerticalAlignment', 'bottom', ...
               'FitBoxToText', 'on');
    
    % Make sure the figure is displayed and refresh
    drawnow;  % Force MATLAB to update the display immediately
    
    % Save figure (optional - keep for record)
    try
        print(fig, '-dpng', '-r300', 'shock_comparison_effects.png');
        fprintf('Plot displayed on screen and saved as shock_comparison_effects.png\n');
    catch
        fprintf('Plot displayed on screen (saving failed but plot is visible)\n');
    end
    
    % Display comparative statistics
    fprintf('\n=== SHOCK COMPARISON SUMMARY ===\n');
    fprintf('Time horizon: %d quarters\n\n', T);
    
    for k = 1:size(shocks, 1)
        shock_name = shocks{k, 1};
        shock_label = shocks{k, 2};
        
        fprintf('%s:\n', shock_label);
        fprintf('================\n');
        
        for i = 1:min(6, n_vars)  % Show first 6 variables
            var_name = available_vars{i};
            irf_name = [var_name '_' shock_name];
            irf_data = oo_.irfs.(irf_name);
            
            peak_impact = 100 * irf_data(abs(irf_data) == max(abs(irf_data)));
            peak_impact = peak_impact(1);  % In case of ties
            peak_quarter = find(abs(irf_data) == max(abs(irf_data)), 1) - 1;
            
            fprintf('  %s: Peak impact = %.3f%% at quarter %d\n', ...
                    available_labels{i}, peak_impact, peak_quarter);
        end
        fprintf('\n');
    end
