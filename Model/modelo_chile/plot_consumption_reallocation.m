% plot_consumption_reallocation.m
% Plots sectoral consumption reallocation: HOME and FOREIGN separately
% Shows CHg/CHs (home) and CFg/CFs (foreign) for each sector

% Close any existing consumption reallocation figures to prevent accumulation
close(findall(0, 'Name', 'Home Consumption Reallocation'));
close(findall(0, 'Name', 'Foreign Consumption (Imports)'));

%% Setup
if ~exist('nsec', 'var')
    nsec = 12;
end

if ~exist('names', 'var')
    warning('Sector names not loaded. Using generic names.');
    names = arrayfun(@(i) sprintf('Sector %d', i), 1:nsec, 'UniformOutput', false);
end

if ~exist('goods', 'var') || ~exist('services', 'var')
    warning('Sector classification not found. Using default.');
    goods = [1 1 1 1 1 0 0 0 0 0 0 0]';
    services = ~goods;
end

% Shorten sector names for plotting
short_names = cell(nsec, 1);
for i = 1:nsec
    if length(names{i}) > 25
        short_names{i} = names{i}(1:22);
    else
        short_names{i} = names{i};
    end
end

T_plot = 40;  % quarters to plot

%% Identify shock
shock_name = 'epsA_3';  % Manufacturing TFP shock (default)
if exist('oo_', 'var') && isfield(oo_, 'irfs')
    irf_fields = fieldnames(oo_.irfs);
    % Try to find which shock is available - prioritize sectoral TFP shocks
    if any(contains(irf_fields, 'epsA_3'))
        shock_name = 'epsA_3';  % Manufacturing TFP shock
    elseif any(contains(irf_fields, 'epsA'))
        % Find any sectoral TFP shock
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
end

fprintf('\n=== Plotting HOME Consumption Reallocation ===\n');

%% ========================================
%% FIGURE 1: HOME CONSUMPTION (CHg + CHs)
%% ========================================

fig1 = figure('Name', 'Home Consumption Reallocation', 'Position', [50, 50, 1600, 900], ...
              'Color', 'white', 'Visible', 'on');
drawnow;

if exist('oo_', 'var') && isfield(oo_, 'irfs')
    irf_names = fieldnames(oo_.irfs);
    irf_length = length(oo_.irfs.(irf_names{1}));
    T_plot = min(T_plot, irf_length);
    
    % Panel 1-12: Sectoral home consumption
    for i = 1:nsec
        subplot(3, 5, i);
        hold on; grid on;
        
        % Choose variable based on sector type
        if goods(i)
            var_name = sprintf('CHg_%d', i);  % Home goods consumption
            line_color = [0.2, 0.4, 0.7];  % Blue for goods
        else
            var_name = sprintf('CHs_%d', i);  % Home services consumption
            line_color = [0.8, 0.3, 0.1];  % Red for services
        end
        
        irf_field = sprintf('%s_%s', var_name, shock_name);
        
        if isfield(oo_.irfs, irf_field)
            irf_data = oo_.irfs.(irf_field);
            plot(0:length(irf_data)-1, 100*irf_data, 'LineWidth', 1.8, 'Color', line_color);
        else
            text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'FontSize', 8);
        end
        
        yline(0, 'k--', 'LineWidth', 0.5);
        title(short_names{i}, 'FontSize', 9, 'FontWeight', 'bold');
        xlabel('Quarters', 'FontSize', 8);
        ylabel('Home Cons. (% dev.)', 'FontSize', 8);
        xlim([0, T_plot-1]);
        set(gca, 'FontSize', 7);
    end
    
    % Panel 13: Model aggregate home goods (C_g)
    subplot(3, 5, 13);
    hold on; grid on;
    if isfield(oo_.irfs, ['C_g_' shock_name])
        C_g_data = oo_.irfs.(['C_g_' shock_name]);
        plot(0:length(C_g_data)-1, 100*C_g_data, 'LineWidth', 3, 'Color', [0.2, 0.4, 0.7]);
    end
    yline(0, 'k--', 'LineWidth', 0.5);
    title('Aggregate: Home Goods (C_g)', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Quarters', 'FontSize', 8);
    ylabel('% deviation', 'FontSize', 8);
    xlim([0, T_plot-1]);
    set(gca, 'FontSize', 8);
    
    % Panel 14: Model aggregate home services (C_s)
    subplot(3, 5, 14);
    hold on; grid on;
    if isfield(oo_.irfs, ['C_s_' shock_name])
        C_s_data = oo_.irfs.(['C_s_' shock_name]);
        plot(0:length(C_s_data)-1, 100*C_s_data, 'LineWidth', 3, 'Color', [0.8, 0.3, 0.1]);
    end
    yline(0, 'k--', 'LineWidth', 0.5);
    title('Aggregate: Home Services (C_s)', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Quarters', 'FontSize', 8);
    ylabel('% deviation', 'FontSize', 8);
    xlim([0, T_plot-1]);
    set(gca, 'FontSize', 8);
    
    % Panel 15: Model aggregate total consumption (C)
    subplot(3, 5, 15);
    hold on; grid on;
    if isfield(oo_.irfs, ['C_' shock_name])
        C_data = oo_.irfs.(['C_' shock_name]);
        plot(0:length(C_data)-1, 100*C_data, 'LineWidth', 3, 'Color', [0.2, 0.6, 0.2]);
    end
    yline(0, 'k--', 'LineWidth', 0.5);
    title('Aggregate: Total Consumption (C)', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Quarters', 'FontSize', 8);
    ylabel('% deviation', 'FontSize', 8);
    xlim([0, T_plot-1]);
    set(gca, 'FontSize', 8);
    
    sgtitle('HOME CONSUMPTION: By Sector (CHg for goods, CHs for services)', ...
            'FontSize', 14, 'FontWeight', 'bold');
end

% Save
print(fig1, '-dpng', '-r300', 'home_consumption_reallocation.png');
fprintf('✓ Home consumption figure saved: home_consumption_reallocation.png\n');

%% ========================================
%% FIGURE 2: FOREIGN CONSUMPTION (CFg + CFs)
%% ========================================

fprintf('\n=== Plotting FOREIGN Consumption (Imports) ===\n');

fig2 = figure('Name', 'Foreign Consumption (Imports)', 'Position', [100, 100, 1600, 900], ...
              'Color', 'white', 'Visible', 'on');
drawnow;

if exist('oo_', 'var') && isfield(oo_, 'irfs')
    
    % Get IRF length
    irf_names = fieldnames(oo_.irfs);
    irf_length = length(oo_.irfs.(irf_names{1}));
    T_plot = min(T_plot, irf_length);
    
    % Compute aggregate foreign goods and services FIRST before plotting
    CFg_total = zeros(1, irf_length);  % Row vector to match IRF structure
    CFs_total = zeros(1, irf_length);  % Row vector to match IRF structure
    
    for i = 1:nsec
        if goods(i)
            CFg_field = sprintf('CFg_%d_%s', i, shock_name);
            if isfield(oo_.irfs, CFg_field)
                CFg_total = CFg_total + oo_.irfs.(CFg_field);
                fprintf('  Added CFg_%d to aggregate (size: %d)\n', i, length(oo_.irfs.(CFg_field)));
            end
        else
            CFs_field = sprintf('CFs_%d_%s', i, shock_name);
            if isfield(oo_.irfs, CFs_field)
                CFs_total = CFs_total + oo_.irfs.(CFs_field);
                fprintf('  Added CFs_%d to aggregate (size: %d)\n', i, length(oo_.irfs.(CFs_field)));
            end
        end
    end
    
    fprintf('CFg_total range: [%.6f, %.6f]\n', min(CFg_total), max(CFg_total));
    fprintf('CFs_total range: [%.6f, %.6f]\n', min(CFs_total), max(CFs_total));
    
    % Panel 1-12: Sectoral foreign consumption
    for i = 1:nsec
        subplot(3, 5, i);
        hold off;  % Ensure clean start
        
        % Choose variable based on sector type
        if goods(i)
            var_name = sprintf('CFg_%d', i);  % Foreign goods consumption (imports)
            line_color = [0.3, 0.6, 0.3];  % Green for goods imports
        else
            var_name = sprintf('CFs_%d', i);  % Foreign services consumption (imports)
            line_color = [0.9, 0.5, 0.1];  % Orange for services imports
        end
        
        irf_field = sprintf('%s_%s', var_name, shock_name);
        
        if isfield(oo_.irfs, irf_field)
            irf_data = oo_.irfs.(irf_field);
            plot(0:length(irf_data)-1, 100*irf_data, 'LineWidth', 1.8, 'Color', line_color);
            hold on;  % Turn on only for adding yline
        else
            text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center', 'FontSize', 8);
        end
        
        yline(0, 'k--', 'LineWidth', 0.5);
        hold off;  % Turn off hold immediately
        grid on;
        title(short_names{i}, 'FontSize', 9, 'FontWeight', 'bold');
        xlabel('Quarters', 'FontSize', 8);
        ylabel('Foreign Cons. (% dev.)', 'FontSize', 8);
        xlim([0, T_plot-1]);
        set(gca, 'FontSize', 7);
    end
    
    % Panel 13: Aggregate foreign goods (sum of CFg across sectors)
    subplot(3, 5, 13);
    hold off;  % Clear any previous plots
    plot(0:length(CFg_total)-1, 100*CFg_total, 'LineWidth', 3, 'Color', [0.3, 0.6, 0.3]);
    hold on;
    yline(0, 'k--', 'LineWidth', 0.5);
    hold off;
    grid on;
    title('Aggregate: Foreign Goods', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Quarters', 'FontSize', 8);
    ylabel('% deviation', 'FontSize', 8);
    xlim([0, T_plot-1]);
    set(gca, 'FontSize', 8);
    
    % Panel 14: Aggregate foreign services (sum of CFs across sectors)
    subplot(3, 5, 14);
    hold off;  % Clear any previous plots
    plot(0:length(CFs_total)-1, 100*CFs_total, 'LineWidth', 3, 'Color', [0.9, 0.5, 0.1]);
    hold on;
    yline(0, 'k--', 'LineWidth', 0.5);
    hold off;
    grid on;
    title('Aggregate: Foreign Services', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Quarters', 'FontSize', 8);
    ylabel('% deviation', 'FontSize', 8);
    xlim([0, T_plot-1]);
    set(gca, 'FontSize', 8);
    
    % Panel 15: Total imports (IMP)
    subplot(3, 5, 15);
    hold on; grid on;
    if isfield(oo_.irfs, ['IMP_' shock_name])
        IMP_data = oo_.irfs.(['IMP_' shock_name]);
        plot(0:length(IMP_data)-1, 100*IMP_data, 'LineWidth', 3, 'Color', [0.4, 0.4, 0.4]);
    end
    yline(0, 'k--', 'LineWidth', 0.5);
    title('Aggregate: Total Imports (IMP)', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Quarters', 'FontSize', 8);
    ylabel('% deviation', 'FontSize', 8);
    xlim([0, T_plot-1]);
    set(gca, 'FontSize', 8);
    
    sgtitle('FOREIGN CONSUMPTION (IMPORTS): By Sector (CFg for goods, CFs for services)', ...
            'FontSize', 14, 'FontWeight', 'bold');
end

% Save
print(fig2, '-dpng', '-r300', 'foreign_consumption_reallocation.png');
fprintf('✓ Foreign consumption figure saved: foreign_consumption_reallocation.png\n');

fprintf('\n=== Consumption reallocation plots completed ===\n');
fprintf('Generated:\n');
fprintf('  1. home_consumption_reallocation.png (CHg/CHs by sector)\n');
fprintf('  2. foreign_consumption_reallocation.png (CFg/CFs by sector)\n\n');
