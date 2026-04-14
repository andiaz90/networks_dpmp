% plot_sectoral_reallocation.m
% Plots sectoral labor and consumption reallocation across different shocks
% Similar to Figure 3 but showing individual sector responses

%% Setup
if ~exist('nsec', 'var')
    nsec = 12;
end

if ~exist('names', 'var')
    warning('Sector names not loaded. Using generic names.');
    names = arrayfun(@(i) sprintf('Sector %d', i), 1:nsec, 'UniformOutput', false);
end

% Shorten sector names for plotting
short_names = cell(nsec, 1);
for i = 1:nsec
    if length(names{i}) > 25
        short_names{i} = names{i}(1:22);
        short_names{i} = strrep(short_names{i}, ' y ', ' ');
    else
        short_names{i} = names{i};
    end
end

% Define time horizon for plots
T_plot = 40;  % quarters to plot

%% Check what data is available
fprintf('\n=== Checking available data for plotting ===\n');

% Check if IRF data exists (preferred for impulse responses)
use_irf_data = false;
if exist('oo_', 'var') && isfield(oo_, 'irfs')
    irf_fields = fieldnames(oo_.irfs);
    fprintf('Found %d IRF series in oo_.irfs\n', length(irf_fields));
    if ~isempty(irf_fields)
        fprintf('Sample IRF fields: %s, %s\n', irf_fields{1}, irf_fields{min(2,end)});
        use_irf_data = true;
    end
else
    fprintf('No IRF data found in oo_.irfs\n');
end

% Check if cell array data exists (from extraction in main_SOE)
use_cell_data = false;
if exist('L', 'var') && iscell(L) && ~isempty(L{1})
    fprintf('Found L cell array with %d sectors\n', length(L));
    fprintf('L{1} has %d time periods\n', length(L{1}));
    use_cell_data = true;
else
    fprintf('No L cell array found\n');
end

if use_irf_data
    fprintf('Using IRF data from oo_.irfs\n');
elseif use_cell_data
    fprintf('Using extracted cell array data\n');
else
    warning('No suitable data found for plotting. Plots may be empty.');
end
fprintf('=========================================\n\n');

% Colors for different shock types
color_tfp = [0.0, 0.4470, 0.7410];      % Blue
color_monetary = [0.8500, 0.3250, 0.0980];  % Red-orange
color_preference = [0.4660, 0.6740, 0.1880];  % Green

%% ========================================
%% CONSOLIDATED SECTORAL REALLOCATION FIGURE
%% ========================================

fig = figure('Name', 'Sectoral Reallocation: Labor IRFs', 'Position', [50, 50, 1600, 900], 'Visible', 'on');
drawnow;

% Check if IRF data exists and extract it
if exist('oo_', 'var') && isfield(oo_, 'irfs')
    fprintf('Using IRF data for plotting...\n');
    
    % Get IRF length
    irf_names = fieldnames(oo_.irfs);
    if ~isempty(irf_names)
        irf_length = length(oo_.irfs.(irf_names{1}));
        T_plot = min(T_plot, irf_length);
        fprintf('IRF length: %d periods\n', irf_length);
    end
    
    % Panel 1-12: Labor reallocation by sector
    for i = 1:nsec
        subplot(3, 5, i);
        hold on; grid on;
        
        % Try different shock names (TFP shock to sector 3 - manufacturing)
        shock_names = {sprintf('epsA_%d', 3), sprintf('epsA_%d', i), 'eps_i', 'eps_om'};
        var_name = sprintf('L_%d', i);
        
        plotted = false;
        for shock_idx = 1:length(shock_names)
            shock = shock_names{shock_idx};
            irf_field = sprintf('%s_%s', var_name, shock);
            
            if isfield(oo_.irfs, irf_field)
                irf_data = oo_.irfs.(irf_field);
                plot(0:length(irf_data)-1, 100*irf_data, 'LineWidth', 1.5, ...
                    'Color', [0.0, 0.4470, 0.7410]);
                plotted = true;
                break;
            end
        end
        
        if ~plotted
            % Fallback: try to find any IRF for this variable
            matching_fields = irf_names(contains(irf_names, var_name));
            if ~isempty(matching_fields)
                irf_data = oo_.irfs.(matching_fields{1});
                plot(0:length(irf_data)-1, 100*irf_data, 'LineWidth', 1.5, ...
                    'Color', [0.0, 0.4470, 0.7410]);
                plotted = true;
            end
        end
        
        if ~plotted
            text(0.5, 0.5, 'No IRF data', 'Units', 'normalized', 'HorizontalAlignment', 'center');
        end
        
        yline(0, 'k--', 'LineWidth', 0.5);
        title(short_names{i}, 'FontSize', 8, 'FontWeight', 'bold');
        xlabel('Quarters', 'FontSize', 7);
        ylabel('Labor (% dev.)', 'FontSize', 7);
        xlim([0, T_plot-1]);
        enforce_min_yaxis(gca);
        set(gca, 'FontSize', 7);
    end
    
    % Panel 13: Aggregate variables
    subplot(3, 5, 13);
    hold on; grid on;
    
    % Plot aggregate consumption IRF
    if isfield(oo_.irfs, 'C_epsA_3')
        plot(0:length(oo_.irfs.C_epsA_3)-1, 100*oo_.irfs.C_epsA_3, 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880]);
    elseif isfield(oo_.irfs, 'C_eps_i')
        plot(0:length(oo_.irfs.C_eps_i)-1, 100*oo_.irfs.C_eps_i, 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880]);
    end
    
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('Quarters', 'FontSize', 8);
    ylabel('Consumption (% dev.)', 'FontSize', 8);
    title('Aggregate Consumption', 'FontSize', 9, 'FontWeight', 'bold');
    xlim([0, T_plot-1]);
    enforce_min_yaxis(gca);
    set(gca, 'FontSize', 8);
    
    % Panel 14: Total labor
    subplot(3, 5, 14);
    hold on; grid on;
    
    if isfield(oo_.irfs, 'N_epsA_3')
        plot(0:length(oo_.irfs.N_epsA_3)-1, 100*oo_.irfs.N_epsA_3, 'LineWidth', 2, 'Color', [0.2, 0.6, 0.3]);
    elseif isfield(oo_.irfs, 'N_eps_i')
        plot(0:length(oo_.irfs.N_eps_i)-1, 100*oo_.irfs.N_eps_i, 'LineWidth', 2, 'Color', [0.2, 0.6, 0.3]);
    end
    
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('Quarters', 'FontSize', 8);
    ylabel('Total Labor (% dev.)', 'FontSize', 8);
    title('Aggregate Labor', 'FontSize', 9, 'FontWeight', 'bold');
    xlim([0, T_plot-1]);
    enforce_min_yaxis(gca);
    set(gca, 'FontSize', 8);
    
    % Panel 15: GDP
    subplot(3, 5, 15);
    hold on; grid on;
    
    if isfield(oo_.irfs, 'GDP_epsA_3')
        plot(0:length(oo_.irfs.GDP_epsA_3)-1, 100*oo_.irfs.GDP_epsA_3, 'LineWidth', 2, 'Color', [0.8, 0.2, 0.2]);
    elseif isfield(oo_.irfs, 'GDP_eps_i')
        plot(0:length(oo_.irfs.GDP_eps_i)-1, 100*oo_.irfs.GDP_eps_i, 'LineWidth', 2, 'Color', [0.8, 0.2, 0.2]);
    end
    
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('Quarters', 'FontSize', 8);
    ylabel('GDP (% dev.)', 'FontSize', 8);
    title('GDP', 'FontSize', 9, 'FontWeight', 'bold');
    xlim([0, T_plot-1]);
    enforce_min_yaxis(gca);
    set(gca, 'FontSize', 8);
    
    sync_ylims(gcf);
    sgtitle('Sectoral Labor Reallocation: IRFs to Manufacturing TFP Shock', 'FontSize', 13, 'FontWeight', 'bold');
    
else
    % Fallback message
    text(0.5, 0.5, 'No IRF data available. Run stoch_simul with irf option.', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14);
    fprintf('ERROR: No IRF data found. Make sure Dynare runs stoch_simul successfully.\n');
end

% Save to output folder if it exists
if ~exist('output', 'dir')
    mkdir('output');
end
saveas(fig, 'output/sectoral_reallocation_IRFs.png');
% Also save in main directory for easier access
print(fig, '-dpng', '-r300', 'sectoral_reallocation_IRFs.png');
fprintf('✓ Sectoral reallocation figure saved: sectoral_reallocation_IRFs.png\n');

% Print available IRF fields for debugging
if exist('oo_', 'var') && isfield(oo_, 'irfs')
    fprintf('\nAvailable IRFs (first 20):\n');
    irf_list = fieldnames(oo_.irfs);
    for i = 1:min(20, length(irf_list))
        fprintf('  %s\n', irf_list{i});
    end
    if length(irf_list) > 20
        fprintf('  ... and %d more\n', length(irf_list) - 20);
    end
end

fprintf('\n=== Sectoral reallocation IRF plot completed ===\n');
