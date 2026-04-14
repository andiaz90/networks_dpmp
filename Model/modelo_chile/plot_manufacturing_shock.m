% plot_manufacturing_shock.m
% Plot IRFs for Manufacturing TFP Shock (epsA_3)

% Define key variables to plot
key_vars = {
    'C',        'Consumption';
    'Y',        'Output';
    'pi',       'Inflation';
    'r',        'Interest Rate';
    'N',        'Labor';
    'w',        'Wage';
    'TB',       'Trade Balance';
    'Q',        'Exchange Rate';
    'X',        'Exports';
    'Y_3',      'Manufacturing Output';
    'P_3',      'Manufacturing Price';
    'L_3',      'Manufacturing Labor'
};

% Manufacturing shock name
shock_name = 'epsA_3';
shock_label = 'Manufacturing TFP Shock';

% Check which IRFs are available
irf_fields = fieldnames(oo_.irfs);
available_vars = {};
available_labels = {};

% Check which variables have IRFs for this shock
for i = 1:size(key_vars, 1)
    var_name = key_vars{i, 1};
    irf_name = [var_name '_' shock_name];
    
    if ismember(irf_name, irf_fields)
        available_vars{end+1} = var_name;
        available_labels{end+1} = key_vars{i, 2};
    else
        fprintf('Warning: IRF for %s not available for manufacturing shock.\n', var_name);
    end
end

if isempty(available_vars)
    error('No IRFs found for manufacturing TFP shock (epsA_3). Check model solution.');
end

% Create figure
n_vars = length(available_vars);
n_rows = ceil(n_vars / 3);
n_cols = min(3, n_vars);

fig = figure('Name', 'Manufacturing TFP Shock - Impulse Response Functions', ...
             'Position', [100 100 1400 900], ...
             'Visible', 'on', ...
             'WindowStyle', 'docked');

% Bring figure to front
figure(fig);
drawnow;

% Time horizon (quarters)
T_full = length(oo_.irfs.([available_vars{1} '_' shock_name]));
T = min(40, T_full);  % Show up to 40 periods for persistent shock
quarters = 0:T-1;

% Plot each variable
for i = 1:n_vars
    subplot(n_rows, n_cols, i);
    
    var_name = available_vars{i};
    irf_name = [var_name '_' shock_name];
    irf_data = oo_.irfs.(irf_name)(1:T);
    
    % Plot IRF
    plot(quarters, 100*irf_data, 'LineWidth', 2.5, 'Color', [0.8 0.2 0.2]);
    hold on;
    
    % Add zero line
    plot(quarters, zeros(size(quarters)), 'k--', 'LineWidth', 0.5);
    
    % Formatting
    title(available_labels{i}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Quarters');
    if contains(var_name, {'pi', 'r'})
        ylabel('Percentage Points');
    else
        ylabel('Percent Deviation from SS');
    end
    
    grid on;
    xlim([0 T-1]);
    
    % Add some visual improvements
    max_abs_val = max(abs(irf_data));
    if max_abs_val > 0
        ylim_range = max_abs_val * 110;  % 110% of max for padding
        ylim([-ylim_range ylim_range]);
    end
    enforce_min_yaxis(gca);
end

% Add main title
sgtitle('Impulse Response Functions: Manufacturing TFP Shock',...
        'FontSize', 16, 'FontWeight', 'bold');

% Add shock information
if exist('manufacturing_sector', 'var')
    sector_info = sprintf('Sector %d: Manufacturing (1%% positive TFP shock)', manufacturing_sector);
else
    sector_info = 'Sector 3: Manufacturing (1% positive TFP shock)';
end

annotation('textbox', [0.70 0.02 0.28 0.07], ...
           'String', sprintf(['%s\n' ...
                              'Persistence: \\rho = %.2f\n' ...
                              'Std Dev: \\sigma = %.2f%%'], ...
                              sector_info, rho_tfp1_val, 1.0), ...
           'Units', 'normalized', ...
           'FontSize', 9, ...
           'EdgeColor', 'black', ...
           'BackgroundColor', 'white', ...
           'HorizontalAlignment', 'left', ...
           'VerticalAlignment', 'bottom', ...
           'FitBoxToText', 'on');

% Force display update
drawnow;

% Save figure
try
    print(fig, '-dpng', '-r300', 'manufacturing_tfp_shock_irfs.png');
    fprintf('\n✓ Manufacturing TFP shock IRFs displayed and saved as manufacturing_tfp_shock_irfs.png\n');
catch
    fprintf('\n✓ Manufacturing TFP shock IRFs displayed (saving failed but plot is visible)\n');
end

% =========================================================================
%% OUTPUT GAP FIGURE 1 – Aggregate gaps: Ygap, GDPgap, Ngap
% =========================================================================
agg_gap_vars   = {'Ygap',   'GDPgap',   'Ngap'};
agg_gap_labels = {'Output Gap  log(Y/Y^{f})', ...
                  'GDP Gap  log(GDP/GDP^{f})', ...
                  'Employment Gap  log(N/N^{f})'};
clr_agg = [0.13 0.47 0.71; 0.84 0.15 0.16; 0.17 0.63 0.17];

fig_agg = figure('Name','Output Gaps – Aggregate (Mfg TFP Shock)', ...
                 'Position',[120 120 1100 320], 'WindowStyle','docked');
for k = 1:3
    vname    = agg_gap_vars{k};
    irf_name = [vname '_' shock_name];
    if isfield(oo_.irfs, irf_name)
        irf_data = 100 * oo_.irfs.(irf_name)(1:T);
    else
        warning('IRF not found: %s. Plotting zeros.', irf_name);
        irf_data = zeros(1, T);
    end
    subplot(1, 3, k);
    plot(quarters, irf_data, 'LineWidth', 2.5, 'Color', clr_agg(k,:));
    hold on;
    plot(quarters, zeros(size(quarters)), 'k--', 'LineWidth', 0.8);
    hold off;
    title(agg_gap_labels{k}, 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');
    xlabel('Quarters', 'FontSize', 9);
    ylabel('% deviation', 'FontSize', 9);
    grid on; box on;
    xlim([0 T-1]);
end
sgtitle(sprintf('Aggregate Output Gaps — %s', shock_label),...
        'FontSize', 13, 'FontWeight', 'bold');
drawnow;

% =========================================================================
%% OUTPUT GAP FIGURE 1b – Goods & Services: output gap vs consumption gap
% =========================================================================
irf_Ygap_g    = get_irf_or_zero(oo_.irfs, 'Ygap_g', shock_name, T);
irf_Ygap_s    = get_irf_or_zero(oo_.irfs, 'Ygap_s', shock_name, T);
irf_Cgap_g    = get_irf_or_zero(oo_.irfs, 'Cgap_g', shock_name, T);
irf_Cgap_s    = get_irf_or_zero(oo_.irfs, 'Cgap_s', shock_name, T);
irf_GDPgap_gs = get_irf_or_zero(oo_.irfs, 'GDPgap', shock_name, T);

clr_Y = [0.84 0.15 0.16];   % output gap  – red
clr_C = [0.13 0.47 0.71];   % consumption gap – blue
clr_GDP = [0.5 0.5 0.5];    % GDP gap – grey

fig_gs = figure('Name','Goods & Services Gaps (Mfg TFP Shock)', ...
                'Position',[130 130 1200 620], 'WindowStyle','docked');

% ---- Column 1: GOODS ----
% (1,1) Goods output gap
subplot(2,3,1);
hold on;
plot(quarters, 100*irf_Ygap_g, 'LineWidth', 2.5, 'Color', clr_Y, ...
     'DisplayName', 'Output gap  log(Y_g/Y_g^f)');
plot(quarters, 100*irf_Cgap_g, 'LineWidth', 2.5, 'Color', clr_C, ...
     'LineStyle', '--', 'DisplayName', 'Consumption gap  log(Ctotg/Ctotg^f)');
plot(quarters, zeros(size(quarters)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold off;
legend('Location', 'best', 'FontSize', 8);
title('Goods: Output gap vs Consumption gap', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Quarters', 'FontSize', 9); ylabel('% deviation from flex-price', 'FontSize', 9);
grid on; box on; xlim([0 T-1]);
enforce_min_yaxis(gca);

% (2,1) Goods output gap alone
subplot(2,3,4);
hold on;
h1 = area(quarters, 100*irf_Ygap_g, 'FaceColor', clr_Y, 'FaceAlpha', 0.25, ...
          'EdgeColor', clr_Y, 'LineWidth', 2.0, 'DisplayName', 'Output gap');
h2 = plot(quarters, 100*irf_Cgap_g, 'LineWidth', 2.5, 'Color', clr_C, ...
          'LineStyle', '--', 'DisplayName', 'Consumption gap');
plot(quarters, zeros(size(quarters)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold off;
legend([h1 h2], 'Location', 'best', 'FontSize', 8);
title('Goods: Gap detail', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Quarters', 'FontSize', 9); ylabel('% deviation', 'FontSize', 9);
grid on; box on; xlim([0 T-1]);
enforce_min_yaxis(gca);

% ---- Column 2: SERVICES ----
% (1,2) Services output gap vs consumption gap
subplot(2,3,2);
hold on;
plot(quarters, 100*irf_Ygap_s, 'LineWidth', 2.5, 'Color', clr_Y, ...
     'DisplayName', 'Output gap  log(Y_s/Y_s^f)');
plot(quarters, 100*irf_Cgap_s, 'LineWidth', 2.5, 'Color', clr_C, ...
     'LineStyle', '--', 'DisplayName', 'Consumption gap  log(Ctots/Ctots^f)');
plot(quarters, zeros(size(quarters)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold off;
legend('Location', 'best', 'FontSize', 8);
title('Services: Output gap vs Consumption gap', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Quarters', 'FontSize', 9); ylabel('% deviation from flex-price', 'FontSize', 9);
grid on; box on; xlim([0 T-1]);
enforce_min_yaxis(gca);

% (2,2) Services gap detail
subplot(2,3,5);
hold on;
h3 = area(quarters, 100*irf_Ygap_s, 'FaceColor', clr_Y, 'FaceAlpha', 0.25, ...
          'EdgeColor', clr_Y, 'LineWidth', 2.0, 'DisplayName', 'Output gap');
h4 = plot(quarters, 100*irf_Cgap_s, 'LineWidth', 2.5, 'Color', clr_C, ...
          'LineStyle', '--', 'DisplayName', 'Consumption gap');
plot(quarters, zeros(size(quarters)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold off;
legend([h3 h4], 'Location', 'best', 'FontSize', 8);
title('Services: Gap detail', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Quarters', 'FontSize', 9); ylabel('% deviation', 'FontSize', 9);
grid on; box on; xlim([0 T-1]);
enforce_min_yaxis(gca);

% ---- Column 3: CROSS-SECTOR OVERLAY ----
% (1,3) Output gaps: goods vs services vs GDP
subplot(2,3,3);
hold on;
plot(quarters, 100*irf_Ygap_g,    'LineWidth', 2.5, 'Color', [0.84 0.15 0.16], ...
     'DisplayName', 'Goods Ygap');
plot(quarters, 100*irf_Ygap_s,    'LineWidth', 2.5, 'Color', [0.13 0.47 0.71], ...
     'DisplayName', 'Services Ygap');
plot(quarters, 100*irf_GDPgap_gs, 'LineWidth', 1.8, 'Color', clr_GDP, ...
     'LineStyle', '--', 'DisplayName', 'GDP gap');
plot(quarters, zeros(size(quarters)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold off;
legend('Location', 'best', 'FontSize', 9);
title('Output Gaps: Goods vs Services', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Quarters', 'FontSize', 9); ylabel('% deviation', 'FontSize', 9);
grid on; box on; xlim([0 T-1]);
enforce_min_yaxis(gca);

% (2,3) Consumption gaps: goods vs services vs GDP
subplot(2,3,6);
hold on;
plot(quarters, 100*irf_Cgap_g,    'LineWidth', 2.5, 'Color', [0.84 0.15 0.16], ...
     'LineStyle', '--', 'DisplayName', 'Goods Cgap');
plot(quarters, 100*irf_Cgap_s,    'LineWidth', 2.5, 'Color', [0.13 0.47 0.71], ...
     'LineStyle', '--', 'DisplayName', 'Services Cgap');
plot(quarters, 100*irf_GDPgap_gs, 'LineWidth', 1.8, 'Color', clr_GDP, ...
     'LineStyle', ':', 'DisplayName', 'GDP gap');
plot(quarters, zeros(size(quarters)), 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold off;
legend('Location', 'best', 'FontSize', 9);
title('Consumption Gaps: Goods vs Services', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Quarters', 'FontSize', 9); ylabel('% deviation', 'FontSize', 9);
grid on; box on; xlim([0 T-1]);
enforce_min_yaxis(gca);

sgtitle(sprintf('Goods & Services — Output Gap vs Consumption Gap  |  %s', shock_label),...
        'FontSize', 13, 'FontWeight', 'bold');
drawnow;

% =========================================================================
%% OUTPUT GAP FIGURE 2 – Sectoral output gaps Ygap_i  (4 × 3 grid)
% =========================================================================
nsec_local = 12;

% Sector labels (use workspace 'names' if available, else generic)
if exist('names', 'var') && numel(names) >= nsec_local
    sec_lbl = cellfun(@(s) s(1:min(15,numel(s))), names(1:nsec_local), ...
                      'UniformOutput', false);
else
    sec_lbl = arrayfun(@(i) sprintf('Sector %d',i), 1:nsec_local, ...
                       'UniformOutput', false);
end

clr_sec = lines(nsec_local);

fig_sec = figure('Name','Output Gaps – Sectoral (Mfg TFP Shock)', ...
                 'Position',[140 140 1300 750], 'WindowStyle','docked');

for i = 1:nsec_local
    vname    = sprintf('Ygap_%d', i);
    irf_name = [vname '_' shock_name];
    if isfield(oo_.irfs, irf_name)
        irf_data = 100 * oo_.irfs.(irf_name)(1:T);
    else
        irf_data = zeros(1, T);
    end

    subplot(3, 4, i);
    plot(quarters, irf_data, 'LineWidth', 2.2, 'Color', clr_sec(i,:));
    hold on;
    plot(quarters, zeros(size(quarters)), 'k--', 'LineWidth', 0.8);
    % Highlight manufacturing sector
    if i == 3
        patch([quarters fliplr(quarters)], ...
              [irf_data fliplr(zeros(size(quarters)))], ...
              clr_sec(i,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    hold off;
    title(sec_lbl{i}, 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlabel('Qtrs', 'FontSize', 7);
    ylabel('%', 'FontSize', 7);
    grid on; box on;
    xlim([0 T-1]);
    enforce_min_yaxis(gca);
    set(gca, 'FontSize', 7);
end

sgtitle(sprintf('Sectoral Output Gaps  Ygap\\_i = log(Y\\_i / Y^f\\_i) — %s', shock_label),...
        'FontSize', 12, 'FontWeight', 'bold');
drawnow;

% =========================================================================
%% OUTPUT GAP FIGURE 3 – Overlay: manufacturing vs aggregate Ygap
% =========================================================================
fig_cmp = figure('Name','Output Gap – Mfg vs Aggregate', ...
                 'Position',[160 160 700 350], 'WindowStyle','docked');

irf_Ygap   = get_irf_or_zero(oo_.irfs, 'Ygap',   shock_name, T);
irf_GDPgap = get_irf_or_zero(oo_.irfs, 'GDPgap', shock_name, T);
irf_Ysec3  = get_irf_or_zero(oo_.irfs, 'Ygap_3', shock_name, T);

subplot(1,2,1);
hold on;
plot(quarters, 100*irf_Ysec3,  'LineWidth', 2.5, 'Color', [0.84 0.15 0.16], ...
     'DisplayName', 'Manufacturing (sec.3)');
plot(quarters, 100*irf_Ygap,   'LineWidth', 2.0, 'Color', [0.13 0.47 0.71], ...
     'DisplayName', 'Aggregate Y gap');
plot(quarters, 100*irf_GDPgap, 'LineWidth', 2.0, 'Color', [0.17 0.63 0.17], ...
     'LineStyle', '--', 'DisplayName', 'GDP gap');
plot(quarters, zeros(size(quarters)), 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
hold off;
legend('Location','best', 'FontSize', 9);
title('Output Gaps: Source vs Aggregate', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Quarters'); ylabel('% deviation');
grid on; box on; xlim([0 T-1]);
enforce_min_yaxis(gca);

subplot(1,2,2);
% Cross-sector dispersion at each horizon
irf_all_sec = zeros(nsec_local, T);
for i = 1:nsec_local
    irf_all_sec(i,:) = 100 * get_irf_or_zero(oo_.irfs, sprintf('Ygap_%d',i), shock_name, T);
end
p10  = prctile(irf_all_sec, 10, 1);
p25  = prctile(irf_all_sec, 25, 1);
p50  = prctile(irf_all_sec, 50, 1);
p75  = prctile(irf_all_sec, 75, 1);
p90  = prctile(irf_all_sec, 90, 1);

hold on;
fill([quarters fliplr(quarters)],[p10 fliplr(p90)], [0.7 0.85 1.0], ...
     'EdgeColor','none','FaceAlpha',0.5,'DisplayName','p10–p90');
fill([quarters fliplr(quarters)],[p25 fliplr(p75)], [0.4 0.65 0.9], ...
     'EdgeColor','none','FaceAlpha',0.5,'DisplayName','p25–p75');
plot(quarters, p50,             'b-',  'LineWidth',1.5, 'DisplayName','Median sector');
plot(quarters, 100*irf_Ygap,    'r--', 'LineWidth',2.0, 'DisplayName','Aggregate Ygap');
plot(quarters, zeros(size(quarters)), 'k:', 'LineWidth',0.8, 'HandleVisibility','off');
hold off;
legend('Location','best','FontSize',8);
title('Cross-Sector Dispersion', 'FontSize',10, 'FontWeight','bold');
xlabel('Quarters'); ylabel('%');
grid on; box on; xlim([0 T-1]);
enforce_min_yaxis(gca);

sgtitle(sprintf('Output Gap Analysis — %s', shock_label),...
        'FontSize', 12, 'FontWeight', 'bold');
drawnow;

try
    print(fig_agg, '-dpng', '-r200', 'gap_aggregate_mfg_shock.png');
    print(fig_gs,  '-dpng', '-r200', 'gap_goods_services_mfg_shock.png');
    print(fig_sec, '-dpng', '-r200', 'gap_sectoral_mfg_shock.png');
    print(fig_cmp, '-dpng', '-r200', 'gap_comparison_mfg_shock.png');
    fprintf('\n✓ Output gap figures saved.\n');
catch
    fprintf('\n✓ Output gap figures plotted (saving failed but plots are visible).\n');
end

% =========================================================================
%% Display impact statistics
% =========================================================================
fprintf('\n=== MANUFACTURING TFP SHOCK IMPACT SUMMARY ===\n');
fprintf('Time horizon: %d quarters\n', T);
fprintf('Shock persistence (rho): %.2f\n\n', rho_tfp1_val);

for i = 1:n_vars
    var_name = available_vars{i};
    irf_name = [var_name '_' shock_name];
    irf_data = oo_.irfs.(irf_name);
    
    % Peak impact
    [max_val, max_idx] = max(abs(irf_data));
    peak_impact = 100 * irf_data(max_idx);
    peak_quarter = max_idx - 1;
    
    % Impact at quarter 4 and 12
    impact_q4 = 100 * irf_data(min(5, length(irf_data)));
    impact_q12 = 100 * irf_data(min(13, length(irf_data)));
    
    fprintf('%-25s: Peak = %7.3f%% (Q%d) | Q4 = %7.3f%% | Q12 = %7.3f%%\n', ...
            available_labels{i}, peak_impact, peak_quarter, impact_q4, impact_q12);
end

% Gap variable summary
fprintf('\n--- Output Gap Impact ---\n');
gap_vars   = [{'Ygap','GDPgap','Ngap','Ygap_g','Ygap_s','Cgap_g','Cgap_s'}, ...
               arrayfun(@(i) sprintf('Ygap_%d',i), 1:nsec_local, 'UniformOutput',false)];
gap_labels = [{'Aggregate Y gap','GDP gap','Employment gap', ...
               'Goods Output gap','Services Output gap', ...
               'Goods Consumption gap','Services Consumption gap'}, sec_lbl'];
for k = 1:numel(gap_vars)
    irf_name = [gap_vars{k} '_' shock_name];
    if isfield(oo_.irfs, irf_name)
        d = oo_.irfs.(irf_name);
        [~, mi] = max(abs(d));
        fprintf('  %-28s: Peak = %+7.3f%% (Q%d) | Q4 = %+7.3f%%\n', ...
            gap_labels{k}, 100*d(mi), mi-1, 100*d(min(5,end)));
    end
end
fprintf('\n');

% =========================================================================
%% LOCAL HELPER
% =========================================================================
function v = get_irf_or_zero(irfs, vname, shock, T)
    fname = [vname '_' shock];
    if isfield(irfs, fname)
        raw = irfs.(fname)(:);
    else
        raw = zeros(T,1);
    end
    if numel(raw) < T; raw(end+1:T) = 0; end
    v = raw(1:T)';
end
