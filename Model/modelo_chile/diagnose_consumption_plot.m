% diagnose_consumption_plot.m
% Check which consumption variables are actually being plotted

fprintf('\n=== DIAGNOSING CONSUMPTION PLOT VARIABLES ===\n\n');

if ~exist('oo_', 'var') || ~isfield(oo_, 'irfs')
    load('model_output_IOSOE_ex2.mat');
    load('params_val.mat');
end

shock_name = 'epsA_3';

fprintf('%-4s %-40s %-12s %-15s %-15s\n', 'Sec', 'Name', 'Type', 'CHg Impact', 'CHs Impact');
fprintf('%s\n', repmat('-', 90, 1));

for i = 1:nsec
    CHg_name = sprintf('CHg_%d_%s', i, shock_name);
    CHs_name = sprintf('CHs_%d_%s', i, shock_name);
    
    CHg_val = 0;
    CHs_val = 0;
    
    if isfield(oo_.irfs, CHg_name)
        CHg_val = oo_.irfs.(CHg_name)(1) * 100;
    end
    
    if isfield(oo_.irfs, CHs_name)
        CHs_val = oo_.irfs.(CHs_name)(1) * 100;
    end
    
    sector_type = 'Mixed';
    if goods(i)
        sector_type = 'GOODS';
    elseif services(i)
        sector_type = 'SERVICES';
    end
    
    % Determine which would be plotted first
    plotted_var = 'CHg';
    plotted_val = CHg_val;
    if abs(CHs_val) > abs(CHg_val)
        plotted_var = 'CHs';
        plotted_val = CHs_val;
    end
    
    fprintf('%-4d %-40s %-12s %+13.6f   %+13.6f', ...
        i, names{i}(1:min(40,length(names{i}))), sector_type, CHg_val, CHs_val);
    
    % Highlight what would actually be plotted
    if strcmp(sector_type, 'SERVICES') && abs(CHg_val) < 1e-10
        fprintf('   << PLOTTING SHOWS ~ZERO (CHg)\n');
    elseif strcmp(sector_type, 'SERVICES')
        fprintf('   << Service sector\n');
    else
        fprintf('\n');
    end
end

fprintf('\n');
fprintf('KEY INSIGHT:\n');
fprintf('- plot_consumption_reallocation.m tries CHg first, then CHs\n');
fprintf('- For SERVICE sectors, CHg ≈ 0 (they produce services, not goods)\n');
fprintf('- So the plot shows ~0 for services when it should show CHs!\n');
fprintf('- Solution: Plot should check sector type and use appropriate variable\n\n');

fprintf('COMPARISON: What aggregate shows vs what sectors show\n');
if isfield(oo_.irfs, ['C_g_' shock_name])
    C_g_impact = oo_.irfs.(['C_g_' shock_name])(1) * 100;
    C_s_impact = oo_.irfs.(['C_s_' shock_name])(1) * 100;
    
    fprintf('Aggregate C_g impact: %+.4f%%\n', C_g_impact);
    fprintf('Aggregate C_s impact: %+.4f%%\n', C_s_impact);
    
    fprintf('\nSum of CHg across all sectors: %+.4f%%\n', sum(arrayfun(@(i) ...
        oo_.irfs.(sprintf('CHg_%d_%s', i, shock_name))(1)*100, 1:nsec)));
    fprintf('Sum of CHs across all sectors: %+.4f%%\n', sum(arrayfun(@(i) ...
        oo_.irfs.(sprintf('CHs_%d_%s', i, shock_name))(1)*100, 1:nsec)));
end

fprintf('\n=== DIAGNOSIS COMPLETE ===\n');
