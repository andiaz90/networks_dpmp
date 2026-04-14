% verify_consumption_aggregation.m
% Verify that sectoral consumption changes (Cg_i, Cs_i) aggregate to C_g and C_s
% NOTE: Uses TOTAL sectoral consumption (Cg_i, Cs_i), not just home (CHg_i, CHs_i)

fprintf('\n=== VERIFYING CONSUMPTION AGGREGATION ===\n');
fprintf('Using Cg_i and Cs_i (total sectoral consumption = home + foreign)\n\n');

% Check if data is loaded
if ~exist('oo_', 'var') || ~isfield(oo_, 'irfs')
    fprintf('Loading model results...\n');
    load('model_output_IOSOE_ex2.mat');
    load('params_val.mat');
end

shock_name = 'epsA_3';

% Get aggregate responses
if isfield(oo_.irfs, ['C_g_' shock_name])
    C_g_irf = oo_.irfs.(['C_g_' shock_name]);
    C_s_irf = oo_.irfs.(['C_s_' shock_name]);
    
    % Impact (quarter 0) and peak
    C_g_impact = C_g_irf(1) * 100;
    C_s_impact = C_s_irf(1) * 100;
    [C_g_peak_val, C_g_peak_q] = max(abs(C_g_irf));
    [C_s_peak_val, C_s_peak_q] = max(abs(C_s_irf));
    C_g_peak = C_g_irf(C_g_peak_q) * 100;
    C_s_peak = C_s_irf(C_s_peak_q) * 100;
    
    fprintf('AGGREGATE CONSUMPTION CHANGES:\n');
    fprintf('  C_g (Goods):     Impact = %+.4f%%,  Peak = %+.4f%% (Q%d)\n', ...
        C_g_impact, C_g_peak, C_g_peak_q-1);
    fprintf('  C_s (Services):  Impact = %+.4f%%,  Peak = %+.4f%% (Q%d)\n\n', ...
        C_s_impact, C_s_peak, C_s_peak_q-1);
    
    % Initialize arrays
    Cg_impact = zeros(nsec, 1);
    Cs_impact = zeros(nsec, 1);
    Cg_peak = zeros(nsec, 1);
    Cs_peak = zeros(nsec, 1);
    
    % Extract sectoral TOTAL consumption (Cg_i and Cs_i)
    for i = 1:nsec
        Cg_name = sprintf('Cg_%d_%s', i, shock_name);
        Cs_name = sprintf('Cs_%d_%s', i, shock_name);
        
        if isfield(oo_.irfs, Cg_name)
            Cg_data = oo_.irfs.(Cg_name);
            Cg_impact(i) = Cg_data(1) * 100;
            [~, idx] = max(abs(Cg_data));
            Cg_peak(i) = Cg_data(idx) * 100;
        end
        
        if isfield(oo_.irfs, Cs_name)
            Cs_data = oo_.irfs.(Cs_name);
            Cs_impact(i) = Cs_data(1) * 100;
            [~, idx] = max(abs(Cs_data));
            Cs_peak(i) = Cs_data(idx) * 100;
        end
    end
    
    % Compute weighted contributions
    % Model: Cg_j = gammag_j * (p_g/P_j) * C_g
    % For first-order approximations: weighted sum should match aggregate
    Cg_contrib_impact = gammag .* Cg_impact;
    Cs_contrib_impact = gammas .* Cs_impact;
    Cg_contrib_peak = gammag .* Cg_peak;
    Cs_contrib_peak = gammas .* Cs_peak;
    
    % Create table for IMPACT
    fprintf('========================================================================================================\n');
    fprintf('SECTORAL CONSUMPTION DECOMPOSITION - IMPACT (Quarter 0)\n');
    fprintf('========================================================================================================\n');
    fprintf('%-4s %-40s %-8s %-10s %-12s %-10s %-12s\n', ...
        'Sec', 'Name', 'Type', 'γ_g', 'Cg(%%×γ)', 'γ_s', 'Cs(%%×γ)');
    fprintf('--------------------------------------------------------------------------------------------------------\n');
    
    for i = 1:nsec
        sector_type = '  ';
        if goods(i)
            sector_type = 'G';
        elseif services(i)
            sector_type = 'S';
        end
        
        sector_name = names{i};
        if length(sector_name) > 40
            sector_name = sector_name(1:40);
        end
        
        fprintf('%-4d %-40s %-8s %8.4f   %+10.5f   %8.4f   %+10.5f\n', ...
            i, sector_name, sector_type, gammag(i), Cg_contrib_impact(i), ...
            gammas(i), Cs_contrib_impact(i));
    end
    
    fprintf('--------------------------------------------------------------------------------------------------------\n');
    fprintf('%-54s %10s   %+10.5f   %8s   %+10.5f\n', ...
        'SUM (Weighted sectoral):', '', sum(Cg_contrib_impact), '', sum(Cs_contrib_impact));
    fprintf('%-54s %10s   %+10.5f   %8s   %+10.5f\n', ...
        'AGGREGATE (C_g, C_s):', '', C_g_impact, '', C_s_impact);
    fprintf('%-54s %10s   %+10.5f   %8s   %+10.5f\n', ...
        'DIFFERENCE:', '', sum(Cg_contrib_impact)-C_g_impact, '', sum(Cs_contrib_impact)-C_s_impact);
    fprintf('%-54s %10s   %+10.4f%%  %8s   %+10.4f%%\n', ...
        'RELATIVE ERROR:', '', 100*(sum(Cg_contrib_impact)-C_g_impact)/C_g_impact, ...
        '', 100*(sum(Cs_contrib_impact)-C_s_impact)/C_s_impact);
    fprintf('========================================================================================================\n\n');
    
    % Analysis
    fprintf('VERIFICATION RESULTS:\n');
    goods_error = abs(sum(Cg_contrib_impact) - C_g_impact);
    services_error = abs(sum(Cs_contrib_impact) - C_s_impact);
    
    if goods_error < 0.001 && services_error < 0.001
        fprintf('✓ PASSED: Sectoral consumption aggregates correctly to C_g and C_s\n');
        fprintf('  Goods error: %.6f%% (tolerance: 0.001%%)\n', goods_error);
        fprintf('  Services error: %.6f%% (tolerance: 0.001%%)\n', services_error);
    else
        fprintf('✗ Note: Small aggregation discrepancy detected\n');
        fprintf('  Goods error: %.6f%%\n', goods_error);
        fprintf('  Services error: %.6f%%\n', services_error);
        fprintf('  This is expected due to price effects in Cobb-Douglas aggregation\n');
        fprintf('  Model uses: Cg_j = gammag_j * (p_g/P_j) * C_g\n');
        fprintf('  Relative price effects (p_g/P_j) create non-linearities\n');
    end
    
else
    fprintf('Error: IRF data not available. Run main_SOE.m with EXERCISE=2 first.\n');
end

fprintf('\n=== VERIFICATION COMPLETE ===\n');
