% analyze_services_demand.m
% Diagnostic script to understand services consumption dynamics

fprintf('\n=== ANALYZING SERVICES DEMAND ===\n\n');

% Display sector classification and weights
fprintf('Sector Classification and Consumption Weights:\n');
fprintf('%-4s %-40s %-10s %-10s %-10s\n', 'Idx', 'Name', 'Goods γ^g', 'Services γ^s', 'Type');
fprintf('%s\n', repmat('-', 80, 1));

for i = 1:nsec
    sector_type = 'Mixed';
    if goods(i)
        sector_type = 'GOODS';
    elseif services(i)
        sector_type = 'SERVICES';
    end
    
    fprintf('%-4d %-40s %8.4f   %8.4f   %-10s\n', ...
        i, names{i}(1:min(40,length(names{i}))), gammag(i), gammas(i), sector_type);
end

fprintf('\n');
fprintf('Sum of goods weights: %.4f\n', sum(gammag));
fprintf('Sum of services weights: %.4f\n', sum(gammas));

% Find main service providers
[sorted_gammas, service_idx] = sort(gammas, 'descend');
fprintf('\nTop 5 Service Providers (by consumption weight γ^s):\n');
for i = 1:min(5, nsec)
    idx = service_idx(i);
    fprintf('%d. Sector %d (%s): γ^s = %.4f\n', ...
        i, idx, names{idx}(1:min(40,length(names{idx}))), gammas(idx));
end

% Check if IRFs exist
if exist('oo_', 'var') && isfield(oo_, 'irfs')
    fprintf('\n=== CHECKING OUTPUT RESPONSES FOR SERVICE SECTORS ===\n');
    shock_name = 'epsA_3';
    
    for i = 1:min(5, nsec)
        idx = service_idx(i);
        Y_irf_name = sprintf('Y_%d_%s', idx, shock_name);
        CHs_irf_name = sprintf('CHs_%d_%s', idx, shock_name);
        
        if isfield(oo_.irfs, Y_irf_name)
            Y_impact = oo_.irfs.(Y_irf_name)(1) * 100;
            Y_peak = max(abs(oo_.irfs.(Y_irf_name))) * 100;
            
            if isfield(oo_.irfs, CHs_irf_name)
                CHs_impact = oo_.irfs.(CHs_irf_name)(1) * 100;
                CHs_peak = max(abs(oo_.irfs.(CHs_irf_name))) * 100;
                
                fprintf('\nSector %d (%s):\n', idx, names{idx}(1:min(30,length(names{idx}))));
                fprintf('  Output (Y):             Impact = %+.3f%%, Peak = %+.3f%%\n', Y_impact, Y_peak);
                fprintf('  Home Cons (CHs):        Impact = %+.3f%%, Peak = %+.3f%%\n', CHs_impact, CHs_peak);
            end
        end
    end
    
    % Check aggregate services consumption
    if isfield(oo_.irfs, ['C_s_' shock_name])
        Cs_irf = oo_.irfs.(['C_s_' shock_name]);
        fprintf('\n=== AGGREGATE SERVICES CONSUMPTION ===\n');
        fprintf('C_s impact: %+.3f%%\n', Cs_irf(1)*100);
        fprintf('C_s peak: %+.3f%% at quarter %d\n', max(abs(Cs_irf))*100, find(abs(Cs_irf)==max(abs(Cs_irf)),1)-1);
    end
    
    % Check which sectors see increased home services consumption
    fprintf('\n=== ALL SECTORS: Home Services Consumption Changes (CHs) ===\n');
    fprintf('%-4s %-40s %-12s %-12s\n', 'Idx', 'Name', 'Impact (%)', 'Peak (%)');
    fprintf('%s\n', repmat('-', 70, 1));
    
    for i = 1:nsec
        CHs_irf_name = sprintf('CHs_%d_%s', i, shock_name);
        if isfield(oo_.irfs, CHs_irf_name)
            CHs_data = oo_.irfs.(CHs_irf_name);
            impact = CHs_data(1) * 100;
            peak = max(abs(CHs_data)) * 100;
            
            fprintf('%-4d %-40s %+10.3f   %+10.3f\n', ...
                i, names{i}(1:min(40,length(names{i}))), impact, peak);
        end
    end
else
    fprintf('\nNo IRF data available. Run main_SOE.m with EXERCISE=2 first.\n');
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');
