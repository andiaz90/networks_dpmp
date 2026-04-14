% Generate Production Structure Tables for Paper
clear; clc;

fprintf('=== GENERATING PRODUCTION STRUCTURE TABLES ===\n\n');

%% Read data
alldata = readtable('Stata_to_excel_few_industries_chile.xls');
names = alldata.Nameshort;
alpha = alldata.alpha;        % Material share (domestic intermediate inputs)
alpha_V = alldata.alpha_V;    % Import share
labor_share = 1 - alpha - alpha_V;  % Labor share (residual)

nsec = length(alpha);

%% Table 1: Production Input Shares
fprintf('TABLE 1: PRODUCTION INPUT SHARES BY SECTOR\n');
fprintf('===========================================\n\n');
fprintf('\\begin{table}[H]\n');
fprintf('\\centering\n');
fprintf('\\caption{Production Input Shares by Sector}\n');
fprintf('\\label{tab:production_shares}\n');
fprintf('\\small\n');
fprintf('\\begin{tabular}{clcccc}\n');
fprintf('\\toprule\n');
fprintf('\\textbf{Sec.} & \\textbf{Sector Name} & $\\alpha_{mi}$ & $\\alpha_{vi}$ & Labor & Check \\\\\n');
fprintf('              &                      & (Materials) & (Imports) & Share & Sum \\\\\n');
fprintf('\\midrule\n');

for i = 1:nsec
    check_sum = alpha(i) + alpha_V(i) + labor_share(i);
    fprintf('%d & %s & %.3f & %.3f & %.3f & %.3f \\\\\n', ...
        i, names{i}, alpha(i), alpha_V(i), labor_share(i), check_sum);
end

fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n\n');

fprintf('\\textbf{Notes:} $\\alpha_{mi}$ is the share of domestic intermediate inputs (materials), $\\alpha_{vi}$ is the share of imported inputs, and the labor share is computed as $1-\\alpha_{mi}-\\alpha_{vi}$. All shares are derived from Chilean Input-Output Tables 2021. The "Check Sum" column verifies that all shares sum to 1 for each sector.\n\n');

%% Table 2: Sector Characteristics
fprintf('\n\nTABLE 2: SECTOR CHARACTERISTICS\n');
fprintf('================================\n\n');
fprintf('\\begin{table}[H]\n');
fprintf('\\centering\n');
fprintf('\\caption{Sector Characteristics}\n');
fprintf('\\label{tab:sector_characteristics}\n');
fprintf('\\footnotesize\n');
fprintf('\\begin{tabular}{clccc}\n');
fprintf('\\toprule\n');
fprintf('\\textbf{Sec.} & \\textbf{Sector Name} & \\textbf{Type} & \\textbf{Import} & \\textbf{Material} \\\\\n');
fprintf('              &                      &                & Intensity & Intensity \\\\\n');
fprintf('\\midrule\n');

% Classify sectors
spend_good = alldata.spend_good;
spend_serv = alldata.spend_serv;

for i = 1:nsec
    if spend_good(i) > spend_serv(i)
        sector_type = 'Goods';
    else
        sector_type = 'Services';
    end
    
    % Classify import intensity
    if alpha_V(i) > 0.2
        import_int = 'High';
    elseif alpha_V(i) > 0.1
        import_int = 'Medium';
    else
        import_int = 'Low';
    end
    
    % Classify material intensity
    if alpha(i) > 0.6
        material_int = 'High';
    elseif alpha(i) > 0.4
        material_int = 'Medium';
    else
        material_int = 'Low';
    end
    
    fprintf('%d & %s & %s & %s & %s \\\\\n', ...
        i, names{i}, sector_type, import_int, material_int);
end

fprintf('\\bottomrule\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n\n');

fprintf('\\textbf{Notes:} Import intensity: High $(>0.20)$, Medium $(0.10-0.20)$, Low $(<0.10)$. Material intensity: High $(>0.60)$, Medium $(0.40-0.60)$, Low $(<0.40)$.\n\n');

%% Summary statistics
fprintf('\n\nSUMMARY STATISTICS\n');
fprintf('==================\n\n');
fprintf('Material Shares (alpha_m):\n');
fprintf('  Mean:   %.3f\n', mean(alpha));
fprintf('  Median: %.3f\n', median(alpha));
fprintf('  Min:    %.3f (Sector %d: %s)\n', min(alpha), find(alpha==min(alpha)), names{alpha==min(alpha)});
fprintf('  Max:    %.3f (Sector %d: %s)\n', max(alpha), find(alpha==max(alpha)), names{alpha==max(alpha)});
fprintf('  Std:    %.3f\n\n', std(alpha));

fprintf('Import Shares (alpha_V):\n');
fprintf('  Mean:   %.3f\n', mean(alpha_V));
fprintf('  Median: %.3f\n', median(alpha_V));
fprintf('  Min:    %.3f (Sector %d: %s)\n', min(alpha_V), find(alpha_V==min(alpha_V)), names{alpha_V==min(alpha_V)});
fprintf('  Max:    %.3f (Sector %d: %s)\n', max(alpha_V), find(alpha_V==max(alpha_V)), names{alpha_V==max(alpha_V)});
fprintf('  Std:    %.3f\n\n', std(alpha_V));

fprintf('Labor Shares:\n');
fprintf('  Mean:   %.3f\n', mean(labor_share));
fprintf('  Median: %.3f\n', median(labor_share));
fprintf('  Min:    %.3f (Sector %d: %s)\n', min(labor_share), find(labor_share==min(labor_share)), names{labor_share==min(labor_share)});
fprintf('  Max:    %.3f (Sector %d: %s)\n', max(labor_share), find(labor_share==max(labor_share)), names{labor_share==max(labor_share)});
fprintf('  Std:    %.3f\n\n', std(labor_share));

%% Correlation analysis
fprintf('CORRELATIONS\n');
fprintf('============\n\n');
fprintf('Correlation between material and import shares: %.3f\n', corr(alpha, alpha_V));
fprintf('Correlation between material and labor shares: %.3f\n', corr(alpha, labor_share));
fprintf('Correlation between import and labor shares: %.3f\n\n', corr(alpha_V, labor_share));

%% Create visualization data
fprintf('=== VISUALIZATION SUGGESTION ===\n');
fprintf('Consider creating a scatter plot showing:\n');
fprintf('  X-axis: Material intensity (alpha_m)\n');
fprintf('  Y-axis: Import intensity (alpha_V)\n');
fprintf('  Size: Labor share\n');
fprintf('  Color: Goods vs Services\n\n');

fprintf('LaTeX tables generated successfully!\n');
