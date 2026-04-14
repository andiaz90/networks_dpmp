% Extract parameter values from params_val.mat for LaTeX tables
clear;
load('params_val.mat');

fprintf('\\n=== SECTORAL PARAMETERS ===\n');
fprintf('Sector | alpha_mi | alpha_vi | kappa | varrho\n');
fprintf('-------|----------|----------|-------|-------\n');
for i = 1:12
    fprintf('%2d     | %.3f    | %.3f    | %.1f  | %.3f\n', ...
        i, modalpha(i), modalphaV(i), modkappa(i), modvarrho(i));
end

fprintf('\n=== CONSUMPTION SHARES (GOODS) ===\n');
fprintf('Sector | omega_G\n');
fprintf('-------|--------\n');
for i = 1:12
    fprintf('%2d     | %.3f\n', i, modgammag(i));
end

fprintf('\n=== CONSUMPTION SHARES (SERVICES) ===\n');
fprintf('Sector | omega_S\n');
fprintf('-------|--------\n');
for i = 1:12
    fprintf('%2d     | %.3f\n', i, modgammas(i));
end

fprintf('\nTotal omega_G: %.4f\n', sum(modgammag));
fprintf('Total omega_S: %.4f\n', sum(modgammas));

%% === SECTORAL CONSUMPTION TABLE (FOR LATEX) ===
fprintf('\n=== SECTORAL CONSUMPTION TABLE (LaTeX Format) ===\n');

% Sector names
sector_names = {
    'Agriculture, Forestry \& Fishing';
    'Mining';
    'Manufacturing';
    'Electricity, Gas, Water \& Waste';
    'Construction';
    'Retail, Hotels \& Restaurants';
    'Transport, Communications \& Info.';
    'Financial Intermediation';
    'Real Estate \& Housing Services';
    'Business Services';
    'Personal Services';
    'Public Administration'
};

% Calculate total consumption by sector (C_gi_ss + C_si_ss)
C_total_ss = C_gi_ss + C_si_ss;

% Calculate shares (as percentage of total consumption)
C_shares = 100 * C_total_ss / sum(C_total_ss);

fprintf('\nLaTeX Table Rows:\n');
fprintf('\\midrule\n');
for i = 1:12
    fprintf('%d & %s & %.3f & %.3f & %.3f & %.1f \\\\\n', ...
        i, sector_names{i}, modgammag(i), modgammas(i), C_total_ss(i), C_shares(i));
end
fprintf('\\midrule\n');
fprintf('\\multicolumn{2}{l}{\\textbf{Total}} & %.3f & %.3f & %.3f & %.1f \\\\\n', ...
    sum(modgammag), sum(modgammas), sum(C_total_ss), sum(C_shares));
fprintf('\\bottomrule\n');

% Also save to a CSV file for easy reference
fprintf('\nSaving consumption table to consumption_table.csv...\n');
T = table((1:12)', sector_names, modgammag, modgammas, C_total_ss, C_shares, ...
    'VariableNames', {'Sector', 'Name', 'omega_G', 'omega_S', 'C_level', 'Share_pct'});
writetable(T, 'consumption_table.csv');
fprintf('Table saved successfully!\n');
