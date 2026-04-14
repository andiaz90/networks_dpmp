% steady_state_table.m
% Creates a comprehensive table of steady state values and economic ratios,
% with an optional comparison column against Chilean long-run data averages.
% Run compute_ss_moments_chile.m first to generate ss_moments_chile.mat.

fprintf('\n========================================\n');
fprintf('STEADY STATE VALUES AND ECONOMIC RATIOS\n');
fprintf('========================================\n\n');

%% Load Chilean data moments for comparison (if available)
ss_file = fullfile(fileparts(mfilename('fullpath')), 'ss_moments_chile.mat');
if isempty(fileparts(mfilename('fullpath'))), ss_file = 'ss_moments_chile.mat'; end
has_data = false;
if exist(ss_file, 'file')
    tmp = load(ss_file);
    if isfield(tmp, 'ss_chile')
        ss_chile  = tmp.ss_chile;
        has_data  = true;
        fprintf('(Comparison data loaded from ss_moments_chile.mat — sample %dQ%d-%dQ%d)\n\n', ...
            ss_chile.sample_start(1), ss_chile.sample_start(2), ...
            ss_chile.sample_end(1),   ss_chile.sample_end(2));
    end
end
if ~has_data
    fprintf('(No ss_moments_chile.mat found — run compute_ss_moments_chile.m for data comparison)\n\n');
end

%% Key Aggregate Variables
fprintf('--- KEY AGGREGATE VARIABLES ---\n');
fprintf('%-30s %12.4f\n', 'Output (Y)', oo_.steady_state(strmatch('Y', M_.endo_names, 'exact')));
fprintf('%-30s %12.4f\n', 'Consumption (C)', oo_.steady_state(strmatch('C', M_.endo_names, 'exact')));
fprintf('%-30s %12.4f\n', 'Total Labor (N)', oo_.steady_state(strmatch('N', M_.endo_names, 'exact')));
fprintf('%-30s %12.4f\n', 'Real Wage (w)', oo_.steady_state(strmatch('w', M_.endo_names, 'exact')));
fprintf('%-30s %12.4f\n', 'Exchange Rate (Q)', oo_.steady_state(strmatch('Q', M_.endo_names, 'exact')));

% Consumption by type
C_idx = strmatch('C', M_.endo_names, 'exact');
Cg_idx = strmatch('C_g', M_.endo_names, 'exact');
Cs_idx = strmatch('C_s', M_.endo_names, 'exact');
if ~isempty(Cg_idx)
    fprintf('%-30s %12.4f\n', 'Consumption Goods (C_g)', oo_.steady_state(Cg_idx));
end
if ~isempty(Cs_idx)
    fprintf('%-30s %12.4f\n', 'Consumption Services (C_s)', oo_.steady_state(Cs_idx));
end

%% Price and Inflation Variables
fprintf('\n--- PRICES AND INFLATION ---\n');
pi_idx = strmatch('pi', M_.endo_names, 'exact');
r_idx = strmatch('r', M_.endo_names, 'exact');
if ~isempty(pi_idx)
    pi_ss = oo_.steady_state(pi_idx);
    fprintf('%-30s %12.4f (%.2f%% ann.)\n', 'Gross Inflation (pi)', pi_ss, 400*(pi_ss-1));
end
if ~isempty(r_idx)
    r_ss = oo_.steady_state(r_idx);
    fprintf('%-30s %12.4f (%.2f%% ann.)\n', 'Interest Rate (r)', r_ss, 400*(r_ss-1));
end

pg_idx = strmatch('p_g', M_.endo_names, 'exact');
ps_idx = strmatch('p_s', M_.endo_names, 'exact');
if ~isempty(pg_idx)
    fprintf('%-30s %12.4f\n', 'Goods Price (p_g)', oo_.steady_state(pg_idx));
end
if ~isempty(ps_idx)
    fprintf('%-30s %12.4f\n', 'Services Price (p_s)', oo_.steady_state(ps_idx));
end

%% External Sector
fprintf('\n--- EXTERNAL SECTOR ---\n');
TB_idx = strmatch('TB', M_.endo_names, 'exact');
X_idx = strmatch('X', M_.endo_names, 'exact');
IMP_idx = strmatch('IMP', M_.endo_names, 'exact');
Bstar_idx = strmatch('Bstar', M_.endo_names, 'exact');

if ~isempty(TB_idx)
    fprintf('%-30s %12.4f\n', 'Trade Balance (TB)', oo_.steady_state(TB_idx));
end
if ~isempty(X_idx)
    fprintf('%-30s %12.4f\n', 'Exports (X)', oo_.steady_state(X_idx));
end
if ~isempty(IMP_idx)
    fprintf('%-30s %12.4f\n', 'Imports (IMP)', oo_.steady_state(IMP_idx));
end
if ~isempty(Bstar_idx)
    fprintf('%-30s %12.4f\n', 'Foreign Debt (Bstar)', oo_.steady_state(Bstar_idx));
end

%% Important Economic Ratios
fprintf('\n--- ECONOMIC RATIOS ---\n');

Y_ss   = oo_.steady_state(strmatch('Y',   M_.endo_names, 'exact'));
C_ss   = oo_.steady_state(strmatch('C',   M_.endo_names, 'exact'));

% Use model GDP (= C + TB, value-added by expenditure identity) as denominator.
% This is directly comparable to national accounts GDP.
% Y (gross output) is NOT used here to avoid the double-counting problem.
GDP_idx = strmatch('GDP', M_.endo_names, 'exact');
if ~isempty(GDP_idx)
    GDP_ss = oo_.steady_state(GDP_idx);
else
    % Fallback: reconstruct from C + TB if GDP not in endo_names
    TB_idx_tmp = strmatch('TB', M_.endo_names, 'exact');
    if ~isempty(TB_idx_tmp)
        GDP_ss = C_ss + oo_.steady_state(TB_idx_tmp);
    else
        GDP_ss = Y_ss;   % last resort: use gross output with a warning
        fprintf('WARNING: GDP variable not found; using gross output Y as denominator.\n');
    end
end

if has_data
    fprintf('%-34s %10s  %10s\n', 'Moment', 'Model', 'Chile Data');
    fprintf('%s\n', repmat('-', 58, 1));
else
    fprintf('%-34s %10s\n', 'Moment', 'Model');
    fprintf('%s\n', repmat('-', 46, 1));
end

% C/GDP
CY_model = 100 * C_ss / GDP_ss;
if has_data
    fprintf('%-34s %9.2f%%  %9.2f%%\n', 'C/GDP (Consumption/GDP)', CY_model, ss_chile.CY_pct);
else
    fprintf('%-34s %9.2f%%\n', 'C/GDP (Consumption/GDP)', CY_model);
end

% Trade balance / GDP
if ~isempty(TB_idx)
    TB_ss    = oo_.steady_state(TB_idx);
    TB_model = 100 * TB_ss / GDP_ss;
    if has_data
        fprintf('%-34s %9.2f%%  %9.2f%%\n', 'TB/GDP (Trade Balance/GDP)', TB_model, ss_chile.TBGDP_pct);
    else
        fprintf('%-34s %9.2f%%\n', 'TB/GDP (Trade Balance/GDP)', TB_model);
    end
end

% Exports / GDP
if ~isempty(X_idx)
    X_ss    = oo_.steady_state(X_idx);
    X_model = 100 * X_ss / GDP_ss;
    if has_data
        fprintf('%-34s %9.2f%%  %9.2f%%\n', 'X/GDP (Exports/GDP)', X_model, ss_chile.XY_pct);
    else
        fprintf('%-34s %9.2f%%\n', 'X/GDP (Exports/GDP)', X_model);
    end
end

% Imports / GDP
if ~isempty(IMP_idx)
    IMP_ss    = oo_.steady_state(IMP_idx);
    IMP_model = 100 * IMP_ss / GDP_ss;
    if has_data
        fprintf('%-34s %9.2f%%  %9.2f%%\n', 'IMP/GDP (Imports/GDP)', IMP_model, ss_chile.IMPY_pct);
    else
        fprintf('%-34s %9.2f%%\n', 'IMP/GDP (Imports/GDP)', IMP_model);
    end
end

% Foreign debt / GDP
if ~isempty(Bstar_idx)
    Bstar_ss  = oo_.steady_state(Bstar_idx);
    Bstar_model = 100 * Bstar_ss / GDP_ss;
    if has_data
        fprintf('%-34s %9.2f%%  %9.2f%%\n', 'Bstar/GDP (For. Debt/GDP)', Bstar_model, 70.0);
    else
        fprintf('%-34s %9.2f%%\n', 'Bstar/GDP (For. Debt/GDP)', Bstar_model);
    end
end

% Gross output / GDP (shows intermediates scale)
YGDP_model = 100 * Y_ss / GDP_ss;
if has_data
    fprintf('%-34s %9.2f%%  %9.2f%%\n', 'Y/GDP (Gross Output/GDP)', YGDP_model, 200.0);
else
    fprintf('%-34s %9.2f%%\n', 'Y/GDP (Gross Output/GDP)', YGDP_model);
end

% Consumption composition (goods vs services)
if ~isempty(Cg_idx) && ~isempty(Cs_idx)
    Cg_ss    = oo_.steady_state(Cg_idx);
    Cs_ss    = oo_.steady_state(Cs_idx);
    Cg_model = 100 * Cg_ss / C_ss;
    Cs_model = 100 * Cs_ss / C_ss;
    if has_data
        fprintf('%-34s %9.2f%%  %9.2f%%\n', 'C_g/C (Goods/Tot. Cons.)',    Cg_model, 100*ss_chile.omG);
        fprintf('%-34s %9.2f%%  %9.2f%%\n', 'C_s/C (Services/Tot. Cons.)', Cs_model, 100*(1-ss_chile.omG));
    else
        fprintf('%-34s %9.2f%%\n', 'C_g/C (Goods/Tot. Cons.)',    Cg_model);
        fprintf('%-34s %9.2f%%\n', 'C_s/C (Services/Tot. Cons.)', Cs_model);
    end
end

%% Sectoral Variables (if nsec_val exists)
if exist('nsec_val', 'var')
    fprintf('\n--- SECTORAL STEADY STATE VALUES ---\n');
    fprintf('%-15s %12s %12s %12s %12s %12s\n', 'Sector', 'Output', 'Labor', 'Price', 'Materials', 'TFP');
    fprintf('%-15s %12s %12s %12s %12s %12s\n', repmat('-', 1, 15), repmat('-', 1, 12), ...
            repmat('-', 1, 12), repmat('-', 1, 12), repmat('-', 1, 12), repmat('-', 1, 12));
    
    for i = 1:min(nsec_val, 12)
        % Get sector name
        if exist('names', 'var') && i <= length(names)
            sector_name = char(names{i});
            if length(sector_name) > 15
                sector_name = sector_name(1:15);
            end
        else
            sector_name = sprintf('Sector %d', i);
        end
        
        % Get steady state values
        Y_idx = strmatch(sprintf('Y_%d', i), M_.endo_names, 'exact');
        L_idx = strmatch(sprintf('L_%d', i), M_.endo_names, 'exact');
        P_idx = strmatch(sprintf('P_%d', i), M_.endo_names, 'exact');
        M_idx = strmatch(sprintf('M_%d', i), M_.endo_names, 'exact');
        A_idx = strmatch(sprintf('A_%d', i), M_.endo_names, 'exact');
        
        Y_val = ~isempty(Y_idx) * oo_.steady_state(Y_idx);
        L_val = ~isempty(L_idx) * oo_.steady_state(L_idx);
        P_val = ~isempty(P_idx) * oo_.steady_state(P_idx);
        M_val = ~isempty(M_idx) * oo_.steady_state(M_idx);
        A_val = ~isempty(A_idx) * exp(oo_.steady_state(A_idx)); % TFP is exp(A)
        
        fprintf('%-15s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
                sector_name, Y_val, L_val, P_val, M_val, A_val);
    end
    
    % Sectoral shares
    fprintf('\n--- SECTORAL OUTPUT SHARES ---\n');
    total_sectoral_Y = 0;
    sectoral_Y = zeros(nsec_val, 1);
    for i = 1:nsec_val
        Y_idx = strmatch(sprintf('Y_%d', i), M_.endo_names, 'exact');
        if ~isempty(Y_idx)
            sectoral_Y(i) = oo_.steady_state(Y_idx);
            total_sectoral_Y = total_sectoral_Y + sectoral_Y(i);
        end
    end

    nsec_print = min(nsec_val, 12);
    if has_data && numel(ss_chile.sec_Y_share) == nsec_print
        fprintf('%-16s %10s  %10s\n', 'Sector', 'Model', 'Chile Data');
        fprintf('%s\n', repmat('-', 40, 1));
    else
        fprintf('%-16s %10s\n', 'Sector', 'Model');
        fprintf('%s\n', repmat('-', 28, 1));
    end

    for i = 1:nsec_print
        if exist('names', 'var') && i <= length(names)
            sector_name = char(names{i});
            if length(sector_name) > 15
                sector_name = sector_name(1:15);
            end
        else
            sector_name = sprintf('Sector %d', i);
        end

        if total_sectoral_Y > 0
            model_share = 100 * sectoral_Y(i) / total_sectoral_Y;
            if has_data && numel(ss_chile.sec_Y_share) == nsec_print
                data_share = ss_chile.sec_Y_share(i);
                fprintf('%-16s %9.2f%%  %9.2f%%\n', sector_name, model_share, data_share);
            else
                fprintf('%-16s %9.2f%%\n', sector_name, model_share);
            end
        end
    end
end

%% Model Parameters (key ones)
fprintf('\n--- KEY MODEL PARAMETERS ---\n');
if exist('beta_val', 'var')
    fprintf('%-30s %12.4f\n', 'Discount Factor (beta)', beta_val);
end
if exist('gamma', 'var')
    fprintf('%-30s %12.4f\n', 'Risk Aversion (gamma)', gamma);
end
if exist('phi_val', 'var')
    fprintf('%-30s %12.4f\n', 'Frisch Elasticity (phi)', phi_val);
end
if exist('epsilon', 'var')
    fprintf('%-30s %12.4f\n', 'Elasticity of Subst. (epsilon)', epsilon);
end
if exist('etastar_val', 'var')
    fprintf('%-30s %12.4f\n', 'Foreign Demand Elast. (etastar)', etastar_val);
end

fprintf('\n========================================\n\n');

%% Export to CSV (optional)
try
    % Create summary table for export
    VarNames = {};
    SteadyState = [];
    
    % Add aggregate variables
    aggregate_vars = {'Y', 'C', 'N', 'w', 'Q', 'pi', 'r', 'TB', 'X'};
    for i = 1:length(aggregate_vars)
        var_idx = strmatch(aggregate_vars{i}, M_.endo_names, 'exact');
        if ~isempty(var_idx)
            VarNames{end+1} = aggregate_vars{i};
            SteadyState(end+1) = oo_.steady_state(var_idx);
        end
    end
    
    % Create table
    T = table(VarNames', SteadyState', 'VariableNames', {'Variable', 'SteadyState'});
    
    % Write to CSV
    writetable(T, 'steady_state_summary.csv');
    fprintf('Steady state values exported to steady_state_summary.csv\n\n');
catch ME
    fprintf('Note: Could not export to CSV (table export failed)\n');
end
