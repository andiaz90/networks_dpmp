function ss_table = extract_steady_state_table(exercise_title, oo_, M_)
% EXTRACT_STEADY_STATE_TABLE Creates a table matching steady_state_table.m format
%
% Inputs:
%   exercise_title - Title/description of the exercise
%   oo_ - Dynare output structure
%   M_ - Dynare model structure
%
% Output:
%   ss_table - Table containing steady state values, ratios, and sectoral data

fprintf('\n========================================\n');
fprintf('EXTRACTING STEADY STATE FOR: %s\n', upper(exercise_title));
fprintf('========================================\n\n');

% Check if Dynare results exist
if nargin < 3
    error('Function requires 3 arguments: exercise_title, oo_, and M_');
end
if isempty(oo_) || isempty(M_)
    error('Dynare output structure (oo_) or model structure (M_) is empty.');
end
if ~isstruct(oo_) || ~isstruct(M_)
    error('oo_ and M_ must be structures.');
end

% Initialize storage
Category = {};
Variable = {};
Value = [];
Description = {};

%% === KEY AGGREGATE VARIABLES ===
fprintf('--- KEY AGGREGATE VARIABLES ---\n');

agg_vars = {'Y', 'C', 'N', 'w', 'Q', 'C_g', 'C_s'};
agg_desc = {'Output', 'Consumption', 'Total Labor', 'Real Wage', 'Exchange Rate', 'Consumption Goods', 'Consumption Services'};

for i = 1:length(agg_vars)
    idx = strmatch(agg_vars{i}, M_.endo_names, 'exact');
    if ~isempty(idx)
        val = oo_.steady_state(idx);
        Category{end+1} = 'Aggregate';
        Variable{end+1} = agg_vars{i};
        Value(end+1) = val;
        Description{end+1} = agg_desc{i};
        fprintf('%-30s %12.4f\n', agg_desc{i}, val);
    end
end

%% === PRICES AND INFLATION ===
fprintf('\n--- PRICES AND INFLATION ---\n');

pi_idx = strmatch('pi', M_.endo_names, 'exact');
r_idx = strmatch('r', M_.endo_names, 'exact');
pg_idx = strmatch('p_g', M_.endo_names, 'exact');
ps_idx = strmatch('p_s', M_.endo_names, 'exact');

if ~isempty(pi_idx)
    pi_ss = oo_.steady_state(pi_idx);
    Category{end+1} = 'Prices';
    Variable{end+1} = 'pi';
    Value(end+1) = pi_ss;
    Description{end+1} = sprintf('Gross Inflation (%.2f%% ann.)', 400*(pi_ss-1));
    fprintf('%-30s %12.4f (%.2f%% ann.)\n', 'Gross Inflation (pi)', pi_ss, 400*(pi_ss-1));
end

if ~isempty(r_idx)
    r_ss = oo_.steady_state(r_idx);
    Category{end+1} = 'Prices';
    Variable{end+1} = 'r';
    Value(end+1) = r_ss;
    Description{end+1} = sprintf('Interest Rate (%.2f%% ann.)', 400*(r_ss-1));
    fprintf('%-30s %12.4f (%.2f%% ann.)\n', 'Interest Rate (r)', r_ss, 400*(r_ss-1));
end

if ~isempty(pg_idx)
    pg_ss = oo_.steady_state(pg_idx);
    Category{end+1} = 'Prices';
    Variable{end+1} = 'p_g';
    Value(end+1) = pg_ss;
    Description{end+1} = 'Goods Price';
    fprintf('%-30s %12.4f\n', 'Goods Price (p_g)', pg_ss);
end

if ~isempty(ps_idx)
    ps_ss = oo_.steady_state(ps_idx);
    Category{end+1} = 'Prices';
    Variable{end+1} = 'p_s';
    Value(end+1) = ps_ss;
    Description{end+1} = 'Services Price';
    fprintf('%-30s %12.4f\n', 'Services Price (p_s)', ps_ss);
end

%% === EXTERNAL SECTOR ===
fprintf('\n--- EXTERNAL SECTOR ---\n');

ext_vars = {'TB', 'X', 'IMP', 'Bstar'};
ext_desc = {'Trade Balance', 'Exports', 'Imports', 'Foreign Debt'};

for i = 1:length(ext_vars)
    idx = strmatch(ext_vars{i}, M_.endo_names, 'exact');
    if ~isempty(idx)
        val = oo_.steady_state(idx);
        Category{end+1} = 'External';
        Variable{end+1} = ext_vars{i};
        Value(end+1) = val;
        Description{end+1} = ext_desc{i};
        fprintf('%-30s %12.4f\n', ext_desc{i}, val);
    end
end

%% === ECONOMIC RATIOS ===
fprintf('\n--- ECONOMIC RATIOS ---\n');

Y_idx     = strmatch('Y',     M_.endo_names, 'exact');
C_idx     = strmatch('C',     M_.endo_names, 'exact');
Cg_idx    = strmatch('C_g',   M_.endo_names, 'exact');
Cs_idx    = strmatch('C_s',   M_.endo_names, 'exact');
TB_idx    = strmatch('TB',    M_.endo_names, 'exact');
X_idx     = strmatch('X',     M_.endo_names, 'exact');
IMP_idx   = strmatch('IMP',   M_.endo_names, 'exact');
Bstar_idx = strmatch('Bstar', M_.endo_names, 'exact');
GDP_idx   = strmatch('GDP',   M_.endo_names, 'exact');

% Use model GDP (= C + TB in .mod file) as denominator for all ratios.
% GDP is value-added by the expenditure identity; Y is gross output (includes
% intermediates) and must NOT be used for ratios compared to national accounts.
Y_ss = oo_.steady_state(Y_idx);
C_ss = oo_.steady_state(C_idx);

if ~isempty(GDP_idx)
    GDP_ss = oo_.steady_state(GDP_idx);
elseif ~isempty(TB_idx)
    GDP_ss = C_ss + oo_.steady_state(TB_idx);   % reconstruct: GDP = C + TB
else
    GDP_ss = Y_ss;
end

if ~isempty(C_idx)
    ratio = 100 * C_ss / GDP_ss;
    Category{end+1} = 'Ratios';
    Variable{end+1} = 'C_GDP_ratio';
    Value(end+1) = ratio;
    Description{end+1} = 'C/GDP (%)';
    fprintf('%-30s %12.2f%%\n', 'C/GDP (Consumption/GDP)', ratio);

    if ~isempty(Cg_idx) && ~isempty(Cs_idx)
        Cg_ss = oo_.steady_state(Cg_idx);
        Cs_ss = oo_.steady_state(Cs_idx);
        Category{end+1} = 'Ratios';
        Variable{end+1} = 'Cg_C_ratio';
        Value(end+1) = 100 * Cg_ss / C_ss;
        Description{end+1} = 'C_g/C (%)';
        fprintf('%-30s %12.2f%%\n', 'C_g/C (Goods/Tot. Cons.)', 100 * Cg_ss / C_ss);

        Category{end+1} = 'Ratios';
        Variable{end+1} = 'Cs_C_ratio';
        Value(end+1) = 100 * Cs_ss / C_ss;
        Description{end+1} = 'C_s/C (%)';
        fprintf('%-30s %12.2f%%\n', 'C_s/C (Services/Tot. Cons.)', 100 * Cs_ss / C_ss);
    end
end

if ~isempty(TB_idx)
    TB_ss = oo_.steady_state(TB_idx);
    ratio = 100 * TB_ss / GDP_ss;
    Category{end+1} = 'Ratios';
    Variable{end+1} = 'TB_GDP_ratio';
    Value(end+1) = ratio;
    Description{end+1} = 'TB/GDP (%)';
    fprintf('%-30s %12.2f%%\n', 'TB/GDP (Trade Balance/GDP)', ratio);
end

if ~isempty(X_idx)
    X_ss = oo_.steady_state(X_idx);
    ratio = 100 * X_ss / GDP_ss;
    Category{end+1} = 'Ratios';
    Variable{end+1} = 'X_GDP_ratio';
    Value(end+1) = ratio;
    Description{end+1} = 'X/GDP (%)';
    fprintf('%-30s %12.2f%%\n', 'X/GDP (Exports/GDP)', ratio);
end

if ~isempty(IMP_idx)
    IMP_ss = oo_.steady_state(IMP_idx);
    ratio = 100 * IMP_ss / GDP_ss;
    Category{end+1} = 'Ratios';
    Variable{end+1} = 'IMP_GDP_ratio';
    Value(end+1) = ratio;
    Description{end+1} = 'IMP/GDP (%)';
    fprintf('%-30s %12.2f%%\n', 'IMP/GDP (Imports/GDP)', ratio);
end

if ~isempty(Bstar_idx)
    Bstar_ss = oo_.steady_state(Bstar_idx);
    ratio = 100 * Bstar_ss / GDP_ss;
    Category{end+1} = 'Ratios';
    Variable{end+1} = 'Bstar_GDP_ratio';
    Value(end+1) = ratio;
    Description{end+1} = 'Bstar/GDP (%)';
    fprintf('%-30s %12.2f%%\n', 'Bstar/GDP (For. Debt/GDP)', ratio);
end

% Show Y/GDP so it's clear how large gross output is relative to value added
ratio_YG = 100 * Y_ss / GDP_ss;
Category{end+1} = 'Ratios';
Variable{end+1} = 'Y_GDP_ratio';
Value(end+1) = ratio_YG;
Description{end+1} = 'Y/GDP (gross output / value added, %)';
fprintf('%-30s %12.2f%%\n', 'Y/GDP (Gross Output/GDP)', ratio_YG);

%% === SECTORAL STEADY STATE VALUES ===
fprintf('\n--- SECTORAL STEADY STATE VALUES ---\n');
fprintf('%-15s %12s %12s %12s %12s %12s\n', 'Sector', 'Output', 'Labor', 'Price', 'Materials', 'TFP');
fprintf('%-15s %12s %12s %12s %12s %12s\n', repmat('-', 1, 15), repmat('-', 1, 12), ...
        repmat('-', 1, 12), repmat('-', 1, 12), repmat('-', 1, 12), repmat('-', 1, 12));

sector_names = {
    'Agropecuario'; 'Minería'; 'Manufactura'; 'Electricidad'; 'Construcción';
    'Comercio'; 'Transporte'; 'Financiero'; 'Inmobiliario'; 'Empresariales';
    'Personales'; 'Adm. Pública'
};

nsec = min(12, length(sector_names));
sectoral_Y = zeros(nsec, 1);
total_Y = 0;

for i = 1:nsec
    Y_idx_s = strmatch(sprintf('Y_%d', i), M_.endo_names, 'exact');
    L_idx_s = strmatch(sprintf('L_%d', i), M_.endo_names, 'exact');
    P_idx_s = strmatch(sprintf('P_%d', i), M_.endo_names, 'exact');
    M_idx_s = strmatch(sprintf('M_%d', i), M_.endo_names, 'exact');
    A_idx_s = strmatch(sprintf('A_%d', i), M_.endo_names, 'exact');
    
    Y_val = (~isempty(Y_idx_s)) * oo_.steady_state(Y_idx_s);
    L_val = (~isempty(L_idx_s)) * oo_.steady_state(L_idx_s);
    P_val = (~isempty(P_idx_s)) * oo_.steady_state(P_idx_s);
    M_val = (~isempty(M_idx_s)) * oo_.steady_state(M_idx_s);
    A_val = (~isempty(A_idx_s)) * exp(oo_.steady_state(A_idx_s));
    
    fprintf('%-15s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
            sector_names{i}, Y_val, L_val, P_val, M_val, A_val);
    
    % Store for table
    if ~isempty(Y_idx_s)
        Category{end+1} = sprintf('Sector %d', i);
        Variable{end+1} = sprintf('Y_%d', i);
        Value(end+1) = Y_val;
        Description{end+1} = sprintf('Output - %s', sector_names{i});
        sectoral_Y(i) = Y_val;
        total_Y = total_Y + Y_val;
    end
    if ~isempty(L_idx_s)
        Category{end+1} = sprintf('Sector %d', i);
        Variable{end+1} = sprintf('L_%d', i);
        Value(end+1) = L_val;
        Description{end+1} = sprintf('Labor - %s', sector_names{i});
    end
    if ~isempty(P_idx_s)
        Category{end+1} = sprintf('Sector %d', i);
        Variable{end+1} = sprintf('P_%d', i);
        Value(end+1) = P_val;
        Description{end+1} = sprintf('Price - %s', sector_names{i});
    end
    if ~isempty(M_idx_s)
        Category{end+1} = sprintf('Sector %d', i);
        Variable{end+1} = sprintf('M_%d', i);
        Value(end+1) = M_val;
        Description{end+1} = sprintf('Materials - %s', sector_names{i});
    end
    if ~isempty(A_idx_s)
        Category{end+1} = sprintf('Sector %d', i);
        Variable{end+1} = sprintf('A_%d', i);
        Value(end+1) = A_val;
        Description{end+1} = sprintf('TFP - %s', sector_names{i});
    end
end

%% === SECTORAL OUTPUT SHARES ===
fprintf('\n--- SECTORAL OUTPUT SHARES ---\n');
for i = 1:nsec
    if total_Y > 0
        share = 100 * sectoral_Y(i) / total_Y;
        fprintf('%-15s %12.2f%%\n', sector_names{i}, share);
        Category{end+1} = 'Sectoral Share';
        Variable{end+1} = sprintf('Y_%d_share', i);
        Value(end+1) = share;
        Description{end+1} = sprintf('Output Share - %s (%%)', sector_names{i});
    end
end

fprintf('\n========================================\n\n');

%% Create final table
ss_table = table(Category', Variable', Value', Description', ...
    'VariableNames', {'Category', 'Variable', 'Value', 'Description'});

fprintf('✓ Extracted %d steady state values and ratios\n\n', height(ss_table));

end
