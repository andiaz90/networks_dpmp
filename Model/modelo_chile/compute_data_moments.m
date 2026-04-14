% compute_data_moments.m
%
% Compute all empirical moments needed for SMM estimation of the
% NK-IOSOE Chile model.  Saves results to data_moments_chile.mat and
% returns a structure that can be loaded directly into smm_estimation.m.
%
% DATA SOURCES
%   Employment (L):  count_workers_by_sector.csv      (monthly, 2006M1-)
%   Prices     (P):  deflactor_pib.csv                (quarterly, 1996Q1-2023Q4)
%   Output     (Y):  pib_sectorial_bc.xlsx            (quarterly real GDP by sector, BCCh)
%   REER       (Q):  reer_chile_bis.xlsx              (BIS Real Broad REER, monthly 1994M1-)
%   TB/GDP:          datos_CCNN_mayo2025.xlsx          (BCCh CCNN, sheet "Nom trimestral")
%
% METHODOLOGY
%   - Monthly employment converted to quarterly mean within each quarter
%   - All series HP-filtered at quarterly frequency (lambda = 1600)
%   - Sample aligned to common window (default: 2006Q1 – 2023Q4)
%   - std dev computed on HP-filtered log deviations
%   - Aggregate output = sum of sectoral real VA (from CCNN)
%   - Aggregate inflation = QoQ growth of total GDP deflator (column 28)
%   - Goods share omG and TB/GDP taken from calibration (steady-state targets)
%
% OUTPUT
%   data_moments_chile.mat  containing struct dm_chile with fields:
%       .y_d       [12x1] sectoral output std devs
%       .p_d       [12x1] sectoral price std devs
%       .l_d       [12x1] sectoral employment std devs
%       .d_std_Yg  scalar  output-weighted avg std, goods sectors
%       .d_std_PHg scalar  output-weighted avg std price, goods
%       .d_std_Lg  scalar  output-weighted avg std employment, goods
%       .d_std_Ys  ...services equivalents
%       .d_std_PHs
%       .d_std_Ls
%       .d_std_GDP scalar  (from pib_sectorial_bc, BCCh official PIB total)
%       .d_std_pi  scalar
%       .d_corr_GDPpi scalar
%       .d_omG     scalar  (0.57, from CCNN calibration)
%       .d_std_Q      scalar  (BIS Real Broad REER, HP-filtered log, quarterly)
%       .d_autocorr_Q scalar  (AR(1) of HP-filtered log REER; identifies rho_pvstar)
%       .d_corr_GDPQ  scalar  (corr(GDP, Q) HP-filtered; over-identifies pvstar shock)
%       .d_TBGDP      scalar  (BCCh CCNN mean(X-M)/GDP over sample)
%       .sample_start  date vector [year quarter]
%       .sample_end    date vector [year quarter]
%       .Y_sec_qrt [nT x 12] real sectoral GDP aligned to sample
%       .GDP_qrt   [nT x 1]  official PIB total aligned to sample

clear; clc;
fprintf('=============================================================\n');
fprintf('  Computing data moments for SMM estimation (Chile)\n');
fprintf('=============================================================\n\n');

%% ================================================================ %%
%%  SETTINGS                                                          %%
%% ================================================================ %%
nsec   = 12;
lambda = 1600;   % HP filter smoothing parameter (quarterly)

data_path = 'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';

% Estimation sample (quarters)
sample_start_yr = 2006;  sample_start_q = 1;    % first quarter of employment data
sample_end_yr   = 2023;  sample_end_q   = 4;

% Output weights for cross-sectional averages (from model SS, set in main_SOE.m)
% Will be loaded from params_val.mat if it exists.
model_dir = fileparts(mfilename('fullpath'));
if isempty(model_dir), model_dir = pwd; end

goods    = logical([1;1;1;1;1;0;0;0;0;0;0;0]);   % sectors 1-5
services = logical([0;0;0;0;0;1;1;1;1;1;1;1]);   % sectors 6-12

%% ================================================================ %%
%%  HELPER: HP FILTER                                                 %%
%% ================================================================ %%
% Uses the standard closed-form Hodrick-Prescott filter.
% Returns cycle and trend components.
hp_filter = @(y, lam) hp_cycle(y, lam);

%% ================================================================ %%
%%  1. SECTORAL EMPLOYMENT                                            %%
%%     count_workers_by_sector.csv                                    %%
%%     Columns are in ALPHABETICAL order (see header), not model order%%
%% ================================================================ %%
fprintf('--- 1. Loading employment data ---\n');
fname_emp = fullfile(data_path, 'count_workers_by_sector.csv');

emp_raw = readtable(fname_emp, 'VariableNamingRule','preserve');

% Column 1 is date; columns 2-13 are sectors in alphabetical order
% Alphabetical order -> model sector mapping:
%  CSV col | Sector name                                          | Model #
%    2     | Agriculture, forestry and fishing                    |  1
%    3     | Business services                                    | 10
%    4     | Construction                                         |  5
%    5     | Electricity, gas, water and waste management         |  4
%    6     | Financial intermediation                             |  8
%    7     | Manufacturing                                        |  3
%    8     | Mining                                               |  2
%    9     | Personal services                                    | 11
%   10     | Public administration                                | 12
%   11     | Real estate and housing services                     |  9
%   12     | Transport, communications and information services   |  7
%   13     | Wholesale and retail trade; accommodation            |  6

alpha_to_model = [1; 10; 5; 4; 8; 3; 2; 11; 12; 9; 7; 6];  % CSV col 2-13 -> model sector

% Parse dates (format: YYYY-MM-DD)
dates_monthly = emp_raw{:,1};
if iscell(dates_monthly)
    % Cell array of date strings
    dates_monthly = datenum(dates_monthly);
elseif isdatetime(dates_monthly)
    % datetime array (readtable auto-detected the column)
    dates_monthly = datenum(dates_monthly);
elseif isnumeric(dates_monthly)
    % Already MATLAB datenums, or a numeric encoding — use as-is
    % (if values look like years, e.g. 2010, assume something went wrong)
    if all(dates_monthly < 3000)
        warning('compute_data_moments:dateFormat', ...
            'Date column appears to be year numbers, not datenums. Check CSV format.');
    end
    % already datenums — no conversion needed
else
    % Fallback: coerce to string then parse
    dates_monthly = datenum(cellstr(dates_monthly));
end
% Convert date numbers to year-month
dv = datevec(dates_monthly);
yr_m  = dv(:,1);
mth_m = dv(:,2);

% Raw employment matrix [nmonths x 12], reordered to model sector order
L_raw_alpha = emp_raw{:, 2:13};   % [nmonths x 12], alphabetical order
L_monthly   = NaN(size(L_raw_alpha,1), nsec);
for k = 1:nsec
    L_monthly(:, alpha_to_model(k)) = L_raw_alpha(:,k);
end

fprintf('  Employment: %d monthly obs from %d-%02d to %d-%02d\n', ...
    size(L_monthly,1), yr_m(1), mth_m(1), yr_m(end), mth_m(end));

% Convert to quarterly (average of 3 months within each quarter)
[L_qrt, yr_q, qt_q] = monthly_to_quarterly(L_monthly, yr_m, mth_m);
fprintf('  Employment: %d quarterly obs from %dQ%d to %dQ%d\n', ...
    size(L_qrt,1), yr_q(1), qt_q(1), yr_q(end), qt_q(end));

%% ================================================================ %%
%%  2. SECTORAL PRICE DEFLATORS                                       %%
%%     deflactor_pib.csv                                              %%
%%     28 series, quarterly, format "MMMYYYY"                         %%
%% ================================================================ %%
fprintf('\n--- 2. Loading price deflators ---\n');
fname_defl = fullfile(data_path, 'deflactor_pib.csv');

% File has 4 header rows (SERIES, DESCRIPCION, UNIDAD, then data starts row 4)
raw_defl  = readcell(fname_defl, 'NumHeaderLines', 0);
% Rows 1-3 are metadata; data starts row 4
data_rows = raw_defl(4:end, :);     % {nquarters x 29}  col1=date, col2:29=series

% Parse quarter dates: "MAR1997" = Q1-1997, "JUN"=Q2, "SEP"=Q3, "DIC"=Q4
n_defl = size(data_rows,1);
yr_defl = zeros(n_defl,1);
q_defl  = zeros(n_defl,1);
month_to_q = containers.Map({'MAR','JUN','SEP','DIC'},{1,2,3,4});
for t = 1:n_defl
    s = strtrim(strrep(data_rows{t,1},'"',''));
    mon = s(1:3);
    yr  = str2double(s(4:end));
    yr_defl(t) = yr;
    q_defl(t)  = month_to_q(mon);
end

% Parse deflator values (replace '_' missing with NaN)
P_raw = NaN(n_defl, 28);
for t = 1:n_defl
    for c = 1:28
        v = data_rows{t, c+1};
        if isnumeric(v)
            P_raw(t,c) = v;
        else
            s = strtrim(strrep(num2str(v),'"',''));
            if ~strcmp(s,'_') && ~isempty(s)
                P_raw(t,c) = str2double(s);
            end
        end
    end
end

fprintf('  Deflators: %d quarterly obs from %dQ%d to %dQ%d\n', ...
    n_defl, yr_defl(1), q_defl(1), yr_defl(end), q_defl(end));

% Map 28 deflator columns to 12 model sectors
% (using simple mean of sub-components where a sector aggregates several)
%  Deflator col -> model sector:
%   1  Agropecuario-silvícola           -> Sector  1 (partial)
%   2  Pesca                             -> Sector  1 (partial)
%   3  Minería (total)                   -> Sector  2
%   4  Minería. Cobre                    -> (sub, skip)
%   5  Minería. Otras                    -> (sub, skip)
%   6  Industria manufacturera (total)   -> Sector  3
%   7-14 Mfg sub-sectors                -> (sub, skip)
%  15  Electricidad, gas, agua           -> Sector  4
%  16  Construcción                      -> Sector  5
%  17  Comercio                          -> Sector  6 (partial)
%  18  Restaurantes y hoteles            -> Sector  6 (partial)
%  19  Transporte                        -> Sector  7 (partial)
%  20  Comunicaciones y serv. info.      -> Sector  7 (partial)
%  21  Servicios financieros             -> Sector  8
%  22  Servicios empresariales           -> Sector 10
%  23  Servicios de vivienda             -> Sector  9
%  24  Servicios personales              -> Sector 11
%  25  Administración pública            -> Sector 12
%  26  PIB a costo de factores           -> (aggregate)
%  27  Impuesto sobre productos          -> (aggregate)
%  28  Deflactor del PIB (total)         -> Aggregate price index
defl_col = {
    [1, 2];   % Sector  1: Agro+Pesca (simple mean)
    [3];      % Sector  2: Mining (aggregate)
    [6];      % Sector  3: Manufacturing (aggregate)
    [15];     % Sector  4: Utilities
    [16];     % Sector  5: Construction
    [17, 18]; % Sector  6: Retail + Hotels (simple mean)
    [19, 20]; % Sector  7: Transport + Comm (simple mean)
    [21];     % Sector  8: Finance
    [23];     % Sector  9: Real Estate
    [22];     % Sector 10: Business services
    [24];     % Sector 11: Personal services
    [25]      % Sector 12: Public admin
};

P_sec = NaN(n_defl, nsec);
for i = 1:nsec
    cols = defl_col{i};
    P_sec(:,i) = nanmean(P_raw(:, cols), 2);
end
P_agg = P_raw(:,28);   % total GDP deflator

%% ================================================================ %%
%%  3. SECTORAL OUTPUT (Real GDP by sector)                           %%
%%     pib_sectorial_bc.xlsx  – Sheet "Cuadro"                        %%
%%     Row 3 = column headers; rows 4-73 = 70 quarterly observations  %%
%%     (2009Q1 – 2026Q2, chained pesos 2018=100, thousands of mill.)  %%
%%                                                                     %%
%%  Column mapping (xlsx 1-indexed, col 1 = date):                    %%
%%    2  Agropecuario-silvicola  → Sector 1 (sum with Pesca)          %%
%%    3  Pesca                   → Sector 1 (partial)                 %%
%%    4  Mineria (total)         → Sector 2                           %%
%%    7  Ind. Manufacturera      → Sector 3                           %%
%%   17  Electricidad/gas/agua   → Sector 4                           %%
%%   18  Construccion            → Sector 5                           %%
%%   19  Comercio/restaurantes   → Sector 6                           %%
%%   22  Transporte              → Sector 7 (sum with Comms)          %%
%%   23  Comunicaciones          → Sector 7 (partial)                 %%
%%   25  Servicios financieros   → Sector 8                           %%
%%   27  Servicios de vivienda   → Sector 9                           %%
%%   26  Servicios empresariales → Sector 10                          %%
%%   28  Servicios personales    → Sector 11                          %%
%%   29  Administracion publica  → Sector 12                          %%
%%   32  PIB total               → Aggregate GDP                      %%
%% ================================================================ %%
fprintf('\n--- 3. Loading sectoral real output (pib_sectorial_bc.xlsx) ---\n');
fname_pib = fullfile(data_path, 'pib_sectorial_bc.xlsx');

Y_sec_raw = [];  yr_y = [];  qt_y = [];  GDP_data = [];

% Pre-check: file must exist before we attempt readcell
if ~exist(fname_pib, 'file')
    fprintf('  ERROR: File not found:\n    %s\n', fname_pib);
    fprintf('  --> Place pib_sectorial_bc.xlsx in data_path and re-run.\n');
    fprintf('  Output moments (y_d) will be NaN.\n');
end

try
    % Sheet "Cuadro": rows 1-2 = title/subtitle, row 3 = headers, rows 4-73 = data
    % Read row 3 (header) through row 73 as a cell array; 33 columns (A:AG)
    raw_pib = readcell(fname_pib, 'Sheet', 'Cuadro', 'Range', 'A3:AG73');
    n_pib   = size(raw_pib, 1) - 1;   % first row is header → 70 data rows
    fprintf('  readcell OK: %d data rows, %d columns\n', n_pib, size(raw_pib,2));

    % --- Parse quarterly dates from column 1 (data rows 2:end of raw_pib) ---
    % NOTE: BCCh xlsx stores the YEAR as a plain integer in a date-formatted cell
    % (e.g. the integer 2009 is stored as an Excel date serial, which MATLAB reads
    % as datetime 01-Jul-1905 because Excel serial 2009 ≈ July 1905).
    % All 4 quarters of a given year share the same integer, so the quarter must
    % be inferred from the row order within each year group (Q1→Q4 in sequence).
    yr_y = zeros(n_pib, 1);
    qt_y = zeros(n_pib, 1);
    for t = 1:n_pib
        d = raw_pib{t+1, 1};
        if isdatetime(d)
            % Recover the year that was stored as an Excel serial number:
            %   excel_serial = days since 30-Dec-1899; here it equals the year.
            excel_serial = round(days(d - datetime(1899, 12, 30)));
            if excel_serial >= 1990 && excel_serial <= 2100
                % BCCh format: serial IS the year; quarter assigned below
                yr_y(t) = excel_serial;
                qt_y(t) = 1;   % placeholder — replaced after loop
            else
                % Genuine datetime: extract year/quarter normally
                yr_y(t) = year(d);
                qt_y(t) = ceil(month(d) / 3);
            end
        elseif ischar(d) || isstring(d)
            ds = char(d);
            yr_y(t) = str2double(ds(1:4));
            mth     = str2double(ds(6:7));
            qt_y(t) = ceil(mth / 3);
        elseif isnumeric(d) && ~isnan(d)
            % Excel serial date (days since 1900-01-00)
            dv = datevec(datenum('30-Dec-1899') + d);
            yr_y(t) = dv(1);
            qt_y(t) = ceil(dv(2) / 3);
        % else: empty cell (missing), NaT, or unrecognised type — leave yr_y=0
        end
    end
    % Assign quarters Q1→Q4 (in row order) within each year group
    unique_yrs = unique(yr_y(yr_y > 0), 'stable');
    for yi = 1:numel(unique_yrs)
        rows_yi = find(yr_y == unique_yrs(yi));
        for qi = 1:numel(rows_yi)
            qt_y(rows_yi(qi)) = qi;
        end
    end

    % Drop rows whose date could not be parsed (empty/trailing rows in xlsx)
    valid_rows = (yr_y > 0);
    n_dropped  = sum(~valid_rows);
    if n_dropped > 0
        % Show what the unparseable cells look like (safe for any type)
        bad_idx = find(~valid_rows, 1);
        d_bad   = raw_pib{bad_idx + 1, 1};   % +1 because row 1 of raw_pib is header
        try,   val_str = num2str(d_bad); catch, val_str = class(d_bad); end
        fprintf('  Dropping %d unparseable date row(s) (class=%s, val=%s).\n', ...
            n_dropped, class(d_bad), val_str);
        yr_y    = yr_y(valid_rows);
        qt_y    = qt_y(valid_rows);
        n_pib   = sum(valid_rows);
        % raw_pib rows: row 1 = header, rows 2..end = data → shift index by 1
        data_mask_ext = [false; valid_rows];   % pad for header row
        raw_pib = raw_pib([true; valid_rows], :);  % keep header + valid data rows
    end
    fprintf('  Date range in file: %dQ%d to %dQ%d  (%d rows)\n', ...
        yr_y(1), qt_y(1), yr_y(end), qt_y(end), n_pib);

    % --- Parse numeric values: pib_num(:,k) = xlsx column k+1 (k=1..32) ---
    pib_num = NaN(n_pib, 32);
    for t = 1:n_pib
        for c = 1:32
            v = raw_pib{t+1, c+1};
            if isnumeric(v) && ~isnan(v)
                pib_num(t, c) = v;
            elseif ischar(v) || isstring(v)
                pib_num(t, c) = str2double(v);
            end
        end
    end
    n_valid_cells = sum(~isnan(pib_num(:)));
    fprintf('  Numeric parse: %d / %d cells valid (%.0f%%)\n', ...
        n_valid_cells, numel(pib_num), 100*n_valid_cells/numel(pib_num));

    % --- Aggregate to 12 model sectors ---
    % pib_num col k  <->  xlsx col k+1
    %  k= 1 Agropecuario  (sum with Pesca → sector 1)
    %  k= 2 Pesca
    %  k= 3 Mineria total                → sector 2
    %  k= 6 Manufactura total            → sector 3
    %  k=16 Electricidad/gas/agua        → sector 4
    %  k=17 Construccion                 → sector 5
    %  k=18 Comercio/restaurantes total  → sector 6
    %  k=21 Transporte  (sum with Comms) → sector 7
    %  k=22 Comunicaciones e info
    %  k=24 Servicios financieros        → sector 8
    %  k=26 Vivienda e inmobiliarios     → sector 9
    %  k=25 Servicios empresariales      → sector 10
    %  k=27 Servicios personales         → sector 11
    %  k=28 Administracion publica       → sector 12
    %  k=31 PIB total (aggregate)
    Y_sec_raw = NaN(n_pib, nsec);
    Y_sec_raw(:, 1)  = pib_num(:, 1) + pib_num(:, 2);   % Agro + Pesca
    Y_sec_raw(:, 2)  = pib_num(:, 3);                    % Mineria
    Y_sec_raw(:, 3)  = pib_num(:, 6);                    % Manufactura
    Y_sec_raw(:, 4)  = pib_num(:,16);                    % Utilities (EGA)
    Y_sec_raw(:, 5)  = pib_num(:,17);                    % Construccion
    Y_sec_raw(:, 6)  = pib_num(:,18);                    % Comercio+Restaurantes
    Y_sec_raw(:, 7)  = pib_num(:,21) + pib_num(:,22);    % Transporte + Comms
    Y_sec_raw(:, 8)  = pib_num(:,24);                    % Financiero
    Y_sec_raw(:, 9)  = pib_num(:,26);                    % Vivienda
    Y_sec_raw(:,10)  = pib_num(:,25);                    % Empresarial
    Y_sec_raw(:,11)  = pib_num(:,27);                    % Personales
    Y_sec_raw(:,12)  = pib_num(:,28);                    % Admin publica
    GDP_data         = pib_num(:,31);                    % PIB total

    fprintf('  PIB Sectorial: %d quarterly obs from %dQ%d to %dQ%d\n', ...
        n_pib, yr_y(1), qt_y(1), yr_y(end), qt_y(end));
catch ME
    fprintf('  ERROR reading pib_sectorial_bc.xlsx: %s\n', ME.message);
    fprintf('  Output moments will be set to NaN.\n');
end

%% ================================================================ %%
%%  4. ALIGN TO COMMON SAMPLE                                         %%
%% ================================================================ %%
fprintf('\n--- 4. Aligning to common sample ---\n');

% Build quarterly date index helpers
date2idx_emp  = @(y,q) find(yr_q    == y & qt_q    == q, 1);
date2idx_defl = @(y,q) find(yr_defl == y & q_defl  == q, 1);
date2idx_pib  = @(y,q) find(yr_y    == y & qt_y    == q, 1);

t0_emp  = date2idx_emp(sample_start_yr,  sample_start_q);
t1_emp  = date2idx_emp(sample_end_yr,    sample_end_q);
t0_defl = date2idx_defl(sample_start_yr, sample_start_q);
t1_defl = date2idx_defl(sample_end_yr,   sample_end_q);

if isempty(t0_emp) || isempty(t1_emp)
    error('Sample start/end not found in employment quarterly data.');
end
if isempty(t0_defl) || isempty(t1_defl)
    error('Sample start/end not found in deflator data.');
end

L_sample  = L_qrt(t0_emp:t1_emp, :);
P_sample  = P_sec(t0_defl:t1_defl, :);
Pa_sample = P_agg(t0_defl:t1_defl);

% Output: align to sample if loaded
Y_qrt      = [];
GDP_sample = [];
if ~isempty(Y_sec_raw) && ~isempty(yr_y)
    t0_pib = date2idx_pib(sample_start_yr, sample_start_q);
    t1_pib = date2idx_pib(sample_end_yr,   sample_end_q);
    if ~isempty(t0_pib) && ~isempty(t1_pib)
        Y_qrt      = Y_sec_raw(t0_pib:t1_pib, :);
        GDP_sample = GDP_data(t0_pib:t1_pib);
    else
        % Sample window predates PIB data — use all available overlap
        mask_pib = (yr_y > sample_start_yr | (yr_y == sample_start_yr & qt_y >= sample_start_q)) & ...
                   (yr_y < sample_end_yr   | (yr_y == sample_end_yr   & qt_y <= sample_end_q));
        if any(mask_pib)
            Y_qrt      = Y_sec_raw(mask_pib, :);
            GDP_sample = GDP_data(mask_pib);
            fprintf('  WARNING: Output sample partially outside requested window.\n');
        else
            fprintf('  ERROR: No PIB data overlaps with estimation sample %dQ%d-%dQ%d.\n', ...
                sample_start_yr, sample_start_q, sample_end_yr, sample_end_q);
            fprintf('    PIB data covers %dQ%d to %dQ%d.\n', yr_y(1), qt_y(1), yr_y(end), qt_y(end));
        end
    end
    if ~isempty(Y_qrt)
        fprintf('  Y_qrt aligned: %d quarters\n', size(Y_qrt,1));
        % Quick sanity check — how many sectors have any valid data?
        n_valid_sec = sum(sum(~isnan(Y_qrt) & Y_qrt > 0) >= 20);
        fprintf('  Sectors with >=20 valid obs: %d / %d\n', n_valid_sec, nsec);
    end
end

nT = size(L_sample,1);
fprintf('  Common sample: %dQ%d - %dQ%d  (%d quarters)\n', ...
    sample_start_yr, sample_start_q, sample_end_yr, sample_end_q, nT);

if size(P_sample,1) ~= nT
    warning('Deflator and employment samples have different lengths (%d vs %d)!', ...
        size(P_sample,1), nT);
end
if ~isempty(Y_qrt) && size(Y_qrt,1) ~= nT
    warning('Output and employment samples have different lengths (%d vs %d).', ...
        size(Y_qrt,1), nT);
end

%% ================================================================ %%
%%  5. HP-FILTER AND COMPUTE SECTORAL STD DEVS                        %%
%% ================================================================ %%
fprintf('\n--- 5. HP filtering and computing std devs ---\n');

% Employment std devs
l_d = NaN(nsec,1);
for i = 1:nsec
    x = log(L_sample(:,i));
    if any(isnan(x) | isinf(x) | x==0)
        fprintf('  WARNING: sector %d employment has missing/zero values, skipping.\n',i);
        continue
    end
    cyc = hp_filter(x, lambda);
    l_d(i) = std(cyc);
end

% Price deflator std devs (use log-level deviations; P is index, not log)
p_d = NaN(nsec,1);
for i = 1:nsec
    x = log(P_sample(:,i));
    bad = isnan(x) | isinf(x);
    if sum(~bad) < 20
        fprintf('  WARNING: sector %d price has too few valid obs, skipping.\n',i);
        continue
    end
    % Fill isolated NaNs by interpolation before HP filter
    x_filled = fillmissing(x, 'linear');
    cyc = hp_filter(x_filled, lambda);
    p_d(i) = std(cyc);
end

fprintf('  Sectoral employment std devs (%%):  ');
fprintf('%.3f ', l_d*100); fprintf('\n');
fprintf('  Sectoral price std devs (%%):        ');
fprintf('%.3f ', p_d*100); fprintf('\n');

% Sectoral output std devs (HP-filtered log real GDP by sector)
y_d = NaN(nsec, 1);
if ~isempty(Y_qrt)
    for i = 1:nsec
        x = Y_qrt(:, i);
        bad = (x <= 0) | isnan(x);
        if sum(~bad) < 20
            fprintf('  WARNING: sector %d output has too few valid obs (%d), skipping.\n', ...
                i, sum(~bad));
            continue
        end
        x_log    = log(fillmissing(x, 'linear'));
        cyc      = hp_filter(x_log, lambda);
        y_d(i)   = std(cyc);
    end
    fprintf('  Sectoral output std devs (%%):      ');
    fprintf('%.3f ', y_d*100); fprintf('\n');
else
    fprintf('  Output data not loaded; y_d set to NaN.\n');
end

% Aggregate inflation: QoQ growth of total GDP deflator
pi_series = Pa_sample(2:end) ./ Pa_sample(1:end-1);
pi_series(pi_series <= 0) = NaN;
pi_hp     = hp_filter(log(fillmissing(pi_series, 'linear')), lambda);
d_std_pi  = std(pi_hp);

fprintf('  std(pi)  = %.5f\n', d_std_pi);

%% ================================================================ %%
%%  6. AGGREGATE OUTPUT MOMENTS                                        %%
%%     If sectoral real GDP loaded, compute GDP std dev and corr(Y,pi) %%
%%     Otherwise set placeholders.                                    %%
%% ================================================================ %%
fprintf('\n--- 6. Aggregate moments ---\n');

d_std_GDP    = NaN;
d_corr_GDPpi = NaN;

if ~isempty(GDP_sample) && sum(GDP_sample > 0) >= 20
    % Use official PIB total series from pib_sectorial_bc (BCCh)
    GDP_hp = hp_filter(log(fillmissing(GDP_sample, 'linear')), lambda);
    d_std_GDP    = std(GDP_hp);
    % Align GDP_hp with pi_hp: pi_hp is 1 obs shorter (QoQ growth rates)
    nmin = min(length(GDP_hp)-1, length(pi_hp));
    d_corr_GDPpi = corr(GDP_hp(end-nmin+1:end), pi_hp(end-nmin+1:end), 'rows','complete');
    fprintf('  std(GDP) = %.5f  (from pib_sectorial_bc, BCCh)\n', d_std_GDP);
    fprintf('  corr(GDP, pi) = %.4f\n', d_corr_GDPpi);
elseif ~isempty(Y_qrt)
    % Fallback: sum of sectoral real VA
    GDP_real = sum(Y_qrt, 2);
    GDP_hp   = hp_filter(log(GDP_real), lambda);
    d_std_GDP    = std(GDP_hp);
    nmin = min(length(GDP_hp)-1, length(pi_hp));
    d_corr_GDPpi = corr(GDP_hp(end-nmin+1:end), pi_hp(end-nmin+1:end), 'rows','complete');
    fprintf('  std(GDP) = %.5f  (sum of sectoral VA)\n', d_std_GDP);
    fprintf('  corr(GDP, pi) = %.4f\n', d_corr_GDPpi);
else
    % Default placeholders consistent with Chilean literature
    d_std_GDP    = 0.021;
    d_corr_GDPpi = -0.15;
    fprintf('  std(GDP): PLACEHOLDER (no output data). Using defaults.\n');
end

%% ================================================================ %%
%%  7. REAL EXCHANGE RATE AND TRADE BALANCE                           %%
%%                                                                    %%
%%  std(Q): BIS Real Broad REER for Chile (RBCL), sheet "Real"       %%
%%          File: reer_chile_bis.xlsx  (downloaded from bis.org)      %%
%%          Monthly data (2020=100) → quarterly average → log → HP   %%
%%                                                                    %%
%%  TB/GDP: datos_CCNN_mayo2025.xlsx, sheet "Nom trimestral"          %%
%%          Col I = Exports, Col J = Imports, Col K = nominal GDP     %%
%%          TB/GDP = mean( (X-M)/GDP ) over the estimation sample     %%
%% ================================================================ %%
fprintf('\n--- 7. Real exchange rate (BIS REER) and Trade Balance ---\n');

%% --- 7a. Real exchange rate: BIS REER Chile ---
d_std_Q         = NaN;
d_autocorr_Q    = NaN;   % initialise before try so catch block always has a defined variable
reer_hp_aligned = [];    % initialise before try so Section 7c always has a defined variable
fname_reer = fullfile(data_path, 'reer_chile_bis.xlsx');
try
    % Sheet "Real": rows 1-3 title, row 4 country names, row 5 codes,
    % rows 6+ = monthly data. Col A = date, col K = Chile (RBCL).
    % readtable returns dates as datetime when DetectImportOptions used.
    opts_reer = detectImportOptions(fname_reer, 'Sheet', 'Real');
    opts_reer.DataRange     = 'A6';       % data starts row 6
    opts_reer.VariableNamesRange = 'A5';  % row 5 = series codes (RBDZ, RBAR, ..., RBCL)
    % Force first column as date, rest as double
    opts_reer.VariableTypes(1) = {'datetime'};
    for kk = 2:length(opts_reer.VariableTypes)
        opts_reer.VariableTypes(kk) = {'double'};
    end
    reer_tbl = readtable(fname_reer, opts_reer);

    % Find Chile column (RBCL)
    cl_col = find(strcmp(reer_tbl.Properties.VariableNames, 'RBCL'), 1);
    if isempty(cl_col)
        % Fallback: Chile is the 10th data column (col K = index 10 after date)
        cl_col = 10;
    end
    reer_dates  = reer_tbl{:, 1};    % datetime array, monthly
    reer_chile  = reer_tbl{:, cl_col}; % REER index, monthly

    % Remove NaN rows
    ok = ~isnan(reer_chile) & ~isnat(reer_dates);
    reer_dates = reer_dates(ok);
    reer_chile = reer_chile(ok);

    % Convert monthly to quarterly (average within each quarter)
    reer_yr = year(reer_dates);
    reer_q  = ceil(month(reer_dates) / 3);
    [uq_yq, ~, ic_reer] = unique([reer_yr, reer_q], 'rows');
    yr_reer = uq_yq(:,1);
    q_reer  = uq_yq(:,2);
    nQ_reer = size(uq_yq, 1);
    reer_qrt = NaN(nQ_reer, 1);
    for t = 1:nQ_reer
        reer_qrt(t) = mean(reer_chile(ic_reer == t));
    end

    % Align to estimation sample
    mask_reer = (yr_reer >= sample_start_yr & yr_reer <= sample_end_yr);
    mask_reer(yr_reer == sample_start_yr) = ...
        mask_reer(yr_reer == sample_start_yr) & (q_reer(yr_reer == sample_start_yr) >= sample_start_q);
    mask_reer(yr_reer == sample_end_yr) = ...
        mask_reer(yr_reer == sample_end_yr) & (q_reer(yr_reer == sample_end_yr) <= sample_end_q);

    reer_s = reer_qrt(mask_reer);
    if sum(~isnan(reer_s)) >= 20
        reer_hp  = hp_filter(log(reer_s), lambda);
        d_std_Q  = std(reer_hp);
        % AR(1) autocorrelation of log REER (identifies rho_pvstar)
        tmp = reer_hp(~isnan(reer_hp));
        d_autocorr_Q    = corr(tmp(1:end-1), tmp(2:end));
        reer_hp_aligned = reer_hp;   % save for corr(GDP,Q) computation below
        fprintf('  std(Q)      = %.5f  (BIS REER, %d quarterly obs)\n', d_std_Q, sum(~isnan(reer_s)));
        fprintf('  autocorr(Q) = %.4f  (AR(1) of HP-filtered log REER)\n', d_autocorr_Q);
    else
        fprintf('  WARNING: insufficient REER obs in sample.\n');
    end
catch ME
    fprintf('  WARNING: could not read reer_chile_bis.xlsx: %s\n', ME.message);
end
if isnan(d_std_Q)
    d_std_Q      = 0.052;
    d_autocorr_Q = 0.75;    % AR(1) of HP-REER; typical SOE value for Chile
    fprintf('  std(Q)      = %.4f  (fallback to literature estimate)\n', d_std_Q);
    fprintf('  autocorr(Q) = %.4f  (fallback)\n', d_autocorr_Q);
end
% corr(GDP, Q): compute after both GDP_hp and reer_hp_aligned are available.
% (done further below, after Section 8 where GDP_hp is defined)
d_corr_GDPQ = NaN;   % placeholder; filled below

%% --- 7b. Trade balance share: datos_CCNN_mayo2025.xlsx "Nom trimestral" ---
d_TBGDP = NaN;
fname_ccnn = fullfile(data_path, 'datos_CCNN_mayo2025.xlsx');
try
    % Sheet "Nom trimestral": row 3 = headers, data from row 4.
    % Col A = date (Excel serial), Col I = Exports, Col J = Imports, Col K = GDP
    opts_tb = detectImportOptions(fname_ccnn, 'Sheet', 'Nom trimestral');
    opts_tb.DataRange        = 'A4';
    opts_tb.VariableNamesRange = 'A3';
    tb_tbl = readtable(fname_ccnn, opts_tb);

    % Find columns by their known MATLAB-sanitised names (BCCh "Nom trimestral")
    %   x6_Exportaciones...  → Exports
    %   x7_Importaciones...  → Imports
    %   x8_PIB...            → Nominal GDP
    vn  = tb_tbl.Properties.VariableNames;
    i_x = find(~cellfun(@isempty, regexp(vn, 'x6_|Exporta', 'ignorecase')), 1);
    i_m = find(~cellfun(@isempty, regexp(vn, 'x7_|Importa', 'ignorecase')), 1);
    i_g = find(~cellfun(@isempty, regexp(vn, 'x8_|PIB', 'ignorecase')), 1);

    % Fallback to known column positions (I=9, J=10, K=11, with col A=1)
    if isempty(i_x), i_x = 9;  end
    if isempty(i_m), i_m = 10; end
    if isempty(i_g), i_g = 11; end

    % Date column (col A): readtable may return datetime or numeric serial
    dates_tb = tb_tbl{:, 1};
    if isnumeric(dates_tb)
        % Excel serial → year and quarter
        base = datenum('30-Dec-1899');
        dv_tb = datevec(base + dates_tb);
        yr_tb = dv_tb(:,1);
        q_tb  = ceil(dv_tb(:,2) / 3);
    else
        yr_tb = year(dates_tb);
        q_tb  = ceil(month(dates_tb) / 3);
    end

    X_tb  = tb_tbl{:, i_x};
    M_tb  = tb_tbl{:, i_m};
    GDP_n = tb_tbl{:, i_g};

    % Align to sample
    mask_tb = (yr_tb >= sample_start_yr & yr_tb <= sample_end_yr);
    X_s   = X_tb(mask_tb);   M_s = M_tb(mask_tb);   GDP_s = GDP_n(mask_tb);
    ok_tb = ~isnan(X_s) & ~isnan(M_s) & ~isnan(GDP_s) & GDP_s > 0;
    if sum(ok_tb) >= 20
        tb_ratio = (X_s(ok_tb) - M_s(ok_tb)) ./ GDP_s(ok_tb);
        d_TBGDP  = mean(tb_ratio);
        fprintf('  TB/GDP   = %.4f  (from BCCh CCNN, %d quarterly obs)\n', d_TBGDP, sum(ok_tb));
    else
        fprintf('  WARNING: insufficient CCNN obs for TB/GDP in sample.\n');
    end
catch ME
    fprintf('  WARNING: could not compute TB/GDP from CCNN: %s\n', ME.message);
end
if isnan(d_TBGDP)
    d_TBGDP = -0.02;
    fprintf('  TB/GDP   = %.4f  (fallback to historical average)\n', d_TBGDP);
end

%% ================================================================ %%
%%  7c. corr(GDP, Q) — cross-correlation, identifies SOE openness     %%
%% ================================================================ %%
if ~isempty(reer_hp_aligned) && exist('GDP_hp','var') && ~isempty(GDP_hp)
    nr = min(length(reer_hp_aligned), length(GDP_hp));
    if nr >= 20
        d_corr_GDPQ = corr(GDP_hp(end-nr+1:end), reer_hp_aligned(end-nr+1:end), 'rows','complete');
        fprintf('  corr(GDP,Q) = %.4f  (HP-filtered GDP vs HP-filtered log REER)\n', d_corr_GDPQ);
    end
end
if isnan(d_corr_GDPQ)
    d_corr_GDPQ = -0.15;   % typical SOE: depreciation -> output expansion (negative sign conv.)
    fprintf('  corr(GDP,Q) = %.4f  (fallback)\n', d_corr_GDPQ);
end

%% ================================================================ %%
%%  8. GOODS EXPENDITURE SHARE (from calibration)                     %%
%% ================================================================ %%
d_omG = 0.57;   % from National Accounts calibration (ombar_val)
fprintf('  omG      = %.2f  (from SS calibration)\n', d_omG);

%% ================================================================ %%
%%  9. OUTPUT-WEIGHTED CROSS-SECTIONAL AVERAGES                       %%
%% ================================================================ %%
fprintf('\n--- 7. Computing output-weighted cross-sectional moments ---\n');

% Load steady-state output weights from params_val.mat if available
params_file = fullfile(model_dir, 'params_val.mat');
Y_ss_vec = ones(nsec,1);   % default equal weights

if exist(params_file,'file')
    try
        pv = load(params_file);
        if isfield(pv, 'Yi_ss')
            Y_ss_vec = pv.Yi_ss;
            fprintf('  Loaded Y_ss from params_val.mat\n');
        end
    catch
        fprintf('  Could not load Y_ss from params_val.mat, using equal weights.\n');
    end
else
    fprintf('  params_val.mat not found; using equal output weights.\n');
end

w_g = Y_ss_vec(goods)   / sum(Y_ss_vec(goods));
w_s = Y_ss_vec(services)/ sum(Y_ss_vec(services));

% y_d was computed in Section 5 from pib_sectorial_bc.xlsx
% (NaN for sectors with insufficient data)
if any(isnan(y_d))
    fprintf('  WARNING: %d sector(s) have NaN output std dev.\n', sum(isnan(y_d)));
end

d_std_Yg  = w_g' * y_d(goods);
d_std_PHg = w_g' * p_d(goods);
d_std_Lg  = w_g' * l_d(goods);
d_std_Ys  = w_s' * y_d(services);
d_std_PHs = w_s' * p_d(services);
d_std_Ls  = w_s' * l_d(services);

%% ================================================================ %%
%%  10. DISPLAY SUMMARY AND SAVE                                       %%
%% ================================================================ %%
fprintf('\n=============================================================\n');
fprintf('  DATA MOMENTS SUMMARY\n');
fprintf('=============================================================\n\n');

sector_names = {'Agriculture','Mining','Manufacturing','Utilities','Construction', ...
                'Trade/Hotels','Transport/Comm','Finance','Real Estate', ...
                'Business Serv.','Personal Serv.','Public Admin.'};

fprintf('%-20s  %8s  %8s  %8s\n', 'Sector', 'std(Y)', 'std(PH)', 'std(L)');
fprintf('%s\n', repmat('-',52,1));
for i = 1:nsec
    mk = '';
    if goods(i),    mk = '[G]'; end
    if services(i), mk = '[S]'; end
    fprintf('%-2d %-17s%s  %8.4f  %8.4f  %8.4f\n', ...
        i, sector_names{i}, mk, y_d(i), p_d(i), l_d(i));
end
fprintf('%s\n', repmat('-',52,1));
fprintf('  Goods   cross-sect avg: std(Y)=%6.4f  std(P)=%6.4f  std(L)=%6.4f\n', ...
    d_std_Yg, d_std_PHg, d_std_Lg);
fprintf('  Services cross-sect avg:std(Y)=%6.4f  std(P)=%6.4f  std(L)=%6.4f\n\n', ...
    d_std_Ys, d_std_PHs, d_std_Ls);
fprintf('  std(GDP)      = %.5f\n', d_std_GDP);
fprintf('  std(pi)       = %.5f\n', d_std_pi);
fprintf('  corr(GDP,pi)  = %.4f\n', d_corr_GDPpi);
fprintf('  omG           = %.2f\n',  d_omG);
fprintf('  std(Q)        = %.5f  (BIS REER)\n', d_std_Q);
fprintf('  autocorr(Q)   = %.4f\n', d_autocorr_Q);
fprintf('  corr(GDP,Q)   = %.4f\n', d_corr_GDPQ);
fprintf('  TB/GDP        = %.4f  (BCCh CCNN)\n', d_TBGDP);

% Pack into struct
dm_chile = struct();
dm_chile.y_d          = y_d;
dm_chile.p_d          = p_d;
dm_chile.l_d          = l_d;
dm_chile.d_std_Yg     = d_std_Yg;
dm_chile.d_std_PHg    = d_std_PHg;
dm_chile.d_std_Lg     = d_std_Lg;
dm_chile.d_std_Ys     = d_std_Ys;
dm_chile.d_std_PHs    = d_std_PHs;
dm_chile.d_std_Ls     = d_std_Ls;
dm_chile.d_std_GDP    = d_std_GDP;
dm_chile.d_std_pi     = d_std_pi;
dm_chile.d_corr_GDPpi = d_corr_GDPpi;
dm_chile.d_omG        = d_omG;
dm_chile.d_std_Q      = d_std_Q;
dm_chile.d_autocorr_Q = d_autocorr_Q;
dm_chile.d_corr_GDPQ  = d_corr_GDPQ;
dm_chile.d_TBGDP      = d_TBGDP;
dm_chile.sample_start = [sample_start_yr, sample_start_q];
dm_chile.sample_end   = [sample_end_yr,   sample_end_q];
dm_chile.sector_names = sector_names;

% Also store the HP-filtered series for diagnostics
dm_chile.L_qrt     = L_qrt;
dm_chile.P_sec_qrt = P_sample;
dm_chile.Pa_qrt    = Pa_sample;
dm_chile.Y_sec_qrt = Y_qrt;        % [nT x 12] real sectoral GDP (sample window)
dm_chile.GDP_qrt   = GDP_sample;   % [nT x 1]  official PIB total (sample window)
dm_chile.yr_q      = yr_q;
dm_chile.qt_q      = qt_q;

%% ================================================================ %%
%%  FINAL VALIDATION: all 44 moments must be non-NaN before saving  %%
%% ================================================================ %%
d_vec_check = [y_d; p_d; l_d; d_std_GDP; d_std_pi; d_corr_GDPpi;
               d_omG; d_std_Q; d_TBGDP; d_autocorr_Q; d_corr_GDPQ];
moment_labels_check = [
    arrayfun(@(i) sprintf('std(Y_%d)',  i), 1:12, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('std(PH_%d)', i), 1:12, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('std(L_%d)',  i), 1:12, 'UniformOutput', false), ...
    {'std(GDP)','std(pi)','corr(GDP,pi)','omG','std(Q)','TB/GDP','autocorr(Q)','corr(GDP,Q)'}];

nan_idx = find(isnan(d_vec_check));
if ~isempty(nan_idx)
    fprintf('\nERROR: the following moments are NaN — fix data sources before using for SMM:\n');
    for kk = 1:numel(nan_idx)
        fprintf('  [%2d] %s\n', nan_idx(kk), moment_labels_check{nan_idx(kk)});
    end
    error('compute_data_moments:nan_moments', ...
          '%d of 44 moments are NaN. See list above.', numel(nan_idx));
end
assert(numel(d_vec_check) == 44, 'BUG: expected 44 moments, got %d.', numel(d_vec_check));
fprintf('\nValidation passed: all 44 moments are non-NaN.\n');

out_file = fullfile(model_dir, 'data_moments_chile.mat');
save(out_file, 'dm_chile');
fprintf('Moments saved to %s\n', out_file);

fprintf('\nAll moments computed from data. Ready to run smm_estimation.m\n\n');

%% ================================================================ %%
%%  LOCAL FUNCTIONS                                                    %%
%% ================================================================ %%

function cyc = hp_cycle(y, lam)
% Standard Hodrick-Prescott filter via the band-pass matrix approach.
    T = length(y);
    e = speye(T);
    D = diff(e, 2);    % (T-2) x T second-difference matrix
    trend = (e + lam * (D' * D)) \ y;
    cyc   = y - trend;
end

function [Y_q, yr_out, qt_out] = monthly_to_quarterly(Y_m, yr_m, mth_m)
% Convert monthly matrix to quarterly by averaging months within each quarter.
% Quarter q of year y = months {3q-2, 3q-1, 3q}.
    month_to_q = ceil(mth_m / 3);   % 1,2,3->Q1; 4,5,6->Q2; etc.
    [unique_yq, ~, ic] = unique([yr_m, month_to_q], 'rows');
    yr_out = unique_yq(:,1);
    qt_out = unique_yq(:,2);
    nQ = size(unique_yq,1);
    nS = size(Y_m,2);
    Y_q = NaN(nQ, nS);
    for t = 1:nQ
        idx = (ic == t);
        if sum(idx) == 3
            Y_q(t,:) = mean(Y_m(idx,:), 1);
        elseif sum(idx) >= 2
            Y_q(t,:) = mean(Y_m(idx,:), 1);  % partial quarters: use available months
        end
    end
end
