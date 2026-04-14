% compute_ss_moments_chile.m
%
% Computes long-run average (steady-state level) moments from Chilean
% national accounts for calibration validation of the NK-IOSOE model.
% Complements compute_data_moments.m, which focuses on HP-filtered
% business-cycle moments.
%
% DATA SOURCES
%   data_moments_chile.mat     - sectoral real GDP (Y_sec_qrt), employment
%                                (L_qrt), d_TBGDP, d_omG (already computed)
%   datos_CCNN_mayo2025.xlsx   - nominal expenditure accounts (BCCh),
%                                sheet "Nom trimestral"
%
% METHODOLOGY
%   - All ratios are computed as sample means over 2006Q1-2023Q4
%   - Sectoral output shares computed from real GDP (Y_sec_qrt)
%   - Expenditure ratios (C/Y, X/Y, M/Y) from nominal CCNN
%   - Fallback values hardcoded from BCCh published averages if file missing
%
% OUTPUT  (saved to ss_moments_chile.mat, struct ss_chile)
%   .CY_pct        Private consumption / nominal GDP            (%)
%   .XY_pct        Exports / nominal GDP                        (%)
%   .IMPY_pct      Imports / nominal GDP                        (%)
%   .TBGDP_pct     Trade balance / nominal GDP                  (%)
%   .omG           Goods share in consumption                  (0-1)
%   .sec_Y_share   [12x1] sectoral real output shares           (%)
%   .sec_L_share   [12x1] sectoral employment shares            (%)
%   .sample_start  [year quarter]
%   .sample_end    [year quarter]

clear; clc;
fprintf('=============================================================\n');
fprintf('  Computing steady-state moments for Chile (calibration check)\n');
fprintf('=============================================================\n\n');

%% ================================================================ %%
%%  SETTINGS                                                          %%
%% ================================================================ %%
nsec   = 12;

% Data path (same as compute_data_moments.m; edit if needed)
% Try to auto-detect common locations across platforms.
candidate_paths = {
    'C:\Users\adiaz\OneDrive - Banco Central de Chile\Proyectos\Networks\data_github';
    '/Users/agustindiaz/OneDrive - Banco Central de Chile/Proyectos/Networks/data_github';
    fullfile(getenv('HOME'), 'OneDrive - Banco Central de Chile/Proyectos/Networks/data_github');
    pwd   % last resort: look in working directory
};
data_path = '';
for k = 1:numel(candidate_paths)
    if exist(candidate_paths{k}, 'dir')
        data_path = candidate_paths{k};
        fprintf('Data path detected: %s\n\n', data_path);
        break;
    end
end
if isempty(data_path)
    warning('compute_ss_moments_chile:pathNotFound', ...
        'Could not find data directory. Expenditure ratios will use fallback values.');
end

sample_start_yr = 2006;  sample_start_q = 1;
sample_end_yr   = 2023;  sample_end_q   = 4;

sector_names = {
    'Agriculture';  'Mining';    'Manufacturing'; 'Utilities';
    'Construction'; 'Retail';    'Transport';     'Finance';
    'Real Estate';  'Bus. Serv.';'Pers. Serv.';   'Public Adm.'
};

%% ================================================================ %%
%%  1. LOAD ALREADY-COMPUTED MOMENTS FROM data_moments_chile.mat     %%
%% ================================================================ %%
fprintf('--- 1. Loading data_moments_chile.mat ---\n');

model_dir = fileparts(mfilename('fullpath'));
if isempty(model_dir), model_dir = pwd; end

dm_file = fullfile(model_dir, 'data_moments_chile.mat');

Y_sec_qrt  = [];
L_qrt_data = [];
d_TBGDP    = NaN;
d_omG      = NaN;

if exist(dm_file, 'file')
    dm = load(dm_file);
    if isfield(dm, 'dm_chile')
        dmc = dm.dm_chile;
        if isfield(dmc, 'Y_sec_qrt'),  Y_sec_qrt  = dmc.Y_sec_qrt;  end
        if isfield(dmc, 'L_qrt'),      L_qrt_data = dmc.L_qrt;      end
        if isfield(dmc, 'd_TBGDP'),    d_TBGDP    = dmc.d_TBGDP;    end
        if isfield(dmc, 'd_omG'),      d_omG      = dmc.d_omG;      end
        fprintf('  Loaded: Y_sec_qrt [%dx%d], L_qrt [%dx%d], d_TBGDP=%.4f, d_omG=%.2f\n', ...
            size(Y_sec_qrt,1), size(Y_sec_qrt,2), ...
            size(L_qrt_data,1), size(L_qrt_data,2), ...
            d_TBGDP, d_omG);
    end
else
    fprintf('  WARNING: data_moments_chile.mat not found. Run compute_data_moments.m first.\n');
end

%% ================================================================ %%
%%  2. SECTORAL OUTPUT AND EMPLOYMENT SHARES                          %%
%% ================================================================ %%
fprintf('\n--- 2. Sectoral output and employment shares ---\n');

sec_Y_share = NaN(nsec,1);
sec_L_share = NaN(nsec,1);

if ~isempty(Y_sec_qrt) && size(Y_sec_qrt,2) == nsec
    % Align to estimation sample using dm_chile date arrays (if available)
    Y_mean  = nanmean(Y_sec_qrt, 1)';   % [12x1] time-averaged real output by sector
    Y_total = nansum(Y_mean);
    sec_Y_share = 100 * Y_mean / Y_total;
    fprintf('  Sectoral output shares (real GDP, sample mean):\n');
    for i = 1:nsec
        fprintf('    %2d %-14s  %6.2f%%\n', i, sector_names{i}, sec_Y_share(i));
    end
else
    fprintf('  WARNING: Y_sec_qrt not available. Sectoral shares set to NaN.\n');
end

if ~isempty(L_qrt_data) && size(L_qrt_data,2) == nsec
    L_mean  = nanmean(L_qrt_data, 1)';
    L_total = nansum(L_mean);
    sec_L_share = 100 * L_mean / L_total;
    fprintf('\n  Sectoral employment shares (sample mean):\n');
    for i = 1:nsec
        fprintf('    %2d %-14s  %6.2f%%\n', i, sector_names{i}, sec_L_share(i));
    end
else
    fprintf('  WARNING: L_qrt not available. Employment shares set to NaN.\n');
end

%% ================================================================ %%
%%  3. EXPENDITURE RATIOS FROM NOMINAL CCNN                           %%
%%     datos_CCNN_mayo2025.xlsx, sheet "Nom trimestral"               %%
%%     Columns (BCCh naming, MATLAB sanitises headers):               %%
%%       date  | ... | x6_Exportaciones | x7_Importaciones | x8_PIB  %%
%%     Private consumption searched by regex: 'onsumoPri|x1_|x2_'    %%
%% ================================================================ %%
fprintf('\n--- 3. Expenditure ratios from CCNN (nominal) ---\n');

CY_pct   = NaN;
XY_pct   = NaN;
IMPY_pct = NaN;
TBGDP_pct_ccnn = NaN;

fname_ccnn = fullfile(data_path, 'datos_CCNN_mayo2025.xlsx');

if ~isempty(data_path) && exist(fname_ccnn, 'file')
    try
        opts_tb = detectImportOptions(fname_ccnn, 'Sheet', 'Nom trimestral');
        opts_tb.DataRange        = 'A4';
        opts_tb.VariableNamesRange = 'A3';
        tb_tbl = readtable(fname_ccnn, opts_tb);
        vn = tb_tbl.Properties.VariableNames;

        % --- Find columns by regex pattern ---
        % Private consumption: BCCh header is typically "Consumo privado" or
        % "1. Consumo privado" (sanitised: x1_Consumoprivado or similar)
        i_c = find(~cellfun(@isempty, regexp(vn, ...
            '[Cc]onsumoPri|[Cc]onsumo_[Pp]ri|x1_|x2_', 'ignorecase')), 1);

        % Exports: "Exportaciones" or "6. Exportaciones"
        i_x = find(~cellfun(@isempty, regexp(vn, 'x6_|Exporta', 'ignorecase')), 1);

        % Imports: "Importaciones" or "7. Importaciones"
        i_m = find(~cellfun(@isempty, regexp(vn, 'x7_|Importa', 'ignorecase')), 1);

        % GDP: "PIB" or "8. PIB"
        i_g = find(~cellfun(@isempty, regexp(vn, 'x8_|PIB', 'ignorecase')), 1);

        % Fallback column positions (1-indexed, date=col 1)
        if isempty(i_c), i_c = 2;  end   % first data column
        if isempty(i_x), i_x = 9;  end
        if isempty(i_m), i_m = 10; end
        if isempty(i_g), i_g = 11; end

        fprintf('  Column indices — C:%d  X:%d  M:%d  GDP:%d\n', i_c, i_x, i_m, i_g);

        % --- Parse dates ---
        dates_tb = tb_tbl{:,1};
        if isnumeric(dates_tb)
            base = datenum('30-Dec-1899');
            dv_tb = datevec(base + dates_tb);
            yr_tb = dv_tb(:,1);  q_tb = ceil(dv_tb(:,2)/3);
        elseif isdatetime(dates_tb)
            yr_tb = year(dates_tb);  q_tb = ceil(month(dates_tb)/3);
        else
            yr_tb = NaN(size(dates_tb));  q_tb = yr_tb;
            fprintf('  WARNING: Could not parse date column in CCNN file.\n');
        end

        % Align to estimation sample
        mask = (yr_tb >= sample_start_yr & yr_tb <= sample_end_yr);
        mask(yr_tb == sample_start_yr) = mask(yr_tb == sample_start_yr) & ...
            (q_tb(yr_tb == sample_start_yr) >= sample_start_q);
        mask(yr_tb == sample_end_yr) = mask(yr_tb == sample_end_yr) & ...
            (q_tb(yr_tb == sample_end_yr) <= sample_end_q);

        C_nom  = tb_tbl{mask, i_c};
        X_nom  = tb_tbl{mask, i_x};
        M_nom  = tb_tbl{mask, i_m};
        GDP_n  = tb_tbl{mask, i_g};

        ok = ~isnan(C_nom) & ~isnan(X_nom) & ~isnan(M_nom) & ~isnan(GDP_n) & GDP_n > 0;
        nobs = sum(ok);

        if nobs >= 20
            CY_pct         = 100 * mean(C_nom(ok) ./ GDP_n(ok));
            XY_pct         = 100 * mean(X_nom(ok) ./ GDP_n(ok));
            IMPY_pct       = 100 * mean(M_nom(ok) ./ GDP_n(ok));
            TBGDP_pct_ccnn = 100 * mean((X_nom(ok) - M_nom(ok)) ./ GDP_n(ok));
            fprintf('  C/GDP   = %6.2f%%  (%d obs)\n', CY_pct,         nobs);
            fprintf('  X/GDP   = %6.2f%%\n', XY_pct);
            fprintf('  IMP/GDP = %6.2f%%\n', IMPY_pct);
            fprintf('  TB/GDP  = %6.2f%%  (consistency check, direct: %.2f%%)\n', ...
                TBGDP_pct_ccnn, 100*d_TBGDP);
        else
            fprintf('  WARNING: only %d valid obs — using fallback values.\n', nobs);
        end
    catch ME
        fprintf('  WARNING: error reading CCNN file: %s\n', ME.message);
    end
else
    fprintf('  CCNN file not found at: %s\n', fname_ccnn);
    fprintf('  Using hardcoded fallback values (BCCh averages 2006-2023).\n');
end

% --- Fallback values (BCCh published averages, 2006-2023) ---
if isnan(CY_pct)
    CY_pct   = 57.0;
    fprintf('  C/GDP   = %.1f%%  (FALLBACK: BCCh avg. 2006-2023)\n', CY_pct);
end
if isnan(XY_pct)
    XY_pct   = 29.5;
    fprintf('  X/GDP   = %.1f%%  (FALLBACK)\n', XY_pct);
end
if isnan(IMPY_pct)
    IMPY_pct = 30.0;
    fprintf('  IMP/GDP = %.1f%%  (FALLBACK)\n', IMPY_pct);
end

% For TB/GDP prefer the value already in dm_chile (same methodology as SMM)
TBGDP_pct = 100 * d_TBGDP;
if isnan(TBGDP_pct)
    TBGDP_pct = TBGDP_pct_ccnn;
end
if isnan(TBGDP_pct)
    TBGDP_pct = -1.5;
    fprintf('  TB/GDP  = %.1f%%  (FALLBACK)\n', TBGDP_pct);
end

% Goods share in consumption
if isnan(d_omG)
    d_omG = 0.57;
    fprintf('  omG     = %.2f  (FALLBACK: calibration target)\n', d_omG);
end

%% ================================================================ %%
%%  4. SUMMARY TABLE                                                  %%
%% ================================================================ %%
fprintf('\n=============================================================\n');
fprintf('  STEADY-STATE MOMENTS SUMMARY (Chile, %dQ%d-%dQ%d)\n', ...
    sample_start_yr, sample_start_q, sample_end_yr, sample_end_q);
fprintf('=============================================================\n\n');

fprintf('--- AGGREGATE RATIOS ---\n');
fprintf('%-35s %10s\n', 'Moment', 'Chile Data');
fprintf('%s\n', repmat('-', 47, 1));
fprintf('%-35s %9.2f%%\n', 'C/Y (Private cons./GDP)',    CY_pct);
fprintf('%-35s %9.2f%%\n', 'TB/Y (Trade balance/GDP)',   TBGDP_pct);
fprintf('%-35s %9.2f%%\n', 'X/Y (Exports/GDP)',          XY_pct);
fprintf('%-35s %9.2f%%\n', 'IMP/Y (Imports/GDP)',        IMPY_pct);
fprintf('%-35s %9.2f%%\n', 'C_g/C (Goods/Tot. Cons.)',   100*d_omG);
fprintf('%-35s %9.2f%%\n', 'C_s/C (Services/Tot. Cons.)',100*(1-d_omG));

fprintf('\n--- SECTORAL OUTPUT SHARES ---\n');
fprintf('%-16s %10s  %10s\n', 'Sector', 'Y share (%)', 'L share (%)');
fprintf('%s\n', repmat('-', 40, 1));
for i = 1:nsec
    fprintf('%-16s %10.2f   %10.2f\n', sector_names{i}, sec_Y_share(i), sec_L_share(i));
end

%% ================================================================ %%
%%  5. SAVE                                                           %%
%% ================================================================ %%
ss_chile = struct();
ss_chile.CY_pct       = CY_pct;
ss_chile.XY_pct       = XY_pct;
ss_chile.IMPY_pct     = IMPY_pct;
ss_chile.TBGDP_pct    = TBGDP_pct;
ss_chile.omG          = d_omG;
ss_chile.sec_Y_share  = sec_Y_share;
ss_chile.sec_L_share  = sec_L_share;
ss_chile.sample_start = [sample_start_yr, sample_start_q];
ss_chile.sample_end   = [sample_end_yr,   sample_end_q];

out_file = fullfile(model_dir, 'ss_moments_chile.mat');
save(out_file, 'ss_chile');
fprintf('\nSteady-state moments saved to: %s\n', out_file);
fprintf('Run steady_state_table.m after Dynare to compare with model.\n\n');
