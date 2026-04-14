% plot_12_figs.m
% Safe plotting script for the 12-sector model.
% It is a script (not a function) so it shares the caller workspace when
% called from another script (like main_SOE.m). If variables are missing
% it will attempt to load 'model_output_chile.mat' and proceed if possible.

% Try to get necessary variables: pmat (prices matrix), Ymat, name_vec
if ~exist('pmat', 'var') || ~exist('Ymat', 'var') || ~exist('name_vec', 'var')
    if exist('model_output_chile.mat', 'file')
        try
            S = load('model_output_chile.mat');
            if ~exist('pmat', 'var') && isfield(S, 'pmat'), pmat = S.pmat; end
            if ~exist('Ymat', 'var') && isfield(S, 'Ymat'), Ymat = S.Ymat; end
            if ~exist('name_vec', 'var') && isfield(S, 'name_vec'), name_vec = S.name_vec; end
        catch ME
            warning('plot_12_figs:failed_to_load','%s', ME.message);
        end
    end
end

% Basic checks
if ~exist('pmat', 'var') || isempty(pmat)
    warning('plot_12_figs: no ''pmat'' found. Nothing to plot.');
    return;
end

nPeriods = size(pmat, 1);
nSectors = size(pmat, 2);

% If name_vec is missing or wrong size, create generic names
if ~exist('name_vec', 'var') || numel(name_vec) ~= nSectors
    name_vec = strings(nSectors,1);
    for ii = 1:nSectors
        name_vec(ii) = sprintf('Sector_%d', ii);
    end
end

% Convert prices to percent deviation from steady state
% Get steady states from oo_.steady_state if available
pmat_pct = zeros(size(pmat));
if exist('oo_', 'var') && exist('M_', 'var')
    for ii = 1:nSectors
        p_idx = strmatch(sprintf('P_%d', ii), M_.endo_names, 'exact');
        if ~isempty(p_idx)
            p_ss = oo_.steady_state(p_idx);
            pmat_pct(:, ii) = 100 * (pmat(:, ii) - p_ss) / p_ss;
        else
            % Fallback: use first period as steady state
            pmat_pct(:, ii) = 100 * (pmat(:, ii) - pmat(1, ii)) / pmat(1, ii);
        end
    end
else
    % Fallback: use first period as steady state
    for ii = 1:nSectors
        pmat_pct(:, ii) = 100 * (pmat(:, ii) - pmat(1, ii)) / pmat(1, ii);
    end
end

% Plot: for readability, plot first 12 sectors (or all if <=12)
fig1 = figure('Name', 'Sectoral price deviations', 'NumberTitle', 'off', 'Visible', 'on');
drawnow;
try
    plot(1:nPeriods, pmat_pct, 'LineWidth', 1.5);
    xlabel('Period'); ylabel('% deviation from SS');
    title('Sectoral Prices');
    if nSectors <= 12
        legend(cellstr(name_vec), 'Location', 'eastoutside');
    end
    grid on;
    catch ME
        warning('plot_12_figs:failed_to_plot_pmat','%s', ME.message);
end

% Additional optional plots if Ymat present
if exist('Ymat', 'var') && ~isempty(Ymat)
    try
        % Convert output to percent deviation from steady state
        Ymat_pct = zeros(size(Ymat));
        if exist('oo_', 'var') && exist('M_', 'var')
            for ii = 1:size(Ymat, 2)
                Y_idx = strmatch(sprintf('Y_%d', ii), M_.endo_names, 'exact');
                if ~isempty(Y_idx)
                    Y_ss = oo_.steady_state(Y_idx);
                    Ymat_pct(:, ii) = 100 * (Ymat(:, ii) - Y_ss) / Y_ss;
                else
                    % Fallback: use first period as steady state
                    Ymat_pct(:, ii) = 100 * (Ymat(:, ii) - Ymat(1, ii)) / Ymat(1, ii);
                end
            end
        else
            % Fallback: use first period as steady state
            for ii = 1:size(Ymat, 2)
                Ymat_pct(:, ii) = 100 * (Ymat(:, ii) - Ymat(1, ii)) / Ymat(1, ii);
            end
        end
        
        fig2 = figure('Name', 'Sectoral output deviations', 'NumberTitle', 'off', 'Visible', 'on');
        drawnow;
        plot(1:size(Ymat,1), Ymat_pct, 'LineWidth', 1.5);
        xlabel('Period'); ylabel('% deviation from SS');
        title('Sectoral Output');
        if size(Ymat,2) <= 12
            legend(cellstr(name_vec), 'Location', 'eastoutside');
        end
        grid on;
    catch ME
        warning('plot_12_figs:failed_to_plot_Ymat','%s', ME.message);
    end
end

% Save figures to files
try
    print(fig1, '-dpng', '-r300', 'sectoral_prices.png');
    fprintf('✓ Sectoral prices figure saved: sectoral_prices.png\n');
    if exist('Ymat','var')
        print(fig2, '-dpng', '-r300', 'sectoral_output.png');
        fprintf('✓ Sectoral output figure saved: sectoral_output.png\n');
    end
catch ME
    warning('Failed to save figures: %s', ME.message);
end
% end
