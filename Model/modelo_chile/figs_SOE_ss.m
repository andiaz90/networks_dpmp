% figs_SOE_ss.m
% Plotting script for main_SOE.m - Shows DEVIATIONS FROM STEADY STATE
% This model does NOT have flexible-price allocations (Y_f), so plots show
% deviations from steady state, not output gaps.

fprintf('Generating plots: Deviations from Steady State...\n');

% Try to get necessary variables: pmat (prices matrix), Ymat, name_vec
if ~exist('pmat', 'var') || ~exist('Ymat', 'var') || ~exist('name_vec', 'var')
    if exist('model_output_IOSOE.mat', 'file') || exist('model_output_IOSOE_ex1.mat', 'file') || ...
       exist('model_output_IOSOE_ex2.mat', 'file') || exist('model_output_IOSOE_ex3.mat', 'file')
        try
            % Try to load the most recent output file
            if exist('EXERCISE', 'var')
                if EXERCISE == 1
                    S = load('model_output_IOSOE_ex1.mat');
                elseif EXERCISE == 2
                    S = load('model_output_IOSOE_ex2.mat');
                elseif EXERCISE == 3
                    S = load('model_output_IOSOE_ex3.mat');
                else
                    S = load('model_output_IOSOE.mat');
                end
            else
                S = load('model_output_IOSOE.mat');
            end
            if ~exist('pmat', 'var') && isfield(S, 'pmat'), pmat = S.pmat; end
            if ~exist('Ymat', 'var') && isfield(S, 'Ymat'), Ymat = S.Ymat; end
            if ~exist('name_vec', 'var') && isfield(S, 'name_vec'), name_vec = S.name_vec; end
        catch ME
            warning('figs_SOE_ss:failed_to_load','%s', ME.message);
        end
    end
end

% Basic checks
if ~exist('pmat', 'var') || isempty(pmat)
    warning('figs_SOE_ss: no ''pmat'' found. Nothing to plot.');
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

% Plot sectoral prices
fig1 = figure('Name', 'Sectoral Price Deviations from SS', 'NumberTitle', 'off', 'Visible', 'on');
drawnow;
try
    plot(1:nPeriods, pmat_pct, 'LineWidth', 1.5);
    xlabel('Period'); ylabel('% Deviation from Steady State');
    title('Sectoral Prices: Deviation from Steady State');
    if nSectors <= 12
        legend(cellstr(name_vec), 'Location', 'eastoutside');
    end
    grid on;
catch ME
    warning('figs_SOE_ss:failed_to_plot_pmat','%s', ME.message);
end

% Plot sectoral output if available
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
        
        fig2 = figure('Name', 'Sectoral Output Deviations from SS', 'NumberTitle', 'off', 'Visible', 'on');
        drawnow;
        plot(1:size(Ymat,1), Ymat_pct, 'LineWidth', 1.5);
        xlabel('Period'); ylabel('% Deviation from Steady State');
        title('Sectoral Output: Deviation from Steady State');
        if size(Ymat,2) <= 12
            legend(cellstr(name_vec), 'Location', 'eastoutside');
        end
        grid on;
    catch ME
        warning('figs_SOE_ss:failed_to_plot_Ymat','%s', ME.message);
    end
end

% Save figures
try
    print(fig1, '-dpng', '-r300', 'sectoral_prices_ss_deviation.png');
    fprintf('✓ Figure saved: sectoral_prices_ss_deviation.png\n');
    if exist('Ymat','var')
        print(fig2, '-dpng', '-r300', 'sectoral_output_ss_deviation.png');
        fprintf('✓ Figure saved: sectoral_output_ss_deviation.png\n');
    end
catch ME
    warning('Failed to save figures: %s', ME.message);
end

fprintf('✓ Steady State Deviation plots completed\n');
