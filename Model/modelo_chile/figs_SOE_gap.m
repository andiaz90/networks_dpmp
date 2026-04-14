% figs_SOE_gap.m
% Plotting script for main_SOE.m - Shows DEVIATIONS FROM STEADY STATE
% This model DOES NOT have flexible-price allocations (Y_f), so plots show
% deviations from steady state (Y_ss), not output gaps relative to flexible price.

fprintf('Generating plots: Steady State Deviations...\n');

% Try to get necessary variables
if ~exist('Ymat', 'var') || ~exist('name_vec', 'var')
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
            if ~exist('Ymat', 'var') && isfield(S, 'Ymat'), Ymat = S.Ymat; end
            if ~exist('name_vec', 'var') && isfield(S, 'name_vec'), name_vec = S.name_vec; end
        catch ME
            warning('figs_SOE_gap:failed_to_load','%s', ME.message);
        end
    end
end

% Basic checks
if ~exist('Ymat', 'var') || isempty(Ymat)
    warning('figs_SOE_gap: no ''Ymat'' found. Nothing to plot.');
    return;
end

nPeriods = size(Ymat, 1);
nSectors = size(Ymat, 2);

% If name_vec is missing or wrong size, create generic names
if ~exist('name_vec', 'var') || numel(name_vec) ~= nSectors
    name_vec = strings(nSectors,1);
    for ii = 1:nSectors
        name_vec(ii) = sprintf('Sector_%d', ii);
    end
end

% Extract aggregate output deviation from steady state if available
if exist('oo_', 'var') && exist('M_', 'var')
    try
        Y_idx = strmatch('Y', M_.endo_names, 'exact');
        if ~isempty(Y_idx)
            % Extract using the correct orientation
            if size(oo_.endo_simul, 1) == M_.endo_nbr
                Y = oo_.endo_simul(Y_idx, :);
            else
                Y = oo_.endo_simul(:, Y_idx)';
            end
            
            % Calculate deviation from steady state
            Y_ss = oo_.steady_state(Y_idx);
            Y_dev = 100 * (log(Y) - log(Y_ss));
            
            % Plot aggregate output deviation from SS
            fig_agg = figure('Name', 'Aggregate Output Deviation from SS', 'NumberTitle', 'off', 'Visible', 'on');
            drawnow;
            plot(1:length(Y_dev), Y_dev, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
            xlabel('Period'); ylabel('Deviation (%)');
            title('Aggregate Output: Deviation from Steady State');
            grid on;
            yline(0, '--k', 'LineWidth', 1);

            % Save figure
            print(fig_agg, '-dpng', '-r300', 'aggregate_output_gap.png');
            fprintf('✓ Figure saved: aggregate_output_gap.png\n');
        else
            warning('figs_SOE_gap: Y variable not found in model');
        end
    catch ME
        warning('figs_SOE_gap:aggregate_gap_failed','%s', ME.message);
    end
end

% Extract sectoral output deviations from steady state if available
sectoral_gaps_found = false;
Y_dev_mat = zeros(nPeriods, nSectors);
if exist('oo_', 'var') && exist('M_', 'var')
    for ii = 1:nSectors
        try
            Y_i_idx = strmatch(sprintf('Y_%d', ii), M_.endo_names, 'exact');
            if ~isempty(Y_i_idx)
                if size(oo_.endo_simul, 1) == M_.endo_nbr
                    Y_i = oo_.endo_simul(Y_i_idx, :)';
                else
                    Y_i = oo_.endo_simul(:, Y_i_idx);
                end
                
                % Calculate deviation from steady state
                Y_i_ss = oo_.steady_state(Y_i_idx);
                Y_dev_mat(:, ii) = 100 * (log(Y_i) - log(Y_i_ss));
                sectoral_gaps_found = true;
            end
        catch
            % Skip this sector
        end
    end
end

if sectoral_gaps_found
    % Plot sectoral output deviations from steady state
    fig_gaps = figure('Name', 'Sectoral Output Deviations from SS', 'NumberTitle', 'off', 'Visible', 'on');
    drawnow;
    plot(1:nPeriods, Y_dev_mat, 'LineWidth', 1.5);
    xlabel('Period'); ylabel('Deviation (%)');
    title('Sectoral Output: Deviation from Steady State');
    if nSectors <= 12
        legend(cellstr(name_vec), 'Location', 'eastoutside');
    end
    grid on;
    yline(0, '--k', 'LineWidth', 0.5);

    % Save figure
    print(fig_gaps, '-dpng', '-r300', 'sectoral_output_gaps.png');
    fprintf('✓ Figure saved: sectoral_output_gaps.png\n');
else
    fprintf('figs_SOE_gap: Using alternative extraction method for Y_%d variables.\n');
    
    % Fallback: Plot deviations from steady state if primary extraction unavailable
    if exist('Ymat', 'var') && ~isempty(Ymat)
        try
            Ymat_pct = zeros(size(Ymat));
            if exist('oo_', 'var') && exist('M_', 'var')
                for ii = 1:size(Ymat, 2)
                    Y_idx = strmatch(sprintf('Y_%d', ii), M_.endo_names, 'exact');
                    if ~isempty(Y_idx)
                        Y_ss = oo_.steady_state(Y_idx);
                        Ymat_pct(:, ii) = 100 * (Ymat(:, ii) - Y_ss) / Y_ss;
                    else
                        Ymat_pct(:, ii) = 100 * (Ymat(:, ii) - Ymat(1, ii)) / Ymat(1, ii);
                    end
                end
            else
                for ii = 1:size(Ymat, 2)
                    Ymat_pct(:, ii) = 100 * (Ymat(:, ii) - Ymat(1, ii)) / Ymat(1, ii);
                end
            end
            
            fig_fallback = figure('Name', 'Sectoral Output (SS Deviation)', 'NumberTitle', 'off', 'Visible', 'on');
            drawnow;
            plot(1:size(Ymat,1), Ymat_pct, 'LineWidth', 1.5);
            xlabel('Period'); ylabel('% Deviation from SS');
            title('Sectoral Output: Deviation from Steady State (Fallback)');
            if size(Ymat,2) <= 12
                legend(cellstr(name_vec), 'Location', 'eastoutside');
            end
            grid on;

            print(fig_fallback, '-dpng', '-r300', 'sectoral_output_ss_deviation_fallback.png');
            fprintf('✓ Figure saved: sectoral_output_ss_deviation_fallback.png\n');
        catch ME
            warning('figs_SOE_gap:fallback_failed','%s', ME.message);
        end
    end
end

fprintf('✓ Steady State Deviation plots completed\n');
