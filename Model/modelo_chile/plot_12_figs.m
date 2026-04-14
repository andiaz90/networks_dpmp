% plot_12_figs.m
% Compatibility wrapper: some versions call plot_12_figs as a script.
% This wrapper delegates to figs_SOE.m if available.
try
    if exist('figs_SOE.m','file')
        figs_SOE;
    elseif exist('plot_12_figs.m','file')
        % Shouldn't happen (this file), but safe-guard
        warning('plot_12_figs: self_reference','plot_12_figs present as wrapper.');
    else
        error('plot_12_figs or figs_SOE not found in folder.');
    end
catch ME
    warning('plot_12_figs:failed','Could not run plot_12_figs compatibility wrapper: %s', ME.message);
end
