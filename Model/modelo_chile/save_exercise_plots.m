function save_exercise_plots(output_dir, exercise_name)
% SAVE_EXERCISE_PLOTS Saves all open figures to the specified directory
%
% Inputs:
%   output_dir - Directory path where plots will be saved
%   exercise_name - Name prefix for the saved plots

fprintf('\nSaving plots to %s...\n', output_dir);

% Get all figure handles
fig_handles = findall(0, 'Type', 'figure');

if isempty(fig_handles)
    warning('No figures found to save.');
    return;
end

fprintf('Found %d figures to save\n', length(fig_handles));

for i = 1:length(fig_handles)
    fig = fig_handles(i);
    
    % Get figure number
    fig_num = get(fig, 'Number');
    
    % Try to get a meaningful name from figure title
    fig_name = get(fig, 'Name');
    if isempty(fig_name)
        % If no name, try to get from title
        ax = findobj(fig, 'Type', 'axes');
        if ~isempty(ax)
            title_obj = get(ax(1), 'Title');
            if ~isempty(title_obj)
                fig_name = get(title_obj, 'String');
                if iscell(fig_name)
                    fig_name = strjoin(fig_name, '_');
                end
            end
        end
    end
    
    % Clean name for filename
    if isempty(fig_name)
        clean_name = sprintf('figure_%02d', fig_num);
    else
        clean_name = regexprep(fig_name, '[^\w\s-]', '');
        clean_name = regexprep(clean_name, '\s+', '_');
        clean_name = sprintf('fig%02d_%s', fig_num, clean_name);
    end
    
    % Limit filename length
    if length(clean_name) > 60
        clean_name = clean_name(1:60);
    end
    
    % Full filename with exercise prefix
    filename_base = sprintf('%s_%s', exercise_name, clean_name);
    
    % Save as both PNG and FIG for flexibility
    png_path = fullfile(output_dir, [filename_base '.png']);
    fig_path = fullfile(output_dir, [filename_base '.fig']);
    
    try
        % Save as PNG (high resolution)
        saveas(fig, png_path);
        fprintf('  Saved: %s\n', [filename_base '.png']);
        
        % Save as FIG (editable MATLAB figure)
        savefig(fig, fig_path);
        fprintf('  Saved: %s\n', [filename_base '.fig']);
    catch ME
        warning('Could not save figure %d: %s', fig_num, ME.message);
    end
end

fprintf('All plots saved successfully!\n\n');

end
