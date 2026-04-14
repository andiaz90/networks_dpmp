function sync_ylims(fig)
% SYNC_YLIMS  Give all subplots in a figure the same y-axis range.
%
%   sync_ylims(fig)   synchronize all axes in figure handle fig
%   sync_ylims()      uses gcf
%
% The shared range is [global_min, global_max] across all axes, expanded
% to enforce the minimum span set by enforce_min_yaxis (default 1e-2).
% Axes that are colorbars, legends, or have no visible lines are skipped.

if nargin < 1 || isempty(fig), fig = gcf; end

% Collect all axes children (excludes colorbars, legends, etc.)
ax_all = findobj(fig, 'Type', 'axes');
if isempty(ax_all), return; end

% Filter: keep only axes that contain at least one line/patch/area
ax_data = gobjects(0);
for k = 1:numel(ax_all)
    kids = get(ax_all(k), 'Children');
    has_data = any(arrayfun(@(h) ...
        ismember(get(h,'Type'), {'line','patch','area','bar','surface'}), kids));
    if has_data
        ax_data(end+1) = ax_all(k); %#ok<AGROW>
    end
end
if isempty(ax_data), return; end

% Find global [ymin, ymax] from current ylim of each data axis
%   (limits already set by the individual subplot code + enforce_min_yaxis)
ylo = arrayfun(@(a) ylim(a), ax_data, 'UniformOutput', false);
ylo = cell2mat(ylo);          % [n x 2]
global_min = min(ylo(:,1));
global_max = max(ylo(:,2));

% Enforce minimum span on the global range too
min_span = 1e-2;
if (global_max - global_min) < min_span
    mid = 0.5 * (global_min + global_max);
    global_min = mid - min_span/2;
    global_max = mid + min_span/2;
end

% Apply to all data axes
for k = 1:numel(ax_data)
    ylim(ax_data(k), [global_min, global_max]);
end
end
