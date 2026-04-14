function enforce_min_yaxis(ax, min_span)
% ENFORCE_MIN_YAXIS  Prevent y-axis from collapsing to numerical-noise scale.
%
%   enforce_min_yaxis(ax)            uses default min_span = 1e-3
%   enforce_min_yaxis(ax, min_span)  uses caller-supplied threshold
%
% If the current y-axis range |ymax - ymin| < min_span, the limits are
% expanded symmetrically around the midpoint to ±(min_span/2).
% This stops MATLAB from showing 1e-16 axes when an IRF is effectively zero.

if nargin < 1 || isempty(ax), ax = gca; end
if nargin < 2 || isempty(min_span), min_span = 1e-2; end

yl = ylim(ax);
if (yl(2) - yl(1)) < min_span
    mid = 0.5 * (yl(1) + yl(2));
    ylim(ax, [mid - min_span/2, mid + min_span/2]);
end
end
