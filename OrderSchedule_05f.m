%% Order unique schedules and assign row colors
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Collapses per-mouse best schedules to unique rows and orders
% them either by how many mice selected each schedule or by dose pattern.
%
% INPUTS:
%   S     - nMice x nDays matrix of selected schedules
%           (0 = none, 1 = OT-1, 2 = gemcitabine)
%   tf    - final simulation time (treatment days are columns 1:tf-9)
%   order - 1: sort by popularity (most mice first)
%           2: sort by number of OT-1 doses, then by the dose sequence
%
% OUTPUTS:
%   nmcs_order    - number of mice associated with each unique row
%   uniqueS       - unique schedule rows (original order from unique)
%   idx           - row permutation of uniqueS used for plotting
%   idx_uniqueS   - map from each mouse (row of R) to a unique schedule
%   rowColors     - RGB colors for the plotted unique rows, in idx order

function [nmcs_order, uniqueS, idx, idx_uniqueS, rowColors] = OrderSchedule_05f(S, tf, order)
Figure_Aesthetics_03

[uniqueS, ~, idx_uniqueS] = unique(S, 'rows', 'stable');
num_mice_choose_sched = zeros(max(idx_uniqueS),1);
for i = 1:max(idx_uniqueS)
    num_mice_choose_sched(i) = sum(idx_uniqueS == i);
end

switch order
    case 1
        [nmcs_order, idx] = sort(num_mice_choose_sched, 'descend');
        ns = size(uniqueS(idx, 1:tf-9), 1);
        rowColors = expand_colormap(palette, ns);

    case 2
        n_u = size(uniqueS, 1);
        S_label = zeros(n_u, 1);
        num_ot1 = zeros(n_u, 1);
        for i = 1:n_u
            vals = nonzeros(uniqueS(i, :));
            % Sort key matches the old coding (OT-1 = 2, gem = 1).
            S_label(i) = str2double(sprintf('%d', 3 - vals));
            num_ot1(i) = sum(vals == 1);
        end

        [~, idx] = sortrows([num_ot1, S_label], [1 2]);

        nmcs_order = num_mice_choose_sched(idx);

        sched_label_order = S_label(idx);
        [unique_label_order, unique_label_idx] = unique(sched_label_order);
        colors = expand_colormap(palette, length(unique_label_order));

        rowColors = zeros(length(idx), 3);
        uli_range = [sort(unique_label_idx); length(idx) + 1];
        for i = 1:length(uli_range) - 1
            for k = uli_range(i):uli_range(i + 1) - 1
                rowColors(k, :) = colors(i, :);
            end
        end

    otherwise
        error('OrderSchedule_05f:InvalidOrder', ...
            'order must be 1 (popularity) or 2 (dose pattern).');
end

end


function pal = expand_colormap(colors, n_rows)
% Stretch a short palette to n_rows colors.
%
% If n_rows is no larger than the palette, return the palette as-is.
% Otherwise keep every original swatch and insert extra colors in the
% gaps between neighbors, spreading those extras as evenly as possible.

n_swatches = size(colors, 1);
n_extra = n_rows - n_swatches;
if n_extra <= 0
    pal = colors;
    return
end

n_gaps = n_swatches - 1;
extras_per_gap = floor(n_extra / n_gaps) * ones(n_gaps, 1);
remainder = rem(n_extra, n_gaps);
extras_per_gap(1:remainder) = extras_per_gap(1:remainder) + 1;

pal = colors(1, :);
for i = 1:n_gaps
    % Start at swatch i (already in pal), add extras_per_gap(i) in-between
    % colors, then the next swatch.
    t = linspace(0, 1, extras_per_gap(i) + 2);
    segment = interp1([0 1], [colors(i, :); colors(i + 1, :)], t);
    pal = [pal; segment(2:end, :)]; 
end

end


