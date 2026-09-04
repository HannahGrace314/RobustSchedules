%% Pie chart of one rank, with schedule numbers on slices and a order legend
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: One pie slice per schedule that attains specificrank for
% any mouse. Slice size is the percent of mice; the label is "S#: p%".
% The legend collapses slices that share a schedule order-style label 
% (GOGG, Original,...) and lists those groups from largest to smallest share.
%
% INPUTS:
%   rank_percent    - (schedule x rank) percent of mice assigning that rank
%   distinguisher   - cost function used to rank the schedules
%   sched_plot_colors - RGB color for each schedule (row-aligned with ranks)
%   labels          - schedule order string for each schedule
%   specificrank    - which rank column to plot (1 = best)
%
% OUTPUT:
%   Pie chart of how the cohort splits among schedules at that rank

function f = SpecificRank_Piechart_06d(rank_percent, distinguisher, sched_plot_colors, ...
    labels, specificrank)

    Figure_Aesthetics_03

    % larger-is-better for survival, smaller-is-better otherwise
    if strcmpi(distinguisher, "Overall Survival (OS)") || ...
            strcmpi(distinguisher, "Disease-Free Survival (DFS)")
        sort_mode = 'Largest';
    else
        sort_mode = 'Smallest';
    end

    % Schedules that actually win this rank for some mice, small → large
    sched_xticks = find(rank_percent(:, specificrank) > 0)';
    [~, sort_size] = sort(rank_percent(sched_xticks, specificrank), 'ascend');
    sched_xticks = sched_xticks(sort_size); 

    % Group slices that share a color (same schedule order), then reverse so
    % the largest groups are drawn first around the pie
    collect_colors = sched_plot_colors(sched_xticks, :);
    [~, uc_idx] = unique(collect_colors, 'rows', 'stable');
    sort_idx = [];
    for i = 1:length(uc_idx)
        check_color = collect_colors(uc_idx(i), :);
        find_idx = find(all(collect_colors == check_color, 2));
        sort_idx = [sort_idx; find_idx];
    end
    sort_idx = flip(sort_idx); 

    order = sched_xticks(sort_idx);
    vals = rank_percent(order, specificrank);
    slice_colors = sched_plot_colors(order, :);
    u_labels = string(labels(order));

    % Whole percents for labels and legend; pie() still uses vals for size
    slice_pct = percents_sum_to_100(vals); 
    pie_text = strings(numel(order), 1);
    for i = 1:numel(order)
        if slice_pct(i) < 1
            pie_text(i) = "";
        else
            pie_text(i) = "S" + string(order(i)) + ": " + string(slice_pct(i)) + "%";
        end
    end

    f = figure;
    h = pie(vals, cellstr(pie_text));
    patches = h(1:2:end); 
    texts = h(2:2:end);
    set(texts, 'FontSize', fs)
    for i = 1:numel(patches) 
        patches(i).FaceColor = slice_colors(i, :);
    end

    % Legend
    % Labels: One legend row per unique schedule order string: sum the labeled slice percents
    [u_label, ~, ic] = unique(u_labels, 'stable');
    n_u = numel(u_label);
    order_pct = zeros(n_u, 1);
    order_col = zeros(n_u, 3);
    for i = 1:n_u
        mask = ic == i;
        order_pct(i) = sum(slice_pct(mask));
        order_col(i, :) = slice_colors(find(mask, 1), :);
    end

    keep = order_pct > 0; 
    u_label = u_label(keep); 
    order_col = order_col(keep, :); 
    order_pct = order_pct(keep);

    [order_pct, sort_leg] = sort(order_pct, 'descend'); 
    u_label = u_label(sort_leg);
    order_col = order_col(sort_leg, :);

    % Dummy patches 
    hold on
    ph = gobjects(numel(u_label), 1);
    for i = 1:numel(u_label)
        ph(i) = patch(NaN, NaN, order_col(i, :), 'EdgeColor', 'none');
    end
    leg_str = u_label + " (" + string(order_pct) + "%)";
    lgd = legend(ph, leg_str, 'Location', 'eastoutside'); 
    lgd.Title.String = 'Schedules';
    set(lgd, 'FontSize', fs - 1, 'Box', 'off')

    title(sprintf('%s %s', sort_mode, distinguisher), 'FontSize', fs + 2)
    set(gcf, 'Position', [2602 242 720 420])

end


function pct = percents_sum_to_100(raw)
    % Integer percents that sum to 100 (largest-remainder method).
    floored = floor(raw);
    frac = raw - floored;
    need = 100 - sum(floored);
    if need > 0
        [~, idx] = sort(frac, 'descend');
        floored(idx(1:need)) = floored(idx(1:need)) + 1;
    end
    pct = floored;
end
