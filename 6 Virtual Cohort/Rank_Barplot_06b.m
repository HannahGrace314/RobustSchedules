%% Stacked bar: percent of the virtual cohort at each rank, by schedule
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% INPUTS:
%   rank_percent    - rows are Schedule numbers and columns are the rank
%   numbers (1 = best, ...)
%   distinguisher   - cost function used to rank the schedules
%   sched_plot_colors - color associated with each schedule
%   labels          - label associated with each schedule
%   sched_for_leg     - schedule label for legend
%   numranks        - greatest number of ranks to plot
%
% OUTPUT:
%   Bar plot figure showing the robustness of schedules across multiple 
%   (or all) ranks for that distinguisher

function f = Rank_Barplot_06b(rank_percent, distinguisher, sched_plot_colors, labels, ...
    sched_for_leg, numranks)

    Figure_Aesthetics_03
    rank_percent = rank_percent(:, 1:numranks);
    nSched = size(rank_percent, 1);
    
    % Group schedules that share a color so stacked segments sit together
    [~, uc_idx] = unique(sched_plot_colors, 'rows', 'stable');
    sort_idx = [];
    for i = 1:length(uc_idx)
        check_color = sched_plot_colors(uc_idx(i), :);
        find_idx = find(all(sched_plot_colors == check_color, 2));
        sort_idx = [sort_idx; find_idx]; 
    end
    sort_idx = flip(sort_idx);
    
    % Removes ranks that no schedule was listed at (hence unnecessary)
    sum_rank_percent = sum(rank_percent, 1);
    find_relevant_ranks = find(sum_rank_percent ~= 0);
    rank_percent = rank_percent(:, find_relevant_ranks);
    
    f = figure;
    tiledlayout(1, 1, 'Padding', 'tight', 'TileSpacing', 'tight');
    nexttile(1)
    h = bar(rank_percent(sort_idx, :)', 'stacked');
    plot_colors = sched_plot_colors(sort_idx, :);
    for i = 1:nSched
        h(i).FaceColor = plot_colors(i, :);
    end
    
    unsort_idx(sort_idx) = 1:length(sort_idx); 
    
    xlabel('Rank')
    ylabel('Cumulative Percentage of Mice')
    title(sprintf('Ranking Virtual Cohort by %s', distinguisher))
    legend(h(unsort_idx(sched_for_leg)), labels{sched_for_leg}, 'Location', 'eastoutside') 
 
    xtick_vals = find_relevant_ranks;
    xticks(xtick_vals)
    xticklabels(string(xtick_vals))
    set(gca, 'FontSize', fs + 2, 'LineWidth', tickWidth)
    set(gcf, 'Position', [136 420 1192 446])
    grid on

end
