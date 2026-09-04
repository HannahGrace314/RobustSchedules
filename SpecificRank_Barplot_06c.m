%% Bar plot of one rank: percent of mice assigning that rank to each schedule
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
%   specificrank    - rank number
%
% OUTPUT:
%   Bar plot figure ordering schedules according to their robustness for
%   that rank and distinguisher
%   Any schedules that are >80% robust for that rank and distinguisher are
%   mentioned in the Command Window

function f = SpecificRank_Barplot_06c(rank_percent, distinguisher, sched_plot_colors, ...
    labels, specificrank)

    Figure_Aesthetics_03
    
    sched_xticks = find(rank_percent(:, specificrank) > 0)';
    [~, sort_size] = sort(rank_percent(sched_xticks, specificrank), 'descend');
    sched_xticks = sched_xticks(sort_size); % largest bars on the left
    
    f = figure;
    tiledlayout(1, 1, 'Padding', 'tight', 'TileSpacing', 'tight');
    nexttile(1)
    hold on
    above80 = [];
    for i = 1:length(sched_xticks)
        j = sched_xticks(i);
        bar(i, rank_percent(j, specificrank), 'FaceColor', sched_plot_colors(j, :));
        if rank_percent(j, specificrank) >= 80
            above80 = [above80, j];
        end
    end
    
    xlabel('Schedule')
    ylabel('Percentage of Mice')
    title(sprintf('Rank %d for %s', specificrank, distinguisher), 'FontSize', fs + 2)
    xticks(1:length(sched_xticks))
    xticklabels(labels(sched_xticks))
    ylim([0 100])
    set(gca, 'FontSize', fs, 'LineWidth', tickWidth)
    set(gcf, 'Position', [1870 -131 1292 293])
    grid on
    
    fprintf('%s: schedules at rank %d for 80%% or more of the cohort: %s\n', ...
        distinguisher, specificrank, mat2str(above80));

end
