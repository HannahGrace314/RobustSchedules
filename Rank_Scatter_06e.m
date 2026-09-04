%% Scatter of which schedule wins a given rank across parameter space
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Each grid point is colored by the schedule that attains which_rank there.
% If several schedules share that rank for a mouse, pick the one that wins
% that rank most often: use cohort win-rates inside the virtual-cohort hull,
% and grid win-rates outside it. The hull is outlined in grey.
%
% Inputs
%   vc_param             - virtual-cohort parameters (for the hull)
%   parameterspace       - full grid of parameter pairs being plotted
%   rank_matrix_space    - (grid parameter space x schedule) rank of each schedule
%   rank_percent_space   - (schedule x rank) % of grid parameter space with that rank
%   rank_percent_cohort  - same, but for the virtual cohort only
%   distinguisher        - ranking metric name (figure title / fprintf)
%   sched_plot_colors      - RGB per schedule
%   labels               - legend strings per schedule
%   which_rank           - ranks to plot (one tile each; e.g. 1 or 1:3)
%   p1, p2               - column indices of the two parameters on the axes
%   pn                   - parameter names used as axis labels
%
% OUTPUT:
%   Scatter plot showing the robustness of schedules across the parameter
%   space and virtual cohort

function g = Rank_Scatter_06e(vc_param, parameterspace, rank_matrix_space, ...
    rank_percent_space, rank_percent_cohort, distinguisher, sched_plot_colors, labels, ...
    which_rank, p1, p2, pn)

    Figure_Aesthetics_03
    n_grid = size(parameterspace, 1);
    n_show = numel(which_rank);
    
    % Per grid mouse: all tied winning schedules, and the one we plot
    cell_collect_sched = cell(n_grid, max(which_rank)); % has multiple schedule per grid mouse per rank if schedules are equally as robust
    collect_sched = zeros(n_grid, max(which_rank)); % even if multiple schedules equally robust, only lists schedule that is most robust in both parameter space and virtual cohort (used for plotting)
    outline_color = [0.7, 0.7, 0.7];

    % Grey outline of the virtual cohort in this 2-D parameter plane
    hull_idx = convhull(vc_param(:, p1), vc_param(:, p2));
    in = inpolygon(parameterspace(:, p1), parameterspace(:, p2), ...
        vc_param(hull_idx, p1), vc_param(hull_idx, p2));

    % Fill collect_sched for every rank before assigning markers or plotting
    for k = 1:n_show
        rnk = which_rank(k);
        for mouse = 1:n_grid
            opt_sched_for_mouse = find(rank_matrix_space(mouse, :) == rnk);
            % Inside the hull, prefer schedules that also win often in the
            % cohort; outside, prefer schedules that win often on the grid
            if in(mouse)
                rank_percent = rank_percent_cohort;
            else
                rank_percent = rank_percent_space;
            end
            most_robust_percent = max(rank_percent(opt_sched_for_mouse, rnk));
            most_robust_sched = find(rank_percent(:, rnk) == most_robust_percent);
            cell_collect_sched{mouse, rnk} = most_robust_sched;
            collect_sched(mouse, rnk) = intersect(most_robust_sched, opt_sched_for_mouse);
        end
    end

    % One marker per unique schedule across all tiles
    marker_symbols = ["square", "o", "^", "x", "diamond", "+", "v", "_", "*", ...
        ">", "<", "pentagram", "hexagram", "|", "."];
    get_sched = unique(collect_sched);
    get_sched = get_sched(get_sched ~= 0);
    nSched = size(rank_matrix_space, 2);
    sched_plot_markers = strings(nSched, 1);
    for ms = 1:numel(get_sched)
        sched_plot_markers(get_sched(ms)) = marker_symbols(ms);
    end

    g = figure;
    if n_show < 4
        scatterplot = tiledlayout(1, n_show, 'TileSpacing', 'compact');
    else
        n_row = ceil(sqrt(n_show));
        n_col = ceil(n_show / n_row);
        scatterplot = tiledlayout(n_row, n_col, 'TileSpacing', 'compact');
    end

    for k = 1:n_show
        rnk = which_rank(k);
        ax = nexttile(k);
        hold on

        h_outline = plot(vc_param(hull_idx, p1), vc_param(hull_idx, p2), ...
            'LineStyle', '-', 'Color', outline_color, 'LineWidth', lineWidth + 0.5);

        % Make these markers filled so easier to see
        filledMarkers = ["o", "diamond", "pentagram", "hexagram", "^", "v"]; 
        for mouse = 1:n_grid
            most_robust = collect_sched(mouse, rnk);
            if ismember(sched_plot_markers(most_robust), filledMarkers)
                scatter(parameterspace(mouse, p1), parameterspace(mouse, p2), 80, ...
                    sched_plot_colors(most_robust, :), sched_plot_markers(most_robust), 'filled')
            else
                scatter(parameterspace(mouse, p1), parameterspace(mouse, p2), 80, ...
                    sched_plot_colors(most_robust, :), sched_plot_markers(most_robust), 'LineWidth', 2)
            end
        end

        % Legend on the last tile only (color + marker per plotted schedule)
        if k == n_show
            unique_collect_sched = unique(collect_sched);
            leg_sched = nonzeros(unique_collect_sched);
            h = gobjects(numel(leg_sched), 1);
            for j = 1:numel(leg_sched)
                if ismember(sched_plot_markers(leg_sched(j)), filledMarkers)
                    h(j) = scatter(nan, nan, 40, sched_plot_colors(leg_sched(j), :), ...
                        sched_plot_markers(leg_sched(j)), 'filled');
                else
                    h(j) = scatter(nan, nan, 40, sched_plot_colors(leg_sched(j), :), ...
                        sched_plot_markers(leg_sched(j)), 'LineWidth', 2);
                end
            end
            if iscell(labels)
                sched_lab = labels(leg_sched);
            else
                sched_lab = cellstr(labels(leg_sched));
            end
            leg = legend(ax, [h;h_outline], [sched_lab(:); {'Virtual cohort'}], ...
                'FontSize', fs - 4);
            leg.Layout.Tile = 'east';
            leg.Title.String = 'Schedules';
        end

        set(gca, 'FontSize', fs - 4, 'LineWidth', tickWidth)
        xlabel(pn(p1), 'FontSize', fs - 2)
        ylabel(pn(p2), 'FontSize', fs - 2)
        subtitle(sprintf('Rank %d', rnk), 'FontSize', fs)

        % Schedules that win this rank somewhere inside the hull
        all_more_robust = [];
        for mouse = 1:n_grid
            if in(mouse)
                all_more_robust = [all_more_robust; cell_collect_sched{mouse, rnk}]; 
            end
        end
        best = unique(all_more_robust);
        fprintf('%s: rank %d: schedules ranking the highest and the most robust were: %s\n', ...
            distinguisher, rnk, mat2str(best'));
    end

    title(scatterplot, sprintf('Parameter Space ranked by %s', distinguisher), ...
        'FontSize', fs + 2)

end
