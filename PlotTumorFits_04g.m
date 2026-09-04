%% Plot tumor fits, mice ranked by R^2 within each cohort
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% One tile per mouse. Columns follow cohort_order; within a column, row 1
% is the best R^2. Optional opts fields (all have defaults):
%   mouse_colors, schedule, xlim, ylim, legend_colors, legend_labels,
%   legend_tile, fs, ax_fs, title_fs, lineWidth, tickWidth

function fig = PlotTumorFits_04g(exp_individual, times_all, T_all, X_all, ...
    Rsquared, MAE, cohort_order, cohort_names, opts)

if nargin < 9 || isempty(opts)
    opts = struct;
end

[nCohorts, ~, nmpc] = size(exp_individual);

if ~isfield(opts, 'mouse_colors'),   opts.mouse_colors = []; end
if ~isfield(opts, 'schedule'),       opts.schedule = []; end
if ~isfield(opts, 'xlim'),           opts.xlim = []; end
if ~isfield(opts, 'ylim'),           opts.ylim = []; end
if ~isfield(opts, 'colors'),         opts.colors = zeros(nCohorts, 3); end
if ~isfield(opts, 'legend_colors'),  opts.legend_colors = []; end
if ~isfield(opts, 'legend_labels'),  opts.legend_labels = {}; end
if ~isfield(opts, 'legend_tile'),    opts.legend_tile = [nmpc nCohorts-1]; end
if ~isfield(opts, 'fs'),             opts.fs = 14; end
if ~isfield(opts, 'ax_fs'),          opts.ax_fs = opts.fs - 2; end
if ~isfield(opts, 'title_fs'),       opts.title_fs = opts.fs + 5; end
if ~isfield(opts, 'lineWidth'),      opts.lineWidth = 2; end
if ~isfield(opts, 'tickWidth'),      opts.tickWidth = 0.8; end


fig = figure;
numsim = tiledlayout(nmpc, nCohorts, 'TileSpacing', 'compact');

for j = 1:nCohorts
    cohort = cohort_order(j);

    % Rank mice in this cohort by R^2 (NaNs dropped, then best = 1)
    r2 = Rsquared(cohort, :);
    r2 = r2(~isnan(r2));
    [~, sortIdx] = sort(r2, 'descend');
    rank_r2 = zeros(size(r2));
    rank_r2(sortIdx) = 1:length(r2);

    for mouse = 1:length(rank_r2)
        i = rank_r2(mouse); % row in this column (1 = best R^2)
        data = exp_individual(cohort, :, mouse);
        times = times_all(cohort, :, mouse);
        non_nan_indices = find(~isnan(data));
        if isempty(non_nan_indices)
            continue
        end

        % Color: per-mouse RGB if supplied, otherwise this cohort's color
        if ~isempty(opts.mouse_colors)
            col = squeeze(opts.mouse_colors(cohort, mouse, :))';
        else
            col = opts.colors(cohort, :);
        end

        nexttile(tilenum(numsim, i, j));
        hold on
        plot(T_all{cohort, mouse}, sum(X_all{cohort, mouse}, 2), ...
            'Color', col, 'LineWidth', opts.lineWidth)
        scatter(times(non_nan_indices), data(non_nan_indices), [], col, 'filled')

        xlab = append('Mouse ', num2str(i), ...
            ': [R^2 = ', num2str(round(Rsquared(cohort, mouse), 2)), ...
            ', MAE = ', num2str(round(MAE(cohort, mouse), 1)), ']');
        if ~isempty(opts.schedule)
            xlab = append(xlab, ', Schedule ', num2str(opts.schedule(cohort, mouse)));
        end
        xlabel(xlab)
        set(gca, 'FontSize', opts.ax_fs, 'LineWidth', opts.tickWidth)

        % Limits: numeric [min max], 'data' (round up from observations), or skip
        if strcmp(opts.xlim, 'data')
            xlim([0 ceil(max(times(non_nan_indices)) / 50) * 50])
        elseif ~isempty(opts.xlim)
            xlim(opts.xlim)
        end
        if strcmp(opts.ylim, 'data')
            ylim([0 ceil(max(data(non_nan_indices)) / 100) * 100])
        elseif ~isempty(opts.ylim)
            ylim(opts.ylim)
        end

        yt = yticks;
        if numel(yt) > 2
            yticks([yt(1), yt(end)]);
        end

        if i == 1
            title(cohort_names{cohort}, 'FontSize', opts.title_fs)
        end
    end
end

% Schedule legend
if ~isempty(opts.legend_colors)
    ax = nexttile(tilenum(numsim, opts.legend_tile(1), opts.legend_tile(2)));
    axis(ax, 'off')
    h = gobjects(size(opts.legend_colors, 1), 1);
    for k = 1:size(opts.legend_colors, 1)
        h(k) = patch(NaN, NaN, opts.legend_colors(k, :));
    end
    lgd = legend(h, opts.legend_labels, 'Orientation', 'horizontal', ...
        'FontSize', opts.fs + 3);
    lgd.Title.String = 'Schedule';
end

ylabel(numsim, 'Total Tumor (mm^3)', 'FontSize', opts.fs + 4)
xlabel(numsim, 'Days', 'FontSize', opts.fs + 4)
set(gcf, 'Position', [1661 -457 1259 1336], 'PaperPositionMode', 'auto')

end
