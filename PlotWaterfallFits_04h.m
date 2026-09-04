%% Waterfall of per-mouse R^2 and MAE, plus summary statistics
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Optional fields in opts:
%   summary_cohorts - which cohort indices enter mean/median/mode (default: all)
%   r2_threshold, mae_threshold, mae_ylim, fs

function fig = PlotWaterfallFits_04h(Rsquared, MAE, cohort_order, cohort_names, opts)

if nargin < 5
    opts = struct;
end
if ~isfield(opts, 'summary_cohorts'), opts.summary_cohorts = sort(cohort_order); end
if ~isfield(opts, 'r2_threshold'), opts.r2_threshold = 0.5; end
if ~isfield(opts, 'mae_threshold'), opts.mae_threshold = round(max(MAE,[],'all')*0.5); end
if ~isfield(opts, 'mae_ylim'), opts.mae_ylim = [0 round(max(MAE,[],'all')*1.05)]; end
if ~isfield(opts, 'fs'), opts.fs = 16; end

nCohorts = size(Rsquared, 1);
y_R2 = cell(nCohorts, 1);
y_MAE = cell(nCohorts, 1);
idx_R2 = cell(nCohorts, 1);
idx_MAE = cell(nCohorts, 1);

for cohort = 1:nCohorts
    % Rsquared waterfall organization
    [keepr,idxr] = sort(Rsquared(cohort,~isnan(Rsquared(cohort,:))),'descend');
    y_R2{cohort} = keepr';
    idx_R2{cohort} = 1:length(idxr);

    % MAE waterfall organization
    pull_MAE = MAE(cohort,~isnan(MAE(cohort,:)));
    [keepm,idxm] = sort(pull_MAE(idxr),'descend');
    y_MAE{cohort} = keepm';
    idx_MAE{cohort} = idxm; 
end


fig = figure;
wfplot = tiledlayout(2, length(cohort_order) * 2 + 1, 'TileSpacing', 'compact');
barWidth = 0.8;

plot_metric_row(y_R2, idx_R2, cohort_order, cohort_names, opts, ...
    true, opts.r2_threshold, [0 1], [10 147 150] / 255, [187 62 3] / 255, barWidth);
plot_metric_row(y_MAE, idx_MAE, cohort_order, cohort_names, opts, ...
    false, opts.mae_threshold, opts.mae_ylim, [187 62 3] / 255, [10 147 150] / 255, barWidth);

xlabel(wfplot, 'Mouse #', 'FontSize', opts.fs + 4)
sgtitle(wfplot, 'Evaluation of Parameter Fitting', 'FontSize', opts.fs + 6)
set(gcf, 'Position', [2064 148 1547 731])

end


function plot_metric_row(y_error, mouse_idx, cohort_order, cohort_names, opts, ...
    isR2, threshold, yl, color_good, color_bad, barWidth)

nCohorts = length(cohort_order);
for j = 1:nCohorts + 1
    if j <= nCohorts
        nexttile([1 2])
        set(gca, 'FontSize', opts.fs + 1)
        cohort = cohort_order(j);
        y = y_error{cohort};
        labels = string(mouse_idx{cohort});
        if isR2
            subtitle(cohort_names{cohort}, 'FontSize', opts.fs + 4)
        end
        is_summary = false;
    else
        nexttile
        set(gca, 'FontSize', opts.fs + 1)
        collect_y = vertcat(y_error{opts.summary_cohorts});
        if isR2
            y = [mean(collect_y); median(collect_y); mode(round(collect_y, 1))];
            subtitle('Summary Stats', 'FontSize', opts.fs + 4)
        else
            y = [mean(collect_y); median(collect_y); mode(round(collect_y))];
        end
        labels = ["mean"; "median"; "mode"];
        is_summary = true;
    end

    hold on
    b = bar(1:length(y), y, 'EdgeColor', 'none', 'BarWidth', barWidth);
    b.FaceColor = 'flat';
    bar_colors = repmat(color_good, length(y), 1);
    bar_colors(y < threshold, :) = repmat(color_bad, sum(y < threshold), 1);
    if is_summary
        bar_colors = repmat([148 210 189] / 255, length(y), 1);
    end
    b.CData = bar_colors;
    yline(threshold, '--', 'Color', [187 62 3] / 255, 'LineWidth', 1.5)
    xticks(1:length(labels))
    xticklabels(labels)
    if j == 1
        if isR2
            ylabel('R^2 Value', 'FontSize', opts.fs + 4)
        else
            ylabel('Mean Absolute Error (MAE)', 'FontSize', opts.fs + 4)
        end
    else
        yticklabels({})
    end
    ylim(yl)
    xlim([0 length(y) + 1])
    hold off
end

end
