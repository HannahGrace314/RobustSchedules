%% Experimental tumor growth and survival on the new schedules
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: 
% Loads BladderData_NewSchedules_07 and plots mean ± SEM tumor volume, 
% individual trajectories, and Kaplan-Meier survival. A mouse's death day is 
% the earlier of first volume >= endpoint_vol and the recorded dayofdeath; 
% otherwise it is censored at the last ultrasound day. 
%
% INPUTS (from BladderData_NewSchedules_07 / Figure_Aesthetics_03)
%   exp_individual, Timedata, nmpc, endpoint_vol
%   km_day_all, km_censored_all
%   colors, cohort_names, cohort_order, fs, lineWidth, markerSize
%
% OUTPUTS
%   mean ± SEM tumor-volume figure 
%   individual tumor-volume figure 
%   Kaplan-Meier figure 

Figure_Aesthetics_03

%% Load data (endpoint volume filled in, for mean curves)

use_endpt_vol = 1;
BladderData_NewSchedules_07

[nCohorts, nTime, nMiceMax] = size(exp_individual);

nCohorts_plot = nCohorts - 1; % removing digital twin cohort from plots

meanTumor = mean(exp_individual, 3, 'omitnan');
nMiceAtTime = sum(~isnan(exp_individual), 3);
stdTumor = std(exp_individual, 0, 3, 'omitnan');
semTumor = stdTumor ./ sqrt(nMiceAtTime);


%% Mean tumor volume ± SEM

lastday = 42;

end_idx = find(Timedata >= lastday, 1, 'first');
if isempty(end_idx)
    end_idx = nTime;
end

figure
hold on
for c = 1:nCohorts_plot
    errorbar(Timedata(1:end_idx), meanTumor(c, 1:end_idx), semTumor(c, 1:end_idx), ...
        '-o', 'Color', colors(c, :), 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, 'MarkerFaceColor', colors(c, :));
end
xlim([Timedata(1), Timedata(end_idx)])
ylim([0, endpoint_vol])
grid on
set(gca, 'FontSize', fs - 2)
xlabel('Days', 'FontSize', fs)
ylabel('Mean Tumor Volume (mm^3)', 'FontSize', fs)
set(gcf, 'Position', [2118 277 709 241])

%% Individual tumor growth (measured volumes, not replaced at endpoint)

use_endpt_vol = 2;
BladderData_NewSchedules_07

figure
tlo = tiledlayout(1, nCohorts_plot, 'TileSpacing', 'compact', 'Padding', 'compact');

for d = 1:nCohorts_plot
    c = cohort_order(d);
    ax = nexttile;
    hold(ax, 'on')
    for m = 1:nMiceMax
        mouseData = squeeze(exp_individual(c, :, m));
        if all(isnan(mouseData))
            continue
        end
        if c == 5
            plot(ax, Timedata, mouseData, '-', 'Color', colors(digtwin_sched(m), :), 'LineWidth', 1);
        else
            plot(ax, Timedata, mouseData, '-', 'Color', colors(c, :), 'LineWidth', 1);
        end
    end
    title(ax, cohort_names{c}, 'FontSize', fs - 2)
    xlim(ax, [0, Timedata(end)])
    ylim(ax, [0, endpoint_vol])
    grid(ax, 'on')
    ax.FontSize = fs - 2;
    if d > 1
        ax.YTickLabel = [];
    end
    if d == nCohorts_plot
        legendHandles = gobjects(nCohorts_plot, 1);
        for k = 1:nCohorts_plot
            ck = cohort_order(k);
            legendHandles(k) = plot(ax, nan, nan, 'Color', colors(ck, :), 'LineWidth', 2);
        end
        legend(ax, legendHandles, cohort_names(cohort_order(1:nCohorts_plot)), ...
            'FontSize', fs, 'Location', 'eastoutside');
    end
end

xlabel(tlo, 'Days', 'FontSize', fs)
ylabel(tlo, 'Tumor Volume (mm^3)', 'FontSize', fs)
set(gcf, 'Position', [2366 92 1028 245])

%% Kaplan-Meier survival (death = min(endpoint day, dayofdeath))

figure
hold on
for c = 1:nCohorts_plot
    t_m = km_day_all(c, :)';
    cens_m = km_censored_all(c, :)';
    keep = isfinite(t_m);
    [x, y] = KaplanMeier_07c(t_m(keep), cens_m(keep), Timedata(end));
    plot(x, 100 * y, '-', 'Color', colors(c, :), 'LineWidth', 2);
end
set(gca, 'FontSize', fs - 2)
xlabel('Days', 'FontSize', fs)
ylabel('% Survival', 'FontSize', fs)
xlim([0, Timedata(end)])
ylim([0, 100])
grid on
set(gcf, 'Position', [1659 277 458 234])
