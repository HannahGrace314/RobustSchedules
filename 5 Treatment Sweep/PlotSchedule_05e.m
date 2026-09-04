%% Plot unique selected schedules as a colored schedule grid
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Displays unique gemcitabine/OT-1 schedules after ordering
% with OrderSchedule_05f. Row color identifies the schedule; squares mark
% gemcitabine and circles mark OT-1. Y-tick labels are mouse counts.
%
% INPUTS:
%   order        - 1: popularity, 2: dose pattern (passed to OrderSchedule_05f)
%   S            - nMice x nDays matrix of selected schedules
%   tf           - final simulation time
%   ofn          - objective functional index (for the figure title)
%   savefigures  - 1 to export a 300 dpi JPEG, otherwise 0
%
% OUTPUT:
%   plot_sched - figure handle

function plot_sched = PlotSchedule_05e(order, S, tf, ofn, savefigures)

[nmcr_order, uniqueS, idx, ~, rowColors] = OrderSchedule_05f(S, tf, order);

plot_sched = figure('Units', 'inches', 'Position', [1 1 9 6]);
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
ax = nexttile;

[nr, nc] = size(uniqueS(idx, 1:tf-9));

% Make colorful grid based on schedules
RGB = ones(nr, nc, 3);
for r = 1:nr
    for c = 1:nc
        val = uniqueS(idx(r), c);
        if val == 1 || val == 2
            RGB(r, c, :) = rowColors(r, :);
        end
    end
end

image(RGB);
axis equal tight;

% Row and Column Labels
labelFontSize = 16;
tickdays = [10:7:tf, tf];
set(ax, 'XTick', tickdays - 9, 'YTick', 1:nr, 'TickLength', [0 0], ...
    'FontSize', labelFontSize, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left');
xticklabels(tickdays);
yticklabels(nmcr_order);

xlabel('Days');
ylabel('Number of Mice Preferred');

hold on

% Markers based on Gem or OT-1 
markerSize = 9;
for r = 1:nr
    for c = 1:nc
        val = uniqueS(idx(r), c);
        if val == 1
            plot(c, r, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
                'MarkerSize', markerSize);
        elseif val == 2
            plot(c, r, 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', [0.3 0.3 0.3], ...
                'MarkerSize', markerSize);
        end
    end
end

% Gridlines
for c = 0:nc
    plot([c + 0.5, c + 0.5], [0.5, nr + 0.5], 'k', 'LineWidth', 1);
end
for r = 0:nr
    plot([0.5, nc + 0.5], [r + 0.5, r + 0.5], 'k', 'LineWidth', 1);
end

% Legend
legendSize = 12;
hGem = plot(nan, nan, 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerSize', legendSize);
hOT1 = plot(nan, nan, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
    'MarkerSize', legendSize);
lgd = legend([hGem, hOT1], {'Gemcitabine', 'OT-1'}, ...
    'Location', 'southoutside', 'Orientation', 'vertical');
set(lgd, 'FontSize', 16, 'Box', 'off');

% Title
obj_name = ["final tumor size", "average tumor size", "minimum tumor size", ...
    "maximum tumor size", "average + final tumor", "average + 2 final tumor"];
title(sprintf('Schedules evaluated by %s', obj_name(ofn)));

hold off
drawnow limitrate

% Save figure
if savefigures == 1
    dateStr = datestr(now, 'ddmmmyyyy');
    filename = ['ts_sched_ofn', num2str(ofn), '_MaxWeeks', num2str((tf - 9) / 7), ...
        '_', dateStr, '.jpg'];
    exportgraphics(plot_sched, filename, 'Resolution', 300);
end

end
