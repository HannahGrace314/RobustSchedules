%% Plotting all Virtual Cohort Response Trajectories versus Original
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% INPUTS:
%   virtualcohort_response - 2x1 cell with simulation output X in {1} and
%                            time T in {2} for each mouse i under each schedule j
%   virtualcohort_analysis - used to grab the volume at the final time (tf)
%                            to calculate the mean and standard deviation
%   sched_plot_colors        - color associated with each schedule
%   tf                     - final time
%
% OUTPUT:
%   Trajectories of virtual cohort under each schedule compared to the
%   original schedule

function f = VC_Trajectories_06f(virtualcohort_response, virtualcohort_analysis, sched_plot_colors, tf)

    nSched = size(sched_plot_colors,1);
    v_n = size(virtualcohort_response{1},1);

    Figure_Aesthetics_03
    
    f = figure;
    numsim = tiledlayout(min(factor(nSched)), max(factor(nSched)), 'TileSpacing', 'compact');
    for s = 1:nSched
        nexttile(s)
        if s < nSched
            title(sprintf('Schedule %d', s), 'FontSize', fs)
        else
            title('Original', 'FontSize', fs)
        end
        hold on
        for i = 1:v_n
            pull_X = virtualcohort_response{1}{i, s};
            pull_T = virtualcohort_response{2}{i, s};
            idx = pull_T <= tf + 1;
            plot(pull_T(idx), sum(pull_X(idx, 1:min(size(pull_X, 2), 3)), 2), ...
                'Color', sched_plot_colors(s, :), 'LineWidth', lineWidth)
        end
        errorbar(tf, mean(virtualcohort_analysis{1, 1}(:, s)), ...
            std(virtualcohort_analysis{1, 1}(:, s)), 'o', ...
            'Color', [174, 32, 18] ./ 255, 'LineWidth', 2)
        xlim([0 tf + 1])
        ylim([0 250])
    end
    set(findall(gcf, 'Type', 'axes'), 'FontSize', fs - 2, 'LineWidth', tickWidth)
    ylabel(numsim, 'Tumor Volume (mm^3)', 'FontSize', fs + 2)
    xlabel(numsim, 'Days', 'FontSize', fs + 2)
    set(gcf, 'Position', [1819, 79, 1430, 611])
end