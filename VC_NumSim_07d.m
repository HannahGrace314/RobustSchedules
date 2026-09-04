%% Numerical simulation of the virtual cohort under predicted schedules
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Simulate each experimentally tested gemcitabine/OT-1
% schedule on every virtual mouse. For each mouse, the digital-twin
% schedule is the one that minimizes final + mean total tumor volume
% through day 23. Trajectories are plotted against ultrasound data for
% each treatment group.
%
% Uses: BladderData_NewSchedules_07.m, BladderFunc_01.m, Evaluate_NumSim_01c.m,
% Figure_Aesthetics_03.m, and vc_param and vc_initialcondition from
% virtualcohort.mat

load('virtualcohort.mat')

%% ODE evaluation

use_endpt_vol = 2;
BladderData_NewSchedules_07

nCohorts = 4; % experimentally tested schedules (excludes digital twin)
timerange = [6, 100];
numsim_times = timerange(1):0.1:timerange(2);

ODEFunc = @BladderFunc_01;
n = length(vc_param);

% cohort x cell type (C, T, M, total) x time x mouse
vc_numsim = nan(nCohorts + 1, 4, length(numsim_times), n);

for cohort = 1:nCohorts
    dose = [cohort_dose_ot1(cohort, :); cohort_dose_gem(cohort, :)];
    parfor mouse = 1:n
        pull_numsim = Evaluate_NumSim_01c(ODEFunc, vc_param(mouse, :), ...
            vc_initialcondition(mouse, :), timerange, dose, dose_times);
        vc_numsim(cohort, :, :, mouse) = [pull_numsim{1}, sum(pull_numsim{1}, 2)]';
    end
end

%% Digital-twin schedule for each virtual mouse

digitaltwin_reg = nan(n, 1);

for mouse = 1:n
    pull_mouseresults = vc_numsim(:, 4, :, mouse);
    cost_function = zeros(nCohorts, 1);
    for cohort = 1:nCohorts
        totaltumor = squeeze(pull_mouseresults(cohort, 1:find(numsim_times==23)));
        cost_function(cohort) = totaltumor(end) + mean(totaltumor);
    end
    [~, digitaltwin_reg(mouse)] = min(cost_function);
    vc_numsim(5, :, :, mouse) = vc_numsim(digitaltwin_reg(mouse), :, :, mouse);
end

%% Figure: virtual-cohort simulations versus ultrasound data

Figure_Aesthetics_03

figure
allgroups_virtualcohort = tiledlayout(4, nCohorts + 1, 'TileSpacing', 'compact');

for j = 1:nCohorts + 1
    cohort = cohort_order(j);

    for celltype = 1:4
        nexttile(tilenum(allgroups_virtualcohort, celltype, j))
        hold on

        for mouse = 1:n
            plot(numsim_times, squeeze(vc_numsim(cohort, celltype, :, mouse)), ...
                'Color', colors(cohort, :) * (mouse / n));
        end

        if strcmp(cell_names{celltype}, 'Total Tumor')
            for mouse = 1:nmpc(cohort)
                scatter(times_all(cohort, :, mouse), exp_individual(cohort, :, mouse), ...
                    markerSize, 'k', 'filled')
            end
        end

        if celltype == 1
            title(cohort_names{cohort}, 'FontWeight', 'bold', 'FontSize', fs - 1)
        end

        if cohort == 1
            ylabel(cell_names{celltype})
        else
            set(gca, 'YTickLabel', {})
        end

        ylim(cell_range(celltype, :))
        xlim([6, Timedata(end)])
        grid on
    end
end

title(allgroups_virtualcohort, 'Virtual Murine Cohort', 'FontSize', fs + 1)
xlabel(allgroups_virtualcohort, 'Time (days)', 'FontSize', fs - 1)
ylabel(allgroups_virtualcohort, 'Volume (mm^3)', 'FontSize', fs - 1)
set(gcf, 'Position', [0 300 1000 800], 'PaperPositionMode', 'auto')
