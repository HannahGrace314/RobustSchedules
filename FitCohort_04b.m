%% Fit one parameter set per mouse using all available ultrasound points
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% OUTPUTS:
%   parametersets_eachmouse - nParam x nCohorts x nMice
%   T_all, X_all            - simulated time and states for the full-data fit
%   convergence, Rsquared, MAE

function [parametersets_eachmouse, T_all, X_all, convergence, Rsquared, MAE] = ...
    FitCohort_04b(ODEFunc, fixed, indexchooseparam, paramrange, ...
    initialcondition_all, initialday_all, time_end, ...
    exp_individual, times_all, cohort_dose_ot1, cohort_dose_gem, ...
    dose_times, digtwin_sched)

if nargin < 13
    digtwin_sched = [];
end

[nCohorts, ~, nMice_per_cohort] = size(exp_individual);

parametersets_eachmouse = NaN(length(fixed), nCohorts, nMice_per_cohort);
convergence = NaN(nCohorts, nMice_per_cohort);
Rsquared = NaN(nCohorts, nMice_per_cohort);
MAE = NaN(nCohorts, nMice_per_cohort);
T_all = cell(nCohorts, nMice_per_cohort);
X_all = cell(nCohorts, nMice_per_cohort);

for cohort = 1:nCohorts
    for mouse = 1:nMice_per_cohort
        non_nan_indices = find(~isnan(squeeze(exp_individual(cohort, :, mouse))));
        if isempty(non_nan_indices)
            continue
        end
        if exp_individual(cohort, non_nan_indices(1), mouse) == 0
            non_nan_indices = non_nan_indices(2:end);
        end

        timerange = [initialday_all(cohort, mouse), time_end];
        ic = initialcondition_all(cohort, mouse, :);
        dose = MouseDose_04c(cohort, mouse, cohort_dose_ot1, cohort_dose_gem, ...
            digtwin_sched, nCohorts);

        estimate_output = ParameterFit_04d(ODEFunc, fixed, indexchooseparam, paramrange, ...
            ic, timerange, exp_individual(cohort, non_nan_indices, mouse), ...
            times_all(cohort, non_nan_indices, mouse), dose, dose_times);

        parametersets_eachmouse(:, cohort, mouse) = estimate_output{1};
        convergence(cohort, mouse) = estimate_output{2};
        Rsquared(cohort, mouse) = estimate_output{3};
        MAE(cohort, mouse) = estimate_output{4};

        numsim = Evaluate_NumSim_01c(ODEFunc, estimate_output{1}, ic, ...
            timerange, dose, dose_times);
        T_all{cohort, mouse} = numsim{2};
        X_all{cohort, mouse} = numsim{1};
    end
end

end
