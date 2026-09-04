%% Fit one parameter set to each mouse on the original experimental schedules
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Fits ntc, rcm, and kgm to ultrasound trajectories from the
% four original arms (untreated, gemcitabine, OT-1, gemcitabine+OT-1).
% Gemcitabine parameters are restored to the nominal values for mice that
% never received gemcitabine. Optional save writes the flattened parameter
% table and the arrays needed by PlotTumorFits_04g.
%
% Set save_output = 1 to write fit_file for Treatment_Sweep_05.

save_output = 0;
fit_file = 'estimate_cohort_fits.mat';

BladderData_OriginalSchedules_02 % exp_individual, times_all, initialcondition_all, initialday_all

Figure_Aesthetics_03


%% Choose parameters to fit

chooseparam = [0;       % 1  - pc
               0;       % 2  - cmax
               0;       % 3  - ktc
               0;       % 4  - kgc
               1;       % 5  - ntc
               0;       % 6  - smt
               0;       % 7  - kgt
               0;       % 8  - dt
               0;       % 9  - T0
               1;       % 10 - rcm
               1;       % 11 - kgm
               0;       % 12 - stm
               0;       % 13 - dm
               0;       % 14 - M0
               0;       % 15 - dg
               0];      % 16 - km

indexchooseparam = find(chooseparam == 1);

%% Starting parameter guess
BladderParam_01b
paramrange = mpr(indexchooseparam,:);
fixed =   param;

time_end = 25;
ODEFunc = @BladderFunc_01;


%% Gradient descent fitting

[parametersets_eachmouse, T_all, X_all, convergence, Rsquared, MAE] = ...
    FitCohort_04b(ODEFunc, fixed, indexchooseparam, paramrange, ...
    initialcondition_all, initialday_all, time_end, ...
    exp_individual, times_all, cohort_dose_ot1, cohort_dose_gem, dose_times, []);


%% For untx and ot-1 treated mice, restore gemcitabine parameters
gem_param = [4,7,11]; 
no_gem_cohorts = [1, 3]; % untreated, OT-1
for c = no_gem_cohorts
    for mouse = 1:nmpc(c)
        if ~any(isnan(parametersets_eachmouse(:,c,mouse)))
            parametersets_eachmouse(gem_param,c,mouse) = param(gem_param);
        end
    end
end

% Flatten to mice x parameters in cohort order (untreated, Gem, OT-1, combo)
allexpmice_parameter_sets = nan(sum(nmpc), length(chooseparam));
k = 1;
for cohort = 1:nCohorts
    for mouse = 1:nmpc(cohort)
        allexpmice_parameter_sets(k,:) = parametersets_eachmouse(:,cohort,mouse)';
        k = k + 1;
    end
end

if save_output
    save(fit_file, 'allexpmice_parameter_sets', 'parametersets_eachmouse', ...
        'exp_individual', 'times_all', 'T_all', 'X_all', 'Rsquared', 'MAE', ...
        'orig_cohort_order', 'orig_cohort_names', 'nmpc', 'nCohorts', '-v7.3');
    fprintf('Wrote %s\n', fit_file);
end


%% Tumor fits (black)

opts_fit = struct;
opts_fit.xlim = [5 25];
opts_fit.fs = fs - 2;
opts_fit.ax_fs = fs - 5;
opts_fit.title_fs = fs + 4;
opts_fit.lineWidth = lineWidth + 0.5;
opts_fit.tickWidth = tickWidth;
opts_fit.colors = zeros(nCohorts, 3);
opts_fit.ylabel_str = 'Total Tumor (mm^3)';

PlotTumorFits_04g(exp_individual, times_all, T_all, X_all, Rsquared, MAE, ...
    orig_cohort_order, orig_cohort_names, opts_fit)


%% Waterfall of R^2 and MAE

opts_wf = struct;
opts_wf.summary_cohorts = sort(orig_cohort_order);
opts_wf.r2_threshold = 0.5;
opts_wf.mae_threshold = 20;
opts_wf.mae_ylim = [0 35];
opts_wf.fs = fs;

PlotWaterfallFits_04h(Rsquared, MAE, orig_cohort_order, orig_cohort_names, opts_wf)
