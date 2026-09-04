%% Sweep gemcitabine/OT-1 schedules on the original experimental mice
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Builds the feasible set of weekday gemcitabine and OT-1
% combinations, evaluates each schedule on every fitted experimental mouse,
% and stores the unique best schedule per mouse for each objective and
% treatment duration. An optional colored tumor-fit figure maps those
% chosen schedules onto the ultrasound fits from Estimate_Cohort_04.
%
% Requires estimate_cohort_fits.mat (or fit_file) from Estimate_Cohort_04.m.

%% User settings

max_weeks_tested = 2;
numDaysOapart    = 3;   % minimum days between OT-1 and a gemcitabine dose
dose_scale       = 0;   % 0 full labeled dose; 1 split across injections
ofn_list         = 1:6; % objective functional_list
%                           1: final tumor (Cf)
%                           2: average tumor
%                           3: minimum tumor
%                           4: maximum tumor
%                           5: final + average
%                           6: 2*final + average
schedule_order   = 2;   % 1 popularity; 2 dose pattern (see OrderSchedule_05f)

savefigures  = 0;       % 1 exports each PlotSchedule_05e figure (300 dpi JPEG)
save_output  = 0;       % 1 writes all_output to output_file
output_file  = 'treatment_sweep_all_output.mat';
fit_file     = 'estimate_cohort_fits.mat';

% Colored tumor-fit figure (schedule colors from all_output)
plot_ofn      = 5;
plot_MaxWeeks = 2;

Figure_Aesthetics_03
tickWidth = 0.8;


%% Candidate Schedules

days = max_weeks_tested * 7;
weeks = reshape(1:days, 7, [])';
numWeeks = size(weeks, 1);

G_schedules = gem_schedules(days, weeks, numWeeks);
O_schedules = ot1_schedules(days, numWeeks, numDaysOapart);
Schedules = combine_schedules(G_schedules, O_schedules, days, numDaysOapart);

nSched = size(Schedules, 1);
fprintf('%d gemcitabine patterns, %d OT-1 patterns, %d combination schedules\n', ...
    size(G_schedules, 1), size(O_schedules, 1), nSched);


%% Experimental mice (fits and initial conditions)

load(fit_file, 'allexpmice_parameter_sets');
originalmice = allexpmice_parameter_sets;
num_expmice = size(originalmice, 1);

BladderData_OriginalSchedules_02
[IC, Iday, FEday, FE_vol] = flatten_original_mice( ...
    initialcondition_all, initialday_all, finalday_all, finalultrasound_all, nmpc);

if size(IC, 1) ~= num_expmice
    error('Treatment_Sweep_05:MouseCount', ...
        ['Flattened initial conditions (%d) do not match fitted mice (%d). ', ...
         'Re-run Estimate_Cohort_04 and reload fit_file.'], size(IC, 1), num_expmice);
end

tf = 9 + days;


%% Evaluate every schedule on every mouse

origC = parallel.pool.Constant(originalmice);
IC_C  = parallel.pool.Constant(IC);
ID_C  = parallel.pool.Constant(Iday);
FeD_C = parallel.pool.Constant(FEday);

endofweek = 16:7:tf;
num_out_eval_sched = 4 * length(endofweek) + 4;
evaluate_sched = nan(nSched, num_expmice, num_out_eval_sched);

fprintf('Evaluating %d schedule-mouse pairs...\n', nSched * num_expmice);
parfor s = 1:nSched
    schedule = Schedules(s, :);
    eval_sched_local = zeros(num_expmice, num_out_eval_sched);
    for mouse = 1:num_expmice
        eval_sched_local(mouse, :) = EvaluateSchedule_05b(schedule, ...
            origC.Value(mouse, :), IC_C.Value(mouse, :), ...
            ID_C.Value(mouse), FeD_C.Value(mouse), tf, dose_scale);
    end
    evaluate_sched(s, :, :) = eval_sched_local;
end


%% Best schedule per mouse, by objective and duration

all_output = cell(max(ofn_list), max_weeks_tested);

for ofn = ofn_list
    for MaxWeeks = 1:max_weeks_tested
        MaxOT1 = 2;
        MaxGem = 2 * MaxWeeks;

        pull_sched = PullSchedule_05c(evaluate_sched, FE_vol, MaxOT1, MaxGem, MaxWeeks);
        idx_best_sched = BestSchedule_05d(pull_sched, ofn);

        S = Schedules(idx_best_sched, :);

        if ofn == plot_ofn && MaxWeeks == plot_MaxWeeks
            PlotSchedule_05e(schedule_order, S, 9 + (7 * MaxWeeks), ofn, savefigures);
        end

        [nmcs_order, uniqueS, idx, idx_uniqueS, rowColors] = ...
            OrderSchedule_05f(S, tf, schedule_order);
        orderedS = uniqueS(idx, :);

        if sum(nmcs_order) < num_expmice
            warning('Treatment_Sweep_05:MissingSchedule', ...
                'Some mouse did not get an optimal schedule (ofn = %d, MaxWeeks = %d).', ...
                ofn, MaxWeeks);
        end

        n_unique = max(idx_uniqueS);
        output = cell(n_unique, 5);
        for s = 1:n_unique
            mouse = find(idx_uniqueS == idx(s));
            pull_sched_row = orderedS(s, :);
            idx_npr = find(pull_sched_row ~= 0);
            nonzero_pull_sched = pull_sched_row(idx_npr);
            resituate_sched = zeros(length(nonzero_pull_sched), 2);
            resituate_sched(:, 1) = 1.426 * (nonzero_pull_sched == 1);
            resituate_sched(:, 2) = 3.8e4 * (nonzero_pull_sched == 2);

            output{s, 1}{1, 1} = resituate_sched;    % schedule
            output{s, 1}{1, 2} = (9 + idx_npr)';   % times
            output{s, 2} = rowColors(s, :); % colors for schedule 
            output{s, 3} = DoseOrderString_05g(resituate_sched, D1, D2); % labels
            output{s, 4} = {originalmice(mouse, :), mouse, IC(mouse, :), Iday(mouse)}; % mouse parameter sets, mouse numbers, initial conditions, and initial days that chose this sched


            NumOT1Taken = evaluate_sched(idx_best_sched(mouse(1)), mouse(1), end - 2);
            NumGemTaken = evaluate_sched(idx_best_sched(mouse(1)), mouse(1), end - 1);
            NumWeeksTaken = evaluate_sched(idx_best_sched(mouse(1)), mouse(1), end);

            output{s, 5}{1,1} = [NumOT1Taken, NumGemTaken, NumWeeksTaken, s]; % doses/weeks chosen
            output{s, 5}{1,2} = [MaxOT1, MaxGem, MaxWeeks, s]; % maximum doses/weeks allowed

        end
        all_output{ofn, MaxWeeks} = output;
    end
end

if save_output
    save(output_file, 'all_output', 'Schedules', 'evaluate_sched', ...
        'ofn_list', 'max_weeks_tested', 'dose_scale', '-v7.3');
    fprintf('Wrote %s\n', output_file);
end


%% Colored tumor fits for the chosen objective and duration 

S = load(fit_file);
need = {'exp_individual', 'times_all', 'T_all', 'X_all', 'Rsquared', 'MAE', ...
    'orig_cohort_order', 'orig_cohort_names', 'nmpc', 'nCohorts'};
if ~all(isfield(S, need))
    warning('Treatment_Sweep_05:SkipTumorPlot', ...
        ['%s is missing PlotTumorFits_04g inputs. Re-run Estimate_Cohort_04 ', ...
         'with save_output = 1.'], fit_file);
    return
end
if isempty(all_output{plot_ofn, plot_MaxWeeks})
    error('Treatment_Sweep_05:NoOutput', ...
        'all_output{%d,%d} is empty. Check plot_ofn and plot_MaxWeeks.', ...
        plot_ofn, plot_MaxWeeks);
end

% Colors and labels for legend
output = all_output{plot_ofn, plot_MaxWeeks};
keep_for_colors = nan(num_expmice, 4);
legend_colors = [];
legend_labels = {};
for s = 1:size(output, 1)
    if isempty(output{s, 4})
        continue
    end
    mice = output{s, 4}{2};
    color = output{s, 2};
    keep_for_colors(mice, 1) = s;
    keep_for_colors(mice, 2:4) = repmat(color, numel(mice), 1);

    if isempty(legend_colors) || ~ismember(color, legend_colors, 'rows')
        legend_colors = [legend_colors; color]; 
        legend_labels{end + 1} = output{s, 3};
    end
end

% Colors for each mouse
mouse_colors = nan(S.nCohorts, max(S.nmpc), 3);
schedule_id = nan(S.nCohorts, max(S.nmpc));
k = 1;
for cohort = 1:S.nCohorts
    for mouse = 1:S.nmpc(cohort)
        schedule_id(cohort, mouse) = keep_for_colors(k, 1);
        mouse_colors(cohort, mouse, :) = keep_for_colors(k, 2:4);
        k = k + 1;
    end
end

opts_fit = struct;
opts_fit.mouse_colors = mouse_colors;
opts_fit.schedule = schedule_id;
opts_fit.xlim = [5 25];
opts_fit.legend_colors = legend_colors;
opts_fit.legend_labels = legend_labels;
opts_fit.legend_tile = [15 3];
opts_fit.fs = fs - 2;
opts_fit.ax_fs = fs - 5;
opts_fit.title_fs = fs + 4;
opts_fit.lineWidth = lineWidth + 0.5;
opts_fit.tickWidth = tickWidth;
opts_fit.ylabel_str = 'Total Tumor (mm^3)';

PlotTumorFits_04g(S.exp_individual, S.times_all, S.T_all, S.X_all, ...
    S.Rsquared, S.MAE, S.orig_cohort_order, S.orig_cohort_names, opts_fit);

if savefigures
    dateStr = datestr(now, 'ddmmmyyyy');
    filename = sprintf('ts_tumorfit_ofn%d_MaxWeeks%d_%s.jpg', ...
        plot_ofn, plot_MaxWeeks, dateStr);
    exportgraphics(gcf, filename, 'Resolution', 300);
end


%% Local functions

function G_schedules = gem_schedules(days, weeks, numWeeks)
G_schedules = [];
    for mask = 0:(2^days - 1)
        G = bitget(mask, 1:days);
        valid = true;
        for w = 1:numWeeks
            if sum(G(weeks(w, :))) > 2
                valid = false;
                break
            end
        end
        if valid
            G_schedules = [G_schedules; G];
        end
    end
    G_schedules = clear_weekends(G_schedules, numWeeks);
    G_schedules = unique(G_schedules, 'rows', 'stable');
end


function O_schedules = ot1_schedules(days, numWeeks, numDaysOapart)
    O_schedules = zeros(1, days);
    for d = 1:days
        O = zeros(1, days);
        O(d) = 1;
        O_schedules = [O_schedules; O];
    end
    for d1 = 1:days
        for d2 = d1 + numDaysOapart:days
            O = zeros(1, days);
            O([d1 d2]) = 1;
            O_schedules = [O_schedules; O];
        end
    end
    O_schedules = clear_weekends(O_schedules, numWeeks);
    O_schedules = unique(O_schedules, 'rows', 'stable');
end


function schedules = clear_weekends(schedules, numWeeks)
    for week = 1:numWeeks
        schedules(:, (week * 7) - 1) = 0;
        schedules(:, week * 7) = 0;
    end
end


function Schedules = combine_schedules(G_schedules, O_schedules, days, numDaysOapart)
    Schedules = [];
    for i = 1:size(G_schedules, 1)
        G = G_schedules(i, :);
        for j = 1:size(O_schedules, 1)
            O = O_schedules(j, :);
            valid = true;
            Odays = find(O);
            for d = Odays
                block = max(1, d - numDaysOapart):min(days, d + numDaysOapart); % blocked days because within numDaysOapart of O dose
                if any(G(block))
                    valid = false;
                    break
                end
            end
            if valid
                S = zeros(1, days);
                S(G == 1) = 2;
                S(O == 1) = 1;
                Schedules(end + 1, :) = S; 
            end
        end
    end
    Schedules = unique(Schedules, 'rows', 'stable'); % keep unique schedules
    Schedules = Schedules(any(Schedules, 2), :); % remove rows that don't have treatment
    Schedules = Schedules(any(Schedules == 1, 2) & any(Schedules == 2, 2), :); % remove rows that don't have both OT-1 and Gem
end


function [IC, Iday, FEday, FE_vol] = flatten_original_mice( ...
    initialcondition_all, initialday_all, finalday_all, finalultrasound_all, nmpc)
    nMice = sum(nmpc);
    IC = nan(nMice, 4);
    Iday = nan(nMice, 1);
    FEday = nan(nMice, 1);
    FE_vol = nan(nMice, 1);
    k = 1;
    for c = 1:length(nmpc)
        for m = 1:nmpc(c)
            IC(k, :) = squeeze(initialcondition_all(c, m, :))';
            Iday(k) = initialday_all(c, m);
            FEday(k) = finalday_all(c, m);
            FE_vol(k) = finalultrasound_all(c, m);
            k = k + 1;
        end
    end
end

