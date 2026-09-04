%% Keep schedules that satisfy experimental and dosing constraints
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: From EvaluateSchedule_05b output, retain weekly tumor summaries
% only for schedules that do not exceed the mouse's final experimental
% volume and that stay within the allowed OT-1, gemcitabine, and week counts.
%
% INPUTS:
%   evaluate_sched - nSched x nMice x nOut array from EvaluateSchedule_05b, with
%                    trailing entries [FC_ode, NumOT1, NumGem, NumWeeks]
%   FE_vol         - nMice x 1 final experimental tumor volume
%   MaxOT1         - maximum allowed OT-1 doses
%   MaxGem         - maximum allowed gemcitabine doses
%   MaxWeeks       - maximum allowed treatment weeks 
%
% OUTPUT:
%   pull_sched - nSched x nMice x 4 array [Cf, average, minimum, maximum]
%                at week MaxWeeks; NaN if the schedule is excluded

function pull_sched = PullSchedule_05c(evaluate_sched, FE_vol, MaxOT1, MaxGem, MaxWeeks)

nSched = size(evaluate_sched, 1);
num_expmice = size(evaluate_sched, 2);
n_out = size(evaluate_sched, 3);

pull_sched = nan(nSched, num_expmice, 4);

% Weekly blocks in evaluate_sched: Cf, average, min, max (length num_eval each)
num_eval = (n_out - 4) / 4; % (n_out - 4 entries for (FC_ode, NumOT1, NumGem, NumWeeks))/ 4 entries for Cf, averagetumor, minimumtumor, maximumtumor (which are the size of the max_weeks_tested))
week_idx = MaxWeeks + (0:3) * num_eval;

for mouse = 1:num_expmice
    for s = 1:nSched
        pull = evaluate_sched(s, mouse, :);
        FC_ode = pull(end - 3);
        NumOT1 = pull(end - 2);
        NumGem = pull(end - 1);
        NumWeeks = pull(end);

        if FC_ode < FE_vol(mouse) && NumOT1 <= MaxOT1 && ...
                NumGem <= MaxGem && NumWeeks <= MaxWeeks
            pull_sched(s, mouse, :) = evaluate_sched(s, mouse, week_idx);
        end
    end
end

end

