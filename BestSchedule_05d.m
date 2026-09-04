%% Select the best remaining schedule for each mouse
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Among schedules retained by PullSchedule_05c, discard those whose
% maximum tumor volume is at least 200 mm^3 (relaxing that cap in 10 mm^3
% steps only if every schedule is excluded). Rank the rest by the chosen
% objective and return the minimizing schedule index.
%
% INPUTS:
%   pull_sched - nSched x nMice x 4 array [Cf, average, minimum, maximum]
%   ofn        - objective functional
%                1: final tumor (Cf)
%                2: average tumor
%                3: minimum tumor
%                4: maximum tumor
%                5: final + average
%                6: 2*final + average
%
% OUTPUT:
%   best_sched - nMice x 1 indices into the schedule list

function best_sched = BestSchedule_05d(pull_sched, ofn)

vol_cap = 200;
[nSched,nMice,~] = size(pull_sched);

best_sched = nan(nMice, 1);

for mouse = 1:nMice
    keep = objective_values(pull_sched, mouse, ofn, vol_cap);

    % If every schedule exceeds the cap, raise it until at least one remains
    if all(isnan(keep))
        for extra = 10:10:100
            keep = objective_values(pull_sched, mouse, ofn, vol_cap + extra);
            if ~all(isnan(keep))
                break
            end
        end
    end

    [~, idx] = min(keep);
    best_sched(mouse) = idx;
end

end


function keep = objective_values(pull_sched, mouse, ofn, max_tumor)

nSched = size(pull_sched, 1);
keep = nan(nSched, 1);

for s = 1:nSched
    pull = squeeze(pull_sched(s, mouse, :));
    if any(isnan(pull))
        continue
    end
    if pull(4) >= max_tumor
        continue
    end

    if ofn <= 4
        keep(s) = pull(ofn);
    elseif ofn == 5
        keep(s) = pull(1) + pull(2);
    elseif ofn == 6
        keep(s) = 2 * pull(1) + pull(2);
    end
end

end
