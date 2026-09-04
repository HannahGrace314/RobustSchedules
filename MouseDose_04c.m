%% Dose matrix for one mouse in a fitting cohort
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"

function dose = MouseDose_04c(cohort, mouse, cohort_dose_ot1, cohort_dose_gem, ...
    digtwin_sched, nCohorts)

if ~isempty(digtwin_sched) && cohort == nCohorts
    s = digtwin_sched(mouse);
    dose = [cohort_dose_ot1(s, :); cohort_dose_gem(s, :)];
else
    dose = [cohort_dose_ot1(cohort, :); cohort_dose_gem(cohort, :)];
end

end
