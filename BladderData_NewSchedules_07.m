%% Bladder Cancer Data for New Schedules
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Ultrasound data over time for the experimental mice on GO,
% GOOG, GGGO, GOGG, and digital-twin schedules (Cage set up sheet).
% Survival day is the earlier of first volume >= endpoint_vol and the
% recorded dayofdeath; mice with neither are censored at the last time.
%
% Set before calling: use_endpt_vol =
%   0 - volumes at/above endpoint_vol replaced with NaN (Violin)
%   1 - those volumes are replaced with endpoint_vol (mean tumor growth)
%   2 - keep the measured volumes (individual growth curves)
%
% OUTPUTS include exp_individual, Timedata, nmpc, endpoint_vol,
% km_day_all (cohort x mouse), km_censored_all (1 = censored, 0 = death).

%% Volume data in terms of mm^3 
excel_sheet = 'New Schedules';
excel_range = [8,29];

% All ultrasound data
Timedata = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_opt/GitHub/ultrasound_data.xlsx','Sheet',excel_sheet,'Range',sprintf('A%s:A%s',num2str(excel_range(1)),num2str(excel_range(2)))));
ultrasound_alldata = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_opt/GitHub/ultrasound_data.xlsx','Sheet',excel_sheet,'Range',sprintf('C%s:AN%s',num2str(excel_range(1)),num2str(excel_range(2)))))';
% GO, GOOG, GGGO, GOGG, and digital twin cohorts--each have at most 8 mice with over a dozen time points

dayofdeath = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_opt/GitHub/ultrasound_data.xlsx','Sheet',excel_sheet,'Range','C31:AN31'))';

mouse_labels = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_opt/GitHub/ultrasound_data.xlsx','Sheet',excel_sheet,'Range','C4:AN4'))';
mddi = [4,7,5,39,12,30]; % mouse died during injection

endpoint_vol = 200;

n_raw = size(ultrasound_alldata, 1);
first_endpt_day = nan(n_raw, 1);

% use_endpt_vol == 0  - keep the original data for individual growth curves
                 % 1  - treat endpoint/death as endpoint_vol

for i = 1:n_raw
    % 250 is a filler value in the excel document indicating that the mouse
    % died (not an actual measurement)
    death_idx = ultrasound_alldata(i,:) == 250;
    mouse = ultrasound_alldata(i,:);
    mouse(death_idx) = NaN;
    ep = find(mouse >= endpoint_vol, 1, 'first');
    if ~isempty(ep)
        first_endpt_day(i) = Timedata(ep);
    end

    ultrasound_alldata(i, death_idx) = NaN;

    if  use_endpt_vol == 0 % set measurements at/above endpoint and death to NaN
        endpoint_idx = ultrasound_alldata(i,:) >= endpoint_vol;
        ultrasound_alldata(i, endpoint_idx | death_idx) = NaN;
    elseif use_endpt_vol == 1 % set measurements at/above endpoint and death to endpoint_vol
        endpoint_idx = ultrasound_alldata(i,:) >= endpoint_vol;
        ultrasound_alldata(i, endpoint_idx | death_idx) = endpoint_vol;
    end
end

% Exclude mice that died on treatment
keep = ~ismember(mouse_labels,mddi);
ultrasound_alldata = ultrasound_alldata(keep, :);
dayofdeath = dayofdeath(keep);
first_endpt_day = first_endpt_day(keep);


%% Sort each mouse into correct cohort

nmpc = [6,6,6,6,8]; % number of mice per cohort

nCohorts = length(nmpc);

exp_individual = nan(nCohorts,length(Timedata),max(nmpc)); % cohorts x times x number of mice
dayofdeath_all = nan(nCohorts, max(nmpc));
first_endpt_all = nan(nCohorts, max(nmpc));
km_day_all = nan(nCohorts, max(nmpc));
km_censored_all = nan(nCohorts, max(nmpc));
for cohort = 1:nCohorts
    lb = sum(nmpc(1:cohort-1)) + 1;
    ub = sum(nmpc(1:cohort));
    exp_individual(cohort,:,1:nmpc(cohort)) = ultrasound_alldata(lb:ub,1:length(Timedata))';
    dayofdeath_all(cohort, 1:nmpc(cohort)) = dayofdeath(lb:ub);
    first_endpt_all(cohort, 1:nmpc(cohort)) = first_endpt_day(lb:ub);
    for mouse = 1:nmpc(cohort)
        t_cands = [first_endpt_all(cohort, mouse), dayofdeath_all(cohort, mouse)];
        t_cands = t_cands(isfinite(t_cands));
        if ~isempty(t_cands)
            km_day_all(cohort, mouse) = min(t_cands);
            km_censored_all(cohort, mouse) = 0; 
        else
            km_day_all(cohort, mouse) = Timedata(end); 
            km_censored_all(cohort, mouse) = 1; % censored because does not die by end of Timedata
        end

    end
end


%% Get times data for each mouse

times_all = repmat(reshape(Timedata,1,[],1), nCohorts, 1, max(nmpc)); % cohorts x times x number of mice
times_all(isnan(exp_individual)) = NaN;


%% Initial and Final Conditions

% Linear interpolation results for C, T, M from flow/histology data from
% Anderson et al. (2026)
linearinterp = [0.971038727058824, 0.00260975294117647, 0.0263515200000000;  % day 6
                0.956558090588235, 0.00391462941176471, 0.0395272800000000]; % day 9

initialday_all = nan(nCohorts,max(nmpc));
initialcondition_all = nan(nCohorts,max(nmpc),4);
finalday_all = nan(nCohorts,max(nmpc));
finalultrasound_all = nan(nCohorts,max(nmpc));

for cohort = 1:nCohorts
    for mouse = 1:max(nmpc)
        vector = exp_individual(cohort,:,mouse);
        if isempty(nonzeros(~isnan(vector)))==0
            % First available measurement
            idx = find(vector ~= 0 & ~isnan(vector), 1, 'first');
            initialday_all(cohort,mouse) = times_all(cohort,idx,mouse);
            interp_idx = find(Timedata == initialday_all(cohort,mouse));
            initialcondition_all(cohort,mouse,:) = [exp_individual(cohort,idx,mouse)*linearinterp(interp_idx,:),0]; 

            % Last available measurement
            idx = find(~isnan(vector), 1, 'last');
            finalday_all(cohort,mouse) = times_all(cohort,idx,mouse);
            finalultrasound_all(cohort,mouse) = exp_individual(cohort,idx,mouse);

        end
    end
end


%% Treatment injection based on treatment cohort
prot1 = 0.75; % percent remaining after OT-1 injection (based on human data - private communication)

dose_times = [10, 14, 17, 18, 21,23];
cohort_dose_ot1 = [0, 1.426*prot1, 0,           0,     0, 0;            % cohort 1 - GO
                   0, 1.426*prot1, 1.426*prot1, 0,     0, 0;            % cohort 2 - GOOG
                   0, 0,           0,           0,     1.426*prot1, 0;  % cohort 3 - GGGO
                   0, 1.426*prot1, 0,           0,     0, 0];           % cohort 4 - GOGG
cohort_dose_gem = [3.8e4, 0,     0,     0,     0, 0;            % cohort 1 - GO
                   3.8e4, 0,     0,     0,     3.8e4, 0;        % cohort 2 - GOOG
                   3.8e4, 3.8e4, 3.8e4, 0,     0, 0;            % cohort 3 - GGGO
                   3.8e4, 0,     0,     3.8e4, 3.8e4, 0];       % cohort 4 - GOGG

digtwin_sched = [4,4,2,4,2,2,4,2];





