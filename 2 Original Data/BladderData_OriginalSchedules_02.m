%% Bladder Cancer Data for Original Schedules
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Ultrasound data over time for each of the 55 experimental 
% mice from 4 experimental cohorts: untreated, Gem, OT-1, Gem+OT-1.

%% Volume data in terms of mm^3 

% All ultrasound data
Timedata = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_opt/GitHub/ultrasound_data.xlsx','Sheet','Original Schedules','Range','A2:A7'));
ultrasound_alldata = table2array(readtable('/Users/4481924/Documents/TIL_ODE/Code/TIL_opt/GitHub/ultrasound_data.xlsx','Sheet','Original Schedules','Range','B2:BI7'))';
% untreated, gem mono, ot1 mono, gem+ot1 combo--each have at most 15 mice with 5 time points
     
for i = 1:size(ultrasound_alldata,1)
    % 250 is a filler value in the excel document indicating that the mouse
    % died (not an actual measurement)
    death_idx = ultrasound_alldata(i,:) == 250;
    ultrasound_alldata(i,death_idx) = NaN;
end

% reshape ultrasound_alldata to remove NaN columns
ultrasound_alldata = ultrasound_alldata(~all(isnan(ultrasound_alldata),2),:);


%% Sort each mouse into correct cohort

nmpc = [14, 14, 12, 15]; % number of mice per cohort
nCohorts = length(nmpc);

exp_individual = nan(nCohorts,length(Timedata),max(nmpc)); % cohorts x times x number of mice
for cohort = 1:nCohorts
    lb = sum(nmpc(1:cohort-1)) +1;
    ub = sum(nmpc(1:cohort));
    exp_individual(cohort,:,1:nmpc(cohort)) = ultrasound_alldata(lb:ub,1:length(Timedata))';
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

dose_times = [10, 14];
cohort_dose_ot1 = [0,     0;            % cohort 1 - untreated
                   0,     0;            % cohort 2 - Gem
                   0,     1.426*prot1;  % cohort 3 - OT1
                   0,     1.426*prot1]; % cohort 4 - Gem+OT1
cohort_dose_gem = [0,     0;            % cohort 1 - untreated
                   3.8e4, 0;            % cohort 2 - Gem
                   0,     0;            % cohort 3 - OT1
                   3.8e4, 0];           % cohort 4 - Gem+OT1
