%% Rank optimized schedules on the virtual murine cohort
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Takes all_output from Treatment_Sweep_05, simulates each
% unique schedule (plus the original schedule) on the virtual cohort and
% a 2-D parameter grid, ranks schedules by tumor and survival metrics, and
% makes bar / pie / scatter figures.
%
% Needs all_output in the workspace or treatment_sweep_all_output.mat, and
% virtualcohort.mat (vc_param, vc_initialcondition, death_vols).

%% Decide on Objective Functional (ofn) and MaxWeeks

ofn = 5;
MaxWeeks = 2;

if ~exist('all_output', 'var')
    load('treatment_sweep_all_output.mat', 'all_output');
end
output = all_output{ofn, MaxWeeks};

[optimized_sched_max, ~] = size(output);
nSched = optimized_sched_max + 1; % optimized schedules + original

%% Call Virtual Cohort and Set up Parameter Space  

% Virtual Cohort
load('virtualcohort.mat') % vc_param, vc_initialcondition, death_vols

% Parameter Space
BladderParam_01b % mpr, param
p1 = 5;   % ntc
p2 = 10;  % rcm

x = linspace(mpr(p1, 1), mpr(p1, 2), 11);
y = linspace(mpr(p2, 1), mpr(p2, 2), 11);
[X, Y] = meshgrid(x, y);
coords = [X(:), Y(:)];
parameterspace = repmat(param', size(coords, 1), 1);
parameterspace(:, p1) = coords(:, 1);
parameterspace(:, p2) = coords(:, 2);


%% Call Treatment Schedules to Test (along with Colors and Labels)  
Figure_Aesthetics_03

ot1_orig = 1.426; % retention is applied in cohort_response
gem_orig = 3.8e4;

schedules = cell(nSched, 1);
times = cell(nSched, 1);
sched_colors = zeros(nSched, 3);
order_labels = strings(nSched, 1);
sched_labels = strings(nSched, 1);
for s = 1:nSched - 1
    % Schedules and times
    schedules{s} = output{s, 1}{1, 1};
    times{s} = output{s, 1}{1, 2};

    % Colors
    sched_colors(s, :) = output{s, 2};

    % Labels
    order_labels(s) = output{s,3};
    sched_labels(s) = strcat('S',string(s));
end

% Original Schedule
schedules{end} = [0, gem_orig; ot1_orig, 0];
times{end} = [10; 14];
sched_colors(end, :) = control_color; % original
order_labels(end) = "Original";
sched_labels(end) = "Original";

[~,sched_for_leg] = unique(order_labels);
sched_for_leg = sort(sched_for_leg);


%% Apply Treatment Schedules to Virtual Cohort and Parameter Space  

t_initial = 6; % initial day
t_endpt = 800; % day mouse dies naturally
ot1retention = 0.75;
ODEFunc = @BladderFunc_01;

% Virtual Cohort 
virtualcohort_response = cohort_response(ODEFunc, vc_param, vc_initialcondition, ...
    times, schedules, t_initial, t_endpt, ot1retention);

% Parameter Space
% Day-6 IC from linear interpolation (same for every grid mouse)
initialcondition = [26.8545112900343, 0.0721738874774849, 0.728763098350346, 0];
initialcondition_ps = repmat(initialcondition, size(parameterspace, 1), 1);
parameterspace_response = cohort_response(ODEFunc, parameterspace, initialcondition_ps, ...
    times, schedules, t_initial, t_endpt, ot1retention);


%% Analyze Response of Virtual Cohort and Parameter Space  

smallest_ultrasound = 3.683; % smallest nonzero ultrasound (mm^3)
exp_tumor_endpt = 200; % experimental endpoint volume (mm^3)

start_times = zeros(nSched, 1);
end_times = zeros(nSched, 1);
for s = 1:nSched
    start_times(s) = times{s}(1);
    end_times(s) = times{s}(end);
end
t_start = min(start_times); % first day of treatment
endofweek = 16:7:37;
candidates = endofweek(endofweek > max(end_times));
tf = min(candidates); % end of week after last day of treatment

virtualcohort_analysis = analyze_response(virtualcohort_response, ...
    death_vols, smallest_ultrasound, t_start, tf, t_endpt);
parameterspace_analysis = analyze_response(parameterspace_response, ...
    exp_tumor_endpt*ones(size(parameterspace, 1),1), smallest_ultrasound, t_start, tf, t_endpt);


%% Rank Schedules according to Different Objectives  

obj_titles = ["Final Tumor Size", "Final Cancer Size", "Cumulative Tumor Size", ...
    "Maximum Tumor Size", "Minimum Tumor Size", "Average Tumor Size", ...
    "Overall Survival (OS)", "First Day Disease-Free", "Disease-Free Survival (DFS)"];
num_obj = length(obj_titles);
rank_cohort = cell(num_obj, 1);
rank_space = cell(num_obj, 1);

for obj = 1:num_obj
    if strcmpi(obj_titles(obj), "Overall Survival (OS)") || ...
            strcmpi(obj_titles(obj), "Disease-Free Survival (DFS)")
        sort_mode = 'large best';
    else
        sort_mode = 'small best';
    end

    rank_cohort{obj} = rank_schedules(virtualcohort_analysis{1, obj}, obj_titles(obj), sort_mode);
    rank_space{obj}  = rank_schedules(parameterspace_analysis{1, obj}, obj_titles(obj), sort_mode);
end


%% Figures - bar graphs, pie chart, scatter plots, and trajectories  

savefigure = [0, 0, 0, 0, 0]; % bar all, specific bar, pie, scatter, trajectories
dateStr = datestr(now, 'ddmmmyyyy');

specificrank = 1;
which_obj = 1; % see obj_titles
numranks = nSched;

jpg_string = ["fts", "fcs", "cumulativetumor", "maxtumor", "mintumor", ...
    "averagetumor", "OS", "firstdaydf", "DFS"];

for k = 1:length(which_obj)
    obj = which_obj(k);

    Rank_Barplot_06b(rank_cohort{obj}{1}, obj_titles(obj), sched_colors, order_labels, ...
        sched_for_leg, numranks);
    save_rank_fig(savefigure(1), jpg_string(obj), "barall", ofn, MaxWeeks, dateStr);

    SpecificRank_Barplot_06c(rank_cohort{obj}{1}, obj_titles(obj), sched_colors, ...
        sched_labels, specificrank);
    save_rank_fig(savefigure(2), jpg_string(obj), "bar", ofn, MaxWeeks, dateStr);

    SpecificRank_Piechart_06d(rank_cohort{obj}{1}, obj_titles(obj), sched_colors, ...
        order_labels, specificrank);
    save_rank_fig(savefigure(3), jpg_string(obj), "pie", ofn, MaxWeeks, dateStr);

    Rank_Scatter_06e(vc_param, parameterspace, rank_space{obj}{2}, ...
        rank_space{obj}{1}, rank_cohort{obj}{1}, obj_titles(obj), sched_colors, ...
        sched_labels, specificrank, p1, p2, pn); % can choose to do vector with multiple ranks here instead of just specificrank
    save_rank_fig(savefigure(4), jpg_string(obj), "scatter_nCTrCM", ofn, MaxWeeks, dateStr);
end

VC_Trajectories_06f(virtualcohort_response, virtualcohort_analysis, sched_colors, tf)
save_rank_fig(savefigure(5),'vc','trajectories',ofn, MaxWeeks, dateStr);


%% Applying Treatment to Cohort  

function cohort_response = cohort_response(ODEFunc, params, IC, dose_times, schedules, t_initial, t_endpt, ot1retention)
    nMice = size(params, 1);
    nSched = numel(schedules);
    response = cell(nMice, nSched);
    times_response = cell(nMice, nSched);
    
    for s = 1:nSched
        dose = [schedules{s}(:, 1)' * ot1retention; schedules{s}(:, 2)'];
        t_dose = dose_times{s}(:)';
        parfor mouse = 1:nMice
            numsim = Evaluate_NumSim_01c(ODEFunc, params(mouse, :), IC(mouse, :), ...
                [t_initial, t_endpt], dose, t_dose, 1);
            response{mouse, s} = numsim{1};
            times_response{mouse, s} = numsim{2};
        end
    end
    cohort_response = {response, times_response};
end


%% Analyze Response  

function analyze_response = analyze_response(cohort_response, ...
    largest_ultrasound, smallest_ultrasound, t_start, tf, t_endpt)

    response = cohort_response{1};
    times_response = cohort_response{2};
    [nMice, nSched] = size(response);
    
    finaltumorsize = zeros(nMice, nSched);
    finalcancersize = zeros(nMice, nSched);
    totaltumorsize = zeros(nMice, nSched);
    maxtumorsize = zeros(nMice, nSched);
    mintumorsize = zeros(nMice, nSched);
    averagetumorsize = zeros(nMice, nSched);
    dayofdeath = t_endpt * ones(nMice, nSched);
    firstdaydiseasefree = t_endpt * ones(nMice, nSched);
    dfs = zeros(nMice, nSched);
    
    for s = 1:nSched
        df_schedule = zeros(nMice, t_endpt); 
        for mouse = 1:nMice
            t = times_response{mouse, s};
            totalsize_eachtime = sum(response{mouse, s}(:, 1:3), 2);
    
            i_ts = find(t <= t_start, 1, 'last'); % start index
            i_tf = find(t >= tf, 1, 'first'); % end index
    
            finaltumorsize(mouse, s) = totalsize_eachtime(i_tf);
            finalcancersize(mouse, s) = response{mouse, s}(i_tf, 1); 
            totaltumorsize(mouse, s) = trapz(t(i_ts:i_tf), totalsize_eachtime(i_ts:i_tf));
            maxtumorsize(mouse, s) = max(totalsize_eachtime(i_ts:i_tf));
            mintumorsize(mouse, s) = min(totalsize_eachtime(i_ts:i_tf));
            averagetumorsize(mouse, s) = mean(totalsize_eachtime(i_ts:i_tf));
    
            dead = find(totalsize_eachtime >= largest_ultrasound(mouse)); % indices
            free = totalsize_eachtime <= floor(smallest_ultrasound); % logical
            if ~isempty(dead)
                dayofdeath(mouse, s) = ceil(t(dead(1)));
                free(dayofdeath(mouse, s):end) = 0;
            end
    
            daysdf = unique(ceil(t(logical(free))));
            if ~isempty(daysdf)
                df_schedule(mouse, daysdf) = 1;
                firstdaydiseasefree(mouse, s) = daysdf(1);
                dfs(mouse,s) = sum(df_schedule(mouse,:));
            end
        end
    end

    analyze_response = {finaltumorsize, finalcancersize, totaltumorsize, ...
        maxtumorsize, mintumorsize, averagetumorsize, dayofdeath, ...
        firstdaydiseasefree, dfs};
end


%% Ranking Schedules  

function rank = rank_schedules(response, distinguisher, sort_mode)

    % Determining sort direction
    if strcmpi(sort_mode, 'small best')
        sort_direction = 'ascend';
    else
        sort_direction = 'descend';
    end
    
    [nMice, nSched] = size(response);
    rank_matrix = zeros(nMice, nSched);
    for mouse = 1:nMice
        [sort_response, sort_idx] = sort(response(mouse,:),sort_direction);
        [~, ~, m_rank_idx] = unique(sort_response); 
        rank_matrix(mouse, sort_idx) = m_rank_idx;
    end

    % Counting how often each schedule gets each rank
    rank_counts = zeros(nSched, nSched);
    for s = 1:nSched
        for rnk = 1:nSched
            rank_counts(s, rnk) = sum(rank_matrix(:, s) == rnk);
        end
    end

    % Convert counts to percentages 
    rank_percent = 100 * rank_counts / nMice; % rows: schedule x cols: percent that rank it that col #

    rank = {rank_percent, rank_matrix};

    % Display the result
    fprintf('Percentage of mice that gave each schedule a specific rank according to %s:\n', distinguisher);
    disp(array2table(rank_percent, ...
        'VariableNames', compose("Rank_%d", 1:nSched), ...
        'RowNames', compose("Schedule_%d", 1:nSched)));
end


%% Save Figure Function  

function save_rank_fig(do_save, metric, kind, ofn, MaxWeeks, dateStr)
    if do_save == 1
        fileName = sprintf('%s_%s_ofn%d_MaxWeeks%d_%s.jpg', metric, kind, ofn, MaxWeeks, dateStr);
        exportgraphics(gcf, fileName, 'Resolution', 300);
    end
end