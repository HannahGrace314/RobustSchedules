%% Validate virtual cohort tumor volumes against experimental ultrasound
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Aligns virtual cohort total tumor (from VC_NumSim_07d) to
% the ultrasound days in BladderData_NewSchedules_07, censors volumes at each
% mouse's death/endpoint, and compares experimental vs virtual cohorts.
% Reports mean/median/max/min/IQR/SD, violin plots of the four tested
% schedules, and Kolmogorov-Smirnov plus two-sample t-tests at each day.
%
% INPUTS (workspace / files)
%   vc_numsim, numsim_times - from VC_NumSim_07d.m
%                            (cohort x cell type x time x virtual mouse;
%                             cell type 4 is total tumor)
%   death_vols              - from virtualcohort.mat (per virtual mouse)
%   BladderData_NewSchedules_07.m - exp_individual, Timedata, nmpc (uses
%                               use_endpt_vol set below)
%   Figure_Aesthetics_03.m     - colors, cohort_order, cohort_names_leg, fs
%
% OUTPUTS
%   datastats, vcstats - structs of summary statistics (cohort x time)
%                        for experimental and virtual tumor volumes
%   pull_vcnumsim      - virtual total tumor at ultrasound days, NaN after
%                        min(death_vols, 200)
%   violin figure      - experimental violins (black) vs virtual (schedule
%                        colors)
%   h_ks, p_ks         - KS test, reject-null flag and p-value
%                        (cohort x ultrasound day)
%   h_ttest, p_ttest   - two-sample t-test on means, same size
%
% A p-value >= 0.05 (h = 0) means the virtual and experimental samples
% are not significantly different at that day.

nCohorts = 5;

%% Summary Statistics of Experimental Cohort
use_endpt_vol = 0;
BladderData_NewSchedules_07

% Experimental Cohort Stats
datastats = summary_stats(exp_individual);


%% Summary Statistics of Virtual Cohort 

[~, idx] = ismember(Timedata, numsim_times);

pull_vcnumsim = squeeze(vc_numsim(:,4,idx,:)); % index 4 to grab total tumor size
for mouse = 1:size(vc_numsim,4)
    grab_matrix = pull_vcnumsim(:,:,mouse);
    grab_matrix(grab_matrix>=min(death_vols(mouse),200))=NaN;
    pull_vcnumsim(:,:,mouse) = grab_matrix;
end
    
% Virtual Cohort Stats
vcstats = summary_stats(pull_vcnumsim); % input: cohorts x times x number of mice


%% Violin Plots
Figure_Aesthetics_03
ec_color = [0,0,0];
nCohorts_plot = nCohorts-1;

figure
outer = tiledlayout(1, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

vp = tiledlayout(outer, 1, nCohorts_plot, 'TileSpacing', 'compact', 'Padding', 'compact');
vp.Layout.Tile = 1;
vp.Layout.TileSpan = [1 5];
for d = 1:nCohorts_plot
    nexttile(vp)
    grid on
    hold on

    cohort = cohort_order(d);
    violinplot(squeeze(exp_individual(cohort,1:6,:))',num2cell(Timedata(1:6)),'ViolinColor',ec_color,'ViolinAlpha',0.2)
    violinplot(squeeze(pull_vcnumsim(cohort,1:6,:))',num2cell(Timedata(1:6)),'ViolinColor',colors(cohort,:),'ShowData',false)
    if d>1
        set(gca,'YTickLabel',[])
    end
    box off
    set(gca,'FontSize',fs-2,'ylim',[0, 200])
end
xlabel(vp,'Days','FontSize',fs)
ylabel(vp,'Tumor Volume (mm^3)','FontSize',fs)

% One legend per axes: a second column so both legends can exist at once
legcol = tiledlayout(outer, 2, 1, 'TileSpacing', 'compact', 'Padding', 'tight');
legcol.Layout.Tile = 6;

nexttile(legcol)
hold on
axis off
ec = fill(nan(4, 1), nan(4, 1), ec_color, 'FaceAlpha', 0.2, 'EdgeColor', ec_color);
data = scatter(nan, nan, [], 'filled', 'MarkerFaceColor', ec_color, 'MarkerEdgeColor', 'k');
lgd_ec = legend([ec, data], [" Violin", "Data Point"]);
lgd_ec.Title.String = 'Experiments';
lgd_ec.FontSize = fs;
lgd_ec.Box = 'on';
lgd_ec.Location = 'west';
lgd_ec.AutoUpdate = 'off';


nexttile(legcol)
hold on
axis off
h_vc = gobjects(nCohorts_plot, 1);
for c = 1:nCohorts_plot
    cohort = cohort_order(c);
    h_vc(c) = fill(nan(4, 1), nan(4, 1), colors(cohort, :), ...
        'FaceAlpha', 0.3, 'EdgeColor', colors(cohort, :));
end
lgd_vc = legend(h_vc, cohort_names_leg(cohort_order(1:nCohorts_plot)));
lgd_vc.Title.String = 'Virtual Predictions';
lgd_vc.FontSize = fs;
lgd_vc.Box = 'on';
lgd_vc.Location = 'west';
lgd_vc.AutoUpdate = 'off';

set(gcf, 'Position',[1744 583 1180 260])


%% Kolmogorov-Smirnov (compares distributions) and Student's t test (compares means)
% for statistically comparing virtual and experimental cohorts
h_ks = nan(nCohorts, length(Timedata));
p_ks = nan(nCohorts, length(Timedata));
h_ttest = nan(nCohorts, length(Timedata));
p_ttest = nan(nCohorts, length(Timedata));

for cohort = 1:nCohorts
    virtual_data = squeeze(pull_vcnumsim(cohort,:,:));
    exp_data = squeeze(exp_individual(cohort,:,:));
    exp_data(all(isnan(exp_data),2),:) = []; % removes rows with no mice data

    for day = 1:length(Timedata)
        data1 = virtual_data(day,:); 
        data2 = exp_data(day,~isnan(exp_data(day,:)));

        % Compare full distributions (Kolmogorov-Smirnov)
        [h_ks(cohort,day),p_ks(cohort,day)] = kstest2(data1, data2); 
        % Want p_ks >= 0.05 (so no significant difference found), and 
        % h_ks = 0 (meaning fail to reject the null hypothesis, ie the data comes from similar population)

        % Compare means (Student's t-test)
        [h_ttest(cohort,day),p_ttest(cohort,day)] = ttest2(data1, data2);
        % again, p_ttest >= 0.05 and h_ttest = 0
    end

end


%% Calculate Stats
function stats = summary_stats(input)
% input is a matrix of nCohorts x tsize (num of timepts) x number of mice
    nCohorts = length(input(:,1,1));
    tsize = length(input(1,:,1));
    
    stats.mean   = zeros(nCohorts, tsize);
    stats.median = zeros(nCohorts, tsize);
    stats.max    = zeros(nCohorts, tsize);
    stats.min    = zeros(nCohorts, tsize);
    stats.iqr    = zeros(nCohorts, tsize);
    stats.std    = zeros(nCohorts, tsize);

    for c = 1:nCohorts
        for t = 1:tsize
            grab_data = squeeze(input(c,t,:));
            grab_data = grab_data(~isnan(grab_data));
            
            stats.mean(c,t)   = mean(grab_data);
            stats.median(c,t) = median(grab_data);
            stats.max(c,t)    = max(grab_data);
            stats.min(c,t)    = min(grab_data);
            stats.iqr(c,t)    = iqr(grab_data);
            stats.std(c,t)    = std(grab_data);
        end
    end
end