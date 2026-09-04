%% Shared figure fonts, line weights, colors, and cohort labels
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Run it so plotting scripts share the same font size, 
% line widths, schedule colors, and new-experiment cohort names. 

fs = 16;
lineWidth = 1.5;
tickWidth = 0.8;
markerSize = 6;
meanMarkerSize = 8;

%% Therapy labels
D1 = 'O'; % OT-1
D2 = 'G'; % Gemcitabine

%% Cells
cell_names = {'Cancer', 'T cells', 'MDSCs', 'Total Tumor'};
cell_range = [0 300; 0 5; 0 25; 0 300];


%% Original Cohorts
orig_cohort_names = {'Untreated','Gem','OT-1','Gem+OT-1'};

orig_cohort_order = [1, 3, 2, 4]; % untreated, OT-1, Gem, Gem+OT-1


%% Cohorts for New Schedules
cohort_names = {'Original (GO)','Schedule 10 (GOOG)','Schedule 1 (GGGO)','Schedule 5 (GOGG)','Digital Twin'};

cohort_names_leg = {' Original (GO)',' S10 (GOOG)',' S1 (GGGO)',' S5 (GOGG)',' Digital Twin'};

cohort_order = [1,3,4,2,5];

% Cohort Colors
control_color = [116,188,31]./255; % green
sched_colors = [238, 155, 0; % marigold
    0, 65 85; % darker blue
    10, 147, 150]./255; % light blue

digitaltwin_color = [187, 62, 3]./255; % orange

colors = [control_color; 
          sched_colors;
          digitaltwin_color];


%% Experiment Colors
ec_color = [0.6, 0.6, 0.6]; % experimental color for comparison to virtual predictions


%% Color Palette (for OrderSchedule_05f.m)
palette = [ ...
    0, 65, 85;     % darker blue
    0, 95, 115;    % dark blue
    10, 147, 150;  % medium blue
    148, 210, 189; % light blue
    233, 216, 166; % beige
    238, 155, 0;   % marigold
    202, 103, 2;   % light orange
    174, 32, 18;   % medium red
    155, 34, 38    % dark red
    ] ./ 255;
