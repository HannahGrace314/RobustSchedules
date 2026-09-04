%% Kaplan-Meier staircase from event times and censoring flags
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Product-limit estimator. Survival drops only at death times;
% censored mice leave the risk set but do not drop the curve. Returns
% staircase coordinates so plot(x, 100*y) is a KM curve from day 0.
%
% INPUTS:
%   times      - event or censoring day for each mouse
%   censored   - 1 if the mouse was censored (still alive), 0 if death
%   t_end      - last day drawn on the plot
%
% OUTPUTS:
%   x, y       - staircase time and survival in [0, 1]

function [x, y] = KaplanMeier_07c(times, censored, t_end)
    
    x = 0;
    y = 1;
    
    death_times = unique(times(censored == 0));
    n_total = length(times);
    for t = death_times'
        n_alive = sum(times > t);
        x = [x; t; t];
        y = [y; y(end); n_alive/n_total]; 
    end
    
    x = [x; t_end];
    y = [y; y(end)];

end
