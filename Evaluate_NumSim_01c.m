%% ODE evaluation function 
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Simulates the bladder cancer model (ODEFunc) with
% gemcitabine and OT-1 doses.
%
%   step - time grid for ode45 (default 0.1). 

function numsim = Evaluate_NumSim_01c(ODEFunc, params, initialcondition, ...
    timerange, dose, dose_times, step)

if nargin < 7 || isempty(step)
    step = 0.1;
end

params = abs(params);
times = [timerange(1), dose_times, timerange(2)];

[T, X] = ode45(@(t, x) ODEFunc(t, x, params), times(1):step:times(2), initialcondition);

for i = 3:length(times)
    [T1, X1] = ode45(@(t, x) ODEFunc(t, x, params), times(i-1):step:times(i), ...
        [X(end, 1), X(end, 2) + dose(1, i-2), X(end, 3), X(end, 4) + dose(2, i-2)]);
    T = [T; T1(2:end)];
    X = [X; X1(2:end, :)];
end

X = X(1:end, 1:3);
numsim = {X, T};

end
