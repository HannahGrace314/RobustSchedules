%% Model Cost
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"

% Description: Cost function comparing ultrasound data to bladder cancer
% model using least squares error.

function error = ModelCost_04f(ODEFunc,params,initialcondition,timerange,data,data_times,dose,dose_times)

% Evaluate ODE:
numsim = Evaluate_NumSim_01c(ODEFunc,params,initialcondition,timerange,dose,dose_times);

T = numsim{2};
X = sum(numsim{1},2);

% Pull out data times:
x = X(ismember(T, data_times), :)';
t = nonzeros(ismember(T, data_times).*T);

% Least Squares Error:
error = sum((data - x).^2); 


end
