%% Fit chosen ODE parameters to an ultrasound trajectory
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: fminsearchbnd minimizes ModelCost_04f (least squares error 
% between total tumor and ultrasound). Only entries in indexchooseparam 
% are estimated; the rest of fixed stay at their given values. 
% After the fit, the ODE is resimulated for R^2 and MAE.
%
% INPUTS:
%   ODEFunc           - ODE right-hand side (e.g. @BladderFunc_01)
%   fixed             - full nominal parameter vector
%   indexchooseparam  - indices of parameters to estimate
%   paramrange        - lower/upper bounds for those parameters (n x 2)
%   initialcondition  - ODE state at timerange(1)
%   timerange         - [t0 tf] for the simulation
%   data, data_times  - ultrasound volumes and their times
%   dose, dose_times  - treatment matrix and injection days
%
% OUTPUT:
%   g - {final_param, exitflag, Rsquared, MAE}
%       final_param is the full parameter vector with fitted entries filled in

function g = ParameterFit_04d(ODEFunc,fixed,indexchooseparam,paramrange, initialcondition,timerange,data,data_times,dose,dose_times)

    params = fixed(indexchooseparam);
    [paramests, fval, exitflag] = fminsearchbnd(@(p) ModelCost_04f(ODEFunc, ChooseParamFunc_04e(fixed,p,indexchooseparam), initialcondition,timerange,data,data_times,dose,dose_times),params,paramrange(:,1),paramrange(:,2),optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000,'TolX',1e-7,'TolFun',1e-7));
    final_param = ChooseParamFunc_04e(fixed,paramests,indexchooseparam);

    % Calculating R^2
    SS_res = fval;
    SS_tot = sum((data-mean(data)).^2);
    Rsquared = 1-(SS_res/SS_tot);

    % Re-Simulating ODE with paramests to determine MAE
    pull_numsim = Evaluate_NumSim_01c(ODEFunc,final_param,initialcondition,timerange,dose,dose_times);
    X = sum(pull_numsim{1},2); 
    T = pull_numsim{2};
    x = X(ismember(T, data_times), :)';

    % Mean Absolute Error (MAE)
    MAE = sum(abs(data - x))/length(data);
    
    g = {final_param,exitflag,Rsquared, MAE};
end