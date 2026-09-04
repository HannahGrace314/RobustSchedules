%% Evaluate a treatment schedule on the bladder cancer ODE model
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"
%
% Description: Simulates one gemcitabine/OT-1 schedule and returns tumor
% summaries at weekly endpoints, the volume on the final experimental day,
% and dose counts.
%
% INPUTS:
%   schedule - 1 x nDays vector of daily treatment codes after day 9
%             (0 = none, 1 = OT-1, 2 = gemcitabine)
%   pars    - model parameter vector
%   IC      - 1 x 4 initial condition [C, T, M, G]
%   Iday    - start day of the simulation
%   FEday   - experimental day used for FC_ode (must lie on the 0.1-day grid)
%   tf      - final simulation time
%   scale   - 0: full labeled dose at each injection
%             1: split the labeled dose equally across injections of that drug
%
% OUTPUT:
%   evaluate_sched = [Cf, averagetumor, minimumtumor, maximumtumor, ...
%                   FC_ode, NumOT1, NumGem, NumWeeks]
%   Weekly summaries (Cf, average, min, max) are concatenated in that order,
%   one entry per week ending on 16:7:tf. Unused later weeks remain NaN.

function evaluate_sched = EvaluateSchedule_05b(schedule, pars, IC, Iday, FEday, tf, scale)

order = nonzeros(schedule);
num_dose = length(order);

Treatment = zeros(num_dose, 2);
Treatment(order == 1, 1) = 1; % OT-1
Treatment(order == 2, 2) = 1; % gemcitabine

% Calendar days of injection (treatment days begin on day 10)
T = 9 + find(schedule ~= 0);

gem_dose = 3.8e4;
ot1_dose = 1.426;
prot1 = 0.75; % OT-1 remaining in bladder after intravesical injection (private communication)

if scale == 0
    dose_scale = [1, 1];
elseif scale == 1
    dose_scale(1) = 1 / sum(Treatment(:, 1));
    dose_scale(2) = 1 / sum(Treatment(:, 2));
end

% Evaluate_NumSim_01c convention: row 1 = OT-1, row 2 = gemcitabine
dose = zeros(2, num_dose);
dose(1, :) = prot1 * ot1_dose * Treatment(:, 1)' * dose_scale(1);
dose(2, :) = gem_dose * Treatment(:, 2)' * dose_scale(2);

numsim = Evaluate_NumSim_01c(@BladderFunc_01, pars, IC, [Iday, tf], dose, T(:)');
y = numsim{1};
t = numsim{2};

NumOT1 = nnz(Treatment(:, 1));
NumGem = nnz(Treatment(:, 2));

endofweek = 16:7:tf;
final_tday = arrayfun(@(tday) min(endofweek(endofweek >= tday)), T(end));
NumWeeks = find(endofweek == final_tday);

nweek = length(endofweek);
Cf = nan(1, nweek);
averagetumor = nan(1, nweek);
minimumtumor = nan(1, nweek);
maximumtumor = nan(1, nweek);

for k = 1:nweek
    if round(t(end)) >= endofweek(k)
        [~, idx] = min(abs(t - endofweek(k)));
        totaltumor = sum(y(1:idx, 1:3), 2);
        Cf(k) = totaltumor(end);
        averagetumor(k) = mean(totaltumor);
        minimumtumor(k) = min(totaltumor);
        maximumtumor(k) = max(totaltumor);
    end
end

FC_ode = sum(y(t == FEday, 1:3), 2);

evaluate_sched = [Cf, averagetumor, minimumtumor, maximumtumor, ...
    FC_ode, NumOT1, NumGem, NumWeeks];

end
