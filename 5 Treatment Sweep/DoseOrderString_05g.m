%% Dose labels by order
% Original Source: Anderson et al. (2026) - "Robust schedules of adoptive
% cell therapy in a virtual murine cohort of bladder cancer with
% experimental validation"

function s = DoseOrderString_05g(dose_mat,D1,D2)
    s = '';
    if isempty(dose_mat)
        return
    end
    for i = 1:size(dose_mat, 1)
        if dose_mat(i, 2) ~= 0
            s = [s, D2]; 
        elseif dose_mat(i, 1) ~= 0
            s = [s, D1];
        end
    end
end
