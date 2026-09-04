function mp = ChooseParamFunc_04e(fixed,params,indexchooseparam)
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: This code allows certain parameters to be fixed while others 
% are varied in order to be fitted to data. This helps with parameter 
% fitting and profile likelihoods if you only want to fit certain parameters.


% EXTRACT PARAMETERS: 
mp = fixed;
mp(indexchooseparam(:)) = params(:);  
    
end