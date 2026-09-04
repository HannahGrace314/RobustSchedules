%% Parameters for Bladder Cancer Model
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

% Description: Parameter set, ranges, and labels for bladder cancer model.
    
ODEFunc = @BladderFunc_01;

%% PARAMETER SET AND RANGES

% SET:
param =   [0.451;       % 1  - pc
           375;         % 2  - cmax
           0.509;       % 3  - ktc 
           0.351;       % 4  - kgc 
           0.229;       % 5  - ntc 
           5.4;         % 6  - smt 
           1.27;        % 7  - kgt 
           0.01;        % 8  - dT
           0.227;       % 9  - T0
           0.132;       % 10 - rcm
           0.398;       % 11 - kgm  
           0.0249;      % 12 - stm
           3.07;        % 13 - dm
           0.0485;      % 14 - M0
           1.54;        % 15 - dg 
           160];        % 16 - km


% RANGES: (Min/Max) 
mpr = [0            1;          % 1  - pc
       100          500;        % 2  - cmax
       0            1.5;        % 3  - ktc 
       0            0.25;       % 4  - kgc 
       0            0.5;        % 5  - nct 
       0            20;         % 6  - smt 
       0            1;          % 7  - kgt 
       0            0.5;        % 8  - dT 
       0            1;          % 9  - T0
       0            0.5;        % 10 - rcm 
       0            0.5;        % 11 - kgm 
       0            5;          % 12 - stm
       0.25         5;          % 13 - dm 
       0            1;          % 14 - M0 
       0            5;          % 15 - dg
       0            1000];      % 16 - km


%% LABELS: 

% PARAMETER NAMES:
pn = ["p_C";
     "C_{max}";
      "k_{TC}";
      "k_{GC}";
      "n_{CT}";
      "s_{MT}";
      "k_{GT}";
      "d_T";
      "T_0";
      "r_{CM}";
      "k_{GM}";
      "s_{TM}";
      "d_M";
      "M_0";
      "d_G";
      "k_M"];


% Specify axis labels for figures
axislabels = strings(size(mpr));
for i = 1:length(mpr(:,1))
    for j = 1:length(mpr(1,:))
        if mpr(i,j) == 0 || mpr(i,j) >=1e-2
            axislabels(i,j)=string(mpr(i,j));
        end
        if mpr(i,j) <1e-2 && mpr(i,j)>0
            first = extractBefore(compose("%5.0e",mpr(i,j)),"e");
            last = extractAfter(compose("%5.0e",mpr(i,j)),"0");
            axislabels(i,j) = append(first,'0^{-',last,'}');
        end
    end
end
