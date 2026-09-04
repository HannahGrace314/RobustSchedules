%% ODE Model of bladder cancer dynamics including treatment with gemcitabine
% Original Source: Anderson et al. (2026) - "A virtual cohort framework
% with applications to adoptive cell therapy in bladder cancer"

function dxdt = BladderFunc_01(t,x,p)

% INPUTS: 
C = x(1);
T = x(2);
M = x(3);
G = x(4);
dxdt = zeros(4,1);

%  x is 4x1 vector with:
%  C = number of tumor cells
%  T = number of T cells
%  M = number of MDSCs 
%  G = gemcitabine


% --------------- Parameter values for cancer equation --------------------
pc = p(1);
cmax = p(2);
ktc = p(3); 
kgc = p(4); 

%  pc   = cancer growth rate
%  cmax = carrying capacity
%  ktc  = T cell kill rate of tumor cells
%  kgc  = gemcitabine kill rate of cancer cells


% --------------- Parameter values for T cell equation --------------------
nct = p(5);
smt = p(6);
kgt = p(7);
dt = p(8);
T0 = p(9);

%  nct = cancer-mediated expansion rate of T cells
%  smt = inhibition rate of T cells by MDSCs
%  kgt = gemcitabine kill rate of T cells
%  dt  = T cell natural death rate
%  T0  = homeostatic T cell population


% --------------- Parameter values for MDSC equation ---------------------- 
rcm = p(10); 
kgm = p(11); 
stm = p(12);
dm = p(13);
M0 = p(14);

%  rcm = cancer-mediated recruitment of MDSCs
%  kgm = gemcitabine kill rate of MDSCs
%  stm = MDSC suppression by T cells
%  dm  = MDSC natural death rate
%  M0  = homeostatic MDSC population


% --------------- Parameter values for gemcitabine equation --------------- 
dg = p(15);

%  dg = gemcitabine decay rate


% --------------- Parameter values for all cell types ---------------------
km = p(16); 

%  km = Michaelis-Menten constant



% -------------------------------------------------------------------------

%  OUTPUTS:
dxdt(1) = pc*C*(1-(C/cmax)) - ktc*T*C - kgc*C*(G/(G+km));
dxdt(2) = nct*C*T -smt*M*T -kgt*T*(G/(G+km)) -dt*T +(smt*M0 +dt)*T0;
dxdt(3) = rcm*C - kgm*M*(G/(G+km)) -stm*T*M -dm*M +(stm*T0 +dm)*M0;
dxdt(4) = -dg*G;
