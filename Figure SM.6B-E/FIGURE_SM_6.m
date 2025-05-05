%% CODE TO GENERATE DATA TO REALIZE FIGURE SM.6B-E
% This code generates input-dependent distributions of our variables of
% interested from stochastic simulations. These simulations will be used to generate Figure SM.6B and as initial conditions
% for simulations in FIGURE SM.6C-E

% Simulations are conducted using the SSA_FIG_SM_6 function.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions

a=0; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

TTT=1000; % Total number N of simulations per distribution

NumberOfDays=30; % Time interval (in days) between two distributions

TTTTT=1; % Total number of distributions to generate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set general biological parameters and initial conditions


Dt=15; % Total number of nucleosomes in the gene of interest

p.V=1; % Reaction volume (used to normalize concentrations)

tfinn=24*NumberOfDays*TTTTT; % Total simulation time (in hours)


%%% Define input range for simulations
% Input (DNMT3) range = [0.6, 2.1]

l1=0.6;
l2=2.1;

rrrrtotal = unifrnd(l1,l2,1,TTT);

DAPG=0; % DAPG = 0 since the external input in this case is DNMT3


% Initialize matrices to collect final simulation values (for each distribution)

ddaaf=zeros(1,TTTTT);% DA
dd11f=zeros(1,TTTTT);% DR1
dd22f=zeros(1,TTTTT);% DR2
dd1122f=zeros(1,TTTTT);% DR12
mmxxf=zeros(1,TTTTT);% mRNA
xf=zeros(1,TTTTT);% Protein


 for i=1:1:TTT;

% Initial condition: Full active state
DR120=0;
DA0=15;
DR10=0;
DR20=0;
m0=0.5326*(DA0+(DR10+DR20+DR120)*0.0417);
p0=144.3299*m0;


% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_6_1(DA0,DR10,DR20,DR120,m0,p0,tfinn,DAPG,rrrrtotal(i));

save('reactions0o0','X0o0','T0o0');  

% Optional: Display current simulation index
    display('-------------------------------------');
        disp(i);
    display('-------------------------------------');


% Loop over each distribution (at each target time)
for iji=1:1:TTTTT;

tend=24*iji*NumberOfDays-1; % Time point of interest (in hours)

load('Tfin.mat')

TTfinn=round(Tfin);

iindex = find(TTfinn>tend); % Find first time exceeding tend

index=min(iindex);


ddaaf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,5)/p.V;
dd11f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,2)/p.V;
dd22f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,3)/p.V;
dd1122f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,4)/p.V;
mmxxf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,6)/p.V;
xf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,7)/p.V;

end

 end


% Store output variables of interest


xf2=xf;



max_value_da = max(ddaaf(:));
max_value_d1 = max(dd11f(:));
max_value_d2 = max(dd22f(:));
max_value_d12 = max(dd1122f(:));
max_value_m  = max(mmxxf(:));


 ddaa=[ddaaf;...
    zeros(1,TTTTT);...
    Dt*ones(1,TTTTT)];


 dd11=[dd11f;...
    zeros(1,TTTTT);...
    Dt*ones(1,TTTTT)/2];


 dd22=[dd22f;...
    zeros(1,TTTTT);...
    Dt*ones(1,TTTTT)];


 dd1122=[dd1122f;...
    zeros(1,TTTTT);...
    Dt*ones(1,TTTTT)/2];


 mmxx=[mmxxf;...
    zeros(1,TTTTT);...
    max_value_m*ones(1,TTTTT)];


x2=[xf2;...
    zeros(1,TTTTT);...
    (10^(4))*ones(1,TTTTT)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Save variables of interest
% 
% var1=dd11(1:end-2,:);
% var2=dd22(1:end-2,:);
% var3=dd1122(1:end-2,:);
% var4=ddaa(1:end-2,:);
% var5=mmxx(1:end-2,:);
% var6=x2(1:end-2,:);
% 
% 
% save('initialdistributionD1','var1');
% save('initialdistributionD2','var2');
% save('initialdistributionD12','var3');
% save('initialdistributionDA','var4');
% save('initialdistributionmx','var5');
% save('initialdistributionx2','var6');


