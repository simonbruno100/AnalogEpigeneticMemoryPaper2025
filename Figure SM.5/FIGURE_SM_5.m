
%% CODE TO GENERATE FIGURE SM.5
% This code generates input-dependent distributions of the variable n_x 
% (protein output) from stochastic simulations under different chromatin input types (e.g., KRAB/DNMT3).
% Simulations are conducted using the SSA_FIG_SM_5 function.
% The results are visualized as logicle-scaled histograms.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions

a=0; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

TTT=1000; % Total number N of simulations per distribution

NumberOfDays=1; % Time interval (in days) between two distributions

TTTTT=13; % Total number of distributions to generate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set general biological parameters and initial conditions


Dt=15; % Total number of nucleosomes in the gene of interest

p.V=1; % Reaction volume (used to normalize concentrations)

tfinn=24*NumberOfDays*TTTTT; % Total simulation time (in hours)

%%% Define input range for simulations

% Select Input ranges for the simulation:
%KRAB: [1.5, 7.5], [20, 100], [45, 225], [4000, 20000]
%DNMT3: [0.2, 0.7], [0.35, 1.23], [0.6, 2.1], [0.9, 3.1]


l1=20;
l2=100;

rrrrtotal = unifrnd(l1,l2,1,TTT);

% Select either one of the two values below
% DAPG=0;  % In case the external input is KRAB
DAPG=50; % In case the external input is DNMT3


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

% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_5(DA0,DR10,DR20,DR120,tfinn,DAPG,rrrrtotal(i));

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

iindex = find(TTfinn>tend);  % Find first time exceeding tend

index=min(iindex);


% Store output variables at selected time point
ddaaf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,5)/p.V;
dd11f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,2)/p.V;
dd22f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,3)/p.V;
dd1122f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,4)/p.V;
mmxxf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,6)/p.V;
xf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,7)/p.V;

end

 end


% Store final variable of interest

xf2=xf; %n_x

% Pad for logicle transform bounds
x2=[xf2;...
    zeros(1,TTTTT);...
    (10^(4))*ones(1,TTTTT)];


%% Generate plots
% Generate logicle-scaled histograms of n_x for each timepoint

% Figure size

width =  101;
height = 31;

% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 



 for iji=1:1:TTTTT;



% Apply logicle transform (cytometry-style scaling)
obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
y = obj.transform(x2(:, iji));
figure(iji)
hold on
box on
% Plot histogram of logicle-transformed values
histogram(y, 60, 'Normalization', 'probability');
axis ([0.0592 0.5282 0 0.4]);

% Set figure size
fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);

% Format axes
ax = gca;
ax.XTick = obj.Tick;
ax.XTickLabel = []; % Hide X labels

set(ax, 'FontName', fontName, 'FontSize', fontSize);

ay = gca;
ax.YTickLabel = []; % Hide Y labels

 end

