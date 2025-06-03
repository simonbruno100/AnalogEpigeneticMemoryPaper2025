
%% CODE TO GENERATE FIGURE SM.3C, SM.3E, SM.1E and 2J 
% This code generates input-dependent distributions of the variable n_x 
% (protein output) and n_{Y1} (Total CpGme) from stochastic simulations under different chromatin input types (e.g., KRAB/DNMT3).
% Simulations are conducted using the SSA_FIG_SM_3 function.
% The results are visualized as logicle-scaled histograms.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions

a=0; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

TTT=1000; % Total number N of simulations per distribution

NumberOfDays=28; % Time interval (in days) between two distributions

TTTTT=1; % Total number of distributions to generate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set general biological parameters and initial conditions


Dt=15; % Total number of nucleosomes in the gene of interest

p.V=1; % Reaction volume (used to normalize concentrations)

tfinn=24*NumberOfDays*TTTTT; % Total simulation time (in hours)


%%% Define input range for simulations

% Choose the appropriate input interval [l1, l2] for the experiment of interest.

% Values used for simulations in the paper:
% DNMT3 in SM.3C (left panel):     [0, 0.7], [0.5, 1.2], [1, 1.7]
% DNMT3 in SM.3C (right panel):    [0, 0.4], [2.4, 2.8], [3.4, 3.8]
% DNMT3 in SM.1E:                  [0, 0.4], [2.4, 2.8], [3.4, 3.8]
% DNMT3 in 2J:                     [0.5, 1.2], [1, 1.7], [1.5, 2.2]
% TET1 in SM.3E:                   [0, 0.05], [0.24, 0.28]


l1=1;
l2=1.7;

rrrrtotal = unifrnd(l1,l2,1,TTT);


% Initialize matrices to collect final simulation values (for each distribution)

ddaaf=zeros(1,TTTTT);% DA
dd11f=zeros(1,TTTTT);% DR1
dd22f=zeros(1,TTTTT);% DR2
dd1122f=zeros(1,TTTTT);% DR12
mmxxf=zeros(1,TTTTT);% mRNA
xf=zeros(1,TTTTT);% Protein


 for i=1:1:TTT;

% Select initial condition: uncomment only one of the two blocks below

% For DNMT3A experiment (Figures SM.3C, SM.1E, and 2.J)
DR120=0;
DA0=15;
DR10=0;
DR20=0;

% % For TET experiment (Figure SM.3E)
% DR120=15;
% DA0=0;
% DR10=0;
% DR20=0;



% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_3(DA0,DR10,DR20,DR120,tfinn,rrrrtotal(i));

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


ddaaf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,5)/p.V;
dd11f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,2)/p.V;
dd22f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,3)/p.V;
dd1122f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,4)/p.V;
mmxxf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,6)/p.V;
xf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,7)/p.V;

end

 end


% Store final variables of interest

xf2=xf; %n_x



% Pad for logicle transform bounds
 dd11=[dd11f;...
    zeros(1,TTTTT);...
    Dt*ones(1,TTTTT)/2];

 dd1122=[dd1122f;...
    zeros(1,TTTTT);...
    Dt*ones(1,TTTTT)/2];

x2=[xf2;...
    zeros(1,TTTTT);...
    (10^(4))*ones(1,TTTTT)];


%% Generate plots
% (1) Generate logicle-scaled histograms of n_x for final time point


% Figure size
width =  76;
height = 52;

% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 



iji=TTTTT;

% Apply logicle transform (cytometry-style scaling)
obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
y = obj.transform(x2(:, iji));
figure(iji)
hold on
box on
% Plot histogram of logicle-transformed values
histogram(y, 60, 'Normalization', 'probability');
axis ([0.0592 0.5282 0 0.3]);

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



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % (2) Generate histograms of n_{Y1} for final time point (Figure SM.1E -
% % RHS panel)
% 
% 
% 
% iji=TTTTT;
% 
% figure(30)
% hold on
% histogram(dd11(:,iji)+dd1122(:,iji),Dt,'Normalization','probability');
% set(gca, 'FontName', 'Times New Roman','FontSize',18)
% box on
% axis ([0 Dt 0 0.6]);
% 
% % Set figure size
% fig=gcf;
% set(fig,'Units', 'points', 'Position', [0, 0, width, height]);
% 
% % Format axes
% ax = gca;
% ax.XTick = obj.Tick;
% ax.XTickLabel = []; % Hide X labels
% 
% set(ax, 'FontName', fontName, 'FontSize', fontSize);
% 
% ay = gca;
% ax.YTickLabel = []; % Hide Y labels
