%% CODE TO GENERATE FIGURE SM.7
% This code generates input-dependent distributions of the variable
% n_{Y1}/Dtot = Fraction of CpG methylated from stochastic simulations
% under different parameter conditions

% Simulations are conducted using the SSA_FIG_SM_7 function.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions

a=0; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

TTT=1000; % Total number N of simulations per distribution

NumberOfDays=30; % Time interval (in days) between two distributions

TTTTT=4; % Total number of distributions to generate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set general biological parameters and initial conditions


Dt=15; % Total number of nucleosomes in the gene of interest

p.V=1; % Reaction volume (used to normalize concentrations)

tfinn=24*NumberOfDays*TTTTT; % Total simulation time (in hours)


%%% No External input -->l1=l2=0

l1=0;
l2=0;

rrrrtotal = unifrnd(l1,l2,1,TTT);


% Initialize matrices to collect final simulation values (for each distribution)

ddaaf=zeros(1,TTTTT);% DA
dd11f=zeros(1,TTTTT);% DR1
dd22f=zeros(1,TTTTT);% DR2
dd1122f=zeros(1,TTTTT);% DR12
mmxxf=zeros(1,TTTTT);% mRNA
xf=zeros(1,TTTTT);% Protein



 for i=1:1:TTT;



% Select initial condition. In the paper we consider DR120 = 0,4,6,14

DR120=6;
DA0=Dt-DR120;
DR10=0;
DR20=0;


% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_7(DA0,DR10,DR20,DR120,tfinn,rrrrtotal(i));

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
% Generate histograms of n_{Y1} for each time point of interest

% Figure size
width =  105;
height = 52;

% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 


 for iji=1:1:TTTTT;


UUUUU = Dt;                         % Number of bins = total nucleosomes


% Histogram
[counts, binEdges] = histcounts(dd11(:,iji)+dd1122(:,iji), UUUUU, 'Normalization', 'probability');
binWidth = binEdges(2) - binEdges(1);
binMidpoints = binEdges(1:end-1) + binWidth / 2;

% Plot smooth line version of the histogram
figure(iji)
hold on;
plot(binMidpoints, counts, 'LineWidth', 2);
hold off;
box on
axis ([0+0.5 Dt-0.5 0 1]);

fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);


ax = gca;
ax.XTickLabel = [];
set(ax, 'FontName', fontName, 'FontSize', fontSize);


ay = gca;
ax.YTickLabel = [];

end
