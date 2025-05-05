%% CODE TO GENERATE FIGURE 4H, J
% This code generates input-dependent distributions of the variable
% n_{Y1}/Dtot = Fraction of CpG methylated from stochastic simulations
% under different parameter conditions, and evaluate the corresponsing gene
% expression level

% Simulations are conducted using the SSA_FIG_4 function.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions

a=0; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

TTT=1000; % Total number N of simulations per distribution

NumberOfDays=9; % Time interval (in days) between two distributions

TTTTT=9; % Total number of distributions to generate


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


% Select initial condition. In the paper we consider DR120 = 13 (Bin 1), 12
% (Bin 2), 8 (Bin 3), 6 (Bin 4), 5 (Bin 5), 4 (Bin 6), 1 (Bin 7), 0 (Bin 8)

DR120=12;
DA0=Dt-DR120;
DR10=0;
DR20=0;


% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_4(DA0,DR10,DR20,DR120,tfinn,rrrrtotal(i));


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



 for iji=[7 9];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Correlation between gene expression and DNA methylation
% % This section extracts and saves the data required to generate the correlation
% % plots shown in Figure 4J, which link CpG methylation to gene expression.
% %
% % Specifically, for each selected time point (e.g., time point 1 = day 63, 
% % time point 2 = day 81), we:
% % - Compute the mean and standard deviation of DNA methylation levels (n_{Y1})
% %   defined as DR1 + DR12
% % - Compute the mean and standard deviation of gene expression levels (n_x)
% % - Save both summary statistics and full gene expression vectors (xf2) 
% %   to be used in an external plotting script (`Code_FIG4I.m`).
% %
% % IMPORTANT:
% % - You must manually replace `_0_` in the filenames below with the appropriate bin number
% %   corresponding to the value of DR120 used at initialization.
% 
% 
% iji = 7;
% 
% DNAMetTot = dd11f(:,iji)+dd1122f(:,iji);
% 
% meanDNAmeth = mean(DNAMetTot);
% 
% stDNAmeth = std(DNAMetTot);
% 
% 
% meanX= mean(xf2(:,iji));
% 
% stX=std(xf2(:,iji));
% 
% valuex=xf2(:,iji);
% 
% 
% save('TimePoint1_DNAmeth_0_mean','meanDNAmeth');
% save('TimePoint1_DNAmeth_0_st','stDNAmeth');
% save('TimePoint1_DNAmeth_0_X_mean','meanX');
% save('TimePoint1_DNAmeth_0_X_st','stX');
% save('values_TimePoint1_DNAmeth_0_X','valuex');
% 
% 
% clear DNAMetTot 
% 
% iji = 9;
% 
% DNAMetTot = dd11f(:,iji)+dd1122f(:,iji);
% 
% meanDNAmeth = mean(DNAMetTot);
% 
% stDNAmeth = std(DNAMetTot);
% 
% 
% meanX= mean(xf2(:,iji));
% 
% stX=std(xf2(:,iji));
% 
% 
% valuex2=xf2(:,iji);
% 
% 
% save('TimePoint2_DNAmeth_0_mean','meanDNAmeth');
% save('TimePoint2_DNAmeth_0_st','stDNAmeth');
% save('TimePoint2_DNAmeth_0_X_mean','meanX');
% save('TimePoint2_DNAmeth_0_X_st','stX');
% save('values_TimePoint2_DNAmeth_0_X','valuex2');
% 
% 
