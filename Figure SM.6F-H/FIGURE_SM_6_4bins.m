%% CODE TO GENERATE FIGURE SM.6F-H
% This code generates  distributions of the variable n_x (protein output) from stochastic simulations under different binning conditions (e.g., bin1,bin2,bin3,bin4).

% Simulations are conducted using the SSA_FIG_SM_6_2 function.

% Results are visualized as logicle-scaled histograms or violin plots.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part of the code to generate distributions in Figure SM.6F

% Load intermediate distribution from previous simulations

 load('IntermediateD10');      % var1 = DR1
 load('IntermediateD220');     % var2 = DR2
 load('IntermediateD120');     % var3 = DR12
 load('IntermediateDA0');      % var4 = DRA
 load('Intermediatem0');       % var5 = mRNA
 load('IntermediateX0');       % var6 = protein (n_x)


% Assign short names to loaded variables for convenience
 dd11=xyxyxyfdd11f;
 dd22=xyxyxyfdd22f;
 dd1122=xyxyxyfdd1122f;
 ddaa=xyxyxyfddaaf;
 mmxx=xyxyxyfmmxxf;
 x2=xyxyxyfxf2;

clear xxxxf
clear xxxx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions 

FacsDayForBin=size(ddaa,2); % Time index to sample for binning

SSize= size(ddaa,1);        % Total number of samples

j=0;

% Filter to generate distribution for bin1,bin2,bin3,bin4 gene exrpression level distributions

 for i=1:1:SSize-2

% Select the bin range depending on the gene expression level of interest
% Bin1: [0, 80]
% Bin2: [90, 150]
% Bin3: [250, 400]
% Bin4: [500, 10000]

if x2(i,FacsDayForBin)>90
    if x2(i,FacsDayForBin)<150

        j=j+1;
        xxxxf(j,1:6)=[ddaa(i,FacsDayForBin) ...   %Da
                      dd11(i,FacsDayForBin) ...   %D1
                      dd22(i,FacsDayForBin) ...   %D2
                      dd1122(i,FacsDayForBin) ... %D12
                      mmxx(i,FacsDayForBin) ...   %mx
                      x2(i,FacsDayForBin)];       %X
        
    end
end

 end

% Pad for logicle transform bounds
xxxx=[xxxxf;...
    zeros(1,6);...
    (10^(4))*ones(1,6)];


%% Generate plots
% Plot initial distribution of n_x (Figure panel F)

% Figure size

width =  109;
height = 66.5;

% Set the line thickness for the plot 
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc. 
fontName = 'Arial';
fontSize = 7; 

% Apply logicle transform (cytometry-style scaling)
obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
yyyy = obj.transform(xxxx(:,6));
figure(1)
hold on
box on
% Plot histogram of logicle-transformed values
histogram(yyyy, 60, 'Normalization', 'probability');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part of the code to generate distributions and violin plots in Figure SM.6G, H

% Clear previous SSA results
clear xyxyxyfxf2 xyxyxyf xyxyxyfddaaf xyxyxyfdd11f xyxyxyfdd22f xyxyxyfdd1122f xyxyxyfmmxxf xyxyxyfxf
clear xyxyxyfx2 xyxyxy xyxyxyfddaa xyxyxyfdd11 xyxyxyfdd22 xyxyxyfdd1122 xyxyxyfmmxx xyxyxyfx

close all

%%% Optional: Increase the initial pool if needed
xxxxf=[xxxxf;...
       xxxxf;...
       xxxxf];

%%% Cap number of simulations
SSSSize= size(xxxxf,1); 

TTT2=1000;

if SSSSize < 1000

TTT2=SSSSize;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions and  general biological parameters

a=10; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

Dt=15; % Total number of nucleosomes in the gene of interest

p.V=1; % Reaction volume (used to normalize concentrations)


NumberOfDays2=2; % Time interval (in days) between two distributions


TTTTT2=7; % Total number of distributions to generate


tfinn2=24*NumberOfDays2*TTTTT2; % Total simulation time (in hours)


% No external input
rrrrrtotal = zeros(1,SSSSize);


 for i=1:1:TTT2;


 j=randi(SSSSize);
 
% Initialize matrices to collect final simulation values (for each distribution)

% Set initial conditions
DA0=xxxxf(j,1);
DR10=xxxxf(j,2);
DR20=xxxxf(j,3);
DR120=xxxxf(j,4);
m0=xxxxf(j,5);
p0 = xxxxf(j,6);

DAPG=0; % DAPG = 0 since no external input


% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_6_2(DA0,DR10,DR20,DR120,m0,p0,tfinn2,DAPG,rrrrrtotal(i));

save('reactions0o0','X0o0','T0o0');  

% Optional: Display current simulation index
    display('-------------------------------------');
        disp(i);
    display('-------------------------------------');



for iji=1:1:TTTTT2;




tend2=24*iji*NumberOfDays2-1;

load('Tfin.mat')

TTfinn2=round(Tfin);

iindex2 = find(TTfinn2>tend2);

index2=min(iindex2);

% Store output variables at selected time point
xyxyxyfddaaf(a*(i-1)+i:a*i+i,iji)=X0o0(index2-a:index2,5)/p.V;
xyxyxyfdd11f(a*(i-1)+i:a*i+i,iji)=X0o0(index2-a:index2,2)/p.V;
xyxyxyfdd22f(a*(i-1)+i:a*i+i,iji)=X0o0(index2-a:index2,3)/p.V;
xyxyxyfdd1122f(a*(i-1)+i:a*i+i,iji)=X0o0(index2-a:index2,4)/p.V;
xyxyxyfmmxxf(a*(i-1)+i:a*i+i,iji)=X0o0(index2-a:index2,6)/p.V;
xyxyxyfxf(a*(i-1)+i:a*i+i,iji)=X0o0(index2-a:index2,7)/p.V;


end

 end

xyxyxyfxf2=xyxyxyfxf;



xyxyxyfx2=[xyxyxyfxf2;...
    zeros(1,TTTTT2);...
    (10^(4))*ones(1,TTTTT2)];




%% Generate plots
% Figure panel H – Plot distribution of n_x at the last time point

% Figure size

width =  87;
height = 50;

% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 


iji = TTTTT2;                       % Final time point index



obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
% Transform and plot histogram of final gene expression level
xyxyxyfy = obj.transform(xyxyxyfx2(:, iji));

figure(20)
hold on
box on
histogram(xyxyxyfy, 60, 'Normalization', 'probability');
axis ([0.0592 0.5282 0 0.4]);

fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);

ax = gca;
ax.XTick = obj.Tick;
ax.XTickLabel = [];

set(ax, 'FontName', fontName, 'FontSize', fontSize);

ay = gca;
ax.YTickLabel = [];


%% Generate Violin plots
% Figure panel G – Violin plots of n_x over time (logicle scaled)


xtot=xyxyxyfxf2;            % Matrix of protein values across time
ijij = 1:1:TTTTT2;          % Time index vector for violin plots


% Figure size
width =  101*5;
height = 31*5;


% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 

% Apply logicle transform (cytometry-style scaling)
obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);


xyxyxyfyTOT = obj.transform(xtot(1:end-2,:));

figure(35)
% Create violin plots for each time point
violinplot(xyxyxyfyTOT,ijij,'ViolinColor',[0.48,0.84,1]); %blue
box on
axis ([0 TTTTT2+1 0.0592 0.5282]);

fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);




ax = gca;
ax.YTick = obj.Tick;

ax.YTickLabel = []; % Hide Y labels

set(ax, 'FontName', fontName, 'FontSize', fontSize);

ax = gca;
ax.XTickLabel = [];% Hide X labels

