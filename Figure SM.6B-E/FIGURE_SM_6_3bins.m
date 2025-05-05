%% CODE TO GENERATE FIGURE SM.6B-E
% This code generates input-dependent distributions of the variable n_x
% (protein output) from stochastic simulations under different binning conditions (e.g., low, medium, high expression).

% Simulations are conducted using the SSA_FIG_SM_6_1 function.

% Results are visualized as logicle-scaled histograms or violin plots.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part of the code to generate distributions in Figure SM.6B

% Load intermediate distribution from previous simulations

load('initialdistributionD1.mat')     % var1 = DR1
load('initialdistributionD2.mat')     % var2 = DR2
load('initialdistributionD12.mat')    % var3 = DR12
load('initialdistributionDA.mat')     % var4 = DA
load('initialdistributionmx.mat')     % var5 = mRNA
load('initialdistributionx2.mat')     % var6 = protein (n_x)


clear xxxxf
clear xxxx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions 

FacsDayForBin=size(var4,2); % Time index to sample for binning

SSize= size(var4,1);        % Total number of samples

 j=0;

 % Filter to generate distribution for low, intermediate and high gene expression level
 for i=1:1:SSize-2


% Select the bin range depending on the gene expression level of interest
% Low:          [0, 100]
% Intermediate: [200, 800]
% High:         [900, 10000]

if var6(i,FacsDayForBin)>0
    if var6(i,FacsDayForBin)<100

        j=j+1;
        xxxxf(j,1:6)=[var4(i,FacsDayForBin) ... %Da
                      var1(i,FacsDayForBin) ... %D1
                      var2(i,FacsDayForBin) ... %D2
                      var3(i,FacsDayForBin) ... %D12
                      var5(i,FacsDayForBin) ... %mx
                      var6(i,FacsDayForBin)];   %X

    end
end

 end


% Pad for logicle transform bounds
xxxx=[xxxxf;...
    zeros(1,6);...
    (10^(4))*ones(1,6)];


%% Generate plots
% Plot initial distribution of n_x (Figure panel B)

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
%% Part of the code to generate distributions and violin plots in Figure SM.6C-E


% Clear previous SSA results
clear xyxyxyfxf2 xyxyxyf xyxyxyfddaaf xyxyxyfdd11f xyxyxyfdd22f xyxyxyfdd1122f xyxyxyfmmxxf xyxyxyfxf
clear xyxyxyfx2 xyxyxy xyxyxyfddaa xyxyxyfdd11 xyxyxyfdd22 xyxyxyfdd1122 xyxyxyfmmxx xyxyxyfx


%%% Optional: Increase the initial pool if needed
% xxxxf=[xxxxf;...
%        xxxxf;...
%        xxxxf];

%%% Cap number of simulations
SSSSize= size(xxxxf,1); 

TTT2=1000;

if SSSSize < 1000

TTT2=SSSSize;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions and  general biological parameters

a=0; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

Dt=15; % Total number of nucleosomes in the gene of interest

p.V=1; % Reaction volume (used to normalize concentrations)


NumberOfDays2=1; % Time interval (in days) between two distributions


TTTTT2=29; % Total number of distributions to generate


tfinn2=24*NumberOfDays2*TTTTT2; % Total simulation time (in hours)


% No external input
rrrrrtotal = zeros(1,SSSSize);


 for i=1:1:TTT2;


 j=randi(SSSSize);
 
% Initialize matrices to collect final simulation values (for each distribution)

DR120=xxxxf(i,4);% DR12
DA0=xxxxf(i,1);% DA
DR10=xxxxf(i,2);% DR1
DR20=xxxxf(i,3);% DR2
m0=xxxxf(i,5);% mRNA
p0=xxxxf(i,6);% Protein


DAPG=0; % DAPG = 0 since no external input


% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_6_1(DA0,DR10,DR20,DR120,m0,p0,tfinn2,DAPG,rrrrrtotal(i));

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



sizXf=size(xyxyxyfxf,1);



for j = 1:1:TTTTT2;


for i = 1:1:sizXf;

    if xyxyxyfxf(i,j)<1
        xyxyxyfxf(i,j)=1;

    end   

end

end



max_value_xyda = max(xyxyxyfddaaf(:));
max_value_xyd1 = max(xyxyxyfdd11f(:));
max_value_xyd2 = max(xyxyxyfdd22f(:));
max_value_xyd12 = max(xyxyxyfdd1122f(:));
max_value_xym  = max(xyxyxyfmmxxf(:));


 xyxyxyfddaa=[xyxyxyfddaaf;...
    zeros(1,TTTTT2);...
    Dt*ones(1,TTTTT2)];


 xyxyxyfdd11=[xyxyxyfdd11f;...
    zeros(1,TTTTT2);...
    Dt*ones(1,TTTTT2)/2];


 xyxyxyfdd22=[xyxyxyfdd22f;...
    zeros(1,TTTTT2);...
    Dt*ones(1,TTTTT2)];


 xyxyxyfdd1122=[xyxyxyfdd1122f;...
    zeros(1,TTTTT2);...
    Dt*ones(1,TTTTT2)/2];


 xyxyxyfmmxx=[xyxyxyfmmxxf;...
    zeros(1,TTTTT2);...
    max_value_xym*ones(1,TTTTT2)];


xyxyxyfx2=[xyxyxyfxf2;...
    zeros(1,TTTTT2);...
    (10^(4))*ones(1,TTTTT2)];
% 


%% Generate plots
% Plot final distribution of n_x (Figure panel D) and corresponding DNA methylation levels (Figure panel E)

%%% Figure panel E – Fraction of CpG methylated distribution


% Figure size
width =  109;
height = 66.5;

% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 


iji = TTTTT2;                       % Index for final time point
UUUUU = Dt;                         % Number of bins = total nucleosomes

% Histogram
 histObj = histogram(xyxyxyfdd11(:,iji)+xyxyxyfdd1122(:,iji),UUUUU,'Normalization','probability');

 % 
% Compute midpoints of bins for plotting
binEdges = histObj.BinEdges;
binWidth = histObj.BinWidth;
binMidpoints = binEdges(1:end-1) + binWidth / 2;

% Plot smooth line version of the histogram
figure(222)
hold on;
plot(binMidpoints, histObj.Values, 'LineWidth', 1.5);
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



%%% Figure panel D – Final protein (n_x) distribution
% This panel shows the final distribution of protein levels (n_x) at the last time point.


% Figure size

width =  87;
height = 50;

% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 


iji = TTTTT2;                       % Final time point index


obj = logicleTransform(262144.0,0.2781913054120009,4.5,0); % Logicle scale transformation

% Transform and plot histogram of final gene expression level
xyxyxyfy = obj.transform(xyxyxyfx2(:, iji));


figure(2)
hold on
box on
histogram(xyxyxyfy, 60, 'Normalization', 'probability');
axis ([0.0592 0.5282 0 0.4]);

ax = gca;
ax.XTick = obj.Tick;
ax.XTickLabel = [];

set(ax, 'FontName', fontName, 'FontSize', fontSize);



ay = gca;
ax.YTickLabel = [];


%% Generate Violin plots
% Figure panel C – Violin plots of n_x over time (logicle scaled)

xtot = xyxyxyfx2;           % Matrix of protein values across time
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

figure(34)
% Create violin plots for each time point
% (Color can be changed to distinguish between conditions: blue/orange/red)
% violinplot(xyxyxyfyTOT,ijij,'ViolinColor',[0.48,0.84,1]); %blue
violinplot(xyxyxyfyTOT,ijij,'ViolinColor',[1,0.82,0.62]); %orange
% violinplot(xyxyxyfyTOT,ijij,'ViolinColor',[1,0.73,0.72]); %red
box on
axis ([0 TTTTT2+1 0.0592 0.5282]);

fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);


ax = gca;
ax.YTick = obj.Tick;
ax.YTickLabel = []; % Hide Y labels

set(ax, 'FontName', fontName, 'FontSize', fontSize);

ax = gca;
ax.XTickLabel = []; % Hide X labels

