
%% CODE TO GENERATE FIGURE SM.8
% This code generates input-dependent distributions of the variable n_x 
% (protein output) from stochastic simulations under different chromatin input types (e.g., KRAB/DNMT3/TET1).
% Initial conditions (intermediate distribution) are sampled from previous
% simulations
% Simulations are conducted using the SSA_FIG_SM_8 function.
% The results are visualized as logicle-scaled histograms.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load and prepare initial distribution data

load('IntermediateD10');   % DR10
load('IntermediateD220');  % DR20
load('IntermediateD120');  % DR120
load('IntermediateDA0');   % DA0
load('Intermediatem0');    % mRNA
load('IntermediateX0');    % Protein (n_x)



% Rename loaded variables to unified names
dd11   = xyxyxyfdd11f;     % DR1
dd22   = xyxyxyfdd22f;     % DR2
dd1122 = xyxyxyfdd1122f;   % DR12
ddaa   = xyxyxyfddaaf;     % DA
mmxx   = xyxyxyfmmxxf;     % mRNA
x2     = xyxyxyfxf2;       % Protein (n_x)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions

clear xxxxf
clear xxxx

a=0; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

TTT=1000; % Total number N of simulations per distribution

Dt=15;  % Total number of nucleosomes

FacsDayForBin = 7; % Day of interest for selecting data (e.g., FACS time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter initial conditions based on n_x range

 SSize= size(ddaa,1);

 j=0;
 for i=1:1:SSize-2

if x2(i,FacsDayForBin)>0
    if x2(i,FacsDayForBin)<10000

        j=j+1;
        xxxxf(j,1:6)=[ddaa(i,FacsDayForBin) ...
                      dd11(i,FacsDayForBin) ...
                      dd22(i,FacsDayForBin) ...
                      dd1122(i,FacsDayForBin) ...
                      mmxx(i,FacsDayForBin) ...
                      x2(i,FacsDayForBin)];

    end
end




 end



% Pad for logicle transform range
xxxx=[xxxxf;...
    zeros(1,6);...
    (10^(4))*ones(1,6)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Plot histogram of initial n_x distribution (logicle scale)
% 
% 
% obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
% yyyy = obj.transform(xxxx(:,6));
% figure(1)
% hold on
% box on
% histogram(yyyy, 60, 'Normalization', 'probability');
% axis ([0.0592 0.5282 0 0.4]);
% 
% fig=gcf;
% set(fig,'Units', 'points', 'Position', [0, 0, width, height]);
% 
% 
% 
% 
% ax = gca;
% ax.XTick = obj.Tick;
% %ax.XTickLabel = obj.TickLabel;
% ax.XTickLabel = [];
% 
% set(ax, 'FontName', fontName, 'FontSize', fontSize);
% 
% 
% 
% ay = gca;
% ax.YTickLabel = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run stochastic simulations for different inputs (KRAB, DNMT3, TET1)

clear xyxyxyfxf2 xyxyxyf xyxyxyfddaaf xyxyxyfdd11f xyxyxyfdd22f xyxyxyfdd1122f xyxyxyfmmxxf xyxyxyfxf
clear xyxyxyfx2 xyxyxy xyxyxyfddaa xyxyxyfdd11 xyxyxyfdd22 xyxyxyfdd1122 xyxyxyfmmxx xyxyxyfx



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions

SSSSize= size(xxxxf,1);

TTT2=TTT;

if SSSSize < TTT;

TTT2=SSSSize;


end



a=0;

NumberOfDays2=2; % Time interval (in days) between two distributions

TTTTT2=7; % Total number of distributions to generate

tfinn2=24*NumberOfDays2*TTTTT2;  % Total simulation time (in hours)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set biological parameters


Dt=15;

p.V=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define input range for simulations

% Select Input ranges for the simulation:
%50,250       KRAB
%0.2,0.7       DNMT3
%0.004, 0.014  TET1


l1=0.004;
l2=0.014;

rrrrrtotal = unifrnd(l1,l2,1,TTT2);  % Input sampled uniformly


% Main simulation loop
 for i=1:1:TTT2;

 % Randomly sample an initial condition

 j=randi(SSSSize);

DA0=xxxxf(j,1);
DR10=xxxxf(j,2);
DR20=xxxxf(j,3);
DR120=xxxxf(j,4);
m0=xxxxf(j,5);
p0 = xxxxf(j,6);

DAPG=0; % no DAPG in this case


% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_8(DA0,DR10,DR20,DR120,m0,p0,tfinn2,DAPG,rrrrrtotal(i));

save('reactions0o0','X0o0','T0o0');  


% Optional: Display current simulation index
    display('-------------------------------------');
        disp(i);
    display('-------------------------------------');



% Loop over each distribution (at each target time)
for iji=1:1:TTTTT2;




tend2=24*iji*NumberOfDays2-1;  % Time point of interest (in hours)

load('Tfin.mat')

TTfinn2=round(Tfin);

iindex2 = find(TTfinn2>tend2); % Find first time exceeding tend

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



% Store final variable of interest

xyxyxyfxf2=xyxyxyfxf; %n_x


xyxyxyfx2=[xyxyxyfxf2;...
    zeros(1,TTTTT2);...
    (10^(4))*ones(1,TTTTT2)];




%% Generate plots
% Generate logicle-scaled histograms of n_x for each timepoint


% Figure size

width =  101;
height = 60;

% Set the line thickness for the plot 
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc. 
fontName = 'Arial';
fontSize = 7; 



for iji=1:1:TTTTT2;

     


ijk=iji;

% Apply logicle transform (cytometry-style scaling)
obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
xyxyxyfy = obj.transform(xyxyxyfx2(:, iji));
figure(ijk)
hold on
box on
% Plot histogram of logicle-transformed values
histogram(xyxyxyfy, 60, 'Normalization', 'probability','FaceColor',[0.93,0.24,0.13]);
axis ([0.0592 0.5282 0 0.4]);

% Set figure size
fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);


% Format axes
ax = gca;
ax.XTick = obj.Tick;
ax.XTickLabel = [];

set(ax, 'FontName', fontName, 'FontSize', fontSize);

ay = gca;
ax.YTickLabel = [];

 end



% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% (Optional) Plot precomputed distributions for low (and optionally high) gene expression
% 
% 
% 
% load('LowX');   % Load precomputed distributions for low gene expression level state
% % load('HighX'); % Load precomputed distributions for high gene expression level state
% 
% 
% %%% Define figure properties
% 
% % Figure size
% width =  101;
% height = 60;
% 
% % Set the line thickness for the plot 
% lineThickness = 0.5;
% 
% % Set the font size for axis labels, titles, etc. 
% fontName = 'Arial';
% fontSize = 7; 
% 
% 
% % Loop through each time point and generate histograms
% for iji=1:1:TTTTT2;
% 
% ijk=iji;
% 
% % Apply logicle transform (cytometry-style scaling)
% obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
% xyxyxyfy = obj.transform(xyxyxyfx2(:, iji));
% figure(ijk)
% hold on
% box on
% % Plot histogram of logicle-transformed values
% histogram(xyxyxyfy, 60, 'Normalization', 'probability','FaceColor',[0.93,0.24,0.13]);
% axis ([0.0592 0.5282 0 0.4]);
% 
% % Set figure size
% fig=gcf;
% set(fig,'Units', 'points', 'Position', [0, 0, width, height]);
% 
% 
% % Format axes
% ax = gca;
% ax.XTick = obj.Tick;
% ax.XTickLabel = [];
% 
% set(ax, 'FontName', fontName, 'FontSize', fontSize);
% 
% ay = gca;
% ax.YTickLabel = [];
% 
%  end
% 
% 
% 
% 
