%% CODE TO GENERATE FIGURE SM.2
% This code generates distributions of the variable n_x from stochastic simulations 
% under specified initial conditions, and visualizes the results as histograms (logicle-scaled).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for distributions

a=0; % Number of final timepoints to retain per simulation (0 means only 1 point kept per sim)

TTT=1000; % Total number N of simulations per distribution

NumberOfDays=4; % Time interval (in days) between two distributions

TTTTT=2; % Total number of distributions to generate


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set general biological parameters and initial conditions


Dt=15; % Total number of nucleosomes in the gene of interest

p.V=1; % Reaction volume (used to normalize concentrations)

tfinn=24*NumberOfDays*TTTTT; % Total simulation time (in hours)

% Initialize matrices to collect final simulation values (for each distribution)

ddaaf=zeros(1,TTTTT);% DA
dd11f=zeros(1,TTTTT);% DR1
dd22f=zeros(1,TTTTT);% DR2
dd1122f=zeros(1,TTTTT);% DR12
mmxxf=zeros(1,TTTTT);% mRNA
xf=zeros(1,TTTTT);% Protein

 for i=1:1:TTT;

% Select initial conditions for the simulation:
% Blue curve in Figure SM.2:   (DR120, DR10, DR20, DA0) = (14,0,1,0)
% Red curve in Figure SM.2:    (DR120, DR10, DR20, DA0) = (6,0,8,1)
% Yellow curve in Figure SM.2: (DR120, DR10, DR20, DA0) = (4,0,5,6)
% Purple curve in Figure SM.2: (DR120, DR10, DR20, DA0) = (1,0,1,13)

DR120=6;
DR10=0;
DR20=8;
DA0=1;



% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_2(DA0,DR10,DR20,DR120,tfinn);

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

x2=[xf2;...
    zeros(1,TTTTT);...
    (10^(4))*ones(1,TTTTT)];



%% Generate plots
% Generate logicle-scaled histograms of n_x for each timepoint

% Figure size

width = 103;
height = 52;

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
axis ([0.0592 0.5282 0 0.5]);

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
