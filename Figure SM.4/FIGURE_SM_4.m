
%% CODE TO GENERATE FIGURE SM.4 - bottom plots
% This code generates input-dependent distributions of the variable n_x
% (protein output) from stochastic simulations for different values of Y1
% Simulations are conducted using the SSA_FIG_SM_4 function.
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

% Generate array of Y1 values (length TTT) to be used across simulations.
% By setting l1 = l2 = X, all entries in 'rrrrtotal' will be equal to X.
% This allows the user to select the desired constant value of Y1
% (e.g., X = 0, 1, 6, 8, or 14) for the entire simulation batch.

l1=1;
l2=1;

rrrrtotal = unifrnd(l1,l2,1,TTT);



% Initialize matrices to collect final simulation values (for each distribution)

ddaaf=zeros(1,TTTTT); %DA
dd22f=zeros(1,TTTTT); %DR2
mmxxf=zeros(1,TTTTT); %DA
xf=zeros(1,TTTTT);% Protein


 for i=1:1:TTT;

% Initial condition
DR10=rrrrtotal(i);
DA0=15-rrrrtotal(i);
DR20=0;


% Run stochastic simulation
[X0o0,T0o0] = SSA_FIG_SM_4(DA0,DR20,tfinn,DR10);


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
ddaaf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,3)/p.V;
dd22f(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,2)/p.V;
mmxxf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,4)/p.V;
xf(a*(i-1)+i:a*i+i,iji)=X0o0(index-a:index,5)/p.V;

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

width =  76;
height = 52;

% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 


iji=TTTTT; % Set the final time point



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


