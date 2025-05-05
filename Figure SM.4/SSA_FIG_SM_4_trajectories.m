%% CODE TO GENERATE FIGURE SM.4 - top plots
% This code generates trajectoros of the variable n_D^A 
% (amount of nucleosomes with activating modifications) from stochastic simulations for different values of Y1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   System reactions

% PLEASE NOTE: D = n_{D}, Da = n_{DA}; D2 = n_{DR2};

%  1. D                  --    kw20           --> D2
%  2. D                  --    kw2            --> D2
%  3. D                  --    km             --> D2
%  4. D                  --    kmbar          --> D2
%  5. D2                 --    delta          --> D
%  6. D2                 --    ke             --> D
%  7. D2                 --    keact          --> D
%  8. D                  --    kwa0           --> Da
%  9. D                  --    kwa            --> Da
%  10. D                  --    kma            --> Da
%  11. Da                 --    delta          --> D
%  12. Da                 --    kea            --> D
%  13. Da                 --    keacta         --> D
%  14. Da                 --    keacta         --> D


%% Parameter values and Initial conditions


Dt = 15;                       % Total number of nucleosomes


% Set the amount of Y1 (range: 0 to Dt). Values used in the paper are: 0,1,6,7,14
rrrr = 1;                      % Amount of Y1

D2_initial = 0; 
Da_initial = Dt - rrrr;        % Set initial Da based on Y1
tttt = 24 * 28;                % Total simulation time (in hours)



% No external erasers or DAPG
Ea=0;
Er1=0;
Er2=0;
p.DAPG=0;

p.V = 1;                       % Volume (used for normalization)

%%%

p.kkw20=0;          %kW20
p.KRAB0=0;          %No KRAB (-->kW2 = 0)

p.Y1=rrrr;          %Amount of Y1 used in the simulation

p.kkm=0.0347;       %km
p.kkmbar=0.0347;    %\bar km

p.delta=0.035;     %delta
p.ke=0.0032;        %\bar kER 
p.kkeact=0.0868;    %kER


%%%

p.kkwa0=0;            %kWA0
p.kkwa=7.8075;        %kWA

p.kkma=0.0347;        %kmA

p.kea=0.0315;         %\bar kEA

p.kkeacta=0.8675;     %kEA


%%%

p.Ea=Ea;

p.Er1=Er1;

p.Er2=Er2;



%% Initial conditions

D_initial = Dt - rrrr - D2_initial - Da_initial;


IC = [D_initial, D2_initial, Da_initial];

%% Initial state
tinterval = [0, tttt]; %hours --> this value multipled by Dt*p.kma will give the corresponding normalized time tau

x0    = IC;
x0    = round(x0);

%% Specify reaction network
pfun = @propensity_functions;

%species: D,  Dr2,  Da, 3

stoich_matrix = [-1  1  0
                 -1  1  0
                 -1  1  0 
                 -1  1  0 
                  1 -1  0 
                  1 -1  0  
                  1 -1  0 
                 -1  0  1   
                 -1  0  1  
                 -1  0  1  
                  1  0 -1  
                  1  0 -1 
                  1  0 -1  
                  1  0 -1]; 
   

%% Run simulation:Direct method
[t,x] = directMethod2(stoich_matrix, pfun, tinterval, x0, p);

%% Plot Da trajectory (smoothed, if needed) over time

Xfin = [x];     % Store state trajectory
Tfin = [t];     % Store time points

mmeann = 100;   % Window size for moving average smoothing 
                % (set to 1 for no smoothing; increase for smoother curve)


A = movingmean(Xfin(:,3), mmeann, 1);   % Compute Da


% Figure size
width =  76;
height = 52;

% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 


figure(6)
box on
hold on;
plot(Tfin/24,A,'LineWidth',2, 'Color', [0 0 0]);
axis([0  tttt/24 0 Dt]);

% Set figure size
fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);

% Format axes
ax = gca;
ax.XTickLabel = []; % Hide X labels

set(ax, 'FontName', fontName, 'FontSize', fontSize);

ay = gca;
ax.YTickLabel = []; % Hide Y labels




%% 
function a = propensity_functions(x, p, dd,dd2, DDAPG,EEa,EEr1,EEr2);

% Return reaction propensities given current state x
%species: D,  Dr2, Da, 3

D          = x(1);
D2         = x(2);
Da         = x(3);



p.km2=(p.kkm*(D2))/p.V;

p.kmbar1=(p.kkmbar*(p.Y1))/p.V;

p.kma=(p.kkma*(Da))/p.V;



p.keBIS=p.ke*(1+EEr2*dd);

p.kkeactBIS=p.kkeact*(1+EEr2*dd);




p.keaBIS=p.kea*(1+EEa*dd);

p.kkeactaBIS=p.kkeacta*(1+EEa*dd);




p.keactBIS=(p.kkeactBIS*Da)/p.V;


p.keacta11BIS=(p.kkeactaBIS*(p.Y1))/p.V;

p.keacta22BIS=(p.kkeactaBIS*(D2))/p.V;




p.kw20=p.kkw20;

p.kkw2=(p.KRAB0)*dd;


p.kw2=p.kkw2;


p.kwa0=p.kkwa0;
p.kwa=p.kkwa;


a = [p.kw20*D;
     p.kw2*D;
     p.km2*D;
     p.kmbar1*D;
     p.delta*D2;
     p.keBIS*D2;
     p.keactBIS*D2;
     p.kwa0*D;
     p.kwa*D;
     p.kma*D;
     p.delta*Da;
     p.kea*Da;
     p.keacta11BIS*Da;
     p.keacta22BIS*Da];
 
end

function [ t, x ] = directMethod2( stoich_matrix, pfun, tinterval, x0,...
                                  p)
    %% Initialization

max_out_length=1000000;
num_rxns = size(stoich_matrix, 1);
num_species = size(stoich_matrix, 2);
T = zeros(max_out_length, 1);
X = zeros(max_out_length, num_species);
DEC = zeros(max_out_length, 1);
DEC2 = zeros(max_out_length, 1);
T(1)     = tinterval(1);
X(1,:)   = x0;
rxn_count = 1;
DEC(1)=1;
DEC2(1)=1;
DDDAPG(1)=0;
EEEa(1)=0;
EEEr1(1)=0;
EEEr2(1)=0;

%% MAIN LOOP
while T(rxn_count) < tinterval(2);
    
    % Calculate reaction propensities
    a = pfun(X(rxn_count,:), p, DEC(rxn_count),DEC2(rxn_count), DDDAPG(rxn_count), EEEa(rxn_count), EEEr1(rxn_count), EEEr2(rxn_count));
    
    % Sample earliest time-to-fire (tau)
    a0 = sum(a);
    r = rand(1,1);
    
    % Sample identity of earliest reaction channel to fire (mu)
    j=1; s=a(1); r0=r*a0;
    while s < r0
      j = j + 1;
      s = s + a(j);
    end
    
    % Sample earliest time-to-fire (tau)
    
    s2= s -a(j);

    r2 = (r0 - s2)/a(j);
    
    tau = -log(r2)/a0; %(1/a0)*log(1/r2);
    

ddayy=9;
    
    
    if T(rxn_count)<ddayy*24;
%     if T(rxn_count)<0;
        
        DDDAPG(rxn_count+1)=0;

        EEEr1(rxn_count+1)=p.Er1;
%         EEEr1(rxn_count+1)=0;

        EEEr2(rxn_count+1)=0;
        
        EEEa(rxn_count+1)=0;
    else
         DDDAPG(rxn_count+1)=p.DAPG;
         
         EEEr1(rxn_count+1)=p.Er1;
         
         EEEr2(rxn_count+1)=p.Er2;
         
         EEEa(rxn_count+1)=p.Ea;
    end
        
    
    
    
    
    
    % Update time and carry out reaction mu
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(j,:); 
    DEC(rxn_count+1)=exp(-p.delta*T(rxn_count+1));   
    DEC2(rxn_count+1)=exp(-p.delta*(T(rxn_count+1)-ddayy*24));   
    rxn_count = rxn_count + 1;



end  

% Return simulation time course
t = T(1:rxn_count);
x = X(1:rxn_count,:);
dec = DEC(1:rxn_count);
dec2 = DEC2(1:rxn_count);
if t(end) > tinterval(2)
    t(end) = tinterval(2);
    x(end,:) = X(rxn_count-1,:);
end    

end
