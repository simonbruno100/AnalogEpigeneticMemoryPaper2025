

function [Xfin,Tfin] = SSA_FIG_SM_4(Da_initial,D2_initial, tttt,rrrr);
%
% INPUTS:
%   Da_initial  - Initial number of Da nucleosomes 
%   D2_initial  - Initial number of D2 nucleosomes 
%   tttt        - Final simulation time (in hours)
%   rrrr        - external input level (Y1)
%
% OUTPUTS:
%   Xfin        - State trajectories over time (matrix with columns: D, D1, D2, Da, mRNA, protein)
%   Tfin        - Time points corresponding to state transitions


% NOTES:
%   - Total number of nucleosomes: Dt = 15
%   - Parameters like `p.kkma`, `p.kea`, etc., can be tuned based on the values reported in the paper
%   - Propensities are defined in `propensity_functions`
%   - Gillespie simulation is implemented via `directMethod2`

%%   System reactions

% PLEASE NOTE: D = n_{D}, Da = n_{DA}; D2 = n_{DR2}; D12 = n_{DR12}

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
%  15. Da                 --    max            --> Da+mx
%  16. D                  --    mx             --> D+mx
%  17. Y1                 --    m1x            --> DR1+mx
%  18. DR2                --    m2x            --> DR2+mx
%  19. mx                 --    px             --> mx+P
%  20. mx                 --    em             --> 0
%  21. X                  --    ex             --> 0

%% Initial conditions


Dt = 15; % Dtot

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

p.pax=0.2556;   %betaA_m
p.px=0.0021;   %beta_m
p.px1=0.0021;  %beta_m
p.px2=0.0021;  %beta_m
p.px12=0.0021; %beta_m

p.pm=2.52;     %beta

p.em=0.24;     %gamma_m
p.ex=0.035;   %delta



m_initial = (p.pax/p.em)*Da_initial; 
P_initial = (p.pm/p.ex)*m_initial; 

%%%

p.Ea=Ea;

p.Er1=Er1;

p.Er2=Er2;


%% Initial conditions

D_initial = Dt - rrrr - D2_initial - Da_initial;


IC = [D_initial, D2_initial, Da_initial, m_initial, P_initial];

%% Initial state
tinterval = [0, tttt]; %hours --> this value multipled by Dt*p.kma will give the corresponding normalized time tau

x0    = IC;
x0    = round(x0);

%% Specify reaction network
pfun = @propensity_functions;

%species: D,  Dr2,  Da, m, P 5

stoich_matrix = [-1  1  0  0  0
                 -1  1  0  0  0
                 -1  1  0  0  0 
                 -1  1  0  0  0 
                  1 -1  0  0  0 
                  1 -1  0  0  0  
                  1 -1  0  0  0  
                 -1  0  1  0  0   
                 -1  0  1  0  0  
                 -1  0  1  0  0  
                  1  0 -1  0  0  
                  1  0 -1  0  0  
                  1  0 -1  0  0  
                  1  0 -1  0  0  
                  0  0  0  1  0
                  0  0  0  1  0
                  0  0  0  1  0
                  0  0  0  1  0
                  0  0  0  0  1
                  0  0  0 -1  0
                  0  0  0  0 -1]; 
   

%% Run simulation:Direct method
[t,x] = directMethod2(stoich_matrix, pfun, tinterval, x0, p);

%% Save data

Xfin = [x];

Tfin =[t];

save('Xfin')

save('Tfin') 

%% 
function a = propensity_functions(x, p, dd,dd2, DDAPG,EEa,EEr1,EEr2);

% Return reaction propensities given current state x
%species: D, Dr1, Da, m, X 5

D          = x(1);
D2         = x(2);
Da         = x(3);
m          = x(4);
P          = x(5);



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
     p.keacta22BIS*Da;
     p.pax*Da;
     p.px*D;
     p.px1*p.Y1;
     p.px2*D2;
     p.pm*m;
     p.em*m;
     p.ex*P];
 
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
    
    tau = -log(r2)/a0; 

ddayy=9;
    
    
    if T(rxn_count)<ddayy*24;
        
        DDDAPG(rxn_count+1)=0;

        EEEr1(rxn_count+1)=p.Er1;

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
end
