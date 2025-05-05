function [Xfin,Tfin] = SSA_FIG_SM_5(Da_initial,D1_initial,D2_initial,D12_initial,tttt,DAPG,rrrr);
%
% INPUTS:
%   Da_initial  - Initial number of Da nucleosomes 
%   D1_initial  - Initial number of D1 nucleosomes 
%   D2_initial  - Initial number of D2 nucleosomes 
%   D12_initial - Initial number of D12 nucleosomes
%   tttt        - Final simulation time (in hours)
%   DAPG        - DAPG level
%   rrrr        - initial external input level (KRAB,DNMT3)
%
% OUTPUTS:
%   Xfin        - State trajectories over time (matrix with columns: D, D1, D2, D12, Da, mRNA, protein)
%   Tfin        - Time points corresponding to state transitions


% NOTES:
%   - Total number of nucleosomes: Dt = 15
%   - Parameters like `p.kkma`, `p.kea`, etc., can be tuned based on the values reported in the paper
%   - Propensities are defined in `propensity_functions`
%   - Gillespie simulation is implemented via `directMethod2`

%%   System reactions

% PLEASE NOTE: D = n_{D}, Da = n_{DA}; D1 = n_{DR1}; D2 = n_{DR2}; D12 = n_{DR12}

%   1. D                  --    kw10           --> D1
%   2. D                  --    kw1            --> D1%   
%  3. D                  --    kmprime        --> D1
%   4. D                  --    kmprime        --> D1
%   2. D                     --    kDMT1            --> D1
%   2. D +D1                 --    k'DNMT1            --> D1 + D1
%   2. D +D12                 --    k'DNMT1            --> D1 + D12
%   3. D1                 --    deltaprime     --> D
%   3. D1                 --    deltaprime'(1/(1+(D1+D12)))     --> D
%   4. D1                 --    ktprime        --> D 
%   5. D1                 --    ktprimeact     --> D
%   6. D1                 --    kw20           --> D12
%   7. D1                 --    km             --> D12
%   8. D1                 --    kmbar          --> D12
%   9. D1                 --    kmbar          --> D12
%  10. D12                --    delta          --> D1
%  11. D12                --    ke             --> D1
%  12. D12                --    keact          --> D1
%  13. D                  --    kw20           --> D2
%  14. D                  --    kw2            --> D2
%  15. D                  --    km             --> D2
%  16. D                  --    kmbar          --> D2
%  17. D                  --    kmbar          --> D2
%  18. D2                 --    delta          --> D
%  19. D2                 --    ke             --> D
%  20. D2                 --    keact          --> D
%  21. D2                 --    kw10           --> D12
%  26. D2                 --    kmprime        --> D12
%  27. D2                 --    kmprime        --> D12
%  21. D2                 --    kDNMT1           --> D12
%  21. D2 +D1                --    k'DNMT1           --> D12 +D1
%  21. D2 +D12                --    k'DNMT1           --> D12 + D12
%  22. D12                --    deltaprime     --> D2
%  22. D12                --    deltaprime'(1/(1+(D1+D12)))     --> D2
%  23. D12                --    ktprime        --> D2
%  24. D12                --    ktprimeact     --> D2
%  25. D                  --    kwa0           --> Da
%  26. D                  --    kwa            --> Da
%  27. D                  --    kma            --> Da
%  28. Da                 --    delta          --> D
%  29. Da                 --    kea            --> D
%  30. Da                 --    keacta         --> D
%  31. Da                 --    keacta         --> D
%  32. Da                 --    keacta         --> D
%  33. Da                 --    max            --> Da+mx
%  34. D                  --    mx             --> D+mx
%  35. DR1                --    m1x            --> DR1+mx
%  36. DR2                --    m2x            --> DR2+mx
%  37. DR12               --    m12x           --> DR12+mx
%  38. mx                 --    px             --> mx+P
%  39. mx                 --    em             --> 0
%  40. X                  --    ex             --> 0

%% Parameter values and Initial conditions


Dt = 15; % Dtot



% External erasers
Ea=0;
Er1=0;
Er2=0;



p.V = 1; % Volume


%%%

p.kkw20=0;          %kW20

% Pick one of the two options below depending on whether the external input is KRAB:
% p.KRAB0=0;               %No KRAB (-->kW2 = 0)
p.KRAB0=rrrr*0.5205;    %With KRAB input   (-->kW2 \neq 0)


p.kkm=0.0347;       %km
p.kkmbar=0.0347;    %\bar km
p.kkmprime=0;       %km '

p.delta=0.035;     %delta
p.ke=0.0018;        %\bar kER 
p.kkeact=0.0868;    %kER


%%%


p.kkw10=0;              %kW10

% Pick one of the two options below depending on whether the external input is DNMT3:
p.DNMT30=0;             %No DNMT3 (-->kW1 = 0)
% p.DNMT30=rrrr*0.5205;  %With DNMT3 input   (-->kW1 \neq 0)


p.DAPG=DAPG;




p.deltaprime=0;       % delta'

p.ktprime=0;   %kT'

p.kktprime=0;  %kT'*


p.kkmDNMT1=0;         % To be set = 0, since this is a parameter associated with the refined DNA methylation maintenance model (see Fig. SM7).
p.kkmDNMT1act=0;      % To be set = 0, since this is a parameter associated with the refined DNA methylation maintenance model (see Fig. SM7).
p.deltaprimeact= 0;   % To be set = 0, since this is a parameter associated with the refined DNA methylation maintenance model (see Fig. SM7). 


%%%

p.kkwa0=0;            %kWA0
p.kkwa=7.8075;        %kWA

p.kkma=0.0347;        %kmA


p.kea=0.0176;         %\bar kEA

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

D_initial = Dt - D1_initial - D2_initial - D12_initial - Da_initial;


IC = [D_initial, D1_initial, D2_initial, D12_initial, Da_initial, m_initial, P_initial];

%% Initial state
tinterval = [0, tttt]; %hours --> this value multipled by Dt*p.kma will give the corresponding normalized time tau

x0    = IC;
x0    = round(x0);

%% Specify reaction network
pfun = @propensity_functions;

%species: D, Dr1, Dr2, Dr12, Da, m, P 7

stoich_matrix = [-1  1  0  0  0  0  0
                 -1  1  0  0  0  0  0
                 -1  1  0  0  0  0  0
                 -1  1  0  0  0  0  0
                 -1  1  0  0  0  0  0
                 -1  1  0  0  0  0  0
                 -1  1  0  0  0  0  0
                  1 -1  0  0  0  0  0
                  1 -1  0  0  0  0  0
                  1 -1  0  0  0  0  0
                  1 -1  0  0  0  0  0
                  0 -1  0  1  0  0  0
                  0 -1  0  1  0  0  0
                  0 -1  0  1  0  0  0           
                  0 -1  0  1  0  0  0
                  0 -1  0  1  0  0  0               
                  0  1  0 -1  0  0  0
                  0  1  0 -1  0  0  0
                  0  1  0 -1  0  0  0
                 -1  0  1  0  0  0  0
                 -1  0  1  0  0  0  0
                 -1  0  1  0  0  0  0
                 -1  0  1  0  0  0  0           
                 -1  0  1  0  0  0  0
                 -1  0  1  0  0  0  0           
                  1  0 -1  0  0  0  0
                  1  0 -1  0  0  0  0
                  1  0 -1  0  0  0  0
                  0  0 -1  1  0  0  0 
                  0  0 -1  1  0  0  0 
                  0  0 -1  1  0  0  0 
                  0  0 -1  1  0  0  0 
                  0  0 -1  1  0  0  0 
                  0  0 -1  1  0  0  0           
                  0  0  1 -1  0  0  0          
                  0  0  1 -1  0  0  0
                  0  0  1 -1  0  0  0
                  0  0  1 -1  0  0  0
                 -1  0  0  0  1  0  0
                 -1  0  0  0  1  0  0
                 -1  0  0  0  1  0  0
                  1  0  0  0 -1  0  0
                  1  0  0  0 -1  0  0
                  1  0  0  0 -1  0  0
                  1  0  0  0 -1  0  0
                  1  0  0  0 -1  0  0
                  0  0  0  0  0  1  0
                  0  0  0  0  0  1  0
                  0  0  0  0  0  1  0
                  0  0  0  0  0  1  0
                  0  0  0  0  0  1  0
                  0  0  0  0  0  0  1
                  0  0  0  0  0 -1  0
                  0  0  0  0  0  0 -1]; 
   

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
%species: D, Dr1, Dr2, Dr12, Da, m, X 7

D          = x(1);
D1         = x(2);
D2         = x(3);
D12        = x(4);
Da         = x(5);
m          = x(6);
P          = x(7);


p.kmprime2=(p.kkmprime*(D2))/p.V;
p.kmprime12=(p.kkmprime*(D12))/p.V;

p.km2=(p.kkm*(D2))/p.V;
p.km12=(p.kkm*(D12))/p.V;

p.kmbar1=(p.kkmbar*(D1))/p.V;
p.kmbar12=(p.kkmbar*(D12))/p.V;

p.kma=(p.kkma*(Da))/p.V;

% 


p.kmDNMT1act1=(p.kkmDNMT1act*(D2))/p.V;

p.kmDNMT1act12=(p.kkmDNMT1act*(D12))/p.V;


p.deltaprimeact1=(p.deltaprimeact*(1/(1+(D2+D12))))/p.V;


p.ktprimeBIS=p.ktprime*(1+EEr1)*dd;

p.kktprimeBIS=p.kktprime*(1+EEr1)*dd;




p.keBIS=p.ke*(1+EEr2*dd);

p.kkeactBIS=p.kkeact*(1+EEr2*dd);




p.keaBIS=p.kea*(1+EEa*dd);

p.kkeactaBIS=p.kkeacta*(1+EEa*dd);



p.ktprimeactBIS=(p.kktprimeBIS*Da)/p.V;

p.keactBIS=(p.kkeactBIS*Da)/p.V;


p.keacta11BIS=(p.kkeactaBIS*(D1))/p.V;
p.keacta112BIS=(p.kkeactaBIS*(D12))/p.V;

p.keacta22BIS=(p.kkeactaBIS*(D2))/p.V;


p.kw10=p.kkw10;

p.kkw1=(p.DNMT30)*dd;



p.kw1=p.kkw1;



p.kw20=p.kkw20;

p.kkw2=(p.KRAB0/(1+DDAPG))*dd;

p.kw2=p.kkw2;


p.kwa0=p.kkwa0;
p.kwa=p.kkwa;


a = [p.kw10*D;
     p.kw1*D;
     p.kmprime2*D;
     p.kmprime12*D;
     p.kkmDNMT1*D;
     p.kmDNMT1act1*D;
     p.kmDNMT1act12*D;
     p.deltaprime*D1;
     p.deltaprimeact1*D1;
     p.ktprimeBIS*D1;
     p.ktprimeactBIS*D1;
     p.kw20*D1;
     p.km2*D1;
     p.km12*D1;
     p.kmbar1*(D1-1);
     p.kmbar12*D1;
     p.delta*D12;
     p.keBIS*D12;
     p.keactBIS*D12;
     p.kw20*D;
     p.kw2*D;
     p.km2*D;
     p.km12*D;
     p.kmbar1*D;
     p.kmbar12*D;
     p.delta*D2;
     p.keBIS*D2;
     p.keactBIS*D2;
     p.kw10*D2;
     p.kmprime2*(D2-1);
     p.kmprime12*D2;
     p.kkmDNMT1*D2;
     p.kmDNMT1act1*D2;
     p.kmDNMT1act12*D2;
     p.deltaprime*D12;
     p.deltaprimeact1*D12;
     p.ktprimeBIS*D12;
     p.ktprimeactBIS*D12;
     p.kwa0*D;
     p.kwa*D;
     p.kma*D;
     p.delta*Da;
     p.kea*Da;
     p.keacta11BIS*Da;
     p.keacta112BIS*Da;
     p.keacta22BIS*Da;
     p.pax*Da;
     p.px*D;
     p.px1*D1;
     p.px2*D2;
     p.px12*D12;
     p.pm*m;
     p.em*m;
     p.ex*P;];
 
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
    

ddayy=6;
    
    
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
%     DEC2(rxn_count+1)=exp(-p.delta*(T(rxn_count+1)));  
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