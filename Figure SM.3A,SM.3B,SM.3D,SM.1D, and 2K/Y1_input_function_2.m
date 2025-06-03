function [xr2,xa] = Y1_input_function_2(uu,xr20,xa0)

% Used for the decreasing-input simulation in Figure SM_3A.

% INPUT:
%   uu  - Input level (i.e., kmet, that represents bar Y1)
%
% OUTPUTS:
%   xr2 - Steady-state value of bar(Dr2) 
%   xa  - Steady-state value of bar(Da)


%% Parameters

Dt = 15;         % Total number of chromatin sites (not explicitly used in dynamics)
kmet = uu;       % Input


a = 1;     %alpha
ab = 1;    %bar alpha
ap = 0;    %alpha prime


b = 1;
mu =0.1;      

e = 0.0001;     %epsilon
ee = 1;         %zeta



uua = 15;       %uA


N=7400;  %Final simulation time

%% Build the simulation

% Define ODE system
% Variables:
%   y1 = bar(Dr2) 
%   y2 = bar(Da)



sys_dyn = @(t,y) [(a*y(1) + ab*kmet)*(1 - y(1) - y(2) - kmet) - (mu*(b*e + ee*y(2)))*y(1);...
                 (uua + y(2))*(1-y(1)-y(2)- kmet)- (e + ee*(y(1) + kmet))*y(2)];


% Simulate dynamics using ODE45


 [T, Y] = ode45(sys_dyn, [0 N], [1 - kmet 0]');




%% Extract steady-state values
                     
                        
  
  
xr2=Y(end,1);  

xa=Y(end,2);


end
