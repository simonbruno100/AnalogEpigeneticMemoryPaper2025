function [yy1,xr2,xa,u1] = W1_input_function(uu)


% INPUT:
%   uu   - value \bar W1 (DNMT3) at t=0
%
% OUTPUTS (steady-state values at final simulation time):
%   yy1   - bar(Y1)
%   xr2   - bar(Dr2)
%   xa    - bar(Da)
%   u1   - bar(W1)
%

%% Parameter values

k = 0;          
Dt = 15;                  % Total number of nucleosomes
DNAlevel = k / Dt;        % Normalized initial amount of Y1



a = 1;     %alpha
ab = 1;    %bar alpha
ap = 0;    %alpha prime


b = 1;
mu =0.1;      

e = 0.1;      %epsilon
ee = 5;       %zeta

uua = 15;  %uA


N=1400;  %Final simuation time


%% Build the simulation

%% Define ODE system
% y1 = bar(Y1)
% y2 = \bar Dr2
% y3 = \bar Da (active)
% y4 = \bar W1 (external input)



sys_dyn = @(t,y) [y(4)*(1-y(1)-y(2)-y(3));...
                 (a*y(2) + ab*y(1))*(1 - y(1) - y(2) - y(3)) - (mu*(b*e + ee*y(3)))*y(2);...
                 (uua + y(3))*(1-y(1)-y(2)- y(3))- (e + ee*(y(1) + y(2)))*y(3);...
                 uu*(-e)*exp(-e*t)];



%% Simulate dynamics using ODE45

 [T, Y] = ode45(sys_dyn, [0 N], [DNAlevel 0 1-DNAlevel uu]');



%% Extract steady-state values

yy1 = Y(end,1);   % bar(Y1) = Dr1 + Dr12           
  
xr2 = Y(end,2);   % bar(Dr2) 

xa  = Y(end,3);   % bar(Da)       
  
u1  = Y(end,4);   % bar(W1)    


end

