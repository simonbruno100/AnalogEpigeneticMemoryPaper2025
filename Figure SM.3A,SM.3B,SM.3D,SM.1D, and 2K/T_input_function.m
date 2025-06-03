function [yy1,xr12,xr2,xa,mu1] = T_input_function(uu)

% INPUT:
%   uu   - value \bar T (TET) at t=0
%
% OUTPUTS (steady-state values at final simulation time):
%   yy1   - bar(Y1)
%   xr12  - bar(Dr12)
%   xr2   - bar(Dr2)
%   xa    - bar(Da)
%   mu1   - bar(T)
%

%% Parameter values

k = 15;          
Dt = 15;                  % Total number of nucleosomes
DNAlevel = k / Dt;        % Normalized initial amount of DR12( and thus Y1)


a = 1;     %alpha
ab = 1;    %bar alpha
ap = 0;    %alpha prime


b = 1;
mu =0.1;      

bb=1;         %beta

e = 0.1;      %epsilon
ee = 5;       %zeta

uua = 15;  %uA


N=1400;  %Final simuation time

%% Build the simulation

%% Define ODE system
% y1 = bar(Y1)
% y2 = \bar Dr12
% y3 = \bar Dr2
% y4 = \bar Da (active)
% y5 = \bar T (external input)


sys_dyn = @(t,y) [-y(5)*(bb*e + ee*y(4))*y(1);...
                 (a*y(3) + ab*y(1))*(y(1) - y(2)) - (y(5)*(bb*e + ee*y(4))+mu*(b*e + ee*y(4)))*y(2);...
                 (a*y(3) + ab*y(1))*(1 - y(1) - y(3) - y(4)) + y(5)*(bb*e + ee*y(4))*y(2)- (mu*(b*e + ee*y(4)))*y(3);...
                 (uua + y(4))*(1 - y(1) - y(3)- y(4))- (e + ee*(y(1) + y(3)))*y(4);...
                 uu*(-e)*exp(-e*t)];


%% Simulate dynamics using ODE45

 [T, Y] = ode45(sys_dyn, [0 N], [DNAlevel DNAlevel 0 1-DNAlevel uu]');



%% Extract steady-state values

yy1  = Y(end,1);   % bar(Y1) = Dr1 + Dr12
xr12 = Y(end,2);   % bar(Dr12)
xr2  = Y(end,3);   % bar(Dr2)
xa   = Y(end,4);   % bar(Da)
mu1  = Y(end,5);   % bar(T)



end

