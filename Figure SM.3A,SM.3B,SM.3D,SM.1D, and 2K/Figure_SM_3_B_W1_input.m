%% Figure_SM_3_B
% This script computes the dose-response curve for the (bar W1, bar Y1) pair for different values of ε and ζ

%% Input range settings

N=4; % Maximum value of the input (\bar W1)

n=0.01; %increment step for the input

T=N/n; % Total number of input values considered


%% Compute steady-state values for increasing input

q = 0;                              % Initialize index
Xa    = zeros(1, T+1);              % Steady-state value of Da
YY1   = zeros(1, T+1);              % Steady-state value of Dr1 + Dr12
Xr2   = zeros(1, T+1);              % Steady-state value of Dr2
U1    = zeros(1, T+1);              % Steady-state value of W1

for j=0:n:N;
  
    [yy1,xr2,xa,u1] = W1_input_function(j); % Compute steady state at input j

    
 
    q=q+1;
    Xa(q)=xa;         %bar Da
    YY1(q)=yy1;       %bar Dr1 + bar Dr12
    Xr2(q)=xr2;       %bar Dr2
    U1(q)=u1;         %bar W1
    
end



%% Plotting parameters

width = 117;      % Width of figure (in points)
height = 39;      % Height of figure (in points)


% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 


%% Generate plot of \bar Y1 vs input

j=0:n:N;

figure(3)
hold on
plot(j, YY1,'LineWidth',1, 'Color', [0 0 0]);
set(gca, 'FontName', 'Times New Roman')
%ylabel('bar Y_1')
axis([0 N 0 1]);    % Set axis limits
box on


% Set figure size
fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);

% Format axes
ax = gca;
ax.XTickLabel = [];
set(ax, 'FontName', fontName, 'FontSize', fontSize);

ay = gca;
ax.YTickLabel = [];
