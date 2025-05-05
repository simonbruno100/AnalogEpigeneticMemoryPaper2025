%% Figure SM_3A
% The script generates input/output steady-state characteristics for the pair 
% (Y1̄, Dā) across different values of epsilon (ε) and zeta (ζ)

%% Input range settings

N=1; % Maximum value of the input (that is, Dt/Dt = 1)

n=0.005; % increment step for the input

T=N/n; % Total number of input values considered


%% Simulate increasing input


q = 0;                              % Initialize index
Xa = zeros(1, T+1);                 % Steady-state value of bar Da
Xr2 = zeros(1, T+1);                % Steady-state value of bar Dr2

for j=0:n:N;
  
    [xr2,xa] = Y1_input_function_1(j); % Compute steady state at input j
    

    q=q+1;
    Xa(q)=xa;         %bar Da
    Xr2(q)=xr2;       %bar Dr2
    
end



%% Plotting parameters

width = 117;
height = 39;


% Set the line thickness for the plot
lineThickness = 0.5;

% Set the font size for axis labels, titles, etc.
fontName = 'Arial';
fontSize = 7; 



%% Generate plot of (Y1̄, Dā) (increasing input)

j=0:n:N;
figure(100)
hold on
plot(j, Xa,'LineWidth',1, 'Color', [0 0 0]);
set(gca, 'FontName', 'Times New Roman')
%ylabel('bar D^A')
axis([0 N 0 1]);
box on


%% Simulate decreasing input


xa0=Xa(end);
xr20=Xr2(end);

% Reset index and output arrays
q=0;
Xa=zeros(1,T+1);
Xr2=zeros(1,T+1);

for j=0:n:N;

    [xr2,xa] = Y1_input_function_2(j,xr20,xa0);

    q=q+1;
    Xa(q)=xa;         %bar Da
    Xr2(q)=xr2;       %bar Dr2 + bar Dr12


end



%% Generate plot of (Y1̄, Dā) (decreasing input)

j=0:n:N;
figure(100)
hold on
plot(j, Xa,'LineWidth',1, 'Color', [0 0 0]);
set(gca, 'FontName', 'Times New Roman')
%ylabel('bar D_A')
axis([0 N 0 1]);
box on

% Final formatting

fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);

ax = gca;
ax.XTickLabel = [];
set(ax, 'FontName', fontName, 'FontSize', fontSize);

ay = gca;
ax.YTickLabel = [];
