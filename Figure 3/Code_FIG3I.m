%% CODE TO GENERATE FIGURE 3I
% Correlation plot between gene expression and Fraction of CpG methylated

%%%% STEP 1: Load and visualize experimental data with error bars

% Mean values of fraction of methylated CpG
YFRACTIONLMHmean = [0.0501    0.3521    0.8029];

% Standard deviation for the CpG methylation
YFRACTIONLMHsd =[0.0109    0.0841    0.0768];  

% Mean values of gene expression
XexpressionLMHmean = [0.5321 0.3585 0.1487];

% Standard deviation for gene expression
XexpressionLMHsd = [0.0321 0.0769 0.0582];


x = XexpressionLMHmean;         % X-coordinates
y = YFRACTIONLMHmean;   % Y-coordinates
std_x = XexpressionLMHsd;  % Standard deviations for X
std_y = YFRACTIONLMHsd;  % Standard deviations for Y


% Plot with horizontal and vertical error bars
figure(1);
hold on
errorbar(y,x, std_x, std_x, std_y, std_y, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

ylabel('Gene expression  (A.U.F.)');
xlabel('Fraction of CpG methylated');
box on;
axis([0 1 0 0.6]);


%%%% STEP 2: Fit a linear regression and overlay the fitted line



% Reverse x and y for regression to keep 'methylation → expression'

x = YFRACTIONLMHmean;       % X-coordinates
y = XexpressionLMHmean;        % Y-coordinates

% Linear regression using polyfit
coefficients = polyfit(x, y, 1);
fitted_line = polyval(coefficients, x);

% Overlay fitted line
figure(1)
hold on;
plot(x, fitted_line, '-r', 'DisplayName', 'Fitted Line');
legend('show');

% Displaying the equation of the fitted line
equation = sprintf('y = %.4fx + %.4f', coefficients(1), coefficients(2));
disp('Fitted Line Equation:');
disp(equation);

% Compute and display p-value
mdl = fitlm(x, y);
disp('P-value:');
disp(mdl.Coefficients.pValue(2));  % p-value for the slope



%%%% STEP 3: Format figure for publication and apply logicle scaling

% Set figure dimensions in points
width = 100*5;
height = 75*5;

% Set the line thickness for the plot 
lineThickness = 2;

% Set the font size for axis labels, titles, etc. 
fontName = 'Arial';
fontSize = 14; 


% Create a figure with specific width and height
fig=gcf;
set(fig,'Units', 'points', 'Position', [0, 0, width, height]);

% Set Y-axis range (for logicle transform)
axis([0 0.9 0.0621 0.5482]);


% Format axes font
ax = gca;
set(ax, 'FontName', fontName, 'FontSize', fontSize);

% % Apply logicle transform (cytometry-style scaling)
obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
ax = gca;
ax.YTick = obj.Tick;
ax.YTickLabel = obj.TickLabel;


