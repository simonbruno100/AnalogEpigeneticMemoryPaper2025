%% CODE TO GENERATE FIGURE 4I
% Correlation plot between gene expression and Fraction of CpG methylated

%%%% STEP 1: Load and visualize experimental data with error bars

% Mean values of fraction of methylated CpG
YFRACTIONmonoclonalmean = [0.9214    0.8560    0.5760    0.4101    0.1588    0.1259    0.0096    0.0013];

% Standard deviation for the CpG methylation
YFRACTIONmonoclonalsd =[0.0337    0.0487    0.0526    0.0546    0.0530    0.0396    0.0196    0.0073];  

% Mean values of gene expression
Xexpressionmonoclonalmean = [0.0887, 0.1264, 0.2214, 0.2820, 0.2979, 0.4002, 0.4557, 0.4649];

% Standard deviation for gene expression
Xexpressionmonoclonalsd = [0.0389, 0.0397, 0.0681, 0.0630, 0.0699, 0.0481, 0.0312, 0.0325];


x = Xexpressionmonoclonalmean;
y = YFRACTIONmonoclonalmean;
std_x = Xexpressionmonoclonalsd;
std_y = YFRACTIONmonoclonalsd;


% Plot with horizontal and vertical error bars
figure(1);
hold on
errorbar(y,x, std_x, std_x, std_y, std_y, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

ylabel('Gene expression (A.U.F.)');
xlabel('Fraction of CpG methylated');
box on;
axis([0 1 0.05 0.55]);

%%%% STEP 2: Fit a linear regression and overlay the fitted line



% Reverse x and y for regression to keep 'methylation → expression'

y = Xexpressionmonoclonalmean;         % X-coordinates
x = YFRACTIONmonoclonalmean;   % Y-coordinates

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
width = 260*5;
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
axis([0 1 0.0641 0.5282]);


% Format axes font
ax = gca;
set(ax, 'FontName', fontName, 'FontSize', fontSize);

% % Apply logicle transform (cytometry-style scaling)
obj = logicleTransform(262144.0,0.2781913054120009,4.5,0);
ax = gca;
ax.YTick = obj.Tick;
ax.YTickLabel = obj.TickLabel;


