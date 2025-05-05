%% CODE TO GENERATE FIGURE 4G
% This code generates the probability distribution of the fraction of CpGs methylated.
% It first loads bisulfite sequencing data from a CSV file, selects the region of interest (e.g., promoter), and computes
% the distribution of the number of CpG sites methylated across that
% region.

%% Load and prepare data

% Step 1: Save your Excel file as a CSV (e.g., BisData.csv)

% Step 2: Read the CSV file as a table 
 tableData = readtable('BisData.csv', 'VariableNamingRule', 'preserve');

% Step 3: Convert the table to a numeric matrix for easier access
matrixData = table2array(tableData);

%% Average methylation levels across replicates


TOT=206;  % Total number of regions (e.g., promoter tiles, genomic segments)

totalbin = 8;

% Step 4: Select the time point
% You can switch between first time point/second time point as needed

% % %%% First time point

bin1Tot=matrixData(:,23);
bin2Tot=matrixData(:,25);
bin3Tot=matrixData(:,27);
bin4Tot=matrixData(:,29);
bin5Tot=matrixData(:,31);
bin6Tot=matrixData(:,33);
bin7Tot=matrixData(:,35);
bin8Tot=matrixData(:,37);

% % %%% Second time point
% % 
% bin1Tot=matrixData(:,39);
% bin2Tot=matrixData(:,41);
% bin3Tot=matrixData(:,43);
% bin4Tot=matrixData(:,45);
% bin5Tot=matrixData(:,47);
% bin6Tot=matrixData(:,49);
% bin7Tot=matrixData(:,51);
% bin8Tot=matrixData(:,53);


%% Select region of interest


% Step 5: Select the region (e.g., promoter + rTetR/dCas9 binding sites)
% You can switch between bin1/bin2/.../bin8 as needed


 Vector = bin1Tot(43:70);
% 
% Vector = bin2Tot(43:70);
% 
% Vector = bin3Tot(43:70);
% 
% Vector = bin4Tot(43:70);
% 
% Vector = bin5Tot(43:70);
% 
% Vector = bin6Tot(43:70);
% 
% Vector = bin7Tot(43:70);
% 
% Vector = bin8Tot(43:70);


% Compute unmethylated fraction for each CpG
VectorNeg=ones(size(Vector))-Vector;


% Step 6: Build a 2xM matrix with:
%   Row 1: probabilities of methylation for each CpG
%   Row 2: probabilities of non-methylation
probArrayVector = [Vector';   % Probability of events happening
                   VectorNeg'];  % Probability of events not happening

M = size(probArrayVector, 2);  % Total number of events


%% Compute probability distribution over the number of methylated CpGs

% Step 7: Loop over all possible numbers of methylated CpGs (from 0 to M)
for j=0:1:size(Vector,1)

% Specify the number of events (CpGs to be methylated) to consider 
N = j;

% Calculate the probability of N events happening and M - N events not happening

eventIndices = nchoosek(1:M, N);  % Generate combinations of event indices

probNEventsMMinusNNot = 0;  % Initialize the probability

% Step 8: Loop over all combinations and compute joint probability    
for i = 1:size(eventIndices, 1)
    
    % Multiply probabilities of methylated and unmethylated events
    probCombination = prod(probArrayVector(1, eventIndices(i, :))) * prod(probArrayVector(2, setdiff(1:M, eventIndices(i, :))));
    
    % Accumulate total probability
    probNEventsMMinusNNot = probNEventsMMinusNNot + probCombination;
end


% Store final probability for j methylated CpGs
probArrayVectorFINAL(j+1)=probNEventsMMinusNNot;

end


%% Plot: Probability distribution of methylated CpG fraction
SIZE = size(Vector,1);

j=0:1:SIZE;
figure(4)
hold on
plot(j/SIZE, probArrayVectorFINAL, 'LineWidth', 2)
xlabel('Fraction of CpGs methylated');
ylabel('Probability');
box on