%% CODE TO GENERATE FIGURE 3H
% This code generates the probability distribution of the fraction of CpGs methylated.
% It first loads bisulfite sequencing data from a CSV file, selects the region of interest (e.g., promoter), and computes
% the distribution of the number of CpG sites methylated across that
% region.

%% Load and prepare data

% Step 1: Save your Excel file as a CSV (e.g., bisulfitedata2.csv)

% Step 2: Read the CSV file as a table 
tableData = readtable('bisulfitedata2.csv', 'VariableNamingRule', 'preserve');

% Step 3: Convert the table to a numeric matrix for easier access
matrixData = table2array(tableData);


%% Average methylation levels across replicates

TOT = 175;  % Total number of regions (e.g., promoter tiles, genomic segments)

ActiveTot = zeros(TOT,1);
IntermediateTot = zeros(TOT,1);
RepressedTot = zeros(TOT,1);
% 


% Step 4: Compute the mean across 3 biological replicates for each chromatin state

for i = 1:TOT
    ActiveTot(i) = mean([matrixData(i,3), matrixData(i,5), matrixData(i,7)]);
    IntermediateTot(i) = mean([matrixData(i,15), matrixData(i,17), matrixData(i,19)]);
    RepressedTot(i) = mean([matrixData(i,9), matrixData(i,11), matrixData(i,13)]);
end



%% Select region of interest

% Step 5: Select the region (e.g., promoter + rTetR/dCas9 binding sites)
% You can switch between ActiveTot / IntermediateTot / RepressedTot as needed

Vector = ActiveTot(12:39);  % Promoter + dCas9/rTEtR sites

% Vector = IntermediateTot(12:39); % Promoter + dCas9/rTEtR sites

% Vector = RepressedTot(12:39); % Promoter + dCas9/rTEtR sites



% Compute unmethylated fraction for each CpG
VectorNeg = 1 - Vector;



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
