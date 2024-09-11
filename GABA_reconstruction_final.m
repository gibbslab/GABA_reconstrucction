%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GABAergic Neuron Metabolic Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Introduction:
% GABAergic neurons are a class of neurons that primarily release gamma-aminobutyric acid (GABA), the main inhibitory neurotransmitter in 
% the mammalian central nervous system. These neurons play a crucial role in regulating neuronal excitability and are involved in various 
% physiological and pathological processes, including epilepsy, anxiety, and sleep regulation. The following script focuses on the integration, 
% correction, and refinement of a tissue-specific metabolic model for GABAergic neurons.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Integration Details %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data Import and Conversion:
% First, you import the transcriptome data and the generic metabolic model. Then, you convert Ensembl IDs to Entrez IDs using a mapping file.

% Mapping Expressions:
% You map the gene expressions to reactions using the mapExpressionToReactions function. This step associates gene expression data with reactions in the model.

% Model Reduction:
% After mapping the expressions, you define parameters for model reduction and use the createTissueSpecificModel function to reduce the model to a tissue-specific one.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Corrections and Clarifications %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definition of 'red':
% Ensure that the variable 'red' is defined before using it in the weights option. It could represent weights for reactions in model reduction. If it's not defined, you should define it appropriately.

% Custom Functions:
% Ensure that functions like mapExpressionToReactions and createTissueSpecificModel are correctly implemented and available in your workspace.

% Input Variables:
% Verify that the variables expresion and Recon3D are properly loaded and contain the expected data.

close all force; clear variables; clc

%% Initialize the COBRA Toolbox
initCobraToolbox(false); % Don't update the toolbox
changeCobraSolver('gurobi', 'all', 1); % Use Gurobi solver for large models

%% Import the transcriptome data and the model
% Load expression data
expresion_Gaba = readtable('GSE115565_norm_counts_TPM_GRCh38.p13_NCBI_2.xls');

% Load generic metabolic model
load('Recon3D_301.mat');

% Convert the first column (IDs) to string cells
EnrezIDs = cellstr(num2str(expresion_Gaba{:, 1}));
expresion_2 = addvars(expresion_Gaba, EnrezIDs, 'after', 1, 'NewVariableNames', 'EntrezIDs_2Strig');

%% Data Cleaning
% Ensure that the 'EntrezIDs_2Strig' column is a character cell array
columnName = 'EntrezIDs_2Strig';

if iscell(expresion_2.(columnName))
    expresion_2.(columnName) = cellstr(expresion_2.(columnName));
    expresion_2.(columnName) = strtrim(expresion_2.(columnName));
else
    expresion_2.(columnName) = strtrim(expresion_2.(columnName));
end

% Display the cleaned data
disp(head(expresion_2.(columnName)));
disp(expresion_2(1:5, :));

% Convert all columns except the first to numeric values
for k = 3:width(expresion_2)
    expresion_2.(k) = str2double(expresion_2{:, k});
end

% Remove the first column after conversion
expresion_2(:, 1) = [];

%% Convert Ensembl to Entrez IDs using https://www.ensembl.org/biomart/martview
f = fopen('mart_export.txt');
C = textscan(f, '%s %s', 'HeaderLines', 1);
fclose(f);

Entrez = regexprep(strIDs, "\.[0-9]*", '\s+', '', "");

disp(Entrez(1:5)); % Display first 5 rows to verify

EnsemblIDs_mart_export = C{1};
EntrezIDs_mart_export = C{2};
EntrezIDs_2Strig = C{2};
EntrezIDs_exp = expresion.(2);

% Map IDs
for i = 1:numel(EntrezIDs_2Strig_exp)
    match = strcmp(EntrezIDs_2Strig, EntrezIDs_2Strig_exp{i});
    
    if any(match)
        EnsemblIDs{i} = EnsemblIDs_mart_export{match};
    else
        EnsemblIDs{i} = ''; 
    end
end

% Add Entrez IDs to the expression table
expresion = addvars(expresion_2, EnsemblIDs, 'before', 1);

height(expresion_2.EntrezIDs_2Strig)
EntrezIDs_2Strig = cellstr(expresion.EntrezIDs_2Strig)
expresion.EntrezIDs_2Strig = EntrezIDs_2Strig

%%% Review
disp('First Recon3D IDs:');
disp(Recon3D.genes(1:5));
disp(expresion.EntrezIDs_2Strig(1:5));

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Map expressions to reactions
Recon3D.genes = regexprep(Recon3D.genes, "\.[0-9]*", "");

expresion_columns = expresion.Properties.VariableNames;
expressionRxns1 = cell(length(expresion_columns) - 2, 1);
FinalExpresion = table();

% Loop over expression columns starting from the third column
for k = 3:length(expresion_columns)
    expressionToMap = expresion(:, [2 k]);
    expressionToMap.Properties.VariableNames = {'gene', 'value'};
    expressionToMap.gene(cellfun('isempty', expressionToMap.gene)) = {' '};
    
    if iscell(expressionToMap.value)
        expressionToMap.value = str2double(expressionToMap.value);
    end
    
    [expressionRxns, parsedGPR] = mapExpressionToReactions(Recon3D, expressionToMap);
    FinalExpresion = addvars(FinalExpresion, expressionRxns, 'NewVariableNames', expresion_columns{k});
end

% Validate the number of columns in FinalExpresion
if size(FinalExpresion, 2) ~= (length(expresion_columns) - 2)
    error('The number of columns in FinalExpresion does not match the expected number.');
end

FinalExpresion.Properties.VariableNames = expresion_columns(3:end);

data_struct.tissues = FinalExpresion.Properties.VariableNames(2:end);
data_struct.genes = FinalExpresion.genes;

data_struct.tissues = FinalExpresion.Properties.VariableNames(2:end)
data_struct.genes = expresion.EntrezIDs_2Strig;
data_struct.levels = table2array(FinalExpresion(:, 2:end));
data_struct.threshold = 1;
data_struct

disp('Column names in FinalExpresion:');
disp(FinalExpresion.Properties.VariableNames);

% Remove NaN genes
validGeneIdx = ~ismissing(FinalExpresion.gene);
validFinalExpresion = FinalExpresion(validGeneIdx, :);

defaultWeight = -1;
weights = defaultWeight * ones(length(Recon3D.genes), 1);

geneExpressionMap = containers.Map(validFinalExpresion.gene, validFinalExpresion.SecondColumn);

% Map weights using 'FinalExpresion'
for i = 1:length(Recon3D.genes)
    gene = Recon3D.genes{i};
    if isKey(geneExpressionMap, gene)
        expressionLevel = geneExpressionMap(gene);
        weights(i) = expressionLevel;
    end
end

options.weights = weights;
options.solver = 'INIT'; 
options.runtime = 14400; 

disp('Configured options:');
disp(options);

% Generate the GABAergic neuron model
GABANEURONModel_03 = createTissueSpecificModel(Recon3D, options);

numReactions = length(GABANEURONModel_03.rxns);
disp(['Number of reactions in the model: ', num2str(numReactions)]);

disp('First 5 reactions in the neuron model:');
disp(GABANEURONModel_03.rxns(1:min(5, length(GABANEURONModel_03.rxns))));

cd('C:\Users\u0500\Documents\MATLAB'); 

filename = 'GABANEURONModel_03.mat';
save(filename, 'GABANEURONModel_03');

disp(['Model saved at: ', fullfile(pwd, filename)]);

%%%% Model Curation

% Identify blocked reactions
[consistModel, BlockedRxns] = identifyBlockedRxns(model);

model = load('GABANEURONModel_03_consistent.mat').consistModel;

epsilon = getCobraSolverParams('LP', 'feasTol') * 100;
[consistModel, BlockedRxns] = identifyBlockedRxns(model, epsilon);

save('GABANEURONModel_03_consistent.mat', 'consistModel');

%%% Mass Balance Verification

loadedData = load('GABANEURONModel_03.mat');
model = loadedData.model; 

[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, Elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(model, 1, 'GABANEURONModel_03_');

disp('Mass imbalance by reaction:');
disp(massImbalance);


