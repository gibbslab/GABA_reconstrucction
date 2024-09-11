loaded_model = load('modelo_actualizado_gaba.mat');
model = loaded_model.model;
% Define the metabolite of interest with the suffix [e]
metabolite_of_interest = 'gg4abut[e]';

% Find the index of the metabolite in the list of metabolites
met_idx = find(strcmp(model.mets, metabolite_of_interest));

% Check if the metabolite exists in the model
if isempty(met_idx)
    error('The metabolite %s is not found in the model.', metabolite_of_interest);
end

% Find the reactions that contain the metabolite
% For each reaction, check if the metabolite is present (non-zero coefficient)
reactions_with_metabolite = {};
for i = 1:length(model.rxns)
    % Get the coefficients of the metabolites for reaction i
    reaction_metabolites = model.S(:, i);
    % Check if the coefficient of the metabolite of interest is non-zero
    if reaction_metabolites(met_idx) ~= 0
        reactions_with_metabolite{end+1} = model.rxns{i}; %#ok<AGROW>
    end
end

% Display the reactions that contain the metabolite
if isempty(reactions_with_metabolite)
    disp(['No reactions were found that contain the metabolite ', metabolite_of_interest]);
else
    disp(['Reactions that contain the metabolite ', metabolite_of_interest, ':']);
    disp(reactions_with_metabolite);
end