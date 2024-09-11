% Name of the reaction to optimize
reaction_name = 'EX_gg4abut[e]';  % Adjust this according to the exact name in your model

% Find the index of the reaction in the model
reaction_idx = find(strcmp(model.rxns, reaction_name));

if isempty(reaction_idx)
    error('The reaction "%s" is not found in the model.', reaction_name);
end

% Create a copy of the model to modify the objective function
model4_opt = model;

% Initialize the objective function coefficient vector
model4_opt.c = zeros(length(model4_opt.rxns), 1);  % Initialize with zeros
model4_opt.c(reaction_idx) = 1;  % Maximize the flux of this reaction