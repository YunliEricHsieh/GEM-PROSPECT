changeCobraSolver('gurobi','all');
changeCobraSolverParams('LP', 'feasTol', 1e-6);
changeCobraSolverParams('QP', 'feasTol', 1e-6);

ncpu = 20;
if isempty(gcp('nocreate'))
    parpool(ncpu);
    parfevalOnAll(@maxNumCompThreads, 0, 2);
end

%% load different type of models
auto_model = readCbModel('Data/pciCre1355/NDLadpraw_Autotrophic_Rep1.xml');
mixo_model = readCbModel('Data/pciCre1355/NDLadpraw_Mixotrophic_Rep1.xml');
hetero_model = readCbModel('Data/pciCre1355/NDLadpraw_Heterotrophic_Rep1.xml');

% Find exchange and uptake reactions
ex_rxns = auto_model.rxns(findExcRxns(auto_model));
upt_rxns = {};
for i = 1:numel(ex_rxns)
    col = find(strcmp(auto_model.rxns, ex_rxns{i}));
    if sum(auto_model.S(:,col)) == 1 && ~contains(ex_rxns{i}, 'DM') && contains(ex_rxns{i}, 'EX')
        upt_rxns{end+1, 1} = ex_rxns{i};
    elseif sum(auto_model.S(:,col)) > 1
        fprintf('Need to check the stoichiometric matrix for col: %d\n',col);
        return;
    end
end

% find the index for uptake reactions
auto_upt_id = find(ismember(auto_model.rxns, upt_rxns));
hetero_upt_id = find(ismember(hetero_model.rxns, upt_rxns));
mixo_upt_id = find(ismember(mixo_model.rxns, upt_rxns));

% remove the upper bound constraint of uptake reactions
auto_model.ub(auto_upt_id(auto_model.ub(auto_upt_id) ~= 0)) = 1000;
hetero_model.ub(hetero_upt_id(hetero_model.ub(hetero_upt_id) ~= 0)) = 1000;
mixo_model.ub(mixo_upt_id(mixo_model.ub(mixo_upt_id) ~= 0)) = 1000;

% find the metabolic reactions
filterRxns = @(model) model.rxns(~contains(model.rxns, {'EX_','draw','prot_pool'}));
rxn_list = intersect(intersect(filterRxns(auto_model), filterRxns(mixo_model)), filterRxns(hetero_model));

% remove the block reaction in rxn_list
auto_FVA = readtable('Results/FVA/auto_FVA_10p.csv');
mixo_FVA = readtable('Results/FVA/mixo_FVA_10p.csv');
hetero_FVA = readtable('Results/FVA/hetero_FVA_10p.csv');

% identify blocked reactions in the models
getBlocked = @(FVA) FVA.RxnID(FVA.minFlux == 0 & FVA.maxFlux == 0);
blocked_all = unique([getBlocked(auto_FVA); getBlocked(mixo_FVA); getBlocked(hetero_FVA)]);
rxn_list = setdiff(rxn_list, blocked_all);

% create reaction index cell array
rxn_cell = arrayfun(@(x) sprintf('Rxn%d', x), (1:numel(rxn_list))', 'UniformOutput', false);

% filter out demand and sudo reactions (containing 'No' or 'DM_')
valid_idx = ~contains(rxn_list, {'No','DM_'});
rxn_list = rxn_list(valid_idx);
rxn_cell = rxn_cell(valid_idx);

% find reversible reactions (forward and backward)
is_rev = contains(rxn_list, '_REV');
rerxn  = rxn_list(is_rev);
rerxn_cell = rxn_cell(is_rev);

forward_rxns = strrep(rerxn, '_REV', '');
forward_cells = rxn_cell(ismember(rxn_list, forward_rxns));

% combine the lists
cell_bine = [forward_cells; rerxn_cell];
rxn_bine = [forward_rxns; rerxn];

% find irreversible reactions 
rxn_list_ir = setdiff(rxn_list, rxn_bine, 'stable');
rxn_cell_ir = setdiff(rxn_cell, cell_bine, 'stable');

%% modeling the treatment
mutant_table = readtable('Data/Mutant_phenotypes_table_filtered_final.csv');
GO_table = readtable('Data/GO_table_filtered.txt');

% Remove the rows that contain 'antiporter'
GO_table(contains(GO_table.GO_Name, 'antiporter'), :) = [];

enzymeTable = [];
% Loop over each row in the table
for i = 1:height(mutant_table)
    % Split the 'UniProtID' by 'or'
    IDs = split(mutant_table.UniProtID{i}, ' or ');

    % Duplicate the row for each split ID
    for j = 1:numel(IDs)
        enzymeTable = [enzymeTable; IDs(j)];  % Append the new row
    end
end

go_enzyme = intersect(GO_table.UniProtID, enzymeTable);
f_mutant_table = mutant_table(contains(mutant_table.UniProtID, go_enzyme),:);

all_ratios = [f_mutant_table.auto_mixo, f_mutant_table.auto_hetero, f_mutant_table.mixo_hetero, ...
              f_mutant_table.auto_auto_CO2, f_mutant_table.mixo_mixo_CO2, f_mutant_table.mixo_hypo10_mixo, ...
              f_mutant_table.mixo_hypo25_mixo, f_mutant_table.mixo_hypo75_mixo];

all_ratios = 2.^all_ratios;
all_ratios(:, 6:8) = 1 ./ all_ratios(:, 6:8);

% Keep only rows where all 8 conditions are valid (~isnan)
valid_rows = sum(~isnan(all_ratios), 2) == 8;
f_mutant_table = f_mutant_table(valid_rows, :);
all_ratios = all_ratios(valid_rows, :);

%% build condition models
alpha = 0.1;
% change the nutrient uptake
% model1 -> auto_model
% model2 -> mixo_model
% model3 -> hetero_model
[model1, model2, model3] = changeuptake(auto_model, mixo_model, hetero_model);

% create CO2 model
% model4 -> auto_CO2
% model5 -> mixo_CO2
[model4, model1] = create_CO2_model(model1, alpha);
[model5, model2] = create_CO2_model(model2, alpha);

% create hypo model
% hypo10 -> model6
% hypo25 -> model7
% hypo75 -> model8
[model6, model7, model8, model2] = create_hypo_model(model2, alpha);

% put all models together
models = {model1; model2; model3; model4; model5; model6; model7; model8};

%% Reversible reactions
% File to store results
out_rev = 'Results/screens/Max_flux_screen_8_Re.csv';
out_irr = 'Results/screens/Max_flux_screen_8.csv';

fid_rev = fopen(out_rev, 'w');
fprintf(fid_rev, '%s\n', strjoin(['EnzymeID', reshape(forward_cells, 1, [])], ','));
fclose(fid_rev);

fid_irr = fopen(out_irr, 'w');
fprintf(fid_irr, '%s\n', strjoin(['EnzymeID', reshape(rxn_cell_ir, 1, [])], ','));
fclose(fid_irr);

for i = 1:height(f_mutant_table)
    enzyme = f_mutant_table.UniProtID{i};
    ratios = all_ratios(i, :);
    
    % Run Reversible
    flux_screens_8_rev = GEM_PROSPECT_reversible_rxns(models, ratios, rerxn, alpha);
    writeToFile(out_rev, {enzyme, flux_screens_8_rev{:}});
    
    % Run Irreversible
    flux_screens_8_irr = GEM_PROSPECT(models, ratios, rxn_list_ir, alpha);
    writeToFile(out_irr, {enzyme, flux_screens_8_irr{:}});
end