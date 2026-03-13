changeCobraSolver('gurobi','all');
changeCobraSolverParams('LP', 'feasTol', 1e-6);
changeCobraSolverParams('QP', 'feasTol', 1e-6);

ncpu = 20;
if isempty(gcp('nocreate'))
    parpool(ncpu);
    parfevalOnAll(@maxNumCompThreads, 0, 1); 
end

%% load different type of models
auto_model   = readCbModel('Data/pciCre1355/NDLadpraw_Autotrophic_Rep1.xml');
mixo_model   = readCbModel('Data/pciCre1355/NDLadpraw_Mixotrophic_Rep1.xml');
hetero_model = readCbModel('Data/pciCre1355/NDLadpraw_Heterotrophic_Rep1.xml');

ex_rxns = auto_model.rxns(findExcRxns(auto_model));

upt_rxns = {};
for i = 1:numel(ex_rxns)
    col = find(strcmp(auto_model.rxns, ex_rxns{i}));
    
    % Use nnz (number of non-zero elements) to mathematically check for exchanges
    if nnz(auto_model.S(:,col)) == 1 && startsWith(ex_rxns{i}, 'EX_') && ~startsWith(ex_rxns{i}, 'DM_')
        upt_rxns{end+1, 1} = ex_rxns{i};
    elseif nnz(auto_model.S(:,col)) > 1
        fprintf('Warning: Need to check the stoichiometric matrix for col: %d (%s)\n', col, ex_rxns{i});
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

auto_rxns   = filterRxns(auto_model);
mixo_rxns   = filterRxns(mixo_model);
hetero_rxns = filterRxns(hetero_model);

rxn_list = intersect(intersect(auto_rxns, mixo_rxns), hetero_rxns);

%% run FVA
% FVA analysis
[auto_minFlux, auto_maxFlux] = FVA_analysis (auto_model, 10, rxn_list);
[mixo_minFlux, mixo_maxFlux]  = FVA_analysis (mixo_model, 10, rxn_list);
[hetero_minFlux, hetero_maxFlux] = FVA_analysis (hetero_model, 10, rxn_list);

% save tables
auto_FVA = table(rxn_list, auto_minFlux, auto_maxFlux, 'VariableNames', {'RxnID', 'minFlux', 'maxFlux'});
mixo_FVA = table(rxn_list, mixo_minFlux, mixo_maxFlux, 'VariableNames', {'RxnID', 'minFlux', 'maxFlux'});
hetero_FVA = table(rxn_list, hetero_minFlux, hetero_maxFlux, 'VariableNames', {'RxnID', 'minFlux', 'maxFlux'});

writetable(auto_FVA,'Results/FVA/auto_FVA_10p.csv');
writetable(mixo_FVA,'Results/FVA/mixo_FVA_10p.csv');
writetable(hetero_FVA,'Results/FVA/hetero_FVA_10p.csv');