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

% find exchange and uptake reactions
ex_rxns = auto_model.rxns(findExcRxns(auto_model));
upt_rxns = {};

for i = 1:numel(ex_rxns)
    col = find(strcmp(auto_model.rxns, ex_rxns{i}));
    % Use nnz (number of non-zero metabolites) to safely identify exchanges
    if nnz(auto_model.S(:,col)) == 1 && startsWith(ex_rxns{i}, 'EX_') && ~startsWith(ex_rxns{i}, 'DM_')
        upt_rxns{end+1, 1} = ex_rxns{i};
    elseif nnz(auto_model.S(:,col)) > 1
        fprintf('Warning: Check stoichiometry for col: %d (%s)\n', col, ex_rxns{i});
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

% find optimize growth
auto_opt = optimizeCbModel(auto_model);
hetero_opt = optimizeCbModel(hetero_model);
mixo_opt = optimizeCbModel(mixo_model);

auto_model.lb(auto_index) = 0.9 * auto_opt.f;
auto_model.ub(auto_index) = 0.9 * auto_opt.f;

hetero_model.lb(hetero_index) = 0.9 * hetero_opt.f;
hetero_model.ub(hetero_index) = 0.9 * hetero_opt.f;

mixo_model.lb(mixo_index) = 0.9 * mixo_opt.f;
mixo_model.ub(mixo_index) = 0.9 * mixo_opt.f;

[a_sampleStruct, a_mixedFraction] = gpSampler(auto_model);
[h_sampleStruct, h_mixedFraction] = gpSampler(hetero_model);
[m_sampleStruct, m_mixedFraction] = gpSampler(mixo_model);

%% create the rxn index
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
is_rev = endsWith(rxn_list, '_REV');
rerxn  = rxn_list(is_rev);
rerxn_cell = rxn_cell(is_rev);

forward_rxns = strrep(rerxn, '_REV', '');
forward_cells = rxn_cell(ismember(rxn_list, forward_rxns));

re_rxn_list = [forward_rxns; rerxn];
re_rxn_cells = [forward_cells; rerxn_cell];

% irreversible reactions
rxn_list_ir = setdiff(rxn_list, re_rxn_list, 'stable');
rxn_cell_ir = setdiff(rxn_cell, re_rxn_cells, 'stable');

%% find the flux distribution for irreversible reactions
process_irr_samples(auto_model, a_sampleStruct.points, rxn_cell_ir, rxn_list_ir, 'Results/flux_sampling/auto_sampling.csv');
process_irr_samples(hetero_model, h_sampleStruct.points, rxn_cell_ir, rxn_list_ir, 'Results/flux_sampling/hetero_sampling.csv');
process_irr_samples(mixo_model, m_sampleStruct.points, rxn_cell_ir, rxn_list_ir, 'Results/flux_sampling/mixo_sampling.csv');

%% find the flux distribution for reversible reactions
process_rev_samples(auto_model, a_sampleStruct.points, forward_cells, forward_rxns, rerxn, 'Results/flux_sampling/auto_sampling_re.csv');
process_rev_samples(hetero_model, h_sampleStruct.points, forward_cells, forward_rxns, rerxn, 'Results/flux_sampling/hetero_sampling_re.csv');
process_rev_samples(mixo_model, m_sampleStruct.points, forward_cells, forward_rxns, rerxn, 'Results/flux_sampling/mixo_sampling_re.csv');

%% locat function
function process_irr_samples(model, samples, rxn_ids, rxn_names, out_path)
    [~, idx] = ismember(rxn_names, model.rxns);
    
    % Calculate mean/std across rows (dimension 2) instantly
    meanFlux = mean(samples(idx, :), 2);
    stdFlux  = std(samples(idx, :), 0, 2);
    
    T = table(rxn_ids, rxn_names, meanFlux, stdFlux, 'VariableNames', {'RxnIndex', 'RxnID', 'meanFlux', 'stdFlux'});
    writetable(T, out_path);
end

function process_rev_samples(model, samples, f_ids, f_names, b_names, out_path)
    [~, idx_f] = ismember(f_names, model.rxns);
    [~, idx_b] = ismember(b_names, model.rxns);
    
    % Subtract full matrices and take absolute value instantly
    net_flux = abs(samples(idx_f, :) - samples(idx_b, :));
    
    meanFlux = mean(net_flux, 2);
    stdFlux  = std(net_flux, 0, 2);
    
    T = table(f_ids, f_names, meanFlux, stdFlux, 'VariableNames', {'RxnIndex', 'RxnID', 'meanFlux', 'stdFlux'});
    writetable(T, out_path);
end