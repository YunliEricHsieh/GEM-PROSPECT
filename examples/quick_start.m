function demo = quick_start()
%QUICK_START Run two candidate reactions for one supplied mutant phenotype.
%   DEMO is a two-row RxnID/Direction/MaxFlux table. Uses the existing case-
%   study condition helpers and optimization functions, with alpha=0.10 and
%   their fixed omega=0.05. Writes only Results/examples/quick_start.csv.
%   Run setup_gem_prospect after initializing COBRA/Gurobi before calling.
%   This small demonstration is not a complete association/evaluation run.

previous_dir = pwd;
restore_dir = onCleanup(@() cd(previous_dir));
setup_gem_prospect();
changeCobraSolver('gurobi', 'all');
changeCobraSolverParams('LP', 'feasTol', 1e-6);
changeCobraSolverParams('QP', 'feasTol', 1e-6);

auto_model = readCbModel('Data/pciCre1355/NDLadpraw_Autotrophic_Rep1.xml');
mixo_model = readCbModel('Data/pciCre1355/NDLadpraw_Mixotrophic_Rep1.xml');
hetero_model = readCbModel('Data/pciCre1355/NDLadpraw_Heterotrophic_Rep1.xml');

% Mirror the main driver's original positive-uptake selection and bound step.
ex_rxns = auto_model.rxns(findExcRxns(auto_model));
upt_rxns = {};
for uptake_index = 1:numel(ex_rxns)
    col = find(strcmp(auto_model.rxns, ex_rxns{uptake_index}));
    if sum(auto_model.S(:, col)) == 1 && ...
            ~contains(ex_rxns{uptake_index}, 'DM') && ...
            contains(ex_rxns{uptake_index}, 'EX')
        upt_rxns{end+1, 1} = ex_rxns{uptake_index}; %#ok<AGROW>
    elseif sum(auto_model.S(:, col)) > 1
        error('GEM_PROSPECT:CheckStoichiometry', ...
            'Review exchange reaction %s before continuing.', ex_rxns{uptake_index});
    end
end
base_models = {auto_model; mixo_model; hetero_model};
for model_index = 1:3
    uptake_ids = find(ismember(base_models{model_index}.rxns, upt_rxns));
    open_ids = uptake_ids(base_models{model_index}.ub(uptake_ids) ~= 0);
    base_models{model_index}.ub(open_ids) = 1000;
end

alpha = 0.10;
[model1, model2, model3] = changeuptake(base_models{1}, base_models{2}, base_models{3});
[model4, model1] = create_CO2_model(model1, alpha);
[model5, model2] = create_CO2_model(model2, alpha);
[model6, model7, model8, model2] = create_hypo_model(model2, alpha);
models = {model1; model2; model3; model4; model5; model6; model7; model8};

phenotypes = readtable('Data/Mutant_phenotypes_table_filtered_final.csv');
mutant_id = 'A0A2K3E4Q7';
row = phenotypes(strcmp(phenotypes.UniProtID, mutant_id), :);
assert(height(row) == 1, 'Expected one supplied phenotype row for %s.', mutant_id);
ratios = 2.^[row.auto_mixo, row.auto_hetero, row.mixo_hetero, ...
    row.auto_auto_CO2, row.mixo_mixo_CO2, row.mixo_hypo10_mixo, ...
    row.mixo_hypo25_mixo, row.mixo_hypo75_mixo];
ratios(6:8) = 1 ./ ratios(6:8);
assert(all(isfinite(ratios) & ratios > 0), 'Example ratios must be finite and positive.');

% These IDs occur in the supplied candidate maps and archived screen tables.
rxn_ir = '3SPYRSPh';
rxn_forward = '10FTHFtx';
rxn_backward = '10FTHFtx_REV';
for model_index = 1:8
    assert(all(ismember({rxn_ir, rxn_forward, rxn_backward}, models{model_index}.rxns)), ...
        'Example reaction IDs must occur in every condition model.');
end
flux_ir = GEM_PROSPECT(models, ratios, {rxn_ir}, alpha);
flux_rev = GEM_PROSPECT_reversible_rxns(models, ratios, {rxn_backward}, alpha);
demo = table(repmat({mutant_id}, 2, 1), {rxn_ir; rxn_forward}, ...
    {'irreversible'; 'reversible_net_magnitude'}, [flux_ir{1}; flux_rev{1}], ...
    'VariableNames', {'EnzymeID', 'RxnID', 'Direction', 'MaxFlux'});
if ~isfolder('Results/examples')
    mkdir('Results/examples');
end
writetable(demo, 'Results/examples/quick_start.csv');
if any(~isfinite(demo.MaxFlux))
    warning('GEM_PROSPECT:NonfiniteExampleFlux', ...
        'A solve did not return finite flux. Inspect feasibility and solver setup.');
end
end
