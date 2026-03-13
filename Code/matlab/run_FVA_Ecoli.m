changeCobraSolver('gurobi','all');
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('QP', 'feasTol', 1e-9);

% load E. coli model
load('Data/Ecoli/iML1515.mat')

% find the metabolic reactions
filterRxns = @(model) model.rxns(~contains(model.rxns, {'EX_','DM_'}));
rxn_list   = filterRxns(model);

% carbon sources reactions
CS = {'Glucose'; 'Mannitol'; 'Glucosamine'; 'Glycerol'; 'Maltose'; 
    'Gluconate'; 'Xylose'; 'Sorbitol'; 'Ribose'; 'Succinate'; 
    'Galactose'; 'Lactate'; 'Alanine'; 'Pyruvate'; 'Oxoglutarate'; 'Acetate'};

CSrxns = {'EX_glc__D_e'; 'EX_mnl_e'; 'EX_gam_e'; 'EX_glyc_e'; 'EX_malt_e';
    'EX_glcn_e'; 'EX_xyl__D_e'; 'EX_sbt__D_e'; 'EX_rib__D_e'; 'EX_succ_e';
    'EX_gal_e'; 'EX_lac__D_e'; 'EX_ala__D_e'; 'EX_pyr_e'; 'EX_akg_e'; 'EX_ac_e'};

[~, orig_cs_idx] = ismember(CSrxns, model.rxns);

% run FVA
ncpu = 20;
parpool(ncpu);
parfevalOnAll(@maxNumCompThreads, 0, 1); 

% FVA analysis
for i = 1:numel(orig_cs_idx)
    condModel = model;
    condModel.lb(orig_cs_idx) = 0;
    condModel.lb(orig_cs_idx(i)) = -10; % following the constraint from Glucose

    [minFlux, maxFlux] = FVA_analysis(condModel, 10, rxn_list);

    % save tables
    finalTable = table(rxn_list, minFlux, maxFlux, 'VariableNames', {'RxnID', 'minFlux', 'maxFlux'});
    writetable(finalTable,fullfile('Results/FVA/Ecoli', [CS{i},'.csv']));
end
