changeCobraSolver('gurobi','all');
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('QP', 'feasTol', 1e-9);

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);
parfevalOnAll(@maxNumCompThreads, 0, 2); 
params.Threads = 2;
params.OutputFlag = 0;

% load E. coli model
load('Data/Ecoli/iML1515.mat')
model = creategrRulesField(model);
model = buildRxnGeneMat(model);

% carbon sources reactions
CS = {'Glucose'; 'Mannitol'; 'Glucosamine'; 'Glycerol'; 'Maltose'; 
    'Gluconate'; 'Xylose'; 'Sorbitol'; 'Ribose'; 'Succinate'; 
    'Galactose'; 'Lactate'; 'Alanine'; 'Pyruvate'; 'Oxoglutarate'; 'Acetate'};

CSrxns = {'EX_glc__D_e'; 'EX_mnl_e'; 'EX_gam_e'; 'EX_glyc_e'; 'EX_malt_e';
    'EX_glcn_e'; 'EX_xyl__D_e'; 'EX_sbt__D_e'; 'EX_rib__D_e'; 'EX_succ_e';
    'EX_gal_e'; 'EX_lac__D_e'; 'EX_ala__D_e'; 'EX_pyr_e'; 'EX_akg_e'; 'EX_ac_e'};

[~, orig_cs_idx] = ismember(CSrxns, model.rxns);

% filter the target reactions
filterRxns = @(model) model.rxns(~contains(model.rxns, {'EX_','DM_', 'BIOMASS_Ec'}));
rxn_list = filterRxns(model);

% remove the block reaction in rxn_list
folderPath = 'Results/FVA/Ecoli';
common_zero_rxns = {}; 

% loop through all carbon sources
for i = 1:numel(CS)
    % construct the filename and read the table
    fileName = fullfile(folderPath, sprintf('%s.csv', CS{i}));
    fvaData = readtable(fileName);
    
    % find reactions where both min and max flux are zero
    is_zero_flux = fvaData.minFlux == 0 & fvaData.maxFlux == 0;
    
    % extract the specific RxnIDs for this carbon source
    current_zero_rxns = fvaData.RxnID(is_zero_flux);
    
    if i == 1
        common_zero_rxns = current_zero_rxns;
    else
        common_zero_rxns = intersect(common_zero_rxns, current_zero_rxns);
    end
    
    fprintf('Processed %s: %d zero-flux reactions found.\n', CS{i}, length(current_zero_rxns));
end

rxn_list = setdiff(rxn_list, common_zero_rxns);

% create condition specific model and put all models together
models = cell(numel(orig_cs_idx), 1);

for i = 1:numel(orig_cs_idx)
    condModel = model;
    condModel.lb(orig_cs_idx) = 0;
    condModel.lb(orig_cs_idx(i)) = -10; % following the constraint from Glucose
    models{i} = condModel;
end

% load gene essentiality data
geneEssentiality = readtable('Data/Ecoli/Ecoli_gene_essentiality.csv');
essentialGene = geneEssentiality.Gene_ID(geneEssentiality.Essentiality_0_1 == 1);

% run GEM-PERSPECT to find the max flux value for each targeted rxns
disp('Precomputing structural matrices...');
numModels = numel(models);
nrxn = zeros(numModels, 1);
ub_cell = cell(numModels, 1);
lb_cell = cell(numModels, 1);
S_cell = cell(numModels, 1);
bio_indices_local = zeros(numModels, 1);

for i = 1:numModels
    bio_idx = find(models{i}.c == 1);
    bio_indices_local(i) = bio_idx;

    ub_cell{i} = models{i}.ub;
    lb_cell{i} = models{i}.lb;
    S_cell{i} = models{i}.S;
    nrxn(i) = numel(models{i}.rxns);
end

% upper and lower bound
base_ub = vertcat(ub_cell{:});
base_lb = vertcat(lb_cell{:});

% calculate global column offsets and global biomass indices
col_offsets = [0; cumsum(nrxn(1:end-1))];
global_bio_indices = bio_indices_local + col_offsets;

% build the super-matrix using sparse matrices
A_block = blkdiag(S_cell{:});

% Preallocate the empty rows for the (v_a - v_m = 0) constraints
Aeq1 = sparse(numModels - 1, size(A_block, 2));  
A_base = [A_block; Aeq1];

beq = zeros(size(A_block, 1), 1);
beq1 = zeros(numModels - 1, 1);

% build base Gurobi problem
problem_base.A = A_base;
problem_base.rhs = [beq; beq1];
problem_base.vtype = repmat('C', size(A_base, 2), 1);
problem_base.sense = repmat('=', size(problem_base.rhs, 1), 1);

% pre-map all reaction indices before the parfor loop
numRxnsList = numel(rxn_list);
[~, orig_idx_list] = ismember(rxn_list, model.rxns);

global_rxn_idx = zeros(numModels, numRxnsList);

for j = 1:numRxnsList
    global_rxn_idx(:, j) = orig_idx_list(j) + col_offsets;
end

% test different threshold 
threshold = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
tName = arrayfun(@(x) sprintf('1e%d', log10(x)), threshold, 'UniformOutput', false);

row_offset = size(A_block, 1);
num_cols = size(A_base, 2);

for t = 1:numel(threshold)
    fprintf('Running GEM-PERSPECT for %s (%d/%d)...\n', tName{t}, t, numel(threshold));

    % use a numeric array pre-filled with NaN
    flux_values = NaN(numRxnsList, 1);

    % run the parallel loop
    parfor j = 1:numRxnsList
        % create a fresh copy of the base UB and update only the biomass bounds
        current_ub = base_ub;
        current_lb = base_lb;

        gene_indices = find(model.rxnGeneMat(orig_idx_list(j), :));
        rxn_genes = model.genes(gene_indices);

        [common_genes, ~, ~] = intersect(rxn_genes, essentialGene);

        if ~isempty(common_genes)
            % Essential: Enforce v_bio <= threshold
            current_ub(global_bio_indices) = threshold(t);
        else
            % Non-Essential: Enforce v_bio >= threshold
            current_lb(global_bio_indices) = threshold(t);
        end

        % assign it to the base problem
        problem_local = problem_base;
        problem_local.ub = current_ub;
        problem_local.lb = current_lb;

        A_local = problem_local.A;
        for i = 1:numModels-1
            % v_con1 - v_con2 = 0
            A_local(row_offset + i, global_rxn_idx(i, j)) = 1;
            A_local(row_offset + i, global_rxn_idx(i+1, j)) = -1;
        end
        problem_local.A = A_local;
 
        % --- maximize target reaction ---
        obj = zeros(num_cols, 1);
        obj(global_rxn_idx(1, j)) = 1;
        problem_local.obj = obj;
        
        % check forward direction
        problem_local.modelsense = 'max';              
        sol_max = gurobi(problem_local, params);

        val_max = 0;
        if strcmp(sol_max.status, 'OPTIMAL')
            val_max = sol_max.objval;
        end

        % if the reaction is reversible, check direction to find absolute
        % flux
        if base_lb(global_rxn_idx(1,j)) < 0
            problem_local.modelsense = 'min';
            sol_min = gurobi(problem_local, params);

            val_min = 0;
            if strcmp(sol_min.status, 'OPTIMAL')
                val_min = sol_min.objval;
            end

            flux_values(j) = max(abs(val_max), abs(val_min));
        else
            flux_values(j) = val_max;
        end
    end

    % build final table
    finaltable = table(rxn_list, flux_values, 'VariableNames', {'RxnID', 'MaxFlux'});
    writetable(finaltable, fullfile('Results/screens/Ecoli', ['max_flux_results_', tName{t},'.csv']));
end

%% Test with WT (single) model ##
disp('Precomputing individual structural matrices...');
numModels = numel(orig_cs_idx);
base_probs = cell(numModels, 1);
bio_indices = zeros(numModels, 1);

for i = 1:numModels
    condModel = model;
    % block all uptakes, open only the specific one
    condModel.lb(orig_cs_idx) = 0;
    condModel.lb(orig_cs_idx(i)) = -10; 
    
    % find biomass index
    bio_indices(i) = find(condModel.c == 1);
    
    % build raw Gurobi structure for this condition
    p.A = sparse(condModel.S);
    p.rhs = condModel.b;
    p.sense = repmat('=', size(p.A, 1), 1);
    p.lb = condModel.lb;
    p.ub = condModel.ub;
    p.vtype = repmat('C', size(p.A, 2), 1);
    
    base_probs{i} = p;
end

% run the condition specific model
for t = 1:numel(threshold)
    fprintf('Running condition specific model for %s (%d/%d)...\n', tName{t}, t, numel(threshold));
    
    % initialize a matrix to store [Reaction x Condition] results
    flux_matrix = NaN(numRxnsList, numModels);
    
    parfor j = 1:numRxnsList
        rxn_idx = orig_idx_list(j);
        
        % essentiality check
        gene_indices = find(model.rxnGeneMat(rxn_idx, :));
        rxn_genes = model.genes(gene_indices);
        [common_genes, ~, ~] = intersect(rxn_genes, essentialGene);
        is_essential = ~isempty(common_genes);
        
        % array to hold max fluxes for this reaction across all 13 conditions
        temp_cond_fluxes = zeros(1, numModels);
        
        for c = 1:numModels
            % pull the pre-built problem for condition 'c'
            prob_local = base_probs{c};
            
            % apply thresholds
            if is_essential
                prob_local.ub(bio_indices(c)) = threshold(t);
            else
                prob_local.lb(bio_indices(c)) = threshold(t);
            end
            
            % maximize target reaction
            obj = zeros(size(prob_local.A, 2), 1);
            obj(rxn_idx) = 1;
            prob_local.obj = obj;
            
            % forward check
            prob_local.modelsense = 'max';              
            sol_max = gurobi(prob_local, params);
            val_max = 0;
            if strcmp(sol_max.status, 'OPTIMAL')
                val_max = sol_max.objval;
            end
            
            % reverse check (if allowed)
            if prob_local.lb(rxn_idx) < 0
                prob_local.modelsense = 'min';
                sol_min = gurobi(prob_local, params);
                val_min = 0;
                if strcmp(sol_min.status, 'OPTIMAL')
                    val_min = sol_min.objval;
                end
                temp_cond_fluxes(c) = max(abs(val_max), abs(val_min));
            else
                temp_cond_fluxes(c) = val_max;
            end
        end
        
        % Store the 13 condition results into the main matrix row
        flux_matrix(j, :) = temp_cond_fluxes;
    end
    
    % build a beautiful final table with 13 data columns
    T_data = array2table(flux_matrix, 'VariableNames', CS');
    T_IDs = table(rxn_list, 'VariableNames', {'RxnID'});
    finaltable = [T_IDs, T_data];
    
    writetable(finaltable, fullfile('Results/screens/Ecoli/CS', ['max_flux_results_', tName{t},'.csv']));
end