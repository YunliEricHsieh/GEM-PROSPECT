function flux_values = GEM_PROSPECT(models, ratios, rxn_list, alpha)
    
    % find the testing models
    ratio_pairs = [1 2; 1 3; 2 3; 1 4; 2 5; 2 6; 2 7; 2 8];

    active_ratios = find(~isnan(ratios));
    active_models_logical = false(8, 1);

    for k = 1:numel(active_ratios)
        pair = ratio_pairs(active_ratios(k), :);
        active_models_logical(pair(1)) = true;
        active_models_logical(pair(2)) = true;
    end

    active_model_idx = find(active_models_logical);
    numModels = numel(active_model_idx);
    test_models = models(active_model_idx);

    abs_to_test = zeros(8, 1);
    abs_to_test(active_model_idx) = 1:numModels;

    %% prepare the matrix
    nrxn = zeros(numModels, 1);
    S_cell = cell(numModels, 1);
    ub_cell = cell(numModels, 1);
    lb_cell = cell(numModels, 1);
    bio_indices_local = zeros(numModels, 1);
    
    for i = 1:numModels
        bio_idx = find(test_models{i}.c == 1);
        bio_indices_local(i) = bio_idx(1);
        
        % Optimize for max biomass
        opt = optimizeCbModel(test_models{i}, 'max');
        
        % Apply alpha constraint to growth
        test_models{i}.lb(bio_idx(1)) = alpha * opt.f;
        
        % Store structural data
        S_cell{i}  = test_models{i}.S;
        ub_cell{i} = test_models{i}.ub;
        lb_cell{i} = test_models{i}.lb;
        nrxn(i)    = numel(test_models{i}.rxns);
    end

    col_offsets = [0; cumsum(nrxn(1:end-1))];

    A_block = blkdiag(S_cell{:});
    beq = zeros(size(A_block, 1), 1);
    
    % Biomass Ratio Inequality Constraints
    num_ineq = numel(active_ratios) * 2;
    Aineq = sparse(num_ineq, size(A_block, 2));
    bineq = zeros(num_ineq, 1);
    
    for k = 1:numel(active_ratios)
        r_idx = active_ratios(k);
        pair = ratio_pairs(r_idx, :);
        ratio_val = ratios(r_idx);
        w = ratio_val * 0.05;
        
        % Map to global column indices in the super-matrix
        m1 = abs_to_test(pair(1));
        m2 = abs_to_test(pair(2));
        g_bio1 = col_offsets(m1) + bio_indices_local(m1);
        g_bio2 = col_offsets(m2) + bio_indices_local(m2);
        
        row1 = (k-1)*2 + 1;
        row2 = (k-1)*2 + 2;
        
        % v_bio1 - ratio * v_bio2 <= w
        Aineq(row1, g_bio1) = 1;
        Aineq(row1, g_bio2) = -ratio_val;
        bineq(row1) = w;
        
        % -v_bio1 + ratio * v_bio2 <= w
        Aineq(row2, g_bio1) = -1;
        Aineq(row2, g_bio2) = ratio_val;
        bineq(row2) = w;
    end

    % Preallocate equality rows for (v_m1 - v_m2 = 0)
    Aeq1 = sparse(numModels - 1, size(A_block, 2));
    beq1 = zeros(numModels - 1, 1);
    
    % Assemble Final Base Problem
    problem_base.A = [A_block; Aineq; Aeq1];
    problem_base.rhs = [beq; bineq; beq1];
    problem_base.lb = vertcat(lb_cell{:});
    problem_base.ub = vertcat(ub_cell{:});
    problem_base.vtype = repmat('C', size(problem_base.A, 2), 1);
    problem_base.sense = [repmat('=', size(beq, 1), 1); ...
                          repmat('<', size(bineq, 1), 1); ...
                          repmat('=', size(beq1, 1), 1)];

    % Reaction indices
    numRxnsList = numel(rxn_list);
    global_rxn_idx = zeros(numModels, numRxnsList);
    
    for i = 1:numModels
        [~, local_idx] = ismember(rxn_list, test_models{i}.rxns);
        global_rxn_idx(i, :) = local_idx + col_offsets(i);
    end

    %% Parallel Solver Loop
    params.Threads = 2;
    params.OutputFlag = 0;

    flux_values = NaN(numRxnsList, 1);
    row_offset = size(A_block, 1) + size(Aineq, 1);
    num_cols = size(problem_base.A, 2);

    parfor j = 1:numRxnsList
        
        problem_local = problem_base;
        A_local = problem_local.A;
        
        % Assign cross-condition equality constraints: v_m(i) - v_m(i+1) = 0
        for i = 1:numModels-1
            A_local(row_offset + i, global_rxn_idx(i, j))   = 1;
            A_local(row_offset + i, global_rxn_idx(i+1, j)) = -1;
        end
        problem_local.A = A_local;
        
        % --- Step 1: Minimize total flux sum ---
        problem1 = problem_local;
        problem1.modelsense = 'min';
        problem1.obj = ones(num_cols, 1);
        
        min_flux_sum = gurobi(problem1, params);
        
        if strcmp(min_flux_sum.status, 'OPTIMAL')
            
            % --- Step 2: Maximize/Minimize target reaction ---
            problem2 = problem_local;
            
            % Add Step 1 objective as a hard constraint
            problem2.A = [problem2.A; ones(1, num_cols)];
            problem2.rhs = [problem2.rhs; min_flux_sum.objval];
            problem2.sense = [problem2.sense; '='];
            
            % Set objective to the target reaction of the first model
            obj = zeros(num_cols, 1);
            obj(global_rxn_idx(1, j)) = 1;
            problem2.obj = obj;
            problem2.modelsense = 'max';
            
            flux_solution = gurobi(problem2, params);
            
            if strcmp(flux_solution.status, 'OPTIMAL')
                flux_values(j) = flux_solution.objval;
            end
        end
    end
end