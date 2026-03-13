function flux_values = GEM_PROSPECT_reversible_rxns(models, ratios, rxn_list, alpha)

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

   % Preallocate equality rows for forward and backward reactions(v_m1 - v_m2 = 0)
    Aeq1 = sparse(2 * (numModels - 1), size(A_block, 2));
    beq1 = zeros(2 * (numModels - 1), 1);
    
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
    global_b_idx = zeros(numModels, numRxnsList);
    global_f_idx = zeros(numModels, numRxnsList);
    
    for j = 1:numRxnsList
        b_id = rxn_list{j};
        f_id = strrep(b_id, '_REV', ''); % Extract forward ID
        
        for i = 1:numModels
            [~, b_loc] = ismember(b_id, test_models{i}.rxns);
            [~, f_loc] = ismember(f_id, test_models{i}.rxns);
            
            global_b_idx(i, j) = b_loc + col_offsets(i);
            global_f_idx(i, j) = f_loc + col_offsets(i);
        end
    end

    %% loop all the reactions
    params.Threads = 2;
    params.OutputFlag = 0;

    flux_values = num2cell(NaN(1, numRxnsList));
    row_offset = size(A_block, 1) + size(Aineq, 1);
    num_cols = size(problem_base.A, 2);

    parfor j = 1:numRxnsList
        
        problem_local = problem_base;
        A_local = problem_local.A;
        
        % Assign cross-condition equality constraints for BOTH directions
        for i = 1:numModels-1
            % v_f(i) - v_f(i+1) = 0
            A_local(row_offset + i, global_f_idx(i, j))   = 1;
            A_local(row_offset + i, global_f_idx(i+1, j)) = -1;
            
            % v_b(i) - v_b(i+1) = 0
            A_local(row_offset + (numModels-1) + i, global_b_idx(i, j))   = 1;
            A_local(row_offset + (numModels-1) + i, global_b_idx(i+1, j)) = -1;
        end
        problem_local.A = A_local;
        
        % --- Step 1: Minimize total flux sum (pFBA) ---
        problem1 = problem_local;
        problem1.modelsense = 'min';
        problem1.obj = ones(num_cols, 1);
        
        min_flux_sum = gurobi(problem1, params);
        
        if strcmp(min_flux_sum.status, 'OPTIMAL')
            
            % --- Step 2: Extract Maximum Absolute Net Flux ---
            problem2 = problem_local;
            
            % Add Step 1 objective as a hard constraint
            problem2.A = [problem2.A; ones(1, num_cols)];
            problem2.rhs = [problem2.rhs; min_flux_sum.objval];
            problem2.sense = [problem2.sense; '='];
            
            % Define the objective: Net Flux = v_f - v_b (Using Model 1)
            obj = zeros(num_cols, 1);
            obj(global_f_idx(1, j)) = 1;
            obj(global_b_idx(1, j)) = -1;
            problem2.obj = obj;
            
            % Maximize (pushes flux forward)
            problem2.modelsense = 'max';
            sol_max = gurobi(problem2, params);
            val_max = 0;
            if strcmp(sol_max.status, 'OPTIMAL'), val_max = sol_max.objval; end
            
            % Minimize (pushes flux backward)
            problem2.modelsense = 'min';
            sol_min = gurobi(problem2, params);
            val_min = 0;
            if strcmp(sol_min.status, 'OPTIMAL'), val_min = sol_min.objval; end
            
            % Save the largest magnitude
            flux_values{j} = max(abs(val_max), abs(val_min));
        end
    end
end