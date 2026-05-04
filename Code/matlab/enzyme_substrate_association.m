%% load model
model = readCbModel('Data/pciCre1355/NDLadpraw_Autotrophic_Rep1.xml');

% Load the CSV file containing potential gene-reaction associations
geneAssoc = readtable('Results/potential_gene_reaction_associations_for_CataPro.csv', 'PreserveVariableNames', true);

% Extract lists directly from the table
rxn_list = geneAssoc.RxnID;
enzyme_list = geneAssoc.enzyme_candidates;

% Find all reaction indices at once
[~, rxnIndices] = ismember(rxn_list, model.rxns);

disp('Extracting LHS substrates and mapping enzymes...');
RxnID       = {};
BiGG_ID     = {};
Substrate   = {};
MetFormula  = {};
Pathway     = {};
Enzyme      = {};
Status      = {};

% Regex to remove stoichiometric coefficients (e.g., "2 ", "1.5 ")
coefRE  = '^\d+(\.\d+)?\s+'; 

for i = 1:numel(rxn_list)
    rxnID = rxn_list{i};
    ridx  = rxnIndices(i);
    
    % Safety check: skip if the reaction isn't actually in the model
    if ridx == 0
        fprintf('Warning: %s not found in the model. Skipping...\n', rxnID);
        continue;
    end
    
    % Handle subsystem (pathway)
    pathway = model.subSystems{ridx};
    if iscell(pathway)
        pathway = strjoin(pathway, '; '); 
    end
    
    % Split the string of enzymes into a cell array
    enzymes = strtrim(strsplit(enzyme_list{i}, ';'));
    
    % Extract formula
    formula = printRxnFormula(model, rxnID, false);
    if iscell(formula), formula = formula{1}; end
    
    % Split by the reaction arrow
    sides = strtrim(strsplit(formula, '->'));
    
    % Strictly process the Left-Hand Side (LHS) as Substrates
    % Split by ' + ' (with spaces) to protect 'H+' metabolites
    LHS_raw = strtrim(strsplit(sides{1}, ' + '));
    LHSmets = regexprep(LHS_raw, coefRE, '');
    
    % Remove empty strings if any crept in
    LHSmets(cellfun(@isempty, LHSmets)) = [];
    
    % Map metabolites and duplicate rows for EACH enzyme
    for m = 1:numel(LHSmets)
        metID = LHSmets{m};
        metIdx = find(strcmp(model.mets, metID), 1);
        
        % Lookup data, handling edge cases where metabolite is missing
        if isempty(metIdx)
            metName = 'Unknown';
            metForm = 'Unknown';
        else
            metName = model.metNames{metIdx};
            metForm = model.metFormulas{metIdx};
        end
        
        % Nested loop: create a row for every enzyme candidate
        for e = 1:numel(enzymes)
            if isempty(enzymes{e}), continue; end % Skip trailing empty splits
            
            RxnID{end+1, 1}      = rxnID;
            BiGG_ID{end+1, 1}    = metID;
            Substrate{end+1, 1}  = metName;
            MetFormula{end+1, 1} = metForm;
            Pathway{end+1, 1}    = pathway;
            Enzyme{end+1, 1}     = enzymes{e};
            Status{end+1, 1}     = 'wild'; 
        end
    end
end

disp('Building table and saving to CSV...');

T = table(RxnID, BiGG_ID, Substrate, MetFormula, Pathway, Enzyme, Status, ...
    'VariableNames', {'RxnID', 'BiGG_ID', 'Substrate', 'MetFormula', 'Pathway', 'Enzyme', 'type'});

% writetable(T, 'Results/CataPro/enzymes_and_substrates.csv');
disp('Success! File saved.');

%% find the reaction in the target pathway
% Extract unique pathways dynamically (removes duplicates and empties)
targetsubSystem = unique(T.Pathway);
targetsubSystem(cellfun(@isempty, targetsubSystem)) = []; 

disp('Extracting reference reactions for targeted pathways...');
RxnID       = {};
BiGG_ID     = {};
Enzyme      = {};
MetFormula  = {};
Substrate   = {};
Pathway     = {};
Type        = {};

% Regex patterns
pattern1 = 'x\((\d+)\)';        % To extract gene indices from model.rules
coefRE   = '^\d+(\.\d+)?\s+';   % To safely remove stoichiometric coefficients

% Flatten subSystems in case COBRA nested them in cell arrays
flatSubSystems = cell(numel(model.subSystems), 1);

for k = 1:numel(model.subSystems)
    sys = model.subSystems{k};
    
    if iscell(sys)
        % Force the contents into strings before joining, preventing the error!
        flatSubSystems{k} = char(strjoin(string(sys), '; '));
    elseif ischar(sys) || isstring(sys)
        % If it is already text, just keep it as is
        flatSubSystems{k} = char(sys);
    else
        % If it is a missing NaN or [], gracefully replace it with a blank string
        flatSubSystems{k} = ''; 
    end
end

for i = 1:numel(targetsubSystem)
    
    % Find all reactions belonging to this specific pathway instantly
    idx_logic = contains(flatSubSystems, targetsubSystem{i});
    rxn_list  = model.rxns(idx_logic);
    
    % Filter out pseudo/demand reactions
    rxn_list = rxn_list(~contains(rxn_list, 'No') & ~startsWith(rxn_list, 'DM_'));
    
    for j = 1:numel(rxn_list)
        rxnID = rxn_list{j};
        ridx  = find(strcmp(model.rxns, rxnID)); % Get numeric index
        
        % Process only if the reaction has Gene-Protein-Reaction (GPR) rules
        if ~isempty(model.rules{ridx})
            
            % --- 1. Parse Substrates (LHS) safely ---
            formula = printRxnFormula(model, rxnID, false);
            if iscell(formula), formula = formula{1}; end
            
            sides   = strtrim(strsplit(formula, '->'));
            LHS_raw = strtrim(strsplit(sides{1}, ' + ')); % Protects 'H+'
            LHSmets = regexprep(LHS_raw, coefRE, '');
            LHSmets(cellfun(@isempty, LHSmets)) = [];
            
            allSides = LHSmets(:); % Force column vector
            
            % Map metabolite indices safely
            metIdx = cellfun(@(m) find(strcmp(model.mets, m), 1), allSides, 'UniformOutput', false);
            valid  = ~cellfun(@isempty, metIdx); % Identify found metabolites
            allSides = allSides(valid);
            metIdx   = cell2mat(metIdx(valid));
            
            if isempty(allSides), continue; end % Skip if no valid substrates
            
            % --- 2. Parse Enzymes (Genes) safely ---
            % Extract the numeric gene IDs from the rule string (e.g., "x(123)")
            eidx    = regexp(model.rules{ridx}, pattern1, 'tokens');
            if isempty(eidx), continue; end
            
            numbers = cellfun(@(x) str2double(x{1}), eidx);
            enzyme  = unique(model.genes(numbers));
            enzyme  = enzyme(:); % Force column vector
            
            % --- 3. Append Data Using Matrix Expansion ---
            n = numel(allSides);
            m = numel(enzyme);
            
            RxnID      = [RxnID;       repmat({rxnID}, n*m, 1)];
            BiGG_ID    = [BiGG_ID;     repmat(allSides, m, 1)];
            MetFormula = [MetFormula;  repmat(model.metFormulas(metIdx), m, 1)];
            Substrate  = [Substrate;   repmat(model.metNames(metIdx), m, 1)];
            Pathway    = [Pathway;     repmat(targetsubSystem(i), n*m, 1)];
            Enzyme     = [Enzyme;      repelem(enzyme, n, 1)];
            Type       = [Type;        repmat({'wild'}, n*m, 1)];
        end
    end
end

disp('Building reference table and saving to CSV...');
T = table(RxnID, BiGG_ID, Substrate, MetFormula, Pathway, Enzyme, Type, ...
    'VariableNames', {'RxnID', 'BiGG_ID', 'Substrate', 'MetFormula', 'Pathway', 'Enzyme', 'type'});

% writetable(T, 'Results/CataPro/enzymes_and_substrates_for_reference.csv');
disp('Success! Reference file saved.');
