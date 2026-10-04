function reproduce_ecoli(stage)
%REPRODUCE_ECOLI Call existing 16-condition FVA and essentiality benchmark.
%   STAGE is 'all', 'fva', or 'prospect' (default). Outputs under Results/FVA/
%   Ecoli and Results/screens/Ecoli (including CS/) are overwritten. The
%   original tau values, reaction selection, bounds and LP logic are retained.
%   Initialize COBRA/Gurobi first; run Python evaluation separately.

if nargin < 1
    stage = 'prospect';
end
stage = validatestring(stage, {'all', 'fva', 'prospect'});
previous_dir = pwd;
restore_dir = onCleanup(@() cd(previous_dir));
setup_gem_prospect();

output_dirs = {'Results/FVA/Ecoli', 'Results/screens/Ecoli', ...
    'Results/screens/Ecoli/CS'};
for dir_index = 1:numel(output_dirs)
    if ~isfolder(output_dirs{dir_index})
        mkdir(output_dirs{dir_index});
    end
end

if any(strcmp(stage, {'all', 'fva'}))
    % The existing FVA script creates a pool unconditionally.
    delete(gcp('nocreate'));
    fprintf('Running existing script: run_FVA_Ecoli\n');
    execute_existing_script('run_FVA_Ecoli');
end
if any(strcmp(stage, {'all', 'prospect'}))
    fprintf('Running existing script: GEM_PROSPECT_Ecoli\n');
    execute_existing_script('GEM_PROSPECT_Ecoli');
end
end

function execute_existing_script(script_name)
eval(script_name);
end
