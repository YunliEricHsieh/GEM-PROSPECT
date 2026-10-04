function reproduce_chlamydomonas(stage)
%REPRODUCE_CHLAMYDOMONAS Call existing modeling scripts without changing settings.
%   STAGE: 'all', 'fva', 'reaction-info', 'prospect' (default), or 'sampling'.
%   'all' runs FVA, metadata, inference,
%   and three-mode sampling. Raw phenotype preprocessing and R analysis are
%   separate. Existing output files are overwritten by the original scripts.
%   COBRA/Gurobi must be initialized; original scripts request 20 workers.

if nargin < 1
    stage = 'prospect';
end
stage = validatestring(stage, {'all', 'fva', 'reaction-info', 'prospect', 'sampling'});
previous_dir = pwd;
restore_dir = onCleanup(@() cd(previous_dir));
setup_gem_prospect();
if any(strcmp(stage, {'all', 'sampling'}))
    assert(~isempty(which('gpSampler')), 'GEM_PROSPECT:MissingSampler', ...
        'COBRA gpSampler is required for reference sampling. Initialize COBRA or reuse the supplied references.');
end

output_dirs = {'Results/FVA', 'Data/Reactions', 'Results/screens', ...
    'Results/flux_sampling'};
for dir_index = 1:numel(output_dirs)
    if ~isfolder(output_dirs{dir_index})
        mkdir(output_dirs{dir_index});
    end
end

switch stage
    case 'all'
        stages = {'run_FVA', 'obtain_reaction_info_from_models', ...
            'GEM_PROSPECT_implementation', 'flux_sampling'};
    case 'fva'
        stages = {'run_FVA'};
    case 'reaction-info'
        stages = {'obtain_reaction_info_from_models'};
    case 'prospect'
        stages = {'GEM_PROSPECT_implementation'};
    case 'sampling'
        stages = {'flux_sampling'};
end

for stage_index = 1:numel(stages)
    fprintf('Running existing script: %s\n', stages{stage_index});
    execute_existing_script(stages{stage_index});
end
end

function execute_existing_script(script_name)
% Give each original script its own workspace so its loop variables cannot
% overwrite the wrapper's stage list. Model/result files remain its outputs.
eval(script_name);
end
