function repo_root = setup_gem_prospect()
%SETUP_GEM_PROSPECT Add active project paths and check modeling dependencies.
%   ROOT = SETUP_GEM_PROSPECT() resolves this checkout, changes the working
%   directory to ROOT (existing scripts use relative paths), and adds only
%   Code/matlab and Code/matlab/function. Initialize COBRA and Gurobi first.
%   This function does not install tools or change scientific parameters.

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
addpath(fullfile(repo_root, 'Code', 'matlab'));
addpath(fullfile(repo_root, 'Code', 'matlab', 'function'));

required = {'readCbModel', 'optimizeCbModel', 'findExcRxns', ...
    'changeCobraSolver', 'changeCobraSolverParams', 'gurobi'};
for setup_index = 1:numel(required)
    assert(~isempty(which(required{setup_index})), ...
        'GEM_PROSPECT:MissingDependency', ...
        'Missing %s. Initialize COBRA/Gurobi as described in README.md.', ...
        required{setup_index});
end
assert(license('test', 'Distrib_Computing_Toolbox'), ...
    'GEM_PROSPECT:MissingParallelToolbox', ...
    'The existing implementation requires Parallel Computing Toolbox.');
end
