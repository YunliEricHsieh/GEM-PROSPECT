# Small example and input template

The example uses the supplied model and phenotype files in `Data/` and the core functions in `Code/matlab/function/`.

From the root run `python3 scripts/check_inputs.py --workflow quick-start`. Initialize MATLAB, COBRA Toolbox, Gurobi and Parallel Computing Toolbox as described in the [main README](../README.md), then:

```matlab
addpath(fullfile(pwd, 'scripts'));
setup_gem_prospect();
addpath(fullfile(pwd, 'examples'));
demo = quick_start();
disp(demo);
```

`quick_start.m` selects supplied phenotype ID `A0A2K3E4Q7` (`GeneID=Cre01.g001800`), converts the eight log2 measurements to the main driver's oriented linear ratios, and constructs its eight condition models using the existing helpers. It then calls `GEM_PROSPECT` for `3SPYRSPh` and `GEM_PROSPECT_reversible_rxns` for `10FTHFtx_REV` / `10FTHFtx`. These two reactions are in the supplied candidate maps. All condition-bound changes and α=0.10 mirror the main driver; ω remains the original 0.05 inside the core functions.

Only those two candidates are evaluated. The example does not rerun FVA, sample fluxes, classify associations or estimate performance. `Results/examples/quick_start.csv` has `EnzymeID,RxnID,Direction,MaxFlux` and two rows; existing manuscript files are untouched. The core's `parfor` may start a pool using your MATLAB profile if no pool is already open. The example does not request the manuscript scripts' 20-worker pool.

Inspect `MaxFlux` and solver messages before interpreting a result. `NaN` indicates unsuccessful inference under the existing handling. A zero can be a valid computed flux; reversible second-stage non-optimal solves may also be encoded as zero. The example demonstrates the optimization interface, while association calls require reference sampling and postprocessing as described in the [main workflow](../docs/reproducing_manuscript.md).

`phenotype_template.csv` is a **header-only** copy of the supplied final-table schema. Values should be dimensionless log2 phenotype ratios. It is not an input dataset until populated and validated. Its extra NaCl/P/N columns match the original table, even though main inference uses only the eight documented comparisons. See [input formats](../docs/input_formats.md) for order, orientation, missing data and ID mapping. New conditions that do not fit the current pair graph need a reviewed code adaptation.
