# Reproducing the existing analyses

These instructions call the existing implementation. They distinguish recomputing outputs from inspecting supplied references. Run in a separate clone if you want to retain archived result files, since case-study scripts write to fixed paths. All original algorithmic constants are retained. Install dependencies as described in the [README](../README.md).

## Setup and input availability

In MATLAB, initialize COBRA and Gurobi, navigate to the repository root, then:

```matlab
addpath(fullfile(pwd, 'scripts'));
setup_gem_prospect();
```

The wrappers create required output directories. The original scripts otherwise assume these directories already exist. They request 20 workers, and Gurobi uses one or two threads per solve; this requires an adequate local parallel-pool configuration. The E. coli wrapper closes an existing pool before its FVA stage because the original FVA script opens a pool unconditionally. Its inference script also closes/reopens the pool itself. For a smaller machine, review/change `ncpu` explicitly and record the change; the wrapper does not silently adjust it. The C. reinhardtii `all` and `sampling` stages check that COBRA's `gpSampler` is available before running.

Input availability can be checked without any scientific dependencies:

```bash
python3 scripts/check_inputs.py --workflow quick-start
python3 scripts/check_inputs.py --workflow chlamydomonas
python3 scripts/check_inputs.py --workflow ecoli
```

The checks validate required paths/headers and result-index joins where relevant. Solver feasibility and installed dependencies are checked when running the MATLAB, Python or R stages.

## C. reinhardtii: start from the supplied processed phenotype table

The required inputs are three model XMLs, the final phenotype CSV, and filtered GO TSV. FVA is needed for candidate selection; reference sampling is needed for downstream flux-reduction analysis. The current main code selects 884 complete activity-annotated rows from the 1,832-row final phenotype table. The archived inference tables also have 884 rows; this agreement is an input audit, not a newly reproduced optimization result.

| Order | Existing script / wrapper stage | Main outputs |
|---|---|---|
| 1 | `run_FVA.m` / `reproduce_chlamydomonas('fva')` | Three `Results/FVA/*_FVA_10p.csv` tables |
| 2 | `obtain_reaction_info_from_models.m` / `'reaction-info'` | Five `Data/Reactions/*.csv` mapping/GPR tables |
| 3 | `GEM_PROSPECT_implementation.m` / `'prospect'` | `Results/screens/Max_flux_screen_8.csv`, `Max_flux_screen_8_Re.csv` |
| 4 | `flux_sampling.m` / `'sampling'` | Six Auto/Hetero/Mixo reference sampling tables |
| 5 | `Code/R/identify gene reaction associations.R` | `Results/potential_gene_reaction_associations.csv` |
| 6 | `Code/R/Fig2.R`, `Fig3.R` | `Results/figures/Fig 2.svg`, `Fig 3.svg` |

Steps 3 and 4 independently depend on FVA; both must finish before steps 5–6. `reproduce_chlamydomonas('all')` runs steps 1–4 in that order. Source-profile preprocessing and R steps are separate. The direct equivalent in MATLAB is:

```matlab
run_FVA;
obtain_reaction_info_from_models;
GEM_PROSPECT_implementation;
flux_sampling;
```

From the repository root, export candidates and figures:

```bash
Rscript scripts/reproduce_chlamydomonas.R associations
Rscript scripts/reproduce_chlamydomonas.R figures
```

To analyze the **supplied** results, skip steps 1–4 and run those two R commands directly. That regenerates postprocessing/plots from archived numerical inputs; it does not rerun inference. Newly generated indices, inference columns, and sampling references must be kept together. `Results/max_flux/` contains older/reference copies and is not the active inference output directory; `Max_flux_screens_8.csv` there has a different candidate-column set.

### R missing-data handling

The association and figure scripts use the same `read_and_clean` helper. It converts infinite fluxes to missing values and negative values to zero, then removes rows whose numeric reaction values are all missing. Next, it removes numeric reaction columns whose values are all missing across the retained rows. The first ID column is retained, and zero fluxes count as valid values. Filtering preserves the data-frame structure even when only one column remains. With the supplied inference tables, both helpers retain 864 of 884 rows; irreversible reaction columns decrease from 834 to 831, and reversible columns from 593 to 574. Downstream zero-flux hit counts treat removed reaction columns as having no hits; this does not change the benchmark reaction lists or recall denominators.

This corrects an earlier indexing error that used a column-derived mask to select rows and did not save the row-filter result. Archived R association tables and figures were not regenerated with this correction; rerun the R postprocessing commands above to produce outputs using the corrected cleaning. The MATLAB inference tables themselves are not altered by R cleaning. Check the cleaning regression tests with `Rscript tests/test_read_and_clean.R` from the repository root; these tests require only base R.

### Condition construction and parameters

The main script opens already nonzero positive uptake bounds to 1000, removes reactions blocked in any base-model FVA table, and excludes exchange/protein-usage/demand and `No`-containing reactions. It splits reverse IDs using `_REV`. It exponentiates the eight phenotype columns and reciprocates the last three. It then executes:

1. `changeuptake`: nutrient upper bounds derived from an optimization at 99% of optimal growth.
2. `create_CO2_model` for Auto then Mixo: adds 3% of the difference between uptake at optimum and at α×optimum to the original upper reference bound; also returns modified base models.
3. `create_hypo_model` on the modified Mixo model: returns the three TAP-fraction models and another modified base Mixo model.
4. `GEM_PROSPECT_reversible_rxns` and `GEM_PROSPECT` for each selected protein row, with α=0.10 and fixed ω=0.05.

Use this order, including the helper-returned updated base models. Details and caveats are in [parameters](parameters.md). The minimization is the sum of all flux variables, including protein pseudo-reactions; it is not a general absolute-value reformulation.

### Preprocessing the source phenotype spreadsheet

The source spreadsheet contains already processed log2 abundance ratios from the original study. The following sequence filters and aggregates these profiles:

```bash
python3 scripts/check_inputs.py --workflow preprocessing
python3 Code/python/pairwise_correlation.py
Rscript "Code/R/filter genes with consistent phynotype.R"
```

Then in MATLAB run `UniProtIDs_mapping` after setup. Then from the root:

```bash
Rscript "Code/R/calculate mean of phynotypes from each screens.R"
```

Those steps overwrite supplied intermediate/final tables. The Spearman threshold is 0.5, the minimum number of different insertion-mutant profiles per gene is three, and the R final selection is strictly greater than 0.5. This row-level consistency filter is separate from averaging repeated-screen columns. Screen names reflect automatic import-name conversion; see [input formats](input_formats.md). The supplied final table is the shortest starting point for standalone inference.

`Code/python/goontology.py` separately retrieves annotations from QuickGO and appends to `Data/GO_table.txt` and `GO_table_filtered.txt`. It is not required when using supplied GO data. Repeating it appends duplicate rows and changes annotation provenance over time; use a separate data copy if researching updated annotations. It is not run automatically by the wrappers.

### Optional evaluations using the supplied main results

| Analysis | Modeling prerequisites | R script / output |
|---|---|---|
| Single-protein reaction subset | Supplied main screen/reference tables, existing GPR maps | `FigS2.R` → `Sup Fig 2.svg` |
| Unique gene–reaction subset | Same base tables; the script further excludes repeated enzyme entries in its subset | `FigS3.R` → `Sup Fig 3.svg` |

After installing the main figure packages, run from the repository root:

```bash
Rscript -e 'dir.create("Results/figures", recursive = TRUE, showWarnings = FALSE)'
Rscript Code/R/FigS2.R
Rscript Code/R/FigS3.R
```

The other supplementary figure scripts require additional condition samples or parameter-comparison outputs. They are independent of standalone inference and the main `all` workflows. Optional annotation scripts and their inputs are listed in [repository map](repository_map.md).

## E. coli analysis

Required inputs: `Data/Ecoli/iML1515.mat` with variable `model`, corresponding SBML for Python, and the semicolon-separated essentiality table. FVA and modeling use the 16 source IDs in [parameters](parameters.md); each opens only the selected one of those sources with lower bound −10, leaving other model settings as supplied.

```matlab
reproduce_ecoli('all');  % run_FVA_Ecoli followed by GEM_PROSPECT_Ecoli
```

Separate stages are `'fva'` and `'prospect'`. Each FVA table has `RxnID`, `minFlux`, `maxFlux`. Candidate reactions exclude `EX_`, `DM_`, and `BIOMASS_Ec`; only reactions blocked in **all** source FVA tables are removed. The script treats a reaction as essential when any associated gene belongs to the essential list. It bounds essential biomass above by τ, and nonessential biomass below by τ, across all conditions.

The τ loop runs `1e-7` to `1e-1` and writes seven `RxnID,MaxFlux` tables under `Results/screens/Ecoli/`. It then performs individual-condition controls and writes seven 16-condition tables under `Results/screens/Ecoli/CS/`. Joint models impose equal target flux across sources, handle signed reverse flux, and do not minimize total flux first. Non-optimal solves contribute zero in this implementation.

Evaluate from the root:

```bash
python3 Code/python/Accuracy_Ecoli.py
```

The script compares `abs(MaxFlux) ≤ 1e-7` to two reaction-label mappings: any essential associated gene, and dependence after collectively knocking out essential genes through Boolean GPR logic. It excludes reactions without mapped GPR labels. NaN flux is treated as zero during evaluation. It writes `Results/figures/GEM_PERSPECT_accuracy_comparison.png`; metrics remain in Python dataframes rather than being exported to a separate CSV. The `CS/` controls are not consumed by this evaluator.

Essentiality labels inform the growth constraints as well as the evaluation labels. Interpret accuracy as association recovery under this labeled benchmark setup; the script does not create a held-out test set.

## Troubleshooting and reproducibility limits

- Missing function: initialize COBRA/Gurobi first, then run `setup_gem_prospect`; use `which ... -all` to inspect shadowing.
- Parallel pool failure: check Parallel Computing Toolbox and the profile's 20-worker capacity. Close an old pool for direct `run_FVA_Ecoli` use.
- Missing file/header: run the appropriate input checker; confirm working directory and the exact names/capitalization in [input formats](input_formats.md).
- `NaN`, `Inf`, or unexplained zeros: inspect feasibility and solver status before biological interpretation. Original outputs do not preserve individual solver statuses.
- Failed R aggregation from other input files: verify the mapped CSV's actual headers against the hard-coded screen names and the [documented import conversion](input_formats.md). The supplied mapped table's headers and final averages have been checked against the existing script.
- Different candidate calls after sampling: `gpSampler(model)` uses sampler defaults and the script does not reset the random seed. Preserve archived references for comparisons to archived results; record the random state and any changed options for future analyses.

For numerical reproduction, run all required stages in a working MATLAB/COBRA/Gurobi installation. Input checks establish file and identifier consistency; they do not execute optimization. Retain supplied reference files for comparisons and record solver settings when generating new results.
