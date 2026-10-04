# Script roles and general/case-study boundaries

The project uses `Code/`, `Data/`, and `Results/` for computation, inputs, and outputs. The `scripts/` entry points call the case-study scripts; `examples/` demonstrates a small call; `docs/` supplies operating instructions.

## Modeling functions in Code/matlab/function

| Function | Purpose and interface boundary |
|---|---|
| `GEM_PROSPECT` | Ratio-constrained, cross-condition candidate evaluation, then flux-sum minimization and irreversible target maximization; fixed eight-slot pair graph |
| `GEM_PROSPECT_reversible_rxns` | Same ratio setup; links both split directions, maximizes absolute net flux magnitude |
| `FVA_analysis` | Fixed-biomass per-reaction minima/maxima; direct Gurobi; returns cells |
| `changeuptake` | C. reinhardtii nutrient-cap construction, not a generic media importer |
| `create_CO2_model`, `create_hypo_model` | C. reinhardtii condition definitions; return variants and updated original model |
| `writeToFile`, `convertToString` | Simple CSV row serialization, with no quoting/escaping |

Reusable optimization functionality remains subject to its model representation, objective, and pair-graph assumptions. Protein constraints are inherited from loaded models; the core does not build a pcGEM or predict kinetic constants.

## Active MATLAB scripts

| Script | Category | Inputs → outputs / role |
|---|---|---|
| `GEM_PROSPECT_implementation.m` | C. reinhardtii main execution | Models, final phenotypes, GO and FVA → two `Results/screens/Max_flux_screen_8*.csv` tables |
| `GEM_PROSPECT_Ecoli.m` | E. coli main execution | MAT model, essentiality and 16 FVA tables → joint τ sweep plus `CS/` controls |
| `run_FVA.m` | C. reinhardtii preprocessing | Three model XMLs → three FVA CSVs |
| `run_FVA_Ecoli.m` | E. coli preprocessing | MAT model → 16 FVA CSVs |
| `UniProtIDs_mapping.m` | C. reinhardtii phenotype preprocessing | Correlation-filtered mutant CSV and mapping TSV → mapped mutant CSV |
| `obtain_reaction_info_from_models.m` | Reaction preprocessing/evaluation metadata | Models, FVA, final phenotypes, GO → five reaction/GPR maps |
| `flux_sampling.m` | Reference baseline/evaluation | Three models, FVA → six sampling-summary CSVs, at exactly 90% biomass optimum |
| `enzyme_substrate_association.m` | Supplementary CataPro preparation | Manually curated candidate table and model → workspace table `T`; both CSV-writing calls are commented out |

## Active Python and R scripts

| Script | Category / dependencies | Role |
|---|---|---|
| `Code/python/pairwise_correlation.py` | Phenotype preprocessing: NumPy, pandas, openpyxl, SciPy | Insertion QC, ≥3 profiles, iterative Spearman consistency and two CSV outputs |
| `Code/R/filter genes with consistent phynotype.R` | Phenotype preprocessing: base R | Retain strictly >0.5 mean Spearman |
| `Code/R/calculate mean of phynotypes from each screens.R` | Phenotype preprocessing: dplyr, stringr | Hard-coded screen aggregation and gene means → final CSV |
| `Code/R/identify gene reaction associations.R` | Candidate export: dplyr, tidyr, purrr | Flux reduction vs three sampling references, count filtering, candidate intersection → association CSV |
| `Code/R/Fig2.R`, `Fig3.R` | Main evaluation: packages in README | Phenotype/flux distributions and three recall views → Figs. 2–3 SVGs |
| `Code/python/Accuracy_Ecoli.py` | E. coli evaluation: pandas, COBRApy, NumPy, matplotlib | τ-wise accuracy against two essentiality/GPR mappings → PNG |
| `Code/R/FigS2.R`, `FigS3.R` | Optional supplementary evaluation | Supplied main screen/reference tables and restricted GPR subsets → supplementary SVGs |
| `Code/R/FigS1.R`, `FigS4.R`, `FigS5.R` | Additional research analyses | Require condition-specific samples or parameter-comparison results in addition to the main inputs; outside the documented main workflow |
| `Code/python/goontology.py` | Supplementary annotation: requests | QuickGO retrieval and activity-term filtering; appends GO TSVs |
| `Code/python/get_seq_uniport_catapro.py` | Supplementary annotation: requests, pandas | Sequence and PubChem structure retrieval; updates candidate substrate CSV |
| `Code/python/get_seq_uniport_transporter.py` | Supplementary annotation: requests, pandas, openpyxl | Sequence retrieval; overwrites transporter XLSX |
| `Code/python/get_SMILE.py` | Supplementary annotation: requests, pandas, openpyxl | PubChem lookup; updates transporter XLSX; may substitute InChI when SMILES is absent |
| `Code/python/get_seq_uniprot_and_catrpro_prediction.py` | Supplementary scratch workflow: requests, pandas | Reference sequence retrieval mixed with a shell command; **not valid standalone Python** as written |
| `Code/python/run_catapro.py` | Supplementary inference launcher | Runs local `CataPro/inference/predict.py` on a curated table; requires separate model weights/environment |
| `Code/R/CataPro.R` | Supplementary evaluation: dplyr, ggplot2, patchwork | Reads supplied kcat/Km predictions → `CataPro.svg` |
| `Code/R/SPOT.R` | Supplementary evaluation: dplyr, ggplot2, readxl | Reads supplied external SPOT predictions → `SPOT.svg`; does not execute SPOT |

The SPOT script reads externally generated predictions. The CataPro launcher requires a separate CataPro installation, model weights, and input columns compatible with that tool. Sequence/substrate preparation requires curated input tables. These annotation workflows are independent of GEM-PROSPECT flux inference and are not invoked by the wrappers.

## Operating entry points

| File | Role |
|---|---|
| `scripts/setup_gem_prospect.m` | Resolve checkout, set root working directory, add active paths, check core MATLAB dependencies |
| `scripts/check_inputs.py` | Read-only path/schema/index checks using standard-library Python |
| `scripts/reproduce_chlamydomonas.m` | Stages `fva`, `reaction-info`, `prospect`, `sampling`; `all` runs these four in order |
| `scripts/reproduce_ecoli.m` | FVA and essentiality benchmark stages; creates `CS/` output directory |
| `scripts/reproduce_chlamydomonas.R` | Stages `associations`, `figures`; `all` runs candidate export and main figures |
| `examples/quick_start.m` | One supplied phenotype, two candidate reactions using the original optimization functions |

The wrappers run the documented main workflows. Optional annotation and supplementary evaluation scripts are run separately with their own required inputs.
