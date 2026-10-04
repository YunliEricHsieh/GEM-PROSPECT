# Parameters in the existing implementation

This catalogue records the implemented settings, code locations, and stages affected by a change. Values were selected for the supplied case studies. Calibrate α, ω, phenotype QC and postprocessing thresholds against suitable benchmarks when applying the framework to another dataset.

Paths are relative to the repository root. Find the exact definitions with `rg -n 'alpha|omega|threshold|0\.03|0\.99|0\.9|0\.5' Code`.

| Parameter / meaning | Current value and values present in code | Where defined | Scope / calibration evidence | What to rerun if changed |
|---|---|---|---|---|
| `alpha`: minimum biomass fraction `v_bio ≥ alpha * optimum` | Main 0.10; manuscript comparisons 0.01, 0.05. Fraction conceptually in [0,1], but no input validator | `GEM_PROSPECT_implementation.m`; `alpha` argument of the two core functions in `function/` | Ratio-based framework parameter; empirically compared at 1%, 5%, 10% for the supplied data | Inference and downstream evaluation. If condition-helper α also changes, rebuild conditions and relevant sampling |
| `omega`: additive biomass-relation slack `abs(v_i - ratio*v_j) ≤ ratio*omega` | Fixed 0.05; manuscript comparisons 0.01, 0.10 | `function/GEM_PROSPECT.m`, `GEM_PROSPECT_reversible_rxns.m`: `w = ratio_val * 0.05` | Ratio-based framework parameter; empirical comparisons apply to the supplied data. Changing ω requires editing separate copies of both functions | Inference and evaluation; no phenotype preprocessing or condition rebuilding |
| FVA biomass percentage | 10 (%); biomass **lb and ub** both set to `0.1*optimum` | `run_FVA.m`, `run_FVA_Ecoli.m`; `function/FVA_analysis.m` | Both case studies; distinct from inference α. No FVA percentage sweep supplied | FVA, reaction selection/maps, inference, sampling and evaluation |
| CO₂ uptake range adjustment | `max_co2 + (max_co2-min_co2)*0.03` | `function/create_CO2_model.m` | **C. reinhardtii case-study specific**. No general calibration search implemented | Condition construction, inference, condition sampling and evaluation |
| Growth levels for CO₂/nutrient reference bounds | At optimum and at `alpha*optimum`; condition-builder α=0.10 | `create_CO2_model.m`, `create_hypo_model.m`; caller `GEM_PROSPECT_implementation.m` | **C. reinhardtii case-study specific**. Values come from fluxes in single optimum solutions, not separate uptake FVA extrema | All derived condition models and dependent analyses |
| Nutrient-cap reference growth | Biomass fixed at 0.99*optimum; selected uptake `ub` set from resulting solution | `function/changeuptake.m` | **C. reinhardtii case-study specific**; no sensitivity sweep supplied | Condition models, inference, condition sampling and evaluation |
| Hypoosmotic fractions | 0.10, 0.25, 0.75; interpolate from min to max uptake. If min=max, multiply min by fraction | `function/create_hypo_model.m` | **C. reinhardtii case-study specific**, corresponding to manuscript conditions | Condition models, inference, condition sampling and evaluation |
| Initial open uptake cap | 1000 for already nonzero uptake upper bounds; closed uptakes remain closed | `GEM_PROSPECT_implementation.m`, `run_FVA.m`, `flux_sampling.m` | **C. reinhardtii case-study specific**, positive split uptake variables | Affected FVA/models/sampling, then dependent inference and evaluation |
| Reference sampling growth | Biomass lb=ub=0.9*optimum in each of Auto/Mixo/Hetero | `flux_sampling.m` | **C. reinhardtii reference/evaluation specific**. `gpSampler(model)` uses toolbox defaults without option overrides or an explicit seed reset | Sampling and downstream candidate calls/figures; no need to rerun inference if its inputs are unchanged |
| Substantial flux reduction | `log10(abs(MaxFlux)/abs(meanFlux)) < -1`; recall scripts sweep −1 through −5 | `Code/R/identify gene reaction associations.R` (`1:1`, `-th`); `Fig2.R`, `Fig3.R`, `FigS2.R`–`FigS5.R` | **C. reinhardtii evaluation specific**. Manuscript describes selection via sensitivity analysis. “One-fold” means a log10 threshold here: <0.1 of reference, strict inequality | R postprocessing/evaluation only; no modeling or raw preprocessing |
| Candidate count filter | Between **1 and 19 inclusive**, separately for Auto/Hetero/Mixo, followed by candidate intersection | `identify gene reaction associations.R`, `between(..., 1, 19)` | **C. reinhardtii candidate-export specific** | Candidate export only |
| `threshold` (τ): essential/nonessential biomass cutoff | Actual loop order `1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1`; essential reaction `ub=τ`, other reaction `lb=τ` | `GEM_PROSPECT_Ecoli.m` | **E. coli case-study specific** absolute growth-flux cutoff, not fraction of optimal growth. Manuscript reports stabilization around 1e-5; no automatic selection implemented | E. coli joint and individual controls, then accuracy evaluation |
| Carbon-source uptake | All listed carbon-source `lb=0`; current source `lb=-10` in original model units | `GEM_PROSPECT_Ecoli.m`, `run_FVA_Ecoli.m` | **E. coli case-study specific**; same uptake magnitude used for every source | FVA, inference and controls, then evaluation |
| E. coli evaluation zero tolerance | Script invocation `1e-7`; helper `run_full_analysis` default `1e-5`; positive prediction when `abs(flux) ≤ tolerance` | `Code/python/Accuracy_Ecoli.py` | **E. coli evaluation specific**; not τ and not solver feasibility tolerance | Python evaluation only |
| Profile consistency (Spearman) | Python trims when mean `<0.5` for ≥4 profiles; R final filter requires mean **>0.5** | `pairwise_correlation.py`; `filter genes with consistent phynotype.R` | **C. reinhardtii preprocessing specific**. Used to select consistent profiles; no universal optimality or threshold tuning established | Raw filtering, correlation table, mapping, aggregation and all dependent analyses |
| Minimum mutant profiles per gene | 3; trimming stops at 3 | `pairwise_correlation.py` (`len(x) >= 3`, loop ≥4) | **C. reinhardtii library specific**. Different insertion-mutant profiles grouped by gene; repeated-screen measurements are separate columns | Full phenotype preprocessing and dependent inference/evaluation |
| Insertion QC | Remove `Feature == "3'UTR"`; retain `Confidence level ≤4`; drop all-missing profiles | `pairwise_correlation.py` | **C. reinhardtii library specific**, follows supplied metadata conventions | Full phenotype preprocessing and dependent inference/evaluation |
| GO/activity selection and ratio completeness | Retrieval retains names containing `activity`; main script removes GO names containing `antiporter`, intersects IDs, requires 8 non-NaN ratios | `goontology.py`, `GEM_PROSPECT_implementation.m` | **C. reinhardtii case-study specific**; a different organism needs its own annotation policy | Annotation/phenotype selection, inference/evaluation; raw replicate QC need not change |
| Fixed ratio graph | `[1 2;1 3;2 3;1 4;2 5;2 6;2 7;2 8]`; 8 model slots | `GEM_PROSPECT.m`, `GEM_PROSPECT_reversible_rxns.m` | Current interface assumption, no arbitrary graph argument | Changing topology requires code review and revalidation, not just configuration |
| COBRA LP/QP feasibility tolerance | 1e-6 for C. reinhardtii, 1e-9 for E. coli | Case scripts' `changeCobraSolverParams` calls | Numerical setting. Direct `gurobi` structs **do not** set `params.FeasibilityTol`, so they use the installed solver's setting | Affected optimization stages and evaluation; capture settings and compare results |
| Parallelism | `ncpu=20`; direct inference `params.Threads=2`, FVA=1; per-worker MATLAB thread caps 1 or 2 | Case scripts and core functions | Execution resource settings; not calibrated biology | Rerun only required computation; numerical repeatability should be assessed if resource settings change |

Unless prefixed, MATLAB script locations above are in `Code/matlab/`. R/Python filenames are in their named directories.

## Hard-coded condition IDs

The C. reinhardtii condition helpers use `EX_co2_e_REV` for CO₂ and the nutrient list:

```text
EX_pi_e_REV, EX_nh4_e_REV, EX_so4_e_REV, EX_fe2_e_REV,
EX_mg2_e_REV, EX_na1_e_REV, EX_ac_e_REV
```

The E. coli source name → reaction mapping is identical in FVA and inference:

| Source | Reaction | Source | Reaction |
|---|---|---|---|
| Glucose | EX_glc__D_e | Ribose | EX_rib__D_e |
| Mannitol | EX_mnl_e | Succinate | EX_succ_e |
| Glucosamine | EX_gam_e | Galactose | EX_gal__D_e |
| Glycerol | EX_glyc_e | Lactate | EX_lac__D_e |
| Maltose | EX_malt_e | Alanine | EX_ala__D_e |
| Gluconate | EX_glcn_e | Pyruvate | EX_pyr_e |
| Xylose | EX_xyl__D_e | Oxoglutarate | EX_akg_e |
| Sorbitol | EX_sbt__D_e | Acetate | EX_ac_e |

Missing exchange IDs are not handled gracefully by the original scripts. Validate membership before substituting models. Preserve source names consistently because they determine FVA filenames and control-table columns.

## Selecting settings for another dataset

Distinguish phenotype QC, condition construction, inference α/ω and postprocessing cutoffs when calibrating a new dataset. Change one stage at a time, preserve its reference settings, and rerun the dependent stages listed above. Core functions accept α and use fixed internal ω=0.05; [the adaptation guide](applying_to_new_dataset.md#step-6--run-the-model) shows the supported call signatures.

Use known associations to compare feasibility, coverage and prediction behavior across a recorded parameter grid. Evaluate on independent evidence when available. Record settings, code changes, benchmark splits, and all rerun stages. The repository provides no automatic parameter search or universal optimal values.
