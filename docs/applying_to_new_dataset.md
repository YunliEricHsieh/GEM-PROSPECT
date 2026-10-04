# Applying GEM-PROSPECT to a new organism or dataset

This guide describes an adaptation of the existing implementation, not a validated turnkey workflow for arbitrary species. Start by running the [small example](../examples/README.md) with the supplied files. Preserve the manuscript scripts and inputs, and put a new dataset in a separate directory with its own driver and output paths.

Two source implementations are available: ratio-based C. reinhardtii inference and an essentiality-based E. coli benchmark. Decide which matches your measurements before translating them. Binary essentiality labels cannot be substituted for measured condition ratios without changing the workflow.

## Step 1 — Prepare the metabolic model

Import SBML with COBRA Toolbox or load a MAT file containing a COBRA model. Verify `S`, `rxns`, `lb`, `ub`, `b`, and a single biomass objective with `c == 1`; optimize biomass in each condition and inspect feasibility and growth. See the [complete model contract](input_formats.md#metabolic-models).

Record reaction IDs, gene IDs and GPRs. Gene/protein IDs used in the mutant table and benchmark must map to model genes. Models without existing GPRs can still be used in reaction optimization if structurally compatible, but those missing associations cannot serve as ground truth for recall. Keep exchange and biomass reactions explicitly identifiable.

The ratio-based implementation minimizes the signed sum of reaction variables. It is designed around the supplied models with split, nonnegative forward/backward variables; it does not convert a conventional model's signed reversible variables automatically. Its reversible helper requires backward IDs with `_REV` and the corresponding forward IDs in every active model. Other representations require a reviewed adaptation. The E. coli implementation handles signed reversible reactions but implements a different optimization procedure.

Extra constraints outside `S`, bounds and the documented fields are not automatically included in the joint problems. Check whether a proposed protein-constrained model encodes its protein capacity in the matrix/bounds. Existing code directly calls Gurobi; verify its MATLAB interface and license using the supplied example.

## Step 2 — Define experimental conditions

Construct a COBRA model for each represented medium/environment, with exchange bounds that reflect your experiment. Inspect uptake sign conventions: C. reinhardtii helper functions cap positive split uptake using upper bounds; E. coli sets conventional negative uptake lower bounds. Check each exchange reaction ID and retain an explicit table of condition names, model files, changed bounds and units.

The supplied `changeuptake`, `create_CO2_model`, and `create_hypo_model` functions encode C. reinhardtii nutrient IDs, 3%-of-range CO₂ adjustment and TAP fractions. They are examples of **case-study specific condition construction**, not general recipes for another organism. Derive new bounds from the new media/measurements; there is no automatic medium importer or inference of nutrient uptake from phenotype data here.

## Step 3 — Prepare mutant phenotype data

Retain mutant ID, affected gene ID, mapping confidence, and a numeric profile across measured comparisons for every strain/profile. In the supplied library, current QC compares different insertion-mutant profiles grouped by `Gene`; repeated screens are represented in comparison columns and averaged separately. Keep distinct mutant alleles, independent experimental repeats, and technical repeats identifiable in a new dataset. See the [source-paper and local-table evidence](input_formats.md#source-phenotype-spreadsheet-and-intermediate-files).

For the supplied ratio interface, prepare a per-gene row of **dimensionless log2 abundance/phenotype ratios**, convert to linear ratios with `2.^x`, and invert comparisons when their orientation differs from the model pair. Do not send log2 values directly to `GEM_PROSPECT`. Determine whether your measurements support growth-ratio constraints; that relationship and normalization need experimental justification for the new dataset.

The exact current slots are `[1/2, 1/3, 2/3, 1/4, 2/5, 2/6, 2/7, 2/8]`, with eight ordered model slots. A matching condition topology can reuse the functions after relabeling the models and measurements. An arbitrary topology, additional comparisons, or a different number of slots requires reviewing both ratio functions, their index allocations, and downstream tables. No exposed `ratio_pairs` argument exists today. Do not silently relabel a measured comparison to an unrelated pair.

Use [the header template](../examples/phenotype_template.csv) only for this existing eight-slot contract. For an essentiality dataset, use the [E. coli format](input_formats.md#e-coli-essentiality-inputs), while recognizing that its script assigns the same gene-essentiality list across all conditions.

## Step 4 — Quality control and filtering

The implemented C. reinhardtii QC removes 3′ UTR insertions, retains confidence levels ≤4, removes entirely missing profiles, and requires ≥3 insertion-mutant profiles per gene. It computes the mean pairwise Spearman correlation between these mutants across the full selected screen profile. For groups of ≥4 profiles, when the mean is <0.5 it removes the profile with greatest summed correlation distance and repeats until the threshold is met or three remain. The next R script retains only genes with mean correlation **>0.5** and drops NaN correlations. This is GEM-PROSPECT's consistency filter across mutant alleles; the source paper's screen-replicate experiments are a separate unit.

Definitions are in `Code/python/pairwise_correlation.py` and `Code/R/filter genes with consistent phynotype.R`. The **0.5** threshold selects phenotype consistency; neither an empirical search for it nor universal applicability is established in this repository. Decide a QC policy appropriate to the new library and document its evidence. The insertion confidence scale and 3′ UTR rule are specific to the supplied dataset.

The averaging script takes means of log2 measurements with `na.rm=TRUE`, first over named screens and then over gene-grouped library rows. The main driver then requires all eight ratios and filters the protein set using GO activity annotations while excluding `antiporter` entries. These are separate selection steps. The core functions can skip NaN ratio slots, but the manuscript driver does not use that flexibility. Do not impute absent measurements as zero (zero log2 means a ratio of one).

For the supplied library, input log2 abundance ratios were already processed by the original study. The MATLAB mapping and R aggregation scripts create the intermediate and final tables; see [the preprocessing sequence](input_formats.md#source-phenotype-spreadsheet-and-intermediate-files). For a new library, document its normalization and adapt the mapping/aggregation step to its identifiers and screen names. Validate substring matching against your identifier conventions rather than assuming the supplied mapping policy transfers unchanged.

## Step 5 — Configure GEM-PROSPECT parameters

| Parameter | Meaning | Current value | Used for | Dataset-specific? | Where defined |
|---|---|---|---|---|---|
| α | Minimum fraction of condition-specific optimal biomass | 0.10; manuscript comparison 0.01/0.05 | Ratio-based optimization | Yes: manuscript explicitly says calibration may be needed | `GEM_PROSPECT_implementation.m`, core functions' `alpha` argument |
| ω | Slack `abs(v_i - r*v_j) ≤ r*ω` | Fixed 0.05; manuscript comparison 0.01/0.10 | Growth-ratio relations | Yes: manuscript explicitly says calibration may be needed | `GEM_PROSPECT.m`, `GEM_PROSPECT_reversible_rxns.m`: `w = ratio_val * 0.05` |
| CO₂ adjustment | 3% of uptake flux difference between two growth-constrained solutions, added to upper reference | 0.03 | CO₂-enriched models | **C. reinhardtii case-study specific** | `create_CO2_model.m` |
| Flux-reduction threshold | log10 maximum/reference flux ratio | −1; recall sweep −1…−5 | Candidate/evaluation calls | **C. reinhardtii case-study specific**, empirically compared in manuscript | R association and recall scripts |
| τ | Absolute biomass cutoff for essential/nonessential reactions | 10⁻⁷…10⁻¹ (executed in ascending order) | Essentiality benchmark | **E. coli case-study specific** | `GEM_PROSPECT_Ecoli.m` |
| Mean Spearman threshold | Select consistent library profiles | Python trimming <0.5; R acceptance >0.5 | Phenotype QC | **C. reinhardtii preprocessing specific**; optimality elsewhere not established | Python correlation and R filtering scripts |

See [parameters](parameters.md) for every important setting, hard-coded IDs, tested script values, numerical tolerances, and rerun dependencies. α=0.10 and FVA percentage=10 are distinct settings; the FVA helper fixes growth exactly at 10% and does not impose a 10%-around-optimum interval. “One-fold flux reduction” is implemented as a strict log10 threshold of −1, i.e. below 0.1 of reference flux.

To calibrate new settings, use appropriate known associations, inspect feasibility and coverage, compare candidate/recall behavior across a recorded parameter grid, and keep an independent evaluation set when possible. This is a proposed research procedure, not an automatic calibration feature. Pass each chosen α explicitly to the core functions and record condition-construction settings separately. Decide explicitly whether new condition constraints should change too.

## Step 6 — Run the model

For the supplied dataset, the exact supported execution order is in [reproduction instructions](reproducing_manuscript.md). For a compatible new ratio dataset, your dataset driver must prepare `models`, `ratios`, `rxn_list_ir`, and backward `rxn_list_rev` according to the contract before calling:

```matlab
% These variables must already contain your prepared models, measurements,
% candidate IDs and explicitly selected parameters.
flux_ir = GEM_PROSPECT(models, ratios, rxn_list_ir, alpha);
flux_rev = GEM_PROSPECT_reversible_rxns(models, ratios, rxn_list_rev, alpha);
```

The core functions use ω=0.05 internally and do not expose an ω argument. Exploring another ω requires an explicit adaptation of separate copies of both core functions at the `w = ratio_val * 0.05` assignment, with recorded changes and validation.

Each returns a 1×N cell array, aligned with the input candidate list. Validate every candidate's presence and split reverse/forward pairing first. No mutation deletion or GPR update is performed by these calls: they test each candidate using the phenotype-derived growth constraints. Save mutant/gene IDs, reaction IDs, all settings and run provenance with your outputs. The existing `writeToFile` is a simple unquoted CSV appender and should not be used for identifiers containing commas or newlines.

For new E. coli-like essentiality data, adapt a separate copy of `GEM_PROSPECT_Ecoli.m`: update model, gene labels, condition/source IDs, candidate exclusion and output paths consistently with FVA. This script builds its own joint LP, not the ratio-based helper. An unlabeled new mutant's function cannot be inferred by simply running this labeled benchmark.

## Step 7 — Inspect and validate predictions

Check solver success, model feasibility, nonfinite outputs, identifier joins and reference coverage before thresholding. `NaN` is not evidence that a reaction is blocked. The existing reversible and E. coli scripts sometimes replace non-optimal second-stage values with zero, so solver-status-aware review is necessary before biological interpretation. The original outputs do not store per-solve status.

The case-study R scripts remove entirely missing profiles, then entirely missing reaction columns, while retaining the ID column and valid zero fluxes. Review the [missing-data policy](reproducing_manuscript.md#r-missing-data-handling) and compare profile IDs before and after cleaning when adapting postprocessing to a new dataset. Profiles with no valid reaction values do not reach candidate calling.

Rebuild reaction-index maps whenever the candidate set changes. Reference sampling must use matching reaction IDs and a documented growth constraint; the supplied workflow uses separate Auto/Hetero/Mixo reference means and intersects candidate proteins across their three comparisons. Zero/missing means make the ratio undefined. There is no confidence probability or validated ranking model.

Recall requires known positive gene–reaction associations and matched gene/reaction IDs. Accuracy additionally requires a defined negative class and denominator. Decide whether you evaluate gene–reaction pairs, reactions, or genes; the C. reinhardtii figures report different views. A GPR containing isoenzymes is different from an “any essential gene” mapping. Report evaluated/excluded counts, feasibility and missing-data policies, and assess novel predictions experimentally or against independent evidence.

The supplied E. coli benchmark uses essentiality labels in both construction and evaluation. Its accuracy is not a held-out cross-species performance claim. Only the two included case studies are documented here; adaptation to a new organism needs its own validation.
