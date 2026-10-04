# Input formats and identifiers

Paths in this document describe the existing implementation. There is no general CSV importer for arbitrary phenotype layouts. All commands run from the repository root unless a wrapper explicitly resolves it.

## Metabolic models

C. reinhardtii scripts load the three `Data/pciCre1355/NDLadpraw_*_Rep1.xml` SBML files with `readCbModel`. E. coli modeling loads `Data/Ecoli/iML1515.mat`, which must define `model`; Python evaluation reads `Data/Ecoli/iML1515.xml` with COBRApy. Confirm that the MAT and XML versions describe the same reactions and GPRs before replacing either file.

The supplied E. coli MAT model stores COBRA-indexed `rules`. Its inference script constructs `grRules` with `creategrRulesField` and `rxnGeneMat` with `buildRxnGeneMat` before mapping essential genes to reactions.

| COBRA field | Required meaning |
|---|---|
| `S` | Metabolites × reactions stoichiometric matrix |
| `lb`, `ub` | Reaction bounds in model flux units, aligned with columns of `S` |
| `rxns` | Unique reaction identifiers, aligned with `S` and bounds |
| `c` | Biomass objective; use exactly one entry equal to 1 |
| `b` | Mass-balance right-hand side; FVA and individual E. coli controls read it |
| `genes`, `rules` | Gene identifiers and COBRA-indexed GPR rules; reaction metadata uses these |
| `rxnGeneMat`, `grRules` | Used/built in E. coli essentiality mapping |
| `mets`, `metNames`, `metFormulas`, `subSystems` | Additional fields for optional substrate extraction |

The C. reinhardtii routines and the E. coli joint optimization construct `S*v = 0` themselves; they do not propagate arbitrary nonzero `b` or additional custom constraint fields. Model transformations that require other constraints need explicit review. Condition models need a common identifier set for the tested reactions; the C. reinhardtii functions map IDs separately per model, while the E. coli script assumes identical reaction ordering across its copied models.

The supplied SBML files name the flux unit `mmol_per_gDW_per_hr`. Metabolic flux is conventionally mmol gDW⁻¹ h⁻¹ and biomass flux represents a specific growth rate (h⁻¹). Protein pseudo-reaction scaling is inherited from the supplied models. No unit conversion is performed by GEM-PROSPECT.

C. reinhardtii reversible reactions are represented by nonnegative forward and backward reactions with backward IDs containing `_REV`. `GEM_PROSPECT_reversible_rxns` receives **backward** IDs, removes `_REV` to obtain the forward ID, links both directions across conditions, and optimizes `v_forward - v_backward`. A model with ordinary signed reversible variables cannot be passed unchanged to this routine. The irreversible routine minimizes the signed sum of all variables, not an automatically constructed sum of absolute values.

COBRA imports may strip SBML prefixes such as `R_`; always use identifiers in the imported `model.rxns`, not guessed raw XML identifiers. Reaction names, generated `RxnIndex`, and model `RxnID` are different objects.

## Processed C. reinhardtii phenotypes

`Data/Mutant_phenotypes_table_filtered_final.csv` has one row per `GeneID`. `UniProtID` identifies the protein(s); multiple IDs use the literal separator ` or `. The file contains 1,832 rows in the inspected checkout. The main script selects activity-annotated proteins, excludes GO rows containing `antiporter`, and requires all eight used ratios to be non-NaN; this gives 884 rows with the current supplied inputs.

All eight CSV values below are **dimensionless log2 phenotype ratios**. The optimization functions instead take **linear** ratios in the explicit order below. They do not read the CSV themselves.

| Ratio slot | CSV column | Conversion in the main script | Model pair / linear ratio |
|---|---|---|---|
| 1 | `auto_mixo` | `2^x` | 1 / 2, Auto / Mixo |
| 2 | `auto_hetero` | `2^x` | 1 / 3, Auto / Hetero |
| 3 | `mixo_hetero` | `2^x` | 2 / 3, Mixo / Hetero |
| 4 | `auto_auto_CO2` | `2^x` | 1 / 4, Auto / Auto_CO2 |
| 5 | `mixo_mixo_CO2` | `2^x` | 2 / 5, Mixo / Mixo_CO2 |
| 6 | `mixo_hypo10_mixo` | `1/(2^x)` | 2 / 6, Mixo / Hypo10 |
| 7 | `mixo_hypo25_mixo` | `1/(2^x)` | 2 / 7, Mixo / Hypo25 |
| 8 | `mixo_hypo75_mixo` | `1/(2^x)` | 2 / 8, Mixo / Hypo75 |

The model cell array is `{Auto; Mixo; Hetero; Auto_CO2; Mixo_CO2; Hypo10; Hypo25; Hypo75}`. Hypoosmotic CSV column names refer to the opposite orientation from the optimizer's pairs, hence the reciprocal. The additional CSV columns `mixo_NaCl_mixo`, `mixo_P_mixo`, `mixo_N_mixo` are aggregated but not used by the main optimization.

For example, an illustrative log2 value of 1 in `auto_mixo` becomes a linear ratio of 2. The same illustrative value in `mixo_hypo10_mixo` becomes a ratio of 0.5. These arithmetic examples are not biological measurements. [The header-only template](../examples/phenotype_template.csv) uses the exact processed-table names and contains no invented results.

Blank/`NA` fields must import as numeric `NaN`. The main script removes rows missing any of the eight ratios. The public functions skip NaN ratio slots and activate only models belonging to the remaining pairs; they do not impute missing ratios. At least one pair must be active. Infinite values and all-NaN vectors are not validated by the original functions; check finite, positive linear ratios before calling them.

`Data/Mutant_phenotypes_table.xlsx` contains the original study's already processed log2 abundance ratios. This repository starts with those processed values and performs profile filtering, screen/gene averaging, and conversion to linear ratios for the biomass constraints. Normalization from raw counts is described in the [original study](https://www.nature.com/articles/s41588-022-01052-9).

## Source phenotype spreadsheet and intermediate files

`Data/Mutant_phenotypes_table.xlsx` supplies already processed log2 abundance ratios from the original study and is read by `pairwise_correlation.py` with the first row as the header. There are five metadata columns, in this order:

| Position | Exact header | Meaning |
|---|---|---|
| 1 | `Mutant ID` | Library strain/insertion identifier |
| 2 | `Gene` | C. reinhardtii gene identifier, e.g. `Cre01.g000250` |
| 3 | `Feature` | Insertion feature; exact value `3'UTR` is excluded |
| 4 | `Confidence level` | Numeric annotation confidence; retain ≤4 |
| 5 | `Gene Name` | Descriptive annotation |
| 6 onward | Existing `c_R…` screen names | Original study's processed log2 abundance ratios across screens |

The script treats **all** columns from position 6 onward as phenotype columns. Grouping is by `Gene`; the correlated rows represent different insertion-mutant profiles assigned to that gene. A local audit of the supplied `filtered_by_num` table found no repeated `Gene`–`Mutant ID` pairs. For example, `Cre01.g000250` has three rows with mutant IDs `LMJ.RY0402.153768`, `LMJ.RY0402.215015`, and `LMJ.RY0402.070496`, with CDS/intron insertion annotations.

[Fauser et al. (2022), High-confidence gene–phenotype relationships and Methods](https://www.nature.com/articles/s41588-022-01052-9) distinguish multiple independent mutant alleles of a gene from independent replicates of a screen. Here, the ≥3-profile rule and Spearman comparison concern mutant profiles grouped by gene; repeated-screen measurements occupy separate comparison columns and are averaged separately by the R aggregation script. Calling the ≥3-row rule “three biological replicates of a mutant” would conflate these two units. The Spearman 0.5 filtering procedure is implemented in GEM-PROSPECT and should not be attributed to the original paper's gene–phenotype statistical test.

Missing pairs are handled by pandas Spearman correlation; mean correlations use `np.nanmean`. The script removes profiles missing all screen values, keeps genes with ≥3 mutant profiles, and for groups of ≥4 removes the profile with the largest sum of `1 - correlation` while mean correlation is <0.5. Three-profile groups are scored without further trimming. All-undefined correlations can remain NaN and are removed by the next R filter. The R filter finally retains **strictly >0.5**.

Outputs in sequence:

1. `Data/Mutant_phenotypes_table_filtered_by_num.csv`: retained library rows, original headers.
2. `Results/mean_correlation_of_phynotype/mean_correlations.csv`: `Gene`, `Mean Correlation`, dimensionless.
3. `Data/Mutant_phenotypes_table_filtered_by_num_and_cor.csv`: rows passing the R correlation filter. Base R converts punctuation in column names to dots.
4. `Data/Mutant_phenotypes_table_filtered_with_uniportIDs.csv`: MATLAB mapping output with appended `UniProtID`; retains the existing filename spelling.
5. `Data/Mutant_phenotypes_table_filtered_final.csv`: R averages named screen replicates and then mutant profiles per gene, using `na.rm=TRUE`; averages log2 values, not linear values. All-missing groups yield NaN. Gene IDs containing `&` are removed.

Both mapping and aggregation are implemented in the repository. `UniProtIDs_mapping.m` reads the source tables, matches gene IDs by substring, joins matching protein IDs with ` or `, removes duplicate protein IDs, and writes the mapped CSV. The subsequent R script selects the named screens and calculates the per-gene means. MATLAB's [default text-table import naming rule](https://www.mathworks.com/help/matlab/ref/readtable.html) modifies headers to valid identifiers; its [identifier conversion rules](https://www.mathworks.com/help/matlab/ref/matlab.lang.makevalidname.html) explain `Gene Name` → `GeneName`, `UniProt ID` → `UniProtID`, and punctuation → underscores. The supplied screen headers also reflect truncation to the identifier-length limit used when they were imported.

The supplied mapped table has 7,190 rows, including six with no matching protein ID and blank `UniProtID`. Profiles within each gene share the same mapped ID string. When replacing the input tables, verify the imported headers, mapping matches, and screen-column names used in `calculate mean of phynotypes from each screens.R` before aggregation. The final phenotype table supplied with the repository is the starting point for the quick start and main inference workflow.

## Annotation and reaction tables

| File | Fields and usage |
|---|---|
| `Data/mart_Cre_Uniprot.txt` | Tab-separated `Gene Name`, `UniProt ID`, `Peptide Name`, `ENZYME ID`; maps gene IDs to protein IDs by substring matching in the existing MATLAB script |
| `Data/GO_table.txt`, `Data/GO_table_filtered.txt` | Tab-separated `UniProtID`, `GO_ID`, `GO_Name`, `GO_Aspect`; filtered file includes GO names containing `activity` from QuickGO retrieval |
| `Data/Reactions/list_of_irreversible_rxns.csv` | `RxnIndex`, `RxnID` |
| `Data/Reactions/list_of_reversible_rxns.csv` | `RxnIndex`, `RxnID`, for both forward and backward directions |
| `Data/Reactions/list_of_rxns_without_GPR.csv` | `RxnIndex`, `RxnID`; orphan-candidate subset |
| `Data/Reactions/list_of_rxns_with_proteins_in_both.csv` | `RxnIndex`, `RxnID`, semicolon-separated `Enzymes`; all associated proteins present in the intersected model/library set |
| `Data/Reactions/list_of_rxns_with_at_least_one_protein_in_both.csv` | Same schema; at least one associated protein is present; used for recall evaluation |

`RxnIndex` is generated from the sorted common reaction set after FVA filtering, **before** the `No`/`DM_` exclusion. Index gaps are expected. Changing models, FVA, or candidate filters requires regenerating metadata, flux columns, and sampling tables together. Do not join regenerated results to older index maps without verifying `RxnID` matches.

`Results/FVA/*.csv` and `Results/FVA/Ecoli/*.csv` contain `RxnID`, `minFlux`, `maxFlux`. The helper fixes biomass exactly at the specified percentage of optimal growth. Missing solver `objval` is written as `Inf` regardless of actual status; this is not a reliable proof of unboundedness. Reaction selection tests exact zero minima **and** maxima. C. reinhardtii excludes reactions blocked in **any** of its three FVA tables; E. coli excludes only reactions blocked in **all** 16 conditions.

Sampling CSVs contain `RxnIndex`, `RxnID`, `meanFlux`, `stdFlux`. For reversible pairs, means/std are calculated from `abs(v_forward - v_backward)` samples. For irreversible reactions they are calculated from raw reaction flux. A missing or zero reference mean makes the downstream ratio undefined/infinite; no imputation or pseudocount is implemented.

## E. coli essentiality inputs

The essentiality CSV is **semicolon separated**, even though its extension is `.csv`. Example from the supplied file:

| Gene_ID | Gene | Essentiality_0_1 |
|---|---|---|
| b2836 | aas | 0 |

`Gene_ID` must match imported model gene IDs. `Essentiality_0_1` is binary and dimensionless; 1 denotes essential. `Uniprot_ID`, `Gene`, `Function`, and `Essentiality` are supplied descriptive fields, not required by the main inference script. Missing gene IDs are not imputed; IDs absent from the model do not contribute to its mappings.

The supplied script uses one essential-gene list across every carbon source; it has no table of measured per-condition ratios. A reaction is assigned the essential-growth bound if **any** associated gene is essential, without evaluating the Boolean GPR in that assignment. The Python evaluator provides both this mapping and a GPR-aware mapping using collective essential-gene knockouts. An absent essentiality label is effectively outside the essential-gene set; this is not an explicit unknown class.

The 16 carbon source names/IDs are listed in [parameters](parameters.md). Original oxygen and other exchange constraints remain from iML1515; only the listed carbon-source lower bounds are changed. Do not infer that all unlisted nutrients were closed.
