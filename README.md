# GEM-PROSPECT

**GEM-PROSPECT** integrates condition-specific growth phenotypes from mutant libraries with protein-constrained genome-scale metabolic models (pcGEMs) to predict protein functions. It links unannotated (orphan) metabolic reactions to the corresponding knocked-out genes and can also refine existing gene-protein-reaction (GPR) associations.

### Main capabilities
1. **Refinement & validation of GPR rules** — using known gene-reaction associations from mutant phenotype data.  
2. **Prediction of candidate genes for orphan reactions** — distinguishing transporters and enzymes, with in silico validation support (SPOT for transporters; CataPro for enzyme kinetics).

The framework was validated on:
- 884 *Chlamydomonas reinhardtii* mutants across eight growth conditions (photoautotrophic, mixotrophic, heterotrophic, etc.) using a recently published pcGEM (https://www.nature.com/articles/s41467-023-40498-1).
- Gene essentiality data from *Escherichia coli* across 16 carbon sources using iML1515 models (maximum accuracy 81.2% when accounting for isoenzymes).

---

### Installation

**Dependencies**  
- MATLAB (tested with R2023b / 23.2.0)  
- COBRA Toolbox v3.0 (or later)  
- Gurobi Optimizer (tested with v11+; academic license required)  
- (Optional) R (for statistical analyses and figure generation)  
- (Optional) Python (for external deep-learning validation tools: CataPro)

**How to install**  
1. Clone the repository:  
   ```
   git clone https://github.com/YunliEricHsieh/GEM-PROSPECT.git
   cd GEM-PROSPECT
2. Add the repository folder (and subfolders) to your MATLAB path.
3. Initialize the COBRA Toolbox and set Gurobi as the solver
   ```
   initCobraToolbox(false);
   changeCobraSolver('gurobi','LP');
   ```

**Repository Structure**

The repository is organized as follows:
```
GEM-PROSPECT/
├── Code/          # Core MATLAB functions and main analysis scripts
├── Data/          # Curated & filtered phenotype dataset (884 mutants), model files, and input data
├── Results/       # Reproduction scripts, intermediate results, and figures from the manuscript
└── README.md
```

### Reproducing the Paper Results
All scripts are in Code/. Run them in MATLAB:
```matlab
% run GEM-PROSPECT for C. reinhardtii
GEM_PROSPECT_implementation.m

% run GEM-PROSPECT for E. coli
GEM_PROSPECT_Ecoli.m
```

### Citation
TODO
