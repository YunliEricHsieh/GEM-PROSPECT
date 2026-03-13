import os
import glob
import pandas as pd
import cobra
import numpy as np
import matplotlib.pyplot as plt
import re

def load_data_and_mappinf(model_path, essentiality_path):
    print('Loading empirical gene essentiality data...')
    df_genes = pd.read_csv(essentiality_path,  sep=';')
    essential_genes = set(df_genes.loc[df_genes['Essentiality_0_1'] == 1, 'Gene_ID'].astype(str))

    print(f'Loading metabolic model from {model_path}...')
    model = cobra.io.read_sbml_model(model_path)

    naive_assoc = {}
    gpr_assoc = {}

    print('Building Navie and GPR-Aware mappings...')

    # builf naive mapping
    for rxn in model.reactions:
        if not rxn.genes:
            continue  # skip reactions without gene associations
        naive_assoc[rxn.id] = any(gene.id in essential_genes for gene in rxn.genes)

    # build GPR-aware mapping
    with model as temp_model:
        for rxn in temp_model.reactions:
            if rxn.genes:
                rxn.bounds = (-1000.0, 1000.0)
        
        for gene_id in essential_genes:
            try:
                temp_model.genes.get_by_id(gene_id).knock_out()
            except KeyError:
                pass

        for rxn in temp_model.reactions:
            if not rxn.genes:
                continue # skip reactions without gene associations
            is_strictly_dependent = (rxn.lower_bound == 0.0 and rxn.upper_bound == 0.0)
            gpr_assoc[rxn.id] = is_strictly_dependent

    return naive_assoc, gpr_assoc

def parse_files(directory):
    files = glob.glob(os.path.join(directory, '*.csv'))
    if not files:
        print(f'WARNING: No CSV files found in {directory}')
        return []
    
    parsed_data = []
    for file_path in files:
        try:
            df = pd.read_csv(file_path, sep='\t')
            if len(df.columns) == 1:
                df = pd.read_csv(file_path, sep=',')
            parsed_data.append((os.path.basename(file_path), df))
        except Exception as e:
            print(f'Error reading {file_path}: {e}')
    return parsed_data

def calculate_accuracy_metrics(file_path, target_column, assoc_dict, TOLERANCE):
    results = []

    for filename, df in file_path:
        if target_column not in df.columns:
            print(f'WARNING: Column {target_column} not found in {filename}. Skipping this file.')
            continue

        tp = tn = fp = fn = 0
        for _, row in df.iterrows():
            rxn_id = str(row['RxnID'])
            if rxn_id not in assoc_dict:
                continue

            flux = float(row.get(target_column, 0.0))
            if pd.isna(flux):
                flux = 0.0

            is_pred = abs(flux) <= TOLERANCE
            is_actual = assoc_dict[rxn_id]

            if is_actual and is_pred: tp += 1
            elif not is_actual and not is_pred: tn += 1
            elif not is_actual and is_pred: fp += 1
            elif is_actual and not is_pred: fn += 1

        total = tp + tn + fp + fn
        accuracy = (tp + tn) / total if total > 0 else 0.0
    
        # extract scientific natation string for label
        match = re.search(r'(\d+\.?\d*e[+-]?\d+)', filename, re.IGNORECASE)
        threshold_label = match.group(1) if match else 'Unknown'
        threshold_val = float(threshold_label) if match else 0.0

        results.append({
            'Threshold_Val': threshold_val,
            'Threshold_Label': threshold_label,
            'Accuracy': accuracy
        })

    # sort from largest to smallest threshold
    return pd.DataFrame(results).sort_values(by='Threshold_Val', ascending=False)

def plot_dual_accuracy_figure(custom_files, naive_assoc, gpr_assoc, output_filename, tolerance):
    print('\nGenerating combined accuracy plot...')

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(9, 6))

    # process and plot naive approach
    df_naive = calculate_accuracy_metrics(custom_files, 'MaxFlux', naive_assoc, TOLERANCE=tolerance)
    if not df_naive.empty:
        max_acc_naive = df_naive['Accuracy'].max()
        ax.plot(df_naive['Threshold_Label'], df_naive['Accuracy'],
                color='darkorange', lw = 2.5, marker='o', markersize=6,
                label=f'Without considering isoenzymes (Max Acc = {max_acc_naive:.3f})')
        
    # process and plot GPR-aware approach
    df_gpr = calculate_accuracy_metrics(custom_files, 'MaxFlux', gpr_assoc, TOLERANCE=tolerance)
    if not df_gpr.empty:
        max_acc_gpr = df_gpr['Accuracy'].max()
        ax.plot(df_gpr['Threshold_Label'], df_gpr['Accuracy'],
                color='navy', lw = 2.5, marker='s', markersize=6,
                label=f'Considering isoenzymes (Max Acc = {max_acc_gpr:.3f})')
        
    # Formatting
    ax.set_xlabel('Simulation Tau', fontsize = 14, fontweight='bold')
    ax.set_ylabel('Predictive Accuracy', fontsize = 14, fontweight='bold')
    
    ax.tick_params(axis='both', labelsize=12)
    plt.xticks(rotation=45)
    ax.legend(loc='lower right', frameon=True, fontsize=12)
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f'Combined accuracy plot saved as {output_filename}')
    plt.show()

def run_full_analysis(model_path, essentiality_path, custom_dir, tolerance=1e-5):
    naive_assoc, gpr_assoc = load_data_and_mappinf(model_path, essentiality_path)

    print('\nReading result files...')
    custom_files = parse_files(custom_dir)

    if not custom_files:
        print('No valid files to process. Exiting.')
        return
    
    plot_dual_accuracy_figure(
        custom_files=custom_files, 
        naive_assoc=naive_assoc, 
        gpr_assoc=gpr_assoc, 
        output_filename='Results/figures/GEM_PERSPECT_accuracy_comparison.png', 
        tolerance=tolerance)
    
    print('\nPipeline execution complete!')

if __name__ == "__main__":
    run_full_analysis(
        model_path='Data/Ecoli/iML1515.xml',
        essentiality_path='Data/Ecoli/Ecoli_gene_essentiality.csv',
        custom_dir='Results/screens/Ecoli',
        tolerance=1e-7
    )
