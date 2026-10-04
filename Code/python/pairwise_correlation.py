"""C. reinhardtii library-profile QC; run from the repository root.

Reads the supplied XLSX (first-row header; five metadata columns), groups
profiles by Gene, and calculates iterative pairwise Spearman consistency.
Writes retained profiles and a mean-correlation table; the R filter later
requires strictly >0.5. This script does not reconstruct raw normalization.
"""

import numpy as np
import pandas as pd

# Specify the file path
file_path = "Data/Mutant_phenotypes_table.xlsx"

# Read the supplied Excel file using its first row as the header.
df = pd.read_excel(file_path)

# Filter the data from a specific column (remove 3'UTR and Confidence level > 4)
filtered_data = df[(df['Feature'] != "3'UTR") & (df['Confidence level'] <= 4)]
filtered_data = filtered_data.dropna(subset=filtered_data.columns[5:], how='all')

# filtered_data.to_csv('Data/Mutant_phenotypes_table_raw.csv', index=False)

# Group the data by the gene name and filter out the genes that have less than 3 data points
print(filtered_data.columns)
grouped = filtered_data.groupby('Gene').filter(lambda x: len(x) >= 3).groupby('Gene')
filtered_data = filtered_data.groupby('Gene').filter(lambda x: len(x) >= 3)

# Calculate the pairwise correlation within the genes for each mutant
mean_correlations = {}
for gene, group in grouped:
    numeric_data = group.iloc[:, 5:].transpose()

    # If the number of columns is 3, calculate the correlation and store the mean correlation coefficient
    if len(numeric_data.columns) == 3:
        corr_matrix = numeric_data.corr(method='spearman')
        off_diagonal = corr_matrix.values[np.triu_indices(n=corr_matrix.shape[0], k=1)]
        print(gene, np.nanmean(off_diagonal))
        mean_correlations[gene] = np.nanmean(off_diagonal)

    else: # If the number of columns is greater than 3
        while len(numeric_data.columns) >= 4:
            corr_matrix = numeric_data.corr(method='spearman')
            off_diagonal = corr_matrix.values[np.triu_indices(n=corr_matrix.shape[0], k=1)]

            # If the mean of the correlation coefficient is less than 0.5, remove the column with the maximum sum distance
            if np.nanmean(off_diagonal) < 0.5:
                new_matrix = 1-corr_matrix
                column_sum = new_matrix.sum()
                max_sum_column = column_sum.idxmax()
                numeric_data = numeric_data.drop(max_sum_column, axis=1)
                filtered_data = filtered_data.drop(max_sum_column, axis=0)
            else:
                print(gene, np.nanmean(off_diagonal))
                mean_correlations[gene] = np.nanmean(off_diagonal)
                break
        else:
            # If the number of columns is 3, calculate the correlation and store the mean correlation coefficient
            corr_matrix = numeric_data.corr(method='spearman')
            off_diagonal = corr_matrix.values[np.triu_indices(n=corr_matrix.shape[0], k=1)]
            print(gene, np.nanmean(off_diagonal))
            mean_correlations[gene] = np.nanmean(off_diagonal)

mean_correlations_df = pd.DataFrame(mean_correlations.items(), columns=['Gene', 'Mean Correlation'])
mean_correlations_df.to_csv('Results/mean_correlation_of_phynotype/mean_correlations.csv', index=False)
filtered_data.to_csv('Data/Mutant_phenotypes_table_filtered_by_num.csv', index=False)
