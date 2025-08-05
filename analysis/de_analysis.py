
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress future warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Use a publication-ready style for plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_context("paper", font_scale=1.5)

# --- 1. Load and Transpose Data ---
try:
    # Load data and transpose so that genes are rows and samples are columns
    control_df = pd.read_csv('/Users/vinayak/Documents/gemini-cli/data/BP1_1M (2).csv', index_col=0).T
    treatment_df = pd.read_csv('/Users/vinayak/Documents/gemini-cli/data/BP1_2M (1).csv', index_col=0).T
except FileNotFoundError as e:
    print(f"Error loading files: {e}")
    exit()

# --- 2. DE Analysis ---
# Find common genes between the two datasets
common_genes = control_df.index.intersection(treatment_df.index)
control_df_aligned = control_df.loc[common_genes]
treatment_df_aligned = treatment_df.loc[common_genes]

# Calculate means for each gene across samples
control_means = control_df_aligned.mean(axis=1)
treatment_means = treatment_df_aligned.mean(axis=1)

# Create a dataframe for analysis results
data = pd.DataFrame({'control': control_means, 'treatment': treatment_means})
data = data.dropna()

# Add a small constant to avoid division by zero
data += 0.01

# Calculate Log2 Fold Change
data['log2FC'] = np.log2(data['treatment'] / data['control'])

# Perform t-test to get p-values
# The t-test is performed on the aligned dataframes with multiple samples
p_values = ttest_ind(control_df_aligned.loc[data.index], treatment_df_aligned.loc[data.index], axis=1).pvalue
data['p_value'] = p_values

# Calculate -log10(p-value)
data['-log10_p_value'] = -np.log10(data['p_value'])

# Determine significance
p_value_cutoff = 0.1  # Less stringent p-value cutoff
data['significant'] = (data['p_value'] < p_value_cutoff) & (abs(data['log2FC']) > 1)
data['regulation'] = 'Not Significant'
data.loc[(data['significant']) & (data['log2FC'] > 1), 'regulation'] = 'Up-regulated'
data.loc[(data['significant']) & (data['log2FC'] < -1), 'regulation'] = 'Down-regulated'

# --- 3. Generate Multi-panel Figure ---
# Use a publication-ready style for plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_context("paper", font_scale=1.5)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
fig.suptitle('Differential Expression Analysis', fontsize=20, fontweight='bold')

# Panel A: Volcano Plot
sns.scatterplot(data=data, x='log2FC', y='-log10_p_value', hue='regulation',
                palette={'Up-regulated': '#E41A1C', 'Down-regulated': '#377EB8', 'Not Significant': '#999999'},
                ax=ax1, alpha=0.7, s=50, edgecolor='w', linewidth=0.5)
ax1.set_title('Volcano Plot', fontsize=16)
ax1.set_xlabel('Log2 Fold Change', fontsize=14)
ax1.set_ylabel('-log10(p-value)', fontsize=14)
ax1.axhline(y=-np.log10(p_value_cutoff), color='k', linestyle='--', linewidth=1)
ax1.axvline(x=1, color='k', linestyle='--', linewidth=1)
ax1.axvline(x=-1, color='k', linestyle='--', linewidth=1)
ax1.text(-0.1, 1.05, 'A', transform=ax1.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
ax1.legend(title='Regulation', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Panel B: Heatmap of top 50 significant genes
significant_genes = data[data['significant']].sort_values(by='p_value').head(50)

if not significant_genes.empty:
    # Combine the data for the heatmap from the original (transposed) dataframes
    heatmap_data = pd.concat([control_df_aligned.loc[significant_genes.index], 
                              treatment_df_aligned.loc[significant_genes.index]], axis=1)

    # Manually calculate z-score to ensure correct broadcasting and handling of zero variance
    mean = heatmap_data.mean(axis=1)
    std = heatmap_data.std(axis=1)
    heatmap_data_zscore = heatmap_data.subtract(mean, axis=0).divide(std, axis=0)
    # Replace any NaNs that result from zero variance (std dev = 0) with 0
    heatmap_data_zscore.fillna(0, inplace=True)

    sns.heatmap(heatmap_data_zscore, ax=ax2, cmap='vlag',
                yticklabels=False, xticklabels=True, cbar_kws={'label': 'Z-score'})
    ax2.set_title('Top 50 DE Genes', fontsize=16)
    ax2.set_xlabel('Samples', fontsize=14)
    ax2.set_ylabel('Genes', fontsize=14)
    ax2.text(-0.1, 1.05, 'B', transform=ax2.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
else:
    ax2.text(0.5, 0.5, 'No significant genes found', ha='center', va='center', fontsize=14)
    ax2.set_title('Top 50 DE Genes', fontsize=16)
    ax2.set_xlabel('Samples', fontsize=14)
    ax2.set_ylabel('Genes', fontsize=14)
    ax2.text(-0.1, 1.05, 'B', transform=ax2.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('de_analysis_figure.png', dpi=300)

print("DE analysis complete. Figure saved as de_analysis_figure.png")
print(f"{len(significant_genes)} significant genes found with a p-value cutoff of {p_value_cutoff}")
