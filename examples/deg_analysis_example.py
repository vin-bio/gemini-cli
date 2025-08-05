#!/usr/bin/env python3
"""
Example: Create test data and demonstrate DEG analysis tool usage

This script creates synthetic gene expression data for testing the DEG analysis tool
and shows how to use it through the Gemini CLI.
"""

import pandas as pd
import numpy as np
import os

def create_test_expression_data():
    """Create synthetic gene expression data for testing"""
    np.random.seed(42)  # For reproducible results
    
    # Create sample metadata
    n_control = 8
    n_treatment = 8
    total_samples = n_control + n_treatment
    
    # Sample names and conditions
    samples = [f"Sample_{i+1}" for i in range(total_samples)]
    conditions = ['control'] * n_control + ['treatment'] * n_treatment
    
    # Gene names (mix of real-looking gene symbols)
    genes = [
        'GAPDH', 'ACTB', 'TP53', 'BRCA1', 'EGFR', 'MYC', 'KRAS', 'PIK3CA',
        'PTEN', 'RB1', 'CDKN2A', 'APC', 'VHL', 'SMAD4', 'TGFBR2', 'MLH1',
        'MSH2', 'MSH6', 'PMS2', 'STK11', 'FHIT', 'WWOX', 'PARK2', 'MACROD2'
    ] + [f"Gene_{i+25}" for i in range(276)]  # Total 300 genes
    
    # Create expression matrix (genes x samples)
    expression_data = {}
    
    for gene_idx, gene in enumerate(genes):
        # Base expression level for this gene
        base_expression = np.random.lognormal(mean=5, sigma=1)
        
        # Control samples - normal variation around base
        control_values = np.random.lognormal(
            mean=np.log(base_expression), 
            sigma=0.3, 
            size=n_control
        )
        
        # Treatment samples - some genes are differentially expressed
        if gene_idx < 20:  # First 20 genes are differentially expressed
            if gene_idx < 10:  # First 10 are upregulated
                fold_change = np.random.uniform(2, 5)  # 2-5x upregulation
                treatment_mean = np.log(base_expression * fold_change)
            else:  # Next 10 are downregulated
                fold_change = np.random.uniform(0.2, 0.5)  # 2-5x downregulation
                treatment_mean = np.log(base_expression * fold_change)
        else:  # Remaining genes are not differentially expressed
            treatment_mean = np.log(base_expression)
        
        treatment_values = np.random.lognormal(
            mean=treatment_mean,
            sigma=0.3,
            size=n_treatment
        )
        
        # Combine all values for this gene
        gene_values = np.concatenate([control_values, treatment_values])
        expression_data[gene] = gene_values
    
    # Create DataFrame
    df = pd.DataFrame(expression_data, index=samples)
    
    # Add condition column
    df['condition'] = conditions
    
    # Reset index to make sample names a column
    df = df.reset_index().rename(columns={'index': 'sample_id'})
    
    # Reorder columns so condition is second
    cols = ['sample_id', 'condition'] + [col for col in df.columns if col not in ['sample_id', 'condition']]
    df = df[cols]
    
    return df

def main():
    """Generate test data and show usage example"""
    print("Creating synthetic gene expression data...")
    
    # Create test data
    expression_df = create_test_expression_data()
    
    # Save to CSV
    output_file = 'test_gene_expression.csv'
    expression_df.to_csv(output_file, index=False)
    
    print(f"✅ Created test dataset: {output_file}")
    print(f"   - Samples: {len(expression_df)} ({expression_df['condition'].value_counts()['control']} control, {expression_df['condition'].value_counts()['treatment']} treatment)")
    print(f"   - Genes: {len(expression_df.columns) - 2}")
    print(f"   - Expected significant genes: ~20 (10 up, 10 down)")
    
    # Show data preview
    print(f"\nData preview:")
    print(expression_df.head())
    
    print(f"\nCondition distribution:")
    print(expression_df['condition'].value_counts())
    
    print(f"\n" + "="*60)
    print("HOW TO USE THE DEG ANALYSIS TOOL")
    print("="*60)
    
    # Show example usage
    example_usage = f"""
# Basic usage with Gemini CLI:
npm start -- --prompt "Perform DEG analysis on the file {os.path.abspath(output_file)}. 
Compare 'control' vs 'treatment' groups using the 'condition' column. 
Use default parameters for statistical testing."

# Or with specific parameters:
npm start -- --prompt "Use the deg_analysis tool with these parameters:
- expression_file: {os.path.abspath(output_file)}
- condition_column: condition  
- group1: control
- group2: treatment
- log2fc_threshold: 1.5
- pvalue_threshold: 0.01
- correction_method: fdr
- statistical_test: ttest
- create_plots: true"

# What the tool will generate:
# - deg_results.csv: Complete results for all genes
# - significant_genes.csv: Only significant genes
# - volcano_plot.png: Volcano plot visualization  
# - ma_plot.png: MA plot visualization
# - deg_summary.txt: Analysis summary
# - deg_analysis_script.py: The analysis script used
"""
    
    print(example_usage)
    
    print("\n" + "="*60)
    print("PYTHON DEPENDENCIES")
    print("="*60)
    print("""
The DEG analysis tool requires these Python packages:
- pandas
- numpy  
- scipy
- matplotlib (optional, for plots)
- seaborn (optional, for plots)

Install with: pip install pandas numpy scipy matplotlib seaborn
""")

if __name__ == "__main__":
    main() 