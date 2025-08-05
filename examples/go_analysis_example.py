#!/usr/bin/env python3
"""
Example: Create test data and demonstrate GO analysis tool usage

This script creates test gene lists for GO analysis and shows how to use 
the GO analysis tool through the Gemini CLI.
"""

import pandas as pd
import numpy as np
import os

def create_test_gene_lists():
    """Create test gene lists for GO analysis"""
    
    # Example gene lists for different biological processes
    gene_lists = {
        'cell_cycle_genes': [
            'CCNA1', 'CCNA2', 'CCNB1', 'CCNB2', 'CCND1', 'CCND2', 'CCND3',
            'CCNE1', 'CCNE2', 'CDK1', 'CDK2', 'CDK4', 'CDK6', 'CDKN1A',
            'CDKN1B', 'CDKN2A', 'CDKN2B', 'E2F1', 'E2F3', 'RB1', 'TP53',
            'MDM2', 'ATM', 'CHEK1', 'CHEK2', 'WEE1', 'CDC25A', 'CDC25B'
        ],
        
        'apoptosis_genes': [
            'BAX', 'BAK1', 'BCL2', 'BCL2L1', 'BID', 'BAD', 'BIM', 'PUMA',
            'NOXA', 'CASP3', 'CASP7', 'CASP8', 'CASP9', 'APAF1', 'CYCS',
            'TP53', 'MDM2', 'P21', 'FAS', 'FASLG', 'TNF', 'TNFRSF1A'
        ],
        
        'dna_repair_genes': [
            'BRCA1', 'BRCA2', 'ATM', 'ATR', 'CHEK1', 'CHEK2', 'TP53',
            'XRCC1', 'XRCC2', 'XRCC3', 'XRCC4', 'PARP1', 'PARP2',
            'RAD51', 'RAD52', 'RAD18', 'MSH2', 'MSH6', 'MLH1', 'PMS2'
        ],
        
        'immune_response_genes': [
            'IL1B', 'IL6', 'IL10', 'TNF', 'IFNG', 'IL2', 'IL4', 'IL12A',
            'TLR4', 'TLR2', 'MYD88', 'NFKB1', 'RELA', 'IRF3', 'IRF7',
            'STAT1', 'STAT3', 'JAK1', 'JAK2', 'CD4', 'CD8A', 'CD14'
        ],
        
        'metabolism_genes': [
            'GAPDH', 'ALDOA', 'ENO1', 'PKM', 'LDHA', 'LDHB', 'G6PD',
            'PFKL', 'PFKM', 'PFKP', 'HK1', 'HK2', 'GPI', 'ALDOC',
            'TPI1', 'PGAM1', 'PGK1', 'ACLY', 'FASN', 'SCD'
        ]
    }
    
    # Create individual gene list files
    for list_name, genes in gene_lists.items():
        with open(f'{list_name}.txt', 'w') as f:
            f.write('\n'.join(genes))
        print(f"✅ Created {list_name}.txt ({len(genes)} genes)")
    
    # Create a CSV file with multiple gene lists (simulating DEG results)
    all_genes = []
    for list_name, genes in gene_lists.items():
        for gene in genes:
            all_genes.append({
                'gene': gene,
                'log2_fold_change': np.random.normal(0, 2),
                'pvalue': np.random.beta(0.5, 5),  # Most p-values small
                'adjusted_pvalue': np.random.beta(0.5, 5),
                'category': list_name.replace('_genes', ''),
                'significant': np.random.choice([True, False], p=[0.7, 0.3])
            })
    
    # Add some random genes
    for i in range(100):
        all_genes.append({
            'gene': f'GENE_{i+1}',
            'log2_fold_change': np.random.normal(0, 1),
            'pvalue': np.random.uniform(0.01, 0.5),
            'adjusted_pvalue': np.random.uniform(0.01, 0.5),
            'category': 'random',
            'significant': False
        })
    
    deg_results = pd.DataFrame(all_genes)
    deg_results.to_csv('simulated_deg_results.csv', index=False)
    print(f"✅ Created simulated_deg_results.csv ({len(deg_results)} genes)")
    
    # Create significant genes only file
    significant_genes = deg_results[deg_results['significant'] == True]['gene'].tolist()
    with open('significant_genes_only.txt', 'w') as f:
        f.write('\n'.join(significant_genes))
    print(f"✅ Created significant_genes_only.txt ({len(significant_genes)} genes)")
    
    return gene_lists

def create_background_gene_set():
    """Create a background gene set file"""
    # Common human gene symbols (simplified example)
    background_genes = [
        'ACTB', 'GAPDH', 'TUBB', 'ALB', 'INS', 'TP53', 'MYC', 'EGFR',
        'VEGFA', 'IL6', 'TNF', 'APOE', 'LDLR', 'PCSK9', 'BRCA1', 'BRCA2'
    ]
    
    # Add more random genes
    background_genes.extend([f'GENE_{i}' for i in range(1000, 5000)])
    
    with open('background_genes.txt', 'w') as f:
        f.write('\n'.join(background_genes))
    
    print(f"✅ Created background_genes.txt ({len(background_genes)} genes)")
    return background_genes

def main():
    """Generate test data and show usage examples"""
    print("Creating test gene lists for GO analysis...")
    
    # Create test gene lists
    gene_lists = create_test_gene_lists()
    background_genes = create_background_gene_set()
    
    print(f"\n" + "="*60)
    print("HOW TO USE THE GO ANALYSIS TOOL")
    print("="*60)
    
    # Show basic usage examples
    basic_examples = f"""
# Basic usage with simple gene list:
npm start -- --prompt "Perform GO enrichment analysis on the file {os.path.abspath('cell_cycle_genes.txt')}. 
Use human organism and default parameters."

# Analysis with specific parameters:
npm start -- --prompt "Use go_analysis tool with these parameters:
- gene_list_file: {os.path.abspath('significant_genes_only.txt')}
- organism: hsapiens
- ontology_categories: ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG']
- pvalue_threshold: 0.05
- correction_method: fdr
- analysis_source: gprofiler
- create_plots: true"

# Analysis with background genes:
npm start -- --prompt "Perform GO analysis using these settings:
- gene_list_file: {os.path.abspath('apoptosis_genes.txt')}
- background_file: {os.path.abspath('background_genes.txt')}
- organism: hsapiens
- min_term_size: 5
- max_term_size: 300
- max_terms_plot: 15"

# Analysis from DEG results (CSV file):
npm start -- --prompt "Run GO enrichment analysis on DEG results:
- gene_list_file: {os.path.abspath('simulated_deg_results.csv')}
- gene_column: gene
- organism: hsapiens
- pvalue_threshold: 0.01
- create_plots: true"
"""
    
    print(basic_examples)
    
    print("\n" + "="*60)
    print("ANALYSIS WORKFLOW EXAMPLES")
    print("="*60)
    
    workflow_examples = """
# Complete DEG → GO Analysis Workflow:

1. First, run DEG analysis:
   npm start -- --prompt "Perform DEG analysis on my expression data"
   
2. Then, run GO analysis on significant genes:
   npm start -- --prompt "Perform GO analysis on the significant genes from DEG analysis"

# Multi-organism analysis:
npm start -- --prompt "Run GO analysis for mouse genes:
- gene_list_file: /path/to/mouse_genes.txt
- organism: mmusculus
- ontology_categories: ['GO:BP', 'GO:MF']"

# Pathway-focused analysis:
npm start -- --prompt "Analyze pathways for my gene list:
- gene_list_file: /path/to/genes.txt
- ontology_categories: ['KEGG', 'REAC', 'WP']
- analysis_source: gprofiler"
"""
    
    print(workflow_examples)
    
    print("\n" + "="*60)
    print("SUPPORTED ORGANISMS")
    print("="*60)
    
    organisms_info = """
Common organism codes for the organism parameter:
- hsapiens: Human (Homo sapiens) - default
- mmusculus: Mouse (Mus musculus)  
- rnorvegicus: Rat (Rattus norvegicus)
- drerio: Zebrafish (Danio rerio)
- dmelanogaster: Fruit fly (Drosophila melanogaster)
- celegans: Roundworm (Caenorhabditis elegans)
- scerevisiae: Baker's yeast (Saccharomyces cerevisiae)
- athaliana: Thale cress (Arabidopsis thaliana)

g:Profiler supports 500+ organisms total.
"""
    
    print(organisms_info)
    
    print("\n" + "="*60)
    print("ONTOLOGY CATEGORIES")
    print("="*60)
    
    categories_info = """
Available ontology categories:
- GO:BP: Gene Ontology Biological Process
- GO:MF: Gene Ontology Molecular Function  
- GO:CC: Gene Ontology Cellular Component
- KEGG: KEGG pathways
- REAC: Reactome pathways
- WP: WikiPathways
- TF: Transcription factor targets (TRANSFAC)
- MIRNA: miRNA targets
- HPA: Human Protein Atlas tissue expression
- CORUM: Protein complexes
- HP: Human Phenotype Ontology
"""
    
    print(categories_info)
    
    print("\n" + "="*60)
    print("OUTPUT FILES")
    print("="*60)
    
    output_info = """
The GO analysis tool generates:
- go_results.csv: Complete results for all tested terms
- enriched_terms.csv: Filtered significant terms only  
- go_summary.txt: Analysis summary and statistics
- go_barplot.png: Bar plot of top enriched terms
- go_dotplot.png: Dot plot showing gene ratios and significance
- go_network.png: Network plot of related terms
- go_analysis_script.py: Reproducible analysis script
"""
    
    print(output_info)
    
    print("\n" + "="*60)
    print("PYTHON DEPENDENCIES")
    print("="*60)
    
    dependencies_info = """
The GO analysis tool requires these Python packages:

Core packages:
- pandas: Data manipulation
- numpy: Numerical computing

Analysis packages (install at least one):
- gprofiler-official: g:Profiler API client (recommended)
  pip install gprofiler-official
  
- enrichr: Enrichr API client  
  pip install enrichr

Visualization packages (optional):
- matplotlib: Basic plotting
- seaborn: Enhanced statistical plots
- networkx: Network analysis and plotting
  pip install matplotlib seaborn networkx

Install all with:
pip install pandas numpy gprofiler-official enrichr matplotlib seaborn networkx
"""
    
    print(dependencies_info)
    
    print("\n" + "="*60)
    print("TEST DATA CREATED")
    print("="*60)
    
    print("✅ Ready to test GO analysis with these files:")
    for filename in ['cell_cycle_genes.txt', 'apoptosis_genes.txt', 'dna_repair_genes.txt', 
                     'immune_response_genes.txt', 'metabolism_genes.txt',
                     'significant_genes_only.txt', 'simulated_deg_results.csv', 'background_genes.txt']:
        if os.path.exists(filename):
            print(f"   - {filename}")

if __name__ == "__main__":
    main() 