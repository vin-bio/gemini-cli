# GO Analysis Tool

The **GO Analysis Tool** (`go_analysis`) is a comprehensive bioinformatics tool for performing Gene Ontology (GO) enrichment analysis on gene lists. It identifies overrepresented biological processes, molecular functions, cellular components, and pathways in your gene sets.

## Overview

Gene Ontology enrichment analysis is a critical step in functional genomics research to understand the biological meaning of gene lists derived from experiments like differential expression analysis, GWAS studies, or proteomics screens. This tool provides a streamlined interface to perform robust GO enrichment analysis with multiple databases and visualization options.

## Features

### Analysis Sources
- **g:Profiler**: Comprehensive multi-organism functional profiling (default, recommended)
- **Enrichr**: Popular web-based enrichment analysis platform
- **DAVID**: Database for Annotation, Visualization and Integrated Discovery

### Organism Support
- **500+ organisms** supported via g:Profiler
- **Human** (Homo sapiens) - default
- **Mouse** (Mus musculus)
- **Rat** (Rattus norvegicus)
- **Zebrafish** (Danio rerio)
- **Fruit fly** (Drosophila melanogaster)
- **Roundworm** (C. elegans)
- **Yeast** (S. cerevisiae)
- **Plant models** (A. thaliana) and many more

### Ontology Categories
- **GO:BP**: Gene Ontology Biological Process
- **GO:MF**: Gene Ontology Molecular Function
- **GO:CC**: Gene Ontology Cellular Component
- **KEGG**: KEGG pathways
- **Reactome**: Reactome pathways
- **WikiPathways**: Community-curated pathways
- **Additional sources**: TRANSFAC, miRNA targets, protein complexes, phenotypes

### Gene ID Types
- **Gene symbols** (HGNC, MGI, etc.) - default and recommended
- **Ensembl gene IDs** (ENSG, ENSMUSG, etc.)
- **Entrez gene IDs**
- **UniProt IDs**
- **RefSeq IDs**
- **Custom organism-specific identifiers**

### Statistical Methods
- **Hypergeometric test** for enrichment significance
- **Multiple testing correction**: FDR (Benjamini-Hochberg), Bonferroni, or none
- **Customizable significance thresholds**
- **Term size filtering** to exclude very small or very large terms

### Visualization Features
- **Bar plots**: Top enriched terms with significance levels
- **Dot plots**: Gene ratio vs. significance with term sizes
- **Network plots**: Relationships between enriched terms
- **Interactive outputs** (when supported by analysis source)

## Parameters

### Required Parameters
- `gene_list_file` (string): Absolute path to file containing gene list

### Optional Parameters

#### Input Configuration
- `organism` (string): Target organism code (default: "hsapiens")
- `gene_id_type` (string): Type of gene identifiers (default: "HGNC")
- `gene_column` (string): Column name for genes in CSV/TSV files
- `separator` (string): File separator (auto-detected if not specified)
- `background_file` (string): Path to background gene set file

#### Analysis Parameters
- `ontology_categories` (array): Categories to analyze (default: ["GO:BP", "GO:MF", "GO:CC"])
- `analysis_source` (string): Analysis database ("gprofiler", "enrichr", "david"; default: "gprofiler")
- `pvalue_threshold` (number): Maximum p-value for significance (default: 0.05)
- `correction_method` (string): Multiple testing correction ("fdr", "bonferroni", "none"; default: "fdr")

#### Filtering Parameters
- `min_term_size` (number): Minimum genes per term (default: 5)
- `max_term_size` (number): Maximum genes per term (default: 500)
- `max_terms_plot` (number): Maximum terms in plots (default: 20)

#### Output Parameters
- `output_directory` (string): Directory for results (default: "./go_analysis_results")
- `create_plots` (boolean): Generate visualizations (default: true)

## Input Data Formats

### Simple Gene Lists
```txt
TP53
BRCA1
BRCA2
ATM
CHEK2
```

### CSV/TSV with Headers
```csv
gene,log2fc,pvalue
TP53,2.3,0.001
BRCA1,1.8,0.005
BRCA2,-1.2,0.01
```

### DEG Analysis Results
```csv
gene,mean_group1,mean_group2,log2_fold_change,pvalue,adjusted_pvalue,significant
TP53,45.2,123.4,1.45,0.001,0.01,TRUE
BRCA1,89.1,156.7,0.81,0.005,0.02,TRUE
```

## Usage Examples

### Basic Usage
```bash
# Simple GO analysis with default parameters
npm start -- --prompt "Perform GO enrichment analysis on /path/to/genes.txt using human organism"
```

### Advanced Analysis
```bash
npm start -- --prompt "Use go_analysis tool with these parameters:
- gene_list_file: /absolute/path/to/gene_list.txt
- organism: hsapiens
- ontology_categories: ['GO:BP', 'GO:MF', 'KEGG', 'REAC']
- pvalue_threshold: 0.01
- correction_method: fdr
- analysis_source: gprofiler
- create_plots: true
- max_terms_plot: 25"
```

### Multi-Organism Analysis
```bash
# Mouse gene analysis
npm start -- --prompt "Analyze mouse genes for GO enrichment:
- gene_list_file: /path/to/mouse_genes.txt
- organism: mmusculus
- ontology_categories: ['GO:BP', 'GO:MF', 'GO:CC']"

# Zebrafish analysis
npm start -- --prompt "Run GO analysis for zebrafish:
- gene_list_file: /path/to/zebrafish_genes.txt  
- organism: drerio
- analysis_source: gprofiler"
```

### Pathway-Focused Analysis
```bash
npm start -- --prompt "Analyze pathways only:
- gene_list_file: /path/to/genes.txt
- ontology_categories: ['KEGG', 'REAC', 'WP']
- min_term_size: 10
- max_term_size: 200"
```

### Analysis with Background
```bash
npm start -- --prompt "GO analysis with custom background:
- gene_list_file: /path/to/significant_genes.txt
- background_file: /path/to/all_expressed_genes.txt
- organism: hsapiens
- pvalue_threshold: 0.05"
```

### Integration with DEG Results
```bash
npm start -- --prompt "Analyze DEG results for GO enrichment:
- gene_list_file: /path/to/deg_results.csv
- gene_column: gene
- organism: hsapiens
- ontology_categories: ['GO:BP', 'KEGG']
- create_plots: true"
```

## Output Files

### Results Files
- **`go_results.csv`**: Complete results for all tested terms
  - Columns: term_id, term_name, database, pvalue, adjusted_pvalue, term_size, intersection_size, genes
- **`enriched_terms.csv`**: Filtered results containing only significant terms
- **`go_summary.txt`**: Analysis summary with statistics and parameters

### Visualization Files (if `create_plots: true`)
- **`go_barplot.png`**: Horizontal bar plot of top enriched terms
- **`go_dotplot.png`**: Dot plot showing gene ratios and significance levels  
- **`go_network.png`**: Network plot of related terms (if networkx available)

### Analysis Script
- **`go_analysis_script.py`**: Complete Python script used for analysis (reproducible)

## System Requirements

### Python Dependencies

#### Core Packages (Required)
```bash
pip install pandas numpy
```

#### Analysis Packages (Install at least one)
```bash
# g:Profiler (recommended)
pip install gprofiler-official

# Enrichr
pip install enrichr

# Statistical testing (if not using g:Profiler)
pip install statsmodels
```

#### Visualization Packages (Optional)
```bash
pip install matplotlib seaborn networkx
```

#### Complete Installation
```bash
pip install pandas numpy gprofiler-official enrichr matplotlib seaborn networkx statsmodels
```

### System Requirements
- **Python 3.7+**
- **Internet connection** (for online databases)
- **Writable output directory**

## Organism Codes Reference

### Common Model Organisms
| Organism | Code | Common Name |
|----------|------|-------------|
| Human | `hsapiens` | Homo sapiens |
| Mouse | `mmusculus` | Mus musculus |
| Rat | `rnorvegicus` | Rattus norvegicus |
| Zebrafish | `drerio` | Danio rerio |
| Fruit fly | `dmelanogaster` | Drosophila melanogaster |
| Roundworm | `celegans` | Caenorhabditis elegans |
| Yeast | `scerevisiae` | Saccharomyces cerevisiae |
| Arabidopsis | `athaliana` | Arabidopsis thaliana |

### Other Organisms
g:Profiler supports 500+ organisms. Use the organism's scientific name in lowercase with spaces replaced by underscores (e.g., `gallus_gallus` for chicken).

## Best Practices

### Gene List Preparation
1. **Use official gene symbols** when possible (HGNC for human, MGI for mouse)
2. **Remove duplicates** from your gene list
3. **Ensure gene identifiers match** the selected organism and ID type
4. **Consider gene list size**: 10-500 genes typically work best
5. **Filter for expressed genes** if using RNA-seq data

### Analysis Configuration
1. **Choose appropriate organism** matching your experimental system
2. **Use FDR correction** for multiple testing (recommended default)
3. **Set meaningful thresholds** based on your research context
4. **Include relevant ontology categories** for your research question
5. **Use background genes** when analyzing a subset of a larger dataset

### Statistical Considerations
1. **Multiple testing correction** is essential when testing many terms
2. **Term size filtering** helps focus on interpretable categories
3. **P-value thresholds** should be adjusted based on study goals
4. **Effect sizes matter**: consider the number of genes in each term
5. **Biological significance** may differ from statistical significance

### Result Interpretation
1. **Examine term hierarchy**: specific terms are often more informative
2. **Consider term redundancy**: related terms often co-occur
3. **Validate key findings** with literature or independent methods
4. **Focus on reproducible patterns** across similar experiments
5. **Integrate with pathway databases** for mechanistic insights

## Analysis Workflows

### Standard DEG → GO Workflow
```bash
# 1. Perform differential expression analysis
npm start -- --prompt "Run DEG analysis on expression data"

# 2. Extract significant genes
# (Use significant_genes.csv from DEG results)

# 3. Perform GO enrichment analysis
npm start -- --prompt "Analyze significant genes for GO enrichment:
- gene_list_file: /path/to/significant_genes.csv
- gene_column: gene
- organism: hsapiens"
```

### Comparative Analysis Workflow
```bash
# 1. Analyze upregulated genes
npm start -- --prompt "GO analysis for upregulated genes:
- gene_list_file: /path/to/upregulated_genes.txt"

# 2. Analyze downregulated genes  
npm start -- --prompt "GO analysis for downregulated genes:
- gene_list_file: /path/to/downregulated_genes.txt"

# 3. Compare results between conditions
```

### Cross-Species Analysis
```bash
# 1. Human analysis
npm start -- --prompt "GO analysis for human genes:
- gene_list_file: /path/to/human_genes.txt
- organism: hsapiens"

# 2. Mouse ortholog analysis
npm start -- --prompt "GO analysis for mouse orthologs:
- gene_list_file: /path/to/mouse_orthologs.txt  
- organism: mmusculus"
```

## Troubleshooting

### Common Issues

**No Enriched Terms Found**
- Check gene identifier format and organism match
- Verify gene list size (too small or too large lists may fail)
- Adjust significance thresholds
- Ensure genes are well-annotated in the database

**Gene ID Conversion Errors**
- Use official gene symbols when possible
- Check for deprecated or outdated identifiers
- Verify organism-specific ID formats
- Consider manual curation of problematic genes

**Network/API Errors**
- Ensure stable internet connection
- Check if analysis service is available
- Try alternative analysis sources (gprofiler, enrichr)
- Implement retry logic for transient failures

**Python Package Errors**
- Install required packages: `pip install gprofiler-official enrichr`
- Update packages to latest versions
- Check Python version compatibility (3.7+)
- Verify package dependencies are met

### Error Messages

**"No genes found in input file"**
- Check file format and gene column specification
- Verify file path is absolute and accessible
- Ensure file contains valid gene identifiers

**"Analysis source not available"**
- Install required Python packages for selected source
- Check internet connectivity for online databases
- Try alternative analysis sources

**"Organism not supported"**
- Verify organism code spelling and format
- Check if organism is supported by selected database
- Use alternative organism codes or databases

## Advanced Features

### Custom Background Sets
```bash
npm start -- --prompt "GO analysis with tissue-specific background:
- gene_list_file: /path/to/brain_degs.txt
- background_file: /path/to/brain_expressed_genes.txt"
```

### Batch Processing
```bash
# Process multiple gene lists
for list in list1.txt list2.txt list3.txt; do
  npm start -- --prompt "GO analysis for $list"
done
```

### Integration with Other Tools
- **Pathway visualization**: Export results to Cytoscape, Gephi
- **Gene set databases**: Compare with MSigDB, Reactome
- **Literature mining**: Cross-reference with PubMed, GeneRIF
- **Protein interactions**: Integrate with STRING, BioGRID

## Examples and Tutorials

See `examples/go_analysis_example.py` for:
- Test gene list generation
- Usage examples for different scenarios
- Integration with DEG analysis results
- Multi-organism analysis examples
- Best practices demonstration

## Version History

- **v1.0.0**: Initial release with GO enrichment functionality
  - Multiple analysis sources (g:Profiler, Enrichr, DAVID)
  - 500+ organism support via g:Profiler
  - Comprehensive ontology category support
  - Statistical testing with multiple correction methods
  - Visualization and reporting features
  - Integration with DEG analysis workflow 