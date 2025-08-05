# Bioinformatics Workflow Guide

The Gemini CLI provides comprehensive bioinformatics analysis capabilities through integrated tools for differential expression gene (DEG) analysis and gene ontology (GO) enrichment analysis. This guide explains how to use these tools individually and in integrated workflows.

## Overview

The bioinformatics workflow in Gemini CLI consists of two main tools that work seamlessly together:

1. **DEG Analysis Tool** (`deg_analysis`) - Identifies differentially expressed genes
2. **GO Analysis Tool** (`go_analysis`) - Performs functional enrichment analysis

These tools are automatically suggested when you mention bioinformatics-related terms in your prompts.

## Automatic Tool Recognition

The CLI automatically recognizes bioinformatics requests when you mention:

### DEG Analysis Keywords
- Gene expression, RNA-seq, microarray
- Differential expression, DEG, DGE
- Volcano plot, MA plot, fold change
- Statistical testing, multiple testing correction
- Control vs treatment, condition comparison

### GO Analysis Keywords
- Pathway analysis, functional enrichment
- Gene ontology, GO terms, GO analysis
- KEGG, Reactome, WikiPathways
- Biological processes, molecular functions, cellular components
- Overrepresentation, enrichment analysis

## Quick Start Examples

### Simple Requests (Natural Language)

```bash
# DEG Analysis
npm start -- --prompt "Analyze my gene expression data for differentially expressed genes"
npm start -- --prompt "Find genes that are significantly different between control and treatment"
npm start -- --prompt "I have RNA-seq data and need to identify DEGs"

# GO Analysis  
npm start -- --prompt "Perform pathway analysis on my gene list"
npm start -- --prompt "I have significant genes and want to understand their biological functions"
npm start -- --prompt "Do functional enrichment analysis on these genes"

# Integrated Workflow
npm start -- --prompt "Analyze my expression data and then do pathway analysis on significant genes"
```

## Complete Workflow Examples

### 1. Standard DEG → GO Pipeline

```bash
# Step 1: Differential Expression Analysis
npm start -- --prompt "Perform DEG analysis on /path/to/expression_data.csv comparing control vs treatment using condition column"

# Step 2: Functional Enrichment (automatic integration)
# The system will automatically offer to analyze significant genes for functional enrichment
# Or you can explicitly request it:
npm start -- --prompt "Now perform GO analysis on the significant genes from the DEG results"
```

### 2. Advanced Analysis with Custom Parameters

```bash
npm start -- --prompt "Run DEG analysis with these parameters:
- expression_file: /path/to/data.csv
- condition_column: treatment_group
- group1: control
- group2: drug_treatment
- log2fc_threshold: 1.5
- pvalue_threshold: 0.01
- statistical_test: mannwhitney
- correction_method: fdr
- create_plots: true

Then perform GO analysis on significant genes with:
- organism: hsapiens
- ontology_categories: ['GO:BP', 'GO:MF', 'KEGG']
- analysis_source: gprofiler"
```

### 3. Multi-Organism Analysis

```bash
# Human analysis
npm start -- --prompt "Analyze human RNA-seq data for differential expression and pathway enrichment"

# Mouse analysis
npm start -- --prompt "Perform DEG and GO analysis on mouse gene expression data using mmusculus organism"

# Cross-species comparison
npm start -- --prompt "Compare differential expression results between human and mouse datasets"
```

## Tool Integration Features

### Seamless Data Flow
- DEG analysis results automatically formatted for GO analysis input
- Significant gene lists extracted and passed between tools
- Consistent organism and identifier handling

### Coordinated Parameters
- Organism settings maintained across tools
- Statistical thresholds applied consistently
- Output directories organized hierarchically

### Comprehensive Reporting
- Combined analysis summaries
- Cross-referenced results between tools
- Integrated visualization outputs

## Input Data Formats

### DEG Analysis Input
```csv
sample_id,condition,Gene1,Gene2,Gene3,...
Sample1,control,45.2,123.4,67.8,...
Sample2,control,42.1,119.7,71.2,...
Sample3,treatment,89.3,234.1,45.6,...
Sample4,treatment,91.7,245.8,43.2,...
```

### GO Analysis Input (from DEG results)
```csv
gene,mean_group1,mean_group2,log2_fold_change,pvalue,adjusted_pvalue,significant
TP53,45.2,123.4,1.45,0.001,0.01,TRUE
BRCA1,89.1,156.7,0.81,0.005,0.02,TRUE
```

### GO Analysis Input (simple gene list)
```txt
TP53
BRCA1
BRCA2
ATM
CHEK2
```

## Output Structure

```
project_directory/
├── deg_analysis_results/
│   ├── deg_results.csv
│   ├── significant_genes.csv
│   ├── volcano_plot.png
│   ├── ma_plot.png
│   └── deg_summary.txt
└── go_analysis_results/
    ├── go_results.csv
    ├── enriched_terms.csv
    ├── go_barplot.png
    ├── go_dotplot.png
    ├── go_network.png
    └── go_summary.txt
```

## Supported Organisms

Both tools support the same organism set for seamless integration:

| Organism | Code | DEG Support | GO Support |
|----------|------|-------------|------------|
| Human | `hsapiens` | ✅ | ✅ |
| Mouse | `mmusculus` | ✅ | ✅ |
| Rat | `rnorvegicus` | ✅ | ✅ |
| Zebrafish | `drerio` | ✅ | ✅ |
| Fruit fly | `dmelanogaster` | ✅ | ✅ |
| Roundworm | `celegans` | ✅ | ✅ |
| Yeast | `scerevisiae` | ✅ | ✅ |
| Arabidopsis | `athaliana` | ✅ | ✅ |

Plus 500+ additional organisms via g:Profiler integration.

## Statistical Considerations

### DEG Analysis
- Multiple statistical tests (t-test, Mann-Whitney, Wilcoxon)
- Multiple testing correction (FDR, Bonferroni)
- Effect size thresholds (log2 fold change)
- Comprehensive diagnostic plots

### GO Analysis
- Hypergeometric enrichment testing
- Multiple testing correction across terms
- Term size filtering for interpretability
- Background gene set support

### Integrated Analysis
- Consistent statistical rigor across tools
- Coordinated significance thresholds
- Reproducible analysis scripts
- Publication-ready outputs

## Best Practices

### Data Preparation
1. **Quality Control**: Remove low-quality samples and genes
2. **Normalization**: Ensure proper data normalization
3. **Sample Size**: Use adequate biological replication (≥3 per group)
4. **Annotation**: Use current gene identifiers

### Analysis Strategy
1. **Exploratory Analysis**: Start with default parameters
2. **Parameter Tuning**: Adjust thresholds based on initial results
3. **Validation**: Cross-validate with independent methods
4. **Integration**: Combine expression and functional analysis

### Result Interpretation
1. **Statistical Significance**: Consider both p-values and effect sizes
2. **Biological Relevance**: Focus on biologically meaningful changes
3. **Pathway Context**: Interpret genes in functional context
4. **Literature Validation**: Cross-reference with existing knowledge

## Advanced Workflows

### Batch Processing
```bash
# Process multiple datasets
for dataset in dataset1.csv dataset2.csv dataset3.csv; do
  npm start -- --prompt "Run DEG and GO analysis on $dataset"
done
```

### Comparative Analysis
```bash
# Compare different conditions
npm start -- --prompt "Compare DEG results across multiple treatment conditions and identify common pathways"
```

### Time Course Analysis
```bash
# Analyze temporal expression patterns
npm start -- --prompt "Perform DEG analysis across time points and identify temporal pathway patterns"
```

### Custom Backgrounds
```bash
# Use tissue-specific backgrounds
npm start -- --prompt "Perform GO analysis with brain-specific background genes for neurological study"
```

## Integration with External Tools

### Export Options
- Results compatible with Cytoscape, Gephi for network visualization
- Integration with R/Bioconductor workflows
- Export to pathway visualization tools
- Compatible with manuscript preparation tools

### API Integration
- Programmatic access to analysis results
- Batch processing capabilities
- Integration with laboratory information systems
- Custom reporting and dashboards

## Troubleshooting

### Common Issues
1. **File Format Errors**: Ensure proper CSV/TSV formatting
2. **Gene ID Mismatches**: Use consistent identifier types
3. **Organism Specification**: Verify organism codes match data
4. **Statistical Thresholds**: Adjust for appropriate stringency

### Performance Optimization
1. **Large Datasets**: Use appropriate filtering thresholds
2. **Memory Usage**: Monitor system resources during analysis
3. **Network Dependencies**: Ensure stable internet for GO analysis
4. **Python Dependencies**: Install required packages

## Getting Help

### Documentation
- Individual tool documentation: `docs/tools/deg-analysis.md` and `docs/tools/go-analysis.md`
- Examples and tutorials: `examples/deg_analysis_example.py` and `examples/go_analysis_example.py`
- API reference: Tool parameter documentation

### Support Commands
```bash
# Get help with bioinformatics workflows
npm start -- --prompt "Help me design a differential expression analysis workflow"

# Tool-specific help
npm start -- --prompt "What parameters should I use for DEG analysis of my dataset?"
```

The Gemini CLI bioinformatics tools provide a complete, integrated solution for functional genomics analysis, from raw expression data to biological insights. 