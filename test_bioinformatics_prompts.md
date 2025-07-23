# Bioinformatics Tool Integration Test

This file contains test prompts to verify that the Gemini CLI correctly recognizes bioinformatics requests and suggests the appropriate DEG and GO analysis tools.

## Test Prompts for DEG Analysis

### Basic DEG Analysis Requests
- "Analyze my gene expression data for differentially expressed genes"
- "I have RNA-seq data and need to find DEGs between control and treatment"
- "Perform differential expression analysis on my expression matrix"
- "Compare gene expression between two conditions"
- "Find significantly up and down regulated genes"

### Advanced DEG Analysis Requests
- "Run DEG analysis with FDR correction and log2FC threshold of 1.5"
- "Analyze differential expression with Mann-Whitney test for non-normal data"
- "Compare multiple groups in my expression dataset"
- "Generate volcano plots for my differential expression results"

## Test Prompts for GO Analysis

### Basic GO Enrichment Requests
- "Perform GO enrichment analysis on my gene list"
- "I have significant genes and need pathway analysis"
- "Do functional enrichment analysis on these genes"
- "Find enriched biological processes in my gene set"
- "Analyze KEGG pathways for my gene list"

### Advanced GO Analysis Requests
- "Run GO analysis with g:Profiler for mouse genes"
- "Perform enrichment analysis for biological processes and molecular functions"
- "Analyze pathways with custom background gene set"
- "Generate network plots for enriched terms"

## Test Prompts for Integrated Workflows

### DEG → GO Workflow
- "Analyze my expression data for DEGs and then do pathway analysis"
- "Run complete differential expression workflow with functional analysis"
- "Find differentially expressed genes and analyze their biological functions"
- "Do DEG analysis followed by GO enrichment on significant genes"

### Multi-organism Analysis
- "Analyze mouse gene expression data for differential expression"
- "Perform GO analysis on zebrafish genes"
- "Run DEG analysis on human RNA-seq data"
- "Do pathway analysis for yeast gene list"

## Expected Behavior

When users provide these prompts, the system should:

1. **Recognize bioinformatics context** from keywords like:
   - Gene expression, RNA-seq, microarray
   - DEG, DGE, differential expression
   - GO, pathway, enrichment, functional analysis
   - KEGG, Reactome, biological processes

2. **Suggest appropriate tools**:
   - `deg_analysis` tool for expression data analysis
   - `go_analysis` tool for gene list enrichment
   - Sequential workflow for integrated analysis

3. **Provide comprehensive parameters** including:
   - Statistical tests and correction methods
   - Organism specification
   - Output visualization options
   - Integration between tools

4. **Generate complete results** with:
   - Statistical tables and summaries
   - Visualizations (volcano plots, bar plots)
   - Reproducible analysis scripts

## Test Commands

To test these integrations:

```bash
# Test DEG analysis recognition
npm start -- --prompt "Analyze my gene expression data for differentially expressed genes"

# Test GO analysis recognition  
npm start -- --prompt "I have a list of genes and need pathway analysis"

# Test integrated workflow recognition
npm start -- --prompt "Run DEG analysis on my data then do functional enrichment"
```

The system should automatically suggest and use the `deg_analysis` and `go_analysis` tools without requiring explicit tool names in the prompts. 