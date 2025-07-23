# DEG Analysis Tool

The **DEG Analysis Tool** (`deg_analysis`) is a comprehensive bioinformatics tool for performing Differential Expression Gene (DEG) analysis on gene expression datasets. It identifies significantly differentially expressed genes between two experimental conditions or groups.

## Features

### Statistical Methods
- **T-test**: Parametric test assuming normal distribution (default)
- **Mann-Whitney U test**: Non-parametric test for non-normal data
- **Wilcoxon signed-rank test**: For paired sample comparisons

### Multiple Testing Correction
- **FDR (False Discovery Rate)**: Benjamini-Hochberg method (recommended, default)
- **Bonferroni correction**: Conservative correction method
- **No correction**: For single hypothesis testing (not recommended for multiple genes)

### Filtering and Thresholds
- **Log2 fold change threshold**: Minimum absolute change for significance (default: 1.0)
- **P-value threshold**: Maximum p-value for significance (default: 0.05)
- **Automatic NA removal**: Removes genes with missing values
- **Log transformation**: Optional log2(x+1) transformation for raw count data

### Output Features
- **Complete results table**: Statistics for all analyzed genes
- **Filtered significant genes**: Only genes meeting significance criteria
- **Volcano plot**: Visualization of fold change vs. significance
- **MA plot**: Average expression vs. fold change plot
- **Summary report**: Analysis overview and statistics
- **Reproducible script**: The exact Python script used for analysis

## Parameters

### Required Parameters
- `expression_file` (string): Absolute path to gene expression data file (CSV/TSV)
- `condition_column` (string): Column name containing group/condition labels
- `group1` (string): Name of the first group for comparison (e.g., "control")
- `group2` (string): Name of the second group for comparison (e.g., "treatment")

### Optional Parameters
- `gene_id_column` (string): Column containing gene identifiers (auto-detected if not specified)
- `log2fc_threshold` (number): Minimum absolute log2 fold change for significance (default: 1.0)
- `pvalue_threshold` (number): Maximum p-value for significance (default: 0.05)
- `correction_method` (string): Multiple testing correction ("fdr", "bonferroni", "none"; default: "fdr")
- `statistical_test` (string): Statistical test to use ("ttest", "mannwhitney", "wilcoxon"; default: "ttest")
- `output_directory` (string): Directory for saving results (default: "./deg_analysis_results")
- `create_plots` (boolean): Generate visualization plots (default: true)
- `remove_na` (boolean): Remove genes with missing values (default: true)
- `log_transform` (boolean): Apply log2 transformation (default: false)
- `separator` (string): File separator (",", "\t", ";"; auto-detected if not specified)

## Input Data Format

The tool accepts CSV or TSV files with:
- **Samples as rows OR columns** (automatically detected)
- **Numeric expression values** for genes
- **Clear condition/group labels** in the specified condition column
- **Optional gene identifiers** column

### Example Data Structure
```csv
sample_id,condition,Gene1,Gene2,Gene3,...
Sample1,control,45.2,123.4,67.8,...
Sample2,control,42.1,119.7,71.2,...
Sample3,treatment,89.3,234.1,45.6,...
Sample4,treatment,91.7,245.8,43.2,...
```

## Usage Examples

### Basic Usage
```bash
# Simple DEG analysis with default parameters
npm start -- --prompt "Perform DEG analysis on /path/to/data.csv comparing control vs treatment groups using the condition column"
```

### Advanced Usage with Custom Parameters
```bash
npm start -- --prompt "Use deg_analysis tool with parameters:
- expression_file: /absolute/path/to/expression_data.csv
- condition_column: treatment_group
- group1: control
- group2: drug_treatment
- log2fc_threshold: 1.5
- pvalue_threshold: 0.01
- correction_method: fdr
- statistical_test: mannwhitney
- create_plots: true
- output_directory: /path/to/results"
```

### Programmatic Usage
```javascript
// Example parameters for the tool
const params = {
  expression_file: "/absolute/path/to/gene_expression.csv",
  condition_column: "condition",
  group1: "control",
  group2: "treatment",
  log2fc_threshold: 1.0,
  pvalue_threshold: 0.05,
  correction_method: "fdr",
  statistical_test: "ttest",
  create_plots: true,
  remove_na: true,
  log_transform: false
};
```

## Output Files

The tool generates several output files in the specified directory:

### Results Files
- **`deg_results.csv`**: Complete results for all analyzed genes
  - Columns: gene, mean_group1, mean_group2, log2_fold_change, statistic, pvalue, adjusted_pvalue, significant, regulation
- **`significant_genes.csv`**: Filtered results containing only significant genes
- **`deg_summary.txt`**: Analysis summary with statistics and parameters

### Visualization Files (if `create_plots: true`)
- **`volcano_plot.png`**: Volcano plot showing log2 fold change vs. -log10(p-value)
- **`ma_plot.png`**: MA plot showing average expression vs. log2 fold change

### Analysis Script
- **`deg_analysis_script.py`**: The complete Python script used for the analysis (reproducible)

## System Requirements

### Python Dependencies
The tool requires Python 3 with the following packages:
```bash
pip install pandas numpy scipy matplotlib seaborn
```

### Required Packages
- **pandas**: Data manipulation and analysis
- **numpy**: Numerical computing
- **scipy**: Statistical functions

### Optional Packages (for plotting)
- **matplotlib**: Basic plotting functionality
- **seaborn**: Enhanced statistical visualizations

## Statistical Considerations

### Choosing Statistical Tests
- **T-test**: Use when data is approximately normally distributed
- **Mann-Whitney U**: Use for non-normal data or ordinal data
- **Wilcoxon**: Use for paired samples (same subjects before/after treatment)

### Multiple Testing Correction
- **FDR (recommended)**: Balances Type I and Type II error rates
- **Bonferroni**: Very conservative, use when Type I error must be minimized
- **None**: Only for single hypothesis or exploratory analysis

### Threshold Selection
- **Log2FC threshold**: 
  - 1.0 = 2-fold change (recommended for most analyses)
  - 0.58 = 1.5-fold change (less stringent)
  - 2.0 = 4-fold change (very stringent)
- **P-value threshold**:
  - 0.05 = Standard significance level
  - 0.01 = More stringent
  - 0.1 = Less stringent (exploratory)

## Best Practices

### Data Preparation
1. **Quality Control**: Remove low-quality samples and lowly expressed genes
2. **Normalization**: Ensure data is properly normalized (e.g., TPM, FPKM, or normalized counts)
3. **Log Transformation**: Use `log_transform: true` for raw count data
4. **Sample Size**: Ensure adequate replication (minimum 3 samples per group, ideally 6+)

### Analysis Parameters
1. **Choose appropriate statistical test** based on data distribution
2. **Use FDR correction** for multiple testing unless specific reasons for other methods
3. **Set meaningful thresholds** based on biological significance, not just statistical significance
4. **Consider effect sizes** alongside p-values

### Result Interpretation
1. **Examine both fold change and p-values** for biological relevance
2. **Validate top candidates** with independent methods (qPCR, etc.)
3. **Consider biological context** when interpreting results
4. **Check for batch effects** or confounding variables

## Troubleshooting

### Common Issues

**File Reading Errors**
- Ensure file path is absolute
- Check file encoding (UTF-8 recommended)
- Verify separator is correctly detected or specified

**No Significant Genes Found**
- Check if thresholds are too stringent
- Verify group labels match exactly (case-sensitive)
- Ensure sufficient sample size and biological effect

**Python/Package Errors**
- Install required Python packages: `pip install pandas numpy scipy matplotlib seaborn`
- Ensure Python 3 is available as `python3`
- Check that output directory is writable

**Statistical Warnings**
- Low sample sizes may produce unreliable results
- Consider non-parametric tests for non-normal data
- Check for outliers that may affect results

### Error Messages
- **"Condition column not found"**: Check spelling and case of condition column name
- **"One or both groups have no samples"**: Verify group names match data exactly
- **"No genes could be analyzed"**: Check data format and ensure numeric expression columns exist

## Advanced Usage

### Custom Workflows
The tool can be integrated into larger bioinformatics pipelines:

1. **Quality Control** → DEG Analysis → **Pathway Analysis**
2. **Normalization** → DEG Analysis → **Gene Set Enrichment**
3. **Batch Correction** → DEG Analysis → **Functional Annotation**

### Extending the Analysis
The generated Python script can be modified for:
- Custom visualization themes
- Additional statistical tests
- Integration with other analysis tools
- Batch processing of multiple comparisons

## Examples and Tutorials

See `examples/deg_analysis_example.py` for:
- Synthetic test data generation
- Step-by-step usage examples
- Expected output interpretation
- Best practices demonstration

## Version History

- **v1.0.0**: Initial release with basic DEG analysis functionality
  - Multiple statistical tests (t-test, Mann-Whitney, Wilcoxon)
  - Multiple testing correction methods (FDR, Bonferroni)
  - Visualization support (volcano plots, MA plots)
  - Comprehensive output and reporting 