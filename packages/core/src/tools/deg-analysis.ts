/**
 * @license
 * Copyright 2025 Google LLC
 * SPDX-License-Identifier: Apache-2.0
 */

import fs from 'fs/promises';
import path from 'path';
import { Config } from '../config/config.js';
import {
  BaseTool,
  ToolResult,
  ToolCallConfirmationDetails,
  ToolExecuteConfirmationDetails,
  ToolConfirmationOutcome,
  ToolLocation,
  Icon,
} from './tools.js';
import { Type } from '@google/genai';
import { SchemaValidator } from '../utils/schemaValidator.js';
import { getErrorMessage } from '../utils/errors.js';

export interface DegAnalysisToolParams {
  expression_file: string;
  condition_column: string;
  group1: string;
  group2: string;
  gene_id_column?: string;
  log2fc_threshold?: number;
  pvalue_threshold?: number;
  correction_method?: 'fdr' | 'bonferroni' | 'none';
  statistical_test?: 'ttest' | 'mannwhitney' | 'wilcoxon';
  output_directory?: string;
  create_plots?: boolean;
  remove_na?: boolean;
  log_transform?: boolean;
  separator?: string;
}

export interface DegAnalysisResult extends ToolResult {
  significant_genes?: number;
  upregulated_genes?: number;
  downregulated_genes?: number;
  output_files?: string[];
}

/**
 * A tool for performing differential expression gene analysis on gene expression data.
 * Supports various statistical tests, multiple testing corrections, and visualization.
 */
export class DegAnalysisTool extends BaseTool<DegAnalysisToolParams, DegAnalysisResult> {
  static readonly Name: string = 'deg_analysis';

  constructor(private readonly config: Config) {
    super(
      DegAnalysisTool.Name,
      'DEG Analysis',
      `Performs comprehensive Differential Expression Gene (DEG) analysis on gene expression data. 
      
This tool analyzes gene expression datasets to identify significantly differentially expressed genes between two conditions or groups. It supports:

**Statistical Methods:**
- T-test (parametric, assumes normal distribution)
- Mann-Whitney U test (non-parametric)
- Wilcoxon signed-rank test (paired samples)

**Multiple Testing Correction:**
- False Discovery Rate (FDR/Benjamini-Hochberg) - Recommended
- Bonferroni correction (conservative)
- No correction (not recommended for multiple genes)

**Filtering Options:**
- Log2 fold change thresholds
- P-value/adjusted p-value thresholds
- Automatic removal of missing values
- Optional log transformation

**Output Features:**
- Results table with statistics for all genes
- Filtered significant genes list
- Optional volcano plot visualization
- Summary statistics report

**Input Requirements:**
- Expression data in CSV/TSV format with genes as rows or columns
- Clear condition/group labels for samples
- Numeric expression values

The tool generates comprehensive results including fold changes, p-values, adjusted p-values, and significance classifications for biological interpretation.`,
      Icon.LightBulb,
      {
        type: Type.OBJECT,
        properties: {
          expression_file: {
            type: Type.STRING,
            description: 'Absolute path to the gene expression data file (CSV or TSV format). Should contain numeric expression values with genes and samples clearly labeled.',
          },
          condition_column: {
            type: Type.STRING,
            description: 'Name of the column containing condition/group labels that distinguish the two groups being compared.',
          },
          group1: {
            type: Type.STRING,
            description: 'Name/label of the first group/condition for comparison (e.g., "control", "untreated").',
          },
          group2: {
            type: Type.STRING,
            description: 'Name/label of the second group/condition for comparison (e.g., "treatment", "disease").',
          },
          gene_id_column: {
            type: Type.STRING,
            description: 'Optional: Name of the column containing gene identifiers. If not specified, will attempt to auto-detect or use row indices.',
          },
          log2fc_threshold: {
            type: Type.NUMBER,
            description: 'Optional: Minimum absolute log2 fold change threshold for significance (default: 1.0). Genes with |log2FC| >= this value are considered significant.',
          },
          pvalue_threshold: {
            type: Type.NUMBER,
            description: 'Optional: Maximum p-value threshold for significance (default: 0.05). Applied after multiple testing correction if used.',
          },
          correction_method: {
            type: Type.STRING,
            description: 'Optional: Multiple testing correction method. Options: "fdr" (Benjamini-Hochberg, recommended), "bonferroni" (conservative), "none" (not recommended). Default: "fdr".',
            enum: ['fdr', 'bonferroni', 'none'],
          },
          statistical_test: {
            type: Type.STRING,
            description: 'Optional: Statistical test to use. Options: "ttest" (parametric), "mannwhitney" (non-parametric), "wilcoxon" (paired). Default: "ttest".',
            enum: ['ttest', 'mannwhitney', 'wilcoxon'],
          },
          output_directory: {
            type: Type.STRING,
            description: 'Optional: Directory path where results will be saved. If not specified, creates "deg_analysis_results" in the current directory.',
          },
          create_plots: {
            type: Type.BOOLEAN,
            description: 'Optional: Whether to generate visualization plots (volcano plot, MA plot). Default: true.',
          },
          remove_na: {
            type: Type.BOOLEAN,
            description: 'Optional: Whether to remove genes with missing (NA/NaN) values before analysis. Default: true.',
          },
          log_transform: {
            type: Type.BOOLEAN,
            description: 'Optional: Whether to apply log2 transformation to expression values. Use if data is not already log-transformed. Default: false.',
          },
          separator: {
            type: Type.STRING,
            description: 'Optional: Field separator for the input file. Options: "," (CSV), "\\t" (TSV), ";" (semicolon). Auto-detected if not specified.',
          },
        },
        required: ['expression_file', 'condition_column', 'group1', 'group2'],
      },
    );
  }

  validateToolParams(params: DegAnalysisToolParams): string | null {
    const validationError = SchemaValidator.validate(this.parameterSchema, params);
    if (validationError) {
      return `Invalid parameters: ${validationError}`;
    }

    // Additional validation
    if (params.log2fc_threshold !== undefined && params.log2fc_threshold < 0) {
      return 'log2fc_threshold must be non-negative';
    }

    if (params.pvalue_threshold !== undefined && (params.pvalue_threshold < 0 || params.pvalue_threshold > 1)) {
      return 'pvalue_threshold must be between 0 and 1';
    }

    if (params.group1 === params.group2) {
      return 'group1 and group2 must be different';
    }

    return null;
  }

  getDescription(params: DegAnalysisToolParams): string {
    const test = params.statistical_test || 'ttest';
    const correction = params.correction_method || 'fdr';
    const log2fc = params.log2fc_threshold || 1.0;
    const pval = params.pvalue_threshold || 0.05;
    
    return `Performing DEG analysis comparing "${params.group1}" vs "${params.group2}" using ${test} test with ${correction} correction (|log2FC| ≥ ${log2fc}, p < ${pval})`;
  }

  toolLocations(params: DegAnalysisToolParams): ToolLocation[] {
    const locations = [
      { path: params.expression_file }
    ];

    const outputDir = params.output_directory || path.join(process.cwd(), 'deg_analysis_results');
    locations.push({ path: outputDir });

    return locations;
  }

  async shouldConfirmExecute(
    params: DegAnalysisToolParams,
    abortSignal: AbortSignal,
  ): Promise<ToolCallConfirmationDetails | false> {
    const validationError = this.validateToolParams(params);
    if (validationError) {
      return false;
    }

    // Check if input file exists
    try {
      await fs.access(params.expression_file);
    } catch {
      return false;
    }

    const description = this.getDescription(params);
    const confirmationDetails: ToolExecuteConfirmationDetails = {
      type: 'exec',
      title: 'Confirm DEG Analysis',
      command: description,
      rootCommand: 'deg_analysis',
      onConfirm: async (outcome: ToolConfirmationOutcome) => {
        // No special handling needed for confirmation
      },
    };
    return confirmationDetails;
  }

  async execute(
    params: DegAnalysisToolParams,
    signal: AbortSignal,
    updateOutput?: (output: string) => void,
  ): Promise<DegAnalysisResult> {
    const validationError = this.validateToolParams(params);
    if (validationError) {
      return {
        llmContent: `Error: Invalid parameters. ${validationError}`,
        returnDisplay: validationError,
      };
    }

    try {
      updateOutput?.('Starting DEG analysis...\n');

      // Set default parameters
      const log2fcThreshold = params.log2fc_threshold ?? 1.0;
      const pvalueThreshold = params.pvalue_threshold ?? 0.05;
      const correctionMethod = params.correction_method ?? 'fdr';
      const statisticalTest = params.statistical_test ?? 'ttest';
      const outputDir = params.output_directory ?? path.join(process.cwd(), 'deg_analysis_results');
      const createPlots = params.create_plots ?? true;
      const removeNA = params.remove_na ?? true;
      const logTransform = params.log_transform ?? false;

      // Create output directory
      await fs.mkdir(outputDir, { recursive: true });

      // Generate Python script for DEG analysis
      const pythonScript = this.generatePythonScript(params, {
        log2fcThreshold,
        pvalueThreshold,
        correctionMethod,
        statisticalTest,
        outputDir,
        createPlots,
        removeNA,
        logTransform,
      });

      const scriptPath = path.join(outputDir, 'deg_analysis_script.py');
      await fs.writeFile(scriptPath, pythonScript);

      updateOutput?.('Generated analysis script. Running DEG analysis...\n');

      // Execute the Python script
      const { spawn } = await import('child_process');
      
      return new Promise<DegAnalysisResult>((resolve) => {
        const python = spawn('python3', [scriptPath], {
          cwd: outputDir,
          stdio: ['pipe', 'pipe', 'pipe'],
        });

        let stdout = '';
        let stderr = '';

        python.stdout.on('data', (data) => {
          const chunk = data.toString();
          stdout += chunk;
          updateOutput?.(chunk);
        });

        python.stderr.on('data', (data) => {
          stderr += data.toString();
        });

        python.on('close', async (code) => {
          if (code !== 0) {
            resolve({
              llmContent: `DEG analysis failed with exit code ${code}.\n\nStderr: ${stderr}\n\nStdout: ${stdout}`,
              returnDisplay: `Analysis failed. Check the error output for details.`,
            });
            return;
          }

          try {
            // Read results summary
            const summaryPath = path.join(outputDir, 'deg_summary.txt');
            const summaryContent = await fs.readFile(summaryPath, 'utf-8');
            
            // Parse summary for numbers
            const summaryLines = summaryContent.split('\n');
            let significantGenes = 0;
            let upregulatedGenes = 0;
            let downregulatedGenes = 0;

            for (const line of summaryLines) {
              if (line.includes('Total significant genes:')) {
                significantGenes = parseInt(line.split(':')[1].trim()) || 0;
              } else if (line.includes('Upregulated genes:')) {
                upregulatedGenes = parseInt(line.split(':')[1].trim()) || 0;
              } else if (line.includes('Downregulated genes:')) {
                downregulatedGenes = parseInt(line.split(':')[1].trim()) || 0;
              }
            }

            // List output files
            const outputFiles = [
              'deg_results.csv',
              'significant_genes.csv',
              'deg_summary.txt',
              'deg_analysis_script.py'
            ];

            if (createPlots) {
              outputFiles.push('volcano_plot.png', 'ma_plot.png');
            }

            const existingFiles = [];
            for (const file of outputFiles) {
              try {
                await fs.access(path.join(outputDir, file));
                existingFiles.push(path.join(outputDir, file));
              } catch {
                // File doesn't exist, skip
              }
            }

            const resultContent = `# DEG Analysis Complete

${summaryContent}

## Output Files:
${existingFiles.map(f => `- ${f}`).join('\n')}

## Analysis Parameters:
- Statistical test: ${statisticalTest}
- Multiple testing correction: ${correctionMethod}
- Log2FC threshold: ${log2fcThreshold}
- P-value threshold: ${pvalueThreshold}
- Groups compared: ${params.group1} vs ${params.group2}

The analysis results are saved in: ${outputDir}
`;

            resolve({
              llmContent: resultContent,
              returnDisplay: resultContent,
              significant_genes: significantGenes,
              upregulated_genes: upregulatedGenes,
              downregulated_genes: downregulatedGenes,
              output_files: existingFiles,
            });
          } catch (error) {
            resolve({
              llmContent: `Analysis completed but failed to read results: ${getErrorMessage(error)}`,
              returnDisplay: 'Analysis may have completed but results could not be read.',
            });
          }
        });

        // Handle abort signal
        signal.addEventListener('abort', () => {
          python.kill('SIGTERM');
          resolve({
            llmContent: 'DEG analysis was aborted by user.',
            returnDisplay: 'Analysis aborted.',
          });
        });
      });

    } catch (error) {
      return {
        llmContent: `Error performing DEG analysis: ${getErrorMessage(error)}`,
        returnDisplay: `Analysis failed: ${getErrorMessage(error)}`,
      };
    }
  }

  private generatePythonScript(
    params: DegAnalysisToolParams,
    options: {
      log2fcThreshold: number;
      pvalueThreshold: number;
      correctionMethod: string;
      statisticalTest: string;
      outputDir: string;
      createPlots: boolean;
      removeNA: boolean;
      logTransform: boolean;
    }
  ): string {
    return `#!/usr/bin/env python3
"""
Differential Expression Gene Analysis Script
Generated by Gemini CLI DEG Analysis Tool
"""

import pandas as pd
import numpy as np
import sys
import os
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Optional imports for plotting
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    PLOTTING_AVAILABLE = True
    # Configure matplotlib for non-interactive use
    plt.switch_backend('Agg')
except ImportError:
    PLOTTING_AVAILABLE = False
    print("Warning: matplotlib/seaborn not available. Plots will be skipped.")

def detect_separator(file_path):
    """Auto-detect file separator"""
    with open(file_path, 'r') as f:
        first_line = f.readline()
    
    separators = [',', '\\t', ';', '|']
    counts = {sep: first_line.count(sep) for sep in separators}
    return max(counts, key=counts.get) if max(counts.values()) > 0 else ','

def load_data(file_path, separator=None):
    """Load gene expression data"""
    if separator is None:
        separator = detect_separator(file_path)
    
    print(f"Loading data from: {file_path}")
    print(f"Using separator: {repr(separator)}")
    
    # Try different encodings
    for encoding in ['utf-8', 'latin-1', 'cp1252']:
        try:
            df = pd.read_csv(file_path, sep=separator, encoding=encoding, low_memory=False)
            print(f"Successfully loaded with encoding: {encoding}")
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError("Could not read file with any common encoding")
    
    print(f"Data shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    return df

def preprocess_data(df, gene_id_column, condition_column, remove_na, log_transform):
    """Preprocess the expression data"""
    print("\\nPreprocessing data...")
    
    # Set gene ID as index if specified
    if gene_id_column and gene_id_column in df.columns:
        df = df.set_index(gene_id_column)
        print(f"Set {gene_id_column} as index")
    
    # Identify numeric columns (expression data)
    numeric_columns = df.select_dtypes(include=[np.number]).columns.tolist()
    if condition_column in numeric_columns:
        numeric_columns.remove(condition_column)
    
    print(f"Found {len(numeric_columns)} numeric expression columns")
    
    if remove_na:
        initial_shape = df.shape
        df = df.dropna(subset=numeric_columns)
        print(f"Removed NA values: {initial_shape} -> {df.shape}")
    
    if log_transform:
        print("Applying log2 transformation...")
        # Add small constant to avoid log(0)
        df[numeric_columns] = np.log2(df[numeric_columns] + 1)
    
    return df, numeric_columns

def perform_statistical_test(group1_data, group2_data, test_type):
    """Perform statistical test between two groups"""
    if test_type == 'ttest':
        statistic, pvalue = stats.ttest_ind(group1_data, group2_data, nan_policy='omit')
    elif test_type == 'mannwhitney':
        statistic, pvalue = stats.mannwhitneyu(group1_data, group2_data, alternative='two-sided')
    elif test_type == 'wilcoxon':
        # For paired test, need same number of samples
        min_len = min(len(group1_data), len(group2_data))
        statistic, pvalue = stats.wilcoxon(group1_data[:min_len], group2_data[:min_len])
    else:
        raise ValueError(f"Unknown test type: {test_type}")
    
    return statistic, pvalue

def multiple_testing_correction(pvalues, method):
    """Apply multiple testing correction"""
    pvalues = np.array(pvalues)
    
    if method == 'fdr':
        # Benjamini-Hochberg FDR correction
        sorted_indices = np.argsort(pvalues)
        sorted_pvalues = pvalues[sorted_indices]
        n = len(pvalues)
        
        adjusted_pvalues = np.zeros_like(pvalues)
        for i in range(n):
            adjusted_pvalues[sorted_indices[i]] = min(1.0, sorted_pvalues[i] * n / (i + 1))
        
        # Ensure monotonicity
        for i in range(n-2, -1, -1):
            adjusted_pvalues[sorted_indices[i]] = min(
                adjusted_pvalues[sorted_indices[i]], 
                adjusted_pvalues[sorted_indices[i+1]]
            )
            
    elif method == 'bonferroni':
        adjusted_pvalues = np.minimum(pvalues * len(pvalues), 1.0)
    else:  # method == 'none'
        adjusted_pvalues = pvalues.copy()
    
    return adjusted_pvalues

def analyze_deg():
    """Main DEG analysis function"""
    
    # Parameters
    expression_file = "${params.expression_file}"
    condition_column = "${params.condition_column}"
    group1 = "${params.group1}"
    group2 = "${params.group2}"
    gene_id_column = ${params.gene_id_column ? `"${params.gene_id_column}"` : 'None'}
    separator = ${params.separator ? `"${params.separator}"` : 'None'}
    
    log2fc_threshold = ${options.log2fcThreshold}
    pvalue_threshold = ${options.pvalueThreshold}
    correction_method = "${options.correctionMethod}"
    statistical_test = "${options.statisticalTest}"
    create_plots = ${options.createPlots}
    remove_na = ${options.removeNA}
    log_transform = ${options.logTransform}
    
    print("=== Differential Expression Gene Analysis ===")
    print(f"Expression file: {expression_file}")
    print(f"Condition column: {condition_column}")
    print(f"Groups: {group1} vs {group2}")
    print(f"Statistical test: {statistical_test}")
    print(f"Correction method: {correction_method}")
    print(f"Thresholds: |log2FC| >= {log2fc_threshold}, p < {pvalue_threshold}")
    
    # Load and preprocess data
    df = load_data(expression_file, separator)
    
    if condition_column not in df.columns:
        raise ValueError(f"Condition column '{condition_column}' not found in data")
    
    df, expression_columns = preprocess_data(df, gene_id_column, condition_column, remove_na, log_transform)
    
    # Split data by groups
    group1_mask = df[condition_column] == group1
    group2_mask = df[condition_column] == group2
    
    group1_samples = df[group1_mask]
    group2_samples = df[group2_mask]
    
    print(f"\\nGroup '{group1}': {len(group1_samples)} samples")
    print(f"Group '{group2}': {len(group2_samples)} samples")
    
    if len(group1_samples) == 0 or len(group2_samples) == 0:
        raise ValueError("One or both groups have no samples")
    
    # Perform analysis for each gene/feature
    results = []
    
    print(f"\\nAnalyzing {len(expression_columns)} genes/features...")
    
    for i, gene in enumerate(expression_columns):
        if i % 1000 == 0:
            print(f"Processed {i}/{len(expression_columns)} genes...")
        
        group1_expr = group1_samples[gene].dropna()
        group2_expr = group2_samples[gene].dropna()
        
        if len(group1_expr) < 2 or len(group2_expr) < 2:
            continue
        
        # Calculate fold change
        mean1 = group1_expr.mean()
        mean2 = group2_expr.mean()
        
        if mean1 == 0 and mean2 == 0:
            log2fc = 0
        elif mean1 == 0:
            log2fc = float('inf')
        elif mean2 == 0:
            log2fc = float('-inf')
        else:
            log2fc = np.log2(mean2 / mean1)
        
        # Perform statistical test
        try:
            statistic, pvalue = perform_statistical_test(group1_expr, group2_expr, statistical_test)
        except:
            continue
        
        results.append({
            'gene': gene,
            'mean_group1': mean1,
            'mean_group2': mean2,
            'log2_fold_change': log2fc,
            'statistic': statistic,
            'pvalue': pvalue
        })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    if len(results_df) == 0:
        print("No genes could be analyzed!")
        return
    
    print(f"Successfully analyzed {len(results_df)} genes")
    
    # Apply multiple testing correction
    print(f"\\nApplying {correction_method} correction...")
    results_df['adjusted_pvalue'] = multiple_testing_correction(results_df['pvalue'], correction_method)
    
    # Determine significance
    pval_col = 'adjusted_pvalue' if correction_method != 'none' else 'pvalue'
    results_df['significant'] = (
        (np.abs(results_df['log2_fold_change']) >= log2fc_threshold) & 
        (results_df[pval_col] < pvalue_threshold)
    )
    
    results_df['regulation'] = 'unchanged'
    results_df.loc[
        (results_df['significant']) & (results_df['log2_fold_change'] > 0), 
        'regulation'
    ] = 'upregulated'
    results_df.loc[
        (results_df['significant']) & (results_df['log2_fold_change'] < 0), 
        'regulation'
    ] = 'downregulated'
    
    # Sort by significance
    results_df = results_df.sort_values(pval_col)
    
    # Save results
    print("\\nSaving results...")
    results_df.to_csv('deg_results.csv', index=False)
    
    # Save significant genes
    significant_genes = results_df[results_df['significant']].copy()
    significant_genes.to_csv('significant_genes.csv', index=False)
    
    # Generate summary
    total_genes = len(results_df)
    significant_count = len(significant_genes)
    upregulated_count = len(significant_genes[significant_genes['regulation'] == 'upregulated'])
    downregulated_count = len(significant_genes[significant_genes['regulation'] == 'downregulated'])
    
    summary_text = f\"\"\"DEG Analysis Summary
====================

Dataset: {expression_file}
Groups compared: {group1} vs {group2}
Statistical test: {statistical_test}
Multiple testing correction: {correction_method}
Significance thresholds: |log2FC| >= {log2fc_threshold}, p < {pvalue_threshold}

Results:
- Total genes analyzed: {total_genes}
- Total significant genes: {significant_count}
- Upregulated genes: {upregulated_count}
- Downregulated genes: {downregulated_count}
- Percentage significant: {(significant_count/total_genes*100):.1f}%

Files generated:
- deg_results.csv: Complete results for all genes
- significant_genes.csv: Filtered significant genes only
- deg_summary.txt: This summary
\"\"\"

    with open('deg_summary.txt', 'w') as f:
        f.write(summary_text)
    
    print(summary_text)
    
    # Create plots if requested and matplotlib is available
    if create_plots and PLOTTING_AVAILABLE and len(results_df) > 0:
        print("\\nGenerating plots...")
        
        # Volcano plot
        plt.figure(figsize=(10, 8))
        
        # Plot non-significant points
        non_sig = results_df[~results_df['significant']]
        plt.scatter(non_sig['log2_fold_change'], -np.log10(non_sig[pval_col]), 
                   c='gray', alpha=0.6, s=20, label='Not significant')
        
        # Plot significant points
        sig_up = results_df[results_df['regulation'] == 'upregulated']
        sig_down = results_df[results_df['regulation'] == 'downregulated']
        
        if len(sig_up) > 0:
            plt.scatter(sig_up['log2_fold_change'], -np.log10(sig_up[pval_col]), 
                       c='red', alpha=0.7, s=20, label=f'Upregulated ({len(sig_up)})')
        
        if len(sig_down) > 0:
            plt.scatter(sig_down['log2_fold_change'], -np.log10(sig_down[pval_col]), 
                       c='blue', alpha=0.7, s=20, label=f'Downregulated ({len(sig_down)})')
        
        # Add threshold lines
        plt.axvline(x=log2fc_threshold, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-log2fc_threshold, color='black', linestyle='--', alpha=0.5)
        plt.axhline(y=-np.log10(pvalue_threshold), color='black', linestyle='--', alpha=0.5)
        
        plt.xlabel('Log2 Fold Change')
        plt.ylabel(f'-Log10 {pval_col.replace("_", " ").title()}')
        plt.title(f'Volcano Plot: {group1} vs {group2}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('volcano_plot.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # MA plot
        plt.figure(figsize=(10, 8))
        
        # Calculate average expression
        results_df['average_expression'] = (results_df['mean_group1'] + results_df['mean_group2']) / 2
        
        # Plot non-significant points
        plt.scatter(results_df.loc[~results_df['significant'], 'average_expression'], 
                   results_df.loc[~results_df['significant'], 'log2_fold_change'], 
                   c='gray', alpha=0.6, s=20, label='Not significant')
        
        # Plot significant points
        if len(sig_up) > 0:
            plt.scatter(sig_up['average_expression'], sig_up['log2_fold_change'], 
                       c='red', alpha=0.7, s=20, label=f'Upregulated ({len(sig_up)})')
        
        if len(sig_down) > 0:
            plt.scatter(sig_down['average_expression'], sig_down['log2_fold_change'], 
                       c='blue', alpha=0.7, s=20, label=f'Downregulated ({len(sig_down)})')
        
        # Add threshold lines
        plt.axhline(y=log2fc_threshold, color='black', linestyle='--', alpha=0.5)
        plt.axhline(y=-log2fc_threshold, color='black', linestyle='--', alpha=0.5)
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        plt.xlabel('Average Expression')
        plt.ylabel('Log2 Fold Change')
        plt.title(f'MA Plot: {group1} vs {group2}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('ma_plot.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Plots saved: volcano_plot.png, ma_plot.png")

if __name__ == "__main__":
    try:
        analyze_deg()
        print("\\n=== Analysis completed successfully! ===")
    except Exception as e:
        print(f"\\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
`;
  }
} 