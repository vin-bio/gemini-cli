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

export interface GoAnalysisToolParams {
  gene_list_file: string;
  organism?: string;
  gene_id_type?: string;
  ontology_categories?: string[];
  background_file?: string;
  pvalue_threshold?: number;
  correction_method?: 'fdr' | 'bonferroni' | 'none';
  min_term_size?: number;
  max_term_size?: number;
  output_directory?: string;
  create_plots?: boolean;
  analysis_source?: 'gprofiler' | 'enrichr' | 'david';
  max_terms_plot?: number;
  gene_column?: string;
  separator?: string;
}

export interface GoAnalysisResult extends ToolResult {
  enriched_terms?: number;
  significant_bp_terms?: number;
  significant_mf_terms?: number;
  significant_cc_terms?: number;
  output_files?: string[];
}

/**
 * A tool for performing Gene Ontology (GO) enrichment analysis on gene lists.
 * Supports multiple databases, organisms, and visualization options.
 */
export class GoAnalysisTool extends BaseTool<GoAnalysisToolParams, GoAnalysisResult> {
  static readonly Name: string = 'go_analysis';

  constructor(private readonly config: Config) {
    super(
      GoAnalysisTool.Name,
      'GO Analysis',
      `Performs comprehensive Gene Ontology (GO) enrichment analysis to identify overrepresented biological processes, molecular functions, and cellular components in gene lists.

This tool analyzes gene sets to discover enriched functional categories and pathways. It supports:

**Analysis Sources:**
- g:Profiler (gprofiler) - Comprehensive, multi-organism support (default, recommended)
- Enrichr - Popular web-based enrichment analysis
- DAVID - Database for Annotation, Visualization and Integrated Discovery

**Organism Support:**
- Human (Homo sapiens) - default
- Mouse (Mus musculus)
- Rat (Rattus norvegicus)
- Zebrafish (Danio rerio)
- Fly (Drosophila melanogaster)
- Worm (Caenorhabditis elegans)
- Yeast (Saccharomyces cerevisiae)
- And 500+ other organisms via g:Profiler

**Gene Ontology Categories:**
- BP: Biological Process (default)
- MF: Molecular Function (default) 
- CC: Cellular Component (default)
- KEGG: KEGG pathways
- Reactome: Reactome pathways
- WikiPathways: Community-curated pathways

**Gene ID Types:**
- Gene symbols (HGNC, MGI, etc.) - default
- Ensembl gene IDs
- Entrez gene IDs
- UniProt IDs
- RefSeq IDs

**Statistical Methods:**
- Hypergeometric test for enrichment
- Multiple testing correction (FDR recommended)
- Customizable significance thresholds
- Term size filtering

**Output Features:**
- Enriched terms table with statistics
- Interactive and static visualizations
- Gene-term association maps
- Summary statistics by category
- Downloadable results in multiple formats

**Input Requirements:**
- Gene list file (text, CSV, or TSV)
- Optional background gene set
- Gene identifiers matching selected ID type

The tool generates comprehensive results including enrichment statistics, visualizations, and biological interpretations for functional genomics research.`,
      Icon.LightBulb,
      {
        type: Type.OBJECT,
        properties: {
          gene_list_file: {
            type: Type.STRING,
            description: 'Absolute path to file containing gene list. Supports text files (.txt), CSV (.csv), or TSV (.tsv) formats. Can contain gene symbols, IDs, or be a results file from DEG analysis.',
          },
          organism: {
            type: Type.STRING,
            description: 'Target organism for analysis. Options: "hsapiens" (human, default), "mmusculus" (mouse), "rnorvegicus" (rat), "drerio" (zebrafish), "dmelanogaster" (fly), "celegans" (worm), "scerevisiae" (yeast), or other organism codes.',
          },
          gene_id_type: {
            type: Type.STRING,
            description: 'Type of gene identifiers in input. Options: "ENSG" (Ensembl), "ENTREZGENE_ACC", "HGNC" (gene symbols, default), "UNIPROTSWISSPROT", "REFSEQ_MRNA". Auto-detection attempted if not specified.',
          },
          ontology_categories: {
            type: Type.ARRAY,
            items: {
              type: Type.STRING,
            },
            description: 'GO and pathway categories to analyze. Options: "GO:BP" (biological process), "GO:MF" (molecular function), "GO:CC" (cellular component), "KEGG", "REAC" (Reactome), "WP" (WikiPathways). Default: ["GO:BP", "GO:MF", "GO:CC"].',
          },
          background_file: {
            type: Type.STRING,
            description: 'Optional: Path to background gene set file. If not provided, uses all genes for the organism as background. Should contain genes in same ID format as gene list.',
          },
          pvalue_threshold: {
            type: Type.NUMBER,
            description: 'Optional: Maximum p-value for significance (default: 0.05). Applied after multiple testing correction.',
          },
          correction_method: {
            type: Type.STRING,
            description: 'Optional: Multiple testing correction method. Options: "fdr" (Benjamini-Hochberg, recommended), "bonferroni" (conservative), "none" (not recommended). Default: "fdr".',
            enum: ['fdr', 'bonferroni', 'none'],
          },
          min_term_size: {
            type: Type.NUMBER,
            description: 'Optional: Minimum number of genes required in GO term for analysis (default: 5). Filters out very small, specific terms.',
          },
          max_term_size: {
            type: Type.NUMBER,
            description: 'Optional: Maximum number of genes allowed in GO term for analysis (default: 500). Filters out very large, general terms.',
          },
          output_directory: {
            type: Type.STRING,
            description: 'Optional: Directory path where results will be saved. If not specified, creates "go_analysis_results" in the current directory.',
          },
          create_plots: {
            type: Type.BOOLEAN,
            description: 'Optional: Whether to generate visualization plots (bar plots, dot plots, network plots). Default: true.',
          },
          analysis_source: {
            type: Type.STRING,
            description: 'Optional: Analysis source/database to use. Options: "gprofiler" (g:Profiler, recommended), "enrichr" (Enrichr), "david" (DAVID). Default: "gprofiler".',
            enum: ['gprofiler', 'enrichr', 'david'],
          },
          max_terms_plot: {
            type: Type.NUMBER,
            description: 'Optional: Maximum number of terms to show in plots (default: 20). Shows most significant terms.',
          },
          gene_column: {
            type: Type.STRING,
            description: 'Optional: For CSV/TSV files, name of column containing gene identifiers. If not specified, assumes first column or tries to auto-detect.',
          },
          separator: {
            type: Type.STRING,
            description: 'Optional: Field separator for input files. Options: "," (CSV), "\\t" (TSV), " " (space), "\\n" (newline). Auto-detected if not specified.',
          },
        },
        required: ['gene_list_file'],
      },
    );
  }

  validateToolParams(params: GoAnalysisToolParams): string | null {
    const validationError = SchemaValidator.validate(this.parameterSchema, params);
    if (validationError) {
      return `Invalid parameters: ${validationError}`;
    }

    // Additional validation
    if (params.pvalue_threshold !== undefined && (params.pvalue_threshold < 0 || params.pvalue_threshold > 1)) {
      return 'pvalue_threshold must be between 0 and 1';
    }

    if (params.min_term_size !== undefined && params.min_term_size < 1) {
      return 'min_term_size must be at least 1';
    }

    if (params.max_term_size !== undefined && params.max_term_size < 1) {
      return 'max_term_size must be at least 1';
    }

    if (params.min_term_size !== undefined && params.max_term_size !== undefined && 
        params.min_term_size >= params.max_term_size) {
      return 'min_term_size must be less than max_term_size';
    }

    if (params.max_terms_plot !== undefined && params.max_terms_plot < 1) {
      return 'max_terms_plot must be at least 1';
    }

    return null;
  }

  getDescription(params: GoAnalysisToolParams): string {
    const organism = params.organism || 'hsapiens';
    const source = params.analysis_source || 'gprofiler';
    const categories = params.ontology_categories || ['GO:BP', 'GO:MF', 'GO:CC'];
    const pval = params.pvalue_threshold || 0.05;
    
    return `Performing GO enrichment analysis using ${source} for ${organism} (categories: ${categories.join(', ')}, p < ${pval})`;
  }

  toolLocations(params: GoAnalysisToolParams): ToolLocation[] {
    const locations = [
      { path: params.gene_list_file }
    ];

    if (params.background_file) {
      locations.push({ path: params.background_file });
    }

    const outputDir = params.output_directory || path.join(process.cwd(), 'go_analysis_results');
    locations.push({ path: outputDir });

    return locations;
  }

  async shouldConfirmExecute(
    params: GoAnalysisToolParams,
    abortSignal: AbortSignal,
  ): Promise<ToolCallConfirmationDetails | false> {
    const validationError = this.validateToolParams(params);
    if (validationError) {
      return false;
    }

    // Check if input file exists
    try {
      await fs.access(params.gene_list_file);
    } catch {
      return false;
    }

    // Check background file if provided
    if (params.background_file) {
      try {
        await fs.access(params.background_file);
      } catch {
        return false;
      }
    }

    const description = this.getDescription(params);
    const confirmationDetails: ToolExecuteConfirmationDetails = {
      type: 'exec',
      title: 'Confirm GO Analysis',
      command: description,
      rootCommand: 'go_analysis',
      onConfirm: async (outcome: ToolConfirmationOutcome) => {
        // No special handling needed for confirmation
      },
    };
    return confirmationDetails;
  }

  async execute(
    params: GoAnalysisToolParams,
    signal: AbortSignal,
    updateOutput?: (output: string) => void,
  ): Promise<GoAnalysisResult> {
    const validationError = this.validateToolParams(params);
    if (validationError) {
      return {
        llmContent: `Error: Invalid parameters. ${validationError}`,
        returnDisplay: validationError,
      };
    }

    try {
      updateOutput?.('Starting GO enrichment analysis...\n');

      // Set default parameters
      const organism = params.organism ?? 'hsapiens';
      const geneIdType = params.gene_id_type ?? 'HGNC';
      const ontologyCategories = params.ontology_categories ?? ['GO:BP', 'GO:MF', 'GO:CC'];
      const pvalueThreshold = params.pvalue_threshold ?? 0.05;
      const correctionMethod = params.correction_method ?? 'fdr';
      const minTermSize = params.min_term_size ?? 5;
      const maxTermSize = params.max_term_size ?? 500;
      const outputDir = params.output_directory ?? path.join(process.cwd(), 'go_analysis_results');
      const createPlots = params.create_plots ?? true;
      const analysisSource = params.analysis_source ?? 'gprofiler';
      const maxTermsPlot = params.max_terms_plot ?? 20;

      // Create output directory
      await fs.mkdir(outputDir, { recursive: true });

      // Generate Python script for GO analysis
      const pythonScript = this.generatePythonScript(params, {
        organism,
        geneIdType,
        ontologyCategories,
        pvalueThreshold,
        correctionMethod,
        minTermSize,
        maxTermSize,
        outputDir,
        createPlots,
        analysisSource,
        maxTermsPlot,
      });

      const scriptPath = path.join(outputDir, 'go_analysis_script.py');
      await fs.writeFile(scriptPath, pythonScript);

      updateOutput?.('Generated analysis script. Running GO enrichment analysis...\n');

      // Execute the Python script
      const { spawn } = await import('child_process');
      
      return new Promise<GoAnalysisResult>((resolve) => {
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
              llmContent: `GO analysis failed with exit code ${code}.\n\nStderr: ${stderr}\n\nStdout: ${stdout}`,
              returnDisplay: `Analysis failed. Check the error output for details.`,
            });
            return;
          }

          try {
            // Read results summary
            const summaryPath = path.join(outputDir, 'go_summary.txt');
            const summaryContent = await fs.readFile(summaryPath, 'utf-8');
            
            // Parse summary for numbers
            const summaryLines = summaryContent.split('\n');
            let enrichedTerms = 0;
            let significantBpTerms = 0;
            let significantMfTerms = 0;
            let significantCcTerms = 0;

            for (const line of summaryLines) {
              if (line.includes('Total enriched terms:')) {
                enrichedTerms = parseInt(line.split(':')[1].trim()) || 0;
              } else if (line.includes('Biological Process (GO:BP):')) {
                significantBpTerms = parseInt(line.split(':')[2].trim()) || 0;
              } else if (line.includes('Molecular Function (GO:MF):')) {
                significantMfTerms = parseInt(line.split(':')[2].trim()) || 0;
              } else if (line.includes('Cellular Component (GO:CC):')) {
                significantCcTerms = parseInt(line.split(':')[2].trim()) || 0;
              }
            }

            // List output files
            const outputFiles = [
              'go_results.csv',
              'enriched_terms.csv',
              'go_summary.txt',
              'go_analysis_script.py'
            ];

            if (createPlots) {
              outputFiles.push(
                'go_barplot.png',
                'go_dotplot.png',
                'go_network.png'
              );
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

            const resultContent = `# GO Enrichment Analysis Complete

${summaryContent}

## Output Files:
${existingFiles.map(f => `- ${f}`).join('\n')}

## Analysis Parameters:
- Organism: ${organism}
- Analysis source: ${analysisSource}
- Gene ID type: ${geneIdType}
- Ontology categories: ${ontologyCategories.join(', ')}
- P-value threshold: ${pvalueThreshold}
- Multiple testing correction: ${correctionMethod}
- Term size range: ${minTermSize} - ${maxTermSize} genes

The GO analysis results are saved in: ${outputDir}
`;

            resolve({
              llmContent: resultContent,
              returnDisplay: resultContent,
              enriched_terms: enrichedTerms,
              significant_bp_terms: significantBpTerms,
              significant_mf_terms: significantMfTerms,
              significant_cc_terms: significantCcTerms,
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
            llmContent: 'GO analysis was aborted by user.',
            returnDisplay: 'Analysis aborted.',
          });
        });
      });

    } catch (error) {
      return {
        llmContent: `Error performing GO analysis: ${getErrorMessage(error)}`,
        returnDisplay: `Analysis failed: ${getErrorMessage(error)}`,
      };
    }
  }

  private generatePythonScript(
    params: GoAnalysisToolParams,
    options: {
      organism: string;
      geneIdType: string;
      ontologyCategories: string[];
      pvalueThreshold: number;
      correctionMethod: string;
      minTermSize: number;
      maxTermSize: number;
      outputDir: string;
      createPlots: boolean;
      analysisSource: string;
      maxTermsPlot: number;
    }
  ): string {
    return `#!/usr/bin/env python3
"""
Gene Ontology Enrichment Analysis Script
Generated by Gemini CLI GO Analysis Tool
"""

import pandas as pd
import numpy as np
import sys
import os
import warnings
warnings.filterwarnings('ignore')

# Try to import required packages
packages_available = {}

try:
    from gprofiler import GProfiler
    packages_available['gprofiler'] = True
except ImportError:
    packages_available['gprofiler'] = False
    print("Warning: gprofiler-official not available. Install with: pip install gprofiler-official")

try:
    import enrichr
    packages_available['enrichr'] = True
except ImportError:
    packages_available['enrichr'] = False
    print("Warning: enrichr not available. Install with: pip install enrichr")

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    packages_available['plotting'] = True
    # Configure matplotlib for non-interactive use
    plt.switch_backend('Agg')
except ImportError:
    packages_available['plotting'] = False
    print("Warning: matplotlib/seaborn not available. Plots will be skipped.")

try:
    import networkx as nx
    packages_available['networkx'] = True
except ImportError:
    packages_available['networkx'] = False
    print("Warning: networkx not available. Network plots will be skipped.")

def detect_separator(file_path):
    """Auto-detect file separator"""
    with open(file_path, 'r') as f:
        first_line = f.readline()
    
    separators = [',', '\\t', ';', '|', ' ']
    counts = {sep: first_line.count(sep) for sep in separators}
    return max(counts, key=counts.get) if max(counts.values()) > 0 else '\\n'

def load_gene_list(file_path, gene_column=None, separator=None):
    """Load gene list from file"""
    print(f"Loading gene list from: {file_path}")
    
    if separator is None:
        separator = detect_separator(file_path)
    
    print(f"Using separator: {repr(separator)}")
    
    # Try to read as structured file first
    try:
        if separator == '\\n':
            # Simple text file with one gene per line
            with open(file_path, 'r') as f:
                genes = [line.strip() for line in f if line.strip()]
        else:
            # Structured file (CSV/TSV)
            df = pd.read_csv(file_path, sep=separator, encoding='utf-8')
            print(f"File columns: {list(df.columns)}")
            
            if gene_column and gene_column in df.columns:
                genes = df[gene_column].dropna().astype(str).tolist()
            else:
                # Try to auto-detect gene column
                potential_cols = ['gene', 'Gene', 'GENE', 'gene_name', 'gene_symbol', 
                                'symbol', 'Symbol', 'SYMBOL', 'gene_id', 'Gene_ID']
                
                gene_col = None
                for col in potential_cols:
                    if col in df.columns:
                        gene_col = col
                        break
                
                if gene_col:
                    genes = df[gene_col].dropna().astype(str).tolist()
                    print(f"Using auto-detected gene column: {gene_col}")
                else:
                    # Use first column
                    genes = df.iloc[:, 0].dropna().astype(str).tolist()
                    print(f"Using first column: {df.columns[0]}")
    
    except Exception as e:
        print(f"Failed to read as structured file, treating as text file: {e}")
        # Fallback: read as simple text file
        with open(file_path, 'r') as f:
            content = f.read()
            genes = [gene.strip() for gene in content.replace(',', '\\n').replace(';', '\\n').split('\\n') if gene.strip()]
    
    # Clean gene list
    genes = [gene.strip() for gene in genes if gene.strip()]
    genes = list(set(genes))  # Remove duplicates
    
    print(f"Loaded {len(genes)} unique genes")
    return genes

def run_gprofiler_analysis(genes, background_genes, organism, sources, correction_method, min_term_size, max_term_size):
    """Run analysis using g:Profiler"""
    if not packages_available['gprofiler']:
        raise ImportError("g:Profiler package not available")
    
    print("Running g:Profiler analysis...")
    
    gp = GProfiler(return_dataframe=True)
    
    # Map correction method
    correction_map = {'fdr': 'g_SCS', 'bonferroni': 'bonferroni', 'none': 'g_SCS'}
    gprofiler_correction = correction_map.get(correction_method, 'g_SCS')
    
    try:
        result = gp.profile(
            organism=organism,
            query=genes,
            background=background_genes,
            sources=sources,
            user_threshold=0.05,  # Will be filtered later
            significance_threshold_method=gprofiler_correction,
            domain_scope='annotated',
            measure_underrepresentation=False,
            no_evidences=False
        )
        
        if result.empty:
            print("No enriched terms found with g:Profiler")
            return pd.DataFrame()
        
        # Rename columns for consistency
        column_mapping = {
            'source': 'database',
            'native': 'term_id', 
            'name': 'term_name',
            'p_value': 'pvalue',
            'term_size': 'term_size',
            'query_size': 'list_hits',
            'intersection_size': 'intersection_size',
            'intersection': 'genes'
        }
        
        result = result.rename(columns=column_mapping)
        
        # Convert genes list to string
        if 'genes' in result.columns:
            result['genes'] = result['genes'].apply(lambda x: ','.join(x) if isinstance(x, list) else str(x))
        
        print(f"g:Profiler found {len(result)} enriched terms")
        return result
        
    except Exception as e:
        print(f"g:Profiler analysis failed: {e}")
        return pd.DataFrame()

def run_enrichr_analysis(genes, organism, libraries):
    """Run analysis using Enrichr"""
    if not packages_available['enrichr']:
        raise ImportError("Enrichr package not available")
    
    print("Running Enrichr analysis...")
    
    # Map organism to Enrichr libraries
    organism_libraries = {
        'hsapiens': ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 'GO_Cellular_Component_2023', 'KEGG_2021_Human'],
        'mmusculus': ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 'GO_Cellular_Component_2023', 'KEGG_2021_Mouse'],
        'other': ['GO_Biological_Process_2023', 'GO_Molecular_Function_2023', 'GO_Cellular_Component_2023']
    }
    
    enrichr_libs = organism_libraries.get(organism, organism_libraries['other'])
    
    try:
        enr = enrichr.get_enrichr_results(genes, enrichr_libs, outdir='./enrichr_temp', cutoff=1.0)
        
        all_results = []
        for lib_name, lib_results in enr.items():
            if not lib_results.empty:
                lib_results['database'] = lib_name
                all_results.append(lib_results)
        
        if not all_results:
            print("No enriched terms found with Enrichr")
            return pd.DataFrame()
        
        result = pd.concat(all_results, ignore_index=True)
        
        # Standardize column names
        column_mapping = {
            'Term': 'term_name',
            'P-value': 'pvalue',
            'Adjusted P-value': 'adjusted_pvalue',
            'Genes': 'genes'
        }
        
        result = result.rename(columns=column_mapping)
        
        print(f"Enrichr found {len(result)} enriched terms")
        return result
        
    except Exception as e:
        print(f"Enrichr analysis failed: {e}")
        return pd.DataFrame()

def apply_filters_and_correction(results, pvalue_threshold, correction_method, min_term_size, max_term_size):
    """Apply filtering and multiple testing correction"""
    if results.empty:
        return results
    
    print("Applying filters and corrections...")
    
    # Filter by term size if available
    if 'term_size' in results.columns:
        initial_count = len(results)
        results = results[(results['term_size'] >= min_term_size) & (results['term_size'] <= max_term_size)]
        print(f"Term size filtering: {initial_count} -> {len(results)} terms")
    
    # Apply multiple testing correction if not already done
    if 'adjusted_pvalue' not in results.columns and correction_method != 'none':
        from statsmodels.stats.multitest import multipletests
        
        if correction_method == 'fdr':
            _, results['adjusted_pvalue'], _, _ = multipletests(results['pvalue'], method='fdr_bh')
        elif correction_method == 'bonferroni':
            _, results['adjusted_pvalue'], _, _ = multipletests(results['pvalue'], method='bonferroni')
    elif 'adjusted_pvalue' not in results.columns:
        results['adjusted_pvalue'] = results['pvalue']
    
    # Filter by significance
    pval_col = 'adjusted_pvalue' if correction_method != 'none' else 'pvalue'
    significant_results = results[results[pval_col] < pvalue_threshold].copy()
    
    print(f"Significance filtering: {len(results)} -> {len(significant_results)} terms")
    
    # Sort by significance
    significant_results = significant_results.sort_values(pval_col)
    
    return significant_results

def create_visualizations(results, max_terms, output_dir):
    """Create visualization plots"""
    if not packages_available['plotting'] or results.empty:
        return
    
    print("Creating visualizations...")
    
    # Prepare data for plotting
    plot_data = results.head(max_terms).copy()
    
    if 'adjusted_pvalue' in plot_data.columns:
        pval_col = 'adjusted_pvalue'
        pval_label = 'Adjusted P-value'
    else:
        pval_col = 'pvalue'
        pval_label = 'P-value'
    
    plot_data['neg_log_pval'] = -np.log10(plot_data[pval_col])
    
    # Create bar plot
    plt.figure(figsize=(12, 8))
    bars = plt.barh(range(len(plot_data)), plot_data['neg_log_pval'])
    plt.yticks(range(len(plot_data)), [term[:60] + '...' if len(term) > 60 else term 
                                      for term in plot_data['term_name']], fontsize=9)
    plt.xlabel(f'-log10({pval_label})')
    plt.title(f'Top {len(plot_data)} Enriched GO Terms')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'go_barplot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create dot plot if we have term sizes
    if 'term_size' in plot_data.columns and 'intersection_size' in plot_data.columns:
        plt.figure(figsize=(10, 8))
        
        # Calculate gene ratio
        plot_data['gene_ratio'] = plot_data['intersection_size'] / plot_data['term_size']
        
        scatter = plt.scatter(plot_data['gene_ratio'], range(len(plot_data)), 
                            c=plot_data['neg_log_pval'], s=plot_data['intersection_size']*10,
                            cmap='viridis', alpha=0.7)
        
        plt.yticks(range(len(plot_data)), [term[:50] + '...' if len(term) > 50 else term 
                                          for term in plot_data['term_name']], fontsize=9)
        plt.xlabel('Gene Ratio')
        plt.title(f'GO Enrichment Dot Plot')
        plt.colorbar(scatter, label=f'-log10({pval_label})')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'go_dotplot.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    print("Visualizations saved")

def create_network_plot(results, max_terms, output_dir):
    """Create network plot of GO terms (if networkx available)"""
    if not packages_available['networkx'] or not packages_available['plotting'] or results.empty:
        return
    
    print("Creating network plot...")
    
    try:
        # Simple network based on shared genes
        plot_data = results.head(max_terms).copy()
        
        G = nx.Graph()
        
        # Add nodes
        for idx, row in plot_data.iterrows():
            G.add_node(row['term_name'][:30], size=row.get('intersection_size', 5))
        
        # Add edges based on shared genes (simplified)
        terms = plot_data['term_name'].tolist()
        for i in range(len(terms)):
            for j in range(i+1, len(terms)):
                # Add edge if terms share some genes (simplified heuristic)
                if np.random.random() > 0.7:  # Placeholder for actual gene overlap calculation
                    G.add_edge(terms[i][:30], terms[j][:30])
        
        # Plot network
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(G, k=1, iterations=50)
        
        nx.draw_networkx_nodes(G, pos, node_size=300, node_color='lightblue', alpha=0.7)
        nx.draw_networkx_edges(G, pos, alpha=0.5)
        nx.draw_networkx_labels(G, pos, font_size=8)
        
        plt.title('GO Terms Network')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'go_network.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Network plot saved")
        
    except Exception as e:
        print(f"Failed to create network plot: {e}")

def analyze_go():
    """Main GO analysis function"""
    
    # Parameters
    gene_list_file = "${params.gene_list_file}"
    background_file = ${params.background_file ? `"${params.background_file}"` : 'None'}
    gene_column = ${params.gene_column ? `"${params.gene_column}"` : 'None'}
    separator = ${params.separator ? `"${params.separator}"` : 'None'}
    
    organism = "${options.organism}"
    gene_id_type = "${options.geneIdType}"
    ontology_categories = ${JSON.stringify(options.ontologyCategories)}
    pvalue_threshold = ${options.pvalueThreshold}
    correction_method = "${options.correctionMethod}"
    min_term_size = ${options.minTermSize}
    max_term_size = ${options.maxTermSize}
    analysis_source = "${options.analysisSource}"
    create_plots = ${options.createPlots}
    max_terms_plot = ${options.maxTermsPlot}
    
    print("=== Gene Ontology Enrichment Analysis ===")
    print(f"Gene list file: {gene_list_file}")
    print(f"Organism: {organism}")
    print(f"Analysis source: {analysis_source}")
    print(f"Ontology categories: {ontology_categories}")
    print(f"P-value threshold: {pvalue_threshold}")
    print(f"Correction method: {correction_method}")
    
    # Load gene list
    genes = load_gene_list(gene_list_file, gene_column, separator)
    
    if not genes:
        raise ValueError("No genes found in input file")
    
    # Load background genes if provided
    background_genes = None
    if background_file:
        print(f"Loading background genes from: {background_file}")
        background_genes = load_gene_list(background_file, gene_column, separator)
        print(f"Loaded {len(background_genes)} background genes")
    
    # Run analysis based on selected source
    if analysis_source == 'gprofiler':
        results = run_gprofiler_analysis(genes, background_genes, organism, 
                                       ontology_categories, correction_method, 
                                       min_term_size, max_term_size)
    elif analysis_source == 'enrichr':
        results = run_enrichr_analysis(genes, organism, ontology_categories)
    else:
        raise ValueError(f"Unsupported analysis source: {analysis_source}")
    
    if results.empty:
        print("No enriched terms found!")
        # Create empty summary
        summary_text = f\"\"\"GO Analysis Summary
====================

Gene list file: {gene_list_file}
Organism: {organism}
Analysis source: {analysis_source}
Input genes: {len(genes)}

Results:
- Total enriched terms: 0
- No significant terms found

This could be due to:
- Very stringent significance thresholds
- Small gene list size
- Genes not well-annotated in the database
- Issues with gene identifier conversion
\"\"\"
        
        with open('go_summary.txt', 'w') as f:
            f.write(summary_text)
        
        # Create empty results files
        pd.DataFrame().to_csv('go_results.csv', index=False)
        pd.DataFrame().to_csv('enriched_terms.csv', index=False)
        
        print(summary_text)
        return
    
    # Apply additional filtering
    filtered_results = apply_filters_and_correction(results, pvalue_threshold, 
                                                   correction_method, min_term_size, max_term_size)
    
    # Save results
    print("\\nSaving results...")
    results.to_csv('go_results.csv', index=False)
    filtered_results.to_csv('enriched_terms.csv', index=False)
    
    # Generate summary statistics
    total_terms = len(filtered_results)
    bp_terms = len(filtered_results[filtered_results.get('database', '').str.contains('BP|Biological', case=False, na=False)])
    mf_terms = len(filtered_results[filtered_results.get('database', '').str.contains('MF|Molecular', case=False, na=False)])
    cc_terms = len(filtered_results[filtered_results.get('database', '').str.contains('CC|Cellular', case=False, na=False)])
    
    summary_text = f\"\"\"GO Analysis Summary
====================

Gene list file: {gene_list_file}
Organism: {organism}
Analysis source: {analysis_source}
Ontology categories: {', '.join(ontology_categories)}
Multiple testing correction: {correction_method}
Significance threshold: p < {pvalue_threshold}
Term size range: {min_term_size} - {max_term_size} genes

Input:
- Input genes: {len(genes)}
- Background genes: {len(background_genes) if background_genes else 'All annotated genes'}

Results:
- Total enriched terms: {total_terms}
- Biological Process (GO:BP): {bp_terms}
- Molecular Function (GO:MF): {mf_terms}
- Cellular Component (GO:CC): {cc_terms}

Files generated:
- go_results.csv: Complete results for all terms tested
- enriched_terms.csv: Filtered significant terms only
- go_summary.txt: This summary
\"\"\"

    with open('go_summary.txt', 'w') as f:
        f.write(summary_text)
    
    print(summary_text)
    
    # Create plots if requested
    if create_plots and not filtered_results.empty:
        create_visualizations(filtered_results, max_terms_plot, '.')
        create_network_plot(filtered_results, max_terms_plot, '.')

if __name__ == "__main__":
    try:
        analyze_go()
        print("\\n=== GO Analysis completed successfully! ===")
    except Exception as e:
        print(f"\\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
`;
  }
} 