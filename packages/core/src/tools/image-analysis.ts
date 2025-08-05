/**
 * @license
 * Copyright 2025 Google LLC
 * SPDX-License-Identifier: Apache-2.0
 */

import { promises as fs } from 'node:fs';
import path from 'node:path';
import { BaseTool, Icon, ToolResult } from './tools.js';
import { Type } from '@google/genai';
import { SchemaValidator } from '../utils/schemaValidator.js';
import { Config } from '../config/config.js';
import { getErrorMessage } from '../utils/errors.js';
import { getResponseText } from '../utils/generateContentResponseUtilities.js';
import { isWithinRoot } from '../utils/fileUtils.js';

/**
 * Parameters for the ImageAnalysisTool.
 */
export interface ImageAnalysisToolParams {
  /**
   * Array of absolute paths to image files to analyze
   */
  image_paths: string[];

  /**
   * The text query or prompt for analyzing the images
   */
  query: string;
}

/**
 * A tool to analyze images using Gemini's vision capabilities.
 */
export class ImageAnalysisTool extends BaseTool<
  ImageAnalysisToolParams,
  ToolResult
> {
  static readonly Name: string = 'analyze_images';

  // Supported image formats by Gemini
  private static readonly SUPPORTED_FORMATS = [
    '.jpg', '.jpeg', '.png', '.gif', '.bmp', '.webp','.tiff'
  ];

  constructor(private readonly config: Config) {
    super(
      ImageAnalysisTool.Name,
      'Image Analysis',
      'Analyzes images using Gemini\'s vision capabilities. Provide image file paths and a text query describing what you want to analyze or understand about the images.',
      Icon.LightBulb,
      {
        type: Type.OBJECT,
        properties: {
          image_paths: {
            type: Type.ARRAY,
            items: {
              type: Type.STRING,
            },
            description: 'Array of absolute paths to image files to analyze. Supports common formats like JPG, PNG, GIF, BMP, WEBP.',
          },
          query: {
            type: Type.STRING,
            description: 'The text query or prompt describing what you want to analyze about the images.',
          },
        },
        required: ['image_paths', 'query'],
      },
    );
  }

  /**
   * Validates the parameters for the ImageAnalysisTool.
   * @param params The parameters to validate
   * @returns An error message string if validation fails, null if valid
   */
  validateToolParams(params: ImageAnalysisToolParams): string | null {
    const errors = SchemaValidator.validate(this.schema.parameters, params);
    if (errors) {
      return errors;
    }

    if (!params.image_paths || params.image_paths.length === 0) {
      return "At least one image path must be provided in 'image_paths'.";
    }

    if (params.image_paths.length > 10) {
      return "Too many images provided. Maximum of 10 images allowed.";
    }

    if (!params.query || params.query.trim() === '') {
      return "The 'query' parameter cannot be empty.";
    }

    // Validate file paths
    for (const imagePath of params.image_paths) {
      if (!path.isAbsolute(imagePath)) {
        return `Image path must be absolute: ${imagePath}`;
      }

      if (!isWithinRoot(imagePath, this.config.getWorkingDir())) {
        return `Image path is outside the working directory: ${imagePath}`;
      }

      const ext = path.extname(imagePath).toLowerCase();
      if (!ImageAnalysisTool.SUPPORTED_FORMATS.includes(ext)) {
        return `Unsupported image format: ${ext}. Supported formats: ${ImageAnalysisTool.SUPPORTED_FORMATS.join(', ')}`;
      }
    }

    return null;
  }

  getDescription(params: ImageAnalysisToolParams): string {
    const imageCount = params.image_paths.length;
    const imageNames = params.image_paths.map(p => path.basename(p)).join(', ');
    return `Analyzing ${imageCount} image${imageCount > 1 ? 's' : ''} (${imageNames}) with query: "${params.query}"`;
  }

  /**
   * Reads an image file and converts it to base64 for Gemini API
   */
  private async readImageAsBase64(imagePath: string): Promise<{
    data: string;
    mimeType: string;
  }> {
    try {
      const imageBuffer = await fs.readFile(imagePath);
      const base64Data = imageBuffer.toString('base64');
      
      // Determine MIME type based on file extension
      const ext = path.extname(imagePath).toLowerCase();
      let mimeType: string;
      
      switch (ext) {
        case '.jpg':
        case '.jpeg':
          mimeType = 'image/jpeg';
          break;
        case '.png':
          mimeType = 'image/png';
          break;
        case '.gif':
          mimeType = 'image/gif';
          break;
        case '.bmp':
          mimeType = 'image/bmp';
          break;
        case '.webp':
          mimeType = 'image/webp';
          break;
        default:
          mimeType = 'image/jpeg'; // fallback
      }

      return {
        data: base64Data,
        mimeType,
      };
    } catch (error) {
      throw new Error(`Failed to read image file ${imagePath}: ${getErrorMessage(error)}`);
    }
  }

  async execute(
    params: ImageAnalysisToolParams,
    signal: AbortSignal,
  ): Promise<ToolResult> {
    const validationError = this.validateToolParams(params);
    if (validationError) {
      return {
        llmContent: `Error: Invalid parameters provided. Reason: ${validationError}`,
        returnDisplay: validationError,
      };
    }

    const geminiClient = this.config.getGeminiClient();

    try {
      // Check if all image files exist
      const missingFiles: string[] = [];
      for (const imagePath of params.image_paths) {
        try {
          await fs.access(imagePath);
        } catch {
          missingFiles.push(imagePath);
        }
      }

      if (missingFiles.length > 0) {
        const errorMsg = `Image files not found: ${missingFiles.join(', ')}`;
        return {
          llmContent: `Error: ${errorMsg}`,
          returnDisplay: errorMsg,
        };
      }

      // Read all images and convert to base64
      const imagePromises = params.image_paths.map(imagePath =>
        this.readImageAsBase64(imagePath)
      );
      
      const images = await Promise.all(imagePromises);

      // Create parts array with text query and images
      const parts = [
        { text: params.query },
        ...images.map(img => ({
          inlineData: {
            data: img.data,
            mimeType: img.mimeType,
          }
        }))
      ];

      // Call Gemini API with vision capabilities
      const response = await geminiClient.generateContent(
        [{ role: 'user', parts }],
        {
          temperature: 0.1, // Lower temperature for more consistent analysis
        },
        signal,
      );

      const responseText = getResponseText(response);

      if (!responseText || !responseText.trim()) {
        return {
          llmContent: `No analysis results returned for images: ${params.image_paths.map(p => path.basename(p)).join(', ')}`,
          returnDisplay: 'No analysis results found.',
        };
      }

      const imageList = params.image_paths.map((p, i) => `${i + 1}. ${path.basename(p)} (${p})`).join('\n');

      return {
        llmContent: `Image analysis results for query: "${params.query}"\n\nImages analyzed:\n${imageList}\n\nAnalysis:\n${responseText}`,
        returnDisplay: `Analysis completed for ${params.image_paths.length} image${params.image_paths.length > 1 ? 's' : ''}.`,
      };

    } catch (error: unknown) {
      const errorMessage = `Error during image analysis: ${getErrorMessage(error)}`;
      console.error(errorMessage, error);
      return {
        llmContent: `Error: ${errorMessage}`,
        returnDisplay: `Error performing image analysis.`,
      };
    }
  }
} 