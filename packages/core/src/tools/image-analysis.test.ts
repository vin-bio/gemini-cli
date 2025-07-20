/**
 * @license
 * Copyright 2025 Google LLC
 * SPDX-License-Identifier: Apache-2.0
 */

import { describe, it, expect, vi, beforeEach } from 'vitest';
import { ImageAnalysisTool, ImageAnalysisToolParams } from './image-analysis.js';
import { Config } from '../config/config.js';
import { GeminiClient } from '../core/client.js';
import { promises as fs } from 'node:fs';

// Mock dependencies
vi.mock('node:fs', () => ({
  promises: {
    readFile: vi.fn(),
    access: vi.fn(),
  },
}));

describe('ImageAnalysisTool', () => {
  let tool: ImageAnalysisTool;
  let mockConfig: Config;
  let mockGeminiClient: GeminiClient;

  beforeEach(() => {
    vi.clearAllMocks();

    mockGeminiClient = {
      generateContent: vi.fn(),
    } as unknown as GeminiClient;

    mockConfig = {
      getGeminiClient: vi.fn().mockReturnValue(mockGeminiClient),
      getWorkingDir: vi.fn().mockReturnValue('/test/dir'),
    } as unknown as Config;

    tool = new ImageAnalysisTool(mockConfig);
  });

  describe('validateToolParams', () => {
    it('should return null for valid parameters', () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['/test/dir/image.jpg'],
        query: 'What do you see in this image?',
      };

      const result = tool.validateToolParams(params);
      expect(result).toBeNull();
    });

    it('should return error for empty image_paths', () => {
      const params: ImageAnalysisToolParams = {
        image_paths: [],
        query: 'What do you see?',
      };

      const result = tool.validateToolParams(params);
      expect(result).toContain('At least one image path must be provided');
    });

    it('should return error for empty query', () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['/test/dir/image.jpg'],
        query: '',
      };

      const result = tool.validateToolParams(params);
      expect(result).toContain('query');
    });

    it('should return error for too many images', () => {
      const params: ImageAnalysisToolParams = {
        image_paths: Array(11).fill('/test/dir/image.jpg'),
        query: 'What do you see?',
      };

      const result = tool.validateToolParams(params);
      expect(result).toContain('Too many images');
    });

    it('should return error for unsupported image format', () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['/test/dir/image.txt'],
        query: 'What do you see?',
      };

      const result = tool.validateToolParams(params);
      expect(result).toContain('Unsupported image format');
    });

    it('should return error for non-absolute paths', () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['relative/path/image.jpg'],
        query: 'What do you see?',
      };

      const result = tool.validateToolParams(params);
      expect(result).toContain('Image path must be absolute');
    });
  });

  describe('getDescription', () => {
    it('should return correct description for single image', () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['/test/dir/photo.png'],
        query: 'Describe this image',
      };

      const description = tool.getDescription(params);
      expect(description).toContain('Analyzing 1 image');
      expect(description).toContain('photo.png');
      expect(description).toContain('Describe this image');
    });

    it('should return correct description for multiple images', () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['/test/dir/photo1.jpg', '/test/dir/photo2.png'],
        query: 'Compare these images',
      };

      const description = tool.getDescription(params);
      expect(description).toContain('Analyzing 2 images');
      expect(description).toContain('photo1.jpg');
      expect(description).toContain('photo2.png');
      expect(description).toContain('Compare these images');
    });
  });

  describe('execute', () => {
    it('should return error for invalid parameters', async () => {
      const params: ImageAnalysisToolParams = {
        image_paths: [],
        query: 'Test query',
      };

      const result = await tool.execute(params, new AbortController().signal);
      expect(result.llmContent).toContain('Error: Invalid parameters');
    });

    it('should return error for missing image files', async () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['/test/dir/missing.jpg'],
        query: 'What do you see?',
      };

      // Mock fs.access to throw (file doesn't exist)
      vi.mocked(fs.access).mockRejectedValue(new Error('File not found'));

      const result = await tool.execute(params, new AbortController().signal);
      expect(result.llmContent).toContain('Error: Image files not found');
    });

    it('should successfully analyze images', async () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['/test/dir/test.jpg'],
        query: 'What do you see?',
      };

      // Mock file system operations
      vi.mocked(fs.access).mockResolvedValue(undefined);
      vi.mocked(fs.readFile).mockResolvedValue(Buffer.from('fake image data'));

      // Mock Gemini API response
      const mockResponse = {
        candidates: [{
          content: {
            parts: [{ text: 'I can see a beautiful landscape with mountains and trees.' }]
          }
        }]
      } as any;
      vi.mocked(mockGeminiClient.generateContent).mockResolvedValue(mockResponse);

      const result = await tool.execute(params, new AbortController().signal);
      
      expect(result.llmContent).toContain('Image analysis results');
      expect(result.llmContent).toContain('beautiful landscape');
      expect(result.returnDisplay).toContain('Analysis completed');
    });

    it('should handle API errors gracefully', async () => {
      const params: ImageAnalysisToolParams = {
        image_paths: ['/test/dir/test.jpg'],
        query: 'What do you see?',
      };

      // Mock file system operations
      vi.mocked(fs.access).mockResolvedValue(undefined);
      vi.mocked(fs.readFile).mockResolvedValue(Buffer.from('fake image data'));

      // Mock API error
      vi.mocked(mockGeminiClient.generateContent).mockRejectedValue(new Error('API Error'));

      const result = await tool.execute(params, new AbortController().signal);
      
      expect(result.llmContent).toContain('Error during image analysis');
      expect(result.returnDisplay).toContain('Error performing image analysis');
    });
  });
}); 