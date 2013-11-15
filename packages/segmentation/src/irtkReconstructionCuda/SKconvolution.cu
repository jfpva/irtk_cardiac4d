/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

//! Convolution implementations for separable kernels

#ifndef SK_CONVOLUTION_CU
#define SK_CONVOLUTION_CU

#include <stdio.h>

#include "gaussFilterConvolution.cuh"
#include "SKconvolution_kernels.cu"

////////////////////////////////////////////////////////////////////////////////
//! Performes separable Kernel convolution in 2D
//! Input and output images have to be of same size! Make sure that output memory
//! is already allocated. The input image is automatically padded in constant style.
//! @param input          pointer to input image (device memory)
//! @param output         pointer to output image (device memory, can be the same as input)
//! @param kernel         pointer to the 1D convolution kernel (host memory)
//! @param klength        kernel length [pixels]
//! @param width          width of the image [pixels]
//! @param height         height of the image [pixels]
//! @param pitch          image pitch [pixels]
//! @param num_ch         number of channels in the image
//! @param temp           optional parameter to set own temporary variable.
//!                       has to be of size pitch*height. If not set, memory
//!                       will be allocated.
////////////////////////////////////////////////////////////////////////////////
int SKConvolution(double* input, double* output, double* kernel, int klength,
                  unsigned int width, unsigned int height, unsigned int pitch,
                  unsigned int num_ch, double* temp)
{
  // prepare fragmentation for processing
  dim3 dimBlockX(BLOCK_SIZE_SK_1, BLOCK_SIZE_SK_2, 1);
  dim3 dimGridX(iDivUp(width, dimBlockX.x), iDivUp(height, dimBlockX.y));
  dim3 dimBlockY(BLOCK_SIZE_SK_2, BLOCK_SIZE_SK_1, 1);
  dim3 dimGridY(iDivUp(width, dimBlockY.x), iDivUp(height, dimBlockY.y));
  size_t sharedMem = BLOCK_SIZE_SK_2*(BLOCK_SIZE_SK_1+klength-1)*sizeof(double);

  if (klength > MAX_LENGTH_SK)
  {
    fprintf(stderr, "Error in CUDA SKConvolution(): Kernel size too big\n");
    return -1;
  }

  // If no temporary variable was given to the function, allocate necassary memory
  bool allocated_temp = false;
  if (temp == NULL)
  {
    size_t pt;
    cudaMallocPitch((void**) &temp, &pt, width*sizeof(double), height);
    CHECK_ERROR(SKConvolution);
    allocated_temp = true;
  }

  // Copy kernel to constant memory
  cudaMemcpyToSymbol(small_kernel_const, kernel, sizeof(double)*klength, 0);
  CHECK_ERROR(SKConvolution);

  // Perform separated convolution
  for (unsigned int i=0; i<num_ch; i++)
  {
    SKConvolutionXKernel<<<dimGridX, dimBlockX, sharedMem>>>
        (&input[i*pitch*height], temp, klength, width, height, pitch);
    CHECK_ERROR(SKConvolution);
    SKConvolutionYKernel<<<dimGridY, dimBlockY, sharedMem>>>
        (temp, &output[i*pitch*height], klength, width, height, pitch);
    CHECK_ERROR(SKConvolution);
  }

  // if we hold memory, clean up
  if (allocated_temp)
  {
    cudaFree(temp);
    CHECK_ERROR(SKConvolution);
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! Performes separable Kernel convolution in 3D
//! Input and output images have to be of same size! Make sure that output memory
//! is already allocated. The input image is automatically padded in constant style.
//! @param input          pointer to input image (device memory)
//! @param output         pointer to output image (device memory, can be the same as input)
//! @param kernel         pointer to the 1D convolution kernel (host memory)
//! @param klength        kernel length [pixels]
//! @param width          width of the image [pixels]
//! @param height         height of the image [pixels]
//! @param depth          depth of the image [pixels]
//! @param pitch_w        image pitch (width) [pixels]
//! @param pitch_h        image pitch (height) [pixels]
//! @param num_ch         number of channels in the image
//! @param temp           optional parameter to set own temporary variable.
//!                       has to be of size pitch_w*pitch_h*depth. If not set, memory
//!                       will be allocated.
////////////////////////////////////////////////////////////////////////////////
int SKConvolution3D(double* input, double* output, double* kernel, int klength,
                    unsigned int width, unsigned int height, unsigned int depth,
                    unsigned int pitch_w, unsigned int pitch_h,
                    unsigned int num_ch, double* temp, double* temp2)
{
  // prepare fragmentation for processing
  dim3 dimBlockX(BLOCK_SIZE_SK_1, BLOCK_SIZE_SK_2, 1);
  dim3 dimGridX(iDivUp(width, dimBlockX.x), iDivUp(height, dimBlockX.y));
  dim3 dimBlockY(BLOCK_SIZE_SK_2, BLOCK_SIZE_SK_1, 1);
  dim3 dimGridY(iDivUp(width, dimBlockY.x), iDivUp(height, dimBlockY.y));
  dim3 dimBlockZ(BLOCK_SIZE_SK_2, BLOCK_SIZE_SK_1, 1);
  dim3 dimGridZ(iDivUp(width, dimBlockY.x), iDivUp(height, dimBlockY.y));

  size_t sharedMem = BLOCK_SIZE_SK_2*(BLOCK_SIZE_SK_1+klength-1)*sizeof(double);

  if (klength > MAX_LENGTH_SK)
  {
    printf("Error in CUDA SKConvolution3D(): Kernel size too big\n");
    //return -1;
  }

  // If no temporary variable was given to the function, allocate necassary memory
  bool allocated_temp = false;
  if (temp == NULL)
  {
    size_t pt;
    cudaMallocPitch((void**) &temp, &pt, pitch_w*sizeof(double), pitch_h*depth);
    CHECK_ERROR(SKConvolution);

    cudaMallocPitch((void**) &temp2, &pt, pitch_w*sizeof(double), pitch_h*depth);
    CHECK_ERROR(SKConvolution);
    allocated_temp = true;
  }

  // Copy kernel to constant memory
  cudaMemcpyToSymbol(small_kernel_const, kernel, sizeof(double)*klength, 0);
  CHECK_ERROR(SKConvolution);

  // Perform separated convolution
  for (unsigned int i=0; i<num_ch; i++)
  {
    for (unsigned int z=0; z<depth; z++)
    {
      SKConvolutionXKernel<<<dimGridX, dimBlockX, sharedMem>>>
          (&input[i*pitch_w*pitch_h*depth+z*pitch_w*pitch_h], &temp[z*pitch_w*pitch_h], klength, width, height, pitch_w);
      CHECK_ERROR(SKConvolution);
    }
    for (unsigned int z=0; z<depth; z++)
    {
      SKConvolutionYKernel<<<dimGridY, dimBlockY, sharedMem>>>
          (&temp[z*pitch_w*pitch_h], &temp2[z*pitch_w*pitch_h], klength, width, height, pitch_w);
      CHECK_ERROR(SKConvolution);
    }
    for (unsigned int y=0; y<height; y++)
    {
      SKConvolutionZKernel<<<dimGridZ, dimBlockZ, sharedMem>>>
          (temp2, &output[i*pitch_w*pitch_h*depth], klength, width, height, depth, pitch_w, pitch_h, y);
          CHECK_ERROR(SKConvolution);
    }
  }

  // if we hold memory, clean up
  if (allocated_temp)
  {
    cudaFree(temp);
    cudaFree(temp2);
    CHECK_ERROR(SKConvolution);
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! Convolves image with 1D kernel in row direction
//! Does not work in place !!! Make sure output != input
//! @param input          pointer to input image
//! @param output         pointer to output image (can be the same as input)
//! @param kernel         pointer to the 1D convolution kernels.
//!                       is assumed to be in host memory.
//! @param klength        kernel length
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
//! @param num_ch         number of channels in the image
////////////////////////////////////////////////////////////////////////////////
int SKConvolutionX(double* input, double* output, double* kernel, int klength,
                   unsigned int width, unsigned int height, unsigned int pitch,
                   unsigned int num_ch)
{
  // prepare fragmentation for processing
  dim3 dimBlockX(BLOCK_SIZE_SK_1, BLOCK_SIZE_SK_2, 1);
  dim3 dimGridX(iDivUp(width, dimBlockX.x), iDivUp(height, dimBlockX.y));
  size_t sharedMem = BLOCK_SIZE_SK_2*(BLOCK_SIZE_SK_1+klength-1)*sizeof(double);

  if (klength > MAX_LENGTH_SK)
  {
    fprintf(stderr, "Error in CUDA SKConvolution(): Kernel size too big\n");
    return -1;
  }

  // Copy kernel to constant memory
  cudaMemcpyToSymbol(small_kernel_const, kernel, sizeof(double)*klength, 0);
  CHECK_ERROR(SKConvolutionX);

  // Perform separated convolution
  for (unsigned int i=0; i<num_ch; i++)
  {
    SKConvolutionXKernel<<<dimGridX, dimBlockX, sharedMem>>>
        (&input[i*pitch*height], &output[i*pitch*height], klength, width, height, pitch);
    CHECK_ERROR(SKConvolutionX);
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! Convolves image with 1D kernel in column direction
//! Does not work in place !!! Make sure output != input
//! @param input          pointer to input image
//! @param output         pointer to output image (can be the same as input)
//! @param kernel         pointer to the 1D convolution kernels.
//!                       is assumed to be in host memory.
//! @param klength        kernel length
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
//! @param num_ch         number of channels in the image
////////////////////////////////////////////////////////////////////////////////
int SKConvolutionY(double* input, double* output, double* kernel, int klength,
                   unsigned int width, unsigned int height, unsigned int pitch,
                   unsigned int num_ch)
{
  // prepare fragmentation for processing
  dim3 dimBlockY(BLOCK_SIZE_SK_2, BLOCK_SIZE_SK_1, 1);
  dim3 dimGridY(iDivUp(width, dimBlockY.x), iDivUp(height, dimBlockY.y));
  size_t sharedMem = BLOCK_SIZE_SK_2*(BLOCK_SIZE_SK_1+klength-1)*sizeof(double);

  if (klength > MAX_LENGTH_SK)
  {
    fprintf(stderr, "Error in CUDA SKConvolution(): Kernel size too big\n");
    return -1;
  }

  // Copy kernel to constant memory
  cudaMemcpyToSymbol(small_kernel_const, kernel, sizeof(double)*klength, 0);
  CHECK_ERROR(SKConvolutionY);

  // Perform separated convolution
  for (unsigned int i=0; i<num_ch; i++)
  {
     SKConvolutionYKernel<<<dimGridY, dimBlockY, sharedMem>>>
         (&input[i*pitch*height], &output[i*pitch*height], klength, width, height, pitch);
     CHECK_ERROR(SKConvolutionY);
  }

  return 0;
}

#endif // SK_CONVOLUTION_CU

