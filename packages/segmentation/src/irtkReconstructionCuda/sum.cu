/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

//! Implementation of array summation algorithms

#ifndef SUM_CU
#define SUM_CU

#include "gaussFilterConvolution.cuh"
#include "sum_kernels.cu"

////////////////////////////////////////////////////////////////////////////////
//! Sums up a 2D array, and writes it to the variable sum.
//! doublehis template was tested for "float" and "int" data types.
//! @param sum            sum of the array in input
//! @param input          pointer to input image
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
////////////////////////////////////////////////////////////////////////////////
//template <class T>
double sum2D(double* sum, double* input, unsigned int width, unsigned int height, unsigned int pitch)
{
  // prepare fragmentation for processing
  dim3 dimBlock(SUM_BLOCK, 1, 1);
  dim3 dimGridX(iDivUp(width, dimBlock.x), 1);

  // allocate memory for row on the host
  double* temp = NULL;
  cudaMalloc((void**) &temp, dimGridX.x*dimBlock.x*sizeof(double));
  CHECK_ERROR(sum2D);

  // Sum up collumns on GPU
  SumXKernel<<<dimGridX, dimBlock>>>(input, temp, width, height, pitch);
  CHECK_ERROR(sum2D);

  // Copy the row to the CPU
  double* temp_host = (double*)malloc(width*sizeof(double));
  cudaMemcpy(temp_host, temp, width*sizeof(double), cudaMemcpyDeviceToHost);
  CHECK_ERROR(sum2D);

  // Calculate final sum on the CPU
  (*sum) = 0;
  for (int i=0; i<width; i++)
    (*sum) += temp_host[i];

  // Clean up memory
  cudaFree(temp);
  free(temp_host);
  CHECK_ERROR(sum2D);

  return 0;
}

float sum2D(float* sum, float* input, unsigned int width, unsigned int height, unsigned int pitch)
{
  // prepare fragmentation for processing
  dim3 dimBlock(SUM_BLOCK, 1, 1);
  dim3 dimGridX(iDivUp(width, dimBlock.x), 1);

  // allocate memory for row on the host
  float* temp = NULL;
  cudaMalloc((void**) &temp, dimGridX.x*dimBlock.x*sizeof(float));
  CHECK_ERROR(sum2D);

  // Sum up collumns on GPU
  SumXKernel<<<dimGridX, dimBlock>>>(input, temp, width, height, pitch);
  CHECK_ERROR(sum2D);

  // Copy the row to the CPU
  float* temp_host = (float*)malloc(width*sizeof(float));
  cudaMemcpy(temp_host, temp, width*sizeof(float), cudaMemcpyDeviceToHost);
  CHECK_ERROR(sum2D);

  // Calculate final sum on the CPU
  (*sum) = 0;
  for (int i=0; i<width; i++)
    (*sum) += temp_host[i];

  // Clean up memory
  cudaFree(temp);
  free(temp_host);
  CHECK_ERROR(sum2D);

  return 0;
}


#endif // SUM_CU
