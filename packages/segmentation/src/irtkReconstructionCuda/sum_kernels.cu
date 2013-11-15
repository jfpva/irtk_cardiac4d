/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

//! Implementation of array summation kernels

#ifndef SUM_KERNEL_CU
#define SUM_KERNEL_CU


////////////////////////////////////////////////////////////////////////////////
//! Summation of collumns in an 2D image
//! doublehis template was tested for "float" and "int" data types.
//! @param input          pointer to input image
//! @param output         pointer to output row
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
////////////////////////////////////////////////////////////////////////////////
//template <class T>
__global__ void SumXKernel(double* input, double* output, int width, int height, int pitch)
{
  // determine starting point for single threads
  int start = IMUL(blockIdx.x, blockDim.x) + threadIdx.x;

  // sum up the columns
  if (start < width)
  {
    double sum = 0;
    for (int x=start; x<start+height*pitch; x+=pitch)
      sum += input[x];

    // write back to global memory
    output[start] = sum;
  }
}

__global__ void SumXKernel(float* input, float* output, int width, int height, int pitch)
{
  // determine starting point for single threads
  int start = IMUL(blockIdx.x, blockDim.x) + threadIdx.x;

  // sum up the columns
  if (start < width)
  {
    double sum = 0;
    for (int x=start; x<start+height*pitch; x+=pitch)
      sum += input[x];

    // write back to global memory
    output[start] = sum;
  }
}


#endif // SUM_KERNEL_CU
