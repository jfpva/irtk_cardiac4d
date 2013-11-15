/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

//! Different implementations of separable kernel convolution

#ifndef SK_CONVOLUTION_KERNELS_CU
#define SK_CONVOLUTION_KERNELS_CU

// Check if BLOCK SIZE was set correctly
#if BLOCK_SIZE_SK_1 < BLOCK_SIZE_SK_2
#error "Configuration Error: Make sure that BLOCK_SIZE_SK_1 >= BLOCK_SIZE_SK_2"
#endif

// Determine maximum kernel length for separable kernel convolution
#define MAX_LENGTH_SK BLOCK_SIZE_SK_1-1

//!//////////////////////////////////////////////////////////////////////////////
//! Separable Kernel Convolution
//!//////////////////////////////////////////////////////////////////////////////

// Device constant array to store kernel coefficients
__constant__  double  small_kernel_const[MAX_LENGTH_SK];

////////////////////////////////////////////////////////////////////////////////
//! Perform a separable kernel convolution with an 1D Kernel in x direction
//! @param input          pointer to input image
//! @param output         pointer to output image
//! @param klength        kernel length
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
//! Make sure that the kernel length "klength" is smaller than "blockDim.x"
////////////////////////////////////////////////////////////////////////////////
__global__ void SKConvolutionXKernel(double* input, double* output, int klength,
                                     int width, int height, int pitch)
{
  // calculate absolute coordinates
  int x = IMUL(blockIdx.x, blockDim.x) + threadIdx.x;
  int y = IMUL(blockIdx.y, blockDim.y) + threadIdx.y;

  int hkl = (klength-1)/2;

  int center = IMUL(y,pitch)+x;
  int shared_width = blockDim.x+klength-1;
  int Scenter = IMUL(threadIdx.y,shared_width)+threadIdx.x+hkl;

  extern __shared__ double input_shared_SKconv[];

  if (y<height)
  {
    // Load data into shared memory
    input_shared_SKconv[Scenter] = input[min(y,height-1)*pitch+min(x,width-1)];

    if (threadIdx.x<hkl)
    {
      // load on left and right side
      input_shared_SKconv[Scenter-hkl] = input[IMUL(y, pitch)+max(0, min(width-1, x-hkl))];
      input_shared_SKconv[Scenter+blockDim.x] = input[IMUL(y, pitch)+max(0, min(width-1, x+blockDim.x))];
    }
  }

  __syncthreads();

  // do convolution
  double value = 0;
  if (x<width && y<height)
  {
    for (int i=-hkl; i<=hkl; i++)
      value += input_shared_SKconv[Scenter+i]*small_kernel_const[i+hkl];
    output[center] = value;
	// if(output[center] != output[center])
		 // output[center] = 0;
		//  printf("NAN ");
  }


}

////////////////////////////////////////////////////////////////////////////////
//! Perform a separable kernel convolution with an 1D Kernel in y direction
//! @param input          pointer to input image
//! @param output         pointer to output image
//! @param klength        kernel length
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
//! Make sure that the kernel length "klength" is smaller than "blockDim.y"
////////////////////////////////////////////////////////////////////////////////
__global__ void SKConvolutionYKernel(double* input, double* output, int klength,
                                     int width, int height, int pitch)
{
  // calculate absolute coordinates
  int x = IMUL(blockIdx.x, blockDim.x) + threadIdx.x;
  int y = IMUL(blockIdx.y, blockDim.y) + threadIdx.y;

  int hkl = (klength-1)/2;

  unsigned int center = IMUL(y,pitch)+x;
  unsigned int Scenter = IMUL(threadIdx.y+hkl,blockDim.x)+threadIdx.x;

  extern __shared__ double input_shared_SKconv[];

  if (x<width)
  {
    // Load data into shared memory
    input_shared_SKconv[Scenter] = input[min(y,(height-1))*pitch+min(x,width-1)];

    if (threadIdx.y < hkl)
    {
      // load on top and bottom
      input_shared_SKconv[Scenter-hkl*blockDim.x] = input[IMUL(max(0, min(height-1, y-hkl)), pitch)+x];
      input_shared_SKconv[Scenter+blockDim.y*blockDim.x] = input[IMUL(max(0, min(height-1, y+blockDim.y)), pitch)+x];
    }
  }

  __syncthreads();

  // Convolve with kernel
  double value = 0;
  if (x<width && y<height)
  {
    for (int i=-hkl; i<=hkl; i++)
      value += input_shared_SKconv[Scenter+i*blockDim.x]*small_kernel_const[i+hkl];
    output[center] = value;
  }
}


////////////////////////////////////////////////////////////////////////////////
//! Perform a separable kernel convolution with an 1D Kernel in z direction
//! @param input          pointer to input image
//! @param output         pointer to output image
//! @param klength        kernel length
//! @param width          width of the image
//! @param height         height of the image
//! @param depth          depth of the image
//! @param pitch_w        image pitch (width)
//! @param pitch_h        image pitch (height)
//! @param y              the current y plane
//! Make sure that the kernel length "klength" is smaller than "blockDim.x"
////////////////////////////////////////////////////////////////////////////////
__global__ void SKConvolutionZKernel(double* input, double* output, int klength,
                                     int width, int height, int depth,
                                     int pitch_w, int pitch_h, int y)
{
  // calculate absolute coordinates
  int x = IMUL(blockIdx.x, blockDim.x) + threadIdx.x;
  int z = IMUL(blockIdx.y, blockDim.y) + threadIdx.y;

  int hkl = (klength-1)/2;

  int center = IMUL(IMUL(z,pitch_w),pitch_h)+IMUL(y,pitch_w)+x;
  int Scenter = IMUL(threadIdx.y+hkl,blockDim.x)+threadIdx.x;

  extern __shared__ double input_shared_SKconv[];

  if (x<width)
  {
    // Load data into shared memory
    input_shared_SKconv[Scenter] = input[IMUL(IMUL(min(z,depth-1),pitch_w),pitch_h)+IMUL(min(y,height-1),pitch_w)+min(x,width-1)];

    if (threadIdx.y < hkl)
    {
      // load on top and bottom
      input_shared_SKconv[Scenter-hkl*blockDim.x] =
          input[IMUL(IMUL(max(min(z-hkl,depth-1),0),pitch_w),pitch_h)+IMUL(min(y,height-1),pitch_w)+min(x,width-1)];
      input_shared_SKconv[Scenter+blockDim.y*blockDim.x] =
          input[IMUL(IMUL(max(min(z+blockDim.y,depth-1),0),pitch_w),pitch_h)+IMUL(min(y,height-1),pitch_w)+min(x,width-1)];
    }
  }

  __syncthreads();

  // do convolution
  double value = 0;
  if (x<width && z<depth)
  {
    for (int i=-hkl; i<=hkl; i++)
      value += input_shared_SKconv[Scenter+i*blockDim.x]*small_kernel_const[i+hkl];
    output[center] = value;
  }
}

#endif // SK_CONVOLUTION_KERNELS_CU

