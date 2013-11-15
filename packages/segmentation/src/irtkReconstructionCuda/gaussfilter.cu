/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

//! Functionality for performing gaussian filtering

#ifndef GAUSSFILTER_CU
#define GAUSSFILTER_CU


#include <stdio.h>
#include <npp.h>
#include "gaussFilterConvolution.cuh"
#include "gaussfilter_kernel.cu"


int iDivUp(int a, int b)
{
  return (a % b != 0) ? (a / b + 1) : (a / b);
}



//!/////////////////////////////////////////////////////////////////////////////
//! General Functions
//!/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! Generate 1D Gaussian convolution kernel
//! @param kernel    resulting kernel (necassary memory will be allocated)
//! @param sigma     sigma
//! @param klength   klength of the kernel
////////////////////////////////////////////////////////////////////////////////
int generateGaussianKernel(double** kernel, double sigma, int klength)
{
  // check for valid filter length
  if ((klength % 2) == 0)
  {
    fprintf(stderr, "Error: Convolution Kernel length even\n");
    return -1;
  }

  // allocate memory for kernel
  *kernel = (double*)malloc(sizeof(double) * klength);

  // sum for normalization
  double sum = 0;

  // compute kernel values
  int mid_point = (int)floor(klength/2.0f);
  for( int i = 0; i < klength; i++)
  {
    // generate value
    (*kernel)[i] = exp(-(double)abs(i-mid_point)*(double)abs(i-mid_point)/(2*sigma*sigma));

    // update sum for normalization
    sum += (*kernel)[i];
  }

  // normalize kernel
  for(int i = 0; i < klength; i++)
    (*kernel)[i] /= sum;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! Create generic gaussian convolution kernel with defective values
//! @param kernel    kernel
//! @param height    height of the kernel (must be odd)
//! @param width     width of the kernel (must be odd)
//! @param sigma     sigma
//! @param mask      mask, can be zero, kernel[i] = 0 if mask[i] < confidence
//! @param confidence
////////////////////////////////////////////////////////////////////////////////
int createGenericGaussianKernel(double* kernel, int height, int width,
                                double sigma,
                                double* mask, double confidence)
{
    if (width%2 == 0 || height%2 == 0)
      fprintf(stderr, "Error: convolution kernel dimensions must be odd\n");

    size_t pitch;

    double *d_kernel, *d_mask;
    cudaMallocPitch((void**) &d_kernel, &pitch, width*sizeof(double), height);
    CHECK_ERROR(createGenericGaussianKernel);

    if (mask != NULL) {
        cudaMallocPitch((void**) &d_mask, &pitch, width*sizeof(double), height);
        CHECK_ERROR(createGenericGaussianKernel);
        cudaMemcpy2D(d_mask, pitch,
                     mask, width*sizeof(double),
                     width*sizeof(double), height,
                     cudaMemcpyHostToDevice);
        CHECK_ERROR(createGenericGaussianKernel);
    }
    else
        d_mask = 0;

    // x * y number of blocks beeing lauched
    dim3 size_of_grid(iDivUp(width,BLOCK_SIZE), iDivUp(height,BLOCK_SIZE),1);

    // x * y * z number of threads per block */
    dim3 size_of_block(BLOCK_SIZE,BLOCK_SIZE,1);

    createGenericGaussian_Kernel<<<size_of_grid, size_of_block>>>
        (d_kernel, height, width, pitch, d_mask, confidence, sigma/width/*2*sigma*sigma*/);

    // compute sum
    double sum = 0;

	//double sum;
	/*	cudaHostRegister(&sum,sizeof(double), cudaHostRegisterMapped);
		Npp64f* d_sum;
		cudaHostGetDevicePointer((void **)&d_sum, (void *)&sum, 0);
     NppiSize oSizeROI = {width,height};
     int hpBufferSize;
		nppiSumGetBufferHostSize_32f_C1R(oSizeROI, &hpBufferSize);
		Npp8u *hpDeviceBuffer;
		cudaMalloc((void **)&hpDeviceBuffer, hpBufferSize);
		nppiSum_32f_C1R ((Npp32f*)(tmp.data()), tmp.size.x*sizeof(float), oSizeROI, hpDeviceBuffer, d_sum);*/


    sum2D(&sum, d_kernel, width, height, pitch/sizeof(double));

    cudaMemcpy2D(kernel, width*sizeof(double),
                 d_kernel, pitch,
                 width * sizeof(double), height,
                 cudaMemcpyDeviceToHost);
    CHECK_ERROR(createGenericGaussianKernel);

    for (int i = 0; i < (height * width); i++)
        kernel[i] /= sum;

    cudaFree(d_kernel);
    if (mask != 0)
        cudaFree(d_mask);

    return 0;
}



inline int round__(double x)
{
  return x > 0 ? int(x + 0.5) : int(x - 0.5);
}

////////////////////////////////////////////////////////////////////////////////
//! Performes simple gaussian filtering of an image
//! The size of the kernel is determined automatically
//! @param input          pointer to input image
//! @param output         pointer to output image (can be the same as input)
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
//! @param num_ch         number of channels in the image
//! @param sigma          sigma parameter to construct kernel
//! @param temp           optional parameter to set own temporary variable.
//!                       has to be of size pitch*height. If not set, memory
//!                       will be allocated.
////////////////////////////////////////////////////////////////////////////////
int FilterGauss(double* input, double* output, unsigned int width,
                unsigned int height, unsigned int pitch, unsigned int num_ch,
                double sigma, double* temp)
{
  int ret;

  // irtkScalarGaussian gaussianX(this->_Sigma/xsize, 1, 1, 0, 0, 0);
  // Create filter kernel for 1D Gaussian function in X
  //irtkGenericImage<irtkRealPixel> kernelX(2*round(4*this->_Sigma/xsize)+1, 1, 1);


  // Automatically determine kernel length
 // int klength = max(min((int)(sigma*5),MAX_LENGTH_SK),7);//max(min((int)(2*(4*sigma)+1),MAX_LENGTH_SK),7);//2*round__(4*sigma)+1;//max(min((int)(sigma*5),MAX_LENGTH_SK),7);
  // ... and make sure it is odd
 // klength -= 1-klength%2;
  int klength = 2*round__(4*sigma)+1;
  klength -= 1-klength%2;

  // Generate small Gaussian Kernel
  double* kernel = NULL;
  ret = generateGaussianKernel(&kernel, sigma, klength);
  if (ret)
  {
    fprintf(stderr, "Error in CUDA FilterGauss(): Could not generate Kernel\n");
    return ret;
  }

  // Convolve image with kernel
  ret = SKConvolution(input, output, kernel, klength, width, height, pitch,
                      num_ch, temp);
  if (ret)
  {
    fprintf(stderr, "Error in CUDA FilterGauss(): Convolution failed\n");
    return ret;
  }

  // Clean up memory
  free(kernel);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! Performes simple gaussian filtering of an 3D image
//! The size of the kernel is determined automatically
//! @param input          pointer to input image
//! @param output         pointer to output image (can be the same as input)
//! @param width          width of the image
//! @param height         height of the image
//! @param depth          depth of the image
//! @param pitch_w        image pitch (width)
//! @param pitch_h        image pitch (height)
//! @param num_ch         number of channels in the image
//! @param sigma          sigma parameter to construct kernel
//! @param temp           optional parameter to set own temporary variable.
//!                       has to be of size pitch*height. If not set, memory
//!                       will be allocated.
////////////////////////////////////////////////////////////////////////////////
int FilterGauss(double* input, double* output, unsigned int width,
                unsigned int height, unsigned int depth,
                unsigned int pitch_w, unsigned int pitch_h, unsigned int num_ch,
                double sigma, double* temp)
{
  int ret;

  // Automatically determine kernel length
  //int klength = max(min((int)(sigma*5),MAX_LENGTH_SK),7);
  // ... and make sure it is odd
  //klength -= 1-klength%2;
  int klength = 2*round__(4*sigma)+1;
  klength -= 1-klength%2;
  //printf("klength %d %d\n", klength, MAX_LENGTH_SK);

  // Generate small Gaussian Kernel
  double* kernel = NULL;
  ret = generateGaussianKernel(&kernel, sigma, klength);
  if (ret)
  {
    fprintf(stderr, "Error in CUDA FilterGauss(): Could not generate Kernel\n");
    return ret;
  }

  // Convolve image with kernel
  ret = SKConvolution3D(input, output, kernel, klength, width, height, depth,
                        pitch_w, pitch_h, num_ch, temp);
  if (ret)
  {
    fprintf(stderr, "Error in CUDA FilterGauss(): Convolution failed\n");
    return ret;
  }

  // Clean up memory
  free(kernel);

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//! Performes simple gaussian filtering of an image only in X direction
//! The size of the kernel is determined automatically
//! Does not work in place !!! Make sure output != input
//! @param input          pointer to input image
//! @param output         pointer to output image (can be the same as input)
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
//! @param num_ch         number of channels in the image
//! @param sigma          sigma parameter to construct kernel
////////////////////////////////////////////////////////////////////////////////
int FilterGaussX(double* input, double* output, unsigned int width,
                 unsigned int height, unsigned int pitch, unsigned int num_ch,
                 double sigma)
{
  int ret;

  // Automatically determine kernel length
  int klength = max(min((int)(sigma*5),MAX_LENGTH_SK),7);
  // ... and make sure it is odd
  klength -= 1-klength%2;

  // Generate small Gaussian Kernel
  double* kernel = NULL;
  ret = generateGaussianKernel(&kernel, sigma, klength);
  if (ret)
  {
    fprintf(stderr, "Error in CUDA FilterGaussX(): Could not generate Kernel\n");
    return ret;
  }

  // Convolve image with kernel
  ret = SKConvolutionX(input, output, kernel, klength, width, height, pitch,
                       num_ch);
  if (ret)
  {
    fprintf(stderr, "Error in CUDA FilterGaussX(): Convolution failed\n");
    return ret;
  }

  // Clean up memory
  free(kernel);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! Performes simple gaussian filtering of an image only in Y direction
//! The size of the kernel is determined automatically
//! Does not work in place !!! Make sure output != input
//! @param input          pointer to input image
//! @param output         pointer to output image (can be the same as input)
//! @param width          width of the image
//! @param height         height of the image
//! @param pitch          image pitch
//! @param num_ch         number of channels in the image
//! @param sigma          sigma parameter to construct kernel
////////////////////////////////////////////////////////////////////////////////
int FilterGaussY(double* input, double* output, unsigned int width,
                 unsigned int height, unsigned int pitch, unsigned int num_ch,
                 double sigma)
{
  int ret;

  // Automatically determine kernel length
  int klength = max(min((int)(sigma*5),MAX_LENGTH_SK),7);
  // ... and make sure it is odd
  klength -= 1-klength%2;

  // Generate small Gaussian Kernel
  double* kernel = NULL;
  ret = generateGaussianKernel(&kernel, sigma, klength);
  if (ret)
  {
    fprintf(stderr, "Error in CUDA FilterGaussY(): Could not generate Kernel\n");
    return ret;
  }

  // Convolve image with kernel
  ret = SKConvolutionY(input, output, kernel, klength, width, height, pitch,
                       num_ch);
  if (ret)
  {
    fprintf(stderr, "Error in CUDA FilterGaussY(): Convolution failed\n");
    return ret;
  }

  // Clean up memory
  free(kernel);

  return 0;
}


#endif // GAUSSFILTER_CU
