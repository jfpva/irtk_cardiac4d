/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

#ifndef gaussFilterConvolution_cuh
#define gaussFilterConvolution_cuh


#include <stdio.h>

#define IMUL(a, b) __mul24(a, b)

#define BLOCK_SIZE 8
#define BLOCK_SIZE_SK_1         256
#define BLOCK_SIZE_SK_2         4

#if 1
#define CHECK_ERROR(function) if (cudaError_t err = cudaGetLastError()) \
{ \
  printf("CUDA Error in " #function "(), line %d: %s\n", __LINE__, cudaGetErrorString(err)); \
}
#else
#define CHECK_ERROR(function)
#endif

#define MAX_LENGdoubleH_SK BLOCK_SIZE_SK_1-1
#define MAX_LENGTH_SK BLOCK_SIZE_SK_1-1
#define SUM_BLOCK 64

int SKConvolution(double* input, double* output, double* kernel, int klength,
                  unsigned int width, unsigned int height, unsigned int pitch,
                  unsigned int num_ch, double* temp = NULL);

int SKConvolution3D(double* input, double* output, double* kernel, int klength,
                    unsigned int width, unsigned int height, unsigned int depth,
                    unsigned int pitch_w, unsigned int pitch_h,
                    unsigned int num_ch, double* temp = NULL, double* temp2 = NULL);


int SKConvolutionX(double* input, double* output, double* kernel, int klength,
                   unsigned int width, unsigned int height, unsigned int pitch,
                   unsigned int num_ch);


int SKConvolutionY(double* input, double* output, double* kernel, int klength,
                   unsigned int width, unsigned int height, unsigned int pitch,
                   unsigned int num_ch);


int generateGaussianKernel(double** kernel, double sigma, int klength);

int createGenericGaussianKernel(double* kernel, int height, int width,
                                double sigma,
                                double* mask, double confidence);

int FilterGauss(double* input, double* output, unsigned int width,
                unsigned int height, unsigned int pitch, unsigned int num_ch,
                double sigma, double* temp = NULL);


int FilterGauss(double* input, double* output, unsigned int width,
                unsigned int height, unsigned int depth,
                unsigned int pitch_w, unsigned int pitch_h, unsigned int num_ch,
                double sigma, double* temp = NULL);

int FilterGaussX(double* input, double* output, unsigned int width,
                 unsigned int height, unsigned int pitch, unsigned int num_ch,
                 double sigma);

int FilterGaussY(double* input, double* output, unsigned int width,
                 unsigned int height, unsigned int pitch, unsigned int num_ch,
                 double sigma);
//template <class T>
double sum2D(double* sum, double* input, unsigned int width, unsigned int height, unsigned int pitch);
float sum2D(float* sum, float* input, unsigned int width, unsigned int height, unsigned int pitch);

int iDivUp(int a, int b);












#endif