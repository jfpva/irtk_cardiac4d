/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

//! compute unnormalized gaussian kernel

//!//////////////////////////////////////////////////////////////////////////////
//! functions
//!//////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//! Compute unnormalized gaussian kernel
//! @param d_kernel       pointer to gaussian kernel data
//! @param width          width of the kernel (must be odd)
//! @param height         height of the kernel (must be odd)
//! @param pitch          kernel pitch (in bytes)
//! @param d_mask         pointer to mask data, can be NULL.
//!                       d_kernel[i] = 0 if d_mask[i] < confidence
//! @param confidence
//! @param sigma2         double variance (sigma2 = 2 * sigma�)
////////////////////////////////////////////////////////////////////////////////
__global__ void createGenericGaussian_Kernel(double *d_kernel,
                                             int height, int width, int pitch,
                                             double *d_mask, double confidence,
                                             double sigma2)
{
    int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    unsigned int poi = y * pitch / sizeof(double) + x;   // 1D coordinate of the image point (point of interest)

    if (x < width && y < height) {
        if (d_mask != 0 && d_mask[poi] < confidence)
            d_kernel[poi] = 0;
        else {
            int radius2_x = abs((width - 1) / 2 - x) * abs((width - 1) / 2 - x);
            int radius2_y = abs((height - 1) / 2 - y) * abs((height - 1) / 2 - y);
            d_kernel[poi] = exp( - (radius2_x + radius2_y) / sigma2);
        }
    }
}

