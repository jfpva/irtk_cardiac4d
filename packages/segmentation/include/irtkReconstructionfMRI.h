/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionfMRI.h 998 2013-10-15 15:24:16Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-10-15 16:24:16 +0100 (Tue, 15 Oct 2013) $
  Version   : $Revision: 998 $
  Changes   : $Author: mm3 $

=========================================================================*/

#ifndef _irtkReconstructionfMRI_H

#define _irtkReconstructionfMRI_H

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkGaussianBlurring.h>
#include <irtkReconstruction.h>
#include <irtkReconstructionb0.h>
#include <irtkMultiChannelImage.h>
#include <irtkBSplineReconstruction.h>

#include <vector>
using namespace std;

class irtkReconstructionfMRI : public irtkReconstructionb0
{

protected:
 
	vector<irtkRealImage> _timeserie;
	
public:
  
	irtkReconstructionfMRI();
	
	void InterpolateBSpline(vector<irtkRealImage>& stacks, int iter);
	
	void InterpolateGaussian(vector<irtkRealImage>& stacks, int iter);
	
	void InterpolateBSplineReordered(vector<irtkRealImage>& stacks, vector<int> multiband_vector , int iter);
		
	void InterpolateGaussianReordered(vector<irtkRealImage>& stacks, vector<int> multiband_vector , int iter);
	  
};

#endif
