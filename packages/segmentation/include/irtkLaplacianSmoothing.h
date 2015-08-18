/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkLaplacianSmoothing.h 998 2013-10-15 15:24:16Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-10-15 16:24:16 +0100 (Tue, 15 Oct 2013) $
  Version   : $Revision: 998 $
  Changes   : $Author: mm3 $

=========================================================================*/

#ifndef _irtkLaplacianSmoothing_H

#define _irtkLaplacianSmoothing_H

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkGaussianBlurring.h>
#include <irtkGaussianBlurringWithPadding.h>
#include <irtkResampling.h>
#include <irtkResamplingWithPadding.h>



class irtkLaplacianSmoothing : public irtkObject
{

protected:

  irtkRealImage _image;
  irtkRealImage _fieldmap;
  irtkRealImage _mask;
  irtkRealImage _m;
  irtkRealImage _target;
  bool _have_target;
  int _directions[13][3]; 

  vector<double> _factor;//(13,0);
  
  double _lap_threshold;
  double _rel_diff_threshold;
  double _relax_iter;

  
  void InitializeFactors();
  
  void Blur(irtkRealImage& image, double sigma, double padding);
  void Resample(irtkRealImage& image, double step, double padding);
  void ResampleOnGrid(irtkRealImage& image, irtkRealImage templ);

  void EnlargeImage(irtkRealImage &image);
  void ReduceImage(irtkRealImage &image);
  
  void CalculateBoundary(irtkRealImage& m);
  void CalculateBoundaryWeights(irtkRealImage& w, irtkRealImage m);
  double LaplacianImage(irtkRealImage image, irtkRealImage m, irtkRealImage& laplacian);
  double LaplacianBoundary(irtkRealImage image, irtkRealImage m, irtkRealImage& laplacian);
  void LLaplacian(irtkRealImage& llaplacian, irtkRealImage laplacian, irtkRealImage mask, irtkRealImage weights);
  void UpdateFieldmap(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage multipliers, 
		    irtkRealImage mask, irtkRealImage weights, double alpha);
  void UpdateFieldmapGD(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage laplacian, irtkRealImage mask, irtkRealImage weights, double alpha, double lambda1, double lambda2);
  
  void UpdateFieldmapWithThreshold(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage multipliers, 
		    irtkRealImage mask, irtkRealImage weights, double alpha, double threshold);
  void UpdateFieldmapHuber(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage multipliers, 
		    irtkRealImage mask, irtkRealImage weights, double alpha);
  
  void Smooth(irtkRealImage& im, irtkRealImage m);
  void SmoothGD(irtkRealImage& im, irtkRealImage m);
  void UpsampleFieldmap(irtkRealImage& target, irtkRealImage mask, irtkRealImage newmask);
  
public:
  
  irtkLaplacianSmoothing();
  
  void SetInput(irtkRealImage image);
  void SetMask(irtkRealImage mask);
  void SetParam(double lap_threshold, double rel_diff_threshold, double relax_iter);
  
  irtkRealImage Run();
  irtkRealImage RunGD();
  irtkRealImage Run3levels();
  irtkRealImage Run1level();

};

inline void irtkLaplacianSmoothing::SetInput(irtkRealImage image)
{
  _image = image;
  _mask = _image;
  _mask=1;
}

inline void irtkLaplacianSmoothing::SetMask(irtkRealImage mask)
{
  _mask = mask;
}

inline void irtkLaplacianSmoothing::SetParam(double lap_threshold, double rel_diff_threshold, double relax_iter)
{
  _lap_threshold = lap_threshold;
  _rel_diff_threshold = rel_diff_threshold;
  _relax_iter = relax_iter;
}

#endif
