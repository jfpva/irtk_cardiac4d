/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkDWImage.h 577 2012-03-30 09:54:05Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2012-03-30 10:54:05 +0100 (Fri, 30 Mar 2012) $
  Version   : $Revision: 577 $
  Changes   : $Author: mm3 $

=========================================================================*/

#ifndef _irtkDWImage_H

#define _irtkDWImage_H

#include <irtkImage.h>
#include <irtkTensor.h>

class irtkDWImage : public irtkRealImage
{
protected:
  irtkMatrix _directions;
  double _b;
public:
  irtkDWImage(const irtkRealImage& image);
  void SetDirections(irtkMatrix directions);
  void SetDirection(int index, irtkVector direction);
  void SetB(double b);
  int GetNumberOfDirections();
  irtkTensor CalculateTensor(int x,int y, int z);
};

inline int irtkDWImage::GetNumberOfDirections()
{
  return GetT()-1;
}

inline void irtkDWImage::SetB(double b)
{
  _b=b;
}


#endif