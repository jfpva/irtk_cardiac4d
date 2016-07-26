
#ifndef _irtkSphericalHarmonics_H

#define _irtkSphericalHarmonics_H

#include <irtkImage.h>
#include <irtkVector.h>
#include <irtkMatrix.h>

class irtkSphericalHarmonics : public irtkObject
{
public:
  
  irtkMatrix _SHT;
  irtkMatrix _iSHT;
  
  irtkMatrix shTransform(irtkMatrix directions, int lmax);
  int NforL (int lmax);
  double LegendrePolynomialsHelper (const double x, const double m);
  void LegendrePolynomials(irtkVector& array, const int lmax, const int m, const double x); 
  int shIndex (int l, int m);
  irtkMatrix cartesian2spherical (irtkMatrix xyz);
  irtkVector LSFit(irtkVector signal, irtkMatrix dirs, int lmax);
  void InitSHT(irtkMatrix dirs, int lmax);
  irtkMatrix SHbasis(irtkMatrix dirs, int lmax);
  irtkVector Coeff2Signal(irtkVector c);
  irtkVector Signal2Coeff(irtkVector s);
  irtkRealImage Signal2Coeff(irtkRealImage signal);
  irtkRealImage Coeff2Signal(irtkRealImage coeff);
  irtkMatrix LaplaceBeltramiMatrix(int lmax);
  void InitSHTRegul(irtkMatrix dirs, double lambda, int lmax);
};

inline int irtkSphericalHarmonics::NforL (int lmax)
{
  return (lmax+1) * (lmax+2) /2;
}

inline double irtkSphericalHarmonics::LegendrePolynomialsHelper (const double x, const double m)
{
  return (m < 1.0 ? 1.0 : x * (m-1.0) / m * LegendrePolynomialsHelper(x, m-2.0));
}

inline int irtkSphericalHarmonics::shIndex (int l, int m)
{
  return l * (l+1) /2 + m;
}

#endif