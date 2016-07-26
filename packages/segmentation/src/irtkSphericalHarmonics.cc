//#include <irtkVector.h>
//#include <irtkMatrix.h>
#include <irtkSphericalHarmonics.h>

/*
      template <class MatrixType>
        Eigen::Matrix<typename MatrixType::Scalar,Eigen::Dynamic, Eigen::Dynamic> init_transform (const MatrixType& dirs, int lmax)
        {
          using namespace Eigen;
          typedef typename MatrixType::Scalar value_type;
          if (dirs.cols() != 2)
            throw Exception ("direction matrix should have 2 columns: [ azimuth elevation ]");
          Matrix<value_type,Dynamic,Dynamic> SHT (dirs.rows(), NforL (lmax));
          Matrix<value_type,Dynamic,1,0,64> AL (lmax+1);
          for (ssize_t i = 0; i < dirs.rows(); i++) {
            value_type x = std::cos (dirs (i,1));
            Legendre::Plm_sph (AL, lmax, 0, x);
            for (int l = 0; l <= lmax; l+=2)
              SHT (i,index (l,0)) = AL[l];
            for (int m = 1; m <= lmax; m++) {
              Legendre::Plm_sph (AL, lmax, m, x);
              for (int l = ( (m&1) ? m+1 : m); l <= lmax; l+=2) {
#ifndef USE_NON_ORTHONORMAL_SH_BASIS
                SHT(i, index(l, m)) = Math::sqrt2 * AL[l]*std::cos (m*dirs (i,0));
                SHT(i, index(l,-m)) = Math::sqrt2 * AL[l]*std::sin (m*dirs (i,0));
#else
                SHT(i, index(l, m)) = AL[l]*std::cos (m*dirs (i,0));
                SHT(i, index(l,-m)) = AL[l]*std::sin (m*dirs (i,0));
#endif
              }
            }
          }
          return SHT;
        }
*/

irtkMatrix irtkSphericalHarmonics::shTransform(irtkMatrix dirs, int lmax)
{
  if (dirs.Cols() != 2)
  {
     cout<<"direction matrix should have 2 columns: [ azimuth elevation ]"<<endl; 
     exit(1);
  }
  
  irtkMatrix SHT(dirs.Rows(), NforL (lmax));
  irtkVector AL(lmax+1);
  
  for (int i = 0; i < dirs.Rows(); i++) 
  {
    double x = std::cos (dirs (i,1));
    LegendrePolynomials (AL, lmax, 0, x);
    for (int l = 0; l <= lmax; l+=2)
      SHT (i,shIndex (l,0)) = AL(l);
    for (int m = 1; m <= lmax; m++) 
    {
      LegendrePolynomials (AL, lmax, m, x);
      for (int l = ( (m&1) ? m+1 : m); l <= lmax; l+=2) 
      {
          SHT(i, shIndex(l, m)) = M_SQRT2 * AL(l)*std::cos (m*dirs (i,0));
          SHT(i, shIndex(l,-m)) = M_SQRT2 * AL(l)*std::sin (m*dirs (i,0));
      }      
    }    
  }
  return SHT;
  
}

void irtkSphericalHarmonics::LegendrePolynomials(irtkVector& array, const int lmax, const int m, const double x)
{
  double x2 = x*x;
  if (m && x2 >= 1.0) 
  {
    for (int n = m; n <= lmax; ++n)
    array(n) = 0.0;
    return;
  }
  array(m) = 0.282094791773878;
  //m>0
  if (m) array(m) *= std::sqrt (double (2*m+1) * LegendrePolynomialsHelper (1.0-x2, 2.0*m));
  //m is odd
  if (m & 1) array(m) = -array(m);
  if (lmax == m) return;

  double f = std::sqrt (double (2*m+3));
  array(m+1) = x * f * array(m);

  for (int n = m+2; n <= lmax; n++) 
  {
    array(n) = x*array(n-1) - array(n-2)/f;
    f = std::sqrt (double (4*n*n-1) / double (n*n-m*m));
    array(n) *= f;
  }
}

irtkMatrix irtkSphericalHarmonics::cartesian2spherical (irtkMatrix xyz)
{
  double r,x,y,z;
  irtkMatrix az_el_r(xyz.Rows(),2);
  
  for (int i=0; i<xyz.Rows();i++)
  {
    x=xyz(i,0);
    y=xyz(i,1);
    z=xyz(i,2);
    r = std::sqrt (x*x + y*y + z*z);
    az_el_r(i,0) = std::atan2 (y, x);
    az_el_r(i,1) = std::acos (z/r);
  }
  return az_el_r;
}

irtkVector irtkSphericalHarmonics::LSFit(irtkVector signal, irtkMatrix dirs, int lmax)
{
  // check that sizes correspond
  if(signal.Rows()!=dirs.Rows())
  {
    cerr<<"dimensions of signal and direction number do not match: "<<signal.Rows()<<" "<<dirs.Rows()<<endl;
    exit(1);
  }
  //convert to spherical coordinates
  irtkMatrix dirs_sph;
  dirs_sph = cartesian2spherical(dirs);
  //calculate matrix of basis functions (directions x basis number)
  irtkMatrix sht = shTransform(dirs_sph,lmax);
  //SH coefficients
  irtkVector c(sht.Cols());
  //calculate least square fit
  irtkMatrix temp = sht;
  temp.Transpose();
  c=temp*signal;
  temp=temp*sht;
  temp.Invert();
  c=temp*c;
  return c;
}

void irtkSphericalHarmonics::InitSHT(irtkMatrix dirs, int lmax)
{
  irtkMatrix dirs_sph;
  dirs_sph = cartesian2spherical(dirs);
  //calculate matrix of basis functions (directions x basis number)
  _SHT = shTransform(dirs_sph,lmax);
  //_SHT.Print();
  
  //calculate inverse matrix for least square fit
  irtkMatrix shtt = _SHT;
  shtt.Transpose();
  _iSHT = shtt*_SHT;
  _iSHT.Invert();
  _iSHT=_iSHT*shtt; 
  //_iSHT.Print();
}

irtkMatrix irtkSphericalHarmonics::SHbasis(irtkMatrix dirs, int lmax)
{
  irtkMatrix dirs_sph;
  dirs_sph = cartesian2spherical(dirs);
  //calculate matrix of basis functions (directions x basis number)
  irtkMatrix SHT = shTransform(dirs_sph,lmax);
  return SHT;
}

void irtkSphericalHarmonics::InitSHTRegul(irtkMatrix dirs, double lambda, int lmax)
{
  irtkMatrix dirs_sph;
  dirs_sph = cartesian2spherical(dirs);
  //calculate matrix of basis functions (directions x basis number)
  _SHT = shTransform(dirs_sph,lmax);
  
  //calculate inverse matrix for least square fit
  irtkMatrix shtt = _SHT;
  shtt.Transpose();
  irtkMatrix LB = LaplaceBeltramiMatrix(lmax);
  _iSHT = shtt*_SHT+LB*lambda;
  _iSHT.Invert();
  _iSHT=_iSHT*shtt; 
}


irtkMatrix irtkSphericalHarmonics::LaplaceBeltramiMatrix(int lmax)
{
  int dim = NforL (lmax);
  int value,l,k;
  irtkMatrix LB(dim, dim);
  //we can skip l==0, because coeff is 0 anyway
  for(l=2;l<=lmax;l=l+2)
  {
    value = l*l*(l+1)*(l+1);
    for(k=NforL(l-2); k<NforL(l); k++)
      LB(k,k)=value;
  }
  //LB.Print();
  return LB;
}

irtkVector irtkSphericalHarmonics::Coeff2Signal(irtkVector c)
{
  if(c.Rows()!=_SHT.Cols())
  {
    cerr<<"dimensions of SH coeffs and number of basis do not match: "<<c.Rows()<<" "<<_SHT.Cols()<<endl;
    exit(1);
  }
  return _SHT*c;  
}

irtkVector irtkSphericalHarmonics::Signal2Coeff(irtkVector s)
{
  if(s.Rows()!=_iSHT.Cols())
  {
    cerr<<"dimensions of signal and number of directions do not match: "<<s.Rows()<<" "<<_iSHT.Cols()<<endl;
    exit(1);
  }
  return _iSHT*s;
  
}

irtkRealImage irtkSphericalHarmonics::Signal2Coeff(irtkRealImage signal)
{
  if(signal.GetT()!=_iSHT.Cols())
  {
    cerr<<"dimensions of signal and number of directions do not match: "<<signal.GetT()<<" "<<_iSHT.Cols()<<endl;
    exit(1);
  }
  //create image with SH coeffs
  irtkImageAttributes attr = signal.GetImageAttributes();
  attr._t = _iSHT.Rows();
  irtkRealImage coeffs(attr);
  irtkVector s(signal.GetT());
  irtkVector a;
  int i,j,k,t;
  for(k=0;k<signal.GetZ();k++)
    for(j=0;j<signal.GetY();j++)
      for(i=0;i<signal.GetX();i++)
      {
	if(signal(i,j,k,0)>0)
	{
	  for(t=0;t<signal.GetT();t++)
	    s(t)=signal(i,j,k,t);
	  a=Signal2Coeff(s);
	  for(t=0;t<coeffs.GetT();t++)
	    coeffs(i,j,k,t)=a(t);
	}
      }
  
  
  return coeffs;
}


irtkRealImage irtkSphericalHarmonics::Coeff2Signal(irtkRealImage coeffs)
{
  if(coeffs.GetT()!=_SHT.Cols())
  {
    cerr<<"dimensions of SH coeffs and number of basis do not match: "<<coeffs.GetT()<<" "<<_SHT.Cols()<<endl;
    exit(1);
  }
  //create image with SH coeffs
  irtkImageAttributes attr = coeffs.GetImageAttributes();
  attr._t = _SHT.Rows();
  irtkRealImage signal(attr);
  irtkVector c(coeffs.GetT());
  irtkVector s;
  int i,j,k,t;
  for(k=0;k<coeffs.GetZ();k++)
    for(j=0;j<coeffs.GetY();j++)
      for(i=0;i<coeffs.GetX();i++)
      {
	//it should be ok - the first coeff should not be negative for positive signal
	//if(coeffs(i,j,k,0)>0)
	//{
	  for(t=0;t<coeffs.GetT();t++)
	    c(t)=coeffs(i,j,k,t);
	  s=Coeff2Signal(c);
	  for(t=0;t<signal.GetT();t++)
	    signal(i,j,k,t)=s(t);
	//}
      }
  
  
  return signal;
}