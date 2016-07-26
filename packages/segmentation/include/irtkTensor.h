/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkTensor.h 577 2012-03-30 09:54:05Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2012-03-30 10:54:05 +0100 (Fri, 30 Mar 2012) $
  Version   : $Revision: 577 $
  Changes   : $Author: mm3 $

=========================================================================*/

#ifndef _irtkTensor_H

#define _irtkTensor_H

#include <irtkImage.h>

class irtkTensor : public irtkVector
{
public:
  irtkTensor();
  irtkTensor(const irtkVector& v);
  double Evaluate(double gx, double gy, double gz, double b);
  double PDF(double dx, double dy, double dz);
  irtkMatrix Matrix();
  void Set(const irtkMatrix& m);
  irtkVector Direction();      
  double FA();
  double Volume();
  
  ///Operators
  void Log();  
  void Exp();
  
  irtkTensor  operator+ (irtkTensor&);
  irtkTensor  operator- (irtkTensor&);
  irtkTensor  operator* (double);
  irtkTensor  operator/ (double);
  
  double Norm();
  double Dist(const irtkTensor&);
  double EuclDist(const irtkTensor&);
  bool IsPD();
  void MakePD();
  double Det();

};


inline irtkTensor irtkTensor::operator+ (irtkTensor& v)
{
  //t=exp(log(t1)+log(t2))
  irtkTensor t1,t2;
  //the tensors to add
  t1=*this;
  t2=v;
  //take logarithms
  t1.Log();
  t2.Log();
  //Vector addition
  t1 = t1.irtkVector::operator+(t2);
  //take exponential 
  t1.Exp();
  
  return t1;
}

inline irtkTensor irtkTensor::operator- (irtkTensor& v)
{
  //t=exp(log(t1)-log(t2))
  irtkTensor t1,t2;
  //the tensors to subtract
  t1=*this;
  t2=v;
  //take logarithms
  t1.Log();
  t2.Log();
  //Vector subtraction
  t1 = t1.irtkVector::operator-(t2);
  //take exponential 
  t1.Exp();
  
  return t1;
}

inline irtkTensor irtkTensor::operator* (double a)
{
  //t=exp(a*log(t1))
  irtkTensor t=*this;
  //take logarithm
  t.Log();
  //Vector mutiptication by scalar
  t = t.irtkVector::operator*(a);
  //take exponential 
  t.Exp();
  
  return t;
}

inline irtkTensor irtkTensor::operator/ (double a)
{
  //t=exp(log(t1)/a)
  irtkTensor t=*this;
  //take logarithm
  t.Log();
  //Vector division by scalar
  t = t.irtkVector::operator/(a);
  //take exponential 
  t.Exp();
  
  return t;
}

inline double irtkTensor::Norm()
{
  irtkTensor t=*this;
  //take logarithm
  t.Log();
  
  double norm = sqrt(t(0)*t(0)+t(1)*t(1)+t(2)*t(2)+
                +2*(t(3)*t(3)+t(4)*t(4)+t(5)*t(5)));
  return norm;
}

inline double irtkTensor::Dist(const irtkTensor& v)
{
  irtkTensor t1,t2;
  t1=*this;
  t2=v;
  //take logarithms
  t1.Log();
  t2.Log();
  //Vector subtraction
  t1 = t1.irtkVector::operator-(t2);
  //calculate norm
  double norm = sqrt(t1(0)*t1(0)+t1(1)*t1(1)+t1(2)*t1(2)+
                +2*(t1(3)*t1(3)+t1(4)*t1(4)+t1(5)*t1(5)));
  return norm;
}

inline double irtkTensor::EuclDist(const irtkTensor& v)
{
  irtkTensor t1,t2;
  t1=*this;
  t2=v;
  //Vector subtraction
  t1 = t1.irtkVector::operator-(t2);
  //calculate norm
  double norm = sqrt(t1(0)*t1(0)+t1(1)*t1(1)+t1(2)*t1(2)+
                +2*(t1(3)*t1(3)+t1(4)*t1(4)+t1(5)*t1(5)));
  return norm;
}

inline bool irtkTensor::IsPD()
{
  bool pd;
  if(_vector[0]<=0){
    pd=false;
  }
  else if ((_vector[0]*_vector[1]-_vector[3]*_vector[3])<=0){
    pd=false;
  }
  else if (Det()<=0){
    pd=false;
  }
  else
    pd=true;
  
  return pd;
}


inline double irtkTensor::Det()
{
  double det = _vector[0]*_vector[1]*_vector[2]
             + 2*_vector[3]*_vector[4]*_vector[5]
             - _vector[2]*_vector[3]*_vector[3]
             - _vector[1]*_vector[4]*_vector[4]
             - _vector[0]*_vector[5]*_vector[5];  
  return det;
}

inline double irtkTensor::Evaluate(double gx, double gy, double gz, double b)
{
  double val = gx*gx*_vector[0]
	     + gy*gy*_vector[1]
	     + gz*gz*_vector[2]
	     + gx*gy*_vector[3]*2
	     + gx*gz*_vector[4]*2
	     + gy*gz*_vector[5]*2;
  return exp(-b*val);     
}
#endif