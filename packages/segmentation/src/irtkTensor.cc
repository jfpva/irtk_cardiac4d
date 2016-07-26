#include <irtkTensor.h>

irtkTensor::irtkTensor():irtkVector(6)
{  
}

irtkTensor::irtkTensor(const irtkVector& v):irtkVector(v)
{  
}

irtkMatrix irtkTensor::Matrix()
{
  irtkMatrix m(3,3);
  m(0,0)=_vector[0];
  m(1,1)=_vector[1];
  m(2,2)=_vector[2];
  
  m(0,1)=_vector[3];
  m(1,0)=_vector[3];
  
  m(0,2)=_vector[4];
  m(2,0)=_vector[4];
  
  m(1,2)=_vector[5];
  m(2,1)=_vector[5];
  
  return m;
}

void irtkTensor::Set(const irtkMatrix& m)
{
  _vector[0]=m(0,0); 
  _vector[1]=m(1,1); 
  _vector[2]=m(2,2); 

  _vector[3]=m(0,1); 
  _vector[4]=m(0,2); 
  _vector[5]=m(1,2); 
}

double irtkTensor::PDF(double dx, double dy, double dz)
{
  irtkMatrix m = Matrix();
  
  m.Invert();
  
  irtkVector v(3);
  v(0)=dx; v(1)=dy; v(2)=dz;
   
  v=m*v;
  
  double value = 0.5*(dx*v(0)+dy*v(1)+dz*v(2));
  return sqrt(m.Det()/pow(2*3.1415,3))*exp(-value);
}

irtkVector irtkTensor::Direction()
{
  irtkMatrix m = Matrix();
  irtkMatrix tmp = m;
  irtkVector v(3);
  m.Eigenvalues(m,v,tmp);
  
    
  //Principal direction
  irtkVector d(3);
  
  //find largest eigenvalue
  int ind = 0;
  double max = v(0);  
  for(int i=1;i<3;i++)
    if(v(i)>max)
    {
      max=v(i);
      ind=i;
    }
    
  //cout<<"Max eigenvalue:"<<max<<endl;
  for(int i=0;i<3;i++)
    d(i)=m(i,ind);  
  
  return d;
  
}


double irtkTensor::FA()
{
  if (!IsPD())
    return 0;
  
  //diffusion matrix
  irtkMatrix m = Matrix();
  irtkMatrix tmp = m;
  
  //eigenvalues
  irtkVector v(3);
  m.Eigenvalues(m,v,tmp);
  
  //average eigenvalue
  double av = (v(0)+v(1)+v(2))/3;
  
  irtkVector va=v;
  for(int i=0;i<3;i++)
    va(i)-=av;
  double fa = sqrt(1.5*(va(0)*va(0)+va(1)*va(1)+va(2)*va(2))/(v(0)*v(0)+v(1)*v(1)+v(2)*v(2)));
   return fa;
}

double irtkTensor::Volume()
{
  //diffusion matrix
  irtkMatrix m = Matrix();
  irtkMatrix tmp = m;
  
  //eigenvalues
  irtkVector v(3);
  m.Eigenvalues(m,v,tmp);
  
  //average eigenvalue
  return v(0)*v(1)*v(2);
}

void irtkTensor::Log()
{
  //difusion matrix
  irtkMatrix m = Matrix();
  irtkMatrix tmp = m;
  
  //eigenanalysis m=r*diag(v)*rT
  irtkVector v(3);
  irtkMatrix r(m);
  r.Eigenvalues(r,v,tmp);
  
  //create diagonal matrix d with log-transform eigenvalues
  irtkMatrix d(3,3);
  d.Ident();
  d(0,0)=log(v(0));
  d(1,1)=log(v(1));
  d(2,2)=log(v(2));
  
  //create log(m)=r*d*rT
  irtkMatrix lm=r*d;
  r.Transpose();
  lm=lm*r;
  
  //Set the tensor according to the new matrix
  Set(lm);  
}

void irtkTensor::Exp()
{
  //difusion matrix
  irtkMatrix m = Matrix();
  irtkMatrix tmp = m;
  
  //eigenanalysis m=r*diag(v)*rT
  irtkVector v(3);
  irtkMatrix r(m);
  r.Eigenvalues(r,v,tmp);
  
  //create diagonal matrix d with exp-transformed eigenvalues
  irtkMatrix d(3,3);
  d.Ident();
  d(0,0)=exp(v(0));
  d(1,1)=exp(v(1));
  d(2,2)=exp(v(2));
  
  //create log(m)=r*d*rT
  irtkMatrix lm=r*d;
  r.Transpose();
  lm=lm*r;
  
  //Set the tensor according to the new matrix
  Set(lm);  
}

void irtkTensor::MakePD()
{
  //difusion matrix
  irtkMatrix m = Matrix();
  irtkMatrix tmp = m;
  
  //eigenanalysis m=r*diag(v)*rT
  irtkVector v(3);
  irtkMatrix r(m);
  r.Eigenvalues(r,v,tmp);

  //Check if the tensor is positive definite and if not correct
  double threshold = 0.000001;
  for(int ii=0;ii<3;ii++)
    if(v(ii)<threshold)
      v(ii)=threshold;
  
  //create diagonal matrix d with exp-transformed eigenvalues
  irtkMatrix d(3,3);
  d.Ident();
  d(0,0)=v(0);
  d(1,1)=v(1);
  d(2,2)=v(2);
  
  //create log(m)=r*d*rT
  irtkMatrix lm=r*d;
  r.Transpose();
  lm=lm*r;
  
  //Set the tensor according to the new matrix
  Set(lm);  
}

