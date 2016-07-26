/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkTensorField.cc 577 2012-03-30 09:54:05Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2012-03-30 10:54:05 +0100 (Fri, 30 Mar 2012) $
  Version   : $Revision: 577 $
  Changes   : $Author: mm3 $

=========================================================================*/



#include <irtkTensorField.h>

irtkTensorField::irtkTensorField()
{  
  _x=0;
  _y=0;
  _z=0;
}

irtkTensorField::irtkTensorField(int x, int y, int z)
{  
  vector <irtkTensor> rows(z);
  vector < vector <irtkTensor> > columns(y,rows);
  vector < vector < vector <irtkTensor> > > data (x,columns);
  _tensorField = data;
  _x=x;
  _y=y;
  _z=z;
  cout<<"Created tensor field "<<_x<<" "<<_y<<" "<<_z<<endl;
}


  
irtkTensorField::irtkTensorField(irtkDWImage& image, double threshold)
{ 
  Initialise(image, threshold);
}


void irtkTensorField::Initialise(irtkDWImage& image, double threshold)
{ 
  //initialise member variables
  *this = irtkTensorField(image.GetX(), image.GetY(), image.GetZ());
  //zero tensor for padding
  irtkTensor zeroTensor;
  //Calculate tensors
  int i,j,k;
  for(i=0;i<_x;i++)
    for(j=0;j<_y;j++)
      for(k=0;k<_z;k++)
	if(image(i,j,k,0)>threshold)
	  _tensorField[i][j][k]=image.CalculateTensor(i,j,k);
	else
          _tensorField[i][j][k]=zeroTensor;
    cout<<"Initialised tensor field "<<_x<<" "<<_y<<" "<<_z<<endl;
}

void irtkTensorField::Log()
{ 
  //Calculate tensors
  int i,j,k;
  for(i=0;i<_x;i++)
    for(j=0;j<_y;j++)
      for(k=0;k<_z;k++)
	if(_tensorField[i][j][k].IsPD())
	  _tensorField[i][j][k].Log();;
    cout<<"Performed log of tensor field. "<<endl;
}
/*
void irtkTensorField::Exp()
{ 
  //Calculate tensors
  int i,j,k;
  for(i=0;i<_x;i++)
    for(j=0;j<_y;j++)
      for(k=0;k<_z;k++)
	if(_tensorField[i][j][k].IsPD())
	  _tensorField[i][j][k].Exp();;
    cout<<"Performed exp of tensor field. "<<endl;
}
*/
irtkTensor irtkTensorField::operator()(double x, double y, double z)
{
    //cout<<"tensorfield double"<<endl;

  if ((x >= 0) && (x <= (_x-1)) && (y >= 0) && (y <= (_y-1)) && (z >= 0) && (z <= (_z-1))) 
  {
    int nx,ny,nz,l,m,n;
    irtkTensor t;
    bool first = true;
    double sum=0;

    nx = (int) floor(x);
    ny = (int) floor(y);
    nz = (int) floor(z);    
    
    for (l = nx; l <= nx + 1; l++)
      if ((l >= 0) && (l < _x))
        for (m = ny; m <= ny + 1; m++)
          if ((m >= 0) && (m < _y))
            for (n = nz; n <= nz + 1; n++)
              if ((n >= 0) && (n < _z))
	      {
		//cout<<l<<" "<<m<<" "<<n<<endl;
		//cout<<(1 - fabs(l - x)) * (1 - fabs(m - y)) * (1 - fabs(n - z))<<endl;
		//_tensorField[l][m][n].Print();
		//cout.flush();
		
		double weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) * (1 - fabs(n - z));
		if ((weight > 0)&&(_tensorField[l][m][n].IsPD()))
		{
		  sum+=weight;
		  //cout<<"point 1; ";
		  //cout.flush();
		 irtkTensor temp = _tensorField[l][m][n]*((1 - fabs(l - x)) * (1 - fabs(m - y)) * (1 - fabs(n - z)));
		  //cout<<"point 2; ";
		  //t.Print();
		  //cout.flush();
		  if (first)
		  {
		    t=temp;
		    first=false;
		  }
		  else
		    t=t+temp;
		}
		  //cout<<"point 3; ";
		  //cout.flush();
		  //t.Print();
		  //cout<<"Sum="<<sum<<endl;

	      }
    if (sum>0.0001)
      return t/sum;
    else
      return t;
    
  } 
  else 
  {
    cout << "irtkTensorField::operator(): parameter out of range\n";
    cout<<x<<" "<<y<<" "<<z<<endl;
    cout<<_x<<" "<<_y<<" "<<_z<<endl;
    exit(1);
  }
}

void irtkTensorField::SmoothOld(double lambda)
{
  int i, j;
  int x,y,z,xx,yy,zz;
  double sum;
  bool first;
  
  int directions[13][3] =
    {
      { 1, 0, -1 },
      { 0, 1, -1 },
      { 1, 1, -1 },
      { 1, -1, -1 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 1, 1, 0 },
      { 1, -1, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { 1, 1, 1 },
      { 1, -1, 1 },
      { 0, 0, 1 }
    };

  double factor[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  for (i = 0; i < 13; i++)
  {
    for (j = 0; j < 3; j++)
      factor[i] += fabs(double(directions[i][j]));
    factor[i] = 1 / factor[i];
  }
  
  irtkTensorField old = *this;
  irtkTensor t,temp,D,Dn;
  
  for (x = 0; x < _x; x++)
    for (y = 0; y < _y; y++)
      for (z = 0; z < _z; z++)
	if(old(x,y,z).IsPD())
        {
	  //Calculate the sum of all coefficients
          sum = 0;
          for (i = 0; i < 13; i++)
          {
	    xx = x + directions[i][0];
	    yy = y + directions[i][1];
	    zz = z + directions[i][2];
	    

	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	        sum += factor[i];
	  }
	  
	  for (i = 0; i < 13; i++)
	  {
	    xx = x - directions[i][0];
	    yy = y - directions[i][1];
	    zz = z - directions[i][2];
	    
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	        sum += factor[i];
	  }
	  
	  //Smooth the image
	  first = true;
	  D=old(x,y,z);
	  
          for (i = 0; i < 13; i++)
          {
	    xx = x + directions[i][0];
	    yy = y + directions[i][1];
	    zz = z + directions[i][2];
	  
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		Dn = old(xx,yy,zz);
		temp = (Dn-D)*(factor[i]*lambda/sum);
	        if(first)
		{
		  t=temp;
		  first=false;
		}
		else
		  t=t+temp;
	      }
	  }
	  
	  for (i = 0; i < 13; i++)
	  {
	    xx = x - directions[i][0];
	    yy = y - directions[i][1];
	    zz = z - directions[i][2];
	    
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		Dn = old(xx,yy,zz);		
		temp = (Dn-D)*(factor[i]*lambda/sum);
	        if(first)
		{
		  t=temp;
		  first=false;
		}
		else
		  t=t+temp;
	      }
	  }
	  
	  //adjust the tensor
	  _tensorField[x][y][z] = _tensorField[x][y][z]+t;
	  
	}
}

void irtkTensorField::Smooth(double lambda)
{
  int i, j;
  int x,y,z,xx,yy,zz;
  double sum;
  bool first;
  
  int directions[13][3] =
    {
      { 1, 0, -1 },
      { 0, 1, -1 },
      { 1, 1, -1 },
      { 1, -1, -1 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 1, 1, 0 },
      { 1, -1, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { 1, 1, 1 },
      { 1, -1, 1 },
      { 0, 0, 1 }
    };

  double factor[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  for (i = 0; i < 13; i++)
  {
    for (j = 0; j < 3; j++)
      factor[i] += fabs(double(directions[i][j]));
    factor[i] = 1 / factor[i];
  }
  
  irtkTensorField old = *this;
  irtkTensorField oldLog = *this;
  oldLog.Log();
  irtkVector t,temp,D,Dn;
  
  for (x = 0; x < _x; x++)
    for (y = 0; y < _y; y++)
      for (z = 0; z < _z; z++)
	if(old(x,y,z).IsPD())
        {
	  //Calculate the sum of all coefficients
          sum = 0;
          for (i = 0; i < 13; i++)
          {
	    xx = x + directions[i][0];
	    yy = y + directions[i][1];
	    zz = z + directions[i][2];

	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	        sum += factor[i];
	  }
	  
	  for (i = 0; i < 13; i++)
	  {
	    xx = x - directions[i][0];
	    yy = y - directions[i][1];
	    zz = z - directions[i][2];
	    
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	        sum += factor[i];
	  }
	  
	  //Smooth the image
	  first = true;
	  D=oldLog(x,y,z);
	  //D.Print();
	  
          for (i = 0; i < 13; i++)
          {
	    xx = x + directions[i][0];
	    yy = y + directions[i][1];
	    zz = z + directions[i][2];
	  
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		Dn = oldLog(xx,yy,zz);
		temp = (Dn-D)*(factor[i]*lambda/sum);
	        if(first)
		{
		  t=temp;
		  first=false;
		}
		else
		  t=t+temp;
	      }
	  }

	  for (i = 0; i < 13; i++)
	  {
	    xx = x - directions[i][0];
	    yy = y - directions[i][1];
	    zz = z - directions[i][2];
	    
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		Dn = oldLog(xx,yy,zz);		
		temp = (Dn-D)*(factor[i]*lambda/sum);
	        if(first)
		{
		  t=temp;
		  first=false;
		}
		else
		  t=t+temp;
	      }
	  }
	  D=D+t;
	  //adjust the tensor
	  _tensorField[x][y][z] = D;
	  _tensorField[x][y][z].Exp();
	  
	}
   //We replaced all PD tensor with their new Log, need to Exp.
   //Exp();
   /*
  for (x = 0; x < _x; x++)
    for (y = 0; y < _y; y++)
      for (z = 0; z < _z; z++)
        if(old(x,y,z).IsPD())
          _tensorField[x][y][z].Print();
	*/
}
/*
void irtkTensorField::EdgePreserveSmooth(double lambda, double delta)
{
  int i, j;
  int x,y,z,xx,yy,zz;
  double sum;
  bool first;
  double b[26];
  
  int directions[13][3] =
    {
      { 1, 0, -1 },
      { 0, 1, -1 },
      { 1, 1, -1 },
      { 1, -1, -1 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 1, 1, 0 },
      { 1, -1, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { 1, 1, 1 },
      { 1, -1, 1 },
      { 0, 0, 1 }
    };

  double factor[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  for (i = 0; i < 13; i++)
  {
    for (j = 0; j < 3; j++)
      factor[i] += fabs(double(directions[i][j]));
    factor[i] = 1 / factor[i];
  }
  
  irtkTensorField old = *this;
  irtkTensorField oldLog = *this;
  oldLog.Log();
  irtkVector t,temp,D,Dn;
  
  for (x = 0; x < _x; x++)
    for (y = 0; y < _y; y++)
      for (z = 0; z < _z; z++)
	if(old(x,y,z).IsPD())
        {
	  //Calculate the sum of all coefficients
          sum = 0;
          for (i = 0; i < 13; i++)
          {
	    xx = x + directions[i][0];
	    yy = y + directions[i][1];
	    zz = z + directions[i][2];

	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		double dist = oldLog(x,y,z).EuclDist(oldLog(xx,yy,zz));
	        sum += factor[i];
		//b[i]=
	      }
	  }
	  
	  for (i = 0; i < 13; i++)
	  {
	    xx = x - directions[i][0];
	    yy = y - directions[i][1];
	    zz = z - directions[i][2];
	    
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	        sum += factor[i];
	  }
	  
	  //Smooth the image
	  first = true;
	  D=oldLog(x,y,z);
	  //D.Print();
	  
          for (i = 0; i < 13; i++)
          {
	    xx = x + directions[i][0];
	    yy = y + directions[i][1];
	    zz = z + directions[i][2];
	  
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		Dn = oldLog(xx,yy,zz);
		temp = (Dn-D)*(factor[i]*lambda/sum);
	        if(first)
		{
		  t=temp;
		  first=false;
		}
		else
		  t=t+temp;
	      }
	  }

	  for (i = 0; i < 13; i++)
	  {
	    xx = x - directions[i][0];
	    yy = y - directions[i][1];
	    zz = z - directions[i][2];
	    
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		Dn = oldLog(xx,yy,zz);		
		temp = (Dn-D)*(factor[i]*lambda/sum);
	        if(first)
		{
		  t=temp;
		  first=false;
		}
		else
		  t=t+temp;
	      }
	  }
	  D=D+t;
	  //adjust the tensor
	  _tensorField[x][y][z] = D;
	  _tensorField[x][y][z].Exp();
	  
	}
   //We replaced all PD tensor with their new Log, need to Exp.
   //Exp();
}
*/


void irtkTensorField::EdgePreserveSmooth(double lambda, double delta)
{
  int i, j;
  int x,y,z,xx,yy,zz;
  double sum;
  bool first;
  double b[26];

  
  int directions[13][3] =
    {
      { 1, 0, -1 },
      { 0, 1, -1 },
      { 1, 1, -1 },
      { 1, -1, -1 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { 1, 1, 0 },
      { 1, -1, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { 1, 1, 1 },
      { 1, -1, 1 },
      { 0, 0, 1 }
    };

  double factor[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  for (i = 0; i < 13; i++)
  {
    for (j = 0; j < 3; j++)
      factor[i] += fabs(double(directions[i][j]));
    factor[i] = 1 / factor[i];
  }
  
  irtkTensorField old = *this;
  irtkTensorField oldLog = *this;
  oldLog.Log();
  irtkVector t,temp,D,Dn;
  
  for (x = 0; x < _x; x++)
    for (y = 0; y < _y; y++)
      for (z = 0; z < _z; z++)
	if(old(x,y,z).IsPD())
        {
	  //Calculate the sum of all coefficients
          sum = 0;
          for (i = 0; i < 13; i++)
          {
	    xx = x + directions[i][0];
	    yy = y + directions[i][1];
	    zz = z + directions[i][2];

	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		double dist = oldLog(x,y,z).EuclDist(oldLog(xx,yy,zz));
		double diff = dist/(sqrt(factor[i])*delta);
		b[i]=factor[i]/sqrt(1+diff*diff);
	        sum += factor[i];
	      }
	  }
	  
	  for (i = 0; i < 13; i++)
	  {
	    xx = x - directions[i][0];
	    yy = y - directions[i][1];
	    zz = z - directions[i][2];
	    
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		double dist = oldLog(x,y,z).EuclDist(oldLog(xx,yy,zz));
		double diff = dist/(sqrt(factor[i])*delta);
		b[i+13]=factor[i]/sqrt(1+diff*diff);
	        sum += factor[i];
	      }
	  }
	  
	  //Smooth the image
	  first = true;
	  D=oldLog(x,y,z);
	  //D.Print();
	  
          for (i = 0; i < 13; i++)
          {
	    xx = x + directions[i][0];
	    yy = y + directions[i][1];
	    zz = z + directions[i][2];
	  
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		Dn = oldLog(xx,yy,zz);
		temp = (Dn-D)*(b[i]*lambda/sum);
	        if(first)
		{
		  t=temp;
		  first=false;
		}
		else
		  t=t+temp;
	      }
	  }

	  for (i = 0; i < 13; i++)
	  {
	    xx = x - directions[i][0];
	    yy = y - directions[i][1];
	    zz = z - directions[i][2];
	    
	    if ((xx >= 0) && (xx < _x) && (yy >= 0) && (yy < _y) && (zz >= 0) && (zz < _z))
	      if(old(xx,yy,zz).IsPD())
	      {
		Dn = oldLog(xx,yy,zz);		
		temp = (Dn-D)*(b[i]*lambda/sum);
	        if(first)
		{
		  t=temp;
		  first=false;
		}
		else
		  t=t+temp;
	      }
	  }
	  D=D+t;
	  //adjust the tensor
	  _tensorField[x][y][z] = D;
	  _tensorField[x][y][z].Exp();
	  
	}
}
 void irtkTensorField::SaveFA(irtkRealImage b0, char * name)
 {
   if((b0.GetX()!=_x)||(b0.GetY()!=_y)||(b0.GetZ()!=_z))
   {
     cout<<"irtkTensorField::Save(): Please give b0 of correct directions!"<<endl;
     exit(1);
  }
  
  for(int i=0;i<_x;i++)
    for(int j=0;j<_y;j++)
      for(int k=0;k<_z;k++)
	if(_tensorField[i][j][k].IsPD())
	  b0(i,j,k)=_tensorField[i][j][k].FA();
	else
	  b0(i,j,k)=-1;
  b0.Write(name);
}

 void irtkTensorField::Save(irtkRealImage b0, char * name)
 {
   if((b0.GetX()!=_x)||(b0.GetY()!=_y)||(b0.GetZ()!=_z))
   {
     cout<<"irtkTensorField::Save(): Please give b0 of correct directions!"<<endl;
     exit(1);
  }
  
  irtkImageAttributes attr = b0.GetImageAttributes();
  attr._t=7;
  irtkRealImage dti(attr);
  
  for(int i=0;i<_x;i++)
    for(int j=0;j<_y;j++)
      for(int k=0;k<_z;k++)
      {
	dti(i,j,k,0)=b0(i,j,k);
	for(int m=1;m<7;m++)
	  dti(i,j,k,m)=_tensorField[i][j][k](m);
      }
  dti.Write(name);
}

void irtkTensorField::SaveDirections(irtkRealImage b0, char * name)
 {
   if((b0.GetX()!=_x)||(b0.GetY()!=_y)||(b0.GetZ()!=_z))
   {
     cout<<"irtkTensorField::Save(): Please give b0 of correct directions!"<<endl;
     exit(1);
  }
  
  irtkImageAttributes attr = b0.GetImageAttributes();
  attr._t=3;
  irtkRealImage PD(attr);
  irtkVector v;
  
  for(int i=0;i<_x;i++)
    for(int j=0;j<_y;j++)
      for(int k=0;k<_z;k++)
	if(_tensorField[i][j][k].IsPD())
	{
	  v=_tensorField[i][j][k].Direction();
	  for(int m=0;m<3;m++)
	    PD(i,j,k,m)=v(m);
	}
	else
	{
	  for(int m=0;m<3;m++)
	    PD(i,j,k,m)=0;  
	}
  PD.Write(name);
}