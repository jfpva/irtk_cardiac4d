/*=========================================================================

 Library   : Image Registration Toolkit (IRTK)
 Module    : $Id: irtkLaplacianSmoothing.cc 837 2013-05-07 12:55:31Z mm3 $
 Copyright : Imperial College, Department of Computing
 Visual Information Processing (VIP), 2011 onwards
 Date      : $Date: 2013-05-07 13:55:31 +0100 (Tue, 07 May 2013) $
 Version   : $Revision: 837 $
 Changes   : $Author: mm3 $

 =========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkGaussianBlurring.h>
#include <irtkGaussianBlurringWithPadding.h>
#include <irtkResampling.h>
#include <irtkResamplingWithPadding.h>
#include <irtkLaplacianSmoothing.h>


irtkLaplacianSmoothing::irtkLaplacianSmoothing()
{
  int directions[13][3] = {
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
        { 0, 0, 1 }};
	
  for(int i=0;i<13;i++)
  {
    for(int j=0;j<3;j++)
    {
      _directions[i][j]=directions[i][j];
    }
  }
    
  _lap_threshold=0.000001;
  _rel_diff_threshold=0.0001;
  _relax_iter=5;
  _boundary_weight=0.4;
}

void irtkLaplacianSmoothing::InitializeFactors()
{
    int i,j;
    
    for (i = 0; i < 13; i++)
      _factor.push_back(0);
    
    double sum=0;
    for (i = 0; i < 13; i++) {
        for (j = 0; j < 3; j++)
            _factor[i] += fabs(double(_directions[i][j]));
        _factor[i] = sqrt(1 / _factor[i]);
	sum += _factor[i];
    }
    sum*=2;
    for (i = 0; i < 13; i++) {
      _factor[i]/=sum;
    }

}

void irtkLaplacianSmoothing::Blur(irtkRealImage& image, double sigma, double padding)
{
  irtkGaussianBlurringWithPadding<irtkRealPixel> gb(sigma,padding);
  gb.SetInput(&image);
  gb.SetOutput(&image);
  gb.Run();
}

void irtkLaplacianSmoothing::Resample(irtkRealImage& image, double step, double padding)
{
  irtkResamplingWithPadding<irtkRealPixel> res(step, step, step, padding);
  irtkLinearInterpolateImageFunction interpolator;
  res.SetInput(&image);
  res.SetOutput(&image);
  res.SetInterpolator(&interpolator);
  res.Run();
}

void irtkLaplacianSmoothing::ResampleOnGrid(irtkRealImage& image, irtkRealImage templ)
{
  irtkImageTransformation *imagetransformation = new irtkImageTransformation;
  irtkImageFunction *interpolator = new irtkNearestNeighborInterpolateImageFunction;
  imagetransformation->PutInterpolator(interpolator);
  //imagetransformation->PutTargetPaddingValue(-1);
  //imagetransformation->PutSourcePaddingValue(-1);
  irtkRigidTransformation id;
  templ=0;
  imagetransformation->SetInput(&image, &id);
  imagetransformation->SetOutput(&templ);
  imagetransformation->Run();
  image = templ;
}

void irtkLaplacianSmoothing::SetMask(irtkRealImage mask)
{
  ResampleOnGrid(mask,_image);
  _mask = mask;
  _mask.Write("_mask_LaplacianSmoothing.nii.gz");
}

void irtkLaplacianSmoothing::EnlargeImage(irtkRealImage &image)
{
  int i,j,k,l;
  irtkRealImage original = image;
  irtkImageAttributes attr = image.GetImageAttributes();
  attr._x += 2;
  attr._y += 2;
  attr._z += 2;
  image.Initialize(attr);
  image=0;
  
  for(l=0;l<original.GetT();l++)
    for(i=0; i<original.GetX();i++)
      for(j=0; j<original.GetY();j++)
        for(k=0; k<original.GetZ();k++)
        {
          image.Put(i+1,j+1,k+1,l,original.Get(i,j,k,l));
        }
}


void irtkLaplacianSmoothing::ReduceImage(irtkRealImage &image)
{
  int i,j,k,l;
  irtkRealImage original = image;
  irtkImageAttributes attr = image.GetImageAttributes();
  attr._x -= 2;
  attr._y -= 2;
  attr._z -= 2;
  image.Initialize(attr);
  image=0;
  
  for(l=0;l<original.GetT();l++)
    for(i=1; i<(original.GetX()-1);i++)
      for(j=1; j<(original.GetY()-1);j++)
        for(k=1; k<(original.GetZ()-1);k++)
        {
          image.Put(i-1,j-1,k-1,l,original.Get(i,j,k,l));
        }
}

void irtkLaplacianSmoothing::CalculateBoundary(irtkRealImage& m)
{
    int x, y, z, xx, yy, zz,i;

    int dx = m.GetX();
    int dy = m.GetY();
    int dz = m.GetZ();
    
    for (x = 0; x < dx; x++)
      for (y = 0; y < dy; y++)
        for (z = 0; z < dz; z++) 
	  if(m(x,y,z)==0)
	  {
            for (i = 0; i < 13; i++) {
              xx = x + _directions[i][0];
              yy = y + _directions[i][1];
              zz = z + _directions[i][2];
              if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	        if (m(xx,yy,zz)==1)
		  m(x,y,z)=2;
	    }

            for (i = 0; i < 13; i++) {
              xx = x - _directions[i][0];
              yy = y - _directions[i][1];
              zz = z - _directions[i][2];
              if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	        if (m(xx,yy,zz)==1)
		  m(x,y,z)=2;
	    }
          }
}

void irtkLaplacianSmoothing::CalculateBoundaryWeights(irtkRealImage& w, irtkRealImage m)
{
    int x, y, z, xx, yy, zz,i;
    double sum;

    int dx = m.GetX();
    int dy = m.GetY();
    int dz = m.GetZ();
    
    for (x = 0; x < dx; x++)
      for (y = 0; y < dy; y++)
        for (z = 0; z < dz; z++) 
	  if(m(x,y,z)==2)
	  {
	    sum=0;
            for (i = 0; i < 13; i++) {
              xx = x + _directions[i][0];
              yy = y + _directions[i][1];
              zz = z + _directions[i][2];
              if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	        if (m(xx,yy,zz)==2)
		  sum += _factor[i];
	    }

            for (i = 0; i < 13; i++) {
              xx = x - _directions[i][0];
              yy = y - _directions[i][1];
              zz = z - _directions[i][2];
              if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	        if (m(xx,yy,zz)==2)
		  sum += _factor[i];
	    }
	    w(x,y,z)=sum;
          }
}



double irtkLaplacianSmoothing::LaplacianImage(irtkRealImage image, irtkRealImage m, irtkRealImage& laplacian)
{
    //cout << "Diffusion Regularization" << endl;
    
    irtkRealImage original = image;
    laplacian=0;

    int dx = image.GetX();
    int dy = image.GetY();
    int dz = image.GetZ();

    int x, y, z, xx, yy, zz;
    double val,sum,diff;
    double laplacian_value=0;
    int i,count=0;
            
    for (x = 0; x < dx; x++)
      for (y = 0; y < dy; y++)
        for (z = 0; z < dz; z++) 
	  if(m(x,y,z)==1)
	  {                    
            val = 0;
            sum = 0;
            for (i = 0; i < 13; i++) {
              xx = x + _directions[i][0];
              yy = y + _directions[i][1];
              zz = z + _directions[i][2];
              if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	      {
                val += _factor[i] * original(xx, yy, zz);
	        sum += _factor[i];
              }
            }

            for (i = 0; i < 13; i++) {
              xx = x - _directions[i][0];
              yy = y - _directions[i][1];
              zz = z - _directions[i][2];
              if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	      {
                val += _factor[i] * original(xx, yy, zz);
	        sum += _factor[i];
              }
            }
          
            if(sum>0.99)//we are inside because all neigbours exist
            {
	      diff = original(x,y,z)-val;
	      laplacian(x,y,z)=diff;
	      laplacian_value+=0.5*diff*diff;
	      count++;
	    }
	    else
	    {
	      cerr<<"Warning: LaplacianImage: sum<=0.99, this should not happen!";
	    }
	  }
    cout<<"Average Laplacian is "<<laplacian_value/count<<" ";
    return laplacian_value/count;
}

double irtkLaplacianSmoothing::LaplacianBoundary(irtkRealImage image, irtkRealImage m, irtkRealImage& laplacian)
{
    irtkRealImage original = image;

    int dx = image.GetX();
    int dy = image.GetY();
    int dz = image.GetZ();

    int x, y, z, xx, yy, zz;
    double val,sum,diff;
    double laplacian_value=0;
    int i,count=0;
            
    for (x = 0; x < dx; x++)
      for (y = 0; y < dy; y++)
        for (z = 0; z < dz; z++) 
	  if(m(x,y,z)==2)
	  {                    
            val = 0;
            sum = 0;
            for (i = 0; i < 13; i++) {
              xx = x + _directions[i][0];
              yy = y + _directions[i][1];
              zz = z + _directions[i][2];
              if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
		if(m(xx,yy,zz)==2)
	        {
                  val += _factor[i] * original(xx, yy, zz);
	          sum += _factor[i];
                }
            }

            for (i = 0; i < 13; i++) {
              xx = x - _directions[i][0];
              yy = y - _directions[i][1];
              zz = z - _directions[i][2];
              if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
		if(m(xx,yy,zz)==2)
	        {
                  val += _factor[i] * original(xx, yy, zz);
	          sum += _factor[i];
                }
            }
          
            if(sum>0)//we are inside because all neigbours exist
            {
	      val/=sum;
	      diff = original(x,y,z)-val;
	      laplacian(x,y,z)=diff;
	      laplacian_value+=0.5*diff*diff;
	      count++;
	    }
	    else
	    {
	      cerr<<"Warning: LaplacianBoundary: sum=0, this should not happen!";
	    }
	  }
    cout<<"Average Laplacian on boundary is "<<laplacian_value/count<<" ";
    return laplacian_value/count;
}

void irtkLaplacianSmoothing::UpdateFieldmapGD(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage laplacian, irtkRealImage mask, irtkRealImage weights, double alpha, double lambda1, double lambda2)
{
    
  irtkRealImage original = fieldmap;
  fieldmap = 0;

  int dx = image.GetX();
  int dy = image.GetY();
  int dz = image.GetZ();

  int x, y, z, xx, yy, zz,i;
  double sum;
   
  for (x = 0; x < dx; x++)
    for (y = 0; y < dy; y++)
      for (z = 0; z < dz; z++) 
	if(mask(x,y,z)>0)
	{
          sum = 0;
          for (i = 0; i < 13; i++) {
            xx = x + _directions[i][0];
            yy = y + _directions[i][1];
            zz = z + _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += lambda1*_factor[i]*laplacian(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += lambda2*_factor[i]*laplacian(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          for (i = 0; i < 13; i++) {
            xx = x - _directions[i][0];
            yy = y - _directions[i][1];
            zz = z - _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += lambda1*_factor[i]*laplacian(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += lambda2*_factor[i]*laplacian(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          
          if(mask(x,y,z)==1)
	    fieldmap(x,y,z)=(1-alpha)*original(x,y,z)+alpha*image(x,y,z)-alpha*lambda1*laplacian(x,y,z)+alpha*sum;
	  if(mask(x,y,z)==2)
            fieldmap(x,y,z)=original(x,y,z)-alpha*lambda2*laplacian(x,y,z)+alpha*sum;  
	}    
}


void irtkLaplacianSmoothing::LLaplacian(irtkRealImage& llaplacian, irtkRealImage laplacian, irtkRealImage mask, irtkRealImage weights)
{
    
  llaplacian = 0;

  int dx = llaplacian.GetX();
  int dy = llaplacian.GetY();
  int dz = llaplacian.GetZ();

  int x, y, z, xx, yy, zz,i;
  double sum;
   
  for (x = 0; x < dx; x++)
    for (y = 0; y < dy; y++)
      for (z = 0; z < dz; z++) 
	if(mask(x,y,z)>0)
	{
          sum = 0;
          for (i = 0; i < 13; i++) {
            xx = x + _directions[i][0];
            yy = y + _directions[i][1];
            zz = z + _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += _factor[i]*laplacian(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += _factor[i]*laplacian(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          for (i = 0; i < 13; i++) {
            xx = x - _directions[i][0];
            yy = y - _directions[i][1];
            zz = z - _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += _factor[i]*laplacian(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += _factor[i]*laplacian(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          
	    llaplacian(x,y,z)=laplacian(x,y,z)-sum;
	}    
}


void irtkLaplacianSmoothing::UpdateFieldmap(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage multipliers, irtkRealImage mask, irtkRealImage weights, double alpha)
{
    
  irtkRealImage original = fieldmap;
  fieldmap = 0;

  int dx = image.GetX();
  int dy = image.GetY();
  int dz = image.GetZ();

  int x, y, z, xx, yy, zz,i;
  double sum;
   
  for (x = 0; x < dx; x++)
    for (y = 0; y < dy; y++)
      for (z = 0; z < dz; z++) 
	if(mask(x,y,z)>0)
	{
          sum = 0;
          for (i = 0; i < 13; i++) {
            xx = x + _directions[i][0];
            yy = y + _directions[i][1];
            zz = z + _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          for (i = 0; i < 13; i++) {
            xx = x - _directions[i][0];
            yy = y - _directions[i][1];
            zz = z - _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          
          if(mask(x,y,z)==1)
	    fieldmap(x,y,z)=(1-alpha)*original(x,y,z)+alpha*image(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;
	  if(mask(x,y,z)==2)
            fieldmap(x,y,z)=original(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;  
	}    
}


void irtkLaplacianSmoothing::UpdateFieldmapWithThreshold(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage multipliers, 
		    irtkRealImage mask, irtkRealImage weights, double alpha, double threshold)
{
    
  irtkRealImage original = fieldmap;
  irtkRealImage tr_mask = fieldmap;
  tr_mask = 0;
  fieldmap = 0;

  int dx = image.GetX();
  int dy = image.GetY();
  int dz = image.GetZ();

  int x, y, z, xx, yy, zz,i;
  double sum,diff;
   
  for (x = 0; x < dx; x++)
    for (y = 0; y < dy; y++)
      for (z = 0; z < dz; z++) 
	if(mask(x,y,z)>0)
	{
          sum = 0;
          for (i = 0; i < 13; i++) {
            xx = x + _directions[i][0];
            yy = y + _directions[i][1];
            zz = z + _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          for (i = 0; i < 13; i++) {
            xx = x - _directions[i][0];
            yy = y - _directions[i][1];
            zz = z - _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          
          diff=original(x,y,z)-image(x,y,z);
          if(mask(x,y,z)==1)
	  {
	    if(fabs(diff)<=threshold)
	    {
	      fieldmap(x,y,z)=(1-alpha)*original(x,y,z)+alpha*image(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;
	      tr_mask(x,y,z)=1;
	    }
	    else
	    {
	      if(diff>threshold)
	        fieldmap(x,y,z)=original(x,y,z)-alpha*threshold-alpha*multipliers(x,y,z)+alpha*sum;
	      if(diff<-threshold)
	        fieldmap(x,y,z)=original(x,y,z)+alpha*threshold-alpha*multipliers(x,y,z)+alpha*sum;
	    }
	  }
	  if(mask(x,y,z)==2)
            fieldmap(x,y,z)=original(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;  
	}    
	tr_mask.Write("tr_mask.nii.gz");
}

void irtkLaplacianSmoothing::UpdateFieldmapHuber(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage multipliers, 
		    irtkRealImage mask, irtkRealImage weights, double alpha)
{
    
  irtkRealImage original = fieldmap;
  irtkRealImage tr_mask = fieldmap;
  tr_mask = 0;
  fieldmap = 0;

  int dx = image.GetX();
  int dy = image.GetY();
  int dz = image.GetZ();

  int x, y, z, xx, yy, zz,i;
  double sum,diff,threshold;
  
  std::vector<double> errors;
  
  //calculate median
  for (x = 0; x < dx; x++)
    for (y = 0; y < dy; y++)
      for (z = 0; z < dz; z++) 
	if(mask(x,y,z)==1)
	{
	  errors.push_back(fabs(image(x,y,z)-original(x,y,z)));
	}
  std::sort(errors.begin(),errors.end());
  threshold =1.35*errors[round(errors.size()/2)];
  cout<<"Threshold is "<<threshold<<endl;
  
  //update fieldmap
  for (x = 0; x < dx; x++)
    for (y = 0; y < dy; y++)
      for (z = 0; z < dz; z++) 
	if(mask(x,y,z)>0)
	{
          sum = 0;
          for (i = 0; i < 13; i++) {
            xx = x + _directions[i][0];
            yy = y + _directions[i][1];
            zz = z + _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          for (i = 0; i < 13; i++) {
            xx = x - _directions[i][0];
            yy = y - _directions[i][1];
            zz = z - _directions[i][2];
            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
	    {
	      if (mask(xx,yy,zz)==1)
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz);
              }
	      if ((mask(xx,yy,zz)==2)&&(mask(x,y,z)==2))
	      {
	        sum += _factor[i]*multipliers(xx,yy,zz)/weights(x,y,z);
              }
	    }
          }
          
          diff=original(x,y,z)-image(x,y,z);
          if(mask(x,y,z)==1)
	  {
	    if(fabs(diff)<=threshold)
	    {
	      fieldmap(x,y,z)=(1-alpha)*original(x,y,z)+alpha*image(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;
	      tr_mask(x,y,z)=1;
	    }
	    else
	    {
	      if(diff>threshold)
	        fieldmap(x,y,z)=original(x,y,z)-alpha*threshold-alpha*multipliers(x,y,z)+alpha*sum;
	      if(diff<-threshold)
	        fieldmap(x,y,z)=original(x,y,z)+alpha*threshold-alpha*multipliers(x,y,z)+alpha*sum;
	    }
	  }
	  if(mask(x,y,z)==2)
            fieldmap(x,y,z)=original(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;  
	}    
	tr_mask.Write("tr_mask.nii.gz");
}

void irtkLaplacianSmoothing::Smooth(irtkRealImage& im, irtkRealImage m)
{
  double prev_lap, lap;
  char buffer[255];
  
  im.Write("im.nii.gz");
  m.Write("m.nii.gz");
  
  irtkRealImage image(im);
  EnlargeImage(image);
  
  irtkRealImage mask(m);
  EnlargeImage(mask);
  CalculateBoundary(mask);

  image.Write("image.nii.gz");
  m.Write("mask.nii.gz");

  
  irtkRealImage weights(mask);
  weights=0;
  CalculateBoundaryWeights(weights,mask);
  
  irtkRealImage fieldmap(image);
  irtkRealImage laplacian(image);
  irtkRealImage llaplacian(image);
  irtkRealImage multipliers(image);
  laplacian=0;
  llaplacian=0;
  multipliers=0;
  image.Write("smooth-image.nii.gz");
  mask.Write("smooth-mask.nii.gz");

  double alpha = 0.5;
  double rel_diff;
  
  for (int iter = 1; iter<1000; iter++)
    {	
      prev_lap=lap;
      lap = LaplacianImage(fieldmap, mask,laplacian);
      laplacian.Write("laplacian.nii.gz");
      LaplacianBoundary(fieldmap, mask,laplacian);
      laplacian.Write("laplacian-with-boundary.nii.gz");
      //LLaplacian(llaplacian,laplacian,mask,weights);
      //llaplacian.Write("llaplacian.nii.gz");
      //exit(1);
      multipliers = multipliers + laplacian*alpha;
      multipliers.Write("multipliers.nii.gz");
      
      //UpdateFieldmapGD(fieldmap, image, laplacian, mask, weights, alpha,1,1);
      UpdateFieldmap(fieldmap, image, laplacian, mask, weights, alpha);
      //UpdateFieldmapWithThreshold(fieldmap, image, multipliers, mask, weights, alpha,0.28);
      //UpdateFieldmapHuber(fieldmap, image, multipliers, mask, weights, alpha);

      //sprintf(buffer,"fieldmap%i.nii.gz",iter);
      //fieldmap.Write(buffer);
      fieldmap.Write("fieldmap.nii.gz");
      
      rel_diff = (prev_lap-lap)/prev_lap;
      cout<<"Iter "<<iter<<" alpha "<<alpha<<" prev_lap "<<prev_lap<<" lap "<<lap<<" reldiff "<<rel_diff<<endl;
      if((iter>1)&&(rel_diff<0.0001)) break;
    }
    //reduce fieldmap to original size
    ReduceImage(fieldmap);
    //mask out all the boundary voxels
    im = fieldmap*m;
    
}

void irtkLaplacianSmoothing::SmoothGD(irtkRealImage& im, irtkRealImage m)
{
  double prev_lap, lap;
  char buffer[255];
  
  //im.Write("im.nii.gz");
  //m.Write("m.nii.gz");
  
  irtkRealImage image(im);
  EnlargeImage(image);
  
  irtkRealImage mask(m);
  EnlargeImage(mask);
  CalculateBoundary(mask);

  //image.Write("image.nii.gz");
  //m.Write("mask.nii.gz");

  
  irtkRealImage weights(mask);
  weights=0;
  CalculateBoundaryWeights(weights,mask);
  
  irtkRealImage fieldmap(image);
  irtkRealImage laplacian(image);
  irtkRealImage llaplacian(image);
  irtkRealImage multipliers(image);
  laplacian=0;
  llaplacian=0;
  multipliers=0;
  //image.Write("smooth-image.nii.gz");
  //mask.Write("smooth-mask.nii.gz");

  double alpha = 0.5;
  double rel_diff;
    alpha = 5;
    for(int aiter = 1; aiter<_relax_iter; aiter++)
    {
      alpha/=10;
      for (int iter = 1; iter<10000; iter++)
      {	
        prev_lap=lap;
        lap = LaplacianImage(fieldmap, mask,laplacian);
        //laplacian.Write("laplacian.nii.gz");
        LaplacianBoundary(fieldmap, mask,laplacian);
        //laplacian.Write("laplacian-with-boundary.nii.gz");
	///Change
        UpdateFieldmapGD(fieldmap, image, laplacian, mask, weights, alpha,0.5/alpha,_boundary_weight*0.5/alpha);
        //UpdateFieldmapGD(fieldmap, image, laplacian, mask, weights, alpha,0.5/alpha,1);
        //UpdateFieldmapGD(fieldmap, image, laplacian, mask, weights, alpha,0.5/alpha,sqrt(0.5/alpha));
        //UpdateFieldmapGD(fieldmap, image, laplacian, mask, weights, alpha,0.5/alpha,0.5/alpha);
        //UpdateFieldmapGD(fieldmap, image, laplacian, mask, weights, alpha,0.5/alpha,0.05/alpha);

        //sprintf(buffer,"fieldmap%f-%i.nii.gz",alpha,iter);
        //fieldmap.Write(buffer);
      
        rel_diff = (prev_lap-lap)/prev_lap;
        cout<<"Iter "<<iter<<" alpha "<<alpha<<" prev_lap "<<prev_lap<<" lap "<<lap<<" reldiff "<<rel_diff<<endl;
        if((iter>1)&&((lap<_lap_threshold)||(rel_diff<_rel_diff_threshold))) break;
      }
      //sprintf(buffer,"final-fieldmap%f.nii.gz",alpha);
      //fieldmap.Write(buffer);
    }
    
    //reduce fieldmap to original size
    ReduceImage(fieldmap);
    //mask out all the boundary voxels
    im = fieldmap*m;
    
}


void irtkLaplacianSmoothing::UpsampleFieldmap(irtkRealImage& target, irtkRealImage mask, irtkRealImage &newmask)
{
  irtkRealImage f = _fieldmap;
  _fieldmap = target;
  _fieldmap = 0;

  //add padding to f
  irtkRealPixel mn,mx;
  f.GetMinMax(&mn,&mx);
  double padding = round(mn-2);
  cout<<"padding = "<<padding<<endl;
  
  irtkRealPixel *pf = f.GetPointerToVoxels();
  irtkRealPixel *pm = mask.GetPointerToVoxels();
  for (int i=0; i<f.GetNumberOfVoxels();i++)
  {
    if(*pm!=1)
      *pf = padding;
    pm++;
    pf++;
  }
  
  f.Write("f.nii.gz");
  
  irtkLinearInterpolateImageFunction interpolator;
  interpolator.SetInput(&f);
  interpolator.Initialize();

  double x,y,z;
  double value;
  irtkImageAttributes attr = _fieldmap.GetImageAttributes();
  for (int t = 0; t < _fieldmap.GetT(); t++) {
    for (int k = 0; k < _fieldmap.GetZ(); k++) {
      for (int j = 0; j < _fieldmap.GetY(); j++) {
        for (int i = 0; i < _fieldmap.GetX(); i++) {
	  if(newmask(i,j,k)==1)
	  {
            x = i;
            y = j;
            z = k;
	    _fieldmap.ImageToWorld(x,y,z);
	    f.WorldToImage(x,y,z);
	    /*
	  if ((x > -0.5) && (x < f.GetX()-0.5) && 
	      (y > -0.5) && (y < f.GetY()-0.5) &&
              (z > -0.5) && (z < f.GetZ()-0.5))
              *//*
	  if ((x > 0) && (x < f.GetX()-1) && 
	      (y > 0) && (y < f.GetY()-1) &&
              (z > 0) && (z < f.GetZ()-1))
              */
	    if ((x > -1) && (x < f.GetX()) && 
	        (y > -1) && (y < f.GetY()) &&
                (z > -1) && (z < f.GetZ()))
	    {
	      value = interpolator.EvaluateWithPadding(x,y,z,t,padding);
	      if(value>padding)
	       _fieldmap(i,j,k,t) = value;
	      else
	      {
		_fieldmap(i,j,k,t) = 0;
		newmask(i,j,k,t) = 0;
	      }
	    //_fieldmap(i,j,k,t) = interpolator.Evaluate(x,y,z,t);
	    }
	    else
	      _fieldmap(i,j,k,t) = 0;
	  }
	  else
	    _fieldmap(i,j,k,t) = 0;
        }
      }
    }
  }  
}

irtkRealImage irtkLaplacianSmoothing::Run()
{
     InitializeFactors();
    
    double step;
    irtkImageAttributes attr = _image.GetImageAttributes();
    irtkRealImage image(_image), mask(_mask);
    //step = attr._dz*8;
    step=6;
    image.Write("image-before.nii.gz");
    Blur(image,step/2,-1000);
    image.Write("image-blurred.nii.gz");
    Resample(image,step,-1000);
    image.Write("image-res.nii.gz");
    mask.Write("mask-before.nii.gz");
    ResampleOnGrid(mask,image);
    image.Write("image-after.nii.gz");
    mask.Write("mask-after.nii.gz");


    
   
    Smooth(image,mask);
    
    image.Write("fieldmap-lr.nii.gz");
    
    _fieldmap=image;
    _m=mask;
    
    image = _image;
    mask=_mask;
    //step = attr._dz*4;
    step=3;
    image.Write("image-before-2.nii.gz");
    Blur(image,step/2,-1000);
    image.Write("image-blurred-2.nii.gz");
    Resample(image,step,-1000);
    image.Write("image-res-2.nii.gz");
    mask.Write("mask-before-2.nii.gz");
    ResampleOnGrid(mask,image);
    image.Write("image-after-2.nii.gz");
    mask.Write("mask-after-2.nii.gz");

    
    UpsampleFieldmap(image,_m,mask);
    _fieldmap.Write("upsampled.nii.gz");
    image*=mask;
    image-=_fieldmap;
    
    Smooth(image,mask);
     image.Write("fieldmap-hr.nii.gz");
    _fieldmap += image;
    _m=mask;
    _fieldmap.Write("fieldmap-added.nii.gz");
    _image.Write("_image.nii.gz");
    _m.Write("_m.nii.gz");
    _mask.Write("_mask.nii.gz");
    UpsampleFieldmap(_image,_m,_mask);
    image.Write("upsampled-2.nii.gz");
    return _fieldmap;
}

irtkRealImage irtkLaplacianSmoothing::RunGD()
{
     InitializeFactors();
    
    double step;
    irtkImageAttributes attr = _image.GetImageAttributes();
    irtkRealImage image(_image), mask(_mask);
    //step = attr._dz*8;
    step=3;
    //image.Write("image-before.nii.gz");
    Blur(image,step/2,-1000);
    //image.Write("image-blurred.nii.gz");
    Resample(image,step,-1000);
    //image.Write("image-res.nii.gz");
   // mask.Write("mask-before.nii.gz");
    ResampleOnGrid(mask,image);
    //image.Write("image-after.nii.gz");
    //mask.Write("mask-after.nii.gz");


    
   
    SmoothGD(image,mask);

     //image.Write("fieldmap-hr.nii.gz");
    _fieldmap = image;
    _m=mask;
    //_fieldmap.Write("fieldmap-added.nii.gz");
    //_image.Write("_image.nii.gz");
    //_m.Write("_m.nii.gz");
    //_mask.Write("_mask.nii.gz");
    UpsampleFieldmap(_image,_m,_mask);
    //_mask.Write("final-fieldmap-mask.nii.gz");
    //image.Write("upsampled-2.nii.gz");
    return _fieldmap;
}


irtkRealImage irtkLaplacianSmoothing::Run1level()
{
     InitializeFactors();
    
    double step;
    irtkImageAttributes attr = _image.GetImageAttributes();
    irtkRealImage image(_image), mask(_mask);
    Smooth(image,mask);
    _fieldmap=image;
    return _fieldmap;
}

irtkRealImage irtkLaplacianSmoothing::Run3levels()
{
     InitializeFactors();
    
    double step;
    irtkImageAttributes attr = _image.GetImageAttributes();
    irtkRealImage image(_image), mask(_mask);
    step = attr._dz*8;
    step = 6;
    Blur(image,step/2,-1000);
    Resample(image,step,-1000);
    Resample(mask,step,0);
    
    Smooth(image,mask);
    _fieldmap=image;
    _m=mask;
    //exit(1);
    
    
    image = _image;
    mask=_mask;
    step = attr._dz*4;
    step = 3;
    Blur(image,step/2,-1000);
    Resample(image,step,-1000);
    Resample(mask,step,0);
    UpsampleFieldmap(image,_m,mask);
    
    image*=mask;
    image-=_fieldmap;
    
    Smooth(image,mask);
    _fieldmap += image;
    _m=mask;
    
    
    image = _image;
    mask=_mask;
    step = attr._dz*2;
    step = 1.5;
    Blur(image,step/2,-1000);
    Resample(image,step,-1000);
    Resample(mask,step,0);

    
    UpsampleFieldmap(image,_m,mask);
    image*=mask;
    image-=_fieldmap;
    
    Smooth(image,mask);
    _fieldmap += image;
    _m=mask;
    
    
    
    UpsampleFieldmap(_image,_m,_mask);
     //_fieldmap.Write("final-fieldmap.nii.gz");
    return _fieldmap;
  
}