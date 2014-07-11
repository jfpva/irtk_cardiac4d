/*=========================================================================

 Library   : Image Registration Toolkit (IRTK)
 Module    : $Id: irtkLaplacianSmoothingMissingData.cc 837 2013-05-07 12:55:31Z mm3 $
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
#include <irtkLaplacianSmoothingMissingData.h>


irtkLaplacianSmoothingMissingData::irtkLaplacianSmoothingMissingData()
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
      cout<<_directions[i][j]<<" ";
    }
    cout<<";";
  }
    
    
}

void irtkLaplacianSmoothingMissingData::InitializeFactors()
{
   // cout<<"Factors: ";
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
      //cout<<_factor[i]<<" ";
    }
    //cout<<endl;

}

void irtkLaplacianSmoothingMissingData::Blur(irtkRealImage& image, double sigma, irtkRealImage& mask, double padding)
{
  irtkRealPixel *pi = image.GetPointerToVoxels();
  irtkRealPixel *pm = mask.GetPointerToVoxels();
  for(int i=0;i<image.GetNumberOfVoxels();i++)
  {
    if(*pm==0)
      *pi=-1000;
    pm++;
    pi++;
  }
  image.Write("image-with-padding.nii.gz");
  irtkGaussianBlurringWithPadding<irtkRealPixel> gb(sigma,padding);
  gb.SetInput(&image);
  gb.SetOutput(&image);
  gb.Run();
}

void irtkLaplacianSmoothingMissingData::Resample(irtkRealImage& image, double step, double padding)
{
  irtkResamplingWithPadding<irtkRealPixel> res(step, step, step, padding);
  irtkLinearInterpolateImageFunction interpolator;
  res.SetInput(&image);
  res.SetOutput(&image);
  res.SetInterpolator(&interpolator);
  res.Run();
}

void irtkLaplacianSmoothingMissingData::SetPaddingToZero(irtkRealImage& image, double padding)
{
  irtkRealPixel *pi = image.GetPointerToVoxels();
  for(int i=0;i<image.GetNumberOfVoxels();i++)
  {
    if(*pi==-1000)
      *pi=-0;
    pi++;
  }
  image.Write("image-res-without-padding.nii.gz");
}

void irtkLaplacianSmoothingMissingData::EnlargeImage(irtkRealImage &image)
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


void irtkLaplacianSmoothingMissingData::ReduceImage(irtkRealImage &image)
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

void irtkLaplacianSmoothingMissingData::CalculateBoundary(irtkRealImage& m)
{
    int x, y, z, xx, yy, zz,i;
    double sum;

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

void irtkLaplacianSmoothingMissingData::CalculateBoundaryWeights(irtkRealImage& w, irtkRealImage m)
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



double irtkLaplacianSmoothingMissingData::LaplacianImage(irtkRealImage image, irtkRealImage m, irtkRealImage& laplacian)
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
	      laplacian_value+=fabs(diff);
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

double irtkLaplacianSmoothingMissingData::LaplacianBoundary(irtkRealImage image, irtkRealImage m, irtkRealImage& laplacian)
{
    //cout << "Diffusion Regularization" << endl;
    
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
	      laplacian_value+=fabs(diff);
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


void irtkLaplacianSmoothingMissingData::UpdateFieldmap(irtkRealImage& fieldmap, irtkRealImage image, irtkRealImage multipliers, 
		    irtkRealImage mask, irtkRealImage image_mask, irtkRealImage weights, double alpha)
{
    //cout << "Diffusion Regularization" << endl;
    
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
          
          if((mask(x,y,z)==1)&&(image_mask(x,y,z)==1))
	    fieldmap(x,y,z)=(1-alpha)*original(x,y,z)+alpha*image(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;
          if((mask(x,y,z)==1)&&(image_mask(x,y,z)==0))
	    fieldmap(x,y,z)=original(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;
	  if(mask(x,y,z)==2)
            fieldmap(x,y,z)=original(x,y,z)-alpha*multipliers(x,y,z)+alpha*sum;  
	}    
}


void irtkLaplacianSmoothingMissingData::Smooth(irtkRealImage& im, irtkRealImage m, irtkRealImage imask)
{
  double prev_lap, lap;
  irtkRealImage image(im);
  image.Write("original.nii.gz");
  EnlargeImage(image);
  image.Write("enlarged.nii.gz");
  
  irtkRealImage mask(m);
  EnlargeImage(mask);
  CalculateBoundary(mask);
  mask.Write("mask-with-boundary.nii.gz");
  
  irtkRealImage image_mask(imask);
  EnlargeImage(image_mask);
  image_mask.Write("image_mask.nii.gz");
  
  irtkRealImage weights(mask);
  weights=0;
  CalculateBoundaryWeights(weights,mask);
  weights.Write("weights.nii.gz");
  //mask.Write("enlarged-mask.nii.gz");
  
  irtkRealImage fieldmap(image);
  irtkRealImage laplacian(image);
  irtkRealImage llaplacian(image);
  irtkRealImage multipliers(image);
  laplacian=0;
  llaplacian=0;
  multipliers=0;
  

  double alpha = 0.5;
  char buffer[256];
  
  for (int iter = 1; iter<1000; iter++)
    {	
      prev_lap=lap;
      lap = LaplacianImage(fieldmap, mask,laplacian);
      laplacian.Write("laplacian.nii.gz");
      LaplacianBoundary(fieldmap, mask,laplacian);
      laplacian.Write("laplacian-with-boundary.nii.gz");
      multipliers = multipliers + laplacian*alpha;
      multipliers.Write("multipliers.nii.gz");
      
      UpdateFieldmap(fieldmap, image, multipliers, mask, image_mask, weights, alpha);
      //fieldmap = fieldmap*(1-alpha)+image*alpha -llaplacian*alpha;

      //sprintf(buffer,"fieldmap%i.nii.gz",iter);
      //fieldmap.Write(buffer);
      fieldmap.Write("fieldmap.nii.gz");
      
    
      cout<<"Iter "<<iter<<" alpha "<<alpha<<" prev_lap "<<prev_lap<<" lap "<<lap<<" reldiff "<<(prev_lap-lap)/prev_lap<<endl;
      
      //if(((prev_lap-lap)/prev_lap < 0.001)&&(iter>1))
	 //break;

    }
    fieldmap.Write("fieldmap.nii.gz");
    //reduce fieldmap to original size
    ReduceImage(fieldmap);
    fieldmap.Write("reduced.nii.gz");
    //mask out all the boundary voxels
    m.Write("m.nii.gz");
    
    /*
    bool at = m.GetImageAttributes()._xorigin==fieldmap.GetImageAttributes()._xorigin;
    cout<<"Attributes: "<<at<<endl;
    
    at = m.GetImageAttributes()._yorigin==fieldmap.GetImageAttributes()._yorigin;
    cout<<"Attributes: "<<at<<endl;
    
    at = m.GetImageAttributes()._zorigin==fieldmap.GetImageAttributes()._zorigin;
    cout<<"Attributes: "<<at<<endl;
    
    at = m.GetImageAttributes()._torigin==fieldmap.GetImageAttributes()._torigin;
    cout<<"Attributes: "<<at<<endl;
    
    //m.GetImageAttributes().Print();
    //fieldmap.GetImageAttributes().Print();
    cout.precision(10);
    cout<<m.GetImageAttributes()._yorigin<<endl;
    cout<<fieldmap.GetImageAttributes()._yorigin<<endl;
    cout<<im.GetImageAttributes()._yorigin<<endl;
    */
    /*(_x  == attr._x)  && (_y  == attr._y)  && (_z  == attr._z) && (_t  == attr._t) &&
          (_dx == attr._dx) && (_dy == attr._dy) && (_dz == attr._dz) && (_dt == attr._dt) &&
          (_xaxis[0] == attr._xaxis[0]) && (_xaxis[1] == attr._xaxis[1]) && (_xaxis[2] == attr._xaxis[2]) &&
          (_yaxis[0] == attr._yaxis[0]) && (_yaxis[1] == attr._yaxis[1]) && (_yaxis[2] == attr._yaxis[2]) &&
          (_zaxis[0] == attr._zaxis[0]) && (_zaxis[1] == attr._zaxis[1]) && (_zaxis[2] == attr._zaxis[2]) &&
          (_xorigin == attr._xorigin) && (_yorigin == attr._yorigin) && (_zorigin == attr._zorigin) && 
          (_torigin == attr._torigin))
          */
    im = fieldmap*m;
    
}

void irtkLaplacianSmoothingMissingData::UpsampleFieldmap(irtkRealImage& target, irtkRealImage mask, irtkRealImage newmask)
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

  f.Write("orig-f-pad.nii.gz");
  
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
		_fieldmap(i,j,k,t) = 0;
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

irtkRealImage irtkLaplacianSmoothingMissingData::Run()
{
     InitializeFactors();
    
    double step;
    irtkImageAttributes attr = _image.GetImageAttributes();
    irtkRealImage image(_image), mask(_mask),image_mask(_image_mask);
    //step = attr._dz*8;
    step=6;
    Blur(image,step/2,_mask,-1000);
    image.Write("blurred-image-8.nii.gz");
    Resample(image,step,-1000);
    SetPaddingToZero(image,-1000);
    Resample(mask,step,0);
    Resample(image_mask,step,0);
    image.Write("image-8.nii.gz");
    mask.Write("mask-8.nii.gz");
    image_mask.Write("image-mask-8.nii.gz");
    
    Smooth(image,mask,image_mask);
    image.Write("fieldmap-8.nii.gz");
    _fieldmap=image;
    _m=mask;
    //exit(1);
    
    /* TEST starts*/
    image = _image;
    mask=_mask;
    image_mask = _image_mask;
    //step = attr._dz*4;
    step=3;
    Blur(image,step/2,_mask, -1000);
    image.Write("blurred-image4.nii.gz");
    Resample(image,step,-1000);
    SetPaddingToZero(image,-1000);
    Resample(mask,step,0);
    Resample(image_mask,step,0);
    image.Write("image4.nii.gz");
    mask.Write("mask4.nii.gz");
    image_mask.Write("image_mask4.nii.gz");

    
    //_fieldmap.Write("_fieldmap.nii.gz");
    //_m.Write("_m.nii.gz");
   // mask.Write("mask.nii.gz");
    UpsampleFieldmap(image,_m,mask);
    _fieldmap.Write("upsampled-fieldmap-8.nii.gz");
    image*=mask;
    image.Write("image-masked.nii.gz");
    image-=_fieldmap;
    image.Write("image-masked-res.nii.gz");
    
    Smooth(image,mask,image_mask);
    image.Write("fieldmap-4-res.nii.gz");
    _fieldmap += image;
    _fieldmap.Write("fieldmap-4.nii.gz");
    _m=mask;
    
    UpsampleFieldmap(_image,_m,_mask);
     //_fieldmap.Write("final-fieldmap.nii.gz");
    /*TEST ends */
    
    return _fieldmap;
    /* 
    if(_have_target)
    {
      irtkImageTransformation itr;
      irtkNearestNeighborInterpolateImageFunction interp;
      itr.PutInterpolator(&interp);
      itr.PutTargetPaddingValue(-1000);
      itr.PutSourcePaddingValue(-1000);
      irtkRigidTransformation id;
      itr.SetInput(&_mask,&id);
      mask = _target;
      itr.SetOutput(&mask);
      itr.Run();
      mask.Write("targetmask.nii.gz");
      UpsampleFieldmap(_target,_mask,mask);
      _fieldmap.Write("targetfieldmap.nii.gz");
      _target.Write("target.nii.gz");
      TransformImage(_target, _fieldmap, mask, 1);
      _target.Write("transformed.nii.gz");
       
    }
*/
}

irtkRealImage irtkLaplacianSmoothingMissingData::Run3levels()
{
     InitializeFactors();
    
    double step;
    irtkImageAttributes attr = _image.GetImageAttributes();
    
    //level 3
    
    irtkRealImage image(_image), mask(_mask),image_mask(_image_mask);
    //step = attr._dz*8;
    step=6;
    Blur(image,step/2,_image_mask,-1000);
    image.Write("blurred-image-8.nii.gz");
    Resample(image,step,-1000);
    SetPaddingToZero(image,-1000);
    Resample(mask,step,0);
    Resample(image_mask,step,0);
    image.Write("image-8.nii.gz");
    mask.Write("mask-8.nii.gz");
    image_mask.Write("image-mask-8.nii.gz");
    
    Smooth(image,mask,image_mask);
    image.Write("fieldmap-8.nii.gz");
    _fieldmap=image;
    _m=mask;

    //level 2
    
    image = _image;
    mask=_mask;
    image_mask = _image_mask;
    //step = attr._dz*4;
    step=3;
    Blur(image,step/2,_image_mask, -1000);
    image.Write("blurred-image4.nii.gz");
    Resample(image,step,-1000);
    SetPaddingToZero(image,-1000);
    Resample(mask,step,0);
    Resample(image_mask,step,0);
    image.Write("image4.nii.gz");
    mask.Write("mask4.nii.gz");
    image_mask.Write("image_mask4.nii.gz");

    
    UpsampleFieldmap(image,_m,mask);
    _fieldmap.Write("upsampled-fieldmap-8.nii.gz");
    image*=mask;
    image.Write("image-masked.nii.gz");
    image-=_fieldmap;
    image*=image_mask;
    image.Write("image-masked-res.nii.gz");
    
    Smooth(image,mask,image_mask);
    image.Write("fieldmap-4-res.nii.gz");
    _fieldmap += image;
    _fieldmap.Write("fieldmap-4.nii.gz");
    _m=mask;
    
    //level 3
    
    image = _image;
    mask=_mask;
    image_mask = _image_mask;
    //step = attr._dz*4;
    step=1.5;
    Blur(image,step/2,_image_mask, -1000);
    image.Write("blurred-image2.nii.gz");
    Resample(image,step,-1000);
    SetPaddingToZero(image,-1000);
    Resample(mask,step,0);
    Resample(image_mask,step,0);
    image.Write("image2.nii.gz");
    mask.Write("mask2.nii.gz");
    image_mask.Write("image_mask2.nii.gz");

    
    UpsampleFieldmap(image,_m,mask);
    _fieldmap.Write("upsampled-fieldmap-4.nii.gz");
    image*=mask;
    image.Write("image-masked.nii.gz");
    image-=_fieldmap;
    image*=image_mask;
    image.Write("image-masked-res.nii.gz");
    
    Smooth(image,mask,image_mask);
    image.Write("fieldmap-2-res.nii.gz");
    _fieldmap += image;
    _fieldmap.Write("fieldmap-2.nii.gz");
    _m=mask;
        
    
    //finish
    UpsampleFieldmap(_image,_m,_mask);

    return _fieldmap;
}

/*
irtkRealImage irtkLaplacianSmoothingMissingData::Run3levels()
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
    image.Write("image-8.nii.gz");
    mask.Write("mask-8.nii.gz");
    
    Smooth(image,mask);
    image.Write("fieldmap-8.nii.gz");
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
    image.Write("image4.nii.gz");
    mask.Write("mask4.nii.gz");

    
    //_fieldmap.Write("_fieldmap.nii.gz");
    //_m.Write("_m.nii.gz");
   // mask.Write("mask.nii.gz");
    UpsampleFieldmap(image,_m,mask);
    _fieldmap.Write("upsampled-fieldmap-8.nii.gz");
    
    image*=mask;
    image.Write("image-masked.nii.gz");
    image-=_fieldmap;
    image.Write("image-masked-res.nii.gz");
    
    Smooth(image,mask);
    image.Write("fieldmap-4-res.nii.gz");
    _fieldmap += image;
    _fieldmap.Write("fieldmap-4.nii.gz");
    _m=mask;
    
    
    image = _image;
    mask=_mask;
    step = attr._dz*2;
    step = 1.5;
    Blur(image,step/2,-1000);
    Resample(image,step,-1000);
    Resample(mask,step,0);
    image.Write("image2.nii.gz");
    mask.Write("mask2.nii.gz");

    
    //_fieldmap.Write("_fieldmap.nii.gz");
    //_m.Write("_m.nii.gz");
   // mask.Write("mask.nii.gz");
    UpsampleFieldmap(image,_m,mask);
    _fieldmap.Write("upsampled-fieldmap-4.nii.gz");
    image*=mask;
    image.Write("image-masked2.nii.gz");
    image-=_fieldmap;
    image.Write("image-masked-res2.nii.gz");
    
    Smooth(image,mask);
    image.Write("fieldmap-2-res.nii.gz");
    _fieldmap += image;
    _fieldmap.Write("fieldmap-2.nii.gz");
    _m=mask;
    
    
    
    UpsampleFieldmap(_image,_m,_mask);
     //_fieldmap.Write("final-fieldmap.nii.gz");
    return _fieldmap;
  
}
*/