#include <irtkImage.h>
#include <irtkLaplacianSmoothingMissingData.h>

irtkRealImage _target, _source, _fieldmap, _mask,_original_source, _image_mask;
bool ok, _swap=true;
int _padding = 0;
char buffer[256];
void usage()
{
  cerr<<"distortionbm [target] [source] [output fieldmap] <-mask mask>"<<endl;
  exit(1);
}

double CC3D(int ti, int tj, int tk, int si, int  sj, int sk, int ns, int tt=0)
{
  int i,j,k,n,count;
  double tmean,smean, tvar, svar,t,s;
  
  tmean=0;smean=0;count=0;
  for (i=-ns; i<=ns;i++)
    for (j=-ns; j<=ns;j++)
      for (k=-ns; k<=ns;k++)
      {
	t =_target(ti+i,tj+j,tk+k,tt);
	if (t>0)
	{
	  tmean+=t;
	  count++;
	  
	  s =_source(si+i,sj+j,sk+k,tt);
	  smean+=s;
	}
      }
  if(count>0)    
  {
    tmean/=count;
    smean/=count;
  }
  else
  {
    tmean=0;
    smean=0;
  }
  
  //cout<<"Means: "<<tmean<<", "<<smean<<endl;
  
  tvar=0;svar=0;
  for (i=-ns; i<=ns;i++)
    for (j=-ns; j<=ns;j++)
      for (k=-ns; k<=ns;k++)
      {
	t =_target(ti+i,tj+j,tk+k,tt);
	if (t>0)
	{
	  t=t-tmean;
	  tvar+=t*t;
	  
	  s =_source(si+i,sj+j,sk+k,tt)-smean;
	  svar+=s*s;
	}
      }
  tvar/=count;
  svar/=count;
  tvar=sqrt(tvar);
  svar = sqrt(svar);
  //cout<<"Variances: "<<tvar<<", "<<svar<<endl;
    
  double ncc=0;
  if((tvar>0)&&(svar>0))
  {
    for (i=-ns; i<=ns;i++)
      for (j=-ns; j<=ns;j++)
        for (k=-ns; k<=ns;k++)
        {
	  t =_target(ti+i,tj+j,tk+k,tt);
	  if (t>0)
	  {
	    t=(t-tmean)/tvar;  
	    s =_source(si+i,sj+j,sk+k,tt);
	    s=(s-smean)/svar;
	    ncc+=t*s;
	  }
        }
    ncc/=count;
  }
  else
    ncc=-2;
  //cout<<"NCC is "<<ncc<<endl;
  return ncc;
}

double CC(int ti, int tj, int tk, int si, int  sj, int sk, int ns, int tt=0)
{
  int i,j,k,n,count;
  double tmean,smean, tvar, svar,t,s;
  
  tmean=0;smean=0;count=0;
  for (i=-ns; i<=ns;i++)
    for (j=-ns; j<=ns;j++)
      //for (k=-ns; k<=ns;k++)
      {
	t =_target(ti+i,tj+j,tk,tt);
	if (t>0)
	{
	  tmean+=t;
	  count++;
	  
	  s =_source(si+i,sj+j,sk,tt);
	  smean+=s;
	}
      }
  if(count>0)    
  {
    tmean/=count;
    smean/=count;
  }
  else
  {
    tmean=0;
    smean=0;
  }
  
  //cout<<"Means: "<<tmean<<", "<<smean<<endl;
  
  tvar=0;svar=0;
  for (i=-ns; i<=ns;i++)
    for (j=-ns; j<=ns;j++)
      //for (k=-ns; k<=ns;k++)
      {
	t =_target(ti+i,tj+j,tk,tt);
	if (t>0)
	{
	  t=t-tmean;
	  tvar+=t*t;
	  
	  s =_source(si+i,sj+j,sk,tt)-smean;
	  svar+=s*s;
	}
      }
  tvar/=count;
  svar/=count;
  tvar=sqrt(tvar);
  svar = sqrt(svar);
  //cout<<"Variances: "<<tvar<<", "<<svar<<endl;
    
  double ncc=0;
  if((tvar>0)&&(svar>0))
  {
    for (i=-ns; i<=ns;i++)
      for (j=-ns; j<=ns;j++)
        //for (k=-ns; k<=ns;k++)
        {
	  t =_target(ti+i,tj+j,tk,tt);
	  if (t>0)
	  {
	    t=(t-tmean)/tvar;  
	    s =_source(si+i,sj+j,sk+k,tt);
	    s=(s-smean)/svar;
	    ncc+=t*s;
	  }
        }
    ncc/=count;
  }
  else
    ncc=-2;
  //cout<<"NCC is "<<ncc<<endl;
  return ncc;
}

int FindDisplacement(int i,int j,int k,int t,int _block_size,int ns, double& ncc)
{
  double tncc;
  int d;
  int displacement = 0;

  ncc = -2;
  for (d=-ns;d<=ns;d++)
  { 
    if((j+d>=0)&&(j+d < _source.GetY()))
    {
      tncc = CC(i,j,k,i,j+d,k,_block_size);
      if(tncc>ncc)
      {
	ncc=tncc;
	displacement=d;
      }
    }
  }
  
  return displacement;
}

void CreateMaskFromTarget()
{
  _mask=_target;
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  irtkRealPixel *pt = _target.GetPointerToVoxels();
  for(int i=0;i<_target.GetNumberOfVoxels();i++)
  {
    if(*pt>0)
      *pm=1;
    else
      *pm=0;
    pm++;
    pt++;
  }
}


void CorrectSource()
{
  //padding for fieldmap
  irtkRealPixel mn,mx;
  _fieldmap.GetMinMax(&mn,&mx);
  int padding = round(mn-2);
  cout<<"padding = "<<padding<<endl;
  irtkRealImage fieldmap(_fieldmap);
  irtkRealPixel *pf = fieldmap.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  for (int j=0; j<fieldmap.GetNumberOfVoxels();j++)
  {
    if(*pm!=1)
      *pf = padding;
    pm++;
    pf++;
  }

  fieldmap.Write("fieldmap-padded.nii.gz");

  //prepare fieldmap for interpolation with padding
  irtkLinearInterpolateImageFunction interpolatorFieldmap;
  interpolatorFieldmap.SetInput(&fieldmap);
  interpolatorFieldmap.Initialize();
 
  //undistort image using fieldmap
  double f,x,y,z;
  irtkLinearInterpolateImageFunction interpolator;
  //irtkSincInterpolateImageFunction interpolator;
  interpolator.SetInput(&_original_source);
  interpolator.Initialize();
  irtkRealImage _f(_source);
  _f=0;
        
  //Correcting image
  irtkImageAttributes attr = _source.GetImageAttributes();
  for (int t = 0; t < _source.GetT(); t++) {
    for (int k = 0; k < _source.GetZ(); k++) {
      for (int j = 0; j < _source.GetY(); j++) {
        for (int i = 0; i < _source.GetX(); i++) {
          //find fieldmap value f
          x = i;
          y = j;
          z = k;
          _source.ImageToWorld(x,y,z);
          _fieldmap.WorldToImage(x,y,z);
          if ((x > -0.5) && (x < _fieldmap.GetX()-0.5) && 
	      (y > -0.5) && (y < _fieldmap.GetY()-0.5) &&
              (z > -0.5) && (z < _fieldmap.GetZ()-0.5) )
	  {
	    f = interpolatorFieldmap.EvaluateWithPadding(x,y,z,t,padding);
	  }
	  else
	    f=padding;
	_f(i,j,k,t)=f;	
	//find displacement
	if (f>padding)  
	{
	  x = i;
	  y = j;
	  z = k;

	  if(_swap)
	    y+=f/attr._dy;
	  else
	    x+=f/attr._dx;
		 
	 if ((x > -0.5) && (x < _source.GetX()-0.5) && 
	     (y > -0.5) && (y < _source.GetY()-0.5) &&
             (z > -0.5) && (z < _source.GetZ()-0.5))
	 {
	    _source(i,j,k,t) = interpolator.Evaluate(x,y,z,t);
	 }
	 else
	  _source(i,j,k,t)=0;
	}
	else
	  _source(i,j,k)=0;
	}//i
      }//j
    }//k
  }//t
  _source.Write("corrected-source.nii.gz");
  _f.Write("_f.nii.gz");
}

void EstimateFieldmap(int _block_size, int _neighbourhood_size)
{
  _image_mask=_target;
  _image_mask=0;
  irtkRealImage f(_fieldmap);
  f=0;
  irtkRealImage _ncc(_fieldmap);
  irtkImageAttributes attr = _target.GetImageAttributes();
  double ncc;
  for (int t = 0; t < _target.GetT(); t++) {
    for (int k = 0; k < _target.GetZ(); k++) {
      for (int j = _block_size; j < _target.GetY()-_block_size; j++) {
        for (int i = _block_size; i < _target.GetX()-_block_size; i++) {
	  if(_target(i,j,k,t)>0)
	  {
	    //cout<<i<<" "<<j<<" "<<k<<endl;
	    double displacement = FindDisplacement(i,j,k,t,_block_size,_neighbourhood_size,ncc);
	    _ncc(i,j,k,t)=ncc;
	    f(i,j,k,t)=displacement*attr._dy;
	    if(ncc>0)
	      _image_mask(i,j,k)=1;
	  }
	}
      }
    }
  }
 
 _image_mask.Write("_image_mask.nii.gz");
  f.Write("displacement.nii.gz");
  
  CreateMaskFromTarget();
  _mask.Write("mask.nii.gz");
  irtkLaplacianSmoothingMissingData smoothing;
  smoothing.SetInput(f);
  smoothing.SetMask(_mask);
  smoothing.SetInputMask(_image_mask);
  f = smoothing.Run3levels();
  _fieldmap+=f;
  _ncc.Write("ncc.nii.gz");
 
}

void EstimateFieldmapRS(int _block_size, int _neighbourhood_size)
{
  _image_mask=_target;
  _image_mask=0;
  irtkRealImage f(_fieldmap), displacement(_fieldmap);
  f=0;
  irtkRealImage _ncc(_fieldmap);
  irtkImageAttributes attr = _target.GetImageAttributes();
  double ncc;
  for (int t = 0; t < _target.GetT(); t++) {
    for (int k = 0; k < _target.GetZ(); k++) {
      for (int j = _block_size; j < _target.GetY()-_block_size; j++) {
        for (int i = _block_size; i < _target.GetX()-_block_size; i++) {
	  if(_target(i,j,k,t)>0)
	  {
	    //cout<<i<<" "<<j<<" "<<k<<endl;
	    double displacement = FindDisplacement(i,j,k,t,_block_size,_neighbourhood_size,ncc);
	    _ncc(i,j,k,t)=ncc;
	    f(i,j,k,t)=displacement*attr._dy;
	    if(ncc>0)
	      _image_mask(i,j,k)=1;
	  }
	}
      }
    }
  }
 
 _image_mask.Write("_image_mask.nii.gz");
  f.Write("displacement.nii.gz");
  displacement = f;
    
  CreateMaskFromTarget();
  _mask.Write("mask.nii.gz");
  irtkLaplacianSmoothingMissingData smoothing;
  smoothing.SetInput(displacement);
  smoothing.SetMask(_mask);
  smoothing.SetInputMask(_image_mask);
  f = smoothing.Run3levels();
  f.Write("f-iter1.nii.gz");
  
  char buffer[256];
  for (int iter=2; iter<10;iter++)
  {
    //calculate error
    irtkRealImage error(_fieldmap);
    error=displacement-f;
    sprintf(buffer,"error%i.nii.gz",iter);
    error.Write(buffer);
    irtkRealImage image_mask(_image_mask);
  
    irtkRealPixel *pe = error.GetPointerToVoxels();
    irtkRealPixel *pm = image_mask.GetPointerToVoxels();
    double dy = attr._dy;
    for(int i=0;i<error.GetNumberOfVoxels();i++)
    {
      if((*pe>dy)||(*pe<-dy))
        *pm=0;
      pm++;
      pe++;
    }
    sprintf(buffer,"image_mask_after_RS%i.nii.gz",iter);
    image_mask.Write(buffer);
  
    smoothing.SetInput(displacement);
    smoothing.SetMask(_mask);
    smoothing.SetInputMask(image_mask);
    f = smoothing.Run3levels();
    sprintf(buffer,"f-iter%i.nii.gz",iter);
    f.Write(buffer);
  }

  /////////////////////////////////// 
  _fieldmap+=f;
  _ncc.Write("ncc.nii.gz");

}



int main(int argc, char **argv)
{
    if (argc < 3)
    usage();
    
    _target.Read(argv[1]);
    argc--;
    argv++;

    _source.Read(argv[1]);
    argc--;
    argv++;
    _original_source=_source;

    char *output_name = argv[1];
    argc--;
    argv++;
    
    _mask=_target;
    _mask=1;
    // Parse options.
    while (argc > 1){
      ok = false;
    
      if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      _mask.Read(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
  
  int _block_size = 3;
  int _neighbourhood_size = 7;
  irtkImageAttributes attr = _target.GetImageAttributes();
  if(_target.GetT()>1)
  {
    cerr<<"Time dimesion is "<<_target.GetT()<<", not implemented yet."<<endl;
    exit(1);
  }
  
  
  
  _fieldmap=_target;
  _fieldmap=0;
  
  for (int iter=1; iter<=1;iter++)
  {
    EstimateFieldmapRS(_block_size, _neighbourhood_size);
    CorrectSource();
    sprintf(buffer,"f%i.nii.gz",iter);
    _fieldmap.Write(buffer);
    sprintf(buffer,"c%i.nii.gz",iter);
    _source.Write(buffer);
  }

  _fieldmap.Write(output_name);
  
 
}