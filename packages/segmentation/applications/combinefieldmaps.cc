#include <irtkImage.h>
#include <irtkTransformation.h>


bool ok;
char buffer[256];

void usage()
{
  cerr<<"fit-fieldmap [fieldmap1] [stack1] [fieldmap2] [stack2]"<<endl;
  exit(1);
}

double CalculatePadding(irtkRealImage &f1, irtkRealImage &f2)
{
  irtkRealPixel min,max;
  double padding;
  f1.GetMinMax(&min,&max);
  padding = min;
  f2.GetMinMax(&min,&max);
  if(min<padding)
    padding = min;
  cout<<"Min = "<<padding<<endl;
  padding = floor(padding)-1;
  cout<<"padding = "<<padding<<endl;
  
  irtkRealPixel *p1 = f1.GetPointerToVoxels();
  irtkRealPixel *p2 = f2.GetPointerToVoxels();
  
  for(int i=0; i<f1.GetNumberOfVoxels();i++)
  {
    if(*p1==0)
      *p1=padding;
    if(*p2==0)
      *p2=padding;
    p1++;
    p2++;
  }
  
  f1.Write("f1-padding.nii.gz");
  f2.Write("f2-padding.nii.gz");
  return padding;
}

irtkRealImage TransformFieldmap(irtkRealImage fieldmap, irtkRealImage image, double c, double padding)
{
  irtkRealImage f(fieldmap);
  f=0;

  //extract PE direction
  double dx=0,dy=1,dz=0;
  double ox=0,oy=0,oz=0; //origin
  
  //cout<<"d in image coord: "<<dx<<" "<<dy<<" "<<dz<<endl;
  //cout<<"o in image coord: "<<ox<<" "<<oy<<" "<<oz<<endl;
  image.ImageToWorld(ox,oy,oz);
  image.ImageToWorld(dx,dy,dz);
  //cout<<"d in world coord: "<<dx<<" "<<dy<<" "<<dz<<endl;
  //cout<<"o in world coord: "<<ox<<" "<<oy<<" "<<oz<<endl;
  dx-=ox;
  dy-=oy;
  dz-=oz;
  //cout<<"PE in world coord: "<<dx<<" "<<dy<<" "<<dz<<endl;
  double norm = sqrt(dx*dx+dy*dy+dz*dz);
  //cout<< "norm = "<<norm<<endl;
  dx/=norm;
  dy/=norm;
  dz/=norm;
  //cout<<"norm PE in world coord: "<<dx<<" "<<dy<<" "<<dz<<endl;
  ox=0;oy=0;oz=0;
  f.WorldToImage(dx,dy,dz);
  f.WorldToImage(ox,oy,oz);
  dx-=ox;
  dy-=oy;
  dz-=oz;
  //cout<<"PE in fieldmap image coord: "<<dx<<" "<<dy<<" "<<dz<<endl;
  norm = sqrt(dx*dx+dy*dy+dz*dz);
  //cout<< "norm = "<<norm<<endl;
  
  irtkLinearInterpolateImageFunction interpolator;
  interpolator.SetInput(&fieldmap);
  interpolator.Initialize();
  
  int i,j,k;
  double x,y,z;
  for(i=0;i<f.GetX();i++)
    for(j=0;j<f.GetY();j++)
      for(k=0;k<f.GetZ();k++)
      {
	x=i-dx*c;
	y=j-dy*c;
	z=k-dz*c;
	if((x>=0)&&(x<f.GetX())&&(y>=0)&&(y<f.GetY())&&(z>=0)&&(z<f.GetZ()))
	{
	  double ff=interpolator.EvaluateWithPadding(x,y,z,0,padding);
	  //double ff=interpolator.Evaluate(x,y,z,0);
	  if(ff!=padding)
	    f(i,j,k)=ff-c;
	  else
	    f(i,j,k)=padding;
	}
	else
	  f(i,j,k)=padding;
	  //if(fieldmap(ii,jj,kk)!=0)
	    //f(i,j,k)=fieldmap(ii,jj,kk)-c;
      }
  return f;
}

double Similarity(irtkRealImage f1, irtkRealImage f2, double padding)
{
  irtkRealPixel *p1 = f1.GetPointerToVoxels();
  irtkRealPixel *p2 = f2.GetPointerToVoxels();
  
  if (f1.GetNumberOfVoxels()!=f2.GetNumberOfVoxels())
  {
    cout<<" Combinefieldmaps: similarity: Please give fieldmaps on the same grid!"<<endl;
    exit(1);
  }
  
  double ssd=0;
  int n=0;
  for(int i=0; i<f1.GetNumberOfVoxels();i++)
  {
    if((*p1!=padding)&&(*p2!=padding))
    {
      ssd += (*p1-*p2)*(*p1-*p2);
      n++;
    }
    p1++;
    p2++;
  }
  ssd/=n;
  return ssd;
}

irtkRealImage Subtract(irtkRealImage f1, irtkRealImage f2)
{
  irtkRealImage diff(f1);
  diff=0;
  
  irtkRealPixel *p1 = f1.GetPointerToVoxels();
  irtkRealPixel *p2 = f2.GetPointerToVoxels();
  irtkRealPixel *pd = diff.GetPointerToVoxels();
  
  if (f1.GetNumberOfVoxels()!=f2.GetNumberOfVoxels())
  {
    cout<<" Combinefieldmaps: similarity: Please give fieldmaps on the same grid!"<<endl;
    exit(1);
  }
  
  double ssd=0;
  int n=0;
  for(int i=0; i<f1.GetNumberOfVoxels();i++)
  {
    if((*p1!=0)&&(*p2!=0))
    {
      *pd = *p1-*p2;
    }
    p1++;
    p2++;
    pd++;
  }
  return diff;
}

int main(int argc, char **argv)
{
  
  irtkRealImage _image1, _fieldmap1;
  irtkRealImage _image2, _fieldmap2;
  
    if (argc < 3)
    usage();
    
    _fieldmap1.Read(argv[1]);
    argc--;
    argv++;

    _image1.Read(argv[1]);
    argc--;
    argv++;

    _fieldmap2.Read(argv[1]);
    argc--;
    argv++;

    _image2.Read(argv[1]);
    argc--;
    argv++;

/*
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

      if ((ok == false) && (strcmp(argv[1], "-image_mask") == 0)){
      argc--;
      argv++;
      _image_mask.Read(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

      if ((ok == false) && (strcmp(argv[1], "-target") == 0)){
      argc--;
      argv++;
      _target.Read(argv[1]);
      have_target=true;
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
*/  


  //double c=5;
  //irtkRealImage f = TransformFieldmap(_fieldmap2,_image2,c);
  //sprintf(buffer,"fieldmap2%f.nii.gz",c);    
  //f.Write(buffer);
  
  double _padding = CalculatePadding(_fieldmap1,_fieldmap2);
  
  cout<<"SSD of fieldmaps is "<<Similarity(_fieldmap1,_fieldmap2,_padding)<<endl;
  
  irtkRealImage sim(11,11,1);
  sim.PutPixelSize(10,10,10);
  irtkRealImage f1,f2;
  double scale = 1;
  
  double min_sim = 10;
  double minc1,minc2;
  
  //for (double c1=-5;c1<=5;c1++)
    //for (double c2=-5;c2<=5;c2++)
  for (int i =-5;i<=5;i++)
    for (int j=-5;j<=5;j++)
    {
      //double c1=c2;
      double c1=i*scale;
      double c2=j*scale;
      char buff[255];
      f1=TransformFieldmap(_fieldmap1,_image1,c1,_padding);
      sprintf(buff,"t1-%f.nii.gz",c1);
      f1.Write(buff);
      f2=TransformFieldmap(_fieldmap2,_image2,c2,_padding);
      sprintf(buff,"t2-%f.nii.gz",c2);
      f2.Write(buff);
      
      //sim(c1+5,c2+5,0)=sqrt(Similarity(f1,f2));
      double s = Similarity(f1,f2,_padding);
      sim(c1+5,c2+5,0)=s*10;
      if (s<min_sim)
      {
	min_sim=s;
	minc1=c1;
	minc2=c2;
      }
      cout<<c1<<" "<<c2<<" "<<s<<endl;
      //sim(c1+5,c2+5,0)=Similarity(f1,f2);
      
    }
      
  sim.Write("sim.nii.gz");
  
  irtkRealImage diff = Subtract(_fieldmap1,_fieldmap2);
  //diff.Write("diff.nii.gz");
  
  cout<<"Minimum at "<<minc1<<" "<<minc2<<endl;
  f1=TransformFieldmap(_fieldmap1,_image1,minc1,_padding);
  f1.Write("transformedFieldmap1.nii.gz");
  f2=TransformFieldmap(_fieldmap2,_image2,minc2,_padding);
  f2.Write("transformedFieldmap2.nii.gz");
  

  
  //_fieldmap.Write(output_name);

}
