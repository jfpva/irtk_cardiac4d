/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: reconstructionb0.cc 1000 2013-10-18 16:49:08Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-10-18 17:49:08 +0100 (Fri, 18 Oct 2013) $
  Version   : $Revision: 1000 $
  Changes   : $Author: mm3 $

=========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkReconstruction.h>
#include <irtkReconstructionb0.h>
#include <vector>
#include <string>
using namespace std;

//Application to perform reconstruction of volumetric MRI from thick slices.

void usage()
{
  cerr << "Usage: ssd [image1] [image2] <-mask> <-file> <-method> <-padding> <-minus>" << endl;
  cerr << "<-param param><-first_param param>\n" << endl;
  cerr << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  
  int ok;
  char buffer[256];
  irtkRealImage image1, image2, m;  
  irtkRealImage *mask=NULL;
  char *file_name=NULL;
  char *method=NULL;
  double padding=0;
  bool minus = false;
  double param = 0, first_param=0;
    
  //if not enough arguments print help
  if (argc < 3)
    usage();
  
  image1.Read(argv[1]);  
  argc--;
  argv++;
  
  image2.Read(argv[1]);  
  argc--;
  argv++;
  
  
  // Parse options.
  while (argc > 1){
    ok = false;
    

    //Read binary mask for final volume
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask= new irtkRealImage(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-file") == 0)){
      argc--;
      argv++;

      file_name = argv[1];

      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-method") == 0)){
      argc--;
      argv++;

      method = argv[1];

      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;

      padding = atof(argv[1]);

      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-param") == 0)){
      argc--;
      argv++;

      param = atof(argv[1]);

      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-first_param") == 0)){
      argc--;
      argv++;

      first_param = atof(argv[1]);

      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-minus") == 0)){
      argc--;
      argv++;
      minus=true;
      cout<< "Sign of phase endoding direction is minus."<<endl;
      cout.flush();
      ok = true;
    }
    
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
  
  if(image1.GetNumberOfVoxels()!=image2.GetNumberOfVoxels())
  {
    cout<<"Give images on the same grid!"<<endl;
    exit(1);
  }
  
  if(mask!=NULL)
  {
    if(image1.GetNumberOfVoxels()!=mask->GetNumberOfVoxels())
    {
      cout<<"Give mask on the same grid as images!"<<endl;
      exit(1);
    }   
    else
      m=*mask;
  }
  else
  {
    m=image1;
    m=1;
  }
  
  irtkRealPixel *p1,*p2,*pm;
  p1 = image1.GetPointerToVoxels();
  p2 = image2.GetPointerToVoxels();
  pm = m.GetPointerToVoxels();
  
  double sum=0,num=0,max = 0,diff,mean,stdev;
  for(int i=0;i<image1.GetNumberOfVoxels();i++)
  {
    if((*pm==1)&&(*p1!=padding)&&(*p2!=padding))
    {
      if (minus)
        diff = *p1+*p2;
      else
        diff = *p1-*p2;
      if(fabs(diff)>max) max=fabs(diff);
      sum+=diff*diff;
      num++;
    }
    p1++;
    p2++;
    pm++;
  }
  mean = sum/num;
  //cout<<"sum = "<<sum<<endl;
  //cout<<"num = "<<num<<endl;
  cout<<"SSD = "<<mean<<endl;
  mean=sqrt(mean);
  cout<<"RMSE = "<<mean<<endl;
  cout<<"max diff = "<<max<<endl;
  
  p1 = image1.GetPointerToVoxels();
  p2 = image2.GetPointerToVoxels();
  pm = m.GetPointerToVoxels();
  
  sum=0,num=0;
  for(int i=0;i<image1.GetNumberOfVoxels();i++)
  {
    if((*pm==1)&&(*p1!=padding)&&(*p2!=padding))
    {
      if (minus)
        diff = fabs(*p1+*p2)-mean;
      else
        diff = fabs(*p1-*p2)-mean;
      sum+=diff*diff;
      num++;
    }
    p1++;
    p2++;
    pm++;
  }
  stdev = sqrt(sum/num);
  cout<<"stdev RMSE = "<<stdev<<endl;
  
 
  if (file_name !=NULL)
  {
    ofstream fileOut(file_name, ofstream::out | ofstream::app);

    if(!fileOut)
    {
      cerr << "Can't open file " << file_name << endl;
      exit(1);
    }

    fileOut.precision(3);  
    
    if(method != NULL)
    {
      if(param == first_param)
        fileOut<<endl<<method<<", ";
      fileOut<<mean<<", ";
    }
    else
    {
      if(param == first_param)
        fileOut<<endl;
      fileOut<<mean<<", ";
    }
  }
  //The end of main()
}  
