/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: reconstructionb0.cc 998 2013-10-15 15:24:16Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-10-15 16:24:16 +0100 (Tue, 15 Oct 2013) $
  Version   : $Revision: 998 $
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
  cerr << "Usage: invertfieldmap [template image] [inputfieldmap] [outputfieldmap] <options>\n" << endl;
  cerr << "Options:\n" << endl;
  cerr << "<-x>        phase encoding direction is x [default: y]" << endl;
  cerr << "<-y>        phase encoding direction is y [default: y]" << endl;
  cerr << "<-minus>    change sign of phase encoding direction [default: plus]" << endl;
  cerr << "<-sinc>     Use sinc interpolation [default: linear]" << endl;
  cerr << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  
  //utility variables
  int ok;
  
  //declare variables for input
  /// Name for output volume
  char * output_name = NULL;
  //phase encode is y
  bool swap = true;
  //sign of phase encoding direction
  bool minus = false;
  //use sinc interpolation
  bool sinc = false;
  
  int i,j,k,t;
  double x,y,z,a;


  if (argc < 3)
    usage();

  //read template
  irtkRealImage image(argv[1]);
  cout<<"Image ... "<<argv[1]<<endl;
  cout.flush();
  argc--;
  argv++;

  //read template
  irtkRealImage fieldmap(argv[1]);
  cout<<"Image ... "<<argv[1]<<endl;
  cout.flush();
  argc--;
  argv++;
  
  //read output name
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"output name ... "<<output_name<<endl;
  cout.flush();

  

  // Parse options.
  while (argc > 1){
    ok = false;
    
    if ((ok == false) && (strcmp(argv[1], "-x") == 0)){
      argc--;
      argv++;
      swap=false;
      cout<< "Phase endoding direction is x."<<endl;
      cout.flush();
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-y") == 0)){
      argc--;
      argv++;
      swap=true;
      cout<< "Phase endoding direction is y."<<endl;
      cout.flush();
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
    
    if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)){
      argc--;
      argv++;
      sinc=true;
      cout<< "Use sinc."<<endl;
      cout.flush();
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  
  //resample fieldmap 
  irtkImageAttributes attr = image.GetImageAttributes();
  attr._t=1;
  irtkRealImage resfieldmap(attr);
  resfieldmap=0;
  irtkImageFunction *interpolatorLin = new irtkLinearInterpolateImageFunction;
  interpolatorLin->SetInput(&fieldmap);
  interpolatorLin->Initialize();
  
  //Step 1: for each line calculate how image gets distorted
  std::vector<double> distorted;
  
  if (swap) //phase encoding is y
  {
    for (k = 0; k < resfieldmap.GetZ(); k++) {
      for (i = 0; i < resfieldmap.GetX(); i++) {
	
	//find distorted grid
	distorted.clear();
	for (j = 0; j < resfieldmap.GetY(); j++) {
	   x = i;
           y = j;
           z = k;
	   resfieldmap.ImageToWorld(x, y, z);
           fieldmap.WorldToImage(x,y,z);
	
	  if ((x > -0.5) && (x < fieldmap.GetX()-0.5) && 
	      (y > -0.5) && (y < fieldmap.GetY()-0.5) &&
              (z > -0.5) && (z < fieldmap.GetZ()-0.5))
	    {
	      //location j after distortion in image coordinates (pixels)
	      if (minus)
	        distorted.push_back(j-interpolatorLin->Evaluate(x,y,z)/attr._dy);
	      else
	        distorted.push_back(j+interpolatorLin->Evaluate(x,y,z)/attr._dy);
	    }
        }
        
        //index - between which two distorted voxels current location falls
        t=-1;
	
	//invert fieldmap
        for (j = 0; j < resfieldmap.GetY(); j++) {
	   x = i;
           y = j;
           z = k;   
	   //check t is still inside ROI
	   if((t+1)<resfieldmap.GetY())
	   {  
	     //keep increasing t until you find that next one is larger than j
	     while (j>=distorted[t+1])
	     {
	        t++;
		//We arrived to the last t
	        if((t+1)==resfieldmap.GetY())
		  break;
	     }
	   }

	   
	     if((t+1)<resfieldmap.GetY())
	     {
	       a=(j-distorted[t])/(distorted[t+1]-distorted[t]);
	       y=(1-a)*t+a*(t+1);
	       
	       
	       
	     //i,j,k is location on the grid of inverted fieldmap
	     //now x,y,z represent location in undistorted space 
	     //but still in image coordinate space of inverted fieldmap
	     //convert to world
	     resfieldmap.ImageToWorld(x,y,z);
	     //convert to fieldmap
	     fieldmap.WorldToImage(x,y,z);
	     //get value of inverted fieldmap using linear interpolation
	     resfieldmap(i,j,k)=-interpolatorLin->Evaluate(x,y,z);
	  }
	
	}
      }
    }
  }
  else //phase encoding is x
  {
    for (k = 0; k < resfieldmap.GetZ(); k++) {
      for (j = 0; j < resfieldmap.GetY(); j++) {
	
	//find distorted grid
	distorted.clear();
	for (i = 0; i < resfieldmap.GetX(); i++) {
	   x = i;
           y = j;
           z = k;
	   resfieldmap.ImageToWorld(x, y, z);
           fieldmap.WorldToImage(x,y,z);
	
	  if ((x > -0.5) && (x < fieldmap.GetX()-0.5) && 
	      (y > -0.5) && (y < fieldmap.GetY()-0.5) &&
              (z > -0.5) && (z < fieldmap.GetZ()-0.5))
	    {
	      //location i after distortion in image coordinates (pixels)
	      if (minus)
	        distorted.push_back(i-interpolatorLin->Evaluate(x,y,z)/attr._dx);
	      else
	        distorted.push_back(i+interpolatorLin->Evaluate(x,y,z)/attr._dx);
	    }
        }
        
        //index - between which two distorted voxels current location falls
        t=-1;
	
	//invert fieldmap
        for (i = 0; i < resfieldmap.GetX(); i++) {
	   x = i;
           y = j;
           z = k;   
	   //check t is still inside ROI
	   if((t+1)<resfieldmap.GetX())
	   {  
	     //keep increasing t until you find that next one is larger than i
	     while (i>=distorted[t+1])
	     {
	        t++;
		//We arrived to the last t
	        if((t+1)==resfieldmap.GetX())
		  break;
	     }
	   }

	   
	     if((t+1)<resfieldmap.GetX())
	     {
	       a=(i-distorted[t])/(distorted[t+1]-distorted[t]);
	       x=(1-a)*t+a*(t+1);
	       
	       
	       
	     //i,j,k is location on the grid of inverted fieldmap
	     //now x,y,z represent location in undistorted space 
	     //but still in image coordinate space of inverted fieldmap
	     //convert to world
	     resfieldmap.ImageToWorld(x,y,z);
	     //convert to fieldmap
	     fieldmap.WorldToImage(x,y,z);
	     //get value of inverted fieldmap using linear interpolation
	     resfieldmap(i,j,k)=-interpolatorLin->Evaluate(x,y,z);
	  }
	
	}
      }
    }
  }
  
  
  
  
  
  resfieldmap.Write(output_name);
  
  //The end of main()

}  
