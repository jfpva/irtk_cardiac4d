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
#include <irtkLaplacianSmoothing.h>
#include <vector>
#include <string>
using namespace std;

//Application to perform reconstruction of volumetric MRI from thick slices.

void usage()
{
  cerr << "Usage: estimate-fieldmap [fieldmap] [distorted_stacks] [simulated_stacks] <options>\n" << endl;
  cerr << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  
  //utility variables
  int i, iter, ok;
  char buffer[256];
  irtkRealImage stack; 
  
  //declare variables for input
  /// Name for output volume
  char * output_name = NULL;
  /// Name for output volume
  char * output_fieldmap_name = NULL;
  /// images
  irtkRealImage original, simulated, corrected;
  ///number of stacks
  int nStacks;

    
  // Default values.
  irtkRealImage *mask=NULL;
  bool debug = false;
  double fieldMapSpacing = 5;
  double penalty = 0;
  double step = 2.5;
  int levels = 4;

  //Create reconstruction object
  irtkReconstructionb0 reconstruction;

  
  irtkRealImage average;

  string log_id;
  
  //distortion correction
  vector<int> groups, stack_group;
  vector<bool> swap;
  groups.push_back(1);
  swap.push_back(true);

  
  //if not enough arguments print help
  if (argc < 4)
    usage();
  
  //read output name
  output_fieldmap_name = argv[1];
  argc--;
  argv++;
  cout<<"Estimated fieldmap name ... "<<output_fieldmap_name<<endl;
  cout.flush();


  // Read stacks 
    original.Read(argv[1]);
    cout<<"Reading original stacks ... "<<argv[1]<<endl;
    cout.flush();
    argc--;
    argv++;
  
  // Read simulated stacks 
    simulated.Read(argv[1]);
    cout<<"Reading simulated stacks ... "<<argv[1]<<endl;
    cout.flush();
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
    
    //Variance of Gaussian kernel to smooth the bias field.
    if ((ok == false) && (strcmp(argv[1], "-fieldMapSpacing") == 0)){
      argc--;
      argv++;
      fieldMapSpacing=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Number of resolution levels for B-spline registration
    if ((ok == false) && (strcmp(argv[1], "-levels") == 0)){
      argc--;
      argv++;
      levels=atoi(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Variance of Gaussian kernel to smooth the bias field.
    if ((ok == false) && (strcmp(argv[1], "-smoothnessPenalty") == 0)){
      argc--;
      argv++;
      penalty=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    

    //Debug mode
    if ((ok == false) && (strcmp(argv[1], "-debug") == 0)){
      argc--;
      argv++;
      debug=true;
      ok = true;
    }
    
    //Read transformations from this folder
    if ((ok == false) && (strcmp(argv[1], "-log_prefix") == 0)){
      argc--;
      argv++;
      log_id=argv[1];
      ok = true;
      argc--;
      argv++;
    }


    //Read phase encoding direction
    if ((ok == false) && (strcmp(argv[1], "-phase") == 0)){
      argc--;
      argv++;
      //if(groups.size()==0)
      //{
//	cout<<"Please give groups before phase encoding direction (one per group)"<<endl;
//	exit(1);
  //    }
      cout<< "Phase encoding is ";
      for (i=0;i<1;i++)
      {
        if (strcmp(argv[1], "x") == 0)
	{
          swap[0]=false;
	  cout<<"x"<<" ";
          argc--;
          argv++;
	}
	else 
	{
	  if (strcmp(argv[1], "y") == 0)
	  {
            swap[0]=true;
	    cout<<"y"<<" ";
            argc--;
            argv++;
	  }
	  else
	  {
	    cerr<<"Unexpected phase encoding axis."<<endl;
	    exit(1);
	  }
	}
       }
       cout<<"."<<endl;
       cout.flush();
      ok = true;
    }



    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
   
  //Output volume
  irtkRealImage reconstructed;

  //Set up for distortion
  reconstruction.SetGroups(stack_group, groups, swap);
  //Set debug mode
  if (debug) reconstruction.DebugOn();
  else reconstruction.DebugOff();
  
  //Set B-spline control point spacing for field map
  reconstruction.SetFieldMapSpacing(fieldMapSpacing);
  
  //Set B-spline control point spacing for field map
  reconstruction.SetSmoothnessPenalty(penalty);
  
  //create template
  reconstruction.CreateTemplate(original,0);
  reconstructed = reconstruction.GetReconstructed();
  reconstructed.Write("template.nii.gz");
  //Set mask to reconstruction object. 
  reconstruction.SetMask(mask,0);   

  //to redirect output from screen to text files
  
  //to remember cout and cerr buffer
  streambuf* strm_buffer = cout.rdbuf();
  streambuf* strm_buffer_e = cerr.rdbuf();
  //files for registration output
  string name;
  name = log_id+"log-registration.txt";
  ofstream file(name.c_str());
  name = log_id+"log-registration-error.txt";
  ofstream file_e(name.c_str());
  //files for reconstruction output
  name = log_id+"log-reconstruction.txt";
  ofstream file2(name.c_str());
  name = log_id+"log-evaluation.txt";
  ofstream fileEv(name.c_str());
  //files for distortion output
  name = log_id+"log-distortion.txt";
  ofstream filed(name.c_str());
  name = log_id+"log-distortion-error.txt";
  ofstream filed_e(name.c_str());
  
  //set precision
  cout<<setprecision(3);
  cerr<<setprecision(3);
  
  irtkRealImage fieldmap = reconstruction.GetReconstructed();
  irtkMultiLevelFreeFormTransformation _fieldMap;
  cout<<"entering fieldmap distortion"<<endl;
  cout.flush();
  //cout<<levels<<endl;
  if(levels==3)
    step = 2.5;
  if(levels==4)
    step = 1.25;
  if(levels==5)
    step = 0.625;
  reconstruction.SetSmoothnessPenalty(penalty);
  reconstruction.FieldMapDistortion(original,simulated,_fieldMap,swap[0],step,0, levels);
  
  cout<<"done with fieldmap distortion"<<endl;
  cout.flush();

  _fieldMap.irtkTransformation::Write("distortion_transformation.dof.gz");
  
  
  
  
  double x,y,z,xx,yy,zz;
  irtkImageAttributes attr = original.GetImageAttributes();
  
  for(int k=0; k<fieldmap.GetZ();k++)
  for(int j=0; j<fieldmap.GetY();j++)
  for(int i=0; i<fieldmap.GetX();i++)
  {
    //transforming to the coordinates of the stack group
    //to get displacement in PE direction
    
    //origin
    xx=i;yy=j;zz=k;
    fieldmap.ImageToWorld(xx,yy,zz);
    original.WorldToImage(xx,yy,zz);

    //transformed
    x=i;y=j;z=k;
    fieldmap.ImageToWorld(x,y,z);
    _fieldMap.Transform(x,y,z);
    original.WorldToImage(x,y,z);
    
    if (swap[0])
      fieldmap(i,j,k)=(y-yy)*attr._dy;
    else
      fieldmap(i,j,k)=(x-xx)*attr._dx;
   }
   
   fieldmap.Write("fieldmap-bspline.nii.gz");
   //exit(1);
   irtkLaplacianSmoothing smoothing;
   smoothing.SetInput(fieldmap);
   smoothing.SetMask(*mask);
   irtkRealImage fieldmap_smooth = smoothing.RunGD();
   fieldmap_smooth.Write("fieldmap-smooth.nii.gz");


   fieldmap_smooth.Write(output_fieldmap_name);

  //prepare larger mask for fieldmap
  //reconstruction.CreateLargerMask(*mask);

/*
    //Print iteration number on the screen
    cout.rdbuf (strm_buffer);
    cout<<"Iteration "<<iter<<". "<<endl;
    cout.flush();
    
    /*change 1: do not uncorrect stack before distortion estimation
    corrected_stacks.clear();
    for(i=0;i<stacks.size();i++)
      corrected_stacks.push_back(stacks[i]);
    *//*
    if(corrected_stacks.size()==0)
    for(i=0;i<stacks.size();i++)
      corrected_stacks.push_back(stacks[i]);
    
    sprintf(buffer,"test-%i.nii.gz",iter);
    corrected_stacks[0].Write(buffer);
    
    //distortion
      //redirect output to files
      cerr.rdbuf(filed_e.rdbuf());
      cout.rdbuf (filed.rdbuf());
      

      reconstruction.Shim(corrected_stacks,iter);
        reconstruction.FieldMap(corrected_stacks,step,iter);
	step/=2;
	reconstruction.SaveDistortionTransformations();
        //change 2: Fieldmap is calculated even for affine only
	reconstruction.SmoothFieldmap(iter);

	corrected_stacks.clear();
        for(i=0;i<stacks.size();i++)
          corrected_stacks.push_back(stacks[i]);
        reconstruction.CorrectStacksSmoothFieldmap(corrected_stacks);
  //save final result
  reconstruction.SaveDistortionTransformations();
*/  
  //The end of main()
}  
