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
  cerr << "Usage: reconstruction [reconstructed] [N] [stack_1] .. [stack_N] [dof_1] .. [dof_N] <options>\n" << endl;
  cerr << endl;

  cerr << "\t[reconstructed]         Name for the reconstructed volume. Nifti or Analyze format." << endl;
  cerr << "\t[N]                     Number of stacks." << endl;
  cerr << "\t[stack_1] .. [stack_N]  The input stacks. Nifti or Analyze format." << endl;
  cerr << "\t[dof_1]   .. [dof_N]    The transformations of the input stack to template" << endl;
  cerr << "\t                        in \'dof\' format used in IRTK." <<endl;
  cerr << "\t                        Only rough alignment with correct orienation and " << endl;
  cerr << "\t                        some overlap is needed." << endl;
  cerr << "\t                        Use \'id\' for an identity transformation for at least" << endl;
  cerr << "\t                        one stack. The first stack with \'id\' transformation" << endl;
  cerr << "\t                        will be resampled as template." << endl;
  cerr << "\t" << endl;
  cerr << "Options:" << endl;
  cerr << "\t-thickness [th_1] .. [th_N] Give slice thickness.[Default: twice voxel size in z direction]"<<endl;
  cerr << "\t-mask [mask]              Binary mask to define the region of interest. [Default: whole image]"<<endl;
  cerr << "\t-second_mask [mask]       Binary mask to define the region of interest for second group. [Default: whole image]"<<endl;
  cerr << "\t-fieldmap_mask [mask]     Binary mask to define the region of interest for fieldmap. [Default: whole image]"<<endl;
  cerr << "\t-packages [num_1] .. [num_N] Give number of packages used during acquisition for each stack."<<endl;
  cerr << "\t                          The stacks will be split into packages during registration iteration 1"<<endl;
  cerr << "\t                          and then into odd and even slices within each package during "<<endl;
  cerr << "\t                          registration iteration 2. The method will then continue with slice to"<<endl;
  cerr << "\t                          volume approach. [Default: slice to volume registration only]"<<endl;
  cerr << "\t-iterations [iter]        Number of registration-reconstruction iterations. [Default: 9]"<<endl;
  cerr << "\t-sigma [sigma]            Stdev for bias field. [Default: 12mm]"<<endl;
  cerr << "\t-resolution [res]         Isotropic resolution of the volume. [Default: 0.75mm]"<<endl;
  cerr << "\t-multires [levels]        Multiresolution smooting with given number of levels. [Default: 3]"<<endl;
  cerr << "\t-average [average]        Average intensity value for stacks [Default: 700]"<<endl;
  cerr << "\t-delta [delta]            Parameter to define what is an edge. [Default: 150]"<<endl;
  cerr << "\t-lambda [lambda]          Smoothing parameter. [Default: 0.02]"<<endl;
  cerr << "\t-lastIter [lambda]        Smoothing parameter for last iteration. [Default: 0.01]"<<endl;
  cerr << "\t-smooth_mask [sigma]      Smooth the mask to reduce artefacts of manual segmentation. [Default: 4mm]"<<endl;
  cerr << "\t-global_bias_correction   Correct the bias in reconstructed image against previous estimation."<<endl;
  cerr << "\t-low_intensity_cutoff     Lower intensity threshold for inclusion of voxels in global bias correction."<<endl;
  cerr << "\t-remove_black_background  Create mask from black background."<<endl;
  cerr << "\t-transformations [folder] Use existing slice-to-volume transformations to initialize the reconstruction."<<endl;
  cerr << "\t-force_exclude [number of slices] [ind1] ... [indN]  Force exclusion of slices with these indices."<<endl;
  cerr << "\t-no_intensity_matching    Switch off intensity matching."<<endl;
  cerr << "\t-group [group1] ... [groupN]  Give group numbers for distortion correction. [Default: all in one group]"<<endl;
  cerr << "\t-phase [axis1] ... [axisG]  For each group give phase encoding axis in image coordinates [x or y]"<<endl;
  cerr << "\t-fieldMapSpacing [spacing]  B-spline control point spacing for field map"<<endl;
  cerr << "\t-smoothnessPenalty [penalty]  smoothness penalty for field map"<<endl;
  cerr << "\t-log_prefix [prefix]      Prefix for the log file."<<endl;
  cerr << "\t-debug                    Debug mode - save intermediate results."<<endl;
  cerr << "\t" << endl;
  cerr << "\t" << endl;
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
  /// Slice stacks
  vector<irtkRealImage> stacks, corrected_stacks;
  /// Stack transformation
  vector<irtkRigidTransformation> stack_transformations;
  /// Stack thickness
  vector<double > thickness;
  ///number of stacks
  int nStacks;
  /// number of packages for each stack
  vector<int> packages;
  ///T2 template
  irtkRealImage t2;


    
  // Default values.
  int templateNumber=-1;
  irtkRealImage *mask = NULL;
  irtkRealImage *mask2 = NULL;
  irtkRealImage *fmask = NULL;
  irtkRealImage b0_mask,T2_mask;
  int iterations = 7;
  bool debug = false;
  double sigma=20;
  double resolution = 0.75;
  double lambda = 0.02;
  double delta = 150;
  int levels = 3;
  double lastIterLambda = 0.01;
  int rec_iterations;
  double averageValue = 700;
  double smooth_mask = 0;
  bool global_bias_correction = false;
  double low_intensity_cutoff = 0.01;
  //folder for slice-to-volume registrations, if given
  char * folder=NULL;
  //flag to remove black background, e.g. when neonatal motion correction is performed
  bool remove_black_background = false;
  //flag to swich the intensity matching on and off
  bool intensity_matching = true;
  bool alignT2 = false;
  double fieldMapSpacing = 5;
  double penalty = 0;
  bool combine_fieldmaps=true;
  irtkRealImage templ;
  irtkRealImage fieldmapMask1, fieldmapMask2;
  irtkRealImage correctedMask1, correctedMask2;

  //Create reconstruction object
  irtkReconstructionb0 reconstruction;

  
  irtkRealImage average;

  string log_id;
  
  //distortion correction
  vector<int> groups, stack_group;
  vector<bool> swap;

  
  //forced exclusion of slices
  int number_of_force_excluded_slices = 0;
  vector<int> force_excluded;
  //if not enough arguments print help
  if (argc < 2)
    usage();
  
  //read output name
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"Recontructed volume name ... "<<output_name<<endl;
  cout.flush();

  //read number of stacks
  nStacks = atoi(argv[1]);
  argc--;
  argv++;
  cout<<"Number 0f stacks ... "<<nStacks<<endl;
  cout.flush();

  // Read stacks 
  for (i=0;i<nStacks;i++)
  {
      //if ( i == 0 )
          //log_id = argv[1];
    stack.Read(argv[1]);
    cout<<"Reading stack ... "<<argv[1]<<endl;
    cout.flush();
    argc--;
    argv++;
    stacks.push_back(stack);
  }
  
  //Read transformation
  for (i=0;i<nStacks;i++)
  {
    irtkTransformation *transformation;
    cout<<"Reading transformation ... "<<argv[1]<<" ... ";
    cout.flush();
    if (strcmp(argv[1], "id") == 0)
    {
      transformation = new irtkRigidTransformation;
      if ( templateNumber < 0) templateNumber = i;
    }
    else
    {
      transformation = irtkTransformation::New(argv[1]);
    }
    cout<<" done."<<endl;
    cout.flush();

    argc--;
    argv++;
    irtkRigidTransformation *rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);
    stack_transformations.push_back(*rigidTransf);
    delete rigidTransf;
  }
  reconstruction.InvertStackTransformations(stack_transformations);

  // Parse options.
  while (argc > 1){
    ok = false;
    
    //Read slice thickness
    if ((ok == false) && (strcmp(argv[1], "-thickness") == 0)){
      argc--;
      argv++;
      cout<< "Slice thickness is ";
      for (i=0;i<nStacks;i++)
      {
        thickness.push_back(atof(argv[1]));
	cout<<thickness[i]<<" ";
        argc--;
        argv++;
       }
       cout<<"."<<endl;
       cout.flush();
      ok = true;
    }
    
    //Read number of packages for each stack
    if ((ok == false) && (strcmp(argv[1], "-packages") == 0)){
      argc--;
      argv++;
      cout<< "Package number is ";
      for (i=0;i<nStacks;i++)
      {
        packages.push_back(atoi(argv[1]));
	cout<<packages[i]<<" ";
        argc--;
        argv++;
       }
       cout<<"."<<endl;
       cout.flush();
      ok = true;
    }

    //Read binary mask for final volume
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask= new irtkRealImage(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Read binary mask for final volume
    if ((ok == false) && (strcmp(argv[1], "-second_mask") == 0)){
      argc--;
      argv++;
      mask2= new irtkRealImage(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Read binary mask for final volume
    if ((ok == false) && (strcmp(argv[1], "-fieldmap_mask") == 0)){
      argc--;
      argv++;
      fmask= new irtkRealImage(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Read number of registration-reconstruction iterations
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)){
      argc--;
      argv++;
      iterations=atoi(argv[1]);
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
    
    //Variance of Gaussian kernel to smooth the bias field.
    if ((ok == false) && (strcmp(argv[1], "-smoothnessPenalty") == 0)){
      argc--;
      argv++;
      penalty=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Variance of Gaussian kernel to smooth the bias field.
    if ((ok == false) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;
      argv++;
      sigma=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Smoothing parameter
    if ((ok == false) && (strcmp(argv[1], "-lambda") == 0)){
      argc--;
      argv++;
      lambda=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Smoothing parameter for last iteration
    if ((ok == false) && (strcmp(argv[1], "-lastIter") == 0)){
      argc--;
      argv++;
      lastIterLambda=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Parameter to define what is an edge
    if ((ok == false) && (strcmp(argv[1], "-delta") == 0)){
      argc--;
      argv++;
      delta=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Isotropic resolution for the reconstructed volume
    if ((ok == false) && (strcmp(argv[1], "-resolution") == 0)){
      argc--;
      argv++;
      resolution=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }

    //Number of resolution levels
    if ((ok == false) && (strcmp(argv[1], "-multires") == 0)){
      argc--;
      argv++;
      levels=atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    //Smooth mask to remove effects of manual segmentation
    if ((ok == false) && (strcmp(argv[1], "-smooth_mask") == 0)){
      argc--;
      argv++;
      smooth_mask=atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    //Switch off intensity matching
    if ((ok == false) && (strcmp(argv[1], "-no_intensity_matching") == 0)){
      argc--;
      argv++;
      intensity_matching=false;
      ok = true;
    }

    //Switch off intensity matching
    if ((ok == false) && (strcmp(argv[1], "-combine_fieldmaps") == 0)){
      argc--;
      argv++;
      combine_fieldmaps=true;
      ok = true;
    }

    //Perform bias correction of the reconstructed image agains the GW image in the same motion correction iteration
    if ((ok == false) && (strcmp(argv[1], "-global_bias_correction") == 0)){
      argc--;
      argv++;
      global_bias_correction=true;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-low_intensity_cutoff") == 0)){
      argc--;
      argv++;
      low_intensity_cutoff=atof(argv[1]);
      argc--;
      argv++;
      ok = true;
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

    //Read transformations from this folder
    if ((ok == false) && (strcmp(argv[1], "-transformations") == 0)){
      argc--;
      argv++;
      folder=argv[1];
      ok = true;
      argc--;
      argv++;
    }

    //Remove black background
    if ((ok == false) && (strcmp(argv[1], "-remove_black_background") == 0)){
      argc--;
      argv++;
      remove_black_background=true;
      ok = true;
    }

    //Force removal of certain slices
    if ((ok == false) && (strcmp(argv[1], "-force_exclude") == 0)){
      argc--;
      argv++;
      number_of_force_excluded_slices = atoi(argv[1]);
      argc--;
      argv++;

      cout<< number_of_force_excluded_slices<< " force excluded slices: ";
      for (i=0;i<number_of_force_excluded_slices;i++)
      {
        force_excluded.push_back(atoi(argv[1]));
	cout<<force_excluded[i]<<" ";
        argc--;
        argv++;
       }
       cout<<"."<<endl;
       cout.flush();

      ok = true;
    }
    //Read stack group number
    if ((ok == false) && (strcmp(argv[1], "-group") == 0)){
      argc--;
      argv++;
      cout<< "Stack group number is ";
      for (i=0;i<nStacks;i++)
      {
        stack_group.push_back(atoi(argv[1]));
	cout<<stack_group[i]<<" ";
        argc--;
        argv++;
       }
       cout<<"."<<endl;
       cout.flush();
       
       bool group_exists;
      for (int i=0;i<stack_group.size();i++)
      {
        //check whether we already have this group
        group_exists=false;
        for (int j=0;j<groups.size();j++)
	  if(stack_group[i]==groups[j])
	    group_exists=true;
        if (!group_exists)
	  groups.push_back(stack_group[i]);
      }
      //sort the group numbers
      sort(groups.begin(),groups.end());
      cout<<"We have "<<groups.size()<<" groups:";
      for (int j=0;j<groups.size();j++) cout<<" "<<groups[j];
      cout<<endl;
    
      ok = true;
    }

    //Read phase encoding direction
    if ((ok == false) && (strcmp(argv[1], "-phase") == 0)){
      argc--;
      argv++;
      if(groups.size()==0)
      {
	cout<<"Please give groups before phase encoding direction (one per group)"<<endl;
	exit(1);
      }
      cout<< "Phase encoding is ";
      for (i=0;i<groups.size();i++)
      {
        if (strcmp(argv[1], "x") == 0)
	{
          swap.push_back(false);
	  cout<<"x"<<" ";
          argc--;
          argv++;
	}
	else 
	{
	  if (strcmp(argv[1], "y") == 0)
	  {
            swap.push_back(true);
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
  
  //Initialise 2*slice thickness if not given by user
  if (thickness.size()==0)
  {
    cout<< "Slice thickness is ";
    for (i=0;i<nStacks;i++)
    {
      double dx,dy,dz;
      stacks[i].GetPixelSize(&dx,&dy,&dz);
      thickness.push_back(dz*2);
      cout<<thickness[i]<<" ";
    }
    cout<<"."<<endl;
    cout.flush();
  }

  
  //Output volume
  irtkRealImage reconstructed;

  //Set up for distortion
  reconstruction.SetGroups(stack_group, groups, swap);
  irtkRealImage alignedT2; //= reconstruction.AlignT2Template(t2);
  //alignedT2.Write("alignedT2.nii.gz");


  //Set debug mode
  if (debug) reconstruction.DebugOn();
  else reconstruction.DebugOff();
  
  //Set B-spline control point spacing for field map
  reconstruction.SetFieldMapSpacing(fieldMapSpacing);
  
  //Set B-spline control point spacing for field map
  reconstruction.SetSmoothnessPenalty(penalty);
  
  //Set force excluded slices
  reconstruction.SetForceExcludedSlices(force_excluded);

  //Set low intensity cutoff for bias estimation
  reconstruction.SetLowIntensityCutoff(low_intensity_cutoff)  ;

  
  // Check whether the template stack can be indentified
  if (templateNumber<0)
  {
    cerr<<"Please identify the template by assigning id transformation."<<endl;
    exit(1);
  }  
  //If no mask was given and flag "remove_black_background" is false, try to create mask from the template image in case it was padded
  if ((mask==NULL)&&(!remove_black_background))
  {
    mask = new irtkRealImage(stacks[templateNumber]);
    *mask = reconstruction.CreateMask(*mask);
  }
  
  
  
/* 
  //Before creating the template we will crop template stack according to the given mask
  if (mask !=NULL)
  {
    //first resample the mask to the space of the stack
    //for template stact the transformation is identity
    irtkRealImage m = *mask;
    reconstruction.TransformMask(stacks[templateNumber],m,stack_transformations[templateNumber]);
    //Crop template stack
    reconstruction.CropImage(stacks[templateNumber],m);
    if (debug)
    {
      m.Write("maskTemplate.nii.gz"); 
      stacks[templateNumber].Write("croppedTemplate.nii.gz");
    }
  }
 */ 
  
  //This is a temporary template
  //Create template volume with isotropic resolution 
  //if resolution==0 it will be determined from in-plane resolution of the image
  resolution = reconstruction.CreateTemplate(stacks[templateNumber],resolution);
  templ=stacks[templateNumber];

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
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////Algorithm starts here///////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  //Split stacks into two groups
  /// Slice stacks
  vector<irtkRealImage> stacks1, stacks2;
  /// Stack transformation
  vector<irtkRigidTransformation> stack_transformations1,stack_transformations2;
  ///Group index
  vector<int> group_index;
  for(uint ind=0;ind<stacks.size();ind++)
  {
    if(stack_group[ind]==groups[0])
    {
      stacks1.push_back(stacks[ind]);
      stack_transformations1.push_back(stack_transformations[ind]);
      group_index.push_back(stacks1.size()-1);
    }
    else
    {
      stacks2.push_back(stacks[ind]);
      stack_transformations2.push_back(stack_transformations[ind]);
      group_index.push_back(stacks2.size()-1);
    }
  }
    

  //perform volumetric registration of the stack groups
  //redirect output to files
  cerr.rdbuf(file_e.rdbuf());
  cout.rdbuf (file.rdbuf());
  
  ////group1
  
  //volumetric registration
  reconstruction.SetMask(mask,smooth_mask);
  reconstruction.irtkReconstruction::StackRegistrations(stacks1,stack_transformations1,0);
  average = reconstruction.CreateAverage(stacks1,stack_transformations1);
  if (debug)
    average.Write("average1.nii.gz");
  for (i=0; i<stacks1.size(); i++)
  {
    //transform the mask
    irtkRealImage m= *mask;
    m.Write("mask1.nii.gz");
    reconstruction.TransformMask(stacks1[i],m,stack_transformations1[i]);
    //Crop template stack
    reconstruction.CropImage(stacks1[i],m);
    if (debug)
    {
      sprintf(buffer,"mask1-%i.nii.gz",i);
      m.Write(buffer); 
      sprintf(buffer,"cropped1-%i.nii.gz",i);
      stacks1[i].Write(buffer);
    }
  }
 

  reconstruction.SetMask(mask2,smooth_mask);
  reconstruction.irtkReconstruction::StackRegistrations(stacks2,stack_transformations2,0);
  average = reconstruction.CreateAverage(stacks2,stack_transformations2);
  if (debug)
    average.Write("average2.nii.gz");
  for (i=0; i<stacks2.size(); i++)
  {
    //transform the mask
    irtkRealImage m= *mask2;
    m.Write("mask2.nii.gz");
    reconstruction.TransformMask(stacks2[i],m,stack_transformations2[i]);
    //Crop template stack
    reconstruction.CropImage(stacks2[i],m);
    if (debug)
    {
      sprintf(buffer,"mask2-%i.nii.gz",i);
      m.Write(buffer); 
      sprintf(buffer,"cropped2-%i.nii.gz",i);
      stacks2[i].Write(buffer);
    }
  }
  
  //put registered and cropped stacks into main vector
  int template2Number = -1;
  stacks.clear();
  stack_transformations.clear();
  for(i=0;i<nStacks;i++)
  {
    if(stack_group[i]==groups[0])
    {
      stacks.push_back(stacks1[group_index[i]]);
      stack_transformations.push_back(stack_transformations1[group_index[i]]);
    }
    else
    {
      if(template2Number==-1)
	template2Number=i;
      stacks.push_back(stacks2[group_index[i]]);
      stack_transformations.push_back(stack_transformations2[group_index[i]]);      
    }
  }
 
  cout<<"template2Number = "<<template2Number<<endl;
  cout.flush();
  
  //Repeat volumetric registrations with cropped stacks

  //volumetric registration
  reconstruction.SetMask(mask,smooth_mask);
  reconstruction.irtkReconstruction::StackRegistrations(stacks,stack_transformations,templateNumber);
  cout<<endl;
  cout.flush();
  //redirect output back to screen
  cout.rdbuf (strm_buffer);
  cerr.rdbuf (strm_buffer_e);
  //reconstruction.InvertStackTransformations(stack_transformations);
  
  
  //Create mask1 using both masks in the space of  the first group and the template
  irtkRealImage m = *mask2;
  m.Write("mask2.nii.gz");
  reconstruction.InvertStackTransformations(stack_transformations); //because transforming from stack to template
  reconstruction.TransformMask(templ,m,stack_transformations[template2Number]);
  reconstruction.InvertStackTransformations(stack_transformations);
  m.Write("mask2-transformed.nii.gz");
  vector<irtkRealImage> masks;
  masks.push_back(*mask);
  mask->Write("mask.nii.gz");
  masks.push_back(*mask2);
  vector<irtkRigidTransformation> mask_transformations;
  mask_transformations.push_back(stack_transformations[templateNumber]);
  mask_transformations.push_back(stack_transformations[template2Number]);
  m=reconstruction.CreateAverage(masks,mask_transformations);
  m.Write("average-mask.nii.gz");
  
  //Create mask2 using both masks in the space of the second group
  mask_transformations.clear();
  reconstruction.InvertStackTransformations(stack_transformations);
  mask_transformations.push_back(stack_transformations[template2Number]);
  mask_transformations.push_back(stack_transformations[templateNumber]);
  reconstruction.InvertStackTransformations(stack_transformations);
  irtkRealImage m2=reconstruction.CreateAverage(masks,mask_transformations);
  m2.Write("average-mask2.nii.gz");
  reconstruction.SetMask2(&m2);
  
  //Create combined mask to crop the template
  irtkRigidTransformation id;
  masks.clear();
  masks.push_back(m);
  masks.push_back(m2);
  mask_transformations.clear();
  mask_transformations.push_back(id);
  mask_transformations.push_back(id);
  m=reconstruction.CreateAverage(masks,mask_transformations);
  m.Write("average-mask-all.nii.gz");
  
  //create template using both masks
  reconstruction.TransformMask(templ,m,stack_transformations[templateNumber]);
  reconstruction.CropImage(templ,m);
  if (debug)
  {
    m.Write("maskTemplate.nii.gz"); 
    templ.Write("croppedTemplate.nii.gz");
  }
  resolution = reconstruction.CreateTemplate(templ,resolution);
  reconstruction.SetMask(&m,resolution,0.45);
  b0_mask=reconstruction.GetMask();
  b0_mask.Write("b0mask.nii.gz");

  
  
  //Rescale intensities of the stacks to have the same average
  if (intensity_matching)
    reconstruction.MatchStackIntensities(stacks,stack_transformations,averageValue);
  else
    reconstruction.MatchStackIntensities(stacks,stack_transformations,averageValue,true);
  average = reconstruction.CreateAverage(stacks,stack_transformations);
  if (debug)
    average.Write("average-all.nii.gz");
  reconstruction.SetReconstructed(average);
  fieldmapMask1=average;
  fieldmapMask1=1;
  fieldmapMask2=average;
  fieldmapMask2=1;
  correctedMask1=*mask;
  correctedMask2=*mask2;
  //exit(1);

  //Create slices and slice-dependent transformations
  reconstruction.CreateSlicesAndTransformations(stacks,stack_transformations,thickness);
  
  //Mask all the slices
  reconstruction.MaskSlices();
   
  //Set sigma for the bias field smoothing
  if (sigma>0)
    reconstruction.SetSigma(sigma);
  else
  {
    //cerr<<"Please set sigma larger than zero. Current value: "<<sigma<<endl;
    //exit(1);
    reconstruction.SetSigma(20);
  }
  
  //Set global bias correction flag
  if (global_bias_correction)
    reconstruction.GlobalBiasCorrectionOn();
  else 
    reconstruction.GlobalBiasCorrectionOff();
    
  //if given read slice-to-volume registrations
  //bool coeff_init=false;  
  if (folder!=NULL)
  {
    reconstruction.ReadTransformation(folder);
    //reconstruction.CoeffInit();
    //coeff_init = true;
  }
    
  //Initialise data structures for EM
  reconstruction.InitializeEM();
  reconstruction.InitializeEMValues();
  
  
  ////TEST
  //Select stacks
  int current_group = -1;
  double step = 2.5;
  //if(iterations>1)
    //iterations = 7;
  //registration iterations
  for (iter=0;iter<(iterations);iter++)
  {
    //Print iteration number on the screen
    cout.rdbuf (strm_buffer);
    cout<<"Iteration "<<iter<<". "<<endl;
    cout.flush();
    
    if(corrected_stacks.size()==0)
    for(i=0;i<stacks.size();i++)
      corrected_stacks.push_back(stacks[i]);
    
    
    //reconstruction.SetMask2(mask2);
    
    //distortion correction
      
    if(iter>0)
    {
      //redirect output to files
      cerr.rdbuf(filed_e.rdbuf());
      cout.rdbuf (filed.rdbuf());
      
    //if(current_group==0)
    //  reconstruction.SetMask(mask2,smooth_mask);
    //else
    //  reconstruction.SetMask(mask,smooth_mask);

      
      reconstruction.CoeffInitBSpline();

      //reconstruction.Shim(corrected_stacks,iter);
      if(iter>0)
      {
        //reconstruction.FieldMapGroup(corrected_stacks,1-current_group,step,iter);
        if((iter>2)&&(iter%8 == 1))
	  step/=2;
	cout<<"step="<<step<<endl;
	//reconstruction.SaveDistortionTransformations();
      }
        //change 2: Fieldmap is calculated even for affine only
        //TODO: will need mask 2 in next iteration
        
      cout<<"Current group = "<<current_group<<endl;
      irtkRealImage m;
      if(current_group == 0)
      {
        m= correctedMask1;
	//if(iter>1)
	  //reconstruction.CorrectMaskSmoothFieldmap(m,0);
        templ=reconstruction.GetMask();
        templ.Write("templ-for-mask.nii.gz");
        m.Write("mask-for-fieldmap-before-tra.nii.gz");
        reconstruction.TransformMask(templ,m,stack_transformations[template2Number]);
      //stacks[template2Number]
        m.Write("mask-for-fieldmap.nii.gz");
	fieldmapMask1=m;
      }
      else
      {
        m= correctedMask2;
	//reconstruction.CorrectMaskSmoothFieldmap(m,1);
	//exit(1);
        templ=reconstruction.GetMask();
        templ.Write("templ-for-mask-2.nii.gz");
        m.Write("mask-for-fieldmap-before-tra-2.nii.gz");
        reconstruction.InvertStackTransformations(stack_transformations);
        reconstruction.TransformMask(templ,m,stack_transformations[template2Number]);
        reconstruction.InvertStackTransformations(stack_transformations);
      //stacks[template2Number]
        m.Write("mask-for-fieldmap-2.nii.gz");
	fieldmapMask2=m;
	
      }
      
      //TEST to override varying fieldmap masks
      if(fmask != NULL)
      {
        fieldmapMask1=*fmask;
        fieldmapMask2=*fmask;
      }
      
      //if(current_group == 0)
        reconstruction.BSplineReconstructionGroup(groups[0]);
	reconstructed=reconstruction.GetReconstructed();
	reconstructed.Write("target0.nii.gz");
	reconstruction.FieldMapGroup(corrected_stacks,fieldmapMask1,1,step,iter);
	reconstruction.RememberDistortion();

      //else
        reconstruction.BSplineReconstructionGroup(groups[1]);
        reconstructed=reconstruction.GetReconstructed();
	reconstructed.Write("target1.nii.gz");
	reconstruction.FieldMapGroup(corrected_stacks,fieldmapMask2,0,step,iter);
	reconstruction.RememberDistortion();
	reconstruction.CombineDistortion();

	//reconstruction.SmoothFieldmapGroup(m,1-current_group,iter);
        //TEST to override varying fieldmap masks
	if(fmask!=NULL)
	  reconstruction.SmoothFieldmapGroup(*fmask,1-current_group,iter, combine_fieldmaps);
	else
	  reconstruction.SmoothFieldmapGroup(m,1-current_group,iter, combine_fieldmaps);

	corrected_stacks.clear();
        for(i=0;i<stacks.size();i++)
          corrected_stacks.push_back(stacks[i]);
        //reconstruction.CorrectStacksSmoothFieldmap(corrected_stacks);
        reconstruction.CorrectStacksSmoothFieldmapWithMasks(corrected_stacks,fieldmapMask1,fieldmapMask2);
	if(current_group==0)
	{
	  correctedMask2=*mask2;
	  reconstruction.CorrectMaskSmoothFieldmap(correctedMask2,m,1);
	}
	else
	{
	  correctedMask1=*mask;
	  reconstruction.CorrectMaskSmoothFieldmap(correctedMask1,m,0);
	}
	//if(iter==2)
	  //exit(1);
	
      //set corrected slices
      reconstruction.UpdateSlices(corrected_stacks,thickness);
      cout<<"PutMask b0"<<endl;
      cout.flush();
      reconstruction.PutMask(b0_mask);
      b0_mask.Write("b0mask.nii.gz");
      cout<<"Mask slices"<<endl;
      cout.flush();
      reconstruction.MaskSlices();
      cout<<" done."<<endl;
      cout.flush();
      

      //redirect output back to screen
      cout.rdbuf (strm_buffer);
      cerr.rdbuf (strm_buffer_e);

    }
      

    //end distortion correction
      
      
      
    //perform slice-to-volume registrations - skip the first iteration 
/*
    if (iter>0)
    {
      cerr.rdbuf(file_e.rdbuf());
      cout.rdbuf (file.rdbuf());
      cout<<"Iteration "<<iter<<": "<<endl;
      
      //reconstruction.SetT2Template(blurredT2);
      
      if((packages.size()>0)&&(iter<=5)&&(iter<(iterations-1)))
      {
	if(iter==1)
          reconstruction.PackageToVolume(corrected_stacks,packages,iter);
	else
	{
	  if(iter==2)
            reconstruction.PackageToVolume(corrected_stacks,packages,iter,true);
	  else
	  {
            if(iter==3)
	      reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true);
	    else
	    {
	      if(iter>=4)
                reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true,iter-2);
	      else
	        reconstruction.SliceToVolumeRegistration();
	    }
	  }
	}
      }
      else
        reconstruction.SliceToVolumeRegistration();
    }
   */
    if (iter>0)
    {
      cerr.rdbuf(file_e.rdbuf());
      cout.rdbuf (file.rdbuf());
      cout<<"Iteration "<<iter<<": "<<endl;
      reconstruction.BSplineReconstruction();
      reconstructed=reconstruction.GetReconstructed();
      sprintf(buffer,"target-for-motion%i.nii.gz",iter);
      reconstructed.Write(buffer);
      //reconstruction.SetT2Template(blurredT2);
      
      if((packages.size()>0)&&(iter<=5)&&(iter<(iterations-2)))
      {
	if(iter==1)
          reconstruction.PackageToVolume(corrected_stacks,packages,iter);
	else
	{
	  if(iter==2)
            reconstruction.PackageToVolume(corrected_stacks,packages,iter,true);
	  else
	  {
            if(iter==3)
	      reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true);
	    else
	    {
	      if(iter>=4)
                reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true,(iter-2)/2);
	      else
	        reconstruction.SliceToVolumeRegistration();
	    }
	  }
	}
      }
      else
        reconstruction.SliceToVolumeRegistration();
    }
  
    //print to screen
    cout.rdbuf (strm_buffer);

    if (current_group == -1)
	current_group = 0;
    else
    {
      current_group++;
      if (current_group >= groups.size())
        current_group = 0;
    }
      
    cout<<endl<<"Current_group = "<<current_group<<endl;
    //if(current_group==0)
    //  reconstruction.SetMask(mask,smooth_mask);
    //else
    //  reconstruction.SetMask(mask2,smooth_mask);
    
    //reconstruction.BSplineReconstructionGroup(groups[current_group]);
    reconstruction.BSplineReconstruction();
    //reconstruction.MaskVolume();
    
    
    if (debug)
    {
      reconstructed=reconstruction.GetReconstructed();
      sprintf(buffer,"image%i.nii.gz",iter);
      reconstructed.Write(buffer);
    }
    cout<<endl;
    cout.flush();
    cerr.rdbuf (strm_buffer_e);
    if(i==1) exit(1);
  }
 
  
  
  //reconstruction.BSplineReconstructionGroup(1);
  reconstruction.BSplineReconstruction();
  if (debug)
  {
    reconstructed=reconstruction.GetReconstructed();
    reconstructed.Write("BsplineImage.nii.gz");
        //sprintf(buffer,"image%i.nii.gz",iter);
        //reconstructed.Write(buffer);
   }
}

////End test
/*
 
  //Create T2 template
  alignedT2= reconstruction.AlignT2Template(t2,0);
  if(debug)
    alignedT2.Write("alignedT2.nii.gz");
  reconstruction.SetT2Template(alignedT2);
  T2_mask=reconstruction.CreateMask(alignedT2);
  if(debug)
    T2_mask.Write("T2mask.nii.gz");


  //Create blurred T2 template for registration
  irtkRealImage blurredT2 = alignedT2;
  irtkGaussianBlurringWithPadding<irtkRealPixel> gb(0.75,0);
  gb.SetInput(&blurredT2);
  gb.SetOutput(&blurredT2);
  gb.Run();
  if(debug)
    blurredT2.Write("blurredT2.nii.gz");

  //prepare larger mask for fieldmap
  reconstruction.CreateLargerMask(T2_mask);


  double step = 2.5;
  if(iterations>1)
    iterations = 7;
  //registration iterations
  for (iter=0;iter<(iterations);iter++)
  {
    //Print iteration number on the screen
    cout.rdbuf (strm_buffer);
    cout<<"Iteration "<<iter<<". "<<endl;
    cout.flush();
    
    if(corrected_stacks.size()==0)
    for(i=0;i<stacks.size();i++)
      corrected_stacks.push_back(stacks[i]);
    
    sprintf(buffer,"test-%i.nii.gz",iter);
    corrected_stacks[0].Write(buffer);
    
    //distortion
    if((iter>=2)||(iterations==1))
    {
      //redirect output to files
      cerr.rdbuf(filed_e.rdbuf());
      cout.rdbuf (filed.rdbuf());
      
      reconstruction.SetT2Template(alignedT2);
      reconstruction.PutMask(T2_mask);
      reconstruction.CoeffInit();

      reconstruction.Shim(corrected_stacks,iter);
      if((iter>2)||(iterations==1)||(iter==(iterations-1)))
      {
        reconstruction.FieldMap(corrected_stacks,step,iter);
	step/=2;
	reconstruction.SaveDistortionTransformations();
      }
        //change 2: Fieldmap is calculated even for affine only
	reconstruction.SmoothFieldmap(iter);

	corrected_stacks.clear();
        for(i=0;i<stacks.size();i++)
          corrected_stacks.push_back(stacks[i]);
        reconstruction.CorrectStacksSmoothFieldmap(corrected_stacks);

      //set corrected slices
      reconstruction.UpdateSlices(corrected_stacks,thickness);
      cout<<"PutMask b0"<<endl;
      cout.flush();
      reconstruction.PutMask(b0_mask);
      b0_mask.Write("b0mask.nii.gz");
      cout<<"Mask slices"<<endl;
      cout.flush();
      reconstruction.MaskSlices();
      cout<<" done."<<endl;
      cout.flush();
      

      //redirect output back to screen
      cout.rdbuf (strm_buffer);
      cerr.rdbuf (strm_buffer_e);

    }

    //perform slice-to-volume registrations - skip the first iteration 
    if (iter>0)
    {
      cerr.rdbuf(file_e.rdbuf());
      cout.rdbuf (file.rdbuf());
      cout<<"Iteration "<<iter<<": "<<endl;
      
      reconstruction.SetT2Template(blurredT2);
      
      if((packages.size()>0)&&(iter<=5)&&(iter<(iterations-1)))
      {
	if(iter==1)
          reconstruction.PackageToVolume(corrected_stacks,packages,iter);
	else
	{
	  if(iter==2)
            reconstruction.PackageToVolume(corrected_stacks,packages,iter,true);
	  else
	  {
            if(iter==3)
	      reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true);
	    else
	    {
	      if(iter>=4)
                reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true,iter-2);
	      else
	        reconstruction.SliceToVolumeRegistration();
	    }
	  }
	}
      }
      else
        reconstruction.SliceToVolumeRegistration();
      
      reconstruction.BSplineReconstruction();
      if (debug)
      {
        reconstructed=reconstruction.GetReconstructed();
        sprintf(buffer,"image%i.nii.gz",iter);
        reconstructed.Write(buffer);
      }
      cout<<endl;
      cout.flush();
      cerr.rdbuf (strm_buffer_e);
    }
  }
    
    //Write to file
    cout.rdbuf (file2.rdbuf());
    cout<<endl<<endl<<"Iteration "<<iter<<": "<<endl<<endl;
    
    //Set smoothing parameters 
    //amount of smoothing (given by lambda) is decreased with improving alignment
    //delta (to determine edges) stays constant throughout
    reconstruction.SetSmoothingParameters(delta,lastIterLambda);
    
    //Use faster reconstruction during iterations and slower for final reconstruction
    reconstruction.SpeedupOff();
    
    //Initialise values of weights, scales and bias fields
    reconstruction.InitializeEMValues();
    
    //Calculate matrix of transformation between voxels of slices and volume
    //if(!coeff_init)
    reconstruction.PutMask(b0_mask);
    reconstruction.CoeffInit();
    
    //Initialize reconstructed image with Gaussian weighted reconstruction
    reconstruction.GaussianReconstruction();
    
    //Simulate slices (needs to be done after Gaussian reconstruction)
    reconstruction.SimulateSlices();
    
    //Initialize robust statistics parameters
    reconstruction.InitializeRobustStatistics();
    
    //EStep
    reconstruction.EStep();

    //number of reconstruction iterations
    rec_iterations = 10;      
        
    //reconstruction iterations
    for (i=0;i<rec_iterations;i++)
    {
      cout<<endl<<"  Reconstruction iteration "<<i<<". "<<endl;
      
      if (intensity_matching)
      {
        //calculate bias fields
        if (sigma>0)
          reconstruction.Bias();
        //calculate scales
        reconstruction.Scale();
      }
      
      //MStep and update reconstructed volume
      reconstruction.Superresolution(i+1);
      
      if (intensity_matching)
      {
        if((sigma>0)&&(!global_bias_correction))
          reconstruction.NormaliseBias(i);
      }

      // Simulate slices (needs to be done
      // after the update of the reconstructed volume)
      reconstruction.SimulateSlices();
      
      reconstruction.MStep(i+1);
      
      //E-step
      reconstruction.EStep();
      
    //Save intermediate reconstructed image
    if (debug)
    {
      reconstructed=reconstruction.GetReconstructed();
      sprintf(buffer,"super%i.nii.gz",i);
      reconstructed.Write(buffer);
    }

      
    }//end of reconstruction iterations
    
    //Mask reconstructed image to ROI given by the mask
    reconstruction.PutMask(T2_mask);
    reconstruction.MaskVolume();


   //Evaluate - write number of included/excluded/outside/zero slices in each iteration in the file
   cout.rdbuf (fileEv.rdbuf());
   reconstruction.Evaluate(iter);
   cout<<endl;
   cout.flush();
   cout.rdbuf (strm_buffer); 
   
  //}// end of interleaved registration-reconstruction iterations

  //save final result
  reconstruction.SaveDistortionTransformations();
  reconstruction.RestoreSliceIntensities();
  reconstruction.ScaleVolume();
  reconstructed=reconstruction.GetReconstructed();
  reconstructed.Write(output_name); 
  reconstruction.SaveTransformations();
  reconstruction.SaveSlices();
  
  //The end of main()
}  
*/