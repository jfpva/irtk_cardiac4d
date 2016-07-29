/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/
#include <vector>
#include <string>
#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkReconstruction.h>
#include <irtkReconstructionb0.h>
#include <irtkReconstructionfMRI.h>

using namespace std;

//Application to perform reconstruction of volumetric MRI from thick slices.

void usage()
{
  cerr << "Usage: reconstruction [reconstructed] [Target] fMRI data <options>\n" << endl;
  cerr << endl;

  cerr << "\t[reconstructed]         Name for the reconstructed volume. Nifti or Analyze format." << endl;
  cerr << "\t-target	                Volume to be used as target (starts from 0)."<<endl;
  cerr << "\tfMRI time serie" << endl;
  cerr << "\t" << endl;
  cerr << "Options:" << endl;
  cerr << "\t-dofin1 [dof_1] .. [dof_N] The transformations of the input stack to template" << endl;
  cerr << "\t                          in \'dof\' format used in IRTK." <<endl;
  cerr << "\t                          Only rough alignment with correct orienation and " << endl;
  cerr << "\t                          some overlap is needed." << endl;
  cerr << "\t                          Firsty, introduce the total number of transformations needed." << endl; 
  cerr << "\t    	  	          Then, the frame number they refer to and dof file (frame by frame). " << endl;
  cerr << "\t-dofin2 [dof_1] .. [dof_N] The transformations of the input stack to template" << endl;
  cerr << "\t                          in \'dof\' format used in IRTK." <<endl;
  cerr << "\t                          Only rough alignment with correct orienation and " << endl;
  cerr << "\t                          some overlap is needed." << endl;
  cerr << "\t                          Firsty, introduce the total number of chunks needed." << endl; 
  cerr << "\t    	  	          Then, first and last frames of each chunk followed by the relative dof file. " << endl;
  cerr << "\t-thickness [th_1] .. [th_N]    Give slice thickness.[Default: twice voxel size in z direction]"<<endl;
  cerr << "\t-mask [mask]              Binary mask to define the region od interest. [Default: whole image]"<<endl;
  cerr << "\t-multiband 		  Multiband factor."<<endl;
  cerr << "\t-packages                 Give number of packages used during acquisition for each stack. [Default: 1]"<<endl;
  cerr << "\t                          The stacks will be split into packages during registration iteration 1"<<endl;
  cerr << "\t                          and then following the specific slice ordering "<<endl;
  cerr << "\t                          from iteration 2. The method will then perform slice to"<<endl;
  cerr << "\t                          volume (or multiband registration)."<<endl;
  cerr << "\t-order                    Slice acquisition order used at acquisition. [Default: ascending (1)]"<<endl;
  cerr << "\t                          Possible values: 2 (descending), 3 (default) 4 (interleaved) and 5 (Customized)."<<endl;
  cerr << "\t-step      		  Forward slice jump for customized (C) slice ordering [Default: 1]"<<endl;
  cerr << "\t-rewinder	          Rewinder for customized slice ordering [Default: 1]"<<endl;
  cerr << "\t-iterations [iter]        Number of registration iterations. [Default is calculated internally]"<<endl;
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
  cerr << "\t-no_robust_statistics     Switch off robust statistics."<<endl;
  cerr << "\t-exclude_slices_only      Do not exclude individual voxels."<<endl;
  cerr << "\t-bspline                  Use multi-level bspline interpolation instead of super-resolution."<<endl;
  cerr << "\t-log_prefix [prefix]      Prefix for the log file."<<endl;
  cerr << "\t-info [filename]          Filename for slice information in\
                                       tab-sparated columns."<<endl;
  cerr << "\t-debug                    Debug mode - save intermediate results."<<endl;
  cerr << "\t-no_log                   Do not redirect cout and cerr to log files."<<endl;
  cerr << "\t" << endl;
  cerr << "\t" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  
  //utility variables
  int i, ok;
  char buffer[256];
  irtkRealImage stack; 
  
  //declare variables for input
  /// Name for output volume
  char * output_name = NULL;
  /// Slice stacks
  vector<irtkRealImage> stacks;
  vector<string> stack_files;
  /// Stack transformation
  vector<irtkRigidTransformation> stack_transformations;
  /// user defined transformations
  bool have_stack_transformations = false;
  /// Stack thickness
  vector<double > thickness;
  ///number of stacks
  int nStacks;
  /// number of packages for each stack
  vector<int> packages;
  vector<int> order_vector;
  vector<int> multiband_vector;
  int step = 1;
  int rewinder = 1;

  // Default values.
  int templateNumber=-1;
  irtkRealImage *mask=NULL;
  int iterations = 0;
  bool debug = false;
  double sigma=20;
  double resolution = 2;
  double lambda = 0.02;
  double delta = 150;
  int levels = 3;
  double lastIterLambda = 0.01;
  int rec_iterations;
  double averageValue = 700;
  double smooth_mask = 4;
  bool global_bias_correction = false;
  double low_intensity_cutoff = 0.01;
  //folder for slice-to-volume registrations, if given
  char * folder=NULL;
  //flag to remove black background, e.g. when neonatal motion correction is performed
  bool remove_black_background = false;
  //flag to swich the intensity matching on and off
  bool intensity_matching = false;
  bool rescale_stacks = false;

  //flag to swich the robust statistics on and off
  bool robust_statistics = false;
  bool robust_slices_only = false;
  //flag to replace super-resolution reconstruction by multilevel B-spline interpolation
  bool bspline = false;
  int multiband_factor = 1;
  
  irtkRealImage average;

  string info_filename = "slice_info.tsv";
  string log_id;
  bool no_log = false;

  //forced exclusion of slices
  int number_of_force_excluded_slices = 0;
  vector<int> force_excluded;

  //Create reconstruction object
  irtkReconstructionfMRI reconstruction;
    
  //if not enough arguments print help
  if (argc < 3)
    usage();
  
  //read output name
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"Recontructed volume name ... "<<output_name<<endl;
 
  templateNumber=atof(argv[1]);
  cout<<"Template Number is ... "<<templateNumber<<endl;
  argc--;
  argv++;
  
  //read 4D image
  irtkRealImage image4D; 
  cout<<"Reading stack ... "<<argv[1]<<endl;
  image4D.Read(argv[1]);
  argc--;
  argv++;
  irtkImageAttributes attr = image4D.GetImageAttributes();
  nStacks = attr._t;
  cout<<"Number 0f stacks ... "<<nStacks<<endl;
	
  // Create stacks 
  for (i = 0; i < nStacks;i++)
  {
		stacks.push_back(image4D.GetRegion(0,0,0,i,attr._x, attr._y,attr._z,i+1));
		sprintf(buffer,"stack%i.nii.gz",i);
		stacks[i].Write(buffer);    
  }
  
  // Parse options.
  while (argc > 1){
    ok = false;

    //Read stack transformations
	if ((ok == false) && (strcmp(argv[1], "-dofin1") == 0)){
		
		argc--;
		argv++;
		
		int quantity = atof(argv[1]);
		argc--;
		argv++;
		
		irtkTransformation *transformation;
		
		int q = 0;
		int minimum = 0;
		while(q < quantity) {
			
			int position = atof(argv[1]);
			argc--;
			argv++;
			
			for (int p = minimum; p < position; p++) {
				// id transformations 
				transformation = new irtkRigidTransformation;
				if ( templateNumber < 0) templateNumber = i;
				irtkRigidTransformation *rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);
				stack_transformations.push_back(*rigidTransf);
				delete rigidTransf;
			}
			
			transformation = irtkTransformation::New(argv[1]);
			argc--;
			argv++;
						
			// actual transformation inserted by user
			irtkRigidTransformation *rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);
			stack_transformations.push_back(*rigidTransf);
			delete rigidTransf;
			
			minimum = position+1;
			q++;
			
			if (q == quantity)	{
				for (int p = minimum; p < nStacks; p++) {
					// id transformations 
					transformation = new irtkRigidTransformation;
					if ( templateNumber < 0) templateNumber = i;
					irtkRigidTransformation *rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);
					stack_transformations.push_back(*rigidTransf);
					delete rigidTransf;
				}
			}
			
		}
		have_stack_transformations = true;	
	}
	
	//Read stack transformations
	if ((ok == false) && (strcmp(argv[1], "-dofin2") == 0)) {
		
		// if dofin1 has not been used, initialize all transformations as id
		if (stack_transformations.size() == 0) {
			for (i=0;i<nStacks;i++)
			{
			  irtkRigidTransformation *rigidTransf = new irtkRigidTransformation;
			  stack_transformations.push_back(*rigidTransf);
			  delete rigidTransf;
			}
		}
		
		argc--;
		argv++;
		
		int quantity = atof(argv[1]);
		argc--;
		argv++;

		int q = 0;
		int min;
		int max;
		
		while(q < quantity)	{
			
			min = atof(argv[1]);
			argc--;
			argv++;
			
			max = atof(argv[1]);
			argc--;
			argv++;
				
			irtkRigidTransformation *rigidTransf;
			irtkTransformation *transformation;
			
			for (int p = min; p <= max; p++) {
				transformation = irtkTransformation::New(argv[1]);
				rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);
				stack_transformations[p] = *rigidTransf;
			}
			q++;			
			argc--;
			argv++;
		}
		have_stack_transformations = true;			
	}
	
    //Read slice thickness
    if ((ok == false) && (strcmp(argv[1], "-thickness") == 0)){
      argc--;
      argv++;
      cout<< "Slice thickness is ";
      for (i=0;i<nStacks;i++)
      {
        thickness.push_back(atof(argv[1]));
        cout<<thickness[i]<<" ";        
      }
      argc--;
      argv++;
      cout<<"."<<endl;
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
      }
      argc--;
      argv++;
      cout<<"."<<endl;
      ok = true;
    }
    
    // Input slice ordering
	if ((ok == false) && (strcmp(argv[1], "-order") == 0)) {
		argc--;
		argv++;
		cout<< "Order is ";
		for (i=0;i<nStacks;i++)
		{
			order_vector.push_back(atoi(argv[1]));		
			cout<<order_vector[i]<<" ";
		}
		argc--;
		argv++;
		cout<<"."<<endl;
		ok = true;
	}
	
	//Multiband factor for each stack
	if ((ok == false) && (strcmp(argv[1], "-multiband") == 0)){
		argc--;
		argv++;
		cout<< "Multiband number is ";
		for (i=0;i<nStacks;i++)
		{
		  multiband_vector.push_back(atoi(argv[1]));
		  cout<<multiband_vector[i]<<" ";
		}
		argc--;
		argv++;
		cout<<"."<<endl;
		ok = true;
	}

    // Forward slice jump for arbitrary slice ordering
	if ((ok == false) && (strcmp(argv[1], "-step") == 0)){
	  argc--;
	  argv++;
	  step=atof(argv[1]);
	  ok = true;
	  argc--;
	  argv++;
	}

	// Rewinder slice jump for arbitrary slice ordering
	if ((ok == false) && (strcmp(argv[1], "-rewinder") == 0)){
	  argc--;
	  argv++;
	  rewinder=atof(argv[1]);
	  ok = true;
	  argc--;
	  argv++;
	}

    //Read binary mask for final volume
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask = new irtkRealImage(argv[1]);
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
    
    //Switch off robust statistics
    if ((ok == false) && (strcmp(argv[1], "-no_robust_statistics") == 0)){
      argc--;
      argv++;
      robust_statistics=false;
      ok = true;
    }
    
    //Switch off robust statistics
    if ((ok == false) && (strcmp(argv[1], "-exclude_slices_only") == 0)){
      argc--;
      argv++;
      robust_slices_only=true;
      ok = true;
    }
    
    //Use multilevel B-spline interpolation instead of super-resolution
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)){
      argc--;
      argv++;
      bspline=true;
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
    
    //Prefix for log files
    if ((ok == false) && (strcmp(argv[1], "-log_prefix") == 0)){
      argc--;
      argv++;
      log_id=argv[1];
      ok = true;
      argc--;
      argv++;
    }

    //No log files
    if ((ok == false) && (strcmp(argv[1], "-no_log") == 0)){
      argc--;
      argv++;
      no_log=true;
      ok = true;
    }

    // rescale stacks to avoid error:
    // irtkImageRigidRegistrationWithPadding::Initialize: Dynamic range of source is too large
    if ((ok == false) && (strcmp(argv[1], "-rescale_stacks") == 0)){
      argc--;
      argv++;
      rescale_stacks=true;
      ok = true;
    }       

    // Save slice info
    if ((ok == false) && (strcmp(argv[1], "-info") == 0)) {
        argc--;
        argv++;
        info_filename=argv[1];
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
       ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

 if (have_stack_transformations == true)
	 reconstruction.InvertStackTransformations(stack_transformations);
  
  if (rescale_stacks)
  {
      for (i=0;i<nStacks;i++)
          reconstruction.Rescale(stacks[i],1000);
  }
  
  //If transformations were not defined by user, set them to identity
  if(!have_stack_transformations)
  {
    for (i=0;i<nStacks;i++)
    {
      irtkRigidTransformation *rigidTransf = new irtkRigidTransformation;
      stack_transformations.push_back(*rigidTransf);
      delete rigidTransf;
    }
  }
  
  // set packages to 1 if not given by user
  if (packages.size() == 0) {
	for (i=0;i<nStacks;i++)	{
		packages.push_back(1);
	}
	cout<<"All packages set to 1"<<endl;
  }
  else // feel with package number from the first stack
  {
	for (i=1;i<nStacks;i++)	{
	  packages.push_back(packages[0]);
	}
  }

  // set multiband to 1 if not given by user
  if (multiband_vector.size() == 0) {
	for (i=0;i<nStacks;i++)	{
		multiband_vector.push_back(1);
	}
	cout<<"Multiband set to 1 for all stacks"<<endl;
  }
  else // feel with multiband from the first stack
  {
	for (i=1;i<nStacks;i++)	{
      multiband_vector.push_back(multiband_vector[0]);
	}
  }
  
  // set order to ascending if not given by user
  if (order_vector.size() == 0) {
  	for (i=0;i<nStacks;i++)	{
      order_vector.push_back(1);
  	}
  	cout<<"Slice ordering set to ascending for all stacks"<<endl;
  }
  else // feel with order from the first stack
  {
  	for (i=1;i<nStacks;i++)	{
        order_vector.push_back(order_vector[0]);
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
  }

  //Output volume
  irtkRealImage reconstructed;

  //Set debug mode
  if (debug) reconstruction.DebugOn();
  else reconstruction.DebugOff();
  
  //Set force excluded slices
  reconstruction.SetForceExcludedSlices(force_excluded);

  //Set low intensity cutoff for bias estimation
  //reconstruction.SetLowIntensityCutoff(low_intensity_cutoff)  ;

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
  //Before creating the template we will crop template stack according to the given mask
  if (mask !=NULL)
  {
    //first resample the mask to the space of the stack
    //for template stact the transformation is identity
    irtkRealImage m = *mask;
    reconstruction.TransformMask(stacks[templateNumber],m,stack_transformations[templateNumber]);
    //Crop template stack
    reconstruction.CropImageIgnoreZ(stacks[templateNumber],m);
    if (debug)
    {
      m.Write("maskTemplate.nii.gz"); 
      stacks[templateNumber].Write("croppedTemplate.nii.gz");
    }
  }
  
  //Create template volume with isotropic resolution 
  //if resolution==0 it will be determined from in-plane resolution of the image
  resolution = reconstruction.CreateTemplate(stacks[templateNumber],resolution);
  
  //Set mask to reconstruction object. 
  reconstruction.SetMask(mask,smooth_mask);   

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
  
  //set precision
  cout<<setprecision(3);
  cerr<<setprecision(3);

  //perform volumetric registration of the stacks
  //redirect output to files
  if ( ! no_log ) {
      cerr.rdbuf(file_e.rdbuf());
      cout.rdbuf (file.rdbuf());
  }
  
  //volumetric registration
  reconstruction.irtkReconstruction::StackRegistrations(stacks,stack_transformations,templateNumber);

  //if remove_black_background flag is set, create mask from black background of the stacks
  if (remove_black_background)
    reconstruction.CreateMaskFromBlackBackground(stacks,stack_transformations, smooth_mask);
  
  cout<<endl;
  //redirect output back to screen
  if ( ! no_log ) {
      cout.rdbuf (strm_buffer);
      cerr.rdbuf (strm_buffer_e);
  }
  
  average = reconstruction.CreateAverage(stacks,stack_transformations);
  if (debug)
    average.Write("average1.nii.gz");

  //Mask is transformed to the all other stacks and they are cropped
  for (i=0; i<nStacks; i++)
  {
    //template stack has been cropped already
    if ((i==templateNumber)&&(!remove_black_background)) continue;
    //transform the mask
    irtkRealImage m=reconstruction.GetMask();
    reconstruction.TransformMask(stacks[i],m,stack_transformations[i]);
    //Crop template stack
    reconstruction.CropImageIgnoreZ(stacks[i],m);
    if (debug)
    {
      sprintf(buffer,"mask%i.nii.gz",i);
      m.Write(buffer); 
      sprintf(buffer,"cropped%i.nii.gz",i);
    }
  }

  // we remove stacks of size 1 voxel (no intersection with ROI)
  vector<irtkRealImage> selected_stacks;
  vector<irtkRigidTransformation> selected_stack_transformations;
  int new_nStacks = 0;
  int new_templateNumber = 0;
  for (i=0; i<nStacks; i++)
  {
      if (stacks[i].GetX() == 1) {
          cerr << "stack " << i << " has no intersection with ROI" << endl;
          continue;
      }

      // we keep it
      selected_stacks.push_back(stacks[i]);
      selected_stack_transformations.push_back(stack_transformations[i]);
      
      if (i == templateNumber)
          new_templateNumber = templateNumber - (i-new_nStacks);

      new_nStacks++;
      
  }
  stacks.clear();
  stack_transformations.clear();
  nStacks = new_nStacks;
  templateNumber = new_templateNumber;
  for (i=0; i<nStacks; i++)
  {
      stacks.push_back(selected_stacks[i]);
      stack_transformations.push_back(selected_stack_transformations[i]);
  }
  
  //Repeat volumetric registrations with cropped stacks
  //redirect output to files
  if ( ! no_log ) {
      cerr.rdbuf(file_e.rdbuf());
      cout.rdbuf (file.rdbuf());
  }
  //volumetric registration
  reconstruction.irtkReconstruction::StackRegistrations(stacks,stack_transformations,templateNumber);
  cout<<endl;

  //redirect output back to screen
  if ( ! no_log ) {
      cout.rdbuf (strm_buffer);
      cerr.rdbuf (strm_buffer_e);
  }
  
  //Rescale intensities of the stacks to have the same average
  if (intensity_matching)
    reconstruction.MatchStackIntensitiesWithMasking(stacks,stack_transformations,averageValue);
  else
    reconstruction.MatchStackIntensitiesWithMasking(stacks,stack_transformations,averageValue,true);
  average = reconstruction.CreateAverage(stacks,stack_transformations);
  if (debug)
    average.Write("average2.nii.gz");
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
  if (folder!=NULL)
    reconstruction.ReadTransformation(folder);
    
  //Initialise data structures for EM
  reconstruction.InitializeEM();
  
  //interleaved registration-reconstruction iterations
  int internal = reconstruction.giveMeDepth(stacks, packages, multiband_vector);
  if (iterations == 0)	{
 	iterations = internal+1;
 	cout<<"Number of iterations is calculated internally: "<<iterations<<endl;
  }
  else if (iterations <= internal)	{
 	iterations = internal+1;
 	cout<<"Number of iterations too small. Iterations are set to :"<<iterations<<endl;
  }
  else {
    cout<<"Number of iterations is :"<<iterations<<endl;
  }
 
  for (int iter=0;iter<iterations;iter++)
  {
	  //Print iteration number on the screen
	  	  if ( ! no_log ) {
	  		  cout.rdbuf (strm_buffer);
	  	  }
	  	  
	  	  cout<<"Iteration"<<iter<<". "<<endl;
	  	  
	  	  if (iter>0)
	  	  {
	  		if ( ! no_log ) {
	  		  cerr.rdbuf(file_e.rdbuf());
	  		  cout.rdbuf (file.rdbuf());
	  		}

	  		vector<int> level;
	  		if(iter == 1) {
	  			reconstruction.newPackageToVolume(stacks, packages, multiband_vector, order_vector, step, rewinder,iter);
	  		}

	  		else if((iter > 1) && (iter < internal-1)){
	  			level = reconstruction.giveMeSplittingVector(stacks, packages, multiband_vector, iter, false);
	  			reconstruction.ChunkToVolume(stacks, packages, level, multiband_vector, order_vector, step, rewinder,iter);
	  		}

	  		else {	
	  			level = reconstruction.giveMeSplittingVector(stacks, packages, multiband_vector, iter, true);
	  			reconstruction.ChunkToVolume(stacks, packages, level, multiband_vector, order_vector, step, rewinder,iter);
	  		}

	  		if ( ! no_log ) {
	  			cerr.rdbuf (strm_buffer_e);
	  		}
	  	  }
	  	  
	  	  if ((iter>0) && (debug))
	  		  reconstruction.SaveRegistrationStep(stacks,iter);
	  	  
	  	  if (iter>0)
	  		  reconstruction.InterpolateGaussian(stacks,iter);
	  
	  //Write to file
	  if ( ! no_log ) {
		 cout.rdbuf (file2.rdbuf());
      }

	  //Set smoothing parameters
	  //amount of smoothing (given by lambda) is decreased with improving alignment
	  //delta (to determine edges) stays constant throughout
	  if(iter==(iterations-1))
		 reconstruction.SetSmoothingParameters(delta,lastIterLambda);
	  else
	  {
		double l=lambda;
		for (i=0;i<levels;i++)
		{
			if (iter==iterations*(levels-i-1)/levels)
			  reconstruction.SetSmoothingParameters(delta, l);
			l*=2;
		}
	  }
	  
	  //Use faster reconstruction during iterations and slower for final reconstruction
	  if ( iter<(iterations-1) )
		  reconstruction.SpeedupOn();
	  else 
		  reconstruction.SpeedupOff();
	  
	  if(robust_slices_only)
		  reconstruction.ExcludeWholeSlicesOnly();

	  //Initialise values of weights, scales and bias fields
	  reconstruction.InitializeEMValues();
	  
	  //Calculate matrix of transformation between voxels of slices and volume
	  reconstruction.SetInterpolationRecon();
	  if (bspline) {
		  reconstruction.CoeffInitBSpline();
	  }
	  else
		  reconstruction.CoeffInit();
	
	  //Initialize reconstructed image with Gaussian weighted reconstruction
	  if (bspline) {
		  reconstruction.BSplineReconstruction();
	  }
	  else
		  reconstruction.GaussianReconstruction();
	
	  //Simulate slices (needs to be done after Gaussian reconstruction)
	  reconstruction.SimulateSlices();

	  //Initialize robust statistics parameters
	  reconstruction.InitializeRobustStatistics();
	
	  //EStep
	  if(robust_statistics)
		  reconstruction.EStep();
	
	  //number of reconstruction iterations
	  if ( iter==(iterations-1) ) 
	  {
		  rec_iterations = 30;      
	  }
	  else 
		  rec_iterations = 10;
	
	  if ((bspline)&&(!robust_statistics)&&(!intensity_matching))
		  rec_iterations=0;
	  
	  // super resolution not needed for fMRI
	  rec_iterations = 0;
	  
	  //reconstruction iterations
	  i=0;
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
	  
	  //Update reconstructed volume
	  if (bspline) {
		  reconstruction.BSplineReconstruction();
	  }
	  else
		  reconstruction.Superresolution(i+1);
	  
	  if (intensity_matching)
	  {
		if((sigma>0)&&(!global_bias_correction))
			reconstruction.NormaliseBias(i);
	  }
	
	  // Simulate slices (needs to be done
	  // after the update of the reconstructed volume)
	  reconstruction.SimulateSlices();
			
	  if(robust_statistics)
		  reconstruction.MStep(i+1);
	  
	  //E-step
	  if(robust_statistics)
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
	if(!bspline)
		reconstruction.MaskVolume();
	
	//Save reconstructed image
	//if (debug)
	//{
	reconstructed=reconstruction.GetReconstructed();
	sprintf(buffer,"image%i.nii.gz",iter);
	reconstructed.Write(buffer);
	//reconstruction.SaveConfidenceMap();
	//}
	
	//Evaluate - write number of included/excluded/outside/zero slices in each iteration in the file
	if ( ! no_log ) {
		cout.rdbuf (fileEv.rdbuf());
	}
	reconstruction.Evaluate(iter);
	cout<<endl;
	
	if ( ! no_log ) {
	    cout.rdbuf (strm_buffer);
	}
	
	}// end of interleaved registration-reconstruction iterations
  
	//save final result
	reconstruction.RestoreSliceIntensities();
	reconstruction.ScaleVolume();
	reconstructed=reconstruction.GetReconstructed();
	reconstructed.Write(output_name); 
	//reconstruction.SaveTransformations();
	reconstruction.SaveSlices();

	// Don't know why it fails here
	/*if ( info_filename.length() > 0 )
	  reconstruction.SlicesInfo( info_filename.c_str(),
								 stack_files );*/
	if(debug)
	{
	reconstruction.SaveWeights();
	reconstruction.SaveBiasFields();
	//reconstruction.SaveConfidenceMap();
	reconstruction.SimulateStacks(stacks);
	for (unsigned int i=0;i<stacks.size();i++)
	{
	  sprintf(buffer,"simulated%i.nii.gz",i);
	  stacks[i].Write(buffer);
	}
  }
  //The end of main()
}
