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
#include <irtkReconstructionCardiac4D.h>

using namespace std;

//Application to perform reconstruction of volumetric cardiac cine MRI from thick-slice dynamic 2D MRI

void usage()
{
  cerr << "Usage: reconstructionCardiac [reconstructed] [N] [stack_1] .. [stack_N] <options>\n" << endl;
  cerr << endl;

  cerr << "NOTE: using temporal PSF = sinc(PI*angdiff/dtrad)*win_Tukey(angdiff,0.3)\n" << endl;

  cerr << "\t[reconstructed]         Name for the reconstructed volume. Nifti or Analyze format." << endl;
  cerr << "\t[N]                     Number of stacks." << endl;
  cerr << "\t[stack_1] .. [stack_N]  The input stacks. Nifti or Analyze format (first taken as reference)." << endl;
  cerr << "\t" << endl;
  cerr << "Options:" << endl;
  cerr << "\t-dofin [dof_1]   .. [dof_N]    The transformations of the input stack to template" << endl;
  cerr << "\t                          in \'dof\' format used in IRTK." <<endl;
  cerr << "\t                          Only rough alignment with correct orienation and " << endl;
  cerr << "\t                          some overlap is needed." << endl;
  cerr << "\t                          Use \'id\' for an identity transformation for at least" << endl;
  cerr << "\t                          one stack. The first stack with \'id\' transformation" << endl;
  cerr << "\t                          will be resampled as template." << endl;
  cerr << "\t-thickness [th_1] .. [th_N]    Give slice thickness.[Default: twice voxel size in z direction]"<<endl;
  cerr << "\t-mask [mask]              Binary mask to define the region od interest. [Default: whole image]"<<endl;
  cerr << "\t-multiband [num_1] .. [num_N]  Multiband factor for each stack for each stack. [Default: 1]"<<endl;
  cerr << "\t-packages [num_1] .. [num_N]   Give number of packages used during acquisition for each stack. [Default: 1]"<<endl;
  cerr << "\t                          The stacks will be split into packages during registration iteration 1"<<endl;
  cerr << "\t                          and then following the specific slice ordering "<<endl;
  cerr << "\t                          from iteration 2. The method will then perform slice to"<<endl;
  cerr << "\t                          volume (or multiband registration)."<<endl;
  cerr << "\t-order                    Vector of slice acquisition orders used at acquisition. [Default: (1)]"<<endl;
  cerr << "\t                          Possible values: 1 (ascending), 2 (descending), 3 (default), 4 (interleaved)"<<endl;
  cerr << "\t                          and 5 (Customized)."<<endl;
  cerr << "\t-step       	          Forward slice jump for customized (C) slice ordering [Default: 1]"<<endl;
  cerr << "\t-rewinder	          Rewinder for customized slice ordering [Default: 1]"<<endl;
  cerr << "\t-cardphase [K] [num_1] .. [num_K]  Cardiac phase (0-2PI) for each of K slices. [Default: 0]."<<endl;
  cerr << "\t-numcardphase             Number of cardiac phases to reconstruct. [Default: 10]."<<endl;
  cerr << "\t-rrinterval [rr]          R-R interval. [Default: 1 s]."<<endl;
  cerr << "\t-iterations [iter]        Number of registration-reconstruction iterations. [Default: calc. internally]"<<endl;
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
  cerr << "\tNOTE: work in progress, use of following options is not recommended..." << endl;
  cerr << "\t\tmultiband, packages, order, step, rewinder, rescale_stacks" << endl;
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
  const double PI = 3.14159265358979323846;
  
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
  // Slice R-R Intervals
  vector<double> rr;
  // Slice cardiac phases
  vector<double> cardPhase;
  
  int step = 1;
  int rewinder = 1;
  
  // Default values.
  int templateNumber=-1;
  irtkRealImage *mask=NULL;
  int iterations = 0;
  bool debug = false;
  double sigma=20;
  double resolution = 0.75;
  int numCardPhase = 10;
  double rrInterval = 1;
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
  char *folder=NULL;
  //flag to remove black background, e.g. when neonatal motion correction is performed
  bool remove_black_background = false;
  //flag to swich the intensity matching on and off
  bool intensity_matching = true;
  bool rescale_stacks = false;

  //flag to swich the robust statistics on and off
  bool robust_statistics = true;
  bool robust_slices_only = false;
  //flag to replace super-resolution reconstruction by multilevel B-spline interpolation
  bool bspline = false;
  vector<int> multiband_vector;
  int multiband_factor=1;
  
  irtkRealImage average;

  string info_filename = "slice_info.tsv";
  string log_id;
  bool no_log = false;

  //forced exclusion of slices
  int number_of_force_excluded_slices = 0;
  vector<int> force_excluded;

  //Create reconstruction object
  irtkReconstructionCardiac4D reconstruction;
    
  //if not enough arguments print help
  if (argc < 5)
    usage();
  
  //read output name
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"Recontructed volume name ... "<<output_name<<endl;

  //read number of stacks
  nStacks = atoi(argv[1]);
  argc--;
  argv++;
  cout<<"Number 0f stacks ... "<<nStacks<<endl;

  // Read stacks 
  for (i=0;i<nStacks;i++)
  {
    //if ( i == 0 )
    //log_id = argv[1];
    stack_files.push_back(argv[1]);
    stack.Read(argv[1]);
    cout<<"Reading stack ... "<<argv[1]<<endl;
    argc--;
    argv++;
    stacks.push_back(stack);
  }
  
  // Parse options.
  while (argc > 1){
    ok = false;
    
    //Read stack transformations
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      
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

        argc--;
        argv++;
        irtkRigidTransformation *rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);
        stack_transformations.push_back(*rigidTransf);
        delete rigidTransf;
      }
      reconstruction.InvertStackTransformations(stack_transformations);
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
        argc--;
        argv++;
       }
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
        argc--;
        argv++;
       }
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
			argc--;
			argv++;
	    }
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
		  argc--;
		  argv++;
	    }
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

  //Read cardiac phases
  if ((ok == false) && (strcmp(argv[1], "-cardphase") == 0)){
    argc--;
    argv++;
    int nSlices = atoi(argv[1]);
    cout<<"Reading cardiac phase for "<<nSlices<<" slices/frames"<<endl;
    argc--;
    argv++;
    cout<< "Cardiac phase values are ";
    for (i=0;i<nSlices;i++)
    {
      cardPhase.push_back(atof(argv[1]));
      cout<<i<<":"<<cardPhase[i]<<", ";
      argc--;
      argv++;
    }
    cout<<"\b\b."<<endl;
    ok = true;
  }

  // Number of cardiac phases in reconstructed volume
	if ((ok == false) && (strcmp(argv[1], "-numcardphase") == 0)){
	  argc--;
	  argv++;
	  numCardPhase=atof(argv[1]);
	  argc--;
	  argv++;
    cout<<"Reconstructing "<<numCardPhase<<" cardiac phases."<<endl;
    ok = true;
	}

  // R-R Interval of Reconstructed Volume
	if ((ok == false) && (strcmp(argv[1], "-rrinterval") == 0)){
	  argc--;
	  argv++;
	  rrInterval=atof(argv[1]);
	  argc--;
	  argv++;
    //cout<<"R-R interval of reconstructed volume is "<<rrInterval<<" s."<<endl;
    ok = true;
    reconstruction.SetReconstructedRRInterval(rrInterval);
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
      cout << "No intensity matching."<<endl;
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

  if (rescale_stacks)
  {
      for (i=0;i<nStacks;i++)
          reconstruction.Rescale(stacks[i],1000);
  }
  
  // set packages to 1 if not given by user
  if (packages.size() == 0)
  	  for (i=0;i<nStacks;i++)	{
  		  packages.push_back(1);
  		  cout<<"All packages set to 1"<<endl;
  	  }
  
  // set multiband to 1 if not given by user
  if (multiband_vector.size() == 0)
	  for (i=0;i<nStacks;i++)	{
		  multiband_vector.push_back(1);
		  cout<<"Multiband set to 1 for all stacks"<<endl;
	  }
  
  // set ascending if not given by user
  if (order_vector.size() == 0)
	  for (i=0;i<nStacks;i++)	{
		  order_vector.push_back(1);	
		  cout<<"Slice order set to ascending for all stacks"<<endl;
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
    templateNumber = 0;  
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
  irtkRealImage volumeweights;
  vector<double> reconstructedCardPhase;
  cout<<setprecision(3);
  cout<<"Reconstructing "<<numCardPhase<<" cardiac phases: ";
  for (i=0;i<numCardPhase;i++)
  {
    reconstructedCardPhase.push_back(2*PI*i/numCardPhase);
    cout<<" "<<reconstructedCardPhase[i]/PI<<",";
  }
  cout<<"\b x PI."<<endl;
  reconstruction.SetReconstructedCardiacPhase( reconstructedCardPhase );
  reconstruction.SetReconstructedTemporalResolution( rrInterval/numCardPhase );
  
  //Set debug mode
  if (debug) reconstruction.DebugOn();
  else reconstruction.DebugOff();
  
  // For now, set R-R for each slice/frame to reconstructed R-R Interval
  // TODO: set R-R interval of each slice/frame based on cardiac trigger times
  for (i=0;i<reconstructedCardPhase.size();i++)
  {
    rr.push_back(rrInterval);
  }
  reconstruction.SetSliceRRInterval(rr);
  
  // Investigate Temporal Weight Calculation
  // reconstruction.TestTemporalWeightCalculation();
  
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
  
  // Initialise Reconstructed Volume
  // Check that mask is provided
  if (mask==NULL)
  {
    cerr<<"Reconstruction of volumetric cardiac cine MRI from thick-slice dynamic 2D MRI requires mask to initilise reconstructed volume."<<endl;
    exit(1);
  }  
  /* Unused if mask required
  //If no mask was given and flag "remove_black_background" is false, try to create mask from the template image in case it was padded
  if ((mask==NULL)&&(!remove_black_background))
  {
    mask = new irtkRealImage(stacks[templateNumber]);
    *mask = reconstruction.CreateMask(*mask);
  }
  */
  // Crop mask
  irtkRealImage maskCropped = *mask;
  reconstruction.CropImage(maskCropped,*mask);  // TODO: TBD: use CropImage or CropImageIgnoreZ
  // Initilaise reconstructed volume with isotropic resolution 
  // if resolution==0 it will be determined from in-plane resolution of the image
  if (resolution <= 0) 
  {
    resolution = reconstruction.GetReconstructedResolutionFromTemplateStack( stacks[templateNumber] );
  }
  if (debug)
    cout << "Initialising volume with isotropic voxel size " << resolution << "mm" << endl;
  
  // Create template 4D volume
  reconstruction.CreateTemplateCardiac4DFromStaticMask( maskCropped, resolution );
  
  // Set mask to reconstruction object
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

  //redirect output to files
  if ( ! no_log ) {
      cerr.rdbuf(file_e.rdbuf());
      cout.rdbuf (file.rdbuf());
  }
  
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

  //Mask is transformed to the all stacks and they are cropped
  for (i=0; i<nStacks; i++)
  {
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
      stacks[i].Write(buffer);
    }
  }

  //Rescale intensities of the stacks to have the same average
  if (intensity_matching)
    reconstruction.MatchStackIntensitiesWithMasking(stacks,stack_transformations,averageValue);
  else
    reconstruction.MatchStackIntensitiesWithMasking(stacks,stack_transformations,averageValue,true);
  average = reconstruction.CreateAverage(stacks,stack_transformations);
  if (debug)
    average.Write("average2.nii.gz");

  //Create slices and slice-dependent transformations
  reconstruction.CreateSlicesAndTransformationsCardiac4D(stacks,stack_transformations,thickness);
  
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
  
  // Calculate Cardiac Phase of Each Slice
  if ( cardPhase.size() == 0 )
  {
    cerr<<"Cardiac 4D reconstruction requires cardiac phase for each slice."<<endl;
    exit(1);
  }
  else
  {
    if(debug)
      cout<<"SetSliceCardiacPhase"<<endl;
    reconstruction.SetSliceCardiacPhase( cardPhase );
  }  
  // Calculate Target Cardiac Phase in Reconstructed Volume for Slice-To-Volume Registration
  reconstruction.CalculateSliceToVolumeTargetCardiacPhase();
  // Calculate Temporal Weight for Each Slice
  reconstruction.CalculateSliceTemporalWeights();  
  
    
  //interleaved registration-reconstruction iterations
  
  if (iterations == 0)	
      iterations = 3;
  if(debug)
      cout<<"Number of iterations is :"<<iterations<<endl;

  for (int iter=0;iter<iterations;iter++)
  {
    //Print iteration number on the screen
    if ( ! no_log ) {
      cout.rdbuf (strm_buffer);
    }
    
    cout<<"Iteration"<<iter<<". "<<endl;

    //perform slice-to-volume registrations
    if ( iter > 0 )
    {
      if ( ! no_log ) {
  		  cerr.rdbuf(file_e.rdbuf());
  		  cout.rdbuf (file.rdbuf());
  		}
      cout<<endl<<endl<<"Iteration "<<iter<<": "<<endl<<endl;
      reconstruction.SliceToVolumeRegistrationCardiac4D();
      cout<<endl;
      if ( ! no_log ) {
          cerr.rdbuf (strm_buffer_e);
      }

    // if ((iter>0) && (debug))
  	// 	  reconstruction.SaveRegistrationStep(stacks,iter); 

      if ( ! no_log ) {
  			cerr.rdbuf (strm_buffer_e);
  		}

    }  // if ( iter > 0 ) 
    
    //Write to file
	  if ( ! no_log ) {
      cout.rdbuf (file2.rdbuf());
    }
	  cout<<endl<<endl<<"Iteration "<<iter<<": "<<endl<<endl;
	
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
      
    //Exclude whole slices only
    if(robust_slices_only)
		  reconstruction.ExcludeWholeSlicesOnly();
  	
	  //Initialise values of weights, scales and bias fields
	  reconstruction.InitializeEMValues();
    
    //Calculate matrix of transformation between voxels of slices and volume
    if (bspline)
    {
      cerr<<"Cannot currently initalise b-spline for cardiac 4D reconstruction."<<endl;
      exit(1);
      // TODO: reconstruction.CoeffInitBSplineCardiac4D();
    }
    else
      reconstruction.CoeffInitCardiac4D();

    //Initialize reconstructed image with Gaussian weighted reconstruction
    if (bspline)
    {
      cerr<<"Cannot currently reconstruct b-spline for cardiac 4D reconstruction."<<endl;
      exit(1);
      // TODO: TBD: reconstruction.BSplineReconstructionCardiac4D();
    }
    else
      reconstruction.GaussianReconstructionCardiac4D();
    
    // Save Initialised Volume to File
    if (debug)
    {
      reconstructed = reconstruction.GetReconstructedCardiac4D();
      sprintf(buffer,"init%i.nii.gz",iter);
      reconstructed.Write(buffer);
      volumeweights = reconstruction.GetVolumeWeights();
      sprintf(buffer,"volume_weights%i.nii.gz",iter);
      volumeweights.Write(buffer);
    }
    
    //Simulate slices (needs to be done after Gaussian reconstruction)
    reconstruction.SimulateSlicesCardiac4D();
      
    //Initialize robust statistics parameters
    reconstruction.InitializeRobustStatistics();

    //EStep
    if(robust_statistics)
      reconstruction.EStep();
        
    //number of reconstruction iterations
    if ( iter==(iterations-1) ) 
    {
      // rec_iterations = 30;
      rec_iterations = 30;  // TBD: optimise number of reconstruction iterations    
    }
    else 
      // rec_iterations = 10;  
      rec_iterations = 10;  // TBD: optimise number of reconstruction iterations    
    
    if ((bspline)&&(!robust_statistics)&&(!intensity_matching))
      rec_iterations=0;
    
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
      if (!bspline)
        reconstruction.SuperresolutionCardiac4D(i+1);
      
      if (intensity_matching)
      {
      if((sigma>0)&&(!global_bias_correction))
        reconstruction.NormaliseBiasCardiac4D(i);
      }
      
      // Simulate slices (needs to be done
      // after the update of the reconstructed volume)
      reconstruction.SimulateSlicesCardiac4D();
      
      if(robust_statistics)
        reconstruction.MStep(i+1);
      
      //E-step
      if(robust_statistics)
        reconstruction.EStep();
      
      //Save intermediate reconstructed image
      if (debug)
      {
        reconstructed=reconstruction.GetReconstructedCardiac4D();
        sprintf(buffer,"super%i_%i.nii.gz",iter,i);
        reconstructed.Write(buffer);
      }    
      
    }//end of reconstruction iterations
  
    //Mask reconstructed image to ROI given by the mask
    if(!bspline)
      reconstruction.StaticMaskReconstructedVolume4D();

    //Save reconstructed image
    reconstructed=reconstruction.GetReconstructedCardiac4D();
    sprintf(buffer,"image%i.nii.gz",iter);
    reconstructed.Write(buffer);

    //Evaluate - write number of included/excluded/outside/zero slices in each iteration in the file
    if ( ! no_log )
      cout.rdbuf (fileEv.rdbuf());
    reconstruction.Evaluate(iter);
    cout<<endl;
    if ( ! no_log )
      cout.rdbuf (strm_buffer);

    if(debug)
    {
      sprintf(buffer,"slice_info_%i.tsv",iter);
      reconstruction.SlicesInfoCardiac4D( buffer, stack_files );
    }

  }// end of interleaved registration-reconstruction iterations

  //save final result
  if(debug)
      cout<<"RestoreSliceIntensities"<<endl;
	reconstruction.RestoreSliceIntensities();
  if(debug)
      cout<<"ScaleVolume"<<endl;
  reconstruction.ScaleVolume();
  if(debug)
      cout<<"Saving Reconstructed Volume"<<endl;
	reconstructed=reconstruction.GetReconstructedCardiac4D();
	reconstructed.Write(output_name); 
  if(debug)
      cout<<"SaveTransformations"<<endl;
	reconstruction.SaveTransformations();
  if(debug)
      cout<<"SaveSlices"<<endl;
	reconstruction.SaveSlices();
  
  /*if(debug)
      cout<<"SaveTransformationsWithTiming"<<endl;
	reconstruction.SaveTransformationsWithTiming();*/
  /*if(debug) 
      cout<<"SaveSlicesWithTiming"<<endl;
	reconstruction.SaveSlicesWithTiming();*/

	if ( info_filename.length() > 0 ) 
  {
    if(debug)
        cout<<"SlicesInfoCardiac4D"<<endl;
    reconstruction.SlicesInfoCardiac4D( info_filename.c_str(),
             								 stack_files );
  }

	if(debug)
	{
      if(debug)
          cout<<"SaveWeights"<<endl;
	    reconstruction.SaveWeights();
      if(debug)
          cout<<"SaveBiasFields"<<endl;
      reconstruction.SaveBiasFields();
      if(debug)
          cout<<"SaveSimulatedSlices"<<endl;
      reconstruction.SaveSimulatedSlices();
	    //reconstruction.SaveConfidenceMap();
  	  /* TODO: update stack simulation for 4D if required 
      reconstruction.SimulateStacks(stacks);
  	  for (unsigned int i=0;i<stacks.size();i++)
  	  {
  	      sprintf(buffer,"simulated%i.nii.gz",i);
  	      stacks[i].Write(buffer);
  	  }
      */
      cout<<"Superheart reconstruction complete."<<endl;
  }  
  //The end of main()
}  