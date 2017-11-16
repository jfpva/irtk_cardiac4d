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
// #include <irtkImageFunction.h>
// #include <irtkResampling.h>
// #include <irtkResamplingWithPadding.h>

using namespace std;

//Application to simulate the acquisition of stacks of dynamic 2D MRI of the heart from a cine volume

void usage()
{
  cerr << "Usage: simulateStacksCardiac [cine_volume] [N] [stack_1] .. [stack_N] <options>\n" << endl;
  cerr << endl;

  cerr << "NOTE: using temporal PSF = sinc(PI*angdiff/dtrad)*win_Tukey(angdiff,0.3)\n" << endl;

  cerr << "\t[cine_volume]              Name for the cine volume. Nifti or Analyze format." << endl;
  cerr << "\t[N]                        Number of stacks." << endl;
  cerr << "\t[stack_1] .. [stack_N]     The input stacks. Nifti or Analyze format." << endl;
  cerr << "\t" << endl;
  cerr << "Options:" << endl;
  cerr << "\t-thickness [th_1]..[th_N]  Give slice thickness.[Default: twice voxel size in z direction]"<<endl;
  cerr << "\t-mask [mask]               Binary mask to define the region of interest. [Default: whole image]"<<endl;
  cerr << "\t-slice_transformations [folder]  Use existing slice-location transformations to initialize the reconstruction."<<endl;
  cerr << "\t-transformations [folder]  Use existing image-frame to volume transformations to initialize the reconstruction."<<endl;
  cerr << "\t-motion_sigma [sigma]      Stdev for smoothing transformations. [Default: 0s, no smoothing]"<<endl;
  cerr << "\t-motion_scale [scale]      Scale existing transformations. [Default: 1, no scaling]"<<endl;
  cerr << "\t-1d                        Perform simulation in through-plane direction only." << endl;
  cerr << "\t-cardphase [K] [num_1]..[num_K]  Cardiac phase (0-2PI) for each image-frames 1-K. [Default: 0]."<<endl;
  cerr << "\t-rrinterval [rr]           R-R interval. [Default: read from cine_volume]."<<endl;
  cerr << "\t-rrintervals [L] [rr_1]..[rr_L]  R-R interval for slice-locations 1-L in input stacks. [Default: 1 s]."<<endl;
  cerr << "\t-simulated_resolution [res_1]..[res_N]  In-plane resolution of simulated stacks. [Default: save voxel size as input stacks.]"<<endl;
  cerr << "\t-smooth_mask [sigma]       Smooth the mask to reduce artefacts of manual segmentation. [Default: 4mm]"<<endl;
  cerr << "\t-force_exclude [n] [ind1]..[indN]  Force exclusion of image-frames with these indices."<<endl;
  cerr << "\t-force_exclude_stack [n] [ind1]..[indN]  Force exclusion of stacks with these indices."<<endl;
  cerr << "\t-speedup                   Use faster, but lower quality reconstruction."<<endl;
  cerr << "\t-debug                     Debug mode - save intermediate results."<<endl;
  cerr << "\t" << endl;

  cerr << "\t" << endl;
  cerr << "\t" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  //utility variables
  int i, ok;
  irtkRealImage stack; 
  const double PI = 3.14159265358979323846;
  
  //declare variables for input
  /// Name of cine volume
  char * cine_name = NULL;
  /// Slice stacks
  vector<irtkRealImage> stacks;
  vector<irtkRealImage> current_stacks;
  vector<string> stack_files;
  /// Stack transformation
  vector<irtkRigidTransformation> stack_transformations;
  /// Stack thickness
  vector<double > thickness;
  ///number of stacks
  int nStacks;
  /// number of packages for each stack
  vector<int> packages;
  vector<int> order_vector;
  // Location R-R Intervals;
  // Slice R-R Intervals
  vector<double> rr_loc;
  vector<double> rr;
  // Slice cardiac phases
  vector<double> cardPhase;
  
  vector <double> simulated_resolution;
  
  // Default values.
  int templateNumber = 0;
  irtkRealImage *mask=NULL;
  bool debug = false;
  int numCardPhase = 0;
  double rrDefault = 0;
  double rrInterval = rrDefault;
  double motion_sigma = 0;
  double motion_scale = 1;
  double smooth_mask = 4;

  //folder for slice-location registrations, if given
  char *slice_transformations_folder=NULL;
  //folder for slice-to-volume registrations, if given
  char *folder=NULL;
  //flag to simulate through-plane only
  bool recon_1d = false;
  vector<int> multiband_vector;
  //flag to use faster lower quality reconstruction
  bool speedup = false;
  
  irtkRealImage average;

  //forced exclusion of slices
  int number_of_force_excluded_slices = 0;
  vector<int> force_excluded;
  int number_of_force_excluded_stacks = 0;
  vector<int> force_excluded_stacks;
  vector<bool> stack_excluded;
  vector<int> current_force_excluded_stacks;

  //Create reconstruction object
  irtkReconstructionCardiac4D reconstruction;
  
  //if not enough arguments print help
  if (argc < 5)
    usage();
  
  //read cine volume name
  cine_name = argv[1];
  argc--;
  argv++;
  cout<<"Recontructed volume name ... "<<cine_name<<endl;

  //read number of stacks
  nStacks = atoi(argv[1]);
  argc--;
  argv++;
  cout<<"Number of stacks ... "<<nStacks<<endl;

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
    
  //Read simulated resolution
  if ((ok == false) && (strcmp(argv[1], "-simulated_resolution") == 0)){
    argc--;
    argv++;
    cout<< "Resolution of simulated stacks will be ";
    for (i=0;i<nStacks;i++)
    {
      simulated_resolution.push_back(atof(argv[1]));
      cout<<" "<<simulated_resolution[i];
      argc--;
      argv++;
     }
     cout<<"."<<endl;
     ok = true;
  }

  //Read stack location R-R Intervals
  if ((ok == false) && (strcmp(argv[1], "-rrintervals") == 0)){
    argc--;
    argv++;
    int nLocs = atoi(argv[1]);
    cout<<"Reading R-R intervals for "<<nLocs<<" slice locations"<<endl;
    argc--;
    argv++;
    cout<< "R-R intervals are ";
    for (i=0;i<nLocs;i++)
    {
      rr_loc.push_back(atof(argv[1]));
      cout<<i<<":"<<rr_loc[i]<<", ";
      argc--;
      argv++;
    }
    cout<<"\b\b."<<endl;
    ok = true;
  }

  //Read cardiac phases
  if ((ok == false) && (strcmp(argv[1], "-cardphase") == 0)){
    argc--;
    argv++;
    int nSlices = atoi(argv[1]);
    cout<<"Reading cardiac phase for "<<nSlices<<" images."<<endl;
    argc--;
    argv++;
    for (i=0;i<nSlices;i++)
    {
      cardPhase.push_back(atof(argv[1]));
      argc--;
      argv++;
    }
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
  
    //Perform reconstruction only in z direction
    if ((ok == false) && (strcmp(argv[1], "-1d") == 0)){
      argc--;
      argv++;
      recon_1d = true;
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
    
    //Debug mode
    if ((ok == false) && (strcmp(argv[1], "-debug") == 0)){
      argc--;
      argv++;
      debug=true;
      ok = true;
    }

    //Read slice-location transformations from this folder
    if ((ok == false) && (strcmp(argv[1], "-slice_transformations") == 0)){
      argc--;
      argv++;
      slice_transformations_folder=argv[1];
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

    //Variance of Gaussian kernel to smooth the motion
    if ((ok == false) && (strcmp(argv[1], "-motion_sigma") == 0)){
      argc--;
      argv++;
      motion_sigma=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    } 
    
    //Factor to scale the motion
    if ((ok == false) && (strcmp(argv[1], "-motion_scale") == 0)){
      argc--;
      argv++;
      motion_scale=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
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

    //Force removal of certain stacks
    if ((ok == false) && (strcmp(argv[1], "-force_exclude_stack") == 0)){
      argc--;
      argv++;
      number_of_force_excluded_stacks = atoi(argv[1]);
      argc--;
      argv++;

      cout<< number_of_force_excluded_stacks<< " force excluded stacks: ";
      for (i=0;i<number_of_force_excluded_stacks;i++)
      {
        force_excluded_stacks.push_back(atoi(argv[1]));
        cout<<force_excluded_stacks[i]<<" ";
        argc--;
        argv++;
       }
       cout<<"."<<endl;
       ok = true;
    }
    
    //Use faster, but lower quality reconstruction
    if ((ok == false) && (strcmp(argv[1], "-speedup") == 0)){
      argc--;
      argv++;
      speedup = true;
      ok = true;
    }
    
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // check that conflicting transformation folders haven't been given
  if ((folder!=NULL)&(slice_transformations_folder!=NULL))
  {
      cerr << "Can not use both -transformations and -slice_transformations arguments." << endl;
      exit(1);
  }

  // set packages to 1
  if (packages.size() == 0)
  	  for (i=0;i<nStacks;i++)	{
  		  packages.push_back(1);
  	  }
  
  // set multiband to 1
  if (multiband_vector.size() == 0)
	  for (i=0;i<nStacks;i++)	{
		  multiband_vector.push_back(1);
	  }
  
  // set ascending
  if (order_vector.size() == 0)
	  for (i=0;i<nStacks;i++)	{
		  order_vector.push_back(1);	
	  }
  
  //Set stack transformations to identity
  irtkRigidTransformation id;
  for (i=0;i<nStacks;i++)
    stack_transformations.push_back(id);
   
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

  //Set debug mode
  if (debug) reconstruction.DebugOn();
  else reconstruction.DebugOff();
  
  //Use high quality factor
  if (speedup)
    reconstruction.SpeedupOn();
  else
    reconstruction.SpeedupOff();

  //Cine volume
  irtkRealImage reconstructed;
  irtkRealImage volumeweights;
  cout<<"Reading cine volume ..." << cine_name << endl;
  reconstructed.Read(cine_name);
  reconstruction.SetReconstructedCardiac4D(reconstructed);
  irtkImageAttributes reconstructedAttr = reconstructed.GetImageAttributes();
  numCardPhase = reconstructedAttr._t;
  vector<double> reconstructedCardPhase;
  if ( rrInterval==rrDefault )
    rrInterval = numCardPhase * reconstructedAttr._dt;
  cout<<"Cine volume has R-R interval "<<rrInterval<<" ms." << endl;    
  cout<<setprecision(3);
  cout<<"Cine volume has "<<numCardPhase<<" cardiac phases: ";
  for (i=0;i<numCardPhase;i++)
  {
    reconstructedCardPhase.push_back(2*PI*i/numCardPhase);
    cout<<" "<<reconstructedCardPhase[i]/PI<<",";
  }
  cout<<"\b x PI."<<endl;
  reconstruction.SetReconstructedCardiacPhase( reconstructedCardPhase );
  reconstruction.SetReconstructedTemporalResolution( rrInterval/numCardPhase );
  
  //initialise stack factors
  reconstruction.InitStackFactor(stacks);
  
  //type of recon - default 3D PSF
  if(recon_1d) {
    cout<<"Set1DRecon()" << endl;
    reconstruction.Set1DRecon();
  }
  
  //Set force excluded slices
  if (debug)
    cout<<"SetForceExcludedSlices()" << endl;
  reconstruction.SetForceExcludedSlices(force_excluded);
  
  //Set force excluded stacks
  if (debug)
    cout<<"SetForceExcludedStacks()" << endl;
  reconstruction.SetForceExcludedStacks(force_excluded_stacks);
  // Identify excluded stacks
  bool isStackExcluded;
  for (unsigned int stackIndex = 0; stackIndex<stacks.size(); stackIndex++){
    isStackExcluded=false;
    for (unsigned int excludedStackIndex = 0; excludedStackIndex<force_excluded_stacks.size(); excludedStackIndex++) {
      if (int(stackIndex)==force_excluded_stacks[excludedStackIndex]) {
        isStackExcluded=true;
      }
    }
    stack_excluded.push_back(isStackExcluded);
  } 
  
  // Set mask
  if (!(mask==NULL)) {
    if (debug)
      cout<<"SetMask()" << endl;
    reconstruction.SetMask(mask,smooth_mask);
  }

  //Set precision
  cout<<setprecision(3);
  cerr<<setprecision(3);
  
  //Resize stacks
  if (simulated_resolution.size()!=0) {
    for (unsigned int stackIndex = 0; stackIndex < stacks.size(); stackIndex++) {
      if (!stack_excluded[stackIndex]) {
          if(debug)
            cout<<"Resizing stack "<<stackIndex<<" with in-plane resolution "<<simulated_resolution[stackIndex]<<" mm."<<endl;
          stack = stacks[stackIndex];
          irtkImageAttributes attr = stack.GetImageAttributes();
          irtkResampling<irtkRealPixel> resampling(simulated_resolution[stackIndex],simulated_resolution[stackIndex],attr._dz);
          irtkNearestNeighborInterpolateImageFunction interpolator;
          resampling.SetInput(&stacks[stackIndex]);
          resampling.SetOutput(&stack);
          resampling.SetInterpolator(&interpolator);
          resampling.Run();
      }
    stacks[stackIndex] = stack;
    }    
  }
  
  //Create slices, required before initialising timing values
  reconstruction.CreateSlicesAndTransformationsCardiac4D(stacks,stack_transformations,thickness);
  
  //Read transformations
  if (folder!=NULL) {
    if (debug)
      cout<<"ReadTransformation()"<<endl;
    reconstruction.ReadTransformation(folder);  // image-frame to volume registrations
  }
  else { 
    if (slice_transformations_folder!=NULL) {    // slice-location to volume registrations
      if (debug)
        cout<<"ReadSliceTransformation()"<<endl;
      reconstruction.ReadSliceTransformation(slice_transformations_folder);
    }
  }
  
  // Set R-R for each image
  if (rr_loc.empty()) 
  {
    if (debug)
      cout<<"SetSliceRRInterval()"<<endl;
    reconstruction.SetSliceRRInterval(rrInterval);
    if (debug)
      cout<<"No R-R intervals specified. All R-R intervals set to "<<rrInterval<<" s."<<endl;
  }  
  else {
    if (debug)
      cout<<"SetLocRRInterval()"<<endl;
    reconstruction.SetLocRRInterval(rr_loc);
  }
  
  // Calculate Cardiac Phase of Each Slice
  if ( cardPhase.size() == 0 ) {  // no cardiac phases specified
    if ( numCardPhase != 1 ) {    // simulating static volume
      irtkImageAttributes attr = stacks[templateNumber].GetImageAttributes();
      if (attr._t > 1) {
        cerr<<"Cardiac simulation requires cardiac phase for each slice."<<endl;
        exit(1);
      }
    }
    else {                        // simulating single cardiac phase volume
      reconstruction.SetSliceCardiacPhase();    // set all cardiac phases to zero
    }
  }
  else {
    reconstruction.SetSliceCardiacPhase( cardPhase );   // set all cardiac phases to given values
  }
  
  // Calculate Temporal Weight for Each Slice
  reconstruction.CalculateSliceTemporalWeights();  
    
  //Smooth transformations
  if(motion_sigma>0)
    reconstruction.SmoothTransformations(motion_sigma);
  
  //Scale transformations
  if(motion_scale!=1)
    reconstruction.ScaleTransformations(motion_scale);
 
  //Initialise data structures for EM
  reconstruction.InitializeEM();
  
  //Initialise values of weights, scales and bias fields
  reconstruction.InitializeEMValues();
    
  // Simulate stacks
  reconstruction.SimulateStacksCardiac4D(stack_excluded);

  // Save simulated stacks
  for (unsigned int stackIndex = 0; stackIndex < stacks.size(); stackIndex++) {
    if (!stack_excluded[stackIndex]) {
      if(debug)
        cout<<"Saving simulated stack "<<stackIndex<<"."<<endl;
      reconstruction.SaveSimulatedSlices(stacks,stackIndex);
    }  
  }  
  
  //Calculate displacements
  double mean_displacement = reconstruction.CalculateDisplacement();
  if (debug)
    cout<<"\tmean displacement = "<<mean_displacement<<" mm."<<endl;
  double mean__weighted_displacement = reconstruction.CalculateWeightedDisplacement();
    if (debug)
      cout<<"\tmean weighted displacement = "<<mean__weighted_displacement<<" mm."<<endl;
      
  // Initialise TRE (required to save info)
  reconstruction.InitTRE();
  
  //Save info
  if (debug)
    cout<<"Saving Info"<<endl;
  reconstruction.SlicesInfoCardiac4D( "info.tsv", stack_files );
  
  //Save transformations
  if(debug)
      cout<<"SaveTransformations"<<endl;
  reconstruction.SaveTransformations();
  
  //Save cine volume
  if(debug)
      cout<<"Saving cine volume."<<endl;
  irtkRealImage c = reconstruction.GetReconstructedCardiac4D();
  c.Write("cine_vol.nii.gz");
  
  if(debug)
      cout<<"SaveTransformations"<<endl;
  reconstruction.SaveTransformations();
  
  //Save mask
  if(debug)
      cout<<"Saving mask."<<endl;
  irtkRealImage m = reconstruction.GetMask();
  m.Write("mask.nii.gz");
  
  //Complete
  cout<<"simulateStacksCardiac complete."<<endl;

  //The end of main()
}  