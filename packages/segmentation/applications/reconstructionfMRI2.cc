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
#include <irtkReconstructionfMRI.h>
#include <irtkLaplacianSmoothing.h>
#include <vector>
#include <string>
using namespace std;

//Application to perform reconstruction of volumetric MRI from thick slices.

void usage()
{
  cerr << "Usage: reconstruction [reconstructed] [T2] [N] [stack_1] .. [stack_N] [dof_1] .. [dof_N] <options>\n" << endl;
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
  cerr << "\t-mask [mask]              Binary mask to define the region od interest. [Default: whole image]"<<endl;
  cerr << "\t-packages [num_1] .. [num_N] Give number of packages used during acquisition for each stack."<<endl;
  cerr << "\t                          The stacks will be split into packages during registration iteration 1"<<endl;
  cerr << "\t                          and then into odd and even slices within each package during "<<endl;
  cerr << "\t                          registration iteration 2. The method will then continue with slice to"<<endl;
  cerr << "\t                          volume approach. [Default: slice to volume registration only]"<<endl;
  cerr << "\t-iterations [iter]        Number of registration-reconstruction iterations. [Default: 9]"<<endl;
  cerr << "\t-sigma [sigma]            Stdev for bias field. [Default: 12mm]"<<endl;
  cerr << "\t-resolution [res]         Isotropic resolution of the volume. [Default: 0.75mm]"<<endl;
  cerr << "\t-levels [levels]          Multiresolution levels for B-spline fieldmap. [Default: 3]"<<endl;
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
  cerr << "\t-group [group1] ... [groupN]  Give group numbers for distortion correction. [Default: all in one group]"<<endl;
  cerr << "\t-phase [axis1] ... [axisG]  For each group give phase encoding axis in image coordinates [x or y]"<<endl;
  cerr << "\t-fieldMapSpacing [spacing]  B-spline control point spacing for field map"<<endl;
  cerr << "\t-smoothnessPenalty [penalty]  smoothness penalty for field map"<<endl;
  cerr << "\t-info [filename]          Filename for slice information in\
                                       tab-sparated columns."<<endl;
  cerr << "\t-log_prefix [prefix]      Prefix for the log file."<<endl;
  cerr << "\t-debug                    Debug mode - save intermediate results."<<endl;
  cerr << "\t" << endl;
  cerr << "\t" << endl;
  exit(1);
}

irtkRealImage MakeSequence(vector<irtkRealImage> stacks)
{
  irtkImageAttributes attr = stacks[0].GetImageAttributes();
  attr._t = stacks.size();
  irtkRealImage image(attr);
  for(int t=0;t<attr._t;t++)
    for(int z=0;z<attr._z;z++)
      for(int y=0;y<attr._y;y++)
        for(int x=0;x<attr._x;x++)
	{
	  image(x,y,z,t)=stacks[t](x,y,z);
	}
  return image;
}

irtkRealImage CreateMask(irtkRealImage image)
{
  irtkImageAttributes attr = image.GetImageAttributes();
  int _t = attr._t;
  attr._t = 1;
  irtkRealImage mask(attr);
  mask=0;
  for(int t=0;t<_t;t++)
    for(int z=0;z<attr._z;z++)
      for(int y=0;y<attr._y;y++)
        for(int x=0;x<attr._x;x++)
	{
	  if(image(x,y,z,t)>0)
	    mask(x,y,z)=1;
	}
  return mask;  
}

int main(int argc, char **argv)
{
  
  //utility variables
  int i,j,k,t,index, iter, ok;
  char buffer[256];
  irtkRealImage stack; 
  
  //declare variables for input
  /// Name for output volume
  char * output_name = NULL;
  /// Slice stacks
  vector<irtkRealImage> stacks, corrected_stacks, simulated;
  vector<string> stack_files;
  /// Stack transformation
  vector<irtkRigidTransformation> stack_transformations;
  /// Stack thickness
  vector<double > thickness;
  ///number of stacks
  int nStacks;
  /// number of packages for each stack
  vector<int> packages;
  char * order = NULL;
  int step = 1;
  int rewinder = 1;
  ///T2 template
  irtkRealImage t2;

  // Default values.
  int templateNumber=-1;
  irtkRealImage *mask=NULL;
  irtkRealImage b0_mask,T2_mask;
  int iterations = 9;
  bool debug = false;
  double sigma=20;
  double resolution = 2.5;
  double lambda = 0.02;
  double delta = 150;
  int levels = 3;
  double lastIterLambda = 0.01;
  int rec_iterations;
  double averageValue = 700;
  double smooth_mask = 4;
  bool global_bias_correction = false;
  bool have_stack_transformations = false;
  double low_intensity_cutoff = 0.01;
  //folder for slice-to-volume registrations, if given
  char * folder=NULL;
  //flag to remove black background, e.g. when neonatal motion correction is performed
  bool remove_black_background = false;
  //flag to swich the intensity matching on and off
  bool intensity_matching = true;
  bool alignT2 = false;
  double fieldMapSpacing = 10;
  double penalty = 0.1;
  //flag to swich the robust statistics on and off
  bool robust_statistics = true;
  bool smooth = false;
  int multiband_factor = 1;

  //Create reconstruction object
  irtkReconstructionfMRI reconstruction;

  irtkRealImage average;

  string log_id;
  string info_filename = "slice_info.tsv";

  
  //distortion correction
  vector<int> groups, stack_group;
  vector<bool> swap;

  
  //forced exclusion of slices
  int number_of_force_excluded_slices = 0;
  vector<int> force_excluded;
  //if not enough arguments print help
  if (argc < 5)
    usage();
  
  //read output name
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"Recontructed volume name ... "<<output_name<<endl;
  cout.flush();

  //read output name
  t2.Read(argv[1]);
  cout<<"T2 template ... "<<output_name<<endl;
  cout.flush();
  argc--;
  argv++;

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
		
		char temp;
		order = argv[1];
		temp = *order;
		
		if ((temp == 'A')) {
				cout<<"Slice acquisition order is ascending"<<endl;
		}
		else if ((temp == 'D')) {
				cout<<"Slice acquisition order is descending"<<endl;
		}
		else if ((temp == 'F')) {
				cout<<"Slice acquisition order is default"<<endl;
		}
		else if ((temp == 'I')) {
				cout<<"Slice acquisition order is interleaved"<<endl;
		}
		else if ((temp == 'C'))	 {
				cout<<"Slice acquisition order is customized"<<endl;
		}
		else	{
				cout<<"Slice acquisition order not recognised, set to ascending (A)"<<endl;
		}
	
	   ok = true;
	   argc--;
	   argv++;
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
    
    //Multiband factor
	if ((ok == false) && (strcmp(argv[1], "-multiband") == 0)){
	  argc--;
	  argv++;
	  multiband_factor=atof(argv[1]);
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
    if ((ok == false) && (strcmp(argv[1], "-levels") == 0)){
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
       }
       argc--;
       argv++;
       cout<<"."<<endl;
       cout.flush();
       
       bool group_exists;
      for (i=0;i<stack_group.size();i++)
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
      for (j=0;j<groups.size();j++) cout<<" "<<groups[j];
      cout<<endl;
    
      ok = true;
    }

    //Read phase encoding direction
    if ((ok == false) && (strcmp(argv[1], "-phase") == 0)){
      argc--;
      argv++;
      
      for (i=0;i<nStacks;i++)
	  {
    	  stack_group.push_back(1);
    	  cout<<stack_group[i]<<" ";
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
  
  if (have_stack_transformations == true)
	  reconstruction.InvertStackTransformations(stack_transformations);
  else {
	  for (i=0;i<nStacks;i++)
	  {
	    irtkRigidTransformation *rigidTransf = new irtkRigidTransformation;
	    stack_transformations.push_back(*rigidTransf);
	    delete rigidTransf;
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

  //Set debug mode
  if (debug) reconstruction.DebugOn();
  else reconstruction.DebugOff();
  
  //Set B-spline control point spacing for field map
  reconstruction.SetFieldMapSpacing(fieldMapSpacing/2);
  
  //Set B-spline control point spacing for field map
  reconstruction.SetSmoothnessPenalty(penalty);
  
  //Set force excluded slices
  reconstruction.SetForceExcludedSlices(force_excluded);

  //Set low intensity cutoff for bias estimation
  reconstruction.SetLowIntensityCutoff(low_intensity_cutoff)  ;

  //Create template volume with isotropic resolution 
  //if resolution==0 it will be determined from in-plane resolution of the image
  resolution = reconstruction.CreateTemplate(t2,0);
  cout<<"Resolution is "<<resolution<<endl;
  
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

  //perform volumetric registration of the stacks
  //redirect output to files
  cerr.rdbuf(file_e.rdbuf());
  cout.rdbuf (file.rdbuf());
  //volumetric registration
  reconstruction.SetT2Template(t2);//SetReconstructed(t2);
  reconstruction.StackRegistrations(stacks,stack_transformations);
  
  cout<<endl;
  cout.flush();
  //redirect output back to screen
  cout.rdbuf (strm_buffer);
  cerr.rdbuf (strm_buffer_e);
  
  //Volumetric registrations are stack-to-template while slice-to-volume 
  //registrations are actually performed as volume-to-slice (reasons: technicalities of implementation)
  //Need to invert stack transformations now.  
  //reconstruction.InvertStackTransformations(stack_transformations);
  average = reconstruction.CreateAverage(stacks,stack_transformations);
  if (debug)
    average.Write("average1.nii.gz");

  //Rescale intensities of the stacks to have the same average
  //if (intensity_matching)
  //  reconstruction.MatchStackIntensities(stacks,stack_transformations,averageValue);
  //else
  //  reconstruction.MatchStackIntensities(stacks,stack_transformations,averageValue,true);
  average = reconstruction.CreateAverage(stacks,stack_transformations);
  if (debug)
    average.Write("average2.nii.gz");
  //reconstruction.SetReconstructed(average);
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

  //prepare larger mask for fieldmap
  reconstruction.CreateLargerMask(*mask);
  

  //interleaved registration-reconstruction iterations
 // int internal = reconstruction.giveMeDepth(stacks, packages, multiband_factor);

  if (iterations == 0)	{
	  if (multiband_factor>1)
		  iterations = internal*2;
	  else
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
  
  double step2 = 2.5;
  //if(iterations>1)
    //iterations = 7;
  //registration iterations
  for (iter=0;iter<(iterations);iter++)
  {
    //Print iteration number on the screen
    cout.rdbuf (strm_buffer);
    cout<<"Iteration "<<iter<<". "<<endl;
    cout.flush();
    
    //Step one: simulate stacks
     reconstruction.SetT2Template(t2);
    reconstruction.CoeffInit();
    cout<<"Going to simulate stacks ...";
    cout.flush();
    simulated.clear();
    for(i=0;i<stacks.size();i++)
      simulated.push_back(stacks[i]);
    //reconstruction.MaskSlices();
    reconstruction.SimulateStacks(simulated);
    //reconstruction.UpdateSlices(stacks,thickness);
    //reconstruction.SaveSlices();
    cout<<"done."<<endl;
  
    for (i=0;i<stacks.size();i++)
    {
      sprintf(buffer,"simulated%i.nii.gz",i);
      simulated[i].Write(buffer);
    }
    
    //Step 2: create 4D images from stacks and simulated
    
    irtkRealImage original = MakeSequence(stacks);
    irtkRealImage simul = MakeSequence(simulated);
    original.Write("original.nii.gz");
    simul.Write("simulated.nii.gz");
    
    //Step 3: calculate distortion
    
    irtkRealImage fieldmap = stacks[0]; //reconstruction.GetReconstructed();
    irtkMultiLevelFreeFormTransformation _fieldMap;
    cout<<"entering fieldmap distortion"<<endl;
    cout.flush();
    //cout<<levels<<endl;
    //reconstruction.SetSmoothnessPenalty(penalty);
    //reconstruction.SetFieldMapSpacing(spacing);
    reconstruction.FieldMapDistortion(original,simul,_fieldMap,swap[0],step2,0, levels);
  
    cout<<"done with fieldmap distortion"<<endl;
    cout.flush();

    _fieldMap.irtkTransformation::Write("distortion_transformation.dof.gz");
    
    double x,y,z,xx,yy,zz;
    irtkImageAttributes attr = original.GetImageAttributes();
    //irtkRealImage jacobian(fieldmap);
    //jacobian=0;
  
    for(k=0; k<fieldmap.GetZ();k++)
    for(j=0; j<fieldmap.GetY();j++)
    for(i=0; i<fieldmap.GetX();i++)
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
      //jacobian(i,j,k)=_fieldMap.irtkTransformation::Jacobian(x,y,z);
      _fieldMap.Transform(x,y,z);
      original.WorldToImage(x,y,z);
    
      if (swap[0])
        fieldmap(i,j,k)=(y-yy)*attr._dy;
      else
        fieldmap(i,j,k)=(x-xx)*attr._dx;
     }
   
    sprintf(buffer,"fieldmap-bspline-%i.nii.gz",iter);
    fieldmap.Write(buffer);
    //sprintf(buffer,"jacobian-fieldmap-bspline-%i.nii.gz",iter);
    //jacobian.Write(buffer);

  
    //Step 4: Smooth fieldmap
    double padding;
    if(smooth)
    {   
      irtkRealImage fieldmap_mask = CreateMask(simul) ;
      fieldmap_mask.Write("fieldmap_mask.nii.gz");
      fieldmap_mask = reconstruction.CreateLargerMask(fieldmap_mask);
      fieldmap_mask.Write("fieldmap_mask_larger.nii.gz");
      irtkLaplacianSmoothing smoothing;
      smoothing.SetInput(fieldmap);
      smoothing.SetMask(fieldmap_mask);
      irtkRealImage fieldmap_smooth = smoothing.RunGD();
      sprintf(buffer,"fieldmap-smooth-%i.nii.gz",iter);
      fieldmap_smooth.Write(buffer);
      fieldmap = fieldmap_smooth;
      fieldmap_mask = smoothing.GetMask();
      fieldmap_mask.Write("fieldmap_mask_new.nii.gz");
   
    
      //introduce padding in fieldmap
      irtkRealPixel mn,mx;
      fieldmap.GetMinMax(&mn,&mx);
      padding = round(mn-2);
      cout<<"padding = "<<padding<<endl;
    
      irtkRealPixel *pf = fieldmap.GetPointerToVoxels();
      irtkRealPixel *pm = fieldmap_mask.GetPointerToVoxels();
      for (j=0; j<fieldmap.GetNumberOfVoxels();j++)
      {
        if(*pm!=1)
          *pf = padding;
        pm++;
        pf++;
      }

      fieldmap.Write("fieldmap-pad.nii.gz");
    }
    else
    {
      irtkRealPixel mn,mx;
      fieldmap.GetMinMax(&mn,&mx);
      padding = round(mn-2);
      cout<<"padding = "<<padding<<endl;
    }
    
    
    //Step 5: Correct stacks by smooth fieldmap

    //resample fieldmap 
    irtkRealImage resfieldmap = fieldmap; //not needed - both on the same grid
    
    //correct stacks
    corrected_stacks.clear();
    for (index = 0; index<stacks.size();index++)
    {
      irtkRealImage image = stacks[index];
      irtkRealImage output = stacks[index];
      output=0;
      
      irtkImageFunction *interpolator;
      interpolator = new irtkLinearInterpolateImageFunction;
  
      interpolator->SetInput(&image);
      interpolator->Initialize();
  
      attr = image.GetImageAttributes();
      for (t = 0; t < image.GetT(); t++) {
        for (k = 0; k < image.GetZ(); k++) {
          for (j = 0; j < image.GetY(); j++) {
            for (i = 0; i < image.GetX(); i++) {
	      if (resfieldmap(i,j,k)>padding)
	      {
                x = i;
                y = j;
                z = k;
	        //move it by fieldmap converted to voxels (reconstructed reslution)
	        if(swap[0])
	        {
	          y+=resfieldmap(i,j,k)/attr._dy;
	        }
	        else
	        {
	          x+=resfieldmap(i,j,k)/attr._dx;
	        }
	  
	        if ((x > -0.5) && (x < image.GetX()-0.5) && 
	            (y > -0.5) && (y < image.GetY()-0.5) &&
                    (z > -0.5) && (z < image.GetZ()-0.5))
	        {
	          output(i,j,k,t) = interpolator->Evaluate(x,y,z,t); // * jacobian(i,j,k);
	        }
	        else
		  output(i,j,k,t) = 0;
	      }
	      else
		output(i,j,k,t) = 0;
            }
          }
        }
      }

      sprintf(buffer,"cor%i-%i.nii.gz",iter,index);
      output.Write(buffer);
      corrected_stacks.push_back(output);

    }

    reconstruction.UpdateSlices(corrected_stacks,thickness);
    reconstruction.MaskSlices();
    
   
    if (iter > 0) {
	  //reconstruction.InterpolateBSpline(stacks,iter);
	  //reconstruction.InterpolateBSplineReordered(stacks,multiband_factor,iter);
    	reconstruction.SetT2Template(corrected_stacks[0]);//SetReconstructed(t2);
	  reconstruction.InterpolateGaussianReordered(stacks,multiband_factor,iter);
	  reconstruction.SaveRegistrationStep(stacks,iter);
	  
	  //reconstruction.InterpolateGaussianReordered(stacks,multiband_factor,iter);
	}
  	
    if (iter>0)
	{
		/*if ( ! no_log ) {
		  cerr.rdbuf(file_e.rdbuf());
		  cout.rdbuf (file.rdbuf());
		}*/
    	
    	reconstruction.SetT2Template(t2);//SetReconstructed(t2);
	
    	vector<int> level;
    				if(iter == 1) {
    					reconstruction.newPackageToVolume(stacks, packages, multiband_factor, *order, step, rewinder,iter);
    				}
    		
    				else if((iter > 1) && (iter < internal-1)){
    					//level = reconstruction.giveMeSplittingVector(stacks, packages, multiband_factor, iter, false);
    					reconstruction.ChunkToVolume(stacks, packages, level, multiband_factor, *order, step, rewinder,iter);
    				}
    		
    				else {
    					
    					//level = reconstruction.giveMeSplittingVector(stacks, packages, multiband_factor, iter, true);
    					if (multiband_factor == 1) {
    						reconstruction.ChunkToVolume(stacks, packages, level, 1, *order, step, rewinder,iter);
    					}
    					else {
    						reconstruction.ChunkToVolume(stacks, packages, level, multiband_factor, *order, step, rewinder,iter);
    					}
    					
    				}
		
		if (debug)
		{
			reconstructed=reconstruction.GetReconstructed();
			sprintf(buffer,"image%i.nii.gz",iter);
			reconstructed.Write(buffer);
		}
	
		/*if ( ! no_log ) {
			cerr.rdbuf (strm_buffer_e);
		}*/
	}
    
    
    //perform slice-to-volume registrations - skip the first iteration 
    /*if (iter>=0)
    {
      cerr.rdbuf(file_e.rdbuf());
      cout.rdbuf (file.rdbuf());
      cout<<"Iteration "<<iter<<": "<<endl;
      
      reconstruction.SetT2Template(t2);//SetReconstructed(t2);
      
      if((packages.size()>0)&&(iter<2)&&(iter<(iterations-1)))
      {
	if(iter==0)
          reconstruction.PackageToVolume(corrected_stacks,packages,iter);
	else
	{
	  if(iter==1)
            reconstruction.PackageToVolume(corrected_stacks,packages,iter,true);
	  else
	  {
            if(iter==2)
	      reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true);
	    else
	    {
	      if(iter>=3)
                reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true,iter-1);
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
    }*/
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
    //reconstruction.PutMask(b0_mask);
    reconstruction.CoeffInit();
    
    //Initialize reconstructed image with Gaussian weighted reconstruction
    reconstruction.GaussianReconstruction();
    
    //Simulate slices (needs to be done after Gaussian reconstruction)
    reconstruction.SimulateSlices();
    
    //Initialize robust statistics parameters
    reconstruction.InitializeRobustStatistics();
    
    //EStep
    if(robust_statistics)
      reconstruction.EStep();

    //number of reconstruction iterations
    rec_iterations = 10;      
        
    // super resolution not needed for fMRI
  	rec_iterations = 0;
  	  
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
    //reconstruction.PutMask(T2_mask);
    reconstruction.MaskVolume();


   //Evaluate - write number of included/excluded/outside/zero slices in each iteration in the file
   cout.rdbuf (fileEv.rdbuf());
   reconstruction.Evaluate(iter);
   cout<<endl;
   cout.flush();
   cout.rdbuf (strm_buffer); 
   
  //}// end of interleaved registration-reconstruction iterations

  //save final result
  //reconstruction.SaveDistortionTransformations();
  //reconstruction.RestoreSliceIntensities();
  //reconstruction.ScaleVolume();
  reconstructed=reconstruction.GetReconstructed();
  reconstructed.Write(output_name); 
  reconstruction.SaveTransformations();
  reconstruction.SaveSlices();
  
  if ( info_filename.length() > 0 )
    reconstruction.SlicesInfo( info_filename.c_str(),
                               stack_files );

  
  //The end of main()
}  
