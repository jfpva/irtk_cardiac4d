/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

 =========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkReconstructionCuda.h>
#include <irtkResampling.h>
#include <vector>
#include <string>
#include <cufft.h>
#include "perfstats.h"

using namespace std;

//2762_1_reconstruction_testGPU.nii 2 2762_1.nii 2762_2.nii id id -mask 2762_1_mask.nii -resolution 0.75
//2762_1_reconstruction_testGPU.nii 4 2762_1.nii 2762_2.nii 2762_3.nii 2762_4.nii id id id id -mask 2762_1_mask.nii -resolution 0.75
//2762_1_reconstruction_testGPU.nii 4 2762_1_d4.nii 2762_2_d4.nii 2762_3_d4.nii 2762_4_d4.nii id id id id -mask 2762_1_mask_d4.nii -resolution 0.75
//2762_1_reconstruction_testGPU.nii 8 2762_1.nii 2762_2.nii 2762_3.nii 2762_4.nii 2762_1.nii 2762_2.nii 2762_3.nii 2762_4.nii id id id id id id id id -mask 2762_1_mask.nii -resolution 0.75
//3115_nomasking.nii 8 3115_1.nii 3115_7.nii 3115_4.nii 3115_3.nii 3115_2.nii  3115_5.nii 3115_8.nii 3115_6.nii id id id id id id id id -mask manual_mask_3115_1.nii -log_prefix 3115_nomasking -smooth_mask 4 -resolution 0.75
//-thickness 1.26 1.26 1.26 1.26 1.26 1.26 1.26 1.26
//D:\Work\mysvn\FAUST\code\reconstruction\data\3115_test

//to be beaten full resolution
// inner loop
// full loop
//reconstruction loop:    5982.76 ms      (max = 36011 ms)
//overall:                659556 ms       (max = 659556 ms)
// overall time: 659.556000 s

//now
//reconstruction loop iter:7554.37 ms     (max = 71572 ms)
//overall:                831113 ms       (max = 831113 ms)
//overall time: 831.113000 s
// reconstruction loop iter:7013.08 ms     (max = 72476 ms)
// overall:                771544 ms       (max = 771544 ms)
// overall time: 771.544000 s


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
	cerr << "\t-mask [mask]              Binary mask to define the region od interest. [Default: whole image]"<<endl;
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
	cerr << "\t-log_prefix [prefix]      Prefix for the log file."<<endl;
	cerr << "\t-debug                    Debug mode - save intermediate results."<<endl;
	cerr << "\t-no_log                   Do not redirect cout and cerr to log files."<<endl;
	cerr << "\t" << endl;
	cerr << "\t" << endl;
	exit(1);
}

//TODO
/*
Interation :            663.333 ms      (max = 1663 ms

SliceToVolume :         112.556 ms      (max = 133 ms)
StackRegistrations:     84 ms   (max = 114 ms)
GaussianReconstruction :46.6667 ms      (max = 131 ms)
Superresolution :       21.3909 ms      (max = 25 ms)
CreateSlices :          20 ms   (max = 20 ms)
SimulateSlices :        16.4454 ms      (max = 25 ms)

CreateAverage:          4 ms    (max = 4 ms)
CropImage:              9 ms    (max = 9 ms)
EStep :                 1.45378 ms      (max = 2 ms)
InitializeEM :          8 ms    (max = 8 ms)
InitializeRS :          0.888889 ms     (max = 1 ms)
MaskSlices :            4 ms    (max = 4 ms)
MatchStack :            4 ms    (max = 4 ms)
recEvaluate :           5.11111 ms      (max = 6 ms)

*/

int main(int argc, char **argv)
{
	cudaDeviceReset();
	//2762_1_reconstruction_test.nii 4 2762_1_d4.nii 2762_2_d4.nii 2762_3_d4.nii 2762_4_d4.nii id id id id -mask 2762_1_mask_d4.nii -thickness 5.05 5.05 5.05 5.05 -resolution 4.00
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
	/// Stack thickness
	vector<double > thickness;
	///number of stacks
	int nStacks;
	/// number of packages for each stack
	vector<int> packages;


	// Default values.
	int templateNumber=-1;
	irtkRealImage *mask=NULL;
	int iterations = 7; //after 7 iterations mask breaks
	bool debug = false;
	double sigma=20;
	double resolution = 0.75;
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
	bool intensity_matching = true;

	irtkRealImage average;

	string log_id;
	bool no_log = false;


	//forced exclusion of slices
	int number_of_force_excluded_slices = 0;
	vector<int> force_excluded;

	//Create reconstruction object
	irtkReconstructionCuda reconstruction;

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
		stack.Read(argv[1]);
		stack_files.push_back(argv[1]);
		//stack.Write(argv[1]);
		reconstruction.Rescale(stack,1000);
		/*stack.GetImageToWorldMatrix().Print();
		std::cout << stack.GetOrigin() << std::endl;*/
		cout<<"Reading stack ... "<<argv[1]<<endl;
		
		/*irtkImageAttributes attr = stack.GetImageAttributes();
		for (int j = 0; j < attr._z; j++) {
			//create slice by selecting the appropreate region of the stack
			irtkRealImage slice = stack.GetRegion(0, 0, j, attr._x, attr._y, j + 1);
			slice.GetImageToWorldMatrix().Print();
		}*/

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

		argc--;
		argv++;
		irtkRigidTransformation *rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);
		stack_transformations.push_back(*rigidTransf);
		delete rigidTransf;
	}
#if 1
	reconstruction.InvertStackTransformations(stack_transformations);
#endif
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

		//Read binary mask for final volume
		if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
			argc--;
			argv++;
			mask= new irtkRealImage(argv[1]);
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

	//Create template volume with isotropic resolution 
	//if resolution==0 it will be determined from in-plane resolution of the image
	resolution = reconstruction.CreateTemplate(stacks[templateNumber],resolution);
	//->GPU

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
#if 1
	//volumetric registration
	reconstruction.StackRegistrations(stacks,stack_transformations,templateNumber);

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
		irtkRealImage m = reconstruction.GetMask();
		reconstruction.TransformMask(stacks[i],m,stack_transformations[i]);
		//Crop template stack
		reconstruction.CropImage(stacks[i],m);
		if (debug)
		{
			sprintf(buffer,"mask%i.nii",i);
			m.Write(buffer); 
			sprintf(buffer,"cropped%i.nii",i);
			stacks[i].Write(buffer);
		}
	}

	//Repeat volumetric registrations with cropped stacks
	//redirect output to files
	if ( ! no_log ) {
		cerr.rdbuf(file_e.rdbuf());
		cout.rdbuf (file.rdbuf());
	}
	//volumetric registration
	reconstruction.StackRegistrations(stacks,stack_transformations,templateNumber);
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
#endif
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

	PerfStats gstats;
	gstats.start();
	PerfStats istats;
	istats.start();
	double start = gstats.get_time();

	//ideally sync only slices and necessaries 

//////////////////////////////////////////////////////////////////////////////////////////////////
	//Crop and transform mask makes things wrong
	//need workaround for corrupted stack idxs -- guys, please use iterators
	std::vector<irtkRealImage> correctStacks;
	for(int i = 0; i < stacks.size(); i++)
	{
		//you crazy IRTK guys. look what I have to do because of your crazy iterator less memory alignment
		irtkRealImage tmp(stacks[i].GetX(), stacks[i].GetY(), stacks[i].GetZ());
		for(int z = 0; z < stacks[i].GetZ(); z++)
		{
			for(int y = 0; y < stacks[i].GetY(); y++)
			{
				for(int x = 0; x < stacks[i].GetX(); x++)
				{
					tmp(x,y,z) = stacks[i](x,y,z);
				}
			}
		}
		correctStacks.push_back(tmp);
	}

	//final sync
	//TODO when CoeffInit done, then only this for slices and transforms w/o CoeffInit sync
	reconstruction.SyncGPU();
	gstats.sample("SyncGPU()");

	for (int iter=0;iter<iterations;iter++)
	{
		//Print iteration number on the screen
		if ( ! no_log ) {
			cout.rdbuf (strm_buffer);
		}
		cout<<"Iteration "<<iter<<". "<<endl;
		//perform slice-to-volume registrations - skip the first iteration 
	#if 1
		if (iter>0)
		{
			if ( ! no_log ) {
				cerr.rdbuf(file_e.rdbuf());
				cout.rdbuf (file.rdbuf());
			}
			cout<<"Iteration "<<iter<<": "<<endl;
			//if((packages.size()>0)&&(iter<(iterations-1)))
			if((packages.size()>0)&&(iter<=iterations*(levels-1)/levels)&&(iter<(iterations-1)))
			{
				if(iter==1)
					reconstruction.PackageToVolume(stacks,packages);
				else
				{
					if(iter==2)
						reconstruction.PackageToVolume(stacks,packages,true);
					else
					{
						if(iter==3)
							reconstruction.PackageToVolume(stacks,packages,true,true);
						else
						{
							if(iter>=4)
								reconstruction.PackageToVolume(stacks,packages,true,true,iter-2);
							else
								reconstruction.SliceToVolumeRegistration();
						}
					}
				}
			}
			else
				reconstruction.SliceToVolumeRegistration();

			cout<<endl;
			if ( ! no_log ) {
				cerr.rdbuf (strm_buffer_e);
			}
		}
#endif
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

		//ideally sync only parameters used in CoeffINit from registration

		//Initialise values of weights, scales and bias fields
		reconstruction.InitializeEMValues();

		//Calculate matrix of transformation between voxels of slices and volume
#if !TEX_TEST
		reconstruction.CoeffInit();
		gstats.sample("CoeffInit inc sync");
#endif
		//Initialize reconstructed image with Gaussian weighted reconstruction
		//reconstruction.ResetSlices(stacks,thickness);
#if TEX_TEST
		reconstruction.UpdateGPUTranformationMatrices();
#endif
		reconstruction.GaussianReconstruction();
		gstats.sample("GaussianReconstruction()");
		reconstructed=reconstruction.GetReconstructedGPU();
		reconstructed.Write("GaussianReconstruction_test.nii");

		//Simulate slices (needs to be done after Gaussian reconstruction)
		//reconstruction.SimulateSlices(true);
		reconstruction.SimulateSlices();

		//Initialize robust statistics parameters
		reconstruction.InitializeRobustStatistics();
		gstats.sample("InitializeRobustStatistics()");

		//std::cin.get();
		//EStep
		reconstruction.EStep();

		//number of reconstruction iterations
		if ( iter==(iterations-1) ) 
		{
			rec_iterations = 30;      
		}
		else 
			rec_iterations = 10;

		//reconstruction iterations


//////////////////////////////////////////////////////////////////////////////////////////////////
		i=0;
		for (i=0;i<rec_iterations;i++)
		{
			cout<<endl<<"  Reconstruction iteration "<<i<<". "<<endl;

			if (intensity_matching)
			{
				//calculate bias fields

				// currently cudarizing
				if (sigma>0)
				{
					reconstruction.Bias();
					istats.sample("Bias()");
				}
				//calculate scales
				reconstruction.Scale(); //TODO gpu to save sync
				istats.sample("Scale()");
			}
			//temporary sync
			//reconstruction.SyncGPU();

			//MStep and update reconstructed volume


			// double start = gstats.get_time();
			reconstruction.Superresolution(i+1);
			istats.sample("Superresolution()");

			// double end = gstats.get_time();
			// printf("gpu+cpu time superres: %f s \n", end-start);

			if (intensity_matching)
			{
				if((sigma>0)&&(!global_bias_correction))
				{
					reconstruction.NormaliseBias(i);
					istats.sample("NormaliseBias()");
				}
			}

			// Simulate slices (needs to be done
			// after the update of the reconstructed volume)
			reconstruction.SimulateSlices();
			istats.sample("SimulateSlices()");

			reconstruction.MStep(i+1);
			istats.sample("MStep()");
			//E-step
			reconstruction.EStep();
			istats.sample("EStep()");
			//Save intermediate reconstructed image
			if (debug)
			{
				//reconstructed=reconstruction.GetReconstructed();
				reconstructed=reconstruction.GetReconstructedGPU();
				sprintf(buffer,"super%i.nii",i);
				reconstructed.Write(buffer);
			}
		gstats.sample("reconstruction loop iter");


		}//end of reconstruction iterations
		//gstats.sample("reconstruction loop");
//////////////////////////////////////////////////////////////////////////////////////////////////
		//temporary sync
		reconstruction.SyncCPU();

		//Mask reconstructed image to ROI given by the mask
		reconstruction.MaskVolume();

		//Save reconstructed image
		if (true)//debug)
		{
			reconstructed=reconstruction.GetReconstructed();
			sprintf(buffer,"image%i.nii",iter);
			reconstructed.Write(buffer);
			//reconstruction.SaveConfidenceMap();
		}

		//Evaluate - write number of included/excluded/outside/zero slices in each iteration in the file
		if ( ! no_log ) {
			cout.rdbuf (fileEv.rdbuf());
		}
		reconstruction.Evaluate(iter);
		cout<<endl;

		if ( ! no_log ) {
			cout.rdbuf (strm_buffer);
		}
		// stats.sample("interleaved registration-reconstruction iterations");
	}// end of interleaved registration-reconstruction iterations
//////////////////////////////////////////////////////////////////////////////////////////////////

	//save final result
	reconstruction.RestoreSliceIntensities();
	reconstruction.ScaleVolume();
	reconstructed=reconstruction.GetReconstructed();
	reconstructed.Write(output_name); 

	reconstruction.SlicesInfo( "SlicesInfo.tsv", stack_files );

	if(debug)
	{
		reconstruction.SaveTransformations();
		reconstruction.SaveSlices();

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
	double end = gstats.get_time();
	istats.sample("overall");
	gstats.print();
	istats.print();
	
	printf("\n\noverall time: %f s\n\n", end-start);

	cudaDeviceReset();

	//The end of main()
}  
