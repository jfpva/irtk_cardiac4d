/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#include <irtkReconstructionCardiac4D.h>
#include <irtkResampling.h>
#include <irtkRegistration.h>
#include <irtkImageRigidRegistration.h>
#include <irtkImageRigidRegistrationWithPadding.h>
#include <irtkTransformation.h>
#include <math.h>


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
irtkReconstructionCardiac4D::irtkReconstructionCardiac4D():irtkReconstruction()
{
    _recon_type = _3D;
}

// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
irtkReconstructionCardiac4D::~irtkReconstructionCardiac4D() { }

// -----------------------------------------------------------------------------
// Set Slice R-R Intervals
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SetSliceRRInterval( vector<double> rr )
{
    _slice_rr = rr;
}

void irtkReconstructionCardiac4D::SetSliceRRInterval( double rr )
{
    vector<double> slice_rr;
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) 
    {
        slice_rr.push_back(rr);
    }
    _slice_rr = slice_rr;
}


// -----------------------------------------------------------------------------
// Set Slice R-R Intervals
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SetLocRRInterval( vector<double> rr )
{
    vector<double> slice_rr;
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) 
    {
        slice_rr.push_back(rr[_loc_index[inputIndex]]);
    }
    _slice_rr = slice_rr;
}


// -----------------------------------------------------------------------------
// Set Slice Cardiac Phases
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SetSliceCardiacPhase( vector<double> cardiacphases )
{
    _slice_cardphase = cardiacphases;
}

void irtkReconstructionCardiac4D::SetSliceCardiacPhase()
{
    _slice_cardphase.clear();
    for (unsigned int i=0; i<_slices.size(); i++)
      _slice_cardphase.push_back( 0 );
}

// -----------------------------------------------------------------------------
// Determine Reconstructed Spatial Resolution
// -----------------------------------------------------------------------------
//determine resolution of volume to reconstruct
double irtkReconstructionCardiac4D::GetReconstructedResolutionFromTemplateStack( irtkRealImage stack )
{
    double dx, dy, dz, d;
    stack.GetPixelSize(&dx, &dy, &dz);
    if ((dx <= dy) && (dx <= dz))
        d = dx;
    else if (dy <= dz)
        d = dy;
    else
        d = dz;
    return d;
}


// -----------------------------------------------------------------------------
// Get Slice-Location Transformations
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::ReadSliceTransformation(char* slice_transformations_folder)
{
    if (_slices.size()==0) 
    {
      cerr << "Please create slices before reading transformations!" << endl;
      exit(1);
    }
    
    int nLoc = _loc_index.back() + 1;

    char name[256];
    char path[256];
    vector<irtkRigidTransformation> loc_transformations;
    irtkTransformation *transformation;
    irtkRigidTransformation *rigidTransf;

    // Read transformations from file
    cout << "Reading transformations:" << endl;
    for (int iLoc = 0; iLoc < nLoc; iLoc++) {
        if (slice_transformations_folder != NULL) {
            sprintf(name, "/transformation%05i.dof", iLoc);
            strcpy(path, slice_transformations_folder);
            strcat(path, name);
        }
        else {
            sprintf(path, "transformation%03i.dof", iLoc);
        }
        transformation = irtkTransformation::New(path);
        rigidTransf = dynamic_cast<irtkRigidTransformation*>(transformation);
        loc_transformations.push_back(*rigidTransf);
        delete transformation;
        cout << path << endl;
    }
    
    // Assign transformations to single-frame images
    _transformations.clear();
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) 
    { 
        _transformations.push_back(loc_transformations[_loc_index[inputIndex]]);
    }
    cout << "ReadSliceTransformations complete." << endl;
    
}


// -----------------------------------------------------------------------------
// Set Reconstructed Cardiac Phases
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SetReconstructedCardiacPhase( vector<double> cardiacphases )
{
    _reconstructed_cardiac_phases = cardiacphases;
}


// -----------------------------------------------------------------------------
// Set Reconstructed R-R Interval
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SetReconstructedRRInterval( double rrinterval )
{
    _reconstructed_rr_interval = rrinterval;
    if (_debug)
      cout<<"Reconstructed R-R interval = "<<_reconstructed_rr_interval<<" s."<<endl;
}


// -----------------------------------------------------------------------------
// Set Reconstructed Cardiac Phases
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SetReconstructedTemporalResolution( double temporalresolution )
{
    _reconstructed_temporal_resolution = temporalresolution;
    if (_debug)
      cout<<"Reconstructed temporal resolution = "<<_reconstructed_temporal_resolution<<" s."<<endl;
}


// -----------------------------------------------------------------------------
// Initialise Reconstructed Volume from Static Mask
// -----------------------------------------------------------------------------
// Create zero image as a template for reconstructed volume
void irtkReconstructionCardiac4D::CreateTemplateCardiac4DFromStaticMask( irtkRealImage mask, double resolution )
{
    // Get mask attributes - image size and voxel size
    irtkImageAttributes attr = mask.GetImageAttributes();
    
    // Set temporal dimension
    attr._t = _reconstructed_cardiac_phases.size();
    attr._dt = _reconstructed_temporal_resolution;

    // Create volume
    irtkRealImage volume4D(attr);
    
    // Initialise volume values
    volume4D = 0;
      
    // Resample to specified resolution
    irtkNearestNeighborInterpolateImageFunction interpolator;
    irtkResampling<irtkRealPixel> resampling(resolution,resolution,resolution);
    resampling.SetInput(&volume4D);
    resampling.SetOutput(&volume4D);
    resampling.SetInterpolator(&interpolator);
    resampling.Run();

    // Set recontructed 4D volume
    _reconstructed4D = volume4D;
    _template_created = true;
    
    // Set reconstructed 3D volume for reference by existing functions of irtkReconstruction class
    irtkImageAttributes attr4d = volume4D.GetImageAttributes();
    irtkImageAttributes attr3d = attr4d;
    attr3d._t = 1;
    irtkRealImage volume3D(attr3d);
    _reconstructed = volume3D;
    
    // Debug
    if (_debug)
    {
      cout << "CreateTemplateCardiac4DFromStaticMask: created template 4D volume with "<<attr4d._t<<" time points and "<<attr4d._dt<<" s temporal resolution."<<endl;
    }
}


// -----------------------------------------------------------------------------
// Match Stack Intensities With Masking
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::MatchStackIntensitiesWithMasking(vector<irtkRealImage>& stacks,
                                               vector<irtkRigidTransformation>& stack_transformations, double averageValue, bool together)
{
    if (_debug)
        cout << "Matching intensities of stacks. ";

    //Calculate the averages of intensities for all stacks
    double sum, num;
    char buffer[256];
    unsigned int ind;
    int i, j, k;
    double x, y, z;
    vector<double> stack_average;
    irtkRealImage m;
        
    //remember the set average value
    _average_value = averageValue;
    
    //averages need to be calculated only in ROI
    for (ind = 0; ind < stacks.size(); ind++) {
        irtkImageAttributes attr = stacks[ind].GetImageAttributes();
        attr._t = 1;
        m.Initialize(attr);
        m = 0;
        sum = 0;
        num = 0;
        for (i = 0; i < stacks[ind].GetX(); i++)
            for (j = 0; j < stacks[ind].GetY(); j++)
                for (k = 0; k < stacks[ind].GetZ(); k++) {
                    //image coordinates of the stack voxel
                    x = i;
                    y = j;
                    z = k;
                    //change to world coordinates
                    stacks[ind].ImageToWorld(x, y, z);
                    //transform to template (and also _mask) space
                    stack_transformations[ind].Transform(x, y, z);
                    //change to mask image coordinates - mask is aligned with template
                    _mask.WorldToImage(x, y, z);
                    x = round(x);
                    y = round(y);
                    z = round(z);
                    //if the voxel is inside mask ROI include it
                     if ((x >= 0) && (x < _mask.GetX()) && (y >= 0) && (y < _mask.GetY()) && (z >= 0)
                          && (z < _mask.GetZ()))
                          {
                            if (_mask(x, y, z) == 1)
                            {
      			                  m(i,j,k)=1;
                              for ( int f = 0; f < stacks[ind].GetT(); f++) 
                              {
                                sum += stacks[ind](i, j, k, f);
                                num++;
                              }
                            }
                    }
                }
         if(_debug)
	 {
           sprintf(buffer,"mask-for-matching%i.nii.gz",ind);   
	   m.Write(buffer);
	 }
        //calculate average for the stack
        if (num > 0)
            stack_average.push_back(sum / num);
        else {
            cerr << "Stack " << ind << " has no overlap with ROI" << endl;
            exit(1);
        }
    }
    
    double global_average;
    if (together) {
        global_average = 0;
        for(i=0;i<stack_average.size();i++)
            global_average += stack_average[i];
        global_average/=stack_average.size();
    }

    if (_debug) {
        cout << "Stack average intensities are ";
        for (ind = 0; ind < stack_average.size(); ind++)
            cout << stack_average[ind] << " ";
        cout << endl;
        cout << "The new average value is " << averageValue << endl;
    }

    //Rescale stacks
    irtkRealPixel *ptr;
    double factor;
    for (ind = 0; ind < stacks.size(); ind++) {
        if (together) {
            factor = averageValue / global_average;
            _stack_factor.push_back(factor);
        }
        else {
            factor = averageValue / stack_average[ind];
            _stack_factor.push_back(factor);

        }

        ptr = stacks[ind].GetPointerToVoxels();
        for (i = 0; i < stacks[ind].GetNumberOfVoxels(); i++) {
            if (*ptr > 0)
                *ptr *= factor;
            ptr++;
        }
    }

    if (_debug) {
        for (ind = 0; ind < stacks.size(); ind++) {
            sprintf(buffer, "rescaled-stack%i.nii.gz", ind);
            stacks[ind].Write(buffer);
        }

        cout << "Slice intensity factors are ";
        for (ind = 0; ind < stack_average.size(); ind++)
            cout << _stack_factor[ind] << " ";
        cout << endl;
        cout << "The new average value is " << averageValue << endl;
    }

}


// -----------------------------------------------------------------------------
// Create Slices and Associated Transformations
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::CreateSlicesAndTransformationsCardiac4D( vector<irtkRealImage> &stacks,
                                                         vector<irtkRigidTransformation> &stack_transformations,
                                                         vector<double> &thickness,
                                                         const vector<irtkRealImage> &probability_maps )
{

    double sliceAcqTime;
    int loc_index = 0;

    if (_debug)
        cout << "CreateSlicesAndTransformations" << endl;
    
    //for each stack
    for (unsigned int i = 0; i < stacks.size(); i++) {
        //image attributes contain image and voxel size
        irtkImageAttributes attr = stacks[i].GetImageAttributes();

        //attr._z is number of slices in the stack
        for (int j = 0; j < attr._z; j++) {
          
            //attr._t is number of frames in the stack
            for (int k = 0; k < attr._t; k++) {
              
                //create slice by selecting the appropreate region of the stack
                irtkRealImage slice = stacks[i].GetRegion(0, 0, j, k, attr._x, attr._y, j + 1, k + 1);
                //set correct voxel size in the stack. Z size is equal to slice thickness.
                slice.PutPixelSize(attr._dx, attr._dy, thickness[i], attr._dt);
                //set slice acquisition time
                sliceAcqTime = attr._torigin + k * attr._dt; // TODO: check calculation of _slice_time from input stack
                _slice_time.push_back(sliceAcqTime);  
                //set slice temporal resolution
                _slice_dt.push_back(attr._dt);
                //remember the slice
                _slices.push_back(slice);
                _simulated_slices.push_back(slice);
                _simulated_weights.push_back(slice);
                _simulated_inside.push_back(slice);
                //remeber stack indices for this slice
                _stack_index.push_back(i);
                _loc_index.push_back(loc_index);
                _stack_loc_index.push_back(j);
                _stack_dyn_index.push_back(k);
                //initialize slice transformation with the stack transformation
                _transformations.push_back(stack_transformations[i]);
                if ( probability_maps.size() > 0 ) {
                    irtkRealImage proba = probability_maps[i].GetRegion(0, 0, j, k, attr._x, attr._y, j + 1, k + 1);
                    proba.PutPixelSize(attr._dx, attr._dy, thickness[i], attr._dt);
                    _probability_maps.push_back(proba);
                }
            }
            loc_index++;
        }
    }
    cout << "Number of images: " << _slices.size() << endl;
}


// -----------------------------------------------------------------------------
// TODO: Calculate Cardiac Phase from Trigger Times
// -----------------------------------------------------------------------------
// void irtkReconstructionCardiac4D::CalculateSliceCardiacPhases( vector<double>& trigger_times )


// -----------------------------------------------------------------------------
// Initialise Slice Temporal Weights
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::InitSliceTemporalWeights()
{
    _slice_temporal_weight.clear();
    _slice_temporal_weight.resize(_reconstructed_cardiac_phases.size());
    for (unsigned int outputIndex = 0; outputIndex < _reconstructed_cardiac_phases.size(); outputIndex++) 
    {
        _slice_temporal_weight[outputIndex].resize(_slices.size());
    }
}


// -----------------------------------------------------------------------------
// Calculate Slice Temporal Weights
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::CalculateSliceTemporalWeights()
{
    if (_debug)
        cout << "CalculateSliceTemporalWeights" << endl;
    InitSliceTemporalWeights();
    for (unsigned int outputIndex = 0; outputIndex < _reconstructed_cardiac_phases.size(); outputIndex++) 
    {
        for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) 
        {
            _slice_temporal_weight[outputIndex][inputIndex] = CalculateTemporalWeight( _reconstructed_cardiac_phases[outputIndex], _slice_cardphase[inputIndex], _slice_dt[inputIndex], _slice_rr[inputIndex], _wintukeypct ); 
        }      
    }
    /* OBSOLETE
    if (_debug)
    {
        for (int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) 
        {
            if ( inputIndex <= 0 )
                cout<<"input index | output index"<<endl;
            cout<<inputIndex<<"          | ";
            for (int outputIndex = 0; outputIndex < _reconstructed_cardiac_phases.size(); outputIndex++) 
            {
                cout<<"d("<<_reconstructed_cardiac_phases[outputIndex]<<","<<_slice_cardphase[inputIndex]<<")="<<_slice_temporal_weight[outputIndex][inputIndex]<<"; ";
            }    
            cout<<"\b\b."<<endl;
        }
    }
    */
}


// -----------------------------------------------------------------------------
// Calculate Angular Difference
// -----------------------------------------------------------------------------
// Angular difference between output cardiac phase (cardphase0) and slice cardiac phase (cardphase)
double irtkReconstructionCardiac4D::CalculateAngularDifference( double cardphase0, double cardphase )
{
    double angdiff;
    angdiff = ( cardphase - cardphase0 ) - ( 2 * PI ) * floor( ( cardphase - cardphase0 ) / ( 2 * PI ) );
    angdiff = ( angdiff <= PI ) ? angdiff : - ( 2 * PI - angdiff);
    return angdiff;
}


// -----------------------------------------------------------------------------
// Calculate Temporal Weight
// -----------------------------------------------------------------------------
double irtkReconstructionCardiac4D::CalculateTemporalWeight( double cardphase0, double cardphase, double dt, double rr, double alpha )
{  
    // Angular Difference
    double angdiff = CalculateAngularDifference( cardphase0, cardphase );
    
    // Temporal Resolution in Radians
    double dtrad = 2 * PI * dt / rr;
        
    // Temporal Weight
    return sinc( PI * angdiff / dtrad ) * wintukey( angdiff, alpha );
}


// -----------------------------------------------------------------------------
// Sinc Function
// -----------------------------------------------------------------------------
double irtkReconstructionCardiac4D::sinc(double x)
{
    if (x == 0)
        return 1;
    return sin(x)/x;
}

// -----------------------------------------------------------------------------
// Tukey Window Function
// -----------------------------------------------------------------------------
double irtkReconstructionCardiac4D::wintukey( double angdiff, double alpha )
{  
    // angdiff = angular difference (-PI to +PI)
    // alpha   = amount of window with tapered cosine edges (0 to 1)
    if ( fabs( angdiff ) > PI * ( 1 - alpha ) )
       return ( 1 + cos( ( fabs( angdiff ) - PI * ( 1 - alpha ) ) / alpha ) ) / 2;
    return 1;
}


// -----------------------------------------------------------------------------
// Test Temporal Weight Calculation
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::TestTemporalWeightCalculation()
{
    double dt = 75;
    double rr = 400;
    double alpha = 0.3;
    cout<<endl;
    cout<<"TestTemporalWeightCalculation"<<endl;
    cout<<"cardphase0, cardphase1, angular difference, temporal weight"<<endl;
    double cardphase0;
    double cardphase, angdiff, temporalweight;
    for ( int j = 0; j < 10; j++ )
    {
        cardphase0 = double( 2 * PI * j / 10 );
        for ( int i = 0; i < 100; i++ )
        {
            cardphase = double( 2 * PI * i / 100 );
            angdiff = CalculateAngularDifference( cardphase0, cardphase );
            temporalweight = CalculateTemporalWeight( cardphase0, cardphase, dt, rr, alpha ); 
            cout<<cardphase0<<", "<<cardphase<<", "<<angdiff<<", "<<temporalweight<<endl;
        }
    }
    cout<<endl;
}


// -----------------------------------------------------------------------------
// ParallelCoeffInitCardiac4D
// -----------------------------------------------------------------------------
class ParallelCoeffInitCardiac4D {
public:
    irtkReconstructionCardiac4D *reconstructor;

    ParallelCoeffInitCardiac4D(irtkReconstructionCardiac4D *_reconstructor) : 
    reconstructor(_reconstructor) { }

    void operator() (const blocked_range<size_t> &r) const {
        
        for ( size_t inputIndex = r.begin(); inputIndex != r.end(); ++inputIndex ) {

            bool slice_inside;

            //current slice
            //irtkRealImage slice;

            //get resolution of the volume
            double vx, vy, vz;
            reconstructor->_reconstructed4D.GetPixelSize(&vx, &vy, &vz);
            //volume is always isotropic
            double res = vx;
        
            //start of a loop for a slice inputIndex
            cout << inputIndex << " ";
            cout.flush();
            //read the slice
            irtkRealImage& slice = reconstructor->_slices[inputIndex];

            //prepare structures for storage
            POINT3D p;
            VOXELCOEFFS empty;
            SLICECOEFFS slicecoeffs(slice.GetX(), vector < VOXELCOEFFS > (slice.GetY(), empty));

            //to check whether the slice has an overlap with mask ROI
            slice_inside = false;

            //PSF will be calculated in slice space in higher resolution

            //get slice voxel size to define PSF
            double dx, dy, dz;
            slice.GetPixelSize(&dx, &dy, &dz);

            //sigma of 3D Gaussian (sinc with FWHM=dx or dy in-plane, Gaussian with FWHM = dz through-plane)
	    
			double sigmax, sigmay, sigmaz;
			if(reconstructor->_recon_type == _3D)
			{
				  sigmax = 1.2 * dx / 2.3548;
				  sigmay = 1.2 * dy / 2.3548;
				  sigmaz = dz / 2.3548;
			}
	
			if(reconstructor->_recon_type == _1D)
			{
				  sigmax = 0.5 * dx / 2.3548;
				  sigmay = 0.5 * dy / 2.3548;
				  sigmaz = dz / 2.3548;
			}
	
			if(reconstructor->_recon_type == _interpolate)
			{
				  sigmax = 0.5 * dx / 2.3548;
				  sigmay = 0.5 * dx / 2.3548;
				  sigmaz = 0.5 * dx / 2.3548;
			}
            /*
              cout<<"Original sigma"<<sigmax<<" "<<sigmay<<" "<<sigmaz<<endl;
        
              //readjust for resolution of the volume
              //double sigmax,sigmay,sigmaz;
              double sigmamin = res/(3*2.3548);
        
              if((dx-res)>sigmamin)
              sigmax = 1.2 * sqrt(dx*dx-res*res) / 2.3548;
              else sigmax = sigmamin;

              if ((dy-res)>sigmamin)
              sigmay = 1.2 * sqrt(dy*dy-res*res) / 2.3548;
              else
              sigmay=sigmamin;
              if ((dz-1.2*res)>sigmamin)
              sigmaz = sqrt(dz*dz-1.2*1.2*res*res) / 2.3548;
              else sigmaz=sigmamin;
        
              cout<<"Adjusted sigma:"<<sigmax<<" "<<sigmay<<" "<<sigmaz<<endl;
            */
        
            //calculate discretized PSF

            //isotropic voxel size of PSF - derived from resolution of reconstructed volume
            double size = res / reconstructor->_quality_factor;

            //number of voxels in each direction
            //the ROI is 2*voxel dimension

            int xDim = round(2 * dx / size);
            int yDim = round(2 * dy / size);
            int zDim = round(2 * dz / size);
			///test to make dimension alwways odd
			xDim = xDim/2*2+1;
			yDim = yDim/2*2+1;
			zDim = zDim/2*2+1;
			///end test

            //image corresponding to PSF
            irtkImageAttributes attr;
            attr._x = xDim;
            attr._y = yDim;
            attr._z = zDim;
            attr._dx = size;
            attr._dy = size;
            attr._dz = size;
            irtkRealImage PSF(attr);

            //centre of PSF
            double cx, cy, cz;
            cx = 0.5 * (xDim - 1);
            cy = 0.5 * (yDim - 1);
            cz = 0.5 * (zDim - 1);
            PSF.ImageToWorld(cx, cy, cz);

            double x, y, z;
            double sum = 0;
            int i, j, k;
            for (i = 0; i < xDim; i++)
                for (j = 0; j < yDim; j++)
                    for (k = 0; k < zDim; k++) {
                        x = i;
                        y = j;
                        z = k;
                        PSF.ImageToWorld(x, y, z);
                        x -= cx;
                        y -= cy;
                        z -= cz;
                        //continuous PSF does not need to be normalized as discrete will be
                        PSF(i, j, k) = exp(
                                           -x * x / (2 * sigmax * sigmax) - y * y / (2 * sigmay * sigmay)
                                           - z * z / (2 * sigmaz * sigmaz));
                        sum += PSF(i, j, k);
                    }
            PSF /= sum;

            if (reconstructor->_debug)
                if (inputIndex == 0)
                    PSF.Write("PSF.nii.gz");

            //prepare storage for PSF transformed and resampled to the space of reconstructed volume
            //maximum dim of rotated kernel - the next higher odd integer plus two to accound for rounding error of tx,ty,tz.
            //Note conversion from PSF image coordinates to tPSF image coordinates *size/res
            int dim = (floor(ceil(sqrt(double(xDim * xDim + yDim * yDim + zDim * zDim)) * size / res) / 2))
                * 2 + 1 + 2;
            //prepare image attributes. Voxel dimension will be taken from the reconstructed volume
            attr._x = dim;
            attr._y = dim;
            attr._z = dim;
            attr._dx = res;
            attr._dy = res;
            attr._dz = res;
            //create matrix from transformed PSF
            irtkRealImage tPSF(attr);
            //calculate centre of tPSF in image coordinates
            int centre = (dim - 1) / 2;

            //for each voxel in current slice calculate matrix coefficients
            int ii, jj, kk;
            int tx, ty, tz;
            int nx, ny, nz;
            int l, m, n;
            double weight;
            for (i = 0; i < slice.GetX(); i++)
                for (j = 0; j < slice.GetY(); j++)
                    if (slice(i, j, 0) != -1) {
                        //calculate centrepoint of slice voxel in volume space (tx,ty,tz)
                        x = i;
                        y = j;
                        z = 0;
                        slice.ImageToWorld(x, y, z);
                        reconstructor->_transformations[inputIndex].Transform(x, y, z);
                        reconstructor->_reconstructed4D.WorldToImage(x, y, z);
                        tx = round(x);
                        ty = round(y);
                        tz = round(z);

                        //Clear the transformed PSF
                        for (ii = 0; ii < dim; ii++)
                            for (jj = 0; jj < dim; jj++)
                                for (kk = 0; kk < dim; kk++)
                                    tPSF(ii, jj, kk) = 0;

                        //for each POINT3D of the PSF
                        for (ii = 0; ii < xDim; ii++)
                            for (jj = 0; jj < yDim; jj++)
                                for (kk = 0; kk < zDim; kk++) {
                                    //Calculate the position of the POINT3D of
                                    //PSF centered over current slice voxel                            
                                    //This is a bit complicated because slices
                                    //can be oriented in any direction 

                                    //PSF image coordinates
                                    x = ii;
                                    y = jj;
                                    z = kk;
                                    //change to PSF world coordinates - now real sizes in mm
                                    PSF.ImageToWorld(x, y, z);
                                    //centre around the centrepoint of the PSF
                                    x -= cx;
                                    y -= cy;
                                    z -= cz;

                                    //Need to convert (x,y,z) to slice image
                                    //coordinates because slices can have
                                    //transformations included in them (they are
                                    //nifti)  and those are not reflected in
                                    //PSF. In slice image coordinates we are
                                    //sure that z is through-plane 

                                    //adjust according to voxel size
                                    x /= dx;
                                    y /= dy;
                                    z /= dz;
                                    //center over current voxel
                                    x += i;
                                    y += j;

                                    //convert from slice image coordinates to world coordinates
                                    slice.ImageToWorld(x, y, z);

                                    //x+=(vx-cx); y+=(vy-cy); z+=(vz-cz);
                                    //Transform to space of reconstructed volume
                                    reconstructor->_transformations[inputIndex].Transform(x, y, z);
                                    //Change to image coordinates
                                    reconstructor->_reconstructed4D.WorldToImage(x, y, z);

                                    //determine coefficients of volume voxels for position x,y,z
                                    //using linear interpolation

                                    //Find the 8 closest volume voxels

                                    //lowest corner of the cube
                                    nx = (int) floor(x);
                                    ny = (int) floor(y);
                                    nz = (int) floor(z);

                                    //not all neighbours might be in ROI, thus we need to normalize
                                    //(l,m,n) are image coordinates of 8 neighbours in volume space
                                    //for each we check whether it is in volume
                                    sum = 0;
                                    //to find wether the current slice voxel has overlap with ROI
                                    bool inside = false;
                                    for (l = nx; l <= nx + 1; l++)
                                        if ((l >= 0) && (l < reconstructor->_reconstructed4D.GetX()))
                                            for (m = ny; m <= ny + 1; m++)
                                                if ((m >= 0) && (m < reconstructor->_reconstructed4D.GetY()))
                                                    for (n = nz; n <= nz + 1; n++)
                                                        if ((n >= 0) && (n < reconstructor->_reconstructed4D.GetZ())) {
                                                            weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) * (1 - fabs(n - z));
                                                            sum += weight;
                                                            if (reconstructor->_mask(l, m, n) == 1) {
                                                                inside = true;
                                                                slice_inside = true;
                                                            }
                                                        }
                                    //if there were no voxels do nothing
                                    if ((sum <= 0) || (!inside))
                                        continue;
                                    //now calculate the transformed PSF
                                    for (l = nx; l <= nx + 1; l++)
                                        if ((l >= 0) && (l < reconstructor->_reconstructed4D.GetX()))
                                            for (m = ny; m <= ny + 1; m++)
                                                if ((m >= 0) && (m < reconstructor->_reconstructed4D.GetY()))
                                                    for (n = nz; n <= nz + 1; n++)
                                                        if ((n >= 0) && (n < reconstructor->_reconstructed4D.GetZ())) {
                                                            weight = (1 - fabs(l - x)) * (1 - fabs(m - y)) * (1 - fabs(n - z));

                                                            //image coordinates in tPSF
                                                            //(centre,centre,centre) in tPSF is aligned with (tx,ty,tz)
                                                            int aa, bb, cc;
                                                            aa = l - tx + centre;
                                                            bb = m - ty + centre;
                                                            cc = n - tz + centre;

                                                            //resulting value
                                                            double value = PSF(ii, jj, kk) * weight / sum;

                                                            //Check that we are in tPSF
                                                            if ((aa < 0) || (aa >= dim) || (bb < 0) || (bb >= dim) || (cc < 0)
                                                                || (cc >= dim)) {
                                                                cerr << "Error while trying to populate tPSF. " << aa << " " << bb
                                                                     << " " << cc << endl;
                                                                cerr << l << " " << m << " " << n << endl;
                                                                cerr << tx << " " << ty << " " << tz << endl;
                                                                cerr << centre << endl;
                                                                tPSF.Write("tPSF.nii.gz");
                                                                exit(1);
                                                            }
                                                            else
                                                                //update transformed PSF
                                                                tPSF(aa, bb, cc) += value;
                                                        }

                                } //end of the loop for PSF points

                        //store tPSF values
                        for (ii = 0; ii < dim; ii++)
                            for (jj = 0; jj < dim; jj++)
                                for (kk = 0; kk < dim; kk++)
                                    if (tPSF(ii, jj, kk) > 0) {
                                        p.x = ii + tx - centre;
                                        p.y = jj + ty - centre;
                                        p.z = kk + tz - centre;
                                        p.value = tPSF(ii, jj, kk);
                                        slicecoeffs[i][j].push_back(p);
                                    }
                    } //end of loop for slice voxels

            reconstructor->_volcoeffs[inputIndex] = slicecoeffs;
            reconstructor->_slice_inside[inputIndex] = slice_inside;

        }  //end of loop through the slices                            
        
    }

    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, reconstructor->_slices.size() ),
                      *this );
        init.terminate();
    }

};


// -----------------------------------------------------------------------------
// Calculate Transformation Matrix Between Slices and Voxels
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::CoeffInitCardiac4D()
{
    if (_debug)
        cout << "CoeffInit" << endl;
    
    //clear slice-volume matrix from previous iteration
    _volcoeffs.clear();
    _volcoeffs.resize(_slices.size());

    //clear indicator of slice having and overlap with volumetric mask
    _slice_inside.clear();
    _slice_inside.resize(_slices.size());

    cout << "Initialising matrix coefficients...";
    cout.flush();
    ParallelCoeffInitCardiac4D coeffinit(this);
    coeffinit();
    cout << " ... done." << endl;

    //prepare image for volume weights, will be needed for Gaussian Reconstruction
    cout << "Computing 4D volume weights..." << endl;
    irtkImageAttributes volAttr = _reconstructed4D.GetImageAttributes();
    _volume_weights.Initialize( volAttr );
    _volume_weights = 0;

    // TODO: investigate if this loop is taking a long time to compute, and consider parallelisation
    int i, j, n, k, outputIndex;
    unsigned int inputIndex;
    POINT3D p;
    cout << "    ... for input slice: ";
    for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
        cout << inputIndex << ", ";
        cout.flush();
        for ( i = 0; i < _slices[inputIndex].GetX(); i++)
            for ( j = 0; j < _slices[inputIndex].GetY(); j++) {
                n = _volcoeffs[inputIndex][i][j].size();
                for (k = 0; k < n; k++) {
                  p = _volcoeffs[inputIndex][i][j][k];
                  for (outputIndex=0; outputIndex<_reconstructed4D.GetT(); outputIndex++)
                  {
                      _volume_weights(p.x, p.y, p.z, outputIndex) += _slice_temporal_weight[outputIndex][inputIndex] * p.value;
                  }
                }
            }
    }
    cout << "\b\b." << endl;
    // if (_debug)
    //     _volume_weights.Write("volume_weights.nii.gz");
    
    //find average volume weight to modify alpha parameters accordingly
    double sum = 0;
    int num=0;
    for (i=0; i<_volume_weights.GetX(); i++)
      for (j=0; j<_volume_weights.GetY(); j++)
        for (k=0; k<_volume_weights.GetZ(); k++)
          if (_mask(i,j,k)==1)
            for (int f=0; f<_volume_weights.GetT(); f++) {
              sum += _volume_weights(i,j,k,f);
              num++;
            }
    
    _average_volume_weight = sum/num;
    
    if(_debug) {
        cout<<"Average volume weight is "<<_average_volume_weight<<endl;
    }
    
} 


// -----------------------------------------------------------------------------
// PSF-Weighted Reconstruction
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::GaussianReconstructionCardiac4D()
{
    if(_debug)
    {
      cout << "Gaussian reconstruction ... " << endl;
      cout << "\tinput slice:  ";
      cout.flush();
    }
    unsigned int inputIndex, outputIndex;
    int k, n;
    irtkRealImage slice;
    double scale;
    POINT3D p;
    vector<int> voxel_num;  
    int slice_vox_num;

    //clear _reconstructed image
    _reconstructed4D = 0;

    for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
        
        if(_debug)
        {
          cout << inputIndex << ", ";
          cout.flush();
        }
        //copy the current slice
        slice = _slices[inputIndex];
        //alias the current bias image
        irtkRealImage& b = _bias[inputIndex];
        //read current scale factor
        scale = _scale[inputIndex];
        
        slice_vox_num=0;

        //Distribute slice intensities to the volume
        for (int i = 0; i < slice.GetX(); i++)
            for (int j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    //biascorrect and scale the slice
                    slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;

                    //number of volume voxels with non-zero coefficients
                    //for current slice voxel
                    n = _volcoeffs[inputIndex][i][j].size();

                    //if given voxel is not present in reconstructed volume at all,
                    //pad it
                    
                    //if (n == 0)
                    //_slices[inputIndex].PutAsDouble(i, j, 0, -1);
                    //calculate num of vox in a slice that have overlap with roi
                    if (n>0)
                        slice_vox_num++;

                    //add contribution of current slice voxel to all voxel volumes
                    //to which it contributes
                    for (k = 0; k < n; k++) {
                        for (outputIndex=0; outputIndex<_reconstructed_cardiac_phases.size(); outputIndex++)
                        {
                            p = _volcoeffs[inputIndex][i][j][k];
                            _reconstructed4D(p.x, p.y, p.z, outputIndex) += _slice_temporal_weight[outputIndex][inputIndex] * p.value * slice(i, j, 0);
                        }
                    }
                }
        voxel_num.push_back(slice_vox_num);
        //end of loop for a slice inputIndex
    }

    //normalize the volume by proportion of contributing slice voxels
    //for each volume voxe
    _reconstructed4D /= _volume_weights;
    
    if(_debug)
    {
      cout << inputIndex << "\b\b." << endl;
      cout << "... Gaussian reconstruction done." << endl << endl;
      cout.flush();
    }    
    

    // if (_debug)
    //     _reconstructed4D.Write("init.nii.gz");

    //now find slices with small overlap with ROI and exclude them.
    
    vector<int> voxel_num_tmp;
    for (unsigned int i=0;i<voxel_num.size();i++)
        voxel_num_tmp.push_back(voxel_num[i]);
    
    //find median
    sort(voxel_num_tmp.begin(),voxel_num_tmp.end());
    int median = voxel_num_tmp[round(voxel_num_tmp.size()*0.5)];
    
    //remember slices with small overlap with ROI
    _small_slices.clear();
    for (unsigned int i=0;i<voxel_num.size();i++)
        if (voxel_num[i]<0.1*median)
            _small_slices.push_back(i);
    
    if (_debug) {
        cout<<"Small slices:";
        for (unsigned int i=0;i<_small_slices.size();i++)
            cout<<" "<<_small_slices[i];
        cout<<endl;
    }
}


// -----------------------------------------------------------------------------
// Scale Volume
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::ScaleVolumeCardiac4D()
{
    if (_debug)
        cout << "Scaling volume: ";
    
    unsigned int inputIndex;
    int i, j;
    double scalenum = 0, scaleden = 0;

    for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        // alias for the current slice
        irtkRealImage& slice = _slices[inputIndex];

        //alias for the current weight image
        irtkRealImage& w = _weights[inputIndex];

        // alias for the current simulated slice
        irtkRealImage& sim = _simulated_slices[inputIndex];
        
        for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    //scale - intensity matching
                    if ( _simulated_weights[inputIndex](i,j,0) > 0.99 ) {
                        scalenum += w(i, j, 0) * _slice_weight[inputIndex] * slice(i, j, 0) * sim(i, j, 0);
                        scaleden += w(i, j, 0) * _slice_weight[inputIndex] * sim(i, j, 0) * sim(i, j, 0);
                    }
                }
    } //end of loop for a slice inputIndex
    
    //calculate scale for the volume
    double scale = scalenum / scaleden;
  
    if(_debug)
        cout<<" scale = "<<scale;
  
    irtkRealPixel *ptr = _reconstructed4D.GetPointerToVoxels();
    for(i=0;i<_reconstructed4D.GetNumberOfVoxels();i++) {
        if(*ptr>0) *ptr = *ptr * scale;
        ptr++;
    }
    cout<<endl;
}


// -----------------------------------------------------------------------------
// Parallel Simulate Slices
// -----------------------------------------------------------------------------
class ParallelSimulateSlicesCardiac4D {
    irtkReconstructionCardiac4D *reconstructor;
        
public:
    ParallelSimulateSlicesCardiac4D( irtkReconstructionCardiac4D *_reconstructor ) : 
    reconstructor(_reconstructor) { }

    void operator() (const blocked_range<size_t> &r) const {
        for ( size_t inputIndex = r.begin(); inputIndex != r.end(); ++inputIndex ) {
            //Calculate simulated slice
            reconstructor->_simulated_slices[inputIndex].Initialize( reconstructor->_slices[inputIndex].GetImageAttributes() );
            reconstructor->_simulated_slices[inputIndex] = 0;

            reconstructor->_simulated_weights[inputIndex].Initialize( reconstructor->_slices[inputIndex].GetImageAttributes() );
            reconstructor->_simulated_weights[inputIndex] = 0;

            reconstructor->_simulated_inside[inputIndex].Initialize( reconstructor->_slices[inputIndex].GetImageAttributes() );
            reconstructor->_simulated_inside[inputIndex] = 0;            

            reconstructor->_slice_inside[inputIndex] = false;
            
            POINT3D p;
            for ( int i = 0; i < reconstructor->_slices[inputIndex].GetX(); i++ )
                for ( int j = 0; j < reconstructor->_slices[inputIndex].GetY(); j++ )
                    if ( reconstructor->_slices[inputIndex](i, j, 0) != -1 ) {
                        double weight = 0;
                        int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                        for ( int k = 0; k < n; k++ ) {
                            p = reconstructor->_volcoeffs[inputIndex][i][j][k];
                            for ( int outputIndex = 0; outputIndex < reconstructor->_reconstructed4D.GetT(); outputIndex++ ) {
                                reconstructor->_simulated_slices[inputIndex](i, j, 0) += reconstructor->_slice_temporal_weight[outputIndex][inputIndex] * p.value * reconstructor->_reconstructed4D(p.x, p.y, p.z, outputIndex);
                                weight += reconstructor->_slice_temporal_weight[outputIndex][inputIndex] * p.value;
                            }
                            if (reconstructor->_mask(p.x, p.y, p.z) == 1) {
                                reconstructor->_simulated_inside[inputIndex](i, j, 0) = 1;
                                reconstructor->_slice_inside[inputIndex] = true;
                            }
                        }                    
                        if( weight > 0 ) {
                            reconstructor->_simulated_slices[inputIndex](i,j,0) /= weight;
                            reconstructor->_simulated_weights[inputIndex](i,j,0) = weight;
                        }
                    }
            
        }
    }
    
    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, reconstructor->_slices.size() ),
                      *this );
        init.terminate();
    }

};


// -----------------------------------------------------------------------------
// Simulate Slices
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SimulateSlicesCardiac4D()
{
  if (_debug)
      cout<<"Simulating Slices..."<<endl;

  ParallelSimulateSlicesCardiac4D parallelSimulateSlices( this );
  parallelSimulateSlices();

  if (_debug)
      cout<<"\t...Simulating Slices done."<<endl;   
}


// -----------------------------------------------------------------------------
// Parallel Processing Class for Normalise Bias
// -----------------------------------------------------------------------------
class ParallelNormaliseBiasCardiac4D{
    irtkReconstructionCardiac4D* reconstructor;

public:
    irtkRealImage bias;
    irtkRealImage volweight3d;

void operator()( const blocked_range<size_t>& r ) {
    for ( size_t inputIndex = r.begin(); inputIndex < r.end(); ++inputIndex) {

        if(reconstructor->_debug) {
            cout<<inputIndex<<" ";
        }
        
        // alias the current slice
        irtkRealImage& slice = reconstructor->_slices[inputIndex];
            
        //read the current bias image
        irtkRealImage b = reconstructor->_bias[inputIndex];
            
        //read current scale factor
        double scale = reconstructor->_scale[inputIndex];

        irtkRealPixel *pi = slice.GetPointerToVoxels();
        irtkRealPixel *pb = b.GetPointerToVoxels();
        for(int i = 0; i<slice.GetNumberOfVoxels(); i++) {
            if((*pi>-1)&&(scale>0))
                *pb -= log(scale);
            pb++;
            pi++;
        }
            
        //Distribute slice intensities to the volume
        POINT3D p;
        for (int i = 0; i < slice.GetX(); i++)
            for (int j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    //number of volume voxels with non-zero coefficients for current slice voxel
                    int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                    //add contribution of current slice voxel to all voxel volumes
                    //to which it contributes
                    for (int k = 0; k < n; k++) {
                        p = reconstructor->_volcoeffs[inputIndex][i][j][k];
                        bias(p.x, p.y, p.z) += p.value * b(i, j, 0);
                        volweight3d(p.x,p.y,p.z) += p.value;
                    }
                }
        //end of loop for a slice inputIndex                
    }
}

ParallelNormaliseBiasCardiac4D( ParallelNormaliseBiasCardiac4D& x, split ) :
    reconstructor(x.reconstructor)
{
    irtkImageAttributes attr = reconstructor->_reconstructed4D.GetImageAttributes();
    attr._t = 1;
    bias.Initialize( attr );
    bias = 0;   
    volweight3d.Initialize( attr );
    volweight3d = 0;
}

void join( const ParallelNormaliseBiasCardiac4D& y ) {       
    bias += y.bias;
    volweight3d += y.volweight3d;
}
         
ParallelNormaliseBiasCardiac4D( irtkReconstructionCardiac4D *reconstructor ) :
reconstructor(reconstructor)
{
    irtkImageAttributes attr = reconstructor->_reconstructed4D.GetImageAttributes();
    attr._t = 1;
    bias.Initialize( attr );
    bias = 0; 
    volweight3d.Initialize( attr );
    volweight3d = 0;
}

// execute
void operator() () {
    task_scheduler_init init(tbb_no_threads);
    parallel_reduce( blocked_range<size_t>(0,reconstructor->_slices.size()),
                     *this );
    init.terminate();
}        
};


// -----------------------------------------------------------------------------
// Normalise Bias
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::NormaliseBiasCardiac4D(int iter, int rec_iter)
{
    if(_debug)
        cout << "Normalise Bias ... ";

    ParallelNormaliseBiasCardiac4D parallelNormaliseBias(this);
    parallelNormaliseBias();
    irtkRealImage bias = parallelNormaliseBias.bias;
    irtkRealImage volweight3d = parallelNormaliseBias.volweight3d;

    // normalize the volume by proportion of contributing slice voxels for each volume voxel
    bias /= volweight3d;
        
    if(_debug)
        cout << "done." << endl;
    
    bias = StaticMaskVolume4D(bias,0);
    irtkRealImage m = _mask;
    irtkGaussianBlurring<irtkRealPixel> gb(_sigma_bias);
    gb.SetInput(&bias);
    gb.SetOutput(&bias);
    gb.Run();
    gb.SetInput(&m);
    gb.SetOutput(&m);
    gb.Run();
    
    bias/=m;
    
    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebias_mc%02isr%02i.nii.gz",iter,rec_iter);
        bias.Write(buffer);
    }
    
    for ( int i = 0; i < _reconstructed4D.GetX(); i++)  
      for ( int j = 0; j < _reconstructed4D.GetY(); j++)  
        for ( int k = 0; k < _reconstructed4D.GetZ(); k++)  
          for ( int f = 0; f < _reconstructed4D.GetT(); f++)  
            if(_reconstructed4D(i,j,k,f)!=-1) 
            _reconstructed4D(i,j,k,f) /=exp(-bias(i,j,k));

}


// -----------------------------------------------------------------------------
// Calculate Target Cardiac Phase in Reconstructed Volume for Slice-To-Volume Registration
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::CalculateSliceToVolumeTargetCardiacPhase()
{
  int card_index;
  double angdiff;
  if (_debug)
    cout << "CalculateSliceToVolumeTargetCardiacPhase" << endl;
  _slice_svr_card_index.clear();    
  for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) 
  {
    angdiff = PI + 0.001;  // NOTE: init angdiff larger than any possible calculated angular difference 
    card_index = -1;
    for (unsigned int outputIndex = 0; outputIndex < _reconstructed_cardiac_phases.size(); outputIndex++)
    {
      if ( fabs( CalculateAngularDifference( _reconstructed_cardiac_phases[outputIndex], _slice_cardphase[inputIndex] ) ) < angdiff)
      {
        angdiff = fabs( CalculateAngularDifference( _reconstructed_cardiac_phases[outputIndex], _slice_cardphase[inputIndex] ) );
        card_index = outputIndex;
      }
    }
    _slice_svr_card_index.push_back(card_index); 
    // if (_debug)
    //   cout << inputIndex << ":" << _slice_svr_card_index[inputIndex] << ", ";
  }
  // if (_debug)
  //   cout << "\b\b." << endl;
}


// -----------------------------------------------------------------------------
// Parallel Slice-to-Volume Registration Class
// -----------------------------------------------------------------------------
class ParallelSliceToVolumeRegistrationCardiac4D {
public:
    irtkReconstructionCardiac4D *reconstructor;

    ParallelSliceToVolumeRegistrationCardiac4D(irtkReconstructionCardiac4D *_reconstructor) : 
    reconstructor(_reconstructor) { }

    void operator() (const blocked_range<size_t> &r) const {

        irtkImageAttributes attr = reconstructor->_reconstructed4D.GetImageAttributes();
        
        for ( size_t inputIndex = r.begin(); inputIndex != r.end(); ++inputIndex ) {

            irtkImageRigidRegistrationWithPadding registration;
            irtkGreyPixel smin, smax;
            irtkGreyImage target;
            irtkRealImage slice, w, b, t;
            irtkResamplingWithPadding<irtkRealPixel> resampling(attr._dx,attr._dx,attr._dx,-1);         
            irtkReconstruction dummy_reconstruction;
            
            // TARGET
            // get current slice   
            t = reconstructor->_slices[inputIndex];
            // resample to spatial resolution of reconstructed volume
            resampling.SetInput(&reconstructor->_slices[inputIndex]);
            resampling.SetOutput(&t);
            resampling.Run();
            target=t;
            // get pixel value min and max
            target.GetMinMax(&smin, &smax);
        
            // SOURCE
            if (smax > -1) {
                // put origin to zero
                irtkRigidTransformation offset;
                dummy_reconstruction.ResetOrigin(target,offset);
                irtkMatrix mo = offset.GetMatrix();
                irtkMatrix m = reconstructor->_transformations[inputIndex].GetMatrix();
                m=m*mo;
                reconstructor->_transformations[inputIndex].PutMatrix(m);

                // TODO: extract nearest cardiac phase from reconstructed 4D to use as source
                irtkGreyImage source = reconstructor->_reconstructed4D.GetRegion( 0, 0, 0, reconstructor->_slice_svr_card_index[inputIndex], attr._x, attr._y, attr._z, reconstructor->_slice_svr_card_index[inputIndex]+1 ); 
                registration.SetInput(&target, &source);
                registration.SetOutput(&reconstructor->_transformations[inputIndex]);
                registration.GuessParameterSliceToVolume();
                registration.SetTargetPadding(-1);
                registration.Run();
                //undo the offset
                mo.Invert();
                m = reconstructor->_transformations[inputIndex].GetMatrix();
                m=m*mo;
                reconstructor->_transformations[inputIndex].PutMatrix(m);
            }      
        }
    }

    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, reconstructor->_slices.size() ),
                      *this );
        init.terminate();
    }

};


// -----------------------------------------------------------------------------
// Slice-to-Volume Registration 
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SliceToVolumeRegistrationCardiac4D()
{
  if (_debug)
      cout << "SliceToVolumeRegistrationCardiac4D" << endl;
  ParallelSliceToVolumeRegistrationCardiac4D registration(this);
  registration();
}


// -----------------------------------------------------------------------------
// Mask Reconstructed Volume
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::StaticMaskReconstructedVolume4D()
{
  for ( int i = 0; i < _mask.GetX(); i++) {
    for ( int j = 0; j < _mask.GetY(); j++) {
      for ( int k = 0; k < _mask.GetZ(); k++) {
        if ( _mask(i,j,k) == 0 ) {
          for ( int t = 0; t < _reconstructed4D.GetT(); t++) {
            _reconstructed4D(i,j,k,t) = -1;
          }
        }
      }
    }
  }
}


// -----------------------------------------------------------------------------
// Apply Static Mask to 4D Volume
// -----------------------------------------------------------------------------
irtkRealImage irtkReconstructionCardiac4D::StaticMaskVolume4D(irtkRealImage volume, double padding)
{
  for ( int i = 0; i < volume.GetX(); i++) {
    for ( int j = 0; j < volume.GetY(); j++) {
      for ( int k = 0; k < volume.GetZ(); k++) {
        if ( _mask(i,j,k) == 0 ) {
          for ( int t = 0; t < volume.GetT(); t++) {
            volume(i,j,k,t) = padding;
          }
        }
      }
    }
  }
  return volume;
}


// -----------------------------------------------------------------------------
// Parallel Super-Resolution Class
// -----------------------------------------------------------------------------
class ParallelSuperresolutionCardiac4D {
    irtkReconstructionCardiac4D* reconstructor;
public:
    irtkRealImage confidence_map;
    irtkRealImage addon;
    
    void operator()( const blocked_range<size_t>& r ) {
        for ( size_t inputIndex = r.begin(); inputIndex < r.end(); ++inputIndex) {
            // read the current slice
            irtkRealImage slice = reconstructor->_slices[inputIndex];
                
            //read the current weight image
            irtkRealImage& w = reconstructor->_weights[inputIndex];
                
            //read the current bias image
            irtkRealImage& b = reconstructor->_bias[inputIndex];
                
            //identify scale factor
            double scale = reconstructor->_scale[inputIndex];

            //Update reconstructed volume using current slice

            //Distribute error to the volume
            POINT3D p;
            for ( int i = 0; i < slice.GetX(); i++)
                for ( int j = 0; j < slice.GetY(); j++)
                    if (slice(i, j, 0) != -1) {
                        //bias correct and scale the slice
                        slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;
                        
                        if ( reconstructor->_simulated_slices[inputIndex](i,j,0) > 0 )
                            slice(i,j,0) -= reconstructor->_simulated_slices[inputIndex](i,j,0);
                        else
                            slice(i,j,0) = 0;

                        int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                        for (int k = 0; k < n; k++) {
                            p = reconstructor->_volcoeffs[inputIndex][i][j][k];
                            for (int outputIndex=0; outputIndex<reconstructor->_reconstructed4D.GetT(); outputIndex++) {
			    if(reconstructor->_robust_slices_only)
			    {
                              addon(p.x, p.y, p.z, outputIndex) += reconstructor->_slice_temporal_weight[outputIndex][inputIndex] * p.value * slice(i, j, 0) * reconstructor->_slice_weight[inputIndex];
                              confidence_map(p.x, p.y, p.z, outputIndex) += reconstructor->_slice_temporal_weight[outputIndex][inputIndex] * p.value * reconstructor->_slice_weight[inputIndex];
			      
			    }
			    else
			    {
                              addon(p.x, p.y, p.z, outputIndex) += reconstructor->_slice_temporal_weight[outputIndex][inputIndex] * p.value * slice(i, j, 0) * w(i, j, 0) * reconstructor->_slice_weight[inputIndex];
                              confidence_map(p.x, p.y, p.z, outputIndex) += reconstructor->_slice_temporal_weight[outputIndex][inputIndex] * p.value * w(i, j, 0) * reconstructor->_slice_weight[inputIndex];
			    }
                            }
                        }
                    }
        } //end of loop for a slice inputIndex
    }
 
    ParallelSuperresolutionCardiac4D( ParallelSuperresolutionCardiac4D& x, split ) :
        reconstructor(x.reconstructor)
    {
        //Clear addon
        addon.Initialize( reconstructor->_reconstructed4D.GetImageAttributes() );
        addon = 0;

        //Clear confidence map
        confidence_map.Initialize( reconstructor->_reconstructed4D.GetImageAttributes() );
        confidence_map = 0;
    }
 
    void join( const ParallelSuperresolutionCardiac4D& y ) {
        addon += y.addon;
        confidence_map += y.confidence_map;
    }
             
    ParallelSuperresolutionCardiac4D( irtkReconstructionCardiac4D *reconstructor ) :
    reconstructor(reconstructor)
    {
        //Clear addon
        addon.Initialize( reconstructor->_reconstructed4D.GetImageAttributes() );
        addon = 0;

        //Clear confidence map
        confidence_map.Initialize( reconstructor->_reconstructed4D.GetImageAttributes() );
        confidence_map = 0;
    }

    // execute
    void operator() () {
        task_scheduler_init init(tbb_no_threads);
        parallel_reduce( blocked_range<size_t>(0,reconstructor->_slices.size()),
                         *this );
        init.terminate();
    }         
};


// -----------------------------------------------------------------------------
// Super-Resolution of 4D Volume
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SuperresolutionCardiac4D( int iter )
{
  if (_debug)
      cout << "Superresolution " << iter << endl;
  
  int i, j, k, t;
  irtkRealImage addon, original;
  

  //Remember current reconstruction for edge-preserving smoothing
  original = _reconstructed4D;

  ParallelSuperresolutionCardiac4D parallelSuperresolution(this);
  parallelSuperresolution();

  addon = parallelSuperresolution.addon;
  _confidence_map = parallelSuperresolution.confidence_map;
  
  if(_debug) {
      // char buffer[256];
      //sprintf(buffer,"confidence-map%i.nii.gz",iter);
      //_confidence_map.Write(buffer);
      _confidence_map.Write("confidence-map.nii.gz");
      //sprintf(buffer,"addon%i.nii.gz",iter);
      //addon.Write(buffer);
  }

  if (!_adaptive) 
      for (i = 0; i < addon.GetX(); i++)
          for (j = 0; j < addon.GetY(); j++)
              for (k = 0; k < addon.GetZ(); k++)
                  for (t = 0; t < addon.GetT(); t++)
                    if (_confidence_map(i, j, k, t) > 0) {
                        // ISSUES if _confidence_map(i, j, k, t) is too small leading
                        // to bright pixels
                        addon(i, j, k, t) /= _confidence_map(i, j, k, t);
                        //this is to revert to normal (non-adaptive) regularisation
                        _confidence_map(i,j,k,t) = 1;
                    }

  _reconstructed4D += addon * _alpha; //_average_volume_weight;
  
  //bound the intensities
  for (i = 0; i < _reconstructed4D.GetX(); i++)
      for (j = 0; j < _reconstructed4D.GetY(); j++)
          for (k = 0; k < _reconstructed4D.GetZ(); k++) 
              for (t = 0; t < _reconstructed4D.GetT(); t++)
          {
              if (_reconstructed4D(i, j, k, t) < _min_intensity * 0.9)
                  _reconstructed4D(i, j, k, t) = _min_intensity * 0.9;
              if (_reconstructed4D(i, j, k, t) > _max_intensity * 1.1)
                  _reconstructed4D(i, j, k, t) = _max_intensity * 1.1;
          }

  //Smooth the reconstructed image
  AdaptiveRegularizationCardiac4D(iter, original);
    
  //Remove the bias in the reconstructed volume compared to previous iteration
  /* TODO: update adaptive regularisation for 4d
  if (_global_bias_correction)
      BiasCorrectVolume(original);
  */
}


// -----------------------------------------------------------------------------
// Parallel Adaptive Regularization Class 1: calculate smoothing factor, b
// -----------------------------------------------------------------------------
class ParallelAdaptiveRegularization1Cardiac4D {
    irtkReconstructionCardiac4D *reconstructor;
    vector<irtkRealImage> &b;
    vector<double> &factor;
    irtkRealImage &original;
        
public:
    ParallelAdaptiveRegularization1Cardiac4D( irtkReconstructionCardiac4D *_reconstructor,
                                     vector<irtkRealImage> &_b,
                                     vector<double> &_factor,
                                     irtkRealImage &_original) : 
        reconstructor(_reconstructor),
        b(_b),
        factor(_factor),
        original(_original) { }

    void operator() (const blocked_range<size_t> &r) const {
        int dx = reconstructor->_reconstructed4D.GetX();
        int dy = reconstructor->_reconstructed4D.GetY();
        int dz = reconstructor->_reconstructed4D.GetZ();
        int dt = reconstructor->_reconstructed4D.GetT();
        for ( size_t i = r.begin(); i != r.end(); ++i ) {
            //b[i] = reconstructor->_reconstructed;
            // b[i].Initialize( reconstructor->_reconstructed.GetImageAttributes() );

            int x, y, z, xx, yy, zz, t;
            double diff;
            for (x = 0; x < dx; x++)
                for (y = 0; y < dy; y++)
                    for (z = 0; z < dz; z++) {
                        xx = x + reconstructor->_directions[i][0];
                        yy = y + reconstructor->_directions[i][1];
                        zz = z + reconstructor->_directions[i][2];
                        for (t = 0; t < dt; t++) {
                            if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)
                                && (reconstructor->_confidence_map(x, y, z, t) > 0) && (reconstructor->_confidence_map(xx, yy, zz, t) > 0)) {
                                diff = (original(xx, yy, zz, t) - original(x, y, z, t)) * sqrt(factor[i]) / reconstructor->_delta;
                                b[i](x, y, z, t) = factor[i] / sqrt(1 + diff * diff);

                            }
                            else
                                b[i](x, y, z, t) = 0;
                        }
                    }
        }
    }

    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, 13),
                      *this );
        init.terminate();
    }

};


// -----------------------------------------------------------------------------
// Parallel Adaptive Regularization Class 2: compute regularisation update
// -----------------------------------------------------------------------------
class ParallelAdaptiveRegularization2Cardiac4D {
    irtkReconstructionCardiac4D *reconstructor;
    vector<irtkRealImage> &b;
    vector<double> &factor;
    irtkRealImage &original;
        
public:
    ParallelAdaptiveRegularization2Cardiac4D( irtkReconstructionCardiac4D *_reconstructor,
                                     vector<irtkRealImage> &_b,
                                     vector<double> &_factor,
                                     irtkRealImage &_original) : 
        reconstructor(_reconstructor),
        b(_b),
        factor(_factor),
        original(_original) { }

    void operator() (const blocked_range<size_t> &r) const {
      int dx = reconstructor->_reconstructed4D.GetX();
      int dy = reconstructor->_reconstructed4D.GetY();
      int dz = reconstructor->_reconstructed4D.GetZ();
      int dt = reconstructor->_reconstructed4D.GetT();
      for ( size_t x = r.begin(); x != r.end(); ++x ) {
        int xx, yy, zz;
        for (int y = 0; y < dy; y++)
          for (int z = 0; z < dz; z++) 
            for (int t = 0; t < dt; t++) {
              if(reconstructor->_confidence_map(x,y,z,t)>0)
		          {                    
                double val = 0;
                double sum = 0;
                for (int i = 0; i < 13; i++) 
                {
                  xx = x + reconstructor->_directions[i][0];
                  yy = y + reconstructor->_directions[i][1];
                  zz = z + reconstructor->_directions[i][2];
                  if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
			               if(reconstructor->_confidence_map(xx,yy,zz,t)>0)
                     {
                       val += b[i](x, y, z, t) * original(xx, yy, zz, t);
                       sum += b[i](x, y, z, t);
                     }
                }

                for (int i = 0; i < 13; i++) {
                  xx = x - reconstructor->_directions[i][0];
                  yy = y - reconstructor->_directions[i][1];
                  zz = z - reconstructor->_directions[i][2];
                  if ((xx >= 0) && (xx < dx) && (yy >= 0) && (yy < dy) && (zz >= 0) && (zz < dz)) 
			              if(reconstructor->_confidence_map(xx,yy,zz,t)>0)
			              {
                      val += b[i](x, y, z, t) * original(xx, yy, zz, t);
                      sum += b[i](x, y, z, t);
                    }    
                }

                val -= sum * original(x, y, z, t);
                val = original(x, y, z, t) 
                      + reconstructor->_alpha * reconstructor->_lambda / (reconstructor->_delta * reconstructor->_delta) * val;
                reconstructor->_reconstructed4D(x, y, z, t) = val;
              } 
            } 
      } 
    }

    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, reconstructor->_reconstructed.GetX()),
                      *this );
        init.terminate();
    }

};


// -----------------------------------------------------------------------------
// Adaptive Regularization
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::AdaptiveRegularizationCardiac4D(int iter, irtkRealImage& original)
{
    if (_debug)
          cout << "AdaptiveRegularizationCardiac4D."<< endl;
        //cout << "AdaptiveRegularizationCardiac4D: _delta = "<<_delta<<" _lambda = "<<_lambda <<" _alpha = "<<_alpha<< endl;

    vector<double> factor(13,0);
    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 3; j++)
            factor[i] += fabs(double(_directions[i][j]));
        factor[i] = 1 / factor[i];
    }

    vector<irtkRealImage> b;//(13);
    for (int i = 0; i < 13; i++)
        b.push_back( _reconstructed4D );

    ParallelAdaptiveRegularization1Cardiac4D parallelAdaptiveRegularization1( this,
                                                                     b,
                                                                     factor,
                                                                     original );
    parallelAdaptiveRegularization1();

    irtkRealImage original2 = _reconstructed4D;
    ParallelAdaptiveRegularization2Cardiac4D parallelAdaptiveRegularization2( this,
                                                                     b,
                                                                     factor,
                                                                     original2 );
    parallelAdaptiveRegularization2();

    if (_alpha * _lambda / (_delta * _delta) > 0.068) {
        cerr
            << "Warning: regularization might not have smoothing effect! Ensure that alpha*lambda/delta^2 is below 0.068."
            << endl;
    }
}


// -----------------------------------------------------------------------------
// ReadTransformation
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::ReadTransformation(char* folder)
{
    int n = _slices.size();
    char name[256];
    char path[256];
    irtkTransformation *transformation;
    irtkRigidTransformation *rigidTransf;

    if (n == 0) {
        cerr << "Please create slices before reading transformations!" << endl;
        exit(1);
    }
    cout << "Reading transformations:" << endl;

    _transformations.clear();
    for (int i = 0; i < n; i++) {
        if (folder != NULL) {
            sprintf(name, "/transformation%05i.dof", i);
            strcpy(path, folder);
            strcat(path, name);
        }
        else {
            sprintf(path, "transformation%05i.dof", i);
        }
        transformation = irtkTransformation::New(path);
        rigidTransf = dynamic_cast<irtkRigidTransformation*>(transformation);
        _transformations.push_back(*rigidTransf);
        delete transformation;
        cout << path << endl;
    }
}

// -----------------------------------------------------------------------------
// Calculate Entropy
// -----------------------------------------------------------------------------
double irtkReconstructionCardiac4D::CalculateEntropy()
{
  
  double x;
  double sum_x_sq = 0;
  double x_max = 0;
  double entropy = 0;
  
  for ( int i = 0; i < _reconstructed4D.GetX(); i++)  
    for ( int j = 0; j < _reconstructed4D.GetY(); j++)  
      for ( int k = 0; k < _reconstructed4D.GetZ(); k++)  
        if ( _mask(i,j,k) == 1 )
          for ( int f = 0; f < _reconstructed4D.GetT(); f++)  
            {
              x = _reconstructed4D(i,j,k,f);
              sum_x_sq += x*x;
            }
            
  x_max = sqrt( sum_x_sq );
  
  for ( int i = 0; i < _reconstructed4D.GetX(); i++)  
    for ( int j = 0; j < _reconstructed4D.GetY(); j++)  
      for ( int k = 0; k < _reconstructed4D.GetZ(); k++)  
        if ( _mask(i,j,k) == 1 )
          for ( int f = 0; f < _reconstructed4D.GetT(); f++)  
          {
            x = _reconstructed4D(i,j,k,f);
            if (x>0)
              entropy += x/x_max * log( x/x_max );
          }
  
  entropy = -entropy;
  
  return entropy;
  
}


// -----------------------------------------------------------------------------
// SaveTransformations
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SaveTransformations()
{
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        sprintf(buffer, "transformation%05i.dof", inputIndex);
        _transformations[inputIndex].irtkTransformation::Write(buffer);
    }
}



// -----------------------------------------------------------------------------
// SaveBiasFields
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SaveBiasFields()
{
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        sprintf(buffer, "bias%05i.nii.gz", inputIndex);
        _bias[inputIndex].Write(buffer);
    }
}

void irtkReconstructionCardiac4D::SaveBiasFields( vector<irtkRealImage> &stacks )
{
        
    if (_debug)
        cout << "SaveBiasFields as stacks ...";
    
    char buffer[256];
    irtkRealImage stack;
    vector<irtkRealImage> biasstacks;

    for (unsigned int i = 0; i < stacks.size(); i++) {
        irtkImageAttributes attr = stacks[i].GetImageAttributes();
        stack.Initialize( attr );
        biasstacks.push_back( stack );
    }

    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) { 
        for (int i = 0; i < _slices[inputIndex].GetX(); i++) {
            for (int j = 0; j < _slices[inputIndex].GetY(); j++) {
                biasstacks[_stack_index[inputIndex]](i,j,_stack_loc_index[inputIndex],_stack_dyn_index[inputIndex]) = _bias[inputIndex](i,j,0);
            }
        }
    }

    for (unsigned int i = 0; i < stacks.size(); i++) {
        sprintf(buffer, "bias%03i.nii.gz", i);
        biasstacks[i].Write(buffer);
    }
    
    if (_debug)
        cout << " done." << endl;

}

void irtkReconstructionCardiac4D::SaveBiasFields( vector<irtkRealImage> &stacks, int iter, int rec_iter )
{
        
    if (_debug)
        cout << "SaveBiasFields as stacks ...";
    
    char buffer[256];
    irtkRealImage stack;
    vector<irtkRealImage> biasstacks;

    for (unsigned int i = 0; i < stacks.size(); i++) {
        irtkImageAttributes attr = stacks[i].GetImageAttributes();
        stack.Initialize( attr );
        biasstacks.push_back( stack );
    }

    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) { 
        for (int i = 0; i < _slices[inputIndex].GetX(); i++) {
            for (int j = 0; j < _slices[inputIndex].GetY(); j++) {
                biasstacks[_stack_index[inputIndex]](i,j,_stack_loc_index[inputIndex],_stack_dyn_index[inputIndex]) = _bias[inputIndex](i,j,0);
            }
        }
    }

    for (unsigned int i = 0; i < stacks.size(); i++) {
        sprintf(buffer, "bias%03i_mc%02isr%02i.nii.gz", i, iter, rec_iter);
        biasstacks[i].Write(buffer);
    }
    
    if (_debug)
        cout << " done." << endl;

}


// -----------------------------------------------------------------------------
// SaveSimulatedSlices
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SaveSimulatedSlices()
{
    if (_debug)
        cout<<"Saving simulated images ...";
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
        {
            sprintf(buffer, "simimage%05i.nii.gz", inputIndex);
            _simulated_slices[inputIndex].Write(buffer);
        }
        cout<<"done."<<endl;
}

void irtkReconstructionCardiac4D::SaveSimulatedSlices( vector<irtkRealImage> &stacks )
{      
  
    if (_debug)
        cout << "Saving simulated images as stacks ...";

    char buffer[256];
    irtkRealImage stack;
    vector<irtkRealImage> simstacks;

    for (unsigned int i = 0; i < stacks.size(); i++) {
          irtkImageAttributes attr = stacks[i].GetImageAttributes();
          stack.Initialize( attr );
          simstacks.push_back( stack );
    }

    for (unsigned int inputIndex = 0; inputIndex < _simulated_slices.size(); ++inputIndex) { 
      for (int i = 0; i < _simulated_slices[inputIndex].GetX(); i++) {
          for (int j = 0; j < _simulated_slices[inputIndex].GetY(); j++) {
              simstacks[_stack_index[inputIndex]](i,j,_stack_loc_index[inputIndex],_stack_dyn_index[inputIndex]) = _simulated_slices[inputIndex](i,j,0);
          }
      }
    }

    for (unsigned int i = 0; i < stacks.size(); i++) {
      sprintf(buffer, "simstack%03i.nii.gz", i);
      simstacks[i].Write(buffer);
    }

    if (_debug)
      cout << " done." << endl;

}

void irtkReconstructionCardiac4D::SaveSimulatedSlices( vector<irtkRealImage> &stacks, int iter, int rec_iter )
{      
  
    if (_debug)
        cout << "Saving simulated images as stacks ...";

    char buffer[256];
    irtkRealImage stack;
    vector<irtkRealImage> simstacks;

    for (unsigned int i = 0; i < stacks.size(); i++) {
          irtkImageAttributes attr = stacks[i].GetImageAttributes();
          stack.Initialize( attr );
          simstacks.push_back( stack );
    }

    for (unsigned int inputIndex = 0; inputIndex < _simulated_slices.size(); ++inputIndex) { 
      for (int i = 0; i < _simulated_slices[inputIndex].GetX(); i++) {
          for (int j = 0; j < _simulated_slices[inputIndex].GetY(); j++) {
              simstacks[_stack_index[inputIndex]](i,j,_stack_loc_index[inputIndex],_stack_dyn_index[inputIndex]) = _simulated_slices[inputIndex](i,j,0);
          }
      }
    }

    for (unsigned int i = 0; i < stacks.size(); i++) {
      sprintf(buffer, "simstack%03i_mc%02isr%02i.nii.gz", i, iter, rec_iter);
      simstacks[i].Write(buffer);
    }

    if (_debug)
      cout << " done." << endl;

}


// -----------------------------------------------------------------------------
// SaveSlices
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SaveSlices()
{

    if (_debug)
        cout << "SaveSlices" << endl;
      
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
        {
            sprintf(buffer, "image%05i.nii.gz", inputIndex);
            _slices[inputIndex].Write(buffer);
        }
}

void irtkReconstructionCardiac4D::SaveSlices( vector<irtkRealImage> &stacks )
{
        
    if (_debug)
        cout << "SaveSlices as stacks ...";
    
    char buffer[256];
    irtkRealImage stack;
    vector<irtkRealImage> imagestacks;

    for (unsigned int i = 0; i < stacks.size(); i++) {
        irtkImageAttributes attr = stacks[i].GetImageAttributes();
        stack.Initialize( attr );
        imagestacks.push_back( stack );
    }

    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) { 
        for (int i = 0; i < _slices[inputIndex].GetX(); i++) {
            for (int j = 0; j < _slices[inputIndex].GetY(); j++) {
                imagestacks[_stack_index[inputIndex]](i,j,_stack_loc_index[inputIndex],_stack_dyn_index[inputIndex]) = _slices[inputIndex](i,j,0);
            }
        }
    }

    for (unsigned int i = 0; i < stacks.size(); i++) {
        sprintf(buffer, "stack%03i.nii.gz", i);
        imagestacks[i].Write(buffer);
    }
    
    if (_debug)
        cout << " done." << endl;

}


// -----------------------------------------------------------------------------
// SaveWeights
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SaveWeights()
{
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        sprintf(buffer, "weights%05i.nii.gz", inputIndex);
        _weights[inputIndex].Write(buffer);
    }
}

void irtkReconstructionCardiac4D::SaveWeights( vector<irtkRealImage> &stacks )
{
        
    if (_debug)
        cout << "SaveWeights as stacks ...";
    
    char buffer[256];
    irtkRealImage stack;
    vector<irtkRealImage> weightstacks;

    for (unsigned int i = 0; i < stacks.size(); i++) {
        irtkImageAttributes attr = stacks[i].GetImageAttributes();
        stack.Initialize( attr );
        weightstacks.push_back( stack );
    }

    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) { 
        for (int i = 0; i < _slices[inputIndex].GetX(); i++) {
            for (int j = 0; j < _slices[inputIndex].GetY(); j++) {
                weightstacks[_stack_index[inputIndex]](i,j,_stack_loc_index[inputIndex],_stack_dyn_index[inputIndex]) = _weights[inputIndex](i,j,0);
            }
        }
    }

    for (unsigned int i = 0; i < stacks.size(); i++) {
        sprintf(buffer, "weight%03i.nii.gz", i);
        weightstacks[i].Write(buffer);
    }
    
    if (_debug)
        cout << " done." << endl;

}

void irtkReconstructionCardiac4D::SaveWeights( vector<irtkRealImage> &stacks, int iter, int rec_iter )
{
        
    if (_debug)
        cout << "SaveWeights as stacks ...";
    
    char buffer[256];
    irtkRealImage stack;
    vector<irtkRealImage> weightstacks;

    for (unsigned int i = 0; i < stacks.size(); i++) {
        irtkImageAttributes attr = stacks[i].GetImageAttributes();
        stack.Initialize( attr );
        weightstacks.push_back( stack );
    }

    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) { 
        for (int i = 0; i < _slices[inputIndex].GetX(); i++) {
            for (int j = 0; j < _slices[inputIndex].GetY(); j++) {
                weightstacks[_stack_index[inputIndex]](i,j,_stack_loc_index[inputIndex],_stack_dyn_index[inputIndex]) = _weights[inputIndex](i,j,0);
            }
        }
    }

    for (unsigned int i = 0; i < stacks.size(); i++) {
        sprintf(buffer, "weights%03i_mc%02isr%02i.nii.gz", i, iter, rec_iter);
        weightstacks[i].Write(buffer);
    }
    
    if (_debug)
        cout << " done." << endl;

}


// -----------------------------------------------------------------------------
// SlicesInfo
// -----------------------------------------------------------------------------
void irtkReconstructionCardiac4D::SlicesInfoCardiac4D( const char* filename,
                                    vector<string> &stack_files )
{
    ofstream info;
    info.open( filename );

    info<<setprecision(3);

    // header
    info << "StackIndex" << "\t"
         << "StackLocIndex" << "\t"
         << "StackDynIndex" << "\t"
         << "LocIndex" << "\t"
         << "InputIndex" << "\t"
         << "File" << "\t"
         << "Scale" << "\t"
         << "StackFactor" << "\t"
         << "Time" << "\t"
         << "TemporalResolution" << "\t"
         << "CardiacPhase" << "\t"
         << "ReconCardPhaseIndex" << "\t"
         << "Included" << "\t" // Included slices
         << "Excluded" << "\t"  // Excluded slices
         << "Outside" << "\t"  // Outside slices
         << "Weight" << "\t"
         << "TranslationX" << "\t"
         << "TranslationY" << "\t"
         << "TranslationZ" << "\t"
         << "RotationX" << "\t"
         << "RotationY" << "\t"
         << "RotationZ" << "\t";
         for (int j = 0; j < _reconstructed4D.GetT(); j++)
            info << "TemporalWeightReconCardPhaseIndex" << j << "\t";
         info << "\b" << endl;       
    
    for (unsigned int i = 0; i < _slices.size(); i++) {
        irtkRigidTransformation& t = _transformations[i];
        info << _stack_index[i] << "\t"
             << _stack_loc_index[i] << "\t" 
             << _stack_dyn_index[i] << "\t" 
             << _loc_index[i] << "\t"
             << i << "\t"
             << stack_files[_stack_index[i]] << "\t"
             << _scale[i] << "\t"
             << _stack_factor[_stack_index[i]] << "\t"
             << _slice_time[i] << "\t"
             << _slice_dt[i] << "\t"
             << _slice_cardphase[i] << "\t"
             << _slice_svr_card_index[i] << "\t"
             << (((_slice_weight[i] >= 0.5) && (_slice_inside[i]))?1:0) << "\t" // Included slices
             << (((_slice_weight[i] < 0.5) && (_slice_inside[i]))?1:0) << "\t"  // Excluded slices
             << ((!(_slice_inside[i]))?1:0) << "\t"  // Outside slices
             << _slice_weight[i] << "\t"
             << t.GetTranslationX() << "\t"
             << t.GetTranslationY() << "\t"
             << t.GetTranslationZ() << "\t"
             << t.GetRotationX() << "\t"
             << t.GetRotationY() << "\t"
             << t.GetRotationZ() << "\t";
             for (int j = 0; j < _reconstructed4D.GetT(); j++)
                info << _slice_temporal_weight[j][i] << "\t";
             info << "\b" << endl;
    }
 
    info.close(); 
}
