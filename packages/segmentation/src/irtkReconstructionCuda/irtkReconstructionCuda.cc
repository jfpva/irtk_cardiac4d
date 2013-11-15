/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

#include <irtkReconstructionCuda.h>
#include <irtkResampling.h>
#include <irtkRegistration.h>
#include <irtkImageRigidRegistration.h>
#include <irtkImageRigidRegistrationWithPadding.h>
#include <irtkTransformation.h>
#include <irtkMeanShift.h>
#include <irtkCRF.h>
#include "perfstats.h"

/* Auxiliary functions (not reconstruction specific) */

void bbox_( irtkRealImage &stack,
           irtkRigidTransformation &transformation,
           double &min_x,
           double &min_y,
           double &min_z,
           double &max_x,
           double &max_y,
           double &max_z ) {

    cout << "bbox_" << endl;
    
    min_x = voxel_limits<irtkRealPixel>::max();
    min_y = voxel_limits<irtkRealPixel>::max();
    min_z = voxel_limits<irtkRealPixel>::max();
    max_x = voxel_limits<irtkRealPixel>::min();
    max_y = voxel_limits<irtkRealPixel>::min();
    max_z = voxel_limits<irtkRealPixel>::min();
    double x,y,z;
    for ( int i = 0; i <= stack.GetX(); i += stack.GetX() )
        for ( int j = 0; j <= stack.GetY(); j += stack.GetY() )
            for ( int k = 0; k <= stack.GetZ(); k += stack.GetZ() ) {
                x = i;
                y = j;
                z = k;
                stack.ImageToWorld( x, y, z );
                // FIXME!!!
                transformation.Transform( x, y, z );
                //transformation.Inverse( x, y, z );
                if ( x < min_x )
                    min_x = x;
                if ( y < min_y )
                    min_y = y;
                if ( z < min_z )
                    min_z = z;
                if ( x > max_x )
                    max_x = x;
                if ( y > max_y )
                    max_y = y;
                if ( z > max_z )
                    max_z = z;
            }
}

void bboxCrop_( irtkRealImage &image ) {
    int min_x, min_y, min_z, max_x, max_y, max_z;
    min_x = image.GetX()-1;
    min_y = image.GetY()-1;
    min_z = image.GetZ()-1;
    max_x = 0;
    max_y = 0;
    max_z = 0;
    for ( int i = 0; i < image.GetX(); i++ )
        for ( int j = 0; j < image.GetY(); j++ )
            for ( int k = 0; k < image.GetZ(); k++ ) {
                if ( image.Get(i, j, k) > 0 ) {
                    if ( i < min_x )
                        min_x = i;
                    if ( j < min_y )
                        min_y = j;
                    if ( k < min_z )
                        min_z = k;
                    if ( i > max_x )
                        max_x = i;
                    if ( j > max_y )
                        max_y = j;
                    if ( k > max_z )
                        max_z = k;                    
                }
            }

    //Cut region of interest
    image = image.GetRegion( min_x, min_y, min_z,
                             max_x+1, max_y+1, max_z+1 );                
}

void centroid_( irtkRealImage &image,
               double &x,
               double &y,
               double &z ) {
    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;
    double norm = 0;
    double v;
    for ( int i = 0; i < image.GetX(); i++ )
        for ( int j = 0; j < image.GetY(); j++ )
            for ( int k = 0; k < image.GetZ(); k++ ) {
                v = image.Get(i, j, k);
                if ( v <= 0 )
                    continue;
                sum_x += v*i;
                sum_y += v*j;
                sum_z += v*k;
                norm += v;
            }

    x = sum_x / norm;
    y = sum_y / norm;
    z = sum_z / norm;
            
    image.ImageToWorld( x, y, z );

    std::cout << "CENTROID:" << x << "," << y <<"," <<z<<"\n\n";
}

/*   end of auxiliary functions */

irtkReconstructionCuda::irtkReconstructionCuda()
{
    _step = 0.0001;
    _debug = false;
    _quality_factor = 2;
    _sigma_bias = 12;
    _sigma_s = 0.025;
    _sigma_s2 = 0.025;
    _mix_s = 0.9;
    _mix = 0.9;
    _delta = 1;
    _lambda = 0.1;
    _alpha = (0.05 / _lambda) * _delta * _delta;
    _template_created = false;
    _have_mask = false;
    _low_intensity_cutoff = 0.01;
    _global_bias_correction = false;
    _adaptive = false;

	//TODO device selection 
	cudaSetDevice(0);

    int directions[13][3] = {
        { 1, 0, -1 },
        { 0, 1, -1 },
        { 1, 1, -1 },
        { 1, -1, -1 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 1, 1, 0 },
        { 1, -1, 0 },
        { 1, 0, 1 },
        { 0, 1, 1 },
        { 1, 1, 1 },
        { 1, -1, 1 },
        { 0, 0, 1 }
    };
    for ( int i = 0; i < 13; i++ )
        for ( int j = 0; j < 3; j++ )
		{
            _directions[i][j] = directions[i][j];
		}


}

irtkReconstructionCuda::~irtkReconstructionCuda() { }

void irtkReconstructionCuda::CenterStacks( vector<irtkRealImage>& stacks,
                                       vector<irtkRigidTransformation>& stack_transformations,
                                       int templateNumber ) {
    // template center
    double x0, y0, z0;
    irtkRealImage mask;
    mask = stacks[templateNumber] != -1;
    centroid_( mask, x0, y0, z0 );

    double x, y, z;
    irtkMatrix m1, m2;
    for ( unsigned int i = 0; i < stacks.size(); i++ ) {
        if (i == templateNumber)
            continue;

        mask = stacks[i] != -1;
        centroid_( mask, x, y, z );

        irtkRigidTransformation translation;
        translation.PutTranslationX( x0-x );
        translation.PutTranslationY( y0-y );
        translation.PutTranslationZ( z0-z );

        std::cout << "TRANSLATION:\n";
        translation.Print();
        std::cout << "\n\n\n";
        
        m1 = stack_transformations[i].GetMatrix();
        m2 = translation.GetMatrix();
        stack_transformations[i].PutMatrix( m2*m1 );
    }

}

//TODO could be cudarized via linear texture lookup but is only executed once in the beginning
class ParallelAverage{
    irtkReconstructionCuda* reconstructor;
    vector<irtkRealImage> &stacks;
    vector<irtkRigidTransformation> &stack_transformations;

    /// Padding value in target (voxels in the target image with this
    /// value will be ignored)
    double targetPadding;

    /// Padding value in source (voxels outside the source image will
    /// be set to this value)
    double sourcePadding;

    double background;

    // Volumetric registrations are stack-to-template while slice-to-volume 
    // registrations are actually performed as volume-to-slice
    // (reasons: technicalities of implementation)
    // so transformations need to be inverted beforehand.

    bool linear;
    
public:
    irtkRealImage average;
    irtkRealImage weights;
    
    void operator()( const blocked_range<size_t>& r ) {

		PerfStats stats;
		stats.start();

        for ( size_t i = r.begin(); i < r.end(); ++i) {
            irtkImageTransformation imagetransformation;
            irtkImageFunction *interpolator;
            if (linear)
                interpolator = new irtkLinearInterpolateImageFunction;
            else
                interpolator = new irtkNearestNeighborInterpolateImageFunction;
              
            irtkRealImage s = stacks[i];
            irtkRigidTransformation t = stack_transformations[i];
            imagetransformation.SetInput(&s, &t);
            irtkRealImage image( reconstructor->_reconstructed.GetImageAttributes() );
            image = 0;
               
            imagetransformation.SetOutput(&image);
            imagetransformation.PutTargetPaddingValue(targetPadding);
            imagetransformation.PutSourcePaddingValue(sourcePadding);
            imagetransformation.PutInterpolator(interpolator);
            imagetransformation.Run();
               
            irtkRealPixel *pa = average.GetPointerToVoxels();
            irtkRealPixel *pi = image.GetPointerToVoxels();
            irtkRealPixel *pw = weights.GetPointerToVoxels();
            for (int p = 0; p < average.GetNumberOfVoxels(); p++) {
                if (*pi != background) {
                    *pa += *pi;
                    *pw += 1;
                }
                pa++;
                pi++;
                pw++;
            }
			stats.sample("ParallelAverage ()");
			stats.print();
            delete interpolator;
        }
    }
 
    ParallelAverage( ParallelAverage& x, split ) :
        reconstructor(x.reconstructor),
        stacks(x.stacks),
        stack_transformations(x.stack_transformations)
    {
        average.Initialize( reconstructor->_reconstructed.GetImageAttributes() );
        average = 0;
        weights.Initialize( reconstructor->_reconstructed.GetImageAttributes() );
        weights = 0;
        targetPadding = x.targetPadding;
        sourcePadding = x.sourcePadding;
        background = x.background;
        linear = x.linear;
    }
 
    void join( const ParallelAverage& y ) {        
        average += y.average;
        weights += y.weights;        
    }
             
    ParallelAverage( irtkReconstructionCuda *reconstructor,
                     vector<irtkRealImage>& _stacks,
                     vector<irtkRigidTransformation>& _stack_transformations,
                     double _targetPadding,
                     double _sourcePadding,
                     double _background,
                     bool _linear=false ) :
        reconstructor(reconstructor),
        stacks(_stacks),
        stack_transformations(_stack_transformations)
    {
        average.Initialize( reconstructor->_reconstructed.GetImageAttributes() );
        average = 0;
        weights.Initialize( reconstructor->_reconstructed.GetImageAttributes() );
        weights = 0;
        targetPadding = _targetPadding;
        sourcePadding = _sourcePadding;
        background = _background;
        linear = _linear;
    }

    // execute
    void operator() () {
        task_scheduler_init init(tbb_no_threads);
        parallel_reduce( blocked_range<size_t>(0,stacks.size()),
                         *this );
        init.terminate();
    }    
};

class ParallelSliceAverage{
    irtkReconstructionCuda* reconstructor;
    vector<irtkRealImage> &slices;
    vector<irtkRigidTransformation> &slice_transformations;
    irtkRealImage &average;
    irtkRealImage &weights;
    
public:
    
    void operator()( const blocked_range<size_t>& r ) const {

		PerfStats stats;
		stats.start();

        for ( int k0 = r.begin(); k0 < r.end(); k0++) {
            for ( int j0 = 0; j0 < average.GetY(); j0++) {
                for ( int i0 = 0; i0 < average.GetX(); i0++) {
                        double x0 = i0;
                        double y0 = j0;
                        double z0 = k0;
                        // Transform point into world coordinates
                        average.ImageToWorld(x0, y0, z0);
                        for (int inputIndex = 0; inputIndex < slices.size(); inputIndex++ ) {
                            double x = x0;
                            double y = y0;
                            double z = z0;
                            // Transform point
                            slice_transformations[inputIndex].Transform(x, y, z);
                            // Transform point into image coordinates
                            slices[inputIndex].WorldToImage(x, y, z);
                            int i = round(x);
                            int j = round(y);
                            int k = round(z);
                            // Check whether transformed point is in FOV of input
                            if ( (i >=  0) && (i < slices[inputIndex].GetX()) &&
                                 (j >= 0) && (j <  slices[inputIndex].GetY()) &&
                                 (k >= 0) && (k < slices[inputIndex].GetZ()) ) {
                                if (slices[inputIndex](i,j,k) > 0) {
                                    average.PutAsDouble( i0, j0, k0, 0, average( i0, j0, k0 ) + slices[inputIndex](i,j,k) );
                                    weights.PutAsDouble( i0, j0, k0, 0, weights( i0, j0, k0 ) + 1 );
                                }
                            }
                        }
                }
            }
        }

		stats.sample("ParallelSliceAverage ()");
		stats.print();

	}
              
    ParallelSliceAverage( irtkReconstructionCuda *reconstructor,
                          vector<irtkRealImage>& _slices,
                          vector<irtkRigidTransformation>& _slice_transformations,
                          irtkRealImage &_average,
                          irtkRealImage &_weights ) :
        reconstructor(reconstructor),
        slices(_slices),
        slice_transformations(_slice_transformations),
        average(_average),
        weights(_weights)
    {
        average.Initialize( reconstructor->_reconstructed.GetImageAttributes() );
        average = 0;
        weights.Initialize( reconstructor->_reconstructed.GetImageAttributes() );
        weights = 0;
    }

    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0,average.GetZ()),
                      *this );
        init.terminate();
    }    
};

irtkRealImage irtkReconstructionCuda::CreateAverage(vector<irtkRealImage>& stacks,
                                                vector<irtkRigidTransformation>& stack_transformations)
{
    if (!_template_created) {
        cerr << "Please create the template before calculating the average of the stacks." << endl;
        exit(1);
    }
    
    InvertStackTransformations(stack_transformations);
    ParallelAverage parallelAverage( this,
                                     stacks,
                                     stack_transformations,
                                     -1,0,0, // target/source/background
                                     true ); 
    parallelAverage();
    irtkRealImage average = parallelAverage.average;
    irtkRealImage weights = parallelAverage.weights;
    average /= weights;
    InvertStackTransformations(stack_transformations);
    return average;
}

double irtkReconstructionCuda::CreateTemplate(irtkRealImage stack, double resolution)
{
    double dx, dy, dz, d;

    //Get image attributes - image size and voxel size
    irtkImageAttributes attr = stack.GetImageAttributes();

    //enlarge stack in z-direction in case top of the head is cut off
    attr._z += 2;

    //create enlarged image
    irtkRealImage enlarged(attr);

    //determine resolution of volume to reconstruct
    if (resolution <= 0) {
        //resolution was not given by user
        // set it to min of res in x or y direction
        stack.GetPixelSize(&dx, &dy, &dz);
        if ((dx <= dy) && (dx <= dz))
            d = dx;
        else if (dy <= dz)
            d = dy;
        else
            d = dz;
    }
    else
        d = resolution;

    cout << "Constructing volume with isotropic voxel size " << d << endl;

    //resample "enlarged" to resolution "d"
    irtkNearestNeighborInterpolateImageFunction interpolator;
    irtkResampling<irtkRealPixel> resampling(d, d, d);
    resampling.SetInput(&enlarged);
    resampling.SetOutput(&enlarged);
    resampling.SetInterpolator(&interpolator);
    resampling.Run();

    //initialize recontructed volume
    _reconstructed = enlarged;
    _template_created = true;

    return d;
}

void irtkReconstructionCuda::CreateMaskFromBlackBackground( vector<irtkRealImage>& stacks,
                                                        vector<irtkRigidTransformation>& stack_transformations,
                                                        double smooth_mask )
{
	//Create average of the stack using currect stack transformations
	irtkGreyImage average = CreateAverage(stacks, stack_transformations);

    irtkGreyPixel* ptr = average.GetPointerToVoxels();
    for (int i = 0; i < average.GetNumberOfVoxels(); i++)
        if (*ptr < 0)
            *ptr = 0;
    
	//Create mask of the average from the black background
	irtkMeanShift msh(average, 0, 256);
	msh.GenerateDensity();
	msh.SetTreshold();
	msh.RemoveBackground();
	irtkGreyImage mask = msh.ReturnMask();

	//Calculate LCC of the mask to remove disconected structures
	irtkMeanShift msh2(mask, 0, 256);
	msh2.SetOutput(&mask);
	msh2.Lcc(1);
	irtkRealImage m = mask;
	SetMask(&m, smooth_mask);
}

double irtkReconstructionCuda::CreateLargeTemplate( vector<irtkRealImage>& stacks,
                                                vector<irtkRigidTransformation>& stack_transformations,
                                                double resolution,
                                                double smooth_mask,
                                                double threshold_mask,
                                                double expand )
{
    if (_debug)
        cout << "CreateLargeTemplate" << endl;
    
    //double min_x, min_y, min_z, max_x, max_y, max_z;
    double min_x = voxel_limits<irtkRealPixel>::max();
    double min_y = voxel_limits<irtkRealPixel>::max();
    double min_z = voxel_limits<irtkRealPixel>::max();
    double max_x = voxel_limits<irtkRealPixel>::min();
    double max_y = voxel_limits<irtkRealPixel>::min();
    double max_z = voxel_limits<irtkRealPixel>::min();
    double min_x0, min_y0, min_z0, max_x0, max_y0, max_z0;
    for (unsigned int i = 0; i < stacks.size(); i++) {
        if (_slice_weight.size() > 0)
            if ( ! ( (_slice_weight[i] >= 0.5) && (_slice_inside[i]) ) )
                continue;
        bbox_( stacks[i],
              stack_transformations[i],
              min_x0, min_y0, min_z0, max_x0, max_y0, max_z0 );
        if ( min_x0 < min_x )
            min_x = min_x0;
        if ( min_y0 < min_y )
            min_y = min_y0;
        if ( min_z0 < min_z )
            min_z = min_z0;
        if ( max_x0 > max_x )
            max_x = max_x0;
        if ( max_y0 > max_y )
            max_y = max_y0;
        if ( max_z0 > max_z )
            max_z = max_z0;
    }
    
    // Set image attributes
    irtkImageAttributes attr;

    // dim
    attr._x = (max_x+1 - min_x + expand)/resolution;
    attr._y = (max_y+1 - min_y + expand)/resolution;
    attr._z = (max_z+1 - min_z + expand)/resolution;
    attr._t = 1;

    // voxel size
    attr._dx = resolution;
    attr._dy = resolution;
    attr._dz = resolution;
    attr._dt = 1;
    
    // origin
    attr._xorigin = min_x + (max_x+1 - min_x)/2;
    attr._yorigin = min_y + (max_y+1 - min_y)/2;
    attr._zorigin = min_z + (max_z+1 - min_z)/2;
    attr._torigin = 0;

    // x-axis
    attr._xaxis[0] = 1;
    attr._xaxis[1] = 0;
    attr._xaxis[2] = 0;

    // y-axis
    attr._yaxis[0] = 0;
    attr._yaxis[1] = 1;
    attr._yaxis[2] = 0;

    // z-axis
    attr._zaxis[0] = 0;
    attr._zaxis[1] = 0;
    attr._zaxis[2] = 1;

    _reconstructed.Initialize(attr);
    _template_created = true;

    CreateMaskFromAllMasks( stacks,
                            stack_transformations,
                            smooth_mask, threshold_mask );
    
    bboxCrop_( _mask );
    _reconstructed = _mask;

    _reconstructed.Write("_reconstructed.nii");

	reconstructionGPU.InitReconstructionVolume(make_uint3(_reconstructed.GetX(), _reconstructed.GetY(), _reconstructed.GetZ()), 
		make_float3(_reconstructed.GetXSize(), _reconstructed.GetYSize(), _reconstructed.GetZSize()), _reconstructed.GetPointerToVoxels(),  _sigma_bias/_reconstructed.GetXSize() );
    
    //return resulting resolution of the template image
    return resolution;
}

void irtkReconstructionCuda::CreateMaskFromAllMasks( vector<irtkRealImage>& stacks,
                                                 vector<irtkRigidTransformation>& stack_transformations,
                                                 double smooth_mask,
                                                 double threshold_mask)
{
    if (_debug)
        cout << "CreateMaskFromAllMasks" << endl;

    // FIXME: use a faster way when summing all slices (for each voxel of final
    // volume, go through each slice)
    InvertStackTransformations(stack_transformations);
    irtkRealImage weights;

    if (stacks[0].GetZ() > 1) {
        ParallelAverage parallelAverage( this,
                                         stacks,
                                         stack_transformations,
                                         voxel_limits<irtkRealPixel>::min(), -1, -1 ); // target/source/background
        parallelAverage();
        weights = parallelAverage.weights;
    }
    else {
        irtkRealImage average;
        ParallelSliceAverage parallelSliceAverage( this,
                                                   stacks,
                                                   stack_transformations,
                                                   average,
                                                   weights ); 
        parallelSliceAverage();
    }
    InvertStackTransformations(stack_transformations);
    
    irtkGaussianBlurring<irtkRealPixel> gb(smooth_mask);
    gb.SetInput(&weights);
    gb.SetOutput(&weights);
    gb.Run();

    weights.PutMinMax(0,1);
    weights.Write("weights4mask.nii");
    // FIXME threshold_mask ???
    irtkGreyImage mask = (weights < threshold_mask) != threshold_mask;
    //Calculate LCC of the mask to remove disconected structures
    irtkMeanShift msh2(mask, 0, 256);
    msh2.SetOutput(&mask);
    msh2.Lcc(1);
    irtkRealImage m = mask;
    SetMask(&m, 0);//smooth_mask
    m.Write("mymask.nii");

}

void irtkReconstructionCuda::UpdateMaskFromAllMasks( double smooth_mask,
                                                 double threshold_mask )
{
    if (_debug)
        cout << "UpdateMaskFromBlackBackground" << endl;
    CreateMaskFromAllMasks( _slices, _transformations, smooth_mask, threshold_mask );
}

void irtkReconstructionCuda::UpdateProbabilityMap()
{
	PerfStats stats;
	stats.start();

    if (_debug)
        cout << "UpdateProbabilityMap" << endl;

    _brain_probability.Initialize( _reconstructed.GetImageAttributes() );
    _brain_probability = 0;
    
    InvertStackTransformations(_transformations);
    // ParallelAverage parallelAverage( this,
    //                                  _probability_maps,
    //                                  _transformations,
    //                                  voxel_limits<irtkRealPixel>::min(), 0, 0 ); // target/source/background
    // parallelAverage();

    irtkRealImage average;
    irtkRealImage weights;
    ParallelSliceAverage parallelSliceAverage( this,
                                               _probability_maps,
                                               _transformations,
                                               average,
                                               weights ); 
    parallelSliceAverage();
    InvertStackTransformations(_transformations);
    //    _brain_probability = parallelAverage.average;
    //irtkRealImage weights = parallelAverage.weights;
    _brain_probability = average / weights;
    
    irtkRealPixel *pr = _brain_probability.GetPointerToVoxels();
    irtkRealPixel *pm = _mask.GetPointerToVoxels();
    for (int i = 0; i < _brain_probability.GetNumberOfVoxels(); i++) {
        if (*pm == 0)
            *pr = 0;
        pm++;
        pr++;
    }

	stats.sample("UpdateProbabilityMap()");
	stats.print();
}

void irtkReconstructionCuda::crf3DMask( double smooth_mask,
                                    double threshold_mask ) {
    if (_debug)
        cout << "crf3DMask" << endl;

    // FIXME: need bbox__crop of the refined mask...
    
    _mask.Write("maskbeforeCRF.nii");
    
    irtkGreyImage labels( _reconstructed.GetImageAttributes() );
    irtkGreyPixel *ptr = labels.GetPointerToVoxels();
    irtkRealPixel *pm = _brain_probability.GetPointerToVoxels();
    for ( unsigned int i = 0; i < labels.GetNumberOfVoxels(); i++ ) {
        if (*pm>=0.5)
            *ptr = 2;
        // else if (*pm>=0)
        //     *ptr = 1;
        else
            *ptr = 0;
    }

    cout << "starting crf...\n";
    irtkCRF crf( _reconstructed, labels, _brain_probability );
    crf.Run();
    // crf( _reconstructed.GetPointerToVoxels(),
    //      _reconstructed.GetZ(), _reconstructed.GetY(), _reconstructed.GetX(),
    //      _reconstructed.GetSD(1),
    //      labels.GetPointerToVoxels(),
    //      _brain_probability.GetPointerToVoxels(),
    //      1.0 );

    labels.Write("maskCRF.nii");

    //Calculate LCC of the mask to remove disconected structures
    irtkMeanShift msh2(labels, 0, 256);
    msh2.SetOutput(&labels);
    msh2.Lcc(1);
    irtkRealImage m = labels;
    SetMask(&m, smooth_mask, threshold_mask );//smooth_mask

    bboxCrop_( _mask );
    _reconstructed = _mask;

    _mask.Write("maskafterCRF.nii");
}

irtkRealImage irtkReconstructionCuda::CreateMask(irtkRealImage image)
{
    //binarize mask
    irtkRealPixel* ptr = image.GetPointerToVoxels();
    for (int i = 0; i < image.GetNumberOfVoxels(); i++) {
        if (*ptr > 0.5)
            *ptr = 1;
        else
            *ptr = 0;
        ptr++;
    }
    return image;
}

void irtkReconstructionCuda::SetMask(irtkRealImage * mask, double sigma, double threshold)
{
    if (!_template_created) {
        cerr
            << "Please create the template before setting the mask, so that the mask can be resampled to the correct dimensions."
            << endl;
        exit(1);
    }

    _mask = _reconstructed;

    if (mask != NULL) {
        //if sigma is nonzero first smooth the mask
        if (sigma > 0) {
            //blur mask
            irtkGaussianBlurring<irtkRealPixel> gb(sigma);
            gb.SetInput(mask);
            gb.SetOutput(mask);
            gb.Run();

            //binarize mask
            irtkRealPixel* ptr = mask->GetPointerToVoxels();
            for (int i = 0; i < mask->GetNumberOfVoxels(); i++) {
                if (*ptr > threshold)
                    *ptr = 1;
                else
                    *ptr = 0;
                ptr++;
            }
        }

        //resample the mask according to the template volume using identity transformation
        irtkRigidTransformation transformation;
        irtkImageTransformation imagetransformation;
        irtkNearestNeighborInterpolateImageFunction interpolator;
        imagetransformation.SetInput(mask, &transformation);
        imagetransformation.SetOutput(&_mask);
        //target is zero image, need padding -1
        imagetransformation.PutTargetPaddingValue(-1);
        //need to fill voxels in target where there is no info from source with zeroes
        imagetransformation.PutSourcePaddingValue(0);
        imagetransformation.PutInterpolator(&interpolator);
        imagetransformation.Run();
    }
    else {
        //fill the mask with ones
        _mask = 1;
    }
    //set flag that mask was created
    _have_mask = true;

    if (_debug)
        _mask.Write("mask.nii.gz");
}

void irtkReconstructionCuda::TransformMask(irtkRealImage& image, irtkRealImage& mask,
                                       irtkRigidTransformation& transformation)
{
    //transform mask to the space of image
    irtkImageTransformation imagetransformation;
    irtkNearestNeighborInterpolateImageFunction interpolator;
    imagetransformation.SetInput(&mask, &transformation);
    irtkRealImage m = image;
    imagetransformation.SetOutput(&m);
    //target contains zeros and ones image, need padding -1
    imagetransformation.PutTargetPaddingValue(-1);
    //need to fill voxels in target where there is no info from source with zeroes
    imagetransformation.PutSourcePaddingValue(0);
    imagetransformation.PutInterpolator(&interpolator);
    imagetransformation.Run();
    mask = m;
}

void irtkReconstructionCuda::ResetOrigin(irtkGreyImage &image, irtkRigidTransformation& transformation)
{
    double ox,oy,oz;
    image.GetOrigin(ox,oy,oz);
    image.PutOrigin(0,0,0);
    transformation.PutTranslationX(ox);
    transformation.PutTranslationY(oy);
    transformation.PutTranslationZ(oz);
    transformation.PutRotationX(0);
    transformation.PutRotationY(0);
    transformation.PutRotationZ(0);
}

class ParallelStackRegistrations {
    irtkReconstructionCuda *reconstructor;
    vector<irtkRealImage>& stacks;
    vector<irtkRigidTransformation>& stack_transformations;
    int templateNumber;
    irtkGreyImage& target;
    irtkRigidTransformation& offset;
        
public:
    ParallelStackRegistrations( irtkReconstructionCuda *_reconstructor,
                                vector<irtkRealImage>& _stacks,
                                vector<irtkRigidTransformation>& _stack_transformations,
                                int _templateNumber,
                                irtkGreyImage& _target,
                                irtkRigidTransformation& _offset ) : 
        reconstructor(_reconstructor),
        stacks(_stacks),
        stack_transformations(_stack_transformations),
        target(_target),
        offset(_offset) {
        templateNumber = _templateNumber;
    }

    void operator() (const blocked_range<size_t> &r) const {
		PerfStats stats;
		stats.start();

        for ( size_t i = r.begin(); i != r.end(); ++i ) {            

            //do not perform registration for template
            if (i == templateNumber)
                continue;

            //rigid registration object
            irtkImageRigidRegistrationWithPadding registration;
            //irtkRigidTransformation transformation = stack_transformations[i];
            
            //set target and source (need to be converted to irtkGreyImage)
            irtkGreyImage source = stacks[i];

            //include offset in trasformation   
            irtkMatrix mo = offset.GetMatrix();
            irtkMatrix m = stack_transformations[i].GetMatrix();
            m=m*mo;
            stack_transformations[i].PutMatrix(m);

            //perform rigid registration
            registration.SetInput(&target, &source);
            registration.SetOutput(&stack_transformations[i]);
            registration.GuessParameterThickSlices();
            registration.SetTargetPadding(0);
            registration.Run();
            
            mo.Invert();
            m = stack_transformations[i].GetMatrix();
            m=m*mo;
            stack_transformations[i].PutMatrix(m);

            //stack_transformations[i] = transformation;            

            //save volumetric registrations
            if (reconstructor->_debug) {
                //buffer to create the name
                char buffer[256];
                registration.irtkImageRegistration::Write((char *) "parout-volume.rreg");
                sprintf(buffer, "stack-transformation%i.dof.gz", i);
                stack_transformations[i].irtkTransformation::Write(buffer);
                target.Write("target.nii.gz");
                sprintf(buffer, "stack%i.nii.gz", i);
                stacks[i].Write(buffer);
            }            
        }
		stats.sample("ParallelStackRegistrations");
		stats.print();

    }

    // execute
    void operator() () const {
		PerfStats stats;
		stats.start();
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, stacks.size() ),
                      *this );
        init.terminate();
		stats.sample("ParallelStackRegistrations");
		stats.print();
    }

};

void irtkReconstructionCuda::StackRegistrations(vector<irtkRealImage>& stacks,
                                            vector<irtkRigidTransformation>& stack_transformations, int templateNumber)
{
    if (_debug)
        cout << "StackRegistrations" << endl;

	PerfStats stats;
	stats.start();
    
    InvertStackTransformations(stack_transformations);
    
    //template is set as the target
    irtkGreyImage target = stacks[templateNumber];
    //target needs to be masked before registration
    if (_have_mask) {
        double x, y, z;
        for (int i = 0; i < target.GetX(); i++)
            for (int j = 0; j < target.GetY(); j++)
                for (int k = 0; k < target.GetZ(); k++) {
                    //image coordinates of the target
                    x = i;
                    y = j;
                    z = k;
                    //change to world coordinates
                    target.ImageToWorld(x, y, z);
                    //change to mask image coordinates - mask is aligned with target
                    _mask.WorldToImage(x, y, z);
                    x = round(x);
                    y = round(y);
                    z = round(z);
                    //if the voxel is outside mask ROI set it to -1 (padding value)
                    if ((x >= 0) && (x < _mask.GetX()) && (y >= 0) && (y < _mask.GetY()) && (z >= 0)
                        && (z < _mask.GetZ())) {
                        if (_mask(x, y, z) == 0)
                            target(i, j, k) = 0;
                    }
                    else
                        target(i, j, k) = 0;
                }
    }

    irtkRigidTransformation offset;
    ResetOrigin(target,offset);

    //register all stacks to the target
    ParallelStackRegistrations registration( this,
                                             stacks,
                                             stack_transformations,
                                             templateNumber,
                                             target,
                                             offset );
    registration();
        
    InvertStackTransformations(stack_transformations);

	stats.sample("StackRegistrations");
	stats.print();
}

void irtkReconstructionCuda::RestoreSliceIntensities()
{
    if (_debug)
        cout << "Restoring the intensities of the slices. "<<endl;
    
    unsigned int inputIndex;
    int i;
    double factor;
    irtkRealPixel *p;
  
    for (inputIndex=0;inputIndex<_slices.size();inputIndex++) {
        //calculate scaling factor
        factor = _stack_factor[_stack_index[inputIndex]];//_average_value;
    
        // read the pointer to current slice
        p = _slices[inputIndex].GetPointerToVoxels();
        for(i=0;i<_slices[inputIndex].GetNumberOfVoxels();i++) {
            if(*p>0) *p = *p / factor;
            p++;
        } 
    }
}

void irtkReconstructionCuda::ScaleVolume()
{
    if (_debug)
        cout << "Scaling volume: ";

	//printf("");
    
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
  
    irtkRealPixel *ptr = _reconstructed.GetPointerToVoxels();
    for(i=0;i<_reconstructed.GetNumberOfVoxels();i++) {
        if(*ptr>0) *ptr = *ptr * scale;
        ptr++;
    }
    cout<<endl;
}


void irtkReconstructionCuda::SimulateSlices()
{
	if(_slice_inside.size() == 0)
	{
		_slice_inside.clear();
		_slice_inside.resize(_slices.size());
	}
	reconstructionGPU.SimulateSlices(_slice_inside);

	/*irtkRealImage simslices(reconstructionGPU.v_simulated_slices.size.x, reconstructionGPU.v_simulated_slices.size.y, reconstructionGPU.v_simulated_slices.size.z);
	reconstructionGPU.debugSimslices(simslices.GetPointerToVoxels());
	simslices.Write("simslices.nii");

	irtkRealImage simweights(reconstructionGPU.v_simulated_weights.size.x, reconstructionGPU.v_simulated_weights.size.y, reconstructionGPU.v_simulated_weights.size.z);
	reconstructionGPU.debugSimweights(simweights.GetPointerToVoxels());
	simweights.Write("simweights.nii");*/

}

void irtkReconstructionCuda::SimulateStacks(vector<irtkRealImage>& stacks)
{
    if (_debug)
        cout<<"Simulating stacks."<<endl;
  
    unsigned int inputIndex;
    int i, j, k, n;
    irtkRealImage sim;
    POINT3D p;
    double weight;
  
    int z, current_stack;
    z=-1;//this is the z coordinate of the stack
    current_stack=-1; //we need to know when to start a new stack


    for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
      
	
        // read the current slice
        irtkRealImage& slice = _slices[inputIndex];

        //Calculate simulated slice
        sim.Initialize( slice.GetImageAttributes() );
        sim = 0;

	//do not simulate excluded slice
        if(_slice_weight[inputIndex]>0.5)
	{
          for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    weight=0;
                    n = _volcoeffs[inputIndex][i][j].size();
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
                        sim(i, j, 0) += p.value * _reconstructed(p.x, p.y, p.z);
                        weight += p.value;
                    }
                    if(weight>0)
                        sim(i,j,0)/=weight;
                }
	}

        if (_stack_index[inputIndex]==current_stack)
            z++;
        else {
            current_stack=_stack_index[inputIndex];
            z=0;
        }
        
        for(i=0;i<sim.GetX();i++)
            for(j=0;j<sim.GetY();j++) {
                stacks[_stack_index[inputIndex]](i,j,z)=sim(i,j,0);
            }
        //end of loop for a slice inputIndex
    }   
}

void irtkReconstructionCuda::MatchStackIntensities(vector<irtkRealImage>& stacks,
                                               vector<irtkRigidTransformation>& stack_transformations, double averageValue, bool together)
{
    if (_debug)
        cout << "Matching intensities of stacks. ";

	PerfStats stats;
	stats.start();

    //Calculate the averages of intensities for all stacks
    double sum, num;
    char buffer[256];
    unsigned int ind;
    int i, j, k;
    double x, y, z;
    vector<double> stack_average;
        
    //remember the set average value
    _average_value = averageValue;
    
    //averages need to be calculated only in ROI
    for (ind = 0; ind < stacks.size(); ind++) {
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
                    // if ((x >= 0) && (x < _mask.GetX()) && (y >= 0) && (y < _mask.GetY()) && (z >= 0)
                    //      && (z < _mask.GetZ()))
                    //      {
                    //if (_mask(x, y, z) == 1)
                    if ( stacks[ind](i, j, k) > 0 ) {
                        sum += stacks[ind](i, j, k);
                        num++;
                    }
                    //}
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

	stats.sample("MatchStackIntensities");
	stats.print();

}


void irtkReconstructionCuda::MatchStackIntensitiesWithMasking(vector<irtkRealImage>& stacks,
                                               vector<irtkRigidTransformation>& stack_transformations, double averageValue, bool together)
{
    if (_debug)
        cout << "Matching intensities of stacks. ";

	PerfStats stats;
	stats.start();

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
        m=stacks[ind];
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
                        sum += stacks[ind](i, j, k);
                        num++;
                      }
                      else
			m(i,j,k)=0;
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

	stats.sample("MatchStackIntensities");
	stats.print();

}


void irtkReconstructionCuda::CreateSlicesAndTransformations( vector<irtkRealImage> &stacks,
                                                         vector<irtkRigidTransformation> &stack_transformations,
                                                         vector<double> &thickness,
                                                         const vector<irtkRealImage> &probability_maps )
{
    if (_debug)
        cout << "CreateSlicesAndTransformations" << endl;
    
    //for each stack
    for (unsigned int i = 0; i < stacks.size(); i++) {
        //image attributes contain image and voxel size
        irtkImageAttributes attr = stacks[i].GetImageAttributes();

        //attr._z is number of slices in the stack
        for (int j = 0; j < attr._z; j++) {
            //create slice by selecting the appropreate region of the stack
            irtkRealImage slice = stacks[i].GetRegion(0, 0, j, attr._x, attr._y, j + 1);
            //set correct voxel size in the stack. Z size is equal to slice thickness.
            slice.PutPixelSize(attr._dx, attr._dy, thickness[i]);
			//std::cout << "slice stack " << i  << std::endl;
			//slice.GetImageToWorldMatrix().Print();
            //remember the slice
            _slices.push_back(slice);
            _simulated_slices.push_back(slice);
            _simulated_weights.push_back(slice);
            _simulated_inside.push_back(slice);
            //remeber stack index for this slice
            _stack_index.push_back(i);
            //initialize slice transformation with the stack transformation
			stack_transformations[i].UpdateParameter();
            _transformations.push_back(stack_transformations[i]);
            if ( probability_maps.size() > 0 ) {
                irtkRealImage proba = probability_maps[i].GetRegion(0, 0, j, attr._x, attr._y, j + 1);
                proba.PutPixelSize(attr._dx, attr._dy, thickness[i]);
                _probability_maps.push_back(proba);
            }
        }
    }
    cout << "Number of slices: " << _slices.size() << endl;

	//TEst///////////////////////////////////////////////7
	double dx, dy, dz;
	_slices[0].GetPixelSize(&dx, &dy, &dz);

	//sigma of 3D Gaussian (sinc with FWHM=dx or dy in-plane, Gaussian with FWHM = dz through-plane)
	double sigmax = 1.2 * dx / 2.3548;
	double sigmay = 1.2 * dy / 2.3548;
	double sigmaz = dz / 2.3548;

	double size = _reconstructed.GetXSize();

	//number of voxels in each direction
	//the ROI is 2*voxel dimension

	int xDim = PSF_SIZE;//round(2 * dx / size);
	int yDim = PSF_SIZE;//round(2 * dy / size);
	int zDim = PSF_SIZE;//round(2 * dz / size);

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
	PSF.GetImageToWorldMatrix().Print();
	//printf("%d %d %d \n\n", xDim, yDim, zDim);

	double x, y, z;
	double sum = 0;
	int i, j, k;
	for (i = 0; i < xDim; i++)
		for (j = 0; j < yDim; j++)
			for (k = 0; k < zDim; k++) {
				x = i;
				y = j;
				z = k;
				//printf("%f %f %f \n", x, y, z);
				PSF.ImageToWorld(x, y, z);
				//printf("%f %f %f \n", x, y, z);
				x -= cx;
				y -= cy;
				z -= cz;
				//printf("%f %f %f \n\n", x, y, z);
				//continuous PSF does not need to be normalized as discrete will be
				PSF(i, j, k) = exp(-x * x / (2 * sigmax * sigmax) - y * y / (2 * sigmay * sigmay)
					- z * z / (2 * sigmaz * sigmaz));
				sum += PSF(i, j, k);
			}
			PSF /= sum;

	//assumption for PSF volume -- all voxel dims same
	reconstructionGPU.generatePSFVolume(PSF.GetPointerToVoxels(), make_float3(_slices[0].GetXSize(),_slices[0].GetYSize(),_slices[0].GetZSize()));
}

void irtkReconstructionCuda::ResetSlices(vector<irtkRealImage>& stacks,
                                     vector<double>& thickness)
{
    if (_debug)
        cout << "ResetSlices" << endl;

    _slices.clear();
    
    //for each stack
    for (unsigned int i = 0; i < stacks.size(); i++) {
        //image attributes contain image and voxel size
        irtkImageAttributes attr = stacks[i].GetImageAttributes();

        //attr._z is number of slices in the stack
        for (int j = 0; j < attr._z; j++) {
            //create slice by selecting the appropreate region of the stack
            irtkRealImage slice = stacks[i].GetRegion(0, 0, j, attr._x, attr._y, j + 1);
            //set correct voxel size in the stack. Z size is equal to slice thickness.
            slice.PutPixelSize(attr._dx, attr._dy, thickness[i]);
            //remember the slice
            _slices.push_back(slice);
        }
    }
    cout << "Number of slices: " << _slices.size() << endl;

    for ( int i = 0; i < _slices.size(); i++ ) {
        _bias[i].Initialize( _slices[i].GetImageAttributes() );
        _weights[i].Initialize( _slices[i].GetImageAttributes() );

    }

}

void irtkReconstructionCuda::SetSlicesAndTransformations( vector<irtkRealImage>& slices,
                                                      vector<irtkRigidTransformation>& slice_transformations,
                                                      vector<int>& stack_ids,
                                                      vector<double>& thickness )
{
    _slices.clear();
    _stack_index.clear();
    _transformations.clear();
    
    //for each slice
    for (unsigned int i = 0; i < slices.size(); i++) {
        //get slice
        irtkRealImage slice = slices[i];
        std::cout << "setting slice " << i << "\n";
        slice.Print();
        //set correct voxel size in the stack. Z size is equal to slice thickness.
        slice.PutPixelSize(slice.GetXSize(), slice.GetYSize(), thickness[i]);
        //remember the slice
        _slices.push_back(slice);

		//reconstructionGPU.UpdateSlice(i, make_uint2(slice.GetX(), slice.GetY()), slice.GetPointerToVoxels());
        //remember stack index for this slice
        _stack_index.push_back(stack_ids[i]);
        //get slice transformation
        _transformations.push_back(slice_transformations[i]);
    }
}

void irtkReconstructionCuda::UpdateSlices(vector<irtkRealImage>& stacks, vector<double>& thickness)
{
    _slices.clear();
	//for each stack
	for (unsigned int i = 0; i < stacks.size(); i++) {
		//image attributes contain image and voxel size
		irtkImageAttributes attr = stacks[i].GetImageAttributes();

		//attr._z is number of slices in the stack
		for (int j = 0; j < attr._z; j++) {
			//create slice by selecting the appropreate region of the stack
			irtkRealImage slice = stacks[i].GetRegion(0, 0, j, attr._x, attr._y, j + 1);
			//set correct voxel size in the stack. Z size is equal to slice thickness.
			slice.PutPixelSize(attr._dx, attr._dy, thickness[i]);
			//remember the slice
			_slices.push_back(slice);
			//reconstructionGPU.UpdateSlice(i, make_uint2(slice.GetX(), slice.GetY()), slice.GetPointerToVoxels());
		}
	}
	cout << "Number of slices: " << _slices.size() << endl;
}

void irtkReconstructionCuda::MaskSlices()
{
    cout << "Masking slices ... ";

    double x, y, z;
    int i, j;

    //Check whether we have a mask
    if (!_have_mask) {
        cout << "Could not mask slices because no mask has been set." << endl;
        return;
    }

    //mask slices
    for (int unsigned inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        irtkRealImage& slice = _slices[inputIndex];
        for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++) {
                //if the value is smaller than 1 assume it is padding
                if (slice(i,j,0) < 0.01)
                    slice(i,j,0) = -1;
                //image coordinates of a slice voxel
                x = i;
                y = j;
                z = 0;
                //change to world coordinates in slice space
                slice.ImageToWorld(x, y, z);
                //world coordinates in volume space
                _transformations[inputIndex].Transform(x, y, z);
                //image coordinates in volume space
                _mask.WorldToImage(x, y, z);
                x = round(x);
                y = round(y);
                z = round(z);
                //if the voxel is outside mask ROI set it to -1 (padding value)
                if ((x >= 0) && (x < _mask.GetX()) && (y >= 0) && (y < _mask.GetY()) && (z >= 0)
                    && (z < _mask.GetZ())) {
                    if (_mask(x, y, z) == 0)
                        slice(i, j, 0) = -1;
                }
                else
                    slice(i, j, 0) = -1;
            }
        //remember masked slice
        //_slices[inputIndex] = slice;
    }
    cout << "done." << endl;
}

class ParallelSliceToVolumeRegistration {
public:
    irtkReconstructionCuda *reconstructor;

    ParallelSliceToVolumeRegistration(irtkReconstructionCuda *_reconstructor) : 
    reconstructor(_reconstructor) { }

    void operator() (const blocked_range<size_t> &r) const {

        irtkImageAttributes attr = reconstructor->_reconstructed.GetImageAttributes();
        
        for ( size_t inputIndex = r.begin(); inputIndex != r.end(); ++inputIndex ) {

            irtkImageRigidRegistrationWithPadding registration;
            irtkGreyPixel smin, smax;
            irtkGreyImage target;
            irtkRealImage slice, w, b, t;
            irtkResamplingWithPadding<irtkRealPixel> resampling(attr._dx,attr._dx,attr._dx,-1);         
            irtkReconstructionCuda dummy_reconstruction;
        
            //target = _slices[inputIndex];
            t = reconstructor->_slices[inputIndex];
            resampling.SetInput(&reconstructor->_slices[inputIndex]);
            resampling.SetOutput(&t);
            resampling.Run();
            target=t;

            target.GetMinMax(&smin, &smax);
        
            if (smax > -1) {
                //put origin to zero
                irtkRigidTransformation offset;
                dummy_reconstruction.ResetOrigin(target,offset);
                irtkMatrix mo = offset.GetMatrix();
                irtkMatrix m = reconstructor->_transformations[inputIndex].GetMatrix();
                m=m*mo;
                reconstructor->_transformations[inputIndex].PutMatrix(m);

                irtkGreyImage source = reconstructor->_reconstructed;
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
		PerfStats stats;
		stats.start();
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, reconstructor->_slices.size() ),
                      *this );
        init.terminate();
		stats.sample("ParallelSliceToVolumeRegistration");
		stats.print();
    }

};

void irtkReconstructionCuda::SliceToVolumeRegistration()
{
	PerfStats stats;
	stats.start();

    if (_debug)
        cout << "SliceToVolumeRegistration" << endl;
    ParallelSliceToVolumeRegistration registration(this);
    registration();

	stats.sample("SliceToVolumeRegistration");
	stats.print();
}

class ParallelCoeffInit {
public:
    irtkReconstructionCuda *reconstructor;

    ParallelCoeffInit(irtkReconstructionCuda *_reconstructor) : 
    reconstructor(_reconstructor) { }

    void operator() (const blocked_range<size_t> &r) const {
        
        for ( size_t inputIndex = r.begin(); inputIndex != r.end(); ++inputIndex ) {

			//TODO cuda
            bool slice_inside;

            //current slice
            //irtkRealImage slice;

            //get resolution of the volume
            double vx, vy, vz;
            reconstructor->_reconstructed.GetPixelSize(&vx, &vy, &vz);
            //volume is always isotropic
            double res = vx;
        
            //start of a loop for a slice inputIndex
            cout << inputIndex << " ";

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
            double sigmax = 1.2 * dx / 2.3548;
            double sigmay = 1.2 * dy / 2.3548;
            double sigmaz = dz / 2.3548;
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

			//printf("xDim %d yDim %d zDim %d \n", xDim, yDim, zDim);

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
                        PSF(i, j, k) = exp(-x * x / (2 * sigmax * sigmax) - y * y / (2 * sigmay * sigmay)
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
			//printf("dim %d \n", dim); //7

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
                        reconstructor->_reconstructed.WorldToImage(x, y, z);
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
                                    reconstructor->_reconstructed.WorldToImage(x, y, z);

									///linear interpolation start

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
                                        if ((l >= 0) && (l < reconstructor->_reconstructed.GetX()))
                                            for (m = ny; m <= ny + 1; m++)
                                                if ((m >= 0) && (m < reconstructor->_reconstructed.GetY()))
                                                    for (n = nz; n <= nz + 1; n++)
                                                        if ((n >= 0) && (n < reconstructor->_reconstructed.GetZ())) {
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
                                        if ((l >= 0) && (l < reconstructor->_reconstructed.GetX()))
                                            for (m = ny; m <= ny + 1; m++)
                                                if ((m >= 0) && (m < reconstructor->_reconstructed.GetY()))
                                                    for (n = nz; n <= nz + 1; n++)
                                                        if ((n >= 0) && (n < reconstructor->_reconstructed.GetZ())) {
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
											///linear interpolation end

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

void irtkReconstructionCuda::CoeffInit()
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
    ParallelCoeffInit coeffinit(this);
    coeffinit();
    cout << " ... done." << endl;

    //prepare image for volume weights, will be needed for Gaussian Reconstruction
    _volume_weights.Initialize( _reconstructed.GetImageAttributes() );
    _volume_weights = 0;

    int inputIndex, i, j, n, k;
    POINT3D p;
    for ( inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
        for ( i = 0; i < _slices[inputIndex].GetX(); i++)
            for ( j = 0; j < _slices[inputIndex].GetY(); j++) {
                n = _volcoeffs[inputIndex][i][j].size();
                for (k = 0; k < n; k++) {
                    p = _volcoeffs[inputIndex][i][j][k];
                    _volume_weights(p.x, p.y, p.z) += p.value;
                }
            }
    }
    if (_debug)
        _volume_weights.Write("volume_weights.nii.gz");
    
	//obsolete
    //find average volume weight to modify alpha parameters accordingly
    irtkRealPixel *ptr = _volume_weights.GetPointerToVoxels();
    irtkRealPixel *pm = _mask.GetPointerToVoxels();
    double sum = 0;
    int num=0;
    for (int i=0;i<_volume_weights.GetNumberOfVoxels();i++) {
        if (*pm==1) {
            sum+=*ptr;
            num++;
        }
        ptr++;
        pm++;
    }
    _average_volume_weight = sum/num;

#if !TEX_TEST
	reconstructionGPU.UpdateSliceCoeff(_volcoeffs, _volume_weights.GetPointerToVoxels());
	reconstructionGPU.UpdateScaleVector(_scale, _slice_weight);
#endif

    if(_debug) {
        cout<<"Average volume weight is "<<_average_volume_weight<<endl;
    }
    
}  //end of CoeffInit()



void irtkReconstructionCuda::GaussianReconstruction()
{
	int i, j, k, n;
	vector<int> voxel_num;  
	reconstructionGPU.GaussianReconstruction(voxel_num);

	//now find slices with small overlap with ROI and exclude them.
	vector<int> voxel_num_tmp;
	for (i=0;i<voxel_num.size();i++)
		voxel_num_tmp.push_back(voxel_num[i]);

	//find median
	sort(voxel_num_tmp.begin(),voxel_num_tmp.end());
	int median = voxel_num_tmp[round(voxel_num_tmp.size()*0.5)];
	// printf("\n median %d \n", median);

	//remember slices with small overlap with ROI
	_small_slices.clear();
	for (i=0;i<voxel_num.size();i++)
		if (voxel_num[i]<0.1*median)
			_small_slices.push_back(i);

	if (_debug) {
		cout<<"Small slices:";
		for (i=0;i<_small_slices.size();i++)
			cout<<" "<<_small_slices[i];
		cout<<endl;
	}
}

void irtkReconstructionCuda::InitializeEM()
{
    if (_debug)
        cout << "InitializeEM" << endl;
    
    _weights.clear();
    _bias.clear();
    _scale.clear();
    _slice_weight.clear();
    
    for (unsigned int i = 0; i < _slices.size(); i++) {
        //Create images for voxel weights and bias fields
        _weights.push_back(_slices[i]);
        _bias.push_back(_slices[i]);

        //Create and initialize scales
        _scale.push_back(1);

        //Create and initialize slice weights
        _slice_weight.push_back(1);
    }

    //Find the range of intensities
    _max_intensity = voxel_limits<irtkRealPixel>::min();
    _min_intensity = voxel_limits<irtkRealPixel>::max();
    for (unsigned int i = 0; i < _slices.size(); i++) {
        //to update minimum we need to exclude padding value
        irtkRealPixel *ptr = _slices[i].GetPointerToVoxels();
        for (int ind = 0; ind < _slices[i].GetNumberOfVoxels(); ind++) {
            if (*ptr > 0) {
                if (*ptr > _max_intensity)
                    _max_intensity = *ptr;
                if (*ptr < _min_intensity)
                    _min_intensity = *ptr;
            }
            ptr++;
        }
    }
}

void irtkReconstructionCuda::InitializeEMValues()
{
    if (_debug)
        cout << "InitializeEMValues" << endl;

	_slice_weight.assign(_slices.size(), 1);
    _scale.assign(_slices.size(), 1);

	reconstructionGPU.InitializeEMValues();

	reconstructionGPU.UpdateScaleVector(_scale, _slice_weight);
}

void irtkReconstructionCuda::InitializeRobustStatistics()
{
    if (_debug)
        cout << "InitializeRobustStatistics" << endl;

	reconstructionGPU.InitializeRobustStatistics(_sigma);

	for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
		//if slice does not have an overlap with ROI, set its weight to zero
		if (!_slice_inside[inputIndex])
			_slice_weight[inputIndex] = 0;
	}
	//Force exclusion of slices predefined by user
    for (unsigned int i = 0; i < _force_excluded.size(); i++)
        _slice_weight[_force_excluded[i]] = 0;


    //initialize sigma for slice-wise robust statistics
    _sigma_s = 0.025;
    //initialize mixing proportion for inlier class in voxel-wise robust statistics
    _mix = 0.9;
    //initialize mixing proportion for outlier class in slice-wise robust statistics
    _mix_s = 0.9;
    //Initialise value for uniform distribution according to the range of intensities
    _m = 1 / (2.1 * _max_intensity - 1.9 * _min_intensity);

    if (_debug)
        cout << "Initializing robust statistics: " << "sigma=" << sqrt(_sigma) << " " << "m=" << _m
             << " " << "mix=" << _mix << " " << "mix_s=" << _mix_s << endl;

	reconstructionGPU.UpdateScaleVector(_scale, _slice_weight);

}

void irtkReconstructionCuda::EStep()
{
    //EStep performs calculation of voxel-wise and slice-wise posteriors (weights)
    if (_debug)
        cout << "EStep: " << endl;

	unsigned int inputIndex;
	vector<double> slice_potential(_slices.size(), 0);

	reconstructionGPU.EStep(_m, _sigma, _mix, slice_potential);

/*	irtkRealImage bweights(reconstructionGPU.v_weights.size.x, reconstructionGPU.v_weights.size.y, reconstructionGPU.v_weights.size.z);
	reconstructionGPU.debugWeights(bweights.GetPointerToVoxels());
	bweights.Write("weigthField.nii");*/

	//can stay on CPU

    //To force-exclude slices predefined by a user, set their potentials to -1
    for (unsigned int i = 0; i < _force_excluded.size(); i++)
        slice_potential[_force_excluded[i]] = -1;
    
    //exclude slices identified as having small overlap with ROI, set their potentials to -1
    for (unsigned int i = 0; i < _small_slices.size(); i++)
        slice_potential[_small_slices[i]] = -1;
    
    //these are unrealistic scales pointing at misregistration - exclude the corresponding slices
    for (inputIndex = 0; inputIndex < slice_potential.size(); inputIndex++)
        if ((_scale[inputIndex]<0.2)||(_scale[inputIndex]>5)) {
            slice_potential[inputIndex] = -1;
        }

    // exclude unrealistic transformations
    /*
      int current_stack = 0;
      int nb_stacks = _stack_index[_stack_index.size()-1];
      double tx,ty,tz,rx,ry,rz, nb;
      for ( int i = 0; i < nb_stacks; i++ ) {
      tx = 0;
      ty = 0;
      tz = 0;
      rx = 0;
      ry = 0;
      rz = 0;
      nb = 0;
      for ( int j = 0; j < _slices.size(); j++ ) {
      if ( _stack_index[j] != i )
      continue;
      if ( slice_potential[j] == -1 )
      continue;
      tx += _transformations[j].GetTranslationX();
      ty += _transformations[j].GetTranslationY();
      tz += _transformations[j].GetTranslationZ();
      rx += _transformations[j].GetRotationX();
      ry += _transformations[j].GetRotationY();
      rz += _transformations[j].GetRotationZ();
      nb++;
      }
      tx /= nb;
      ty /= nb;
      tz /= nb;
      rx /= nb;
      ry /= nb;
      rz /= nb;
      for ( int j = 0; j < _slices.size(); j++ ) {
      if ( _stack_index[j] != i )
      continue;
      if ( slice_potential[j] == -1 )
      continue;
      if ( abs( tx - _transformations[j].GetTranslationX() ) > 20
      || abs( ty - _transformations[j].GetTranslationY() ) > 20
      || abs( tz - _transformations[j].GetTranslationZ() ) > 20
      || abs( rx - _transformations[j].GetRotationX() ) > 5
      || abs( ry - _transformations[j].GetRotationY() ) > 5
      || abs( rz - _transformations[j].GetRotationZ() ) > 5 )
      slice_potential[j] = -1;
      }        
      }
    */  
    if(_debug) {
        cout<<endl<<"Slice potentials: ";
        for (inputIndex = 0; inputIndex < slice_potential.size(); inputIndex++)
            cout<<slice_potential[inputIndex]<<" ";
        cout<<endl;
    }


    //Calulation of slice-wise robust statistics parameters.
    //This is theoretically M-step,
    //but we want to use latest estimate of slice potentials
    //to update the parameters

    //Calculate means of the inlier and outlier potentials
    double sum = 0, den = 0, sum2 = 0, den2 = 0, maxs = 0, mins = 1;
    for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
        if (slice_potential[inputIndex] >= 0) {
            //calculate means
            sum += slice_potential[inputIndex] * _slice_weight[inputIndex];
            den += _slice_weight[inputIndex];
            sum2 += slice_potential[inputIndex] * (1 - _slice_weight[inputIndex]);
            den2 += (1 - _slice_weight[inputIndex]);

            //calculate min and max of potentials in case means need to be initalized
            if (slice_potential[inputIndex] > maxs)
                maxs = slice_potential[inputIndex];
            if (slice_potential[inputIndex] < mins)
                mins = slice_potential[inputIndex];
        }

    if (den > 0)
        _mean_s = sum / den;
    else
        _mean_s = mins;

    if (den2 > 0)
        _mean_s2 = sum2 / den2;
    else
        _mean_s2 = (maxs + _mean_s) / 2;

    //Calculate the variances of the potentials
    sum = 0;
    den = 0;
    sum2 = 0;
    den2 = 0;
    for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
        if (slice_potential[inputIndex] >= 0) {
            sum += (slice_potential[inputIndex] - _mean_s) * (slice_potential[inputIndex] - _mean_s)
                * _slice_weight[inputIndex];
            den += _slice_weight[inputIndex];

            sum2 += (slice_potential[inputIndex] - _mean_s2) * (slice_potential[inputIndex] - _mean_s2)
                * (1 - _slice_weight[inputIndex]);
            den2 += (1 - _slice_weight[inputIndex]);

        }

    //_sigma_s
    if ((sum > 0) && (den > 0)) {
        _sigma_s = sum / den;
        //do not allow too small sigma
        if (_sigma_s < _step * _step / 6.28)
            _sigma_s = _step * _step / 6.28;
    }
    else {
        _sigma_s = 0.025;
        if (_debug) {
            if (sum <= 0)
                cout << "All slices are equal. ";
            if (den < 0) //this should not happen
                cout << "All slices are outliers. ";
            cout << "Setting sigma to " << sqrt(_sigma_s) << endl;
        }
    }

    //sigma_s2
    if ((sum2 > 0) && (den2 > 0)) {
        _sigma_s2 = sum2 / den2;
        //do not allow too small sigma
        if (_sigma_s2 < _step * _step / 6.28)
            _sigma_s2 = _step * _step / 6.28;
    }
    else {
        _sigma_s2 = (_mean_s2 - _mean_s) * (_mean_s2 - _mean_s) / 4;
        //do not allow too small sigma
        if (_sigma_s2 < _step * _step / 6.28)
            _sigma_s2 = _step * _step / 6.28;

        if (_debug) {
            if (sum2 <= 0)
                cout << "All slices are equal. ";
            if (den2 <= 0)
                cout << "All slices inliers. ";
            cout << "Setting sigma_s2 to " << sqrt(_sigma_s2) << endl;
        }
    }

    //Calculate slice weights
    double gs1, gs2;
    for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        //Slice does not have any voxels in volumetric ROI
        if (slice_potential[inputIndex] == -1) {
            _slice_weight[inputIndex] = 0;
            continue;
        }

        //All slices are outliers or the means are not valid
        if ((den <= 0) || (_mean_s2 <= _mean_s)) {
            _slice_weight[inputIndex] = 1;
            continue;
        }

        //likelihood for inliers
        if (slice_potential[inputIndex] < _mean_s2)
            gs1 = G(slice_potential[inputIndex] - _mean_s, _sigma_s);
        else
            gs1 = 0;

        //likelihood for outliers
        if (slice_potential[inputIndex] > _mean_s)
            gs2 = G(slice_potential[inputIndex] - _mean_s2, _sigma_s2);
        else
            gs2 = 0;

        //calculate slice weight
        double likelihood = gs1 * _mix_s + gs2 * (1 - _mix_s);
        if (likelihood > 0)
            _slice_weight[inputIndex] = gs1 * _mix_s / likelihood;
        else {
            if (slice_potential[inputIndex] <= _mean_s)
                _slice_weight[inputIndex] = 1;
            if (slice_potential[inputIndex] >= _mean_s2)
                _slice_weight[inputIndex] = 0;
            if ((slice_potential[inputIndex] < _mean_s2) && (slice_potential[inputIndex] > _mean_s)) //should not happen
                _slice_weight[inputIndex] = 1;
        }
    }

    //Update _mix_s this should also be part of MStep
    sum = 0;
	int num = 0;
    for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
        if (slice_potential[inputIndex] >= 0) {
            sum += _slice_weight[inputIndex];
            num++;
        }

    if (num > 0)
        _mix_s = sum / num;
    else {
        cout << "All slices are outliers. Setting _mix_s to 0.9." << endl;
        _mix_s = 0.9;
    }

    if (_debug) {
        cout << setprecision(3);
        cout << "Slice robust statistics parameters: ";
        cout << "means: " << _mean_s << " " << _mean_s2 << "  ";
        cout << "sigmas: " << sqrt(_sigma_s) << " " << sqrt(_sigma_s2) << "  ";
        cout << "proportions: " << _mix_s << " " << 1 - _mix_s << endl;
        cout << "Slice weights: ";
        for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
            cout << _slice_weight[inputIndex] << " ";
        cout << endl;
    }

	//TODO only slice weight
	reconstructionGPU.UpdateSliceWeights(_slice_weight);

}

void irtkReconstructionCuda::Scale()
{
    if (_debug)
        cout << "Scale" << endl;

	reconstructionGPU.CalculateScaleVector();
	_scale = reconstructionGPU.h_scales;
}

irtkRealImage irtkReconstructionCuda::GetReconstructedGPU()
{
	reconstructionGPU.syncCPU(_reconstructed.GetPointerToVoxels());
    return _reconstructed;
}

void irtkReconstructionCuda::SyncCPU()
{
	/*for(int i = 0; i < _slices.size(); i++)
	{
		reconstructionGPU.SyncBSlice(i, _bias[i].GetPointerToVoxels());
	}*/
	reconstructionGPU.syncCPU(_reconstructed.GetPointerToVoxels());

}

Matrix4 irtkReconstructionCuda::toMatrix4(irtkMatrix mat)
{
	Matrix4 mmat;
	mmat.data[0] = make_double4(mat(0,0), mat(0,1), mat(0,2), mat(0,3));
	mmat.data[1] = make_double4(mat(1,0), mat(1,1), mat(1,2), mat(1,3));
	mmat.data[2] = make_double4(mat(2,0), mat(2,1), mat(2,2), mat(2,3));
	mmat.data[3] = make_double4(mat(3,0), mat(3,1), mat(3,2), mat(3,3));
	return mmat;
}

void irtkReconstructionCuda::UpdateGPUTranformationMatrices()
{

	std::vector<Matrix4> sI2W;
	std::vector<Matrix4> sW2I;
	for(int i = 0; i < _slices.size(); i++)
	{
		_transformations[i].UpdateMatrix();
		irtkMatrix m_new = _transformations[i].GetMatrix()*_slices[i].GetImageToWorldMatrix();
		Matrix4 si2w = toMatrix4(m_new);
		sI2W.push_back(si2w);
		irtkMatrix i_m_new = m_new;
		i_m_new.Invert();
		Matrix4 w2i = toMatrix4(i_m_new);
		sW2I.push_back(w2i);
		//std::cout << "m_new: " << std::endl;
		//m_new.Print();
		//std::cout << "i_m_new: " << std::endl;
		//i_m_new.Print();
		//std::cout << "slice_transformations: " << std::endl;

	}
	reconstructionGPU.SetSliceMatrices(sI2W, sW2I, toMatrix4(_reconstructed.GetImageToWorldMatrix()), toMatrix4(_reconstructed.GetWorldToImageMatrix()));
}


void irtkReconstructionCuda::SyncGPU()
{

	cudaDeviceSynchronize();

	reconstructionGPU.InitReconstructionVolume(make_uint3(_reconstructed.GetX(), _reconstructed.GetY(), _reconstructed.GetZ()), 
		make_float3(_reconstructed.GetXSize(), _reconstructed.GetYSize(), _reconstructed.GetZSize()), _reconstructed.GetPointerToVoxels(),  _sigma_bias/_reconstructed.GetXSize() );
	reconstructionGPU.setMask(make_uint3(_mask.GetX(), _mask.GetY(), _mask.GetZ()), make_float3(_mask.GetXSize(), _mask.GetYSize(), _mask.GetZSize()), 
		_mask.GetPointerToVoxels(), _sigma_bias/_reconstructed.GetXSize());

	uint3 reqVDims = make_uint3(INT_MIN, INT_MIN, 0);
	uint3 wasteVDims = make_uint3(INT_MAX, INT_MAX, 0);

	for(int i = 0; i < _slices.size(); i++)
	{
		reqVDims = make_uint3(max(reqVDims.x, _slices[i].GetX()), max(reqVDims.y, _slices[i].GetY()), _slices.size());
		wasteVDims = make_uint3(min(wasteVDims.x, _slices[i].GetX()), min(wasteVDims.y, _slices[i].GetY()), _slices.size());
	}

	double waste = (((reqVDims.x-wasteVDims.x) * (reqVDims.y-wasteVDims.y) * _slices.size())*sizeof(double)*5.0)/1024.0;
	printf("GPU memory waste approx: %f KB with %d %d %d %d\n", waste, reqVDims.x, reqVDims.y, wasteVDims.x, wasteVDims.y);

	irtkRealImage combinedStacks(reqVDims.x, reqVDims.y, reqVDims.z);
	//combinedStacks.Clear();

	//all the same dimensions
	reconstructionGPU.initStorageVolumes(reqVDims, make_float3(_slices[0].GetXSize(), _slices[0].GetYSize(), _slices[0].GetZSize()));
	vector<int> sizesX;
	vector<int> sizesY;
	std::vector<float3> slice_dims;
	combinedStacks = -1.0;

	double* ptr = combinedStacks.GetPointerToVoxels();
	for(int i = 0; i < _slices.size(); i++)
	{
		irtkRealImage slice = _slices[i];
		//We need to do this line wise because of different cropping sizes
		for(int y = 0; y < slice.GetY(); y++)
		{
			memcpy(ptr, slice.GetPointerToVoxels()+y*slice.GetX(), slice.GetX()*sizeof(double));
			ptr += combinedStacks.GetX();
		}
		ptr += abs(combinedStacks.GetY()-slice.GetY())*combinedStacks.GetX();
		//TODO - set slice sizes
		sizesX.push_back(slice.GetX());
		sizesY.push_back(slice.GetY());
		slice_dims.push_back(make_float3(_slices[i].GetXSize(), _slices[i].GetYSize(), _slices[i].GetZSize()));
	}
	//combinedStacks.Write("combinedStacks.nii");

	reconstructionGPU.FillSlices(combinedStacks.GetPointerToVoxels(), sizesX, sizesY);

	reconstructionGPU.setSliceDims(slice_dims);
	reconstructionGPU.reconstructedVoxelSize = _reconstructed.GetXSize();
	reconstructionGPU.UpdateScaleVector(_scale, _slice_weight);

	//debug
	/*irtkRealImage debugSlices(reqVDims.x, reqVDims.y, reqVDims.z);
	reconstructionGPU.getSlicesVol_debug(debugSlices.GetPointerToVoxels());
	debugSlices.Write("debugSlices.nii");*/
}


void irtkReconstructionCuda::Bias()
{
    if (_debug)
        cout << "Correcting bias ...";

	if(_global_bias_correction)
		printf("_global_bias_correction is not yet fully implemented in CUDA\n");

	reconstructionGPU.CorrectBias(_sigma_bias/_slices[0].GetXSize(), _global_bias_correction); //assuming globally constant pixel size

	/*irtkRealImage bimg(reconstructionGPU.v_bias.size.x, reconstructionGPU.v_bias.size.y, reconstructionGPU.v_bias.size.z);
	reconstructionGPU.debugBias(bimg.GetPointerToVoxels());
	bimg.Write("biasField.nii");*/
	
    if (_debug)
        cout << "done. " << endl;
}

void irtkReconstructionCuda::Superresolution(int iter)
{
    if (_debug)
        cout << "Superresolution " << iter << endl;

	int i, j, k;
	irtkRealImage addon, original;
	//Remember current reconstruction for edge-preserving smoothing
	original = _reconstructed;

	reconstructionGPU.Superresolution(iter, _slice_weight, _adaptive, _alpha, _min_intensity, _max_intensity, _delta, 
		_lambda, _global_bias_correction, _sigma_bias/_slices[0].GetXSize(), _low_intensity_cutoff); //assuming isotrop constant voxel size
}

void irtkReconstructionCuda::MStep(int iter)
{
	reconstructionGPU.MStep(iter, _step, _sigma, _mix, _m);

    if (_debug) {
        cout << "Voxel-wise robust statistics parameters: ";
        cout << "sigma = " << sqrt(_sigma) << " mix = " << _mix << " ";
        cout << " m = " << _m << endl;
    }

}
/*
void irtkReconstructionCuda::AdaptiveRegularization(int iter, irtkRealImage& original)
{
    if (_debug)
        cout << "AdaptiveRegularization" << endl;

    vector<double> factor(13,0);
    for (int i = 0; i < 13; i++) {
        for (int j = 0; j < 3; j++)
            factor[i] += fabs(double(_directions[i][j]));
        factor[i] = 1 / factor[i];
			//	printf("cpu factor %f\n", factor[i]);
    }

    vector<irtkRealImage> b;//(13);
    for (int i = 0; i < 13; i++)
        b.push_back( _reconstructed );

    ParallelAdaptiveRegularization1 parallelAdaptiveRegularization1( this,
                                                                     b,
                                                                     factor,
                                                                     original );
    parallelAdaptiveRegularization1();

    irtkRealImage original2 = _reconstructed;
    ParallelAdaptiveRegularization2 parallelAdaptiveRegularization2( this,
                                                                     b,
                                                                     factor,
                                                                     original2 );
    parallelAdaptiveRegularization2();

    if (_alpha * _lambda / (_delta * _delta) > 0.068) {
        cerr
            << "Warning: regularization might not have smoothing effect! Ensure that alpha*lambda/delta^2 is below 0.068."
            << endl;
    }
}*/

void irtkReconstructionCuda::BiasCorrectVolume(irtkRealImage& original)
{
    //remove low-frequancy component in the reconstructed image which might have accured due to overfitting of the biasfield
    irtkRealImage residual = _reconstructed;
    irtkRealImage weights = _mask;

    //_reconstructed.Write("super-notbiascor.nii.gz");

    //calculate weighted residual
    irtkRealPixel *pr = residual.GetPointerToVoxels();
    irtkRealPixel *po = original.GetPointerToVoxels();
    irtkRealPixel *pw = weights.GetPointerToVoxels();
    for (int i = 0; i < _reconstructed.GetNumberOfVoxels(); i++) {
        //second and term to avoid numerical problems
        if ((*pw == 1) && (*po > _low_intensity_cutoff * _max_intensity)
            && (*pr > _low_intensity_cutoff * _max_intensity)) {
            *pr /= *po;
            *pr = log(*pr);
        }
        else {
            *pw = 0;
            *pr = 0;
        }
        pr++;
        po++;
        pw++;
    }
    //residual.Write("residual.nii.gz");
    //blurring needs to be same as for slices
    irtkGaussianBlurring<irtkRealPixel> gb(_sigma_bias);
    //blur weigted residual
    gb.SetInput(&residual);
    gb.SetOutput(&residual);
    gb.Run();
    //blur weight image
    gb.SetInput(&weights);
    gb.SetOutput(&weights);
    gb.Run();

    //calculate the bias field
    pr = residual.GetPointerToVoxels();
    pw = weights.GetPointerToVoxels();
    irtkRealPixel *pm = _mask.GetPointerToVoxels();
    irtkRealPixel *pi = _reconstructed.GetPointerToVoxels();
    for (int i = 0; i < _reconstructed.GetNumberOfVoxels(); i++) {

        if (*pm == 1) {
            //weighted gaussian smoothing
            *pr /= *pw;
            //exponential to recover multiplicative bias field
            *pr = exp(*pr);
            //bias correct reconstructed
            *pi /= *pr;
            //clamp intensities to allowed range
            if (*pi < _min_intensity * 0.9)
                *pi = _min_intensity * 0.9;
            if (*pi > _max_intensity * 1.1)
                *pi = _max_intensity * 1.1;
        }
        else {
            *pr = 0;
        }
        pr++;
        pw++;
        pm++;
        pi++;
    }

    //residual.Write("biasfield.nii.gz");
    //_reconstructed.Write("super-biascor.nii.gz");

}

void irtkReconstructionCuda::Evaluate(int iter)
{
    cout << "Iteration " << iter << ": " << endl;

    cout << "Included slices: ";
    int sum = 0;
    unsigned int i;
    for (i = 0; i < _slices.size(); i++) {
        if ((_slice_weight[i] >= 0.5) && (_slice_inside[i])) {
            cout << i << " ";
            sum++;
        }
    }
    cout << endl << "Total: " << sum << endl;

    cout << "Excluded slices: ";
    sum = 0;
    for (i = 0; i < _slices.size(); i++) {
        if ((_slice_weight[i] < 0.5) && (_slice_inside[i])) {
            cout << i << " ";
            sum++;
        }
    }
    cout << endl << "Total: " << sum << endl;

    cout << "Outside slices: ";
    sum = 0;
    for (i = 0; i < _slices.size(); i++) {
        if (!(_slice_inside[i])) {
            cout << i << " ";
            sum++;
        }
    }
    cout << endl << "Total: " << sum << endl;

}


void irtkReconstructionCuda::NormaliseBias(int iter)
{


	

	reconstructionGPU.NormaliseBias(iter, _sigma_bias/_reconstructed.GetXSize());

	/*char buffer[256];
	irtkRealImage nbias(reconstructionGPU.bias_.size.x, reconstructionGPU.bias_.size.y,  reconstructionGPU.bias_.size.z);
	reconstructionGPU.debugNormalizeBias(nbias.GetPointerToVoxels());
	sprintf(buffer,"nbiasGPU%i.nii",iter);
	nbias.Write(buffer);
	double mmin, mmax;
	nbias.GetMinMax(&mmin, &mmax);
	printf("nbias minmax %f %f \n", mmin, mmax);*/

	/*irtkRealImage smoothMask(reconstructionGPU.bias_.size.x, reconstructionGPU.bias_.size.y,  reconstructionGPU.bias_.size.z);
	reconstructionGPU.debugSmoothMask(smoothMask.GetPointerToVoxels());
	sprintf(buffer,"smoothMaskGPU%i.nii",iter);
	smoothMask.Write(buffer);
	smoothMask.GetMinMax(&mmin, &mmax);
	printf("mask minmax %f %f \n", mmin, mmax);*/
}

/* Set/Get/Save operations */

void irtkReconstructionCuda::ReadTransformation(char* folder)
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
            sprintf(name, "/transformation%i.dof", i);
            strcpy(path, folder);
            strcat(path, name);
        }
        else {
            sprintf(path, "transformation%i.dof", i);
        }
        transformation = irtkTransformation::New(path);
        rigidTransf = dynamic_cast<irtkRigidTransformation*>(transformation);
        _transformations.push_back(*rigidTransf);
        delete transformation;
        cout << path << endl;
    }
}

void irtkReconstructionCuda::SetReconstructed( irtkRealImage &reconstructed )
{
    _reconstructed = reconstructed;
    _template_created = true;
}

void irtkReconstructionCuda::SetTransformations( vector<irtkRigidTransformation>& transformations )
{
    _transformations.clear();
    for ( int i = 0; i < transformations.size(); i++ )
        _transformations.push_back(transformations[i]);

}

void irtkReconstructionCuda::SaveBiasFields()
{
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        sprintf(buffer, "bias%i.nii.gz", inputIndex);
        _bias[inputIndex].Write(buffer);
    }
}

void irtkReconstructionCuda::SaveConfidenceMap()
{
    _confidence_map.Write("confidence-map.nii.gz");
}

void irtkReconstructionCuda::SaveSlices()
{
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
        {
            sprintf(buffer, "slice%i.nii.gz", inputIndex);
            _slices[inputIndex].Write(buffer);
        }
}

void irtkReconstructionCuda::SaveWeights()
{
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        sprintf(buffer, "weights%i.nii.gz", inputIndex);
        _weights[inputIndex].Write(buffer);
    }
}

void irtkReconstructionCuda::SaveTransformations()
{
    char buffer[256];
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
        sprintf(buffer, "transformation%i.dof", inputIndex);
        _transformations[inputIndex].irtkTransformation::Write(buffer);
    }
}

void irtkReconstructionCuda::GetTransformations( vector<irtkRigidTransformation> &transformations )
{
    transformations.clear();
    for (unsigned int inputIndex = 0; inputIndex < _slices.size(); inputIndex++)
        transformations.push_back( _transformations[inputIndex] );
}

void irtkReconstructionCuda::GetSlices( vector<irtkRealImage> &slices )
{
    slices.clear();    
    for (unsigned int i = 0; i < _slices.size(); i++)
        slices.push_back(_slices[i]);
}

void irtkReconstructionCuda::SaveProbabilityMap( int i )
{
    char buffer[256];
    sprintf( buffer, "probability_map%i.nii", i );
    _brain_probability.Write( buffer );
}

void irtkReconstructionCuda::SlicesInfo( const char* filename,
                                     vector<string> &stack_files )
{
    ofstream info;
    info.open( filename );

    // header
    info << "stack_index" << "\t"
         << "stack_name" << "\t"
         << "included" << "\t" // Included slices
         << "excluded" << "\t"  // Excluded slices
         << "outside" << "\t"  // Outside slices
         << "weight" << "\t"
         << "scale" << "\t"
        // << _stack_factor[i] << "\t"
         << "TranslationX" << "\t"
         << "TranslationY" << "\t"
         << "TranslationZ" << "\t"
         << "RotationX" << "\t"
         << "RotationY" << "\t"
         << "RotationZ" << endl;
    
    for (int i = 0; i < _slices.size(); i++) {
        irtkRigidTransformation& t = _transformations[i];
        info << _stack_index[i] << "\t"
             << stack_files[_stack_index[i]] << "\t"
             << (((_slice_weight[i] >= 0.5) && (_slice_inside[i]))?1:0) << "\t" // Included slices
             << (((_slice_weight[i] < 0.5) && (_slice_inside[i]))?1:0) << "\t"  // Excluded slices
             << ((!(_slice_inside[i]))?1:0) << "\t"  // Outside slices
             << _slice_weight[i] << "\t"
             << _scale[i] << "\t"
            // << _stack_factor[i] << "\t"
             << t.GetTranslationX() << "\t"
             << t.GetTranslationY() << "\t"
             << t.GetTranslationZ() << "\t"
             << t.GetRotationX() << "\t"
             << t.GetRotationY() << "\t"
             << t.GetRotationZ() << endl;
    }
 
    info.close(); 
}


/* end Set/Get/Save operations */

/* Package specific functions */
void irtkReconstructionCuda::SplitImage(irtkRealImage image, int packages, vector<irtkRealImage>& stacks)
{
    irtkImageAttributes attr = image.GetImageAttributes();
  
    //slices in package
    int pkg_z = attr._z/packages;
    double pkg_dz = attr._dz*packages;
    //cout<<"packages: "<<packages<<"; slices: "<<attr._z<<"; slices in package: "<<pkg_z<<endl;
    //cout<<"slice thickness "<<attr._dz<<"; slickess thickness in package: "<<pkg_dz<<endl;
  
    char buffer[256];
    int i,j,k,l;
    double x,y,z,sx,sy,sz,ox,oy,oz;
    for(l=0;l<packages;l++) {
        attr = image.GetImageAttributes();
        if((pkg_z*packages+l)<attr._z)
            attr._z = pkg_z+1;
        else
            attr._z = pkg_z;
        attr._dz = pkg_dz;
    
        cout<<"split image "<<l<<" has "<<attr._z<<" slices."<<endl;
  
        //fill values in each stack
        irtkRealImage stack(attr);
        stack.GetOrigin(ox,oy,oz);

        cout<<"Stack "<<l<<":"<<endl;
        for(k=0; k<stack.GetZ();k++)
            for(j=0; j<stack.GetY();j++)
                for(i=0; i<stack.GetX();i++)
                    stack.Put(i,j,k,image(i,j,k*packages+l));
    
        //adjust origin
     
        //original image coordinates
        x=0;y=0;z=l;
        image.ImageToWorld(x,y,z);
        cout<<"image: "<<x<<" "<<y<<" "<<z<<endl;
        //stack coordinates
        sx=0;sy=0;sz=0;
        stack.PutOrigin(ox,oy,oz); //adjust to original value
        stack.ImageToWorld(sx,sy,sz);
        cout<<"stack: "<<sx<<" "<<sy<<" "<<sz<<endl;
        //adjust origin
        cout<<"adjustment needed: "<<x-sx<<" "<<y-sy<<" "<<z-sz<<endl;
        stack.PutOrigin(ox + (x-sx), oy + (y-sy), oz + (z-sz));
        sx=0;sy=0;sz=0;
        stack.ImageToWorld(sx,sy,sz);
        cout<<"adjusted: "<<sx<<" "<<sy<<" "<<sz<<endl;
    
        //sprintf(buffer,"stack%i.nii.gz",l);
        //stack.Write(buffer);
        stacks.push_back(stack);
    }
    cout<<"done."<<endl;

}

void irtkReconstructionCuda::SplitImageEvenOdd(irtkRealImage image, int packages, vector<irtkRealImage>& stacks)
{
    vector<irtkRealImage> packs;
    vector<irtkRealImage> packs2;
    cout<<"Split Image Even Odd: "<<packages<<" packages."<<endl;

    stacks.clear();
    SplitImage(image,packages,packs);
    for (int i=0;i<packs.size();i++) {
        cout<<"Package "<<i<<": "<<endl;
        packs2.clear();
        SplitImage(packs[i],2,packs2);
        stacks.push_back(packs2[0]);
        stacks.push_back(packs2[1]);
    }
   
    cout<<"done."<<endl;
}

void irtkReconstructionCuda::SplitImageEvenOddHalf(irtkRealImage image, int packages, vector<irtkRealImage>& stacks, int iter)
{
    vector<irtkRealImage> packs;
    vector<irtkRealImage> packs2;
   
    cout<<"Split Image Even Odd Half "<<iter<<endl;
    stacks.clear();
    if (iter>1)
        SplitImageEvenOddHalf(image,packages,packs,iter-1);
    else
        SplitImageEvenOdd(image,packages,packs);
    for (int i=0;i<packs.size();i++) {
        packs2.clear();
        HalfImage(packs[i],packs2);
        for (int j=0; j<packs2.size(); j++)
            stacks.push_back(packs2[j]);
    }     
}


void irtkReconstructionCuda::HalfImage(irtkRealImage image, vector<irtkRealImage>& stacks)
{
    irtkRealImage tmp;
    irtkImageAttributes attr = image.GetImageAttributes();
    stacks.clear();
  
    //We would not like single slices - that is reserved for slice-to-volume
    if(attr._z>=4) {
        tmp = image.GetRegion(0,0,0,attr._x,attr._y,attr._z/2);
        stacks.push_back(tmp);
        tmp = image.GetRegion(0,0,attr._z/2,attr._x,attr._y,attr._z);
        stacks.push_back(tmp);
    }
    else
        stacks.push_back(image);
}


void irtkReconstructionCuda::PackageToVolume(vector<irtkRealImage>& stacks, vector<int> &pack_num, bool evenodd, bool half, int half_iter)
{
	PerfStats stats;
	stats.start();

    irtkImageRigidRegistrationWithPadding rigidregistration;
    irtkGreyImage t,s;
    //irtkRigidTransformation transformation;
    vector<irtkRealImage> packages;
    char buffer[256];
  
    int firstSlice = 0;
    cout<<"Package to volume: "<<endl;
    for (unsigned int i = 0; i < stacks.size(); i++) {
        cout<<"Stack "<<i<<": First slice index is "<<firstSlice<<endl;

        packages.clear();
        if (evenodd) {
            if(half)
                SplitImageEvenOddHalf(stacks[i],pack_num[i],packages,half_iter);
            else
                SplitImageEvenOdd(stacks[i],pack_num[i],packages);
        }
        else
            SplitImage(stacks[i],pack_num[i],packages);
    
        for (unsigned int j = 0; j < packages.size(); j++) {
            cout<<"Package "<<j<<" of stack "<<i<<endl;
            if (_debug) {
                sprintf(buffer,"package%i-%i.nii.gz",i,j);
                packages[j].Write(buffer);
            }
      
            t=packages[j];
            s=_reconstructed;
      
            //find existing transformation
            double x,y,z;
            x=0;y=0;z=0;
            packages[j].ImageToWorld(x,y,z);
            stacks[i].WorldToImage(x,y,z);
      
            int firstSliceIndex = round(z)+firstSlice;
            cout<<"First slice index for package "<<j<<" of stack "<<i<<" is "<<firstSliceIndex<<endl;
            //transformation = _transformations[sliceIndex];
      
            //put origin in target to zero
            irtkRigidTransformation offset;
            ResetOrigin(t,offset);
            irtkMatrix mo = offset.GetMatrix();
            irtkMatrix m = _transformations[firstSliceIndex].GetMatrix();
            m=m*mo;
            _transformations[firstSliceIndex].PutMatrix(m);

            rigidregistration.SetInput(&t, &s);
            rigidregistration.SetOutput(&_transformations[firstSliceIndex]);
            rigidregistration.GuessParameterSliceToVolume();
            if(_debug)
                rigidregistration.Write("par-packages.rreg");
            rigidregistration.Run();
      
            //undo the offset
            mo.Invert();
            m = _transformations[firstSliceIndex].GetMatrix();
            m=m*mo;
            _transformations[firstSliceIndex].PutMatrix(m);
      
            if (_debug) {
                sprintf(buffer,"transformation%i-%i.dof",i,j);
                _transformations[firstSliceIndex].irtkTransformation::Write(buffer);
            }

      
            //set the transformation to all slices of the package
            cout<<"Slices of the package "<<j<<" of the stack "<<i<<" are: ";
            for (int k = 0; k < packages[j].GetZ(); k++) {
                x=0;y=0;z=k;
                packages[j].ImageToWorld(x,y,z);
                stacks[i].WorldToImage(x,y,z);
                int sliceIndex = round(z)+firstSlice;
                cout<<sliceIndex<<" "<<endl;
    
                if(sliceIndex>=_transformations.size()) {
                    cerr<<"irtkRecnstruction::PackageToVolume: sliceIndex out of range."<<endl;
                    cerr<<sliceIndex<<" "<<_transformations.size()<<endl;
                    exit(1);
                }
    
                if(sliceIndex!=firstSliceIndex) {
                    _transformations[sliceIndex].PutTranslationX(_transformations[firstSliceIndex].GetTranslationX());
                    _transformations[sliceIndex].PutTranslationY(_transformations[firstSliceIndex].GetTranslationY());
                    _transformations[sliceIndex].PutTranslationZ(_transformations[firstSliceIndex].GetTranslationZ());
                    _transformations[sliceIndex].PutRotationX(_transformations[firstSliceIndex].GetRotationX());
                    _transformations[sliceIndex].PutRotationY(_transformations[firstSliceIndex].GetRotationY());
                    _transformations[sliceIndex].PutRotationZ(_transformations[firstSliceIndex].GetRotationZ());
                    _transformations[sliceIndex].UpdateMatrix();
                }
            }
      
      
        }
        cout<<"End of stack "<<i<<endl<<endl;
    
        firstSlice += stacks[i].GetZ();
    }
	
	stats.sample("MatchStackIntensities");
	stats.print();

}

/* end Package specific functions */

/* Utility functions */

void irtkReconstructionCuda::CropImage(irtkRealImage& image, irtkRealImage& mask)
{
    //Crops the image according to the mask

    int i, j, k;
    //ROI boundaries
    int x1, x2, y1, y2, z1, z2;

    //Original ROI
    x1 = 0;
    y1 = 0;
    z1 = 0;
    x2 = image.GetX();
    y2 = image.GetY();
    z2 = image.GetZ();

    //upper boundary for z coordinate
    int sum = 0;
    for (k = image.GetZ() - 1; k >= 0; k--) {
        sum = 0;
        for (j = image.GetY() - 1; j >= 0; j--)
            for (i = image.GetX() - 1; i >= 0; i--)
                if (mask.Get(i, j, k) > 0)
                    sum++;
        if (sum > 0)
            break;
    }
    z2 = k;

    //lower boundary for z coordinate
    sum = 0;
    for (k = 0; k <= image.GetZ() - 1; k++) {
        sum = 0;
        for (j = image.GetY() - 1; j >= 0; j--)
            for (i = image.GetX() - 1; i >= 0; i--)
                if (mask.Get(i, j, k) > 0)
                    sum++;
        if (sum > 0)
            break;
    }
    z1 = k;

    //upper boundary for y coordinate
    sum = 0;
    for (j = image.GetY() - 1; j >= 0; j--) {
        sum = 0;
        for (k = image.GetZ() - 1; k >= 0; k--)
            for (i = image.GetX() - 1; i >= 0; i--)
                if (mask.Get(i, j, k) > 0)
                    sum++;
        if (sum > 0)
            break;
    }
    y2 = j;

    //lower boundary for y coordinate
    sum = 0;
    for (j = 0; j <= image.GetY() - 1; j++) {
        sum = 0;
        for (k = image.GetZ() - 1; k >= 0; k--)
            for (i = image.GetX() - 1; i >= 0; i--)
                if (mask.Get(i, j, k) > 0)
                    sum++;
        if (sum > 0)
            break;
    }
    y1 = j;

    //upper boundary for x coordinate
    sum = 0;
    for (i = image.GetX() - 1; i >= 0; i--) {
        sum = 0;
        for (k = image.GetZ() - 1; k >= 0; k--)
            for (j = image.GetY() - 1; j >= 0; j--)
                if (mask.Get(i, j, k) > 0)
                    sum++;
        if (sum > 0)
            break;
    }
    x2 = i;

    //lower boundary for x coordinate
    sum = 0;
    for (i = 0; i <= image.GetX() - 1; i++) {
        sum = 0;
        for (k = image.GetZ() - 1; k >= 0; k--)
            for (j = image.GetY() - 1; j >= 0; j--)
                if (mask.Get(i, j, k) > 0)
                    sum++;
        if (sum > 0)
            break;
    }

    x1 = i;

    if (_debug)
        cout << "Region of interest is " << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2
             << " " << z2 << endl;

    //Cut region of interest
    image = image.GetRegion(x1, y1, z1, x2+1, y2+1, z2+1);
}

void irtkReconstructionCuda::InvertStackTransformations( vector<irtkRigidTransformation>& stack_transformations)
{
    //for each stack
    for (unsigned int i = 0; i < stack_transformations.size(); i++) {
        //invert transformation for the stacks
	stack_transformations[i].Invert();
	stack_transformations[i].UpdateParameter();
    }
}

void irtkReconstructionCuda::MaskVolume()
{
    irtkRealPixel *pr = _reconstructed.GetPointerToVoxels();
    irtkRealPixel *pm = _mask.GetPointerToVoxels();
    for (int i = 0; i < _reconstructed.GetNumberOfVoxels(); i++) {
        if (*pm == 0)
            *pr = -1;
        pm++;
        pr++;
    }
}

void irtkReconstructionCuda::MaskImage(irtkRealImage& image, double padding)
{
    if(image.GetNumberOfVoxels()!=_mask.GetNumberOfVoxels()) {
        cerr<<"Cannot mask the image - different dimensions"<<endl;
        exit(1);
    }
    irtkRealPixel *pr = image.GetPointerToVoxels();
    irtkRealPixel *pm = _mask.GetPointerToVoxels();
    for (int i = 0; i < image.GetNumberOfVoxels(); i++) {
        if (*pm == 0)
            *pr = padding;
        pm++;
        pr++;
    }
}

/// Like PutMinMax but ignoring negative values (mask)
void irtkReconstructionCuda::Rescale( irtkRealImage &img, double max)
{
    int i, n;
    irtkRealPixel *ptr, min_val, max_val;

    // Get lower and upper bound
    img.GetMinMax(&min_val, &max_val);

    n   = img.GetNumberOfVoxels();
    ptr = img.GetPointerToVoxels();
    for (i = 0; i < n; i++)
        if ( ptr[i] > 0 )
            ptr[i] = double(ptr[i]) / double(max_val) * max;
}
/* end Utility functions */
