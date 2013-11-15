/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

#ifndef RECONSTRUCTION_CUDA_CUH
#define RECONSTRUCTION_CUDA_CUH

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <vector_types.h>
#include <vector_functions.h>
#include "cutil_math.h"
#include <cufft.h>

#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/functional.h>
#include <thrust/advance.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>
#include <thrust/count.h>

#define __step  0.0001
#define PSF_SIZE 64
#define MAX_SLICES_PER_RUN 200
#define TEX_TEST 1

#if 1
#define CHECK_ERROR(function) if (cudaError_t err = cudaGetLastError()) \
{ \
  printf("CUDA Error in " #function "(), line %d: %s\n", __LINE__, cudaGetErrorString(err)); \
}
#else
#define CHECK_ERROR(function)
#endif


inline int divup(int a, int b) { return (a % b != 0) ? (a / b + 1) : (a / b); }
inline dim3 divup( uint2 a, dim3 b) { return dim3(divup(a.x, b.x), divup(a.y, b.y)); }
inline dim3 divup( dim3 a, dim3 b) { return dim3(divup(a.x, b.x), divup(a.y, b.y), divup(a.z, b.z)); }

struct POINT3D
{
    short x;
    short y;
    short z;
    double value;
};

struct Matrix4 {
    double4 data[4];
 //   inline __host__ __device__ float3 get_translation() const {
 //       return make_float3(data[0].w, data[1].w, data[2].w);
 //   }
};
//std::ostream & operator<<( std::ostream & out, const Matrix4 & m );
//Matrix4 operator*( const Matrix4 & A, const Matrix4 & B);
//Matrix4 inverse( const Matrix4 & A );


inline __host__ __device__ double dot(double4 a, double4 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline __host__ __device__ double dot(double3 a, double3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline __host__ __device__ double4 operator+( const double4 & a, const double4 & b){
	return make_double4( a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);
}

inline __host__ __device__ double3 operator+( const double3 & a, const double3 & b){
	return make_double3( a.x+b.x, a.y+b.y, a.z+b.z);
}

inline __host__ __device__ double4 operator-( const double4 & a, const double4 & b){
	return make_double4( a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
}

inline __host__ __device__ double3 operator-( const double3 & a, const double3 & b){
	return make_double3( a.x-b.x, a.y-b.y, a.z-b.z);
}

inline __host__ __device__ double4 operator*( const Matrix4 & M, const double4 & v){
	return make_double4( dot(M.data[0], v), 
						 dot(M.data[1], v), 
						 dot(M.data[2], v), 
						 dot(M.data[3], v));
}

inline __host__ __device__ double3 operator*( const Matrix4 & M, const double3 & v){
   
	return make_double3(
        dot(make_double3(M.data[0].x, M.data[0].y, M.data[0].z), v) + M.data[0].w,
        dot(make_double3(M.data[1].x, M.data[1].y, M.data[1].z), v) + M.data[1].w,
        dot(make_double3(M.data[2].x, M.data[2].y, M.data[2].z), v) + M.data[2].w);
}

inline __host__ __device__ double3 rotate( const Matrix4 & M, const double3 & v){
    return make_double3(
        dot(make_double3(M.data[0].x, M.data[0].y, M.data[0].z), v),
        dot(make_double3(M.data[1].x, M.data[1].y, M.data[1].z), v),
        dot(make_double3(M.data[2].x, M.data[2].y, M.data[2].z), v));
}


typedef std::vector<POINT3D> VOXELCOEFFS; 
typedef std::vector<std::vector<VOXELCOEFFS> > SLICECOEFFS;

struct Ref {
    Ref( void * d = NULL) : data(d) {}
    void * data;
};

struct Host {
    Host() : data(NULL) {}
    ~Host() { cudaFreeHost( data ); }

    void alloc( uint size ) { cudaHostAlloc( &data, size, cudaHostAllocDefault); }
    void * data;

};

struct Device {
    Device() : data(NULL) {}
    ~Device() { cudaFree( data ); }

    void alloc( uint size ) { cudaMalloc( &data, size ); }
    void * data;

};

struct HostDevice {
    HostDevice() : data(NULL) {}
    ~HostDevice() { cudaFreeHost( data ); }

    void alloc( uint size ) { cudaHostAlloc( &data, size,  cudaHostAllocMapped ); }
    void * getDevice() const {
        void * devicePtr;
        cudaHostGetDevicePointer(&devicePtr, data, 0);
        return devicePtr;
    }
    void * data;

	void release(){
        cudaFree(data);
        data = NULL;
    }
};

inline __device__ uint2 thr2pos2(){
#ifdef __CUDACC__
    return make_uint2( __umul24(blockDim.x, blockIdx.x) + threadIdx.x,
                       __umul24(blockDim.y, blockIdx.y) + threadIdx.y);
#else
    return make_uint2(0);
#endif
}

template <typename T>
struct Volume {
    uint3 size;
    float3 dim;
    T * data;

    Volume() { size = make_uint3(0); dim = make_float3(1); data = NULL; }
#ifdef __CUDACC__
    __device__ T operator[]( const uint3 & pos ) const {
        const T d = data[pos.x + pos.y * size.x + pos.z * size.x * size.y];
        return d; //  / 32766.0f
    }

    __device__ T v(const uint3 & pos) const {
        return operator[](pos);
    }

    __device__ T vs(const uint3 & pos) const {
        return data[pos.x + pos.y * size.x + pos.z * size.x * size.y];
    }

    __device__ void set(const uint3 & pos, const T & d ){
        data[pos.x + pos.y * size.x + pos.z * size.x * size.y] = d;
    }

    __device__ float3 pos( const uint3 & p ) const {
        return make_float3((p.x + 0.5f) * dim.x / size.x, (p.y + 0.5f) * dim.y / size.y, (p.z + 0.5f) * dim.z / size.z);
    }

  /*  __device__ float interp( const float3 & pos ) const {
#if 0   // only for testing without linear interpolation
        const float3 scaled_pos = make_float3((pos.x * size.x / dim.x) , (pos.y * size.y / dim.y) , (pos.z * size.z / dim.z) );
        return v(make_uint3(clamp(make_int3(scaled_pos), make_int3(0), make_int3(size) - make_int3(1))));

#else
        const float3 scaled_pos = make_float3((pos.x * size.x / dim.x) - 0.5f, (pos.y * size.y / dim.y) - 0.5f, (pos.z * size.z / dim.z) - 0.5f);
        const int3 base = make_int3(floorf(scaled_pos));
        const float3 factor = fracf(scaled_pos);
        const int3 lower = max(base, make_int3(0));
        const int3 upper = min(base + make_int3(1), make_int3(size) - make_int3(1));
        return (
              ((vs(make_uint3(lower.x, lower.y, lower.z)) * (1-factor.x) + vs(make_uint3(upper.x, lower.y, lower.z)) * factor.x) * (1-factor.y)
             + (vs(make_uint3(lower.x, upper.y, lower.z)) * (1-factor.x) + vs(make_uint3(upper.x, upper.y, lower.z)) * factor.x) * factor.y) * (1-factor.z)
            + ((vs(make_uint3(lower.x, lower.y, upper.z)) * (1-factor.x) + vs(make_uint3(upper.x, lower.y, upper.z)) * factor.x) * (1-factor.y)
             + (vs(make_uint3(lower.x, upper.y, upper.z)) * (1-factor.x) + vs(make_uint3(upper.x, upper.y, upper.z)) * factor.x) * factor.y) * factor.z
            ) * 0.00003051944088f;
#endif
    }

    __device__ float3 grad( const float3 & pos ) const {
        const float3 scaled_pos = make_float3((pos.x * size.x / dim.x) - 0.5f, (pos.y * size.y / dim.y) - 0.5f, (pos.z * size.z / dim.z) - 0.5f);
        const int3 base = make_int3(floorf(scaled_pos));
        const float3 factor = fracf(scaled_pos);
        const int3 lower_lower = max(base - make_int3(1), make_int3(0));
        const int3 lower_upper = max(base, make_int3(0));
        const int3 upper_lower = min(base + make_int3(1), make_int3(size) - make_int3(1));
        const int3 upper_upper = min(base + make_int3(2), make_int3(size) - make_int3(1));
        const int3 & lower = lower_upper;
        const int3 & upper = upper_lower;

        float3 gradient;

        gradient.x =
              (((vs(make_uint3(upper_lower.x, lower.y, lower.z)) - vs(make_uint3(lower_lower.x, lower.y, lower.z))) * (1-factor.x)
            + (vs(make_uint3(upper_upper.x, lower.y, lower.z)) - vs(make_uint3(lower_upper.x, lower.y, lower.z))) * factor.x) * (1-factor.y)
            + ((vs(make_uint3(upper_lower.x, upper.y, lower.z)) - vs(make_uint3(lower_lower.x, upper.y, lower.z))) * (1-factor.x)
            + (vs(make_uint3(upper_upper.x, upper.y, lower.z)) - vs(make_uint3(lower_upper.x, upper.y, lower.z))) * factor.x) * factor.y) * (1-factor.z)
            + (((vs(make_uint3(upper_lower.x, lower.y, upper.z)) - vs(make_uint3(lower_lower.x, lower.y, upper.z))) * (1-factor.x)
            + (vs(make_uint3(upper_upper.x, lower.y, upper.z)) - vs(make_uint3(lower_upper.x, lower.y, upper.z))) * factor.x) * (1-factor.y)
            + ((vs(make_uint3(upper_lower.x, upper.y, upper.z)) - vs(make_uint3(lower_lower.x, upper.y, upper.z))) * (1-factor.x)
            + (vs(make_uint3(upper_upper.x, upper.y, upper.z)) - vs(make_uint3(lower_upper.x, upper.y, upper.z))) * factor.x) * factor.y) * factor.z;

        gradient.y =
              (((vs(make_uint3(lower.x, upper_lower.y, lower.z)) - vs(make_uint3(lower.x, lower_lower.y, lower.z))) * (1-factor.x)
            + (vs(make_uint3(upper.x, upper_lower.y, lower.z)) - vs(make_uint3(upper.x, lower_lower.y, lower.z))) * factor.x) * (1-factor.y)
            + ((vs(make_uint3(lower.x, upper_upper.y, lower.z)) - vs(make_uint3(lower.x, lower_upper.y, lower.z))) * (1-factor.x)
            + (vs(make_uint3(upper.x, upper_upper.y, lower.z)) - vs(make_uint3(upper.x, lower_upper.y, lower.z))) * factor.x) * factor.y) * (1-factor.z)
            + (((vs(make_uint3(lower.x, upper_lower.y, upper.z)) - vs(make_uint3(lower.x, lower_lower.y, upper.z))) * (1-factor.x)
            + (vs(make_uint3(upper.x, upper_lower.y, upper.z)) - vs(make_uint3(upper.x, lower_lower.y, upper.z))) * factor.x) * (1-factor.y)
            + ((vs(make_uint3(lower.x, upper_upper.y, upper.z)) - vs(make_uint3(lower.x, lower_upper.y, upper.z))) * (1-factor.x)
            + (vs(make_uint3(upper.x, upper_upper.y, upper.z)) - vs(make_uint3(upper.x, lower_upper.y, upper.z))) * factor.x) * factor.y) * factor.z;

        gradient.z =
              (((vs(make_uint3(lower.x, lower.y, upper_lower.z)) - vs(make_uint3(lower.x, lower.y, lower_lower.z))) * (1-factor.x)
            + (vs(make_uint3(upper.x, lower.y, upper_lower.z)) - vs(make_uint3(upper.x, lower.y, lower_lower.z))) * factor.x) * (1-factor.y)
            + ((vs(make_uint3(lower.x, upper.y, upper_lower.z)) - vs(make_uint3(lower.x, upper.y, lower_lower.z))) * (1-factor.x)
            + (vs(make_uint3(upper.x, upper.y, upper_lower.z)) - vs(make_uint3(upper.x, upper.y, lower_lower.z))) * factor.x) * factor.y) * (1-factor.z)
            + (((vs(make_uint3(lower.x, lower.y, upper_upper.z)) - vs(make_uint3(lower.x, lower.y, lower_upper.z))) * (1-factor.x)
            + (vs(make_uint3(upper.x, lower.y, upper_upper.z)) - vs(make_uint3(upper.x, lower.y, lower_upper.z))) * factor.x) * (1-factor.y)
            + ((vs(make_uint3(lower.x, upper.y, upper_upper.z)) - vs(make_uint3(lower.x, upper.y, lower_upper.z))) * (1-factor.x)
            + (vs(make_uint3(upper.x, upper.y, upper_upper.z)) - vs(make_uint3(upper.x, upper.y, lower_upper.z))) * factor.x) * factor.y) * factor.z;

        return gradient * make_float3(dim.x/size.x, dim.y/size.y, dim.z/size.z) * (0.5f * 0.00003051944088f);
    }*/
#endif
    void init(uint3 s, float3 d){
        size = s;
        dim = d;
        cudaMalloc((void **)&data, size.x * size.y * size.z * sizeof(T));
    }

    void release(){
        cudaFree(data);
        data = NULL;
    }
};



template <typename T, typename Allocator = Ref>
struct Slice : public Allocator {
    typedef T PIXEL_TYPE;
    uint2 size;

    Slice() : Allocator() { size = make_uint2(0);  }
    Slice( const uint2 & s ) { alloc(s); }

    void alloc( const uint2 & s ){
        if(s.x == size.x && s.y == size.y)
            return;
        Allocator::alloc( s.x * s.y * sizeof(T) );
        size = s;
    }
#ifdef __CUDACC__
    __device__ T & el(){
        return operator[](thr2pos2());
    }

    __device__ const T & el() const {
        return operator[](thr2pos2());
    }

    __device__ T & operator[](const uint2 & pos ){
        return static_cast<T *>(Allocator::data)[pos.x + size.x * pos.y];
    }

    __device__ const T & operator[](const uint2 & pos ) const {
        return static_cast<const T *>(Allocator::data)[pos.x + size.x * pos.y];
    }
#endif
    Slice<T> getDeviceSlice() {
        return Slice<T>(size, Allocator::getDevice());
    }

    operator Slice<T>() {
        return Slice<T>(size, Allocator::data);
    }

    template <typename A1>
    Slice<T, Allocator> & operator=( const Slice<T, A1> & other ){
        Slice_copy(*this, other, size.x * size.y * sizeof(T));
        return *this;
    }

    T * data() {
        return static_cast<T *>(Allocator::data);
    }

    const T * data() const {
        return static_cast<const T *>(Allocator::data);
    }
};

template <typename T>
struct Slice<T, Ref> : public Ref {
    typedef T PIXEL_TYPE;
    uint2 size;

    Slice() { size = make_uint2(0,0); }
    Slice( const uint2 & s, void * d ) : Ref(d), size(s) {}
#ifdef __CUDACC__
    __device__ T & el(){
        return operator[](thr2pos2());
    }

    __device__ const T & el() const {
        return operator[](thr2pos2());
    }

    __device__ T & operator[](const uint2 & pos ){
        return static_cast<T *>(Ref::data)[pos.x + size.x * pos.y];
    }

    __device__ const T & operator[](const uint2 & pos ) const {
        return static_cast<const T *>(Ref::data)[pos.x + size.x * pos.y];
    }
#endif

    T * data() {
        return static_cast<T *>(Ref::data);
    }

    const T * data() const {
        return static_cast<const T *>(Ref::data);
    }

};



struct Reconstruction {

	Reconstruction();
	~Reconstruction();
    Volume<double> reconstructed_;
	Volume<double> reconstructed_volWeigths;
	Volume<double> mask_;
	Volume<double> maskC_;
	Volume<double> addon_;
	Volume<double> bias_;
	Volume<double> confidence_map_;
	Volume<double> volume_weights_;
	Volume<int> sliceVoxel_count_;

	Volume<double> v_PSF_sums_;
	Volume<double> v_slices;
	Volume<double> v_bias;
	Volume<double> v_weights;
	Volume<double> v_simulated_weights;
	Volume<double> v_simulated_slices;
	Volume<char> v_simulated_inside;
	Volume<double> v_wresidual;
	Volume<double> v_wb;
	Volume<double> v_buffer;

	Volume<float> PSFlut;

	int* d_slice_sicesX;
	int* d_slice_sicesY;
	double* d_scales;
	double* d_slice_weights;

	float3* d_sliceDims;
	float reconstructedVoxelSize;

	Matrix4* d_slicesI2W; //unfortunately too large for constant mem
	Matrix4* d_slicesW2I;
	Matrix4* d_slicesTransformation;
	Matrix4* d_slicesInvTransformation;

	Matrix4 _ri2w;
	Matrix4 _rw2i;

	thrust::device_vector<double> scale_;
	std::vector<double> h_scales;
	std::vector<double> h_slices_weights;
	thrust::device_vector<double> d_scale;
	thrust::device_vector<double> d_slices_weights;

	std::vector<Matrix4> sliceI2W;
	std::vector<Matrix4> sliceW2I;
	std::vector<Matrix4> sliceTransformation;
	std::vector<Matrix4> sliceInvTransformation;
	Matrix4 reconstructedW2I;
	Matrix4 reconstructedI2W;

#if !TEX_TEST
	void UpdateSliceCoeff(std::vector<SLICECOEFFS> _volcoeffs, double* volWeights);
	std::vector<Volume<POINT3D>> volcoeffs_; //this is not a good idea for full res, todo improve, e.g. per slice
	unsigned int maxVolcoDepth;
	Volume<int> v_volcoeffsOfsSlices_X;
	Volume<int> v_volcoeffsOfsSlices_Y;
	POINT3D* volcoeffsarray_;
#endif

	uint2 maxSliceDimension;
	bool* slicesInit;
	unsigned int num_slices_;

	//debugFuctions
	void debugWeights(double* weights);
	void debugBias(double* bias);
	void debugNormalizeBias(double* nbias);
	void debugSmoothMask(double* smoothMask);
	void debugSimslices(double* simslices);
	void debugSimweights(double* simweights);
	void getSlicesVol_debug(double* h_imdata);
	void syncCPU(double* reconstructed);
	void InitializeEMValues();

	//init/update functions
	void initStorageVolumes(uint3 size, float3 dim);
	void FillSlices(double* sdata,  std::vector<int> sizesX, std::vector<int> sizesY);
	void generatePSFVolume(double* CPUPSF, float3 sliceVoxelDim);
	void setSliceDims(std::vector<float3> slice_dims);
	void SetSliceMatrices(std::vector<Matrix4>& matsI2W, std::vector<Matrix4>& matsW2I, Matrix4 reconI2W, Matrix4 reconW2I);
	void UpdateSliceWeights(std::vector<double> slices_weights);
    void InitReconstructionVolume(uint3 s, float3 dim, double* data, double sigma_bias); // allocates the volume and image data on the device
	void InitSlices(int num_slices);
	void UpdateReconstructed(const uint3 vsize, double* data); //temporary for CPU GPU sync
	void SyncConfidenceMapAddon(double* cmdata, double* addondata);
	void UpdateScaleVector(std::vector<double> scales, std::vector<double> slices_weights);//{h_scales = scales; /*thrust::copy(scales.begin(), scales.end(), scale_.begin());*/};
	void CalculateScaleVector();
	void setMask(uint3 s, float3 dim, double* data, double sigma_bias);

	//calculation functions
	void NormaliseBias(int iter, double sigma_bias);
	void EStep(double _m, double _sigma, double _mix, std::vector<double>& slice_potential);
	void MStep(int iter,  double _step, double& _sigma, double& _mix, double& _m);
	void SimulateSlices(std::vector<bool>& slice_inside);
	void SyncSimulatedCPU(int slice_idx, double* ssdata, double* swdata, double* sidata);
	void InitializeRobustStatistics(double& _sigma);
	void CorrectBias(double sigma_bias, bool _global_bias_correction);
	void Superresolution(int iter, std::vector<double> _slice_weight, bool _adaptive, double alpha, 
		double _min_intensity, double _max_intensity, double delta, double lambda,  bool _global_bias_correction, double sigma_bias,
		double _low_intensity_cutoff);

	void GaussianReconstruction(std::vector<int>& voxel_num);
	
};



#endif //RECONSTRUCTION_CUDA_CUH