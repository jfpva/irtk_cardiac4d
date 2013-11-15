/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

#include "reconstruction_cuda.cuh"
#include "gaussFilterConvolution.cuh"
#define _USE_MATH_DEFINES
#include <math.h>
#include <thrust/inner_product.h>

using namespace thrust;

__constant__ int d_directions[13][3];
__constant__ double d_factor[13];
__constant__ Matrix4 d_reconstructedW2I;
__constant__ Matrix4 d_reconstructedI2W;

texture<float, 3, cudaReadModeElementType > psfTex;

__forceinline__ __device__ float sq( const float x ){
	return x*x;
}

inline __host__ __device__ inline double G(double x,double s)
{
	return __step*exp(-x*x/(2*s))/(sqrt(6.28*s)); //todo define step in symbol
}

inline __host__ __device__ inline double M(double m)
{
	return m*__step;
}

void checkGPUMemory()
{
	size_t free_byte;
	size_t total_byte ;
	cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;

	if ( cudaSuccess != cuda_status ){
		printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
		//exit(1);
	}

	double free_db = (double)free_byte ;
	double total_db = (double)total_byte ;
	double used_db = total_db - free_db ;
	printf("GPU memory usage: \nused = %f, free = %f MB, total = %f MB\n",
		used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
}

__host__ __device__ int round_(double x)
{
	return (x > 0) ? (int)(x + 0.5) : (int)(x - 0.5);
}

__device__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
		(unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
			__double_as_longlong(val +
			__longlong_as_double(assumed)));
	} while (assumed != old);
	return __longlong_as_double(old);
}


///////////////////////////////test area PSF texture

#define EXTENT 4 

#if TEX_TEST
//in slice coords
__device__ double calcPSF(float3 sPos, float3 dim)
{

	//z must be normal to the transfomed (registered) slice
	// center_slice = sliceImageToWorld*transformation * (x,y,0)'
	// dx, dy, dz = slice pixel sizes
	// need to be centered 
	// for all number of voxels in each direction
	//the ROI is 2*voxel dimension
	//int xDim = round(2 * dx / size);
	// int yDim = round(2 * dy / size);
	//int zDim = round(2 * dz / size);
	//reconstructor->_reconstructed.WorldToImage
	const double sigmax = 1.2 * dim.x / 2.3548;
	const double sigmay = 1.2 * dim.y / 2.3548;
	const double sigmaz = dim.z / 2.3548;

	return exp(-sPos.x * sPos.x / (2 * sigmax * sigmax) - sPos.y * sPos.y / (2 * sigmay * sigmay)
		- sPos.z * sPos.z / (2 * sigmaz * sigmaz));

}

//TODO add extend as quality_factor parameter

__global__ void gaussianReconstructionKernel3D_tex(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict slices, double* __restrict bias2D, double* __restrict scales, Volume<double> reconstructed, 
	Volume<double> reconstructed_volWeigths, Volume<double> v_PSF_sums,  Volume<int> sliceVoxel_count, double* __restrict mask, uint2 vSize,
	Matrix4* __restrict sliceI2W, Matrix4* __restrict sliceW2I, float3* __restrict d_slicedim, float reconVSize, int step)
{

	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z + step);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;
	unsigned int idx = pos.x + pos.y*vSize.x + pos.z*vSize.x*vSize.y;
	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];
	double s = slices[idx];
	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0))
		return;

	double b = bias2D[idx];
	double scale = scales[pos.z];

	float3 sliceDim = d_slicedim[pos.z];

	double sliceVal = s * exp(-b) * scale;

	//first slice->trans->volume = centerpos
	//for 8 neighbors
	// centerpos+xyz ->world->invTran->slicepos
	//(slicepos - centerpos)*slicedim = pos in PSF relative to center
	//sample PSFtex = psfval
	//insert sliceVal*psfval at centerpos+xyz
	double3 slicePos = make_double3((float)pos.x+0.5, (float)pos.y+0.5, 0.5);
	double3 wpos = sliceI2W[pos.z]*slicePos;
	double3 volumePos = d_reconstructedW2I * wpos;
	double3 PSFcenter = make_double3(PSF_SIZE/2.0 + 0.5,PSF_SIZE/2.0 + 0.5,PSF_SIZE/2.0 + 0.5);
	//volumePos = make_double3(round_(volumePos.x), round_(volumePos.y), round_(volumePos.z)); 

	double sume = 0;
	for(int z = -EXTENT; z < EXTENT; z++)
	{
		for(int y = -EXTENT; y < EXTENT; y++)
		{
			for(int x = -EXTENT; x < EXTENT; x++)
			{
				double3 ofsPos = make_double3(volumePos.x+x, volumePos.y+y, volumePos.z+z);
				wpos = d_reconstructedI2W*ofsPos;
				double3 nslicePos = sliceW2I[pos.z]*wpos;
				double3 psfCOfs = nslicePos-slicePos;
				psfCOfs = PSFcenter + make_double3(psfCOfs.x*sliceDim.x,psfCOfs.y*sliceDim.y,psfCOfs.z*sliceDim.z);
				//psfCOfs = PSFcenter + make_double3(psfCOfs.x,psfCOfs.y,psfCOfs.z);
				double psfval = tex3D(psfTex, psfCOfs.x/(float)PSF_SIZE, psfCOfs.y/(float)PSF_SIZE, psfCOfs.z/(float)PSF_SIZE);
				//double psfval = calcPSF(make_float3(psfCOfs.x,psfCOfs.y,psfCOfs.z), sliceDim);

				uint3 apos = make_uint3(round_(ofsPos.x),round_(ofsPos.y),round_(ofsPos.z)); //NN
				if(apos.x < reconstructed.size.x && apos.y < reconstructed.size.y && apos.z < reconstructed.size.z && 
					apos.x >= 0 && apos.y >= 0 && apos.z >= 0 /*&& mask[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y] == 1*/)
				{
					sume += psfval;
				}
			}
		}
	}
	//storage for PSF norm
	v_PSF_sums.set(pos, sume);

	if(sume <= 0)
		return;

	bool addvoxNum = false;
	for(int z = -EXTENT; z < EXTENT; z++)
	{
		for(int y = -EXTENT; y < EXTENT; y++)
		{
			for(int x = -EXTENT; x < EXTENT; x++)
			{
				double3 ofsPos = make_double3(volumePos.x+x, volumePos.y+y, volumePos.z+z);
				wpos = d_reconstructedI2W*ofsPos;
				double3 nslicePos = sliceW2I[pos.z]*wpos;
				double3 psfCOfs = nslicePos-slicePos;
				psfCOfs = PSFcenter + make_double3(psfCOfs.x*sliceDim.x,psfCOfs.y*sliceDim.y,psfCOfs.z*sliceDim.z);
				//psfCOfs = PSFcenter + make_double3(psfCOfs.x,psfCOfs.y,psfCOfs.z);
				double psfval = tex3D(psfTex, psfCOfs.x/(float)PSF_SIZE, psfCOfs.y/(float)PSF_SIZE, psfCOfs.z/(float)PSF_SIZE)/sume;
				//double psfval = calcPSF(make_float3(psfCOfs.x,psfCOfs.y,psfCOfs.z), sliceDim)/sume;
				
				uint3 apos = make_uint3(round_(ofsPos.x),round_(ofsPos.y),round_(ofsPos.z)); //NN
				if(apos.x < reconstructed.size.x && apos.y < reconstructed.size.y && apos.z < reconstructed.size.z && 
					apos.x >= 0 && apos.y >= 0 && apos.z >= 0 /*&& mask[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y] == 1*/)
				{
					atomicAdd(&(reconstructed.data[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y]), 
						psfval*sliceVal);
					atomicAdd(&(reconstructed_volWeigths.data[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y]), 
						psfval);//sume);
					addvoxNum = true;
				}

				//linear extrapolation works but no huge difference
				/*int3 ll = make_int3(floor(ofsPos.x/reconVSize),floor(ofsPos.y/reconVSize),floor(ofsPos.z/reconVSize));
				for (int l = ll.x; l <= ll.x + 1; l++)
				{
				for (int m = ll.y; m <= ll.y + 1; m++)
				{
				for (int n = ll.z; n <= ll.z + 1; n++)
				{
				uint3 apos = make_uint3(l,m,n); //linear extrapolate
				//check if inside 
				if(apos.x < reconstructed.size.x && apos.y < reconstructed.size.y && apos.z < reconstructed.size.z && 
				apos.x >= 0 && apos.y >= 0 && apos.z >= 0 && mask[apos] == 1)
				{
				double weight = (fabs((float)l - ofsPos.x)) * (fabs((float)m - ofsPos.y)) * 
				(fabs((float)n - ofsPos.z));
				//double psfvalw = psfval * weight / sume;
				atomicAdd(&(reconstructed.data[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y]), 
				weight*psfval*sliceVal);
				atomicAdd(&(reconstructed_count.data[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y]), 
				psfval*weight);
				}
				}
				}
				}*/

			}
		}
	}

	if(addvoxNum)
	{
		atomicAdd(&(sliceVoxel_count.data[idx]), 1);
	}
}


__global__ void simulateSlicesKernel3D_tex(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict slices,  Volume<double> simslices, Volume<double> simweights, Volume<char> siminside,
	Volume<double> reconstructed, double* __restrict v_PSF_sums, double* __restrict mask, uint2 vSize,
	Matrix4* __restrict sliceI2W, Matrix4* __restrict sliceW2I, float3* __restrict d_slicedim, float reconVSize, int step)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z + step);

	if(pos.z >= numSlices || pos.z < 0)
	return;
	unsigned int idx = pos.x + pos.y*vSize.x + pos.z*vSize.x*vSize.y;
	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];
	double s = slices[idx];
	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0))
	return;

	double simulated_sliceV = 0;
	bool slice_inside = false;
	double weight = 0;
	float3 sliceDim = d_slicedim[pos.z];

	double3 slicePos = make_double3((float)pos.x+0.5, (float)pos.y+0.5, 0.5);
	double3 wpos = sliceI2W[pos.z]*slicePos;
	double3 volumePos = d_reconstructedW2I * wpos;
	double3 PSFcenter = make_double3(PSF_SIZE/2.0 + 0.5,PSF_SIZE/2.0 + 0.5,PSF_SIZE/2.0 + 0.5);
	//volumePos = make_double3(round_(volumePos.x), round_(volumePos.y), round_(volumePos.z)); 

	double sume = v_PSF_sums[idx];

	if(sume <= 0)
		return;

	for(int z = -EXTENT; z < EXTENT; z++)
	{
		for(int y = -EXTENT; y < EXTENT; y++)
		{
			for(int x = -EXTENT; x < EXTENT; x++)
			{
				double3 ofsPos = make_double3(volumePos.x+x, volumePos.y+y, volumePos.z+z);
				wpos = d_reconstructedI2W*ofsPos;
				double3 nslicePos = sliceW2I[pos.z]*wpos;
				double3 psfCOfs = nslicePos-slicePos;
				psfCOfs = PSFcenter + make_double3(psfCOfs.x*sliceDim.x,psfCOfs.y*sliceDim.y,psfCOfs.z*sliceDim.z);
				//psfCOfs = PSFcenter + make_double3(psfCOfs.x,psfCOfs.y,psfCOfs.z);
				double psfval = tex3D(psfTex, psfCOfs.x/(float)PSF_SIZE, psfCOfs.y/(float)PSF_SIZE, psfCOfs.z/(float)PSF_SIZE)/sume;
				//double psfval = calcPSF(make_float3(psfCOfs.x,psfCOfs.y,psfCOfs.z), sliceDim)/sume;

				uint3 apos = make_uint3(round_(ofsPos.x),round_(ofsPos.y),round_(ofsPos.z)); //NN
				if(apos.x < reconstructed.size.x && apos.y < reconstructed.size.y && apos.z < reconstructed.size.z && 
					apos.x >= 0 && apos.y >= 0 && apos.z >= 0  /*&& mask[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y] == 1*/)
				{
					simulated_sliceV += psfval * reconstructed[apos];
					weight += psfval;

					if(mask[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y] == 1)
					{
						slice_inside = true; 
					}
				}
			}
		}
	}

	siminside.set(pos, (char)slice_inside);

	if( weight > 0 ) 
	{
		simslices.set(pos,simulated_sliceV/weight);
		simweights.set(pos,weight);
	}
	else
	{
		simslices.set(pos,simulated_sliceV);
	}

}


__global__ void SuperresolutionKernel3D_tex(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	Volume<double> slices, Volume<double> bias, Volume<double> weights, Volume<double> simslices, 
	double* __restrict slice_weights, double* __restrict scales, double* __restrict v_PSF_sums, double* __restrict mask,
	Volume<double> addon, double* confidence_map, 
	Matrix4* __restrict sliceI2W, Matrix4* __restrict sliceW2I, float3* __restrict d_slicedim, float reconVSize, int step)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z + step);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;
	unsigned int idx = pos.x + pos.y*slices.size.x + pos.z*slices.size.x*slices.size.y;
	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];

	double s = slices[pos];
	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0))
		return;

	double b = bias[pos];
	double w = weights[pos];
	double ss = simslices[pos];
	double slice_weight = slice_weights[pos.z];
	double scale = scales[pos.z];

	double sliceVal = s * exp(-b) * scale;
	if ( ss > 0 )
		sliceVal -= ss;
	else
		sliceVal = 0;

	float3 sliceDim = d_slicedim[pos.z];

	double3 slicePos = make_double3((float)pos.x+0.5, (float)pos.y+0.5, 0.5);
	double3 wpos = sliceI2W[pos.z]*slicePos;
	double3 volumePos = d_reconstructedW2I * wpos;
	double3 PSFcenter = make_double3(PSF_SIZE/2.0 + 0.5,PSF_SIZE/2.0 + 0.5,PSF_SIZE/2.0 + 0.5);
	//volumePos = make_double3(round_(volumePos.x), round_(volumePos.y), round_(volumePos.z)); 

	double sume = v_PSF_sums[idx];
	
	if(sume <= 0)
		return;

	for(int z = -EXTENT; z < EXTENT; z++)
	{
		for(int y = -EXTENT; y < EXTENT; y++)
		{
			for(int x = -EXTENT; x < EXTENT; x++)
			{
				double3 ofsPos = make_double3(volumePos.x+x, volumePos.y+y, volumePos.z+z);
				wpos = d_reconstructedI2W*ofsPos;
				double3 nslicePos = sliceW2I[pos.z]*wpos;
				double3 psfCOfs = nslicePos-slicePos;
				psfCOfs = PSFcenter + make_double3(psfCOfs.x*sliceDim.x,psfCOfs.y*sliceDim.y,psfCOfs.z*sliceDim.z);
				//psfCOfs = PSFcenter + make_double3(psfCOfs.x,psfCOfs.y,psfCOfs.z);
				double psfval = tex3D(psfTex, psfCOfs.x/(float)PSF_SIZE, psfCOfs.y/(float)PSF_SIZE, psfCOfs.z/(float)PSF_SIZE)/sume;
				//double psfval = calcPSF(make_float3(psfCOfs.x,psfCOfs.y,psfCOfs.z), sliceDim)/sume;

				uint3 apos = make_uint3(round_(ofsPos.x),round_(ofsPos.y),round_(ofsPos.z)); //NN
				if(apos.x < addon.size.x && apos.y < addon.size.y && apos.z < addon.size.z && 
					apos.x >= 0 && apos.y >= 0 && apos.z >= 0 /*&& mask[apos.x + apos.y*addon.size.x + apos.z*addon.size.x*addon.size.y] == 1*/)
				{
					atomicAdd(&(addon.data[apos.x + apos.y*addon.size.x + apos.z*addon.size.x*addon.size.y]), 
						psfval * w * slice_weight * sliceVal);
					atomicAdd(&(confidence_map[apos.x + apos.y*addon.size.x + apos.z*addon.size.x*addon.size.y]), 
						psfval * w * slice_weight);
				}
			}
		}
	}

}


__global__ void normalizeBiasKernel3D_tex(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	Volume<double> slices, Volume<double> bias2D, double* __restrict scales, double* __restrict v_PSF_sums, double* __restrict mask,
	Matrix4* __restrict sliceI2W, Matrix4* __restrict sliceW2I, float3* __restrict d_slicedim, float reconVSize, 
	Volume<double> bias, Volume<double> reconVolWeights, int step)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z + step);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;
	unsigned int idx = pos.x + pos.y*slices.size.x + pos.z*slices.size.x*slices.size.y;
	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];

	double s = slices[pos];
	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0))
		return;

	double nbias = bias2D[pos];
	double scale = scales[pos.z];
	if(s > -1.0 && scale > 0)
	{
		nbias -= log(scale);
	}

	float3 sliceDim = d_slicedim[pos.z];

	double3 slicePos = make_double3((float)pos.x+0.5, (float)pos.y+0.5, 0.5);
	double3 wpos = sliceI2W[pos.z]*slicePos;
	double3 volumePos = d_reconstructedW2I * wpos;
	double3 PSFcenter = make_double3(PSF_SIZE/2.0 + 0.5,PSF_SIZE/2.0 + 0.5,PSF_SIZE/2.0 + 0.5);
	//volumePos = make_double3(round_(volumePos.x), round_(volumePos.y), round_(volumePos.z)); 

	double sume = v_PSF_sums[idx];
	
	if(sume <= 0)
		return;

	for(int z = -EXTENT; z < EXTENT; z++)
	{
		for(int y = -EXTENT; y < EXTENT; y++)
		{
			for(int x = -EXTENT; x < EXTENT; x++)
			{
				double3 ofsPos = make_double3(volumePos.x+x, volumePos.y+y, volumePos.z+z);
				wpos = d_reconstructedI2W*ofsPos;
				double3 nslicePos = sliceW2I[pos.z]*wpos;
				double3 psfCOfs = nslicePos-slicePos;
				psfCOfs = PSFcenter + make_double3(psfCOfs.x*sliceDim.x,psfCOfs.y*sliceDim.y,psfCOfs.z*sliceDim.z);
				double psfval = tex3D(psfTex, psfCOfs.x/(float)PSF_SIZE, psfCOfs.y/(float)PSF_SIZE, psfCOfs.z/(float)PSF_SIZE)/sume;
				//double psfval = calcPSF(make_float3(psfCOfs.x,psfCOfs.y,psfCOfs.z), sliceDim)/sume;

				uint3 apos = make_uint3(round_(ofsPos.x),round_(ofsPos.y),round_(ofsPos.z)); //NN
				if(apos.x < bias.size.x && apos.y < bias.size.y && apos.z < bias.size.z && 
					apos.x >= 0 && apos.y >= 0 && apos.z >= 0 /*&& mask[apos.x + apos.y*bias.size.x + apos.z*bias.size.x*bias.size.y] == 1*/)
				{
					atomicAdd(&(bias.data[apos.x + apos.y*bias.size.x + apos.z*bias.size.x*bias.size.y]), 
						psfval*nbias);
				}
			}
		}
	}

}

#endif //TEX_TEST
///////////////////////////////test area PSF texture end

Reconstruction::Reconstruction()
{
	maxSliceDimension.x = 0; 
	maxSliceDimension.y = 0;

	int h_directions[][3] = {
		{ 1, 0,-1 },
		{ 0, 1,-1 },
		{ 1, 1,-1 },
		{ 1,-1,-1 },
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 1, 1, 0 },
		{ 1,-1, 0 },
		{ 1, 0, 1 },
		{ 0, 1, 1 },
		{ 1, 1, 1 },
		{ 1,-1, 1 },
		{ 0, 0, 1 }
	};
	// this seems to be super constant -> constant mem
	double factor[13];
	for (int i = 0; i < 13; i++) {
		factor[i] = 0;
	}
	for (int i = 0; i < 13; i++) {
		for (int j = 0; j < 3; j++)
		{
			factor[i] += fabs((double)(h_directions[i][j]));
		}
		factor[i] = 1.0 / factor[i];
	}

#if ! TEX_TEST
	volcoeffsarray_ = NULL;
#endif

	checkCudaErrors(cudaMemcpyToSymbol(d_factor, (void*)factor, 13*sizeof(double)));
	checkCudaErrors(cudaMemcpyToSymbol(d_directions, (void*)h_directions, 3*13*sizeof(int)));
	
	d_slice_sicesX = NULL;
	d_slice_sicesY = NULL;
	d_scales = NULL;
	d_slice_weights = NULL;
	
	d_slicesI2W = NULL;
	d_sliceDims = NULL;
	d_slicesW2I = NULL;
}

void Reconstruction::generatePSFVolume(double* CPUPSF, float3 sliceVoxelDim)
{
	printf("PSF lut tex %f %f %f \n", sliceVoxelDim.x, sliceVoxelDim.y, sliceVoxelDim.z);
	float *d_psf = NULL;
	checkCudaErrors(cudaMalloc((void **) &d_psf, PSF_SIZE*PSF_SIZE*PSF_SIZE*sizeof(float)));

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaArray *cuArray;

	cudaExtent asize;
	asize.height = PSF_SIZE;
	asize.width = PSF_SIZE;
	asize.depth = PSF_SIZE;
	checkCudaErrors(cudaMalloc3DArray(&cuArray, &channelDesc, asize));

	dim3 blockSize3 = dim3(8,8,8);
	dim3 gridSize3 = divup(dim3(PSF_SIZE, PSF_SIZE, PSF_SIZE), blockSize3);

	//TODO get slice dim -- must be euqal for all slices
	PSFlut.init(make_uint3(PSF_SIZE,PSF_SIZE,PSF_SIZE), sliceVoxelDim);
	float* h_PSF = new float[PSF_SIZE*PSF_SIZE*PSF_SIZE];
	for(int i = 0; i < PSF_SIZE*PSF_SIZE*PSF_SIZE; i++)
		h_PSF[i] = CPUPSF[i];

	cudaMemcpy3DParms copyParams = {0};
	copyParams.srcPtr   = make_cudaPitchedPtr((void *)h_PSF, PSF_SIZE*sizeof(float), PSF_SIZE, PSF_SIZE);//make_cudaPitchedPtr((void *)PSFlut.data, PSF_SIZE*sizeof(float), PSF_SIZE, PSF_SIZE);
	copyParams.dstArray = cuArray;
	copyParams.extent   = asize;
	copyParams.kind     = cudaMemcpyHostToDevice;
	checkCudaErrors(cudaMemcpy3D(&copyParams));

	psfTex.addressMode[0] = cudaAddressModeClamp;
	psfTex.addressMode[1] = cudaAddressModeClamp;
	psfTex.addressMode[3] = cudaAddressModeClamp;
	psfTex.filterMode = cudaFilterModeLinear;
	psfTex.normalized = true;

	checkCudaErrors(cudaBindTextureToArray(psfTex, cuArray, channelDesc));
	CHECK_ERROR(cudaBindTextureToArray);
}


template <typename T>
__global__ void initVolume( Volume<T> volume, const T val ){
	uint3 pos = make_uint3(thr2pos2());
	for(pos.z = 0; pos.z < volume.size.z; ++pos.z)
	{
		if(pos.x < volume.size.x && pos.y < volume.size.y && pos.z < volume.size.z &&
			pos.x >= 0 && pos.y >= 0 && pos.z >= 0)
			volume.set(pos, val);
	}
}

template <class T>
void initMem(T* mem, unsigned int N, T val)
{
	thrust::device_ptr<T> dev_p(mem);
	thrust::fill(dev_p, dev_p + N, val);
}


void Reconstruction::setSliceDims(std::vector<float3> slice_dims)
{
	if(d_sliceDims == NULL)  checkCudaErrors(cudaMalloc((void**)&d_sliceDims, slice_dims.size()*sizeof(float3))); 
	checkCudaErrors(cudaMemcpy(d_sliceDims, &slice_dims[0], slice_dims.size()*sizeof(float3), cudaMemcpyHostToDevice));
}


void Reconstruction::SetSliceMatrices(std::vector<Matrix4>& matsI2W, std::vector<Matrix4>& matsW2I, Matrix4 reconI2W, Matrix4 reconW2I)
{
	sliceW2I = matsW2I;
	sliceI2W = matsI2W;
	_rw2i = reconW2I;
	_ri2w = reconI2W;

	printf("mat slices %d %d \n", matsI2W.size(), matsW2I.size());

	checkCudaErrors(cudaMemcpyToSymbol(d_reconstructedW2I, &(reconW2I), sizeof(Matrix4)));
	checkCudaErrors(cudaMemcpyToSymbol(d_reconstructedI2W, &(reconI2W), sizeof(Matrix4)));
	CHECK_ERROR(cudaMemcpyToSymbol);

	if(d_slicesI2W == NULL) checkCudaErrors(cudaMalloc((void**)&d_slicesI2W, sliceI2W.size()*sizeof(Matrix4))); 
	checkCudaErrors(cudaMemcpy(d_slicesI2W, &sliceI2W[0], sliceI2W.size()*sizeof(Matrix4), cudaMemcpyHostToDevice));
	CHECK_ERROR(d_slicesI2W);

	//Matrix4* test = new Matrix4[sliceI2W.size()];
	//checkCudaErrors(cudaMemcpy(test, d_slicesI2W, sliceI2W.size()*sizeof(Matrix4), cudaMemcpyDeviceToHost));
	//CHECK_ERROR(SetSliceMatrices);
	/*for(int j = 0; j < sliceI2W.size(); j++)
	{
		Matrix4 si2w = test[j];
		for(int i = 0; i < 4; i++)
		{
			printf("%f %f %f %f \n", si2w.data[i].x, si2w.data[i].y, si2w.data[i].z, si2w.data[i].w); //000?
		}
	}*/

	if(d_slicesW2I == NULL) checkCudaErrors(cudaMalloc((void**)&d_slicesW2I, sliceW2I.size()*sizeof(Matrix4))); 
	checkCudaErrors(cudaMemcpy(d_slicesW2I, &sliceW2I[0], sliceW2I.size()*sizeof(Matrix4), cudaMemcpyHostToDevice));
	CHECK_ERROR(d_slicesW2I);
}

void Reconstruction::setMask(uint3 s, float3 dim, double* data, double sigma_bias)
{
	mask_.init(s,dim);
	checkCudaErrors(cudaMemcpy(mask_.data, data, s.x*s.y*s.z*sizeof(double), cudaMemcpyHostToDevice));		
	maskC_.init(s,dim);
	checkCudaErrors(cudaMemcpy(maskC_.data, data, s.x*s.y*s.z*sizeof(double), cudaMemcpyHostToDevice));	

	/*dim3 blockSize3 = dim3(8,8,8);
	dim3 gridSize3 = divup(dim3(maskC_.size.x, maskC_.size.y, maskC_.size.z), blockSize3);

	Volume<double> mbuf;
	mbuf.init(maskC_.size, maskC_.dim);

	//int klength = max(5, (unsigned int)ceil(sigma_bias * 5) * 2 + 1);
	int klength = 2*round_(4*sigma_bias)+1;
	GaussianConvolutionKernel3D<<<gridSize3,blockSize3>>>(mask_.data, maskC_, klength, sigma_bias, 0);
	GaussianConvolutionKernel3D<<<gridSize3,blockSize3>>>(maskC_.data, mbuf, klength, sigma_bias, 1);
	GaussianConvolutionKernel3D<<<gridSize3,blockSize3>>>(mbuf.data, maskC_, klength, sigma_bias, 2);

	mbuf.release();*/

	//filterGaussFFT(maskC_);
	FilterGauss(maskC_.data, maskC_.data, maskC_.size.x, maskC_.size.y, maskC_.size.z, maskC_.size.x, maskC_.size.y, 1, sigma_bias );
}


void Reconstruction::InitReconstructionVolume(uint3 s, float3 dim, double *data, double sigma_bias) 
{
	//printf("InitReconstructionVolume %d %d %d \n", s.x, s.y, s.z);
	reconstructed_.init(s,dim);
	bias_.init(s,dim);
	reconstructed_volWeigths.init(s,dim);
	addon_.init(s,dim);
	confidence_map_.init(s,dim);

	dim3 block(32,16);
	dim3 grid = divup(dim3(reconstructed_.size.x, reconstructed_.size.y), block);
	unsigned int N = addon_.size.x*addon_.size.y*addon_.size.z;
	cudaMemset(reconstructed_.data, 0, N*sizeof(double));

	if(data != NULL) checkCudaErrors(cudaMemcpy(reconstructed_.data, data, s.x*s.y*s.z*sizeof(double), cudaMemcpyHostToDevice));

	cudaMemset(reconstructed_volWeigths.data, 0, N*sizeof(double));
	cudaMemset(addon_.data, 0, N*sizeof(double));
	cudaMemset(confidence_map_.data, 0, N*sizeof(double));
	cudaMemset(bias_.data, 0, N*sizeof(double));

}

Reconstruction::~Reconstruction()
{
	reconstructed_.release();
	bias_.release();
	reconstructed_volWeigths.release();
	addon_.release();
	confidence_map_.release();

	sliceVoxel_count_.release();
	v_slices.release();
	v_bias.release();
	v_weights.release();
	v_simulated_weights.release();
	v_simulated_slices.release();
	v_wresidual.release();
	v_wb.release();
	v_buffer.release();
#if !TEX_TEST
	v_volcoeffsOfsSlices_X.release();
	v_volcoeffsOfsSlices_Y.release();
#endif
	v_PSF_sums_.release();
}

void Reconstruction::UpdateSliceWeights(std::vector<double> slices_weights)
{
	h_slices_weights = slices_weights;
	checkCudaErrors(cudaMemcpy(d_slice_weights, &slices_weights[0], slices_weights.size()*sizeof(double), cudaMemcpyHostToDevice));
}

void Reconstruction::UpdateScaleVector(std::vector<double> scales, std::vector<double> slices_weights)
{
	//obsolete -- scales are calcualted here
	h_scales = scales;
	h_slices_weights = slices_weights;
	if(d_scales == NULL) checkCudaErrors(cudaMalloc((void **)&d_scales, num_slices_*sizeof(double))); 
	if(d_slice_weights == NULL) checkCudaErrors(cudaMalloc((void **)&d_slice_weights, num_slices_*sizeof(double))); 

	checkCudaErrors(cudaMemcpy(d_scales, &scales[0], scales.size()*sizeof(double), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_slice_weights, &slices_weights[0], slices_weights.size()*sizeof(double), cudaMemcpyHostToDevice));
}


void Reconstruction::initStorageVolumes(uint3 size, float3 dim)
{
	v_slices.init(size, dim);
	v_bias.init(size, dim);
	v_weights.init(size, dim);
	v_simulated_weights.init(size, dim);
	v_simulated_slices.init(size, dim);
	v_wresidual.init(size, dim);
	v_wb.init(size, dim);
	v_buffer.init(size, dim);
	sliceVoxel_count_.init(size, dim);
	v_simulated_inside.init(size, dim);
#if !TEX_TEST
	v_volcoeffsOfsSlices_X.init(size, dim);
	v_volcoeffsOfsSlices_Y.init(size, dim);
#endif
	v_PSF_sums_.init(size, dim);
	unsigned int N = size.x*size.y*size.z;

	cudaMemset(v_slices.data, 0, N*sizeof(double));
	//initMem<double>(v_slices.data, N, -1.0); //not necessary, init later to padding
	cudaMemset(v_bias.data, 0, N*sizeof(double));
	cudaMemset(v_weights.data, 0, N*sizeof(double));
	cudaMemset(v_simulated_weights.data, 0, N*sizeof(double));
	cudaMemset(v_simulated_slices.data, 0, N*sizeof(double));
	cudaMemset(v_wresidual.data, 0, N*sizeof(double));
	cudaMemset(v_wb.data, 0, N*sizeof(double));
	cudaMemset(v_buffer.data, 0, N*sizeof(double));
	cudaMemset(v_simulated_inside.data, 0, N*sizeof(char));
	cudaMemset(v_PSF_sums_.data, 0, N*sizeof(double));
	
#if !TEX_TEST
	cudaMemset(v_volcoeffsOfsSlices_X.data, 0, N*sizeof(int));
	cudaMemset(v_volcoeffsOfsSlices_Y.data, 0, N*sizeof(int));
#endif
	cudaMemset(sliceVoxel_count_.data, 0, N*sizeof(int));
	
	maxSliceDimension.x = size.x;
	maxSliceDimension.y = size.y;
}

void Reconstruction::FillSlices(double* sdata, std::vector<int> sizesX, std::vector<int> sizesY)
{
	num_slices_ = sizesX.size();
	checkCudaErrors(cudaMemcpy(v_slices.data, sdata, v_slices.size.x*v_slices.size.y*v_slices.size.z*sizeof(double), cudaMemcpyHostToDevice));
	CHECK_ERROR(FillSlices);

	if(d_slice_sicesX == NULL) checkCudaErrors(cudaMalloc((void **)&d_slice_sicesX, num_slices_*sizeof(int))); 
	if(d_slice_sicesY == NULL) checkCudaErrors(cudaMalloc((void **)&d_slice_sicesY, num_slices_*sizeof(int))); 

	checkCudaErrors(cudaMemcpy(d_slice_sicesX, &sizesX[0], sizesX.size()*sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_slice_sicesY, &sizesX[1], sizesY.size()*sizeof(int), cudaMemcpyHostToDevice));
}

void Reconstruction::getSlicesVol_debug(double* h_imdata)
{
	checkCudaErrors(cudaMemcpy(h_imdata, v_slices.data, v_slices.size.x*v_slices.size.y*v_slices.size.z*sizeof(double), cudaMemcpyDeviceToHost));
}


void Reconstruction::UpdateReconstructed(const uint3 vsize, double* data)
{
	cudaMemcpy(reconstructed_.data, data, vsize.x*vsize.y*vsize.z*sizeof(double), cudaMemcpyHostToDevice);
}



__global__ void GaussianConvolutionKernel(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict input, Volume<double> output,
	int kernel_size, double sigma, bool horizontal)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z);

	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];

	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0)
		return;

	unsigned int idx = pos.x + __umul24(pos.y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));

	//double s = input[idx];
	//////////////////

	double sum = 0;
	int half_kernel_elements = (kernel_size - 1) / 2;

	/*if(src[pos] != src[pos])
	printf("NAN %d %f \n ", horizontal, src[pos]);*/

	if (horizontal) {
		// convolve horizontally
		double g0 = 1.0 / (sqrt(2.0 * M_PI) * sigma); //norm
		double g1 = exp(-0.5 / (sigma * sigma));
		double g2 = g1 * g1;
		sum = g0 * input[idx];
		double sum_coeff = g0;
		//#pragma unroll 10
		for (unsigned int i = 1; i <= half_kernel_elements; i++){
			g0 *= g1;
			g1 *= g2;
			//border repeat
			unsigned int src_x = max(0, min(ssizeX - 1, (int) pos.x + (int) i));
			unsigned int idx2 = src_x + __umul24(pos.y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(src_x, pos2.y)]; 
			src_x = max(0, min(ssizeX - 1, (int) pos.x - (int) i));
			idx2 = src_x + __umul24(pos.y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(src_x, pos2.y)];
			sum_coeff += 2*g0;
		}
		output.set(pos, sum/sum_coeff);
	}
	else {
		// convolve vertically
		double g0 = 1.0 / (sqrt(2.0 * M_PI) * sigma);
		double g1 = exp(-0.5 / (sigma * sigma));
		double g2 = g1 * g1;
		sum = g0 * input[idx];
		double sum_coeff = g0;
		//#pragma unroll 10
		for (unsigned int j = 1; j <= half_kernel_elements; j++){
			g0 *= g1;
			g1 *= g2;
			//border repeat
			unsigned int src_y = max(0, min(ssizeY - 1, (int) pos.y + (int) j));
			unsigned int idx2 = pos.x + __umul24(src_y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(pos2.x, src_y)];
			src_y = max(0, min(ssizeY - 1, (int) pos.y - (int) j));
			idx2 = pos.x + __umul24(src_y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(pos2.x, src_y)];
			sum_coeff += 2*g0;
		}
		output.set(pos, sum/sum_coeff);
	}

}



__global__ void GaussianConvolutionKernel3D(double* __restrict input, Volume<double> output, int kernel_size, double sigma, short dir)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z);


	if(pos.x >= output.size.x || pos.y >= output.size.y || pos.z >= output.size.z || pos.x < 0 || pos.y < 0 || pos.z < 0)
		return;

	unsigned int idx = pos.x + __umul24(pos.y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));

	//double s = input[idx];
	//////////////////

	double sum = 0;
	int half_kernel_elements = (kernel_size - 1) / 2;

	/*if(src[pos] != src[pos])
	printf("NAN %d %f \n ", horizontal, src[pos]);*/

	if (dir == 0) {
		// convolve horizontally
		double g0 = 1.0 / (sqrt(2.0 * M_PI) * sigma); //norm
		double g1 = exp(-0.5 / (sigma * sigma));
		double g2 = g1 * g1;
		sum = g0 * input[idx];
		double sum_coeff = g0;
		//#pragma unroll 10
		for (unsigned int i = 1; i <= half_kernel_elements; i++){
			g0 *= g1;
			g1 *= g2;
			//border repeat
			unsigned int src_x = max(0, min(output.size.x - 1, (int) pos.x + (int) i));
			unsigned int idx2 = src_x + __umul24(pos.y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(src_x, pos2.y)]; 
			src_x = max(0, min(output.size.x - 1, (int) pos.x - (int) i));
			idx2 = src_x + __umul24(pos.y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(src_x, pos2.y)];
			sum_coeff += 2*g0;
		}
		output.set(pos, sum/sum_coeff);
	}
	else if(dir == 1)
	{
		// convolve vertically
		double g0 = 1.0 / (sqrt(2.0 * M_PI) * sigma);
		double g1 = exp(-0.5 / (sigma * sigma));
		double g2 = g1 * g1;
		sum = g0 * input[idx];
		double sum_coeff = g0;
		//#pragma unroll 10
		for (unsigned int j = 1; j <= half_kernel_elements; j++){
			g0 *= g1;
			g1 *= g2;
			//border repeat
			unsigned int src_y = max(0, min(output.size.y - 1, (int) pos.y + (int) j));
			unsigned int idx2 = pos.x + __umul24(src_y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(pos2.x, src_y)];
			src_y = max(0, min(output.size.y - 1, (int) pos.y - (int) j));
			idx2 = pos.x + __umul24(src_y,output.size.x) + __umul24(pos.z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(pos2.x, src_y)];
			sum_coeff += 2*g0;
		}
		output.set(pos, sum/sum_coeff);
	}
	else if(dir == 2)
	{
		// convolve depth
		double g0 = 1.0 / (sqrt(2.0 * M_PI) * sigma);
		double g1 = exp(-0.5 / (sigma * sigma));
		double g2 = g1 * g1;
		sum = g0 * input[idx];
		double sum_coeff = g0;
		//#pragma unroll 10
		for (unsigned int j = 1; j <= half_kernel_elements; j++){
			g0 *= g1;
			g1 *= g2;
			//border repeat
			unsigned int src_z = max(0, min(output.size.z - 1, (int) pos.z + (int) j));
			unsigned int idx2 = pos.x + __umul24(pos.y,output.size.x) + __umul24(src_z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(pos2.x, src_y)];
			src_z = max(0, min(output.size.z - 1, (int) pos.z - (int) j));
			idx2 = pos.x + __umul24(pos.y,output.size.x) + __umul24(src_z,__umul24(output.size.x,output.size.y));
			sum += g0 * input[idx2];//src[make_uint2(pos2.x, src_y)];
			sum_coeff += 2*g0;
		}
		output.set(pos, sum/sum_coeff);
	}

}

__global__ void calculateResidual3D_adv(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict slices, double* __restrict bias, double* __restrict weights, double* __restrict simweights, 
	double* __restrict simslices, Volume<double> wb_, Volume<double> wr_, double* scales)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];
	unsigned int idx = pos.x + pos.y*wb_.size.x + pos.z*wb_.size.x*wb_.size.y;
	double s = slices[idx];
	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0))
		return;

	double b = bias[idx];
	double w = weights[idx];
	double sw = simweights[idx];
	double ss = simslices[idx];
	double wb = wb_[pos];
	double wr = wr_[pos];
	double scale = scales[pos.z];
	double wbo = 0.0;
	double wro = 0.0;

	if( sw > 0.99)
	{
		double eb = exp(-b);
		double sliceVal = s*(eb * scale);
		wbo =  w * sliceVal;

		if( (ss > 1.0) && (sliceVal > 1.0)) {
			wro = log(sliceVal / ss) * wbo;
		}
	}

	wb_.set(pos, wbo);
	wr_.set(pos, wro);
}


__global__ void updateBiasField3D_adv(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict slices, Volume<double> bias, double* __restrict wb_, double* __restrict wr_) 
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z);

	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];

	unsigned int idx = pos.x + pos.y*bias.size.x + pos.z*bias.size.x*bias.size.y;

	double s = slices[idx];

	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0))
		return;

	//double b = bias[pos];
	double wb = wb_[idx];
	double wr = wr_[idx];

	if(wb > 0)
	{
		//b[idx2] += wr[idx2] / wb[idx2];
		bias.set(pos, bias[pos] + wr/wb);
	}
}

struct not_equal_than
{
	double val;
	not_equal_than(double t) { val = t; }
	__host__ __device__
		bool operator()(double x) { return x != val; }
};

void Reconstruction::debugWeights(double* weights)
{
	cudaDeviceSynchronize();
	checkCudaErrors(cudaMemcpy(weights, v_weights.data, v_weights.size.x*v_weights.size.y*v_weights.size.z*sizeof(double), cudaMemcpyDeviceToHost));	
}

void Reconstruction::debugBias(double* bias)
{
	cudaDeviceSynchronize();
	checkCudaErrors(cudaMemcpy(bias, v_bias.data, v_bias.size.x*v_bias.size.y*v_bias.size.z*sizeof(double), cudaMemcpyDeviceToHost));		
}

void Reconstruction::debugSimslices(double* simslices)
{
	cudaDeviceSynchronize();
	checkCudaErrors(cudaMemcpy(simslices, v_simulated_slices.data, v_simulated_slices.size.x*v_simulated_slices.size.y*v_simulated_slices.size.z*sizeof(double), cudaMemcpyDeviceToHost));		
}

void Reconstruction::debugSimweights(double* simweights)
{
	cudaDeviceSynchronize();
	checkCudaErrors(cudaMemcpy(simweights, v_simulated_weights.data, v_simulated_weights.size.x*v_simulated_weights.size.y*v_simulated_weights.size.z*sizeof(double), cudaMemcpyDeviceToHost));		
}


void Reconstruction::debugNormalizeBias(double* nbias)
{
	cudaDeviceSynchronize();
	checkCudaErrors(cudaMemcpy(nbias, bias_.data, bias_.size.x*bias_.size.y*bias_.size.z*sizeof(double), cudaMemcpyDeviceToHost));		
}

void Reconstruction::debugSmoothMask(double* smoothMask)
{
	cudaDeviceSynchronize();
	checkCudaErrors(cudaMemcpy(smoothMask, maskC_.data, maskC_.size.x*maskC_.size.y*maskC_.size.z*sizeof(double), cudaMemcpyDeviceToHost));		
}


void Reconstruction::CorrectBias(double sigma_bias, bool _global_bias_correction)
{
	dim3 blockSize3 = dim3(8,8,8); //check this invalid config
	dim3 gridSize3 = divup(dim3(v_slices.size.x, v_slices.size.y, num_slices_), blockSize3);
	unsigned int N = v_wb.size.x*v_wb.size.y*v_wb.size.z;
	cudaMemset(v_wb.data, 0, N*sizeof(double));
	cudaMemset(v_wresidual.data, 0, N*sizeof(double));


	calculateResidual3D_adv<<<gridSize3, blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_slices.data, v_bias.data, v_weights.data, 
		v_simulated_weights.data, v_simulated_slices.data, v_wb, v_wresidual, d_scales);

#if 0// this is 20% faster
	//TODO buf in gauss 3D kernel
	//int kernel_size, double sigma, bool horizontal
	//int klength = max(5, (unsigned int)ceil(sigma_bias * 5) * 2 + 1);//max(min((int)(sigma_bias*5),MAX_LENGTH_SK),7); // ... and make sure it is odd
	int klength = 2*round_(4*sigma_bias)+1;
	//printf("klength %d\n", klength);
	//klength -= 1-klength%2;
	//printf("klength %d\n ", klength);
	//only correct with buffer
	GaussianConvolutionKernel<<<gridSize3, blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_wb.data, v_buffer, klength, sigma_bias, true);
	GaussianConvolutionKernel<<<gridSize3, blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_buffer.data, v_wb, klength, sigma_bias, false);
	GaussianConvolutionKernel<<<gridSize3, blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_wresidual.data, v_buffer, klength, sigma_bias, true);
	GaussianConvolutionKernel<<<gridSize3, blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_buffer.data, v_wresidual, klength, sigma_bias, false);
	//printf("%d %f \n", klength, sigma_bias, 0.39*sigma_bias); //63*/

#else
	for(int i = 0; i < num_slices_; i++)
	{
		//Slice<double, Device>* s = slices_[i];
		//dim3 blockSize = dim3(8,8,1); //check this invalid config
	//	dim3 gridSize = divup(dim3(v_wb.size.x, v_wb.size.y, 1), blockSize);

		FilterGauss(v_wb.data + (i*v_wb.size.x*v_wb.size.y), v_wb.data + (i*v_wb.size.x*v_wb.size.y),
			v_wb.size.x, v_wb.size.y, v_wb.size.x, 1,  sigma_bias);
		FilterGauss(v_wresidual.data + (i*v_wresidual.size.x*v_wresidual.size.y), v_wresidual.data + (i*v_wresidual.size.x*v_wresidual.size.y),
			v_wresidual.size.x, v_wresidual.size.y, v_wresidual.size.x, 1,  sigma_bias);
	}
#endif

	updateBiasField3D_adv<<<gridSize3, blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_slices.data, v_bias, v_wb.data, v_wresidual.data);

	//untested
	if(_global_bias_correction)
	{
		/*	for(int i = 0; i < num_slices_; i++)
		{
		Slice<double, Device>* s = slices_[i];
		dim3 blockSize = dim3(8,8,1); //check this invalid config
		dim3 gridSize = divup(dim3(s->size.x, s->size.y, 1), blockSize);

		unsigned int N = slices_[i]->size.x*slices_[i]->size.y;
		not_equal_than pred(-1.0);
		thrust::device_ptr<double> dev_slice(slices_[i]->data());
		int num = thrust::count_if(dev_slice, dev_slice+N, pred);

		double sum2 = 0;
		thrust::device_ptr<double> dev_ptr(bias_[i]->data());
		sum2 = thrust::reduce(dev_ptr, dev_ptr + bias_[i]->size.x*bias_[i]->size.y, 0.0, thrust::plus<double>());

		/*Npp8u num2 = 0;
		thrust::device_ptr<Npp8u> dev_ptr2(mask.data());
		num2 = thrust::reduce(dev_ptr2, dev_ptr2 + mask.size.x*mask.size.y, 0.0, thrust::plus<Npp8u>());*/

		/*	double mean = sum2/(double)num; //weired mean caluculation

		CHECK_ERROR(updateBiasField);

		//printf("GPU %d %d \n", num2, num);

		normalizeBiasToZeroMean<<<gridSize, blockSize>>>(*(slices_[i]), *(bias_[i]), mean);
		}*/
	}
	//CHECK_ERROR(cudaDeviceSynchronize());

}


__global__ void AdaptiveRegularizationPrep(bool _adaptive, double _alpha, Volume<double> reconstructed, Volume<double> addon, Volume<double> confidence_map,
	double _min_intensity, double _max_intensity)
{
	const uint3 pos = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y,
		blockIdx.z * blockDim.z + threadIdx.z);

	if(pos.x >= reconstructed.size.x || pos.y >= reconstructed.size.y || pos.z >= reconstructed.size.z
		|| pos.x < 0 || pos.y < 0 || pos.z < 0 /*|| confidence_map[pos] <= 0*/)
		return;

	if(!_adaptive)
	{
		if(confidence_map[pos] != 0)
		{
			addon.set(pos, addon[pos] / confidence_map[pos]);
			confidence_map.set(pos, 1.0);
		}
	}

	reconstructed.set(pos, reconstructed[pos] + addon[pos] * _alpha);

	if (reconstructed[pos] < _min_intensity * 0.9)
		reconstructed.set(pos, _min_intensity * 0.9);
	if (reconstructed[pos] > _max_intensity * 1.1)
		reconstructed.set(pos, _max_intensity * 1.1);

}


__global__ void	weightedResidulaKernel(Volume<double> original, Volume<double> residual, Volume<double> weights, double _low_intensity_cutoff,
	double _max_intensity)
{

	const uint3 pos = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y,
		blockIdx.z * blockDim.z + threadIdx.z);

	if(pos.x >= original.size.x || pos.y >= original.size.y || pos.z >= original.size.z
		|| pos.x < 0 || pos.y < 0 || pos.z < 0)
		return;

	if ((weights[pos] == 1) && (original[pos] > _low_intensity_cutoff * _max_intensity) 
		&& (residual[pos] > _low_intensity_cutoff * _max_intensity)) 
	{
		if(original[pos] != 0)
		{
			double val = residual[pos]/original[pos];
			residual.set(pos, log(val));
		}
	}
	else
	{
		residual.set(pos, 0.0);
		weights.set(pos, 0.0);
	}

}

__global__ void	calcBiasFieldKernel(Volume<double> reconstructed, Volume<double> residual, Volume<double> weights, Volume<double> mask, 
	double _min_intensity, double _max_intensity)
{
	const uint3 pos = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y,
		blockIdx.z * blockDim.z + threadIdx.z);

	if(pos.x >= reconstructed.size.x || pos.y >= reconstructed.size.y || pos.z >= reconstructed.size.z
		|| pos.x < 0 || pos.y < 0 || pos.z < 0)
		return;

	if (mask[pos] == 1.0 && weights[pos] != 0) 
	{
		double val = residual[pos]/weights[pos];
		val = exp(val);
		residual.set(pos, val);
		double rval = reconstructed[pos]/val;
		if(rval  < _min_intensity * 0.9)
			rval = _min_intensity * 0.9;
		if (rval > _max_intensity * 1.1)
			rval = _max_intensity * 1.1;

		reconstructed.set(pos, rval);
	}
	else
	{
		residual.set(pos,0.0);
	}
}

void Reconstruction::syncCPU(double* reconstructed)
{
	cudaMemcpy(reconstructed, reconstructed_.data, reconstructed_.size.x*reconstructed_.size.y*reconstructed_.size.z*sizeof(double), cudaMemcpyDeviceToHost);
}

void Reconstruction::SyncConfidenceMapAddon(double* cmdata, double* addondata)
{
	cudaMemcpy(cmdata, confidence_map_.data, confidence_map_.size.x*confidence_map_.size.y*confidence_map_.size.z*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(addondata, addon_.data, addon_.size.x*addon_.size.y*addon_.size.z*sizeof(double), cudaMemcpyDeviceToHost);
}

//original without prep
__device__ double AdaptiveRegularization1(int i, uint3 pos, uint3 pos2, Volume<double> original, Volume<double> confidence_map, double delta)
{
	if(pos.x >= original.size.x || pos.y >= original.size.y || pos.z >= original.size.z
		|| pos.x < 0 || pos.y < 0 || pos.z < 0 || confidence_map[pos] <= 0 || confidence_map[pos2] <= 0 ||
		pos2.x >= original.size.x || pos2.y >= original.size.y || pos2.z >= original.size.z
		|| pos2.x < 0 || pos2.y < 0 || pos2.z < 0)
		return 0.0;

	//central differences would be better... improve with texture linear interpolation
	double diff = (original[pos2] - original[pos]) * sqrt(d_factor[i]) / delta;
	return d_factor[i] / sqrt(1.0 + diff * diff);
}

// 50% occupancy, 3.5 ms -- improve
//original with prep
__global__ void AdaptiveRegularizationKernel(Volume<double> reconstructed, Volume<double> original,
	Volume<double> confidence_map, double delta, double alpha, double lambda)
{
	uint3 pos = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y,
		blockIdx.z * blockDim.z + threadIdx.z);

	if(pos.x >= reconstructed.size.x || pos.y >= reconstructed.size.y || pos.z >= reconstructed.size.z
		|| pos.x < 0 || pos.y < 0 || pos.z < 0 || confidence_map[pos] <= 0.0)
		return;

	double val = 0;
	double valW = 0;
	double sum = 0;

	for (int i = 0; i < 13; i++) 
	{
		uint3 pos2 = make_uint3(pos.x + d_directions[i][0],pos.y + d_directions[i][1],pos.z + d_directions[i][2]);

		if ((pos2.x >= 0) && (pos2.x  < original.size.x) && (pos2.y >= 0) && (pos2.y < original.size.y) 
			&& (pos2.z >= 0) && (pos2.z < original.size.z)) 
		{
			double bi =  AdaptiveRegularization1(i, pos, pos2, original, confidence_map, delta);
			val += bi * reconstructed[pos2] * confidence_map[pos2]; //reconstructed == original2
			valW += bi * confidence_map[pos2];
			sum += bi;
		}

		uint3 pos3 = make_uint3(pos.x - d_directions[i][0],pos.y - d_directions[i][1],pos.z - d_directions[i][2]); //recycle pos register

		if ((pos3.x >= 0) && (pos3.x < original.size.x) && (pos3.y >= 0) && (pos3.y < original.size.y)
			&& (pos3.z >= 0) && (pos3.z < original.size.z) &&
			(pos2.x >= 0) && (pos2.x  < original.size.x) && (pos2.y >= 0) && (pos2.y < original.size.y) 
			&& (pos2.z >= 0) && (pos2.z < original.size.z)
			) 
		{
			double bi =  AdaptiveRegularization1(i, pos3, pos2, original, confidence_map, delta);
			val += bi * reconstructed[pos2] * confidence_map[pos2];
			valW += bi * confidence_map[pos2];
			sum += bi;
		}

	}

	val -= sum * reconstructed[pos] * confidence_map[pos];
	valW -= sum * confidence_map[pos];

	val = reconstructed[pos] * confidence_map[pos] + alpha * lambda / (delta * delta) * val;
	valW = confidence_map[pos] + alpha * lambda / (delta * delta) * valW;

	if (valW > 0) {
		reconstructed.set(pos, val / valW);
		//printf("%f ", val);
	}
	else
	{
		reconstructed.set(pos, 0.0);
	}

}

__global__ void SuperresolutionKernel3D_adv2(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	Volume<double> slices, Volume<double> bias, Volume<double> weights, Volume<double> simslices, 
	double* __restrict slice_weights, double* __restrict scales,
	Volume<int> ofssliceX, Volume<int> ofssliceY,
	Volume<double> addon, Volume<double> confidence_map, POINT3D* volcoeffsarray_)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];

	double s = slices[pos];
	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0))
		return;

	int ofsY = ofssliceY[pos];

	if((s == -1.0) || ofsY < 0)
		return;

	double b = bias[pos];
	double w = weights[pos];
	double ss = simslices[pos];
	double slice_weight = slice_weights[pos.z];
	double scale = scales[pos.z];

	if(pos.x == 25 && pos.y == 25 && pos.z == 25)
		printf("%f %f %f %f %f %f\n", s,b,w,ss,slice_weight,scale);

	double sliceVal = s * exp(-b) * scale;
	if ( ss > 0 )
		sliceVal -= ss;
	else
		sliceVal = 0;

	int n = ofssliceX[pos];
	int ofy = ofsY;
	//if(threadIdx.x == 0)
	//  printf("%f %f %f %d \n", b, w, ss, n);
	POINT3D	p;
	for (int k = 0; k < n; k++) {
		p = volcoeffsarray_[ofy + k];

		uint3 apos = make_uint3(p.x, p.y, p.z);
		if(apos.x < addon.size.x && apos.y < addon.size.y && apos.z < addon.size.z && 
			apos.x >= 0 && apos.y >= 0 && apos.z >= 0 )
		{
			//addon.set(apos, addon[apos] + p.value * sliceVal * w * slice_weight); //does this need an atomic??
			//confidence_map.set(apos, confidence_map[apos] + p.value * w * slice_weight); //does this need an atomic??
			atomicAdd(&(addon.data[apos.x + apos.y*addon.size.x + apos.z*addon.size.x*addon.size.y]), 
				p.value * sliceVal * w * slice_weight);
			atomicAdd(&(confidence_map.data[apos.x + apos.y*addon.size.x + apos.z*addon.size.x*addon.size.y]), 
				p.value * w * slice_weight);
		}
	}

}

void Reconstruction::Superresolution(int iter, std::vector<double> _slice_weight, bool _adaptive, double alpha, 
	double _min_intensity, double _max_intensity, double delta, double lambda, bool _global_bias_correction, double sigma_bias,
	double _low_intensity_cutoff)
{
	UpdateSliceWeights(_slice_weight);
	if (alpha * lambda / (delta * delta) > 0.068) 
	{
		printf("Warning: regularization might not have smoothing effect! Ensure that alpha*lambda/delta^2 is below 0.068.");
	}

	Volume<double> original;
	original.init(reconstructed_.size, reconstructed_.dim);
	checkCudaErrors(cudaMemcpy(original.data, reconstructed_.data, 
		reconstructed_.size.x*reconstructed_.size.y*reconstructed_.size.z*sizeof(double), cudaMemcpyDeviceToDevice));

	unsigned int N = addon_.size.x*addon_.size.y*addon_.size.z;
	cudaMemset(addon_.data, 0, N*sizeof(double));
	cudaMemset(confidence_map_.data, 0, N*sizeof(double));

	dim3 blockSize3D = dim3(8,8,8);
	dim3 gridSize3D = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, num_slices_), blockSize3D);

#if TEX_TEST
	int rest = num_slices_;
	for(int i = 0; i < ceil(num_slices_/(float)MAX_SLICES_PER_RUN); i++)
	{
		int thisrun = (rest >= MAX_SLICES_PER_RUN) ? MAX_SLICES_PER_RUN : rest;
		//printf("run %d %d %d \n", i, thisrun, num_slices_);

		dim3 blockSize3 = dim3(8,8,8);
		dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, MAX_SLICES_PER_RUN), blockSize3);

		SuperresolutionKernel3D_tex<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, 
			v_slices, v_bias, v_weights, v_simulated_slices, d_slice_weights, d_scales, v_PSF_sums_.data, mask_.data,
			addon_, confidence_map_.data, d_slicesI2W, d_slicesW2I, d_sliceDims, reconstructedVoxelSize, i*MAX_SLICES_PER_RUN);

		rest -= MAX_SLICES_PER_RUN;
	}

	CHECK_ERROR(SuperresolutionKernel3D_tex);
#else
	//for(int i = 0; i < h_slices_weights.size(); i++)
	//	printf("sw: %f \n", h_slices_weights[i]);

	SuperresolutionKernel3D_adv2<<<gridSize3D,blockSize3D>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, 
		v_slices, v_bias, v_weights, v_simulated_slices, d_slice_weights, d_scales, 
		v_volcoeffsOfsSlices_X, v_volcoeffsOfsSlices_Y,
		addon_, confidence_map_, volcoeffsarray_);
#endif //TEX_TEST

	CHECK_ERROR(SuperresolutionKernel);

	dim3 blockSize = dim3(8,8,8);
	dim3 gridSize = divup(dim3(addon_.size.x, addon_.size.y, addon_.size.z), blockSize);

	AdaptiveRegularizationPrep<<<gridSize,blockSize>>>(_adaptive, alpha, reconstructed_, 
		addon_, confidence_map_, _min_intensity, _max_intensity);

	AdaptiveRegularizationKernel<<<gridSize,blockSize>>>(reconstructed_, original, confidence_map_, delta, alpha, lambda);

	//untested!
	if(_global_bias_correction)
	{
		//TODO untested
		printf("_global_bias_correction \n");
		//need mask
		//BiasCorrectVolume

		// irtkRealImage residual = _reconstructed;
		//irtkRealImage weights = _mask;
		Volume<double> residual;
		residual.init(reconstructed_.size, reconstructed_.dim);
		checkCudaErrors(cudaMemcpy(residual.data, reconstructed_.data, 
			reconstructed_.size.x*reconstructed_.size.y*reconstructed_.size.z*sizeof(double), cudaMemcpyDeviceToDevice));

		Volume<double> weights;
		weights.init(mask_.size, mask_.dim);
		checkCudaErrors(cudaMemcpy(weights.data, mask_.data, mask_.size.x*mask_.size.y*mask_.size.z*sizeof(double), cudaMemcpyDeviceToDevice));

		//weighted residual kernel
		weightedResidulaKernel<<<gridSize,blockSize>>>(original, residual, weights, _low_intensity_cutoff, _max_intensity);
		//CHECK_ERROR(cudaDeviceSynchronize());
		//blur residual
		FilterGauss(residual.data,residual.data, residual.size.x, residual.size.y, residual.size.z, residual.size.x, residual.size.y, 1, sigma_bias );
		//blur weights
		FilterGauss(weights.data, weights.data, weights.size.x, weights.size.y, weights.size.z, weights.size.x, weights.size.y, 1,  sigma_bias );

		// calculate the bias field kenrel
		calcBiasFieldKernel<<<gridSize,blockSize>>>(reconstructed_, residual, weights, mask_, _min_intensity, _max_intensity);
		//CHECK_ERROR(cudaDeviceSynchronize());
	}

	original.release();
}


struct transformRS
{
	transformRS(){}

	__host__ __device__
		thrust::tuple<double, unsigned int> operator()(const thrust::tuple<double, char, double, double>& v)
	{
		const double s_ = thrust::get<0>(v);
		const char si_ = thrust::get<1>(v);
		const double ss_ = thrust::get<2>(v);
		const double sw_ = thrust::get<3>(v);

		if(s_ != -1 && si_ == 1 && sw_ > 0.99)
		{
			double sval = s_ - ss_;
			return make_tuple(sval*sval, 1);
		}
		else
		{
			return make_tuple(0.0, 0);
		}
	}

};

struct reduceRS
{
	reduceRS(){}

	__host__ __device__
		thrust::tuple<double, unsigned int> operator()(const thrust::tuple<double, unsigned int>& a,
		const thrust::tuple<double, unsigned int>& b)
	{
		return make_tuple(thrust::get<0>(a)+thrust::get<0>(b), thrust::get<1>(a)+thrust::get<1>(b));
	}
};

void Reconstruction::InitializeRobustStatistics(double& _sigma)
{
#if 1
	thrust::device_ptr<double> d_s(v_slices.data);//s->data());
	thrust::device_ptr<char> d_si(v_simulated_inside.data);//si->data());
	thrust::device_ptr<double> d_ss(v_simulated_slices.data);//ss->data());
	thrust::device_ptr<double> d_sw(v_simulated_weights.data);//sw->data());

	unsigned int N = v_slices.size.x*v_slices.size.y*v_slices.size.z;

	thrust::tuple<double, unsigned int> out =  transform_reduce(make_zip_iterator(make_tuple(d_s, d_si, d_ss, d_sw)), 
		make_zip_iterator(make_tuple(d_s+N, d_si+N, d_ss+N, d_sw+N)), transformRS(), 
		make_tuple(0.0,0), reduceRS());

	_sigma = get<0>(out)/(double)get<1>(out);
#else

	//TODO arrange slices in Volume and store offset pointer
	double sig = 0;
	unsigned int num = 0;
	for(int i = 0; i < num_slices_; i++)
	{
		Slice<double, Device>* s = slices_[i];
		/*		Slice<char, Device>* si = simulated_inside_[i];
		Slice<double, Device>* ss = simulated_slices_[i];
		Slice<double, Device>* sw = simulated_weights_[i];*/
		unsigned int N = slices_[i]->size.x*slices_[i]->size.y;

		thrust::device_ptr<double> d_s(v_slices.data + (i*v_slices.size.x*v_slices.size.y));//s->data());
		thrust::device_ptr<char> d_si(v_simulated_inside.data + (i*v_slices.size.x*v_slices.size.y));//si->data());
		thrust::device_ptr<double> d_ss(v_simulated_slices.data + (i*v_slices.size.x*v_slices.size.y));//ss->data());
		thrust::device_ptr<double> d_sw(v_simulated_weights.data + (i*v_slices.size.x*v_slices.size.y));//sw->data());

		thrust::tuple<double, unsigned int> out =  transform_reduce(make_zip_iterator(make_tuple(d_s, d_si, d_ss, d_sw)), 
			make_zip_iterator(make_tuple(d_s+N, d_si+N, d_ss+N, d_sw+N)), transformRS(), 
			make_tuple(0.0,0), reduceRS());

		sig += get<0>(out);
		num += get<1>(out);
		//if(get<1>(out) != 0)
		//_sigma +=  get<0>(out)/get<1>(out);
	}

	_sigma = sig/(double)num;

	//printf("GPU: %f \n", _sigma);
#endif
}


__global__ void gaussianReconstructionKernel3D(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict slices, double* __restrict bias2D, double* __restrict scales, int* __restrict ofssliceX, int* __restrict ofssliceY, 
	POINT3D* volcoeffsarray_, Volume<double> reconstructed, uint2 vSize)
{

	const uint3 pos = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y,
		blockIdx.z * blockDim.z + threadIdx.z);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;

	unsigned int idx = pos.x + pos.y*vSize.x + pos.z*vSize.x*vSize.y;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];

	double s = slices[idx];
	int ofy = ofssliceY[idx];

	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0) || ofy < 0)
		return;

	double b = bias2D[idx];
	int ofsX = ofssliceX[idx];

	double scale = scales[pos.z];

	int n = ofssliceX[idx];
	POINT3D p;

	double sliceVal = s * exp(-b) * scale;

	for (int k = 0; k < n; k++) {
		p = volcoeffsarray_[ofy + k];

		uint3 apos = make_uint3(p.x, p.y, p.z);
		if(apos.x < reconstructed.size.x && apos.y < reconstructed.size.y && apos.z < reconstructed.size.z && 
			apos.x >= 0 && apos.y >= 0 && apos.z >= 0 )
		{
			//bias.set(apos, bias[apos] + p.value*nbias); //atomic add if necessary
			atomicAdd(&(reconstructed.data[apos.x + apos.y*reconstructed.size.x + apos.z*reconstructed.size.x*reconstructed.size.y]), 
				p.value*sliceVal);
		}
	}


}

template< typename T >
class divS
{
public:
	T operator()(T a, T b)
	{
		return (b != 0) ? a / b : 0;
	}
};

template< typename T >
class divSame
{
public:
	T operator()(T a, T b)
	{
		return (b != 0) ? a / b : a;
	}
};


class divMin
{
public:
	double operator()(double a, double b)
	{
		return (b != 0) ? a / b : DBL_MIN;
	}
};

struct is_larger_zero
{
	__host__ __device__
		bool operator()(const int &x)
	{
		return (x > 0);
	}
};

void Reconstruction::GaussianReconstruction(std::vector<int>& voxel_num)
{
	//Test
	unsigned int N2 = v_weights.size.x*v_weights.size.y*v_weights.size.z;
	//cudaMemset(v_bias.data, 0, N*sizeof(double));
	cudaMemset(v_weights.data, 0, N2*sizeof(double));
	cudaMemset(v_simulated_weights.data, 0, N2*sizeof(double));
	cudaMemset(v_simulated_slices.data, 0, N2*sizeof(double));
	cudaMemset(v_wresidual.data, 0, N2*sizeof(double));
	cudaMemset(v_wb.data, 0, N2*sizeof(double));
	cudaMemset(v_buffer.data, 0, N2*sizeof(double));
	cudaMemset(v_simulated_inside.data, 0, N2*sizeof(char));

	//clear _reconstructed image
	//initMem<double>(reconstructed_.data, reconstructed_.size.x*reconstructed_.size.y*reconstructed_.size.z, 0.0);
	cudaMemset(reconstructed_.data, 0, reconstructed_.size.x*reconstructed_.size.y*reconstructed_.size.z*sizeof(double));

	dim3 blockSize3 = dim3(8,8,8);
	dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, num_slices_), blockSize3);
	/*	gaussianReconstructionKernel3D<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, 
	d_slices, d_bias, d_scales, 
	d_volcoeffsOfsSlices_XPtr, d_volcoeffsOfsSlices_YPtr, volcoeffsarray_, reconstructed_);*/
	unsigned int N = reconstructed_.size.x*reconstructed_.size.y*reconstructed_.size.z;

#if TEX_TEST

	cudaMemset(reconstructed_volWeigths.data, 0, N*sizeof(double));

	unsigned int N1 = sliceVoxel_count_.size.x*sliceVoxel_count_.size.y*sliceVoxel_count_.size.z;
	cudaMemset(sliceVoxel_count_.data, 0, N1*sizeof(int));


	int rest = num_slices_;
	for(int i = 0; i < ceil(num_slices_/(float)MAX_SLICES_PER_RUN); i++)
	{
		int thisrun = (rest >= MAX_SLICES_PER_RUN) ? MAX_SLICES_PER_RUN : rest;
		//printf("run %d %d %d \n", i, thisrun, num_slices_);

		dim3 blockSize3 = dim3(8,8,8);
		dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, MAX_SLICES_PER_RUN), blockSize3);
		gaussianReconstructionKernel3D_tex<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, 
			v_slices.data, v_bias.data, d_scales, 
			reconstructed_, reconstructed_volWeigths, v_PSF_sums_, sliceVoxel_count_, mask_.data,
			make_uint2(v_slices.size.x, v_slices.size.y), 
			d_slicesI2W, d_slicesW2I, d_sliceDims, reconstructedVoxelSize, i*MAX_SLICES_PER_RUN);

		rest -= MAX_SLICES_PER_RUN;
	}

	thrust::device_ptr<double> ptr_recons(reconstructed_.data);
	thrust::device_ptr<double> ptr_count(reconstructed_volWeigths.data);
	thrust::transform(ptr_recons, ptr_recons+N, ptr_count, ptr_recons, divSame<double>());

	for(int i = 0; i < num_slices_; i++)
	{
		unsigned int N2 = sliceVoxel_count_.size.x*sliceVoxel_count_.size.y;
		thrust::device_ptr<int> ptr_count(sliceVoxel_count_.data + (i*sliceVoxel_count_.size.x*sliceVoxel_count_.size.y)  );
		//only ++ if n > 0 == set pixel
		int vnum = thrust::count_if(ptr_count, ptr_count+N2, is_larger_zero()); 
		voxel_num.push_back(vnum);
		//printf("%d \n", vnum);
	}
#else

	gaussianReconstructionKernel3D<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, 
		v_slices.data, v_bias.data, d_scales, 
		v_volcoeffsOfsSlices_X.data, v_volcoeffsOfsSlices_Y.data, volcoeffsarray_, 
		reconstructed_, make_uint2(v_slices.size.x, v_slices.size.y));
	thrust::device_ptr<double> ptr_recon(reconstructed_.data);
	thrust::device_ptr<double> ptr_volw(volume_weights_.data);
	thrust::transform(ptr_recon, ptr_recon+N, ptr_volw, ptr_recon, divS<double>());
	CHECK_ERROR(gaussianReconstructionKernel3D);

	//TODO make volume for speedup!!
	//TODO this has to be different with inplace...
	for(int i = 0; i < num_slices_; i++)
	{
		unsigned int N1 = v_volcoeffsOfsSlices_X.size.x*v_volcoeffsOfsSlices_X.size.y;
		thrust::device_ptr<int> ptr_volcoOfsX(v_volcoeffsOfsSlices_X.data + (i*v_volcoeffsOfsSlices_X.size.x*v_volcoeffsOfsSlices_X.size.y)  );
		//only ++ if n > 0 == set pixel
		int vnum = thrust::count_if(ptr_volcoOfsX, ptr_volcoOfsX+N1, is_larger_zero()); 
		voxel_num.push_back(vnum);
	}
#endif //TEX_TEST
	CHECK_ERROR(gaussianReconstructionKernel3Dcount_if);
}


__global__ void normalizeBiasKernel3D_adv(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	Volume<double> slices, Volume<double> bias2D, double* __restrict scales,
	Volume<int> ofssliceX, Volume<int> ofssliceY,
	POINT3D* volcoeffsarray_, Volume<double> bias)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];

	double s = slices[pos];
	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0))
		return;

	int ofsY = ofssliceY[pos];

	if((s == -1.0) || ofsY < 0)
		return;

	double nbias = bias2D[pos];
	double scale = scales[pos.z];

	if(s > -1.0 && scale > 0)
	{
		nbias -= log(scale);
	}

	int n = ofssliceX[pos];
	int ofy = ofsY;

	POINT3D	p;
	for (int k = 0; k < n; k++) {
		p = volcoeffsarray_[ofy + k];

		uint3 apos = make_uint3(p.x, p.y, p.z);
		if(apos.x < bias.size.x && apos.y < bias.size.y && apos.z < bias.size.z && 
			apos.x >= 0 && apos.y >= 0 && apos.z >= 0 )
		{
			atomicAdd(&(bias.data[apos.x + apos.y*bias.size.x + apos.z*bias.size.x*bias.size.y]), p.value*nbias);
		}
	}

}


template< typename T >
class maskS
{
public:
	T operator()(T a, T b)
	{
		return (b != 0) ? a : 0;
	}
};

template< typename T >
class divexp
{
public:
	T operator()(T a, T b)
	{
		return  (a!=-1) ? (a / exp(-b)) : a;
	}
};




void Reconstruction::NormaliseBias(int iter, double sigma_bias)
{

	cudaMemset(bias_.data, 0.0, bias_.size.x*bias_.size.y*bias_.size.z*sizeof(double));

	dim3 blockSize3 = dim3(8,8,8);
	dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, num_slices_), blockSize3);


#if TEX_TEST
	int rest = num_slices_;
	for(int i = 0; i < ceil(num_slices_/(float)MAX_SLICES_PER_RUN); i++)
	{
		int thisrun = (rest >= MAX_SLICES_PER_RUN) ? MAX_SLICES_PER_RUN : rest;
		//printf("run %d %d %d \n", i, thisrun, num_slices_);

		dim3 blockSize3 = dim3(8,8,8);
		dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, MAX_SLICES_PER_RUN), blockSize3);

		normalizeBiasKernel3D_tex<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_slices, v_bias, d_scales,  v_PSF_sums_.data, mask_.data,
		d_slicesI2W, d_slicesW2I, d_sliceDims, reconstructedVoxelSize, bias_, reconstructed_volWeigths, i*MAX_SLICES_PER_RUN);

		rest -= MAX_SLICES_PER_RUN;
	}

	unsigned int N = bias_.size.x*bias_.size.y*bias_.size.z;
	thrust::device_ptr<double> ptr_bias(bias_.data);

	//zip iterator 
	//bias /= _volume_weights; -- this could be obsolete now
	thrust::device_ptr<double> ptr_volWeights(reconstructed_volWeigths.data);
	thrust::transform(ptr_bias, ptr_bias+N, ptr_volWeights, ptr_bias, divS<double>()); //TODO check if divS or divSame

	//MaskImage bias
	thrust::device_ptr<double> ptr_mask(mask_.data); //mask does not need to be a double
	thrust::transform(ptr_bias, ptr_bias+N, ptr_mask, ptr_bias, maskS<double>()); //test this

	thrust::device_ptr<double> ptr_maskC(maskC_.data); //mask does not need to be a double

	dim3 blockSize3k = dim3(8,8,8);
	dim3 gridSize3k = divup(dim3(bias_.size.x, bias_.size.y, bias_.size.z), blockSize3k);

	Volume<double> mbuf;
	mbuf.init(maskC_.size, maskC_.dim);

	//gives us still to high values
	//int klength = max(5, (unsigned int)ceil(sigma_bias * 5) * 2 + 1);
	/*int klength = 2*round_(4*sigma_bias)+1;
	GaussianConvolutionKernel3D<<<gridSize3k,blockSize3k>>>(bias_.data, mbuf, klength, sigma_bias, 0);
	GaussianConvolutionKernel3D<<<gridSize3k,blockSize3k>>>(mbuf.data, bias_, klength, sigma_bias, 1);
	GaussianConvolutionKernel3D<<<gridSize3k,blockSize3k>>>(bias_.data, mbuf, klength, sigma_bias, 2);
	checkCudaErrors(cudaMemcpy(bias_.data, mbuf.data, bias_.size.x*bias_.size.y*bias_.size.z*sizeof(double), cudaMemcpyDeviceToDevice));*/
	//normalize??
	/*double norm = sqrt((double)thrust::inner_product(ptr_bias, ptr_bias+N, ptr_bias, 0));
	using namespace thrust::placeholders;
    thrust::transform(ptr_bias, ptr_bias+N, ptr_bias, _1 /= norm);*/

	//FilterGauss mask and bias
	FilterGauss(bias_.data, bias_.data, bias_.size.x, bias_.size.y, bias_.size.z, bias_.size.x, bias_.size.y, 1,  sigma_bias );
	//FilterGauss(maskC.data, maskC.data, maskC.size.x, maskC.size.y, maskC.size.z, maskC.size.x, maskC.size.y, 1,  sigma_bias );

	//bias/=m;
	thrust::transform(ptr_bias, ptr_bias+N, ptr_maskC, ptr_bias, divS<double>());

	thrust::device_ptr<double> ptr_reconstructed(reconstructed_.data);
	// *pi /=exp(-(*pb));
	thrust::transform(ptr_reconstructed, ptr_reconstructed+N, ptr_bias, ptr_reconstructed, divexp<double>());

	//thrust::transform(ptr_reconstructed, ptr_reconstructed+N, ptr_mask, ptr_reconstructed, maskS<double>()); //test this

#else
	normalizeBiasKernel3D_adv<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_slices, v_bias, d_scales,
		v_volcoeffsOfsSlices_X, v_volcoeffsOfsSlices_Y, volcoeffsarray_, bias_);


	unsigned int N = bias_.size.x*bias_.size.y*bias_.size.z;
	thrust::device_ptr<double> ptr_bias(bias_.data);

	//zip iterator 
	//bias /= _volume_weights;
	thrust::device_ptr<double> ptr_volWeights(volume_weights_.data);
	thrust::transform(ptr_bias, ptr_bias+N, ptr_volWeights, ptr_bias, divS<double>());

	//MaskImage bias
	thrust::device_ptr<double> ptr_mask(mask_.data); //mask does not need to be a double
	thrust::transform(ptr_bias, ptr_bias+N, ptr_mask, ptr_bias, maskS<double>()); //test this

	thrust::device_ptr<double> ptr_maskC(maskC_.data);
	//FilterGauss mask and bias
	FilterGauss(bias_.data, bias_.data, bias_.size.x, bias_.size.y, bias_.size.z, bias_.size.x, bias_.size.y, 1,  sigma_bias );
	//FilterGauss(maskC.data, maskC.data, maskC.size.x, maskC.size.y, maskC.size.z, maskC.size.x, maskC.size.y, 1,  sigma_bias );

	//bias/=m;
	thrust::transform(ptr_bias, ptr_bias+N, ptr_maskC, ptr_bias, divS<double>());

	thrust::device_ptr<double> ptr_reconstructed(reconstructed_.data);
	// *pi /=exp(-(*pb));
	thrust::transform(ptr_reconstructed, ptr_reconstructed+N, ptr_bias, ptr_reconstructed, divexp<double>());

#endif

}

__global__ void simulateSlicesKernel3D_adv(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	Volume<double> slices,  Volume<double> simslices, Volume<double> simweights, Volume<char> siminside,
	Volume<double> reconstructed, Volume<double> mask,
	Volume<int> ofssliceX, Volume<int> ofssliceY, POINT3D* volcoeffsarray_)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];

	double s = slices[pos];
	int n = ofssliceX[pos];
	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1.0) || n < 0)
		return;

	double simulated_sliceV = 0;
	bool slice_inside = false;
	double weight = 0;

	int ofy = ofssliceY[pos];
	POINT3D p;
	for (int k = 0; k < n; k++) {
		p = volcoeffsarray_[ofy + k];

		uint3 apos = make_uint3(p.x, p.y, p.z);
		if(apos.x < reconstructed.size.x && apos.y < reconstructed.size.y && apos.z < reconstructed.size.z && 
			apos.x >= 0 && apos.y >= 0 && apos.z >= 0 )
		{
			simulated_sliceV += p.value * reconstructed[apos];
			weight += p.value;

			if(mask[apos] != 0)
			{
				slice_inside = true; 
			}
		}
	}

	siminside.set(pos, (char)slice_inside);

	if( weight > 0 ) 
	{
		simslices.set(pos,simulated_sliceV/weight);
		simweights.set(pos,weight);
	}
	else
	{
		simslices.set(pos,simulated_sliceV);
	}

}

void Reconstruction::SimulateSlices(std::vector<bool>& slice_inside)
{
	dim3 blockSize3 = dim3(8,8,8);
	dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, num_slices_), blockSize3);

#if  TEX_TEST
	int rest = num_slices_;
	for(int i = 0; i < ceil(num_slices_/(float)MAX_SLICES_PER_RUN); i++)
	{
		int thisrun = (rest >= MAX_SLICES_PER_RUN) ? MAX_SLICES_PER_RUN : rest;
		//printf("run %d %d %d \n", i, thisrun, num_slices_);

		dim3 blockSize3 = dim3(8,8,8);
		dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, MAX_SLICES_PER_RUN), blockSize3);

		simulateSlicesKernel3D_tex<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, 
			v_slices.data, v_simulated_slices, v_simulated_weights, v_simulated_inside, 
			reconstructed_, v_PSF_sums_.data, mask_.data, make_uint2(v_slices.size.x, v_slices.size.y), 
			d_slicesI2W, d_slicesW2I, d_sliceDims, reconstructedVoxelSize, i*MAX_SLICES_PER_RUN);

		rest -= MAX_SLICES_PER_RUN;
	}

#else
	simulateSlicesKernel3D_adv<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, 
		v_slices, v_simulated_slices, v_simulated_weights, v_simulated_inside, 
		reconstructed_, mask_, v_volcoeffsOfsSlices_X, v_volcoeffsOfsSlices_Y, volcoeffsarray_);
#endif //TEX_TEST

	//checkCudaErrors(cudaDeviceSynchronize());
	CHECK_ERROR(simulateSlicesKernel3D());

	for(int i = 0; i < num_slices_; i++)
	{
		unsigned int N = v_simulated_inside.size.x*v_simulated_inside.size.y;
		thrust::device_ptr<char> d_si(v_simulated_inside.data + (i*v_simulated_inside.size.x*v_simulated_inside.size.y));
		int h_sliceInside = thrust::reduce(d_si, d_si + N, 0, thrust::plus<char>());
		if(h_sliceInside > 0)
			slice_inside[i] = true;
		else
			slice_inside[i] = false;
	}

	unsigned int N = bias_.size.x*bias_.size.y*bias_.size.z;
	thrust::device_ptr<double> ptr_mask(mask_.data);
	thrust::device_ptr<double> ptr_reconstructed(reconstructed_.data);
	thrust::transform(ptr_reconstructed, ptr_reconstructed+N, ptr_mask, ptr_reconstructed, maskS<double>());
}



struct transformSlicePotential
{
	__host__ __device__
		tuple<double,double> operator()(const tuple<double,double>& a)
		//tuple<d_w,d_sw> out tuple<potenial,num>
	{

		if(thrust::get<1>(a) > 0.99 && thrust::get<0>(a) != 0.0)
		{
			return thrust::make_tuple(((1.0 - thrust::get<0>(a)) * (1.0 - thrust::get<0>(a))), 1.0);
		}
		else
		{
			return thrust::make_tuple(0.0,0.0);
		}

	}
};


__global__ void EStepKernel3D_tex(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict slices, double* __restrict bias, double* weights, 
	double* __restrict simslices, double* __restrict simweights, double* __restrict scales, double _m, double _sigma, double _mix, uint2 vSize)
{

	const uint3 pos = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y,
		blockIdx.z * blockDim.z + threadIdx.z);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];
	unsigned int idx = pos.x + pos.y*vSize.x + pos.z*vSize.x*vSize.y;

	double s = slices[idx];
	double sw = simweights[idx];

	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1) || sw <= 0)
		return;

	double b = bias[idx];
	double w = weights[idx];
	double ss = simslices[idx];
	double scale = scales[pos.z];

	double sliceVal = s * exp(-b) * scale;

	sliceVal -= ss;

	//Gaussian distribution for inliers (likelihood)
	double g = G(sliceVal, _sigma);
	//Uniform distribution for outliers (likelihood)
	double m = M(_m);

	double weight = g * _mix / (g *_mix + m * (1.0 - _mix));
	weights[idx] = weight;
	//w[idx2] = weight;
}


__global__ void EStepKernel3D(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict slices, double* __restrict bias, double* weights, 
	double* __restrict simslices, double* __restrict simweights, double* __restrict scales,
	int* __restrict ofssliceX, double _m, double _sigma, double _mix, uint2 vSize)
{

	const uint3 pos = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y,
		blockIdx.z * blockDim.z + threadIdx.z);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];
	unsigned int idx = pos.x + pos.y*vSize.x + pos.z*vSize.x*vSize.y;

	double s = slices[idx];
	int ofsX = ofssliceX[idx];
	double sw = simweights[idx];

	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0 || (s == -1) || sw <= 0 || ofsX <= 0)
		return;

	double b = bias[idx];
	double w = weights[idx];
	double ss = simslices[idx];
	double scale = scales[pos.z];

	double sliceVal = s * exp(-b) * scale;

	sliceVal -= ss;

	//Gaussian distribution for inliers (likelihood)
	double g = G(sliceVal, _sigma);
	//Uniform distribution for outliers (likelihood)
	double m = M(_m);

	double weight = g * _mix / (g *_mix + m * (1.0 - _mix));
	weights[idx] = weight;
}


struct reduceSlicePotential
{
	__host__ __device__
		tuple<double,double> operator()(const tuple<double,double>& a, const tuple<double,double>& b)
		//tuple<d_w,d_sw> out tuple<potenial,num>
	{
		return thrust::make_tuple(thrust::get<0>(a) + thrust::get<0>(b), thrust::get<1>(a) + thrust::get<1>(b));
	}
};


void Reconstruction::EStep(double _m, double _sigma, double _mix, std::vector<double>& slice_potential)
{
	dim3 blockSize3 = dim3(8,8,8);
	dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, num_slices_), blockSize3);

	unsigned int N = v_weights.size.x*v_weights.size.y*v_weights.size.z; 	
	cudaMemset(v_weights.data, 0, N*sizeof(double));

#if TEX_TEST
	EStepKernel3D_tex<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_slices.data, v_bias.data, v_weights.data,
		v_simulated_slices.data, v_simulated_weights.data, d_scales, _m, _sigma, _mix, 
		make_uint2(v_slices.size.x, v_slices.size.y));
#else
	EStepKernel3D<<<gridSize3,blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, v_slices.data, v_bias.data, v_weights.data,
		v_simulated_slices.data, v_simulated_weights.data, d_scales, v_volcoeffsOfsSlices_X.data, _m, _sigma, _mix, 
		make_uint2(v_slices.size.x, v_slices.size.y));
#endif

	for(int i = 0; i < num_slices_; i++)
	{

		unsigned int N = v_weights.size.x*v_weights.size.y;

		thrust::device_ptr<double> d_w(v_weights.data + (i*v_weights.size.x*v_weights.size.y));//w->data());
		thrust::device_ptr<double> d_sw(v_simulated_weights.data + (i*v_weights.size.x*v_weights.size.y));//sw->data());
		tuple<double,double> out = transform_reduce(make_zip_iterator(make_tuple(d_w, d_sw)), make_zip_iterator(make_tuple(d_w+N, d_sw+N)), 
			transformSlicePotential(), make_tuple(0.0, 0.0),  reduceSlicePotential());

		if (thrust::get<1>(out) > 0)
		{
			slice_potential[i] = sqrt(thrust::get<0>(out) / thrust::get<1>(out));
			//printf("%f \n", slice_potential[i]);
		}
		else
		{
			slice_potential[i] = -1; // slice has no unpadded voxels
		}
	//	printf("slice_potential %f \n", slice_potential[i]);	
	}
}


struct transformMStep
{
	double scale;
	transformMStep(double scale_){scale = scale_;}

	__host__ __device__
		thrust::tuple<double, double, double, double, double> operator()(const thrust::tuple<double, double, double, double, double>& v)
		//thrust::tuple<sigma_, mix_, count, e, e> //this order is very important for thrust -- why?
	{
		double s_ = thrust::get<0>(v);
		const double b_ = thrust::get<1>(v);
		const double w_ = thrust::get<2>(v);
		const double ss_ = thrust::get<3>(v);
		const double sw_ = thrust::get<4>(v);

		double sigma_ = 0.0;
		double mix_ = 0.0;
		double count = 0.0;
		double e = 0.0;

		if(s_ != -1.0)
		{
			s_ *= exp(-b_) * scale;

			if ( sw_ > 0.99 ) {
				s_ -= ss_;
				e = s_;			
				sigma_ = e * e * w_;
				mix_ = w_;
				count = 1.0;
			}
		}
		else
		{
			s_ = -1.0;
		}
		return thrust::make_tuple(sigma_, mix_, count, e, e);
	}
};

struct transformMStep3D
{
	__host__ __device__
		thrust::tuple<double, double, double, double, double> operator()(const thrust::tuple<double, double, double, double, double, double>& v)
		//thrust::tuple<sigma_, mix_, count, e, e> //this order is very important for the thrust optization
	{
		double s_ = thrust::get<0>(v);
		const double b_ = thrust::get<1>(v);
		const double w_ = thrust::get<2>(v);
		const double ss_ = thrust::get<3>(v);
		const double sw_ = thrust::get<4>(v);
		const double scale = thrust::get<5>(v);

		double sigma_ = 0.0;
		double mix_ = 0.0;
		double count = 0.0;
		double e = 0.0;

		if(s_ != -1.0)
		{
			s_ *= exp(-b_) * scale; // *= exp(-b(i, j, 0)) * scale;
			if ( sw_ > 0.99 ) {
				s_ -= ss_;
				e = s_;			
				sigma_ = e * e * w_;
				mix_ = w_;
				count = 1.0;
			}
		}
		return thrust::make_tuple(sigma_, mix_, count, e, e);
	}
};



struct reduceMStep 
{
	__host__ __device__
		thrust::tuple<double, double, double, double, double> operator()(const thrust::tuple<double, double, double, double, double>& a,
		const thrust::tuple<double, double, double, double, double>& b)
	{
		return thrust::make_tuple(thrust::get<0>(a)+thrust::get<0>(b), thrust::get<1>(a)+thrust::get<1>(b),
			thrust::get<2>(a)+thrust::get<2>(b), min(thrust::get<3>(a),thrust::get<3>(b)),
			max(thrust::get<4>(a),thrust::get<4>(b)));
	}
};

void Reconstruction::MStep(int iter, double _step, double& _sigma, double& _mix, double& _m)
{

	thrust::device_ptr<double> d_s(v_slices.data);
	thrust::device_ptr<double> d_b(v_bias.data);
	thrust::device_ptr<double> d_w(v_weights.data);
	thrust::device_ptr<double> d_ss(v_simulated_slices.data);
	thrust::device_ptr<double> d_sw(v_simulated_weights.data);
	thrust::device_ptr<double> d_buf(v_buffer.data);

	//we use buffer for scales -- TODO improved with const_iterator
	for(int i = 0; i < h_scales.size(); i++)
	{
		unsigned int N1 = v_buffer.size.x*v_buffer.size.y;
		initMem<double>(v_buffer.data,  N1, h_scales[i]);
	}

	unsigned int N3 = v_buffer.size.x*v_buffer.size.y*v_buffer.size.z;
	tuple<double,double,double,double,double> out =  transform_reduce(make_zip_iterator(make_tuple(d_s, d_b, d_w, d_ss, d_sw, d_buf)), 
		make_zip_iterator(make_tuple(d_s+N3, d_b+N3, d_w+N3, d_ss+N3, d_sw+N3, d_buf+N3)), transformMStep3D(), 
		make_tuple<double,double,double,double,double>(0.0, 0.0, 0.0, 0.0, 0.0), reduceMStep());

	double sigma = get<0>(out);
	double mix = get<1>(out);
	double num = get<2>(out);
	double min_ = get<3>(out);
	double max_ = get<4>(out);
	//	printf("GPU min %f, max %f, sigma %f, mix %f, num %f\n", _min, _max, _sigma, _mix, _num);

	if (mix > 0) {
		_sigma = sigma / mix;
	}
	else {
		printf("Something went wrong: sigma= %f mix= %f\n", sigma, mix);
		//exit(1);
	}
	if (_sigma < _step * _step / 6.28) 
		_sigma = _step * _step / 6.28;
	if (iter > 1)
		_mix = mix / num;

	//Calculate m
	_m = 1 / (max_ - min_);
}


struct transformScale
{
	transformScale(){}

	__host__ __device__
		thrust::tuple<double, double> operator()(const thrust::tuple<double, double, double, double, double>& v)
	{
		double s_ = thrust::get<0>(v);
		const double b_ = thrust::get<1>(v);
		const double w_ = thrust::get<2>(v);
		const double ss_ = thrust::get<3>(v);
		const double sw_ = thrust::get<4>(v);

		if( (s_ == -1) || sw_ <= 0.99)
		{
			return make_tuple(0.0,0.0);
		}
		else
		{
			double eb = exp(-(b_));
			double scalenum = w_ * s_ * eb * ss_;
			double scaleden = w_ * s_ * eb * s_ * eb;
			return make_tuple(scalenum, scaleden);
		}
	}
};


struct reduceScale
{
	reduceScale(){}

	__host__ __device__
		thrust::tuple<double, double> operator()(const thrust::tuple<double, double>& a, const thrust::tuple<double, double>& b)
	{
		return make_tuple(thrust::get<0>(a)+thrust::get<0>(b), thrust::get<1>(a)+thrust::get<1>(b));
	}
};

void Reconstruction::CalculateScaleVector()
{
	h_scales.clear();

	//this is a transform_reduce with 5er tuple
	for(int i = 0; i < num_slices_; i++)
	{
		unsigned int N = v_slices.size.x*v_slices.size.y;
#if 1
		thrust::device_ptr<double> d_s(v_slices.data + (i*v_slices.size.x*v_slices.size.y));
		thrust::device_ptr<double> d_b(v_bias.data + (i*v_slices.size.x*v_slices.size.y));
		thrust::device_ptr<double> d_w(v_weights.data + (i*v_slices.size.x*v_slices.size.y));
		thrust::device_ptr<double> d_ss(v_simulated_slices.data + (i*v_slices.size.x*v_slices.size.y));
		thrust::device_ptr<double> d_sw(v_simulated_weights.data + (i*v_slices.size.x*v_slices.size.y));

		thrust::tuple<double, double> out =  transform_reduce(make_zip_iterator(make_tuple(d_s, d_b, d_w, d_ss, d_sw)), 
			make_zip_iterator(make_tuple(d_s+N, d_b+N, d_w+N, d_ss+N, d_sw+N)), transformScale(), 
			make_tuple(0.0,0.0), reduceScale());

		if (thrust::get<1>(out) > 0)
		{
			h_scales.push_back(thrust::get<0>(out)/thrust::get<1>(out));
		}
		else
		{
			h_scales.push_back(1.0);
		}

		//printf("scale %f \n", h_scales[i]);		
#else
		Slice<double, Device> scalenum;
		scalenum.alloc(slices_[i]->size);
		Slice<double, Device> scaleden;
		scaleden.alloc(slices_[i]->size);
		//kernel bias, weights, slice, simluated slices
		/*Slice<double, Device>* s = slices_[i];
		Slice<double, Device>* w = weights_[i];
		Slice<double, Device>* b = bias_[i];
		Slice<double, Device>* ss = simulated_slices_[i];
		Slice<double, Device>* sw = simulated_weights_[i];*/
		dim3 blockSize = dim3(BLOCK_SIZE,BLOCK_SIZE,1);
		dim3 gridSize = divup(dim3(s->size.x, s->size.y, 1), blockSize);

		CalculateScaleSlices<<<gridSize,blockSize>>>(*s, *w, *b, *ss, *sw, scalenum, scaleden);
		//CHECK_ERROR(cudaDeviceSynchronize());

		//thrust reduction of slices
		thrust::device_ptr<double> scalenum_dev_ptr(scalenum.data());
		thrust::device_ptr<double> scaleden_dev_ptr(scaleden.data());

		double h_scalenum = thrust::reduce(scalenum_dev_ptr, scalenum_dev_ptr + N, 0.0, thrust::plus<double>());
		double h_scaleden = thrust::reduce(scaleden_dev_ptr, scaleden_dev_ptr + N, 0.0, thrust::plus<double>());

		//scale[i] = sum scalenum / sum scaleden;
		if (h_scaleden > 0)
			h_scales.push_back(h_scalenum/h_scaleden);
		else
			h_scales.push_back(1.0);

		//double test = (thrust::get<1>(out) > 0) ? thrust::get<0>(out)/thrust::get<1>(out) : 1.0;
		//printf("%f %f \n", test , h_scales[i]);
#endif
	}

	checkCudaErrors(cudaMemcpy(d_scales, &h_scales[0], num_slices_*sizeof(double), cudaMemcpyHostToDevice));
}

__global__ void InitializeEMValuesKernel(int numSlices, int* __restrict d_sliceSicesX, int* __restrict d_sliceSicesY, 
	double* __restrict slices, double* bias, double* weights, uint2 vSize)
{
	const uint3 pos = make_uint3(__umul24(blockIdx.x, blockDim.x) + threadIdx.x,
		__umul24(blockIdx.y, blockDim.y) + threadIdx.y,
		__umul24(blockIdx.z, blockDim.z) + threadIdx.z);

	//z is slice index
	if(pos.z >= numSlices || pos.z < 0)
		return;

	int ssizeX = d_sliceSicesX[pos.z];
	int ssizeY = d_sliceSicesY[pos.z];
	unsigned int idx = pos.x + pos.y*vSize.x + pos.z*vSize.x*vSize.y;

	double s = slices[idx];

	if(pos.x >= ssizeX || pos.y >= ssizeY || pos.x < 0 || pos.y < 0)
		return;

	if(s != -1)
	{
		weights[idx] = 1;
	}
	else
	{
		weights[idx] = 0;
	}
	bias[idx] = 0;
}

void Reconstruction::InitializeEMValues()
{
	dim3 blockSize3 = dim3(8,8,8);
	dim3 gridSize3 = divup(dim3(maxSliceDimension.x, maxSliceDimension.y, num_slices_), blockSize3);

	InitializeEMValuesKernel<<<gridSize3, blockSize3>>>(num_slices_, d_slice_sicesX, d_slice_sicesY, 
		v_slices.data, v_bias.data, v_weights.data, make_uint2(v_slices.size.x, v_slices.size.y));
}


#if !TEX_TEST
void Reconstruction::UpdateSliceCoeff(std::vector<SLICECOEFFS> _volcoeffs, double* volWeights)
{
	if(volume_weights_.data == NULL)
		volume_weights_.init(reconstructed_.size, reconstructed_.dim);
	checkCudaErrors(cudaMemcpy(volume_weights_.data, volWeights, 
		volume_weights_.size.x*volume_weights_.size.y*volume_weights_.size.z*sizeof(double), cudaMemcpyHostToDevice));

	unsigned int maxPdepth = 0;
	unsigned int numElems = 0;

	int* hv_volcoeffsOfsSlices_X = new int[v_volcoeffsOfsSlices_X.size.x*v_volcoeffsOfsSlices_X.size.y*v_volcoeffsOfsSlices_X.size.z];
	int* hv_volcoeffsOfsSlices_Y = new int[v_volcoeffsOfsSlices_Y.size.x*v_volcoeffsOfsSlices_Y.size.y*v_volcoeffsOfsSlices_Y.size.z];

	for(int idx = 0; idx < _volcoeffs.size(); idx++)
	{
		for(int i = 0; i < _volcoeffs[idx].size(); i++) 
		{
			for(int j = 0; j < _volcoeffs[idx][i].size(); j++) 
			{
				if(!_volcoeffs[idx][i][j].empty() && _volcoeffs[idx][i][j].size() < 100000) //avoid not inizialized vectors
				{
					hv_volcoeffsOfsSlices_X[i + j*v_volcoeffsOfsSlices_X.size.x + idx*v_volcoeffsOfsSlices_X.size.x*v_volcoeffsOfsSlices_X.size.y] = _volcoeffs[idx][i][j].size();
					hv_volcoeffsOfsSlices_Y[i + j*v_volcoeffsOfsSlices_Y.size.x + idx*v_volcoeffsOfsSlices_Y.size.x*v_volcoeffsOfsSlices_Y.size.y] = numElems;
					numElems += (unsigned int)_volcoeffs[idx][i][j].size();
				}
				else
				{
					hv_volcoeffsOfsSlices_X[i + j*v_volcoeffsOfsSlices_X.size.x] = 0;
					hv_volcoeffsOfsSlices_Y[i + j*v_volcoeffsOfsSlices_Y.size.x] = -1;
				}
			}
		}
	}

	checkCudaErrors(cudaMemcpy(v_volcoeffsOfsSlices_X.data, hv_volcoeffsOfsSlices_X,
		v_volcoeffsOfsSlices_X.size.x*v_volcoeffsOfsSlices_X.size.y*v_volcoeffsOfsSlices_X.size.z*sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(v_volcoeffsOfsSlices_Y.data, hv_volcoeffsOfsSlices_Y,
		v_volcoeffsOfsSlices_Y.size.x*v_volcoeffsOfsSlices_Y.size.y*v_volcoeffsOfsSlices_Y.size.z*sizeof(int), cudaMemcpyHostToDevice));
	CHECK_ERROR(UpdateSliceCoeffcudaMemcpy);

	delete hv_volcoeffsOfsSlices_X;
	delete hv_volcoeffsOfsSlices_Y;

	printf("num elems %llu %Lf MB estimated\n", numElems, (double)(numElems*(sizeof(POINT3D)))/1024.0/1024.0);

	if(volcoeffsarray_ != NULL)
	{
		cudaFree(volcoeffsarray_);
		cudaDeviceSynchronize();
	}
	CHECK_ERROR(UpdateSliceCoeffcudaFree);

	checkCudaErrors(cudaMalloc((void **)&volcoeffsarray_, numElems*sizeof(POINT3D))); 
	CHECK_ERROR(UpdateSliceCoeffcudaMalloc);

	POINT3D* ptr = volcoeffsarray_;
	for(int idx = 0; idx < _volcoeffs.size(); idx++) //_volcoeffs 
	{
		for(int i = 0; i < _volcoeffs[idx].size(); i++) 
		{
			for(int j = 0; j < _volcoeffs[idx][i].size(); j++) 
			{
				if(!_volcoeffs[idx][i][j].empty() && _volcoeffs[idx][i][j].size() < 100000) //avoid not inizialized vectors
				{
					checkCudaErrors(cudaMemcpy(ptr, &(_volcoeffs[idx][i][j][0]), _volcoeffs[idx][i][j].size()*sizeof(POINT3D), cudaMemcpyHostToDevice));
					ptr += _volcoeffs[idx][i][j].size();
					maxPdepth = max((int)maxPdepth, (int)_volcoeffs[idx][i][j].size());
				}
			}
		}
	}
	CHECK_ERROR(UpdateSliceCoeff);
	//cudaDeviceSynchronize();
	maxVolcoDepth = maxPdepth;
	
	//printf("maxPdepth %d\n", maxPdepth);
	checkGPUMemory();
	printf("CoeffInit done\n");
	CHECK_ERROR(UpdateSliceCoeff);
}
#endif
