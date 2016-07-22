/*=========================================================================

 Library   : Image Registration Toolkit (IRTK)
 Module    : $Id: irtkReconstructionfMRI.cc 837 2013-05-07 12:55:31Z mm3 $
 Copyright : Imperial College, Department of Computing
 Visual Information Processing (VIP), 2011 onwards
 Date      : $Date: 2013-05-07 13:55:31 +0100 (Tue, 07 May 2013) $
 Version   : $Revision: 837 $
 Changes   : $Author: mm3 $

 =========================================================================*/

#include <irtkReconstruction.h>
#include <irtkReconstructionb0.h>
#include <irtkReconstructionfMRI.h>

irtkReconstructionfMRI::irtkReconstructionfMRI()
{
	irtkReconstructionb0();
}

void irtkReconstructionfMRI::InterpolateBSpline(vector<irtkRealImage>& stacks, int iter) {

	// clear timeserie from previous iterations
	_timeserie.clear();
	
	vector<irtkRigidTransformation> currentTransformations;
	vector<irtkRealImage> currentSlices;
	irtkRealImage interpolated;
	irtkImageAttributes attr;
	
	int counter = 0;
	for (int dyn = 0; dyn < stacks.size(); dyn++)  {
	
		attr = stacks[dyn].GetImageAttributes();
		interpolated = stacks[dyn];
		
		for (int slice = 0; slice < attr._z; slice ++) {
			currentTransformations.push_back(_transformations[counter + slice]);
			currentSlices.push_back(_slices[counter + slice]);
		}
		
		_bSplineReconstruction.Reconstruct(6,1,interpolated,currentSlices,currentTransformations); // to be tuned 	
		_timeserie.push_back(interpolated);
		counter = counter + attr._z;
		currentTransformations.clear();
		currentSlices.clear();
	}
	
	if (true) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "BsplineS%04iVolume%04i.nii.gz",iter,dyn);
			_timeserie[dyn].Write(buffer);	
		}
	}
}

void irtkReconstructionfMRI::InterpolateGaussian(vector<irtkRealImage>& stacks, int iter) {

	// clear timeserie from previous iterations
	_timeserie.clear();
	
	vector<irtkRigidTransformation> currentTransformations;
	vector<irtkRealImage> currentSlices;
	vector<double> currentScales;
	vector<irtkRealImage> currentBiases;
	irtkRealImage volumeWeights;
	
	irtkRealImage interpolated;
	irtkImageAttributes attr, attr2;
	
	irtkRealImage slice;
	double scale;
	int n;
	POINT3D p;
	int slice_vox_num = 0;
    
	int counter = 0;
	interpolated  = _reconstructed;
	volumeWeights = _reconstructed;
	
    for (int dyn = 0; dyn < stacks.size(); dyn++)  {
	
		attr = stacks[dyn].GetImageAttributes();
		attr2 = interpolated.GetImageAttributes();
		
		for (int s = 0; s < attr._z; s ++) {
			currentTransformations.push_back(_transformations[counter + s]);
			currentSlices.push_back(_slices[counter + s]);
			currentScales.push_back(_scale[counter + s]);
			currentBiases.push_back(_bias[counter + s]);
		}

		// cleaning interpolated and volumeWeights
		for (int k = 0; k < attr2._z; k++) {
			for (int j = 0; j < attr2._y; j++) {
				for (int i = 0; i < attr2._x; i++) {
					interpolated(i,j,k) = 0;
					volumeWeights(i,j,k) = 0;
				}
			}
		}
		
		// weights
		for ( int s = 0; s < attr._z; s ++ ) {
			for ( int i = 0; i < _slices[counter + s].GetX(); i++)
				for ( int j = 0; j < _slices[counter + s].GetY(); j++) {
					n = _volcoeffs[counter + s][i][j].size();
					for (int k = 0; k < n; k++) {			
						p = _volcoeffs[counter + s][i][j][k];
						volumeWeights(p.x, p.y, p.z) += p.value;
					}
			}
		}
		
		for (int s = 0; s < currentSlices.size(); s++) {
			
			//copy the current slice
			slice = currentSlices[s];
			//alias the current bias image
			irtkRealImage& b = currentBiases[s];
			//read current scale factor
			scale = currentScales[s];
			
			for (int i = 0; i < slice.GetX(); i++)
				for (int j = 0; j < slice.GetY(); j++)
					if (slice(i, j, 0) != -1) {
						//biascorrect and scale the slice
						slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;
			
						//number of volume voxels with non-zero coefficients
						//for current slice voxel
						n = _volcoeffs[counter + s][i][j].size();
			
						//if given voxel is not present in reconstructed volume at all,
						//pad it
						
						//if (n == 0)
						//_slices[inputIndex].PutAsDouble(i, j, 0, -1);
						//calculate num of vox in a slice that have overlap with roi
						if (n>0)
							slice_vox_num++;
			
						//add contribution of current slice voxel to all voxel volumes
						//to which it contributes
						for (int k = 0; k < n; k++) {
							p = _volcoeffs[counter + s][i][j][k];
							interpolated(p.x, p.y, p.z) += p.value * slice(i, j, 0);
						}
					}
		}
		counter = counter + attr._z;
		currentSlices.clear();
		currentBiases.clear();
		currentTransformations.clear();
		currentScales.clear();
		interpolated /= volumeWeights;
		_timeserie.push_back(interpolated);
    }
    
    if (true) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "GaussianS%04iVolume%04i.nii.gz",iter,dyn);
			_timeserie[dyn].Write(buffer);	
		}
	}
}

void irtkReconstructionfMRI::InterpolateBSplineReordered(vector<irtkRealImage>& stacks, int multiband, int iter) {
	
	// clear timeserie from previous iterations
	_timeserie.clear();
	
	vector<irtkRigidTransformation> currentTransformations;
	vector<irtkRealImage> currentSlices;
	irtkRealImage interpolated;
	irtkImageAttributes attr;
	
	int counter = 0;
	int position;
	int grouping;
	for (int dyn = 0; dyn < stacks.size()-multiband+1; dyn++)  {
	
		attr = stacks[dyn].GetImageAttributes();
		interpolated = stacks[dyn];
		grouping = attr._z/multiband;
		
		for (int m = 0; m < multiband; m++) {
			for (int g = 0; g < grouping; g++) {
				position = counter + g + attr._z*m + grouping*m;
				currentTransformations.push_back(_transformations[position]);
				currentSlices.push_back(_slices[position]);
			}
		}
		
		_bSplineReconstruction.Reconstruct(6,1,interpolated,currentSlices,currentTransformations); // to be tuned 	
		_timeserie.push_back(interpolated);
		counter = counter + attr._z;
		currentTransformations.clear();
		currentSlices.clear();
	}
	
	if (true) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "BSplineR%04iVolume%04i.nii.gz",iter,dyn);
			_timeserie[dyn].Write(buffer);	
		}
	}
}

void irtkReconstructionfMRI::InterpolateGaussianReordered(vector<irtkRealImage>& stacks, int multiband, int iter) {

	// clear timeserie from previous iterations
	_timeserie.clear();
	
	vector<irtkRigidTransformation> currentTransformations;
	vector<irtkRealImage> currentSlices;
	vector<double> currentScales;
	vector<irtkRealImage> currentBiases;
	irtkRealImage volumeWeights;
	
	irtkRealImage interpolated;
	irtkImageAttributes attr, attr2;
	
	irtkRealImage slice;
	double scale;
	int n;
	POINT3D p;
	int slice_vox_num = 0;
    
	int counter = 0;
	int sliceIndex = 0;
	interpolated  = _reconstructed;
	volumeWeights = _reconstructed;

	char buffer[256];
	int grouping;
	int position;
    for (int dyn = 0; dyn < stacks.size()-multiband+1; dyn++)  {
	
		attr = stacks[dyn].GetImageAttributes();
		attr2 = interpolated.GetImageAttributes();
		grouping = attr._z/multiband;
		
		for (int m = 0; m < multiband; m++) {
			for (int g = 0; g < grouping; g++) {	
				position = counter + g + attr._z*m + grouping*m;
				currentTransformations.push_back(_transformations[position]);
				currentSlices.push_back(_slices[position]);
				currentScales.push_back(_scale[position]);
				currentBiases.push_back(_bias[position]);
			}
		}

		// cleaning interpolated and volumeWeights
		for (int k = 0; k < attr2._z; k++) {
			for (int j = 0; j < attr2._y; j++) {
				for (int i = 0; i < attr2._x; i++) {
					interpolated(i,j,k) = 0;
					volumeWeights(i,j,k) = 0;
				}
			}
		}
		
		for (int m = 0; m < multiband; m++) {
			for (int g = 0; g < grouping; g++) {	
				position = counter + g + attr._z*m + grouping*m;
				for ( int i = 0; i < _slices[position].GetX(); i++)
					for ( int j = 0; j < _slices[position].GetY(); j++) {
						n = _volcoeffs[position][i][j].size();
						for (int k = 0; k < n; k++) {			
							p = _volcoeffs[position][i][j][k];
							volumeWeights(p.x, p.y, p.z) += p.value;
						}
				}
			}
		}

		for (int m = 0; m < multiband; m++) {
			for (int g = 0; g < grouping; g++) {	
				
				position = counter + g + attr._z*m + grouping*m;
				
				//copy the current slice
				slice = currentSlices[sliceIndex];
				//alias the current bias image
				irtkRealImage& b = currentBiases[sliceIndex];
				//read current scale factor
				scale = currentScales[sliceIndex];
				
				for (int i = 0; i < slice.GetX(); i++)
					for (int j = 0; j < slice.GetY(); j++)
						if (slice(i, j, 0) != -1) {
							//biascorrect and scale the slice
							slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;
				
							//number of volume voxels with non-zero coefficients
							//for current slice voxel
							n = _volcoeffs[position][i][j].size();
				
							//if given voxel is not present in reconstructed volume at all,
							//pad it
							
							//if (n == 0)
							//_slices[inputIndex].PutAsDouble(i, j, 0, -1);
							//calculate num of vox in a slice that have overlap with roi
							if (n>0)
								slice_vox_num++;
				
							//add contribution of current slice voxel to all voxel volumes
							//to which it contributes
							for (int k = 0; k < n; k++) {
								p = _volcoeffs[position][i][j][k];
								interpolated(p.x, p.y, p.z) += p.value * slice(i, j, 0);
							}
						}
				sliceIndex++;
				}
		}
		sliceIndex = 0;
		counter = counter + attr._z;
		currentSlices.clear();
		currentBiases.clear();
		currentTransformations.clear();
		currentScales.clear();
		interpolated /= volumeWeights;
		_timeserie.push_back(interpolated);		
    }
    if (true) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "GuassianR%04iVolume%04i.nii.gz",iter,dyn);
			_timeserie[dyn].Write(buffer);	
		}
	}
}