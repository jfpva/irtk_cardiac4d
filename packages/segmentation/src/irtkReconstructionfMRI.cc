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
	
	if (_debug) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "BsplineS%04iVolume%04i.nii.gz",iter,dyn);
			_timeserie[dyn].Write(buffer);	
		}
	}
}

void irtkReconstructionfMRI::InterpolateBSplineReordered(vector<irtkRealImage>& stacks, vector<int> multiband_vector, int iter) {
	
	// clear timeserie from previous iterations
	_timeserie.clear();
	
	vector<irtkRigidTransformation> currentTransformations;
	vector<irtkRealImage> currentSlices;
	irtkRealImage interpolated;
	irtkImageAttributes attr;
	
	int multiband;
	int counter = 0;
	int position;
	int grouping;
	for (int dyn = 0; dyn < stacks.size()-multiband+1; dyn++)  {
	
		attr = stacks[dyn].GetImageAttributes();
		interpolated = stacks[dyn];
		multiband = multiband_vector[dyn];
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
	
	if (_debug) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "BSplineR%04iVolume%04i.nii.gz",iter,dyn);
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
	
	irtkRealImage interpolated;
	irtkImageAttributes attr, attr2;
	
	int multiband;
	irtkRealImage slice;
	double scale;
	int n;
	POINT3D p;
    
	int counter = 0;
	interpolated  = _reconstructed;
	attr2 = interpolated.GetImageAttributes();
	
    for (int dyn = 0; dyn < stacks.size(); dyn++)  {
	
		attr = stacks[dyn].GetImageAttributes();
			
		CoeffInitSF(counter,counter+attr._z);
		
		// cleaning interpolated
		for (int k = 0; k < attr2._z; k++) {
			for (int j = 0; j < attr2._y; j++) {
				for (int i = 0; i < attr2._x; i++) {
					interpolated(i,j,k) = 0;
				}
			}
		}
		
		for (int s = 0; s < attr._z; s ++) {
			currentTransformations.push_back(_transformations[counter + s]);
			currentSlices.push_back(_slices[counter + s]);
			currentScales.push_back(_scale[counter + s]);
			currentBiases.push_back(_bias[counter + s]);		
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
						n = _volcoeffsSF[s][i][j].size();
						
						//add contribution of current slice voxel to all voxel volumes
						//to which it contributes
						for (int k = 0; k < n; k++) {
							p = _volcoeffsSF[s][i][j][k];
							interpolated(p.x, p.y, p.z) += p.value * slice(i, j, 0);
						}
					}
		}
		
		counter = counter+attr._z;
		currentSlices.clear();
		currentBiases.clear();
		currentTransformations.clear();
		currentScales.clear();
		interpolated /= _volume_weightsSF;
		_timeserie.push_back(interpolated);	
    }
    
	if (_debug) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "GaussianS%04iVolume%04i.nii.gz",iter,dyn);
			_timeserie[dyn].Write(buffer);	
		}
	}
	
}

void irtkReconstructionfMRI::InterpolateGaussianReordered(vector<irtkRealImage>& stacks, vector<int> multiband_vector, int iter) {

	// clear timeserie from previous iterations
	_timeserie.clear();
	
	vector<irtkRigidTransformation> currentTransformations;
	vector<irtkRealImage> currentSlices;
	vector<double> currentScales;
	vector<irtkRealImage> currentBiases;
	
	irtkRealImage interpolated;
	irtkImageAttributes attr, attr2;
	
	int multiband;
	irtkRealImage slice;
	double scale;
	int n;
	POINT3D p;
    
	int counter = 0;
	int stackCounter = 0;
	int sliceIndex = 0;
	interpolated  = _reconstructed;
	attr2 = interpolated.GetImageAttributes();
	
	char buffer[256];
	int grouping;
	int position;
	
	// clear from previous iterations stages
    _slicesRwithMB.clear();
    _transformationsRwithMB.clear();
	for (int dyn = 0; dyn < stacks.size()-multiband+1; dyn++)  {
		
		attr = stacks[dyn].GetImageAttributes();
		multiband = multiband_vector[dyn];
		grouping = attr._z/multiband;
		for (int m = 0; m < multiband; m++) {
			for (int g = 0; g < grouping; g++) {	
				position = counter + g + attr._z*m + grouping*m;
				_slicesRwithMB.push_back(_slices[position]); 
				_transformationsRwithMB.push_back(_transformations[position]);					
			}
		}
		counter = counter+attr._z;
	}
	
	counter = 0;	
    for (int dyn = 0; dyn < stacks.size()-multiband+1; dyn++)  {
    	
    	attr = stacks[dyn].GetImageAttributes();
    	CoeffInitSF(counter,counter+attr._z);
    	
    	// cleaning interpolated
		for (int k = 0; k < attr2._z; k++) {
			for (int j = 0; j < attr2._y; j++) {
				for (int i = 0; i < attr2._x; i++) {
					interpolated(i,j,k) = 0;
				}
			}
		}
		
    	for (int m = 0; m < multiband; m++) {
    		for (int g = 0; g < grouping; g++) {	
    			position = counter + g + attr._z*m + grouping*m;
				currentTransformations.push_back(_transformations[position]);
				currentSlices.push_back(_slices[position]);
				currentScales.push_back(_scale[position]);
				currentBiases.push_back(_bias[position]);						
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
						n = _volcoeffsSF[s][i][j].size();
						
						//add contribution of current slice voxel to all voxel volumes
						//to which it contributes
						for (int k = 0; k < n; k++) {
							p = _volcoeffsSF[s][i][j][k];
							interpolated(p.x, p.y, p.z) += p.value * slice(i, j, 0);
						}
					}
		}
		counter = counter+attr._z;
		currentSlices.clear();
		currentBiases.clear();
		currentTransformations.clear();
		currentScales.clear();
		interpolated /= _volume_weightsSF;
		_timeserie.push_back(interpolated);
    }
    
	if (_debug) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "GaussianR%04iVolume%04i.nii.gz",iter,dyn);
			_timeserie[dyn].Write(buffer);	
		}
	}   
}

void irtkReconstructionfMRI::writefMRI() {
	
	char buffer[256];
	
	irtkImageAttributes attr = _timeserie[0].GetImageAttributes();
	attr._t = _timeserie.size();
	irtkRealImage image(attr);
	for(int t=0;t<attr._t;t++)
	for(int z=0;z<attr._z;z++)
	  for(int y=0;y<attr._y;y++)
		for(int x=0;x<attr._x;x++)
		{
			image(x,y,z,t) = _timeserie[t](x,y,z);
		}
	
	sprintf(buffer,"output.nii.gz");
	image.Write(buffer);
	
}