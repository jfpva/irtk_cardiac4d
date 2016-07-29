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
	
		CoeffInitfMRI(counter,counter + attr._z -1);
		
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
	
	if (true) {
		char buffer[256];
		for (int dyn = 0; dyn < _timeserie.size(); dyn++) {
			sprintf(buffer, "BSplineR%04iVolume%04i.nii.gz",iter,dyn);
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
	irtkRealImage volumeWeights;
	
	irtkRealImage interpolated;
	irtkImageAttributes attr, attr2;
	
	int multiband;
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
		multiband = multiband_vector[dyn];
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
			sprintf(buffer, "GaussianR%04iVolume%04i.nii.gz",iter,dyn);
			_timeserie[dyn].Write(buffer);	
		}
	}
}

class ParallelCoeffInitfMRI {
public:
    irtkReconstruction *reconstructor;
    size_t _begin;
    size_t _end;

    ParallelCoeffInitfMRI(irtkReconstructionfMRI *_reconstructor, size_t begin, size_t end) : 
    reconstructor(_reconstructor) 
    { 
    	_begin = begin;
    	_end = end;
    }

    void operator() (const blocked_range<size_t> &r) const {
        
        for ( size_t inputIndex = r.begin(); inputIndex != r.end(); ++inputIndex ) {

            bool slice_inside;
            cerr << inputIndex << " ";
            //current slice
            //irtkRealImage slice;

            //get resolution of the volume
            double vx, vy, vz;
            reconstructor->_reconstructed.GetPixelSize(&vx, &vy, &vz);
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
            cerr<<" Done "<<inputIndex<<endl;
        }  //end of loop through the slices                            

    }

    // execute
    void operator() () const {
    	cerr<<"Bounds: "<<_begin<<" "<<_end<<endl;
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(_begin, _end +1 ),
                      *this );
        init.terminate();
    }

};

void irtkReconstructionfMRI::CoeffInitfMRI(int begin, int end)
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
	
    ParallelCoeffInitfMRI coeffinit(this,begin,end);
    coeffinit();
    cerr << " ... done." << endl;
    /*
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
    
    if(_debug) {
        cout<<"Average volume weight is "<<_average_volume_weight<<endl;
    }
    */
}  //end of CoeffInit()