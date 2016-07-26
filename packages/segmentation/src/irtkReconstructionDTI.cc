/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#include <irtkReconstruction.h>
#include <irtkReconstructionDTI.h>


irtkReconstructionDTI::irtkReconstructionDTI():irtkReconstruction()
{
    _recon_type = _interpolate;

}

irtkReconstructionDTI::~irtkReconstructionDTI() { }

void irtkReconstructionDTI::GaussianReconstruction4D(int nStacks)
{
    cout << "Gaussian reconstruction ... ";
    unsigned int inputIndex;
    int i, j, k, n;
    irtkRealImage slice;
    double scale;
    POINT3D p;
    int dirIndex,origDir;
    double gx,gy,gz,dx,dy,dz,dotp,sigma=0.2,w,tw;

    //clear _reconstructed image
    //_reconstructed = 0;
    irtkImageAttributes attr = _reconstructed.GetImageAttributes();
    attr._t = nStacks;
    cout<<nStacks<<" stacks."<<endl;
    irtkRealImage recon4D(attr);
    irtkRealImage weights(attr);

    for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
      cout<<"Processing slice "<<inputIndex<<" stack "<<_stack_index[inputIndex]<<". ";
        //copy the current slice
        slice = _slices[inputIndex];
        //alias the current bias image
        irtkRealImage& b = _bias[inputIndex];
        //read current scale factor
        scale = _scale[inputIndex];
	
	//direction for current slice
        dirIndex = _stack_index[inputIndex]+1;
	gx=_directions[0][dirIndex];
	gy=_directions[1][dirIndex];
	gz=_directions[2][dirIndex];
	RotateDirections(gx,gy,gz,inputIndex);
	b=_bvalues[dirIndex];
	
	/*
	//cout<<"Dot product is ";
        //Loop over directions
	for (origDir=1; origDir<_bvalues.size();origDir++)
	{
	  
	  dx=_directions[0][origDir];
	  dy=_directions[1][origDir];
	  dz=_directions[2][origDir];
	  
	  dotp = (dx*gx+dy*gy+dz*gz)/sqrt((dx*dx+dy*dy+dz*dz)*(gx*gx+gy*gy+gz*gz));
	  //cout<<"origDir="<<origDir<<": "<<dotp<<"; ";
	  tw = (fabs(dotp)-1)/sigma;
	  w=exp(-tw*tw)/(6.28*sigma);
	  //cout<<"weight = "<<w<<"; ";
	  */
	
	  //loop over slice voxels
          for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    //biascorrect and scale the slice
                    slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;

                    //number of volume voxels with non-zero coefficients
                    //for current slice voxel
                    n = _volcoeffs[inputIndex][i][j].size();

                    //add contribution of current slice voxel to all voxel volumes
                    //to which it contributes
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
                        //recon4D(p.x, p.y, p.z,origDir-1) += p.value * slice(i, j, 0);
                        //weights(p.x, p.y, p.z,origDir-1) += p.value;
                        recon4D(p.x, p.y, p.z,_stack_index[inputIndex]) += p.value * slice(i, j, 0);
                        weights(p.x, p.y, p.z,_stack_index[inputIndex]) += p.value;
			
                    }
                }

	//} //end of loop for direction origDir  
	cout<<endl;
    }//end of loop for a slice inputIndex

    //normalize the volume by proportion of contributing slice voxels
    //for each volume voxe
    recon4D /= weights;
        
    cout << "done." << endl;

    if (_debug)
    {
        recon4D.Write("recon4D.nii.gz");
	weights.Write("weights.nii.gz");
    }

}

void irtkReconstructionDTI::GaussianReconstruction4D2(int nStacks)
{
    cout << "Gaussian reconstruction ... ";
    unsigned int inputIndex;
    int i, j, k, n;
    irtkRealImage slice;
    double scale;
    POINT3D p;
    int dirIndex, origDir;
    double bval,gx,gy,gz,dx,dy,dz,dotp,sigma=0.02,w,tw;

    irtkImageAttributes attr = _reconstructed.GetImageAttributes();
    attr._t = nStacks;
    cout<<nStacks<<" stacks."<<endl;
    irtkRealImage recon4D(attr);
    irtkRealImage weights(attr);

    for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
      //copy the current slice
      slice = _slices[inputIndex];
      //alias the current bias image
      irtkRealImage& b = _bias[inputIndex];
      //read current scale factor
      scale = _scale[inputIndex];
        
      //direction for current slice
      dirIndex = _stack_index[inputIndex]+1;
      gx=_directions[0][dirIndex];
      gy=_directions[1][dirIndex];
      gz=_directions[2][dirIndex];
      RotateDirections(gx,gy,gz,inputIndex);
      bval=_bvalues[dirIndex];
	
      for (origDir=1; origDir<_bvalues.size();origDir++)
      {
	//cout<<origDir<<" ";
	
	  dx=_directions[0][origDir];
	  dy=_directions[1][origDir];
	  dz=_directions[2][origDir];
	  
	  dotp = (dx*gx+dy*gy+dz*gz)/sqrt((dx*dx+dy*dy+dz*dz)*(gx*gx+gy*gy+gz*gz));
	  //cout<<"origDir="<<origDir<<": "<<dotp<<"; ";
	  tw = (fabs(dotp)-1)/sigma;
	  w=exp(-tw*tw)/(6.28*sigma);
	  //cout<<"weight = "<<w<<"; ";

	
        //Distribute slice intensities to the volume
        for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    //biascorrect and scale the slice
                    slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;

                    //number of volume voxels with non-zero coefficients
                    //for current slice voxel
                    n = _volcoeffs[inputIndex][i][j].size();

                    //add contribution of current slice voxel to all voxel volumes
                    //to which it contributes
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
                        recon4D(p.x, p.y, p.z,origDir-1) += _slice_weight[inputIndex] * w * p.value * slice(i, j, 0);
                        weights(p.x, p.y, p.z,origDir-1) += _slice_weight[inputIndex] * w * p.value;
                    }
                }
      } //end of loop for origDir
      //cout<<endl;
        //end of loop for a slice inputIndex
    }

    //normalize the volume by proportion of contributing slice voxels
    //for each volume voxe
   recon4D /= weights;
   _simulated_signal = recon4D;
        
    cout << "done." << endl;

    if (_debug)
    {
        recon4D.Write("recon4Dgaussian.nii.gz");
        weights.Write("weights.nii.gz");
    }

}

class ParallelSimulateSlicesDTI {
    irtkReconstructionDTI *reconstructor;
        
public:
    ParallelSimulateSlicesDTI( irtkReconstructionDTI *_reconstructor ) : 
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
	    
	    
	    //direction for current slice
            int dirIndex = reconstructor->_stack_index[inputIndex]+1;
            double gx=reconstructor->_directions[0][dirIndex];
            double gy=reconstructor->_directions[1][dirIndex];
            double gz=reconstructor->_directions[2][dirIndex];
            reconstructor->RotateDirections(gx,gy,gz,inputIndex);
            double bval=reconstructor->_bvalues[dirIndex];
	    irtkSphericalHarmonics sh;
	    irtkMatrix dir(1,3);
	    dir(0,0)=gx;
	    dir(0,1)=gy;
	    dir(0,2)=gz;
	    irtkMatrix basis = sh.SHbasis(dir,2);
	    if(basis.Cols() != reconstructor->_SH_coeffs.GetT())
	    {
	      cerr<<"ParallelSimulateSlicesDTI:basis numbers does not match SH coefficients number."<<endl;
	      exit(1);
	    }
	    double sim_signal;
            
            POINT3D p;
            for ( unsigned int i = 0; i < reconstructor->_slices[inputIndex].GetX(); i++ )
                for ( unsigned int j = 0; j < reconstructor->_slices[inputIndex].GetY(); j++ )
                    if ( reconstructor->_slices[inputIndex](i, j, 0) != -1 ) {
                        double weight = 0;
                        int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                        for ( unsigned int k = 0; k < n; k++ ) {
			     //PSF
                            p = reconstructor->_volcoeffs[inputIndex][i][j][k];
			    //signal simulated from SH
			     sim_signal = 0;
			     for(unsigned int l = 0; l < basis.Cols(); l++ )
			       sim_signal += reconstructor->_SH_coeffs(p.x, p.y, p.z,l)*basis(0,l);
			     //update slice
                            reconstructor->_simulated_slices[inputIndex](i, j, 0) += p.value * sim_signal;
                            weight += p.value;
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

void irtkReconstructionDTI::SimulateSlicesDTI()
{
    if (_debug)
        cout<<"Simulating slices DTI."<<endl;

    ParallelSimulateSlicesDTI parallelSimulateSlicesDTI( this );
    parallelSimulateSlicesDTI();

    if (_debug)
        cout<<"done."<<endl;    
}



class ParallelSuperresolutionDTI {
    irtkReconstructionDTI* reconstructor;
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
	    
	    //direction for current slice
            int dirIndex = reconstructor->_stack_index[inputIndex]+1;
            double gx=reconstructor->_directions[0][dirIndex];
            double gy=reconstructor->_directions[1][dirIndex];
            double gz=reconstructor->_directions[2][dirIndex];
            reconstructor->RotateDirections(gx,gy,gz,inputIndex);
            double bval=reconstructor->_bvalues[dirIndex];
	    irtkSphericalHarmonics sh;
	    irtkMatrix dir(1,3);
	    dir(0,0)=gx;
	    dir(0,1)=gy;
	    dir(0,2)=gz;
	    irtkMatrix basis = sh.SHbasis(dir,2);
	    if(basis.Cols() != reconstructor->_SH_coeffs.GetT())
	    {
	      cerr<<"ParallelSimulateSlicesDTI:basis numbers does not match SH coefficients number."<<endl;
	      exit(1);
	    }

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
			    if(reconstructor->_robust_slices_only)
			    {
			      for(unsigned int l = 0; l < basis.Cols(); l++ )
			      {
                                addon(p.x, p.y, p.z,l) += p.value * basis(0,l) * slice(i, j, 0) * reconstructor->_slice_weight[inputIndex];
                                confidence_map(p.x, p.y, p.z,l) += p.value *reconstructor->_slice_weight[inputIndex];
			      }
			      
			    }
			    else
			    {
                              for(unsigned int l = 0; l < basis.Cols(); l++ )
			       {
				  addon(p.x, p.y, p.z, l) += p.value * basis(0,l) * slice(i, j, 0) * w(i, j, 0) * reconstructor->_slice_weight[inputIndex];
                                confidence_map(p.x, p.y, p.z, l) += p.value * w(i, j, 0) * reconstructor->_slice_weight[inputIndex];
                                //p.value * basis(0,l) * w(i, j, 0) * reconstructor->_slice_weight[inputIndex];
			       }
			    }
                        }
                    }
        } //end of loop for a slice inputIndex
    }
 
    ParallelSuperresolutionDTI( ParallelSuperresolutionDTI& x, split ) :
        reconstructor(x.reconstructor)
    {
        //Clear addon
        addon.Initialize( reconstructor->_SH_coeffs.GetImageAttributes() );
        addon = 0;

        //Clear confidence map
        confidence_map.Initialize( reconstructor->_SH_coeffs.GetImageAttributes() );
        confidence_map = 0;
    }
 
    void join( const ParallelSuperresolutionDTI& y ) {
        addon += y.addon;
        confidence_map += y.confidence_map;
    }
             
    ParallelSuperresolutionDTI( irtkReconstructionDTI *reconstructor ) :
    reconstructor(reconstructor)
    {
        //Clear addon
        addon.Initialize( reconstructor->_SH_coeffs.GetImageAttributes() );
        addon = 0;

        //Clear confidence map
        confidence_map.Initialize( reconstructor->_SH_coeffs.GetImageAttributes() );
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

void irtkReconstructionDTI::SuperresolutionDTI(int iter)
{
    if (_debug)
        cout << "SuperresolutionDTI " << iter << endl;
    
    int i, j, k, l;
    irtkRealImage addon, original;
    

    //Remember current reconstruction for edge-preserving smoothing
    //original = _reconstructed;

    ParallelSuperresolutionDTI parallelSuperresolutionDTI(this);
    parallelSuperresolutionDTI();
    addon = parallelSuperresolutionDTI.addon;
    _confidence_map = parallelSuperresolutionDTI.confidence_map;
    //_confidence4mask = _confidence_map;
    
    if(_debug) {
        char buffer[256];
        //sprintf(buffer,"confidence-map%i.nii.gz",iter);
        //_confidence_map.Write(buffer);
	_confidence_map.Write("confidence-map-superDTI.nii.gz");
        sprintf(buffer,"addon%i.nii.gz",iter);
        addon.Write(buffer);
    }

    if (!_adaptive) 
        for (i = 0; i < addon.GetX(); i++)
            for (j = 0; j < addon.GetY(); j++)
                for (k = 0; k < addon.GetZ(); k++)
                     for (l = 0; l < addon.GetT(); l++)
                     if (_confidence_map(i, j, k, l) > 0) {
                        // ISSUES if _confidence_map(i, j, k) is too small leading
                        // to bright pixels
                        addon(i, j, k,l) /= _confidence_map(i, j, k, l);
                        //this is to revert to normal (non-adaptive) regularisation
                        _confidence_map(i,j,k,l) = 1;
                    }

    _SH_coeffs += addon * _alpha * _average_volume_weight; //_average_volume_weight;
  /*  
    //bound the intensities
    for (i = 0; i < _reconstructed.GetX(); i++)
        for (j = 0; j < _reconstructed.GetY(); j++)
            for (k = 0; k < _reconstructed.GetZ(); k++) {
                if (_reconstructed(i, j, k) < _min_intensity * 0.9)
                    _reconstructed(i, j, k) = _min_intensity * 0.9;
                if (_reconstructed(i, j, k) > _max_intensity * 1.1)
                    _reconstructed(i, j, k) = _max_intensity * 1.1;
            }

    //Smooth the reconstructed image
    AdaptiveRegularization(iter, original);
    //Remove the bias in the reconstructed volume compared to previous iteration
    if (_global_bias_correction)
        BiasCorrectVolume(original);
*/
}



void irtkReconstructionDTI::RotateDirections(double &dx, double &dy, double &dz, int i)
{

  //vector end-point
  double x,y,z;
  //origin
  double ox,oy,oz;
  
  if (_debug)
  {
    //cout<<"Original direction "<<i<<"(dir"<<_stack_index[i]+1<<"): ";
    //cout<<dx<<", "<<dy<<", "<<dz<<". ";
    //cout<<endl;
  }

  //origin
  ox=0;oy=0;oz=0;
  _transformations[i].Transform(ox,oy,oz);
    
  //end-point
  x=dx;
  y=dy;
  z=dz;
  _transformations[i].Transform(x,y,z);
    
  dx=x-ox;
  dy=y-oy;
  dz=z-oz;
    
  if (_debug)
  {
    //cout<<"Rotated direction "<<i<<"(dir"<<_stack_index[i]+1<<"): ";
    //cout<<dx<<", "<<dy<<", "<<dz<<". ";
    //cout<<endl;
  }
}

void irtkReconstructionDTI::CreateSliceDirections(vector< vector<double> >& directions, vector<double>& bvalues)
{
  _directions = directions;
  _bvalues = bvalues;
  
  cout<<"B-values: ";
  for(uint i=0;i<_bvalues.size(); i++)
    cout<<_bvalues[i]<<" ";
  cout<<endl;

  cout<<"B-vectors: ";
  for(uint j=0; j<_directions.size(); j++)
  {
    for(uint i=0;i<_directions[j].size(); i++)
      cout<<_directions[j][i]<<" ";
    cout<<endl;
  }
  cout<<endl;
}

void irtkReconstructionDTI::InitSH(irtkMatrix dirs, int order)
{
  irtkSphericalHarmonics sh;
  sh.InitSHT(dirs,order);//lambda,order);
  _SH_coeffs = sh.Signal2Coeff(_simulated_signal);
  _SH_coeffs.Write("_SH_coeffs.nii.gz");
}


void irtkReconstructionDTI::SaveSHcoeffs(int iteration)
{
  char buffer[256];
  sprintf(buffer, "shCoeff%i.nii.gz", iteration);
  _SH_coeffs.Write(buffer);
}

