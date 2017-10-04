/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#ifndef _irtkReconstructionCardiac4D_H

#define _irtkReconstructionCardiac4D_H

#include <irtkReconstruction.h>

#include <vector>

using namespace std;

/*

  Reconstruction of 4D cardiac cine volume from 2D slices

*/

class irtkReconstructionCardiac4D : public irtkReconstruction
{

protected:
   
  // Stacks
  vector<int> _loc_index;        // running index of all 2D slice locations
  vector<int> _stack_loc_index;  // index of 2D slice location in M2D stack
  vector<int> _stack_dyn_index;  // index of dynamic in M2D stack

   // PI
   const double PI = 3.14159265358979323846;
   
   // Reconstructed 4D Cardiac Cine Image
   irtkRealImage _reconstructed4D;  // TODO: replace _reconstructed4D with _reconstructed and fix conflicts between irtkReconstruction and irtkReconstructionCardiac4D use of _reconstruction and_reconstruction4D
    
   // Reconstructed Cardiac Phases
   vector<double> _reconstructed_cardiac_phases;
 
   // Reconstructe R-R Interval
   double _reconstructed_rr_interval;
 
   // Reconstructed Temporal Resolution
   double _reconstructed_temporal_resolution;
 
   // Slice Acquisition Time
   vector<double> _slice_time;
   
   // Slice Temporal Resolution
   vector<double> _slice_dt;
   
   // Slice R-R interval
   vector<double> _slice_rr; 
   
   // Slice Cardiac Phase
   vector<double> _slice_cardphase;
   
   // Slice Temporal Weight
   // access as: _slice_temporal_weight[iReconstructedCardiacPhase][iSlice]
   vector< vector<double> > _slice_temporal_weight; 
   
   // Slice SVR Target Cardiac Phase
   vector<int> _slice_svr_card_index;
  
   // Initialise Slice Temporal Weights
   void InitSliceTemporalWeights();
   
   // Calculate Angular Difference
   double CalculateAngularDifference( double cardphase0, double cardphase );

   // Temporal Weighting Tukey Window Edge Percent
   double _wintukeypct = 0.3;
   
   // Sinc Function
   double sinc( double x );

   // Tukey Window Function
   double wintukey( double angdiff, double alpha );
   
   // Calculate Temporal Weight
   double CalculateTemporalWeight( double cardphase0, double cardphase, double dt, double rr, double alpha );
  
 public:
   
   // Constructor
   irtkReconstructionCardiac4D();
   
   // Destructor
   ~irtkReconstructionCardiac4D();
   
   // Get Reconstructed 4D Volume
   inline irtkRealImage GetReconstructedCardiac4D();

   // Get Volume Weights
   inline irtkRealImage GetVolumeWeights();

   // Get Slice-Location transformations
   void ReadSliceTransformation(char* slice_transformations_folder);

   // Set Slice R-R Intervals
   void SetSliceRRInterval( vector<double> rr );
   void SetSliceRRInterval( double rr );
   
   // Set Loc R-R Intervals
   void SetLocRRInterval( vector<double> rr );

   // Set Slice Cardiac Phases
   void SetSliceCardiacPhase( vector<double> cardiacphases );
   void SetSliceCardiacPhase();
   
   // Set Reconstructed R-R Interval
   void SetReconstructedRRInterval( double rrinterval );
   
   // Set Reconstructed Temporal Resolution
   void SetReconstructedTemporalResolution( double temporalresolution );
   
   // Set Reconstructed Cardiac Phases
   void SetReconstructedCardiacPhase( vector<double> cardiacphases );

   // Set Reconstructed Volume Spatial Resolution
   double GetReconstructedResolutionFromTemplateStack( irtkRealImage stack );
   
   // Initialise Reconstructed Volume from Static Mask
   void CreateTemplateCardiac4DFromStaticMask( irtkRealImage mask, double resolution );
  
   ///Match stack intensities with masking
   void MatchStackIntensitiesWithMasking (vector<irtkRealImage>& stacks,
                               vector<irtkRigidTransformation>& stack_transformations,
                               double averageValue,
                               bool together=false);
                               
   // Create slices from the stacks and 
   // slice-dependent transformations from stack transformations
   void CreateSlicesAndTransformationsCardiac4D( vector<irtkRealImage>& stacks,
                                        vector<irtkRigidTransformation>& stack_transformations,
                                        vector<double>& thickness,
                                        const vector<irtkRealImage> &probability_maps=vector<irtkRealImage>() );
  
   // TODO: Calculate Cardiac Phase from Trigger Times
   // void CalculateSliceCardiacPhases( vector<double>& trigger_times );
   
   // Calculate Slice Temporal Weights
   void CalculateSliceTemporalWeights();

   void TestTemporalWeightCalculation();

   
   // Calculate Transformation Matrix Between Slices and Voxels
   void CoeffInitCardiac4D();

   // PSF-Weighted Reconstruction
   void GaussianReconstructionCardiac4D();
   
   ///Scale volume to match the slice intensities
   void ScaleVolumeCardiac4D();
   
   // Simulate Slices
   void SimulateSlicesCardiac4D();
   
   // Normalise Bias
   void NormaliseBiasCardiac4D(int iter, int rec_inter);
  
   // Slice-to-Volume Registration
   void CalculateSliceToVolumeTargetCardiacPhase();
   void SliceToVolumeRegistrationCardiac4D();
  
   // Apply Static Mask to Reconstructed 4D Volume
   void StaticMaskReconstructedVolume4D();
   
   // Apply Static Mask 4D Volume
   irtkRealImage StaticMaskVolume4D(irtkRealImage volume, double padding);
   
   // Superresolution 
   void SuperresolutionCardiac4D( int iter );

   /// Edge-Preserving Regularization with Confidence Map
   void AdaptiveRegularizationCardiac4D(int iter, irtkRealImage& original);

   
   /// Read Transformations
   void ReadTransformation( char* folder );

   /// Calculate Entropy
   double CalculateEntropy();

   ///Save transformations
   void SaveTransformations();

   // Save Bias Fields
   void SaveBiasFields();
   void SaveBiasFields(vector<irtkRealImage>& stacks);
   void SaveBiasFields(vector<irtkRealImage>& stacks, int iter, int rec_iter);

   // Save Slices
   void SaveSimulatedSlices();
   void SaveSimulatedSlices(vector<irtkRealImage>& stacks);
   void SaveSimulatedSlices(vector<irtkRealImage>& stacks, int iter, int rec_iter);
   
   // Save Slices
   void SaveSlices();
   void SaveSlices(vector<irtkRealImage>& stacks);
   
   // Save Weights
   void SaveWeights();
   void SaveWeights(vector<irtkRealImage>& stacks);
   void SaveWeights(vector<irtkRealImage>& stacks, int iter, int rec_iter);

   // Slice Info 
   void SlicesInfoCardiac4D( const char* filename, vector<string> &stack_filenames );

   // Access to Parallel Processing Classes
   friend class ParallelCoeffInitCardiac4D;
   friend class ParallelSliceToVolumeRegistrationCardiac4D;
   friend class ParallelSimulateSlicesCardiac4D;
   friend class ParallelNormaliseBiasCardiac4D;
   friend class ParallelSuperresolutionCardiac4D;   
   friend class ParallelAdaptiveRegularization1Cardiac4D;
   friend class ParallelAdaptiveRegularization2Cardiac4D;
   
};  // end of irtReconstructionCardiac4D class definition


// -----------------------------------------------------------------------------
// Get Reconstructed 4D Volume
// -----------------------------------------------------------------------------
inline irtkRealImage irtkReconstructionCardiac4D::GetReconstructedCardiac4D()
{
    return _reconstructed4D;
}


// -----------------------------------------------------------------------------
// Get Volume Weights
// -----------------------------------------------------------------------------
inline irtkRealImage irtkReconstructionCardiac4D::GetVolumeWeights()
{
    return _volume_weights;
}


#endif
