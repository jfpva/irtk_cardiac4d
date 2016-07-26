/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#ifndef _irtkReconstructionDTI_H

#define _irtkReconstructionDTI_H

#include <irtkReconstruction.h>
#include <irtkSphericalHarmonics.h>


#include <vector>
using namespace std;

/*

  Reconstruction of volume from 2D slices

*/

class irtkReconstructionDTI : public irtkReconstruction
{

 private:
 vector<vector<double> > _directions;
 vector<double> _bvalues;
 irtkRealImage _simulated_signal;
 irtkRealImage _SH_coeffs;
  
 public:
   irtkReconstructionDTI();
   ~irtkReconstructionDTI();
   void GaussianReconstruction4D(int nStacks);
   void GaussianReconstruction4D2(int nStacks);
   void RotateDirections(double &dx, double &dy, double &dz, int i);
   void CreateSliceDirections(vector< vector<double> >& directions, vector<double>& bvalues);
   void InitSH(irtkMatrix dirs,int order);
   void SimulateSlicesDTI();
   void SuperresolutionDTI(int iter);
   void SaveSHcoeffs(int iteration);

   
   friend class ParallelSimulateSlicesDTI;
   friend class ParallelSuperresolutionDTI;
};

#endif
