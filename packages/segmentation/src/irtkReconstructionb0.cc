/*=========================================================================

 Library   : Image Registration Toolkit (IRTK)
 Module    : $Id: irtkReconstructionb0.cc 837 2013-05-07 12:55:31Z mm3 $
 Copyright : Imperial College, Department of Computing
 Visual Information Processing (VIP), 2011 onwards
 Date      : $Date: 2013-05-07 13:55:31 +0100 (Tue, 07 May 2013) $
 Version   : $Revision: 837 $
 Changes   : $Author: mm3 $

 =========================================================================*/

#include <irtkReconstruction.h>
#include <irtkReconstructionb0.h>
#include <irtkResampling.h>
#include <irtkRegistration.h>
#include <irtkImageRigidRegistration.h>
#include <irtkImageRigidRegistrationWithPadding.h>
#include <irtkImageAffineRegistrationWithPadding.h>
#include <irtkImageFreeFormRegistrationWithPadding.h>
#include <irtkTransformation.h>
#include <irtkLaplacianSmoothing.h>

irtkReconstructionb0::irtkReconstructionb0()
{
  _have_larger_mask = false;
}

void irtkReconstructionb0::StackRegistrations(vector<irtkRealImage>& stacks,
		vector<irtkRigidTransformation>& stack_transformations)
{
      InvertStackTransformations(stack_transformations);
	//rigid registration object
	irtkImageRigidRegistrationWithPadding registration;
	//buffer to create the name
	char buffer[256];

	//template is set as the target
	irtkGreyImage target = _reconstructed;
	//target needs to be masked before registration
	if(_debug)
	  target.Write("target-nomask.nii.gz");
	if (_have_mask)
	{
		double x, y, z;
		for (int i = 0; i < target.GetX(); i++)
			for (int j = 0; j < target.GetY(); j++)
				for (int k = 0; k < target.GetZ(); k++)
				{
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
							&& (z < _mask.GetZ()))
							{
						if (_mask(x, y, z) == 0)
							target(i, j, k) = 0;
					}
					else
						target(i, j, k) = 0;
				}
	}

        if(_debug)
          target.Write("target.nii.gz");
        irtkRigidTransformation offset;
	ResetOrigin(target,offset);

	//register all stacks to the target
	for (int i = 0; i < (int)stacks.size(); i++)
	{
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


		//save volumetric registrations
		if (_debug)
		{
			registration.irtkImageRegistration::Write((char *) "parout-volume.rreg");
			sprintf(buffer, "stack-transformation%i.dof.gz", i);
			stack_transformations[i].irtkTransformation::Write(buffer);
			sprintf(buffer, "stack%i.nii.gz", i);
			stacks[i].Write(buffer);
		}
	}

      InvertStackTransformations(stack_transformations);
}

void irtkReconstructionb0::SetT2Template(irtkRealImage T2)
{
  irtkRealImage t2template = _reconstructed;
  t2template=0;
  irtkRigidTransformation tr;
  
  irtkImageTransformation *imagetransformation = new irtkImageTransformation;
  irtkImageFunction *interpolator = new irtkLinearInterpolateImageFunction;
  imagetransformation->SetInput(&T2, &tr);
  imagetransformation->SetOutput(&t2template);
  //target contains zeros, need padding -1
  imagetransformation->PutTargetPaddingValue(-1);
  //need to fill voxels in target where there is no info from source with zeroes
  imagetransformation->PutSourcePaddingValue(0);
  imagetransformation->PutInterpolator(interpolator);
  imagetransformation->Run();
  
  _reconstructed = t2template;
  _reconstructed.Write("t2template.nii.gz");

}

irtkRealImage irtkReconstructionb0::AlignT2Template(irtkRealImage T2, double sigma)
{
  irtkImageRigidRegistrationWithPadding registration;
  irtkRigidTransformation offset,tr;
  irtkGreyImage target = _reconstructed;
  ResetOrigin(target,offset);

  if (sigma>0)
  {
    irtkGaussianBlurringWithPadding<irtkRealPixel> gb(sigma,0);
    gb.SetInput(&T2);
    gb.SetOutput(&T2);
    gb.Run();
  }
  irtkGreyImage source = T2;
	
  //include offset in trasformation	
  irtkMatrix mo = offset.GetMatrix();
  irtkMatrix m = tr.GetMatrix();
  m=m*mo;
  tr.PutMatrix(m);
  registration.SetInput(&target, &source);
  registration.SetOutput(&tr);
  registration.GuessParameter();
  registration.SetTargetPadding(0);
  registration.Run();
  //undo the offset
  mo.Invert();
  m = tr.GetMatrix();
  m=m*mo;
  tr.PutMatrix(m);
  
  if (_debug)
    tr.irtkTransformation::Write("T2-to-b0.dof");
  
  //transform T2
  irtkRealImage t2template = _reconstructed;
  t2template=0;
  
  irtkImageTransformation *imagetransformation = new irtkImageTransformation;
  irtkImageFunction *interpolator = new irtkLinearInterpolateImageFunction;
  imagetransformation->SetInput(&T2, &tr);
  imagetransformation->SetOutput(&t2template);
  //target contains zeros, need padding -1
  imagetransformation->PutTargetPaddingValue(-1);
  //need to fill voxels in target where there is no info from source with zeroes
  imagetransformation->PutSourcePaddingValue(0);
  imagetransformation->PutInterpolator(interpolator);
  imagetransformation->Run();
  //t2template.Write("at2.nii.gz");
  
  irtkRealImage mask = T2, alignedmask = _reconstructed;
  irtkRealPixel *pm,*pt;
  pm = mask.GetPointerToVoxels();
  for (int i=0;i<mask.GetNumberOfVoxels();i++)
  {
    if(*pm>0) 
      *pm=1;
    pm++;  
  }
  //mask.Write("T2mask.nii.gz");
  imagetransformation->SetInput(&mask, &tr);
  imagetransformation->SetOutput(&alignedmask);
  imagetransformation->Run();
  //alignedmask.Write("alignedT2mask.nii.gz");
  pm = alignedmask.GetPointerToVoxels();
  pt = t2template.GetPointerToVoxels();
  for (int i=0;i<alignedmask.GetNumberOfVoxels();i++)
  {
    if(*pm>=0.5)
      *pt/=*pm;
    else       
      *pt=0;
    pm++;  
    pt++;  
  }
  //t2template.Write("at2masked.nii.gz");
  
  
  return t2template;

}


irtkRealImage irtkReconstructionb0::AdjustOrientation(irtkRealImage &image, bool swap)
{
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //change orientation of the image so that we can use only sy, syx and syz affine parameters   //
  //as given by image coordinate system. This means extracting rotation (and reflection)        //
  //from image header and swapping x and y axis                                                 //
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  //image attributes
  irtkImageAttributes attr = image.GetImageAttributes();
    
  //set orientation in preparation for distortion correction
  
  if (swap)
  {
    attr._xaxis[0] = 0;
    attr._xaxis[1] = 1;
    attr._yaxis[0] = 1;
    attr._yaxis[1] = 0;
  }
  else
  {
    attr._xaxis[0] = 1;
    attr._xaxis[1] = 0;
    attr._yaxis[0] = 0;
    attr._yaxis[1] = 1;    
  }
  
  attr._xaxis[2] = 0;
  attr._yaxis[2] = 0;
  attr._zaxis[0] = 0;
  attr._zaxis[1] = 0;
  attr._zaxis[2] = 1;
    
  //reset origin
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  
  //create image with identity orientation and zero origin
  irtkRealImage newimage(attr);
    
  for (int i=0; i<attr._x; i++)
    for (int j=0; j<attr._y; j++)
      for (int k=0; k<attr._z; k++)
        for (int l=0; l<attr._t; l++)
          newimage(i,j,k,l)=image(i,j,k,l);
  
  return newimage;
}

irtkAffineTransformation irtkReconstructionb0::AdjustOrientationTransformation(irtkRealImage &image, bool swap)
{
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //give transformation of the image with changed orientation (for purpose of distortion        //
  //correction, so that we can use only sy, syx and syz affine parameters) to the original      //
  //image                                                                                       //
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  //image attributes
  irtkImageAttributes attr = image.GetImageAttributes();
  
  irtkMatrix orient(4,4);
  orient.Ident();
  
  //find orientation axis
  //transposing equals inverting because it is rotation matrix
  orient(0, 0) = attr._xaxis[0];
  orient(0, 1) = attr._xaxis[1];
  orient(0, 2) = attr._xaxis[2];
  orient(1, 0) = attr._yaxis[0];
  orient(1, 1) = attr._yaxis[1];
  orient(1, 2) = attr._yaxis[2];
  orient(2, 0) = attr._zaxis[0];
  orient(2, 1) = attr._zaxis[1];
  orient(2, 2) = attr._zaxis[2];
  
  //offset vector
  irtkVector offset(4);
  offset(0)= attr._xorigin;
  offset(1)= attr._yorigin;
  offset(2)= attr._zorigin;
  offset(3)= 1;

  //needs to be rotated
  offset = orient*offset;

  //adjust translation 
  orient(0,3)=-offset(0);
  orient(1,3)=-offset(1);
  orient(2,3)=-offset(2);

  if (swap)
  {
    //matrix to swap x and y axis
    irtkMatrix axis_swap(4,4);
    axis_swap.Ident();
    axis_swap(0,0)=0;
    axis_swap(1,1)=0;
    axis_swap(0,1)=1;
    axis_swap(1,0)=1;
  
    //multiply the matrices
    orient = axis_swap*orient;
  }
  
  //create affine transformation
  irtkAffineTransformation tr;
  tr.PutMatrix(orient);
  //tr.irtkTransformation::Write("adjusted-orient.dof");
  return tr;

}


void irtkReconstructionb0::ShimDistortion(irtkRealImage &acquired, irtkRealImage &simulated, irtkAffineTransformation &shim, bool swap)
{ 
  
  irtkRealImage acq,sim;
  acq = AdjustOrientation(acquired,swap);
  sim = AdjustOrientation(simulated,swap);
  
  irtkAffineTransformation orient = AdjustOrientationTransformation(acquired,swap);
  //orient.irtkTransformation::Write("orient.dof");
  //acq.Write("adjusted.nii.gz");
  //sim.Write("simulated.nii.gz");
  
  //constrain distortion transformation
  irtkAffineTransformation distortion;
  //distortion.PutStatus(TX,  _Passive);
  distortion.PutStatus(TY,  _Passive);
  distortion.PutStatus(TZ,  _Passive);
  distortion.PutStatus(RX,  _Passive);
  distortion.PutStatus(RY,  _Passive);
  distortion.PutStatus(RZ,  _Passive);
  distortion.PutStatus(SY,  _Passive);
  distortion.PutStatus(SZ,  _Passive);
  distortion.PutStatus(SYZ, _Passive);

  double xsize,ysize,zsize;
  _reconstructed.GetPixelSize(&xsize, &ysize, &zsize);
  irtkImageAffineRegistrationWithPadding registration;
  irtkImageTransformation imagetransformation; 
  irtkLinearInterpolateImageFunction interpolator;
  irtkGreyImage t,s;
  t=acq;
  s=sim;
  registration.SetInput(&t,&s);
  registration.SetOutput(&distortion);
  registration.GuessParameterDistortion(xsize);
  registration.SetTargetPadding(0);
  registration.Run();
  distortion.irtkTransformation::Write("d.dof");
  distortion.Print();
  if(_debug)
    registration.Write("par-shim.areg");
  
  irtkMatrix mo = orient.GetMatrix();
  irtkMatrix md = distortion.GetMatrix();
  md = md*mo;
  mo.Invert();
  md = mo*md;
  distortion.PutMatrix(md);
  shim.PutMatrix(md);

  //distortion.irtkTransformation::Write("shim.dof");
  
  //return distortion;
  //irtkMatrix ms = _shim.GetMatrix();
  //ms=ms*md;
  //_shim.PutMatrix(ms);
  //_shim.irtkTransformation::Write("shim.dof");

}

void irtkReconstructionb0::Shim(vector<irtkRealImage> &stacks, int iter)
{
  cout<<"Shimming."<<endl;
  cout.flush();
  if(stacks.size()==0)
  {
    cerr<<"irtkReconstructionb0: Please set the stacks!"<<endl;
    exit(1);
  }
  vector<irtkRealImage> simulated;
  vector<irtkRealImage> stacks2;
  
  unsigned int ind;
  int i,j,k;
  char buffer[256];
  irtkRealImage image;
  double valst,valsim;
  
  //simulate stacks
  for(ind = 0; ind<stacks.size();ind++)
  {
    simulated.push_back(stacks[ind]);
    stacks2.push_back(stacks[ind]);
  }
  SimulateStacks(simulated);
  if(!_have_larger_mask)
    CreateStackMask(simulated);
  if(_debug)
  {
    for(ind = 0; ind<stacks.size();ind++)
    {
      sprintf(buffer,"st%i.nii.gz",ind);
      stacks[ind].Write(buffer);
      sprintf(buffer,"sim%i.nii.gz",ind);
      simulated[ind].Write(buffer);
    }
  }
  
  
  irtkImageTransformation *imagetransformation = new irtkImageTransformation;
  irtkImageFunction *interpolator = new irtkNearestNeighborInterpolateImageFunction;
  irtkImageFunction *interpolatorLin = new irtkLinearInterpolateImageFunction;
  irtkImageFunction *interpolatorSinc = new irtkSincInterpolateImageFunction;
  imagetransformation->PutInterpolator(interpolator);
  //target probably has padding 0, need padding -1
  imagetransformation->PutTargetPaddingValue(-1);
  //need to fill voxels in target where there is no info from source with zeroes
  imagetransformation->PutSourcePaddingValue(0);

  cout<<"done."<<endl;
  cout.flush();
  cout<<"Groups: ";
  cout.flush();
  _shim.clear();
  for(unsigned int g=0; g<_groups.size();g++)
  {
    cout<<g<<" ";
    cout.flush();
    //find all the stacks in the group g
    vector<int> stacknum;
    for(ind = 0; ind<stacks.size();ind++)
      if(_stack_group[ind]==_groups[g])
	stacknum.push_back(ind);

    irtkRealImage templateStack = stacks[stacknum[0]];
    irtkImageAttributes attr = templateStack.GetImageAttributes();
    attr._t = stacknum.size();
    cout<<endl<<"We have "<<attr._t<<" images in group "<<_groups[g]<<endl;
    irtkRealImage stack(attr), simul(attr);
  
    irtkRigidTransformation id;
    
    irtkRealImage st(templateStack),sim(templateStack);
    //resample all on the same grid
    for(ind=0; ind<stacknum.size(); ind++)
    {
      imagetransformation->SetInput(&stacks[stacknum[ind]], &id);
      imagetransformation->SetOutput(&st);
      imagetransformation->Run();
    
      imagetransformation->SetInput(&simulated[stacknum[ind]], &id);
      imagetransformation->SetOutput(&sim);
      imagetransformation->Run();
    
      for(i=0;i<templateStack.GetX();i++)
        for(j=0;j<templateStack.GetY();j++)
          for(k=0;k<templateStack.GetZ();k++)
	  {
	    valst = st(i,j,k);
	    valsim = sim(i,j,k);
	    if(valst>0.01)
	      stack(i,j,k,ind)=valst;
	    else
	      stack(i,j,k,ind)=0;
	    if(valsim>0.01)
	      simul(i,j,k,ind)=valsim; 
	    else
	      simul(i,j,k,ind)=0;
	    /*
	    if ((valst>0)&&(valsim>0))
	    {
	      stack(i,j,k,ind)=valst;
	      simul(i,j,k,ind)=valsim; 
	    }
	    else
	    {
	      stack(i,j,k,ind)=0;
	      simul(i,j,k,ind)=0;
	    }
	    */
	  }
    }
    
    if(_debug)
    {
      sprintf(buffer,"stacks%i-%i.nii.gz",iter,g);
      stack.Write(buffer);
      sprintf(buffer,"sims%i-%i.nii.gz",iter,g);
      simul.Write(buffer);
    }
    
    //calculate shim
    irtkAffineTransformation shim;
    //ShimDistortion(stack,simul,shim,false);
    
    //if(g==0)
      ShimDistortion(stack,simul,shim,_swap[g]);
    //else
      //ShimDistortion(stack,simul,shim,true);
    //sprintf(buffer,"ishim%i.nii.gz",g);
    //shim.irtkTransformation::Write(buffer);
    shim.Invert();
    shim.UpdateParameter();
    shim.Print();
    if(_debug)
    {
      sprintf(buffer,"shim%i-%i.dof",iter,g);
      shim.irtkTransformation::Write(buffer);
    }
    _shim.push_back(shim);
    
    //TODO: correct sinc artefacts - need padding
     imagetransformation->PutInterpolator(interpolatorSinc);
    
    for(ind=0; ind<stacknum.size(); ind++)
    {
      //sprintf(buffer,"s1-%i.nii.gz",stacknum[ind]);
      //stacks[stacknum[ind]].Write(buffer);
      //sprintf(buffer,"s2-%i.nii.gz",stacknum[ind]);
      //stacks2[stacknum[ind]].Write(buffer);
      cout<<"Correcting stack "<<stacknum[ind]<<endl;
      imagetransformation->SetInput(&stacks2[stacknum[ind]], &shim);
      imagetransformation->SetOutput(&stacks[stacknum[ind]]);
      imagetransformation->Run();
      
      //sprintf(buffer,"ss1-%i.nii.gz",stacknum[ind]);
      //stacks[stacknum[ind]].Write(buffer);
      //sprintf(buffer,"ss2-%i.nii.gz",stacknum[ind]);
      //stacks2[stacknum[ind]].Write(buffer);
    }
  }
  delete imagetransformation;
  delete interpolator;
  delete interpolatorLin;
}

void irtkReconstructionb0::FieldMap(vector<irtkRealImage> &stacks, double step, int iter)
{
  cout<<"Field map correction."<<endl;
  cout<<"FieldMap correction: Warning: only implemented for stacks of same geometry at present!."<<endl;
  cout.flush();
  if(stacks.size()==0)
  {
    cerr<<"irtkReconstructionb0: Please set the stacks!"<<endl;
    exit(1);
  }
  vector<irtkRealImage> simulated;
  vector<irtkRealImage> stacks2;
  
  unsigned int ind;
  int i,j,k;
  char buffer[256];
  irtkRealImage image;
  double valst,valsim;
  
  //simulate stacks
  for(ind = 0; ind<stacks.size();ind++)
  {
    simulated.push_back(stacks[ind]);
    stacks2.push_back(stacks[ind]);
  }
  SimulateStacks(simulated);
  if(_debug)
  {
    for(ind = 0; ind<stacks.size();ind++)
    {
      sprintf(buffer,"st%i.nii.gz",ind);
      stacks[ind].Write(buffer);
      sprintf(buffer,"sim%i.nii.gz",ind);
      simulated[ind].Write(buffer);
    }
  }
  
  //Resample stacks and simulated to the same geometry and create 4D nifti
  irtkImageTransformation *imagetransformation = new irtkImageTransformation;
  irtkImageFunction *interpolator = new irtkNearestNeighborInterpolateImageFunction;
  irtkImageFunction *interpolatorLin = new irtkLinearInterpolateImageFunction;
  imagetransformation->PutInterpolator(interpolator);
  //target probably has padding 0, need padding -1
  imagetransformation->PutTargetPaddingValue(-1);
  //need to fill voxels in target where there is no info from source with zeroes
  imagetransformation->PutSourcePaddingValue(0);


  irtkRealImage templateStack = stacks[0];
  irtkImageAttributes attr = templateStack.GetImageAttributes();
  attr._t = stacks.size();
  cout<<endl<<"We have "<<attr._t<<" images "<<endl;
  irtkRealImage stack(attr), simul(attr);
  
  irtkRigidTransformation id;
    
  irtkRealImage st(templateStack),sim(templateStack);
  //resample all on the same grid
  for(ind=0; ind<stacks.size(); ind++)
  {
    imagetransformation->SetInput(&stacks[ind], &id);
    imagetransformation->SetOutput(&st);
    imagetransformation->Run();
    
    imagetransformation->SetInput(&simulated[ind], &id);
    imagetransformation->SetOutput(&sim);
    imagetransformation->Run();
    
    for(i=0;i<templateStack.GetX();i++)
      for(j=0;j<templateStack.GetY();j++)
        for(k=0;k<templateStack.GetZ();k++)
	{
	  valst = st(i,j,k);
	  valsim = sim(i,j,k);
	  if(valst>0.01)
	    stack(i,j,k,ind)=valst;
	  else
	    stack(i,j,k,ind)=0;
	  if(valsim>0.01)
	    simul(i,j,k,ind)=valsim; 
	  else
	    simul(i,j,k,ind)=0;
	}
  }
  
  if (_debug)
  {
    sprintf(buffer,"fmstacks%i.nii.gz",iter);
    stack.Write(buffer);
    sprintf(buffer,"fmsims%i.nii.gz",iter);
    simul.Write(buffer);
  }
    
  //calculate b0 field distortiom
  //irtkMultiLevelFreeFormTransformation dist;
  //FieldMapDistortion(stack,simul,dist,_swap[0]);
  
  //clear fieldmap
  if(_fieldMap.NumberOfLevels()>0)
  {
    irtkFreeFormTransformation *tr = _fieldMap.PopLocalTransformation();
    delete tr;
  }
  _fieldMap.irtkTransformation::Write("fieldMap_clean.dof");
  cout<<"entering FieldMapDistortion."<<endl;
  cout<<_swap.size()<<endl;
  cout.flush();
  FieldMapDistortion(stack,simul,_fieldMap,_swap[0],step,0,3);
  
  if(_debug)
  {
    sprintf(buffer,"fmdist%i.dof",iter);
    _fieldMap.irtkTransformation::Write(buffer);
  }
    
  //Corect the stacks
  imagetransformation->PutInterpolator(interpolatorLin); 
  for(ind=0; ind<stacks.size(); ind++)
  {
    cout<<"Correcting stack "<<ind<<endl;
    imagetransformation->SetInput(&stacks2[ind], &_fieldMap);//&dist);
    imagetransformation->SetOutput(&stacks[ind]);
    imagetransformation->Run();
  }
  
  delete imagetransformation;
  delete interpolator;
  delete interpolatorLin;
}

void irtkReconstructionb0::SimulateStacksWithMask(vector<irtkRealImage>& stacks, irtkRealImage mask)
{
    if (_debug)
        cout<<"Simulating stacks."<<endl;
  
    unsigned int inputIndex;
    int i, j, k, n;
    irtkRealImage sim;
    POINT3D p;
    double weight;
    double xx,yy,zz;
  
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
	    {
	      //check whether the slice voxel is in ROI defined by mask
	      xx=i;yy=j;zz=0;
	      slice.ImageToWorld(xx,yy,zz);
	      mask.WorldToImage(xx,yy,zz);
	      xx=round(xx);yy=round(yy);zz=round(zz);
	      if((xx>=0)&&(xx<mask.GetX())&&(yy>=0)&&(yy<mask.GetY())&&(zz>=0)&&(zz<mask.GetZ()))
		if(mask(xx,yy,zz)>0)
		{   
		  weight=0;
                  n = _volcoeffs[inputIndex][i][j].size();
                  for (k = 0; k < n; k++) {
                    p = _volcoeffs[inputIndex][i][j][k];
                    sim(i, j, 0) += p.value * _reconstructed(p.x, p.y, p.z);
                    weight += p.value;
                  }
                  if(weight>0.99)
                    sim(i,j,0)/=weight;
		   else
		     sim(i,j,0)=0;
		}
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

void irtkReconstructionb0::FieldMapGroup(vector<irtkRealImage> &stacks, irtkRealImage stackMask, int group, double step, int iter)
{
  cout<<"Field map correction."<<endl;
  cout<<"FieldMap correction: Warning: only implemented for stacks of same geometry at present!."<<endl;
  cout.flush();
  if(stacks.size()==0)
  {
    cerr<<"irtkReconstructionb0: Please set the stacks!"<<endl;
    exit(1);
  }
  vector<irtkRealImage> simulated;
  vector<irtkRealImage> simulatedgroup;
  vector<irtkRealImage> stacksgroup;
  vector<irtkRealImage> stacks2;
  
  unsigned int ind;
  int i,j,k;
  char buffer[256];
  irtkRealImage image;
  double valst,valsim;
  
  //simulate stacks
  for(ind = 0; ind<stacks.size();ind++)
  {
    simulated.push_back(stacks[ind]);
    stacks2.push_back(stacks[ind]);
  }
  
  cout<<"Group = "<<group<<endl;
  //if(group==0)
    SimulateStacksWithMask(simulated,stackMask);
  //else
    //SimulateStacksWithMask(simulated,_mask);
  
  cout<<"group: "<<group<<" _groups[group]"<<_groups[group]<<endl;
  cout<<"_stack_group: ";
  //choose stacks for current group 
  for(ind = 0; ind<stacks.size();ind++)
  {
    if(_stack_group[ind]==_groups[group])
    {
      cout<<ind<<":"<<_stack_group[ind]<<" ";
      simulatedgroup.push_back(simulated[ind]);
      stacksgroup.push_back(stacks[ind]);
    }
  }
  
  if(_debug)
  {
    for(ind = 0; ind<stacksgroup.size();ind++)
    {
      sprintf(buffer,"st%i.nii.gz",ind);
      stacksgroup[ind].Write(buffer);
      sprintf(buffer,"sim%i.nii.gz",ind);
      simulatedgroup[ind].Write(buffer);
    }
  }
  
  //Resample stacks and simulated to the same geometry and create 4D nifti
  irtkImageTransformation *imagetransformation = new irtkImageTransformation;
  irtkImageFunction *interpolator = new irtkNearestNeighborInterpolateImageFunction;
  irtkImageFunction *interpolatorLin = new irtkLinearInterpolateImageFunction;
  imagetransformation->PutInterpolator(interpolator);
  //target probably has padding 0, need padding -1
  imagetransformation->PutTargetPaddingValue(-1);
  //need to fill voxels in target where there is no info from source with zeroes
  imagetransformation->PutSourcePaddingValue(0);


  irtkRealImage templateStack = stacksgroup[0];
  irtkImageAttributes attr = templateStack.GetImageAttributes();
  attr._t = stacksgroup.size();
  cout<<endl<<"We have "<<attr._t<<" images "<<endl;
  irtkRealImage stack(attr), simul(attr);
  
  irtkRigidTransformation id;
    
  irtkRealImage st(templateStack),sim(templateStack);
  //resample all on the same grid
  for(ind=0; ind<stacksgroup.size(); ind++)
  {
    imagetransformation->SetInput(&stacksgroup[ind], &id);
    imagetransformation->SetOutput(&st);
    imagetransformation->Run();
    
    imagetransformation->SetInput(&simulatedgroup[ind], &id);
    imagetransformation->SetOutput(&sim);
    imagetransformation->Run();
    
    for(i=0;i<templateStack.GetX();i++)
      for(j=0;j<templateStack.GetY();j++)
        for(k=0;k<templateStack.GetZ();k++)
	{
	  valst = st(i,j,k);
	  valsim = sim(i,j,k);
	  if(valst>0.01)
	    stack(i,j,k,ind)=valst;
	  else
	    stack(i,j,k,ind)=0;
	  if(valsim>0.01)
	    simul(i,j,k,ind)=valsim; 
	  else
	    simul(i,j,k,ind)=0;
	}
  }
  
  if (_debug)
  {
    sprintf(buffer,"fmstacks%i-%i.nii.gz",iter,group);
    stack.Write(buffer);
    sprintf(buffer,"fmsims%i-%i.nii.gz",iter,group);
    simul.Write(buffer);
  }
    
  //calculate b0 field distortiom
  //irtkMultiLevelFreeFormTransformation dist;
  //FieldMapDistortion(stack,simul,dist,_swap[0]);
  
  //clear fieldmap
  if(_fieldMap.NumberOfLevels()>0)
  {
    irtkFreeFormTransformation *tr = _fieldMap.PopLocalTransformation();
    delete tr;
  }
  _fieldMap.irtkTransformation::Write("fieldMap_clean.dof");
  FieldMapDistortion(stack,simul,_fieldMap,_swap[group],step,iter);
  
  if(_debug)
  {
    sprintf(buffer,"fmdist%i-%i.dof",iter,group);
    _fieldMap.irtkTransformation::Write(buffer);
  }
    /*
  //Corect the stacks
  imagetransformation->PutInterpolator(interpolatorLin); 
  for(ind=0; ind<stacks.size(); ind++)
  {
    cout<<"Correcting stack "<<ind<<endl;
    imagetransformation->SetInput(&stacks2[ind], &_fieldMap);//&dist);
    imagetransformation->SetOutput(&stacks[ind]);
    imagetransformation->Run();
  }
  */
  //resample distortion on a template
  irtkRealImage tempDist = templateStack;
  _distortion = _reconstructed;
  _distortion=0;
  tempDist = 0;
  irtkMultiLevelFreeFormTransformation dist(_fieldMap);
  irtkMatrix m;
  double x,y,z,xx,yy,zz;
  attr = tempDist.GetImageAttributes();

  //TODO: add shim
  //Calculate distortion displacements from the transformation
 // m = _shim[group].GetMatrix();
  //dist.PutMatrix(m);

  tempDist=0;
  for(int k=0; k<_distortion.GetZ();k++)
  for(int j=0; j<_distortion.GetY();j++)
  for(int i=0; i<_distortion.GetX();i++)
  {
    //transforming to the coordinates of the stack group
    //to get displacement in PE direction
    
    //origin
    xx=i;yy=j;zz=k;
    _distortion.ImageToWorld(xx,yy,zz);
    tempDist.WorldToImage(xx,yy,zz);

    //transformed
    x=i;y=j;z=k;
    _distortion.ImageToWorld(x,y,z);
    dist.Transform(x,y,z);
    tempDist.WorldToImage(x,y,z);
    
    if (_swap[group])
      _distortion(i,j,k)=(y-yy)*attr._dy;
    else
      _distortion(i,j,k)=(x-xx)*attr._dx;
   }

    sprintf(buffer,"fmdist%i-%i.nii.gz",iter,group);
    _distortion.Write(buffer);
  
  delete imagetransformation;
  delete interpolator;
  delete interpolatorLin;
}


void  irtkReconstructionb0::FieldMapDistortion(irtkRealImage &stack,irtkRealImage &simul, irtkMultiLevelFreeFormTransformation &distortion, bool swap, double step, int iter, int levels)
{
  cout<<"Entered FieldMapDistortion."<<endl;
  cout.flush();
  //Adjust orientation
  irtkRealImage st,sm;
  sm = AdjustOrientation(simul,false);
  st = AdjustOrientation(stack,false);
  if (_debug)
  {
    sm.Write("sm.nii.gz");
    st.Write("st.nii.gz");
  }
  irtkAffineTransformation orient = AdjustOrientationTransformation(simul,false);
  //orient.irtkTransformation::Write("orient.dof");
  
  //register acquired stacks to simulated
  irtkImageFreeFormRegistrationWithPadding registration;
  if(swap)
    registration.SetMode(RegisterY);
  else
    registration.SetMode(RegisterX);
  
  irtkGreyImage t,s;
  t=sm;
  s=st;
  registration.SetInput(&t,&s);
  registration.SetOutput(&distortion);
  //int levels;
  double res;
  
  /*
  if(iter<=4)
  {
    levels=3;
    res = 1.5;
  }
  else
  {
    levels=4;
    res = 0.75;
  }
  */
  //levels=1;
  irtkImageAttributes attrt = t.GetImageAttributes();
  res = attrt._dx;
  
  //registration.GuessParameterDistortion(res,_fieldMapSpacing*8,levels,step,_smoothnessPenalty);
  //registration.GuessParameterDistortion(res,_fieldMapSpacing*2,levels,step,_smoothnessPenalty);
  registration.GuessParameterDistortion(res,_fieldMapSpacing*8,levels,step*2,_smoothnessPenalty);
  if(_debug)
  {
    char buffer[256];
    sprintf(buffer,"par-dist-%f.nreg",step*4);
    registration.irtkImageRegistration::Write(buffer);
  }
  registration.SetTargetPadding(0);
  //if (_smoothnessPenalty>0)
    registration.RunRelax();
  //else
    //registration.Run();
  //distortion.irtkTransformation::Write("fmd.dof");
  
  //adjust lattice of Bspine transformation according to the original images
  irtkImageAttributes attr = simul.GetImageAttributes();
  irtkFreeFormTransformation3D *bspline = dynamic_cast<irtkFreeFormTransformation3D *>(distortion.GetLocalTransformation(0));
  bspline->PutOrigin(attr._xorigin, attr._yorigin, attr._zorigin);
  bspline->PutOrientation(attr._xaxis, attr._yaxis, attr._zaxis);
  
  //orient Bspline control points
  irtkMatrix mo = orient.GetMatrix();
  mo.Invert();

  irtkVector cp(4);
  cp(3)=0;
  for (int k = 0; k < bspline->GetZ(); ++k) {
    for (int j = 0; j < bspline->GetY(); ++j) {
      for (int i = 0; i < bspline->GetX(); ++i) {
        bspline->Get(i, j, k, cp(0), cp(1), cp(2));
	cp=mo*cp;
        bspline->Put(i, j, k, cp(0), cp(1), cp(2));
      }
    }  
  }
  //distortion.irtkTransformation::Write("fmdist.dof");
}


irtkRealImage irtkReconstructionb0::Create4DImage(vector<irtkRealImage> &stacks)
{
  int i,j,k;
  double val;
  irtkImageTransformation *imagetransformation = new irtkImageTransformation;
  irtkImageFunction *interpolator = new irtkNearestNeighborInterpolateImageFunction;
  //irtkImageFunction *interpolatorLin = new irtkLinearInterpolateImageFunction;
  imagetransformation->PutInterpolator(interpolator);
  //target probably has padding 0, need padding -1
  imagetransformation->PutTargetPaddingValue(-1);
  //need to fill voxels in target where there is no info from source with zeroes
  imagetransformation->PutSourcePaddingValue(0);
  
  irtkRealImage templateStack = stacks[0];
  irtkImageAttributes attr = templateStack.GetImageAttributes();
  attr._t = stacks.size();
  irtkRealImage stack(attr);
  
  irtkRigidTransformation id;
  irtkRealImage st(templateStack);
  //resample all on the same grid
  for(unsigned int ind=0; ind<stacks.size(); ind++)
  {
    imagetransformation->SetInput(&stacks[ind], &id);
    imagetransformation->SetOutput(&st);
    imagetransformation->Run();
    
    for(i=0;i<templateStack.GetX();i++)
      for(j=0;j<templateStack.GetY();j++)
        for(k=0;k<templateStack.GetZ();k++)
        {
	  val = st(i,j,k);
	  if (val>0)
	    stack(i,j,k,ind)=val;
	  else
	  stack(i,j,k,ind)=0;
	}
    }
    //stack.Write("4D.nii.gz");
    return stack;
  
}

void irtkReconstructionb0::CreateSimulated(vector<irtkRealImage> &stacks)
{
  _simulated.clear();
  for (unsigned int i=0;i<stacks.size();i++)
  {
    _simulated.push_back(stacks[i]);
    _simulated[i]=0;
  }
  SimulateStacks(_simulated);
}

void irtkReconstructionb0::CreateStackMask(vector<irtkRealImage> &simulated)
{
  _larger_mask = _reconstructed;
  _larger_mask = 0;
  _have_larger_mask=true;
  
  /*
  irtkRealPixel *pi,*pm;
  
  for(uint t=0; t<simulated.size();t++)
  {
    pi = simulated[t].GetPointerToVoxels();
    pm = _larger_mask.GetPointerToVoxels();
    for(int i=0;i<_larger_mask.GetNumberOfVoxels();i++)
    {
      if(*pi>0)
	*pm=1;
      pi++;
      pm++;
    }
  }
  */
 
  double x,y,z,value;
  char buffer[256];
  irtkRealImage testImage =_larger_mask;
  testImage=0;
  irtkLinearInterpolateImageFunction interpolator;
  cout<<"Creating stack mask ... ";
  cout.flush();
  for(uint t=0; t<simulated.size();t++)
  {
    cout<<t<<" ";
    cout.flush();
    interpolator.SetInput(&(simulated[t]));
    interpolator.Initialize();
    testImage=0;
    
    for(int i=0;i<_larger_mask.GetX();i++)
      for(int j=0;j<_larger_mask.GetY();j++)
        for(int k=0;k<_larger_mask.GetZ();k++)
        {
	  x=i;y=j;z=k;
	  _larger_mask.ImageToWorld(x,y,z);
	  simulated[t].WorldToImage(x,y,z);
	  value = interpolator.Evaluate(x,y,z);
	  testImage(i,j,k)=value;
	  if(value>0)
	    _larger_mask(i,j,k)=1;
        }
        
     sprintf(buffer,"testImage%i.nii.gz",t);
    testImage.Write(buffer);
  }
  cout<<"done."<<endl;
  cout.flush();
  _larger_mask.Write("stackmask.nii.gz");
  CreateLargerMask(_larger_mask);
  _larger_mask.Write("stackmask-larger.nii.gz");
  //exit(1);
}

void irtkReconstructionb0::WriteSimulated()
{
  char buffer[256];
  for(unsigned int i=0;i<_simulated.size();i++)
  {
    sprintf(buffer,"simulated%i.nii.gz",i);
    _simulated[i].Write(buffer);
  }
}

void irtkReconstructionb0::SaveDistortionTransformations()
{
  char buffer[256];
  irtkMultiLevelFreeFormTransformation dist(_fieldMap);
  irtkMatrix m;
  irtkRealImage distortion = _reconstructed;
  double x,y,z;
  irtkImageAttributes attr = distortion.GetImageAttributes();
  double res = attr._dx;

  for(unsigned int ind=0;ind<_shim.size();ind++)
  {
    m = _shim[ind].GetMatrix();
    dist.PutMatrix(m);

    distortion=0;
    for(int k=0; k<distortion.GetZ();k++)
    for(int j=0; j<distortion.GetY();j++)
    for(int i=0; i<distortion.GetX();i++)
    {
      x=i;y=j;z=k;
      distortion.ImageToWorld(x,y,z);
      dist.Transform(x,y,z);
      distortion.WorldToImage(x,y,z);
      if (_swap[ind])
	distortion(i,j,k)=(y-j)*res;
      else
	distortion(i,j,k)=(x-i)*res;
    }
      
    if(_shim.size()>1)
    {
      sprintf(buffer,"distortion%i.dof",ind);
      dist.irtkTransformation::Write(buffer);
      sprintf(buffer,"distortion%i.nii.gz",ind);
      distortion.Write(buffer);
    }
    else
    {
      dist.irtkTransformation::Write("distortion.dof");
      distortion.Write("distortion.nii.gz");
    }
  }
  
  
}

irtkRealImage irtkReconstructionb0::CreateLargerMask(irtkRealImage mask)
{
    _larger_mask = mask;
    irtkImageAttributes attr = _mask.GetImageAttributes();
    irtkGaussianBlurring<irtkRealPixel> gb(attr._dx*2);
    gb.SetInput(&_larger_mask);
    gb.SetOutput(&_larger_mask);
    gb.Run();
    //binarize mask
    irtkRealPixel* ptr = _larger_mask.GetPointerToVoxels();
    for (int i = 0; i < _larger_mask.GetNumberOfVoxels(); i++) {
      if (*ptr > 0.05)
	*ptr = 1;
      else
	*ptr = 0;
      ptr++;
    }
    _larger_mask.Write("larger_mask.nii.gz");
    return _larger_mask;
 
}

void irtkReconstructionb0::SmoothFieldmap(int iter)
{
  
  char buffer[256];
  irtkMultiLevelFreeFormTransformation dist(_fieldMap);
  irtkMatrix m;
  irtkRealImage distortion = _reconstructed;
  double x,y,z;
  irtkImageAttributes attr = distortion.GetImageAttributes();
  double res = attr._dx;
  //change 3: do not clear fieldmap, but add to it
  //_smoothFieldMap.clear();

  bool first_iter = false;
  if (_smoothFieldMap.size()==0)
    first_iter=true;
  
  for(unsigned int ind=0;ind<_shim.size();ind++)
  {
    //Calculate distortion displacements from the transformation
    m = _shim[ind].GetMatrix();
    dist.PutMatrix(m);

    distortion=0;
    for(int k=0; k<distortion.GetZ();k++)
    for(int j=0; j<distortion.GetY();j++)
    for(int i=0; i<distortion.GetX();i++)
    {
      x=i;y=j;z=k;
      distortion.ImageToWorld(x,y,z);
      dist.Transform(x,y,z);
      distortion.WorldToImage(x,y,z);
      if (_swap[ind])
	distortion(i,j,k)=(y-j)*res;
      else
	distortion(i,j,k)=(x-i)*res;
    }

    //Smooth the distortion displacements
    irtkLaplacianSmoothing smoothing;
    smoothing.SetInput(distortion);
    smoothing.SetMask(_larger_mask);
    irtkRealImage fieldmap = smoothing.RunGD();
    distortion.Write("d.nii.gz");
    fieldmap.Write("f.nii.gz");
    
    //change 4: add new only at first iter otherwise add to existing
    if(first_iter)
      _smoothFieldMap.push_back(fieldmap);
    else
      _smoothFieldMap[ind]+=fieldmap;
    
     
    if(_shim.size()>1)
    {
      sprintf(buffer,"fieldmap%i-%i.nii.gz",iter,ind);
      fieldmap.Write(buffer);
    }
    else
    {
      sprintf(buffer,"incr-fieldmap%i.nii.gz",iter);
      fieldmap.Write(buffer);
      sprintf(buffer,"fieldmap%i.nii.gz",iter);
      _smoothFieldMap[ind].Write(buffer);
      sprintf(buffer,"distortion%i.nii.gz",iter);
      distortion.Write(buffer);
    }
  }
  
  
}


void irtkReconstructionb0::SmoothFieldmapGroup(irtkRealImage mask, int group, int iter, bool combine_fieldmap, double boundary_weight,bool minus)
{
  
  char buffer[256];

  cout<<"SmoothFieldmapGroup"<<endl;
  cout.flush();
  
  irtkRealImage fieldmap = _distortion;
  
  if (_smoothFieldMap.size()==0)
  {
      fieldmap=0;
      for(unsigned int ind=0;ind<_groups.size();ind++)
	_smoothFieldMap.push_back(fieldmap);
  }
  
  //Smooth the distortion displacements
  cout<<"test"<<endl;
  cout.flush();
  irtkLaplacianSmoothing smoothing;
  sprintf(buffer,"_distortion%i-%i.nii.gz",iter,group);
  _distortion.Write(buffer);
  smoothing.SetInput(_distortion);
  sprintf(buffer,"fieldmapmask%i-%i.nii.gz",iter,group);
  mask.Write(buffer);
  smoothing.SetMask(mask);
  _larger_mask = _mask;
  smoothing.SetBoundaryWeight(boundary_weight);
  fieldmap = smoothing.RunGD();
  //_distortion.Write("d.nii.gz");
  //fieldmap.Write("f.nii.gz");
    
  //change: add to both
  if (combine_fieldmap)
  {
     _smoothFieldMap[0]+=fieldmap;
     if(minus)
       _smoothFieldMap[1]-=fieldmap;
     else
       _smoothFieldMap[1]+=fieldmap;
  }
  else
  {
    if((minus)&&(group==1))
     _smoothFieldMap[group]-=fieldmap;
    else
     _smoothFieldMap[group]+=fieldmap;
  }
    
     
    if(_groups.size()>1)
    {
      sprintf(buffer,"fieldmap%i-%i.nii.gz",iter,group);
      fieldmap.Write(buffer);
      for(int g=0;g<_groups.size();g++)
      {
        sprintf(buffer,"addedfieldmap%i-%i.nii.gz",iter,g);
        _smoothFieldMap[g].Write(buffer);
      }
    }
    else
    {
      sprintf(buffer,"incr-fieldmap%i.nii.gz",iter);
      fieldmap.Write(buffer);
      sprintf(buffer,"fieldmap%i.nii.gz",iter);
      _smoothFieldMap[0].Write(buffer);
      sprintf(buffer,"distortion%i.nii.gz",iter);
      _distortion.Write(buffer);
    }
    
}


void irtkReconstructionb0::CorrectStacks(vector<irtkRealImage> &stacks)
{
  char buffer[256];
  irtkMatrix m;
  unsigned int ind, i;
  
  //make a copy of stacks
  vector<irtkRealImage> stacks2;
  for(unsigned int ind = 0; ind<stacks.size();ind++)
  {
    stacks2.push_back(stacks[ind]);
  }

  
  for(i=0;i<_shim.size();i++)
  {
    //compose shim and fieldmap
    irtkMultiLevelFreeFormTransformation dist(_fieldMap);
    m = _shim[i].GetMatrix();
    dist.PutMatrix(m);
    
    //apply distortion to relevant stacks
    irtkImageTransformation *imagetransformation = new irtkImageTransformation;
    irtkImageFunction *interpolatorLin = new irtkLinearInterpolateImageFunction;
    irtkImageFunction *interpolatorSinc = new irtkSincInterpolateImageFunction;
    imagetransformation->PutInterpolator(interpolatorSinc);
    //target probably has padding 0, need padding -1
    imagetransformation->PutTargetPaddingValue(-1);
    //need to fill voxels in target where there is no info from source with zeroes
    imagetransformation->PutSourcePaddingValue(0);
    
    for(ind = 0; ind<stacks.size();ind++)
      if(_stack_group[ind]==_groups[i])
      {
	imagetransformation->SetInput(&stacks2[ind], &dist);
        imagetransformation->SetOutput(&stacks[ind]);
        imagetransformation->Run();
        if (_debug)
        {
          sprintf(buffer,"cor%i.nii.gz",ind);
	  stacks[ind].Write(buffer);
        }
      }

    
/*    if(_debug)
    {
      if(_shim.size()>1)
      {
        sprintf(buffer,"distortion%i.dof",i);
        dist.irtkTransformation::Write(buffer);
      }
      else
        dist.irtkTransformation::Write("distortion.dof");
    }
    */
  }
}

void irtkReconstructionb0::CorrectStacksSmoothFieldmap(vector<irtkRealImage> &stacks)
{
  char buffer[256];
  irtkMatrix m;
  unsigned int index, ind;
  int i,j,k,t;
  irtkRealImage fieldmap, fieldmap2, image;
  int padding;
  
  for(index=0;index<_smoothFieldMap.size();index++)
  {
    cout<<"Index = "<<index<<endl;
    //padding for fieldmap - careful vector = only copies reference
    fieldmap2 = _smoothFieldMap[index];
    fieldmap = fieldmap2;
    irtkRealPixel mn,mx;
    fieldmap.GetMinMax(&mn,&mx);
    padding = round(mn-2);
    cout<<"padding = "<<padding<<endl;
    irtkRealPixel *pf = fieldmap.GetPointerToVoxels();
    irtkRealPixel *pm = _larger_mask.GetPointerToVoxels();
    for (j=0; j<fieldmap.GetNumberOfVoxels();j++)
    {
      if(*pm!=1)
        *pf = padding;
      pm++;
      pf++;
    }

    fieldmap.Write("orig-f-pad.nii.gz");
    
    //prepare fieldmap for interpolation with padding
    irtkLinearInterpolateImageFunction interpolatorFieldmap;
    interpolatorFieldmap.SetInput(&fieldmap);
    interpolatorFieldmap.Initialize();


    //irtkImageFunction *interpolatorSinc = new irtkSincInterpolateImageFunction;
    
    //undistort stack using fieldmap
    double f,x,y,z;
    for(ind = 0; ind<stacks.size();ind++)
    {
      cout<<"ind = "<<ind<<endl;
      cout<<"_stack_group[ind] = "<<_stack_group[ind]<<endl;
      cout<<"_groups[index] = "<<_groups[index]<<endl;
      
      if(_stack_group[ind]==_groups[index])
      {
	cout<<"I am inside!"<<endl;
	image = stacks[ind];
        //irtkLinearInterpolateImageFunction interpolator;
	irtkSincInterpolateImageFunction interpolator;
        interpolator.SetInput(&image);
        interpolator.Initialize();
        
	cout<<"Correcting stack "<<ind<<endl;
	irtkImageAttributes attr = image.GetImageAttributes();
        for (t = 0; t < image.GetT(); t++) {
          for (k = 0; k < image.GetZ(); k++) {
            for (j = 0; j < image.GetY(); j++) {
              for (i = 0; i < image.GetX(); i++) {
		//find fieldmap value f
                x = i;
                y = j;
                z = k;
		image.ImageToWorld(x,y,z);
		fieldmap.WorldToImage(x,y,z);
		if ((x > -0.5) && (x < fieldmap.GetX()-0.5) && 
	            (y > -0.5) && (y < fieldmap.GetY()-0.5) &&
                    (z > -0.5) && (z < fieldmap.GetZ()-0.5) )
	        {
	         f = interpolatorFieldmap.EvaluateWithPadding(x,y,z,t,padding);
	        }
	        else
		  f=padding;
		
		//find displacement
		if (f>padding)  
		{
		  x = i;
		  y = j;
		  z = k;

		  if(_swap[index])
	            y+=f/attr._dy;
	          else
	            x+=f/attr._dx;
		 
		 if ((x > -0.5) && (x < image.GetX()-0.5) && 
	             (y > -0.5) && (y < image.GetY()-0.5) &&
                     (z > -0.5) && (z < image.GetZ()-0.5))
		 {
	             stacks[ind](i,j,k,t) = interpolator.Evaluate(x,y,z,t);
		 }
	         else
		   stacks[ind](i,j,k)=0;
		}
		else
		  stacks[ind](i,j,k)=0;
	      }//i
	    }//j
	  }//k
	}//t

	
        if (_debug)
        {
          sprintf(buffer,"cor%i.nii.gz",ind);
	  stacks[ind].Write(buffer);
        }
      }//if   
    }//ind
  }//index
}//function

void irtkReconstructionb0::CorrectStacksSmoothFieldmapWithMasks(vector<irtkRealImage> &stacks, irtkRealImage fieldmapMask1, irtkRealImage fieldmapMask2, bool minus)
{
  char buffer[256];
  irtkMatrix m;
  unsigned int index, ind;
  int i,j,k,t;
  irtkRealImage fieldmap, fieldmap2, image;
  int padding;
  
  
  fieldmapMask1.Write("fieldmapMask1.nii.gz");
  fieldmapMask2.Write("fieldmapMask2.nii.gz");
  
  for(index=0;index<_smoothFieldMap.size();index++)
  {
    cout<<"Index = "<<index<<endl;
    //padding for fieldmap - careful vector = only copies reference
    fieldmap2 = _smoothFieldMap[index];
    fieldmap = fieldmap2;
    irtkRealPixel mn,mx;
    fieldmap.GetMinMax(&mn,&mx);
    padding = round(mn-2);
    cout<<"padding = "<<padding<<endl;
    irtkRealPixel *pf = fieldmap.GetPointerToVoxels();
    irtkRealPixel *pm;
    if(index==1)
      pm = fieldmapMask1.GetPointerToVoxels();
    else
      pm = fieldmapMask2.GetPointerToVoxels();

    for (j=0; j<fieldmap.GetNumberOfVoxels();j++)
    {
      if(*pm!=1)
        *pf = padding;
      pm++;
      pf++;
    }

    if(index==0)
      fieldmap.Write("orig-f-pad-1.nii.gz");
    else
      fieldmap.Write("orig-f-pad-2.nii.gz");
    
    //prepare fieldmap for interpolation with padding
    irtkLinearInterpolateImageFunction interpolatorFieldmap;
    interpolatorFieldmap.SetInput(&fieldmap);
    interpolatorFieldmap.Initialize();


    //irtkImageFunction *interpolatorSinc = new irtkSincInterpolateImageFunction;
    
    //undistort stack using fieldmap
    double f,x,y,z;
    for(ind = 0; ind<stacks.size();ind++)
    {
      cout<<"ind = "<<ind<<endl;
      cout<<"_stack_group[ind] = "<<_stack_group[ind]<<endl;
      cout<<"_groups[index] = "<<_groups[index]<<endl;
      
      if(_stack_group[ind]==_groups[index])
      {
	cout<<"I am inside!"<<endl;
	image = stacks[ind];
        //irtkLinearInterpolateImageFunction interpolator;
	irtkSincInterpolateImageFunction interpolator;
        interpolator.SetInput(&image);
        interpolator.Initialize();
        
	cout<<"Correcting stack "<<ind<<endl;
	irtkImageAttributes attr = image.GetImageAttributes();
        for (t = 0; t < image.GetT(); t++) {
          for (k = 0; k < image.GetZ(); k++) {
            for (j = 0; j < image.GetY(); j++) {
              for (i = 0; i < image.GetX(); i++) {
		//find fieldmap value f
                x = i;
                y = j;
                z = k;
		image.ImageToWorld(x,y,z);
		fieldmap.WorldToImage(x,y,z);
		if ((x > -0.5) && (x < fieldmap.GetX()-0.5) && 
	            (y > -0.5) && (y < fieldmap.GetY()-0.5) &&
                    (z > -0.5) && (z < fieldmap.GetZ()-0.5) )
	        {
	         f = interpolatorFieldmap.EvaluateWithPadding(x,y,z,t,padding);
	        }
	        else
		  f=padding;
		
		//find displacement
		if (f>padding)  
		{
		  x = i;
		  y = j;
		  z = k;
                 if((index==1)&&(minus))
		  {
		    if(_swap[index])
	              y-=f/attr._dy;
	            else
	              x-=f/attr._dx;
		  }
		  else
		  {
		    if(_swap[index])
	              y+=f/attr._dy;
	             else
	               x+=f/attr._dx;
		   }
		 
		 if ((x > -0.5) && (x < image.GetX()-0.5) && 
	             (y > -0.5) && (y < image.GetY()-0.5) &&
                     (z > -0.5) && (z < image.GetZ()-0.5))
		 {
	             stacks[ind](i,j,k,t) = interpolator.Evaluate(x,y,z,t);
		 }
	         else
		   stacks[ind](i,j,k)=0;
		}
		else
		  stacks[ind](i,j,k)=0;
	      }//i
	    }//j
	  }//k
	}//t

	
        if (_debug)
        {
          sprintf(buffer,"cor%i.nii.gz",ind);
	  stacks[ind].Write(buffer);
        }
      }//if   
    }//ind
  }//index
}//function


void irtkReconstructionb0::CorrectMaskSmoothFieldmap(irtkRealImage& mask, irtkRealImage fieldmapMask, int group, bool minus )
{
  char buffer[256];
  //irtkMatrix m;
  unsigned int index, ind;
  int i,j,k,t;
  irtkRealImage fieldmap, fieldmap2, image;
  int padding;
  
  cout<<"Group = "<<group<<endl;
  //padding for fieldmap - careful vector = only copies reference
  fieldmap2 = _smoothFieldMap[group];
  fieldmap = fieldmap2;
  irtkRealPixel mn,mx;
  fieldmap.GetMinMax(&mn,&mx);
  padding = round(mn-2);
  cout<<"padding = "<<padding<<endl;
  irtkRealPixel *pf = fieldmap.GetPointerToVoxels();
  irtkRealPixel *pm = fieldmapMask.GetPointerToVoxels();
  for (j=0; j<fieldmap.GetNumberOfVoxels();j++)
  {
    if(*pm!=1)
      *pf = padding;
    pm++;
    pf++;
  }

  fieldmap.Write("orig-f-pad-mask.nii.gz");
    
  //prepare fieldmap for interpolation with padding
  irtkLinearInterpolateImageFunction interpolatorFieldmap;
  interpolatorFieldmap.SetInput(&fieldmap);
  interpolatorFieldmap.Initialize();


  //irtkImageFunction *interpolatorSinc = new irtkSincInterpolateImageFunction;
  //undistort stack using fieldmap
  double f,x,y,z;
  irtkLinearInterpolateImageFunction interpolator;
  //irtkSincInterpolateImageFunction interpolator;
  interpolator.SetInput(&mask);
  interpolator.Initialize();
        
   cout<<"Correcting mask. "<<endl;
   image = mask;
   irtkImageAttributes attr = image.GetImageAttributes();
     for (t = 0; t < image.GetT(); t++) {
       for (k = 0; k < image.GetZ(); k++) {
         for (j = 0; j < image.GetY(); j++) {
           for (i = 0; i < image.GetX(); i++) {
	    //find fieldmap value f
            x = i;
            y = j;
            z = k;
            image.ImageToWorld(x,y,z);
	    fieldmap.WorldToImage(x,y,z);
	    if ((x > -0.5) && (x < fieldmap.GetX()-0.5) && 
	        (y > -0.5) && (y < fieldmap.GetY()-0.5) &&
                (z > -0.5) && (z < fieldmap.GetZ()-0.5) )
	     {
	       f = interpolatorFieldmap.EvaluateWithPadding(x,y,z,t,padding);
	     }
	     else
		f=padding;
		
	     //find displacement
	     if (f>padding)  
	     {
		x = i;
		y = j;
		z = k;
               if((group==1)&&(minus))
	       {
		  if(_swap[group])
	            y-=f/attr._dy;
	          else
	            x-=f/attr._dx;
	       }
	       else
	       {
		  if(_swap[group])
	            y+=f/attr._dy;
	          else
	            x+=f/attr._dx;
	       }
		 
	      if ((x > -0.5) && (x < image.GetX()-0.5) && 
	          (y > -0.5) && (y < image.GetY()-0.5) &&
                  (z > -0.5) && (z < image.GetZ()-0.5))
	      {
	          image(i,j,k,t) = interpolator.Evaluate(x,y,z,t);
	      }
	      else
		 image(i,j,k)=0;
	      }
	      else
		image(i,j,k)=0;
	      }//i
	    }//j
	  }//k
	}//t

	mask=image;
        if (_debug)
        {
          sprintf(buffer,"mask%i-cor.nii.gz",group+1);
	  mask.Write(buffer);
        }
}//function



void irtkReconstructionb0::BSplineReconstructionGroup(int g)
{
  vector<irtkRealImage> slices;
  vector<irtkRigidTransformation> transformations;
  
  irtkRealImage slice,b;
  int i,j;
  
  double scale;
  
  int inputIndex;
  for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex)
  {
    //correct and exclude slices  
    if ((_slice_weight[inputIndex]>=0.5)&&(_stack_group[_stack_index[inputIndex]]==g))
    {
      cout<<inputIndex<<" ";
      // read the current slice
      slice=_slices[inputIndex];
      //read the current bias image
      b=_bias[inputIndex];
      //identify scale factor
        scale = _scale[inputIndex];
    
      //correct the slice      
      for (i=0;i<slice.GetX();i++)
        for (j=0;j<slice.GetY();j++)
          if (slice(i,j,0)!=-1)
	    slice(i,j,0)*=exp(-b(i,j,0))*scale;
       //prepare slices for BSpline reconstruction
       slices.push_back(slice);
       transformations.push_back(_transformations[inputIndex]);
    }
  }

  _bSplineReconstruction.Reconstruct(6,1,_reconstructed,slices,transformations);
  _reconstructed.Write("reconBSpline.nii.gz");
}

void irtkReconstructionb0::RememberDistortion()
{
  _BSplineField.AddImage(_distortion);
}

void irtkReconstructionb0::CombineDistortion()
{
  _distortion=_BSplineField.Average(); 
  _BSplineField.Clear();
  _distortion.Write("average-distortion.nii.gz");
}


void irtkReconstructionb0::CombineFieldmaps()
{
  irtkMultiChannelImage mch;
  mch.AddImage(_smoothFieldMap[0]);
  mch.AddImage(_smoothFieldMap[1]);
  _smoothFieldMap[0]=0;
  _smoothFieldMap[1]=0;
  _smoothFieldMap[0]+=mch.Average();
  _smoothFieldMap[1]+=mch.Average();
  _smoothFieldMap[0].Write("average-fieldmap-0-it4.nii.gz");
  _smoothFieldMap[1].Write("average-fieldmap-1-it4.nii.gz");
}
