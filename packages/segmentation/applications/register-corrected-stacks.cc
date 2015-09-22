#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkReconstruction.h>
#include <irtkReconstructionb0.h>
#include <irtkRegistration.h>
#include <irtkImageAffineRegistrationWithPadding.h>



irtkRealImage _target, _source, _fieldmap;
vector<irtkRigidTransformation> _dofs;
vector<bool> _have_dof;
vector<double> _translations;
vector<int> _excluded;
bool ok;
char buffer[256];
double x,y,z;
int i,j,k,l;

bool _swap = true;
bool _minus = false;


void usage()
{
  cerr<<"register-corrected-stacks [target] [source] [fieldmap] [registered_source] [corrected_fieldmap]"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 3)
    usage();
    
    _target.Read(argv[1]);
    argc--;
    argv++;

    _source.Read(argv[1]);
    argc--;
    argv++;

    _fieldmap.Read(argv[1]);
    argc--;
    argv++;

    char *output_name = argv[1];
    argc--;
    argv++;
    
    char *output_fieldmap_name = argv[1];
    argc--;
    argv++;
    
    // Parse options.
    
    while (argc > 1){
      ok = false;
    
      if ((ok == false) && (strcmp(argv[1], "-exclude") == 0)){
        argc--;
        argv++;
        int num_excluded = atoi(argv[1]);
        argc--;
        argv++;
	for (int i=0;i<num_excluded; i++)
	{
	  int slice_num = atoi(argv[1]);
	  _excluded.push_back(slice_num);
	  argc--;
          argv++;
	}
        ok = true;
      }

      if ((ok == false) && (strcmp(argv[1], "-x") == 0)){
        argc--;
        argv++;
        _swap=false;
        ok = true;
      }

      if ((ok == false) && (strcmp(argv[1], "-minus") == 0)){
        argc--;
        argv++;
        _minus=true;
        ok = true;
      }

      if (ok == false){
        cerr << "Can not parse argument " << argv[1] << endl;
        usage();
      }
  }
  
 
  //resample source on the grid of target
  irtkNearestNeighborInterpolateImageFunction interpolator;
  irtkRealImage source(_target);
  source=0;
  irtkRigidTransformation id;
  irtkImageTransformation imagetransformation;
  imagetransformation.SetInput(&_source,&id);
  imagetransformation.SetOutput(&source);
  imagetransformation.PutInterpolator(&interpolator);
  imagetransformation.Run();
  source.Write("source.nii.gz");
  

  //adjust orientation
  irtkReconstructionb0 reconstruction;
  irtkRealImage t,s;

  t=reconstruction.AdjustOrientation(_target,_swap);
  s=reconstruction.AdjustOrientation(source,_swap);
  irtkAffineTransformation orient = reconstruction.AdjustOrientationTransformation(_target,_swap);
  t.Write("t.nii.gz");
  s.Write("s.nii.gz");
  orient.irtkTransformation::Write("o.dof");
  
  //choose slice
  irtkImageAttributes attr = t.GetImageAttributes();
  //for(i=62; i<63;i++)
  for(i=0; i<attr._z;i++)
  {
    cout<<i<<endl;
    //get the slice
    irtkRealImage slice = t.GetRegion(0,0,i,attr._x,attr._y,i+1);
    //sprintf(buffer,"slice%i.nii.gz",i);
    //slice.Write(buffer);

    //check whether the slice in empty
    irtkRealImage slice2 = s.GetRegion(0,0,i,attr._x,attr._y,i+1);
    //sprintf(buffer,"slice2-%i.nii.gz",i);
    //slice2.Write(buffer);
    irtkRealPixel smin,smax;
    slice2.GetMinMax(&smin,&smax);
    bool have_dof;
    if(smax==0)
    {
      _have_dof.push_back(false);
       have_dof=false;
    }
    else
    {
      _have_dof.push_back(true);
      have_dof=true;
    }
    
    //create dof
    irtkAffineTransformation dof;
    dof.PutStatus(TX,  _Active);
    dof.PutStatus(TY,  _Passive);
    dof.PutStatus(TZ,  _Passive);
    dof.PutStatus(RX,  _Passive);
    dof.PutStatus(RY,  _Passive);
    dof.PutStatus(RZ,  _Passive);
    dof.PutStatus(SX,  _Passive);
    dof.PutStatus(SY,  _Passive);
    dof.PutStatus(SZ,  _Passive);
    dof.PutStatus(SXY, _Passive);
    dof.PutStatus(SXZ, _Passive);
    dof.PutStatus(SYZ, _Passive);
  
    //register source to the target slice
    irtkImageAffineRegistrationWithPadding registration;
    irtkGreyImage gt(slice),gs(s);
    //gt.Write("gt.nii.gz");
    //gs.Write("gs.nii.gz");
    registration.SetInput(&gt,&gs);
    registration.SetOutput(&dof);
    registration.GuessParameterDistortion(1);
    registration.SetTargetPadding(0);
    registration.Run();
    registration.Write("par.areg");
  
    //sprintf(buffer,"dof%i.dof",i);
    //dof.irtkTransformation::Write(buffer);
    _translations.push_back(dof.GetTranslationX());

    //re-orient transformation
    irtkMatrix mo = orient.GetMatrix();
    irtkMatrix md = dof.GetMatrix();
    md = md*mo;
    mo.Invert();
    md = mo*md;
    dof.PutMatrix(md);
    
    //sprintf(buffer,"translation%i.dof",i);
    //dof.irtkTransformation::Write(buffer);
    
    //remember slice trasformation
    _dofs.push_back(dof);
  }
  
  cout<<"Translations: ";
  for(uint index=0; index<_translations.size();index++)
  {
    cout<<"("<<_translations[index]<<", "<<_have_dof[index]<<"); ";
  }

  //create registered source
  irtkLinearInterpolateImageFunction int_lin;
  int_lin.SetInput(&_source);
  int_lin.Initialize();
  
  irtkRealImage output(_target);
  output=0;
  
  for (l = 0; l < output.GetT(); l++) {
    for (k = 0; k < output.GetZ(); k++) {
      for (j = 0; j < output.GetY(); j++) {
        for (i = 0; i < output.GetX(); i++) {
          x = i;
          y = j;
          z = k;
	  output.ImageToWorld(x,y,z);
	  _dofs[k].Transform(x,y,z);
	  _source.WorldToImage(x,y,z);
	  if ((x > -0.5) && (x < output.GetX()-0.5) && 
	      (y > -0.5) && (y < output.GetY()-0.5) &&
              (z > -0.5) && (z < output.GetZ()-0.5))
	  {
	    //output(i,j,k,l) = int_lin.Evaluate(x,y,z,l);
	    output(i,j,k,l) = int_lin.EvaluateWithPadding2D(x,y,z,l,0);
	  }
        }
      }
    }
  }
  output.Write(output_name);
  
 //create registered fieldmap
  //irtkLinearInterpolateImageFunction int_lin;
  int_lin.SetInput(&_fieldmap);
  int_lin.Initialize();
  
  
  irtkRealImage output_fieldmap(_target);
  output_fieldmap = 0;
  attr = output_fieldmap.GetImageAttributes();
  double res = attr._dy;
  double xx,yy,zz;
  double adjustment,value;
  
  //make sure excluded slices apply
  for(i=0;i<_excluded.size();i++)
    _have_dof[_excluded[i]]=false;
  
  cout<<endl<<"Adjustment: ";
  for (l = 0; l < output_fieldmap.GetT(); l++) {
    for (k = 0; k < output_fieldmap.GetZ(); k++) {
      if(_have_dof[k])
      {
      for (j = 0; j < output_fieldmap.GetY(); j++) {
        for (i = 0; i < output_fieldmap.GetX(); i++) {
          //find fieldmap value
	  x = i;
          y = j;
          z = k;
	  output_fieldmap.ImageToWorld(x,y,z);
	  _dofs[k].Transform(x,y,z);
	  _fieldmap.WorldToImage(x,y,z);
	  //find adjustment
	  xx = i;
          yy = j;
          zz = k;
	  output_fieldmap.ImageToWorld(xx,yy,zz);
	  _dofs[k].Transform(xx,yy,zz);
	  output_fieldmap.WorldToImage(xx,yy,zz);
	  
	  if(_swap)
	    adjustment = (yy-j)*attr._dy;
	  else
	    adjustment = (xx-i)*attr._dx;
	  if((i==0)&&(j==0))
	    cout<<adjustment<<" ";
	  
	  if ((x > -0.5) && (x < _fieldmap.GetX()-0.5) && 
	      (y > -0.5) && (y < _fieldmap.GetY()-0.5) &&
              (z > -0.5) && (z < _fieldmap.GetZ()-0.5))
	  {
	    value = int_lin.EvaluateWithPadding2D(x,y,z,l,0);
	    if (_minus)
	      output_fieldmap(i,j,k,l) = -value;
	    else
	      output_fieldmap(i,j,k,l) = value;
	    if(output_fieldmap(i,j,k,l)!=0) 
	      if(_minus)
	        output_fieldmap(i,j,k,l)-=adjustment;
	      else
	       output_fieldmap(i,j,k,l)+=adjustment;
	  }
        }
      }
      }
    }
  }
  //if(_minus)
    //output_fieldmap = output_fieldmap*(-1);
  
  output.Write(output_name);
  output_fieldmap.Write(output_fieldmap_name);

}