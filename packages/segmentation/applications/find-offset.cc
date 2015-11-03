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
  cerr<<"find-offset [target] [source] [registered_source] [dofout]"<<endl;
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

    char *output_name = argv[1];
    argc--;
    argv++;
    
    char *output_dof_name = argv[1];
    argc--;
    argv++;
    
    // Parse options.
    
    while (argc > 1){
      ok = false;
    

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
    irtkGreyImage gt(t),gs(s);
    registration.SetInput(&gt,&gs);
    registration.SetOutput(&dof);
    registration.GuessParameterDistortion(1);
    registration.SetTargetPadding(0);
    registration.Run();
    registration.Write("par.areg");
    
    dof.Print();
    
    cout<<"offset="<<dof.GetTranslationX()<<endl;
  
    //re-orient transformation
    irtkMatrix mo = orient.GetMatrix();
    irtkMatrix md = dof.GetMatrix();
    md = md*mo;
    mo.Invert();
    md = mo*md;
    dof.PutMatrix(md);
    
    dof.irtkTransformation::Write("translation.dof");
    
}