#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkReconstruction.h>

irtkRealImage _image;
double sigma, timestep;
bool ok;
char buffer[256];

void usage()
{
  cerr<<"superhart [image] [output] [sigma] [timestep]"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 3)
    usage();
    
    _image.Read(argv[1]);
    argc--;
    argv++;

    char *output_name = argv[1];
    argc--;
    argv++;
    
    sigma = atof(argv[1]);
    argc--;
    argv++;

    timestep = atof(argv[1]);
    argc--;
    argv++;

    /*
    // Parse options.
    while (argc > 1){
      ok = false;
    
      if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      _mask.Read(argv[1]);
      argc--;
      argv++;
      ok = true;
    }


    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
*/
    
  irtkReconstruction reconstruction;
  reconstruction.DebugOn();
  
  irtkImageAttributes attr = _image.GetImageAttributes();
  
  attr._z = attr._t*2;
  attr._t=1;
  attr._dz = timestep;
  
  
  irtkRealImage temp(attr);
  
  for(int t=0;t<_image.GetT();t++)
    for(int j=0;j<_image.GetY();j++)
      for(int i=0;i<_image.GetX();i++)
      {
	temp(i,j,t)=_image(i,j,0,t);
	temp(i,j,t+_image.GetT())=_image(i,j,0,t);
      }
          
   temp.Write("temp.nii.gz");
 
     double res = reconstruction.CreateTemplate(temp,0);
     irtkRealImage reconstructed = reconstruction.GetReconstructed();
 
   reconstructed.Write("template.nii.gz");
   irtkRealImage mask = reconstructed;
   mask=1;
   
   reconstruction.SetMask(&mask,0);
     
  /// Slice stacks
  vector<irtkRealImage> slices;
  vector<irtkRigidTransformation> transformations;
  vector<double> thickness;

 for(int t=0; t<_image.GetT();t++)
 {
     slices.push_back(_image.GetRegion(0, 0, 0, t, attr._x, attr._y, 1,t+1));
     irtkRigidTransformation tr;
     tr.PutTranslationZ((t-_image.GetT()+0.5)*timestep);
     transformations.push_back(tr);
     thickness.push_back(sigma);
 }
 
 for(int t=0; t<_image.GetT();t++)
 {
     slices.push_back(_image.GetRegion(0, 0, 0, t, attr._x, attr._y, 1,t+1));
     irtkRigidTransformation tr;
     tr.PutTranslationZ((t+0.5)*timestep);
     transformations.push_back(tr);
     thickness.push_back(sigma);
 }

 
 ///test
 //slices.push_back(temp);
 //irtkRigidTransformation tr;
 //transformations.push_back(tr);
 //thickness.push_back(sigma);
 ///end test
 
 reconstruction.CreateSlicesAndTransformations(slices,transformations,thickness);
 reconstruction.SaveSlices();
 reconstruction.SaveTransformations();
 
 reconstruction.SpeedupOff();
 cout<<"InitializeEM ...";
 reconstruction.InitializeEM();
 reconstruction.InitializeEMValues();
 cout<<"done."<<endl;
 cout.flush();
 
 
 cout<<"Sizes:"<<slices.size()<<" "<<transformations.size()<<" "<<thickness.size()<<endl;
 cout<<"CoeffInit ...";
  cout.flush();
 reconstruction.CoeffInit();
 cout<<"done."<<endl;
 cout.flush();

 cout<<"GaussianReconstruction ...";
 reconstruction.GaussianReconstruction();
 cout<<"done."<<endl;
 cout.flush();
 reconstructed = reconstruction.GetReconstructed();
 reconstructed.Write("init.nii.gz"); 
 
 reconstruction.SimulateSlices();
 reconstruction.InitializeRobustStatistics();
 
 for(int i = 0; i<10; i++)
 {
   cout<<endl<<"  Reconstruction iteration "<<i<<". "<<endl;
   reconstruction.Superresolution(i+1);
   reconstruction.SimulateSlices();
   
   reconstructed=reconstruction.GetReconstructed();
   sprintf(buffer,"super%i.nii.gz",i);
   reconstructed.Write(buffer);
 }
 
 reconstructed = reconstruction.GetReconstructed();
 reconstructed.Write(output_name); 
  

}
