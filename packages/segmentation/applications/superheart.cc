#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkReconstruction.h>
#include <vector>

irtkRealImage _image, *_mask=NULL;
double frame_duration, cardiac_cycle;
vector<double > times;
double delta = 50, lambda = 0.01;
double resolution = 0;

bool ok;
char buffer[256];

void usage()
{
  cerr<<"superhart [image] [output] [frame duration] [cardiac cycle]"<<endl;
  cerr<<"         <-mask mask> <-times n time_1 ... time_n> <-delta delta> <-lambda lambda>"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 3)
    usage();
    
    _image.Read(argv[1]);
    cout<<"Reading image "<<argv[1]<<endl;
    argc--;
    argv++;

    char *output_name = argv[1];
    argc--;
    argv++;
    cout<<"Output name "<<output_name<<endl;
    
    frame_duration = atof(argv[1]);
    argc--;
    argv++;
    cout<<"Frame duration is "<<frame_duration<<"ms."<<endl;

    cardiac_cycle = atof(argv[1]);
    argc--;
    argv++;
    cout<<"Cardiac cycle is "<<cardiac_cycle<<"ms."<<endl;

    
    // Parse options.
    while (argc > 1){
      ok = false;
    
      if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      cout<<"Reading mask "<<argv[1]<<" ...";
      cout.flush();
      _mask= new irtkRealImage(argv[1]);
      _mask->Write("_mask.nii.gz");
      cout<<" done."<<endl;
      cout.flush();
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-times") == 0)){
      argc--;
      argv++;
      
      int n = atoi(argv[1]);
      argc--;
      argv++;
      cout<<"Number of frames ... "<<n<<endl;

      cout<< "Times are ";
      for (int i=0;i<n;i++)
      {
        times.push_back(atof(argv[1]));
	cout<<times[i]<<" ";
        argc--;
        argv++;
       }
       cout<<"."<<endl;
      ok = true;
    }

    //Parameter to define what is an edge
    if ((ok == false) && (strcmp(argv[1], "-delta") == 0)){
      argc--;
      argv++;
      delta=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }

    if ((ok == false) && (strcmp(argv[1], "-lambda") == 0)){
      argc--;
      argv++;
      lambda=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }

   //Isotropic resolution for the reconstructed volume
    if ((ok == false) && (strcmp(argv[1], "-resolution") == 0)){
      argc--;
      argv++;
      resolution=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

    
  irtkReconstruction reconstruction;
  reconstruction.DebugOn();
  
  irtkImageAttributes attr = _image.GetImageAttributes();
  
  double no_frames =  attr._t;
  cout<<"Number of frames is "<<no_frames<<endl;
  attr._z = attr._t*2;
  attr._t=1;
  
  double timestep;
  if (resolution == 0)
  {
    timestep = attr._dx;
    if(attr._dy<timestep)
      timestep=attr._dy;
  }
  else
    timestep=resolution;
  attr._dz = timestep;
  
  cout<<"Virtual time step is "<<timestep<<endl;
  cout<<"Virtual cardiac cycle is "<<timestep*no_frames<<endl;
  cout<<"Real time step is "<<cardiac_cycle/no_frames<<endl;

  double sigma = timestep*no_frames*frame_duration/cardiac_cycle;
  cout<<"Virtual frame duration is "<<sigma<<endl;
  //sigma=0.75*sigma;
  cout<<"FWHM of time frame is "<<sigma<<endl;
  
  if(times.size()==0)
  {
    for(int i=0;i<no_frames;i++)
      times.push_back(i);
  }
  else
  {
    cout<<"Time frames: ";
    for(int i=0;i<times.size();i++)
    {
      times[i]*=(no_frames/cardiac_cycle);
      cout<<times[i]<<" ";
    }
    cout<<endl;
  }
  
  reconstruction.SetSmoothingParameters(delta,lambda);
  
  irtkRealImage temp(attr);
  irtkRealImage mask(attr);
  
  for(int t=0;t<_image.GetT();t++)
    for(int j=0;j<_image.GetY();j++)
      for(int i=0;i<_image.GetX();i++)
      {
	temp(i,j,t)=_image(i,j,0,t);
	temp(i,j,t+_image.GetT())=_image(i,j,0,t);
        if(_mask!=NULL)
        {
	  mask(i,j,t)=_mask->Get(i,j,0);
	  mask(i,j,t+_image.GetT())=_mask->Get(i,j,0);
        }
      }
          
   temp.Write("temp.nii.gz");
   mask.Write("mask.nii.gz");
 
   if (resolution == 0)
     resolution = reconstruction.CreateTemplate(temp,0);
   else
     resolution = reconstruction.CreateTemplate(temp,resolution);
   
   irtkRealImage reconstructed = reconstruction.GetReconstructed();
   reconstructed.Write("template.nii.gz");
   
   if(_mask==NULL)
   {
     mask = reconstructed;
     mask=1;
     reconstruction.SetMask(&mask,0);
   }
   else
     reconstruction.SetMask(&mask,0);
     
  /// Slice stacks
  vector<irtkRealImage> slices;
  vector<irtkRigidTransformation> transformations;
  vector<double> thickness;

 for(int t=0; t<_image.GetT();t++)
 {
     slices.push_back(_image.GetRegion(0, 0, 0, t, attr._x, attr._y, 1,t+1));
     irtkRigidTransformation tr;
     tr.PutTranslationZ((times[t]-_image.GetT()+0.5)*timestep);
     transformations.push_back(tr);
     thickness.push_back(sigma);
 }
 
 for(int t=0; t<_image.GetT();t++)
 {
     slices.push_back(_image.GetRegion(0, 0, 0, t, attr._x, attr._y, 1,t+1));
     irtkRigidTransformation tr;
     tr.PutTranslationZ((times[t]+0.5)*timestep);
     transformations.push_back(tr);
     thickness.push_back(sigma);
 }

 
 reconstruction.CreateSlicesAndTransformations(slices,transformations,thickness);
 reconstruction.MaskSlices();
 //reconstruction.SaveSlices();
 //reconstruction.SaveTransformations();
 
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
 
 reconstruction.MaskVolume();
 reconstructed = reconstruction.GetReconstructed();
 reconstructed.Write("reconstructed.nii.gz"); 
 
 
 double x=0,y=0,z;
 //time 0
 cout<<"time 0: ";
 z=no_frames;
 temp.ImageToWorld(x,y,z);
 cout<<z<<" in world coord."<<endl;
 reconstructed.WorldToImage(x,y,z);
 cout<<z<<" in recon coord."<<endl;
 int z0=round(z);
 
 /*
 cout<<"All times:"<<endl;
 
 for (int i=0; i<reconstructed.GetZ();i++)
 {
   cout<<"Recon frame "<<i<<": ";
   z=i;
   reconstructed.ImageToWorld(x,y,z);
   temp.WorldToImage(x,y,z);
   cout<<z<<" in temp; ";
   if(z>=no_frames)
     z=z-no_frames;
   cout<<z<<" in no frames; ";
   double time = z*cardiac_cycle/no_frames;
   cout<<time<<"ms in cardiac_cycle"<<endl;
   
 }
 */
 
 irtkImageAttributes attr2 = reconstructed.GetImageAttributes();
 attr = _image.GetImageAttributes();
 attr2._t=no_frames;
 attr2._z=1;
 attr2._dz=attr._dz;
 attr2._dt=cardiac_cycle/no_frames;
 
 irtkRealImage result(attr2);
 
 for(int ind=z0-(no_frames/2);ind<z0+(no_frames/2)-1;ind++)
 {
   cout<<"recon"<<ind<<"; ";
   z=ind;
   reconstructed.ImageToWorld(x,y,z);
   temp.WorldToImage(x,y,z);
   if(z>=no_frames)
     z=z-no_frames;
   int t = round(z);
   cout<<"frame "<<t<<endl;
   for(int i=0;i<reconstructed.GetX();i++)
     for(int j=0;j<reconstructed.GetY();j++)
     {
       result(i,j,0,t)=reconstructed(i,j,ind);  
     }
 }
 
 cout<<"Timing of frames in ms: ";
 for(int i=0;i<result.GetT();i++)
 {
   cout<<i*cardiac_cycle/no_frames<<" ";
 }
 
 result.Write(output_name);  

}
