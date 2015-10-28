#include <irtkImage.h>
#include <irtkLaplacianSmoothing.h>

irtkRealImage _image, _fieldmap, _mask, _target,_image_mask;
bool ok, have_target=false;
double relax_iter=5;
double lap_threshold=0.000001;
double rel_diff_threshold=0.0001;
double boundary_weight=1;

void usage()
{
  cerr<<"fit-fieldmap [image] [output] <-mask mask> <-relax_iter iter> <-lap_threshold threshold> <-rel_diff rel_diff>"<<endl;
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
    
    _mask=_image;
    _mask=1;
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

      if ((ok == false) && (strcmp(argv[1], "-image_mask") == 0)){
      argc--;
      argv++;
      _image_mask.Read(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

      if ((ok == false) && (strcmp(argv[1], "-relax_iter") == 0)){
      argc--;
      argv++;
      relax_iter =  atof(argv[1]);     
      argc--;
      argv++;
      ok = true;
    }

      if ((ok == false) && (strcmp(argv[1], "-lap_threshold") == 0)){
      argc--;
      argv++;
      lap_threshold = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

      if ((ok == false) && (strcmp(argv[1], "-rel_diff") == 0)){
      argc--;
      argv++;
      rel_diff_threshold = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

      if ((ok == false) && (strcmp(argv[1], "-boundary_weight") == 0)){
      argc--;
      argv++;
      boundary_weight = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
   
  irtkLaplacianSmoothing smoothing;
  smoothing.SetInput(_image);
  smoothing.SetMask(_mask);  
  smoothing.SetParam(lap_threshold,rel_diff_threshold,relax_iter,boundary_weight);
  _fieldmap = smoothing.RunGD();
  //_fieldmap = smoothing.Run();
  //_fieldmap = smoothing.Run1level();
  _fieldmap.Write(output_name);

}
