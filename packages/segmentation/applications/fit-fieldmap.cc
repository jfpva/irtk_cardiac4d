#include <irtkImage.h>
#include <irtkLaplacianSmoothing.h>

irtkRealImage _image, _fieldmap, _mask, _target,_image_mask;
bool ok, have_target=false;

void usage()
{
  cerr<<"fit-fieldmap [image] [output] <-mask mask> <-target target>"<<endl;
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

      if ((ok == false) && (strcmp(argv[1], "-target") == 0)){
      argc--;
      argv++;
      _target.Read(argv[1]);
      have_target=true;
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
  //_fieldmap = smoothing.Run();
  _fieldmap = smoothing.Run1level();
  _fieldmap.Write(output_name);

}
