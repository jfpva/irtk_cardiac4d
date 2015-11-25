/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: jacobian.cc 337 2011-06-04 18:41:48Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2011-06-04 19:41:48 +0100 (Sat, 04 Jun 2011) $
  Version   : $Revision: 337 $
  Changes   : $Author: ws207 $

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

// Default filenames
char *input_name = NULL, *output_name, *dof_name  = NULL;

void usage()
{
  cerr << "Usage: displacement [input] [output] [ffd]\n" << endl;
  cerr << "<-padding value>     Padding value" << endl;
  cerr << "<-world>             Use world coordinates [default: image]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, ok, padding;
  double x, y, z, xx, yy, zz;
  bool world;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  dof_name = argv[1];
  argc--;
  argv++;

  // Initialize padding value
  padding = MIN_GREY;
  world = false;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-world") == 0)) {
      argc--;
      argv++;
      world = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image
  cout << "Reading image ... "; cout.flush();
  irtkRealImage image = irtkRealImage(input_name);
  irtkImageAttributes attr = image.GetImageAttributes();
  attr._t=3;
  irtkRealImage displacement(attr);
  cout << "done" << endl;

  // Read transformation
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffd;
  mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dof_name);
  cout << "done" << endl;

  for (k = 0; k < image.GetZ(); k++) {
    for (j = 0; j < image.GetY(); j++) {
      for (i = 0; i < image.GetX(); i++) {
        if (image.Get(i, j, k) > padding) {
          x = i;
          y = j;
          z = k;
	  xx = x;
	  yy = y;
	  zz = z;
          image.ImageToWorld(x, y, z);
          image.ImageToWorld(xx, yy, zz);
          mffd->Transform(xx, yy, zz);
	  if(world)
	  {
	    displacement.Put(i,j,k,0,xx-x);
	    displacement.Put(i,j,k,1,yy-y);
	    displacement.Put(i,j,k,2,zz-z);
	  }
	  else
	  {
	    image.WorldToImage(xx,yy,zz);
	    displacement.Put(i,j,k,0,(xx-i)*attr._dx);
	    displacement.Put(i,j,k,1,(yy-j)*attr._dy);
	    displacement.Put(i,j,k,2,(zz-k)*attr._dz);
	  }
        }
      }
    }
  }

  // Write the final transformation estimate
  displacement.Write(output_name);
}
