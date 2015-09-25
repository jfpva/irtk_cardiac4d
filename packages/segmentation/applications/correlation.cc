/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: blur.cc 772 2013-03-15 14:46:38Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2013-03-15 14:46:38 +0000 (Fri, 15 Mar 2013) $
  Version   : $Revision: 772 $
  Changes   : $Author: ws207 $

=========================================================================*/

#include <irtkImage.h>
#include <irtkFileToImage.h>

#include <irtkGaussianBlurring.h>
#include <irtkMultiChannelImage.h>


char *input_name = NULL, *output_name = NULL, *file_name = NULL, *par_name = NULL, *first_par_name = NULL, *par2_name = NULL;

void usage()
{
  cerr << "Usage: correlation [im1] [im2] <options>" << endl;
  cerr << "where <options> are one or more of the following:\n";
  cerr << "\t<-padding p>       Set padding." << endl;
  cerr << "\t<-file_name name>  Name of output file." << endl;
  cerr << "\t<-par_name value>  Value of parameter for inner loop." << endl;
  cerr << "\t<-first_par_name value>  First value of parameter for inner loop." << endl;
  cerr << "\t<-par2_name value>  Value of parameter for outer loop." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  int padding = 0;

  if (argc < 3) {
    usage();
  }

  // Parse parameters
  input_name  = argv[1];
  argc--;
  argv++;
  cout<<"Input 1 image: "<<input_name<<endl;
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"Input 2 image: "<<output_name<<endl;

  while (argc > 1) {
    ok=false;
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;

      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-file_name") == 0)) {
      argc--;
      argv++;
      file_name = argv[1];
      argc--;
      argv++;

      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-par_name") == 0)) {
      argc--;
      argv++;
      par_name = argv[1];
      argc--;
      argv++;

      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-first_par_name") == 0)) {
      argc--;
      argv++;
      first_par_name = argv[1];
      argc--;
      argv++;

      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-par2_name") == 0)) {
      argc--;
      argv++;
      par2_name = argv[1];
      argc--;
      argv++;

      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }
  cout<<"Padding: "<<padding<<endl;

  //read input images
  irtkRealImage x,y;
  x.Read(input_name);
  y.Read(output_name);
  irtkRealPixel *px = x.GetPointerToVoxels();
  irtkRealPixel *py = y.GetPointerToVoxels();

  double sum = 0;
  int num = 0;
  double value;
  
  for (int i=0; i<x.GetNumberOfVoxels();i++)
  {
    if(*px!=padding)
    {
      value = *px * *py; 
      //*px = value;
      sum += value;
      num++;
    }
    px++;
    py++;
  }

  //x.Write("lc.nii.gz");
  double cor = sum/num;
  cout<<"correlation is "<<cor<<endl;
  
  if (file_name != NULL)
  {
    ofstream fileOut(file_name, ofstream::out | ofstream::app);

    if(!fileOut)
    {
      cerr << "Can't open file " << output_name << endl;
      exit(1);
    }

    fileOut.precision(6);

    if (par_name == NULL)
    {
      if ( par2_name != NULL )
      {
        fileOut << par2_name << "," ;
      }

      fileOut << cor <<endl;
    }
    else
    {
      if ( first_par_name != NULL)
      {
        if ( strcmp(par_name, first_par_name) == 0 )
        {
          fileOut << endl;
          if ( par2_name != NULL )
          {
            fileOut << par2_name << "," ;
          }
        }
      }

      fileOut << cor <<",";

    }
    
  }
  
  return 0;
}
