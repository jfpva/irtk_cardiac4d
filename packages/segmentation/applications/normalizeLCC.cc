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


char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: normalizeLCC [in] [out] [sigma] <options>" << endl;
  cerr << "where <options> are one or more of the following:\n";
  cerr << "\t<-padding p>       Set padding." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  double sigma;
  int padding = 0;
  bool have_padding = false;

  if (argc < 4) {
    usage();
  }

  // Parse parameters
  input_name  = argv[1];
  argc--;
  argv++;
  cout<<"Input image: "<<input_name<<endl;
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"Output image: "<<output_name<<endl;
  sigma = atof(argv[1]);
  argc--;
  argv++;
  cout<<"Sigma: "<<sigma<<endl;

  while (argc > 1) {
    ok=false;
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      have_padding=true;
      padding = atoi(argv[1]);
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

  irtkRealImage input;
  input.Read(input_name);
  
  // Blur image
  //if(have_padding)
    irtkGaussianBlurringWithPadding<irtkRealPixel> gaussianBlurring(sigma, padding);
  //else
    //irtkGaussianBlurring<irtkRealPixel> gaussianBlurring(sigma);
  
  irtkImageAttributes attr = input.GetImageAttributes();
  cout<<"Calculating image mean ...";
  irtkRealImage mean(attr);
  gaussianBlurring.SetInput (&input);
  gaussianBlurring.SetOutput(&mean);
  gaussianBlurring.Run();
  mean.Write("mean.nii.gz");
  cout<<"done.";
  
  irtkMultiChannelImage mch;
  mch.AddImage(input);
  mch.AddImage(mean);
  mch.SetPadding(padding);
  mch.CreateMask();
  mch.Brainmask();
  irtkRealImage var = mch.Subtract();
  var.Write("dif.nii.gz");
  mch.SetImage(0,var);
  mch.SetImage(1,var);
  var = mch.Multiply();
  var.Write("dif2.nii.gz");

  gaussianBlurring.SetInput (&var);
  gaussianBlurring.SetOutput(&var);
  gaussianBlurring.Run();
  var.Write("var.nii.gz");
  
  mch.SetImage(0,var);
  mch.Sqrt(0);
  var = mch.GetImage(0);
  var.Write("std.nii.gz");
  
  mch.SetImage(0,input);
  mch.SetImage(1,mean);
  irtkRealImage norm(attr);
  
  norm = mch.Subtract();
  norm.Write("norm-mean.nii.gz");
  
  mch.SetImage(0,norm);
  mch.SetImage(1,var);
  norm=mch.Divide(1);
  
  mch.SetPadding(-10);
  mch.SetImage(0,norm);
  mch.Brainmask();
  mch.Write(0,output_name);
  

  return 0;
}
