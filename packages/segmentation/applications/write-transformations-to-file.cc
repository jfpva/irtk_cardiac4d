/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: reconstructionb0.cc 1000 2013-10-18 16:49:08Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-10-18 17:49:08 +0100 (Fri, 18 Oct 2013) $
  Version   : $Revision: 1000 $
  Changes   : $Author: mm3 $

=========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <string>
using namespace std;

//Application to perform reconstruction of volumetric MRI from thick slices.

void usage()
{
  cerr << "Usage: write-transformations-to-file [file] [tr1 ... trn]\n" << endl;
  cerr << endl;
  exit(1);
}

int main(int argc, char **argv)
{
      
  //if not enough arguments print help
  if (argc < 3)
    usage();
  
  //read file name
  char* file_name=argv[1];
  argc--;
  argv++;
  
  cout<<"Writing to file "<<file_name<<endl;

  //create file
  ofstream fileOut(file_name, ofstream::out | ofstream::app);
  if(!fileOut)
  {
    cerr << "Can't open file " << file_name << endl;
    exit(1);
  }
  fileOut.precision(3);  
 
  //number of transformations
  int n = argc-1;
  cout<<"We have "<<n<<"transformations."<<endl;
  
  //loop through transformations
  
  irtkRigidTransformation t;
  for(int i=0;i<n;i++)
  {
    t.irtkTransformation::Read(argv[1]);
    argc--;
    argv++;

    fileOut<<t.GetTranslationX()<<" ";
    fileOut<<t.GetTranslationY()<<" ";
    fileOut<<t.GetTranslationZ()<<" ";
    fileOut<<t.GetRotationX()<<" ";
    fileOut<<t.GetRotationY()<<" ";
    fileOut<<t.GetRotationZ()<<endl;
  }  
  //The end of main()
}  
