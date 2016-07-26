#include <irtkImage.h>
#include <irtkTensor.h>
#include <irtkDWImage.h>
#include <irtkTensorField.h>
#include <irtkResampling.h>
#include <irtkVector.h>
#include <irtkMatrix.h>


#include <irtkSphericalHarmonics.h>

#include <vector>
using namespace std;

double directions[3][64];
double bvalues[64];
int nDir;

int i,j,k;
double x,y,z;
double num;
int order = 2;
double lambda = 0;

void usage()
{
  cerr<<"signal2SH [4Dimage] [bvalues] [directions] [SH order] [lambda] "<<endl;
}

int main(int argc, char **argv)
{
    if (argc < 4)
    {
      usage();
      exit(1);
    }
    
    
    //read the image
    irtkRealImage image(argv[1]);
    argc--;
    argv++;
    
    char *btextfile = argv[1];
    argc--;
    argv++;

    char *textfile = argv[1];
    argc--;
    argv++;
    
    order = atoi(argv[1]);
    argc--;
    argv++;
    cout<<"order is "<<order<<endl;
    
    lambda = atof(argv[1]);
    argc--;
    argv++;
    cout<<"lambda is "<<lambda<<endl;
    
    nDir = image.GetT();
    cout<<"Number of directions is "<<nDir-1<<endl;
    
    //Read the b-values from the text file
    ifstream inb(btextfile);
    i = 0;
    cout<<"Reading b-values: "<<endl;
    if (inb.is_open()) 
    {
      while (!inb.eof()) 
      {

	inb >> num;
	if (i<=nDir)
	  bvalues[i]=num;
	cout<<num<<" ";
	i++;
      }
      cout<<endl;
      inb.close();
    } 
    else 
    {
      cout << "signal2SH: Unable to open file " << btextfile << endl;
      exit(1);
    }
    
    //Read directions from the text file
    ifstream in(textfile);
    int coord = 0;
    int dir = 0;
    cout<<"Reading directions: "<<endl;
    if (in.is_open()) 
    {
      while (!in.eof()) 
      {
	in >> num;
	if ((coord<3)&&(dir<=nDir))
	  directions[coord][dir]=num;
	cout<<num<<" ";
	dir++;
	if (dir>nDir)
	{
	  dir=0;
	  coord++;
	  cout<<endl;
	}
      }
      in.close();
    } 
    else 
    {
      cout << "signal2SH: Unable to open file " << textfile << endl;
      exit(1);
    }
    
    ////////////////real directions
    irtkMatrix dirs_xyz(nDir,3);
    for(int i=0;i<nDir;i++)
      for(int j=0;j<3;j++)
	dirs_xyz(i,j) = directions[j][i+1];
    
   dirs_xyz.Print();
   
   irtkSphericalHarmonics sh;
   //sh.LaplaceBeltramiMatrix(order);
   
   sh.InitSHTRegul(dirs_xyz,lambda,order);
   //sh.InitSHT(dirs_xyz,order);
   irtkRealImage shc = sh.Signal2Coeff(image);
   shc.Write("shCoeffs.nii.gz");
   
   irtkRealImage sim = sh.Coeff2Signal(shc);
   sim.Write("simsignal.nii.gz");
   
}    
