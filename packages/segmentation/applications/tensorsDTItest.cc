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

vector<irtkRealImage> images;
double directions[3][64];
double bvalues[64];
int nDir;
int nImages,noEntries;

int i,j,k;
double x,y,z;
irtkTensor D;
double num;

void usage()
{
  cerr<<"tensorDTI [noImages][4Dimage 1] ... [4Dimage n] [bvalues] [directions] "<<endl;
}

int main(int argc, char **argv)
{
    if (argc < 4)
    {
      usage();
      exit(1);
    }
    
    //read number of images
    nImages = atoi(argv[1]);
    argc--;
    argv++;
    
    //read the images
    for (i=0;i<nImages; i++)
    {
      irtkRealImage im(argv[1]);
      images.push_back(im);
      argc--;
      argv++;
    }
    
    char *btextfile = argv[1];
    argc--;
    argv++;

    char *textfile = argv[1];
    argc--;
    argv++;
    
    
    nDir = images[0].GetT();
    cout<<"Number of directions is "<<nDir-1<<endl;
    
    noEntries = nDir*nImages;
    cout<<"Number of entries is "<<noEntries<<endl;

    //Read the b-values from the text file
    ifstream inb(btextfile);
    i = 0;
    cout<<"Reading b-values: "<<endl;
    if (inb.is_open()) 
    {
      while (!inb.eof()) 
      {

	inb >> num;
	if (i<noEntries)
	  bvalues[i]=num;
	cout<<num<<" ";
	i++;
      }
      cout<<endl;
      inb.close();
    } 
    else 
    {
      cout << "tensorDTI: Unable to open file " << btextfile << endl;
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
	if ((coord<3)&&(dir<noEntries))
	  directions[coord][dir]=num;
	cout<<num<<" ";
	dir++;
	if (dir>=noEntries)
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
      cout << "tensorDTI: Unable to open file " << textfile << endl;
      exit(1);
    }
    /*
    cout<<"test."<<endl;
    irtkMatrix test(6,2);

    test(0,1)=0;
    test(0,0)=0;
    
    test(1,1)=0.166*3.1415;
    test(1,0)=0.333*3.1415;
    
    test(2,1)=0.333*3.1415;
    test(2,0)=0.666*3.1415;
    
    test(3,1)=0.5*3.1415;
    test(3,0)=1*3.1415;
    
    test(4,1)=0.666*3.1415;
    test(4,0)=1.333*3.1415;
    
    test(5,1)=0.833*3.1415;
    test(5,0)=1.666*3.1415;
    
    
    test.Print();
    
    irtkSphericalHarmonics sh;
    irtkMatrix sht = sh.shTransform(test,2);
    sht.Print();
    sht.Invert();
    sht.Print();
    
    irtkVector c(6),s(6);
    
    for(int i=0;i<6;i++)
      s(i)=1;
    s.Print();
    
    c=sht*s;
    c.Print();
    */
    
    ////////////////real directions
    irtkMatrix dirs_xyz(nDir-1,3);
    for(int i=0;i<nDir-1;i++)
      for(int j=0;j<3;j++)
	dirs_xyz(i,j) = directions[j][i+1];
    
   dirs_xyz.Print();
   
   irtkSphericalHarmonics sh;
   /*
   irtkMatrix dirs_sph;
   dirs_sph = sh.cartesian2spherical(dirs_xyz);
   dirs_sph.Print();
   irtkMatrix sht = sh.shTransform(dirs_sph,2);
   sht.Print();
   
   
   irtkVector c(nDir-1),s(nDir-1);
   for(int i=0;i<(nDir-1);i++)
      s(i)=1;
    s.Print();
   
   irtkMatrix temp = sht;
   temp.Transpose();
   c=temp*s;
   temp=temp*sht;
   temp.Invert();
   c=temp*c;
   */

   irtkVector s(nDir-1);
   for(int i=0;i<(nDir-1);i++)
   s(i)=1;
   s.Print();

   irtkVector c = sh.LSFit(s,dirs_xyz,2);
   c.Print();
   
   sh.InitSHT(dirs_xyz,2);
   irtkVector cc = sh.Signal2Coeff(s);
   cc.Print();
   
   irtkVector ss = sh.Coeff2Signal(c);
   ss.Print();
   
   /*
   irtkMatrix d(1,3);
   irtkMatrix sph(1,2);
   d(0,0)=-1;
   d(0,1)=0;
   d(0,2)=0;
   sph = sh.cartesian2spherical(d);
   sph.Print();
   */

    
    
/*
    
    //create structures for DW images
    vector<irtkDWImage> DWimages;
    cout<<endl;
    for(i=0;i<images.size();i++)
    {
      cout<<"Image "<<i<<": "<<endl;
      irtkDWImage dw = images[i];
      dw.SetB(bvalues[i*nDir+1]);
      cout<<"b-value "<<bvalues[i*nDir+1]<<endl;
      irtkMatrix d(3,nDir-1);
      for (j=0;j<3;j++)
        for (k=0;k<(nDir-1);k++)
	  d(j,k)=directions[j][k+i*nDir+1];
      cout<<"Directions: "<<endl;
      d.Print();
      dw.SetDirections(d);
      DWimages.push_back(dw);
    }
    
    images.clear();
    
    //calculate directions and FA maps
    
    irtkImageAttributes attr;
    
    for (uint index=0; index<DWimages.size();index++)//DWimages.size()
    {
      char buffer[256];
      sprintf(buffer,"image%i.nii.gz",index);
      DWimages[index].Write(buffer);
      attr= DWimages[index].GetImageAttributes();
    
      //principal direction
      attr._t = 3;
      irtkRealImage PD(attr);
      
      //FA map
      attr._t = 1;
      irtkRealImage FA(attr);
      
      //Tensor map
      attr._t = 6;
      irtkRealImage DTI(attr);
      
      //irtkTensorField field;
      //field.Initialise(DWimages[index]);
      irtkTensorField field(DWimages[index]);
    
      irtkTensor D;
      irtkVector dr;
      for(i=0;i<DWimages[index].GetX();i++)
        for(j=0;j<DWimages[index].GetY();j++)
          for(k=0;k<DWimages[index].GetZ();k++)
	    if(DWimages[index](i,j,k,0)>250)
	    {
	      D = field(i,j,k);//DWimages[index].CalculateTensor(i,j,k);
	      dr = D.Direction();
	      PD(i,j,k,0)=dr(0);
	      PD(i,j,k,1)=dr(1);
	      PD(i,j,k,2)=dr(2);
	      FA(i,j,k)=D.FA();
	      DTI(i,j,k,0)=D(0);
	      DTI(i,j,k,1)=D(1);
	      DTI(i,j,k,2)=D(2);
	      DTI(i,j,k,3)=D(3);
	      DTI(i,j,k,4)=D(4);
	      DTI(i,j,k,5)=D(5);
	    }
	  
      
      sprintf(buffer,"PD%i.nii.gz",index);
      PD.Write(buffer);
      sprintf(buffer,"FA%i.nii.gz",index);
      FA.Write(buffer);
      sprintf(buffer,"DTI%i.nii.gz",index);
      DTI.Write(buffer);
    }
    */
}    
