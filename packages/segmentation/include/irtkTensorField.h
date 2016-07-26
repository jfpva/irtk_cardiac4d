/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkTensorField.h 577 2012-03-30 09:54:05Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2012-03-30 10:54:05 +0100 (Fri, 30 Mar 2012) $
  Version   : $Revision: 577 $
  Changes   : $Author: mm3 $

=========================================================================*/

#ifndef _irtkTensorField_H

#define _irtkTensorField_H

#include <irtkImage.h>
#include <irtkTensor.h>
#include <irtkDWImage.h>

#include <vector>
using namespace std;


class irtkTensorField : public irtkObject
{
  vector<vector<vector<irtkTensor> > > _tensorField;
  int _x;
  int _y;
  int _z
  ;
public:
  irtkTensorField();
  irtkTensorField(int x, int y, int z);
  irtkTensorField(irtkDWImage& image, double threshold = 0);
  void Initialise(irtkDWImage& image, double threshold = 0);
  irtkTensor operator () (int,int,int);
  irtkTensor operator () (double,double,double);
  irtkTensor Get(int,int,int);
  bool InROI(double x, double y, double z);
  void SmoothOld(double lambda = 0.25);
  void Smooth(double lambda = 0.25);
  void EdgePreserveSmooth(double lambda = 0.25, double delta = 0);
  void Log();
  //void Exp();
  
  void Save(irtkRealImage b0, char * name);
  void SaveFA(irtkRealImage b0, char * name);
  void SaveDirections(irtkRealImage b0, char * name);
  
  void SetValue(int i, int j, int k, int l, double value);
  bool IsPD(int i, int j, int k);
  void MakePD(int i, int j, int k);

};

inline irtkTensor irtkTensorField::operator()(int i, int j, int k)
{
  if ((i >= 0) && (i < _x) && (j >= 0) && (j < _y) && (k >= 0) && (k < _z)) 
  {
    return _tensorField[i][j][k];
  } 
  else 
  {
    cout << "irtkTensorField::operator(): parameter out of range\n";
    cout<<i<<" "<<j<<" "<<k<<endl;
    cout<<_x<<" "<<_y<<" "<<_z<<endl;
    exit(1);
  }
}

inline irtkTensor irtkTensorField::Get(int i, int j, int k)
{
  if ((i >= 0) && (i < _x) && (j >= 0) && (j < _y) && (k >= 0) && (k < _z)) 
  {
    return _tensorField[i][j][k];
  } 
  else 
  {
    cout << "irtkTensorField::operator(): parameter out of range\n";
    cout<<i<<" "<<j<<" "<<k<<endl;
    cout<<_x<<" "<<_y<<" "<<_z<<endl;
    exit(1);
  }
}

inline void irtkTensorField::SetValue(int i, int j, int k, int  l, double value)
{
  if ((i >= 0) && (i < _x) && (j >= 0) && (j < _y) && (k >= 0) && (k < _z) && (l >= 0) && (l < 7)) 
  {
    _tensorField[i][j][k](l)=value;
  } 
  else 
  {
    cout << "irtkTensorField::operator(): parameter out of range\n";
    cout<<i<<" "<<j<<" "<<k<<endl;
    cout<<_x<<" "<<_y<<" "<<_z<<endl;
    exit(1);
  }
}

inline bool irtkTensorField::InROI(double x, double y, double z)
{
  if ((x >= 0) && (x <= (_x-1)) && (y >= 0) && (y <= (_y-1)) && (z >= 0) && (z <= (_z-1))) 
    return true;
  else
    return false;
}

inline bool irtkTensorField::IsPD(int x, int y, int z)
{
  if ((x >= 0) && (x <= (_x-1)) && (y >= 0) && (y <= (_y-1)) && (z >= 0) && (z <= (_z-1))) 
    if(_tensorField[x][y][z].IsPD())
      return true;
    else
      return false;
  else
    return false;
}

inline void irtkTensorField::MakePD(int x, int y, int z)
{
  if ((x >= 0) && (x <= (_x-1)) && (y >= 0) && (y <= (_y-1)) && (z >= 0) && (z <= (_z-1))) 
    _tensorField[x][y][z].MakePD();
}

#endif