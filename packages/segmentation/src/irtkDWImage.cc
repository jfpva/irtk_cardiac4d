#include <irtkDWImage.h>

irtkDWImage::irtkDWImage(const irtkRealImage& image):irtkRealImage(image)
{  
}

void irtkDWImage::SetDirections(irtkMatrix directions)
{
  //check there is right number of directions
  if (directions.Cols() != (GetT()-1))
  {
    cout<<"Please give "<<GetT()-1<<"directions."<<endl;
    exit(1);
  }

  //check directions have 3 codinates
  if (directions.Rows() != 3)
  {
    cout<<"Please give 3 coordinates of direction vectors."<<endl;
    exit(1);
  }
  
  //make sure the direction vectors have unit lengths
  double norm;
  for (int i=0; i< directions.Cols();i++)
  {
    norm=0;
    for (int j=0; j<3; j++)
      norm += directions(j,i)*directions(j,i);
    norm = sqrt(norm);
    for (int j=0; j<3; j++)
      directions(j,i)/=norm;
  }
  
  _directions = directions;

}

void irtkDWImage::SetDirection(int index, irtkVector direction)
{
  //check the index is witin the bounds
  if (index >= (GetT()-1))
  {
    cout<<"Direction index "<<index<<"out of range."<<endl;
    exit(1);
  }

  //check the directions have been set
  if (_directions.Cols() != (GetT()-1))
  {
    cout<<"Direction have not been set yet."<<endl;
    exit(1);
  }

  //check directions have 3 coordinates
  if (direction.Rows() != 3)
  {
    cout<<"Please give 3 coordinates of direction vectors."<<endl;
    exit(1);
  }
  
  //make sure the direction vectors have unit lengths
  double norm=0;
  for (int j=0; j<3; j++)
    norm += direction(j)*direction(j);
  norm = sqrt(norm);
  for (int j=0; j<3; j++)
    direction(j)/=norm;
  
  for (int j=0; j<3; j++)
   _directions(j,index) = direction(j);

}

irtkTensor irtkDWImage::CalculateTensor(int x,int y, int z)
{
  int nDir = GetNumberOfDirections();
  irtkVector C(nDir);
  irtkMatrix b(nDir,6);
  irtkTensor D;
  int i,j,k;
  
  //tensor does not exist, return zero tensor
  for (i=0;i<=nDir;i++)
    if (this->GetAsDouble(x,y,z,i)<=0.0001)
    {
      //D.MakePD();
      return D;
    }

  for (i=0;i<nDir;i++)
  {
    C(i)=-log(this->GetAsDouble(x,y,z,i+1)/this->GetAsDouble(x,y,z,0))/_b;
    for (j=0;j<3;j++)
    {
      b(i,j)=_directions(j,i)*_directions(j,i);
      for(k=j+1;k<3;k++)
	b(i,j+k+2)=2*_directions(j,i)*_directions(k,i);
    }
  }
  
  irtkMatrix bt = b;
  bt.Transpose();
  
  irtkMatrix btb = bt*b;
  btb.Invert();
  
  irtkVector d = bt*C;
  d=btb*d;
  D=d;
  if (!D.IsPD())
    D.MakePD();
  return D;
}

