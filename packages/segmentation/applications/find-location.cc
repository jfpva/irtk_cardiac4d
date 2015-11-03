#include <irtkImage.h>
#include <irtkReconstruction.h>
#include <irtkReconstructionb0.h>
#include <irtkRegistration.h>
#include <irtkImageAffineRegistrationWithPadding.h>
irtkRealImage target, source;
bool ok, _swap=true;

void usage()
{
  cerr<<"find-location [target] [source] [output_transformation] <-x> <-y>"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 3)
    usage();
    
    target.Read(argv[1]);
    argc--;
    argv++;

    source.Read(argv[1]);
    argc--;
    argv++;

    char *output_name = argv[1];
    argc--;
    argv++;
    
    // Parse options.
    while (argc > 1){
      ok = false;
    
      if ((ok == false) && (strcmp(argv[1], "-x") == 0)){
      argc--;
      argv++;
      _swap = false;
      ok = true;
    }

      if ((ok == false) && (strcmp(argv[1], "-y") == 0)){
      argc--;
      argv++;
      _swap=true;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
   
  irtkReconstructionb0 reconstruction;
  
  double x,y,z;
  double x0,y0,z0,x1,y1,z1;
  double xt,yt,zt,xs,ys,zs;
  double ptx,pty,ptz,psx,psy,psz;

  
  irtkAffineTransformation tr;
   
  tr.PutStatus(TX,  _Active);
  tr.PutStatus(TY,  _Active);
  tr.PutStatus(TZ,  _Active);
  tr.PutStatus(RX,  _Passive);
  tr.PutStatus(RY,  _Passive);
  tr.PutStatus(RZ,  _Passive);
  tr.PutStatus(SX,  _Passive);
  tr.PutStatus(SY,  _Passive);
  tr.PutStatus(SZ,  _Passive);
  tr.PutStatus(SXY, _Passive);
  tr.PutStatus(SXZ, _Passive);
  tr.PutStatus(SYZ, _Passive);

  
  irtkImageAffineRegistrationWithPadding registration;
  irtkGreyImage t = target, s=source;
  t.Write("t.nii.gz");
  s.Write("s.nii.gz");
  registration.SetInput(&t,&s);
  registration.SetOutput(&tr);
  registration.GuessParameter();
  registration.Run();
  
  tr.irtkTransformation::Write("translation.dof");
  
  cout<<"tx = "<<tr.GetTranslationX()<<" ty = "<<tr.GetTranslationY()<<" tz = "<<tr.GetTranslationZ()<<endl;
  

  
  x=0;y=0;z=0;
  tr.Transform(x,y,z);
  
  cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<endl; //x,y,z is not the same as x,ty,tz;
  
  //find PE for target 
  x0=0;y0=0;z0=0;
  x1=0;y1=1;z1=0;
  
  target.ImageToWorld(x0,y0,z0);
  target.ImageToWorld(x1,y1,z1);
  
  xt=x1-x0;
  yt=y1-y0;
  zt=z1-z0;
  
  cout<<"PE direction for target is "<<xt<<" "<<yt<<" "<<zt<<" "<<endl;
  
  x0=0;y0=0;z0=0;
  x1=0;y1=1;z1=0;
  
  source.ImageToWorld(x0,y0,z0);
  source.ImageToWorld(x1,y1,z1);
  
  xs=x1-x0;
  ys=y1-y0;
  zs=z1-z0;
  
  cout<<"PE direction for source is "<<xs<<" "<<ys<<" "<<zs<<" "<<endl;
  
  irtkImageAttributes attr=target.GetImageAttributes();
  cout<<attr._dy<<endl;
  double dy = attr._dy;

  x0=0;y0=0;z0=0;
  x1=x;y1=y;z1=z;
  
  target.WorldToImage(x0,y0,z0);
  target.WorldToImage(x1,y1,z1);
  
  
  ptx=(x1-x0);
  pty=(y1-y0);
  ptz=(z1-z0);
  
  cout<<"translation in image coord of target is "<<ptx/dy<<" "<<pty/dy<<" "<<ptz/dy<<" "<<endl;
 
  x0=0;y0=0;z0=0;
  x1=x;y1=y;z1=z;
  
  source.WorldToImage(x0,y0,z0);
  source.WorldToImage(x1,y1,z1);
  
  psx=(x1-x0);
  psy=(y1-y0);
  psz=(z1-z0);
  
  cout<<"translation in image coord of source is "<<psx/dy<<" "<<psy/dy<<" "<<psz/dy<<" "<<endl;
  exit(1);
  double dist_transl = (psy-pty)/2;
  
  cout<<"distortion translation is: "<<dist_transl/dy<<endl;
  
  ptx=0;ptz=0;pty=dist_transl;
  target.ImageToWorld(ptx,pty,ptz);
  x0=0;y0=0;z0=0;
  target.ImageToWorld(x0,y0,z0);
  irtkAffineTransformation trt;
  trt.PutTranslationX(ptx-x0);
  trt.PutTranslationY(pty-y0);
  trt.PutTranslationZ(ptz-z0);
  trt.irtkTransformation::Write("trt.dof.gz");
  
  psx=0;psz=0;psy=dist_transl;
  source.ImageToWorld(psx,psy,psz);
  x0=0;y0=0;z0=0;
  source.ImageToWorld(x0,y0,z0);
  irtkAffineTransformation trs;
  trs.PutTranslationX(psx-x0);
  trs.PutTranslationY(psy-y0);
  trs.PutTranslationZ(psz-z0);
  trs.irtkTransformation::Write("trs.dof.gz");
  
  //irtkAffineTransformation distortion();
  
  irtkRigidTransformation motion_t, motion_s;
  irtkRealImage cor_target(target), cor_source(source);
  
  irtkImageTransformation imagetransformation;
  irtkLinearInterpolateImageFunction interpolator;
  imagetransformation.PutInterpolator(&interpolator);
  imagetransformation.PutTargetPaddingValue(-1);
  imagetransformation.PutSourcePaddingValue(0);
  
  imagetransformation.SetInput(&target, &trt);
  imagetransformation.SetOutput(&cor_target);
  imagetransformation.Run();

  imagetransformation.SetInput(&source, &trs);
  imagetransformation.SetOutput(&cor_source);
  imagetransformation.Run();
  
  cor_target.Write("cor_target.nii.gz");
  cor_source.Write("cor_source.nii.gz");
  
  t=cor_target;
  s=cor_source;
  
  irtkImageRigidRegistrationWithPadding rigidreg;
  
  rigidreg.SetInput(&t,&s);
  rigidreg.SetOutput(&motion_s);
  rigidreg.GuessParameter();
  rigidreg.Run();
  
  motion_s.irtkTransformation::Write("motion_s.dof");

  rigidreg.SetInput(&s,&t);
  rigidreg.SetOutput(&motion_t);
  rigidreg.GuessParameter();
  rigidreg.Run();
  
  motion_t.irtkTransformation::Write("motion_t.dof");
  
  irtkRealImage aligned_target(cor_source), aligned_source(cor_target);

  imagetransformation.SetInput(&cor_source, &motion_s);
  imagetransformation.SetOutput(&aligned_source);
  imagetransformation.Run();
  
  imagetransformation.SetInput(&cor_target, &motion_t);
  imagetransformation.SetOutput(&aligned_target);
  imagetransformation.Run();
  
  aligned_target.Write("aligned_target.nii.gz");
  aligned_source.Write("aligned_source.nii.gz");
  
  irtkAffineTransformation shim_t, shim_s;
  
  reconstruction.ShimDistortion(aligned_source,cor_target,shim_t,true);
  shim_t.irtkTransformation::Write("shim_t.dof");
  
  reconstruction.ShimDistortion(aligned_target,cor_source,shim_s,true);
  shim_s.irtkTransformation::Write("shim_s.dof");
  
  
  
  
  /////////////////////////////////////
  //ROTATION
  /////////////////////////////////////
  /*
  irtkRigidTransformation offset;
  reconstruction.ResetOrigin(t,offset);
  
  offset.irtkTransformation::Write("offset.dof");
  
  irtkMatrix mo = offset.GetMatrix();
  irtkMatrix m = tr.GetMatrix();
  m=m*mo;
  tr.PutMatrix(m);
  
  tr.PutStatus(RX,  _Active);
  tr.PutStatus(RY,  _Active);
  tr.PutStatus(RZ,  _Active);

  t.Write("t.nii.gz");
  s.Write("s.nii.gz");
  registration.SetInput(&t,&s);
  registration.SetOutput(&tr);
  registration.GuessParameter();
  registration.Run();
  tr.irtkTransformation::Write("rot.dof.gz");

  mo.Invert();
  m = tr.GetMatrix();
  m=m*mo;
  tr.PutMatrix(m);
  tr.irtkTransformation::Write(output_name);
  
  x=0;y=0;z=0;
  tr.Transform(x,y,z);
  
  cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<endl; //x,y,z is not the same as x,ty,tz;
 
  
  x0=0;y0=0;z0=0;
  x1=x;y1=y;z1=z;
  
  target.WorldToImage(x0,y0,z0);
  target.WorldToImage(x1,y1,z1);
  
  
  ptx=(x1-x0);
  pty=(y1-y0);
  ptz=(z1-z0);
  
  cout<<"translation in image coord of target is "<<ptx/dy<<" "<<pty/dy<<" "<<ptz/dy<<" "<<endl;
 
  x0=0;y0=0;z0=0;
  x1=x;y1=y;z1=z;
  
  source.WorldToImage(x0,y0,z0);
  source.WorldToImage(x1,y1,z1);
  
  psx=(x1-x0);
  psy=(y1-y0);
  psz=(z1-z0);
  
  cout<<"translation in image coord of source is "<<psx/dy<<" "<<psy/dy<<" "<<psz/dy<<" "<<endl;
  
  ptx=0;ptz=0;pty=-pty;
  target.ImageToWorld(ptx,pty,ptz);
  x0=0;y0=0;z0=0;
  target.ImageToWorld(x0,y0,z0);
  trt.PutTranslationX(ptx-x0);
  trt.PutTranslationY(pty-y0);
  trt.PutTranslationZ(ptz-z0);
  trt.irtkTransformation::Write("trt.dof.gz");
  
  psx=0;psz=0;
  source.ImageToWorld(psx,psy,psz);
  x0=0;y0=0;z0=0;
  source.ImageToWorld(x0,y0,z0);
  trs.PutTranslationX(psx-x0);
  trs.PutTranslationY(psy-y0);
  trs.PutTranslationZ(psz-z0);
  trs.irtkTransformation::Write("trs.dof.gz");
  
  //test
  //tr.PutRotationX(0);
  //tr.PutRotationY(0);
  //tr.PutRotationZ(0);
  //tr.irtkTransformation::Write("tr-translation.dof.gz");
  irtkAffineTransformation rot(tr);
  rot.PutTranslationX(0);
  rot.PutTranslationY(0);
  rot.PutTranslationZ(0);
  rot.irtkTransformation::Write("tr-rotation.dof.gz");
  
  m = rot.GetMatrix();
  mo = offset.GetMatrix();
  irtkMatrix mi(mo);
  mi.Invert();
  
  m=mo*m*mi;
  rot.PutMatrix(m);
  rot.irtkTransformation::Write("tr-rotation-off.dof.gz");
  
  irtkAffineTransformation tra;
  tra.PutTranslationX(tr.GetTranslationX()-rot.GetTranslationX());
  tra.PutTranslationY(tr.GetTranslationY()-rot.GetTranslationY());
  tra.PutTranslationZ(tr.GetTranslationZ()-rot.GetTranslationZ());
  tra.irtkTransformation::Write("tr-tra-off.dof.gz");
  
  x=0;y=0;z=0;
  tra.Transform(x,y,z);
  
  cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<endl; //x,y,z is not the same as x,ty,tz;
 
  
  x0=0;y0=0;z0=0;
  x1=x;y1=y;z1=z;
  
  target.WorldToImage(x0,y0,z0);
  target.WorldToImage(x1,y1,z1);
  
  
  ptx=(x1-x0);
  pty=(y1-y0);
  ptz=(z1-z0);
  
  cout<<"translation in image coord of target is "<<ptx/dy<<" "<<pty/dy<<" "<<ptz/dy<<" "<<endl;
 
  x0=0;y0=0;z0=0;
  x1=x;y1=y;z1=z;
  
  source.WorldToImage(x0,y0,z0);
  source.WorldToImage(x1,y1,z1);
  
  psx=(x1-x0);
  psy=(y1-y0);
  psz=(z1-z0);
  
  cout<<"translation in image coord of source is "<<psx/dy<<" "<<psy/dy<<" "<<psz/dy<<" "<<endl;
  
  */

  
  //todo: normalize both vectors and express translation through these two directions plus the one orthogonal to both?
  //result should be same scale times each vector to get from common location to the position in distorted space

}
