// Test interpolation

#include "Lineshapes.h"
#include "simpsonio.h"
#include "ttyio.h"

using namespace libcmatrix;

#define USETYPE complex

int main(int argc, const char* argv[])
{
  int count=1;
  const size_t oldnp=100;
  
  LorentzianGaussian shapegen(50.0);
  List<USETYPE> data(oldnp);
  double sw=1000.0;
  LineshapeHistogram<USETYPE> hist(data,shapegen,sw);
  hist.add(1.0,0.0);
  
  const size_t newnp=getint(argc,argv,count,"New points? ",100);
  const bool dofold=getlogical(argc,argv,count,"Fold? ");
  const double oldstart=-sw/2;
  const double olddx=sw/oldnp;
  double start=oldstart;
  double dx=olddx;
  if (!dofold) {
    start=getfloat(argc,argv,count,"Start frequency? ",oldstart);
    dx=getfloat(argc,argv,count,"Frequency step? ",(olddx*oldnp)/newnp);
  }
  char fnamebase[256];
  char fname[256];
  getstring(argc,argv,count,"Output name? ",fnamebase,sizeof(fnamebase),"testinterp");
  
  sprintf(fname,"%s_orig.spe",fnamebase);
  write_simpson(fname,data,sw,true);

  List<USETYPE> newdata;
  if (dofold)
    newdata=cubic_interpolate(newnp,data);
  else
    newdata=cubic_interpolate(newnp,start,dx,data,oldstart,olddx);
  sprintf(fname,"%s.spe",fnamebase);
  write_simpson(fname,newdata,newnp*dx,true);

  return 0;
}
