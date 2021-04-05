#include "optim.h"
#include "cmatrix_utils.h"
#include "ttyio.h"

enum { AMP, ISO, J, LW };
const int NPARS=4;
const int FILE_MAXLEN=64;

/* fitdoublet - 200 0.01 400.0 200 100 20 testspec
 */

using namespace libcmatrix;
using namespace std;

int npts;
double noise,sw;

void add_signal(BaseList<double> &dest,double amp,double iso,double lw)
{
  int i;
  double a,lw2,x;

  /*  printf("%g  %g  %g\n",amp,iso,lw);*/

  lw2=lw*lw/4.0;
  for (i=npts;i--;) {
    x=((double)i/npts)-iso;
    a=lw2/(x*x+lw2);
    dest(i)+=a*amp;
  }
}

void trial_doublet(BaseList<double> &dest,const BaseList<double> &paras,void *)
{
  //int i;

  const double amp=paras(AMP);
  const double lw=paras(LW)/sw;
  const double iso=paras(ISO)/sw;
  const double j=paras(J)/sw;

  dest=0.0;
  //for (i=npts;i--;) dest(i)=0.0;
  add_signal(dest,0.5*amp,iso-j/2.0,lw);
  add_signal(dest,0.5*amp,iso+j/2.0,lw);
  // add_noise_ip(dest,noise);
}

main(int argc,const char *argv[])
{
  int count=1;
  char fname[FILE_MAXLEN];
  List<double> paras(NPARS);
  List<double> data;
  rmatrix covar;

  getstring(argc,argv,count,"Data file (<CR> for simple spectrum) ? ",fname,FILE_MAXLEN);
  int isfitting=(fname[0]!='\0');
  
  if (isfitting) {
    if (error_filter(read_vector(data,fname))) return(1);
    npts=data.length();
    cout << "Number of data points: " << npts << endl;
  }
  else
    npts=getint(argc,argv,count,"Number of points in spectrum? ",512);

  noise=getfloat(argc,argv,count,"Estimate of noise? ",0.0);
  sw=getfloat(argc,argv,count,"Spectral width (Hz) ? ");
  paras(ISO)=getfloat(argc,argv,count,"Centre (Hz) ? ");
  paras(J)=getfloat(argc,argv,count,"Splitting (Hz) ? ");
  paras(LW)=getfloat(argc,argv,count,"Linewidth (Hz) ? ");
  paras(AMP)=1.0;
  getstring(argc,argv,count,"Output file? ",fname,FILE_MAXLEN);

  if (isfitting) {
    List<double> errs(NPARS);

    const int varpars=4;
    List<size_t> actord(varpars);
    actord(0)=AMP;
    actord(1)=ISO;
    actord(2)=J;
    actord(3)=LW;
  
    double s=sum(data);
    printf("Integral: %g\n",s);

    paras(AMP)=2;

    errs(AMP)=0.5*s;
    errs(ISO)=0.1*sw;
    errs(J)=0.1*sw;
    errs(LW)=0.5*sw;

    double chisqr=fitdata(covar,paras,trial_doublet,(void *)NULL,data,actord,errs,noise,2);

    cout << "\nFitted   \tError:\n";
    for (int j=0;j<covar.rows();j++) {
      const int which=actord(j);
      cout << paras(which) << " \t" << sqrt(covar(j,j)) << endl;  
    } 
  }

  data.create(npts);
  trial_doublet(data,paras,(void *)NULL);
  if (!isfitting && noise)
    for (int i=npts;i--;) data(i)+=gauss(noise);

  return error_filter(write_vector(fname,data));
}
