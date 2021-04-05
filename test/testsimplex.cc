#include "optim.h"
#include "ttyio.h"
#include "cmatrix_utils.h"

using namespace std;
using namespace libcmatrix;

void testfunc(BaseList<double> dest, const BaseList<double>& paras, void* =NULL)
{
  const int npars=dest.length();
  const double a0=paras(2);
  const double a1=paras(1);
  const double a2=paras(0);

  for (int n=npars;n--;) {
    double x=n+1;
    dest(n)=a2*x*x+a1*x+a0;
  }
}

int main()
{
  const double noise=0.02;
  const int npars=3;

  List<double> apars(3);
  apars(0)=1.0;
  apars(1)=0.0;
  apars(2)=1.0;

  const int npts=10;
  List<double> data(npts);

  testfunc(data,apars);
  cout << "Initial: " << data << endl;
  
  int i;
  for (i=npts;i--;) data(i)+=gauss(noise);
  cout << "Noisy: " << data << endl;

  List<double> tpars(npars);
  tpars(0)=4.0; tpars(1)=1.0; tpars(2)=1.0;

  const int spar=0;
  const int epar=npars-1;
  List<size_t> actord(epar-spar+1);

  set_seed();
  for (i=spar;i<=epar;i++) {
    tpars(i)+=random(0.5);
    actord(i-spar)=i;
  }

  cout << "Initial parameters: " << tpars << endl;

  List<double> errs(npars,0.5);

  rmatrix covar;

  do {
    (void)simplex_fitdata(covar,tpars,FitFunction(testfunc),data,actord,errs,noise,3,1e-4);
  } while(getlogical("Continue? ",1));

  cout << "\nFitted   \tError:\n";
  for (int j=0;j<covar.rows();j++) {
    const int k=actord(j);
    cout << tpars(k) << " \t" << sqrt(covar(k,k)) << endl;   
  }

  List<double> newd(npts);
  testfunc(newd,tpars);
  cout << newd << endl;

  return 0;
}
