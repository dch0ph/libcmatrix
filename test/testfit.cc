#include "optim.h"
#include "cmatrix_utils.h"

using namespace std;
using namespace libcmatrix;

void testfunc(BaseList<double> dest, const BaseList<double>& paras,void* =NULL)
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

  cout << "Correct parameters: " << apars << endl;

  const int npts=10;
  List<double> data(npts);

  testfunc(data,apars);
  cout << "Initial data: " << data << endl;
  
  int i;
  for (i=npts;i--;) data(i)+=gauss(noise);
  cout << "Noisy data: " << data << endl;

  List<double> opars(npars);
  opars(0)=2.0; opars(1)=3.0; opars(2)=4.0;

  set_seed();
  for (i=3;i--;) 
    opars(i)+=random(0.5);
  List<double> tpars(opars);

  const List<size_t> actord(slice(0,npars));

  cout << "Initial parameters: " << tpars << endl;
  cout << "Fitting parameters: " << actord << endl;

  List<double> errs(npars,0.02);

  rmatrix covar;
  cout.precision(8);
  const double chisqr=fitdata(covar,tpars,FitFunction(testfunc),data,actord,errs,noise,2,1e-4);

  List<double> newd(npts);
  testfunc(newd,tpars);
  cout << "\nFitted data: " << newd << endl;
  cout << "\nFitted   \tError:\n";
  for (int j=0;j<covar.rows();j++)
    cout << tpars(j) << " \t" << sqrt(covar(j,j)) << endl;   

  cout << "\nExtended fitting model\n";
  List<Parameter> newpars;
  for (size_t i=0;i<npars;i++)
    newpars.push_back(Parameter(opars(i),errs(i)));
  SimpleBoundsState constraint(2.0,4.0); //limit parameter between 0 and 2
  cout << "Constraint: " << constraint << '\n';
  newpars(2U).set(2.5); //!< move within constraint
  newpars(2U).constrain(*(constraint.function())); //!< apply constrain function
  const double chisqr2=fitdata(covar,newpars,FitFunction(testfunc),data,noise,2,1e-4);

  cout << "\nFitted   \tError:\n";
  for (int j=0;j<covar.rows();j++)
    cout << newpars(j).get() << " \t" << sqrt(covar(j,j)) << endl;   

  return 0;
}
