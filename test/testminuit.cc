#include <iostream>
#include "cmatrix_complex.h"
#include "List.h"
#include "optim.h"
#include "ttyio.h"
#include <vector>
#include "Minuit/VariableMetricMinimizer.h"
#include "Minuit/FunctionMinimum.h"

using namespace std;
using namespace libcmatrix;

// testminuit <PASS pitch> <output file>  < testminuit.in

const double TWOPI=2.0*M_PI;

class PassValues : public BaseMinFunction {
public:
  PassValues(double pitch_) : pitch(pitch_) {}
  double operator()(const BaseList<double>&) const;
private:
  double pitch;
  complex passval(const BaseList<double>& vals, double thetan, int m) const;
};

int main(int argc,const char *argv[])
{
  int count=1;
  const double pitch=getfloat(argc,argv,count,"Pitch? ",0.0);

  //char fname[64];
  //getstring(argc,argv,count,"Output file? ",fname,64);

  vector<double> init_par;
  init_par.push_back(0.2);
  init_par.push_back(0.4);
  init_par.push_back(0.6);
  init_par.push_back(0.8);

  vector<double> init_err(4,0.05);
  
  VariableMetricMinimizer theMinimizer;
  PassValues theFCN(pitch);
  const FunctionMinimum min(theMinimizer.minimize(MinuitAdaptor<PassValues>(theFCN),init_par,init_err));
  //  cout << "Minimum: " << min << '\n';
  if (!min.isValid()) {
    std::cerr << "Minuit failed to converge\n";
    return 1;
  }
  const MnUserParameterState& parstate(min.userState());
  cout << "Optimised parameters:\n";
  for (size_t i=0;i<4;i++) 
    cout << i << ": " << parstate.value(i) << " +/- " << parstate.error(i) << '\n';

  //fpout=fopen(fname,"a");

  //minuitinit("testminuit.in");
  //minuit_(minuitlink,passfunc);

  //if (fpout) fclose(fpout);

  //f_exit(); // cleanup FORTRAN I/O
//     List<double> avals(n);
//     int j;
//     for (j=n-1;j--;) avals(j)=vals(j);
//     avals(n-1)=thetan;
//     avals.sort();
//     fprintf(fpout,"%.5f",pitch);
//     for (j=0;j<n;j++)
//       fprintf(fpout,"  %.5f",avals(j));
//     fprintf(fpout,"\n");

  return 0;
}

inline int powm1(int n) { return (n & 1) ? -1 : 1; }

complex PassValues::passval(const BaseList<double>& vals, double thetan, int m) const
{
  complex sum(0,0);
  const int n=vals.size()+1;
  for (int i=n;i--;) {
    const double p=TWOPI*((i==n-1) ? thetan : vals(i));
    if (i & 1)
      sum+=expi(m*p);
    else
      sum-=expi(m*p);
  }
  return 2.0*sum+1.0-powm1(n)*expi(TWOPI*m*pitch);
}

double PassValues::operator()(const BaseList<double>& vals) const
{
  
  const int n=vals.size()+1;

  const int sig=powm1(n);
  
  double sumfirst=0.0;

  for (int i=n-1;i--;)
    sumfirst+=powm1(i+1)*vals(i);
  const double thetan=-sig*sumfirst+0.5;
  
  const complex m1=passval(vals,thetan,1);
  const complex m2=passval(vals,thetan,2);

  //cout << "m1: " << m1 << " \t m2: " << m2 << endl;

  return sqrt(norm(m1)+norm(m2));
}
