/* Test matrix exponential */

#include "cmatrix.h"
#include "rmatrix.h"
#include "timer.h"
#include "cmatrix_utils.h"
#include "ttyio.h"

using namespace libcmatrix;
using namespace std;

// Adapted from expokit padm.m 

int main(int argc, const char* argv[])
{
  int count=1;
  const int n=getint(argc,argv,count,"Matrix size? ",3);
  const int ntimes=getint(argc,argv,count,"Repetitions? ",1);
  const double scalef=1.0;//getfloat(argc,argv,count,"Scale factor? ",1.0);

  cmatrix A(n,n);
  for (size_t i=n;i--;) {
    for (size_t j=n;j--;) {
      if (i==j)
	A(i,i)=complex(0,scalef);
      else
	A(i,j)=random(1.0);
    }
  }

  if (n<10)
    cout << "Input matrix:\n" << A << '\n';

  cmatrix eB1,eB2,eB3;

  cout.precision(12);
  timer<> stopwatch;
  for (size_t r=ntimes;r--;)
    eB1=exp(A,1.0);
  cout << "exp(A) [via eigenvalues]  time=" << stopwatch() << " s\n";
  if (n<10)
    cout << eB1;
 
  stopwatch.reset();
  for (size_t r=ntimes;r--;) {
    eB2=A;
    exppade_ip(eB2,6U);
  }

  cout << "exp(A) [via Pade n=6]  time=" << stopwatch() << " s\n";
  if (n<10)
    cout << eB2;

  stopwatch.reset();
  for (size_t r=ntimes;r--;) {
    eB3=A;
    exppade_ip(eB3,5U);
  }
  
  cout << "exp(A) [via Pade n=5]  time=" << stopwatch() << " s\n";
  if (n<10)
    cout << eB3;
   
  return 0;
}
