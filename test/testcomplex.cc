
#include <cstdlib>
#include <iostream>
#include "timer.h"

//using namespace std;

#ifdef NEWSTYLE
#include <complex>
#define COMPLEX complex<double>

//const complex<float> junk;
#else
#include "cmatrix_complex.h"
#define COMPLEX complex

// ostream& operator<< (ostream& ostr, const complex &z)
// {
//   return ostr << '(' << real(z) << "," << imag(z) << ')';
// }

#endif

using namespace libcmatrix;

int main()
{
#ifdef NEWSTYLE
  std::cout << "New style complex\n";
#else
  std::cout << "Old style complex\n";
#endif

  const COMPLEX a(2.0,3.0);
  const COMPLEX b(0.7071,0.7071);
  COMPLEX c1,c2,c3,c4,c5;

  const int times=200000000;

  timer<WallTimer> stopwatch;
  int i;

  c1=c2=c3=c4=c5=COMPLEX(0,0);

  stopwatch.reset();
  for (i=times/5;i--;) {
    c1=c2+b;
    c2=c3+b;
    c3=c4+b;
    c4=c5+b;
    c5=c1+b;
  }
  std::cout << "Simple add: " << c1 << std::endl;
  std::cout << "Time taken per operation " << (1e6*stopwatch()/times) << " us\n\n";

  stopwatch.reset();

  c1=c2=c3=c4=c5=a;

  for (i=times/5;i--;) {
    c1=b*c2;
    c2=b*c3;
    c3=b*c4;
    c4=b*c5;
    c5=b*c1;
  }
  std::cout << "multiply: " << c1 << std::endl;
  std::cout << "Time taken per operation " << (1e6*stopwatch()/times) << " us\n\n";

//   c1=c2=c3=c4=c5=COMPLEX(0,0);

//   stopwatch.reset();
//   for (i=times/5;i--;) {
//     mla(c1,a,b);
//     mla(c2,a,b);
//     mla(c3,a,b);
//     mla(c4,a,b);
//     mla(c5,a,b);
//   }
//   std::cout << "mla: " << c1 << std::endl;
//   std::cout << "Time taken per operation " << (1e6*stopwatch()/times) << " us\n\n";

  c1=c2=c3=c4=c5=COMPLEX(0,0);
  stopwatch.reset();
  for (i=times/5;i--;) {
    c1=c2;
    c1+=a;
    c2=c3;
    c2+=a;
    c3=c4;
    c3+=a;
    c4=c5;
    c4+=a;
    c5=c1;
    c5+=a;
  }
  std::cout << "copy add: " << c1 << std::endl;
  std::cout << "Time taken per operation " << (1e6*stopwatch()/times) << " us\n\n";

  const complex complexi(0,1);
  const complex sqrti=sqrt(complexi);

  std::cout << "sqrt -1: " << sqrt(complex(-1)) << std::endl;  
  std::cout << "sqrt i: " << sqrti << "   (sqrt i)^2: " << (sqrti*sqrti) << std::endl;
  
  const double ang=M_PI/3;
  std::cout << "exp(i pi/3): " << exp(complexi*ang) << "   expi(pi/3): " << expi(ang) << std::endl;
  std::cout << "cos pi/3, sin pi/3: " << cos(ang) << "," << sin(ang) << std::endl;

  const double mynan=sqrt(-1);
  const double myinf=sqrt(mynan);
  std::cout << "NAN: " << mynan << "   Infty: " << myinf << std::endl;
  std::cout << "sqrt(NAN): " << sqrt(complex(mynan,0)) << "   sqrt(complex(NAN,NAN)): " << sqrt(complex(mynan,mynan)) << std::endl;
  std::cout << "log(complex(0)): " << log(complex(0)) << std::endl;
  std::cout << "log(i): " << log(complexi) << "  exp(log(i)): " << exp(log(complexi)) << std::endl;
  std::cout << "log(complex(myinf)): " << log(complex(myinf)) << std::endl;
  
  return 0;
}

