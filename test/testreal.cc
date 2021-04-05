
#include "NMR.h"
#include "rmatrix.h"
#include "timer.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

const double tip=M_PI/2;
const double phase=M_PI/2;

int main()
{
  spin_system ax(2);

  rmatrix Fx=real(F(ax,'x'));
  List<double> Fz=diag_Fz(ax);

  const int reptimes=getint("Repeat times? ");
  int i;

  timer<> stopwatch;

  cmatrix U;

  for (i=reptimes;i--;) U=Upulse(ax,NULL_NUCLEUS,tip,phase);
  cout << "Default: " << stopwatch()*1e6/reptimes << " us\n" << U << endl;

  stopwatch.reset();
  for (i=reptimes;i--;) U=Upulse(Fx,Fz,tip,phase);
  cout << "With Fx,z: " << stopwatch()*1e6/reptimes << " us\n" << U << endl;

  return 0;
}
