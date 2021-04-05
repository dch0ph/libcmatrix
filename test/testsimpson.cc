#include "simpsonio.h"
#include "NMR.h"

using namespace std;
using namespace libcmatrix;

int main()
{
  Matrix<complex> a;
  read_simpson(a,"simptest.spe");
  cout << a << endl;
  write_simpson("simptestn.spe",a,10e3,0.0);
  
  spin_system sys(3,"1H");
  sys(0).isotope("13C");
  
  cout << "+1 coherence (all spins)\n" << coherencematrix(sys,1) << endl;
  cout << "-1,+1 coherence (all spins)\n" << coherencematrix(sys,ExplicitList<2,int>(1,-1)) << endl;
  const Matrix<bool> mask2=coherencematrix(sys,"1H",1);
  cout << "+1 coherence (1H)\n" << mask2 << endl;
  const Matrix<bool> mask=coherencematrix(sys,"13C",ExplicitList<2,int>(1,-1));
  cout << "-1,+1 coherence (13C)\n" << mask << endl;

  cout << "Or masks\n" << (mask | mask2) << endl;
  
  const int n=sys.size();
  rmatrix H(n,n,4.0);
  H.emultiply(mask);
  cout << "After masking\n" << H << endl;
  return 0;
}
