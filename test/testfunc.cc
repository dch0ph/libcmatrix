/* Test table/function objects */

#include "FunctionObject.h"

using namespace std;
using namespace libcmatrix;

  
int main()
{
  double durs[]={ 1.0, 2.0, 3.0, 4.0 };

  SequentialIndex myind(BaseList<double>(4,durs));
  SequentialCyclicIndex myindc(BaseList<double>(4,durs));

  TableIndex tindex(2,10.0);
  CyclicIndex cindex(2,10.0);

  cout << "Cyclic\n";
  for (int i=0;i<20;i++)
    cout << myindc(i+0.5) << "  " << cindex(i+0.5) << "\n";

  cout << "Table\n";
  for (int i=20;i--;) {
    double t=(i/2.0)+0.25;
    cout << myind(t) << "  " << tindex(t) << "\n";
  }

  return 0;
}
