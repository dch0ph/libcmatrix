#include "ListList.h"
#include <cstdarg>

using namespace libcmatrix;
using namespace std;

int main()
{
  size_t rsizes[]={1,3,3,1};
  double rdata[]={ 1.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -1.5};
  const BaseList<size_t> str(4,rsizes);
  BaseList<double> data(8,rdata);
  cout << ExplicitList<3,double>(1.0,1.5,2.0) << endl;

  const int ridata[]= { 1, 1, 1, 1, -1, -1, -1, -1};
  BaseList<const int> idata(8,ridata);

  //data+=idata;
  
  ListList<double> A(str,data);
  cout << "A: " << A << endl;

  A+=1.0;
  cout << "After add 1.0: " << A << endl;

  const ListList<int> B(str,idata);
  cout << "B: " << B << endl;

  A+=B;
  cout << "A+=B: " << A << endl;

  ListList<double> C(2,10);//re-allocation forced on 2nd push
  C.push_back(data);
  C.push_back(idata);
  cout << "C: " << C << endl;

  return 0;
}
