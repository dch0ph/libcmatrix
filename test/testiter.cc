/* Try using Standard C++ Library iteration with List's */

#include <algorithm>
#include <numeric>
#include "cmatrix_complex.h"
#include "cmatrix.h"

char errmesg[256];

using namespace std;
using namespace libcmatrix;

template<typename T> T times2(T a) { return 2*a; }

template<typename T> void transform_ip(T* begin,T* end,T (*func)(T))
{
  while (begin!=end) {
    *begin=func(*begin);
    begin++;
  }
}

template<typename T> void pluseq(T& a,const T& b) { a+=b; }

template<typename T> void transform_ip2(T* begin,T* end,T& (T::*pmf)())
{
  while (begin!=end) ((begin++)->*pmf)();
}

template<class Iter,class T> void transform_ip3(Iter begin,Iter end,void (*pf)(T&,const T&),T arg)
{
  while (begin!=end) pf(*begin++,arg);
}

int main()
{
  int items[]={1,9,5,6,8,9,9,9,5};
  
  try {

  BaseList<int> mylist(8,items);
  cout << "Before: " << mylist << endl;

  replace(mylist.begin(),mylist.end(),5,10);
  cout << "After replace 5's by 10's: " << mylist << endl;

  Matrix<int> mymat(3,3,items);
  cout << "Before:\n" << mymat << endl;
  
  replace(mymat.begin(),mymat.end(),9,-9);
  cout << "After replacing 9 by -9:\n" << mymat << endl;

  //sort(mymat.begin(),mymat.end());
  //cout << "After sorting:\n" << mymat << endl;

  //transform_ip(mymat.begin(),mymat.end(),times2<int>);
  //cout << "After times 2:\n" << mymat << endl;

  //transform_ip2(mymat.begin(),mymat.end(),&complex::conj);
  //cout << "After conj:\n" << mymat << endl;

  //  transform_ip3(mymat.begin(),mymat.end(),&pluseq,2);
  //cout << "After +=2:\n" << mymat << endl;

  const List<int> c_mylist(mylist);

  //SliceList<int> mysel=mylist(slice(1,3,2));
  //const SliceList<int> c_mysel=c_mylist(slice(1,3,2));
  size_t indices[]={1,2,4};
  const BaseList<size_t> selinds(3,indices);

  IndirectList<int> mysel(mylist(selinds));
  const IndirectList<int> c_mysel=c_mylist(selinds);
  cout << "Selected: " << mysel << endl;

  mysel+=List<int>(mysel);
  cout << "After adding selection to self: " << mylist << endl;

  cout << "Sum: " << accumulate(mysel.begin(),mysel.end(),0) << endl;
  
  replace(mysel.begin(),mysel.end(),18,14);
  cout << "After replace 18 by 14: " << mylist << endl;

  //transform_ip3(mysel.begin(),mysel.end(),&pluseq,int(2));
  //cout << "After +=2 on selected: " << mylist << endl;
  
  mysel=0;
  cout << "=0: " << mylist << endl;

  mymat(1,range())=1;
  cout << "After set row 1 to 1:\n" << mymat << endl;

  //const Matrix<int> c_mymat=mymat;
  
  // mymat(slice(),1)=mymat(slice(),2);
  //cout << "After set col 1:\n" << mymat << endl;

  cout << "(slice(0,2),1): " << mymat(slice(0,2),1) << endl;
  cout << "Row 2: " << mymat(2,range()) << endl;
  cout << "Column 2: " << mymat(range(),2) << endl;

  cout << "Rows 0-1, cols 0,2:\n";
  IndirectMatrix<int,range,slice> matsel=mymat(range(0,1),slice(0,2,2));

  cout << matsel << endl;

  cout << "diagonal block 1-2:\n";
  //  size_t inds[]={1,2};
  const ExplicitList<2,size_t> myinds(1,2);
  cout << mymat(myinds,myinds) << endl;

  cout << "(1,{1,2}): " << mymat(1,myinds) << endl;
  cout << "({1,2},1): " << mymat(myinds,1) << endl;

  matsel+=4;
  cout << "After +=4\n" << mymat << endl;

  add_ip(matsel,-4);
  cout << "After -=4\n" << mymat << endl;

  Matrix<int> res;
  add(res,matsel,mymat(myinds,myinds));
  cout << "After add together\n" << res << endl;

  matsel=53;
  cout << "After set 53\n" << mymat << endl;

#ifdef LCM_COMPLEX_CHEAT
  List<complex> clist(4,5.0);
  imags(clist)=1.0;
  cout << "After set imag to 1: " << clist << endl;
#endif
 
  }
  catch (MatrixException& exc) {
    cerr << exc << endl;
  }
  return 0;
}
 
