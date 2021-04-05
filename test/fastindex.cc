#include <functional>
#include "spin_system.h"

using namespace libcmatrix;
using namespace std;

template<class T> List<size_t> makeindex(BaseList<int>& rindex,const BaseList<T>& a,size_t (T::*func)() const)
{
  const size_t max=rindex.length();
  rindex=-1;

  List<size_t> index(mxflag::temporary);
  index.create(max);

  size_t types=0;

  for (size_t n=a.length();n--;) {
    const size_t val=(a(n).*func)();
    if (val>max) throw Failed("makeindex: index out of range");
    if (rindex(val)<0) {
      rindex(val)=types;
      index(types)=val;
      types++;
    }
  }
  index.resize(types);
  return index;
}
  
  
int main()
{
  spin_system sys(4,"1H");
  sys.isotope(1,"13C");
  sys.isotope(2,"15N");

  List<int> rindex(MAX_NUCLEUS);
  List<size_t> index= makeindex(rindex,sys,&spin::nucleus);

  cout << rindex << endl;
  cout << index << endl;

  return 0;
}
