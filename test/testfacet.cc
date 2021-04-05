/* test use of iostream manipulators */

#include <iostream>
#include <iomanip>
#include <map>
#include "cmatrix.h"

using namespace libcmatrix;

struct tcomplex {
  tcomplex(double rv, double iv) 
    : r_(rv), i_(iv) {}
  double r_,i_;
};

enum lcmio_t { LCM_PREC=1 };
typedef std::map<int,int> lcmiomap_t;
lcmiomap_t lcmiomap;

struct lcmioword { 
  lcmioword(int whichv,int precv) : which_(whichv), prec_(precv) {}
  int which_;
  int prec_;
};

//inline lcmioword_ lcmio(lcmio_t which,int prec) { return lcmioword_(which,prec); }

std::ostream& operator<< (std::ostream& ostr, const lcmioword& w)
{
  const lcmiomap_t::iterator iter(lcmiomap.find(w.which_));
  int allocind;
  if (iter==lcmiomap.end()) {    
    allocind=ostr.xalloc();
    lcmiomap[w.which_]=allocind;
  }
  else
    allocind=iter->second;
  ostr.iword(allocind)=w.prec_;
  return ostr;
}

int getlcmword(const std::ostream& ostr, lcmio_t which)
{
  const lcmiomap_t::iterator iter(lcmiomap.find(which));
  if (iter==lcmiomap.end())
    throw std::bad_cast();
  return const_cast<std::ostream& >(ostr).iword(iter->second);
}

int getlcmprec(const std::ostream& ostr)
{
  const lcmiomap_t::iterator iter(lcmiomap.find(LCM_PREC));
  return (iter==lcmiomap.end()) ? ostr.precision() : const_cast<std::ostream& >(ostr).iword(iter->second);
}

std::ostream& operator<< (std::ostream& ostr, const tcomplex& z)
{
  const int width=getlcmprec(ostr);
  ostr << width << ' ';
  return ostr << '(' << std::setw(width) << z.r_ << ',' << z.i_ << ')';
}

int main()
{
  cmatrix a(2,2,complex(2.333322,3.141345));
  std::cout << a << '\n';
  std::cout << complex(1.234567890) << '\n';
  std::cout << complex(0,1.234567890) << '\n';
  std::cout << setmatrixprecision(0) << a << '\n';

  std::cout << setmatrixprecision(2) << a << '\n';
  cmatrix_ostream_controller(std::cout).complexview=ostream_controller::withi;
  std::cerr << a << '\n';
  std::cout << complex(1.234567890) << '\n';
  std::cout << complex(0,1.234567890) << '\n';
  //  std::cerr << lcmioword(LCM_PREC,3) << a << '\n';

  return 0;
}
