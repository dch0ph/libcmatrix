#include "BaseMetaPropagation.h"
//#include "spin_system.h"
#include "cmatrix.h"
#include "MultiMatrix.h"
#include <numeric>

using namespace libcmatrix;

// void dumppattern(const SpinOpGenerator& opgen, const char* nuc, char op)
// {
//   const operator_spec spec(nuc,op);
//   const block_pattern blkspec(opgen,spec);
//   std::cout << blkspec << '\n';
//   Matrix<bool> present;
//   blkspec.makepresent(present);
//   std::cout << present << '\n';
// }


//template<class T> inline LCM_VAL(T) product(const T& a) { return std::accumulate(a.begin(),a.end(),LCM_VAL(T)(1),std::multiplies< LCM_VAL(T) >()); }

void dumppattern(const blocking_information& blkinfo, size_t which, int coher)
{
  std::cout << "Index: " << which << "  Coherence: " << coher << '\n';
  block_pattern blkspec(blkinfo,which,coher,false);
  std::cout << blkspec;
  block_pattern::iterator iter(blkspec);
  Matrix<bool> present(blkinfo.totallevels,blkinfo.totallevels,false);
  size_t r,c;
  while (iter.next(r,c)) {
    std::cout << r << ',' << c << '\n';
    present(r,c)=true;
  }  
  std::cout << '\n';
  spy(std::cout,present);
}

int main()
{
//   spin_system sys(3,"1H");
//   sys.isotope(1,"19F");
//   sys.isotope(2,"13C");

//   List<nuclei_spec> blockedby;
//   blockedby.push_back("1H");
//   blockedby.push_back("19F");
//   blockedby.push_back("13C");

//   const SpinOpGenerator opgen(sys,CrystalStructure(),blockedby);
  
//   dumppattern(opgen,"1H",'x');
//   dumppattern(opgen,"13C",'-');
//   dumppattern(opgen,"19F",'z');

  const blocking_information blockinfo(ScratchList<size_t>(2,3,2));
  std::cout << blockinfo;
  
  dumppattern(blockinfo,0,1);
  dumppattern(blockinfo,1,1);
  dumppattern(blockinfo,2,1);
  dumppattern(blockinfo,1,-2);
  dumppattern(blockinfo,0,0);
  dumppattern(blockinfo,1,0);
  dumppattern(blockinfo,2,0);
  
  return 0;
}
