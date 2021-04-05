// test meta blocking

#include "MetaPropagation.h"

using namespace libcmatrix;
using namespace std;

#define OPGEN SpinOpGenerator

int main()
{
  enum { XONE=1, XTHREEHALF=2, XTWOTHREEHALF=3, YNONE=4, YONE=8, YTWOHALF=12 };
  enum { XBLOCKED=1, YBLOCKED=2 };
  const char* ops="+-xz";

  size_t ndone=0;
  size_t nfailed=0;
  BlockedOperator sigma,sigmanormal;

  //  for (int system=0;system<16;system++) {
    int system=7; {
    const char* Xisotope="1H";
    size_t nspinsX=1;
    switch (system & 3) {
    case XONE:
      Xisotope="2H";
      break;
    case XTHREEHALF:
      Xisotope="23Na";
      break;
      //    case XTWOHALF:
    case XTWOTHREEHALF:
      Xisotope="11B";
      nspinsX=2;
      break;
    }
    const char* Yisotope="13C";
    size_t nspinsY=1;
    switch (system & 12) {
    case YNONE:
      nspinsY=0;
      break;
    case YONE:
      Yisotope="14N";
      break;
    case YTWOHALF:
      nspinsY=2;
      break;
    }

    spin_system sys(nspinsX+nspinsY,Xisotope);
    for (size_t i=nspinsY;i--;)
      sys(i+nspinsX).isotope(Yisotope);

    std::cout << "System: " << system << ": " << sys << '\n';

    const OPGEN opgen_notblocked(sys);

    for (size_t blocktype=0;blocktype<4;blocktype++) {
    //  size_t blocktype=1; {
      const bool yblocked=(blocktype & YBLOCKED);
      if (yblocked && (nspinsY==0))
	continue;
      std::cout << "Blocking type: " << blocktype << "   ";      

      const bool xblocked=(blocktype & XBLOCKED);      
      cout << "Blocked on X: " << (xblocked ? "Yes  " : "No  ");
      cout << "Blocked on Y: " << (yblocked ? "Yes\n" : "No\n");

      List<nuclei_spec> blocking;
      if (xblocked)
	blocking.push_back(nuclei_spec(Xisotope));
      if (yblocked)
	blocking.push_back(nuclei_spec(Yisotope));
    
      const OPGEN opgen_blocked(sys,blocking,0,2);

      cmatrix matrix_notblocked,matrix_blocked,difference;

      const int isoloop=nspinsY ? 2 : 1;
      for (int isotype=0;isotype<isoloop;isotype++) {
	for (int optype=0;optype<strlen("ops");optype++) {
	  //	  const operator_spec sigma0spec((isotype==0) ? Xisotope : Yisotope,ops[optype]);
	  const operator_spec sigma0spec(size_t(0),ops[optype]);
	  std::cout << "Operator: " << sigma0spec << '\n';

	  const BlockedOperator sigma0notblocked(opgen_notblocked,sigma0spec);
	  const BlockedOperator sigma0blocked(opgen_blocked,sigma0spec);

	  sigma0notblocked.full(matrix_notblocked,opgen_notblocked);
	  sigma0blocked.full(matrix_blocked,opgen_blocked);
	  
	  cmatrix difference(matrix_notblocked);
	  difference-=matrix_blocked;
	  const bool failed= (norm(difference)>1e-5);
	  std::cout << "Passed: " << (failed ? "No\n" : "Yes\n");
	  if (failed) {
	    std::cout << "Without blocking:\n" << matrix_notblocked;
	    std::cout << "Blocked:\n" << matrix_blocked << '\n';
	    nfailed++;
	  }
	  ndone++;
	}
      }
    }
  }
  std::cout << "\nFailures: " << nfailed << '/' << ndone << '\n';

  return nfailed;
}
