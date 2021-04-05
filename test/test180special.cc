// test objects for creating ideal 180 inversion pulses

#include "lcm_CommonSequence.hpp"
#include "MetaPropagation.h"
#include "InversionGenerator.h"

using namespace libcmatrix;
using namespace std;


bool xblocked=true;
bool randomise=false;//true;

#define OPGEN SpinOpGenerator

void cleanup(BaseList<complex> a, double tol =1e-10)
{
  for (size_t i=a.size();i--;) {
    complex& z(a(i));
    if (fabs(z.real())<tol)
      z.real(0.0);
    if (fabs(z.imag())<tol)
      z.imag(0.0);
  }
}

int main(int argc, const char* argv[])
{
  enum { XONE=2, XTHREEHALF=4, XTWOTHREEHALF=6, YNONE=8, YTHREEHALF=16, YTWOHALF=24, YBLOCKED=1 };
  enum { STARTX=1, APPLYy=2 };

  if (argc>1) {
    if (strcmp(argv[1],"-noxblocking")==0)
      xblocked=false;
    else {
      std::cerr << "Unrecognised flag: " << argv[1] << '\n';
      return 1;
    }
  }

  size_t ndone=0;
  size_t nfailed=0;
  BlockedOperator sigma,sigmanormal;
  BlockedMatrix<complex> Upulse;
  cmatrix sigmaexpand,sigmanormalexpand;

  for (int system=0;system<32;system++) {
    //   int system=14; {
    const char* Xisotope="1H";
    size_t nspinsX=1;
    switch (system & 6) {
    case XONE:
      Xisotope="2H";
      break;
    case XTHREEHALF:
      Xisotope="23Na";
      break;
      //case XTWOHALF:
    case XTWOTHREEHALF:
      Xisotope="11B";
      nspinsX=2;
      break;
    }
    const char* Yisotope="13C";
    size_t nspinsY=1;
    switch (system & 24) {
    case YNONE:
      nspinsY=0;
      break;
//     case YONE:
//       Yisotope="14N";
//       break;
    case YTHREEHALF:
      Yisotope="35Cl";
      break;
    case YTWOHALF:
      nspinsY=2;
      break;
    }

    spin_system sys(nspinsX+nspinsY,Xisotope);
    for (size_t i=nspinsY;i--;)
      sys(i+nspinsX).isotope(Yisotope);

    std::cout << "System: " << system << ": " << sys << '\n';
    const bool yblocked=(system & YBLOCKED);
    cout << "Blocked on spin Y: " << (yblocked ? "Yes\n" : "No\n");

    const nuclei_spec Xspec(Xisotope);
    List<nuclei_spec> blocking;
    if (xblocked)
      blocking.push_back(Xspec);
    List<nuclei_spec> RFblocking;
    if (yblocked) {
      blocking.push_back(nuclei_spec(Yisotope));
      RFblocking.push_back(nuclei_spec(Yisotope));
    }
    
    const OPGEN opgen(sys,blocking,0,2);
    const OPGEN RFopgen(sys,RFblocking,0,2);

    const InversionGenerator invgen(opgen,blocking,Xspec);
    cout << invgen << '\n';
     const PulseGenerator pgen(RFopgen,Xspec);

     //     for (int pulsetype=0;pulsetype<4;pulsetype++) {
       int pulsetype=1; {
       std::cout << "pulsetype: " << pulsetype << '\n';
       //       const operator_spec sigma0spec((pulsetype & STARTX) ? Xisotope : Yisotope,'x');
       const operator_spec sigma0spec(size_t(0),'y');
       const double phase=(pulsetype & APPLYy) ? M_PI/2 : 0;
       std::cout << "Starting operator: " << sigma0spec << '\n';
       std::cout << "Phase: " << ((pulsetype & APPLYy) ? "y\n" : "x\n");
      
       const BlockedOperator sigma0(opgen,sigma0spec);
       BaseList<complex> rawop(sigma0.row().row());
       if (randomise) {
	 for (size_t i=rawop.size();i--;)
	   rawop(i)=complex(i,i+1);
       }
       cout << "sigma0\n" << sigma0;

       invgen(sigma,sigma0,phase);
       cout << "sigma (post inversion)\n" << sigma;

       BlockedOperator sigma0normal(RFopgen,sigma0spec);
       if (randomise) {
	 sigma0.full(sigmaexpand,opgen);
	 const BaseList<complex> source(sigmaexpand.row());
	 BaseList<complex> dest(sigma0normal.row().row());
	 dest=source;
       }

       pgen(Upulse,M_PI,phase);
       unitary_simtrans(sigmanormal,sigma0normal,Upulse);
       cleanup(Upulse.row());
       cout << "Upulse\n" << Upulse << '\n';
       cleanup(sigmanormal.row().row());
       cout << "sigma (normal)\n" << sigmanormal;

       sigmanormal.full(sigmanormalexpand,RFopgen);
       sigma.full(sigmaexpand,opgen);
       cmatrix difference(sigmaexpand);
       difference-=sigmanormalexpand;
       const bool failed= (norm(difference)>1e-5);
       std::cout << "Passed: " << (failed ? "No\n" : "Yes\n");
       if (failed) {
	 std::cout << "After normal inversion:\n" << sigmanormalexpand;
	 std::cout << "After special inversion:\n" << sigmaexpand << '\n';
// 	 BlockedMatrix<complex> tmp,tmp2;
// 	 multiply(tmp,Upulse,sigma0normal);
// 	 std::cout << "Upulse*sigma:\n" << tmp << '\n';
// 	 multiply_conj_transpose(tmp2,tmp,Upulse);
// 	 std::cout << "Upulse*sigma*Upulse':\n" << tmp2 << '\n';
	 nfailed++;
       }
       ndone++;
     }
    }
    std::cout << "\nFailures: " << nfailed << '/' << ndone << '\n';

    return nfailed;
}
