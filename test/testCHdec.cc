
#include "MAS.h"
#include "geometry.h"
#include "ttyio.h"

using namespace libcmatrix;

int main(int argc,const char *argv[])
{
  int count=1;  
  //  const size_t nbeta=getint(argc,argv,count,"Number of beta steps ",30);
  const double DCH=1e3*getfloat(argc,argv,count,"D (kHz)? ",23);
  const double nurf=1e3*getfloat(argc,argv,count,"vRF (kHz)? ",102.06);
  const double charf=DCH*DCH*DCH/(nurf*nurf);
  std::cout << "D^3/vrf^2: " << charf << " kHz\n";
  const double nueff=nurf*sqrt(1.5);
  const double theta=MAGIC_ANGLE;
  //std::cout << "Theta: " << (theta*180.0/M_PI) << " degrees\n";
  const double C=6.0*sin(theta)*sin(theta)*cos(theta)/(nueff*nueff);
  //  std::cout << "nueff: " << nueff << " Hz";
  //std::cout << "C: " << C << '\n';

  const space_T D_PAS(spatial_tensor(DCH));
  std::cout << "Tensor in PAS: " << D_PAS << '\n';
  Euler PAS_to_RF(0.0,0.0,0.0);
  RotorInfo rinfo; //!< default is MAS
  const Tensor<double>& dvalues(rinfo.dvalues());

  double tot=0.0;
  double tweight=0.0;

  const double dwig1=dvalues(2,1);
  const double dwig2=dvalues(2,2);
  std::cout << "Wigner factors: " << dwig1 << "  " << dwig2 << '\n';

  //  for (size_t i=nbeta;i--;) {
  //const double beta=(M_PI*(i+0.5))/nbeta;    
  const double beta=MAGIC_ANGLE;
    PAS_to_RF.beta=beta;
    const space_T D_RF(rotate(D_PAS,PAS_to_RF));
    std::cout << D_RF << '\n';
//     std::cout << ((180.0/M_PI)*PAS_to_RF.beta) << ": ";
//     for (int j=-2;j<=2;j++) {
//       const complex dfac=D_RF(2,j)*dvalues(2,j);
//       std::cout << "  " << dfac;
//     }
//    std::cout << '\n';
    const double dfac2=real(D_RF(2,2))*dwig2;
    const double dfac1=real(D_RF(2,1))*dwig1;
    const double weight=sin(beta);
    const double contrib=C*dfac2*dfac1*dfac1;
    const double betadeg=(180.0/M_PI)*beta;
    std::cout << betadeg << ": " << dfac1 << ' ' << dfac2 << ' ' << contrib <<  '\n';
    tot+=contrib*weight;
    tweight+=weight;
    //  }
    //  std::cout << tweight << '\n';
  const double final=(tot/tweight);
  std::cout << "Constant: " << (final/charf) << '\n';
  //  std::cout << "Estimated mean LB: " << final << " Hz  Constant: " << (final*normfac) << '\n';
  return 0;
}
