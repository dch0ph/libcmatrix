//#include "CrystalSystem.h"
#include "spinhalf_system.h"
#include "MetaPropagation.h"
#include "ttyio.h"

using namespace libcmatrix;
using namespace std;

#define OPGEN SpinOpGenerator

const size_t npts=8;
const double dt=1.0/300.0;

//functions are dubious because they ignore "middle block"

//  double hermitian_trace_multiply(const BlockedMatrix<complex>& a, const BlockedMatrix<complex>& b)
//  {
//    double sum=0.0;
//    if (!arematching(a,b))
//      throw Mismatch("hermitian_trace_multiply");

//    for (size_t i=a.size();i--;) {
//      if (!!b(i))
//        sum+=hermitian_trace_multiply(a(i),b(i));
//    }
//    return 2.0*sum;
//  }

//  complex trace_multiply(const BlockedMatrix<complex>& a, const BlockedMatrix<complex>& b)
//  {
//    complex sum(0.0);
//    if (!arematching(a,b))
//      throw Mismatch("trace_multiply");
//    for (size_t i=a.size();i--;) {
//      if (!!b(i))
//        sum+=trace_multiply(a(i),b(i));
//    }
//    return sum;
//  }

template<class OpGen> void docalc(const OpGen& opgen, const HamiltonianStore<double>& Hstore, char sigma0_op, char detect_op)
{
  const operator_spec sigma0spec(0U,sigma0_op);
  const operator_spec detectspec(1U,detect_op);  

  try {
    const BlockedOperator sigma0(opgen,sigma0spec);
    cout << "sigma0\n" << sigma0;
    const BlockedOperator detect(opgen,detectspec);
    cout << "detect\n" << detect;
    
    BlockedStaticHamiltonian<typename OpGen::base_type> Ham(opgen,Hstore);
    cout << "Hamiltonian\n" << Ham;
    
    BlockedMatrix<complex> U;
    propagator(U,Ham,dt);
    cout << "Propagator\n" << U << '\n';
    
    BlockedOperator sigma(sigma0);
    List<complex> FID(npts,complex(0.0));
    for (size_t n=0;n<npts;n++) {
      FID(n)=trace_multiply(sigma,detect);
      sigma.unitary_simtrans(U);
    }
    cout << "FID: " << FID << '\n';
  }
  catch (MatrixException& exc) {
    cerr << exc << '\n';
  }
}

int main(int argc, const char* argv[])
{
  HamiltonianStore<double> Hstore(2);
  Hstore.set_coupling(I_DIPOLE,0,1,0,1000);
  
  const spinhalf_system sys(2,"1H");

  int count=1;
  const bool usemz=getlogical(argc,argv,count,"Use mz symmetry? ",true);
  const int flags=usemz ? MetaFlags::UseMzSymmetry : 0;
  //OPGEN opgen_unblocked(sys,flags,1);
  const List<nuclei_spec> blocking(1,nuclei_spec("1H"));
  OPGEN opgen_blocked(sys,blocking,flags,1);

  for (size_t which=0;which<=3;which++) {
    const char sigma0_op= (which & 1) ? 'x' : '+';
    const char detect_op= (which & 2) ? 'x' : '-';
    cout << "sigma0: " << sigma0_op << '\n';
    cout << "detect: " << detect_op << '\n';

    //cout << "Unblocked\n";
    //docalc(opgen_unblocked,Hstore,sigma0_op,detect_op);
    cout << "Blocked\n";
    docalc(opgen_blocked,Hstore,sigma0_op,detect_op);
  }

  return 0;
}
