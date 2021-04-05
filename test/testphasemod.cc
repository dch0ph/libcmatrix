#include "PhaseModulatedPropagator.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

int main(int argc, const char* argv[])
{
  int count=1;
  spin_system sys(2,"1H");
  SpinOpGenerator opgen(sys);
  PulseGenerator pgen(opgen,"1H",1);
  
  HamiltonianStore<space_T> Hstore(2);
  Hstore.set_coupling(I_DIPOLE,0,1,rotate(spatial_tensor(10000.0),Euler(0,M_PI/3,0)));

  cout.precision(5);
   const double vrf=150e3;
   const double rotor_speed=33.333333e3;
   const double maxdt=1e-6;
   const double dt=1.0/rotor_speed;
   const double period=0.8*dt;
   const size_t nperiods=3;
   const double tp=0.5*period;
   Sequence XiX;
   XiX.push_back(new CWPulse(pgen,tp,vrf,0.0),0.0);
   XiX.push_back(new CWPulse(pgen,tp,vrf,M_PI),tp);

   cout << XiX << '\n';

   PhaseModulation pmod(XiX,period,1e-9,2);
   pmod.build(2);
   cout << pmod << endl;

   BlockedSpinningHamiltonian<double> Ham(opgen,Hstore,rotor_speed,0.0);
   
   const int explicitn=getint(argc,argv,count,"Rotor cycle steps (0 for normal PM)? ",20);
   if (explicitn)
     std::cout << "Cache step: " << (1e6*dt/explicitn) << " us\n";

   PhaseModulatedPropagator propgen(Ham,NULL,maxdt,0.0,pmod,0.0,NULL,1e-9,2,0,explicitn,maxdt);
   BlockedMatrix<complex> U;
   for (size_t n=0;n<nperiods;n++) {
     const double t=n*dt;
     propgen(U,t,t+dt);
     std::cout << "Propagator for rotor period " << n << '\n' << U << '\n';
   }  
   return 0;
}
