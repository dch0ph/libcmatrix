/* basic simulation of RFDR */

#include "Sequence.h"
#include "MAS.h"
#include "Propagation.h"
#include "MetaPropagation.h"
#include "ttyio.h"

using namespace libcmatrix;
using namespace std;

int main(int argc,const char* argv[])
{
  try {

    //BlockedSpinningHamiltonian<SpinOpGenerator> tmp;
  int count=1;

  const double rotor_phase=0.0;
  const bool twospins=getlogical(argc,argv,count,"Two spins? ");
  const double dHH=twospins? 1e3*getfloat(argc,argv,count,"HH coupling (kHz)? ",10.0) : 0.0;
  const double rotor_speed=1e3*getfloat(argc,argv,count,"Rotor speed (kHz)? ",20.0);
  if (rotor_speed==0.0) throw InvalidParameter("rotor_speed");
  //  const int maxn=getint(argc,argv,count,"Maximum rotor cycles? ");
  //  const double dur=getfloat(argc,argv,count,"Nominal pulse duration (us)? ",5.0)*1e-6;
  const double vRF=1e3*getfloat(argc,argv,count,"RF field strength (kHz, 0 for ideal pulses)? ",0.0);
  const int int_steps=getint(argc,argv,count,"Integration steps per rotor cycle? ",200);

  spin_system sys(twospins ? 2 : 1,"1H");

  const cmatrix Ix(I(sys,0,'x'));
  const cmatrix Iy(I(sys,0,'y'));
  cout << "Iy\n" << Iy;

  SpinOpGenerator opgen(sys);

  //Test PulseGenerator
  PulseGenerator Ugen(opgen);
  cout << "90x\n" << Ugen(M_PI/2,'x');
  cout << "90y\n" << Ugen(M_PI/2,'y');
  cout << "180\n" << Ugen(M_PI,'x');

  const cmatrix simple90x(hermitian_exp(Ix,complex(0,-M_PI/2)));
  cout << "exp(-i*(pi/2)Ix)\n" << simple90x;

  const cmatrix simple90y(hermitian_exp(Iy,complex(0,-M_PI/2)));
  cout << "exp(-i*(pi/2)Iy)\n" << simple90y;
  
  RFEvent* eventp;
  if (vRF)
    eventp=new SoftPulse(Ugen,0.5/vRF,vRF,'x');
  else
    eventp=new HardPulse(Ugen,10e-6,50e3,'x'); //use 50 kHz (purely nominal);

  //Output pulse used
  cout << "Pulse object: " << (*eventp) << '\n';

  const double rotor_period=1.0/fabs(rotor_speed);

  //Create sequence
  Sequence RFDR(rotor_period);
  RFDR.push_back(eventp,0.5);
  
  //  cout << "Sequence\n" << RFDR << endl;
  cout << "Sequence\n" << RFDR << endl;

  //Create Hamiltonian
  Euler PAS_to_RF(0.0,M_PI/3,0.0);
  const space_T A_RF=rotate(spatial_tensor(dHH),PAS_to_RF);
  cout << "Dipolar tensor\n" << A_RF << endl;

  HamiltonianStore<space_T> interactions(sys.nspins());
  if (twospins)
    interactions.set_coupling(I_DIPOLE,0,1,0.0,A_RF);

  BlockedSpinningHamiltonian<double> Ham(opgen,interactions,rotor_speed,rotor_phase);
  //const cmatrix Hdip=spin_dipolar(sys,0,1);
  //Ham.add(A_RF,Hdip);
  cout << "Spinning Hamiltonian:\n" << Ham << endl;

  BlockedMatrix<double> H0;
  BlockedMatrix<complex> Usoft;
  Ham(H0,0.5*rotor_period);
  if (vRF) {
    Ugen(Usoft,H0,0.5/vRF,0.0,M_PI,'x');
    cout << "Usoft:\n" << Usoft << '\n';
  }

  //Create propagator over rotor cycle
  SequencePropagator Ugen_rf(Ham,rotor_period/int_steps,RFDR,0.0,rotor_period,1e-9,2);
  BlockedMatrix<complex> Ucycle;
  Ugen_rf(Ucycle,0.0,rotor_period);
  cout << "Cycle propagator\n" << Ucycle << endl;

  //SequencePropagator Ugen_norf(Ham,rotor_speed,rotor_phase,rotor_period/int_steps);
  //cout << "Cycle propagator without RF\n" << Ugen_norf(0.0,rotor_period) << endl;

  //Calculate FID
  //StaticFID_U obj;
  //obj.set_U(Ucycle);
  //const cmatrix sigma0=I(sys,0,'z');
  //const cmatrix detect=I(sys,1,'z');
  //cout << "'FID': " << obj.FID_hermitian(maxn,sigma0,detect);

  //FSLG
  if (vRF) {
    const HamiltonianStore<double> sinteractions(interactions,Euler(0,0,0));
    BlockedStaticHamiltonian<complex> sHam(opgen,sinteractions);
    cout << "Static Hamiltonian:\n" << sHam << endl;
    
    PulseGenerator Ugen_cache(opgen,sHam());
    cout << "90x\n" << Ugen_cache(M_PI/2,'x');
    cout << "90y\n" << Ugen_cache(M_PI/2,'y');
 
    const double offset=vRF/sqrt(2.0);
    const double fslg_period=2.0/(sqrt(3.0/2)*vRF);

    CWPulse LGplus(Ugen_cache,0.5,vRF,0.0,offset);
    CWPulse LGminus(Ugen_cache,0.5,vRF,180.0,-offset);

    Sequence FSLG(fslg_period);
    FSLG.push_back(&LGplus,0.0);
    FSLG.push_back(&LGminus,0.5);
    cout << "FSLG Sequence\n" << FSLG << endl;

    SequencePropagator Ugen_FSLG(sHam,FSLG,fslg_period,1e-9,2);
    Ugen_FSLG(Ucycle,0.0,fslg_period);
    cout << "FSLG propagator\n" << Ucycle << endl;
  }
  
  } catch (MatrixException& exc) {
    cerr << exc;
    return 1;
  }
  return 0;
}

