/* test effect of phase transients */

#include "Sequence.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

const double rad_to_deg=180.0/M_PI;

struct trans_state : public complex {

  trans_state() {}
  trans_state(double inv,double outv) : complex(inv,outv) {}
  complex Brf(const complex& preBrf, const complex& postBrf) const;

  void Hrf(cmatrix& Hrf, const cmatrix& preHrf_in, const cmatrix& preHrf_out, const cmatrix& postHrf_in, const cmatrix& postHrf_out) const;
};

void trans_state::Hrf(cmatrix& Hrf, const cmatrix& preHrf_in, const cmatrix& preHrf_out, const cmatrix& postHrf_in, const cmatrix& postHrf_out) const 
{
  multiply(Hrf,re,preHrf_in);
  mla(Hrf,im,preHrf_out);
  mla(Hrf,1.0-re,postHrf_in);
  mla(Hrf,-im,postHrf_out);
}

complex trans_state::Brf(const complex& preBrf, const complex& postBrf) const
{
  const complex deltaBrf(preBrf-postBrf);
  return postBrf+re*deltaBrf+im*complex(-deltaBrf.im,deltaBrf.re);
}

inline std::ostream& operator<< (std::ostream& ostr, const trans_state& a)
{
  return ostr << "In: " << a.re << "  Out: " << a.im;
}

class transient_model {
public:
  transient_model(double trans_amp =0.0, double tauQ =0.0, double resv =0.0);

  void reset() {
    facQ_=0.0;
    facT_= lambdaT_ ? 1.0 : 0.0;
  }

  void initialise(double preampv, double prephasev, double postampv, double postphasev, double resv =0.0);
  trans_state operator()();
  trans_state operator()(double t) const;
private:
  double lambdaT_,lambdaQ_;
  const double res_;
  double facT_step_,facQ_step_;
  double facT_,facQ_;
};

transient_model::transient_model(double trans_amp, double tauQ, double resv)
  : lambdaT_(tauQ ? 1.0/tauQ : 0.0), res_(resv)
{
  lambdaQ_=lambdaT_*lambdaT_*(trans_amp/360e3);
  std::cout << "lambdaT: " << lambdaT_ << " Hz\nlambdaQ: " << lambdaQ_ << " Hz\n";
  if (res_) {
    facT_step_=exp(-lambdaT_*res_);
    facQ_step_=lambdaQ_*res_;
  }
  reset();
}

trans_state transient_model::operator()(double t) const
{
  if (t<0.0)
    return trans_state(1.0,0.0);
  if (lambdaT_) {
    const double in=exp(-lambdaT_*t);
    return trans_state(in,in*lambdaQ_*t);
  }
  return trans_state(0.0,0.0);
}

trans_state transient_model::operator()()
{
  if (res_==0.0)
    throw Failed("transient_model: default time resolution unset");
  const trans_state tmp(trans_state(facT_,facT_*facQ_));
  facT_*=facT_step_;
  facQ_+=facQ_step_;
  return tmp;
}

int main(int argc,const char* argv[])
{
  try {

    cout.precision(4);
  const spin_system sys(1,"1H");
  const List<double> sigma0(diag_Fz(sys));
  const rmatrix detect(real(F(sys,'+')));
  const double dt=50e-9; //time resolution after pulse
  const size_t N=25; //points to sample after pulse

  const SpinOpGenerator opgen(sys);
  const BlockedStaticHamiltonian<double> H(opgen);//empty Hamiltonian

  const double vRF=50e3;
  const double tip=M_PI/2;
  const double tau=(tip/(2*M_PI))*(1/vRF);
  cout << "Tip duration: " << (tau*1e6) << " us\n";

  enum { TRANS_NONE=0, TRANS_SIMPLE, TRANS_FULL };
  int count=1;
  cout << "N - None\nS - Simple\nF - full\n";
  int transtype=getoption(argc,argv,count,"Transient model? ","NSF");

  PulseGeneratorBase* pgenp=NULL;
  const double trans_amp=(transtype!=TRANS_NONE) ? getfloat(argc,argv,count,"Transient amplitude (degrees / kHz)? ",0.01) : 0.0;

  double tauT=0.0;

  switch (transtype) {
  case TRANS_NONE:
    pgenp=new PulseGenerator(opgen,NULL_NUCLEUS,2);
    break;
  case TRANS_SIMPLE: {
    PulseGenerator_SimpleTransients* transgenp=new PulseGenerator_SimpleTransients(opgen,NULL_NUCLEUS,2);
    transgenp->transient_amplitude(trans_amp);
    pgenp=transgenp;
  }
    break;
  case TRANS_FULL:
    pgenp=new PulseGenerator(opgen,NULL_NUCLEUS,2);
    tauT=getfloat(argc,argv,count,"Probe (Q) time constant (us)? ",0.3)*1e-6;
    break;
  default:
    throw Failed("Unsupported model");
  }

  if (transtype==TRANS_FULL) {
    complex Boff(0.0,0.0);
    complex Bon(vRF,0.0);

    transient_model fulltrans(trans_amp,tauT,dt);
    cout << "Pulse on\n";
    double inttrans=0.0;
    for (size_t n=0;n<N;n++) {
      const trans_state state=fulltrans();
      const complex Brf(state.Brf(Boff,Bon));
      const double t=n*dt;
      cout << "At t=" << (1e6*t) << " us: \t" << state << " \tBrf: " << Brf << " \tAmp: " << (abs(Brf)*1e-3) << " kHz  Phase: " << (arg(Brf)*rad_to_deg) << '\n';
      inttrans+=imag(Brf)*dt;
    }
    cout << "Integrated out-of-phase tip: " << (inttrans*360.0) << " degrees  (Ideal: " << (trans_amp*vRF*1e-3) << ")\n";
    cout << "\nPulse off\n";
    fulltrans.reset();
    for (size_t n=0;n<N;n++) {
      const trans_state state=fulltrans();
      const complex Brf(state.Brf(Bon,Boff));
      const double t=n*dt;
      cout << "At t=" << (1e6*t) << " us: \t" << state << " \tBrf: " << Brf << " \tAmp: " << (abs(Brf)*1e-3) << " kHz  Phase: " << (arg(Brf)*rad_to_deg) << '\n';
    }    
    return 0;
  }

  Sequence seq;  
  
  seq.push_back(new CWPulse(*pgenp,tau,vRF,'y'),0.0);
  cout << "Sequence:\n" << seq;

  SequencePropagator Ugen(H,seq,0.0,0.0,1e-9,2);
  BlockedMatrix<complex> U;
  Ugen(U,0.0,tau);
  cout << "Pulse propagator\n" << U.front();

  cmatrix sigma;
  unitary_simtrans(sigma,sigma0,U.front());
  cout << "Density matrix after pulse\n" << sigma;
  for (size_t n=0;n<N;n++) {
    const double t=tau+n*dt;
    cout << "xy magnetisation at " << (t*1e6) << " us: " << trace_multiply(detect,sigma) << '\n';
    Ugen(U,t,t+dt);
    const cmatrix& U0(U.front());
    if (!!U0)
      sigma.unitary_simtrans(U0);
  }
  } catch (MatrixException& exc) {
    cerr << exc << '\n';
    return 1;
  }
  return 0;
}
