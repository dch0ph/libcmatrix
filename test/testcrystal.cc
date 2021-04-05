#include "CrystalSystem.h"
#include "Propagation.h"
#include "Histogram.h"
#include "MetaPropagation.h"
#include "Sequence.h"
#include "ttyio.h"

using namespace libcmatrix;
using namespace std;

FILE* fp=NULL;
const char nuctype[]="1H";
int verbose=0;
size_t gammasteps=0;
double rotor_speed=0.0;
size_t nobs=0;
double intdt=0.0;
double rf=0.0;
double sw=0.0;
bool isherm=false;

void readout(BaseHistogram<double>& hist, SpectrumIterator<double>& obj)
{
  double amp,freq;
  while (obj(amp,freq)) {
    if (fabs(amp)>1e-8) {
      if (verbose) {
	cout << amp << " of " << freq << '\n';
	if (fp)
	  fprintf(fp,"%g %g\n",amp,freq);
      }
      hist.add(amp,freq);
    }
  }
}

template<class OpGen> void docalc(BaseHistogram<double>& hist,const OpGen& opgen, const HamiltonianStore<space_T>& interactions)
{
  if (verbose)
    cout << "OpGen\n" << opgen << '\n';

  const operator_spec ospec("1H",isherm ? 'x' : '+');
  BlockedOperator sigma0det(opgen,ospec);
  if (verbose)
    cout << "sigma0/detect\n" << sigma0det << '\n';

  Sequence rfseq;
  static PulseGenerator* pgenp=NULL;
  if (rf) {
    if (!pgenp)
      pgenp=new PulseGenerator(opgen,"1H",verbose);
    rfseq.push_back(new CWPulse(*pgenp,1e6,rf,0.0),0.0);
    if (verbose)
      cout << "Sequence: " << rfseq << '\n';
  }
  
  if (rotor_speed) {
    BlockedSpinningHamiltonian<OpGen> H(opgen,interactions,rotor_speed,0.0);
    if (verbose)
      cout << "Hamiltonian\n" << H << '\n';
    const double rotor_period=1.0/fabs(rotor_speed);
    if (rf) {
      SequencePropagator propgen(H,intdt,rfseq,0.0,1e-9,verbose); //GammaPeriodic OK as rfseq is time independent
      MetaSpectrumED<GammaPeriodicSpectrumED> obj(propgen,gammasteps,0.0,rotor_period,sigma0det,GammaPeriodicSpectrumED(nobs,verbose),verbose);
      readout(hist,obj);
    }
    else {
      MetaPropagator propgen(H,intdt);
      MetaSpectrumED<GammaPeriodicSpectrumED> obj(propgen,gammasteps,0.0,rotor_period,sigma0det,GammaPeriodicSpectrumED(nobs,verbose),verbose);
      readout(hist,obj);
    }
  }
  else {
    const HamiltonianStore<double> interactions_LF(interactions,Euler(0,0,0));
    BlockedHamiltonian<OpGen> H(opgen,interactions_LF);
    if (verbose)
      cout << "Hamiltonian\n" << H << '\n';

    if (rf) {
      SequencePropagator Ugen(H,rfseq,0.0,1e-9,verbose);
      MetaSpectrumED<StaticSpectrumED> obj(Ugen,1,0.0,1.0/sw,sigma0det,verbose);
      readout(hist,obj);
    }
    else {
      MetaSpectrumED<StaticSpectrumED> obj(H,sigma0det,verbose);
      readout(hist,obj);
    }
  }
}

int main(int argc,const char *argv[])
{
  int count=1;

  try {
    const int M=getint(argc,argv,count,"Number of spins per cell? ",1);
  const int N=getint(argc,argv,count,"Number of cells? ",2);
  if (M*N<2)
    throw Failed("Must have at least 2 spins!");

  const double theta=(M>1) ? getfloat(argc,argv,count,"Angle between spin pair and translation axis? ",0.0)*M_PI/180.0 : 0.0;
  const double ratio=(M>1) ? getfloat(argc,argv,count,"Ratio of inter/intra cell distances? ",0.3) : 1.0;
  rotor_speed=1e3*getfloat(argc,argv,count,"Rotor speed (kHz)? ",0.0);
  if (rotor_speed) {
    nobs=getint(argc,argv,count,"Observations per rotor cycle? ",1);
    sw=nobs*rotor_speed;
    cout << "Spectral width: " << (sw/1e3) << " kHz\n";
    gammasteps=getint(argc,argv,count,"Gamma integration steps? ",4*nobs);
    intdt=1e-6*getfloat(argc,argv,count,"Integration time step (us)? ",2);
  }
  else 
    sw=1e3*getfloat(argc,argv,count,"Spectral width (kHz)? ");

  const int npts=getint(argc,argv,count,"Histogram points? ",1000);
  const bool testmode=getlogical(argc,argv,count,"Test? ",false);
  bool useksym=true;
  rf=1e3*getfloat(argc,argv,count,"RF nutation rate (kHz)? ",0.0);
  const bool isblocked=rf ? false : getlogical(argc,argv,count,"Use mz blocking? ",true);
  const bool usemzsym=false; //isblocked ? getlogical(argc,argv,count,"Exploit mz symmetry? ",true) : false;
  if (!testmode)
    useksym=(N>2) ? getlogical(argc,argv,count,"Exploit +/-k symmetry? ",true) : false;
  isherm=getlogical(argc,argv,count,"Hermitian sigma0/detect? ",false);

  char fname[128];
  getstring(argc,argv,count,"Output file? ",fname,128);
  verbose=getint(argc,argv,count,"Verbosity? ",1);

  const double r=1.5;
  const double R=r/ratio;

  List<vector3> coords(M);
  switch (M) {
  case 2: {
    const double x=(r/2)*cos(theta);
    const double y=(r/2)*sin(theta);
    coords(0)=vector3(x,y,0);
    coords(1)=vector3(-x,-y,0);
  }
    break;
  case 1:
    coords(0)=vector3(0,0,0);
    break;
  default:
    throw Failed("Not supported");
  }

  spinhalf_system lsys(M,nuctype);
  CrystalGeometry<1> geomspec(N,vector3(R,0,0),coords);
  if (verbose)
    cout << geomspec;

  HamiltonianStore<space_T> interactions_MF(M,N);
  const size_t total=M*N;
  //build coupling network
  for (int j=0;j<M;j++)
    for (size_t sk=j+1;sk<total;sk++) 
      interactions_MF.set_dipole(j,sk,geomspec.dipolar_tensor(lsys,j,sk));

  if (verbose)
    cout << "Interactions\n" << interactions_MF << '\n';

  List<double> spec(npts,0.0);
  Histogram<double> spechist(spec,sw);

  int flags= (useksym ? MetaFlags::UseEigSymmetry : 0) | (usemzsym ? MetaFlags::UseMzSymmetry : 0);

  fp=fopen(fname,"w");

  if (testmode) {
    SpinOpGenerator* opgenp=isblocked ? new SpinOpGenerator(lsys,N,"1H",flags) : new SpinOpGenerator(lsys,N);
    docalc(spechist,*opgenp,interactions_MF);
  }
  else {
    if (useksym)
      flags|=MetaFlags::UseEigSymmetry;
    CrystalOpGenerator* opgenp= isblocked ? new CrystalOpGenerator(lsys,N,"1H",flags,verbose) : new CrystalOpGenerator(lsys,N,flags,verbose);
    docalc(spechist,*opgenp,interactions_MF);
  }
  if (fp)
    fclose(fp);

  } catch (MatrixException& exc) {
    cerr << exc;
    return 1;
  }
  return 0;
}
