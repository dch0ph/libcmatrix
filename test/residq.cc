/* Calculate spectrum due to relayed quadrupolar coupling under MAS */

#include "NMR.h"
#include "MAS.h"
#include "powder.h"
#include "Propagation.h"
#include "tensorop.h"
#include "ttyio.h"
#include <cstring>
#include "cmatrix_utils.h"
#include "optim.h"
#include "Histogram.h"
#include "simpsonio.h"

using namespace std;
using namespace libcmatrix;

double df=0.0;
double rphase=0.0;

List<cmatrix> NT2Q(5);
//List<cmatrix> MT2Q(5);
cmatrix MQ;
List<cmatrix> HdipM(3);
List<double> Hzeeman,Hcs;

cmatrix sigma0det;

int intsteps,nprops,nobs,zcw;

//bool fullH=false;
bool haveM,secularM;
bool manualgamma=false;

double larmorM;

double sw=0.0;
double specleft=0.0;
size_t npts=0;

bool tdomain=false;

const char Ntype[]="14N";
const char Mtype[]="1H";

Euler default_orient(0,M_PI/6,0);

enum { AMP, MISO, QPOLE, QETA, DIPOLE, DTOQBETA, DTOQGAMMA, JMX, MCSA, CSAETA, CSATOQALPHA, CSATOQBETA, CSATOQGAMMA, QPOLEM, QETAM, QMTOQALPHA, QMTOQBETA, QMTOQGAMMA, SPEED, LW };
const size_t NPARS=LW+1;

const char *par_names[NPARS]={"amplitude","isotropic",
  "Q constant","Q asymmetry","NM D coupling",
    "D->Q beta","D->Q gamma",
    "JMX",
    "CSA","CSA asymmetry",
    "CSA->Q alpha","CSA->Q beta","CSA->Q gamma",
    "M Q constant","M Q asymmetry",
    "QM->Q alpha","QM->Q beta","QM->Q gamma",
    "rotor speed", "linewidth"};

//errors are in *internal* units (not input)
double rerrs[NPARS]={1.0, 1.0,
			     1e5, 0.2, 100,
			     0.2,0.2, //angles in radians
			     5.0,//J's in Hz
			     100,0.1,
			     0.2,0.2,0.2,//angles in radians
			     1e5,0.1,
			     0.2,0.2,0.2,//angles in radians
			     10.0,2.0};

BaseList<double> errs(NPARS,rerrs);

inline double force_eta(double eta) {
  if (eta<-1.0) return -1.0;
  if (eta>1.0) return 1.0;
  return eta;
}
  
inline bool isfixed(istream& istr)
{
  return (istr.peek()=='F') ? (istr.get(),true) : false;
}

void dump_diag(List<double>& eigs,const cmatrix& Ucycle,double period)
{
  cmatrix V;

  cout << "Cycle propagator\n" << Ucycle << endl;
  spy(Ucycle);
  diag_propagator(V,eigs,Ucycle,period);
  cout << "Eigenvalues: " << eigs << endl;
}

void dump_diag(List<double>& eigs,const cmatrix& H)
{
  cmatrix V;

  hermitian_eigensystem(V,eigs,H);
  cout << "Eigenvalues: " << eigs << endl;
  //  cout << "Values:\n" << V << endl;
}

void calc_spec(BaseList<double> spec,const BaseList<double>& pars, void* =NULL)
{
  const size_t inpts=spec.length();

  if (!ispowerof2(inpts)) throw Failed("Data length must be power of 2");

  FoldingInterpHistogram<double> histo(spec,sw);
  spec=0.0;

  static space_T CSA_PAS, Q_PAS, D_PAS, CSA_QPAS, D_QPAS, DM_QPAS, QM_QPAS;

  const Euler D_to_Q(0,pars(int(DTOQBETA)),pars(int(DTOQGAMMA)));
  const Euler CSA_to_Q(pars(int(CSATOQALPHA)),pars(int(CSATOQBETA)),pars(int(CSATOQGAMMA)));
  const Euler QM_to_Q(pars(int(QMTOQALPHA)),pars(int(QMTOQBETA)),pars(int(QMTOQGAMMA)));

  if (zcw<3) {
    cout << "D_to_Q: " << D_to_Q << endl;
    cout << "CSA_to_Q: " << CSA_to_Q << endl;
    if (haveM) cout << "QM_to_Q: " << QM_to_Q << endl;
  }

  double xx,yy,zz;

  const size_t mloop=1;
  //  const int mloop=(haveM && secularM) ? HdipMsec.length() : 1;

  const double scale=pars(int(AMP))/(mloop*NT2Q(2).rows());

  const double CSA=larmorM*pars(int(MCSA))/1e6;
  const double iso=0.0;//larmorM*pars(int(MISO))/1e6-specleft;

  const bool haveCS=((CSA!=0.0) || (iso!=0.0));

  if (haveCS) {
    //Use spatial tensor here since Hamiltonian is high-field
    CSA_QPAS=rotate(spatial_tensor(iso,CSA,force_eta(pars(int(CSAETA)))),CSA_to_Q);
    if (zcw<3) cout << "CSA_QPAS: " << CSA_QPAS << endl;
  }

  anisotropy_to_cartesian(xx,yy,zz,pars(int(QPOLE)),force_eta(pars(int(QETA))));  
  Q_PAS=A2(xx,yy,zz,2);
  Q_PAS*=0.666666; //fudge factor

  anisotropy_to_cartesian(xx,yy,zz,0.5*pars(size_t(JMX)),pars(size_t(DIPOLE)),0.0);  
  D_QPAS=rotate(A2(xx,yy,zz),D_to_Q);

  if (zcw<3) {
    cout << "Q_PAS: " << Q_PAS;
    cout << "D_QPAS: " << D_QPAS;
  }

  //  const double QpoleM=pars(size_t(QPOLEM));
  const double QetaM=force_eta(pars(size_t(QETAM)));
  anisotropy_to_cartesian(xx,yy,zz,pars(size_t(QPOLEM)),QetaM);
  QM_QPAS=rotate(A2(xx,yy,zz,2),QM_to_Q);
  QM_QPAS*=0.666666; //fudge factor
  
  if (zcw<3)
    cout << "QM_QPAS: " << QM_QPAS;

  List<cmatrix> Ualphas(nprops);
  List<cmatrix> Ubetas(nprops);
  List<cmatrix> Us(nprops);
  
  List<double> eigs,eigsa,eigsb,eigsaH,eigsbH;

  const double rotor_speed=pars(size_t(SPEED));
  const double period=1.0/fabs(rotor_speed);
  
  SpinningHamiltonian Ham(rotor_speed);
  //  SpinningHamiltonian Hcouple(rotor_speed);
  //  SpinningHamiltonian Htmp;

  const double dt=1.0/sw;

  Euler powder=default_orient;
  cmatrix Dvals;
  double weight=1.0;

  List<double> damping(inpts);

  const double lw=pars(size_t(LW));
  if (lw) {
    const double dampf=exp(-dt*M_PI*lw);

    double dv=1.0/inpts;
    for (int i=0;i<inpts;i++) {
      damping(i)=dv;
      dv*=dampf;
    }
    damping(0)*=0.5;
    if (zcw<3) cout << "Damped to: " << dv << endl;
  }

  //cout << "Decay constants: " << dcons << endl;

  PlanarZCW powdmeth(zcw);

  const double intdt=period/intsteps;

  if (zcw<3)
    cout << "Integration step: " << (intdt*1e6) << " us\n";

  List<complex> FID(inpts,0.0);
  space_T CSA_RF;

  while (powdmeth.next(powder,weight)) {
    if (zcw<4) 
      cout << powder << " (" << weight << ")" << endl;
    else {
      cout << '.';
      cout.flush();
    }

    weight*=scale;
    
    Dvals=D(2,powder);

    const space_T Q_RF(rotate(Q_PAS,Dvals));
    const space_T QM_RF(rotate(QM_QPAS,Dvals));
    const space_T D_RF(rotate(D_QPAS,Dvals));
    if (haveCS)
      CSA_RF=rotate(CSA_QPAS,Dvals);

    //    for (int m=mloop;m--;) {
    //  if (zcw<3) cout << "M state: " << m << endl;

      Ham.reset();
      //Ham.add(Q_RF,NT2Q);
      Ham.add(Q_RF,NT2Q(2));
      if (!!MQ)
	Ham.add(QM_RF,MQ);
       Ham+=Hzeeman;
      
      //      if (fullH) {
      Ham.add(D_RF,HdipM);
	if (haveCS)
	  Ham.add(CSA_RF,Hcs);
		
	if (zcw<3) {
	  cout << "H\n" << Ham << endl;
	  cout << "H(0)\n" << Ham(0) << endl;
	  cout << "H(tau_r/5)\n" << Ham(period/5.0) << endl;
	  cout << "H(tau_r)\n" << Ham(period) << endl;
	}
	cmatrix Utmp;
	propagator(Utmp,Ham(0),period);
	cout << "Fake U\n" << Utmp;
	
	HomogeneousPropagator<cmatrix> Ugen(Ham,intdt);
	propagators(Us,Ugen,0,period);
	if (zcw<3) dump_diag(eigs,Us(nprops-1),period);
//       }
//       else {
// 	Hcouple.reset();
// 	Hcouple.add(D_RF,Hdip);
// 	if (haveCS)
// 	  Hcouple.add(CSA_RF,0.5);

// 	if (zcw<3) {
// 	  cout << "Hbase\n" << Ham << endl;
// 	  cout << "Hcouple\n" << Hcouple << endl;
// 	}
	
// 	Htmp=Ham;
// 	Htmp+=Hcouple;
	
// 	if (zcw<3) {
// 	  cout << "H+\n" << Htmp << endl;
// 	  dump_diag(eigsaH,Htmp.component(0));
	  
// 	  cout << "H(0)\n" << Htmp(0) << endl;
// 	  cout << "H(pi/3)\n" << Htmp(M_PI/3) << endl;
// 	  dump_diag(eigs,Htmp(M_PI/3));
// 	}
// 	StepMASPropagator<cmatrix> Ugena(Htmp,rotor_speed);
// 	Ugena.rotor_phase(rphase);
// 	propagators(Ualphas,Ugena,0,period,intdt);
// 	//if (zcw<3) cout << "alpha Propagator(0)\n" << Ualphas << endl;
// 	if (zcw<3) dump_diag(eigsa,Ualphas(nprops-1),period);
	
// 	Ham-=Hcouple;
	
// 	if (zcw<3) {
// 	  cout << "H-\n" << Ham << endl;
// 	  dump_diag(eigsbH,Ham.component(0));
// 	  cout << "H eig diff: " << (eigsaH-eigsbH) << endl;
// 	}
// 	StepMASPropagator<cmatrix> Ugenb(Ham,rotor_speed);
// 	Ugenb.rotor_phase(rphase);
// 	propagators(Ubetas,Ugenb,0,period,intdt);
// 	if (zcw<3) {
// 	  dump_diag(eigsb,Ubetas(nprops-1),period);
// 	  cout << "Eig diff: " << (eigsa-eigsb) << endl;
// 	}
//       }
      
     //SimpleSpectrumED obj;
      
      //	add_GammaPeriodicFID(FID,weight,nobs,Ualphas,Ubetas);	
      GammaPeriodicSpectrumED obj(nobs);
      
      //      if (fullH) {
	obj.set_Us(Us,period);
	obj.observe(sigma0det);
 //      }
//       else {
// 	obj.set_Us('R',Ualphas,period);
// 	obj.set_Us('C',Ubetas,period);
// 	obj.observe();
//       }
      if (zcw<3) cout << obj;
      
      double a,f;
      while (obj(a,f)) {
	if (fabs(a)>1e-10) {
	  a*=weight;
	  f-=sw/2.0;
	  if (f<0.0)
	    f+=sw;
	  if (zcw<3) cout << a << " of " << f << endl;
	  if (tdomain)
	    add_FID(FID,a,f/sw);
	  else
	    histo.add(a,f);
	}     
      }
    }
  //} //mloop

  if (zcw>=4) cout << '\n';

  if (lw) {
    if (!tdomain) FID=fft(spec,FT_BACKWARD);
    FID*=damping;
    fft_ip(FID,FT_FORWARD);
    real(spec,FID);
  }
}

double moment(const BaseList<double>& a)
{
  const double s=sum(a);
  double m=0.0;
  for (size_t i=a.length();i--;) m+=i*a(i);
  return m/s;
}

int main(int argc,const char *argv[])
{
  try {

    //cout << dipolar_coupling(gamma("14N"),gamma("13C"),1.081e-10) << endl;
    //cout << dipolar_coupling(gamma("65Cu"),gamma("13C"),1.876e-10) << endl;

  List<double> pars(NPARS,0.0);
  List<size_t> varywhich(NPARS);
  size_t point=0;

  args_iter getargs(1,argc,argv,cin,cout);

  double larmorH;
  getargs(larmorH,"1H Larmor frequency (MHz)? "); larmorH*=1e6;
  const double gamma1H=gamma("1H");

  const double larmorN=larmorH*gamma(Ntype)/gamma1H;
  
  cout << "Larmor frequency of 14N: " << larmorN/1e6 << " MHz\n";

  char fname[128];
  FILE *fpin;
  getargs(fpin,"Spectrum to fit <CR if none>? ","r",fname,128,mxflag::allownull);

  size_t inpts=0;
  List<double> spec;
  const bool havespec=(fpin!=NULL);

  double midp=0.0;

  rmatrix rspec;

  if (havespec) {
    read_matrix(rspec,fpin);
    spec=rspec(slice(),1);
    inpts=spec.length();
    cout << "Points in spectrum: " << inpts << endl;

    specleft=rspec(0,0)*(larmorM/1e6);
    df=(larmorM/1e6)*(rspec(inpts-1,0)-rspec(0,0))/(inpts-1);
    //specright=specleft+npts*df;

    cout << "Frequency step: " << df << " Hz\n";

    pars(size_t(AMP))=sum(spec);
    errs(size_t(AMP))=pars(size_t(AMP))/100.0;

    double mom=moment(spec);
    cout << "Midpoint: " << mom << endl;
    midp=(mom*df+specleft)*(1e6/larmorM);
  }
  else
    pars(size_t(AMP))=1.0;
  varywhich(point++)=AMP;
    
  //getargs(fullH,"Use full Hilbert space? ");

  double& qpole=pars(size_t(QPOLE));
  {
    istream& istr=getargs(qpole,"14N quadrupolar coupling (MHz)? ");
    if (qpole && !isfixed(istr)) varywhich(point++)=QPOLE; 
    qpole*=1e6;

    if (qpole && !isfixed(getargs(pars(size_t(QETA)),"Quadrupolar asymmetry? ",Within<double>(0.0,1.0))))
      varywhich(point++)=QETA;
  }

  bool arecol;
  getargs(arecol,"Assume interactions are collinear? ",true);

  {
    double& dipole=pars(size_t(DIPOLE));
    istream& istr=getargs(dipole,"XN dipolar coupling (inc. J anisotropy) (kHz)? ");
    if (dipole && !isfixed(istr)) varywhich(point++)=DIPOLE;
    dipole*=1e3;
  }

  if (!arecol && qpole) {
    if (!isfixed(getargs(pars(size_t(DTOQBETA)),"dipole->quad PAS beta angle? "))) varywhich(point++)=DTOQBETA;
    if (!isfixed(getargs(pars(size_t(DTOQGAMMA)),"dipole->quad PAS gamma angle? "))) varywhich(point++)=DTOQGAMMA;
    pars(size_t(DTOQBETA))*=M_PI/180.0;
    pars(size_t(DTOQGAMMA))*=M_PI/180.0;
  }

  larmorM=larmorH*gamma(Mtype)/gamma1H;
  cout << "Larmor frequency of 2H: " << larmorM/1e6 << " MHz\n";

  double& qpoleM=pars(size_t(QPOLEM));
  if (!isfixed(getargs(qpoleM,"M quadrupolar coupling (MHz)? "))) varywhich(point++)=QPOLEM; 
  qpoleM*=1e6;
      
  if (!isfixed(getargs(pars(size_t(QETAM)),"Quadrupolar asymmetry? ",Within<double>(0.0,1.0))))
    varywhich(point++)=QETAM;
      
  if (!arecol) {
    if (pars(size_t(QETAM))) {
      if (!isfixed(getargs(pars(size_t(QMTOQALPHA)),"quad M->quad PAS alpha angle? ")))
	varywhich(point++)=QMTOQALPHA;
    }
    if (!isfixed(getargs(pars(size_t(QMTOQBETA)),"quad M->quad PAS beta angle? ")))
      varywhich(point++)=QMTOQBETA;
    if (!isfixed(getargs(pars(size_t(QMTOQGAMMA)),"quad M->quad PAS gamma angle? ")))
      varywhich(point++)=QMTOQGAMMA;
    pars(size_t(QMTOQALPHA))*=M_PI/180.0;
    pars(size_t(QMTOQBETA))*=M_PI/180.0;
    pars(size_t(QMTOQGAMMA))*=M_PI/180.0;
  }
    
  if (!isfixed(getargs(pars(size_t(JMX)),"J_MX coupling (Hz)? ")))
    varywhich(point++)=JMX;

  if (havespec) {
    if (!isfixed(getargs(pars(size_t(MISO)),"Isotropic frequency (ppm)? ",midp)))
      varywhich(point++)=MISO;
  }
  else
    pars(size_t(MISO))=0.0; //ignore shift if not fitting
  
  double& CSA=pars(size_t(MCSA));
  //std::istream& istr=getargs(CSA,"Shift anisotropy (ppm)? ");
  CSA=0.0; //neglect 2H CSA
//   if (CSA) {
//     if (!isfixed(istr)) varywhich(point++)=MCSA;
    
//     if (!isfixed(getargs(pars(size_t(CSAETA)),"CSA asymmetry? ",Within<double>(0.0,1.0))))
//       varywhich(point++)=CSAETA;
//     if (!arecol) {
//       if (pars(size_t(CSAETA))) {
// 	if (!isfixed(getargs(pars(size_t(CSATOQALPHA)),"CSA->quad PAS alpha angle? ")))
// 	  varywhich(point++)=CSATOQALPHA;
//       }
//       if (!isfixed(getargs(pars(size_t(CSATOQBETA)),"CSA->quad PAS beta angle? ")))
// 	varywhich(point++)=CSATOQBETA;
//       if (!isfixed(getargs(pars(size_t(CSATOQGAMMA)),"CSA->quad PAS gamma angle? ")))
// 	varywhich(point++)=CSATOQGAMMA;
//       pars(size_t(CSATOQALPHA))*=M_PI/180.0;
//       pars(size_t(CSATOQBETA))*=M_PI/180.0;
//       pars(size_t(CSATOQGAMMA))*=M_PI/180.0;
//     }
//   }

  double& rotor_speed=pars(size_t(SPEED));
  getargs(rotor_speed,"Rotor speed (kHz)? "); rotor_speed*=1e3;

  if (havespec) {
    const int mins=int(fabs(inpts*df/rotor_speed)+1.0);
    getargs(nobs,"Sidebands to calculate? ",mins);
    sw=fabs(inpts*df);
  }
  else {
    getargs(nobs,"Sidebands to calculate? ");
    sw=fabs(nobs*rotor_speed);
  }

  getargs(nprops,"Gamma integration steps (multiple of sidebands)? ",2*nobs);
  if (nprops % nobs) throw Failed("Gamma steps must be multiple of sidebands");
  getargs(intsteps,"Integration steps in rotor cycle? ",100);

  getargs(zcw,"ZCW parameter? ");
  if (zcw==0) {
    getargs(default_orient.alpha,"Alpha? ");
    getargs(default_orient.beta,"Beta? ");
    getargs(default_orient.gamma,"Gamma? ");
    default_orient*=M_PI/180.0;
  }

  getargs(rphase,"Rotor phase (degrees)? ",0.0); rphase*=M_PI/180.0;

  if (!isfixed(getargs(pars(size_t(LW)),"(Lorentzian) linewidth (Hz)? ")))
    varywhich(point++)=LW;

  if (!havespec) {
    getargs(npts,"Number of points in spectrum (ideally 2^N)? ",128U);
    spec.create(npts);
  }
  else {
    const double dt=(1.0/rotor_speed)/nobs;
    const double tmax=1.0/fabs(df);
    //cout << "tmax: " << tmax*1e3 << " ms\n";
    
    double actpts=tmax/dt;
    //cout << "Points needed: " << actpts << endl;
    npts=int(0.5+actpts);
  }

  if (spec.length()<=point) {
    cerr << "Fewer data points than parameters!\n";
    return FAILED;
  }

  //getargs(tdomain,"Time domain calculation? ",false);

  bool dofit=false;
  if (havespec) getargs(dofit,"Perform fitting? ",true);
  bool calcbounds=dofit;
  //  if (!dofit) getargs(calcbounds,"Determine errors? ",false);

  double noise= 0.0;
  if (calcbounds) getargs(noise,"Noise level? ");

  const size_t NN=0;
  const size_t NM=1;

  spin_system sys(2,Ntype);
  sys(NM).isotope(Mtype);
  const bool isquad=(sys(NM).deg()>2);

  cout.precision(3);

  //cout << "T00\n" << T2(sys,NN,0,0) << endl;
  for (int m=-2;m<=2;m++) {
    NT2Q(m+2)=T2(sys,NN,2,m);
    cout << "m=" << m << endl << NT2Q(m+2) << endl;
    //    MT2Q(m+2)=T2(sys,NM,2,m);
  }
  if (isquad)
  MQ=spin_quadrupolar(sys,NM);
  
  Hzeeman=larmorN*diag_Iz(sys,NN);
  //  Hzeeman+=larmorM*diag_Iz(sys,NM);
  cout << "Zeeman Hamiltonian: " << Hzeeman << '\n';

  HdipM(0)=0.5*I(sys,NM,'-');
  HdipM(1)=sqrt(2.0/3.0)*I(sys,NM,'z');
  HdipM(2)=-0.5*I(sys,NM,'+');
 
  // if (fullH) {
  Hcs=diag_Iz(sys,NM);
  sigma0det=I(sys,NM,'+');
  cout << "sigma0:\n" << sigma0det << endl;
  //}

  char outname[140];
  getargs(outname,1,128,"Output filename? ",havespec ? fname : NULL);  
  if (havespec) strcat(outname,"Fit");

  rmatrix covar;

  if (calcbounds) {
    varywhich.resize(point);
    //cout << "Parameters being optimised: " << varywhich << endl;

    List<bool> isfixed(NPARS,true);
    isfixed(varywhich)=false;

    cout << "\n";
    for (int i=0;i<NPARS;i++) {
      cout << par_names[i] << ":\t" << pars(i) << "\t";
      if (isfixed(i))
	cout << "fixed\n";
      else
	cout << "+/- " << errs(i) << "\n";
    }
    cout << "\nVariable parameters: " << point << endl;
  }

  if (dofit) {
    double chisqr,ochisqr;

    LevMarqFit fitobj(ochisqr,pars,FitFunction(calc_spec),spec,varywhich,errs,noise);
    fitobj.asymmetric(true);
    fitobj.verbose(2);

    cout << "Initial chi^2: " << ochisqr << endl;
    const double chistop=1e-4*spec.length();
    cout << "Terminates when chi^2 varies by less than: " << chistop << endl;

    int iter_count=0;

    const char contops[]="crsa"; 

    bool nopause=false;

    const double tol=chistop/1e6;

    for (;;) {

      if (!nopause) {
	rspec(slice(),1)=fitobj.current_trial();
	write_matrix(outname,rspec,NULL,mxflag::block);
	//write_vector(outname,fitobj.current_trial());
	
	size_t opt;
	getargs(opt,"(c)ontinue, (r)un free, (s)top, (a)bort? ",contops,0);
	const char chosen=contops[opt];
	if (chosen=='s') break;
	switch (chosen) {
	case 'a': return 1;
	case 'r': nopause=true; break;
	}
     }

      cout << (iter_count++) << ": ";
      chisqr=fitobj.next();
      if (chisqr<ochisqr+tol) {
	if (ochisqr-chisqr<chistop) break;
	ochisqr=chisqr;
      }
    }

    chisqr=fitobj.final(covar,pars);

    cout << "Final chi^2: " << chisqr << endl;

    rspec(slice(),1)=fitobj.current_trial();
    error_filter(write_matrix(outname,rspec));
  }
  else {
    if (havespec) {
      List<double> tspec(spec.length());

      calc_spec(tspec,pars);
      error_filter(write_vector(outname,tspec));
      cout << "(Non-normalised) chi^2: " << norm(tspec-spec) << endl;
    }
    else {
      calc_spec(spec,pars);
      //      error_filter(write_vector(outname,spec));
      char fname[128];
      List<complex> cspec(spec);
      sprintf(fname,"%s.spe",outname);
      write_simpson(fname,cspec,sw,!tdomain);      
    }

    if (calcbounds) {
      List<double> weights(npts,noise);
      covariance(covar,FitFunction(calc_spec),pars,varywhich,errs,weights);    
    }
  }
  
  if (isdefined(covar)) {
    int i;
    List<double> bounds(point);

    cout << endl;
    for (i=0;i<point;i++) {
      const size_t apar=varywhich(i);
      bounds(i)=sqrt(covar(i,i));
      cout << i << "  " << par_names[apar] << ": " << pars(apar) << " +/- " << bounds(i) << "\n";
    }        

    cout << endl;
    for (i=0;i<point;i++) {
      const size_t apari=varywhich(i);
      for (int j=i+1;j<point;j++) {
	cout << "cor(" << par_names[apari] << "," << par_names[varywhich(j)] << ")=" << covar(i,j)/(bounds(i)*bounds(j)) << "\n";
      }
    }
  }
  }
  catch (MatrixException& exc) {
    cerr << exc << endl;
    return 1;
  }
  return 0;
}
