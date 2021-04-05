#include <cstdlib>
#include "MAS.h"
#include "NMR.h"
#include "Propagation.h"
#include "ttyio.h"
#include "timer.h"

/* This is the program used to test inhomogeneous propagation in PNMRS article */

// testinhomo 1 - 500 n - 24 n 8 1 100
// testinhomo 1 - 500 n - 24 y 256 8 1

using namespace std;
using namespace libcmatrix;

timer<> stopwatch;

int errcount=0;
List<double> refreal;
List<complex> refcomplex;
double period=0.0;

template<class T> void print_result_(const char *header,const BaseList<T>& FID,List<T>& ref)
{
  if (ref.empty()) {
    cout << header << FID << endl;
    ref=FID;
  }
  else {
    const bool isbad=(sqrt(norm(FID-ref))>0.01*FID.length());
    if (isbad) errcount++;
    cout << header << (isbad ? "(BAD) " : "(OK) ") << FID << endl;
  }
}

void print_result(const char *header,const BaseList<complex>& FID)
{
  print_result_(header,FID,refcomplex);
}

void print_result(const char *header,const BaseList<double>& FID)
{
  print_result_(header,FID,refreal);
}

template<class T> void split(List<Matrix<T> >& dest,const Matrix<T>& source,const ListList<size_t>& block_str,int offset,int verbose =0)
{
  const int n=block_str.size()-abs(offset);
  dest.create(n);
  if (offset<0) {
    for (int i=0;i<n;i++) {
      dest(i)=source(block_str(i-offset),block_str(i));
      if (verbose) cout << "Block " << i << '\n' << dest(i) << '\n';
    }
  }
  else {
    for (int i=0;i<n;i++) {
      dest(i)=source(block_str(i),block_str(i+offset));
      if (verbose) cout << "Block " << i << '\n' << dest(i) << '\n';
    }
  }
}

void (*addfunc)(BaseList<complex> &,complex,double);

void print_usage(double overhead,int rep_times)
{
  double ttake=stopwatch();
  cout << "Time taken: " << ttake << " s\n";
  overhead*=1e3/rep_times;
  ttake*=1e3/rep_times;
  cout << "Per calculation: (" << overhead << ") " << (overhead+ttake) << " ms\n";
}

int n;
int spins;
double sw;
bool verbose=false;

template<class T> void readout(BaseList<complex>& Output, SpectrumIterator<T>& obj)
{
  double freq;
  T amp;
  while (obj(amp,freq)) {
    addfunc(Output,amp,freq/sw);
    if (verbose) cout << amp << " of " << freq << endl;
  }
}

template<class F> void docalc2(List<complex>& Output,F& obj)
{
  if (verbose) cout << obj << endl;
  readout(Output,obj);
}

template<class F> void dump_add(F& obj)
{
  const Matrix< typename F::transition_type>& amps(obj.amplitudes());

  const size_t dim=amps.rows();
  for (size_t i=0;i<dim;i++) {
    for (size_t j=0;j<dim;j++) {
      if (amps(i,j)!=0.0) {
	cout << "Transition " << i << "," << j << '\n';
	cout << "Frequencies: ";
	for (int k=0;k<n;k++)
	  cout << obj.frequency(i,j,k) << ' ';
	cout << endl;
	List< typename F::amplitude_type > samps;
	obj.add(samps,i,j);
	cout << "Amplitudes: " << samps << endl;
      }
    }
  }
}

template<class F> void docalc(List<complex>& Output,F& obj,const rmatrix& H, DynamicPhase& phasecalc, const cmatrix& sigma0, const cmatrix& detect)
{
  obj.set_H(H);
  obj.set_phases(phasecalc);
  obj.observe(sigma0,detect);
  docalc2(Output,obj);

  if (verbose)
    dump_add(obj);
}

template<class F> void dosubcalc(List<complex>& Output,F& obj,const BaseList<rmatrix>& Hs, DynamicPhase& phasecalc, const BaseList<cmatrix>& sigma0s, const BaseList<cmatrix>& detects)
{
  obj.set_phases(phasecalc);
  obj.set_H('C',Hs(0));
  const int blks=Hs.length();
  for (int blk=1;blk<blks;blk++) {
    obj.shuffle('C');
    obj.set_H('C',Hs(blk));
    obj.observe(sigma0s(blk-1),detects(blk-1));
    docalc2(Output,obj);
  }
}

template<class F> void dosubcalc(List<complex>& Output,F& obj,const BaseList<rmatrix>& Hs, DynamicPhase& phasecalc, const BaseList<cmatrix>& sigma0dets)
{
  obj.set_H('C',Hs(0));
  obj.set_phases(phasecalc);
  const int blks=Hs.length();
  for (int blk=1;blk<blks;blk++) {
    obj.shuffle('C');
    obj.set_H('C',Hs(blk));
    obj.observe(sigma0dets(blk-1));
    docalc2(Output,obj);
  }
}

template<class F> void docalc(List<complex>& Output, F& obj,const rmatrix& H, DynamicPhase& phasecalc, const cmatrix& sigma0det)
{
  obj.set_H(H);
  obj.set_phases(phasecalc);
  obj.observe(sigma0det);
  docalc2(Output,obj);

  if (verbose)
    dump_add(obj);
}

template<class F> void docalc(List<complex>& Output, F& obj,const cmatrix& Us, const cmatrix& sigma0det)
{
  obj.set_Us(Us,period);
  obj.observe(sigma0det);
  docalc2(Output,obj);

  if (verbose)
    dump_add(obj);
}

template<class F> void docalc(List<complex>& Output, F& obj,const cmatrix& Us, const cmatrix& sigma0, const cmatrix& detect)
{
  obj.set_Us(Us,period);
  obj.observe(sigma0,detect);
  docalc2(Output,obj);

  if (verbose)
    dump_add(obj);
}

template<class F,class M> void dohomocalc(BaseList<complex>& Output,F& obj,BaseList<cmatrix>& shiftedUs,double period,const M& sigma0, bool verbose =true)
{
  obj.set_Us(shiftedUs,period);
  obj.observe(sigma0);
  if (verbose) cout << obj << endl; 
  readout(Output,obj);
}

template<class F,class M> void dohomocalc(BaseList<complex>& Output,F& obj,BaseList<cmatrix>& shiftedUs,double period,const M& sigma0, const M& det,bool verbose =true)
{
  obj.set_Us(shiftedUs,period);
  obj.observe(sigma0,det);
  if (verbose) cout << obj << endl; 
  readout(Output,obj);
}

void add_signal(BaseList<complex> &FID,complex a,double w)
{
  if (w==0) {
    FID+=a;
    return;
  }

  const complex fac=expi(2*M_PI*w);

  const int npts=FID.length();
  complex* data=FID.vector();
  for (int i=0;i<npts;i++) {
    data[i]+=a;
    a*=fac;
  }
}

int speclength=0;

inline void add_spectrum(BaseList<complex> &spec,complex a,double w)
{
  const int bin=int(0.5+w*speclength) % speclength;
  spec(bin)+=a;
}

int main(int argc,const char *argv[])
{
  int gn,r,s,j,k,cyc;
  int count=1;

  try {

  spins=getint(argc,argv,count,"Number of spins? ",1,2,1);

  const double approx_sw=getfloat(argc,argv,count,"Spectral width? ",4000.0);

  const double rotor_speed=getfloat(argc,argv,count,"Rotor speed (Hz) ? ");
  n=int(approx_sw/rotor_speed);
  if (n<1) n=1;

  sw=n*rotor_speed;
  
  cout << "Observations: " << n << endl;
  cout << "Spectral width: " << sw << endl;

  const bool isherm=getlogical(argc,argv,count,"Hermitian sigma0 and detect? ");
  const int maxdosub= isherm ? 0 : 1;
  const bool samedetinit=getlogical(argc,argv,count,"Matching excitation and detection? ");

  const double iso=getfloat(argc,argv,count,"Isotropic frequency? ",200.0);

  period=1.0/rotor_speed;
  const double stept=period/n;

  const int gamma_steps=getint(argc,argv,count,"Gamma steps? ");

  Euler PAS_to_RF(0,M_PI/3,0);
  if (gamma_steps==1)
    PAS_to_RF.gamma=getfloat(argc,argv,count,"Gamma? ",0.0);
  else {
    if (gamma_steps % n) {
      cerr << "Gamma steps must be a multiple of observations per rotor period.\n";
      exit(1);
    }
  }
  
  const bool asfreq=getlogical(argc,argv,count,"Histogram? ");
  if (asfreq)
    speclength=getint(argc,argv,count,"Number of points in spectrum? ",256);

  const int cycles=getint(argc,argv,count,"Number of periods? ",1);

  addfunc= asfreq ? add_spectrum : add_signal;

  const int rep_times=getint(argc,argv,count,"Repeat times? ",1);
  verbose=(rep_times==1);
  int loop;

  const int intsteps=gamma_steps/n;

  //Do homogeneous calculation as check

  const int isteps=getint(argc,argv,count,"Number of integration steps (0 for no homo check)? ",100);

  spin_system sys(spins);

  rmatrix H,HCS;
  space_T A_PAS;
  if (spins==1) {
    H=real(I(sys,0,'z'));
    A_PAS=spatial_tensor(iso,2000,0.0);
  }
  else {
    H=real(spin_dipolar(sys,0,1));
    A_PAS=spatial_tensor(2000);
    cout.precision(4);
  }

  space_T A_RF=rotate(A_PAS,PAS_to_RF);
  cout << "Interaction tensor\n" << A_RF << endl;

  const int dim=H.rows();

  const cmatrix sigma0=I(sys,0,(isherm || !samedetinit)? 'x' : '+');
  const rmatrix rsigma0=real(sigma0);
  cmatrix detect;
  if (samedetinit)
    detect=conj_transpose(sigma0);
  else
    detect=I(sys,spins-1,isherm ? 'x' : '-');
  const rmatrix rdetect=real(detect);
  
  cout << "sigma0\n" << rsigma0 << endl;
  cout << "detect\n" << rdetect << endl;

  const ListList<size_t> block_str=find_blocks(diag_Fz(sys),0.0);
  const int blks=block_str.size();
  cout << "Block structure: " << block_str << '\n';

  List<cmatrix> subsigma0s,subdetects;
  List<rmatrix> subHs;
  if (maxdosub) {
    split(subHs,H,block_str,0,verbose);
    split(subsigma0s,sigma0,block_str,1,verbose);
    split(subdetects,detect,block_str,-1,verbose);
  }

  const int nphases= (gamma_steps>1) ? gamma_steps : n;
  DynamicPhase phasecalc(rotor_speed,0.0,nphases,A_RF);
  cout << "Obs: " << phasecalc.observations() << '\n';

  List<double> phase(nphases);
  List<double> phaseshift(nphases);

  cout << "B values\n";
  for (j=-2;j<=2;j++) cout << j << ": " << phasecalc.component(j) << endl;
  //  const double B0=real(phasecalc.component(0));

  cmatrix U,S,R,sigma;

  const int total_steps=n*cycles;

  List<complex> FID(total_steps);
  List<double> FIDr(total_steps);
  List<complex> Spec(speclength);

  List<complex>& Output = asfreq ? Spec : FID;

  stopwatch.reset();
  int red_loop=rep_times/40;
  if (red_loop<1) red_loop=1;

  for (loop=red_loop;loop--;) {
    FID=complex(0.0);
    
    for (gn=0;gn<gamma_steps;gn++) {
      if (gamma_steps>1)
	phasecalc.rotor_phase((2*M_PI*gn)/gamma_steps);
  
      cmatrix U;

      FID(0)+=trace_multiply(sigma0,detect);
      for (j=1;j<total_steps;j++) {
	const double dphase=phasecalc((j-1)*stept,j*stept);
	U&=hermitian_expi(H,dphase);
	//	if (rep_times==1) 
	//  cout << j*stept << " (" << gn << "): " << dphase << '\n' << U;
	unitary_simtrans(sigma,sigma0,U);
	FID(j)+=trace_multiply(sigma,detect);
      }
    }
    if (gamma_steps>1) FID/=gamma_steps;
  }
  print_result("Very simple propagation: ",FID);
  print_usage(0,red_loop);  

  phasecalc.rotor_phase(0.0);

  cmatrix detectt;
  rmatrix D;

  List<double> eigs;

  stopwatch.reset();
  
  // all this is independent of alpha, beta (usually)
  // therefore *not* included in timings

  hermitian_eigensystem(D,eigs,H);
  cout << "Eigenvalues: " << eigs << endl;
  cout << "Transformation matrix:\n" << D << endl;

  const double period_iso=phasecalc.isotropic(nphases);
    
  unitary_isimtrans(S,sigma0,D);
  cout << "sigma0':\n" << S << endl;
  unitary_isimtrans(detectt,detect,D);
  cout << "detect':\n" << detectt << endl;
  detectt.transpose();
  S.emultiply(detectt);
  cout << "Transition probs:\n" << S << endl;
   
  for (loop=rep_times;loop--;)
    for (gn=nphases;gn--;) phase(gn)=phasecalc(gn);

  for (gn=nphases;gn--;) phaseshift(gn)=phasecalc(gn+1);

  const double phase_calc2=stopwatch();
  cout << "Preliminaries: "; print_usage(0,rep_times);
  cout << "Phase table: " << phase << endl;

  List<complex> phasebuf2(gamma_steps);
  complex *const vecbuf=phasebuf2.vector();
  
  // in basis propagation
  stopwatch.reset();
  for (loop=rep_times;loop--;) {
    FID=complex(0);
    for (int r=dim;r--;) {
      for (int s=dim;s--;) {
	const complex &amp=S(r,s);
	if (norm(amp)>1e-20) {
	  //const complex v=conj_multiply(evolph(s),evolph(r));
	  const double diff=eigs(r)-eigs(s);
	  const complex v=expi(diff*period_iso);
	  cout << "cycle evol: " << v << endl;
	  for (gn=gamma_steps;gn--;)
	    vecbuf[gn]=expi(diff*phase(gn));
	  //vecbuf[gn+gamma_steps]=isofac*(vecbuf[gn]=expi(diff*phase(gn)));

	  for (j=0;j<n;j++) {
	    complex p=0;
	    if (gamma_steps>1) {
	      const int off=j*intsteps;
	      complex *vecoff=vecbuf+off;
	      for (gn=gamma_steps-off;gn--;)
	    	mla_conj(p,vecoff[gn],vecbuf[gn]);
	      vecoff-=gamma_steps;
	      complex p2=0;
	      for (gn=gamma_steps-off;gn<gamma_steps;gn++)
	    	mla_conj(p2,vecoff[gn],vecbuf[gn]);
	      mla(p,v,p2);
	    }
	    else
	      p=expi(diff*phase(j));
	    p*=amp;
	    k=j;
	    for (cyc=cycles;cyc--;) {
	      FID(k)+=p;
	      p*=v;
	      k+=n;
	    }
	  }
	}
      }
    }
    FID/=gamma_steps;
  }

  //rmatrix rS=real(S);

  print_result("Inhomogeneous in basis: ",FID);
  print_usage(phase_calc2,rep_times);

  cmatrix Udiags(nphases,H.cols());

  int maxdotype=0;
  if (spins==1) {
    maxdotype=1;
    DiagonalSpinningHamiltonian Hspin(rotor_speed,0.0);
    Hspin.add(A_RF,diag_Iz(sys,0));
    cout << "Hspin\n" << Hspin;
    cout << "Hspin(tr/5): " << Hspin(period/5.0) << '\n';

    DiagonalInhomogeneousPropagator Ugen(Hspin);
    const double dt=period/nphases;
    
    for (j=0;j<nphases;j++) {
      BaseList<complex> Urowj(Udiags.row(j));
      Ugen(Urowj,0.0,(j+1)*dt);
    }
    cout << "Propagators\n" << Udiags;
  }

  for (int dotype=0;dotype<=maxdotype;dotype++) {
    for (int dosub=0;dosub<=maxdosub;dosub++) {
      if (dosub && dotype) continue; //don't bother with blocking 
      
      if (gamma_steps>1) {
	GammaInhomogeneousFID obj(n);
	
	stopwatch.reset();
	for (loop=rep_times;loop--;) {
	  if (dotype==0)
	    obj.set_phases(phasecalc);
	  
	  if (dosub) {
	    FID=complex(0.0);
	    obj.set_H('C',subHs(0));
	    for (int blk=1;blk<blks;blk++) {
	      obj.shuffle('C');
	      obj.set_H('C',subHs(blk));
	      if (samedetinit)
		obj.add_FID(FID,1.0,subsigma0s(blk-1));
	      else
		obj.add_FID(FID,1.0,subsigma0s(blk-1),subdetects(blk-1));
	    }
	  }
	  else {
	    if (dotype==0)
	      obj.set_H(H);
	    else
	      obj.set_Us(Udiags);

	    if (verbose) 
	      cout << obj << endl;
	    
	    if (isherm) {
	      if (samedetinit)
		FIDr=obj.FID_hermitian(total_steps,rsigma0);
	      else
		FIDr=obj.FID_hermitian(total_steps,rsigma0,rdetect);
	    }
	    else {
	      if (samedetinit)
		FID=obj.FID(total_steps,rsigma0);
	      else
		FID=obj.FID(total_steps,rsigma0,rdetect);
	    }
	  }
	}
      }
      else {
	
	InhomogeneousFID obj;
	
	stopwatch.reset();
	for (loop=rep_times;loop--;) {

	  if (isherm)
	    FIDr=0.0;
	  else
	    FID=complex(0.0);

	  if (dotype==0) {
	    obj.set_H(H);
	    obj.set_phases(phasecalc);
	  }
	  else
	    obj.set_Us(Udiags);

	  if (verbose)
	    cout << obj << endl;

	  if (isherm) {
	    if (samedetinit)
	      FIDr=obj.FID_hermitian(total_steps,rsigma0);
	    else
	      FIDr=obj.FID_hermitian(total_steps,rsigma0,rdetect);
	  }
	  else {
	    if (samedetinit)
	      FID=obj.FID(total_steps,rsigma0);
	    else
	      FID=obj.FID(total_steps,rsigma0,rdetect);
	  }
	}
      }

      cout << (dotype ? "via DiagonalSpinningHamiltonian\n" : "via DynamicPhase\n");
      const char *lstr=dosub ? "(Gamma)InhomogeneousFID (sub): " : "(Gamma)InhomogenousFID: ";
      if (isherm)
	print_result(lstr,FIDr);
      else
	print_result(lstr,FID);
      print_usage(phase_calc2,rep_times);
    }
  }

  if (gamma_steps>1) {

    stopwatch.reset();
    for (loop=rep_times;loop--;) {
      for (int gn=nphases;gn--;) {
	phase(gn)=phasecalc.anisotropic(gn);
	phaseshift(gn)=phasecalc.anisotropic(gn+1);
      }
    }
    const double phase_calc=stopwatch();
    cout << "Preliminaries: "; print_usage(0,rep_times);

    List<complex> ft_facs(gamma_steps);
    ft_create_table(ft_facs,FT_BACKWARD);

    const int usefft=ispowerof2(gamma_steps);
    cout << "Fast FFT: " << (usefft ? "Yes" : "No") << endl;

    List<complex> sideband(gamma_steps);
    List<complex> scr(gamma_steps);

    List<complex> & store = usefft ? sideband : scr;

    stopwatch.reset();
    const double scale=1.0/(gamma_steps*gamma_steps);
    for (loop=rep_times;loop--;) {
      FID=complex(0,0);
            
      for (r=0;r<dim;r++) {
	for (s=0;s<dim;s++) {
	  complex amp=S(r,s);
	  if (norm(amp)>1e-20) {
	    const double diff=eigs(r)-eigs(s);
	    const double viso=phasecalc.component0()*diff;
	    if (diff) {
	      amp*=scale;
	      for (j=gamma_steps;j--;) store(j)=expi(-diff*phase(j));
	      if (usefft)
		fft_ip(sideband,FT_BACKWARD);
	      else
		ft(sideband,store,ft_facs);
	      for (j=0;j<n;j++) {
		double totamp=0;
		for (k=j;k<gamma_steps;k+=n) totamp+=norm(sideband(k));
		double f=-viso+j*rotor_speed;
		if (verbose) 
		  cout << "Adding " << (amp*totamp) << " of " << f << '\n';
		addfunc(Output,amp*totamp,f/sw);
	      }
	    }
	    else
	      addfunc(Output,amp,viso/sw);
	  }
	}
      }
    }
    print_result("Inhomogeneous FT: ",Output);
    print_usage(phase_calc,rep_times);

    GammaInhomogeneousSpectrumED objED(n);
    GammaInhomogeneousSpectrum obj(n);

    for (int dotype=0;dotype<=maxdotype;dotype++) {
      for (int dosub=0;dosub<=maxdosub;dosub++) {
	if (dosub && dotype)
	  continue;

	stopwatch.reset();
	for (loop=rep_times;loop--;) {
	  Output=complex(0.0);
	  
	  if (dosub) {
	    if (samedetinit)
	      dosubcalc(Output,objED,subHs,phasecalc,subsigma0s);
	    else
	      dosubcalc(Output,obj,subHs,phasecalc,subsigma0s,subdetects);
	  }
	  else {
	    if (dotype) {
	      if (samedetinit)
		docalc(Output,objED,Udiags,sigma0);
	      else
		docalc(Output,obj,Udiags,sigma0,detect);
	    }
	    else {
	      if (samedetinit)
		docalc(Output,objED,H,phasecalc,sigma0);
	      else
		docalc(Output,obj,H,phasecalc,sigma0,detect);
	    }
	  }
	}
	cout << (dotype ? "via DiagonalSpinningHamiltonian\n" : "via DynamicPhase\n");
	print_result(dosub ? "GammaInhomogeneousSpectrum(ED) (sub): " : "GammaInhomogeneousSpectrum(ED): ",Output);
	print_usage(phase_calc,rep_times);
      }
    }
  }
  if (gamma_steps==1) {

    InhomogeneousSpectrum obj;
    InhomogeneousSpectrumED objED;

    for (int dotype=0;dotype<=maxdotype;dotype++) {
      for (int dosub=0;dosub<=maxdosub;dosub++) {
	if (dosub && dotype)
	  continue;

	stopwatch.reset();
	for (loop=rep_times;loop--;) {
	  
	  Output=complex(0.0);
	  
	  if (dosub) {
	    if (samedetinit)
	      dosubcalc(Output,objED,subHs,phasecalc,subsigma0s);
	    else
	      dosubcalc(Output,obj,subHs,phasecalc,subsigma0s,subdetects);
	  }
	  else {
	    if (dotype) {
	      if (samedetinit)
		docalc(Output,objED,Udiags,sigma0);
	      else
		docalc(Output,obj,Udiags,sigma0,detect);
	    }
	    else {
	      if (samedetinit)
		docalc(Output,objED,H,phasecalc,sigma0);
	      else
		docalc(Output,obj,H,phasecalc,sigma0,detect);
	    }
	  }
	}
	cout << (dotype ? "via DiagonalSpinningHamiltonian\n" : "via DynamicPhase\n");
	print_result(dosub ? "InhomogeneousSpectrum (sub): " : "InhomogeneousSpectrum: ",Output);
	print_usage(phase_calc2,rep_times);
      }
    }
  }

  if (isteps) {

    const int usteps=(gamma_steps==1) ? n : gamma_steps;
    const double dt= period/usteps;
    
    RealSpinningHamiltonian H_cycle(rotor_speed,0.0);
    H_cycle.add(A_RF,H);
    cout << H_cycle;
    cout << "Hspin(tr/5):\n" << H_cycle(period/5.0);

    const double intdt=period/isteps;

    //rotor_phase always zero (gamma included in A_RF)
    HomogeneousPropagator<rmatrix> Ugen(H_cycle,intdt);
    List<cmatrix> shiftedUs(usteps);

    cmatrix U;

    for (j=0;j<usteps;j++) {
      Ugen(U,j*dt,(j+1)*dt);
      if (j)
	multiply(shiftedUs(j),U,shiftedUs(j-1));
      else
	shiftedUs(j)=U;
    }
    cout << "U(0)\n" << shiftedUs(0);

    if (asfreq) { 
      Output=complex(0.0);
      
      //don't bother with samedetinit - this is only a check
      if (gamma_steps>1) {
	GammaPeriodicSpectrum pobj(n);
	dohomocalc(Output,pobj,shiftedUs,period,rsigma0,rdetect);
      }
      else {
	PeriodicSpectrum pobj;
	dohomocalc(Output,pobj,shiftedUs,period,rsigma0,rdetect);
      }
      print_result("(Gamma)PeriodicSpectrum: ",Output);
    }
    else {
      if (gamma_steps>1) {
	GammaPeriodicFID obj(n);
	
	obj.set_Us(shiftedUs);
	if (verbose) cout << obj << '\n';
	if (isherm)
	  Output=obj.FID_hermitian(total_steps,sigma0,detect);
	else
	  Output=obj.FID(total_steps,sigma0,detect);
      }
      else {
	PeriodicFID obj;

	obj.set_Us(shiftedUs);
	if (verbose) cout << obj << '\n';
	if (isherm)
	  Output=obj.FID_hermitian(total_steps,sigma0,detect);
	else
	  Output=obj.FID(total_steps,sigma0,detect);
      }
      matrix_traits<double>::allocator::print(std::cout);
      print_result("(Gamma)PeriodicFID: ",Output);
    }
  }//end of homo check



  } catch (MatrixException& exc) {
    cerr << exc << endl;
    return 1;
  }

  cout << "Failures: " << errcount << '\n';
  return errcount;
}
