#include <cstdlib>
#include "MAS.h"
#include "NMR.h"
#include "ttyio.h"
#include "timer.h"
#include "Propagation.h"
#include "AsynchronousPropagation.h"
//#include "Floquet.h"

/* This is the program used to test homogeneous propagation for PNMRS article */
// testhomo <spins> <SW> <rspeed> <isotropic> <gsteps> <isteps> <samedi> 
// testhomo 1 - 500 - 24 - n - n 8 1
// testhomo 1 - 500 - 1 0 - y - n 2 1
// testhomo 2 - 500 - 1 0 - y - n n 2 1
// testhomo 1 - 500 - 1 3 - y - n 3 1

using namespace std;
using namespace libcmatrix;

timer<> stopwatch;

int verbose_level=1;
int errcount=0;
List<double> refreal;
List<complex> refcomplex;

template<class T> void print_result_(const char *header,const BaseList<T>& FID,List<T>& ref,bool asfreq)
{
  if (asfreq) {
    cout << header << FID << endl;
    return;
  }
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

void print_result(const char *header,const BaseList<complex>& FID,bool asfreq =false)
{
  print_result_(header,FID,refcomplex,asfreq);
}

void print_result(const char *header,const BaseList<double>& FID, bool asfreq =false)
{
  print_result_(header,FID,refreal,asfreq);
}

inline void mla(cmatrix &d,const cmatrix &a,const cmatrix &b)
{
  BaseList<complex> dr=d.row();
  mla(dr,a.row(),b.row());
}

void split(List<cmatrix>& dest,const cmatrix& source,const ListList<size_t>& block_str,int offset,int verbose =0)
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

template<class T> void dump_table(const BaseList<T>& tab)
{
  const int n=tab.length();
  for (int j=0;j<n;j++)
    cout << j << '\n' << tab(j);
  cout << '\n';
}

void print_usage(double overhead,int rep_times)
{
  double ttake=stopwatch();
  cout << "Time taken: " << ttake << " s\n";
  overhead*=1e3/rep_times;
  ttake*=1e3/rep_times;
  cout << "Per calculation: (" << overhead << ") " << (overhead+ttake) << " ms\n";
}

int spins;
int n;

double sw;
void (*addfunc)(BaseList<complex>,complex,double);

inline double norm(double x) { return x*x; }

template<class T> void readout(BaseList<complex>& Output,SpectrumIterator<T>& obj,bool verbose)
{
  T famp;
  double ffreq;
  
  while (obj(famp,ffreq)) {
    if (norm(famp)>1e-12) {
      if (verbose) cout << famp << " of " << ffreq << endl;
      addfunc(Output,famp,ffreq/sw);
    }
  }
}

template<class F> void docalc(BaseList<complex>& Output,F& obj,BaseList<cmatrix>& shiftedUs,double period,const cmatrix& sigma0, bool verbose)
{
  obj.set_Us(shiftedUs,period);
  obj.observe(sigma0);
  if (verbose) cout << obj << endl; 
  readout(Output,obj,verbose);

  if (spins==1) {
    cout << "(0,1) Frequencies: ";
    for (int i=0;i<n;i++) cout << obj.frequency(0,1,i) << ' ';
    cout << endl;
    typedef typename F::amplitude_type amp_t;
    List<amp_t> samps;
    obj.add(samps,amp_t(1.0),0,1);
    cout << "(0,1) Amplitudes: " << samps << endl;
  }
}

template<class F> void docalc(BaseList<complex>& Output,F& obj,BaseList< List<cmatrix> >& subUs,double period,const BaseList<cmatrix>& sigma0s, bool verbose)
{
  const int blks=subUs.length();

  for (size_t blk=0;blk<blks;blk++) {
    obj.shuffle('C');
    obj.set_Us('C',subUs(blk),period);

    if (blk) {
      obj.observe(sigma0s(blk-1));
      if (verbose) cout << obj << endl; 
      readout(Output,obj,verbose);
    }
  }
}

template<class F> void docalc(BaseList<complex>& Output,F& obj,BaseList< List<cmatrix> >& subUs,double period,const BaseList<cmatrix>& sigma0s, const BaseList<cmatrix>& detects,bool verbose)
{
  const int blks=subUs.length();

  for (size_t blk=0;blk<blks;blk++) {
    obj.shuffle('C');
    obj.set_Us('C',subUs(blk),period);

    if (blk) {
      obj.observe(sigma0s(blk-1),detects(blk-1));
      if (verbose) cout << obj << endl; 
      readout(Output,obj,verbose);
    }
  }
}

// template<class F,class TH,class TM> void dofloquetspectrum(BaseList<complex>& Output,F& obj,const TH& FHam,int N,const TM& sigma0,const TM& detect,bool verbose)
// {
//   obj.set_H(FHam,N);
//   obj.observe(sigma0,detect);
//   if (verbose) cout << obj << endl; 
//   readout(Output,obj,verbose);
// }

template<class F> void docalc(BaseList<complex>& Output,F& obj,BaseList<cmatrix>& shiftedUs,double period,const cmatrix& sigma0, const cmatrix& det,bool verbose)
{
  obj.set_Us(shiftedUs,period);
  obj.observe(sigma0,det);
  if (verbose) cout << obj << endl; 
  readout(Output,obj,verbose);
}

inline void add_spectrum(BaseList<complex> spec,complex a,double w)
{
  static const int speclength(spec.length());
  if (w<0.0)
    w+=1.0;
  const int bin=int(0.5+w*speclength) % speclength;
  if ((bin<0) || (bin>=speclength))
    throw BadIndex("add_spectrum");
  spec(bin)+=a;
}

inline void add_FID_complex(BaseList<complex> FID,complex a,double w) 
{ add_FID(FID,a,w); }

int main(int argc,const char *argv[])
{
  int i,gn,r,s,j,k,cyc;
  int count=1;

  try {

  spins=getint(argc,argv,count,"Number of spins? ",1,3,1);

  const double approx_sw=getfloat(argc,argv,count,"Spectral width? ",4000.0);

  const double rotor_speed=getfloat(argc,argv,count,"Rotor speed (Hz) ? ");
  n=int(approx_sw/rotor_speed);
  if (n<1) n=1;

  sw=n*rotor_speed;
  
  cout << "Observations: " << n << endl;
  cout << "Spectral width: " << sw << endl;

  spin_system sys(spins);

  cmatrix sigma;
  rmatrix H,Hdip12;
  List<double> HCS;
  space_T A_PAS;
  
  const double iso=getfloat(argc,argv,count,"Isotropic frequency? ",200.0);

  if (spins==1) {
    H=real(I(sys,0,'z'));
    A_PAS=spatial_tensor(iso,2000,0.5);
    //A_PAS=spatial_tensor(iso,0.0,0.0);
  }
  else {
    H=real(spin_dipolar(sys,0,1));
    HCS=iso*diag_Iz(sys,0);
    A_PAS=spatial_tensor(2000);
    if (spins>2) Hdip12=real(spin_dipolar(sys,1,2));
    cout.precision(4);
  }
  const List<double> Fz=diag_Fz(sys);
  const int dim=H.rows();
  const ListList<size_t> block_str=find_blocks(Fz,0.0);
  const int blks=block_str.size();
  cout << "Block structure: " << block_str << '\n';

  const double period=1.0/rotor_speed;

  const double stept=period/n;

  const int gamma_steps=getint(argc,argv,count,"Gamma steps? ");
  if (gamma_steps<1) throw InvalidParameter("Gamma steps must be >0");

  Euler PAS_to_RF(0,M_PI/4,0);
  Euler PAS_to_RF2(M_PI/3,M_PI/3,M_PI/3);

  if (gamma_steps==1)
    PAS_to_RF.gamma=getfloat(argc,argv,count,"Gamma? ",0.0)*M_PI/180.0;
  else {
    if (gamma_steps % n) {
      cerr << "Gamma steps must be a multiple of observations per rotor period.\n";
      exit(1);
    }
  }
  
  space_T A_RF=rotate(A_PAS,PAS_to_RF);
  cout << "A_RF\n" << A_RF << endl;

  space_T A_RF2=rotate(A_PAS,PAS_to_RF2);
  if (spins>2) cout << "A_RF2" << A_RF2 << endl;

  const int isteps=getint(argc,argv,count,"Number of integration steps? ",100);

  const double jitter=getfloat(argc,argv,count,"Timing jitter (ns)? ",0.0)*1e-9;
  IntervalSampler sampler(jitter);

  const bool samedetinit=(gamma_steps>1) ? getlogical(argc,argv,count,"Same detection and initial? ") : false;
  const bool isdiag=(spins==1) ? false : getlogical(argc,argv,count,"Diagonal operators? ");
  const bool isherm=isdiag ? true : getlogical(argc,argv,count,"Hermitian detection? ");
  //  const int floquet_N=getint(argc,argv,count,"Floquet order (0 to skip)? ",n/2);

  const int asfreq=getlogical(argc,argv,count,"Histogram? ");
  const int speclength = asfreq ? getint(argc,argv,count,"Number of points in spectrum? ",1) : 0;

  const int cycles=getint(argc,argv,count,"Number of periods? ",1);

  addfunc = asfreq ? add_spectrum : add_FID_complex;

  const int total_steps=n*cycles;

  List<complex> FID(total_steps);
  List<double> FIDr(total_steps);
  List<complex> Spec(speclength);

  List<complex>& Output = asfreq ? Spec : FID;

  const int rep_times=getint(argc,argv,count,"Repeat times? ",1);
  const bool verbose=1;//(rep_times==1);
  int loop;
  
  List<double> dsigma0,ddetect;
  cmatrix sigma0,detect;

  if (isdiag) {
    dsigma0=diag_Iz(sys,0);
    full(sigma0,dsigma0);
    ddetect=diag_Iz(sys,samedetinit ? 0 : 1);
    full(detect,ddetect);
  }
  else {
    sigma0=I(sys,0,(isherm || !samedetinit) ? 'x' : '+');
    if (samedetinit)
      detect=conj_transpose(sigma0);
    else
      detect=I(sys,0,isherm ? 'y' : '-');
  }
  
  const bool dosub=!isherm;
  const int maxdosub= dosub ? 1 : 0;

  cout << "Initial density matrix\n" << sigma0 << endl;
  cout << "Detection operator\n" << detect << endl;

  List<cmatrix> subsigma0s,subdetects;
  if (dosub) {
    split(subsigma0s,sigma0,block_str,1,verbose);
    split(subdetects,detect,block_str,-1,verbose);
  }

  const double intdt=period/isteps;

  const double rotor_phase=gamma_steps==1 ? PAS_to_RF.gamma : 0.0;

  RealSpinningHamiltonian H_cycle(rotor_speed,rotor_phase,MASRotorInfo,sampler);
  H_cycle.add(A_RF,H);
  if (spins>1)
    H_cycle+=HCS;
  if (spins>2)
    H_cycle.add(A_RF2,Hdip12);

  cout << H_cycle;

  const double rotor_period=1.0/fabs(rotor_speed);
  cout << "Hamiltonian at 0:\n" << H_cycle(0.0) << endl;
  cout << "Hamiltonian at t_r/2:\n" << H_cycle(rotor_period/2) << endl;
  cout << "Hamiltonian at t_r:\n" << H_cycle(rotor_period) << endl;

  HomogeneousPropagator<rmatrix> Ugen(H_cycle,intdt);
  Ugen.verbose(verbose_level);

  cmatrix U;
  Ugen(U,0.0,period/5);
  cout << "Tr/5:\n" << U << endl;
   
  Ugen(U,0.0,period/2);
  cout << "Tr/2:\n" << U << endl;

  Ugen(U,0.0,period);
  cout << "Tr:\n" << U << endl;

  Ugen(U,period,2*period);
  cout << "Tr-2Tr:\n" << U << endl;

  stopwatch.reset();
  int red_loop=rep_times/100;
  if (red_loop<1) red_loop=1;

  for (loop=red_loop;loop--;) {
    FID=complex(0.0);
    
    for (gn=0;gn<gamma_steps;gn++) {
      if (gamma_steps>1)
	H_cycle.rotor_phase(2*M_PI*gn/gamma_steps);

      cmatrix U;
      cmatrix Utot;

      FID(0)+=trace_multiply(sigma0,detect);
      for (j=1;j<total_steps;j++) {
	Ugen(U,(j-1)*stept,j*stept);
	Utot&=U;
	unitary_simtrans(sigma,sigma0,Utot);
	FID(j)+=trace_multiply(detect,sigma);
      }
    }
    if (gamma_steps>1) FID/=gamma_steps;
  }
  print_result("Very simple propagation: ",FID);
  print_usage(0,red_loop);

  if (gamma_steps==1) {
    AsynchronousFID obj(0.0,intdt);
    H_cycle.rotor_phase(rotor_phase);
    cmatrix Utmp;

//     List<cmatrix> Ualls(total_steps-1);
//     List< List<cmatrix> > subUalls(blks);
//     for (int k=blks;k--;) subUalls(k).create(total_steps-1);
//     for (j=0;j<total_steps-1;j++) {
//       Ugen(Utmp,j*stept,(j+1)*stept);
//       cmatrix& U(Ualls(j));
//       if (j)
// 	multiply(U,Utmp,Ualls(j-1));
//       else
// 	U=Utmp;
//       for (int blk=blks;blk--;)
// 	subUalls(blk)(j)=U(block_str(blk),block_str(blk));
//     }
//    if (verbose)
    //     cout << "Asynchronous propagators:\n" << Ualls << '\n';

    stopwatch.reset();
    for (loop=red_loop;loop--;) {
      if (isherm)
	FID=0.0;
      else
	FID=complex(0.0);
      
      obj.set_Ugenerator(Ugen);
      if (isherm)
	obj.add_FID_hermitian(FIDr,1.0,sigma0,detect);
      else
	obj.add_FID(FID,1.0,sigma0,detect);
    }
    const char* lstr="AsynchronousFID: ";
    if (isherm) 
      print_result(lstr,FIDr);
    else
      print_result(lstr,FID);
    print_usage(0,red_loop);
  }
    
// Simple propagation
  cmatrix *detect_tab=new cmatrix[n];

  cmatrix Ucycle;
  const int nprops=(gamma_steps>1) ? gamma_steps : n;
  List<cmatrix> shiftedUs(nprops);
  List< List<cmatrix> > subUs(blks);
  for (int blk=blks;blk--;) subUs(blk).create(nprops);

  stopwatch.reset();
  for (loop=red_loop;loop--;) {
    FID=complex(0.0);
    for (gn=0;gn<gamma_steps;gn++) {
      if (gamma_steps>1)
	H_cycle.rotor_phase(2*M_PI*gn/gamma_steps);

      cmatrix Utmp,U;

      for (j=0;j<n;j++) {
	if (j==0)
	  detect_tab[j]=detect;
	else
	  unitary_isimtrans(detect_tab[j],detect,U);
	Ugen(Utmp,j*stept,(j+1)*stept);
	U&=Utmp;
	if (gamma_steps==1) {
	  shiftedUs(j)=U;
	  for (int blk=blks;blk--;)
	    subUs(blk)(j)=U(block_str(blk),block_str(blk));
	}
      }
      Ucycle=U;

      k=0;
      sigma=sigma0;
      for (cyc=0;cyc<cycles;cyc++) {
	for (j=0;j<n;j++)
	  FID(k++)+=trace_multiply(detect_tab[j],sigma);
	sigma.unitary_simtrans(Ucycle);
      }
    }
    if (gamma_steps>1) FID/=gamma_steps;
  }
  print_result("Simple propagation (explicit): ",FID);
  print_usage(0,red_loop);
  if (verbose) 
    cout << "shiftedUs:\n" << shiftedUs;

#ifndef NDEBUG
  cout.precision(3);
#endif

  if (gamma_steps==1) {
    PeriodicFID obj;
    const int blks=subUs.length();

    for (int ldosub=0;ldosub<=maxdosub;ldosub++) {

      stopwatch.reset();

      for (loop=red_loop;loop--;) {
	if (isherm || isdiag)
	  FIDr=0.0;
	else
	  FID=complex(0.0);

	if (ldosub) {
	  obj.set_Us('C',subUs(0));

	  for (size_t blk=1;blk<blks;blk++) {
	    obj.shuffle('C');
	    obj.set_Us('C',subUs(blk));
	    obj.add_FID(FID,1.0,subsigma0s(blk-1),subdetects(blk-1));
	  }
	}
	else {
	  obj.set_Us(shiftedUs);
	  
	  if (isdiag)
	    obj.add_FID(FIDr,1.0,dsigma0,ddetect);
	  else {
	    if (isherm)
	      obj.add_FID_hermitian(FIDr,1.0,sigma0,detect);
	    else
	      obj.add_FID(FID,1.0,sigma0,detect);
	  }
	}
      }
      const char* lstr=ldosub ? "PeriodicFID (sub): " : "PeriodicFID: ";
      if (isherm || isdiag) 
	print_result(lstr,FIDr);
      else
	print_result(lstr,FID);
      print_usage(0,red_loop);
    }
  }

  const rmatrix rsigma0(real(sigma0));
  const rmatrix rdetect(real(detect));
  
//   if (floquet_N) {
//     rmatrix FHam=real(floquet_hamiltonian(H_cycle,floquet_N,rotor_speed));

//     cout << "Floquet Hamiltonian:\n"; 
//   //set_frequency(FHam,floquet_N,rotor_speed);
//   //cout << "Floquet Hamiltonian (set_frequency):\n"; 

// #ifdef NDEBUG
//     spy(FHam);
// #else
//     cout << FHam << '\n';
// #endif

//     FloquetFID ffid;
//     GammaFloquetFID gffid;
    
//     stopwatch.reset();
//     for (loop=red_loop;loop--;) {
//       if (gamma_steps==1) {
// 	ffid.set_H(FHam,floquet_N,1.0/sw);
// 	if (samedetinit) 
// 	  FID=ffid.FID(total_steps,rsigma0,rdetect);
// 	else
// 	  FID=ffid.FID(total_steps,sigma0,detect);
//       }
//       else {
// 	gffid.set_H(FHam,floquet_N,1.0/sw);
// 	if (samedetinit) 
// 	  FIDr=gffid.FID_hermitian(total_steps,rsigma0,rdetect);
// 	else
// 	  FID=gffid.FID(total_steps,sigma0,detect);
//       }      
//     }

//     if (samedetinit && (gamma_steps>1)) 
//       print_result("Floquet propagation (TD): ",FIDr);
//     else
//       print_result("Floquet propagation (TD): ",FID);
    
//     print_usage(0,red_loop);
    
//     FloquetSpectrum fobj;
//     GammaFloquetSpectrum gfobj;
    
//     stopwatch.reset();
//     for (loop=red_loop;loop--;) {
//       if (gamma_steps==1)
// 	dofloquetspectrum(FID,fobj,FHam,floquet_N,rsigma0,rdetect,verbose);
//       else {
// 	if (samedetinit)
// 	  dofloquetspectrum(FID,gfobj,FHam,floquet_N,rsigma0,rdetect,verbose);
// 	else
// 	  dofloquetspectrum(FID,gfobj,FHam,floquet_N,sigma0,detect,verbose);
//       }
//     }
//     print_result("Floquet propagation (FD): ",FID);
//     print_usage(0,red_loop);
//   }
    
//   else { //gamma_steps!=1
//     cout << "Simple propagation (add_GammaPeriodicFID): ";

//     stopwatch.reset();
//     List<complex> tFID(total_steps);
//     for (loop=red_loop;loop--;) {
//       tFID=complex(0.0);  
//       add_GammaPeriodicFID(tFID,1.0,n,shiftedUs,sigma0,detect);
//     }
//     cout << tFID << endl;
//     print_usage(0,red_loop);
//   }
  
  double normdiff=-1;
  bool isdiagonal;
  List<complex> cycle_eigs;

  {
    cmatrix R,UT,S,D,sigma0t;
    cmatrix *Ustore=new cmatrix[n];

    stopwatch.reset();
    for (loop=red_loop;loop--;) {
      FID=complex(0.0);
      for (gn=0;gn<gamma_steps;gn++) {
	if (gamma_steps>1)
	  H_cycle.rotor_phase(2*M_PI*gn/gamma_steps);
	
	cmatrix U,Utmp;

	for (j=0;j<n;j++) {
	  if (j) Ustore[j]=U;
	  Ugen(Utmp,j*stept,(j+1)*stept);
	  U&=Utmp;
	}
	Ucycle=U;
	
	if (normdiff<0) {
	  normdiff=sqrt(norm(sum(Ucycle)-trace(Ucycle)));
	  cout << "Norm diff: " << normdiff << endl;
	  isdiagonal=(normdiff<1e-8);
	  cout << "Diagonal: " << (isdiagonal ? "Yes" : "No") << endl;
	}
	
	if (isdiagonal) {
	  cycle_eigs=diag(Ucycle);
	  D.identity(dim);
	}
	else
	  eigensystem(D,cycle_eigs,Ucycle);	
	unitary_isimtrans(sigma0t,sigma0,D);

	evolution_matrix(R,cycle_eigs);
	
	for (j=0;j<n;j++) {
	  if (isdefined(Ustore[j])) {
	    multiply(UT,Ustore[j],D);
	    unitary_isimtrans(S,detect,UT);
	  }
	  else
	    unitary_isimtrans(S,detect,D);
	  S.transpose();
	  S.emultiply(sigma0t);
	  
	  k=j;
	  for (cyc=cycles;cyc--;) {
	    FID(k)+=sum(S);
	    k+=n;
	    S.emultiply(R);
	  }
	}
      }
      if (gamma_steps>1) FID/=gamma_steps;
    }
    print_result("Simple in-basis propagation: ",FID);
    print_usage(0,red_loop);
  }

  const int intsteps=gamma_steps/n;
  
  const double microst=period/gamma_steps;

  cmatrix *Ustore=NULL;

  double prop_calc=0;
  if (gamma_steps>1) {
    Ustore=new cmatrix[gamma_steps];

    H_cycle.rotor_phase(0.0);
    stopwatch.reset();
    cmatrix Utmp;

    for (loop=rep_times;loop--;) {
      Ucycle.kill();
      for (i=0;i<gamma_steps;i++) {
	if (i==0)
	  Ustore[i].identity(dim);
	else
	  Ustore[i]=Ucycle;
	Ugen(Utmp,i*microst,(i+1)*microst);
	Ucycle&=Utmp;
	shiftedUs(i)=Ucycle;
	for (int blk=blks;blk--;)
	  subUs(blk)(i)=Ucycle(block_str(blk),block_str(blk));
      }
    }
    prop_calc=stopwatch();
    cout << "Calculation of propagators: "; print_usage(0,rep_times);
  }
  else {
    Ucycle.kill();
    Ugen(Ucycle,0,period);
  }

  cout << "Cycle propagator:\n" << Ucycle << endl;

  if (gamma_steps>1) {
    cmatrix T,detect_cycle,sigmanext;
    cmatrix *detect_table = new cmatrix[gamma_steps];
    
    stopwatch.reset();

    for (loop=rep_times;loop--;) {
      for (gn=0;gn<gamma_steps;gn++) {
	if (gn==0)
	  detect_table[gn]=detect;
	else
	  unitary_isimtrans(detect_table[gn],detect,Ustore[gn]);
      }
    }
    const double table_calc=stopwatch();
    cout << "Table: "; print_usage(0,rep_times);

    // Over 50% of the time is spent in trace in this algorithm!
    stopwatch.reset();
    for (loop=rep_times;loop--;) {
      FID=complex(0,0);
    
      for (gn=0;gn<gamma_steps;gn++) {
	unitary_isimtrans(sigma,sigma0,Ustore[gn]);

	k=0;
	for (cyc=0;cyc<cycles;cyc++) {
	  unitary_simtrans(sigmanext,sigma,Ucycle);
	  for (j=0;j<n;j++) {
	    const int which=j*intsteps+gn;
	    if (which>=gamma_steps)
	      FID(k++)+=trace_multiply(detect_table[which-gamma_steps],sigmanext);
	    else
	      FID(k++)+=trace_multiply(detect_table[which],sigma);
	  }
	  sigma=sigmanext;
	}
      }
      FID/=gamma_steps;
    }
    print_result("Stored propagators: ",FID);
    print_usage(prop_calc+table_calc,rep_times);
  }

  // Homogenous in basis set

  cmatrix *sigma_table=NULL;
  cmatrix *detect_table=NULL;

  normdiff=sqrt(norm(sum(Ucycle)-trace(Ucycle)));
  isdiagonal=(normdiff<1e-8);
  cout << "Diagonal: " << (isdiagonal ? "Yes" : "No") << endl;

  double transfops=-1e30;

  cmatrix R,D;

  if (gamma_steps>1) {
    cmatrix tmp,UT;
  
    sigma_table=new cmatrix[gamma_steps];
    detect_table=new cmatrix[2*gamma_steps];

    stopwatch.reset();
    for (loop=rep_times;loop--;) {
      if (isdiagonal)
	cycle_eigs=diag(Ucycle);
      else
	eigensystem(D,cycle_eigs,Ucycle);

      evolution_matrix(R,cycle_eigs);
      
      for (gn=0;gn<gamma_steps;gn++) {
	if (isdefined(D))
	  multiply(UT,Ustore[gn],D);
	else
	  UT=Ustore[gn];
	unitary_isimtrans(detect_table[gn],detect,UT);
	unitary_isimtrans(sigma_table[gn],sigma0,UT);
	//cout << "S" << gn << "\n" << sigma_table[gn] << endl;
	// only needed by some of the methods
	detect_table[gn].transpose();
	//cout << "det" << gn << "\n" << detect_table[gn] << endl;
	emultiply(detect_table[gn+gamma_steps],detect_table[gn],R);
      }
    }
    transfops=stopwatch();
  }
  else {
    if (isdiagonal)
      cycle_eigs=diag(Ucycle);
    else
      eigensystem(D,cycle_eigs,Ucycle);
  }

  cout << "Eigenvalues: " << cycle_eigs << endl;
  cout << "R:\n" << R << endl;
  if (isdefined(D)) {
    cmatrix partial;
    multiply(partial,D,cycle_eigs);
    cout << "Recon:\n" << partial*conj_transpose(D) << endl;
  }
  cout << "Transformed operators: "; print_usage(0,rep_times);
  cout << "Eigenvalues: " << cycle_eigs << endl;
  if (isdefined(D))
    cout << "Eigenvector orthogonality test:\n" << conj_transpose(D)*D << endl;
  
  // in basis propagation
  if (gamma_steps>1) {
    cmatrix *Sstore=new cmatrix[n];
    cmatrix tmp;

    stopwatch.reset();
    for (loop=rep_times;loop--;) {
      for (j=0;j<n;j++) {
	cmatrix &S=Sstore[j];
	cmatrix S2;
	const int offset=j*intsteps;
	emultiply(S,detect_table[offset],sigma_table[0]);
	for (gn=1;gn<gamma_steps;gn++) {
	  const int which=offset+gn;
	  if (which>=gamma_steps) {
	    if (!isdefined(S2))
	      emultiply(S2,detect_table[which-gamma_steps],sigma_table[gn]);
	    else
	      mla(S2,detect_table[which-gamma_steps],sigma_table[gn] );
	  }
	  else
	    mla(S,detect_table[which],sigma_table[gn]);
	}
	if (isdefined(S2)) mla(S,R,S2);
      }
    }    
    const double prep_time=stopwatch();
    cout << "In basis prep: ";
    print_usage(0,rep_times);

    stopwatch.reset();
    for (loop=rep_times;loop--;) {
      FID=complex(0);
      for (j=0;j<n;j++) {
	const cmatrix &S=Sstore[j];
	for (r=dim;r--;) {
	  for (s=dim;s--;) {
	    complex a=S(r,s);
	    if (norm(a)>1e-12) {
	      const complex &p=R(r,s);
	      k=j;
	      for (cyc=cycles;cyc--;) {
		FID(k)+=a;
		a*=p;
		k+=n;
	      }
	    }
	  }
	}
      }
      FID/=double(gamma_steps);
    }
    print_result("Gamma Homogenous in basis: ",FID);
    print_usage(transfops+prop_calc+prep_time,rep_times);

    GammaPeriodicFID obj(n);

    for (int ldosub=0;ldosub<=maxdosub;ldosub++) {

      stopwatch.reset();
      
      for (loop=rep_times;loop--;) {
	if (isdiag || isherm)
	  FIDr=0.0;
	else
	  FID=complex(0.0);

	if (ldosub) {
	  obj.set_Us('C',subUs(0));

	  for (int blk=1;blk<blks;blk++) {
	    obj.shuffle('C');
	    obj.set_Us('C',subUs(blk));
	    if (samedetinit)
	      obj.add_FID(FID,1.0,subsigma0s(blk-1));
	    else
	      obj.add_FID(FID,1.0,subsigma0s(blk-1),subdetects(blk-1));
	  }
	}
	else {
	  obj.set_Us(shiftedUs);
	  
	  if (samedetinit) {
	    if (isdiag)
	      obj.add_FID(FIDr,1.0,dsigma0);
	    else {
	      if (isherm)
		obj.add_FID_hermitian(FIDr,1.0,sigma0);
	      else
		obj.add_FID(FID,1.0,sigma0);
	    }
	  }
	  else {
	    if (isdiag)
	      obj.add_FID(FIDr,1.0,dsigma0,ddetect);
	    else {
	      if (isherm)
		obj.add_FID_hermitian(FIDr,1.0,sigma0,detect);
	      else
		obj.add_FID(FID,1.0,sigma0,detect);
	    }
	  }
	}
      }
      const char* lstr=ldosub ? "GammaPeriodicFID (sub): " : "GammaPeriodicFID: ";
      if (isdiag || isherm)
	print_result(lstr,FIDr);
      else
	print_result(lstr,FID);
      print_usage(prop_calc,rep_times);    
    }
  }
  else { // single orientation
    List<cmatrix> Sstore(n);
    List<cmatrix> lSstore(n);
    cmatrix tmp,sigma0t;

    for (loop=rep_times;loop--;) {
    
      evolution_matrix(R,cycle_eigs);
      R.transpose();

      if (!isdiagonal) {
	unitary_isimtrans(sigma0t,sigma0,D);
	sigma0t.transpose();
      }
      else
	transpose(sigma0t,sigma0);

      cmatrix U,Utmp;

      for (j=0;j<n;j++) {
	if (isdefined(U))  {
	  if (isdefined(D)) {
	    multiply(tmp,U,D);
	    unitary_isimtrans(Sstore(j),detect,tmp);
	  }
	  else
	    unitary_isimtrans(Sstore(j),detect,U);
	}
	else {
	  if (isdefined(D))
	    unitary_isimtrans(Sstore(j),detect,D);
	  else
	    Sstore(j)=detect;
	}
	Sstore(j).emultiply(sigma0t);
	Ugen(Utmp,j*stept,(j+1)*stept);
	U&=Utmp;
      }
    }
    const double overhead=stopwatch();

    stopwatch.reset();
    for (loop=rep_times;loop--;) {
      k=0;
      lSstore=Sstore;
      for (cyc=cycles;cyc--;) {
	for (j=0;j<n;j++) {
	  FID(k++)=sum(lSstore(j));
	  lSstore(j).emultiply(R);
	}
      }
    }
    
    print_result("Homogeneous (SO) in basis: ",FID);
    print_usage(overhead,rep_times);
  }
  
  if (gamma_steps>1) {
    const int usefft=ispowerof2(gamma_steps);
    cout << "Fast FFT: " << (usefft ? "Yes" : "No") << endl;

    List<complex> ft_facs(gamma_steps);
    ft_create_table(ft_facs,FT_FORWARD);
    
    List<double> eff_freq(dim);

    for (j=0;j<dim;j++)
      eff_freq(j)=-rotor_speed*arg(cycle_eigs(j))/(2*M_PI);
    cout << "Effective frequencies: " << eff_freq << endl;

    // Homogenous in basis FT
    
    List<complex> sideband(gamma_steps);
    List<complex> sigma0ft(gamma_steps);
    List<complex> detectft(gamma_steps);
    List<complex> scr(gamma_steps);
    
    List<complex> & store=usefft ? sideband : scr;

    complex** sigma_addr=new complex* [gamma_steps];
    complex** detect_addr=new complex* [gamma_steps];

    stopwatch.reset();
    for (loop=rep_times;loop--;) {
      Output=complex(0,0);
            
      for (r=0;r<dim;r++) {
	for (gn=gamma_steps;gn--;) {
	  sigma_addr[gn]=sigma_table[gn].vector(r);
	  detect_addr[gn]=detect_table[gn].vector(r);
	}

	for (s=0;s<dim;s++) {
	  const double diff=eff_freq(r)-eff_freq(s);
	  
	  scr=complex(0,0);

	  complex p=1.0;
	  const complex pfac=expi(2*M_PI*diff*microst);

	  for (j=0;j<gamma_steps;j++) {	 
	    complex sum(0,0);

	    for (int which=j;which<gamma_steps;which++)
	      mla(sum,sigma_addr[which-j][s],detect_addr[which][s]);
	    if (j) {
	      const int offset=gamma_steps-j;
	      complex s2(0,0);
	      for (int which2=j;which2--;) 
		mla(s2,sigma_addr[which2+offset][s],detect_addr[which2][s]);
	      mla(sum,R(r,s),s2);
	    }
	    store(j)=sum*p; 
	    p*=pfac;
	  }
	  if (verbose) {
	    cout << "Evolution phase factor for " << r << ", " << s << " transition: " << pfac << '\n';
	    cout << "Raw factors: " << store << '\n';
	  }

	  if (usefft)
	    fft_ip(sideband,FT_FORWARD);
	  else
	    ft(sideband,scr,ft_facs);

	  for (j=n;j--;) {
	    complex totamp=0;
	    for (k=j;k<gamma_steps;k+=n) totamp+=sideband(k);
	    if (norm(totamp)>1e-12) {
	      double freq=-diff+j*rotor_speed;
	      if (verbose) cout << r << ", " << s << ", " << j << ": " << freq << "  " << totamp << endl;
	      addfunc(Output,totamp,freq/sw);
	    }
	  }
	}
      }
      Output/=double(gamma_steps*gamma_steps);
    }
    print_result("Homogeneous sum FT: ",Output,asfreq);
    print_usage(transfops+prop_calc,rep_times);

    List<complex> scr2(gamma_steps);

    List<complex> & store1=usefft ? sigma0ft : scr;
    List<complex> & store2=usefft ? detectft : scr2;
    
    stopwatch.reset();
    const double sfac=1.0/(gamma_steps*gamma_steps);

    if (verbose) {
      cout << "sigma_table\n";
      dump_table(BaseList<cmatrix>(gamma_steps,sigma_table));
      cout << "detect_table\n";
      dump_table(BaseList<cmatrix>(gamma_steps,detect_table));
    }

    for (loop=rep_times;loop--;) {
      Output=complex(0.0);

      for (r=0;r<dim;r++) {
	for (s=0;s<dim;s++) {
	  const double diff=eff_freq(r)-eff_freq(s);

	  const complex pfac=expi(diff*2*M_PI*microst);
	  complex phas(1.0);
	  
	  for (j=0;j<gamma_steps;j++) {
	    store1(j)=conj_multiply(sigma_table[j](r,s),phas);
	    store2(j)=detect_table[j](r,s)*phas;
	    phas*=pfac;
	  }

	  if (verbose) {
	    cout << "Evolution phase factor for " << r << ", " << s << " transition: " << pfac << '\n';
	    cout << "Raw factors (sigma0): " << store1 << '\n';
	    cout << "Raw factors (detect): " << store2 << '\n';
	  }

	  // Note that many of these vectors are in fact zero; it would be more efficient to skip over such r,s pairs
	  // Nearly half the time is spent in the FTs
	  if (usefft) {
	    fft_ip(sigma0ft,FT_FORWARD);
	    fft_ip(detectft,FT_FORWARD);
	  }
	  else {
	    ft(sigma0ft,store1,ft_facs);
	    ft(detectft,store2,ft_facs);
	  }
	  //cout << r << ", " << s << ": " << diff << "  " << sigma0ft << endl;
	  
	  for (j=0;j<n;j++) {
	    complex amp(0,0);
	    for (k=j;k<gamma_steps;k+=n)
	      mla_conj(amp,detectft(k),sigma0ft(k));
	    const double freq=-diff+j*rotor_speed;
	    if (verbose) cout << r << ", " << s << ", " << j << ": " << freq << "  " << amp*sfac << endl;
	    if (norm(amp)>1e-12)
	      addfunc(Output,amp,freq/sw);
	  }
	}
      }
      Output*=sfac;
    }
    print_result("Homogeneous FT: ",Output,asfreq);
    print_usage(transfops+prop_calc,rep_times);

    GammaPeriodicSpectrumED pobjED(n);
    GammaPeriodicSpectrum pobj(n);

    if (dosub) {
      stopwatch.reset();

      for (loop=rep_times;loop--;) {
	Output=complex(0.0);
	
	if (samedetinit)
	  docalc(Output,pobjED,subUs,period,subsigma0s,verbose);
	else
	  docalc(Output,pobj,subUs,period,subsigma0s,subdetects,verbose);
      }
      print_result("GammaPeriodicSpectrum (sub): ",Output,asfreq);
      print_usage(prop_calc,rep_times);
    }

    stopwatch.reset();

    for (loop=rep_times;loop--;) {
      Output=complex(0.0);

      if (samedetinit)
	docalc(Output,pobjED,shiftedUs,period,sigma0,verbose);
      else
	docalc(Output,pobj,shiftedUs,period,sigma0,detect,verbose);
    }
  }
  else {
    PeriodicSpectrum pobj;
    
    if (dosub) {
      stopwatch.reset();

      for (loop=rep_times;loop--;) {
	Output=complex(0.0);
	docalc(Output,pobj,subUs,period,subsigma0s,subdetects,verbose);
      }
      print_result("PeriodicSpectrum (sub): ",Output,asfreq);
      print_usage(prop_calc,rep_times);
    }

    stopwatch.reset();

    for (loop=rep_times;loop--;) {
      Output=complex(0.0);      
      if (samedetinit)
	docalc(Output,pobj,shiftedUs,period,sigma0,verbose);
      else
	docalc(Output,pobj,shiftedUs,period,sigma0,detect,verbose);
    }
  }
  print_result("Homogenous FT 2/(Gamma)PeriodicSpectrum: ",Output,asfreq);
  print_usage(prop_calc,rep_times);

  } catch (MatrixException& exc) {
    cerr << exc << endl;
    return 1;
  }

  cout << "Failures: " << errcount << '\n';
  return errcount;
}
