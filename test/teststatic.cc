#include "spinhalf_system.h"
#include "Propagation.h"
#include "NMR.h"
#include "ListList.h"
#include "ttyio.h"
#include "timer.h"

using namespace std;
using namespace libcmatrix;

//teststatic 2 y n y 1
timer<> stopwatch;

int errcount=0;
List<double> refreal;
List<complex> refcomplex;

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

double SW=0.0;

inline double norm(double x) { return x*x; }

template<class F> void docalc(BaseList<complex> FID,F obj,bool verbose)
{
  typename F::amplitude_type amp;
  double freq;
  if (verbose) cout << obj << endl;
  while (obj(amp,freq)) {
    if (norm(amp)>1e-12) {
      if (verbose) cout << amp << " of " << freq << endl;
      add_FID(FID,amp,freq/SW);
    }
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

void print_usage(int rep_times)
{
  double ttake=stopwatch();
  cout << "Time taken: " << ttake << "    " << "Per calculation: "<< (ttake*1e3/rep_times) << " ms\n";
}

template<class T> void split(List< Matrix<T> >& dest,const Matrix<T>& source,const ListList<size_t>& block_str,int offset,int verbose =0)
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

int main(int argc,const char *argv[])
{
  try {
  int count=1;
  const int SPINS=getint(argc,argv,count,"Number of spins? ",2,3);
  const bool useU=getlogical(argc,argv,count,"Use U? ");
  const bool useED=getlogical(argc,argv,count,"Use ED forms? ");
  const bool isherm=getlogical(argc,argv,count,"Hermitian detection? ");
  const int reptimes=getint(argc,argv,count,"Repeat times? ",1);
  const bool verbose=(reptimes==1);
  const int maxdosub=isherm ? 0 : 1;

  spinhalf_system ax(SPINS);

  const List<double> Fzop=diag_Fz(ax);
  const ListList<size_t> str=find_blocks(Fzop,0.0);

  const double dips[3][3]={1000,2000,3000,2000,4000,5000,3000,4000,1500};

  cmatrix Hdipc;

  for(int i=0;i<SPINS;i++)
    for (int j=0;j<i;j++)
      mla(Hdipc,dips[i][j],spin_dipolar(ax,i,j));
  rmatrix Hdip=real(Hdipc);

  cout << "H\n" << Hdip << endl;

  int loop;

  const cmatrix sigma0=F(ax,isherm ? 'x' : '+');
  const cmatrix detect=conj_transpose(sigma0);
  
  cmatrix sigma0T,detectT;

  rmatrix D;
  List<double> eigs;

  SW=10000;
  const double dt=1/SW;
  const int npts=10;
  
  List<complex> FID(npts);
  List<double> FIDr(npts);

  cmatrix A;

  stopwatch.reset();
  for (loop=reptimes;loop--;) {
    FID=0.0;
    hermitian_eigensystem(D,eigs,Hdip);

    if (reptimes==1) cout << "D\n" << D << endl;

    unitary_isimtrans(sigma0T,sigma0,D);
    if (reptimes==1) cout << "Transformed sigma0\n" << sigma0T << endl;
       
    if (!useED) {
      unitary_isimtrans(detectT,detect,D);
      if (reptimes==1) cout << "Transformed detect\n" << detectT << endl;
      emultiply(A,transpose(detectT),sigma0T);
    }
    else
      A=enorm(sigma0T);
    if (reptimes==1) cout << "Transition amplitudes:\n" << A << endl;
    
    const int dim=eigs.length();
    for (int r=0;r<dim;r++) {
      for (int c=0;c<dim;c++) {
	const complex amp=A(r,c);
	const double freq=eigs(r)-eigs(c);
	if (norm(amp)>1e-12) {
	  // cout << amp << " of " << freq << endl;
	  add_FID(FID,amp,freq/SW);
	}
      }
    }
  }
  print_usage(reptimes);
  print_result("Explicit propagation: ",FID);

  const cmatrix U(propagator(Hdip,dt));

  //create blocks
  const int nblks=str.size();
  List<cmatrix> UBs(nblks);
  List<rmatrix> HBs(nblks);
  List<cmatrix> sigma0s(nblks-1);
  List<cmatrix> detects(nblks-1);

  if (maxdosub) {
    split(HBs,Hdip,str,0,verbose);
    split(UBs,U,str,0,verbose);
    split(sigma0s,sigma0,str,1,verbose);
    split(detects,detect,str,-1,verbose);
  }

  StaticFID_H fgen_H;
  StaticFID_U fgen_U;

  int dosub;
  for (dosub=0;dosub<=maxdosub;dosub++) {
    stopwatch.reset();
    for (loop=reptimes;loop--;) {
      if (isherm)
	FIDr=0.0;
      else
	FID=complex(0.0);
      
      if (dosub) {
	for (int n=0;n<str.size();n++) {

	  if (useU)
	    fgen_U.set_U('C',UBs(n));
	  else
	    fgen_H.set_H('C',HBs(n),dt);
	  
	  if (n) {
	    const cmatrix& lsigma0=sigma0s(n-1);
	    const cmatrix& ldetect=detects(n-1);
	    if (useU)
	      fgen_U.add_FID(FID,1.0,lsigma0,ldetect);
	    else {
	      if (useED)
		fgen_H.add_FID(FID,1.0,lsigma0);
	      else
		fgen_H.add_FID(FID,1.0,lsigma0,ldetect);	
	    }
	  }
	  if (useU)
	    fgen_U.shuffle('C');
	  else
	    fgen_H.shuffle('C');
	}
      }
      else {
	if (useU) {
	  fgen_U.set_U(U);
	  if (isherm)
	    fgen_U.add_FID_hermitian(FIDr,1.0,sigma0,detect);
	  else
	    fgen_U.add_FID(FID,1.0,sigma0,detect);
	}
	else {
	  fgen_H.set_H(Hdip,dt);
	  
	  if (isherm) {
	    if (useED)
	      fgen_H.add_FID_hermitian(FIDr,1.0,sigma0);
	    else
	      fgen_H.add_FID_hermitian(FIDr,1.0,sigma0,detect);	
	  }
	  else {
	    if (useED)
	      fgen_H.add_FID(FID,1.0,sigma0);
	    else
	      fgen_H.add_FID(FID,1.0,sigma0,detect);	
	  }
	}
      }
    }
    print_usage(reptimes);
    const char * name = dosub ? "StaticFID (sub): " : "StaticFID: ";
    if (isherm)
      print_result(name,FIDr);
    else
      print_result(name,FID);
  }

  StaticSpectrum sgen2;
  StaticSpectrumED sgen2ED;

  for (dosub=0;dosub<=maxdosub;dosub++) {
    stopwatch.reset();
    for (loop=reptimes;loop--;) {
      FID=0.0;
      
      if (dosub) {
	for (int n=0;n<str.size();n++) {
	  
	  if (useED) {
	    if (!useU)
	      sgen2ED.set_H('C',HBs(n));
	    else
	      sgen2ED.set_U('C',UBs(n),dt);
	  }
	  else {
	    if (!useU)
	      sgen2.set_H('C',HBs(n));
	    else
	      sgen2.set_U('C',UBs(n),dt);
	  }
	  
	  if (n) {
	    if (useED) {
	      sgen2ED.observe(sigma0s(n-1));
	      docalc(FID,sgen2ED,verbose);
	    }
	    else {
	      sgen2.observe(sigma0s(n-1),detects(n-1));
	      docalc(FID,sgen2,verbose);
	    }
	  }
	  if (useED)
	    sgen2ED.shuffle('C');
	  else
	    sgen2.shuffle('C');
	}
      }
      else {
	if (useED) {
	  if (!useU)
	    sgen2ED.set_H(Hdip);
	  else
	    sgen2ED.set_U(U,dt);
	  sgen2ED.observe(sigma0);
	  docalc(FID,sgen2ED,verbose);	  
	}
	else {
	  if (!useU)
	    sgen2.set_H(Hdip);
	  else
	    sgen2.set_U(U,dt);
	  sgen2.observe(sigma0,detect);
	  docalc(FID,sgen2,verbose);
	}
      }
    }  
    print_usage(reptimes);
    print_result(dosub ? "StaticSpectrum(ED) (sub): " : "StaticSpectrum(ED): ",FID);
  }
  }
  catch (MatrixException &exc) {
    cerr << exc << endl;
    return 1;
  }

  cout << "Failures: " << errcount << '\n';
  return errcount;
}
