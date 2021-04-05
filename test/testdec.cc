// test homogeneous propagation for heteronuclear decoupling

#include "spinhalf_system.h"
#include "Propagation.h"
#include "FunctionObject.h"
#include "NMR.h"
#include "ttyio.h"
#include "timer.h"
#include "rmatrix.h"

timer<> stopwatch;

void print_usage(int rep_times)
{
  double ttake=stopwatch();
  cout << "Time taken: " << ttake << "    " << "Per calculation: "<< (ttake*1e3/rep_times) << " ms\n";
}

int main(int argc,const char *argv[])
{
  try {
  int count=1;
  const int SPINS=getint(argc,argv,count,"Number of spins? ",1,2);
  const bool useED=getlogical(argc,argv,count,"Use ED forms? ");

  int form;
  enum { ASC, ASR, ASID };
  if (useED) {
    cout << "C - complex full matrix\nR - real full matrix\nI - identity\n";
    form=getoption(argc,argv,count,"Form of sigma0/detection operators? ","CRI");
  }
  else
    form=getlogical(argc,argv,count,"Complex matrix? ") ? ASC : ASR;

  const int reptimes=getint(argc,argv,count,"Repeat times? ",1);

  spinhalf_system ax(SPINS);

  List<complex> DUMBO1_coeffs(7);
  DUMBO1_coeffs(0)=0.1056;
  DUMBO1_coeffs(1)=complex(0.0325,0.1310);
  DUMBO1_coeffs(2)=complex(0.0189,0.1947);
  DUMBO1_coeffs(3)=complex(0.0238,0.0194);
  DUMBO1_coeffs(4)=complex(0.0107,0.1124);
  DUMBO1_coeffs(5)=complex(0.0038,-0.0456);
  DUMBO1_coeffs(6)=complex(-0.0013,0.0869);

  const double dips[3][3]={{0,2000,3000},{3000,0,5000},{3000,5000,0}};
  const double wRF=getfloat(argc,argv,count,"RF strength? ",10000);

  enum { CW, LG, FSLG, DUMBO };
  int dectype=CW;

  if (wRF) {
    cout << "C - CW\nL - LG\nF - FSLG\nD - DUMBO1\n";
    dectype=getoption(argc,argv,count,"Decouple method? ","CLFD");
  }

  cmatrix Hbasec;
  List<double> Hhet;

  for(int i=0;i<SPINS;i++) {
    mla(Hhet,dips[0][i+1],diag_Iz(ax,i));
    for (int j=0;j<i;j++)
      mla(Hbasec,dips[i+1][j+1],spin_dipolar(ax,i,j));
  }
  const List<double> Fz=diag_Fz(ax);
  const int dim=Hhet.length();
  
  rmatrix Hrf=wRF*real(F(ax,'x'));

  rmatrix Halpha,Hbeta;
  if (isdefined(Hbasec)) {
    const rmatrix Hbase=real(Hbasec);
    Halpha=Hbase+Hhet;
    Hbeta=Hbase-Hhet;
  }
  else {
    Halpha=full(Hhet);
    Hbeta=full(-Hhet);
  }

  cmatrix Ualpha,Ubeta;
  double cycle_period;
  
  switch (dectype) {
  case CW:
    Halpha+=Hrf;
    Hbeta+=Hrf;
    break;
  case DUMBO:
    {
      const double mod_freq=wRF/3;
      cycle_period=1.0/mod_freq;
      SymmetricModFunc DUMBO(DUMBO1_coeffs,mod_freq);
      const double dt=cycle_period/64;
      Halpha+=Hrf;
      Hbeta+=Hrf;
      const cmatrix Ualpha0=propagator(Halpha,dt);
      const cmatrix Ubeta0=propagator(Hbeta,dt);

      cout << "DUMBO phases: ";
      for (int j=0;j<64;j++) {
	const double phase=DUMBO(j/(64*mod_freq));
	cout << phase << ' ';
	Ualpha&=rotatez(Ualpha0,Fz,phase);
	Ubeta&=rotatez(Ubeta0,Fz,phase);
      }
      cout << endl;
    } 
    break;
  case FSLG: 
    {
      mla(Hrf,wRF/sqrt(2.0),Fz);
      const double cycletime=1.0/(sqrt(3.0/2)*wRF);
      Ualpha=propagator(Halpha+Hrf,cycletime);
      Ualpha&=propagator(Halpha-Hrf,cycletime);
      Ubeta=propagator(Hbeta+Hrf,cycletime);
      Ubeta&=propagator(Hbeta-Hrf,cycletime);
      cycle_period=2*cycletime;
    }
    break;
  case LG:
    Hrf+=wRF/sqrt(2.0)*Fz;
    Halpha+=Hrf;
    Hbeta+=Hrf;
    break;
  }

  if (isdefined(Ualpha)) {
    cout << "Ualpha\n" << Ualpha << std::endl;
    cout << "Ubeta\n" << Ubeta << std::endl;
  }
  else {
    cout << "Halpha\n" << Halpha << std::endl;
    cout << "Hbeta\n" << Hbeta << std::endl;
  }

  const rmatrix rident=identity(dim);
  const cmatrix cident=rident;
  
  //  List<double> FIDref(npts);

  stopwatch.reset();
  for (int loop=reptimes;loop--;) {

    if (useED) {
      SimpleSpectrumED specobj;
      if (isdefined(Ualpha)) {
	specobj.set_U('r',Ualpha,cycle_period);
	specobj.set_U('c',Ubeta,cycle_period);
      }
      else {
	specobj.set_H('r',Halpha);
	specobj.set_H('c',Hbeta);
      }
            
      switch (form) {
      case ASC: specobj.observe(cident); break;
      case ASR: specobj.observe(rident); break;
      case ASID: specobj.observe(); break;
      }
     
      if (reptimes==1) cout << specobj << endl;

      double amp,freq;
      while (specobj(amp,freq)) {
	if (reptimes==1)
	  cout << amp << " of " << freq << endl;
      }
    }
    else {
      SimpleSpectrum specobj;
      if (isdefined(Ualpha)) {
	specobj.set_U('r',Ualpha,cycle_period);
	specobj.set_U('c',Ubeta,cycle_period);
      }
      else {
	specobj.set_H('r',Halpha);
	specobj.set_H('c',Hbeta);
      }
      
      if (form==ASC) 
	specobj.observe(cident,cident);
      else
	specobj.observe(rident,rident);
      
      if (reptimes==1) cout << specobj << endl;

      double freq;
      complex amp;
      while (specobj(amp,freq)) {
	if (reptimes==1) cout << amp << " of " << freq << endl;
      }
    }
  }
  if (reptimes!=1) print_usage(reptimes);
  }
  catch (MatrixException &exc) {
    cerr << exc << endl;
  }
  return 0;
}
