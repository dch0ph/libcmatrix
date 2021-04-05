#include "spinhalf_system.h"
#include "NMR.h"
#include "tensorop.h"
#include "cmatrix.h"
#include "Histogram.h"
#include "simpsonio.h"
#include "ttyio.h"
#include "powder.h"
#include "MAS.h"
#include "Propagation.h"

using namespace libcmatrix;
using namespace std;

//**************************************************************************Change left and right part of spectrum*************************************************************************************
void halfshift(List<double> & spectrum) {

int L=spectrum.length();
List<double> newspectrum(L);

  for (int i=0; i<L/2; i++){
     newspectrum(i)=spectrum(i+L/2);
     newspectrum(i+L/2)=spectrum(i);
      }
    spectrum=newspectrum;
    }//*********************************************************************************************************************************************************************************************************



int main(int argc,const char *argv[]) {


//Define spin system
try {

cout.precision(3);
char fname[128]; //filename
int count=1;
spinhalf_system lsys(2, "1H"); //define two protons system
getstring(argc,argv,count,"Output file? ",fname,128);
rmatrix H;
cmatrix sigma0=F(lsys,'+'); //define start density matrix F+
//cmatrix detect_operator=F(lsys,'+');

cout<<"F+=\n"<< sigma0<<endl;

const double gam=gamma("1H"); //magnetogyric ratio for protons


//Define Hamiltonian

const double d=dipolar_coupling(gam, gam, 1.5*1e-10); //calculate dipolar coupling for 2 angstrom

cout<<"d="<<d<<endl;

//const rmatrix spin_part=real(spin_dipolar(lsys,0,1)); //spin part of homonuclear dipolar hamiltonian

//const RotorInfo MASinfo(2,MAGIC_ANGLE);

//List<double> cshifts(lsys.nspins());
//cshifts(0)=0;
//cshifts(1)=4000;
//const List<double> Hcs=diag_HCS(lsys, cshifts);
//cout<<"HCS: \n"<<Hcs<<endl;


//const cmatrix spin=I(lsys,0,'z')*I(lsys,1,'z');
//cout<<"spin=\n"<<spin_part<<endl;



//Define MAS parameters

const int nobs=16; //number of observation per rotor period
const int gammasteps=32; //number of gamma steps per rotor period
//const int step_ratio=gammasteps/nobs; //number of gammasteps per 1 observation interval
const int n_rotor_cycles=32; //number of full rotor cycles during FID recording
const double nu_r=10000; //MAS rate in Hz
const int spec_length=nobs*n_rotor_cycles; //number of points in FID
const double sw=nu_r*nobs; //Spectrum width in Hz
//const double intdt=1/(nu_r*gammasteps*16); //integration step = 16 times gammastep
const double rotor_period=1/nu_r;

H=real(spin_dipolar(lsys,0,1));
cout<<endl<<"Hamiltonian: "<<endl<<H<<endl;
space_T spat_MF=spatial_tensor(d); //spatial part of Hamiltonian in principal axis
cout<<"spat=\n"<<spat_MF<<endl;
const int nphases=  gammasteps;


//Powder avaraging cycle

PlanarZCW zcw(7); //Powder method - ZCW
Euler powder; //Euler - angles of single powder orientation
double weight; //Statistical weight of the orientation

List<double> spectrum(spec_length, 0.0); //free space for spectrum storage
FoldingHistogram<double> spechist(spectrum,sw);

DynamicPhase phasecalc(nu_r,0.0,nphases);
GammaInhomogeneousSpectrumED obj(nobs);
obj.set_H(H);
List<double> sb_amp;

//Powder cycle

  cout << "sigma=\n"<<sigma0<<endl;


while (zcw.next(powder, weight)){
cout<<"POWDER Angles: \n"<<powder<<endl;
space_T  spat_RF=rotate(spat_MF, powder);//spatial tensor in rotor frame
cout<<spat_RF<<endl;
phasecalc.tensor(spat_RF);
cout<<"Dynamic phases: "<<phasecalc<<endl;
obj.set_phases(phasecalc);
obj.observe(sigma0);

    //output
   double amp,freq;
   while (obj(amp,freq)) {
   if (fabs(amp)>1e-8) {
      amp*=weight;
      cout << amp << " of " << freq << " Hz \n";
      spechist.add(amp,freq);
    }
    }

	obj.add(sb_amp,0,1);
	obj.add(sb_amp,2,3);

	cout << "Intensity pattern: " << sb_amp << '\n';
 }


//cout << "sideband pattern "<<sb_amp;
//while (zcw.next(powder,weight)) {

  // cout<<"POWDER Angles: \n"<<powder<<endl;

    //space_T  spat_RF=rotate(spat_MF, powder);//spatial tensor in rotor frame
//    cout<<spat_RF<<endl;

//    RealSpinningHamiltonian H_MAS;  //(MASinfo);
//    	cout << "spin part\n" << spin_part << '\n';
//    H_MAS.add(spat_RF,spin_part);  //Entire hamiltonian
   // H_MAS.add(Hcs);

//    cout<<"H=\n"<<H_MAS<<endl;
//    MASPropagator<rmatrix> Ugen(H_MAS,nu_r,0.0,intdt);

//    List<cmatrix> U(gammasteps);
//    propagators(U, Ugen, 0.0, rotor_period);
    //cout << "U=\n"<<U<<endl;
//    GammaPeriodicSpectrumED obj(nobs);
//    obj.set_Us(U, rotor_period);
//    obj.observe(sigma0);
//    cout << "sigma=\n"<<sigma0<<endl;

    //output
//   double amp,freq;
 //  while (obj(amp,freq)) {
 //   if (fabs(amp)>1e-8) {
//      amp*=weight;
//      cout << amp << " of " << freq << '\n';
//      spechist.add(amp,freq);
 //   }
//    }
// }

 char tmpname[80];

  sprintf(tmpname,"%s.spe",fname);


  halfshift(spectrum);

  //cmatrix spec_com(1,spectrum.length());
  //for (int i=0; i<spectrum.length(); i++) spec_com(0,i)=complex(spectrum(i));
  //write_simpson(tmpname, spec_com, sw, 1);
  write_simpson(tmpname,List<complex>(spectrum), sw, true);

 } catch(MatrixException& exc) {
	cerr << exc << '\n';
	return 1;
}
	return 0;
 }


