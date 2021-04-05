/* Check 1st and second order quadrupoles */

#include "ttyio.h"
#include "NMR.h"
#include "MAS.h"
#include "tensorop.h"
#include "Propagation.h"
#include "MetaPropagation.h"
#include "powder.h"
#include "simpsonio.h"
//#include "args_iter.h"

using namespace libcmatrix;
using namespace std;

const double deg_to_rad=M_PI/180.0;

bool verbose=false;
int propverbose=0;
size_t gamma_steps=0;
double rotor_speed=0.0;
double rotor_phase=0.0;
double intdt=0.0;
size_t nobs=0;
List<double> Hzeeman;
List<double> Hzeemant;

inline double norm(double x) { return x*x; }

void calchomo_(BaseList<complex> FID,double scale, const PropGen_t& Ugen, const cmatrix& sigma0det)
{
  static List<cmatrix> Us(gamma_steps);
  propagators(Us,Ugen,0.0,1.0/fabs(rotor_speed));
  if (verbose)
    cout << "Propagators\n" << Us;
  GammaPeriodicFID obj(nobs,propverbose);
  obj.set_Us(Us);
  obj.add_FID(FID,scale,sigma0det);
}

template<class HType> void calchomo(BaseList<complex> FID, double scale, const HType& Ham, const cmatrix& sigma0det)
{
  if (verbose) {
    cout << "Hamiltonian\n" << Ham;
    cout << "H(0)\n" << Ham(0.0);
  }
  HomogeneousPropagator<typename HType::result_type> Ugen(Ham,intdt,Hzeeman);
  calchomo_(FID,scale,Ugen,sigma0det);
}
  
template<class HType> void calchomo2nd(BaseList<complex> FID, double scale, const HType& Ham, const ListList<size_t>& blkstr, const cmatrix& sigma0det)
{
  if (verbose) {
    cout << "Hamiltonian\n" << Ham;
    cout << "H(0)\n" << Ham(0.0);
  }
  HomogeneousPropagator_second<typename HType::result_type> Ugen(Ham,intdt,blkstr,Hzeemant);
  calchomo_(FID,scale,Ugen,sigma0det);
}
  
void calcinhomo(BaseList<complex> FID, double scale, const DiagonalSpinningHamiltonian& dHam, const rmatrix& sigma0det)
{
  if (verbose) {
    cout << "Hamiltonian:\n" << dHam;
    cout << "H(0): " << dHam(0.0) << '\n';
  }
  DiagonalInhomogeneousPropagator Ugen(dHam);
  static cmatrix dUs;
  propagators(dUs,gamma_steps,Ugen,0.0,1.0/fabs(rotor_speed));
  if (verbose)
    cout << "Propagators\n" << dUs;
  GammaInhomogeneousFID obj(nobs,propverbose);
  obj.set_Us(dUs);
  obj.add_FID(FID,scale,sigma0det);
}

double larmor;
		    
int main(int argc,const char *argv[])
{
  try {
    //  args_iter getargs(1,argc,argv,cin,cout);
    int count=1;

  char name[6],name2[6];
  getstring(argc,argv,count,"Nucleus type 1? ",name,sizeof(name),"14N");
  getstring(argc,argv,count,"Nucleus type 2 (<CR> for none)? ",name2,sizeof(name2));
  const bool twospins=(name2[0]!='\0');
  const size_t nspins=twospins ? 2 : 1;
  HamiltonianStore<space_T> Hstore(nspins);
  spin_system sys(nspins,name);
  if (twospins)
    sys(1).isotope(name2);
  if (sys(0).deg()<3) {
    cerr << "Spin 1 is not quadrupolar!\n";
    return 1;
  }
  const bool haveq2=twospins ? (sys(1).deg()>2) : false;

  const size_t order=getint(argc,argv,count,"Order to use (0 for exact)? ");
  //  const bool classic = (order==2) ? getlogical(argc,argv,count,"Classic second order? ");

  if (order!=1)
    larmor=1e6*getfloat(argc,argv,count,"Larmor frequency 1 (MHz) (or 0)? ",30.0);

  const double quad=1e6*getfloat(argc,argv,count,"Quadrupole coupling 1 (MHz)? ",10.0);
  const double eta=getfloat(argc,argv,count,"Asymmetry 1? ");

  const bool nonsec=(larmor!=0.0);  
  const space_T Q_MF(spatial_tensor(quad,eta));
  Hstore.set_quadrupole(0,Q_MF,order);
    
  const double csa=1e3*getfloat(argc,argv,count,"CSA (kHz)? ",0.0);

  double beta;
  space_T CSA_MF;
  if (csa) {
    beta=deg_to_rad*getfloat(argc,argv,count,"Q -> CSA beta? ",0.0);
    CSA_MF=rotate(spatial_tensor(csa,0.0),Euler(0,beta,0));
    Hstore.set_shift(0,CSA_MF);
  }

  double larmor2=0.0;
  double quad2=0.0;
  double eta2=0.0;
  double dipole=0.0;
  double J=0.0;

  space_T Q2_MF,D_MF;
     
  if (twospins) {
    if (haveq2) {
      if (order!=1)
	larmor2=1e6*getfloat(argc,argv,count,"Larmor frequency 2 (MHz) (or 0)? ",30.0);

      quad2=1e6*getfloat(argc,argv,count,"Quadrupole coupling 2 (MHz)? ",0.2);
      beta=deg_to_rad*getfloat(argc,argv,count,"Q -> Q2 beta? ",0.0);
      Q2_MF=rotate(spatial_tensor(quad2,eta2),Euler(0,beta,0));
      Hstore.set_quadrupole(1,Q2_MF,larmor2 ? order : 1);
    }
    dipole=1e3*getfloat(argc,argv,count,"Dipole coupling (kHz)? ",1.0);
    if (dipole) {
      beta=deg_to_rad*getfloat(argc,argv,count,"Q -> D beta? ",0.0);
      D_MF=rotate(spatial_tensor(dipole,0.0,nonsec),Euler(0,beta,0));
      Hstore.set_coupling(I_DIPOLE,0,1,D_MF);
    }
    J=getfloat(argc,argv,count,"J coupling (Hz)? ",1.0);
    if (J)
      Hstore.set_coupling(I_J,0,1,space_T(0,J));
  }

  cout << "HamiltonianStore\n" << Hstore;
  HamiltonianStructure Hstruct(sys,Hstore);
  cout << "Structure\n" << Hstruct;
  const SpinOpGenerator opgen(Hstruct);
  cout << "OpGen\n" << opgen;

  rotor_speed=1e3*getfloat(argc,argv,count,"Rotor speed (kHz)? ",0.0);

  double sw=0.0;
  double misset=0.0;
  if (rotor_speed) {
    misset=deg_to_rad*getfloat(argc,argv,count,"Angle misset (degrees)? ",0.0);
    nobs=getint(argc,argv,count,"Observations per rotor period? ");
    sw=nobs*rotor_speed;
    cout << "Spectral width: " << (sw/1000.0) << " kHz\n";
    gamma_steps=getint(argc,argv,count,"Gamma steps? ",2*nobs);
    intdt=1e-6*getfloat(argc,argv,count,"Integration time step (us)? ",2.0);
  }
  else
    sw=1e3*getfloat(argc,argv,count,"Spectral width (kHz)? ");

  const size_t npts=getint(argc,argv,count,"Points in FID? ");

  const int zcw=getint(argc,argv,count,"ZCW parameter? ",1);
  verbose=(zcw<4);
  propverbose=verbose ? 2 : 0;

  PlanarZCW powdmeth(zcw);

  cout << "Orients: " << powdmeth.orientations() << endl;

  Euler powder(0,0,0);  
  if (zcw==0) {
    if (eta)
      powder.alpha=deg_to_rad*getfloat(argc,argv,count,"Q->LAB alpha? ");
    powder.beta=deg_to_rad*getfloat(argc,argv,count,"Q->LAB beta? ");
    powder.gamma=deg_to_rad*getfloat(argc,argv,count,"Q->LAB gamma? ");
  }

  char fname[128];
  getstring(argc,argv,count,"Filename? ",fname,sizeof(fname));
      
  double weight;

  //  const bool ishet=!sys.ishomonuclear();
  const bool isdiag=Hstruct.isdiagonal();
  cout << "isdiag: " << isdiag << '\n';

  RotorInfo rinfo(2,MAGIC_ANGLE+misset);

  const size_t obsspin=twospins ? 1 : 0;
  const double obslarmor=twospins ? larmor2 : larmor;
  const rmatrix sigma0det=real(I(sys,obsspin,'+'));
  cout << "sigma0det\n" << sigma0det;
  const cmatrix csigma0det(sigma0det);
  
  ListList<size_t> blkstr;
  if (larmor) {
    mla(Hzeeman,larmor,diag_Iz(sys,0));
    if (larmor2)
      mla(Hzeeman,larmor2,diag_Iz(sys,1));
    cout << "Hzeeman: " << Hzeeman << '\n';
    blkstr=find_blocks(Hzeeman,1e-5);
    cout << "Zeeman structure: " << blkstr << '\n';
    Hzeemant=Hzeeman(blkstr.row()); //re-order
    cout << "Hzeeman (blocked): " << Hzeemant << '\n';
  }

  const double dt=1.0/sw;

  List<complex> FID(npts,0.0);

  cmatrix V,Ht,sigma0dett;
  List<double> Hzeemant,seigs;

  List<complex> Us(gamma_steps);
  
  while (powdmeth.next(powder,weight)) {

    if (verbose)
      cout << "Orientation: " << powder << '\n';

    const HamiltonianStore<space_T> Hstore_MF(Hstore,powder);      
    if (verbose)
      cout << Hstore_MF;
    const HamiltonianStore<double> Hstore_MFs(Hstore,powder);

    if (order==1) {
      if (rotor_speed) {
	if (isdiag) {
	  BlockedDiagonalSpinningHamiltonian dHam(opgen,Hstore_MF,rotor_speed,rotor_phase,rinfo);
	  calcinhomo(FID,weight,dHam(0,0),sigma0det);
	}
	else {
	  BlockedSpinningHamiltonian<double> Ham(opgen,Hstore_MF,rotor_speed,rotor_phase,rinfo);
	  calchomo(FID,weight,Ham(0,0),csigma0det);
	}
	continue;
      }
      //static
      BlockedStaticHamiltonian<double> fHtot(opgen,Hstore_MFs);
      rmatrix Htot(fHtot(0,0));
      if (isdiag) {
	diag(seigs,Htot);
	Htot.clear();
      }

      StaticFID_H obj(verbose ? 2 : 0);
      
      if (!Htot) {
	if (verbose)
	  cout << "Diagonal Hamiltonian: " << seigs << endl;
	obj.set_H(seigs,dt);
      }
      else {
	if (verbose) {
	  cout << "Hamiltonian:\n" << Htot << endl;
	  spy(cout,Htot);
	}
	obj.set_H(Htot,dt);
      }
      obj.add_FID(FID,weight,sigma0det);
      continue;
    }

    // order 0 and 2
    if (rotor_speed) {
      BlockedSpinningHamiltonian<complex> Ham(opgen,Hstore_MF,rotor_speed,rotor_phase,rinfo);
      if (verbose)
	cout << "H (full)\n" << Ham;

      if (order==0) { //exact
	if (verbose) {
	  const BlockedMatrix<complex> H0(Ham(0.0));
	  cmatrix vectors;
	  List<double> eigs;
	  hermitian_eigensystem_ns(vectors,eigs,H0.front(),Hzeeman);
	  cout << "H(0) eigenvalues: " << eigs << '\n';
	}
	calchomo(FID,weight,Ham(0,0),csigma0det);
	continue;
      }
      //2nd order
      if (isdiag) { //construct DiagonalSpinningHamiltonian
	if (verbose) {
	  const BlockedMatrix<complex> H0(Ham(0.0));
	  List<double> seigs;
	  hermitian_eigenvalues_second(seigs,H0.front(),Hzeeman); //apply 2nd order corrections
	  cout << "H(0) 2nd order eigs: " << seigs << '\n';
	}
	DiagonalSpinningHamiltonian dHam(Ham(0,0),Hzeeman);
	calcinhomo(FID,weight,dHam,sigma0det);
      }
      else
	calchomo2nd(FID,weight,Ham(0,0),blkstr,csigma0det);
      continue;
    }

    //static order 0,2
    BlockedStaticHamiltonian<complex> fH(opgen,Hstore_MFs);
    cmatrix H(fH(0,0));
    if (verbose) {
      cout << "H (full)\n" << H;
      spy(cout,H);
    }
   
    StaticFID_H obj(propverbose);
    switch (order) {
    case 0: {
      if (obslarmor)
	obj.larmor(obslarmor); //select only frequencies compatible with Larmor frequency
      H+=Hzeeman;
      obj.set_H(H,dt);
      obj.add_FID(FID,weight,sigma0det);
      break;
    }
    case 2: {
      if (isdiag) {
	hermitian_eigenvalues_second(seigs,H,Hzeeman); //apply 2nd order corrections
	obj.set_H(seigs,dt);
	obj.add_FID(FID,weight,sigma0det);
      }
      else {
	hermitian_eigensystem_second(V,seigs,H,blkstr,Hzeemant);
	unitary_isimtrans(sigma0dett,sigma0det,V);
	if (verbose)
	  cout << "sigma0/detect in eigenbasis\n" << sigma0dett;
	obj.set_H(seigs,dt);
	obj.add_FID(FID,weight,sigma0dett);
      }
    }
      break;
    default:
      cerr << "Bad order: " << order << endl;
      return 1;
    }
  }

    //    if (!tdomain) FID=fft(spec,FT_BACKWARD);
  char scratch[sizeof(fname)+4];
  sprintf(scratch,"%s.fid",fname);
  write_simpson(scratch,FID,sw,false);
  exponential_multiply_ip(FID,2.5);
  FID(slice(1,FID.size()/2,2))*=-1.0;
  List<complex> eFID(FID.size()*2,0.0);
  eFID(range(0,FID.size()-1))=FID;
  fft_ip(eFID,FT_FORWARD);
  sprintf(scratch,"%s.spe",fname);
  write_simpson(scratch,eFID,sw,true);  

  } catch (MatrixException& exc) {
    cerr << exc << endl;
    return 1;
  }
  return 0;
}

//     cmatrix tmp;
//     hermitian_eigensystem(V,eigs,H);
//     cout << "Eigenvectors:\n" << V << endl;
//     multiply_conj_transpose(tmp,V,V);
//     cout << "V V'\n" << tmp;
//     cout << "Eigenvalues: " << eigs << endl;
//     cout << "Eigenvalues-Hzeeman: " << (eigs-Hzeeman) << endl;
//     unitary_simtrans(Hre,eigs,V);
//     cout << "Reconstructed H\n" << Hre;

//     cmatrix H2;
//     mla(H2,A_RF(2,0),HQs(2));
//     if (twospins) {
//       mla(H2,A2_LAB(2,0),HQ2s(2));
//       mla(H2,D_LAB(2,0),HDs(2));
//     }
//     cout << "H(1)\n" << diag(H2) << endl;
//     for (m=-2;m<=2;m++) {
//       if (m==0)
// 	continue;
//       const complex scaleQ(A_LAB(2,m)*A_LAB(2,-m)/(m*larmor));
//       multiply(tmp,HQs(m+2),HQs(2-m));
//       cout << "HQ(" << m << ")HQ(" << -m << "): " << scaleQ << "\n" << tmp;
//       mla(H2,scaleQ,tmp);

//       if (twospins) {
// 	const complex scaleQ2(A2_LAB(2,m)*A2_LAB(2,-m)/(m*larmor2));
#// 	multiply(tmp,HQ2s(m+2),HQ2s(2-m));
// 	cout << "HQ2(" << m << ")HQ2(" << -m << "): " << scaleQ2 << "\n" << tmp;
// 	mla(H2,scaleQ2,tmp);
	
// 	if (abs(m)<=1) {
// 	  const complex scaleD(A_LAB(2,m)*D_LAB(2,-m)/(m*larmor));
// 	  multiply(tmp,HQs(m+2),HDs(2-m));
// 	  cout << "HQ(" << m << ")HD(" << -m << "): " << scaleD << "\n" << tmp;
// 	  mla(H2,scaleD,tmp);

// 	  const complex scaleD2(A2_LAB(2,m)*D_LAB(2,-m)/(m*larmor2));
// 	  multiply(tmp,HQ2s(m+2),HDs(2-m));
// 	  cout << "HQ2(" << m << ")HD(" << -m << "): " << scaleD2 << "\n" << tmp;
// 	  mla(H2,scaleD2,tmp);
// 	}
//       }
//     }

//     cout << "H(2)\n" << diag(H2) << endl;

//     H-=Hzeeman;
//     cout << "H\n" << H << endl;
//     hermitian_eigenvalues1(V2,eigs1,H,Hzeeman);
//     cout << "Eigenvectors (pert. theory):\n" << V2;
//     multiply_conj_transpose(tmp,V2,V2);
//     cout << "V V'\n" << tmp;
//     cout << "Eigenvalues (first order pert. theory): " << eigs1 << endl;
//     eigs1+=Hzeeman;
//     unitary_simtrans(Hre,eigs1,V2);
//     cout << "Reconstructed H\n" << Hre;
