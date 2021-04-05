/* Two site exchange under MAS */

#include "ttyio.h"
#include "MAS.h"
#include "powder.h"
#include "simpsonio.h"
#include "MetaPropagation.h"
#include "superop.h"

using namespace std;
using namespace libcmatrix;

static const double deg_to_rad=M_PI/180.0;

template<class HType> class BlockedLiouvillian {
public:
  BlockedLiouvillian(size_t nsites =1U) : hams_(nsites) {} //!< construct for n-sites (HType has default constructor)
  template<class CType> BlockedLiouvillian(const CType& base, size_t nsites =1U) : hams_(int(nsites),base) {} //!< construct for n-sites (from object)
  BlockedLiouvillian(const rmatrix& X) : hams_(X.rows()) { exchange(X); } //!< build Liouvillian based on site exchange matrix \a X
  template<class CType> BlockedLiouvillian(const CType& base, const rmatrix& X) : hams_(X.rows(),base) { exchange(X); } //!< build Liouvillian based on site exchange matrix \a X

  double period() const;

  HType& hamiltonian(size_t which =0U) {
    Ltmp_.clear();
    return hams_(which);
  }
  const HType& hamiltonian(size_t which =0U) const { return hams_(which); }
  void clear() { exchange_.clear(); hams_.clear(); Ltmp_.clear(); }
  void propagator(BlockedMatrix<complex>&, double, double, double intdt =0) const;
  void liouvillian(BlockedMatrix<complex>&, double t =0.0) const; //!< return Liouvillian for given time point
  void liouvillian(cmatrix& dest, double t, size_t blk =0U) const; //!< return Liouvillian block for given time point
  // void size() const { return hams_.size()*hams_.front().row().size(); }  //!< return dimensionality (0 if undefined)
  void sanitycheck() const; //!< throw exception if hamiltonian sizes etc. are inconsistent
  void exchange(const rmatrix&);

private:
  List<HType> hams_;
  rmatrix exchange_;
  mutable List<size_t> bsizes_; //!< block sizes
  BlockedMatrix<complex> Ltmp_;
  void propagator_(BlockedMatrix<complex>&, double t1, double t2, double, Bool2Type<false>) const; //!< propagator generator for time-dependent problems
  void propagator_(BlockedMatrix<complex>&, double t1, double t2, double, Bool2Type<true>) const; //!< propagator generator for constant Hamiltonians
  void ensure(BlockedMatrix<complex>&) const;
  void getperiod() const; //!< check period

  //void liouvillian_(BlockedMatrix<complex>&, double t, Bool2Type<false>) const; //!< return Liouvillian for given time point
  //void liouvillian_(BlockedMatrix<complex>&, double t, Bool2Type<true>) const; //!< return Liouvillian for given time point
};

//!< helper class to extract Liouvillian sub block
template<class HType> class LiouvillianBlock : public UnaryFunction<cmatrix,double> {
public:
  LiouvillianBlock(const BlockedLiouvillian<HType>& L, size_t blkv =0U)
    : L_(L), blk_(blkv), period_(L.period()) {}

  double period() const { return period_; }
  void operator()(cmatrix& dest, double t) const { L_.liouvillian(dest,t,blk_); }
  cmatrix operator()(double t) const {
    cmatrix d(mxflag::temporary);
    (*this)(d,t);
    return d;
  }
private:
  const BlockedLiouvillian<HType>& L_;
  size_t blk_;
  double period_;
};

template<class HType> double BlockedLiouvillian<HType>::period() const
{
  if (hams_.empty())
    throw Undefined("BlockedLiouvillian::period");
  double p=-1.0;
  for (size_t i=hams_.size();i--;) {
    if (p<0.0)
      p=hams_(i).period();
    else {
      if (hams_(i).period()!=p)
	throw Failed("BlockedLiouvillian::period: different Hamiltonians have different periods!");
    }
  }
  return p;
}


// template<class HType> size_t BlockedLiouvillian<HType>::size() const 
// { 
//   const size_t nbasis=hams_.front().size();
//   return hams_.size()*nbasis*nbasis;
// }

template<class HType> void BlockedLiouvillian<HType>::sanitycheck() const
{
  if (!exchange_.empty() && (exchange_.rows()!=hams_.size()))
    throw Mismatch("BlockedLiouvillian: mismatch between exchange matrix and number of site Hamiltonians");
  if (hams_.empty())
    throw Undefined("BlockedLiouvillian: Hamiltonians undefined");
  size_t s=hams_.back().size();
  for (size_t n=hams_.size()-1;n--;) {
    if (hams_(n).size()!=s)
      throw Mismatch("BlockedLiouvillian: Hamiltonians have inconsistent dimensions");
  }
}

template<class HType> void BlockedLiouvillian<HType>::ensure(BlockedMatrix<complex>& dest) const
{
  const HType& Href(hams_.front()); //!< must assume all Hamiltonians have same structure
  bsizes_.create(Href.totalblocks());
  const size_t nsites=hams_.size();
  for (size_t i=bsizes_.size();i--;) {
    const size_t s=Href(i).size();
    bsizes_(i)=s*s*nsites;
  }
  dest.create(bsizes_);
}
  
template<class HType> void BlockedLiouvillian<HType>::propagator_(BlockedMatrix<complex>& dest, double t1, double t2, double, Bool2Type<true>) const //!< propagator generator for time-independent problems
{
  if (!Ltmp_)
    liouvillian(Ltmp_); //!< time is irrelevant
  dest.duplicate_structure(Ltmp_);
  const double dt=t2-t1;
  for (size_t i=Ltmp_.size();i--;) {
    cmatrix cdest(dest(i));
    propagatorL(cdest,Ltmp_(i),dt);
  }
}

template<class HType> void BlockedLiouvillian<HType>::propagator_(BlockedMatrix<complex>& dest, double t1, double t2, double intdt, Bool2Type<false>) const //!< propagator generator for time-dependent problems
{
  if (intdt<=0.0)
    throw InvalidParameter("BlockedLiouvillian: integration time step must be >0");
  ensure(dest);
  for (size_t i=dest.size();i--;) {
    LiouvillianBlock<HType> Lblk(*this,i);
    HomogeneousPropagator<> propgen(Lblk,intdt,HomogeneousPropagator<>::LiouvilleSpace);
    cmatrix cdest(dest(i));
    propgen(cdest,t1,t2);
  }
}

template<class HrType, class HType> void getH(HrType& dest, const HType& H, size_t blk, double t,Bool2Type<false>)
{
  size_t mzblk,eigblk;
  H.reverse(mzblk,eigblk,blk);
  H(mzblk,eigblk)(dest,t);
}

template<class HrType, class HType> void getH(HrType& dest, const HType& H, size_t blk, double,Bool2Type<true>)
{
  dest=H(blk);
}
	
template<class HType> void getsuperH(rmatrix& superHtmp, const HType& Hcur, size_t blk, double t, Bool2Type<false>)
{
  rmatrix Htmp;
  getH(Htmp,Hcur,blk,t,Bool2Type< Ham_traits<HType>::isconstant >());
  commutator(superHtmp,Htmp);
}

template<class HType> void getsuperH(rmatrix& superHtmp, const HType& Hcur, size_t blk, double t, Bool2Type<true>)
{
  List<double> Htmp,superHtmpd;
  getH(Htmp,Hcur,blk,t,Bool2Type< Ham_traits<HType>::isconstant >());
  commutator(superHtmpd,Htmp);
  full(superHtmp,superHtmpd);
}
				 
template<class HType> void BlockedLiouvillian<HType>::liouvillian(cmatrix& dest, double t, size_t blk) const
{
  rmatrix superHtmp; //!< temporary for H superoperator
  cmatrix Ltmp;
  static const complex MINUS_TWOPI_I(0.0,-2.0*M_PI);
  dest=0.0;
  size_t ptr=0U;
  const size_t nsites=hams_.size();
  const bool haveX=!!exchange_;
  for (size_t site=0;site<nsites;site++) {
    const HType& Hcur(hams_(site));
    getsuperH(superHtmp,Hcur,blk,t,Bool2Type< Ham_traits<HType>::isdiagonal >());
    multiply(Ltmp,MINUS_TWOPI_I,superHtmp);
    const size_t Ldim=Ltmp.rows();
    const range rsel(ptr,ptr+Ldim-1);
    dest(rsel,rsel)=Ltmp;
    if (haveX) {
      for (size_t csite=0;csite<nsites;csite++) {
	const double kval=exchange_(site,csite);
	const size_t rptr=site*Ldim;
	const size_t cptr=csite*Ldim;
	for (size_t i=Ldim;i--;)
	  dest(rptr+i,cptr+i)+=kval;
      }
    }
    ptr+=Ldim;
  }
}

template<class HType> void BlockedLiouvillian<HType>::propagator(BlockedMatrix<complex>& dest, double t1, double t2, double intdt) const
{
  if (t2<t1)
    throw InvalidParameter("BlockedLiouvillian: propagator interval is <0");
  if (t1==t2) {
    dest.clear();
    return;
  }
  propagator_(dest,t1,t2,intdt,Bool2Type< Ham_traits<HType>::isconstant >());
}

template<class HType> void BlockedLiouvillian<HType>::liouvillian(BlockedMatrix<complex>& dest, double t) const
{
  sanitycheck();
  ensure(dest);
  for (size_t blk=dest.size();blk--;) {
    cmatrix& cdest(dest(blk));
    liouvillian(dest,t,blk);
  }
}

template<class HType> void BlockedLiouvillian<HType>::exchange(const rmatrix& X)
{
  if (!issquare(X))
    throw InvalidParameter("BlockedLiouvillian: exchange matrix must be non-empty and square");
  // check validity?
  if (X.rows()!=hams_.size()) {
    //    if (hams_.empty())
    //   hams_.create(X.rows());
    //else
    throw Mismatch("BlockedLiouvillian: exchange matrix doesn't match number of sites");
  }
  exchange_=X;
  Ltmp_.clear();
}

int main(int argc, const char* argv[])
{
  int count=1;

  const double aniso=(200.0*2.0/3)*470; //!< CSA (frequency units) - factor of 2/3 "corrects" CASTEP-derived CSA
  const double angle=getfloat(argc,argv,count,"Jump angle? ",40.0)*deg_to_rad;

  const double rotor_speed=8e3;
  //const size_t steps=64;

  const space_T A_PAS(spatial_tensor(aniso));
  const space_T A_MF1(rotate(A_PAS,Euler(0.0,M_PI/2,-angle/2)));
  space_T A_MF2;

  if (getlogical(argc,argv,count,"Include fast averaging? ")) {
    const space_T A_MF2a=rotate(A_PAS,Euler(0.0,M_PI/2,angle));
    const space_T A_MF2b=rotate(A_PAS,Euler(0.0,M_PI/2,0.0));
    A_MF2=A_MF2a+A_MF2b;
    //A_MF2+=A_MF2b;
    A_MF2*=0.5; //!< average over two orientations
  }
  else
    A_MF2=rotate(A_PAS,Euler(0.0,M_PI/2,angle/2));

  cout << "Tensor 1 (MF):\n" << A_MF1;
  cout << "Tensor 2 (MF):\n" << A_MF2;

  Euler RF_to_LF(0.0,MAGIC_ANGLE,0.0);

  //  const double R=0.0;
  const size_t ksteps=getint(argc,argv,count,"Exchange rate steps? ",20);  
  assert(ksteps>0);
  double kstart,kend;
  size_t nacq=1;
  //  double refLW=0.0;
  if (ksteps<2) {
    kstart=kend=getfloat(argc,argv,count,"Exchange rate (Hz) ? ");
    nacq=16;
  }
  else {
    kstart=getfloat(argc,argv,count,"Starting exchange rate (Hz)? ");
    kend=getfloat(argc,argv,count,"Finish exchange rate (Hz)? ");
    //    refLW=getfloat(argc,argv,count,"Reference linewidth (Hz)? ",100);
  }   
  cout.precision(5);
  //  const double refR=M_PI*refLW;
    
  //const size_t nsteps=16; //!< steps per rotor period

  double k=kstart;
  List<complex> S(nacq);
  Matrix<double> Results(ksteps,2,0.0);
  const size_t nzcw=getint(argc,argv,count,"Powder quality? ",1);
  PlanarZCW powdmeth(nzcw);

  //  const double TWOPI=2.0*M_PI;
  lcm_pade_use=getlogical(argc,argv,count,"Use Pade approximants? ");

  const double kfac=(ksteps>1) ? pow(kend/kstart,1.0/(ksteps-1)) : 1.0;

  HamiltonianStore<space_T> HstoreA(1);
  HamiltonianStore<space_T> HstoreB(1);
  const spin_system sys(1,"2H");
  SpinOpGenerator opgen(sys);
  
  for (size_t kstep=0;kstep<ksteps;kstep++,k*=kfac) {
    cout << "Jump rate: " << k << " Hz\n";
    S=complex(0.0);
    
    rmatrix X(2,2);
    X(0U,0U)=-k;
    X(0U,1U)=k;
    X(1U,0U)=k;
    X(1U,1U)=-k;

    //    const double dt=1.0/(rotor_speed*nsteps);
    
    cmatrix Utmp;
    BlockedMatrix<complex> UB;
    
    Euler powder(0.0,0.0,0.0);
    
    double weight;
    powdmeth.reset();
    typedef BlockedDiagonalSpinningHamiltonian htype;
    const htype baseH(opgen,rotor_speed,0.0);
    BlockedLiouvillian<htype> Liouv(baseH,X);

    //double scalefac=0.0;
    while (powdmeth.next(powder,weight)) {
      
      //scalefac+=weight;
      const space_T A_RF1(rotate(A_MF1,powder));
      const space_T A_RF2(rotate(A_MF2,powder));
      
      cmatrix U;

 //      for (size_t i=0;i<steps;i++) {
// 	RF_to_LF.alpha=TWOPI*(0.5+i)/steps;
// 	const double w1=-M_PI*real(rotate(A_RF1,2,0,RF_to_LF));
// 	const double w2=-M_PI*real(rotate(A_RF2,2,0,RF_to_LF));
// 	//if (nzcw==1)
// 	//  std::cout << w1 << '\t' << w2 << '\n';
	
// 	L(0,0)=complex(-k-R,w1);
// 	L(1,1)=complex(-k-R,w2);

// 	propagatorL(Utmp,L,dt);
// 	//	Utmp=exp(L,dt);
// 	U&=Utmp;
//       }
      
      HstoreA.set_shift(0U,A_RF1);
      HstoreB.set_shift(0U,A_RF2);
      
      Liouv.hamiltonian(0U).interactions(HstoreA);
      Liouv.hamiltonian(1U).interactions(HstoreB);
      Liouv.propagator(UB,0.0,1.0/rotor_speed);

      const ScratchList<double,2> M0(2,weight*0.5);
      List<complex> M;
    
      if (nacq==1) {
	multiply(M,U,M0);
	const complex val(sum(M));
	S.front()+=val;
	//	cout << "Ratio: " << (val/weight) << '\n';
      }
      else {
	if (nzcw==1)
	  cout << "Rotor period propagator:\n" << U << '\n';
	
	M=M0;
	List<complex> tmp;
	S.front()+=sum(M);
	for (size_t i=1;i<nacq;i++) {
	  multiply(tmp,U,M);
	  M.swap(tmp);
	  S(i)+=sum(M);
	}
      }
    }
    const double rawfac=real(S(nacq==1 ? 0U : 1U));
    const double LW=M_PI*(-log(rawfac)*rotor_speed);

    cout << "Raw reduction factor: " << rawfac << "   Additional LW: " << LW << " Hz\n";
    if (nacq==1) {
      Results(kstep,0U)=k;
      Results(kstep,1U)=LW;
    }
  }
  if (S.size()>1)
    write_simpson("siteexchange.fid",S,rotor_speed,false);
  else
    write_matrix("siteexchange_rs",Results);
  return 0;
}


  
