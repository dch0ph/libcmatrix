#undef LCM_SUPPRESS_VIEWS
#include <cstdlib>
#include "MAS.h"
#include "rmatrix.h"


namespace libcmatrix {

#include "lcm_accumulate.hpp"

  static const double deg_to_rad=M_PI/180.0;
  static const double symmetrytol=1e-6;

  void DynamicPhase::dump() const
  {
    std::cout << (*this);
  }

  DynamicPhase::DynamicPhase(double speed,double lrotor_phase,int nobs_,const RotorInfo& info_)
    : RotorSpec(info_,speed,lrotor_phase), nobs(nobs_), Bvalues(infop_->rank(),complex(0.0))
{ 
  B0=0.0;
}

 DynamicPhase::DynamicPhase(double speed,double lrotor_phase,int nobs_,const space_T& A,const RotorInfo& info_)
   : RotorSpec(info_,speed,lrotor_phase), nobs(nobs_), Bvalues(A.rank())
{ 
  tensor(A);
}

  DynamicPhase::DynamicPhase(const DiagonalSpinningHamiltonian& Hspin, size_t which,int nobs_)
  : RotorSpec(Hspin), nobs(nobs_)
{
   if (which>=Hspin.size())
    throw BadIndex("DiagonalSpinningHamiltonian->DynamicPhase");
  B0=Hspin.H0(which);
  Bvalues=Hspin.tensor_values(which,range());  
  update_tensor();
}

  DynamicPhase::DynamicPhase(const RealSpinningHamiltonian& Hspin, size_t which, int nobs_)
  : RotorSpec(Hspin), nobs(nobs_)
{
   if (which>=Hspin.size())
    throw BadIndex("RealSpinningHamiltonian->DynamicPhase");
   if (!isuniform())
     throw Failed("DynamicPhase: can't be initialised if rotor spin rate is unstable");

  B0=Hspin.H0(which,which);
  size_t m=Hspin.tensor_values.size();
  Bvalues.create(m);
  for (;m--;)
    Bvalues(m)=Hspin.tensor_values(m)(which,which);
  update_tensor();
}

  DynamicPhase::DynamicPhase(const SpinningHamiltonian& Hspin, size_t which, int nobs_)
    : RotorSpec(Hspin), nobs(nobs_), Bvalues(Hspin.rank())
{
  if (!isuniform())
    throw Failed("DynamicPhase: can't be initialised if rotor spin rate is unstable");
  Hspin.squeeze(B0,Bvalues,which,symmetrytol);
}

static void expi(BaseList<complex> U, const BaseList<double>& H)
{
  size_t n=H.size();
  if (U.size()!=n)
    throw Mismatch("expi");
  for (;n--;)
    U(n)=expi(H(n));
}

size_t DiagonalInhomogeneousPropagator::size() const
{
  if (Hs.empty()) {
    const size_t s=Hspinp->size();
    assert(s!=0);
    return s;
  }
  return Hs.cols();
}

void DiagonalInhomogeneousPropagator::operator()(List<complex>& U, double t1, double t2) const
{
  U.create(size());
  (*this)(static_cast<BaseList<complex> >(U),t1,t2);
}

List<complex> DiagonalInhomogeneousPropagator::operator() (double t1,double t2) const 
{
  List<complex> U(mxflag::temporary);
  operator()(U,t1,t2);
  return U;
};

void DiagonalInhomogeneousPropagator::operator() (BaseList<complex> U,double t1,double t2) const
{
  if (!!Hspinp) {
    ScratchList<double> Hint(Hspinp->component0());
    Hint*=MINUS_TWO_PI*(t2-t1);
    const double phase1=Hspinp->phase(t1);
    const double phase2=Hspinp->phase(t2);
    const double twodivrotorspeed=-2.0/Hspinp->rotor_speed();
    for (int m=Hspinp->rank();m>=1;m--) {
      const complex fac=(expi(-m*phase2)-expi(-m*phase1))*complex(0,twodivrotorspeed/m);
      real_mla(Hint,fac,Hspinp->component(m));
    }
    expi(U,Hint);
    return;
  }

  ScratchList<complex> Utmp(Hs.cols());
  bool haveflag=false;
  for (size_t j=dphases.length();j--;) {
    dphases(j).propagator(Utmp,Hs.row(j),dphases(j)(t1,t2));
    accumulate_(U,Utmp,haveflag);
  }
  if (Hcon.length()) {
    propagator(Utmp,Hcon,t2-t1);
    accumulate_(U,Utmp,haveflag);
  }
}

void InhomogeneousPropagator::operator() (cmatrix& U,double t1,double t2) const
{
  List<complex> Utot;
  dobj(Utot,t1,t2);
  if (havetransform()) {
    if (!VC)
      unitary_simtrans(U,Utot,VR);
    else
      unitary_simtrans(U,Utot,VC);
  }
  else 
    full(U,Utot);
}

// reject tensors with anti-symmetric components
void DynamicPhase::tensor(const space_T& A)
{
  const int ran=A.max_rank();
  if (ran > rank())
    throw Mismatch("DynamicPhase::tensor");

  B0=0.0;
  if (ran==0)
    Bvalues.clear();
  else
    Bvalues.create(rank(),complex(0.0));

  const Tensor<double> &dvals=infop_->dvalues();

  for (int l=0;l<=ran;l++) {
    if (A.have_rank(l)) {
      if (l & 1)
		throw Failed("DynamicPhase::tensor");
      else
		B0+=dvals(l,0)*real(A(l,0));
	  for (int m=1;m<=l;m++)
		mla(Bvalues(m-1),dvals(l,m),A(l,m));
    }
  }
  update_tensor();
}

double DynamicPhase::anisotropic(double start,double end) const
{
  double dphase=0;

  const double diff_phase=omegar_*(end-start)/2;
  const double mid_phase(phase((start+end)/2));
  //std::cout << (mid_phase*180/M_PI) << '\n';
  const double scalefac=-4/rotor_speed_;

  //std::cout << "Mid: " << (mid_phase*180.0/M_PI) << "  Diff: " << (diff_phase*180.0/M_PI) << '\n';

  for (int m=Bvalues.length();m>0;m--) {
    const complex fac=(::std::sin(m*diff_phase)/m)*expi(-m*mid_phase);
    //std::cout << scalefac << "  " << (sin(m*diff_phase)/m) << "  " << expi(-m*mid_phase) << '\n';
    //std::cout << m << "  " << Bvalues(m-1) << "  " << (scalefac*fac) << '\n';
    real_mla(dphase,Bvalues(m-1),fac);
    //real_mla_conj(phase,B[-m],fac);
  }
  dphase*=scalefac;
  return dphase;
}

void DynamicPhase::observations(int n)
{
  if (n==observations())
    return;

  if (n==0) {
    factor_table.clear();
    phases_table.clear();
  }
  else {
    nobs=n;
    update();
  }
}

void DynamicPhase::make_factortable()
{
  const size_t r=Bvalues.size(); 
  if ((nobs==0) || (r==0)) {
    factor_table.clear();
    return;
  }

  factor_table.create(nobs,r);
  
  const double cfac=-4/rotor_speed_;
  
  for (size_t i=0;i<nobs;i++) {
    const double diff_phase=(M_PI*i)/nobs;
    const double mid_phase=diff_phase+rotor_phase_;

    //std::cout << i << "  Mid: " << (mid_phase*180.0/M_PI) << "  Diff: " << (diff_phase*180.0/M_PI) << '\n';

    for (int m=1;m<=r;m++) {
      //std::cout << cfac << "  " << (sin(m*diff_phase)/m) << "  " << expi(-m*mid_phase) << '\n';
      factor_table(i,m-1)=(cfac*::std::sin(m*diff_phase)/m)*expi(-m*mid_phase);
      //std::cout << m << "  " << factor_table(i,m-1) << '\n';
    }
  }
}

void DynamicPhase::update() // update after change of speed, rotor_phase
{
  if (nobs==0) return;
  make_factortable();
  update_tensor();
}

void DynamicPhase::update_tensor() // update after change of tensor
{
  const size_t n=observations();
  if (n==0) return;

  const int ran=Bvalues.size();
  if (ran==0) {
    factor_table.clear();
    phases_table.create(n,0.0);
    return;
  }

  phases_table.create(n);
  phase_step=MINUS_TWO_PI/(n*rotor_speed_);

  //ensure factor_table exists
  if (!factor_table)
    make_factortable();

  for (size_t i=n;i--;) {
    double dphase=0;
    //    std::cout << i << '\n';
    for (int mm1=0;mm1<ran;mm1++) {
      //std::cout << (mm1+1) << "  " << Bvalues(mm1) << "  " << factor_table(i,mm1) << '\n';
      real_mla(dphase,Bvalues(mm1),factor_table(i,mm1) );
    }
    phases_table(i)=dphase;
  }
}

complex DynamicPhase::component(int m) const
{ 
  if (m==0)
    return complex(B0);
  if (m<0)
    return ((m & 1) ? -conj(Bvalues(-m-1)) : conj(Bvalues(-m-1)));
  else
    return Bvalues(m-1);
}

double DynamicPhase::instant_phase(double t) const
{
  const int ran=Bvalues.size();
  if (ran==0) 
    return B0;

  const double rphase(phase(t));
  double lphase=B0;
  const complex cfac=expi(rphase);
  complex fac=cfac;
  for (int mm1=0;mm1<ran;mm1++) {
    real_mla(lphase,Bvalues(mm1),fac);
    fac*=cfac;
  }
  return lphase;
}

void DynamicPhase::propagator(BaseList<complex> U,const BaseList<double>& eigs,double lphase) const
{
  size_t m=eigs.length();
  if (m!=U.length())
    throw Mismatch("DynamicPhase::propagator");
  for (;m--;)
    U(m)=expi(eigs(m)*lphase);
}

void DynamicPhase::propagators(cmatrix& Us,const BaseList<double>& eigs) const
{
  if (!nobs)
    throw Failed("DynamicPhase: no set number of observations"); 
  Us.create(nobs,eigs.size());
  for (size_t n=nobs;n--;)
    propagator(Us.row(n),eigs,(*this)(n+1));
}

void DynamicPhase::propagators(BaseList<complex> Us) const
{
  if (Us.size()!=nobs)
    throw Mismatch("DynamicPhase::propagators");
  for (size_t n=nobs;n--;)
    Us(n)=expi((*this)(n+1));
}

std::ostream& operator << (std::ostream& ostr,const DynamicPhase& a)
{
  ostr << "Angle: " << a.infop_->angle()/deg_to_rad << " degrees\n";
  ostr << "Speed: " << a.rotor_speed() << " Hz\n";
  ostr << "Initial phase: " << a.rotor_phase()/deg_to_rad << " degrees\n";
  ostr << "B0: " << a.B0;
  for (int m=1;m<=a.Bvalues.size();m++)
    ostr << "   B" << m << ": " << a.Bvalues(m-1);
  ostr << std::endl;
  ostr << "Observations per rotor cycle: ";
  if (a.observations())
    ostr << a.observations() << std::endl;
  else
    ostr << "unset\n";
  return ostr;
}
  
  void DiagonalSpinningHamiltonian::dump() const
  {
    std::cout << (*this);
  }

void DiagonalSpinningHamiltonian::clear()
{
  H0.clear();
  tensor_values.clear();
}

size_t DiagonalSpinningHamiltonian::size() const
{
  if (!(*this))
    throw Undefined("DiagonalSpinningHamiltonian");
  return H0.empty() ? tensor_values.cols() : H0.size();
}

DiagonalSpinningHamiltonian::DiagonalSpinningHamiltonian(const RealSpinningHamiltonian& Hspin,size_t which)
  : SpinningHamiltonianBase< List<double> >(Hspin)
{
  if (!isuniform())
    throw Failed("DiagonalSpinningHamiltonian: can't be initialised if rotor spin rate is unstable");

  if (which>=Hspin.size())
    throw BadIndex("DiagonalSpinningHamiltonian");
  H0.create(1,Hspin.H0(which,which));
  size_t m=Hspin.tensor_values.size();
  tensor_values.create(1,m);
  BaseList<complex> asrow(tensor_values.row());
  for (;m--;)
    asrow(m)=Hspin.tensor_values(m)(which,which);
}

// void DiagonalSpinningHamiltonian::addQ2_(size_t n, double scalef, const space_T& V, const BaseList<double>& ops, bool verbose)
// {
//   const cmatrix& d2s(d2values());
//   for (int m=5;m--;)
//     tmp_(m)=d2s(m,n+2)*V(2,m-2);

//   coeffs_=complex(0.0);
//   for (int m=5;m--;) {
//     for (int mp=5;m--;) {
//       int diffm=m-mp;
//       if (diffm>=0)
// 	conj_mla(coeffs_(diffm),tmp_(mp),tmp_(m));
//     }
//   }
//   coeffs_*=scalef;
//   for (size_t m=5;m--;) {
//     if (verbose)
//       std::cout << "Second order correction from n=" << n << " to rank " << m << ": " << coeffs(m) << " H\n";
//     if (m==0)
//       real_mla(component0(),coeffs_.front(),ops);
//     else
//       mla(component(m),coeffs_(m),ops);
//   }
// }
      
// void DiagonalSpinningHamiltonian::addQ2(double scalef,const space_T& V, const BaseList<double>& ops1, const BaseList<double>& ops2, bool verbose)
// {
//   if (rank()!=4)
//     throw Failed("DiagonalSpinningHamiltonian::addQ2: Hamiltonian is not rank 4");
//   addQ2_(1,scalef,V,ops1,verbose);
//   addQ2_(2,-scalef,V,ops2,verbose);
// }

DiagonalSpinningHamiltonian::DiagonalSpinningHamiltonian(const SpinningHamiltonian& Hspin,size_t which)
  : SpinningHamiltonianBase< List<double> >(Hspin), H0(1), tensor_values(1,Hspin.rank())
{
  if (!isuniform())
    throw Failed("DiagonalSpinningHamiltonian: can't be initialised if rotor spin rate is unstable");
  BaseList<complex> asrow(tensor_values.row());
  Hspin.squeeze(H0.front(),asrow,which,symmetrytol);
}

DiagonalSpinningHamiltonian::DiagonalSpinningHamiltonian(const DiagonalSpinningHamiltonian& Hspin,const BaseList<size_t>& which)
  : SpinningHamiltonianBase< List<double> >(Hspin)
{
  if (!isvalid_indexlist(which,Hspin.size()))
    throw BadIndex("DiagonalSpinningHamiltonian");
  H0=Hspin.H0(which);
  tensor_values=Hspin.tensor_values(which,range());
}

DiagonalSpinningHamiltonian::DiagonalSpinningHamiltonian(const RotorSpec& spec_)
  : SpinningHamiltonianBase< List<double> >(spec_) { 
  if (!isuniform())
    throw Failed("DiagonalSpinningHamiltonian: can't be initialised if rotor spin rate is unstable");
  clear(); 
}

DiagonalSpinningHamiltonian::DiagonalSpinningHamiltonian(const DiagonalSpinningHamiltonian& Hspin,size_t which)
  : SpinningHamiltonianBase< List<double> >(Hspin)
{
  if (which>=Hspin.size())
    throw BadIndex("DiagonalSpinningHamiltonian");
  H0.create(1,Hspin.H0(which));
  tensor_values.create(1,Hspin.tensor_values.cols());
  BaseList<complex> asrow(tensor_values.row());
  asrow=Hspin.tensor_values(which,range());
}

void DiagonalSpinningHamiltonian::add(const space_T& A, const BaseList<double>& H)
{
  if (hasodd(A))
    throw InvalidParameter(NOODD);
  const int userank=A.max_rank();
  if (userank>RotorSpec::rank())
    throw Mismatch("DiagonalSpinningHamiltonian::add");

  const Tensor<double> &dvals=infop_->dvalues();
  
  if (!isdefined(tensor_values))
    tensor_values.create(RotorSpec::rank(),H.length(),complex(0.0));

  for (int m=0;m<=userank;m++) {
    complex tot(0.0);
    int l=m;
    if (l & 1)
      l++;
    for (;l<=userank;l+=2) {
      if (A.have_rank(l))
	mla(tot,dvals(l,m),A(l,m));
    }
    //std::cout << "Adding " << tot << " of " << m << '\n';
    if (m) {
      BaseList<complex> crow=component(m);
      mla(crow,tot,H);
    }
    else
      mla(H0,real(tot),H);
  }
}

std::ostream& operator<< (std::ostream& ostr, const DiagonalSpinningHamiltonian &a)
{
  if (!a) 
    return ostr << "<empty>\n";

  ostr << "H_0: " << a.component0() << '\n';
  for (int m=1;m<=a.rank();m++)
    ostr << "H_" << m << ": " << a.component(m) << '\n';
  return ostr;
}

void DiagonalSpinningHamiltonian::operator() (BaseList<double> Htot,double t) const
{
  if (!(*this))
    throw Undefined("DiagonalSpinningHamiltonian()");
  if (Htot.size()!=size())
    throw Mismatch("DiagonalSpinningHamiltonian()");

  if (H0.length())
    Htot=H0;
  else
    Htot=0.0;

  const double iphase(phase(t));
  for (int m=1;m<=rank();m++)
    real_mla(Htot,2.0*expi(-m*iphase),tensor_values.row(m-1));
}

  DiagonalSpinningHamiltonian& DiagonalSpinningHamiltonian::operator+= (const DiagonalSpinningHamiltonian& a) {
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("DiagonalSpinningHamiltonian::add");
    H0+=a.H0;
    tensor_values+=a.tensor_values;
    return *this;
  }
  DiagonalSpinningHamiltonian& DiagonalSpinningHamiltonian::operator-= (const DiagonalSpinningHamiltonian& a) {
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("DiagonalSpinningHamiltonian::add");
    H0-=a.H0;
    tensor_values-=a.tensor_values;
    return *this;
  }

template<class HType> void DiagonalSpinningHamiltonian::interactions_(const cmatrix Hm[5], const HType& Ham, const BaseList<double>& Hzeeman, const BaseList<size_t>& order, double degfac)
{
  const size_t dim(Hm[0].rows());
  const int L=Ham.rank();
  tensor_values.create(4,dim); 
  for (size_t m=1;m<=4;m++) {
    BaseList<complex> dest(tensor_values.row(m-1));
    if (m<=L)
      dest=diag(Ham.component(m));
    else
      dest=complex(0.0);
  }

#ifndef NDEBUG
  for (int m=0;m<=4;m++)
    ::std::cout << "Hm(" << m << ")\n" << Hm[m];

  ::std::cout << "First order Hamiltonian:\n" << (*this);
#endif
  
  if (degfac<=0.0)
    throw InvalidParameter("DiagonalSpinningHamiltonian: degeneracy criterion cannot be <=0");
  
  //Apply second order correction
  for (size_t i=dim;i--;) {
    const size_t acti=order(i);
    for (size_t j=dim;j--;) {
      if (i==j)
	continue;
      
      const size_t actj=order(j);
      const double diff=Hzeeman(i)-Hzeeman(j);

      if (fabs(diff)<degfac) {
	for (size_t m=0;m<=4;m++) {
	  if (norm(Hm[m](acti,actj))>degfac) //fail if cross-term 
	    throw Failed("Degenerate Zeeman eigenvalues");
	}
      }
      else {
	const double fac=1.0/diff;
	H0(acti)+=real(Hm[0](acti,actj))*fac;
	for (size_t m=1;m<=4;m++) 
	  tensor_values(m-1,acti)+=Hm[m](acti,actj)*fac;
      }
    }
  }
}

void DiagonalSpinningHamiltonian::interactions(const SpinningHamiltonian& Ham, const BaseList<double>& Hzeeman, const BaseList<size_t>& order, double degfac)
{
  const size_t dim(Hzeeman.size());
  if (Ham.size()!=dim)
    throw Mismatch("DiagonalSpinningHamiltonian()");
  const int L=Ham.rank();
  if ((L!=0) && (L!=2))
    throw Failed("DiagonalSpinningHamiltonian: can only construct 2nd order Hamiltonian from rank 2 SpinningHamiltonian");
  cmatrix Hm[5];
  //construct cross-terms
  cmatrix tmp;
  for (int m1=-L;m1<=L;m1++) {
    for (int m2=m1;m2<=L;m2++) {
      conj_emultiply(tmp,Ham.component(m1),Ham.component(m2));
      Hm[m2-m1]+=tmp;
    }
  }

  //Copy diagonals (rank 1)
  H0.create(dim);
  for (size_t j=dim;j--;)
    H0(j)=real(Ham.component(0)(j,j));

  interactions_(Hm,Ham,Hzeeman,order,degfac);
}

void DiagonalSpinningHamiltonian::interactions(const RealSpinningHamiltonian& Ham, const BaseList<double>& Hzeeman, const BaseList<size_t>& order, double degfac)
{
  const size_t dim(Hzeeman.size());
  if (Ham.size()!=dim)
    throw Mismatch("DiagonalSpinningHamiltonian()");
  const int L=Ham.rank();
  if ((L!=0) && (L!=2))
    throw Failed("DiagonalSpinningHamiltonian: can only construct 2nd order Hamiltonian from rank 2 SpinningHamiltonian");
  cmatrix Hm[5];
  //construct cross-terms
  cmatrix tmp;
  rmatrix tmpr;

  // Hm[0] = H0.H0 + 2( H1.H1* + H2.H2*)
  // Hm[1] = 2( H0.H1 + H1*.H2)
  // Hm[2] = H1.H1 + 2 H0.H2
  // Hm[3] = 2 H1.H2
  // Hm[4] = H2.H2
  if (L) {
    enorm(tmpr,Ham.component(1));
    Hm[0]=tmpr;
    enorm(tmpr,Ham.component(2));
    Hm[0]+=tmpr;
    Hm[0]*=2.0;
  }
  enorm(tmpr,Ham.component0());
  Hm[0]+=tmpr;

  if (L) {
    emultiply(Hm[1],Ham.component0(),Ham.component(1));
    conj_emultiply(tmp,Ham.component(1),Ham.component(2));
    Hm[1]+=tmp;
    Hm[1]*=2.0;
    
    emultiply(Hm[2],Ham.component0(),Ham.component(2));
    Hm[2]*=2.0;
    emultiply(tmp,Ham.component(1),Ham.component(1));
    Hm[2]+=tmp;
    
    emultiply(Hm[3],Ham.component(1),Ham.component(2));
    Hm[3]*=2.0;
    
    emultiply(Hm[4],Ham.component(2),Ham.component(2));
  }

  //Copy diagonals (rank 1)
  H0.create(dim);
  for (size_t j=dim;j--;)
    H0(j)=Ham.component0()(j,j);

  interactions_(Hm,Ham,Hzeeman,order,degfac);
}

}
