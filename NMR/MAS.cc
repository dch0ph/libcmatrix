#undef LCM_SUPPRESS_VIEWS
#include <cstdlib>
#include "MAS.h"
#include "NMR.h"
#include "rmatrix.h"
#include "cmatrix.h"
#include "superop.h"
#include "timer.h" //!< for prettyprint_time

namespace libcmatrix {

  const double MAGIC_ANGLE=double(54.7356103172453)*M_PI/180.0;
  static const double deg_to_rad=M_PI/180.0;
  static const double symmetrytol=1e-6;
  double lcm_MAS_timingtol=1e-9; //!< time rounding error is ~1 ns
  static const double magicangtol=1e-9; //!< max error allowed for detection of magic angle (radians)

const char NOODD[]="RealSpinningHamiltonian: odd-rank components not allowed";

  RotorInfo::RotorInfo(int rankv,double angle_,double gamma_) : _rotor_gamma(gamma_)
{
    if (rankv<1)
      throw InvalidParameter("RotorInfo: rank must be >0");
    _dvalues.create(rankv);
    angle(angle_);
}  

  RotorInfo::RotorInfo(const RotorInfo& base,int rankv)
    : _rotor_gamma(base._rotor_gamma)
  {
    if (rankv<1)
      throw InvalidParameter("RotorInfo: rank must be >0");
    _dvalues.create(rankv);
    angle(base._rotor_angle);
  }  
    
  void RotorInfo::angle(double angle_)
  {
    _rotor_angle=angle_;
    ismagic=(fabs(_rotor_angle-MAGIC_ANGLE)<magicangtol);
    for (int l=rank();l>=0;l--) {
      _dvalues.ensure_rank(l);
      for (int m=-l;m<=l;m++)
	_dvalues(l,m)=(l==2 && m==0 && ismagic) ? 0.0 : d(l,m,0,_rotor_angle);
    }
    d2values_.clear(); //invalidate d2values
  }
  
  RotorSpec::RotorSpec(const RotorSpec& spec, int rankv)
     : rotor_phase_(spec.rotor_phase_)
   {
     const RotorInfo& curinfo=*(spec.infop_);
     if (rankv>curinfo.rank())
       infop_.reset(new RotorInfo(curinfo,rankv));
     else
       infop_.reset(&curinfo,mxflag::nondynamic);
     rotor_speed(spec.rotor_speed_);
   }

  void RotorSpec::rotor_speed(double rotor_speedv)
  {
    if (rotor_speedv==0.0)
      throw InvalidParameter("Rotor speed cannot be zero");
    rotor_speed_=rotor_speedv;
    omegar_=TWO_PI*rotor_speedv;
    rotor_period_=1.0/fabs(rotor_speed_);
  }

  size_t SpinningHamiltonian::usage() const
  {
    size_t tot=0;
    for (size_t i=tensor_values.size();i--;)
      tot+=tensor_values(i).size();
    return tot;
  }

void SpinningHamiltonian::squeeze(double& H0,BaseList<complex> Hm,size_t which,double tol) const
{
  if (tol<0.0)
    throw InvalidParameter("SpinningHamiltonian: tolerance must be >=0");

  const int L=rank();
  if (L<0)
    throw Failed("SpinningHamiltonian is empty");
  const size_t dim=size();
  if ((which>=dim) || (L<0) || (L!=Hm.size()))
    throw BadIndex("SpinningHamiltonian");
  H0=real(tensor_values(L)(which,which));
  const double toltol=tol*tol;
  for (int m=L;m>0;m--) {
    const complex valp=tensor_values(L+m)(which,which);
    const complex valm=tensor_values(L-m)(which,which);    
    const double err=(m & 1) ? norm(valp+conj(valm)) : norm(valp-conj(valm));
    if (err>toltol)
      throw Failed("SpinningHamiltonian is not valid");
    Hm(m-1)=valp;
  }
}

const RotorInfo MASRotorInfo(2,MAGIC_ANGLE);

const cmatrix& RotorInfo::d2values() const
{
  if (!d2values_) {
    cmatrix& nonconst_d2values=const_cast<RotorInfo*>(this)->d2values_;
    nonconst_d2values=D(2,Euler(0,_rotor_angle,_rotor_gamma));
    if (ismagic)
      nonconst_d2values(2,2)=0.0; //!< ensure identically zero for magic angle
  }
  return d2values_;
}


void SpinningHamiltonian::clear()
{
  if (tensor_values.empty())
    tensor_values.create(2*RotorSpec::rank()+1);
  else {
    for (size_t m=tensor_values.size();m--;)
      tensor_values(m).clear();
  }
}

  void SpinningHamiltonian::dump() const
  {
    std::cout << (*this);
  }

  void RealSpinningHamiltonian::dump() const
  {
    std::cout << (*this);
  }

void RealSpinningHamiltonian::clear()
{
  H0.clear();
  tensor_values.clear();
  //int m=rank();
  //tensor_values.create(m);
  //for (;m--;)
  //  tensor_values(m).clear();
}

void RealSpinningHamiltonian::tensor(const Tensor<cmatrix>& A)
{
  const int ran=A.max_rank();
  if (ran>RotorSpec::rank())
    throw Mismatch("RealSpinningHamiltonian::tensor");
  if (hasodd(A))
    throw InvalidParameter(NOODD);
  
  clear();
  const Tensor<double>& dvals=infop_->dvalues();
  rmatrix tmpr;

  for (int l=0;l<=ran;l+=2) {
    if (!A.have_rank(l))
      continue;

    if (isdefined(A(l,0))) {
      real(tmpr,A(l,0));
      mla(H0,dvals(l,0),tmpr);
    }
    for (int m=1;m<=l;m++) {
      const cmatrix& Alm=A(l,m);
      if (isdefined(Alm)) {
	const bool madenew=ensure_tensor(Alm.rows());
	cmatrix tmp(tensor_values(size_t(m-1)));
	if (madenew)
	  multiply(tmp,dvals(l,m),Alm);
	else
	  mla(tmp,dvals(l,m),Alm);
      }
    }
  }
}

void RealSpinningHamiltonian::tensor(const space_T& A)
{
  const int ran=A.max_rank();
  if (ran>RotorSpec::rank())
    throw Mismatch("RealSpinningHamiltonian::tensor");
  if (hasodd(A))
    throw InvalidParameter(NOODD);
  
  clear();
  const Tensor<double>& dvals=infop_->dvalues();

  for (int l=0;l<=ran;l+=2) {
    if (!A.have_rank(l))
      continue;

    if (!H0) {
      H0.create(1,1);
      H0=0.0;
    }
    H0(0U,0U)+=dvals(l,0)*real(A(l,0));

    (void)ensure_tensor(1);
    for (int m=1;m<=l;m++)
      mla(tensor_values(size_t(m-1),0U,0U),dvals(l,m),A(l,m));
  }
}

bool RealSpinningHamiltonian::ensure_tensor(size_t r)
{
  if (tensor_values.empty()) {
    tensor_values.create(RotorSpec::rank(),r,r,complex(0.0));
    return true;
  }
  return false;
}

bool SpinningHamiltonian::ensure_tensor(size_t)
{
  if (tensor_values.empty()) {
    //    tensor_values.create(2*rank()+1,r,r,complex(0.0));
    tensor_values.create(2*RotorSpec::rank()+1);
    return true;
  }
  return false;
}

RealSpinningHamiltonian::RealSpinningHamiltonian(const RealSpinningHamiltonian& Hspin,const BaseList<size_t>& which)
  : SpinningHamiltonianBase<rmatrix>(Hspin)
{
  if (!isvalid_indexlist(which,Hspin.size()))
    throw BadIndex("RealSpinningHamiltonian");
  H0=Hspin.H0(which,which);
  size_t m=Hspin.tensor_values.dimension(0);
  const size_t r=which.size();
  tensor_values.create(m,r,r);
  for (;m--;)
    tensor_values(m)=Hspin.tensor_values(m)(which,which);
}

SpinningHamiltonian::SpinningHamiltonian(const SpinningHamiltonian& Hspin, const BaseList<size_t>& which)
  : SpinningHamiltonianBase<cmatrix>(Hspin)
{
  if (!isvalid_indexlist(which,Hspin.size()))
    throw BadIndex("SpinningHamiltonian()");
  size_t m=Hspin.tensor_values.size();
  tensor_values.create(m);
  for (;m--;) {
    const cmatrix& Hcur=Hspin.tensor_values(m);
    if (isdefined(Hcur))
      tensor_values(m)=Hcur(which,which);
  }
}

void SpinningHamiltonian::tensor(const Tensor<cmatrix>& A)
{
  const int ran=A.max_rank();
  if (ran>RotorSpec::rank())
    throw Mismatch("SpinningHamiltonian::tensor");
  
  clear();

  const Tensor<double>& dvals=infop_->dvalues();
  for (int l=0;l<=ran;l++) {
    if (!A.have_rank(l))
      continue;
    ensure_tensor(0); //size argument is dummy
    for (int m=-l;m<=l;m++) {
      if (isdefined(A(l,m)))
	mla(component(m),dvals(l,m),A(l,m));
    }
  }
}

void SpinningHamiltonian::tensor(const space_T& A)
{
  const int ran=A.max_rank();
  if (ran>RotorSpec::rank())
    throw Mismatch("SpinningHamiltonian::tensor");
  
  clear();

  const Tensor<double>& dvals=infop_->dvalues();
  for (int l=0;l<=ran;l++) {
    if (!A.have_rank(l))
      continue;

    ensure_tensor(1);
    for (int m=-l;m<=l;m++) {
      cmatrix& dest(component(m));
      const complex val(dvals(l,m)*A(l,m));
      if (!dest)
	dest.create(1,1,val);
      else
	dest(0U,0U)+=val;
    }
  }
}

static const double sqrtfac=-std::sqrt(0.5);

template<class M> void SpinningHamiltonian::add_ns(const space_T& A, const M& H)
{
  const int useL=getL(H.dimension(0U)); //validateL(A,H.dimension(0U));
  if (useL>2)
    throw InvalidParameter("add_ns");
  const cmatrix& d2s=infop_->d2values();

  if (tensor_values.empty())
    clear();

  //  if (useL) {
  assert(useL<=2);
  if (A.have_rank(2)) {
    for (int n=-useL;n<=useL;n++) {
      const complex anglefac=expi(-n*infop_->orientation());
      for (int m=-2;m<=2;m++) {
	const complex coeff=anglefac*d2s(m+2,n+2)*A(2,m);
#ifndef NDEBUG
	//      std::cout << "Adding " << coeff << " of H(" << n <<") to component " << m << "\n";
#endif
	mla(component(m),coeff,H(n+useL));
      }
    }
  }
  //}
  if (A.have_rank(0))
    mla(component(0),sqrtfac*real(A(0,0)),H(useL));
}

void SpinningHamiltonian::add(const space_T& A, const BaseList<cmatrix>& H)
{
  add_ns(A,H);
}

void SpinningHamiltonian::add(const space_T& A, const BaseList<rmatrix>& H)
{
  add_ns(A,H);
}

void SpinningHamiltonian::add(const space_T& A, const MultiMatrix<double,3>& H)
{
  add_ns(A,H);
}

std::ostream& operator<< (std::ostream& ostr, const SpinningHamiltonian &a)
{
  if (!a)
    return ostr << "<empty>\n";

  for (int m=-a.rank();m<=a.rank();m++) {
    const Matrix<complex>& am(a.component(m));
    switch (am.rows()) {
    case 0: break;
    case 1:
      ostr << "H_" << m << ": " << am(0,0) << '\n';
      break;
    default:
      ostr << "H_" << m << "\n" << am << '\n';      
    }
  }
  return ostr;
}

std::ostream& operator<< (std::ostream& ostr, const RealSpinningHamiltonian &a)
{
  if (!a)
    return ostr << "<empty>\n";
  
  const Matrix<double>& am0(a.component0());
  switch (am0.rows()) {
  case 0: break;
  case 1:
    ostr << "H_0: " << am0(0,0) << '\n';
    break;
  default:
    ostr << "H_0\n" << am0 << '\n';
  }
    
  for (int m=1;m<=a.rank();m++) {
    const Matrix<complex>& am(a.component(m));
    switch (am.rows()) {
    case 0: break;
    case 1:
      ostr << "H_" << m << ": " << am(0,0) << '\n';
      break;
    default:
      ostr << "H_" << m << "\n" << am << '\n';      
    }
  }
  return ostr;
}

void SpinningHamiltonian::operator() (cmatrix &Htot,double t) const
{
  if (rank()<0)
    throw Undefined("SpinningHamiltonian()");
  if (isdefined(component(0)))
    Htot=component(0);
  else {
    if (Htot.isdynamic())
      Htot.clear();
    else
      Htot=complex(0.0);
  }
  
  const double iphase(phase(t));
  for (int m=1;m<=rank();m++) {
    const complex fac=expi(-m*iphase);
    const cmatrix& Hm=component(m);
    if (isdefined(Hm)) {
      mla(Htot,fac,Hm);
      mla(Htot,conj(fac),component(-m));
    }
  }
}

template<class T> void RealSpinningHamiltonian::operator() (Matrix<T>& Htot,double t) const
{
  if (isdefined(component0()))
    Htot=component0();
  else {
    if (Htot.isdynamic())
      Htot.clear();
    else
      Htot=T(0.0);
  }

  if (!tensor_values.empty()) {
    const double iphase(phase(t));
    for (int m=1;m<=rank();m++) {
      const complex fac=2.0*expi(-m*iphase);
      const cmatrix& Hm=component(m);
      if (isdefined(Hm))
	real_mla(Htot,fac,Hm);
    }
  }
}

template void RealSpinningHamiltonian::operator() (Matrix<double>&, double) const;
template void RealSpinningHamiltonian::operator() (Matrix<complex>&, double) const;

void unitary_simtrans(SpinningHamiltonian& dest, const SpinningHamiltonian& source,const cmatrix& U)
{
  const int r=dest.rank();
  if (source.rank()!=r)
    throw Mismatch("unitary_simtrans");
  for (int m=-r;m<=r;m++) {
    if (isdefined(source.component(m)))
      unitary_simtrans(dest.component(m),source.component(m),U);
    else
      dest.component(m).clear();
  }
}

void unitary_simtrans(SpinningHamiltonian& dest, const RealSpinningHamiltonian& source,const cmatrix& U)
{
  const int r=dest.rank();
  if (source.rank()!=r)
    throw Mismatch("unitary_simtrans");
  if (isdefined(source.component0()))
    unitary_simtrans(dest.component(0),source.component0(),U);
  else
    dest.component(0).clear();

  cmatrix tmp;
  for (int m=1;m<=r;m++) {
    if (isdefined(source.component(m))) {
      tmp=source.component(m);
      unitary_simtrans(dest.component(m),tmp,U);
      tmp.conj();
      unitary_simtrans(dest.component(-m),tmp,U);
    }
    else {
      dest.component(m).clear();
      dest.component(-m).clear();
    }
  }
}

void unitary_simtrans(SpinningHamiltonian& dest, const DiagonalSpinningHamiltonian& source,const cmatrix& U)
{
  const int r=dest.rank();
  if (source.rank()!=r)
    throw Mismatch("unitary_simtrans");
  if (source.component0().length())
    unitary_simtrans(dest.component(0),source.component0(),U);
  else
    dest.component(0).clear();

  List<complex> tmp;
  for (int m=1;m<=r;m++) {
    tmp=source.component(m);
    unitary_simtrans(dest.component(m),tmp,U);
    conj_ip(tmp);
    unitary_simtrans(dest.component(-m),tmp,U);
  }
}

  SpinningHamiltonian& SpinningHamiltonian::operator+= (const SpinningHamiltonian& a) {
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("SpinningHamiltonian::add");

    for (size_t m=tensor_values.length();m--;) {
      const cmatrix& am=a.tensor_values(m);
      if (isdefined(am))
	tensor_values(m)+=am;
    }
    return *this;
  }

  SpinningHamiltonian& SpinningHamiltonian::operator-= (const SpinningHamiltonian& a) {
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("SpinningHamiltonian::add");
    for (size_t m=tensor_values.length();m--;) {
      const cmatrix& am=a.tensor_values(m);
      if (isdefined(am))
	tensor_values(m)-=am;
    }
    return *this;
  }

  SpinningHamiltonian& SpinningHamiltonian::operator+= (const DiagonalSpinningHamiltonian& a) {
    int m=a.rank();
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("SpinningHamiltonian::add");
    component(0)+=a.component0();
    for (;m>0;m--) {
      component(m)+=a.component(m);
      component(-m)+=a.component(m);
    }
    return *this;
  }

  SpinningHamiltonian& SpinningHamiltonian::operator-= (const DiagonalSpinningHamiltonian& a) {
    int m=a.rank();
    if ((m>RotorSpec::rank()) || (static_cast<const RotorSpec&>(a)!=*this))
      throw Mismatch("SpinningHamiltonian::add");
    component(0)-=a.component0();
    for (;m>0;m--) {
      component(m)-=a.component(m);
      component(-m)-=a.component(m);
    }
    return *this;
  }

  RealSpinningHamiltonian& RealSpinningHamiltonian::operator+= (const RealSpinningHamiltonian& a) {
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("RealSpinningHamiltonian::add");
    H0+=a.H0;
    tensor_values+=a.tensor_values;
    return *this;
  }

  RealSpinningHamiltonian& RealSpinningHamiltonian::operator-= (const RealSpinningHamiltonian& a) {
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("RealSpinningHamiltonian::add");
    H0-=a.H0;
    tensor_values-=a.tensor_values;
    return *this;
  }

  RealSpinningHamiltonian& RealSpinningHamiltonian::operator+= (const DiagonalSpinningHamiltonian& a) {
    int m=a.rank();
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("RealSpinningHamiltonian::add");
    H0+=a.component0();
    if (m) {
      ensure_tensor(a.component(1).size());
      for (;m>0;m--)
	component(m)+=a.component(m);
    }
    return *this;
  }

  RealSpinningHamiltonian& RealSpinningHamiltonian::operator-= (const DiagonalSpinningHamiltonian& a) {
    int m=a.rank();
    if (static_cast<const RotorSpec&>(a)!=*this)
      throw Mismatch("RealSpinningHamiltonian::add");
    H0-=a.component0();
    if (m) {
      ensure_tensor(a.component(1).size());
      for (;m>0;m--)
	component(m)-=a.component(m);
    }
    return *this;
  }

namespace {
  template<class T> double real_(const T& a) { return real(a); }
  template<> double real_(const double& a) { return a; }
  
  template<class T> double norm_(const T& a) { return norm(a); }
  template<> double norm_(const double& a) { return a*a; }
}

template<class T> void hermitian_eigenvalues1(List<double>& eigs,const Matrix<T>& H,const BaseList<double>& Hzeeman)
{
  if (!issquare(H))
    throw NotSquare("hermitian_eigenvalues1");
  const int n=H.rows();
  if (n!=Hzeeman.length())
    throw Mismatch("hermitian_eigenvalues1");
  eigs.create(n);

  for (size_t i=0;i<n;i++) {
    double val(real_(H(i,i)));
    for (size_t j=0;j<n;j++) {
      if (i!=j) {
 	const double diff=Hzeeman(i)-Hzeeman(j);
 	if (diff==0.0)
	  throw Failed("hermitian_eigenvalues1: degenerate eigenvalues");
 	val+=norm_(H(i,j))/diff;
      }
    }
    eigs(i)=val;
  }
}

#include "lcm_accumulate.hpp"

 template<class M> void HomogeneousPropagator<M>::operator() (cmatrix& U,double t1,double t2) const
   {       
     if (!Hamp)
       throw InternalError("HomogeneousPropagator()");

     if (period && (t2-t1+lcm_MAS_timingtol>2.0*period)) {
       const int reptimes=int((t2-t1+lcm_MAS_timingtol)/period);
       cmatrix Utmp2; //!< must use different temporary
       prop_U(Utmp2,t1,t1+period);
       if (verbose_)
	 std::cout << "HomogeneousPropagator: accumulating propagator " << reptimes << " times\n";
       pow(U,Utmp2,reptimes);
       if (verbose_>1)
	 std::cout << U;
       t1+=reptimes*period;
     }
     prop_U(U,t1,t2);
   }

inline void propagatorL_(cmatrix& Utmp, const cmatrix& H, double dt, const matrix_partition* partp)
{
  propagatorL(Utmp,H,dt,partp);
}

inline void propagatorL_(cmatrix&, const rmatrix&, double, const matrix_partition*)
{
  throw InternalError("propagatorL: Liouvillian must be complex");
}

IntervalSampler lcm_defsampler; //!< default sampler - can be shared as carries no state

Warning<> rotordrift_warning("timing jitter has resulted in significant drift between timing used for propagator integration and nominal timing. This may cause unphysical results.",&lcm_base_warning);

 template<class M> void HomogeneousPropagator<M>::prop_U(cmatrix& U,double t1,double t2) const
 {  
   if (fabs(t2-t1)<lcm_MAS_timingtol)
     return;

   if (t2<t1)
     throw InvalidParameter("HomogeneousPropagator: start time is after end time!");

   if (verbose_) {     
     std::cout << "HomogeneousPropagator: accumulating propagator from ";
     prettyprint_time(t1) << " to ";
     prettyprint_time(t2) << '\n';
   }
   cmatrix Usupertmp;
   bool needwarn=false;   
   bool haveflag=false;
   const bool usens=(Zeeman.length()!=0);
   if ((mode_==LiouvilleSpace) && usens)
     throw Failed("Can't (currently) combine Liouville space propagation and non-secular Hamiltonians");

   IntervalSamplerBase* usesamplerp=&lcm_defsampler;
   if (Hamp->samplerp() ) {
     if (!Samplerp)
       Samplerp.reset( (Hamp->samplerp() )->clone());       
     usesamplerp=Samplerp.get();
   }

   const double maxdttol=intdt+lcm_MAS_timingtol;
   double midt=-1e30; //!< keep midt so we can check for significant drift
   for (;;) {
     double dt=t2-t1;
     if (fabs(dt)<lcm_MAS_timingtol)
       break;
     if (dt>maxdttol)
       dt=intdt;
     midt=(*usesamplerp)(t1,t1+dt);

     (*Hamp)(H,midt);
     if (verbose_) {
       std::cout << "Sampling H at 'midpoint' of " << (1e6*t1) << " and " << (1e6*(t1+dt)) << ", t=" << (1e6*midt) << " us\n";
       if (verbose_>1)
	 std::cout << H;
     }
     if (H.empty())
       needwarn=true;
     else {
       switch (mode_) {
       case LiouvilleSpace:
	 propagatorL_(Utmp,H,dt,partp);
	 break;
       case Chebyshev:
	 chebyshev_propagator(Utmp,H,dt,partp);
	 break;
       case Diagonalisation:
	 if (usens)
	   propagator_ns(Utmp,H,dt,Zeeman);
	 else
	   propagator(Utmp,H,dt);
	 break;
       default:
	 throw InternalError("HomogeneousPropagator");
       }
       if (verbose_>1)
	 std::cout << "U at time t=" << (1e6*midt) << " us\n" << U;
       accumulate_(U,Utmp,Usupertmp,haveflag);
     }
     t1+=dt;
   }
   if (fabs(midt-t2)>intdt) {
     char scratch[256];
     snprintf(scratch,sizeof(scratch),"Expected end time: %g us  Last sampling time: %g us",t2*1e6,midt*1e6);
     rotordrift_warning.raise(scratch);
   }     
   if (needwarn)
     propagation_emptyhamiltonian_warning.raise();
   if (!haveflag) {
     if (U.isdynamic())
       U.clear();
     else
       U.identity();
   }
   if (verbose_>1)
     std::cout << "Result\n" << U;
 }
 
template void HomogeneousPropagator<rmatrix>::operator() (cmatrix& U,double t1,double t2) const;
template void HomogeneousPropagator<cmatrix>::operator() (cmatrix& U,double t1,double t2) const;

template<class T> void SpinningHamiltonian::_add(const space_T& A, const T& H, size_t sizev)
{
  const int userank=A.max_rank();
  if (userank==0) {
    add(real(A(0,0)),H);
    return;
  }
  if (userank>RotorSpec::rank())
    throw Mismatch("SpinningHamiltonian::add");

  const Tensor<double>& dvals=infop_->dvalues();
  (void)ensure_tensor(sizev);

  for (int m=-userank;m<=userank;m++) {
    complex tot(0.0);
    for (int l=::std::abs(m);l<=userank;l++) {
      if (A.have_rank(l))
	mla(tot,dvals(l,m),A(l,m));
    }
    //std::cout << "Adding " << tot << " of " << m << '\n';
    mlad(component(m),tot,H);
  }
}

template void SpinningHamiltonian::_add(const space_T&, const cmatrix&, size_t);
template void SpinningHamiltonian::_add(const space_T&, const rmatrix&, size_t);
template void SpinningHamiltonian::_add(const space_T&, const BaseList<double>&, size_t);

template<class T> void RealSpinningHamiltonian::_add(const space_T& A, const T& H, const size_t sizev)
{
  const int userank=A.max_rank();
  if (userank==0) {
    add(real(A(0,0)),H);
    return;
  }
  if (hasodd(A))
    throw InvalidParameter(NOODD);
  if (userank>RotorSpec::rank()) 
    throw Mismatch("RealSpinningHamiltonian::add");

  const Tensor<double> &dvals=infop_->dvalues();
  (void)ensure_tensor(sizev);

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
      cmatrix tmp(component(m));
      mlad(tmp,tot,H);
    }
    else
      mlad(H0,real(tot),H);
  }
}

template void RealSpinningHamiltonian::_add(const space_T& A, const rmatrix& H, size_t);
template void RealSpinningHamiltonian::_add(const space_T& A, const BaseList<double>& H, size_t);

#include "lcm_CommonNMR.cc"

cmatrix gammareduce(const BaseList<cmatrix>& source, size_t n)
{
  cmatrix U(mxflag::temporary);
  gammareduce_(U,source,source.size(),n);
  return U;
}

void gammareduce(BaseList<cmatrix> dest, const BaseList<cmatrix>& source, size_t n)
{
  gammareduce_(dest,source,n);
}

IntervalSampler::IntervalSampler(double jitterv, double cortimv, bool randomisev)
  : IntervalSamplerBase(randomisev),
    needsreset_(true)
  { 
    jitter(jitterv,cortimv);
  }

  double IntervalSampler::operator() (double t1,double t2)
  {
    double midt=0.5*(t1+t2);
    if (jitter_==0.0)
      return midt;

    const double interval=t2-t1;
    if (needsreset_) {
      lastresettime_=midt;
      reset();
      lastjitter_=0.0;
      needsreset_=false;
    }
    else {
      if (midt==lastmidt_)
	return midt+lastjitter_;           

      if (midt<lastmidt_) {
	// re-winding
	if (fabs(midt-lastresettime_)>1e-10) {
	  char scratch[256];
	  snprintf(scratch,sizeof(scratch)," %g vs %g",midt,lastresettime_);
	  invalidreset_warning.raise(scratch);
	  //!< can't do anything sensible here
	}
	else {// reset
	  reset();
	  lastjitter_=0.0;
	}
      }
    }	  
    lastmidt_=midt;
    lastjitter_+=jitter_*ggen_(interval);
    midt+=lastjitter_;
//     if (midt<t1) {
//       excessivejitter_warning.raise();
//       return t1;
//     }
//     if (midt>t2) {
//       excessivejitter_warning.raise();
//       return t2;
//     }
    return midt;
  }

void IntervalSampler::jitter(double jitterv, double cortimv)
  {
    if (jitterv<0.0)
      throw InvalidParameter("IntervalSampler: jitter must be >=0");
    jitter_=jitterv;
    ggen_.set_correlationtime(cortimv);
  }
  
//  Warning<> IntervalSampler::excessivejitter_warning("Jitter is excessive compared to interval required",&lcm_base_warning,BaseWarning::FirstOnly);
Warning<> IntervalSampler::invalidreset_warning("Re-winding IntervalSampler to time which doesn't match initial reset time",&lcm_base_warning);

IntervalSamplerBase::synchronise_t IntervalSamplerBase::get_synchronisation() const
{
  return randomise_ ?
    RandomGenerator::get_seed_random()
    : RandomGenerator::get_seed_default();
}
  
  void IntervalSamplerBase::reset()
  {
    set_synchronisation(IntervalSamplerBase::get_synchronisation());
  }
  
void IntervalSampler::set_synchronisation(synchronise_t seedval)
{
#ifndef NDEBUG
    std::cout << "Resetting RNG with seed " << seedval << '\n';
#endif
  ggen_.set_seed(seedval);
}

}//namespace libcmatrix
