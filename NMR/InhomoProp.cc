#include "NMR.h"
#include "MAS.h"
#include "Propagation.h"
#include "ScratchList.h"
#include "HomoProp_shared.h"
#include "InhomoProp_shared.h"

namespace libcmatrix {

  template<class T> void unitary_isimtrans_In(cmatrix& dest,const RCmatrix& R,const Matrix<T>& a,const RCmatrix& C,Matrix<T>* tmp)
{
  if (R.type()==RCmatrix::NONE) {
    switch (C.type()) {
    case RCmatrix::NONE:
      dest=a;
      return;
    case RCmatrix::REAL:
      multiply(dest,a,C.get_real());
      return;
    case RCmatrix::COMPLEX:
      multiply(dest,a,C.get_complex());
      return;
    }
  }
  else { 
    Matrix<T> Ttmp_;
    Matrix<T>& Ttmp = tmp ? *tmp : Ttmp_;

    if (R.type()==RCmatrix::REAL) {
      switch (C.type()) {
      case RCmatrix::NONE:
	transpose_multiply(dest,R.get_real(),a);
	return;
      case RCmatrix::REAL:
	unitary_isimtransLR(dest,R.get_real(),a,C.get_real(),&Ttmp);
	return;
      case RCmatrix::COMPLEX:
	unitary_isimtransLR(dest,R.get_real(),a,C.get_complex(),&Ttmp);
	return;
      }
    }
    else {
      if (R.type()==RCmatrix::COMPLEX) {
	cmatrix ctmp;

	switch (C.type()) {
	case RCmatrix::NONE:
	  conj_transpose_multiply(dest,R.get_complex(),a);
	  return;
	case RCmatrix::REAL:
	  unitary_isimtransLR(dest,R.get_complex(),a,C.get_real(),&Ttmp);
	  return;
	case RCmatrix::COMPLEX:
	  unitary_isimtransLR(dest,R.get_complex(),a,C.get_complex(),&ctmp);
	  return;
	}
      }
    }
  }
  throw Failed("unitary_isimtrans: transformation undefined");
}


template void unitary_isimtrans_In(cmatrix&, const RCmatrix&, const Matrix<complex>&, const RCmatrix&, Matrix<complex>*);
template void unitary_isimtrans_In(cmatrix&, const RCmatrix&, const Matrix<double>&, const RCmatrix&, Matrix<double>*);

template<class T> void unitary_isimtrans_allownull(Matrix<T>& dest,const rmatrix& R,const Matrix<T>& a,const rmatrix& C,Matrix<T>* tmp)
{
  if (!R) {
    if (!C) 
      dest=a;
    else
      multiply(dest,a,C);
  }
  else { 
    if (!C)
      transpose_multiply(dest,R,a);
    else {
      Matrix<T> Ttmp_;
      unitary_isimtransLR(dest,R,a,C,tmp ? tmp : &Ttmp_);
    }
  }
}

template void unitary_isimtrans_allownull(Matrix<double>&, const rmatrix&, const Matrix<double>&, const rmatrix&, Matrix<double>*);
template void unitary_isimtrans_allownull(Matrix<complex>&, const rmatrix&, const Matrix<complex>&, const rmatrix&, Matrix<complex>*);

void unitary_isimtrans_identity(cmatrix &UT,const RCmatrix& R,const RCmatrix& C)
{
  switch (R.type()) {
  case RCmatrix::REAL:
    switch (C.type()) {
    case RCmatrix::REAL:
      transpose_multiply(UT,R.get_real(),C.get_real());
      return;
    case RCmatrix::COMPLEX:
      transpose_multiply(UT,R.get_real(),C.get_complex());
      return;
    }
    break;
  case RCmatrix::COMPLEX:
    switch (C.type()) {
    case RCmatrix::REAL:
      conj_transpose_multiply(UT,R.get_complex(),C.get_real());
      return;
    case RCmatrix::COMPLEX:
      conj_transpose_multiply(UT,R.get_complex(),C.get_complex());
      return;
    }
    break;
  }
  throw Failed("unitary_isimtrans_identity");
}

void unitary_isimtrans_In(cmatrix& dest,const BaseList<double>& a,const RCmatrix& stash)
{
  switch (stash.type()) {
  case RCmatrix::NONE:
    full(dest,a);
    return;
  case RCmatrix::REAL: {
    rmatrix tmpr;
    transpose_multiply(tmpr,stash.get_real(),a);
    multiply(dest,tmpr,stash.get_real());
  }
    return;
  case RCmatrix::COMPLEX:
    unitary_isimtrans(dest,a,stash.get_complex());
    return;
  }
  throw Failed("unitary_isimtrans: transformation undefined");
}

void InhomoStash::set_H(const rmatrix& H)
{ 
  invalidatephases();
  const size_t dim=H.rows();
  switch (dim) {
  case 1:
    if (!issquare(H))
      throw NotSquare("InhomoStash::set_H");
    set_H(H(0,0));
    break;
  case 0:
    clear();
    eff_freqs.clear();
    phasestore.clear();
    break;
  default:
    eff_freqs.create(1,dim);
    BaseList<double> asr=eff_freqs.row();
    hermitian_eigensystem(set_real(),asr,H);    
    break;
  }
}

void InhomoStash::set_H(const rmatrix& V, const BaseList<double>& eigs)
{ 
  invalidatephases();
  if (V.rows()!=eigs.size())
    throw Mismatch("InhomoStash::set_H");
  eff_freqs.create(1,eigs.size(),eigs.vector());
  set_real()=V;
}

void InhomoStash::set_H(const cmatrix& V, const BaseList<double>& eigs)
{ 
  invalidatephases();
  if (V.rows()!=eigs.size())
    throw Mismatch("InhomoStash::set_H");
  eff_freqs.create(1,eigs.size(),eigs.vector());
  set_complex()=V;
}

void InhomoStash::set_H(const cmatrix& H)
{ 
  invalidatephases();
  const size_t dim=H.rows();
  switch (dim) {
  case 1:
    if (!issquare(H)) throw NotSquare("InhomoStash::set_H");
    set_H(real(H(0,0)));
    break;
  case 0:
    clear();
    eff_freqs.clear();
    phasestore.clear();
    break;
  default:
    eff_freqs.create(1,dim);
    BaseList<double> asr=eff_freqs.row();
    hermitian_eigensystem(set_complex(),asr,H);
    break;
  }
}

void InhomoStash::set_H(const BaseList<double>& H)
{
  invalidatephases();
  clear();
  if (H.empty()) {
    eff_freqs.clear();
    phasestore.clear();
  }
  else {
    //    std::cout << "set_H6: " << H << " " << H.size() << std::endl;
    eff_freqs.create(1,H.size(),H.vector());
  }
  //  std::cout << "set_H7" << std::endl;
}

void InhomoStash::set_H(double H)
{
  invalidatephases();
  clear();
  eff_freqs.create(1,1,H);
}

void InhomoStash::set_U(const cmatrix& U, double period)
{
  const size_t dim=U.rows();
  if (dim==1) {
    if (!issquare(U))
      throw NotSquare("InhomoStash::set_U");
    set_U(U(0,0),period);
    return;
  }
  invalidatephases();
  eff_freqs.create(1,dim);
  BaseList<double> asr=eff_freqs.row();
  diag_propagator(set_complex(),asr,U,period);
}
  
void InhomoStash::set_U(const complex& U, double period)
{
  invalidatephases();
  clear();
  eff_freqs.create(1,1,diag_propagator(U,period));
}

void InhomoStash::set_U(const BaseList<complex>& Uv, double periodv)
{
  invalidatephases();
  clear();
  eff_freqs.create(1,Uv.size());
  BaseList<double> asr=eff_freqs.row();
  diag_propagator(asr,Uv,periodv);
}

bool InhomoIter::advance()
{
  if (++r==dimR) {
    if (offdiagonly) {
      s++;
      r=s+1;
      if (r==dimR) {
	finished_=true;
	return false;
      }
    }
    else {
      r=0;
      if (++s==dimC) {
	finished_=true;
	return false;
      }
    }
  }
  return true;
}

std::ostream& operator<< (std::ostream& ostr,const InhomoIter& a)
{
  if (a.finished_)
    return ostr << "Position: finished\n";
  else
    return ostr << "Position: " << a.r << "," << a.s << '\n';
}

std::ostream& operator<< (std::ostream& ostr,const InhomoStash& a)
{
  if (!a)
    return ostr << "<unset>\n";
  if (!(a.eff_freqs.empty()))
    ostr << "Frequencies:\n" << a.eff_freqs << '\n';
  ostr << "Transformation matrix:\n" << static_cast<const RCmatrix&>(a);
  ostr << "Propagators:";
  if (a.isvalid)
    ostr << '\n' << a.phasestore << '\n';
  else
    ostr << " <not valid>\n";
  return ostr;
}

std::ostream& operator<< (std::ostream& ostr,const StashH& a)
{
  if (!a)
    return ostr << "<unset>\n";
  ostr << "Frequencies: " << a.eff_freqs << '\n';
  ostr << "Transformation matrix:\n" << static_cast<const RCmatrix&>(a);
  return ostr;
}

 std::ostream& operator<< (std::ostream& ostr,const BaseMASInhomoSpectrum& a)
{ 
  ostr << static_cast<const MASInhomoObj& >(a) << a.iter_;
  if (!a.finished)
    ostr << "Sideband position: " << a.ncount << '\n';
  return ostr;
}

  void InhomoIter::reset(size_t dimRv,size_t dimCv,bool offdiagv)
  {
    dimR=dimRv; dimC=dimCv; offdiagonly=offdiagv;
    s=0;
    r=offdiagonly ? 1 : 0;
    finished_=false;
  }

  void InhomoStash::updatephases(const rmatrix& phases,bool istimedom)
  {
    if (isvalid)
      return;
    if (eff_freqs.empty()) {
      if (phasestore.empty())
	return; //assume phasestore set explicitly
      throw Failed("InhomoStash: can't updatephases if set_H not used");
    }
    const size_t nints=this->interactions();
    const size_t neigs=eff_freqs.cols();
    if (nints!=phases.rows())
      throw Mismatch("updatephases");
    const size_t nobs=phases.cols();
    phasestore.create(neigs,nobs);
    for (size_t eig=neigs;eig--;) {
      for (size_t obs=nobs;obs--;) {
	double s=0.0;
	for (size_t intn=nints;intn--;)
	  s+=eff_freqs(intn,eig)*phases(intn,obs);
	if (obs || istimedom)	  
	  phasestore(eig,obs)=expi(s);
	else
	  phasestore(eig,0)=s;
      }
    }
#ifndef NDEBUG
    std::cout << "Propagators:\n" << phasestore;
#endif
    isvalid=true;
  }


void MASInhomoObj::sideband(int val)
{
  if (tdomain)
    throw Failed("Cannot restrict sideband for time-domain propagation");  
  restrictsideband=true;
  valsideband=val;
  valrestrict=-1; //!< invalid value
}

  void MASInhomoObj::updatephases() const
    //combine dynamic phases with Hamiltonian eigenvalues
  {
    if (!*this)
      throw Failed("MASInhomoObj: Hamiltonian(s)/phases not set");
    MASInhomoObj* mutable_this= const_cast< MASInhomoObj* >(this);
    if (!(row().isvalid))
      mutable_this->row().updatephases(phase,tdomain);
    if (!(col().isvalid))
      mutable_this->col().updatephases(phase,tdomain);
  }

  complex MASInhomoObj::calcps(size_t r,size_t s) const
  {
    const BaseList<complex> Rr=row().phasestore.row(r);
    const BaseList<complex> Cs=col().phasestore.row(s);

    sideband_.create(this->gammasteps());
    sideband_(0)=1.0;
    for (size_t gn=this->gammasteps();--gn;)
      //      sideband_(gn)=multiply_conj(Rr(gn),Cs(gn));
      sideband_(gn)=conj_multiply(Cs(gn),Rr(gn));

    //    return multiply_conj(Rr(0),Cs(0));
    return conj_multiply(Cs(0),Rr(0));
  }

void GammaInhomogeneousFID::add_FID_(BaseList<complex> FIDv,complex scale) const
{
  if (!*this)
    throw Failed("GammaInhomogeneousFID: Hamiltonian(s)/phases not set");

  scale/=this->gammasteps();
  const size_t npoints=FIDv.size();
  if (npoints==0)
    throw Undefined("add_FID: zero-length FID!");

  const bool isdiag=isdiagonal();

  if (isdiag) {
    const complex intzero=scale*trace(A);
    if (intzero!=0.0) FIDv+=intzero;
  }

  const size_t nobs_=this->observations();
  const size_t dimR=rows();
  const size_t dimC=cols();

  if (isdiag) {
    for (size_t lr=1;lr<dimR;lr++) {
      for (size_t ls=lr;ls--;) {
	complex amprs=scale*A(lr,ls);
	complex ampsr=scale*A(ls,lr);
	const bool dors=(norm(amprs)>1e-16);
	const bool dosr=(norm(ampsr)>1e-16);
	if (dosr || dors) {
	  const complex v(calcps(lr,ls));
	  const complex vc(conj(v));

	  for (size_t j=0;j<nobs_;j++) {
	    const complex p=getp(j,v);

	    if (dors) {
	      complex ploc=p*amprs;
	      if (verbose_>1)
		std::cout << lr << "," << ls << "," << j << ": adding " << ploc << " of phase " << v << '\n';
	      for (size_t k=j;k<npoints;k+=nobs_) {
		FIDv(k)+=ploc;
		ploc*=v;
	      }
	    }
	    if (dosr) {
	      complex ploc=conj_multiply(p,ampsr);
	      if (verbose_>1) 
		std::cout << ls << ',' << lr << ',' << j << ": adding " << ploc << " of phase " << v << '\n';
	      for (size_t k=j;k<npoints;k+=nobs_) {
		FIDv(k)+=ploc;
		ploc*=vc;
	      }
	    }
	  }
	}
      }
    }
  }
  else {
    for (size_t lr=dimR;lr--;) {
      for (size_t ls=dimC;ls--;) {
      
	complex amp=A(lr,ls);
	if (norm(amp)>1e-16) {
	  amp*=scale;

	  const complex v=calcps(lr,ls);
	  
	  for (size_t j=0;j<nobs_;j++) {
	    complex p=amp*getp(j,v);
	    if (verbose_>1)
	      std::cout << lr << ',' << ls << ',' << j << ": adding " << p << " of phase " << v << '\n';
	    for (size_t k=j;k<npoints;k+=nobs_) {
	      FIDv(k)+=p;
	      p*=v;
	    }
	  }
	}
      }
    }
  }
}

void InhomogeneousFID::add_FID_(BaseList<complex> FIDv,complex scale) const
{
  if (!*this)
    throw Failed("InhomogeneousFID: Hamiltonian(s)/phases not set");

  const size_t npoints=FIDv.size();
  if (npoints==0)
    throw Undefined("add_FID: zero-length FID!");

  const bool isdiag=isdiagonal();
  if (isdiag) {
    const complex intzero=scale*trace(A);
    if (intzero!=0.0)
      FIDv+=intzero;
  }

  const cmatrix& C=col().phasestore;

  const size_t dimR=rows();
  const size_t dimC=cols();
  const size_t nobs_=this->observations();

  for (size_t r=dimR;r--;) {
    const BaseList<complex> Rr=row().phasestore.row(r);

    for (size_t s=dimC;s--;) {
      if (isdiag && (r==s))
	continue;

      const complex amp=scale*A(r,s);
      if (norm(amp)>1e-16) {
	const BaseList<complex> Cs=C.row(s);
	const complex v=multiply_conj(Rr(0),Cs(0));
	for (size_t j=0;j<nobs_;j++) {
	  complex p= j ? amp*multiply_conj(Rr(j),Cs(j)) : amp;
	  if (verbose_>1)
	    std::cout << r << ',' << s << ',' << j << ": adding " << p << " of phase " << v << '\n';
	  for (size_t k=j;k<npoints;k+=nobs_) {
	    FIDv(k)+=p;
	    p*=v;
	  }
	}
      }
    }
  }
}

void MASInhomoObj::set_Hs(const rmatrix& H)
{
  invalidatephases();
  InhomoStash& stash=setdiagonal();
  stash.clear();
  stash.eff_freqs=H;
}

void MASInhomoObj::set_Hs(char sel,const rmatrix& H)
{
  invalidatephases();
  InhomoStash& which=setRC(sel);
  which.clear();
  which.eff_freqs=H;
}

double BaseMASInhomoSpectrum::frequency(size_t lr,size_t ls) const
{
  updatephases();
  return real(col().phasestore(ls,0))-real(row().phasestore(lr,0));
}

void InhomogeneousFID::add_FID_hermitian_(BaseList<double> FIDv,double scale) const
{
  if (!*this)
    throw Failed("InhomogeneousFID: Hamiltonian(s)/phases not set");

  if (!isdiagonal()) 
    throw Failed("add_FID_hermitian: diagonal blocks only");

  const size_t npoints=FIDv.size();
  if (npoints==0)
    throw Undefined("add_FID_hermitian: zero-length FID!");

  const double intzero=real_trace(A);
  if (intzero)
    FIDv+=scale*intzero;
  scale*=2.0;

  const cmatrix& F=col().phasestore;

  const size_t dim=rows();
  const size_t nobs_=this->observations();

  for (size_t r=1;r<dim;r++) {
    for (size_t s=r;s--;) {
      const complex amp=scale*A(r,s);
      if (norm(amp)>1e-16) {
	const complex v=multiply_conj(F(r,0),F(s,0));
	
	for (size_t j=0;j<nobs_;j++) {
	  complex p= j ? amp*multiply_conj(F(r,j),F(s,j)) : amp;

	  if (verbose_>1)
	    std::cout << r << "," << s << "," << j << ": adding " << p << " of phase " << v << '\n';
	  
	  for (size_t k=j;k<npoints;k+=nobs_) {
	    FIDv(k)+=real(p);
	    p*=v;
	  }
	}
      }
    }
  }
}

void InhomogeneousFID::add_FID_hermitian_(BaseList<complex> FIDv,complex scale) const
{
  if (!*this)
    throw Failed("InhomogeneousFID: Hamiltonian(s)/phases not set");

  if (!isdiagonal()) 
    throw Failed("add_FID_hermitian: diagonal blocks only");

  const size_t npoints=FIDv.size();
  if (npoints==0)
    throw Undefined("add_FID_hermitian: zero-length FID!");

  const double intzero=real_trace(A);
  if (intzero)
    FIDv+=scale*intzero;
  scale*=2.0;

  const cmatrix& F=col().phasestore;

  const size_t dim=rows();
  const size_t nobs_=this->observations();

  for (size_t r=1;r<dim;r++) {
    for (size_t s=r;s--;) {
      const complex amp=A(r,s);
      if (norm(amp)>1e-16) {
	const complex v=multiply_conj(F(r,0),F(s,0));
	
	for (size_t j=0;j<nobs_;j++) {
	  complex p= j ? amp*multiply_conj(F(r,j),F(s,j)) : amp;

	  if (verbose_>1)
	    std::cout << r << "," << s << "," << j << ": adding " << p << " of phase " << v << '\n';
	  
	  for (size_t k=j;k<npoints;k+=nobs_) {
	    mla(FIDv(k),real(p),scale);
	    p*=v;
	  }
	}
      }
    }
  }
}

void MASInhomoObj::invalidatephases()
{ 
  if (explicit_Us())
    throw Failed("Can't set H/phases if propagators are being set explicitly");

  row().invalidatephases();
  col().invalidatephases();
}

complex MASInhomoObj::getp(size_t j,const complex& v) const
{
  //  const size_t gamma_steps=this->gammasteps();
  size_t off=j*gamma_steps/this->observations();
  complex p(0.0);
  size_t gn;
  const complex* vecbuf=sideband_.vector();
  const complex* vecoff=vecbuf+off;
  for (gn=gamma_steps-off;gn--;)
    mla_conj(p,vecoff[gn],vecbuf[gn]);
  vecoff-=gamma_steps;
  complex p2(0.0);
  for (gn=gamma_steps-off;gn<gamma_steps;gn++)
    mla_conj(p2,vecoff[gn],vecbuf[gn]);
  mla(p,v,p2);

  return p;
}

void GammaInhomogeneousFID::add_FID_hermitian_(BaseList<complex> FIDv, complex scale) const
{
  if (!*this)
    throw Failed("GammaInhomogeneousFID: Hamiltonian(s)/phases not set");

  if (!isdiagonal())
    throw Failed("add_FID_hermitian: diagonal blocks only");
  scale/=this->gammasteps();

  const size_t npoints=FIDv.size();
  if (npoints==0)
    throw Undefined("add_FID_hermitian: zero-length FID!");

  const double intzero=real_trace(A);
  if (intzero)
    FIDv+=scale*intzero;
  scale*=2.0;

  const size_t dim=rows();
  const size_t nobs_=this->observations();

  for (size_t r=1;r<dim;r++) {
    for (size_t s=r;s--;) {
      complex amp=A(r,s);
      if (norm(amp)>1e-16) {

	const complex v=calcps(r,s);
	for (size_t j=0;j<nobs_;j++) {
	  complex p=amp*getp(j,v);
	  if (verbose_>1)
	    std::cout << r << "," << s <<  "," << j << ": adding " << p << " of phase " << v << '\n';
	  for (size_t k=j;k<npoints;k+=nobs_) {
	    mla(FIDv(k),scale,real(p));
	    p*=v;
	  }
	}
      }
    }
  }
}

void GammaInhomogeneousFID::add_FID_hermitian_(BaseList<double> FIDv, double scale) const
{
  if (!*this)
    throw Failed("GammaInhomogeneousFID: Hamiltonian(s)/phases not set");

  if (!isdiagonal()) 
    throw Failed("add_FID_hermitian: diagonal blocks only");
  scale/=this->gammasteps();

  const size_t npoints=FIDv.size();
  if (npoints==0)
    throw Undefined("add_FID_hermitian: zero-length FID!");

  const double intzero=real_trace(A);
  if (intzero)
    FIDv+=scale*intzero;
  scale*=2.0;

  const size_t dim=rows();
  const size_t nobs_=this->observations();

  for (size_t r=1;r<dim;r++) {
    for (size_t s=r;s--;) {
      complex amp=scale*A(r,s);
      if (norm(amp)>1e-16) {

	const complex v=calcps(r,s);
	for (size_t j=0;j<nobs_;j++) {
	  complex p=amp*getp(j,v);
	  if (verbose_>1)
	    std::cout << r << ',' << s <<  ',' << j << ": adding " << p << " of phase " << v << '\n';
	  for (size_t k=j;k<npoints;k+=nobs_) {
	    FIDv(k)+=real(p);
	    p*=v;
	  }
	}
      }
    }
  }
}

void
BaseMASInhomoSpectrum::reset_()
{ 
  updatephases();
  if (!*this)
    lock();
  else {
    iter_.reset(rows(),cols(),isdiagonal()); 
    havenext=false; ncount=0;
    finished=iter_.isfinished();
  }
  setclean();
}

template<class T> void
BaseSingleInhomo<T>::reset()
{
  BaseMASInhomoSpectrum::reset_();
  scale_factor=1.0/this->observations();
  if (isdiagonal() && ::libcmatrix::isdefined(A) ) {
    const T zeroint=trace(A); // will catch non-square S
    if (zeroint!=0.0) {
      havenext=true;
      nextamp=zeroint;
      nextfreq=0.0;
    }
  }
}

template<class T> void
BaseGammaInhomo<T>::reset()
{
  BaseMASInhomoSpectrum::reset_();
  //const size_t gamma_steps=this->gammasteps();
  scale_factor=1.0/(gamma_steps*gamma_steps);
  if (isdiagonal() && ::libcmatrix::isdefined(A) ) {
    const T zeroint=trace(A); // will catch non-square S
    if (zeroint!=0.0) {
      havenext=true;
      nextamp=zeroint;
      nextfreq=0.0;
    }
  }
}

void MASInhomoObj::set_steps(int stepsv) 
{ 
  if (PeriodicObj::set_steps(stepsv))
    invalidatephases();
  if (restrictsideband)
    throw InternalError("set_steps"); //!< method should be over-ridden
}

bool PeriodicObj::set_steps(int stepsv)
{
  if (stepsv<1)
    throw InvalidParameter("PeriodicObj: number of steps must be non-zero");

  size_t& steps= isgamma ? gamma_steps : nobs;

  if (steps==stepsv)
    return false;

  //  if (steps) throw Failed("PeriodicObj: can't change number of steps");  
  
  if (isgamma) {
    if (stepsv % nobs)
	throw Failed("Gamma steps must be integer multiple of observations");
    //    intsteps=stepsv/nobs;
  }
  else
    gamma_steps=1;
  
  steps=stepsv;
  return true;
}

void BaseMASInhomoSpectrum::set_steps(int stepsv)
{
  if (PeriodicObj::set_steps(stepsv)) {
    usefft=ispowerof2(stepsv);
    if (!usefft) {
      ft_facs.create(stepsv);
      ft_create_table(ft_facs,FT_FORWARD);
      scr.create(stepsv);
    }
  }
  if (restrictsideband) {
    const int nobs2=nobs/2;
    
    if ((valsideband>nobs2) || (valsideband<-nobs2))
      throw InvalidParameter("BaseMASInhomoSpectrum: sideband out of range");
    valrestrict= (valsideband<0) ? valsideband+nobs : valsideband;
  }
}

Warning<> MASInhomoObj::ignoring_period_warning("set_Us: period specified for time domain propagation is ignored",&lcm_base_warning);

void MASInhomoObj::set_Us(const cmatrix& Us, double lperiod)
{
  if (Us.empty())
    throw Undefined("set_Us");
  if (tdomain) {
    if (lperiod)
      ignoring_period_warning.raise();
  }
  else {
    if (!lperiod)
      throw Failed("set_Us: must specify period for frequency domain propagation");
    this->period(lperiod);
  }
  set_steps(Us.rows());
  setdiagonal().set_Us(Us,lperiod);
}

void InhomoStash::set_Us(const cmatrix& Us, double period)
{
  clear();
  const size_t dim=Us.cols();
  const size_t steps=Us.rows();
  phasestore.create(dim,steps);
  for (size_t k=dim;k--;) {
    BaseList<complex> storek(phasestore.row(k));
    if (period==0.0) {
      size_t j=steps-1;
      storek(0U)=Us(j,k);
      for (;j--;) 
	storek(j+1)=Us(j,k);
    }
    else {  
      const double freq=diag_propagator(Us(steps-1,k),period);
      storek(0U)=freq;
      const complex fac=propagator(-freq,period/steps);
      complex accfac(fac);
      for (size_t j=1;j<steps;j++) {
	storek(j)=Us(j-1,k)*accfac;
	accfac*=fac;
      }
    }
  }
  eff_freqs.clear();
  isvalid=true;
}

void MASInhomoObj::set_Us(char sel,const cmatrix& Us, double lperiod)
{
  if (Us.empty())
    throw Undefined("set_Us");
  set_steps(Us.rows());
  setRC(sel).set_Us(Us,lperiod);
}

void MASInhomoObj::set_phases(const BaseList<DynamicPhase>& cycle)
{
  const size_t ints=cycle.size();
  const size_t steps=cycle(0).observations();
  if (steps==0)
    throw Failed("MASInhomoObj: number of sampling steps unset");

  invalidatephases();

  phase.create(ints,steps);
  double speed=0.0;

  for (size_t k=0;k<ints;k++) {
    const DynamicPhase& ccycle=cycle(k);
    if (ccycle.observations()!=steps) throw Mismatch("MASInhomObj: number of phase steps is not consistent");
    if (tdomain) {
      phase(k,0)=ccycle.isotropic(steps);
      for (size_t n=1;n<steps;n++) phase(k,n)=ccycle(n);
    }
    else {
      phase(k,0)=ccycle.component0();
      for (size_t n=1;n<steps;n++) phase(k,n)=ccycle.anisotropic(n);
    }
    if (k==0)
      speed=fabs(ccycle.rotor_speed());
    else {
      if (speed!=fabs(ccycle.rotor_speed()))
	throw Mismatch("phases: rotor speeds must be equal!");
    }
    this->period(1.0/fabs(speed));
  }
  set_steps(steps);
}

void MASInhomoObj::set_phases(const BaseList<double>& phasesv,double periodv)
{
  const size_t steps=phasesv.size();
  if (steps<1)
    throw InvalidParameter("MASInhomoObj: phases vector has zero size!");
  invalidatephases();

  phase.create(1,steps);
  BaseList<double> phaserow=phase.row();

  if (tdomain)
    phaserow=phasesv;
  else {
    this->period(periodv);
    const double finphase=phasesv(steps-1);
    
    phaserow(0)=finphase/(-2*M_PI*periodv);
    for (size_t n=1;n<steps;n++)
      phaserow(n)=phasesv(n-1)-(n*finphase)/steps;
  }
  set_steps(steps);
}

void MASInhomoObj::set_phases(const rmatrix& phasesv,double periodv)
{
  const size_t steps=phasesv.cols();
  const size_t ints=phasesv.rows();

  invalidatephases();

  if (tdomain)
    phase=phasesv;
  else {
    this->period(periodv);
    phase.create(ints,steps);
    
    for (size_t k=ints;k--;) {
      const double finphase=phasesv(k,steps-1);
      phase(k,0)=finphase/(-2*M_PI*periodv);
      for (size_t n=1;n<steps;n++)
	phase(k,n)=phasesv(k,n-1)-(n*finphase)/steps;
    }
  }
  set_steps(steps);
}

std::ostream& operator<< (std::ostream &ostr,const MASInhomoObj& a)
{
  ostr << static_cast<const SwapStore<InhomoStash>& >(a);
  if (a.observations()==0)
    return ostr << "<unset>\n";
  ostr << "Dynamic phases: " << a.phase << '\n';
  return ostr;
}

template<typename T> std::ostream& operator<< (std::ostream &ostr,const BaseSingleInhomo<T>& a)
{
  if (a.observations()==0)
    return ostr << "<unset>\n";
  ostr << static_cast< const BaseMASInhomoSpectrum& >(a);
  ostr << "Transition amplitudes:\n" << a.A << '\n';
  return ostr;
}

template<typename T> std::ostream& operator<< (std::ostream &ostr,const BaseGammaInhomo<T>& a)
{
  if (a.observations()==0)
    return ostr << "<unset>\n";
  ostr << static_cast< const BaseMASInhomoSpectrum& >(a);
  ostr << "Gamma integration steps: " << a.gammasteps() << '\n';
  ostr << "Transition amplitudes:\n" << a.A << '\n';
  return ostr;
}

void BaseMASInhomoSpectrum::calcamps_(size_t lr,size_t ls)
{
  viso=frequency(lr,ls);

  const BaseList<complex> R=row().phasestore.row(lr);
  const BaseList<complex> C=col().phasestore.row(ls);

  List<complex>& store=usefft ? sideband_ : scr;
  store.create(R.size());
  store(0)=1.0;

  for (size_t j=R.size()-1;j;j--)
    store(j)=multiply_conj(R(j),C(j));

  if (usefft)
    fft_ip(sideband_,FT_FORWARD);
  else
    ft(sideband_,store,ft_facs);
}

template<typename T> void BaseGammaInhomo<T>::calcamps(size_t lr,size_t ls)
{
  calcamps_(lr,ls);
  
  amps.create(nobs);
  for (size_t n=nobs;n--;) {
    double totamp=0.0;
    for (size_t k=n;k<gamma_steps;k+=nobs)
      totamp+=norm(sideband_(k));
    amps(n)=totamp;
  }
}

template<typename T> bool BaseSingleInhomo<T>::operator() (complex& amp,double &freq)
{
  if (finished)
    return false;
  verifyclean();

  if (havenext) {
    amp=nextamp;
    freq=nextfreq;
    havenext=false;
    return true;
  }

  const bool isdiag=isdiagonal();

  if (ncount==0) {
    if (iter_.isfinished()) {
      finished=true;
      return false;
    }

    if (isdiag) {
      while ( (A(iter_.r,iter_.s)==0.0) && (A(iter_.s,iter_.r)==0.0)) { // skip over null transition probabilities
	if (!iter_.advance())
	  return false;
      }
      ampfac=A(iter_.r,iter_.s);
    }
    else {
      while ( (ampfac=A(iter_.r,iter_.s))==0.0) {
	if (!iter_.advance())
	  return false;
      }
    }

    calcamps(iter_.r,iter_.s);
    ampfac*=scale_factor;
    if (isdiag)
      nextampfac=A(iter_.s,iter_.r)*scale_factor;

    iter_.advance();
  }

  size_t usencount=ncount;
  if (restrictsideband)
    usencount=valrestrict;
  else {
    if (++ncount==this->observations())
      ncount=0;
  }

  const complex& totamp=sideband_(usencount);
  amp=ampfac*totamp;
  const double negfreq=(this->observations()-usencount)*this->speed()-viso;
  
  if (isdiag) {
    const double posfreq=usencount*this->speed()+viso;
    nextfreq=posfreq;
    nextamp=conj_multiply(totamp,nextampfac);

    if (ampfac!=0.0) {
      freq=negfreq;
      havenext=true;
    }
    else {
      amp=nextamp;
      freq=nextfreq;
    }
  }
  else
    freq=negfreq;
     
  return true;
}

template<typename T> bool BaseGammaInhomo<T>::operator() (T& amp,double &freq)
{
  if (finished)
    return false;
  verifyclean();

  if (havenext) {
    amp=nextamp;
    freq=nextfreq;
    havenext=false;
    return true;
  }

  const bool isdiag=isdiagonal();
  
  if (ncount==0) {
    if (iter_.isfinished()) {
      finished=true;
      return false;
    }

    if (isdiag) {
      while ( (A(iter_.r,iter_.s)==0.0) && (A(iter_.s,iter_.r)==0.0)) { // skip over null transition probabilities
	if (!iter_.advance())
	  return false;
      }
      ampfac=A(iter_.r,iter_.s);
    }
    else {
      while ( (ampfac=A(iter_.r,iter_.s))==0.0) {
	if (!iter_.advance())
	  return false;
      }
    }

    calcamps(iter_.r,iter_.s);
    ampfac*=scale_factor;
    if (isdiag)
      nextampfac=A(iter_.s,iter_.r)*scale_factor;

    iter_.advance();
  }//ncount==0

  size_t usencount=ncount;
  if (restrictsideband)
    usencount=valrestrict;
  else {
    if (++ncount==this->observations())
      ncount=0;
  }

  amp=ampfac*amps(usencount);

  const double negfreq=(this->observations()-usencount)*this->speed()-viso;  
  if (isdiag) {
    const double posfreq=usencount*this->speed()+viso;
    nextfreq=posfreq;
    nextamp=nextampfac*amps(usencount);

    if (ampfac!=0.0) {
      freq=negfreq;
      havenext=true;
    }
    else {
      amp=nextamp;
      freq=nextfreq;
    }
  }
  else
    freq=negfreq;
     
  return true;
}

template<class T> void BaseSingleInhomo<T>::add(List<complex>& damps,size_t lr,size_t ls)
{
  if (damps.size()==0)
    damps.create(this->observations(),complex(0.0));
  add( static_cast< BaseList<complex>& >(damps),lr,ls);
}

template<class T> void BaseGammaInhomo<T>::add(List<T>& damps,size_t lr,size_t ls)
{
  if (damps.size()==0)
    damps.create(this->observations(),T(0.0));
  add( static_cast< BaseList<T>& >(damps),lr,ls);
}

template<typename T> void BaseGammaInhomo<T>::add(BaseList<T>& ampout,size_t lr,size_t ls)
{
  if (lr>=iter_.dimR || ls>=iter_.dimC)
    throw BadIndex("InhomogeneousSpectrum::add");
  if (ampout.size()!=this->observations())
    throw Mismatch("InhomogeneousSpectrum::add");

  iter_.lock(); //calcamps invalidates iterator
  const T amp(A(lr,ls));
  if (amp!=0.0) {
    calcamps(lr,ls);
    mla(ampout,amp*scale_factor,amps);
  }
}

template<typename T> void BaseSingleInhomo<T>::add(BaseList<complex>& ampout,size_t lr,size_t ls)
{
  if (lr>=iter_.dimR || ls>=iter_.dimC)
    throw BadIndex("InhomogeneousSpectrum::add");
  if (ampout.size()!=this->observations())
    throw Mismatch("InhomogeneousSpectrum::add");

  iter_.lock(); //calcamps invalidates iterator
  const T amp(A(lr,ls));
  if (amp!=0.0) {
    calcamps(lr,ls);
    mla(ampout,amp*scale_factor,sideband_);
  }
}

void InhomogeneousSpectrum::observe(const cmatrix& sigma0,const cmatrix& det) {
  InhomoHelper_::observe(*this,A,sigma0,det);
}
void InhomogeneousSpectrum::observe(const rmatrix& sigma0,const rmatrix& det) {
  InhomoHelper_::observe(*this,A,sigma0,det);
}
void InhomogeneousSpectrum::observe(const BaseList<double>& sigma0,const BaseList<double>& det) {
  InhomoHelper_::observe(*this,A,sigma0,det);
}

void InhomogeneousFID::observe() { 
  InhomoHelper_::observe(*this,A);
} 
void InhomogeneousFID::observe(const rmatrix& sigma0det) {
  InhomoHelper_::observe(*this,A,sigma0det);
} 
void InhomogeneousFID::observe(const cmatrix& sigma0det) { 
  InhomoHelper_::observe(*this,A,sigma0det);
} 
void InhomogeneousFID::observe(const BaseList<double>& sigma0det) {
  InhomoHelper_::observe(*this,A,sigma0det);
} 
void InhomogeneousFID::observe(const rmatrix& sigma0, const rmatrix& detect) {
  InhomoHelper_::observe(*this,A,sigma0,detect);
} 
void InhomogeneousFID::observe(const cmatrix& sigma0, const cmatrix& detect) {
  InhomoHelper_::observe(*this,A,sigma0,detect);
} 
void InhomogeneousFID::observe(const BaseList<double>& sigma0, const BaseList<double>& detect) {
  InhomoHelper_::observe(*this,A,sigma0,detect);
} 

void InhomogeneousSpectrumED::observe(const cmatrix& sigma0det)
{ 
  InhomoHelper_::observe(*this,A,sigma0det);
}
void InhomogeneousSpectrumED::observe(const rmatrix& sigma0det) {
  InhomoHelper_::observe(*this,A,sigma0det);
}
void InhomogeneousSpectrumED::observe(const BaseList<double>& sigma0det) {
  InhomoHelper_::observe(*this,A,sigma0det); 
}
void InhomogeneousSpectrumED::observe() { 
  InhomoHelper_::observe(*this,A);
}

void GammaInhomogeneousSpectrumED::observe(const cmatrix& sigma0det) { 
  InhomoHelper_::observe(*this,A,sigma0det);
}
void GammaInhomogeneousSpectrumED::observe(const rmatrix& sigma0det) {
  InhomoHelper_::observe(*this,A,sigma0det); 
}
void GammaInhomogeneousSpectrumED::observe(const BaseList<double>& sigma0det) {
  InhomoHelper_::observe(*this,A,sigma0det);
}
void GammaInhomogeneousSpectrumED::observe() { 
  InhomoHelper_::observe(*this,A);
}

void GammaInhomogeneousSpectrum::observe(const cmatrix& sigma0, const cmatrix& det) { 
  InhomoHelper_::observe(*this,A,sigma0,det);
}
void GammaInhomogeneousSpectrum::observe(const rmatrix& sigma0, const rmatrix& det) { 
  InhomoHelper_::observe(*this,A,sigma0,det);
}
void GammaInhomogeneousSpectrum::observe(const BaseList<double>& sigma0,const BaseList<double>& det) {
  InhomoHelper_::observe(*this,A,sigma0,det); 
}

void GammaInhomogeneousFID::observe() { 
  InhomoHelper_::observe(*this,A);
} 
void GammaInhomogeneousFID::observe(const rmatrix& sigma0det) { 
  InhomoHelper_::observe(*this,A,sigma0det);
} 
void GammaInhomogeneousFID::observe(const cmatrix& sigma0det) { 
  InhomoHelper_::observe(*this,A,sigma0det);
} 
void GammaInhomogeneousFID::observe(const BaseList<double>& sigma0det) {
  InhomoHelper_::observe(*this,A,sigma0det); 
} 
void GammaInhomogeneousFID::observe(const rmatrix& sigma0, const rmatrix& detect) { 
  InhomoHelper_::observe(*this,A,sigma0,detect); 
} 
void GammaInhomogeneousFID::observe(const cmatrix& sigma0, const cmatrix& detect) { 
  InhomoHelper_::observe(*this,A,sigma0,detect);
} 
void GammaInhomogeneousFID::observe(const BaseList<double>& sigma0, const BaseList<double>& detect) { 
  InhomoHelper_::observe(*this,A,sigma0,detect);
} 

std::ostream& operator<< (std::ostream& ostr,const GammaInhomogeneousFID& a) {
  return ostr << static_cast<const MASInhomoObj&>(a) << "Transition matrix\n" << a.A;
}

std::ostream& operator<< (std::ostream& ostr,const InhomogeneousFID& a) {
  return ostr << static_cast<const MASInhomoObj&>(a) << "Transition matrix\n" << a.A;
}
   
//template class SwapStore<InhomoStash>;

template class BaseSingleInhomo<double>;
template class BaseSingleInhomo<complex>;
template class BaseGammaInhomo<complex>;
template class BaseGammaInhomo<double>;

template std::ostream& operator<< (std::ostream &ostr,const BaseSingleInhomo<complex>&);
template std::ostream& operator<< (std::ostream &ostr,const BaseSingleInhomo<double>&);
template std::ostream& operator<< (std::ostream &ostr,const BaseGammaInhomo<double>&);
template std::ostream& operator<< (std::ostream &ostr,const BaseGammaInhomo<complex>&);

}//namespace libcmatrix

// void InhomoStash::write(FILE *fp,const char *name) const
// {
//   ScratchList<char> tmpstr_loc(strlen(name)+12);
//   char *tmpstr=tmpstr_loc.vector();
//   sprintf(tmpstr,"%sFreqs",name);
//   if (WriteMATLAB(fp,eff_freq,tmpstr)) throw Failed("WriteMATLAB");
//   sprintf(tmpstr,"%sV",name);
//   switch (type) {
//   case NONE:
//     {
//       rmatrix dummy;// can't assume that DR is undefined
//       if (!WriteMATLAB(fp,dummy,tmpstr)) return;
//     }
//     break;
//   case REAL:
//     if (!WriteMATLAB(fp,DR,tmpstr)) return;
//     break;
//   case COMPLEX:
//     if (!WriteMATLAB(fp,DC,tmpstr)) return;
//   default: break;
//   }
//   throw Failed("WriteMATLAB");
// }

// int InhomoStash::read(FILE *fp)
// {
//   int fail,rowcol;
//   if ((fail=ReadMATLAB(eff_freq,fp))) return fail;
//   bool iscomplex;
//   if ((fail=ReadMATLABHeader(rowcol,rowcol,iscomplex,fp))) return fail;
//   if (rowcol==0) {
//     type=NONE;
//     return 0;
//   }
//   if (iscomplex) {
//     type=COMPLEX;
//     return ReadMATLAB(DC,fp);
//   }
//   else {
//     type=REAL;
//     return ReadMATLAB(DR,fp);
//   }
// }
