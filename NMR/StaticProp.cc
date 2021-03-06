#include "BasePropagation.h"
#include "lcm_StaticProp.h"
#include "HomoProp_shared.h"
#include "InhomoProp_shared.h"
#include "NMR.h"

namespace libcmatrix {

inline void doscale(cmatrix& sigma,complex scale)
{
  if (imag(scale))
    sigma*=scale;
  else {
    double rscale=real(scale);
    if (rscale!=1.0)
      sigma*=rscale;
  }
}

  //Not terribly efficient, but only used by StaticFID_U, which is not terribly efficient anyway
void unitary_simtrans_gen_ip(cmatrix& d,const cmatrix& VR,const cmatrix& VC,cmatrix& tmp)
{
  if (VR.rows()==1) {
    if (d.rows()!=1)
      throw Mismatch("unitary_simtrans_gen_ip");
    if (VC.rows()==1) {
      if (d.cols()!=1)
	throw Mismatch("unitary_simtrans_gen_ip");
      unitary_simtransLR(d(0,0),VR(0,0),d(0,0),VC(0,0));
    }
    else {
      tmp.create(1,VC.rows());
      unitary_simtransLR(tmp.row(),VR(0,0),d.row(),VC);
      d=tmp;
    }
  }
  else {
    if (VC.rows()==1) {
      if (d.cols()!=1)
	throw Mismatch("unitary_simtrans_gen_ip");
      tmp.create(VR.rows(),1);
      unitary_simtransLR(tmp.row(),VR,d.row(),VC(0,0));
      d=tmp;
    }
    else
      unitary_simtransLR(d,VR,d,VC,&tmp);
  }
}

void 
StaticSpectrum::reset()
{ 
  if (!*this)
    iter_.lock();
  else {
    iter_.reset(rows(),cols(),false);
    eff_freqsR.create(row().eff_freqs.row());
    eff_freqsC.create(col().eff_freqs.row());
    setclean();
  }
}

void
StaticSpectrumED::reset()
{
  if (!*this)
    iter_.lock();
  else {
    iter_.reset(rows(),cols(),false);
    eff_freqsR.create(row().eff_freqs.row());
    eff_freqsC.create(col().eff_freqs.row());
    setclean();
  }
}

template<class T> bool StaticBase_::get_trans(T& amp,double& freq,const Matrix<T>& A)
{
  if (iter_.isfinished())
    return false;
  verifyclean();
  const bool nominmax(min_f==max_f);
  for (;;) {
    while ( (amp=A(iter_.r,iter_.s))==0.0) {
      if (!iter_.advance())
	return false;
    }
    freq=eff_freqsR(iter_.r)-eff_freqsC(iter_.s);
    const bool isOK(iter_.advance());
    if (nominmax)
      return true;
    if ((freq>min_f) && (freq<max_f)) {
      freq-=larmor_f;
      return true;
    }
    if (!isOK)
      return false;
  }
}

bool StaticSpectrum::operator() (complex& amp,double& freq)
{
  return get_trans(amp,freq,AC);
}

bool StaticSpectrumED::operator() (double& amp,double& freq)
{
  return get_trans(amp,freq,AR);
}

void StaticSpectrum::print(std::ostream& ostr) const
{
  ostr << static_cast<const SwapStore<InhomoStash>& >(*this);
  ostr << iter_;
  ostr << "Transition amplitudes:\n" << AC << std::endl;
}

void StaticSpectrumED::print(std::ostream& ostr) const
{
  ostr << static_cast<const SwapStore<InhomoStash>& >(*this);
  ostr << iter_;
  ostr << "Transition amplitudes:\n" << AR << std::endl;
}

static const char FIDU_NOPROPS[]="StaticFID_U: propagators not set";

void StaticFID_U::add_FID(BaseList<complex> FIDv,complex scale,const cmatrix& sigma0,const cmatrix& det)
{
  if (!*this)
    throw Failed(FIDU_NOPROPS);
  if (!aretranspose(sigma0,det))
    throw Mismatch("add_FID: sigma0 and detect operator incompatible");

  sigma=sigma0;
  doscale(sigma,scale);
  const size_t n=FIDv.size();
  if (n==0)
    throw Undefined("add_FID: zero-length FID!");
  const cmatrix& VR=row();
  const cmatrix& VC=col();

  for (size_t i=0;i<n;i++) {
    FIDv(i)+=trace_multiply(det,sigma);
    unitary_simtrans_gen_ip(sigma,VR,VC,tmp);
  }
}

void StaticFID_U::add_FID(BaseList<complex> FIDv,complex scale)
{
  if (!*this)
    throw Failed(FIDU_NOPROPS);
  if (isdiagonal())
    throw Failed("add_FID: only valid for off-diagonal propagation");

  const size_t n=FIDv.size();
  if (n==0)
    throw Undefined("add_FID: zero-length FID!");

  const cmatrix& VR=row();
  const cmatrix& VC=col();

  FIDv(0)+=scale*float_t(VR.rows());
  multiply_conj_transpose(sigma,VR,VC);
  doscale(sigma,scale);

  for (size_t i=1;i<n;i++) {
    FIDv(i)+=trace(sigma);
    if (i<n-1)
      unitary_simtransLR(sigma,VR,sigma,VC,&tmp);
  }
}

void StaticFID_H::larmor(double larmor_)
{
  larmor(larmor_,larmor_*0.5,larmor_*1.5);
}

void StaticBase_::larmor(double larmor_)
{
  larmor(larmor_,larmor_*0.5,larmor_*1.5);
}

void StaticFID_H::larmor(double larmor_,double min_, double max_)
{
  if (min_==max_)
    throw InvalidParameter("larmor: width cannot be zero");
    
  if (min_>max_) 
    ::std::swap(min_,max_);

  if ((larmor_<=min_) || (larmor_>=max_))
    throw InvalidParameter("larmor: larmor frequency outside range");

  larmor_f=larmor_;
  max_f=max_;
  min_f=min_;
}
 
void StaticBase_::larmor(double larmor_,double min_, double max_)
{
  if (min_==max_)
    throw InvalidParameter("larmor: width cannot be zero");
    
  if (min_>max_) 
    ::std::swap(min_,max_);

  if ((larmor_<min_f) || (larmor_>max_f))
    throw InvalidParameter("larmor: larmor frequency outside range");

  larmor_f=larmor_;
  max_f=max_;
  min_f=min_;
}
 
template<class T> void StaticFID_U::add_FID_hermitian(BaseList<T> FIDv,T scale,const cmatrix& sigma0,const cmatrix& det)
{
  if (!*this)
    throw Failed(FIDU_NOPROPS);
  const bool forceherm=!isdiagonal();
  const size_t n=FIDv.size();
  if (n==0)
    throw Undefined("add_FID_hermitian: zero-length FID!");
  sigma=sigma0;
  const cmatrix& UR(row());
  if (forceherm) {
    const cmatrix& UC(col());
    for (size_t i=0;i<n;i++) {
      mla(FIDv(i),scale,real(trace_multiply(det,sigma)));
      unitary_simtrans_gen_ip(sigma,UR,UC,tmp);
    }
  }
  else {
    for (size_t i=0;i<n;i++) {
      mla(FIDv(i),scale,hermitian_trace_multiply(det,sigma));
      unitary_simtrans_gen_ip(sigma,UR,UR,tmp);
    }
  }
}

template void StaticFID_U::add_FID_hermitian(BaseList<double>, double, const cmatrix&, const cmatrix&);
template void StaticFID_U::add_FID_hermitian(BaseList<complex>, complex, const cmatrix&, const cmatrix&);

void StashH::set_H(const cmatrix& H,double dt)
{
  const size_t dim=H.rows();
  if (dim==1) {
    if (!issquare(H))
      throw NotSquare("StashH::set_H");
    set_H(real(H(0,0)),dt);
    return;
  }
  eff_freqs.create(dim);
  hermitian_eigensystem(set_complex(),eff_freqs,H);
  finish_set_(dt);
}

void StashH::set_H(const rmatrix& V,const BaseList<double>& eigs,double dt)
{ 
  if (V.rows()!=eigs.size())
    throw Mismatch("StashH::set_H");
  eff_freqs=eigs;
  set_real()=V;
  finish_set_(dt);
}

void StashH::set_H(const cmatrix& V,const BaseList<double>& eigs,double dt)
{ 
  if (V.rows()!=eigs.size())
    throw Mismatch("StashH::set_H");
  eff_freqs=eigs;
  set_complex()=V;
  finish_set_(dt);
}

void StashH::set_U(const cmatrix& U)
{
  evol_facs.create(U.rows());
  eigensystem(set_complex(),evol_facs,U);
  dtstore=0.0;
  eff_freqs.clear();
}

void StashH::set_H(const rmatrix& H,double dt)
{ 
  const size_t dim=H.rows();
  if (dim==1) {
    if (!issquare(H))
      throw NotSquare("StashH::set_H");
    set_H(H(0,0),dt);
    return;
  }
  eff_freqs.create(dim);
  hermitian_eigensystem(set_real(),eff_freqs,H);
  finish_set_(dt);
}

void StashH::finish_set_(double dt)
{
  evol_facs.create(eff_freqs.size());
  propagator(evol_facs,eff_freqs,dt);
  dtstore=dt;
}

void StashH::set_H(const BaseList<double>& H,double dt)
{
  clear();
  eff_freqs=H;
  finish_set_(dt);
}

void StashH::set_H(double H,double dt)
{
  clear();
  eff_freqs.create(1);
  eff_freqs=H;
  finish_set_(dt);
}

double StaticFID_H::get_dwell() const
{
  double dt=row().dtstore;
  if (dt!=col().dtstore)
    throw Mismatch("StaticFID_H: bra and ket spaces have different dwell times!");
  if (dt==0.0)
    throw Failed("StaticFID_H: dwell time not set (set_U rather than set_H?)");
  return dt;
}

void StaticFID_H::add_FID_(BaseList<complex> FIDv,complex scale) const
{
  const bool checkminmax=(min_f!=max_f);
  const double dt=checkminmax ? get_dwell() : 0.0;
  const BaseList<double>& ER(row().eff_freqs);
  const BaseList<double>& es(col().eff_freqs); //!< Nasty symbol clash on ES with some systems
  const BaseList<complex>& PR(row().evol_facs);
  const BaseList<complex>& PS(col().evol_facs);

  if (verbose_)
    std::cout << *this;
  
  const size_t n=FIDv.size();
  complex v;
  if (n==0)
    throw Undefined("add_FID: zero-length FID!");
  const size_t Acols=A.cols();
  for (size_t r=A.rows();r--;) {
    for (size_t s=Acols;s--;) {
      complex amp=A(r,s);
      if (norm(amp)>1e-16) {
	amp*=scale;
	if (checkminmax) {
	  double freq=ER(r)-es(s);
	  if ((freq<=min_f) || (freq>=max_f)) {
	    if (verbose_>1)
	      std::cout << "Rejecting " << amp << " of " << freq << '\n';
	    continue; //frequency outside detect width
	  }
	  freq-=larmor_f;
	  v=propagator(freq,dt); //we don't use pre-calculated values as these are subject to high numerical error
	  if (verbose_>1)
	    std::cout << "Adding " << amp << " of " << freq << '\n';
	}
	else
	  v=multiply_conj(PR(r),PS(s));

	for (size_t i=0;i<n;i++) {
	  FIDv(i)+=amp;
	  amp*=v;
	}
      }
    }
  }
}

template<class T> void StaticFID_H::add_FID_hermitian__(BaseList<T> FIDv, T scale) const
{
  const bool forceherm=!isdiagonal();

  const BaseList<complex>& PR(row().evol_facs);
  const BaseList<complex>& PC(col().evol_facs);
  const size_t n=FIDv.size();
  if (n==0)
    throw Undefined("add_FID_hermitian: zero-length FID!");
  if (!forceherm) {
    FIDv+=scale*hermitian_trace(A);
    scale*=2.0;
  }
  double initsum=0.0;
  const size_t startr= forceherm ? 0 : 1;
  for (size_t r=startr;r<A.rows();r++) {
    const size_t start=forceherm ? A.cols() : r;
    for (size_t s=start;s--;) {
      complex amp=A(r,s);
      if (norm(amp)>1e-16) {
	amp*=scale;
	const complex v=multiply_conj(PR(r),PC(s));
	initsum+=real(amp);
	for (size_t i=1;i<n;i++) {
	  amp*=v;
	  FIDv(i)+=real(amp);
	}
      }
    }
  }
  FIDv.front()+=initsum;
}

template void StaticFID_H::add_FID_hermitian__(BaseList<double>,double) const;
template void StaticFID_H::add_FID_hermitian__(BaseList<complex>,complex) const;

void StaticFID_H::observe() 
{ InhomoHelper_::observe(*this,A); } 

void StaticFID_H::observe(const rmatrix& sigma0det)
{ InhomoHelper_::observe(*this,A,sigma0det); } 

void StaticFID_H::observe(const cmatrix& sigma0det)
{ InhomoHelper_::observe(*this,A,sigma0det); } 

void StaticFID_H::observe(const BaseList<double>& sigma0det)
{ InhomoHelper_::observe(*this,A,sigma0det); } 

void StaticFID_H::observe(const rmatrix& sigma0, const rmatrix& detect)
{ InhomoHelper_::observe(*this,A,sigma0,detect); } 

void StaticFID_H::observe(const cmatrix& sigma0, const cmatrix& detect)
{ InhomoHelper_::observe(*this,A,sigma0,detect); } 

void StaticFID_H::observe(const BaseList<double>& sigma0, const BaseList<double>& detect)
{ InhomoHelper_::observe(*this,A,sigma0,detect); } 

void StaticSpectrum::observe(const cmatrix& sigma0, const cmatrix& det)
{ InhomoHelper_::observe(*this,AC,sigma0,det); }

void StaticSpectrum::observe(const rmatrix& sigma0, const rmatrix& det)
{ InhomoHelper_::observe(*this,AC,sigma0,det); }

void StaticSpectrum::observe(const BaseList<double>& sigma0,const BaseList<double>& det) 
{ InhomoHelper_::observe(*this,AC,sigma0,det); }
  
void StaticSpectrumED::observe(const cmatrix& sigma0det)
{ InhomoHelper_::observe(*this,AR,sigma0det); }

void StaticSpectrumED::observe(const rmatrix& sigma0det)
{ InhomoHelper_::observe(*this,AR,sigma0det); }

void StaticSpectrumED::observe(const BaseList<double> &sigma0det)
{ InhomoHelper_::observe(*this,AR,sigma0det); }

void StaticSpectrumED::observe()
{ InhomoHelper_::observe(*this,AR); }

} //namespace libcmatrix
