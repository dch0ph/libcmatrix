#include "NMR.h"
#include "HomoProp_shared.h"
#include "lcm_HomoProp.h"

namespace libcmatrix {

  static const double zerotol=1e-12; //!< ignore transitions with *norm* below this value

  Warning<InvalidParameter> propagation_closetodiagonal_warning("close-to-diagonal propagator detected.  System is inhomogeneous?",&lcm_base_warning);

static void mla(cmatrix& a,const cmatrix& b,const cmatrix& c)
{
  if (isdefined(a)) {
    BaseList<complex> arow=a.row();
    if (b.rows()!=c.rows())
      throw Mismatch("mla"); // only need check one since length is checked later
    mla(arow,b.row(),c.row());
  }
  else
    emultiply(a,b,c);
}

void doadd_PeriodicFID(BaseList<complex>& FIDv,complex sfac,const cmatrix& detect_tab,const complex& UcycleL,const BaseList<complex>& sigmav,const cmatrix& UcycleR)
{
  ScratchList<complex> sigma(sigmav);
  ScratchList<complex> tmp(sigma.size());
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_PeriodicFID: zero-length FID!");
  const size_t nprops=detect_tab.rows();

  size_t which=0;
  for (size_t k=0;k<totpoints;k++) {
    FIDv(k)+=sfac*trace_multiply(detect_tab.row(which),sigma);
    if (++which==nprops) {
      which=0;
      multiply_conj_transpose(tmp,sigma,UcycleR);
      sigma=tmp;
      sfac*=UcycleL;
    }
  }
}

void doadd_PeriodicFID(BaseList<complex>& FIDv,complex sfac,const cmatrix& detect_tab,const cmatrix& UcycleL,const BaseList<complex>& sigmav,const complex& UcycleR)
{
  ScratchList<complex> sigma(sigmav);
  ScratchList<complex> tmp(sigma.size());
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_PeriodicFID: zero-length FID!");
  const size_t nprops=detect_tab.rows();

  size_t which=0;
  for (size_t k=0;k<totpoints;k++) {
    FIDv(k)+=sfac*trace_multiply(detect_tab.row(which),sigma);
    if (++which==nprops) {
      which=0;
      multiply(tmp,UcycleL,sigma);
      sigma=tmp;
      sfac*=UcycleR;
    }
  }
}

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

void doadd_PeriodicFID(BaseList<complex>& FIDv,complex scale,const BaseList<complex>& detect_tab,const complex& UcycleL,complex sigma,const complex& UcycleR)
{
  sigma*=scale;
  //  const complex evolfac=multiply_conj(UcycleL,UcycleR);
  const complex evolfac=conj_multiply(UcycleR,UcycleL);
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_PeriodicFID: zero-length FID!");
  const size_t nprops=detect_tab.size();

  size_t which=0;
  for (size_t k=0;k<totpoints;k++) {
    FIDv(k)+=detect_tab(which)*sigma;
    if (++which==nprops) {
      which=0;
      sigma*=evolfac;
    }
  }
}

void doadd_PeriodicFID(BaseList<complex>& FIDv,complex scale,const MultiMatrix<complex,3>& detect_tab,const cmatrix& UcycleL,const cmatrix& sigmaref, const cmatrix& UcycleR,cmatrix& tmp)
{
  cmatrix sigma(sigmaref,mxflag::normal);
  doscale(sigma,scale);
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_PeriodicFID: zero-length FID!");
  const size_t nprops=detect_tab.dimension(0);

  size_t which=0;
  for (size_t k=0;k<totpoints;k++) {
    FIDv(k)+=trace_multiply(detect_tab(which),sigma);
    if (++which==nprops) {
      which=0;
      multiply_conj_transpose(tmp,sigma,UcycleR);
      multiply(sigma,UcycleL,tmp);
    }
  }
}

void doadd_PeriodicFID(BaseList<complex>& FIDv,complex scale,const MultiMatrix<complex,3>& detect_tab,const cmatrix& UcycleL,const cmatrix& UcycleR, cmatrix& tmp)
{
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_PeriodicFID: zero-length FID!");
  const size_t nprops=detect_tab.dimension(0);
  const size_t todo= (totpoints<nprops) ? totpoints : nprops;

  FIDv(0)+=scale*double(UcycleL.rows());
  size_t k=1;
  for (;k<todo;k++)
    FIDv(k)+=scale*trace(detect_tab(k));
  if (k==totpoints) return;
  
  size_t which=0;
  cmatrix sigma;
  multiply_conj_transpose(sigma,UcycleL,UcycleR);
  doscale(sigma,scale);

  for (;k<totpoints;k++) {
    FIDv(k)+=(which ? trace_multiply(detect_tab(which),sigma) : trace(sigma));
    if (++which==nprops) {
      which=0;
      multiply_conj_transpose(tmp,sigma,UcycleR);
      multiply(sigma,UcycleL,tmp);
    }
  }
}

template<class T> inline void fillfull(cmatrix& d,const Matrix<T>& a) { d=a; }
template<class T> inline void fillfull(cmatrix& d,const BaseList<T>& a) { full(d,a); }

template<class T,size_t N> struct getdim_ {
  static inline size_t func(const T& a) { return a.dimension(N); }
};
template<> struct getdim_< BaseList<cmatrix>, 1> {
  static inline size_t func(const BaseList<cmatrix>& a) { return a(0).rows(); }
};
template<> struct getdim_< BaseList<cmatrix>, 2> {
  static inline size_t func(const BaseList<cmatrix>& a) { return a(0).cols(); }
};

template<class TR,class TC> void table_create(MultiMatrix<complex,3>& table,const TR& Rs,const TC& Cs) {
  const size_t n=Rs.dimension(0);
  if (n!=Cs.dimension(0))
    throw Mismatch("table_create");  
  table.create(n,getdim_<TR,1>::func(Rs),getdim_<TC,2>::func(Cs));
}

void table_create(cmatrix& table,const BaseList<cmatrix>& Rs,const BaseList<complex>& Cs) {
  const size_t n=Rs.size();
  if (n!=Cs.size())
    throw Mismatch("table_create");  
  table.create(n,Rs(0).rows());
}

void table_create(cmatrix& table,const MultiMatrix<complex,3>& Rs,const BaseList<complex>& Cs) {
  const size_t n=Rs.dimension(0);
  if (n!=Cs.size())
    throw Mismatch("table_create");  
  table.create(n,Rs.dimension(1));
}

void table_create(cmatrix& table,const BaseList<complex>& Rs,const BaseList<cmatrix>& Cs) {
  const size_t n=Cs.size();
  if (n!=Rs.size())
    throw Mismatch("table_create");  
  table.create(n,Cs(0).cols());
}

void table_create(cmatrix& table,const BaseList<complex>& Rs,const BaseList<complex>& Cs) {
  const size_t n=Cs.size();
  if (n!=Rs.size())
    throw Mismatch("table_create");  
  table.create(n,1);
}

void table_create(cmatrix& table,const BaseList<complex>& Rs,const MultiMatrix<complex,3>& Cs) {
  const size_t n=Cs.dimension(0);
  if (n!=Rs.size())
    throw Mismatch("table_create");  
  table.create(n,Cs.dimension(2));
}

struct transform_TD_ {

  template<typename T,class M> inline void operator() (MultiMatrix<complex,3>& table,const M& UTRs,const M& UTCs,const T& a) const
  {
    table_create(table,UTRs,UTCs);
    cmatrix tmp;

    for (size_t j=UTRs.dimension(0);j--;) {
      cmatrix tablej(table(j));
      if (j)
	unitary_isimtransLR(tablej,UTRs(j-1),a,UTCs(j-1),&tmp);
      else
	fillfull(tablej,a);
    }
  }
    
  template<typename TR,typename TC> inline void operator() (cmatrix& table,const TR& UTRs,const TC& UTCs,const BaseList<complex>& a) const
  {
    table_create(table,UTRs,UTCs);

    for (size_t j=UTRs.dimension(0);j--;) {
      if (j)
	unitary_isimtransLR(table.row(j),UTRs(j-1),a,UTCs(j-1));
      else
	table.row(0)=a;
    }
  }
  
  template<class M> inline void operator() (MultiMatrix<complex,3>& table,const M& UTRs,const M& UTCs) const
{
   const size_t nprops=UTCs.dimension(0);  
   table_create(table,UTRs,UTCs);
   //don't require table(0)
   for (size_t j=1;j<nprops;j++) {
     cmatrix tablej(table(j));
     unitary_isimtrans_identity(tablej,UTRs(j-1),UTCs(j-1));
   }

}

  template<typename T> inline void operator() (List<complex>& table,const BaseList<complex>& UTRs,const BaseList<complex>& UTCs,T a) const
 {
   const size_t nprops=UTCs.size();
   if (nprops!=UTRs.size())
     throw Mismatch("transform_TD");
   table.create(nprops);
   table(0)=a;
   for (size_t j=1;j<nprops;j++)
     table(j)=a*conj_multiply(UTRs(j-1),UTCs(j-1));
 }
};

  template<class T> void doadd_PeriodicFID_hermitian(BaseList<T>& FIDv, T scale, const MultiMatrix<complex,3>& detect_tab, const cmatrix& UcycleL, const cmatrix& sigmaref, const cmatrix& UcycleR, cmatrix& tmp, bool forceherm =false)
{
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_PeriodicFID_hermitian: zero-length FID!");
  const size_t nprops=detect_tab.dimension(0);
  if (forceherm)
    scale*=2.0;	      
  cmatrix sigma(sigmaref,mxflag::normal);
  doscale(sigma,complex(scale));

  size_t which=0;
  for (size_t k=0;k<totpoints;k++) {
    if (forceherm)
      FIDv(k)+=real(trace_multiply(detect_tab(which),sigma)); //real(trace_multiply) != hermitian_trace_multiply!
    else
      FIDv(k)+=hermitian_trace_multiply(detect_tab(which),sigma);
    if (++which==nprops) {
      which=0;
      unitary_simtransLR(sigma,UcycleL,sigma,UcycleR,&tmp);
    }
  }
}

template<class T> void transpose_all(MultiMatrix<T,3>& a,const MultiMatrix<T,3>& b)
{
  size_t n=b.dimension(0);
  a.create(n,b.dimension(2),b.dimension(1));
  for (;n--;) {
    Matrix<T> an(a(n));
    transpose(an,b(n));
  }
}

void make_S(MultiMatrix<complex,3>& Sstore,const MultiMatrix<complex,3>& detect_table,const MultiMatrix<complex,3>& sigma_table,const cmatrix& R,size_t nobs)
{
  const size_t gamma_steps=sigma_table.dimension(0);
  Sstore.create(nobs,R.rows(),R.cols());
  MultiMatrix<complex,3> tdetect_table;
  transpose_all(tdetect_table,detect_table);

  if (gamma_steps % nobs)
    throw Failed("add_GammaPeriodicFID: gamma steps must be a multiple of observations");
  const size_t intsteps=gamma_steps/nobs;

  for (size_t j=0;j<nobs;j++) {
    cmatrix S(Sstore(j));
    cmatrix S2;
    const size_t offset=j*intsteps;
    emultiply(S,tdetect_table(offset),sigma_table(0));
    for (size_t gn=1;gn<gamma_steps;gn++) {
      const size_t which=offset+gn;
      if (which>=gamma_steps)
	mla(S2,tdetect_table(which-gamma_steps),sigma_table(gn));
      else
	mla(S,tdetect_table(which),sigma_table(gn));
    }
    if (isdefined(S2))
      mla(S,R,S2);
  }
  Sstore*=1.0/gamma_steps;
}

void make_S(cmatrix& Sstore,const cmatrix& detect_table,const cmatrix& sigma_table,const BaseList<complex>& R,size_t nobs)
{
  Sstore.create(nobs,R.size());
  const size_t gamma_steps=sigma_table.rows();

  if (gamma_steps % nobs)
    throw Failed("add_GammaPeriodicFID: gamma steps must be a multiple of observations");
  const size_t intsteps=gamma_steps/nobs;

  for (size_t j=0;j<nobs;j++) {
    BaseList<complex> Sj=Sstore.row(j);
    List<complex> S2;
    const size_t offset=j*intsteps;
    multiply(Sj,detect_table.row(offset),sigma_table.row(0));
    for (size_t gn=1;gn<gamma_steps;gn++) {
      const size_t which=offset+gn;
      if (which>=gamma_steps)
	mla(S2,detect_table.row(which-gamma_steps),sigma_table.row(gn));
      else
	mla(Sj,detect_table.row(which),sigma_table.row(gn));
    }
    if (j)
      mla(Sj,R,S2);
  }
  Sstore*=1.0/gamma_steps;
}

// case where sigma0 is conjugate transpose of detect  
void make_S(MultiMatrix<complex,3>& Sstore,const MultiMatrix<complex,3>& table,const cmatrix& R,size_t nobs)
{
  const size_t gamma_steps=table.dimension(0);
  Sstore.create(nobs,R.rows(),R.cols());

  if (gamma_steps % nobs)
    throw Failed("add_GammaPeriodicFID: gamma steps must be a multiple of observations");
  const size_t intsteps=gamma_steps/nobs; 

  for (size_t j=0;j<nobs;j++) {
    cmatrix Sj(Sstore(j));
    if (j) {
      cmatrix S2;
      const size_t offset=j*intsteps;
      conj_emultiply(Sj,table(offset),table(0));
      for (size_t gn=1;gn<gamma_steps;gn++) {
	const size_t which=offset+gn;
	if (which>=gamma_steps)
	  conj_mla(S2,table(which-gamma_steps),table(gn));
	else
	  conj_mla(Sj,table(which),table(gn));
      }
      mla(Sj,R,S2);
    }
    else {
      conj_emultiply(Sj,table(0),table(0));
      for (size_t gn=1;gn<gamma_steps;gn++)
	conj_mla(Sj,table(gn),table(gn));
    }
  }
  Sstore*=1.0/gamma_steps;
}
  
void make_S(cmatrix& Sstore,const cmatrix& table,const BaseList<complex>& R,size_t nobs)
{
  Sstore.create(nobs,R.size());
  const size_t gamma_steps=table.rows();

  if (gamma_steps % nobs)
    throw Failed("add_GammaPeriodicFID: gamma steps must be a multiple of observations");
  const size_t intsteps=gamma_steps/nobs; 

  { // minor optimisation for nobs=1 case
    BaseList<complex> S0=Sstore.row(0);
    conj_multiply(S0,table.row(0),table.row(0));
    for (size_t gn=1;gn<gamma_steps;gn++)
      conj_mla(S0,table.row(gn),table.row(gn));
  }

  for (size_t j=1;j<nobs;j++) {
    BaseList<complex> Sj=Sstore.row(j);
    List<complex> S2;
    const size_t offset=j*intsteps;
    conj_multiply(Sj,table.row(offset),table.row(0));
    for (size_t gn=1;gn<gamma_steps;gn++) {
      const size_t which=offset+gn;
      if (which>=gamma_steps)
	conj_mla(S2,table.row(which-gamma_steps),table.row(gn));
      else
	conj_mla(Sj,table.row(which),table.row(gn));
    }
    mla(Sj,R,S2);
  }
  Sstore*=1.0/gamma_steps;
}

void make_S_mh(MultiHolder<complex>& S,const MultiHolder<complex>& table,const cmatrix& R,size_t nobs)
{
  switch (table.dimensions()) {
  case 3:
    S.dimensions(3);
    make_S(S.multimatrix3(),table.multimatrix3(),R,nobs);
    return;
  case 2:
    S.dimensions(2);
    make_S(S.matrix(),table.matrix(),R.row(),nobs);
    return;
  }
  throw Failed("make_S_mh: not supported");
}

void make_S_mh(MultiHolder<complex>& S,const MultiHolder<complex>& detect_table,const MultiHolder<complex>& sigma_table,const cmatrix& R,size_t nobs)
{
  switch (detect_table.dimensions()) {
  case 3:
    S.dimensions(3);
    make_S(S.multimatrix3(),detect_table.multimatrix3(),sigma_table.multimatrix3(),R,nobs);
    return;
  case 2:
    S.dimensions(2);
    make_S(S.matrix(),detect_table.matrix(),sigma_table.matrix(),R.row(),nobs);
    return;
  }
  throw Failed("make_S_mh: not supported");
}

  void GammaPeriodicFID::doadd(BaseList<complex> FIDv,complex scale) const
  {
    const MultiMatrix<complex,3>& Sstore(S.multimatrix3());
    if (verbose_>1) {
      std::cout << "R matrix\n" << R << std::endl;
      std::cout << "S matrices\n" << Sstore << std::endl;
    }
    
  const size_t nobs=Sstore.dimension(0);
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_GammaPeriodicFID: zero-length FID!");
  const size_t rows=R.rows();
  const size_t cols=R.cols();

  for (size_t j=0;j<nobs;j++) {
    const cmatrix& Sj=Sstore(j);
    for (size_t r=rows;r--;) {
      for (size_t s=cols;s--;) {
	complex a=Sj(r,s)*scale;
	if (norm(a)>zerotol) {
	  const complex p(getR(r,s));
	  for (size_t k=j;k<totpoints;k+=nobs) {
	    FIDv(k)+=a;
	    a*=p;
	  }
	}
      }
    }
  }
}

complex GammaPeriodicFID::getR(size_t r, size_t s) const
{
  return quickpow(R(r,s),reduction_factor_);
}

  template<class T> void GammaPeriodicFID::doadd_hermitian(BaseList<T> FIDv,T scale,bool forceherm) const
{
  const MultiMatrix<complex,3>& Sstore(S.multimatrix3());
  scale/=this->gammasteps();
  const size_t nobs=Sstore.dimension(0);
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_GammaPeriodicFID_hermitian: zero-length FID!");
  const size_t rows=R.rows();
  const size_t cols=R.cols();
  const T scale2(2.0*scale);

  for (size_t j=0;j<nobs;j++) {
    const cmatrix& Sj=Sstore(j);
    if (forceherm) {
      for (size_t r=rows;r--;) {
	for (size_t s=cols;s--;) {
	  complex a=Sj(r,s)*scale2;
	  if (norm(a)>zerotol) {
	    const complex p(getR(r,s));
	    for (size_t k=j;k<totpoints;k+=nobs) {
	      FIDv(k)+=real(a);
	      a*=p;
	    }
	  }
	}
      }
    }
    else {
      const double amp(hermitian_trace(Sj));
      if (amp) {
	const T con(amp*scale);
	for (size_t k=j;k<totpoints;k+=nobs)
	  FIDv(k)+=con;
      }
      for (size_t r=rows;r--;) {
	for (size_t s=r;s--;) {
	  complex a=Sj(r,s)*scale2;
	  if (norm(a)>zerotol) {
	    const complex p(getR(r,s));
	    for (size_t k=j;k<totpoints;k+=nobs) {
	      FIDv(k)+=real(a);
	      a*=p;
	    }
	  }
	}
      }
    }
  }
}

  void GammaPeriodicFID::doadd_compressed(BaseList<complex> FIDv,complex scale) const
{
  const cmatrix& Sstore(S.matrix());
  const BaseList<complex>& Rrow(R.row());
  scale/=this->gammasteps();
  const size_t nobs=Sstore.rows();
  const size_t totpoints=FIDv.size();
  if (totpoints==0)
    throw Undefined("add_GammaPeriodicFID: zero-length FID!");
  const size_t dim=Rrow.size();

  for (size_t j=nobs;j--;) {
    const BaseList<complex> Sj=Sstore.row(j);
    for (size_t s=dim;s--;) {
      complex a=Sj(s)*scale;
      if (norm(a)>zerotol) {
	const complex p=quickpow(Rrow(s),reduction_factor_);
	for (size_t k=j;k<totpoints;k+=nobs) {
	  FIDv(k)+=a;
	  a*=p;
	}
      }
    }
  }
}

void make_tables(cmatrix& detect_table,cmatrix& sigma_table,const BaseList<complex>& URs,const BaseList<cmatrix>& UCs,const BaseList<complex>& sigma0,const BaseList<complex>& detect,const cmatrix& DC)
{
  const size_t gamma_steps=UCs.size();
  const size_t dim=DC.rows();
 
  detect_table.create(gamma_steps,dim);
  sigma_table.create(gamma_steps,dim);

  for (size_t gn=0;gn<gamma_steps;gn++) {
    BaseList<complex> det_row=detect_table.row(gn);
    BaseList<complex> sig_row=sigma_table.row(gn);
    if (gn==0) {
      conj_transpose_multiply(det_row,DC,detect);
      multiply(sig_row,sigma0,DC);
    }
    else {
      unitary_isimtransLR(det_row,UCs(gn-1),detect,URs(gn-1));
      unitary_isimtransLR(sig_row,URs(gn-1),sigma0,UCs(gn-1));
    }
  }
}

void make_tables(cmatrix& detect_table,cmatrix& sigma_table,const BaseList<cmatrix>& URs,const BaseList<complex>& UCs,const BaseList<complex>& sigma0,const BaseList<complex>& detect,const cmatrix& DR)
{
  const size_t gamma_steps=UCs.size();
  const size_t dim=DR.rows();
 
  detect_table.create(gamma_steps,dim);
  sigma_table.create(gamma_steps,dim);

  for (size_t gn=0;gn<gamma_steps;gn++) {
    BaseList<complex> det_row=detect_table.row(gn);
    BaseList<complex> sig_row=sigma_table.row(gn);
    if (gn==0) {
      multiply(det_row,detect,DR);
      conj_transpose_multiply(sig_row,DR,sigma0);
    }
    else {
      unitary_isimtransLR(det_row,UCs(gn-1),detect,URs(gn-1));
      unitary_isimtransLR(sig_row,URs(gn-1),sigma0,UCs(gn-1));
    }
  }
}

struct transform_GammaTD_ {
  template<class T,class MR,class MC> void operator() (MultiMatrix<complex,3>& table,const MR& URs,const MC& UCs,const T& a) const
  {
    table_create(table,URs,UCs);
    for (size_t gn=URs.dimension(0);gn--;) {
      cmatrix tablegn(table(gn));
      unitary_isimtransLR(tablegn,URs(gn),a,UCs(gn));
    }
  }

  template<class TR,class TC> void operator() (cmatrix& table,const TR& URs,const TC& UCs,const BaseList<complex>& a) const
  {
    table_create(table,URs,UCs);
    for (size_t gn=URs.dimension(0);gn--;)
       unitary_isimtransLR(table.row(gn),URs(gn),a,UCs(gn));
  }
  template<class T> void operator() (cmatrix& table,const BaseList<complex>& URs,const BaseList<complex>& UCs,T a) const
  {
    const size_t gamma_steps=URs.size();
    if (UCs.size()!=gamma_steps)
      throw Mismatch("transform_GammaTD_");
    table.create(gamma_steps,1);
    BaseList<complex> tabler=table.row();
    for (size_t gn=0;gn<gamma_steps;gn++)
      unitary_isimtransLR(tabler(gn),URs(gn),a,UCs(gn));
  }

  void operator() (cmatrix& table,const BaseList<complex>& URs,const BaseList<complex>& UCs) const
  {
    table_create(table,URs,UCs);
    BaseList<complex> tabler(table.row());
    for (size_t gn=URs.size();gn--;)
      tabler(gn)=unitary_isimtrans_identity(URs(gn),UCs(gn));
  }
  void operator() (MultiMatrix<complex,3>& table,const MultiMatrix<complex,3>& URs,const MultiMatrix<complex,3>& UCs) const
  {
    table_create(table,URs,UCs);
    for (size_t gn=URs.dimension(0);gn--;) {
      cmatrix tablegn(table(gn));
      unitary_isimtrans_identity(tablegn,URs(gn),UCs(gn));
    }
  }
};

template<class T,class M> void make_tables(MultiMatrix<complex,3>& detect_table,MultiMatrix<complex,3>& sigma_table,const M& URs,const M& UCs,const T& sigma0,const T& detect,const cmatrix& DR,const cmatrix& DC)
{
  table_create(detect_table,URs,UCs);
  table_create(sigma_table,URs,UCs);

  cmatrix UTL,UTR,tmp;
  const bool sameLR=(&URs==&UCs);
  cmatrix& use_UTR=sameLR ? UTL : UTR;

  for (size_t gn=URs.dimension(0);gn--;) {
    cmatrix detectgn(detect_table(gn));
    cmatrix sigmagn(sigma_table(gn));

    if (gn==0) {
      unitary_isimtrans_using(detectgn,DC,detect,DR,tmp);
      unitary_isimtrans_using(sigmagn,DR,sigma0,DC,tmp);
    }
    else {
      multiply(UTL,URs(gn-1),DR);
      if (!sameLR) multiply(UTR,UCs(gn-1),DC);
      unitary_isimtrans_using(detectgn,use_UTR,detect,UTL,tmp);
      unitary_isimtrans_using(sigmagn,UTL,sigma0,use_UTR,tmp);
    }
    detectgn.transpose();
  }
}

void make_tables(MultiMatrix<complex,3>& detect_table,MultiMatrix<complex,3>& sigma_table,const cmatrix& URs,const cmatrix& UCs,const cmatrix& sigma0,const cmatrix& detect)
{
  const size_t gamma_steps=URs.size();
  if (gamma_steps!=UCs.size())
    throw Mismatch("make_tables");
 
  detect_table.create(gamma_steps,detect.rows(),detect.cols());
  sigma_table.create(gamma_steps,sigma0.rows(),sigma0.cols());

  for (size_t gn=0;gn<gamma_steps;gn++) {
    if (gn==0) {
      detect_table(0)=detect;
      sigma_table(0)=sigma0;
    }
    else {
      cmatrix detectgn(detect_table(gn));
      cmatrix sigmagn(sigma_table(gn));
      unitary_isimtransLR(detectgn,UCs.row(gn-1),detect,URs.row(gn-1));
      unitary_isimtransLR(sigmagn,URs.row(gn-1),sigma0,UCs.row(gn-1));
    }
    detect_table(gn).transpose();
  }
}

template<class T,class M> void make_table(MultiMatrix<complex,3>& table,const M& URs,const M& UCs,const T& sigma0det,const cmatrix& DR,const cmatrix& DC)
{
  table_create(table,URs,UCs);

  cmatrix UTL,UTR,tmp;
  const bool sameLR=(&DC==&DR);
  cmatrix& use_UTR=sameLR ? UTL : UTR;

  for (size_t gn=URs.dimension(0);gn--;) {
    cmatrix tablegn(table(gn));
    if (gn==0)
      unitary_isimtrans_using(tablegn,DR,sigma0det,DC,tmp);
    else {
      multiply(UTL,URs(gn-1),DR);
      if (!sameLR) multiply(UTR,UCs(gn-1),DC);
      unitary_isimtrans_using(tablegn,UTL,sigma0det,use_UTR,tmp);
    }
  }
}

template<class M> void make_table(cmatrix& table,const BaseList<complex>& URs,const M& UCs,const BaseList<complex>& sigma0det,const cmatrix& DC)
{
  const size_t gamma_steps=URs.size();
  const size_t dim=DC.rows();
  table.create(gamma_steps,dim);
 
   for (size_t gn=0;gn<gamma_steps;gn++) {
     BaseList<complex> tab_row=table.row(gn);
     if (gn==0)
       multiply(tab_row,sigma0det,DC);
     else
       unitary_isimtransLR(tab_row,URs(gn-1),sigma0det,UCs(gn-1));
  }
}

void make_table(cmatrix& table,const BaseList<cmatrix>& URs,const BaseList<complex>& UCs,const BaseList<complex>& sigma0det,const cmatrix& DR)
{
  const size_t gamma_steps=URs.size();

  table.create(gamma_steps,DR.rows());

  for (size_t gn=0;gn<gamma_steps;gn++) {
    BaseList<complex> tab_row=table.row(gn);
    if (gn==0)
      conj_transpose_multiply(tab_row,DR,sigma0det);
    else
      unitary_isimtransLR(tab_row,URs(gn-1),sigma0det,UCs(gn-1));
  }
}

template<class M> void make_table(MultiMatrix<complex,3>& table,const M& URs,const M& UCs,const cmatrix& DR,const cmatrix& DC)
{
  table_create(table,URs,UCs);

  cmatrix UTR,UTC;
  for (size_t gn=URs.dimension(0);gn--;) {
    cmatrix tablegn(table(gn));
    if (gn==0) 
      unitary_isimtrans_identity(tablegn,DR,DC);
    else {
      multiply(UTR,URs(gn-1),DR);
      multiply(UTC,UCs(gn-1),DC);
      unitary_isimtrans_identity(tablegn,UTR,UTC);
    }
  }
}

//static const double TWO_PI=2.0*M_PI;

void add_FID(BaseList<complex> FIDv,complex a,double f)
{
  complex p=expi(TWO_PI*f);
  const size_t n=FIDv.size();
  if (n==0)
    throw Undefined("add_FID: zero-length FID!");
  for (size_t i=0;i<n;i++) {
    FIDv(i)+=a;
    a*=p;
  }
}

void add_FID(BaseList<double> FIDv,complex a,double f)
{
  complex p=expi(TWO_PI*f);
  const size_t n=FIDv.size();
  if (n==0)
    throw Undefined("add_FID: zero-length FID!");
  for (size_t i=0;i<n;i++) {
    FIDv(i)+=real(a);
    a*=p;
  }
}

template<class F> void transform_TD(MultiHolder<complex>& detect_tab,const F& obj,const MultiHolder<complex>& URs,const MultiHolder<complex>& UCs,const cmatrix& detect)
{
  if (URs.dimensions()==3) {
    if (UCs.dimensions()==3) {
      detect_tab.dimensions(3);
      obj(detect_tab.multimatrix3(),URs.multimatrix3(),UCs.multimatrix3(),detect);
    }
    else {
      detect_tab.dimensions(2);
      obj(detect_tab.matrix(),URs.multimatrix3(),UCs.list(),detect.row());
    }
  }
  else {
    if (UCs.dimensions()==3) {
      detect_tab.dimensions(2);
      obj(detect_tab.matrix(),URs.list(),UCs.multimatrix3(),detect.row());
    }
    else {
      if (detect.rows()!=1 || detect.cols()!=1)
	throw Mismatch("PeriodicFID::add_FID");
      detect_tab.dimensions(1);
      obj(detect_tab.list(),URs.list(),UCs.list(),detect(0,0));
    }
  }
}

void PeriodicFID::add_FID(BaseList<complex> FIDv, complex scale, const cmatrix& sigma0,const cmatrix& detect) const
{
  if (!*this)
    throw Failed("Propagators not set");

  const MultiHolder<complex>& URs=row().Us;
  const MultiHolder<complex>& UCs=col().Us;

  transform_TD(detect_tab,transform_TD_(),UCs,URs,detect);

  if (URs.dimensions()==3) {
    const MultiMatrix<complex,3>& Rs(URs.multimatrix3());
    const size_t nprops=Rs.dimension(0);
    if (UCs.dimensions()==3)
      doadd_PeriodicFID(FIDv,scale,detect_tab.multimatrix3(),Rs(nprops-1),sigma0,UCs.multimatrix3()(nprops-1),tmp);
    else
      doadd_PeriodicFID(FIDv,scale,detect_tab.matrix(),Rs(nprops-1),sigma0.row(),UCs.list()(nprops-1));
  }
  else {
    const List<complex>& Rs=URs.list();
    const size_t nprops=Rs.size();
    if (UCs.dimensions()==3)
      doadd_PeriodicFID(FIDv,scale,detect_tab.matrix(),Rs(nprops-1),sigma0.row(),UCs.multimatrix3()(nprops-1));
    else {
      if (sigma0.rows()!=1 || sigma0.cols()!=1)
	throw Mismatch("PeriodicFID::add_FID");
      doadd_PeriodicFID(FIDv,scale,detect_tab.list(),Rs(nprops-1),sigma0(0,0),UCs.list()(nprops-1));
    }
  }
}

template<class T> void PeriodicFID::add_FID_(BaseList<T> FIDv, T scale, const BaseList<double>& sigma0,const BaseList<double>& detect) const
{
  if (!isdiagonal())
    throw Failed("PeriodicFID: diagonal operators only valid for diagonal blocks");
  if (!*this)
    throw Failed("Propagators not set");

  const MultiMatrix<complex,3>& Us(row().Us.multimatrix3());
  transform_TD_ obj;
  detect_tab.dimensions(3);
  obj(detect_tab.multimatrix3(),Us,Us,detect);

  const cmatrix Ucycle(Us.back());
  cmatrix sigma0f;
  full(sigma0f,sigma0);
  doadd_PeriodicFID_hermitian(FIDv,scale,detect_tab.multimatrix3(),Ucycle,sigma0f,Ucycle,tmp);
}

template void PeriodicFID::add_FID_(BaseList<double>, double, const BaseList<double>&, const BaseList<double>&) const;
template void PeriodicFID::add_FID_(BaseList<complex>, complex, const BaseList<double>&, const BaseList<double>&) const;

template<class T> void PeriodicFID::add_FID_hermitian(BaseList<T> FIDv, T scale, const cmatrix& sigma0,const cmatrix& detect) const
{
  if (!*this)
    throw Failed("Propagators not set");

  detect_tab.dimensions(3);
  MultiMatrix<complex,3>& detect_tab_ref(detect_tab.multimatrix3());

  transform_TD_ obj;
  const MultiMatrix<complex,3>& URs(row().Us.multimatrix3());
  const MultiMatrix<complex,3>& UCs(col().Us.multimatrix3());
  obj(detect_tab_ref,UCs,URs,detect);

  doadd_PeriodicFID_hermitian(FIDv,scale,detect_tab_ref,URs.back(),sigma0,UCs.back(),tmp,!isdiagonal());
}

template void PeriodicFID::add_FID_hermitian(BaseList<double>, double, const cmatrix&, const cmatrix&) const;
template void PeriodicFID::add_FID_hermitian(BaseList<complex>, complex, const cmatrix&, const cmatrix&) const;

void PeriodicFID::add_FID(BaseList<complex> FIDv, complex scale) const
{
//   if (isdiagonal())
//     throw Failed("add_FID: only valid for off-diagonal blocks");
  if (!*this) 
    throw Failed("Propagators not set");

  if ( (row().Us.dimensions()!=3) || (col().Us.dimensions()!=3) )
    throw InternalError("add_FID: illegal combination");
  const MultiMatrix<complex,3>& URs(row().Us.multimatrix3());
  if (isdiagonal()) {
    scale*=URs.back().rows();
    if (verbose_)
      std::cout << "PeriodicFID::add_FID: identity matrix optimisation: adding constant " << scale << '\n';
    FIDv+=scale;
  }
  else {
    detect_tab.dimensions(3);
    MultiMatrix<complex,3>& dest(detect_tab.multimatrix3());
    
    transform_TD_ obj;
    const MultiMatrix<complex,3>& UCs(col().Us.multimatrix3());
    obj(dest,UCs,URs); //Transforming detection operator, hence swap R/C
    doadd_PeriodicFID(FIDv,scale,dest,URs.back(),UCs.back(),tmp);
  }
}

void HomoStashGammaTD::operator= (const BaseList<complex>& Usv)
{
  size_t n=Usv.size();
  cycle_eigs.create(1,Usv(n-1));
  Us.dimensions(1);
  BaseList<complex> Urow=Us.list();
  Urow.create(n);
  Urow(0)=1.0;
  for (size_t gn=1;gn<n;gn++)
    Urow(gn)=Usv(gn-1);
}

void HomoStashGammaTD::operator= (const BaseList<cmatrix>& Usv)
{
  const size_t n=Usv.size();
  const cmatrix& Ucycle=Usv(n-1);
  const size_t dim=Ucycle.rows();
  cycle_eigs.create(dim);
  Us.dimensions(3);
  MultiMatrix<complex,3>& UT(Us.multimatrix3());
  UT.create(n,dim,dim);
  cmatrix D(UT(0));
  if (eigensystem(D,cycle_eigs,Usv(n-1),offdiagonal_tolerance)) {
    propagation_closetodiagonal_warning.raise();
    for (size_t gn=1;gn<n;gn++)
      UT(gn)=Usv(gn-1);
  }
  else {
    for (size_t gn=1;gn<n;gn++) {
      cmatrix UTgn(UT(gn));
      multiply(UTgn,Usv(gn-1),D);
    }
  }
}

void GammaPeriodicFID::reduction_factor(int factorv)
{
  if (reduction_factor_==factorv)
    return;
  if (factorv<1)
    throw InvalidParameter("GammaPeriodicFID::reduction_factor must be >0");
  reduction_factor_=factorv;
  clear(); //!< invalid contents
}

void GammaPeriodicFID::make_R() const
{
  if (isdiagonal())
    evolution_matrix(R,row().cycle_eigs);
  else
    evolution_matrix(R,row().cycle_eigs,col().cycle_eigs);
}

template<class T> void GammaPeriodicFID::add_FID_(BaseList<T> FIDv, T scale, const BaseList<double>& sigma0,const BaseList<double>& detect) const
{
  if (!isdiagonal())
    throw Failed("GammaPeriodicFID: diagonal operators only valid for diagonal blocks");
  if (!*this)
    throw Failed("Propagators not set");

  transform_mh(detect_tab,transform_GammaTD_(),row().Us,detect);
  transform_mh(sigma0_tab,transform_GammaTD_(),row().Us,sigma0);
  make_R();
  make_S_mh(S,detect_tab,sigma0_tab,R,this->observations());
  doadd_hermitian(FIDv,scale);
}

template void GammaPeriodicFID::add_FID_(BaseList<double>, double, const BaseList<double>&, const BaseList<double>&) const;
template void GammaPeriodicFID::add_FID_(BaseList<complex>, complex, const BaseList<double>&, const BaseList<double>&) const;

void GammaPeriodicFID::add_FID(BaseList<complex> FIDv, complex scale, const cmatrix& sigma0,const cmatrix& detect) const
{
  if (!*this)
    throw Failed("Propagators not set");
  transform_mh(detect_tab,transform_GammaTD_(),col().Us,row().Us,detect);
  transform_mh(sigma0_tab,transform_GammaTD_(),row().Us,col().Us,sigma0);
  make_R();
  make_S_mh(S,detect_tab,sigma0_tab,R,this->observations());
  if (verbose_>1) {
    std::cout << "Transformed sigma0\n" << sigma0_tab << '\n';
    std::cout << "Transformed detect\n" << detect_tab << '\n';
  }
  switch (S.dimensions()) {
  case 3:
    doadd(FIDv,scale);
    return;
  case 2:
    doadd_compressed(FIDv,scale);
    return;
  }
  throw InternalError("GammaPeriodicFID::add_FID");
}

template<class T> void GammaPeriodicFID::add_FID_hermitian(BaseList<T> FIDv, T scale, const cmatrix& sigma0,const cmatrix& detect) const
{
  if (!*this)
    throw Failed("Propagators not set");
  transform_mh(detect_tab,transform_GammaTD_(),col().Us,row().Us,detect);
  transform_mh(sigma0_tab,transform_GammaTD_(),row().Us,col().Us,sigma0);
  make_R();
  make_S_mh(S,detect_tab,sigma0_tab,R,this->observations());
  doadd_hermitian(FIDv,scale,!isdiagonal());
}

template void GammaPeriodicFID::add_FID_hermitian(BaseList<double>, double, const cmatrix&, const cmatrix&) const;
template void GammaPeriodicFID::add_FID_hermitian(BaseList<complex>, complex, const cmatrix&, const cmatrix&) const;

template<class T> void GammaPeriodicFID::add_FID_(BaseList<T> FIDv, T scale, const BaseList<double>& sigma0det) const
{
  if (!*this)
    throw Failed("Propagators not set");
  if (!isdiagonal())
    throw Failed("add_FID: only valid for diagonal blocks");
  transform_mh(detect_tab,transform_GammaTD_(),row().Us,sigma0det);
  sigma0_tab.clear();
  make_R();
  make_S_mh(S,detect_tab,R,this->observations());
  doadd_hermitian(FIDv,scale);
}

template void GammaPeriodicFID::add_FID_(BaseList<double>, double, const BaseList<double>&) const;
template void GammaPeriodicFID::add_FID_(BaseList<complex>, complex, const BaseList<double>&) const;

void GammaPeriodicFID::add_FID(BaseList<complex> FIDv, complex scale, const cmatrix& sigma0det) const
{
  if (!*this) 
    throw Failed("Propagators not set");
  transform_mh(detect_tab,transform_GammaTD_(),row().Us,col().Us,sigma0det);
  sigma0_tab.clear();
  make_R();
  make_S_mh(S,detect_tab,R,this->observations());
  if (verbose_>1)
    std::cout << "Transformed Qs\n" << detect_tab << '\n';
  switch (S.dimensions()) {
  case 3:
    doadd(FIDv,scale);
    return;
  case 2:
    doadd_compressed(FIDv,scale);
    return;
  }
  throw InternalError("GammaPeriodicFID::add_FID");
}

template<class T> void GammaPeriodicFID::add_FID_hermitian(BaseList<T> FIDv, T scale, const cmatrix& sigma0det) const
{
  if (!*this)
    throw Failed("Propagators not set");
  //  if (!isdiagonal())
  //throw Failed("add_FID_hermitian: only valid for diagonal blocks");
  transform_mh(detect_tab,transform_GammaTD_(),row().Us,col().Us,sigma0det);
  sigma0_tab.clear();
  make_R();
  make_S_mh(S,detect_tab,R,this->observations());
  doadd_hermitian(FIDv,scale,!isdiagonal());
}

template void GammaPeriodicFID::add_FID_hermitian(BaseList<double> FID, double, const cmatrix&) const;
template void GammaPeriodicFID::add_FID_hermitian(BaseList<complex> FID, complex, const cmatrix&) const;

void GammaPeriodicFID::add_FID(BaseList<complex> FIDv, complex scale) const
{
  if (!*this) 
    throw Failed("Propagators not set");
  if (isdiagonal()) {
    //    throw Failed("add_FID: only valid for off-diagonal blocks");
    scale*=col().size();
    if (verbose_)
      std::cout << "GammaPeriodicFID::add_FID: identity matrix optimisation - adding constant " << scale << '\n';
    FIDv+=scale;
  }
  else {  
    transform_mh(detect_tab,transform_GammaTD_(),row().Us,col().Us);
    sigma0_tab.clear();
    make_R();
    make_S_mh(S,detect_tab,R,this->observations());
    doadd(FIDv,scale);
  }
}

// void GammaPeriodicFID::add_FID(BaseList<complex> FID, complex scale) const
// {
//   if (isdiagonal())
//     evolution_matrix(R,row().cycle_eigs);    
//   else
//     evolution_matrix(R,row().cycle_eigs,col().cycle_eigs);

//   scale*=scale_factor;
//   switch (S.dimensions()) {
//   case 3:
//     doadd_GammaPeriodicFID(FID,scale/gamma_steps,S.multimatrix3(),R);
//     return;
//   case 2:
//     doadd_GammaPeriodicFID(FID,scale/gamma_steps,S.matrix(),R.row());
//     return;
//   }
//   throw Failed("GammaPeriodicFID: not supported");
// }

std::ostream& operator<< (std::ostream& ostr, const GammaPeriodicFID& a)
{
  if (a.verbose_)
    ostr << static_cast<const SwapStore<HomoStashGammaTD>&>(a);
  return ostr << static_cast<const PeriodicObj&>(a);
}
   
std::ostream& operator<< (std::ostream& ostr, const PeriodicFID& a)
{
  if (a.verbose_) 
    ostr << static_cast<const SwapStore<HomoStashTD>&>(a);
  return ostr << static_cast<const PeriodicObj&>(a);
}

void PeriodicFID::set_Us(const BaseList<cmatrix>& Usv) {
  setdiagonal()=Usv;
  this->set_steps(Usv.size());
}
void PeriodicFID::set_Us(const BaseList<complex>& Usv) {
  this->set_steps(Usv.size());
  setdiagonal()=Usv;
}
void PeriodicFID::set_Us(char sel,const BaseList<cmatrix>& Usv) { 
  setRC(sel)=Usv;
  this->set_steps(Usv.size());
}
void PeriodicFID::set_Us(char sel,const BaseList<complex>& Usv) {
  setRC(sel)=Usv;
  this->set_steps(Usv.size()); 
}

GammaPeriodicFID::GammaPeriodicFID(int nobsv, int verbosev)
  : PeriodicObj(nobsv,verbosev), reduction_factor_(1)
{
  //  row().showwarning=col().showwarning=(verbosev>0);
}

void GammaPeriodicFID::set_Us(const BaseList<cmatrix>& Usv)
{
  setdiagonal()=Usv;
  this->set_steps(Usv.size());
}
void GammaPeriodicFID::set_Us(const BaseList<complex>& Usv) {
  this->set_steps(Usv.size());
  setdiagonal()=Usv;
}
void GammaPeriodicFID::set_Us(char sel,const BaseList<cmatrix>& Usv)
{ 
  setRC(sel)=Usv;
  this->set_steps(Usv.size());
}
void GammaPeriodicFID::set_Us(char sel,const BaseList<complex>& Usv)
{
  setRC(sel)=Usv;
  this->set_steps(Usv.size());
}
   
}//namespace libcmatrix
