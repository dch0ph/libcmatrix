#include "NMR.h"
#include "HomoProp_shared.h"
#include "lcm_HomoProp.h"
#include "ScratchList.h"
#include "MultiMatrix.h"
#include <cstring>
#include "timer.h"

namespace libcmatrix {

static const char NOPROPS[]="observe: propagators not defined";

void multiply(cmatrix& dest,const rmatrix& a,const rmatrix& b)
{
  rmatrix tmpr;
  multiply(tmpr,a,b);
  dest=tmpr;
}

void unitary_isimtrans_identity(cmatrix& d, const BaseList<complex>& VL, const BaseList<complex>& VR)
{
  size_t n=VL.size();
  if (n!=VR.size())
    throw Mismatch("unitary_isimtrans_identity");
  d.create(n,n);
  d=0.0;
  for (;n--;)
    d(n,n)=conj_multiply(VL(n),VR(n));
}

void conj_mla(cmatrix& a, const cmatrix& b, const cmatrix& c)
{
  if (!a)
    conj_emultiply(a,b,c);
  else {
    if (!arematching(b,c) || !arematching(a,b))
      throw Mismatch("conj_mla");
//     const complex *bd=b.vector();
//     const complex *cd=c.vector();
//     complex *ad=a.vector();
//     for (size_t i=a.size();i--;) 
//       conj_mla(ad[i],bd[i],cd[i]);
    apply_ip3(doesconj_mla<complex>(),a,b,c);
  }
}

struct transform_homo_ {
  template<class T> inline void operator() (cmatrix& table,const BaseList<complex>& UTRs,const MultiMatrix<complex,3>& UTCs,const BaseList<T>& a) const
  {
    size_t gn=UTRs.size();
    table.create(gn,a.size());
    for (;gn--;)
      unitary_isimtransLR(table.row(gn),UTRs(gn),a,UTCs(gn));
  }
  
  template<class T> inline void operator() (cmatrix& table,const MultiMatrix<complex,3>& UTRs,const BaseList<complex>& UTCs,const BaseList<T>& a) const
  {
    size_t gn=UTCs.size();
    table.create(gn,a.size());
    for (;gn--;)
      unitary_isimtransLR(table.row(gn),UTRs(gn),a,UTCs(gn));
  }
  
  inline void operator() (cmatrix& table,const BaseList<complex>& UTRs,const BaseList<complex>& UTCs) const
  {
    size_t gn=UTCs.size();
    table.create(gn,1);
    BaseList<complex> tabler=table.row();
    for (;gn--;)
      tabler(gn)=unitary_isimtrans_identity(UTRs(gn),UTCs(gn));
  }
  template<typename T> inline void operator() (cmatrix& table,const BaseList<complex>& UTRs,const BaseList<complex>& UTCs,T a) const
  {
    size_t gn=UTCs.size();
    table.create(gn,1);
    BaseList<complex> tabler=table.row();
    for (;gn--;)
      unitary_isimtransLR(tabler(gn),UTRs(gn),a,UTCs(gn));
  }
  
  template<typename T> inline void operator() (MultiMatrix<complex,3>& table,const MultiMatrix<complex,3>& UTRs,const MultiMatrix<complex,3>& UTCs,const T& a) const
  {
    size_t gn=UTCs.dimension(0);
    table.create(gn,UTRs.dimension(1),UTCs.dimension(2));
    cmatrix tmp;
    for (;gn--;) {
      Matrix<complex> tablegn(table(gn));
      unitary_isimtransLR(tablegn,UTRs(gn),a,UTCs(gn),&tmp);
    }
  }

  inline void operator() (MultiMatrix<complex,3>& table,const MultiMatrix<complex,3>& UTRs,const MultiMatrix<complex,3>& UTCs) const
  {
    size_t gn=UTCs.dimension(0);
    table.set_dimensions(UTRs);
    for (;gn--;) {
      Matrix<complex> tablegn(table(gn));
      unitary_isimtrans_identity(tablegn,UTRs(gn),UTCs(gn));
    }
  }
};
  
template<typename T> void do_transform_mh0(MultiHolder<complex>& table,const MultiHolder<complex>& UTRs,const MultiHolder<complex>& UTCs,const Matrix<T>& a)
{ 
  table.dimensions(2);
  cmatrix& dest=table.matrix();

  if (UTRs.dimensions()==3) {
    const MultiMatrix<complex,3>& Rs=UTRs.multimatrix3();
    if (UTCs.dimensions()==3)
      unitary_isimtransLR(dest,Rs(0),a,UTCs.multimatrix3()(0));
    else {
      dest.create(a.rows(),1);
      unitary_isimtransLR(dest.row(),Rs(0),a.row(),UTCs.list()(0));
    }
  }
  else {
    const List<complex>& Rs=UTRs.list();
    if (UTCs.dimensions()==3) {
      dest.create(1,a.cols());
      unitary_isimtransLR(dest.row(),Rs(0),a.row(),UTCs.multimatrix3()(0));
    }
    else {
      if (a.rows()!=1 || a.cols()!=1)
	throw Mismatch("do_transform_mh0");
      dest.create(1,1);
      unitary_isimtransLR(dest(0,0),Rs(0),a(0,0),UTCs.list()(0));
    }
  }
}

void do_transform_mh0(MultiHolder<complex>& table, const MultiHolder<complex>& UTs, const BaseList<double>& a)
{
  table.dimensions(2);
  cmatrix& dest=table.matrix();

  if (UTs.dimensions()==3)
    unitary_isimtrans(dest,a,UTs.multimatrix3()(0));
  else {
    if (a.size()!=1)
      throw Mismatch("do_transform");
    dest.create(1,1,complex(a.front()));  //OK, this is wasteful, but it's a rare case!
  }
  throw InternalError("do_transform");
}

template<class T> void
BasePeriodicSpectrum<T>::reset()
{
  if (!*this)
    lock();
  else {
    dimR=rows();
    dimC=cols();
    r=s=ncount=0;
    finished=false;
    fscale=1.0/this->period(); //conversion from fractions of cycle period to Hz
  }
  setclean();
}

template<class T> void
BasePeriodicSpectrum<T>::sideband(int val)
{
  eigrestrict=true;
  valsideband=val;
}
  
const char LCM_INCOMPAT[]="sigma0 and detection operators incompatible (transpose relationship) or undefined";

void GammaPeriodicSpectrum::observe(const cmatrix& sigma0,const cmatrix& det)
{
  if (!aretranspose(sigma0,det))
    throw Failed(LCM_INCOMPAT);
  if (!*this)
    throw Failed(NOPROPS);
  transform_mh(sigma_table,transform_homo_(),row().Us,col().Us,sigma0);
  transform_mh(detect_table,transform_homo_(),col().Us,row().Us,det);
  reset();
}

void PeriodicSpectrum::observe(const cmatrix& sigma0,const cmatrix& det)
{
  if (!aretranspose(sigma0,det))
    throw Failed(LCM_INCOMPAT);
  if (!*this)
    throw Failed(NOPROPS);
  do_transform_mh0(sigma_table,row().Us,col().Us,sigma0);
  transform_mh(detect_table,transform_homo_(),col().Us,row().Us,det);
  reset();
}

void PeriodicSpectrum::observe()
{
  if (!*this)
    throw Failed(NOPROPS);
  transform_mh(detect_table,transform_homo_(),row().Us,col().Us);
  reset();
}

void PeriodicSpectrum::observe(const cmatrix& sigma0det)
{
  if (!*this)
    throw Failed(NOPROPS);
  transform_mh(detect_table,transform_homo_(),col().Us,row().Us,sigma0det);
  sigma_table.clear();
  reset();
}

void PeriodicSpectrum::observe(const rmatrix& sigma0det)
{
  if (!*this)
    throw Failed(NOPROPS);
  transform_mh(detect_table,transform_homo_(),col().Us,row().Us,sigma0det);
  sigma_table.clear();
  reset();
}

void PeriodicSpectrum::observe(const BaseList<double> &sigma0,const BaseList<double> &det)
{
  if (!isdiagonal())
    throw Failed("Diagonal operators only valid for diagonal blocks");
  if (!arematching(sigma0,det))
    throw Failed(LCM_INCOMPAT);
  if (!*this)
    throw Failed(NOPROPS);
  do_transform_mh0(sigma_table,row().Us,sigma0);
  transform_mh(detect_table,transform_homo_(),row().Us,det);
  reset();
}

void GammaPeriodicSpectrum::observe(const BaseList<double> &sigma0,const BaseList<double> &det)
{
  if (!isdiagonal()) 
    throw Failed("Diagonal operators only valid for diagonal blocks");
  if (!arematching(sigma0,det))
    throw Failed(LCM_INCOMPAT);
  if (!*this) 
    throw Failed(NOPROPS);
  transform_mh(sigma_table,transform_homo_(),row().Us,sigma0);
  transform_mh(detect_table,transform_homo_(),row().Us,det);
  reset();
}

void GammaPeriodicSpectrum::observe(const rmatrix& sigma0,const rmatrix& det)
{
  if (!aretranspose(sigma0,det))
    throw Failed(LCM_INCOMPAT);
  if (!*this)
    throw Failed(NOPROPS);
  transform_mh(sigma_table,transform_homo_(),row().Us,col().Us,sigma0);
  transform_mh(detect_table,transform_homo_(),col().Us,row().Us,det);
  reset();
}

void PeriodicSpectrum::observe(const rmatrix& sigma0,const rmatrix& det)
{
  if (!aretranspose(sigma0,det))
    throw Failed(LCM_INCOMPAT);
  if (!*this)
    throw Failed(NOPROPS);
  do_transform_mh0(sigma_table,row().Us,col().Us,sigma0);
  transform_mh(detect_table,transform_homo_(),col().Us,row().Us,det);
  reset();
}

void PeriodicSpectrum::observe(const BaseList<double>& sigma0det)
{
  if (!isdiagonal())
    throw Failed("Diagonal operators only valid for diagonal blocks");
  if (!*this) 
    throw Failed(NOPROPS);
  transform_mh(detect_table,transform_homo_(),row().Us,sigma0det);
  sigma_table.clear();
  reset();
}

template<class T> void 
BasePeriodicSpectrum<T>::set_Us(const BaseList<cmatrix>& Us,double periodv)
{
  set_steps(Us.size());
  this->period(periodv);
  setdiagonal()=Us;
}

template<class T> void
BasePeriodicSpectrum<T>::set_Us(const BaseList<complex>& Us,double periodv)
{
  set_steps(Us.size());
  this->period(periodv);
  setdiagonal()=Us;
}

template<class T> void 
BasePeriodicSpectrum<T>::set_Us(char sel,const BaseList<cmatrix>& Us,double periodv)
{
  set_steps(Us.size());
  this->period(periodv);
  setRC(sel)=Us;
}

template<class T> void
BasePeriodicSpectrum<T>::set_Us(char sel,const BaseList<complex>& Us,double periodv)
{
  set_steps(Us.size());
  this->period(periodv);
  setRC(sel)=Us;
}

void HomoStash::operator= (const BaseList<complex>& Usv)
{
  Us.dimensions(1);
  List<complex>& UT=Us.list();
  const size_t n=Usv.size();
  UT.create(n);
  UT(0)=1.0;
  for (size_t gn=1;gn<n;gn++) UT(gn)=Usv(gn-1);
  eff_freq.create(1,diag_propagator(Usv.back(),1.0));
}

void HomoStash::operator= (const BaseList<cmatrix>& Usv)
{
  const size_t n=Usv.size();
  const cmatrix& Utop=Usv(n-1);
  const size_t dim=Utop.rows();
  switch (dim) {
    case 0: throw Undefined("Homogeneous propagation: propagators undefined");
  case 1: { //handle special case of single element propagators for convenience
    if (Utop.cols()!=1)
      throw NotSquare("Homogeneous propagation: invalid propagator");
    Us.dimensions(1);
    List<complex>& UT=Us.list();
    UT.create(n);
    UT(0)=1.0;
    for (size_t gn=1;gn<n;gn++) {
      const cmatrix& Utmp = Usv(gn-1);
      if ((Utmp.rows()!=1) || (Utmp.cols()!=1))
	throw Mismatch("Homogeneous propagation: propagators vary in dimension!");
      UT(gn)=Utmp(0,0);
    }
    eff_freq.create(1,diag_propagator(Utop(0,0),1.0));
  }
    return;
  }
  ScratchList<complex> ceigs(dim);
  MultiMatrix<complex,3>& UT=Us.create(n,dim,dim);
  cmatrix D(UT(0));
  if (eigensystem(D,ceigs,Utop,offdiagonal_tolerance)) {
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
  diag_propagator(eff_freq,ceigs,1.0);
}
   
template<class T> void BasePeriodicSpectrum<T>::set_steps(int stepsv)
{
  if (PeriodicObj::set_steps(stepsv)) {
    const size_t ftsteps = usegammasteps_ ? stepsv : nobs;
    usefft=ispowerof2(ftsteps);        
    if (!usefft) {
      scr.create(ftsteps);
      ft_facs.create(ftsteps);
      ft_create_table(ft_facs,FT_FORWARD);
    }

    if (eigrestrict) {
      const int nobs2=nobs/2;
      
      if ((valsideband>nobs2) || (valsideband<-nobs2))
	throw InvalidParameter("BasePeriodicSpectrum: sideband out of range");
      valrestrict= (valsideband<0) ? valsideband+nobs : valsideband;
    }
  }
  reset();
}

bool PeriodicSpectrum::calcamps(size_t r,size_t s)
{
  const bool isel=iselement();
  const bool isED=!sigma_table;

  List<complex>& store=usefft ? sideband_ : scr;
  store.create(nobs);
  sideband_.create(nobs);
 
  diff=-frequency(r,s);

  const double phasestep=2*M_PI*this->period()/nobs;
  const complex pfac=expi(diff*phasestep);
  
  const size_t ind=whichrs();
  complex p;

  if (isED)
    p= conj(isel ? detect_table.matrix()(0,ind) : detect_table.multimatrix3()(0,s,r));
  else
    p= sigma_table.matrix()(r,s);
  p/=nobs;
 
  if (norm(p)<1e-12) {
    sideband_=0.0;
    return false;
  }
  
  if (isel) {
    const cmatrix& detect_matrix=detect_table.matrix();
    for (size_t j=0;j<nobs;j++) { 
      store(j)=p*detect_matrix(j,ind);
      p*=pfac;
    }
  }
  else {
    const MultiMatrix<complex,3>& detect_mm=detect_table.multimatrix3();
    for (size_t j=0;j<nobs;j++) {	 
      store(j)=p*detect_mm(j,s,r);
      p*=pfac;
    }
  }

  if (usefft)
    fft_ip(sideband_,FT_FORWARD);
  else
    ft(sideband_,scr,ft_facs);

  return true;
}

bool hasnz(const BaseList<complex>& a,double tol)
{
  tol*=tol;
  const size_t n=a.size();
  double s=0.0;
  for (size_t gn=0;gn<n;gn++) {
    s+=norm(a(gn));
    if (s>tol)
      return true;
  }
  return false;
}

bool GammaPeriodicSpectrum::calcamps(size_t r,size_t s)
{
  sideband_.create(nobs);

  const bool isel=iselement();
  const size_t ind=whichrs();
  if (nobs==1) {
    complex sum(0,0);
    if (isel) {
      const cmatrix& sigma_matrix(sigma_table.matrix());
      const cmatrix& detect_matrix(detect_table.matrix());

      for (size_t which=gamma_steps;which--;)
	mla(sum,sigma_matrix(which,ind),detect_matrix(which,ind));
    }
    else {
      const MultiMatrix<complex,3>& sigma_mm(sigma_table.multimatrix3());
      const MultiMatrix<complex,3>& detect_mm(detect_table.multimatrix3());
      
      for (size_t which=gamma_steps;which--;)
	mla(sum,sigma_mm(which,r,s),detect_mm(which,s,r));
    }
    return (norm(sideband_.front()=sum/double(gamma_steps))>1e-12);
  }

  const double gamma_scale=1.0/(nobs*gamma_steps);
  const size_t rat=gamma_steps/nobs;

  usefft=ispowerof2(nobs);
  List<complex> &store=usefft ? sideband_ : scr;
  store.create(nobs);

  diff=-frequency(r,s);

  const double phasestep=2*M_PI*this->period()/nobs;
  const complex pfac=expi(diff*phasestep);
  const complex R=expi(-2*M_PI*diff*this->period());

  if (verbose_)
    std::cout << "Evolution phase factor for " << r << ", " << s << " transition: " << pfac << '\n';

  complex p(gamma_scale);
      
  int ptr=0;
  for (size_t j=0;j<gamma_steps;j+=rat,ptr++) {	 
    complex sum(0,0);
    if (isel) {
      const cmatrix& sigma_matrix =sigma_table.matrix();
      const cmatrix& detect_matrix =detect_table.matrix();

      for (size_t which=j;which<gamma_steps;which++)
	mla(sum,sigma_matrix(which-j,ind),detect_matrix(which,ind));
      if (j) {
	const int offset=gamma_steps-j;
	complex s2(0,0);
	for (size_t which2=j;which2--;) 
	  mla(s2,sigma_matrix(which2+offset,ind),detect_matrix(which2,ind));
	mla(sum,R,s2);
      }
    }
    else {
      const MultiMatrix<complex,3>& sigma_mm=sigma_table.multimatrix3();
      const MultiMatrix<complex,3>& detect_mm=detect_table.multimatrix3();
      
      for (size_t which=j;which<gamma_steps;which++)
	mla(sum,sigma_mm(which-j,r,s),detect_mm(which,s,r));
      if (j) {
	const int offset=gamma_steps-j;
	complex s2(0,0);
	for (size_t which2=j;which2--;) 
	  mla(s2,sigma_mm(which2+offset,r,s),detect_mm(which2,s,r));
	mla(sum,R,s2);
      }
    }
    store(ptr)=sum*p; 
    p*=pfac;
  }

  if (verbose_)
    std::cout << "Raw factors: " << store << '\n';

  if (!hasnz(store,1e-7)) {
    sideband_=0.0;
    return false;
  }

  if (usefft)
    fft_ip(sideband_,FT_FORWARD);
  else
    ft(sideband_,scr,ft_facs);
  
  return true;
}

bool GammaPeriodicSpectrumED::calcamps(size_t r,size_t s)
{
  if (!!sigma_table)
    throw Failed("sigma_table shouldn't be defined!");
  const bool isel=iselement();
  const size_t ind=whichrs();

  sideband_.create(nobs);

  if (nobs==1) {
    double nsum=0.0;
    if (isel) {
      const cmatrix& rtable=detect_table.matrix();
      for (size_t j=0;j<gamma_steps;j++)
	nsum+=norm(rtable(j,ind));
    }
    else {
      const MultiMatrix<complex,3>& detect_mm=detect_table.multimatrix3();
      for (size_t j=0;j<gamma_steps;j++)
	nsum+=norm(detect_mm(j,r,s));
    }
    sideband_(0)=nsum/gamma_steps;
    return true;
  }

  const double gamma_scale=1.0/(gamma_steps*gamma_steps);

  List<complex> &store=usefft ? gft : scr;
  store.create(gamma_steps);
  diff=-frequency(r,s);

  const double phasestep=2*M_PI*this->period()/gamma_steps;

  if (isel) {
    const cmatrix &rtable=detect_table.matrix();
    
    const complex pfac=expi(diff*phasestep);
    if (verbose_>1)
      std::cout << "Evolution phase factor for " << r << ", " << s << " transition: " << pfac << '\n';

    complex phas(1.0);
    
    for (size_t j=0;j<gamma_steps;j++) {
      store(j)=conj_multiply(rtable(j,ind),phas);
      phas*=pfac;
    }
  }
  else {      
    const MultiMatrix<complex,3>& detect_mm=detect_table.multimatrix3();

    if (diff==0.0) {
      for (size_t j=gamma_steps;j--;)
	store(j)=conj(detect_mm(j,r,s));
    }
    else {
      const complex pfac=expi(diff*phasestep);
      if (verbose_>1)
	std::cout << "Evolution phase factor for " << r << ", " << s << " transition: " << pfac << '\n';

      complex phas(1.0);
      
      for (size_t j=0;j<gamma_steps;j++) {
	store(j)=conj_multiply(detect_mm(j,r,s),phas);
	phas*=pfac;
      }
    }
  }
  if (verbose_>1)
    std::cout << "Raw factors: " << store << '\n';

  if (!hasnz(store,1e-7)) {
    sideband_=0.0;
    return false;
  }

  if (usefft)
    fft_ip(gft,FT_FORWARD);
  else
    ft(gft,store,ft_facs);
  
  for (size_t n=nobs;n--;) {
    double amp=0.0;
    for (size_t k=n;k<gamma_steps;k+=nobs) amp+=norm(gft(k));
    sideband_(n)=amp*gamma_scale;
  }
  return true;
}

template<class T> bool BasePeriodicSpectrum<T>::advance()
{
  if (++r==dimR) {
    r=0;
    s++;
    if (s==dimC) {
      finished=true;
      return false;
    }
  }
  return true;
}

template bool BasePeriodicSpectrum<double>::advance();
template bool BasePeriodicSpectrum<complex>::advance();

void GammaPeriodicSpectrumED::observe(const cmatrix& sigma0det)
{
  if (!*this)
    throw Failed(NOPROPS);
  transform_mh(detect_table,transform_homo_(),row().Us,col().Us,sigma0det);
  reset();
}

void GammaPeriodicSpectrumED::observe(const rmatrix& sigma0det)
{
  if (!*this) 
    throw Failed(NOPROPS);
  transform_mh(detect_table,transform_homo_(),row().Us,col().Us,sigma0det);
  reset();
}

void GammaPeriodicSpectrumED::observe()
{
  if (!*this) 
    throw Failed(NOPROPS);
  transform_mh(detect_table,transform_homo_(),row().Us,col().Us);
  reset();
}

void GammaPeriodicSpectrumED::observe(const BaseList<double>& sigma0det)
{
  if (!isdiagonal())
    throw Failed("Diagonal operators only valid for diagonal blocks");
  if (!*this) 
    throw Failed(NOPROPS);
  transform_mh(detect_table,transform_homo_(),row().Us,sigma0det);
  reset();
}

std::ostream& operator<< (std::ostream& ostr,const PeriodicObj& a)
{
  ostr << "Cycle frequency: " << a.speed() << " Hz    (period: ";
  prettyprint_time(a.period(),ostr) << ")\n";
  ostr << "Observations: ";
  if (a.nobs==0)
    ostr << "<unset>\n";
  else
    ostr << a.nobs << '\n';
  ostr << "Gamma integration: ";
  if (a.isgamma) {    
    if (a.gammasteps()==0)
      return ostr << "<unset>\n";
    return ostr << a.gammasteps() << " steps\n";
  }
  return ostr << "none\n";
}

std::ostream& operator<< (std::ostream& ostr,const HomoStash& a)
{
  ostr << "Frequencies (/cycle frequency): " << a.eff_freq << '\n';
  ostr << "Transformation matrices: " << a.Us << '\n';
  return ostr;
}

std::ostream& operator<< (std::ostream& ostr,const HomoStashGammaTD& a)
{
  ostr << "Cycle evolution factors: " << a.cycle_eigs << '\n';
  ostr << "Transformation matrix: " << a.D << '\n';
  return ostr;
}

// void HomoStash::write(FILE* fp,const char* name) const
// {
//   ScratchList<char> tmpstr_loc(strlen(name)+12);
//   char *tmpstr=tmpstr_loc.vector();
//   sprintf(tmpstr,"%sFreqs",name);
//   WriteMATLAB(fp,eff_freq,tmpstr);
//   sprintf(tmpstr,"%sTransforms",name);
//   WriteMATLAB(fp,Us,tmpstr);
// }

// int HomoStash::read(FILE *fp)
// {
//   int fail=ReadMATLAB(eff_freq,fp);
//   return fail ? fail : ReadMATLAB(Us,fp);
// }

template class SwapStore<HomoStash>;

template<class T> std::ostream& operator<< (std::ostream& ostr, const BasePeriodicSpectrum<T>& a) 
{
  if (!a)
    return ostr << "<undefined>\n";

  ostr << static_cast<const PeriodicObj&>(a);
  ostr << "Block dimensions: " << a.rows() << " x " << a.cols() << '\n';
  ostr << "Block treated as vector: " << (a.iselement() ? "yes" : "no") << '\n';
  if (a.isdiagonal())
    ostr << "diagonal block: " << a.row() << '\n';
  else {
    ostr << "bra states: " << a.row() << '\n';
    ostr << "ket states: " << a.col() << '\n';
  }

  ostr << "Using FFT: " << (a.usefft ? "yes" : "no") << '\n';
  ostr << "Position: " << a.ncount << "/" << a.observations() << '\n';
  if (a.ncount) {
    ostr << a.r << "," << a.s << ": " << a.diff << " Hz" << '\n';
    ostr << "Current intensities: " << a.sideband_ << '\n';
  }
  if (a.verbose_>1) {
    ostr << "Sigma table:\n" << a.sigma_table << '\n';
    ostr << "Detect table:\n" << a.detect_table << '\n';
  }
  return ostr;
}

template<class T> bool
BasePeriodicSpectrum<T>::operator() (T& amp, double& freq)
{
  if (finished)
    return false;
  verifyclean();
  if (ncount==0) {
    while (!calcamps(r,s)) {
      if (!advance())
	return false;
    }
  }
  if (eigrestrict) {
    amp=sideband_(valrestrict);
    freq=diff+valrestrict*this->speed();
    advance();
  }
  else {
    amp=sideband_(ncount);
    freq=diff+(this->observations()-ncount)*this->speed();
    
    if (++ncount==this->observations()) {
      ncount=0;
      advance();
    }
  }
  return true;
}

template<class T,class U> inline void supermla(List<U>& amps,U scale,const BaseList<T>& sideband_)
{
  mla(amps,scale,sideband_);
}
template<class T> inline void supermla(List<complex>& amps,complex scale,const BaseList<T>& sideband_)
{
  if (imag(scale))
    mla(amps,scale,sideband_);
  else
    mla(amps,real(scale),sideband_);
}

template<class T,class U> inline void supermla(BaseList<U> amps,U scale,const BaseList<T>& sideband_)
{
  mla(amps,scale,sideband_);
}
template<class T> inline void supermla(BaseList<complex> amps,complex scale,const BaseList<T>& sideband_)
{
  if (imag(scale))
    mla(amps,scale,sideband_);
  else
    mla(amps,real(scale),sideband_);
}

   template<class T>
   template<class U>
   void BasePeriodicSpectrum<T>::add(BaseList<U> amps,U scale,size_t lr,size_t ls)
   {
     if (calcamps(lr,ls))
       supermla(amps,scale,sideband_);
   }

template void BasePeriodicSpectrum<double>::add(BaseList<double> amps,double,size_t,size_t);
template void BasePeriodicSpectrum<double>::add(BaseList<complex> amps,complex,size_t,size_t);
template void BasePeriodicSpectrum<complex>::add(BaseList<complex> amps,complex,size_t,size_t);

   template<class T>
   template<class U>
   void BasePeriodicSpectrum<T>::add(List<U>& amps,U scale,size_t lr,size_t ls)
   {
     if (calcamps(lr,ls))
       supermla(amps,scale,sideband_);
   }

template void BasePeriodicSpectrum<double>::add(List<double>& amps,double,size_t,size_t);
template void BasePeriodicSpectrum<double>::add(List<complex>& amps,complex,size_t,size_t);
template void BasePeriodicSpectrum<complex>::add(List<complex>& amps,complex,size_t,size_t);

template class BasePeriodicSpectrum<double>;
template class BasePeriodicSpectrum<complex>;

template std::ostream& operator<< (std::ostream&, const BasePeriodicSpectrum<double>&);
template std::ostream& operator<< (std::ostream&, const BasePeriodicSpectrum<complex>&);

}//namespace libcmatrix


