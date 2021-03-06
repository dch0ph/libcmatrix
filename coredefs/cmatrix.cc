#undef LCM_SUPPRESS_VIEWS
// flag that Level 2 (matrix) data types should be instantiated
#define LCM_LEVEL2_INSTANT_COMPLEX
#include <iomanip>
#include "cmatrix.h"
#include "List.h"
#ifdef LCM_USE_EXTERNAL
#include "cmatrix_external.h"
#endif

namespace libcmatrix {

  namespace {
    const complex one(1.0);
    const complex zero(0.0);

    int ensure_pword(std::ostream& ostr) {
      static bool made_pword=false;
      static int pword_ind;
      if (!made_pword) {
	pword_ind=ostr.xalloc();
	made_pword=true;
      }
      return pword_ind;
    }

    ostream_controller cmatrix_ostream_controller_;
  }

  ostream_controller::ostream_controller() :
    matrixprecision(LCM_DEFAULT_MATRIX_PRECISION),
    complexview(pair),
    timeprecision(LCM_DEFAULT_TIME_PRECISION),
    spytolerance(1e-10),
    indexbase(0)
  {}

  std::ostream& operator<< (std::ostream& ostr, const setmatrixprecision& a)
  {
    cmatrix_ostream_controller(ostr).matrixprecision=a.prec_;
    return ostr;
  }

  std::ostream& operator<< (std::ostream& ostr, const settimeprecision& a)
  {
    if (a.prec_==0)
      throw InvalidParameter("settimeprecision");
    cmatrix_ostream_controller(ostr).timeprecision=a.prec_;
    return ostr;
  }

  void cmatrix_ostream_controller(std::ostream& ostr, ostream_controller& a)
  {
    ostr.pword(ensure_pword(ostr))=&a;
  }

  ostream_controller& cmatrix_ostream_controller(std::ostream& ostr)
  {
    void* ap(ostr.pword(ensure_pword(ostr)));
    return ap ? *static_cast<ostream_controller*>(ap) : cmatrix_ostream_controller_;
  }
     
void print_cmatrix(std::ostream& ostr, const Matrix<complex>& a)
{
  if (!a) {
    ostr << "<empty> [" << a.rows() << " x " << a.cols() << "]\n";
    return;
  }
  ostream_controller& ctrl(cmatrix_ostream_controller(ostr));  
  if (ctrl.matrixprecision==0) {
    spy(ostr,a,ctrl.spytolerance);
    return;
  }
  int width;
  int oldprec=-1;
  if (ctrl.matrixprecision<0)
    width=ostr.precision();
  else {
    width=ctrl.matrixprecision;
    oldprec=ostr.precision(ctrl.matrixprecision);
  }
  // Output is unreadable without fixed format
  //const LCM_IOS::fmtflags oldflags=ostr.setf(LCM_IOS::floatfield,LCM_IOS::fixed);
  const bool isfixed=ctrl.complexfixedwidth();
  width+=2; //!< allow for . and sign
  int pad=2*(ostr.width()-width);
  if (pad<1)
    pad=1;

  for (size_t i=0;i<a.rows();i++) {
    for (size_t j=0;j<a.cols();j++) {
      ctrl.print(ostr,a(i,j),width);
      if (isfixed) {
	//      ostr << '(' << std::setw(width) << real(z) << ',' << std::setw(width) << imag(z) << ")";
	for (size_t k=pad;k--;)
	  ostr << ' ';
      }
      else
	ostr << " \t"; //!< can't meaningfully pad, tab and hope for the best
    }
    ostr << '\n';
  }
  //! reset stream - not strictly exception safe but who cares some flags are not reset?
  //  ostr.flags(oldflags);
  if (oldprec>=0)
    ostr.precision(oldprec);
}

#ifndef LCM_USE_EXTERNTEMPLATE
  template void Matrix<complex>::dump() const;
  template void BaseList<complex>::dump() const;
  template void BaseList<size_t>::dump() const;
#endif

  void conj_ip(BaseList<complex>& a)
  { 
    doesconj_ip<complex> obj;
    for (size_t n=a.size();n--;)
    obj(a(n));
  }

  template<> void Matrix<complex_t>::conj()
  {
    BaseList<complex> asrow(row());
    conj_ip(asrow);
  }
 
double hermitian_trace(const cmatrix& a)
{
  if (!issquare(a))
    throw NotSquare("trace");

  size_t i=a.rows()-1;  
  double sum=real(a(i,i));
  for (;i--;)
    sum+=real(a(i,i));
  return sum;
}

/* Simple maths operators for matrices */

void conj_transpose(cmatrix &d,const cmatrix &a)
{
  if (!a)
    throw Undefined("conj_transpose");
  if (!!d && issame(d,a))
    throw ArgumentClash("conj_transpose");

  size_t ar=a.rows();
  size_t ac=a.cols();
  d.create(ac,ar);

  for (;ar--;) {
    const complex *source=a.vector(ar);
    for (size_t c=ac;c--;)
      d(c,ar)=conj(source[c]);
  }
}

template<> void Matrix<complex>::conj_transpose()
{
  if (!issquare(*this))
    throw NotSquare("conj_transpose");
  Matrix<complex>& data=*this;

  complex tmp;
  doesconj<complex> conjobj;
  doesconj_ip<complex> conjipobj;
  for (size_t r=rows()-1;r>0;r--) {
    complex *vec=vector(r);
    conjipobj(vec[r]);
    for (size_t c=r;c--;) {
      tmp=conjobj(vec[c]);
      complex& datacr(data(c,r));
      vec[c]=conjobj(datacr);
      datacr=tmp;
    }
  }
}

#include "../coredefs/shared.cc"

void kronecker(cmatrix& d, const cmatrix& a, int n)
{
  kronecker_(d,a,n);
}

void kronecker(cmatrix& d, int n,const cmatrix &a)
{
  kronecker_(d,n,a);
}
	
cmatrix kronecker(const cmatrix& a, const cmatrix& b)
{
  cmatrix d(mxflag::temporary);
  kronecker_(d,a,b);
  return d;
}

void kronecker(cmatrix& d, const cmatrix& a, const cmatrix& b)
{
  kronecker_(d,a,b,false);
}

void kronecker_transpose(cmatrix& d, int n, const cmatrix& b)
{
  kronecker_(d,n,b,true);
}
	
void kronecker_transpose(cmatrix& d, const cmatrix& a, const cmatrix& b)
{
  kronecker_(d,a,b,true);
}
	
/* Functions that are specific to hermitian matrices */

double hermitian_sum(const cmatrix &a)
{
  if (!issquare(a))
    throw NotSquare("hermitian_sum");

  double sumd=0.0;
  double sumo=0.0;
  
  const size_t n=a.rows();

  for (size_t i=n;i--;) {
    const complex *v=a.vector(i)+i;
    sumd+=real(*v++);
    for (size_t j=i+1;j<n;j++) sumo+=real(*v++);
  }
  return sumd+2.0*sumo;
}

double hermitian_trace_multiply(const cmatrix& b,const cmatrix& c)
{
  if (!issquare(b) || !issquare(c))
    throw NotSquare("hermitian_trace_multiply");

  const size_t n=c.cols();
  if (b.cols()!=n)
    throw Mismatch("hermitian_trace_multiply");

  double sumd=0.0;
  double sumo=0.0;

  for (size_t i=0;i<n;i++) {
    const complex *bv=b.vector(i);

    sumd+=real(bv[i])*real(c(i,i));
    for (size_t j=i+1;j<n;j++)
      real_mla(sumo,bv[j],c(j,i));
  }

  return sumd+2.0*sumo;
}

double hermitian_trace_multiply(const BaseList<double>& b,const cmatrix& c)
{
  if (!issquare(c))
    throw NotSquare("hermitian_trace_multiply");
  const size_t n=b.size();
  if (n!=c.rows())
    throw Mismatch("hermitian_trace_multiply");

  double sum=0.0;

  for (size_t i=n;i--;)
    sum+=b(i)*real(c(i,i));

  return sum;
}

void spy(std::ostream& ostr, const BaseList<complex>& a, double tol)
{
  checksnonzero_tolerance<complex> nonzerocheck(tol);
  for (size_t c=0;c<a.size();c++) {
    const complex& val(a(c));
#ifdef HAVE_ISFINITE
    if (!std::isfinite(val.real()) || !std::isfinite(val.imag())) {
      ostr << '!';
      continue;
    }
#endif
    ostr << (nonzerocheck(val) ? 'X' : '.');
  }
}

void spy(std::ostream& ostr,const Matrix<complex>& a,double tol)
{
  if (!a) {
    print(a,ostr);
    return;
  }
  for (size_t r=0;r<a.rows();r++) {
    spy(ostr,a.row(r),tol);
    if (r==a.rows()-1)
      ostr << std::endl;
    else
      ostr << '\n';
  }
}

void spy(std::ostream& ostr,const Matrix<bool>& a)
{
  if (!a) {
    print(a,ostr);
    return;
  }
  for (size_t r=0;r<a.rows();r++) {
    const BaseList<bool> crow=a.row(r);
    for (size_t c=0;c<crow.size();c++)
      ostr << (crow(c) ? 'T' : '.');
    if (r==a.rows()-1)
      ostr << std::endl;
    else
      ostr << '\n';
  }
}

void conj_multiply(BaseList<complex> a,const BaseList<complex>& b,const BaseList<complex>& c)
{
  const size_t n=b.size();
  if ( (n!=c.size()) || (n!=a.size()))
    throw Mismatch("conj_multiply");
  const complex *bd=b.vector();
  const complex *cd=c.vector();
  complex *ad=a.vector();  
  for (size_t i=n;i--;)
    ad[i]=conj_multiply(bd[i],cd[i]);
}
 
void conj_multiply(BaseList<complex> a,const BaseList<complex>& b,const BaseList<double>& c)
{
  const size_t n=b.size();
  if ( (n!=c.size()) || (n!=a.size()))
    throw Mismatch("conj_multiply");
  const complex* bd=b.vector();
  const double* cd=c.vector();
  complex* ad=a.vector();  
  for (size_t i=n;i--;)
    ad[i]=conj_multiply(bd[i],cd[i]);
}
 
void real_conj_multiply(BaseList<double> a,const BaseList<complex>& b,const BaseList<complex>& c)
{
  const size_t n=b.size();
  if ( (n!=c.size()) || (n!=a.size()))
    throw Mismatch("conj_multiply");
  const complex* bd=b.vector();
  const complex* cd=c.vector();
  double* ad=a.vector();  
  for (size_t i=n;i--;)
    ad[i]=real_conj_multiply(bd[i],cd[i]);
}
 
void conj_emultiply(cmatrix& a,const cmatrix& b,const cmatrix& c)
{
  if (!arematching(b,c))
    throw Mismatch("conj_emultiply",b.rows(),b.cols(),c.rows(),c.cols());
  a.create(b.rows(),b.cols());
  conj_multiply(a.row(),b.row(),c.row());
}
 
void real_conj_emultiply(Matrix<double>& a,const cmatrix& b,const cmatrix& c)
{
  if (!arematching(b,c))
    throw Mismatch("real_conj_emultiply",b.rows(),b.cols(),c.rows(),c.cols());
  a.create(b.rows(),b.cols());
  real_conj_multiply(a.row(),b.row(),c.row());
}
 
void conj_mla(BaseList<complex> a,const BaseList<complex>& b,const BaseList<complex>& c)
{
  const size_t n=b.size();
  if (n!=c.size() || a.size()!=n)
    throw Mismatch("conj_mla");
  const complex *bd=b.vector();
  const complex *cd=c.vector();
  complex *ad=a.vector();
  for (size_t i=n;i--;)
    conj_mla(ad[i],bd[i],cd[i]);
}

void conj_mla(List<complex>& a,const BaseList<complex>& b,const BaseList<complex>& c)
{
  if (a.size())
    conj_mla( static_cast< BaseList<complex> >(a),b,c);
  else
    conj_multiply(a,b,c);
}

bool hasoffdiagonal(const cmatrix& U, double ptol)
{
  if (!issquare(U))
    throw NotSquare("isdiagonal");
  if (ptol<0)
    throw InvalidParameter("hasoffdiagonal");
  ptol*=ptol;

  for (size_t j=U.rows();j--;) {
    for (size_t i=j;i--;) {
      if ( (norm(U(i,j))>ptol) || (norm(U(j,i))>ptol)) return true;
    }
  }
  return false;
}

#ifdef LCM_COMPLEX_CHEAT

void copyascomplex(List<complex>& d,const BaseList<double>& s) {
  if (s.size() % 2)
    throw Failed("copyascomplex: odd number of data items");
  d.create(s.size()/2);
  d=(complex *)(s.vector());
}

#else

void copyascomplex(List<complex>& d,const BaseList<double>& s) {
  if (s.size() % 2)
    throw Failed("copyascomplex: odd number of data items");
  const size_t n=s.size()/2;
  d.create(n);
  size_t i,j;
  for (i=j=0;i<n;j+=2)
    d(i++)=complex(s(j),s(j+1));
}

#endif

#ifdef LCM_USE_EXTERNAL

void lapack_multiply(RawMatrix<complex>& d, const RawMatrix<complex>& a, const RawMatrix<complex>& b, bool acc, bool Atrans,bool Btrans)
{
  size_t ar,ac,br,bc;
  getrowcols(ar,ac,a,Atrans);
  getrowcols(br,bc,b,Btrans);

  if (ac!=br)
    throw Mismatch("lapack_multiply",ar,ac,br,bc);

#ifdef USE_SUNPERFACML
  zgemm(Atrans ? 'Y' : 'N',
	Btrans ? 'Y' : 'N',
	ar,bc,ac,
	lapack_pass(one),
	lapack_pass(a),a.step(),
	lapack_pass(b),b.step(),
	lapack_pass(acc ? one : zero),
	lapack_pass(d),d.step());
#else
#ifdef USE_ATLAS
  cblas_zgemm(LCM_ORDER,
	      Atrans ? CblasConjTrans : CblasNoTrans,
	      Btrans ? CblasConjTrans : CblasNoTrans,
	ar,bc,ac,
	lapack_pass(one),
	      lapack_pass(a),a.step(),
	      lapack_pass(b),b.step(),
	      lapack_pass(acc ? one : zero),
	      lapack_pass(d),d.step());
#else
#error "Unsupported external library"
//   zgemm_("N","N",&ar,&bc,&ac,
// 	 lapack_pass(&one),
// 	 lapack_pass(a),&iac,
// 	 lapack_pass(b),&ibc,
// 	 lapack_pass(&zero),
// 	 lapack_pass(d),&bc);  
#endif
#endif
}

//NB effect of transpose on is to return ba rather than ab

void lapack_multiply(BaseList<complex>& d, const Matrix<complex>& a, const BaseList<complex>& b, bool Atrans)
{
  size_t ar,ac;
  getrowcols(ar,ac,a,Atrans);
  const size_t iac=a.cols();

  if ( (ac!=b.length()) || (ar!=d.length()))
    throw Mismatch("lapack_multiply");

#ifdef USE_ATLAS
  cblas_zgemv(LCM_ORDER,
	      Atrans ? CblasTrans : CblasNoTrans,
	ar,ac,
	      lapack_pass(one),
	      lapack_pass(a),iac,
	lapack_pass(b),1,
	      lapack_pass(zero),
	lapack_pass(d),1);
#else
#ifdef USE_SUNPERFACML
  zgemv(Atrans ? 'Y' : 'N',
	ar,ac,
	      lapack_pass(one),
	lapack_pass(a),iac,
	lapack_pass(b),1,
	      lapack_pass(zero),
	lapack_pass(d),1);
#else
#error "Unsupported external library"
#endif
#endif
}

#endif //USE_EXTERNAL

void cmatrix_multiply(RawMatrix<complex>& dest, const RawMatrix<complex>& a, const RawMatrix<complex>& b, bool acc)
{
  //#if (LCM_COMPLEX_BYREF)
  if (a.rows()<LCM_NAIVE_BELOW)
    multiply_naive_ref_(dest,a,b,acc);
  else
    multiply_ref_(dest,a,b,acc);
// #else
//   if (a.rows()<LCM_NAIVE_BELOW)
//     multiply_naive_direct_(dest,a,b,acc);
//   else
//     multiply_direct_(dest,a,b,acc);
// #endif
}

template<class T> void lcm_matrix_ensure_(Matrix<T>& d, size_t r, size_t c, bool acc)
{
  if (!d) {
    if (acc)
      d.create(r,c,T(0.0));
    else
      d.create(r,c);
  }
}

void cmatrix_multiply(Matrix<complex>& d, const RawMatrix<complex>& a, const RawMatrix<complex>& b, bool acc)
{
  lcm_matrix_ensure_(d,a.rows(),b.cols(),acc);
  RawMatrix<complex> rawd(d);
  cmatrix_multiply(rawd,a,b,acc);
}

void cmatrix_multiply(Matrix<double>& d, const RawMatrix<double>& a, const RawMatrix<double>& b, bool acc)
{
  lcm_matrix_ensure_(d,a.rows(),b.cols(),acc);
  RawMatrix<double> rawd(d);
  cmatrix_multiply(rawd,a,b,acc);
}

template<class InputIter, class T, class Func> T accumulate_ip(InputIter first, InputIter last, T init, Func fobj)
{
  for (;first!=last;++first)
    fobj(init,*first);
  return init;
}

double real_trace(const cmatrix& a)
{
  const IndirectList<double,slice> aslice=real_diag(a);
  return accumulate_ip(aslice.begin(),aslice.end(),0.0,doesadd_ip<double,double>());
}

  bool ishermitian(const cmatrix& a, double tol) {
    if (!issquare(a))
      throw NotSquare("ishermitian");
    if (tol) {
      if (tol<0)
	throw InvalidParameter("ishermitian: tolerance must be >=0");
      const double toltol=tol*tol;
      for (size_t r=a.rows();r--;) {
	if (imag(a(r,r))>tol)
	  return false;
	for (size_t c=r;c--;) {
	  if (norm(a(r,c)-conj(a(c,r)))>=toltol)
	    return false;
	}
      }
    }
    else {
      for (size_t r=a.rows();r--;) {
	if (imag(a(r,r)))
	  return false;
	for (size_t c=r;c--;) {
	  if (a(r,c)!=conj(a(c,r)))
	    return false;
	}
      }
    }
    return true;
  }

cmatrix conj_transpose(const cmatrix& a)
{ 
  cmatrix d(mxflag::temporary); 
  conj_transpose(d,a);
  return d;
}

double
imagfrac(const BaseList<complex>& a)
{
  double isum=0.0;
  double nsum=0.0;
  for (size_t i=a.length();i--;) {
    const complex& v=a(i);
    isum+=norm(imag(v));
    nsum+=norm(v);
  }
  return ::std::sqrt(isum/nsum);
}

double
imagfrac(const cmatrix& a) { return imagfrac(a.row()); }

// Withdrawn in ABI 3.5.0 in favour of more general dot
// complex trace_multiply_conj(const BaseList<complex>& b,const BaseList<complex>& c)
// {
//   size_t n=b.length();
//   if (n!=c.length())
//     throw Mismatch("trace_multiply",n,c.length());
//   complex sum(0.0);
//   for (;n--;)
//     mla_conj(sum,b(n),c(n));
//   return sum;
// }

} //namespace libcmatrix
