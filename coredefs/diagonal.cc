#include "rmatrix.h"
#include "cmatrix.h"
#include "ScratchList.h"

namespace libcmatrix {

void full(cmatrix &a,const BaseList<float_t> &v)
{
  size_t n=v.length();
  a.create(n,n);
  a=0.0;
  for (;n--;) a(n,n)=v(n);
}

double norm(const BaseList<complex> &a)
{
  double ret=0.0;
  for (size_t i=a.length();i--;) ret+=norm(a(i));
  return ret;
}

inline double square(double a) { return a*a; }

double norm(const BaseList<float_t> &a)
{
  double ret=0.0;
  for (size_t i=a.length();i--;) ret+=square(a(i));
  return ret;
}

template<class T> void real_mla_(BaseList<T>& a, complex b, const BaseList<complex>& c)
{
  const size_t n=c.length();
  if (a.length()!=n)
    throw Mismatch("real_mla",a.length(),n);
  const complex *source=c.vector();
  T *dest=a.vector();
  for (size_t i=n;i--;)
    real_mla(dest[i],b,source[i]);
}

void real_mla(BaseList<float_t> a, complex b, const BaseList<complex>& c)
{
  real_mla_(a,b,c);
}

void real_mla(BaseList<complex> a, complex b, const BaseList<complex>& c)
{
  real_mla_(a,b,c);
}

complex trace(const BaseList<complex>& b, const cmatrix& c)
{
  if (!issquare(c))
    throw NotSquare("trace");
  const size_t n=b.length();
  if (n!=c.rows())
    throw Mismatch("trace");

  complex sum(0.0,0.0);
  for (size_t i=n;i--;) mla(sum,b(i),c(i,i));

  return sum;
}

template<class T1,class T2> inline void _mla(cmatrix& a, T1 scale, const BaseList<T2>& b)
{
  if (!a) {
    full(a,b);
    a*=scale;
  }
  const size_t n=a.rows();
  if (a.cols()!=n)
    throw NotSquare("mla");
  if (b.length()!=n)
    throw Mismatch("mla");
  
  for (size_t i=n;i--;) mla(a(i,i),scale,b(i));
}

void mla(cmatrix &a,double scale,const BaseList<float_t> &b)
{
  _mla(a,scale,b);
}

void mla(cmatrix &a,complex scale,const BaseList<float_t> &b)
{
  _mla(a,scale,b);
}

static const char CM[]="conj_multiply";

template<class T> inline void _conj_multiply(cmatrix& dest,const BaseList<complex>& b,const Matrix<T>& a)
{
  if (!a)
    throw Undefined(CM);

  const size_t m=b.length();
  const size_t n=a.cols();
  if (n!=a.rows())
    throw Mismatch(CM);

  dest.create(m,n);

  complex *destp=dest.vector();
  const T* source=a.vector();

  if (issame(destp,source)) {
    for (size_t i=0;i<m;i++) {
      const complex scale=conj(b(i));
      for (size_t j=n;j--;)
	(*destp++)*=scale;
    }
  }
  else {
    for (size_t i=0;i<m;i++) {
      const complex scale=conj(b(i));
      for (size_t j=n;j--;)
	*destp++=(*source++)*scale;
    }
  }
}

void conj_multiply(cmatrix& dest,const BaseList<complex>& b,const cmatrix& a)
{
  _conj_multiply(dest,b,a);
}

void conj_multiply(cmatrix& dest,const BaseList<complex>& b,const rmatrix& a)
{
  _conj_multiply(dest,b,a);
}

inline void multiply_conj_ip(complex &z1,const complex &z2)
//{ const double r=z1.real(); z1.re=r*z2.re+z1.im*z2.im; z1.im=z1.im*z2.re-r*z2.im; }
{ z1=complex(z1.real()*z2.real()+z1.imag()*z2.imag(), z1.imag()*z2.real()-z1.real()*z2.imag()); }

static const char MC[]="multiply_conj";

inline void multiply_conj_ip(double,const complex&) { throw Failed("Internal error"); }

template<class T> inline void _multiply_conj(cmatrix& dest,const Matrix<T>& a,const BaseList<complex>& b)
{
  if (!a)
    throw Undefined(MC);

  const size_t n=b.length();
  if (n!=a.cols())
    throw Mismatch(MC);

  const size_t m=a.rows();

  dest.create(m,n);
  complex *destp=dest.vector();
  const T* source=a.vector();
  const complex *svals=b.vector();

  if (issame(destp,source)) {
    for (size_t i=m;i--;) {
      for (size_t j=0;j<n;j++)
	multiply_conj_ip(*destp++,svals[j]);
    }
  }
  else {
    for (size_t i=m;i--;) {
      for (size_t j=0;j<n;j++)
	*destp++=conj_multiply(svals[j],*source++);
    }
  }
}

void multiply_conj(cmatrix& dest,const cmatrix& a,const BaseList<complex>& b)
{
  _multiply_conj(dest,a,b);
}

void multiply_conj(cmatrix& dest,const rmatrix& a,const BaseList<complex>& b)
{
  _multiply_conj(dest,a,b);
}


}//namespace libcmatrix
