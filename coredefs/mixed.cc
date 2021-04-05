#ifdef LCM_SUPPRESS_VIEWS
#undef LCM_SUPPRESS_VIEWS
#endif

/* Operators that convert between rmatrices and cmatrices */

#include "cmatrix.h"
#include "rmatrix.h"
#include "List.h"

namespace libcmatrix {
 
template<class T> void
real_multiply(const BaseList<T>& a, complex b, const BaseList<complex>& c)
{
  const size_t n=c.length();
  if (a.length()!=n)
    throw Mismatch("real_multiply",a.length(),n);
  const complex *source=c.vector();
  T *dest=a.vector();
  for (size_t i=n;i--;)
    dest[i]=real_multiply(b,source[i]);
}

template<class T> void real_mla_(Matrix<T>& a,complex b,const cmatrix& c)
{
  if (!c) throw Undefined("real_mla");
  if (!a) {
    a.create(c.rows(),c.cols());
    real_multiply(a.row(),b,c.row());
  }
  else {
    if (!arematching(a,c))
      throw Mismatch("real_mla",a.rows(),a.cols(),c.rows(),c.cols());
    real_mla(a.row(),b,c.row());
  }
}

void real_mla(rmatrix& a,complex b,const cmatrix& c)
{
  real_mla_(a,b,c);
}

void real_mla(cmatrix& a,complex b,const cmatrix& c)
{
  real_mla_(a,b,c);
}

//R3 - only need to define "unnatural" orderings

cmatrix operator* (const Matrix<double>& a,complex b)
{
  cmatrix d(mxflag::temporary);
  emultiply(d,a,b);
  return d;
}

// cmatrix operator* (complex b,const Matrix<double>& a)
// {
//   cmatrix d(mxflag::temporary);
//   apply_sd(d,doesmultiply<complex,double,complex>(),a,b);
//   return d;
// }

// cmatrix operator+ (const Matrix<double>& a,complex b)
// {
//   cmatrix d(mxflag::temporary);
//   apply_sd(d,doesadd<complex,double,complex>(),a,b);
//   return d;
// }

cmatrix operator+ (complex b,const Matrix<double>& a)
{
  cmatrix d(mxflag::temporary);
  add(d,b,a);
  return d;
}

// cmatrix operator- (const Matrix<double>& a,complex b)
// {
//   cmatrix d(mxflag::temporary);
//   apply(d,doessubtract<complex,double,complex>(),a,b);
//   return d;
// }

cmatrix operator- (complex b,const Matrix<double>& a)
{
  cmatrix d(mxflag::temporary);
  subtract(d,b,a);
  return d;
}

cmatrix operator* (const cmatrix& a,double b)
{
  cmatrix d(mxflag::temporary);
  emultiply(d,a,b);
  return d;
}

void conj_multiply(List<complex_t> &a,const BaseList<complex_t> &b,const BaseList<complex_t> &c)
{
  a.create(b.size()); 
  conj_multiply( static_cast< BaseList<complex_t>& >(a),b,c);
}

void conj_multiply(List<complex_t> &a,const BaseList<complex_t> &b,const BaseList<double> &c)
{
  a.create(b.size());
  conj_multiply( static_cast< BaseList<complex_t>& >(a),b,c);
}


#ifdef LCM_COMPLEX_CHEAT

const BaseList<float_t> asdoubles(const BaseList<complex_t>& a)
{ return BaseList<float_t>(a.size()*2,(float_t*)(a.vector())); }

BaseList<float_t> asdoubles(BaseList<complex_t>& a)
{ return BaseList<float_t>(a.size()*2,(float_t*)(a.vector())); }

const BaseList<complex_t> ascomplex(const BaseList<float_t>& a) {
  if (a.size() % 2)
    throw InvalidParameter("ascomplex: input not of even length");
  else
    return BaseList<complex_t>(a.size()/2,(complex_t*)(a.vector()));
}
BaseList<complex_t> ascomplex(BaseList<float_t>& a) {
  if (a.size() % 2)
    throw InvalidParameter("ascomplex: input not of even length");
  else
    return BaseList<complex_t>(a.size()/2,(complex_t*)(a.vector()));
}

//Matrix<float_t> asdoubles(Matrix<complex_t>& a) { return Matrix<float_t>(a.cols()*2,a.rows(),(float_t*)(a.vector()),mxflag::nondynamic); }
//const Matrix<float_t> asdoubles(const Matrix<complex_t>& a) { return Matrix<float_t>(a.cols()*2,a.rows(),(float_t*)(a.vector()),mxflag::nondynamic); }
IndirectList<float_t,slice> imags(BaseList<complex_t>& a)
{ return asdoubles(a)(slice(1,a.size(),2)); }

const IndirectList<float_t,slice> imags(const BaseList<complex_t>& a)
{ return asdoubles(a)(slice(1,a.size(),2)); }

IndirectList<float_t,slice> reals(BaseList<complex_t>& a)
{ return asdoubles(a)(slice(0,a.size(),2)); }

const IndirectList<float_t,slice> reals(const BaseList<complex_t>& a)
{ return asdoubles(a)(slice(0,a.size(),2)); }

IndirectList<float_t,slice> imags(Matrix<complex_t>& a)
{ return asdoubles(a.row())(slice(1,a.size(),2)); }

const IndirectList<float_t,slice> imags(const Matrix<complex_t>& a)
{ return asdoubles(a.row())(slice(1,a.size(),2)); }

IndirectList<float_t,slice> reals(Matrix<complex_t>& a)
{ return asdoubles(a.row())(slice(0,a.size(),2)); }

const IndirectList<float_t,slice> reals(const Matrix<complex_t>& a)
{ return asdoubles(a.row())(slice(0,a.size(),2)); }

IndirectList<float_t,slice> real_diag(Matrix<complex_t>& a)
 { if (!issquare(a))
     throw NotSquare("real_diag");
   return asdoubles(a.row())(slice(0,a.rows(),2*(a.rows()+1))); }

const IndirectList<float_t,slice> real_diag(const Matrix<complex_t>& a) {
  if (!issquare(a))
    throw NotSquare("real_diag");
  return asdoubles(a.row())(slice(0,a.rows(),2*(a.rows()+1))); }

//COMPLEX_CHEAT
#endif


// cmatrix operator* (double b,const cmatrix& a)
// {
//   cmatrix d(mxflag::temporary);
//   apply2(d,doesmultiply<complex,double,complex>(),b,a);
//   return d;
// }

}//namespace libcmatrix
