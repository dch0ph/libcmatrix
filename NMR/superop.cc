#include "superop.h"

namespace libcmatrix {

  //  static const double MINUS_TWO_PI=-2.0*M_PI;

template<class T> void superop_multiply_(Matrix<T> &d,const Matrix<T> &a,const Matrix<T> &b)
{
  if (!issquare(a) || !issquare(b))
    throw NotSquare("superop_multiply");

  const size_t na=a.rows();
  const size_t nb=b.rows();

  if (na==nb*nb) { // a is super-op
    d.create(nb,nb);
    BaseList<T> dr(d.row());
    multiply(dr,a,b.row());
  }
  else {
    if (na*na==nb) { // b is super-op
      d.create(na,na);
      BaseList<T> dr(d.row());
      multiply(dr,a.row(),b);
    }
    else
      throw Mismatch("superop_multiply",na*na,nb);
  }
}

void superop_multiply(cmatrix &d,const cmatrix &a,const cmatrix &b)
{
  superop_multiply_(d,a,b);
}

void superop_multiply(rmatrix &d,const rmatrix &a,const rmatrix &b)
{
  superop_multiply_(d,a,b);
}

cmatrix superop_multiply(const cmatrix &a,const cmatrix &b)
{
  cmatrix d(mxflag::temporary);
  superop_multiply_(d,a,b);
  return d;
}

rmatrix superop_multiply(const rmatrix &a,const rmatrix &b)
{
  rmatrix d(mxflag::temporary); 
  superop_multiply_(d,a,b);
  return d;
}

void commutator(rmatrix& d, const rmatrix& a)
{
  if (!issquare(a))
    throw NotSquare("commutator");
  const int n=a.rows();
  kronecker(d,a,n);
  rmatrix tmp;
  kronecker_transpose(tmp,n,a);
  d-=tmp;
}

rmatrix commutator(const rmatrix &a)
{
  rmatrix d(mxflag::temporary);
  commutator(d,a);
  return d;
}

void commutator(cmatrix& d, const cmatrix &a)
{
  if (!issquare(a))
    throw NotSquare("commutator");
  const int n=a.rows();
  kronecker(d,a,n);
  cmatrix tmp;
  kronecker_transpose(tmp,n,a);
  d-=tmp;
}

template<typename T> void commutator_(BaseList<T> d, const BaseList<T>& a)
{
  const size_t n=a.size();
  if (d.size()!=n*n)
    throw Mismatch("commutator: destination and source don't match",d.size(),n*n);
  for (size_t i=n;i--;) {
    const T ai(a(i));    
    const size_t base=i*n;
    for (size_t j=n;j--;)
      d(base+j)=ai-a(j);    
  }
}

inline size_t square(size_t n) { return n*n; }

List<double> commutator(const BaseList<double>& a)
{
  List<double> d(square(a.size()),mxflag::temporary);
  commutator_(d,a);
  return d;
}

List<complex> commutator(const BaseList<complex>& a)
{
  List<complex> d(square(a.size()),mxflag::temporary);
  commutator_(d,a);
  return d;
}

void commutator(BaseList<double> d, const BaseList<double>& a)
{
  commutator_(d,a);
}

void commutator(BaseList<complex> d, const BaseList<complex>& a)
{
  commutator_(d,a);
}

cmatrix commutator(const cmatrix &a)
{
  cmatrix d(mxflag::temporary);
  commutator(d,a);
  return d;
}

template<class T> inline void double_commutator_(Matrix<T> &d,const Matrix<T> &a)
{
  if (!issquare(a))
    throw NotSquare("double_commutator");

  const int n=a.rows();
 
  Matrix<T> a2,tmp;
  multiply(a2,a,a);

  kronecker_transpose(d,a,a);
  d*=-2.0;
  kronecker(tmp,a2,n);
  d+=tmp;
  a2.transpose();
  kronecker(tmp,n,a2);
  d+=tmp;
}

cmatrix double_commutator(const cmatrix &a)
{
  cmatrix d(mxflag::temporary);
  double_commutator_(d,a);
  return d;
}

rmatrix double_commutator(const rmatrix &a)
{
  rmatrix d(mxflag::temporary);
  double_commutator_(d,a);
  return d;
}

//commutator for selected sub-Hilbert space

template<class T> void commutator_(Matrix<T>& dest, const Matrix<T>& H, const BaseList<size_t>& bras, const BaseList<size_t>& kets)
{
  if (!issquare(H))
    throw NotSquare("commutator");

  const size_t nbras=bras.length();
  const size_t nkets=kets.length();
  const size_t nL=nbras*nkets;

  dest.create(nL,nL,T(0));

  size_t r,c;

  for (r=0;r<nbras;r++) {
    const size_t rn=r*nkets;
    for (c=0;c<nbras;c++) {
      const size_t cn=c*nkets;
      const T Hval=H(bras(r),bras(c));
      for (size_t k=nkets;k--;)
	dest(rn+k,cn+k)=Hval;
    }
  }
  for (r=0;r<nkets;r++) {
    for (c=0;c<nkets;c++) {
      const T Hval=H(kets(c),kets(r));
      for (size_t k=nbras;k--;)
	dest(k*nkets+r,k*nkets+c)-=Hval;
    }
  }
}

cmatrix commutator(const cmatrix& a, const BaseList<size_t>& brasel, const BaseList<size_t>& ketsel)
{
  cmatrix d(mxflag::temporary);
  commutator_(d,a,brasel,ketsel);
  return d;
}

void commutator(cmatrix& d,const cmatrix& a, const BaseList<size_t>& brasel, const BaseList<size_t>& ketsel)
{
  commutator_(d,a,brasel,ketsel);
}

rmatrix commutator(const rmatrix& a, const BaseList<size_t>& brasel, const BaseList<size_t>& ketsel)
{
  rmatrix d(mxflag::temporary);
  commutator_(d,a,brasel,ketsel);
  return d;
}

void commutator(rmatrix& d,const rmatrix& a, const BaseList<size_t>& brasel, const BaseList<size_t>& ketsel)
{
  commutator_(d,a,brasel,ketsel);
}

  bool lcm_pade_use=false; //!< diagonalisation is default

  Warning<> partitioning_ignored_warning("Partitioning ignored as diagonalisation is being used",&lcm_base_warning,LCM_DEBUG_WARNING);

  cmatrix propagatorL(const cmatrix& L,double t, const matrix_partition* partp)
  {
    cmatrix dest(mxflag::temporary);
    propagatorL(dest,L,t,partp);
    return dest;
  }
  
  void propagatorL(cmatrix& dest, const cmatrix& L, double t, const matrix_partition* partp)
  {
    if (lcm_pade_use) 
      pade_propagatorL(dest,L,t,partp);
    else {
      if (partp)
	partitioning_ignored_warning.raise();
      dest=exp(L,complex(t));
    }
  }
    
}//namespace libcmatrix
