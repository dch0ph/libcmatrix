#ifndef _space_T_h_
#define _space_T_h_

#include "Tensor.h"
#include "wigner.h"
#include "cmatrix.h"
#include "rmatrix.h"

namespace libcmatrix {

typedef Tensor<complex> space_T;

complex A2(const rmatrix &,int,int);
double A2(double,double,double,int,int);
space_T A2(const rmatrix &);
space_T A2(const rmatrix &,int);
space_T A2(double,double,double);
space_T A2(double,double,double,int);

complex A1(double,double,double,int);
space_T A1(double,double,double);

 space_T rotate(const space_T&, const Euler&);
 space_T rotate(const space_T&, const cmatrix&);
 void rotate(space_T&, const space_T&, const cmatrix&);
 space_T rotate(const space_T&, const BaseList<cmatrix>&);
 void rotate(space_T&, const space_T&, const BaseList<cmatrix>&);
   
  bool areequal(const space_T&, const space_T&, double);

complex rotate(const space_T&, int l,int m,const Euler&);
complex rotate(const space_T&, int m,const cmatrix&);

  double sumzero(const space_T&);
  double sumzero_rotate(const space_T&, const Matrix<complex>&);

  void tensor_to_matrix(Matrix<double>&, const space_T&); //!< convert (symmetric) tensor to matrix representation

template<class T> void multiply(Tensor<T> &X,const space_T& a,const T &b)
{
  X.create(a.rank());

  for (int i=a.rank();i>=0;i--) {
    if (a.have_rank(i)) {
      X.ensure_rank(i);
      for (int j=-i;j<=i;j++) {
	X(i,j)=b;
	X(i,j)*=a(i,j);
      }
    }
  }
}

template<class T> void mla(Tensor<T> &X,const space_T& a,const T &b)
{
  const int r=X.rank();
  if (r<0) {
    multiply(X,a,b);
    return;
  }
  if (a.max_rank()>r) throw Mismatch();

  for (int i=0;i<=r;i++) {
    if (a.have_rank(i)) {
      X.ensure_rank(i);
      for (int j=-i;j<=i;j++) mla(X(i,j),a(i,j),b);
    }
  }
}

template<class T> void multiply(Tensor<T> &c,const T &a,const space_T& b) { multiply(c,b,a); }
template<class T> void mla(Tensor<T> &c,const T &a,const space_T& b) { mla(c,b,a); }

} //namespace libcmatrix

#endif
