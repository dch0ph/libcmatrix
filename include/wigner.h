#ifndef _wigner_h_
#define _wigner_h_

#include "cmatrix_complex.h"
#include "Euler.h"
#include "Matrix.h"

namespace libcmatrix {

  // reduced Wigner matrix elements
  double d(int l,int m,int n,double beta); // d(1)_m,n(beta)
  double d1(int m,int l,double beta); // special case d(1)
  double d2(int m,int l,double beta); // special case d(2)
  double legendre(int l,double); // special case d(l)_0,0(beta)
  
  Matrix<complex> D(int l,const Euler&);
  complex D(int l,int m,int n,const Euler&);

  inline complex wigner_factor(int n,int m,double alpha,double gamma) { return expi(-alpha*n-gamma*m); }

}

#endif
