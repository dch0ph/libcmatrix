#include "rmatrix.h"
#include "cmatrix.h"

namespace libcmatrix {

  void unitary_simtrans(rmatrix& d,const BaseList<double>& D,const rmatrix& V, rmatrix* tmp)
{
  rmatrix to_;
  rmatrix& to(tmp ? (*tmp) : to_);
  multiply_transpose(to,D,V);
  multiply(d,V,to);
}

void unitary_simtrans(cmatrix& d,const BaseList<complex>& D,const rmatrix& V, cmatrix* tmp)
{
  cmatrix to_;
  cmatrix& to= tmp ? (*tmp) : to_;
  multiply_transpose(to,D,V);
  multiply(d,V,to);
}

  void unitary_isimtrans(rmatrix &d,const BaseList<double> &D,const rmatrix &V, rmatrix* tmp)
{
  rmatrix to_;
  rmatrix& to(tmp ? (*tmp) : to_);
  transpose_multiply(to,V,D);
  multiply(d,to,V);
}

void unitary_isimtrans(cmatrix& d, const BaseList<complex>& D, const rmatrix& V, cmatrix* tmp)
{
  cmatrix to_;
  cmatrix& to(tmp ? (*tmp) : to_);
  transpose_multiply(to,V,D);
  multiply(d,to,V);
}

void unitary_simtrans(rmatrix &dest,const rmatrix &a,const rmatrix &b, rmatrix* tmp)
{
  rmatrix d_;
  rmatrix& d(tmp ? (*tmp) : d_);
  multiply_transpose(d,a,b);
  multiply(dest,b,d);
}

void unitary_simtrans(rmatrix &dest,const rmatrix &VR,const rmatrix &a,const rmatrix &VC, rmatrix* tmp)
{
  rmatrix d_;
  rmatrix& d(tmp ? (*tmp) : d_);
  multiply_transpose(d,a,VC);
  multiply(dest,VR,d);
}

void unitary_isimtrans(rmatrix &dest,const rmatrix &a,const rmatrix &b, rmatrix* tmp)
{
  rmatrix d_;
  rmatrix& d(tmp ? (*tmp) : d_);
  transpose_multiply(d,b,a);
  multiply(dest,d,b);
}

void unitary_isimtrans(rmatrix &dest,const rmatrix &VR,const rmatrix &a,const rmatrix &VC, rmatrix* tmp)
{
  rmatrix d_;
  rmatrix& d(tmp ? (*tmp) : d_);
  transpose_multiply(d,VR,a);
  multiply(dest,d,VC);
}

}//namespace libcmatrix
