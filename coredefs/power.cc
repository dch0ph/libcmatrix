#include "cmatrix.h"
#include "rmatrix.h"
#include "ScratchList.h"

namespace libcmatrix {

cmatrix exp(const cmatrix& a,complex scale)
{
  const size_t n=a.rows();
  cmatrix vectors;
  ScratchList<complex> eigs(n);
  
  eigensystem(vectors,eigs,a);
  for (size_t i=n;i--;)
    eigs(i)=exp(scale*eigs(i));

  cmatrix d(mxflag::temporary);
  simtrans(d,eigs,vectors); // this must be a "real" similarity transform since the eigenvector matrix is not guarenteed to be unitary
  return d;
}

void hermitian_exp(cmatrix& d,const cmatrix& a,complex scale)
{
  const size_t n=a.rows();
  cmatrix vectors;
  ScratchList<double> reigs(n);
  ScratchList<complex> eigs(n);

  hermitian_eigensystem(vectors,reigs,a);
  if (real(scale)==0.0) { // special case for purely imaginary
    const double b=imag(scale);
    for (size_t i=n;i--;)
      eigs(i)=expi(b*reigs(i));
  }
  else {
    for (size_t i=n;i--;)
      eigs(i)=exp(scale*reigs(i));
  }

  unitary_simtrans(d,eigs,vectors);
  
//   cmatrix tmpv;
//   std::cout << "H\n" << a;
//   conj_transpose_multiply(tmpv,vectors,vectors);
//   std::cout << "Eigs: " << eigs << std::endl;
//   std::cout << tmpv << std::endl;
//   conj_transpose_multiply(tmpv,d,d);
//   std::cout << tmpv << std::endl;
}

void hermitian_expi(cmatrix& d,const rmatrix& a,double scale)
{
  const size_t n=a.rows();
  rmatrix vectors;
  ScratchList<double> eigs(n);
  ScratchList<complex> ceigs(n);

  hermitian_eigensystem(vectors,eigs,a);
  for (size_t i=n;i--;)
    ceigs(i)=expi(scale*eigs(i));
  unitary_simtrans(d,ceigs,vectors);
}
        
  namespace {
    template<class T> void pow_(Matrix<T>& d, const Matrix<T>& a, int n)
    {
      if (!issquare(a))
	throw NotSquare("pow");
      if (n<0)
	throw InvalidParameter("pow: exponent cannot be <0");

      switch (n) {
      case 0:
	d.identity(a.rows());
	return;
      case 1:
	d=a;
	return;
      case 2:
	multiply(d,a,a);
	return;
      }
      Matrix<T>* lastcachep(const_cast<Matrix<T>* >(&a)); //!< harmless const cast
      Matrix<T> cache,tmp;
      bool done=false;
      for (;;) {
	if (n & 1) {
	  if (done) {
	    multiply(tmp,d,*lastcachep);
	    d.swap(tmp);	    
	  }
	  else {
	    d=*lastcachep;
	    done=true;
	  }
	}
	n>>=1;
	if (n) {
	  multiply(tmp,*lastcachep,*lastcachep);
	  lastcachep=&cache;
	  lastcachep->swap(tmp);
	}
	else {
	  assert(done);
	  return;
	}
      }
    }
	
  }

	
 void pow(cmatrix& d,const cmatrix& a,int n)
 {
   pow_(d,a,n);
 }

cmatrix pow(const cmatrix& a,double rz)
{
  const size_t n=a.rows();
  cmatrix vectors;
  ScratchList<complex> eigs(n);
  
  eigensystem(vectors,eigs,a);
  for (size_t i=n;i--;)
    eigs(i)=pow(eigs(i),rz);

  cmatrix d(mxflag::temporary);
  simtrans(d,eigs,vectors);
  return d;
}

void hermitian_pow(cmatrix& d,const cmatrix& a,double rz)
{
  const size_t n=a.rows();
  cmatrix vectors;
  ScratchList<double> reigs(n);
  
  hermitian_eigensystem(vectors,reigs,a);
  for (size_t i=n;i--;)
    reigs(i)=std::pow(reigs(i),rz);
  unitary_simtrans(d,reigs,vectors);
}

cmatrix hermitian_pow(const cmatrix& a,double z)
{
  cmatrix d(mxflag::temporary);
  hermitian_pow(d,a,z);
  return d;
}

cmatrix hermitian_exp(const cmatrix& a,complex scale)
{
  cmatrix d(mxflag::temporary);
  hermitian_exp(d,a,scale);
  return d;
}

cmatrix hermitian_expi(const rmatrix& a,double scale)
{
  cmatrix d(mxflag::temporary);
  hermitian_expi(d,a,scale);
  return d;
}

}
