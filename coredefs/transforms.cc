/* Similarity transformations for general matrices */

#include "cmatrix.h"
#include "rmatrix.h"

namespace libcmatrix {

const char TM[]="transpose_multiply";
const char MT[]="multiply_transpose";

static inline void three_multiply(cmatrix &dest,const cmatrix &A,const cmatrix&B,const cmatrix &C)
{
  cmatrix tmp;
  multiply(tmp,B,C);
  multiply(dest,A,tmp);
}

static inline void three_multiply(cmatrix &dest,const cmatrix &A,const BaseList<complex>&B,const cmatrix &C)
{
  cmatrix tmp;
  multiply(tmp,B,C);
  multiply(dest,A,tmp);
}

void simtrans(cmatrix &d,const cmatrix &a,const cmatrix &b)
{
  three_multiply(d,b,a,inv(b));
}

void isimtrans(cmatrix &d,const cmatrix &a,const cmatrix &b)
{
  three_multiply(d,inv(b),a,b);
}

void simtrans(cmatrix &d,const BaseList<complex> &a,const cmatrix &b)
{
  three_multiply(d,b,a,inv(b));
}

void isimtrans(cmatrix &d,const BaseList<complex> &a,const cmatrix &b)
{
  three_multiply(d,inv(b),a,b);
}

/* Similarity transformations for unitary matrices */

#define MCT "multiply_conj_transpose"

template<class T> void __multiply_conj_transpose(cmatrix &to,const BaseList<T> &D,const cmatrix &V)
{
  if (!!to && issame(to,V))
    throw ArgumentClash(MCT);

  const size_t n=D.length();
  if (V.cols()!=n)
    throw Mismatch(MCT);
  const size_t m=V.rows();

  to.create(n,m);
  complex *destp=to.vector();

  for (size_t i=0;i<n;i++) {
    const T& v=D(i);
    for (size_t k=0;k<m;k++)
      *destp++=conj_multiply(V(k,i),v);
  }
}

void multiply_conj_transpose(cmatrix &to,const BaseList<complex> &D,const cmatrix &V)
{
  __multiply_conj_transpose(to,D,V);
}

void multiply_conj_transpose(cmatrix &to,const BaseList<double> &D,const cmatrix &V)
{
  __multiply_conj_transpose(to,D,V);
}

void unitary_simtrans(cmatrix &d,const BaseList<double> &D,const cmatrix &V, cmatrix* tmp)
{
  cmatrix to_;
  cmatrix& to= tmp ? (*tmp) : to_;
  __multiply_conj_transpose(to,D,V);
  multiply(d,V,to);
}

void unitary_simtrans(cmatrix &d,const BaseList<complex> &D,const cmatrix &V, cmatrix* tmp)
{
  cmatrix to_;
  cmatrix& to= tmp ? (*tmp) : to_;
  __multiply_conj_transpose(to,D,V);
  multiply(d,V,to);
}

#define CTM "conj_transpose_multiply"

template<class T> void _conj_transpose_multiply(cmatrix& to, const cmatrix& V, const BaseList<T>& D)
{
  if (!!to && issame(to,V))
    throw ArgumentClash(CTM);

  const size_t n=D.length();
  if (V.rows()!=n)
    throw Mismatch(CTM);
  const size_t m=V.cols();

  to.create(m,n);

  for (size_t k=n;k--;) {
    const T& v=D(k);
    const complex *Vr=V.vector(k);
    for (size_t i=m;i--;)
      to(i,k)=conj_multiply(Vr[i],v);
  }
}

void conj_transpose_multiply(cmatrix &to,const cmatrix &V,const BaseList<double> &D)
{
  _conj_transpose_multiply(to,V,D);
}

void conj_transpose_multiply(cmatrix &to,const cmatrix &V,const BaseList<complex> &D)
{
  _conj_transpose_multiply(to,V,D);
}

void unitary_isimtrans(cmatrix& d,const BaseList<double>& D,const cmatrix& V, cmatrix* tmp)
{
  cmatrix to_;
  cmatrix& to= tmp ? (*tmp) : to_;
  _conj_transpose_multiply(to,V,D);
  multiply(d,to,V);
}

void unitary_isimtrans(cmatrix& d,const BaseList<complex>& D,const cmatrix& V, cmatrix* tmp)
{
  cmatrix to_;
  cmatrix& to= tmp ? (*tmp) : to_;
  _conj_transpose_multiply(to,V,D);
  multiply(d,to,V);
}

//inline void mla_conj(complex &z,const complex &a,double b) { z.re+=a.re*b; z.im+=a.im*b; }
inline void mla_conj(complex &z,const complex &a,double b) { mla(z,a,b); }

template<class T1,class T2> inline void _multiply_conj_transpose(cmatrix& to,const Matrix<T1>& a,const Matrix<T2>& V)
{
  if (!!to && (issame(to,a) || issame(to,V)))
    throw ArgumentClash(MCT);
  if (!V)
    throw Undefined(MCT);
  const size_t n=a.cols();
  if (n!=V.cols())
    throw Mismatch(MCT,a.rows(),a.cols(),V.rows(),V.cols());

  const size_t ra=a.rows();
  const size_t rt=V.rows();
  to.create(ra,rt);
  to=complex(0.0);

  size_t ib=0;
  while (ib<ra) {
    size_t fini=ib+LCM_BLOCK_FACTOR;
    if (fini>ra)
      fini=ra;
    
    size_t kb=0;
    while (kb<rt) {
      size_t fink=kb+LCM_BLOCK_FACTOR;
      if (fink>rt)
	fink=rt;

      size_t jb=0;
      while (jb<n) {
	size_t finj=jb+LCM_BLOCK_FACTOR;
	if (finj>n)
	  finj=n;

	complex* di=to.vector(ib);
	const T1* ai=a.vector(ib);
	for (size_t i=ib;i<fini;i++) {
	  const T2* Vk=V.vector(kb);
	  for (size_t k=kb;k<fink;k++) {
	    //#if (LCM_COMPLEX_BYREF)
	    complex& dik = di[k];
	    //#else
	    //complex dik=di[k];
	    //#endif
	    for (size_t j=jb;j<finj;j++)
	      mla_conj(dik,ai[j],Vk[j]);
	    //#if (LCM_COMPLEX_BYREF==0)
	    di[k]=dik;
	    //#endif
	    Vk+=n;
	  }
	  di+=rt;
	  ai+=n;
	}
	jb=finj;
      }
      kb=fink;
    }
    ib=fini;
  }
}

void multiply_conj_transpose(BaseList<complex>& to, const BaseList<complex>& a, const cmatrix& T)
{
  if (issame(to,a))
    throw ArgumentClash(MCT);
  const size_t n=a.length();
  if (n!=T.cols())
    throw Mismatch(MCT);

  const size_t rt=T.rows();
  if (to.length()!=rt)
    throw Mismatch(MCT);
  to=complex(0.0);

  const complex *ap=a.vector();
  complex *destp=to.vector();
  const complex *Tp=T.vector();

  for (size_t k=0;k<rt;k++) {
    complex &sum=*destp++;
    for (size_t j=n;j--;) 
      mla_conj(sum,ap[j],Tp[j]);
    Tp+=n;
  }
}

void cmatrix_MCT(cmatrix &to,const cmatrix &a,const cmatrix &T)
{
  _multiply_conj_transpose(to,a,T);
}

void multiply_transpose(cmatrix &to,const cmatrix &a,const rmatrix &T)
{
  _multiply_conj_transpose(to,a,T);
}

template<class M> void _conj_transpose_multiply(cmatrix& to,const cmatrix& T,const Matrix<M>& a)
{
  if (!T)
    throw Undefined(CTM);
  if (!!to && (issame(to,T) || issame(to,a)))
    throw ArgumentClash(CTM);

  const size_t n=a.rows();
  if (n!=T.rows())
    throw Mismatch(CTM,T.rows(),T.cols(),a.rows(),a.cols());

  const size_t ca=a.cols();
  const size_t ct=T.cols();
  to.create(ct,ca);
  to=complex(0.0);

  size_t jb=0;
  while (jb<n) {
    size_t finj=jb+LCM_BLOCK_FACTOR;
    if (finj>n)
      finj=n;
    
    size_t ib=0;
    while (ib<ct) {
      size_t fini=ib+LCM_BLOCK_FACTOR;
      if (fini>ct)
	fini=ct;
      
      size_t kb=0;
      while (kb<ca) {
	size_t fink=kb+LCM_BLOCK_FACTOR;
	if (fink>ca)
	  fink=ca;
	
	for (size_t j=jb;j<finj;j++) {
	  const M* aj=a.vector(j);
	  const complex* Tj=T.vector(j);
	  complex* desti=to.vector(ib);
	  for (size_t i=ib;i<fini;i++) {
	    const complex Tji=Tj[i];

	    for (size_t k=kb;k<fink;k++)
	      mla_conj(desti[k],aj[k],Tji);
	    desti+=ca;
	  }
	}
	kb=fink;
      }
      ib=fini;
    }
    jb=finj;
  }
}

void cmatrix_CTM(cmatrix& to, const cmatrix& T, const cmatrix& a)
{
  _conj_transpose_multiply(to,T,a);
}

void conj_transpose_multiply(cmatrix& to, const cmatrix& T, const rmatrix& a)
{
  _conj_transpose_multiply(to,T,a);
}

template<typename M> void _conj_transpose_multiply(BaseList<complex>& to, const cmatrix& T, const BaseList<M>& a)
{
  if (issame(to,a))
    throw ArgumentClash(CTM);
  const size_t n=a.length();
  if (n!=T.rows())
    throw Mismatch(CTM);

  const size_t ct=T.cols();
  if (ct!=to.length())
    throw Mismatch(CTM);
  to=complex(0.0);
  complex *destp=to.vector();

  for (size_t i=0;i<ct;i++) {
    complex& sum=*destp++;
    for (size_t j=n;j--;) 
      mla_conj(sum,a(j),T(j,i));
  }
}

void conj_transpose_multiply(BaseList<complex> &to,const cmatrix &T,const BaseList<complex> &a)
{
  _conj_transpose_multiply(to,T,a);
}

void conj_transpose_multiply(BaseList<complex> &to,const cmatrix &T,const BaseList<double> &a)
{
  _conj_transpose_multiply(to,T,a);
}

void unitary_simtrans(cmatrix& dest,const cmatrix& a,const cmatrix& b, cmatrix* tmp)
{
  cmatrix d_;
  cmatrix& d= tmp ? (*tmp) : d_;
  _multiply_conj_transpose(d,a,b);
  multiply(dest,b,d);
}

void unitary_simtrans(cmatrix& dest, const rmatrix& a, const cmatrix& b, cmatrix* tmp)
{
  cmatrix d_;
  cmatrix& d= tmp ? (*tmp) : d_;
  _multiply_conj_transpose(d,a,b);
  multiply(dest,b,d);
}

void unitary_simtrans(cmatrix& dest,const cmatrix& a,const rmatrix& b, cmatrix* tmp)
{
  cmatrix d_;
  cmatrix& d= tmp ? (*tmp) : d_;
  multiply(d,b,a);
  multiply_transpose(dest,d,b);
}

void unitary_isimtrans(cmatrix& dest, const cmatrix& a, const rmatrix& b, cmatrix* tmp)
{
  cmatrix d_;
  cmatrix& d= tmp ? (*tmp) : d_;
  transpose_multiply(d,b,a);
  multiply(dest,d,b);
}

void unitary_isimtrans(cmatrix& dest,const cmatrix& a,const cmatrix& b, cmatrix* tmp)
{
  cmatrix d_;
  cmatrix& d= tmp ? (*tmp) : d_;
  conj_transpose_multiply(d,b,a);
  multiply(dest,d,b);
}

void unitary_simtrans(cmatrix& dest,const cmatrix& a, const BaseList<complex>& D)
{
  multiply(dest,D,a);
  multiply_conj(dest,dest,D);
}

void unitary_isimtrans(cmatrix& dest,const rmatrix& a,const cmatrix& b, cmatrix* tmp)
{
  cmatrix d_;
  cmatrix& d= tmp ? (*tmp) : d_;
  conj_transpose_multiply(d,b,a);
  multiply(dest,d,b);
}

void unitary_isimtrans(cmatrix &dest,const cmatrix &a,const BaseList<complex> &D)
{
  multiply(dest,a,D);
  conj_multiply(dest,D,dest);
}

}//namespace libcmatrix
