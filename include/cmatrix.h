/* Operations particular to complex matrices */

#ifndef _cmatrix_h_
#define _cmatrix_h_

#include "cmatrix_complex.h"
#include "Matrix.h"
#include "List.h"
#include <cstdio>
#include "Warnings.h"

namespace libcmatrix {
  typedef COMPLEX_T complex_t;
  typedef Matrix<complex_t> cmatrix;

  cmatrix operator* (const cmatrix&, double);
  //cmatrix operator* (double ,const cmatrix&);
  cmatrix operator* (const Matrix<double>&,complex);
  //cmatrix operator* (complex,const Matrix<double>&);

  //cmatrix operator+ (const cmatrix&, double);
  cmatrix operator+ (double, const cmatrix&);
  //cmatrix operator- (const cmatrix&, double);
  cmatrix operator- (double, const cmatrix&);
  inline cmatrix operator/ (const cmatrix& a, double b) { return a*(1.0/b); }

bool hasoffdiagonal(const cmatrix&, double =0.0);
double imagfrac(const cmatrix&);

 extern Warning<> invalidmaxiter_warning; //!< failed to parse LCM_MAXITERATIONS_PER_EIGENVALUE
  
  struct eigensystem_controller {
    double tolerance;
    double effectivezero;
    bool throwexception;
    int verbose;
    size_t chebyshev_iterations;
    static const size_t max_chebyshev_iterations;
    bool useexternal;
    Warning<> nonorthogonal_warning;

    eigensystem_controller();
  };

  extern eigensystem_controller cmatrix_eigensystem_controller;

  void eigensystem(cmatrix&, List<complex_t>&, const cmatrix&, cmatrix* =NULL);
  void eigensystem(cmatrix&, BaseList<complex_t>, const cmatrix&, cmatrix* =NULL);
//returns true if matrix is diagonal with tol parameter
  bool eigensystem(cmatrix&, BaseList<complex_t>, const cmatrix&, double tol, cmatrix* =NULL);
  void hermitian_eigensystem(cmatrix&, List<double>&, const cmatrix&, cmatrix* =NULL);
  void hermitian_eigensystem(cmatrix&, BaseList<double>, const cmatrix&, cmatrix* =NULL);

 bool ishermitian(const cmatrix&, double tol);
 bool isunitary(const cmatrix&, const char* ="U");
 bool isunitary(const Matrix<double>&, const char* ="V");

void inv(cmatrix &,const cmatrix&);
cmatrix inv(const cmatrix&);
 void multiply_inv_ip2(cmatrix& a,cmatrix& eqns);
 void multiply_inv_ip2(Matrix<double>& a,cmatrix& eqns);
 void multiply_inv_ip2(cmatrix& a,BaseList<complex>& eqns);
 void multiply_inv_ip2(Matrix<double>& a,BaseList<complex>& eqns);
 List<complex> multiply_inv(const cmatrix& a,const BaseList<complex>& eqns);
 List<complex> multiply_inv(const Matrix<double>& a,const BaseList<complex>& eqns);

  void conj_transpose(cmatrix &,const cmatrix &);
  cmatrix conj_transpose(const cmatrix&);

void conj_mla(BaseList<complex_t>, const BaseList<complex_t>&, const BaseList<complex_t>&);
void conj_mla(List<complex_t>&, const BaseList<complex_t>&, const BaseList<complex_t>&);

void real_mla(Matrix<double>&, complex,const cmatrix&);
void real_mla(cmatrix&, complex,const cmatrix&);
 double real_trace(const cmatrix&);

  //! default number of Pade iterations
#define LCM_DEFAULT_PADE_ITERATIONS 6U
  //! raw matrix exponential via Pade approximants
  void exppade_ip(Matrix<double>&, size_t =LCM_DEFAULT_PADE_ITERATIONS);
  void exppade_ip(Matrix<complex>&, size_t =LCM_DEFAULT_PADE_ITERATIONS);

  extern size_t lcm_pade_iterations; //!< value used by higher level functions

void simtrans(cmatrix&, const cmatrix&, const cmatrix&);
void unitary_simtrans(cmatrix&, const cmatrix&, const cmatrix&, cmatrix* =NULL);
void unitary_simtrans(cmatrix&, const cmatrix&, const cmatrix&, const cmatrix&, cmatrix* =NULL);
void unitary_simtrans(cmatrix&, const cmatrix&, const BaseList<complex_t>&);
template<> inline void Matrix<complex_t>::unitary_simtrans(const cmatrix& b,cmatrix* tmp) { ::libcmatrix::unitary_simtrans(*this,*this,b,tmp); }
template<> inline void Matrix<complex_t>::unitary_simtrans(const BaseList<complex_t>& b) { ::libcmatrix::unitary_simtrans(*this,*this,b); }

void isimtrans(cmatrix&, const cmatrix&, const cmatrix&);
void unitary_isimtrans(cmatrix&, const cmatrix&, const cmatrix&, cmatrix* =NULL);
void unitary_isimtrans(cmatrix&, const cmatrix&, const BaseList<complex_t>&);
template<> inline void Matrix<complex_t>::unitary_isimtrans(const cmatrix& b, cmatrix* tmp) { ::libcmatrix::unitary_isimtrans(*this,*this,b,tmp); }
template<> inline void Matrix<complex_t>::unitary_isimtrans(const BaseList<complex_t>& b) { ::libcmatrix::unitary_simtrans(*this,*this,b); }  
// Diagonal matrices
void simtrans(cmatrix &,const BaseList<complex_t> &,const cmatrix &);
void isimtrans(cmatrix &,const BaseList<complex_t> &,const cmatrix &);

void unitary_simtrans(cmatrix&, const BaseList<double>&, const cmatrix&, cmatrix* =NULL);
void unitary_simtrans(cmatrix&, const BaseList<complex_t>&, const cmatrix&, cmatrix* =NULL);
void unitary_simtrans(cmatrix&, const Matrix<float_t>&, const cmatrix&, cmatrix* =NULL);
void unitary_simtrans(cmatrix&, const cmatrix&, const Matrix<float_t>&, cmatrix* =NULL);
void unitary_simtrans(cmatrix&, const BaseList<complex_t>&, const Matrix<double>&, cmatrix* =NULL);

void unitary_isimtrans(cmatrix&, const BaseList<double>&, const cmatrix&, cmatrix* =NULL);
void unitary_isimtrans(cmatrix&, const BaseList<complex_t>&, const cmatrix&, cmatrix* =NULL);
void unitary_isimtrans(cmatrix&, const Matrix<float_t>&, const cmatrix&, cmatrix* =NULL);
void unitary_isimtrans(cmatrix&, const BaseList<complex_t>&, const Matrix<double>&, cmatrix* =NULL);
void unitary_isimtrans(cmatrix&, const cmatrix&, const Matrix<float_t>&, cmatrix* =NULL);
 
void cmatrix_CTM(cmatrix&, const cmatrix&, const cmatrix&);
void conj_transpose_multiply(cmatrix&, const cmatrix&, const Matrix<double>&);
void cmatrix_MCT(cmatrix&, const cmatrix&, const cmatrix&);
void conj_transpose_multiply(cmatrix&, const cmatrix&, const BaseList<double>&);
void multiply_conj_transpose(cmatrix&, const BaseList<double>&, const cmatrix&);
void conj_transpose_multiply(cmatrix&, const cmatrix&, const BaseList<complex_t>&);
void multiply_conj_transpose(cmatrix&, const BaseList<complex_t>&, const cmatrix&);

#ifdef LCM_USE_EXTERNAL
 inline void lapack_CTM(cmatrix& d, const cmatrix& a, const cmatrix &b)
   {
     d.create(a.cols(),b.cols());
     RawMatrix<complex> draw(d);
     lapack_multiply(draw,a,b,false,true,false);
   }
 inline void lapack_MCT(cmatrix& d, const cmatrix& a, const cmatrix &b)
   {
     d.create(a.rows(),b.rows());
     RawMatrix<complex> draw(d);
     lapack_multiply(draw,a,b,false,false,true);
   }
 
 inline void multiply_conj_transpose(cmatrix& d, const cmatrix& a, const cmatrix& b)
   {
     if (issimple(a,b,LCM_INTERNAL_ZMT))
       cmatrix_MCT(d,a,b);
     else
       lapack_MCT(d,a,b);
   }
 inline void conj_transpose_multiply(cmatrix& d, const cmatrix& a, const cmatrix& b)
   {
     if (issimple(a,b,LCM_INTERNAL_ZTM))
       cmatrix_CTM(d,a,b);
     else
       lapack_CTM(d,a,b);
   }
#else
 inline void multiply_conj_transpose(cmatrix& d, const cmatrix& a, const cmatrix& b)
   { cmatrix_MCT(d,a,b); }
 inline void conj_transpose_multiply(cmatrix& d, const cmatrix& a, const cmatrix& b)
   { cmatrix_CTM(d,a,b); }
#endif

void multiply_conj(cmatrix&, const cmatrix&, const BaseList<complex_t>&);
void multiply_conj(cmatrix&, const Matrix<double>&, const BaseList<complex_t>&);
void conj_multiply(cmatrix&, const BaseList<complex_t>&, const cmatrix&);
void conj_multiply(cmatrix&, const BaseList<complex_t>&, const Matrix<double>&);

void conj_multiply(BaseList<complex_t>, const BaseList<complex_t>&, const BaseList<complex_t>&);
void conj_multiply(BaseList<complex_t>, const BaseList<complex_t>&, const BaseList<double>&);

  void conj_multiply(List<complex_t>&, const BaseList<complex_t>&, const BaseList<complex_t>&);
  void conj_multiply(List<complex_t>&, const BaseList<complex_t>&, const BaseList<double>&);

void conj_emultiply(cmatrix&, const cmatrix&, const cmatrix&); // A(i,j)=conj(B(i,j)) x C(i,j)
void real_conj_emultiply(Matrix<double>&, const cmatrix&, const cmatrix&); // A(i,j)=real(conj(B(i,j)) x C(i,j))

 template<typename T1,typename T2,typename T3> inline void unitary_isimtransLR(Matrix<T1>& d, const Matrix<double>& VL, const Matrix<T2>& a, const Matrix<T3>& VR, Matrix< LCM_NEWTYPE(double,T2,true) >* tmp_ =NULL)
    {
      Matrix< LCM_NEWTYPE(double,T2,true) > to_;
      Matrix< LCM_NEWTYPE(double,T2,true) >& tmp = tmp_ ? (*tmp_) : to_;
      transpose_multiply(tmp,VL,a);
      multiply(d,tmp,VR);
    }

 template<typename T1,typename T2,typename T3> inline void unitary_isimtransLR(Matrix<T1>& d, const BaseList<double>& VL, const Matrix<T2>& a, const Matrix<T3>& VR, Matrix< LCM_NEWTYPE(double,T2,true) >* tmp_ =NULL)
    {
      Matrix< LCM_NEWTYPE(double,T2,true) > to_;
      Matrix< LCM_NEWTYPE(double,T2,true) >& tmp = tmp_ ? (*tmp_) : to_;
      transpose_multiply(tmp,VL,a);
      multiply(d,tmp,VR);
    }

  template<typename T1,typename T2> inline void unitary_isimtransLR(cmatrix& d, const cmatrix& VL, const Matrix<T1>& a, const Matrix<T2>& VR, Matrix< LCM_NEWTYPE(T1,T2,true) >* tmp_ =NULL)
  {
      Matrix< LCM_NEWTYPE(T1,T2,true) > to_;
      Matrix< LCM_NEWTYPE(T1,T2,true) >& tmp = tmp_ ? (*tmp_) : to_;
      multiply(tmp,a,VR);
      conj_transpose_multiply(d,VL,tmp);
    }

  template<typename T1,typename T2> inline void unitary_isimtransLR(cmatrix& d, const cmatrix& VL, const BaseList<T1>& a, const Matrix<T2>& VR, Matrix< LCM_NEWTYPE(T1,T2,true) >* tmp_ =NULL)
  {
      Matrix< LCM_NEWTYPE(T1,T2,true) > to_;
      Matrix< LCM_NEWTYPE(T1,T2,true) >& tmp = tmp_ ? (*tmp_) : to_;
      multiply(tmp,a,VR);
      conj_transpose_multiply(d,VL,tmp);
    }
 
/*   template<typename T> inline void unitary_isimtransLR(cmatrix& d, const cmatrix& VL, const T& a, const cmatrix& VR, cmatrix* tmp_ =NULL) */
/*     { */
/*       cmatrix to_; */
/*       cmatrix& tmp = tmp_ ? (*tmp_) : to_; */
/*       multiply(tmp,a,VR); */
/*       conj_transpose_multiply(d,VL,tmp); */
/*     } */

/*   template<typename T> inline void unitary_simtransLR(cmatrix& d, const cmatrix& VL, const T& a, const cmatrix& VR, cmatrix* tmp_ =NULL) */
/*     { */
/*       cmatrix to_; */
/*       cmatrix& tmp = tmp_ ? (*tmp_) : to_; */
/*       multiply(tmp,VL,a); */
/*       multiply_conj_transpose(d,tmp,VR); */
/*     } */

  template<typename T1,typename T2> inline void unitary_simtransLR(cmatrix& d, const Matrix<T2>& VL, const Matrix<T1>& a, const cmatrix& VR, cmatrix* tmp_ =NULL)
    {
      cmatrix to_;
      cmatrix& tmp(tmp_ ? (*tmp_) : to_);
      multiply_conj_transpose(tmp,a,VR);
      multiply(d,VL,tmp);
    }
 
  template<typename T1> inline void unitary_isimtransLR(BaseList<complex> d, const complex& VL, const BaseList<T1>& a, const cmatrix& VR)
    {
      multiply(d,a,VR);
      d*=conj(VL);
    }

  template<typename T1> inline void unitary_simtransLR(BaseList<complex> d, const complex& VL, const BaseList<T1>& a, const cmatrix& VR)
    {
      multiply_conj_transpose(d,a,VR);
      d*=VL;
    }

  template<typename T> inline void unitary_isimtransLR(BaseList<complex> d, const cmatrix& VL, const BaseList<T>& a, const complex& VR)
    {
      conj_transpose_multiply(d,VL,a);
      d*=VR;
    }

  template<typename T> inline void unitary_isimtransLR(cmatrix& d, const BaseList<complex>& VL, const Matrix<T>& a, const BaseList<complex>& VR)
    {
      conj_multiply(d,VL,a);
      d*=VR;
    }
  
  template<typename T> inline void unitary_isimtransLR(BaseList<complex> d, const complex& VL, const BaseList<T>& a, const BaseList<complex>& VR)
    {
      multiply(d,a,VR);
      d*=conj(VL);
    }

  template<typename T> inline void unitary_isimtransLR(BaseList<complex> d, const BaseList<complex>& VL, const BaseList<T>& a, const complex& VR)
    {
      conj_multiply(d,VL,a);
      d*=VR;
    }

  template<typename T> inline void unitary_simtransLR(BaseList<complex> d, const cmatrix& VL, const BaseList<T>& a, const complex& VR)
    {
      multiply(d,VL,a);
      d*=conj(VR);
    }

  template<typename T> inline void unitary_isimtransLR(complex& d, const complex& VL, const T& a, const complex& VR)
    {
      d=conj_multiply(VL,a)*VR;
    }

  template<typename T> inline void unitary_simtransLR(complex& d, const complex& VL, const T& a, const complex& VR)
    {
      d=VL*multiply_conj(a,VR);
    }

  template<typename T> inline void unitary_isimtransLR(cmatrix& d, const BaseList<complex>& VL, const Matrix<T>& a, const cmatrix& VR)
    {
      multiply(d,a,VR);
      conj_multiply(d,VL,d);
    }
  
  template<typename T> inline void unitary_isimtransLR(cmatrix& d, const cmatrix& VL, const Matrix<T>& a, const BaseList<complex>& VR) // should do biggest operation first in case input is real
    {
      conj_transpose_multiply(d,VL,a);
      d*=VR;
    }


// Row vectors
void conj_transpose_multiply(BaseList<complex_t>&, const cmatrix&, const BaseList<complex_t>&);
void conj_transpose_multiply(BaseList<complex_t>&, const cmatrix&, const BaseList<double>&);
void multiply_conj_transpose(BaseList<complex_t>&, const BaseList<complex_t>&, const cmatrix&);
void multiply_conj_transpose(BaseList<complex_t>&, const BaseList<double>&, const cmatrix&);

  cmatrix exp(const cmatrix &a,complex scale =complex(1.0)); // only new-object form for "inefficient" operations
  cmatrix hermitian_exp(const cmatrix &a,complex scale =complex(1.0));
  void hermitian_exp(cmatrix &,const cmatrix &,complex scale =complex(1.0));

cmatrix pow(const cmatrix& a,double);
  void pow(cmatrix&, const cmatrix&, int);
void hermitian_pow(cmatrix &,const cmatrix &,double);
cmatrix hermitian_pow(const cmatrix &,double);

cmatrix hermitian_expi(const Matrix<double>&, double);
void hermitian_expi(cmatrix&, const Matrix<double>&, double);

complex trace_multiply(const BaseList<complex_t>&, const cmatrix&);
inline complex trace_multiply(const cmatrix& a, const BaseList<complex_t>& b) { return trace_multiply(b,a); }
double hermitian_trace(const cmatrix&);
double hermitian_trace_multiply(const cmatrix&, const cmatrix&);
double hermitian_sum(const cmatrix&);
double hermitian_trace_multiply(const BaseList<double>&, const cmatrix&);
inline double hermitian_trace_multiply(const cmatrix& a,const BaseList<double>& b) { return hermitian_trace_multiply(b,a); }

inline Matrix<double> real(const cmatrix& a) { Matrix<double> d(mxflag::temporary); real(d,a); return d; }
inline Matrix<double> imag(const cmatrix& a) { Matrix<double> d(mxflag::temporary); imag(d,a); return d; }
inline Matrix<double> enorm(const cmatrix& a) { Matrix<double> d(mxflag::temporary); enorm(d,a); return d; }
inline Matrix<complex_t> conj(const cmatrix& a) { Matrix<complex_t> d(mxflag::temporary); conj(d,a); return d; }

  template<> struct doesconj_ip<cmatrix> { void operator()(cmatrix& v) const { v.conj(); } };

  cmatrix kronecker(const cmatrix&, const cmatrix&);
  void kronecker(cmatrix&, const cmatrix&, int);
  void kronecker(cmatrix&, int, const cmatrix&);
  void kronecker(cmatrix&, const cmatrix&, const cmatrix&);
  void kronecker_transpose(cmatrix&, int, const cmatrix&);
  void kronecker_transpose(cmatrix&, const cmatrix&, const cmatrix&);

  std::pair<double,double> hermitian_eigenvalue_range(const cmatrix&); //!< return estimate of maximum eigenvalue
  std::pair<double,double> eigenvalue_range(const Matrix<double>&); //!< return estimate of maximum eigenvalue
  
  //inline void spy(const cmatrix& a,double tol =1e-10) { spy(std::cout,a,tol); }
  //void spy(std::ostream&, const Matrix<bool>&);
  //inline void spy(const Matrix<bool>& a) { spy(std::cout,a); }

// I/O

void read_vector(List<complex_t>&, const char*);
void read_vector(List<complex_t>&, FILE*);
void read_vector_ascii(List<complex_t>&, FILE*);
void read_matrix(cmatrix&, const char*);
void read_matrix(cmatrix&, FILE*);

 void write_comment(FILE *, const char *c, const char *intro);
void write_vector(FILE *,const BaseList<complex_t> &,int =mxflag::doublep);
void write_vector(const char *,const BaseList<complex_t> &,int =mxflag::doublep);
void write_matrix(FILE *,const cmatrix &,const char * =NULL,int =(mxflag::doublep | mxflag::block));
void write_matrix(const char *,const cmatrix &,const char * =NULL,int =(mxflag::doublep | mxflag::block));

void cmatrix_new_thread();
void cmatrix_resize_stack(size_t);
void reset_stats();
void print_stats(std::ostream& =std::cout);
void copyascomplex(List<complex_t>&, const BaseList<double>&);

#ifdef LCM_COMPLEX_CHEAT
const BaseList<float_t> asdoubles(const BaseList<complex_t>&);
BaseList<float_t> asdoubles(BaseList<complex_t>&);
const BaseList<complex_t> ascomplex(const BaseList<float_t>&);
BaseList<complex_t> ascomplex(BaseList<float_t>&);
#ifndef LCM_SUPPRESS_VIEWS
IndirectList<float_t,slice> imags(BaseList<complex_t>&);
const IndirectList<float_t,slice> imags(const BaseList<complex_t>&);
IndirectList<float_t,slice> reals(BaseList<complex_t>&);
const IndirectList<float_t,slice> reals(const BaseList<complex_t>&);
IndirectList<float_t,slice> imags(Matrix<complex_t>&);
const IndirectList<float_t,slice> imags(const Matrix<complex_t>&);
IndirectList<float_t,slice> reals(Matrix<complex_t>&);
const IndirectList<float_t,slice> reals(const Matrix<complex_t>&);
IndirectList<float_t,slice> real_diag(Matrix<complex_t>&);
const IndirectList<float_t,slice> real_diag(const Matrix<complex_t>&);
#endif
#endif

} //namespace libcmatrix

#endif
