#if !defined(lcm_external_h) && defined(LCM_USE_EXTERNAL)
#define lcm_external_h

#define DREAL_PTR double*

#if defined(HAVE_LIBACML)

#define complex singlecomplex
#define LCM_EXTERNAL_NAME ACML
namespace acml {
#include "acml.h"
}
#define LCM_EXT_PREFIX acml
#define LCM_UPLO 'U'
#define DCMPLX_PTR LCM_EXT_PREFIX::doublecomplex*
#define DCMPLX_CPTR LCM_EXT_PREFIX::doublecomplex*
#define DREAL_CPTR double*

#define USE_SUNPERFACML

#else 

#if defined(HAVE_LIBATLAS) || defined(HAVE_LIBMKL) || defined(HAVE_LIBOPENBLAS)
#define USE_ATLAS
extern "C" {
#ifdef HAVE_LIBMKL
#define LCM_EXTERNAL_NAME MKL
#include "mkl_cblas.h"
#else
#ifdef HAVE_LIBOPENBLAS
#define LCM_EXTERNAL_NAME OpenBLAS
#ifdef HAVE_CBLAS_H
#include "cblas.h"
#else
#ifdef HAVE_CBLAS64_H
#include "cblas64.h"
#else
#include "OpenBLAS/cblas.h"
#endif
#endif
#else
#define LCM_EXTERNAL_NAME ATLAS
#include "cblas.h"
#endif
#endif
}
#define LCM_ORDER CblasRowMajor
#define LCM_UPLO CblasUpper
#define DCMPLX_PTR void*
#define DCMPLX_CPTR const void*
#define DREAL_CPTR const double* 

/* #define real singlereal */
/* extern "C" { */
/* #include "cmatrix_f2c.h" */
/* #include "blaswrap.h" */
/* #include "fblaswr.h" */
/* #include "clapack.h" */
/* } */
/* #undef real */
/* #define CMATRIX_UPLO "U" */

#else
#error "Unknown external library"
#endif
#endif

#include "cmatrix_complex.h"
#include "Matrix.h"

namespace libcmatrix {
  inline DCMPLX_PTR lapack_pass(Matrix<complex>& a) { return reinterpret_cast<DCMPLX_PTR>(a.vector()); }
  inline DCMPLX_PTR lapack_pass(RawMatrix<complex>& a) { return reinterpret_cast<DCMPLX_PTR>(a.vector()); }
  inline DCMPLX_CPTR lapack_pass(const Matrix<complex>& a) { return (DCMPLX_CPTR)(a.vector()); }
  inline DCMPLX_CPTR lapack_pass(const RawMatrix<complex>& a) { return (DCMPLX_CPTR)(a.vector()); }
  inline DREAL_PTR lapack_pass(Matrix<double>& a) { return const_cast<DREAL_PTR>(a.vector()); }
  inline DREAL_CPTR lapack_pass(const Matrix<double>& a) { return const_cast<DREAL_CPTR>(a.vector()); }
  inline DREAL_PTR lapack_pass(RawMatrix<double>& a) { return const_cast<DREAL_PTR>(a.vector()); }
  inline DREAL_CPTR lapack_pass(const RawMatrix<double>& a) { return const_cast<DREAL_CPTR>(a.vector()); }
  inline DREAL_PTR lapack_pass(BaseList<double>& a) { return const_cast<DREAL_PTR>(a.vector()); }
  inline DREAL_CPTR lapack_pass(const BaseList<double>& a) { return const_cast<DREAL_CPTR>(a.vector()); }
  inline DCMPLX_PTR lapack_pass(BaseList<complex>& a) { return reinterpret_cast<DCMPLX_PTR>(a.vector()); }
  inline DCMPLX_CPTR lapack_pass(const BaseList<complex>& a) { return reinterpret_cast<DCMPLX_CPTR>(a.vector()); }
  inline DCMPLX_PTR lapack_pass(complex* a) { return reinterpret_cast<DCMPLX_PTR>(a); }
  inline DCMPLX_CPTR lapack_pass(const complex& a) { return (DCMPLX_CPTR)(&a); }

 template<typename T> void getrowcols(size_t& r, size_t& c,const T& a, bool trans)
   {
     if (trans) {
       r=a.cols();
       c=a.rows();
     }
     else {
       r=a.rows();
       c=a.cols();
     }
   }

}

#ifdef HAVE_LIBLAPACK
typedef libcmatrix::complex lcmcomplex; //!< need to rename to avoid clash with complex used by f2c
//!< flag that f2c declarations have been made
#define LCM_F2C_SETUP
extern "C" {
#define complex scomplex
#define doublecomplex dcomplex
#include "f2c.h"
#undef doublecomplex
#define doublecomplex lcmcomplex
  //! LAPACK declarations.  File is clapack.h but renamed to prevent clash with ATLAS clapack.h
#include "lapack.h"
#undef complex
#undef doublecomplex
}
#endif

#endif
