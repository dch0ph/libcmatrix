#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include "basedefs.h" 
#include <cmath>
#include <iosfwd>
#include <functional>

//! if LCM_DISABLE_IMPICIT_COMPLEX set disable implicit construction of complex type from double
#ifndef LCM_IMPLICIT_COMPLEX
#ifdef LCM_DISABLE_IMPLICIT_COMPLEX
#define LCM_IMPLICIT_COMPLEX explicit
#else
#define LCM_IMPLICIT_COMPLEX
#endif
#endif

#ifdef HAVE_SINCOS
# ifdef HAVE_SUNMATH_H
#  include <sunmath.h>
# else
extern "C" void sincos(double, double*, double*); // emergency definition in case found library but not header!
# endif
namespace libcmatrix {
  inline void cmatrix_sincos(double th, double& s, double& c) { 
    sincos(th,&s,&c); 
  }
}
#else
namespace libcmatrix {
  inline void cmatrix_sincos(double th, double& s, double& c) {
    s=::std::sin(th);
    c=::std::cos(th);
  }
}
#endif

namespace libcmatrix {

void sinhcosh(double,double &,double &);
  complex expi(float_t);  // quick calculation of e^i theta

// #if (LCM_COMPLEX_BYREF==0)
// #define CMPLX_ARG complex
// #else
// #define CMPLX_ARG complex&
// #endif

 template<class T1 =double,class T2 =complex> struct doesnorm {};
 template<class T> struct doesnorm<T,T> : public  ::std::unary_function<T,T> {
   T operator()(const T& x) const { return x*x; }
 };

 template<typename T> inline double norm(const T& a)
   { return Sum_<LCM_DIM(T)>::func(a,doesnorm<double,LCM_VAL(T)>()); }

  template<> inline double norm(const complex&);
}

#ifdef LCM_USE_STDCOMPLEX
#define LCM_HAS_CONJ_MEMBER 0
#define LCM_HAS_NEGATE_MEMBER 0
#define LCM_HAS_PUBLICREIM 0
#include <complex>

namespace libcmatrix {
  typedef ::std::complex<double> complex;
  inline void real(complex& Z, float_t R) { Z.real()=R; }
  inline void imag(complex& Z, float_t I) { Z.imag()=I; }
}
#else

#ifdef LCM_USE_SSECOMPLEX
#define LCM_HAS_CONJ_MEMBER 0
#define LCM_HAS_NEGATE_MEMBER 0
#include "lcm_ssecomplex.h"
#else
#define LCM_HAS_CONJ_MEMBER 1
#define LCM_HAS_NEGATE_MEMBER 1
#include "lcm_cmatrix_complex.h"
#endif

#define LCM_HAS_PUBLICREIM 0

namespace libcmatrix {

#ifdef LCM_NEED_128BIT_ALIGN
  template<> struct memory_traits<complex> {
    static const size_t alignment=16;
  };
#endif

  inline void real(complex& Z, float_t R) { Z.real(R); }
  inline void imag(complex& Z, float_t I) { Z.imag(I); }

  complex cos(const complex&);
  complex cosh(const complex&);
  complex exp(const complex&);
  complex log(const complex&); 
  complex log10(const complex&);
  complex pow(float_t, const complex&);
  complex pow(const complex&, float_t);
  complex pow(const complex&, const complex&);
  complex pow(const complex&, int);
  complex sin(const complex&);
  complex sinh(const complex&);
  complex sqrt(const complex&);
  complex tan(const complex&);
  complex tanh(const complex&);
  
  // these are standard functions in std::complex, therefore defined here with class
  inline float_t real(const complex& z) { return z.real(); }
  inline float_t imag(const complex& z) { return z.imag(); }
  inline float_t arg(const complex& z) { return std::atan2(z.imag(), z.real()); }
  inline float_t abs(const complex& z) { return std::sqrt(norm(z)); }
  // Create a complex object given polar coordinates
  inline complex polar(float_t magv, float_t argv) { return magv*expi(argv); }

/* Complex stream I/O */
  std::ostream& operator<< (std::ostream&, const complex&);
  std::istream& operator>> (std::istream&, complex&);
}
#endif

namespace libcmatrix {

//   bool operator< (const complex& z1, const complex& z2)
//   { return norm(z1)<norm(z2); }

//   bool operator> (const complex& z1, const complex& z2)
//   { return norm(z1)>norm(z2); }

  complex quickpow(const complex& z, size_t n); //!< z^n calculated by repeated squaring (fast for modest n)

 template<> struct type_traits<complex> {
   static const bool trivialconstructor=true;
   static const size_t dimensionality=0;
   static const size_t rank=800;
   typedef complex value_type;
#ifdef LCM_NEED_128BIT_ALIGN
   static const size_t alignment=16;
#endif
 };
 
/*   template<class T> inline void */
/*   real_multiply(T& d, const complex& b, const complex& c) */
/*   { d=b.real()*c.real()-b.imag()*c.imag(); } */

//  inline void multiply_conj(complex& zd, const complex& z1,const complex& z2)
//    { zd.re=z1.re*z2.re+z1.im*z2.im; zd.im=z1.im*z2.re-z1.re*z2.im; }

 // (z1*)*z2;
#ifndef LCM_USE_SSECOMPLEX
 inline complex multiply_conj(const complex& z1,const complex& z2)
 { return complex(z1.real()*z2.real()+z1.imag()*z2.imag(), z2.real()*z1.imag()-z2.imag()*z1.real()); }

 inline complex multiply_conj(const double v,const complex& z)
 { return complex(v*z.real(),-v*z.imag()); }

  inline double real_multiply(const complex& b, const complex& c)
  { return b.real()*c.real()-b.imag()*c.imag(); }
  
  inline double real_conj_multiply(const complex& z1,const complex& z2)
  { return z1.real()*z2.real()+z1.imag()*z2.imag(); }
#endif
  inline complex conj_multiply(const complex& z1, const complex& z2) { return multiply_conj(z2,z1); }
  
//  inline void conj_multiply(complex& zd, const complex& z1,const complex& z2)
//  { zd=conj_multiply(z1,z2); }

 // z1*(z2*)
//  inline complex multiply_conj(const complex& z1,const complex& z2)
//    { return complex(z1.re*z2.re+z1.im*z2.im, z1.im*z2.re-z1.re*z2.im); }
 
 // (z1*)*v;
  inline complex conj_multiply(const complex& z1,float_t v)
  //  { return complex(z1.real()*v,-z1.imag()*v); }
  { return multiply_conj(v,z1); }
 // v*(z2*);
//  inline complex multiply_conj(float_t v,const complex& z2)
//    { return complex(z2.re*v,-z2.im*v); }

//  inline void conj_multiply(complex& zd, const complex& z1,const complex& z2)
//    { zd.re=z1.re*z2.re+z1.im*z2.im; zd.im=z1.re*z2.im-z1.im*z2.re; }
 
//  inline double real_conj_multiply(const complex& z1,const complex& z2)
//    { return z1.re*z2.re+z1.im*z2.im; }
 
//  inline void multiply_conj(complex& zd, const complex& z1,const complex& z2)
//    { zd.re=z1.re*z2.re+z1.im*z2.im; zd.im=z1.im*z2.re-z1.re*z2.im; }

//  // (z1*)*z2;
//  inline complex conj_multiply(const complex& z1,const complex& z2)
//    { return complex(z1.re*z2.re+z1.im*z2.im, z1.re*z2.im-z1.im*z2.re); }
//  // z1*(z2*)
//  inline complex multiply_conj(const complex& z1,const complex& z2)
//    { return complex(z1.re*z2.re+z1.im*z2.im, z1.im*z2.re-z1.re*z2.im); }
 
//  // (z1*)*v;
//  inline complex conj_multiply(const complex& z1,float_t v)
//    { return complex(z1.re*v,-z1.im*v); }
//  // v*(z2*);
//  inline complex multiply_conj(float_t v,const complex& z2)
//    { return complex(z2.re*v,-z2.im*v); }
 
 template<> struct doesnegate_ip<complex> {
   typedef complex argument_type;
   void operator()(complex& a) const {
   //   { a.re=-a.re; a.im=-a.im; }
#if LCM_HAS_NEGATE_MEMBER
     a.negate();
#else
     a=-a;
#endif
   }
 };
 template<> struct doesnegate<complex,double> : public ::std::unary_function<double,complex> {
   complex operator()(double a) const { 
     return complex(-a); }
 };

 // a+=b*c;
/*  template<class T> struct doesreal_mla : public binary_function_ip<T,complex,complex> { */
/*    inline void operator()(T& a, const complex& b, const complex& c) const  */
/*    { a=b.re*c.re-b.im*c.im; } */
/*  }; */

/*  template<> struct doesmla<complex,complex,complex> : public binary_function_ip<complex,complex,complex>  { */
/*    inline void operator()(complex& a,const complex& b,const complex& c) const */
/*    //   { a.re+=b.re*c.re-b.im*c.im; a.im+=b.im*c.re+b.re*c.im; } */
/*    { a+=b*c; } */
/*    //   { mla(a,b,c); } */
/*  }; */
 
/*  template<> struct doesmla<complex,float_t,complex> : public binary_function_ip<complex,float_t,complex> { */
/*    inline void operator()(complex& a, float_t b, const complex& c) const */
/*    //  { a.re+=b*c.re; a.im+=b*c.im; } */
/*    { a+=b*c; } */
/*    //   { throw Failed(""); } */
/*  }; */
 
/*  template<> struct doesmla<complex,float_t,float_t> : public binary_function_ip<complex,float_t,float_t>  { */
/*    inline void operator()(complex& a, float_t b, float_t c) const */
/*    { a+=b*c; } */
/*    //   { mla(a,b,c); } */
/*  }; */
 
/*  template<> struct doesmla<complex,complex,float_t> : public binary_function_ip<complex,complex,float_t> { */
/*    inline void operator()(complex& a, const complex& c, float_t b) const */
/*    { a+=b*c; } */
/*    //  { a.re+=b*c.re; a.im+=b*c.im; } */
/*    //   { mla(a,b,c); } */
/*  }; */

template<typename T> struct doesconj_mla : public binary_function_ip<complex,complex,T> {
  inline void operator()(complex& d, const complex& a, const T& b) const {
    conj_mla(d,a,b);
  }
};

  template<> struct checksnonzero_tolerance<complex> : public ::std::unary_function<complex,bool> {
    checksnonzero_tolerance(double tolv) : tol_(tolv*tolv) {}
    bool operator()(const complex& z) const { return (norm(z)>tol_); }
    double tol_;
  };
 
//  template<class T> inline void real_mla(T& a, const complex& b, const complex& c)
//    { a+=b.re*c.re-b.im*c.im; }
 
//  // a+=b*(c*)
//  inline void mla_conj(complex &a,const complex& b,const complex& c)
//    { a.re+=b.re*c.re+b.im*c.im; a.im+=b.im*c.re-b.re*c.im; }
//  inline void mla_conj(complex &a,float_t b,const complex& c)
//    { a.re+=b*c.re; a.im-=b*c.im; }

//  // a+=(b*)*c
//  inline void conj_mla(complex &a,const complex& b,const complex& c)
//    { a.re+=b.re*c.re+b.im*c.im; a.im+=b.re*c.im-b.im*c.re; }
//  inline void conj_mla(complex &a,const complex& b,float_t c)
//    { a.re+=b.re*c; a.im-=b.im*c; }
 
//  inline void real_mla_conj(float_t &a,const complex& b,const complex& c)
//    { a+=b.re*c.re+b.im*c.im; }

  template<class T> inline void real_mla(T& a, const complex& b, const complex& c)
   // { a+=b.real()*c.real()-b.imag()*c.imag(); }
  // inline void real_mla(float_t& a, const complex& b, const complex& c)
   { a+=real_multiply(b,c); }
 
 // a+=b*(c*)
  inline void mla_conj(complex &a,const complex& b,const complex& c)
   //{ a+=complex(b.real()*c.real()+b.imag()*c.imag(),b.imag()*c.real()-b.real()*c.imag()); }
    { a+=multiply_conj(b,c); }

  inline void mla_conj(complex &a, const float_t b, const complex& c)
    { a+=multiply_conj(b,c); }

 // a+=(b*)*c

  inline void conj_mla(complex &a,const complex& b,const complex& c)
    //  { a+=complex(b.real()*c.real()+b.imag()*c.imag(),b.real()*c.imag()-b.imag()*c.real()); }
    { a+=multiply_conj(c,b); }
  inline void conj_mla(complex &a,const complex& b,float_t c)
    //    { a+=complex(b.real()*c,-c*b.imag()); }
    { a+=multiply_conj(c,b); }
  
  inline void real_mla_conj(float_t &a,const complex& b,const complex& c)
    //  { a+=b.real()*c.real()+b.imag()*c.imag(); }
    { a+=real_conj_multiply(c,b); }
 
 template<> struct doesnorm<double,complex> : public  ::std::unary_function<complex,double> {
   //   double operator()(const complex& z) const { return z.real()*z.real()+z.imag()*z.imag(); } };
   double operator()(const complex& z) const { return norm(z); } };

 template<class T1 =double,class T2 =complex> struct doesreal {};
 template<> struct doesreal<double,complex> : public  ::std::unary_function<complex,double> {
   double operator()(const complex& z) const { return real(z); }
 };

 template<class T1 =double,class T2 =complex> struct doesimag {};
 template<> struct doesimag<double,complex> : public  ::std::unary_function<complex,double> {
   double operator()(const complex& z) const { return imag(z); }
 };

  inline complex expi(float_t _theta) {
    double _sina,_cosa;
    cmatrix_sincos(_theta,_sina,_cosa);
    return complex(float_t(_cosa),float_t(_sina));
  }

  template<class T> bool isreal(const T& a) {
    if (a.empty())
      throw Undefined("isreal");
    typename T::const_iterator aiter=a.begin();
    for (size_t n=a.size();n--;++aiter) {
      if (imag(*aiter))
	return false;
    }
    return true;
  }

  template<class T> bool isreal(const T& a,double tol) {
    if (a.empty())
      throw Undefined("isreal");
    typename T::const_iterator aiter=a.begin();
    if (tol<0.0)
      throw InvalidParameter("isreal");
    for (size_t n=a.size();n--;++aiter) {
      if (fabs(imag(*aiter))>tol)
	return false;
    }
    return true;
  }

 template<class T =complex> struct doesconj {};
 template<> struct doesconj<complex> : public  ::std::unary_function<complex,complex> { 
   complex operator()(const complex& z) const { return complex(z.real(),-z.imag()); }
 };

 template<class T =complex> struct doesconj_ip {};
  template<> struct doesconj_ip<complex> { void operator()(complex& v) const { imag(v,-(v.imag())); } };

 template<class T1,class T2> inline void real(T1& d,const T2& a)
   { apply(d,doesreal<LCM_VAL(T1),LCM_VAL(T2)>(),a); }
 template<class T1,class T2> inline void imag(T1& d,const T2& a)
   { apply(d,doesimag<LCM_VAL(T1),LCM_VAL(T2)>(),a); }
  template<class T1,class T2> inline void conj(T1& d,const T2& a)
   { apply(d,doesconj<LCM_VAL(T2)>(),a); }
 template<class T1,class T2> inline void enorm(T1& d,const T2& a)
   { apply(d,doesnorm<LCM_VAL(T1),LCM_VAL(T2)>(),a); }

 template<class T> void conj_ip(T& a)
   { apply_ip(doesconj_ip<LCM_VAL(T)>(),a); }
 
//  template<> struct doesadd_sd<complex,complex,complex> : public binary_function_ip<complex,complex,complex> {
//    inline void operator()(complex& d,const complex& a, const complex& b) const
//    { d.re=a.re+b.re; d.im=a.im+b.im; }
//  };
 
//  template<> struct doesadd_sd<complex,complex,double> : public binary_function_ip<complex,complex,double>  {
//    inline void operator()(complex& d, const complex& a, double b) const
//    { d.re=a.re+b; d.im=a.im; }
//  };
 
//  template<> struct doesadd_sd<complex,double,complex> : public binary_function_ip<complex,double,complex>  {
//    inline void operator()(complex& d, double b, const complex& a) const
//    { d.re=a.re+b; d.im=a.im; }
//  };
 
//  template<> struct doessubtract_sd<complex,complex,complex> : public binary_function_ip<complex,complex,complex> {
//    inline void operator()(complex& d,const complex& a, const complex& b) const
//    { d.re=a.re-b.re; d.im=a.im-b.im; }
//  };
 
//  template<> struct doessubtract_sd<complex,complex,double> : public binary_function_ip<complex,complex,double>  {
//    inline void operator()(complex& d, const complex& a, double b) const
//    { d.re=a.re-b; d.im=a.im; }
//  };
 
//  template<> struct doessubtract_sd<complex,double,complex> : public binary_function_ip<complex,double,complex>  {
//    inline void operator()(complex& d, double b, const complex& a) const
//    { d.re=b-a.re; d.im=a.im; }
//  };
 
//  template<> struct doesmultiply_sd<complex,complex,complex> : public binary_function_ip<complex,complex,complex> {
//    inline void operator()(complex& d, const complex& z1, const complex& z2) const
//    { d.re=z1.re*z2.re-z1.im*z2.im, d.im=z1.re*z2.im+z1.im*z2.re; }
//  };
 
//  template<> struct doesmultiply_sd<complex,complex,double> : public binary_function_ip<complex,complex,double> {
//    inline void operator()(complex& d, const complex& z, double v) const
//    { d.re=z.re*v, d.im=z.im*v; }
//  };
 
//  template<> struct doesmultiply_sd<complex,double,complex> : public binary_function_ip<complex,complex,double> {
//    inline void operator()(complex& d, double v,const complex& z) const
//    { d.re=z.re*v, d.im=z.im*v; }
//  };
  
} //namespace libcmatrix

#endif /* _COMPLEX_H_ */
