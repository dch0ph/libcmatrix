#include <pmmintrin.h>

#ifndef LCM_NEED_128BIT_ALIGN
#error SSE complex requires 128 bit alignment to be enabled
#endif

namespace libcmatrix {

  typedef __m128d ssecomplex_t __attribute__ (( __aligned__ (16)));

  class complex {
public:
  complex() {}
  
  complex(double re, double im) 
    : z_(_mm_set_pd(im,re)) {}

  LCM_IMPLICIT_COMPLEX complex(double x) 
    : z_(_mm_set_sd(x)) {}

  complex& operator= (const double x) { z_=_mm_set_sd(x); return *this; }

  explicit complex(const ssecomplex_t& z) : z_(z) {}
  operator ssecomplex_t() const { return z_; }

  double real() { return reinterpret_cast<double&>(z_); }
  const double& real() const { return reinterpret_cast<const double& >(z_); }
  void real(const double);
  void imag(const double);

  double imag() const {
    ssecomplex_t swapped(_mm_shuffle_pd(z_,z_,1));
    return reinterpret_cast<const double&>(swapped); }

  friend complex operator+ (const complex&, const double);
  friend complex operator- (const complex&, const double);
  friend complex operator* (const complex&, const double);
  friend complex operator/ (const complex&, const double);
  friend complex operator+ (const complex&, const complex&);
  friend complex operator- (const complex&, const complex&);
  friend complex operator* (const complex&, const complex&);
  friend complex operator/ (const complex&, const complex&);
  friend complex operator+ (const double, const complex&);
  friend complex operator- (const double, const complex&);
  friend complex operator* (const double, const complex&);
  friend complex operator/ (const double, const complex&);
  friend complex conj(const complex& z);
  friend complex multiply_conj(double, const complex&);

  complex& operator+= (const complex& z2) { z_=_mm_add_pd(z_,z2.z_); return *this; }
  complex& operator-= (const complex& z2) { z_=_mm_sub_pd(z_,z2.z_); return *this; }
  complex& operator*= (const complex& z2) { z_=(*this)*z2; return *this; }
  complex& operator/= (const complex& z2) { z_=(*this)/z2; return *this; }
  complex& operator+= (const double b) { *this=(*this)+b; return *this; }
  complex& operator-= (const double b) { *this=(*this)-b; return *this; }
  complex& operator*= (const double b) { *this=(*this)*b; return *this; }
  complex& operator/= (const double b) { *this=(*this)/b; return *this; }

  // Usual caveats about comparing floating point numbers!
  friend bool operator== (const complex& z1, const complex& z2) { 
    const ssecomplex_t cmpres(_mm_cmpneq_pd(z1.z_,z2.z_));
    return (_mm_movemask_pd(cmpres)==0);
  }

  friend bool operator== (const complex& z, float_t v)
  { return (z==complex(v)); }
  friend bool operator== (float_t v, const complex& z)
  { return (z==complex(v)); }
  
  friend bool operator!= (const complex& z1, const complex& z2) { 
    const ssecomplex_t cmpres(_mm_cmpneq_pd(z1.z_,z2.z_));
    return (_mm_movemask_pd(cmpres)!=0);
  }
    
  friend bool operator!= (const complex& z, float_t v)
  { return (z!=complex(v)); }
  friend bool operator!= (float_t v, const complex& z)
  { return (z!=complex(v)); }
  
private:
  ssecomplex_t z_;
} __attribute__ (( __aligned__ (16)));

inline void complex::real(const double re)
{ z_= _mm_loadl_pd(z_, &re); }

inline void complex::imag(const double im)
{ z_= _mm_loadh_pd(z_, &im); }

inline complex operator+ (const complex& z1, const complex& z2) 
{ return complex(_mm_add_pd(z1.z_,z2.z_)); }

inline complex operator- (const complex& z1, const complex& z2)
{ return complex(_mm_sub_pd(z1.z_,z2.z_)); }

inline complex operator+ (const complex& z1, const double b)
{ return complex(_mm_add_sd(z1.z_, _mm_load_sd(&b))); }

inline complex operator+ (const double b, const complex& z1)
{ return complex(_mm_add_sd(z1.z_, _mm_load_sd(&b))); }

inline complex operator- (const complex& z1, const double b)
{ return complex(_mm_sub_sd(z1.z_, _mm_load_sd(&b))); }

inline complex operator- (const double b, const complex& z1)
{ return complex(_mm_sub_pd(_mm_load_sd(&b), z1.z_)); }

inline complex operator* (const complex& a, const double b) 
{ return complex(_mm_mul_pd(a.z_, _mm_set1_pd(b))); }

inline complex operator* (const double b, const complex& a) 
{ return complex(_mm_mul_pd(a.z_, _mm_set1_pd(b))); }

inline complex operator/ (const complex& a, const double b) 
{ return complex(_mm_mul_pd(a.z_, _mm_set1_pd(1.0/b))); }

inline complex operator* (const complex& a, const complex& b) {
  const ssecomplex_t b_im = _mm_shuffle_pd(b,b,3); // Imag. part of b in both
  const ssecomplex_t b_re = _mm_shuffle_pd(b,b,0); // Real part of b in both
  const ssecomplex_t tmp=_mm_mul_pd(a,b_re);
  const ssecomplex_t a_flipped = _mm_shuffle_pd(a,a,1); // Swap real and imag parts of a
  const ssecomplex_t tmp2=_mm_mul_pd(a_flipped, b_im);
  return complex(_mm_addsub_pd(tmp,tmp2));
}

inline double real_multiply(const complex& a, const complex& b) {
  ssecomplex_t tmp=_mm_mul_pd(a,b);
  tmp=_mm_hsub_pd(tmp,tmp);
  return reinterpret_cast<const double& >(tmp);
}

inline double real_conj_multiply(const complex& a, const complex& b) {
  ssecomplex_t tmp=_mm_mul_pd(a,b);
  tmp=_mm_hadd_pd(tmp,tmp);
  return reinterpret_cast<const double& >(tmp);
}
 
template<> double norm(const complex& b)
{
  const ssecomplex_t b_square(_mm_mul_pd(b,b));
  const ssecomplex_t res(_mm_hadd_pd(b_square,b_square));
  return reinterpret_cast<const double& >(res);
}

template<typename T> void mla(complex& a, const T& b, const complex& c)
{
  a+=b*c;
}

inline complex multiply_conj(const complex& a, const complex& b)
{
  static const union {
    unsigned int i[4]; ssecomplex_t v; //!< 9/2/2016  int -> unsigned int 
  } signbithigh = {{0,0,0,0x80000000}};
  ssecomplex_t b_im = _mm_shuffle_pd(b,b,3); // Imag. part of b in both
  const ssecomplex_t b_re = _mm_shuffle_pd(b,b,0); // Real part of b in both
  const ssecomplex_t tmp=_mm_mul_pd(a,b_re);
  b_im = _mm_xor_pd(b_im, signbithigh.v); // Change sign of high
  const ssecomplex_t a_flipped = _mm_shuffle_pd(a,a,1); // Swap real and imag parts of a
  //return complex(_mm_mul_pd(a, b_re)) +
  return complex(tmp+_mm_mul_pd(a_flipped, b_im));
}

inline complex multiply_conj(double b, const complex& a)
{
  ssecomplex_t b_re = _mm_set1_pd(b);
  static const union {
    unsigned int i[4]; ssecomplex_t v;
  } signbithigh = {{0,0,0,0x80000000}};
  b_re = _mm_xor_pd(b_re, signbithigh.v); // Change sign of high
  return complex(_mm_mul_pd(a.z_, b_re));
}

inline complex operator / (const complex& a, const complex& b)
{
  return multiply_conj(a,b)/norm(b);
}

inline complex operator / (const double a, const complex& b)
{
  return multiply_conj(a,b)/norm(b);
}

inline complex operator- (const complex& a) {
  static const union { // (signbit,signbit)
    unsigned int i[4]; ssecomplex_t v;
  } signbits = {{0,0x80000000,0,0x80000000}};
  return complex(_mm_xor_pd(a, signbits.v)); // Change sign of both elements
}

inline complex conj(const complex& a) {
  static const union { // (signbit,signbit)
    unsigned int i[4]; ssecomplex_t v;
  } signbithigh = {{0,0,0,0x80000000}};
  return complex(_mm_xor_pd(a.z_, signbithigh.v)); // Change sign of imag. part
}

} //namespace libcmatrix
