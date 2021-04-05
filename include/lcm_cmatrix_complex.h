#ifndef LCM_COMPLEX_H_
#define LCM_COMPLEX_H_

// internal file defining normal complex type

namespace libcmatrix {

class complex {
public:
  typedef float_t value_type;

  complex(float_t v, float_t i) : re(v), im(i) {}
  LCM_IMPLICIT_COMPLEX complex(float_t v) : re(v), im(0.0) {}
  complex() {} // NB. Unlike other complex implementations, not initialised to zero in the interests of speed

  complex& operator= (float_t v) { re=v; im=0.0; return *this; }
    complex& operator+= (const complex& z) { re+=z.re; im+=z.im; return *this; } 
    complex& operator+= (float_t v) { re+=v; return *this; }
    complex& operator-= (const complex& z) { re-=z.re; im-=z.im; return *this; }
    complex& operator-= (float_t v) { re-=v; return *this; }
    complex& operator*= (const complex& z)
      { register float_t r=re; re=r*z.re-im*z.im; im=im*z.re+r*z.im; return *this; } 
    complex& operator*= (float_t v) { re*=v; im*=v; return *this; }
    complex& operator/= (const complex& z) { return (*this=(*this)/z); }
    complex& operator/= (float_t v) { re/=v; im/=v; return *this; }

  friend float_t real(const complex&);
  friend float_t imag(const complex&);
  friend complex conj(const complex&);

    // member functions also defined by better conformance with std::complex
  const float_t& real() const { return re; }
  const float_t& imag() const { return im; }
    void real(float_t rev) { re=rev; }
    void imag(float_t imv) { im=imv; }
    
// The conj and negate member functions are additional
    complex& conj() { im=-im; return *this; }
    complex& negate() { re=-re; im=-im; return *this; } 

// Note norm returns the square of magnitude; technically incorrect, but
// kept for compatibility other complex implementations

// Most of the following have been changed to call-by-reference
// Unary Operator Functions
    friend inline complex operator- (const complex& z) { return complex(-z.re,-z.im);}

// Binary Operator Functions
    friend complex operator+ (const complex& z1, const complex& z2)
      { return complex(z1.re + z2.re, z1.im + z2.im); }

    friend complex operator+ (float_t v, const complex& z)
        { return complex(v + z.re, z.im); }

    friend complex operator+ (const complex& z, float_t v)
        { return complex(z.re+v,z.im); }

    friend complex operator- (const complex& z1, const complex& z2)
        { return complex(z1.re-z2.re, z1.im-z2.im); }

    friend complex operator- (float_t v1, const complex& z2)
	{ return complex(v1-z2.re,-z2.im); }

    friend complex operator- (const complex& z1, float_t v2)
        { return complex(z1.re-v2, z1.im); }

    friend complex operator* (const complex& z1, const complex& z2)
        { return complex(z1.re*z2.re-z1.im*z2.im, z1.re*z2.im+z1.im*z2.re); }

    friend complex operator* (const complex& z, float_t v)
        { return complex(z.re*v, z.im*v); }

    friend complex operator* (float_t v, const complex& z)
        { return complex(z.re*v, z.im*v); }

    friend complex operator/ (const complex& z1, const complex& z2) {
      const float_t _n=norm(z2);
      return complex( (z1.re*z2.re+z1.im*z2.im)/_n, (z1.im*z2.re-z1.re*z2.im)/_n);
    }

    friend complex operator/ (const complex& z1, float_t v2)
        { return complex(z1.re/v2, z1.im/v2); }

    friend complex operator/ (float_t v1, const complex& z2)
      { float_t scale=v1/(z2.re*z2.re+z2.im*z2.im); return complex(scale*z2.re,-scale*z2.im); }

    // Usual caveats about comparing floating point numbers!
    friend bool operator== (const complex& z1, const complex& z2) 
        { return (z1.re==z2.re) && (z1.im==z2.im); }
    friend bool operator== (const complex& z, float_t v)
      { return (z.re==v) && (z.im==0.0); }
    friend bool operator== (float_t v, const complex& z)
      { return (z.re==v) && (z.im==0.0); }

    friend bool operator!= (const complex& z1, const complex& z2)
        { return (z1.re!=z2.re) || (z1.im!=z2.im); }
    friend bool operator!= (const complex& z, float_t v)
        { return (z.re!=v) || (z.im!=0.0); }
    friend bool operator!= (float_t v, const complex& z)
        { return (z.re!=v) || (z.im!=0.0); }

#if !LCM_HAS_PUBLICREIM
private:
#endif
  float_t re,im;
};

  inline complex conj(const complex& z) { return complex(z.real(), -z.imag()); }
  template<> inline double norm(const complex& z) { return z.real()*z.real()+z.imag()*z.imag(); }

} //namespace libcmatrix

#endif
