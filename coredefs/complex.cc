#include <cmath>
#include "cmatrix_complex.h"
#include <cassert>

namespace libcmatrix {

  using ::std::pow;

void sinhcosh(double x,double &sinhx,double &coshx)
{
  register double eix=std::exp(x);
  register double eimx=1.0/eix;
  coshx=(eix+eimx)/2.0;
  sinhx=(eix-eimx)/2.0;
}

complex cos(const complex& x)
{
  double cosre,sinre, sinhim, coshim;
  cmatrix_sincos(x.real(),cosre, sinre);
  sinhcosh(x.imag(),sinhim,coshim);

  return complex(cosre*coshim, -sinre*sinhim);
}


complex cosh(const complex& x)
{
  double sinh_val,cosh_val,sinix,cosix;
  sinhcosh(x.real(),sinh_val,cosh_val);
  cmatrix_sincos(x.imag(),sinix,cosix);
  
  return complex(cosh_val*cosix, sinh_val*sinix);
}

complex sin(const complex& x)
{
  double sinh_val,cosh_val,sinix,cosix;
  sinhcosh(x.imag(),sinh_val,cosh_val);
  cmatrix_sincos(x.real(),sinix,cosix);
  
  return complex( cosh_val * sinix, sinh_val * cosix);
}


complex sinh(const complex& x)
{
  double sinh_val,cosh_val,sinix,cosix;
  sinhcosh(x.real(),sinh_val,cosh_val);
  cmatrix_sincos(x.imag(),sinix,cosix);
  
  return complex( sinh_val * cosix, cosh_val * sinix);
}


complex tan(const complex& x)
{
  double sin2re,cos2re,cosh2im,sinh2im;
  cmatrix_sincos(2*x.real(),sin2re,cos2re);
  sinhcosh(2*x.imag(),sinh2im,cosh2im);

  const double den = cos2re+cosh2im;
  
  return complex( sin2re/den, sinh2im/den);
}

complex tanh(const complex& x)
{
  double sinh2re,cosh2re,cos2im,sin2im;
  cmatrix_sincos(2*x.imag(),sin2im,cos2im);
  sinhcosh(2*x.real(),sinh2re,cosh2re);

  const double den = cosh2re+cos2im;
  
  return complex( sinh2re/den, sin2im/den);
}

complex pow(const complex& x,const complex& y)
{
  return exp(y*log(x));
}

complex pow(const complex& x,int n)
{
  if (x.imag()==0) {
    if (x.real()<0)
      return pow(x,float_t(n));
    return complex(pow(x.real(),float_t(n))); // this will catch x==0
  }
  return exp(float_t(n)*log(x));
}

complex pow(const complex& x,float_t y)
{
  return exp(y*log(x));
}

complex pow(float_t x,const complex& z)
{
  if (x<0)
    return pow(complex(x),z);
  if (z.imag()==0)
    return complex(pow(x,z.real()));
  return exp(z*std::log(x));
}

  static inline complex square(const complex& z) { return complex(z.real()*z.real()-z.imag()*z.imag(), 2.0*z.real()*z.imag()); }

complex quickpow(const complex& a, size_t n)
{
  if (n==0)
    return complex(1.0);
  bool done=false;
  complex result;
  complex lastpower(a);
  for (;;) {
    if (n & 1) {
      if (done)
	result*=lastpower;
      else {
	result=lastpower;
	done=true;
      }
    }    
    n>>=1;
    if (n)
      lastpower=square(lastpower);
    else {
      assert(done);
      return result;
    }
  }
}
    
complex log10(const complex& x)
{
  return complex( std::log10(abs(x)), arg(x));
}

complex log(const complex& x)
{
  return complex( std::log(abs(x)), arg(x));
}

complex exp(const complex& x)
{
  double sinix,cosix;
  cmatrix_sincos(x.imag(),sinix,cosix);
  
  if (x.real()==0.0)
    return complex(cosix, sinix);
  register float_t exp_val = std::exp(x.real());
  return complex( exp_val * cosix, exp_val*sinix);
}

complex sqrt(const complex& x)
{
  const float_t re=real(x);
  const float_t im=imag(x);

  if (im == 0.0) {
    if (re < 0.0)
      return complex( 0.0, LCM_COPYSIGN(std::sqrt(-re),im));
    else
      return complex( fabs(std::sqrt(re)), LCM_COPYSIGN(0.0,im));
  }
 
  if (re == 0.0) {
    const double r = std::sqrt (0.5 * fabs (im));
    return complex( LCM_COPYSIGN(r,im), r);
  }

  const double r=std::sqrt(hypot(re,im));
  const double p=arg(x)/2.0;

  double pc,ps;
  cmatrix_sincos(p,ps,pc);
  return complex(r*pc,r*ps);
}

namespace {
  inline std::ostream& dump_(std::ostream& ostr, double a, size_t width) {
    if (width)
      ostr << std::setw(width);
    return ostr << a;
  }

  std::ostream& pretty_print(std::ostream& ostr, double val, bool needintro)
  {
    if (needintro && (val>=0.0))
      ostr << '+';
    if (val!=1.0)
      ostr << val;
    return ostr;
  }
  
}

// int ostream_controller::complexpadding() const
// {
//   switch (complexview) {
//     case pair:
//     return 2;
//   case withi:
//     return 2;
//   }
//   return -1; //!< don't know
// }

std::ostream& ostream_controller::print(std::ostream& ostr, const complex& z, size_t width) const
{
  switch (complexview) {
  case pair:
    ostr << '(';
    dump_(ostr,z.real(),width) << ',';
    dump_(ostr,z.imag(),width) << ')';
    break;
  case withi: case compact: {
    const LCM_IOS::fmtflags oldflags=ostr.flags();
    bool isdone=false;
    const bool iscompact(complexview==compact);
    if (!iscompact || real(z)) {
      dump_(ostr,z.real(),width);
      isdone=true;
    }    
    if (!iscompact || imag(z)) {
      if (iscompact && (imag(z)==1.0)) {
	if (isdone || (oldflags & LCM_IOS::showpos))
	  ostr << '+';
	ostr << 'i';
      }
      else {
	if (isdone)
	  ostr.setf(LCM_IOS::showpos); //!< force leading + if previous output      
	dump_(ostr,z.imag(),width) << 'i';
      }
      isdone=true;
    }
    if (!isdone)
      dump_(ostr,0.0,width);
    ostr.flags(oldflags);
  }
    break;    
  default:
    throw InternalError("ostream_controller::print");
  }
  return ostr;
}

std::ostream& operator<< (std::ostream& ostr, const complex& z)
{
  cmatrix_ostream_controller(ostr).print(ostr,z);
  return ostr;
}

std::istream& operator >> (std::istream& istr, complex& z)
{
  char ch;
  float_t re,im;

  istr >> ch;

  if (ch == '(') {
    istr >> re >> ch;
    switch(ch) {
    case ',':
      istr >> im >> ch;
      if (ch == ')') {
	z=complex(re,im);
	return istr;
      }
      break;
    case ')':
      z = complex(re);
      return istr;
    }
  }
  else {
    istr.putback(ch);
    istr >> re;
    z= complex(re);
    return istr;
  }
  istr.setstate(LCM_IOS::failbit);
  return istr;
}

}//namespace libcmatrix

