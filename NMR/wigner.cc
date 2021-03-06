#include <cstdlib>
#include "basedefs.h"
#include "wigner.h"

namespace libcmatrix {

static const double deg_to_rad=M_PI/180.0;

static double fact(int a)
{
  double r=1.0;
  for (;a>0;a--)
    r*=a;
  return r;
}

std::ostream &operator << (std::ostream &ostr,const Euler &EA)
{
  return ostr << '(' << EA.alpha/deg_to_rad << ", " << EA.beta/deg_to_rad << ", " << EA.gamma/deg_to_rad << ')';
}

complex D(int rank,int n,int m,const Euler& EA)
{
  return wigner_factor(n,m,EA.alpha,EA.gamma)*d(rank,n,m,EA.beta); //!< Feb 13 swapped n,m for consistency with space_T::rotate
}

Matrix<complex> D(int L,const Euler& EA)
{
  Matrix<complex> d(mxflag::temporary);
  d.create(2*L+1,2*L+1);

  for (int m=-L;m<=L;m++)
    for (int n=-L;n<=L;n++)
      d(m+L,n+L)=D(L,m,n,EA);
  return d;
}

double legendre(int n,double beta)
{
  switch (n) {
  case 2:
    return 0.75*std::cos(2*beta)+0.25;
  case 1:
    return std::cos(beta);
  case 0:
    return 1.0;
  }
  return d(n,0,0,beta);
}

double d(int J,int m,int n,double beta)
{
  switch (J) {
  case 0:
    if (m==0 && n==0)
      return 1.0;
    break;
  case 1: 
    return d1(m,n,beta);
  case 2:
    return d2(m,n,beta);
  }

  if ((J < 0) || (::std::abs(m) > J) || (::std::abs(n) >J))
    throw InvalidParameter("Bad Wigner matrix element");

  const double betaover2 = beta/2;
  double COSF,SINF,COSB,mSINB;
  cmatrix_sincos(betaover2,COSB,mSINB); mSINB=-mSINB;
  double djval = 0.0;					// Value of Wigner element
  int Jpm=J+m, Jmm=J-m, Jpn=J+n, Jmn=J-n;		// We need these
  double prefact =					// Constant prefactor
        std::sqrt(fact(Jpm)*fact(Jmm)*fact(Jpn)*fact(Jmn));
  double signf=1.0, denom=0.0;
  int dd1,dd2,dd3;
  for (int k=0; k<=2*J; k++) {
    dd1 = Jmm-k;
    dd2 = Jpn-k;
    dd3 = k+m-n;
    if ((dd1>=0) && (dd2>=0) && (dd3>=0)) {
      COSF = std::pow(COSB, double(Jpn + Jmm - 2*k));		// Cosine to power
      SINF = std::pow(mSINB, double(2*k + m - n));			// Sine to power
      denom = fact(dd1)*fact(dd2)*fact(dd3)*fact(k);		// Get demoninator
      djval += signf*COSF*SINF/denom;			// Unscaled dJ contribution
    }
    signf = -signf;					// Switch sign 1 <-> -1
  }
  return djval*prefact;
}
    
  //!< more-or-less checks out by Table B.2 of Klaus SR (p451) - must have come from somewhere else

double d2(int m,int n,double beta)
{
  //NB shouldn't be global otherwise creates initialisation order problems
  // static const double sqrthalf=std::sqrt(0.5);
  static const double sqrt1half=std::sqrt(1.5);
  static const double sqrt3o8=std::sqrt(3.0/8);
  
  double x;

  if (::std::abs(n)>::std::abs(m))
    return (((m-n) % 2) ? -1 : 1)*d2(n,m,beta);

  switch (m) {
  case -2: case -1:
    return (((m-n) % 2) ? -1 : 1)*d2(-m,-n,beta);
  case 2:

    switch (n) {
    case -2:
      x=std::sin(beta/2);
      x*=x;
      return x*x;
    case -1:
      return 0.5*std::sin(beta)*(std::cos(beta)-1);
    case 0:
      x=std::sin(beta);
      return sqrt3o8*x*x;
    case 1:
      return -0.5*std::sin(beta)*(1+std::cos(beta));
    case 2:
      x=std::cos(beta/2);
      x*=x;
      return x*x;
    }
  case 1:
    x=std::cos(beta);
    switch (n) {
    case -1:
      return (x+0.5)*(1-x);
    case 0:
      return -sqrt1half*std::sin(beta)*x;
    case 1:
      return (x-0.5)*(1+x);
    }
  case 0:
    x=std::cos(beta);
    if (n==0) return 0.5*(3*x*x-1);

  default:
    throw InvalidParameter("d(2)");
  }
}

double d1(int m,int n,double beta )
{
  static const double sqrthalf=std::sqrt(0.5);

  double r;

  switch(m) {
  case -1:
    switch (n) {
    case -1: return  d1(1,1,beta);
    case  0: return -d1(1,0,beta);
    case  1: return  d1(1,-1,beta);
    }
    break;

  case 0:
    switch(n) {
    case -1: return -sqrthalf*std::sin(beta);
    case  0: return std::cos(beta);
    case  1: return sqrthalf*std::sin(beta);
    }
    break;

  case 1:
    switch(n) {
    case -1:
      r=std::sin(beta/2);
      return r*r;
    case 0:
      return -sqrthalf*std::sin(beta);
    case 1:
      r=std::cos(beta/2);
      return r*r;
    }
    break;
  }
  throw InvalidParameter("d(1)");
} 

}//namespace libcmatrix
