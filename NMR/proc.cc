#include "List.h"
#include "cmatrix.h"
#include "rmatrix.h"
#include <cmath>
#include "ScratchList.h"
#include "cmatrix_utils.h"

namespace libcmatrix {

  using ::std::exp;

static const double gaussian_factor=M_PI/(2.0*std::sqrt(std::log(2.0)));

void make_gaussian(BaseList<double> d, double lb, double dt)
{
  const double gfac=lb*dt*gaussian_factor;
  const size_t n=d.size();
  d.front()=1.0;
  for (size_t c=1;c<n;c++) {
    const double x=gfac*c;
    d(c)=exp(-x*x);
  }
}

template<class T> inline void 
_exponential_multiply_ip(BaseList<T>& d, double cons)
{
  const size_t rs=d.length();
  const double rfac=exp(-cons/rs);
  double val=1.0;
  for (size_t r=0;r<rs;r++) {
    d(r)*=val;
    val*=rfac;
  }
}

void 
exponential_multiply_ip(BaseList<double> d, double cons)
{
  _exponential_multiply_ip(d,cons);
}

void
exponential_multiply_ip(BaseList<complex> d, double cons)
{
  _exponential_multiply_ip(d,cons);
}

// Note this is optimised for exponentials applied along rows
// This row/column indexing is not optimal for memory access

template<class T> void
_exponential_multiply_ip(Matrix<T>& d, double rcons, double ccons, size_t skip)
{
  if (!d)
    throw Undefined("exponential_multiply_ip");

  const size_t rs=d.rows();
  const size_t cs=d.cols();

  switch (skip) {
  case 0:
    throw InvalidParameter("exponential_multiply: skip factor must be >0");
  case 1:
    break;
  default:
    if (rs % skip)
      throw Failed("exponential_multiply: no. of rows is not multiple of skip factor");
  }
  if ((ccons==0) && (rcons==0))
    return;

  const size_t rsteps=rs/skip;
  const double rrfac=rcons ? std::exp(-rcons/rsteps) : 1.0;
  const double ccfac=ccons ? std::exp(-ccons/cs) : 1.0;    

  double cfac=1.0;
  size_t r;
    
  for (size_t c=0;c<cs;c++) {
    if (rcons) {
      r=0;
      double fac=cfac;
      for (size_t ro=0;ro<rsteps;ro++) {
	for (size_t i=skip;i--;)
	  d(r++,c)*=fac;
	fac*=rrfac;
      }
    }
    else {
      for (r=rs;r--;)
	d(r,c)*=cfac;
    }
    cfac*=ccfac;
  }
}

void 
exponential_multiply_ip(cmatrix& a, double rtc, double ctc, size_t skip)
{
  _exponential_multiply_ip( static_cast< Matrix<complex>& >(a),rtc,ctc,skip);
}
	  
void
exponential_multiply_ip(rmatrix& a, double rtc, double ctc, size_t skip)
{
  _exponential_multiply_ip( static_cast< Matrix<float_t>& >(a),rtc,ctc,skip);
}
	  
List<complex> 
exponential_multiply(const BaseList<complex>& a, double tc)
{
  List<complex> d(a,mxflag::temporary);
  exponential_multiply_ip(d,tc);
  return d;
}

List<double>
exponential_multiply(const BaseList<double>& a, double tc)
{
  List<double> d(a,mxflag::temporary);
  exponential_multiply_ip(d,tc);
  return d;
}

rmatrix 
exponential_multiply(const rmatrix& a, double rtc, double ctc, size_t skip) 
{
  rmatrix d(a,mxflag::temporary);
  exponential_multiply_ip(d,rtc,ctc,skip);
  return d;
}

cmatrix 
exponential_multiply(const cmatrix& a, double rtc, double ctc, size_t skip)
{
  cmatrix d(a,mxflag::temporary);
  exponential_multiply_ip(d,rtc,ctc,skip);
  return d;
}

void
phase_correct_ip(BaseList<complex> a, double zero, double first, double pivot)
{
  size_t n=a.size();
  if (n==0)
    return;
  if (first==0.0) {
    if (zero)
      a*=expi(zero);
    return;
  }
  const double scale=1.0/n;
  for (;n--;) {
    const double phasec=phase_correction(n*scale,zero,first,pivot);
    a(n)*=expi(phasec);
  }
}

//only calculate expi once per column
void
phase_correct_ip(cmatrix& a, double zero, double first, double pivot)
{
  if (!a)
    throw Undefined("phase_correct_ip");
  if (first==0.0) {
    if (zero)
      a*=expi(zero);
    return;
  }
  const size_t nr=a.rows();
  size_t nc=a.cols();
  const double scale=1.0/nc;
  for (;nc--;) {
    const double phasec=zero+first*(nc*scale-pivot);
    const complex p=expi(phasec);
    for (size_t r=nr;r--;)
      a(r,nc)*=p;
  }
}
  
void 
phase_correct_ip(cmatrix& a, double zero2, double first2, double pivot2, double zero1, double first1, double pivot1)
{
  phase_correct_ip(a,zero2+zero1,first2,pivot2);
  if (first1) {
    size_t nr=a.rows();
    const double scale=1.0/nr;
    for (;nr--;) {
      const double phasec=first1*(nr*scale-pivot1);
      a.row(nr)*=expi(phasec);
    }
  }
}

inline void 
multiply_ip_pair(complex& C, complex& S, const complex& z)
{
  const complex correctedre=complex(C.real(),S.real())*z;
  const complex correctedim=complex(C.imag(),S.imag())*z;
  C=complex(correctedre.real(),correctedim.real());
  S=complex(correctedre.imag(),correctedim.imag());
//   const complex corrected=complex(re,im)*z;
//   re=real(corrected);
//   im=imag(corrected);
}

void
ampphase_correct_ip(cmatrix& ca, cmatrix& sa, double zero2, double first2, double pivot2, double zero1, double first1, double pivot1)
{
  if (!arematching(ca,sa))
    throw Mismatch("ampphase_correct_ip");

  phase_correct_ip(ca,zero2,first2,pivot2);
  phase_correct_ip(sa,zero2,first2,pivot2);

  if ((zero1==0.0) && (first1==0.0))
    return;

  size_t nr=ca.rows();
  const double scale=1.0/nr;
  for (;nr--;) {
    BaseList<complex> crow=ca.row(nr);
    BaseList<complex> srow=sa.row(nr);
    const double phasec=zero1+first1*(nr*scale-pivot1);
    const complex z=expi(phasec);
    for (size_t c=ca.cols();c--;) {
//       multiply_ip_pair(crow(c).re,srow(c).re,z);
//       multiply_ip_pair(crow(c).im,srow(c).im,z);
      multiply_ip_pair(crow(c),srow(c),z);
    }
  }
}

template<typename T> void 
gaussian_multiply_ip_(Matrix<T>& d, double lbr, double dtr, double lbc, double dtc, size_t skip)
{
  if (!d)
    throw Undefined("gaussian_multiply_ip");

  const size_t rs=d.rows();
  const size_t cs=d.cols();

  switch (skip) {
  case 0:
    throw InvalidParameter("gaussian_multiply: skip factor must be >0");
  case 1:
    break;
  default:
    if (rs % skip)
      throw Failed("gaussian_multiply: no. of rows is not multiple of skip factor");
  }
  if ((lbr==0) && (lbc==0))
    return;
  if ( (lbr && (dtr<=0.0)) || (lbc && (dtc<=0.0)))
    throw InvalidParameter("gaussian_multiply: required dwell time is <=0");

  const size_t slowrows=rs/skip;
  ScratchList<double> rfacs(slowrows);
  if (lbr)
    make_gaussian(rfacs,lbr,dtr);
  ScratchList<double> cfacs(cs);
  if (lbc)
    make_gaussian(cfacs,lbc,dtc);

  size_t r=0;
  for (size_t slowrow=0;slowrow<slowrows;slowrow++) {
    for (size_t loop=skip;loop--;) {
      BaseList<T> currow(d.row(r));
      if (lbr) {
	const double rscale=rfacs(slowrow);
	if (lbc) {
	  for (size_t c=cs;c--;)
	    currow(c)*=rscale*cfacs(c);
	}
	else
	  currow*=rscale;
      }
      else
	multiply_ip(currow,cfacs);
      r++;
    }
  }
}

template<typename T> void 
gaussian_multiply_ip_(BaseList<T> d, double lb, double dt)
{
  if (d.empty())
    throw Undefined("gaussian_multiply_ip");
  
  const size_t n=d.size();

  if (lb==0.0)
    return;
  if (dt<=0.0)
    throw InvalidParameter("gaussian_multiply: dwell time is <=0");

  const double gfac=lb*dt*gaussian_factor;
  for (size_t c=1;c<n;c++) {
    const double x=gfac*c;
    d(c)*=exp(-x*x);
  }
}

void 
gaussian_multiply_ip(cmatrix& a, double lbr, double dtr, double lbc, double dtc, size_t skip)
{
  gaussian_multiply_ip_(a,lbr,dtr,lbc,dtc,skip);
}
	  
void 
gaussian_multiply_ip(BaseList<complex> a, double lb, double dt)
{
  gaussian_multiply_ip_(a,lb,dt);
}
	  
}//namespace libcmatrix
