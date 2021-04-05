#include <cstdlib>
#include "cmatrix_utils.h"
#include "ScratchList.h"

namespace libcmatrix {

void ft_create_table(BaseList<complex> &w,int isign)
{
  const size_t n=w.length();
  if (n<1) throw InvalidParameter("ft_create_table: number of points must be >0");

  const complex sfac=expi(2*isign*M_PI/n);

  complex fac(1.0);

  for (size_t i=0;i<n;i++) {
    w(i)=fac;
    fac*=sfac;
  }
}

void ft_row(complex *dest,const complex *source,const complex *facs,int n,double scale)
{
  complex start=scale*source[0];

  for (size_t i=0;i<n;i++) {
    size_t k=i;
    complex sum(start);

    for (size_t j=1;j<n;j++) {
      mla(sum,source[j],facs[k]);
      k+=i;
      if (k>=n) k-=n;
    }
    dest[i]=sum;
  }
}

void ft_row(double *dest,const complex *source,const complex *facs,int n,double scale)
{
  double start=scale*real(source[0]);

  for (size_t i=0;i<n;i++) {
    size_t k=i;
    double sum(start);

    for (size_t j=1;j<n;j++) {
      real_mla(sum,source[j],facs[k]);
      k+=i;
      if (k>=n) k-=n;
    }
    dest[i]=sum;
  }
}

template<class T> void _ft(BaseList<T> &r,const BaseList<complex> &d,int isign,double scale)
{
  const size_t n=d.length();
  if (r.length()!=n)
    throw Mismatch("ft",r.length(),n);

  ScratchList<complex> w(n);
  ft_create_table(w,isign);
  ft_row(r.vector(),d.vector(),w.vector(),n,scale);
}


void ft(BaseList<complex>& dest, const BaseList<complex>& source, int dir, double scale)
{
  _ft(dest,source,dir,scale);
}

void real_ft(BaseList<double>& dest, const BaseList<complex>& source, int dir, double scale)
{
  _ft(dest,source,dir,scale);
}

template<class T> void _ft(BaseList<T>& r,const BaseList<complex>& d,const BaseList<complex>& facs, double scale)
{
  const size_t n=d.length();
  if (n!=facs.length() || n!=r.length())
    throw Mismatch("ft");
  ft_row(r.vector(),d.vector(),facs.vector(),n,scale);
}

void real_ft(BaseList<double>& d, const BaseList<complex>& s, const BaseList<complex>& facs, double scale)
{
  _ft(d,s,facs,scale);
}

void ft(BaseList<complex>& d, const BaseList<complex>& s, const BaseList<complex>& facs, double scale)
{
  _ft(d,s,facs,scale);
}

void ft(List<complex>& d, const BaseList<complex>& s, const BaseList<complex>& facs, double scale)
{
  d.create(s.length());
  _ft(d,s,facs,scale);
}

void real_ft(BaseList<complex>& d, const BaseList<complex>& s, const BaseList<complex>& facs, double scale)
{
  _ft(d,s,facs,scale);
}

void ft(cmatrix& b, const cmatrix& a, int isign, double scalefac)
{
  if (!a)
    throw Undefined("ft");

  const size_t n=a.cols();
  b.create(a.rows(),n);
  
  ScratchList<complex> w(n);
  ft_create_table(w,isign);
  
  for (size_t r=a.rows();r--;)
    ft_row((b.row(r)).vector(),(a.row(r)).vector(),w.vector(),n,scalefac);
}

List<complex> fft(const BaseList<complex>& a,int dir,double scale)
{
  List<complex> d(a,mxflag::temporary);
  fft_ip(d,dir,scale);
  return d;
}

List<complex> fft(const BaseList<double>&a,int dir,double scale)
{
  List<complex> d(a,mxflag::temporary);
  fft_ip(d,dir,scale);
  return d;
}

List<complex> ft(const BaseList<complex> &a,int dir,double scale)
{
  List<complex> d(mxflag::temporary);
  ft(d,a,dir,scale);
  return d;
}

cmatrix ft(const cmatrix &a,int dir,double scale)
{
  cmatrix d(mxflag::temporary);
  ft(d,a,dir,scale);
  return d;
}

bool ispowerof2(int n)
{
  if (n<1) 
    return false;
  return ((n & (n-1))==0);
  //   while (!(n & 1)) n>>=1;
//   return (n==1);
}

inline void validate_ftargs(int n,int sign)
{
  if (!ispowerof2(n))
    throw InvalidParameter("n must be power of 2");
  if (sign!=FT_FORWARD && sign!=FT_BACKWARD)
    throw InvalidParameter("FT sign must be +/- 1");
}

void _fft(complex *data,int nn,int isign,double scale)
{
  int i,j,istep,m,mmax;
  double wtemp,theta;
  complex tmp,wp,w;
  
  data[0]*=scale;
  
  j=0;
  for (i=0;i<nn;i++) {
    if (j>i) {
      tmp=data[j]; data[j]=data[i]; data[i]=tmp;
    }
    m=nn>>1;
    while (m>=1 && j>=m) {
      j-=m;
      m>>=1;
    }
    j+=m;
  }
  mmax=1;
  while (nn>mmax) {
    istep=2*mmax;
    theta=(M_PI)/(isign*mmax);
    wtemp=std::sin(0.5*theta);
    wp=complex(-2.0*wtemp*wtemp,std::sin(theta));
    
    w=1.0;
    for (m=0;m<mmax;m++) {
      for (i=m;i<nn;i+=istep) {
	j=i+mmax;
	tmp=w*data[j];
	data[j]=data[i]-tmp;
	data[i]+=tmp;
      }
      w+=w*wp;
    }
    mmax=istep;
  }
}

void fft_ip(BaseList<complex> &dlist,int isign,double scale)
{
  const size_t n=dlist.length();
  validate_ftargs(n,isign);
  _fft(dlist.vector(),n,isign,scale);
}

void fft_ip(cmatrix &a,int isign,double scalefac)
{
  if (!a) throw Undefined("fft_ip");
  for (size_t i=a.rows();i--;) {
    BaseList<complex> row=a.row(i);
    fft_ip(row,isign,scalefac);
  }
}

cmatrix fft(const cmatrix& a,int isign,double scalefac)
{
  cmatrix b(a,mxflag::temporary);
  fft_ip(b,isign,scalefac);
  return b;
}

void phasefft_ip(cmatrix &a,int isign,double rscale,double cscale)
{
  fft_ip(a,isign,rscale);

  size_t i,j;
  const size_t r=a.rows();

  ScratchList<complex> tmp(r);

  for (i=a.cols();i--;) {
    for (j=r;j--;)
      tmp(j)=a(j,i);
    tmp.front()*=cscale;
    fft_ip(tmp,isign,cscale);
    for (j=r;j--;)
      a(j,i)=tmp(j);
  }
}

cmatrix phasefft(const cmatrix &a,int isign,double rscale,double cscale)
{
  cmatrix b(a,mxflag::temporary);
  phasefft_ip(b,isign,rscale,cscale);
  return b;
}

/* Assume null sin fids, throw away resultant hypercomplex part */
void ampfft_ip(cmatrix &a,int isign,double rscale,double cscale)
{
  fft_ip(a,isign,rscale);

  size_t i,j;
  const size_t r=a.rows();
  const size_t c=a.cols();

  ScratchList<complex> tmpi(r);
  ScratchList<complex> tmpr(r);

  for (i=c;i--;) {
    for (j=r;j--;) {
      tmpr(j)=real(a(j,i));
      tmpi(j)=imag(a(j,i));
    }
    tmpr.front()*=cscale;
    tmpi.front()*=cscale;
    fft_ip(tmpr,isign,cscale);
    fft_ip(tmpi,isign,cscale);
    for (j=r;j--;) a(j,i)=complex(real(tmpr(j)),real(tmpi(j)));
  }
}

cmatrix ampfft(const cmatrix& a, int isign, double rscale, double cscale)
{
  cmatrix b(a,mxflag::temporary);
  ampfft_ip(b,isign,rscale,cscale);
  return b;
}

void ampfft_ip(cmatrix& cf, cmatrix& sf, int isign, double rscale, double cscale)
{
  const size_t r=cf.rows();
  const size_t c=cf.cols();

  if (r!=sf.rows() || c!=sf.cols())
    throw Mismatch("ampfft_ip",r,c,sf.rows(),sf.cols());

  fft_ip(cf,isign,rscale);
  fft_ip(sf,isign,rscale);

  size_t i,j;

  ScratchList<complex> tmpi(r);
  ScratchList<complex> tmpr(r);

  for (i=c;i--;) {
    for (j=r;j--;) {
      tmpr(j)=complex(real(cf(j,i)),real(sf(j,i)));
      tmpi(j)=complex(imag(cf(j,i)),imag(sf(j,i)));
    }
    tmpr.front()*=cscale;
    tmpi.front()*=cscale;
    fft_ip(tmpr,isign,cscale);
    fft_ip(tmpi,isign,cscale);
    for (j=r;j--;) {
      cf(j,i)=complex(real(tmpr(j)),real(tmpi(j)));
      sf(j,i)=complex(imag(tmpr(j)),imag(tmpi(j)));
    }
  }
}

}//namespace libcmatrix
