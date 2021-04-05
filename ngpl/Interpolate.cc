		#include "Lineshapes.h"

// cubic interpolation from Numerical Recipes via SIMPSON 

/* Given four data points y1,y2,y3,y4 with x-values starting from x1
   and incremented with dx1, return the function value
   at x-value x using a polynomial of third degree as the interpolating
   function. Calculated using Lagrange's classical formula e.g. described in
   Numerical Recipes. */

namespace libcmatrix {

template<typename T> T pol3int(double x,double x1,double dx,const T y[4])
{   
   const double t1=x-x1;
   const double t2=t1-dx;
   const double t3=t2-dx;
   const double t4=t3-dx;

   const double dx3=dx*dx*dx;
   const double tedge=t2*t3/(6*dx3);
   const double tmid=t1*t4/(2*dx3);

   T r((-t4*tedge)*y[0]);
   mla(r, t3*tmid,y[1]);
   mla(r,-t2*tmid,y[2]);
   mla(r, t1*tedge,y[3]);
  
   if (verbose_interpolate)
     std::cout << "Interpolated (" << x1 << ',' << y[0] << ") (" << (x1+dx) << ',' << y[1] << ") (" << (x1+2*dx) << ',' << y[2] << ") (" << (x1+3*dx) << ',' << y[3] << ") at " << x << " to give " << r << '\n'; 
   return r;   
}

//! interpolate with values outside range clamped to \a padding
template<typename T> void cubic_interpolate_(BaseList<T>& dest, double newstart, double newdx, const BaseList<T>& source, double oldstart, double olddx, T padding =T(0))
{
  const size_t newnp=dest.size();
  const size_t oldnp=source.size();
  if ((oldnp<3) || (newnp<1))
    throw Failed("cubic_interpolate: insufficient data points for interpolation");
  if ((olddx<=0.0) || (newdx<=0.0))
    throw InvalidParameter("cubic_interpolate: frequency step must be >=0");

  //const double olddx=(oldend-oldstart)/oldnp;
  //const double newdx=(newend-newstart)/newnp;
  //  const size_t oldnp2=oldnp-2;
  T y[4];
  for (size_t j=0;j<newnp;j++) {
    const double newf=j*newdx+newstart;
    const int oldi=int((newf-oldstart)/olddx+0.5); //!< nearest old point
    if (verbose_interpolate)
      std::cout << "Interpolating new point " << j << " (" << newf << ") starting from point " << oldi << '\n';
    for (size_t k=4;k--;) {
      int curi=oldi+k-1;
      y[k]=((curi<0) || (curi>=oldnp))
	? padding
	: source(size_t(curi));
    }
    dest(j)=pol3int(newf,oldstart+(oldi-1)*olddx,olddx,y);
//     if (oldi<1) {
//       dest(j)=pol3int(newf,oldstart,olddx,source(0U),source(1U),source(2U),source(3U));
//     }
//     else {
//       if (oldi>=oldnp2)
// 	dest(j)=pol3int(newf,oldstart+(oldnp-4)*olddx,olddx,source(oldnp2-2),source(oldnp2-1),source(oldnp2),source(oldnp2+1));
//       else
// 	dest(j)=pol3int(newf,oldstart+(oldi-1)*olddx,olddx,source(oldi-1),source(oldi),source(oldi+1),source(oldi+2));
//     }   
  }
}

//! interpolate with folding (can't change window)
template<typename T> void cubic_interpolate_(BaseList<T> dest, const BaseList<T>& source)
{
  const size_t newnp=dest.size();
  const size_t oldnp=source.size();
  if ((oldnp<3) || (newnp<1))
    throw Failed("cubic_interpolate: insufficient data points for interpolation");
  const double r=double(oldnp)/newnp;
  T y[4];
  for (size_t j=0;j<newnp;j++) {
    const double oldf=j*r;
    const size_t oldi=int(oldf+0.5); //!< nearest old point
    size_t k=0;
    if (j==0) {
      y[0]=source.back();
      k++;
    }    
    for (;k<4;k++) {
      size_t curi=oldi+k-1;
      if (curi>=oldnp)
	y[k]=source(size_t(curi-int(oldnp)));
      else
	y[k]=source(curi);
    }
    dest(j)=pol3int(oldf,oldi-1.0,1.0,y);
  }
}

void cubic_interpolate(BaseList<double> dest, double newstart, double newdx, const BaseList<double>& source, double oldstart, double olddx, double padding)
{
  cubic_interpolate_(dest,newstart,newdx,source,oldstart,olddx,padding);
}

List<double> cubic_interpolate(size_t newnp, double newstart, double newdx, const BaseList<double>& source, double oldstart, double olddx, double padding)
{
  List<double> dest(newnp);
  cubic_interpolate_(dest,newstart,newdx,source,oldstart,olddx,padding);
  return dest;
}

void cubic_interpolate(BaseList<double> dest, const BaseList<double>& source)
{
  cubic_interpolate_(dest,source);
}

List<double> cubic_interpolate(size_t newnp, const BaseList<double>& source)
{
  List<double> dest(newnp);
  cubic_interpolate_(dest,source);
  return dest;
}

void cubic_interpolate(BaseList<complex> dest, double newstart, double newdx, const BaseList<complex>& source, double oldstart, double olddx, complex padding)
{
  cubic_interpolate_(dest,newstart,newdx,source,oldstart,olddx,padding);
}

List<complex> cubic_interpolate(size_t newnp, double newstart, double newdx, const BaseList<complex>& source, double oldstart, double olddx, complex padding)
{
  List<complex> dest(newnp);
  cubic_interpolate_(dest,newstart,newdx,source,oldstart,olddx,padding);
  return dest;
}

void cubic_interpolate(BaseList<complex> dest, const BaseList<complex>& source)
{
  cubic_interpolate_(dest,source);
}

List<complex> cubic_interpolate(size_t newnp, const BaseList<complex>& source)
{
  List<complex> dest(newnp);
  cubic_interpolate_(dest,source);
  return dest;
}

} //namespace libcmatrix
