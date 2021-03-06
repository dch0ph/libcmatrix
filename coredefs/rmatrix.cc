// flag that Level 2 (matrix) data types should be instantiated
#include "basedefs.h"
#ifdef LCM_USE_EXTERNTEMPLATE
#define LCM_LEVEL2_INSTANT_REAL
#undef LCM_SUPPRESS_VIEWS
#endif
#include "rmatrix.h"
#include "cmatrix_external.h"

namespace libcmatrix {

#ifndef LCM_USE_EXTERNTEMPLATE
  template void Matrix<float_t>::dump() const;
  template void BaseList<float_t>::dump() const;
#endif

rmatrix identity(int n)
{
  rmatrix d(mxflag::temporary);
  d.identity(n);
  return d;
}

template<> bool tryspy_(const Matrix<double>& a, std::ostream& ostr, double tol) {
  spy(ostr,a,tol);
  return true;
}

template<> bool tryspy_(const Matrix<bool>& a, std::ostream& ostr, double) {
  spy(ostr,a);
  return true;
}

void spy(std::ostream& ostr,const BaseList<float_t>& a, double tol)
{ 
  checksnonzero_tolerance<double> nonzerocheck(tol);
  for (size_t c=0;c<a.length();c++) {
    const double val(a(c));
#ifdef HAVE_ISFINITE
    if (!std::isfinite(val)) {
      ostr << '!';
      continue;
    }
#endif
    ostr << (nonzerocheck(val)? 'X' : '.');
  }
}

void spy(std::ostream& ostr,const Matrix<double>& a,double tol)
{
  if (!a) {
    print(a,ostr);
    return;
  }
  for (size_t r=0;r<a.rows();r++) {
    spy(ostr,a.row(r),tol);
    if (r==a.rows()-1)
      ostr << std::endl;
    else
      ostr << '\n';
  }
}

#include "../coredefs/shared.cc"

void kronecker(rmatrix& d, const rmatrix& a, int n)
{
  kronecker_(d,a,n);
}

void kronecker(rmatrix& d, int n, const rmatrix& a)
{
  kronecker_(d,n,a);
}
	
void kronecker(rmatrix& d, const rmatrix& a, const rmatrix& b)
{
  kronecker_(d,a,b,false);
}

void kronecker_transpose(rmatrix& d, int n, const rmatrix& b)
{
  kronecker_(d,n,b,true);
}

void kronecker_transpose(rmatrix& d, const rmatrix& a, const rmatrix& b)
{
  kronecker_(d,a,b,true);
}

double stddev(const BaseList<double>& a)
{
  if (a.size()<3)
    throw Failed("stddev: called with data set of <3 values!");
  const size_t n=a.size();
  const double mean=sum(a)/n;
  double tot=0.0;
  for (size_t i=n;i--;) {
    const double diff=a(i)-mean;
    tot+=diff*diff;
  }
  return std::sqrt(tot/n); // fix Jan 07
}

bool hasoffdiagonal(const rmatrix& U, double ptol)
{
  if (!issquare(U))
    throw NotSquare("isdiagonal");
  if (ptol<0)
    throw InvalidParameter("hasoffdiagonal");

  for (size_t j=U.rows();j--;) {
    for (size_t i=j;i--;) {
      if ( (fabs(U(i,j))>ptol) || (fabs(U(j,i))>ptol))
	return true;
    }
  }
  return false;
}

void cmatrix_multiply(RawMatrix<double>& dest, const RawMatrix<double>& a, const RawMatrix<double>& b, bool acc)
{
  if (a.rows()<LCM_NAIVE_BELOW)
    multiply_naive_direct_(dest,a,b,acc);
  else
    multiply_direct_(dest,a,b,acc);
}

#ifdef LCM_USE_EXTERNAL

void lapack_multiply(RawMatrix<double>& d, const RawMatrix<double>& a, const RawMatrix<double>& b, bool acc, bool Atrans, bool Btrans)
{
  size_t ar,ac,br,bc;
  getrowcols(ar,ac,a,Atrans);
  getrowcols(br,bc,b,Btrans);
   
  if (ac!=br)
    throw Mismatch("lapack_multiply",ar,ac,br,bc);

#ifdef USE_SUNPERFACML
   LCM_EXT_PREFIX::dgemm(
	 Atrans ? 'Y' : 'N',
	 Btrans ? 'Y' : 'N',
	 ar,bc,ac,
	 1.0,
	 lapack_pass(a),a.step(),
	 lapack_pass(b),b.step(),
	 acc ? 1.0 : 0.0,
	 lapack_pass(d),d.step());
#else
#ifdef USE_ATLAS

  cblas_dgemm(LCM_ORDER,
	      Atrans ? CblasTrans : CblasNoTrans,
	      Btrans ? CblasTrans : CblasNoTrans,
	ar,bc,ac,
	1.0,
	      lapack_pass(a),a.step(),
	      lapack_pass(b),b.step(),
	      acc ? 1.0 : 0.0,
	      lapack_pass(d),d.step());
// #else
//   double one=1.0;
//   double zero=0.0;

//   dgemm_("N","N",&ar,&bc,&ac,
// 	 &one,
// 	 lapack_pass(a),&iac,
// 	 lapack_pass(b),&ibc,
// 	 &zero,
// 	 lapack_pass(d),&bc);  
#else
#error "Unsupported external library"
#endif
#endif
}

void lapack_multiply(BaseList<double>& d, const Matrix<double>& a, const BaseList<double>& b, bool Atrans)
{
  size_t ar,ac;
  getrowcols(ar,ac,a,Atrans);
  const size_t iac=a.cols();

  if ( (ac!=b.length()) || (ar!=d.length()))
    throw Mismatch("lapack_multiply");

#ifdef USE_SUNPERFACML
   LCM_EXT_PREFIX::dgemv(
	      Atrans ? 'Y' : 'N',
	 ar,ac,
	 1.0,
	 lapack_pass(a),iac,
	 lapack_pass(b),1,
	 0.0,
	 lapack_pass(d),1);
#else
#ifdef USE_ATLAS
  cblas_dgemv(LCM_ORDER,
	      Atrans ? CblasTrans : CblasNoTrans,
	ar,ac,
	1.0,
	lapack_pass(a),iac,
	lapack_pass(b),1,
	0.0,
	lapack_pass(d),1);
#else
#error "Unsupported external library"
#endif
#endif
}

#endif  //USE_EXTERNAL

}//namespace libcmatrix


