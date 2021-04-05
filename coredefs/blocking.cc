#define LCM_LEVEL3_INSTANT
#include "cmatrix.h"
#include "ListList.h"
#include "ScratchList.h"
#include "BlockedMatrix.h"

namespace libcmatrix {

#ifndef LCM_USE_EXTERNTEMPLATE
  template void ListList<size_t>::dump() const;
  template void ListList<float_t>::dump() const;
#endif
  
static void recursive_inc(const Matrix<bool> &connected,BaseList<bool> &inout,size_t n,size_t start)
{
  inout(start)=true;
  
  for (size_t i=n;i--;) {
    if (!inout(i) && connected(i,start))
      recursive_inc(connected,inout,n,i);
  }
}

ListList<size_t> find_blocks(const Matrix<complex> &a,double cutoff)
{
  if (!issquare(a))
    throw NotSquare("find_blocks");

  size_t i;
  const size_t n=a.rows();

  Matrix<bool> connected(n,n,false);

  cutoff*=cutoff;

  for (i=0;i<n;i++) {
    for (size_t j=i+1;j<n;j++) {
      if (norm(a(i,j))>cutoff || norm(a(j,i))>cutoff)
	connected(i,j)=connected(j,i)=true;
    }
  }

  ScratchList<bool,2*SCRATCH_SIZE> boolbase(2*n);
  BaseList<bool> used(n,boolbase.vector());
  BaseList<bool> inout(n,boolbase.vector()+n);
  ScratchList<size_t,2*SCRATCH_SIZE> sizebase(2*n);
  BaseList<size_t> tmpsize(n,sizebase.vector());
  BaseList<size_t> asvec(n,sizebase.vector()+n);

  used=false;

  size_t blockcount=0;
  size_t totcount=0;
  
  for (;totcount<n;) {
    size_t start;

    for (start=0;used(start);start++);

    if (start>=n)
      throw InternalError("find_blocks");

    inout=false;

    recursive_inc(connected,inout,n,start);

    size_t count=0;
    for (i=0;i<n;i++) {
      if (inout(i)) {
	asvec(totcount++)=i;
	used(i)=true;
	count++;
      }
    }

    if (!(count && blockcount<n))
      throw InternalError("find_blocks");
    
    tmpsize(blockcount++)=count;
  }
  return ListList<size_t>(tmpsize.truncate(blockcount),asvec);
}

   
// Construct block list from identical diagonal elements

namespace {
  struct SameWithinTol : public std::binary_function<double,double,bool> {
    SameWithinTol(double tol_) : tol(tol_) {
      if (tol_<0.0)
	throw InvalidParameter("SameWithinTol");
    }
    bool operator()(double a, double b) const { 
      return (std::fabs(a-b)<tol);
    }
    double tol;
  };
}

ListList<size_t> find_blocks(const BaseList<double>& a,double tolerance)
{
  return find_blocks(a,SameWithinTol(tolerance));
}

double hermitian_trace(BlockedMatrix<complex>& a)
{
  if (!a)
    throw Undefined("trace");
  double sum(0.0);
  for (size_t i=a.size();i--;)
    sum+=hermitian_trace(a(i));
  return sum;
}

template<typename T> void pow(BlockedMatrix<T>& d, const BlockedMatrix<T>& a, int n)
{
  if (!a)
    throw Undefined("pow");
  d.duplicate_structure(a);
  for (size_t i=a.size();i--;) {
    cmatrix& di(d(i));
    if (!!di)
      pow(di,a(i),n);
  }
}

template void pow(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, int);

}//namespace libcmatrix
