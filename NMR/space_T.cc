#include "space_T.h"
#include "Matrix.h"
#include "wigner.h"
#include <cstdlib>
#include <sstream>

using std::fabs;

namespace libcmatrix {

template<> std::ostream& operator << (std::ostream& ostr,const Tensor<double>& a)
{ 
  a.printsimple(ostr);
  return ostr;
}

template<> std::ostream& operator << (std::ostream& ostr,const Tensor<complex>& a)
{ 
  a.printsimple(ostr);
  return ostr;
}

space_T A1(double x,double y,double z)
{
  space_T a(1,mxflag::maximum);

  for (int m=-1;m<=1;m++)
    a(1,m)=A1(x,y,z,m);
  return a;
}
 
BadRank::BadRank(const char* errm_, int rank_, int maxrank_)
  : MatrixException("Bad rank"),
    std::invalid_argument(errm_)
{
  std::ostringstream ostr(errm_);
  if (rank_<maxrank_)
    ostr << " (requested rank " << rank_ << " is missing)";
  else {
    ostr << " (rank " << rank_ << " requested when ";
    if (maxrank_>=0)
      ostr << "maximum rank is " << maxrank_ << ')';
    else
      ostr << "tensor is undefined)";
  }
  errm=ostr.str();
}
 
static const double sqrthalf=std::sqrt(0.5);
static const double sqrt6th=std::sqrt(1.0/6);
static const double sqrt3rd=std::sqrt(1.0/3);

complex A1(double x,double y,double z,int m)
{
  switch (m) {
  case -1:
    return sqrthalf*complex(x,-y);
  case 0:
    return complex(z);
  case 1:
    return -sqrthalf*complex(x,y);
  }
  throw BadIndex("A1");
}

bool areequal(const BaseList<complex>& a, const BaseList<complex>& b, double tol)
{
  const double tol2=tol*tol;
  for (size_t j=a.size();j--;) {
    const complex& aj(a(j));
    const complex& bj(b(j));
    const double sum=norm(aj)+norm(bj);
    if ((norm(aj-bj)>tol*sum) && (sum>tol2))
      return false;
  }
  return true;
}

bool areequal(double a,double b,double tol) { 
  const double sum=fabs(a)+fabs(b);
  return (std::fabs(a-b)<(tol*sum)) || (sum<tol);
}
      
bool areequal(const space_T& a, const space_T& b, double tol)
  {
    if ((tol>=1.0) || (tol<=0.0))
      throw InvalidParameter("areequal");
    int l=a.max_rank();
    if (l!=b.max_rank())
      return false;
    for (;l>=0;l--) {
      const bool isactive(a.have_rank(l));
      if (isactive ^ b.have_rank(l))
	return false;
      if (isactive && !areequal(a(l),b(l),tol))
	return false;
    }
    return true;
  }

space_T A2(const Matrix<double>& a)
{
  if (a.cols()!=3 || a.rows()!=3) 
    throw Mismatch("A2");

  space_T T(2,mxflag::all);

  for (int rank=0;rank<=2;rank++) {
    for (int order=-rank;order<=rank;order++)
      T(rank,order)=A2(a,rank,order);
  }
  return T;
}

space_T A2(double xx,double yy,double zz)
{
  space_T T(2,mxflag::maximum);

  T.ensure_rank(0);
  T(0,0)=A2(xx,yy,zz,0,0);

  for (int order=-2;order<=2;order++)
    T(2,order)=A2(xx,yy,zz,2,order);
  return T;
}

space_T A2(double xx,double yy,double zz,int rank)
{
  space_T T(2);

  T.ensure_rank(rank);
  for (int order=-rank;order<=rank;order++)
    T(rank,order)=A2(xx,yy,zz,rank,order);

  return T;
}

space_T A2(const Matrix<double> &a,int rank)
{
  if (a.cols()!=3 || a.rows()!=3)
    throw Mismatch("A2");
  if (rank<0 || rank>2)
    throw BadRank("A2");

  space_T T(2);
  T.ensure_rank(rank);

  for (int order=-rank;order<=rank;order++)
    T(rank,order)=A2(a,rank,order);

  return T;
}

double A2(double xx,double yy,double zz,int rank,int order)
{
  switch (rank) {
  case 0:
    if (order==0)
      return -sqrt3rd*(xx+yy+zz);
    throw BadIndex("A2");
    
  case 1:
    if (::std::abs(order)>1)
      throw BadIndex("A2");
    return 0.0;

  case 2:
    switch (std::abs(order)) {
    case 2:
      return 0.5*(xx-yy);
    case 1:
      return 0.0;
    case 0:
      return (2*zz-yy-xx)*sqrt6th;
    }
    throw BadIndex("A2");
  }    
  throw BadRank("A2");
}

// Checked against Spiess & Schmidt-Rohr Table B.1
complex A2(const Matrix<double>& a, int rank, int order)
{
  if (a.cols()!=3 || a.rows()!=3)
    throw Mismatch("A2");

  switch (rank) {
  case 0:
    if (order==0)
      return complex(-sqrt3rd*trace(a));
    throw BadIndex("A2");

  case 1:
    switch (order) {
    case -1:
      return -0.5*(a(2,0)-a(0,2)-complex(0,a(2,1)-a(1,2)));
    case 0:
      return complex(0,-sqrthalf*(a(0,1)-a(1,0)));
    case 1:
      return -0.5*(a(2,0)-a(0,2)+complex(0,a(2,1)-a(1,2)));
    }
    throw BadIndex("A2");

  case 2:
    switch (order) {
    case -2:
      return 0.5*(a(0,0)-a(1,1)-complex(0,a(0,1)+a(1,0)));
    case -1:
      return 0.5*(a(0,2)+a(2,0)-complex(0,a(1,2)+a(2,1)));
    case 0:
      return complex(sqrt6th*(2.0*a(2,2)-a(1,1)-a(0,0)));
    case 1:
      return -0.5*(a(0,2)+a(2,0)+complex(0,a(1,2)+a(2,1)));
    case 2:
      return 0.5*(a(0,0)-a(1,1)+complex(0,a(0,1)+a(1,0)));
    }
    throw BadIndex("A2");

  }    
  throw BadRank("A2");
}

space_T rotate(const space_T& a, const Euler& EA)
{
  space_T T(a.rank());

  const double& beta=EA.beta;

  for (int l=a.rank();l>=0;l--) {
    if (!a.have_rank(l))
      continue;

    T.ensure_rank(l);

    if (l==0) {
      T(0,0)=a(0,0);
      continue;
    }

    if (EA.alpha || EA.gamma) {
      for (int m=-l;m<=l;m++) {
	complex sum(0.0,0.0);
	for (int n=-l;n<=l;n++)
	  mla(sum,D(l,n,m,EA),a(l,n));
	T(l,m)=sum;
      }
    }
    else {
      for (int m=-l;m<=l;m++) {
	complex sum(0.0,0.0);
	for (int n=-l;n<=l;n++)
	  mla(sum,d(l,n,m,beta),a(l,n)); //!< form agrees with B.25 of Klaus SR (p450)
	T(l,m)=sum;
      }
    }
  }
  return T;
}

inline int whichrank(int span) {
  if (span % 2)
    return (span-1)/2;
  else
    throw InvalidParameter("Dimension must be 2L+1"); 
}

inline int whichrank(const Matrix<complex> &DM) {
  if (!issquare(DM))
    throw NotSquare("Bad Wigner rotation matrix");
  else
    return whichrank(DM.rows());
}

// static void do_mult(BaseList<complex> dest, const BaseList<complex>& a, int l, const cmatrix& DM)
// {
//   for (int m=-l;m<=l;m++) {
//     complex sum(0.0,0.0);
//     const int ml=m+l;
//     for (int n=-l;n<=l;n++)
//       mla(sum,DM(n+l,ml),a[n]);
//     dest[m]=sum;
//   }
// }

complex rotate(const space_T& a,int m,const Matrix<complex>& DM)
{
  const int l=whichrank(DM);
  if ((m<-l) || (m>l))
    throw BadIndex("rotate",m,l);

  if (!a.have_rank(l))
    return complex(0.0);
  
  complex sum(0.0,0.0);
  const BaseList<complex> al(a(l));
  m+=l;
  for (size_t n=DM.rows();n--;)
    mla(sum,al(n),DM(n,m));
    
  return sum;
}

double sumzero(const space_T& a)
{
  int r=a.max_rank();
  if (r<0)
    throw Undefined("Tensor<T>::sumzero");
  double result=real(a(r,0));
  for (;r--;) {
    if (a.have_rank(r))
      result+=real(a(r,0));
  }
  return result;
}

double sumzero_rotate(const space_T& a,const Matrix<complex>& DM)
{
  const int r=a.max_rank();
  if (r<0)
    throw Undefined("sumzero_rotate");
  if (r==0)
    return real(a(0,0));

  if ((whichrank(DM)!=2) || (r!=2) || !a.have_rank(2))
    throw InvalidParameter("sumzero_rotate: only applicable to rank 2 tensors");

  double res=real(rotate(a,0,DM));
  if (a.have_rank(0))
    res+=real(a(0,0));
  
  return res;
}

void rotate(space_T& dest, const space_T& a, const Matrix<complex>& DM)
{
  dest.create(a.rank(),mxflag::none);
  const int modrank=whichrank(DM);

  for (int l=a.rank();l>=0;l--) {
    if (!a.have_rank(l))
      continue;
    dest.ensure_rank(l);

    if (l==0)
      dest(0,0)=a(0,0);
    else {
      if (l!=modrank)
	throw Mismatch("Ranks don't match");      
      multiply(dest(l),a(l),DM);
    }
  }
}

space_T rotate(const space_T& a,const Matrix<complex>& DM)
{
  space_T dest;
  rotate(dest,a,DM);
  return dest;
}

void rotate(space_T& dest, const space_T& a,const BaseList<Matrix<complex> >& DMS)
{
  if (DMS.length()<a.max_rank())
    throw Mismatch("Insufficient rank");

  dest.create(a.rank(),mxflag::none);

  for (int l=a.rank();l>=0;l--) {
    if (!a.have_rank(l))
      continue;
    
    dest.ensure_rank(l);

    if (l==0)
      dest(0,0)=a(0,0);
    else {
      const Matrix<complex> &DM=DMS(l);
      if (whichrank(DM)!=l)
	throw Mismatch("Ranks don't match");
      multiply(dest(l),a(l),DM);
    }
  }
}

space_T rotate(const space_T& a,const BaseList<Matrix<complex> >& DMS)
{
  space_T dest;
  rotate(dest,a,DMS);
  return dest;
}

complex rotate(const space_T& a,int l,int m,const Euler& EA)
{
  if (l>a.rank() || !a.have_rank(l))
    return complex(0.0);

  complex sum(0.0,0.0);

  if (EA.alpha || m*EA.gamma) {
    for (int n=-l;n<=l;n++)
      mla(sum,a(l,n),D(l,n,m,EA));
  }
  else {
    const double& beta=EA.beta;

    for (int n=-l;n<=l;n++)
      mla(sum,a(l,n),d(l,n,m,beta));
  }
    
  return sum;
}

void tensor_to_matrix(Matrix<double>& Am, const space_T& A)
{
  static const double minussqrt3=-std::sqrt(3.0);
  static const double sqrt6=std::sqrt(6.0);
  if (A.max_rank()!=2)
    throw InvalidParameter("tensor_to_matrix: expecting rank 2 tensor");
  if (A.have_rank(1U))
    throw InvalidParameter("tensor_to_matrix: tensor should not have antisymmetric component");

  Am.create(3,3,0.0);
  double trA=0.0;
  if (A.have_rank(0U))
    trA=minussqrt3*real(A(0,0));
  const complex S22=A(2,2)+A(2,-2);
  const complex D22=A(2,2)-A(2,-2);
  const complex S21=A(2,1)+A(2,-1);
  const complex D21=A(2,1)-A(2,-1);
  //    cout << S22 << "   " << D22 << "  " << S21 << "  " << D21 << '\n';
  Am(2,0)=Am(0,2)=-0.5*real(D21);
  Am(1,2)=Am(2,1)=-0.5*imag(S21);
  Am(1,0)=Am(0,1)=0.5*imag(D22);
  Am(2,2)=(sqrt6*real(A(2,0))+trA)/3.0;
  const double S0011=trA-Am(2,2);
  Am(0,0)=0.5*(real(S22)+S0011);
  Am(1,1)=0.5*(S0011-real(S22));
}

}//namespace libcmatrix
