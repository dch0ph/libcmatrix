#include "cmatrix.h"
#include "List.h"
#include "ScratchList.h"
#include <cmath>

using namespace libcmatrix;

namespace libcmatrix {

const double TINY=1e-20;

/* LU Decomposition *************************************-*-c++-*-
**							 	**
**  This routine is taken from "Numerical Recipies", W.H.Press,	**
**  B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge	**
**  University Press, 1989.  See pages 31 through 36.  The code **
**  has been modified for GAMMA into C++.			**
**							 	**
**  The algorithm takes an input matrix A & overwrites it with  **
**  the LU decompositions of A' where A' is A with any needed	**
**  permutations. Any permutations used	are stored in the	**
**  integer array indx.  The function returns the integer d	**
**  which will be +1 if an even	number of permutations was used **
**  or -1 if an odd number of permutations was necessary.  Sub- **
**  sequently, d can be used in taking the determinant of the	**
**  input matrix A (or LU).					** 
**								**  
**  The LU decomposition formulates the equation		**
**								**  
**		     	        A = LU				**  
**								**  
**  where L is lower triangular and U is upper triangular.	**
**								**
**  [A11 A12 A13 A14]   [L11  0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21 L22  0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32 L33  0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43 L44] [ 0   0   0  U44]	**
**								**
**  Both L and U have elements on the diagonal, but those of L	**
**  are always set to be 1: <i|L|i> = 1.			**
**								**
**  [A11 A12 A13 A14]   [ 1   0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21  1   0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32  1   0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43  1 ] [ 0   0   0  U44]	**
**								**
**  Neglecting the diagonal of L, both L and U can be overlaid	**
**  for storage in a single array (here overwriting A).		**
**								**
**                      [U11 U12 U12 U14]			**
**  			[L21 U22 U23 U24|			**
**  			[L31 L32 U33 U34|			**
**  			[L41 L42 L43 U44]			**
**								**
**  Finally, the algorithm uses Crouts method to determine the	**
**  elements of L and U.  This sets the diagonal elements of L	**
**  to 1 and follows a specific order in computation of the L	**
**  and U elements which allows A to be overwritten as L and U	**
**  are determined.  For Crouts method to be stable, partial	**
**  pivoting is used (row interchanges).			**
**								**
*****************************************************************/

static bool LU_decomp(cmatrix& a,BaseList<size_t> indx)

  // Input		A     : Input normal matrix (this)
  //			indx  : Integer array for permutations
  // Return		A     : LU Decomposition of input A
  // 			indx  : Index of Permutations
  // Note			      : Input matrix is overwritten!
  
{
  if (!issquare(a))
    throw NotSquare("LU_decomp");

  size_t n = a.rows();			// Get dimension of A
  if (n<indx.length())
    throw Mismatch("LU_decomp");
  static const complex complex_zero(0.0,0.0);

  int d = 1;				// Permutation counter (1 even, -1 odd)
  size_t i;
  static const double epsilon=1e-12;
  static const double tiny = 1.0e-20;		// Set this to a small number
  ScratchList<double> vv(n); // Vector for row scaling factors
  double big, temp;

  for(i=n;i--;) {			// For each row of A find the					// the largest element for pivoting 
    big=0.0;				// and store corresponding scaling
    for(int j=n; j--;) {		// factor in array vv
      temp = norm(a(i,j));
      if (temp > big)
	big = temp;
    }
    big=std::sqrt(big);
    if(big<epsilon)
      throw SingularMatrix<complex>("LU_decomp",a);
    vv(i)=1.0/big;
  }

  complex Uij, Lij;			// Element of U and L
  size_t imax=0;
  double dum;
  complex dumz;
  for (size_t j=0; j<n; j++) {		// Begin going through each A column

// First implement equation below to compute the values of U (except <i|U|i>)
//
//		       i-1
//		       ---
// <i|U|j> = <i|A|j> - \   <i|L|k><k|U|j>
// 		       /
//		       ---
//		       k=0
//
// This is equation 2.3.12 on page 34 of the referenced text and its use is
// restricted to terms above the diagonal.

    for(i=0; i<j; i++) {
      Uij = a(i,j);
      for (size_t k=0; k<i; k++)
        Uij -= a(i,k)*a(k,j);
      a(i,j)=Uij;
    }

// Next implement the equation below to compute the values of U (i>j, i<n)
// For the case when i=j, this equation is the same as thh previous equation
// when the scaling factor is neglected (and this is done).
//
//			       j-1
//	        1	       ---
// <i|L|j> = ------- <i|A|j> - \   <i|L|k><k|U|j>
// 	     <j|U|j>	       /
//			       ---
//			       k=0

    big=0.0;				// Used to look for largest pivot
    for (i=j; i<n; i++) {
      Lij = a(i,j);
      for (size_t k=0; k<j; k++)
	Lij -= a(i,k)*a(k,j);
      a(i,j)=Lij;
      dum = vv(i) * std::sqrt(norm(Lij));   	// For looking at the pivot
      if (dum >= big) {			// Use it if better that big
	big=dum;
	imax=i;
      }
    }

    if (j != imax) {			// This part interchanges the
                			// rows of A via permutation
      for (size_t k=0; k<n; k++) {		// The sign of d is switched
        				// and the permutation stored
        dumz = a(imax,k);		// in the array indx
        a(imax,k) = a(j,k);
	a(j,k)=dumz;
      }
      d = -d;
      vv(imax) = vv(j);
    }

    indx(j) = imax;
    if (a(j,j)==complex_zero)		// If <j|A|j> is zero make it
      a(j,j)=tiny;			// small instead!
    if (j != n-1) {			// Divide by the pivot element
      					// If <j|A|j> was zero it has 
      dumz = 1.0/a(j,j);		// since been set to tiny so it
      for (i=j+1; i<n; i++)		// won't blow up!
        a(i,j)*= dumz;
    }
  }					// Return to treat next column
  return true;
}

static bool LU_decomp(Matrix<double>& a,BaseList<size_t> indx)
{
  if (!issquare(a))
    throw NotSquare("LU_decomp");

  int i,imax,j,k;
  double big,dum,temp,sum;
  int d=1; /* No row interchanges */

  size_t n=a.rows();
  if (n<indx.length())
    throw Mismatch("LU_decomp");

  ScratchList<double> vv(n);

  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++) 
      if ((temp=fabs(a(i,j)))>big)
	big=temp;

    if (big==0.0)
      return false;

    vv(i)=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a(i,j);
      for (k=0;k<i;k++)
	sum-=a(i,k)*a(k,j);
      a(i,j)=sum;
    }
    big=0.0;
    imax=-1;
    for (i=j;i<n;i++) {
      sum=a(i,j);
      for (k=0;k<j;k++)
	sum-=a(i,k)*a(k,j);
      a(i,j)=sum;
      /* Find best pivot */
      if ((dum=vv(i)*fabs(sum))>=big) {
	big=dum;
	imax=i;
      }		
    }
    if (imax<0)
      throw SingularMatrix<double>("LU_decomp",a);

    if (j!=imax) {
      for (k=0;k<n;k++) {
	dum=a(imax,k);
	a(imax,k)=a(j,k);
	a(j,k)=dum;
      }
      d=-d;
      vv(imax)=vv(j);
    }
    indx(j)=imax;
    if (a(j,j)==0.0)
      a(j,j)=TINY;
    if (j!=n-1) {
      dum=1.0/a(j,j);
      for (i=j+1;i<n;i++)
	a(i,j)*=dum;
    }
  }
  return true;
}

/* LU Back Substitution *********************************-*-c++-*-
**							 	**
**  This routine is taken from "Numerical Recipies", W.H.Press,	**
**  B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge	**
**  University Press, 1989.  See pages 31 through 38.  The code **
**  has been modified for GAMMA into C++.			**
**							 	**
**  The algorithm takes an input matrix A' in its LU decomp-	**
**  osition format where A' is some original matrix A with any	**
**  needed row permutations to attain the LU form. The row	**
**  permutations used to relate A to A' are stored in the	**
**  integer array indx.  The function the proceeds to solve	**
**								**  
**			A|x> = |b>				**
**								**  
**  for |x> by considering the problem as			** 
**								**  
**              A|x> = (LU)|x> = L(U|x>) |b> 			**
**								**  
**  and first solving for the vector |y> where |y> = U|x>	**
**								**  
**			L|y> = |b>				**
**								**
**  followed by solving for |x> by				**
**								**
**			U|x> = |y>				**
**								**
**  Due to the triagular nature of the arrays L and U it is	**
**  relatively easy to solve the vector equations involving	**
**  them.  The only further complexity in this algorithm is	**
**  that the input |b> must be first permuted to |b'> so that	**
**  its element ordering match that of A'.  Following the	**
**  solution in which |x'> is obtained it is then un-permuted	**
**  to match the original A ordering.  The first equation	**
**  to be solved (actually involving L' not L) appears as	**
**								**
** 		 [ 1   0   0   0 ] [y1] = [b1]			**
** 		 [L21  1   0   0 ] [y2] = |b2|			**
** 		 [L31 L32  1   0 ] [y3] = |b3|			**
** 		 [L41 L42 L43  1 ] [y4] = [b4]			**
**								**
**  because the diagonal elements of the matrix L are always 1.	**
**  The first element of |y> (actually |y'>) is given by	**
**								**
**			y1 = b1					**
**								**
**  and then subseqent elements of this vector are found by	**
**								**
**			      [	     ---	  ]		**
**		         1    |      \	 	  |		**
**	         yi = ------- | b  - /  <i|L|j>y  |		**
**		      <i|L|i> [  i   ---	j ]		**
**								**
*****************************************************************/

template<class T1,class T2> void lubksb_ip(Matrix<T1>& eqns,const Matrix<T2>& in,const BaseList<size_t> &indx)
{
  if (!issquare(in))
    throw NotSquare("lubksb_ip");
  if (issame(&eqns,&in))
    throw ArgumentClash("lubksb_ip");
  const size_t n=in.rows();  
  if (eqns.rows()!=n)
    throw Mismatch("lubksb_ip");
 
  size_t i,j,k;

  for (k=eqns.cols();k--;) {
    int ii=-1;
    for (i=0;i<n;i++) {
      const int ip=indx(i);
      T1 sum=eqns(ip,k);
      eqns(ip,k)=eqns(i,k);
      if (ii>=0)
	for (j=ii;j<i;j++)
	  sum-=in(i,j)*eqns(j,k);
      else
	if (sum!=0.0)
	  ii=i;
      eqns(i,k)=sum;
    }
    for (i=n;i--;) {
      T1 sum=eqns(i,k);
      for (j=i+1;j<n;j++)
	sum-=in(i,j)*eqns(j,k);
      eqns(i,k)=sum/in(i,i);
    }
  }
}

template<class T1,class T2> void lubksb_ip(BaseList<T1>& eqns,const Matrix<T2>& in,const BaseList<size_t>& indx)
{
  if (!issquare(in))
    throw NotSquare("lubksb_ip");
  const size_t n=in.rows();  
  if (eqns.length()!=n)
    throw Mismatch("lubksb_ip");
  
  size_t i,j;
  int ii=-1;
  for (i=0;i<n;i++) {
    const size_t ip=indx(i);
    T1 sum=eqns(ip);
    eqns(ip)=eqns(i);
    if (ii>=0) {
      for (j=ii;j<i;j++)
	sum-=in(i,j)*eqns(j);
    }
    else {
      if (sum!=0.0)
	ii=i;
    }
    eqns(i)=sum;
  }
  for (i=n;i--;) {
    T1 sum=eqns(i);
    for (j=i+1;j<n;j++)
      sum-=in(i,j)*eqns(j);
    eqns(i)=sum/in(i,i);
  }
}

template<class T> void inv_(Matrix<T>& d,const Matrix<T>& bv)
{
  Matrix<T> b(bv);
  if (!issquare(b))
    throw NotSquare("inv");
  const int n=b.rows();
  ScratchList<size_t> indx(n);
  if (!LU_decomp(b,indx))	// LU Decomposition of b
    throw SingularMatrix<T>("inv: singular matrix",b);

  d.identity(n);
  lubksb_ip(d,b,indx);
}

void inv(cmatrix& d,const cmatrix& b) { inv_(d,b); }  
void inv(Matrix<double>& d,const Matrix<double>& b) { inv_(d,b); }  

cmatrix inv(const cmatrix& a)
{
  cmatrix b(mxflag::temporary);
  inv_(b,a);
  return b;
}

Matrix<double> inv(const Matrix<double>& a)
{
  Matrix<double> b(mxflag::temporary);
  inv_(b,a);
  return b;
}

template<class T1,class T2> inline void multiply_inv_ip_(Matrix<T1>& a,Matrix<T2>& eqns)
{
  ScratchList<size_t> indx(a.rows());
  if (!LU_decomp(a,indx))
    throw SingularMatrix<T1>("multiply_inv_ip: matrix singular?",a);
  lubksb_ip(eqns,a,indx);
}

void multiply_inv_ip2(cmatrix& a,cmatrix& eqns)
{
  multiply_inv_ip_(a,eqns);
}

void multiply_inv_ip2(Matrix<double>& a,cmatrix& eqns)
{
  multiply_inv_ip_(a,eqns);
}

void multiply_inv_ip2(Matrix<double>& a,Matrix<double>& eqns)
{
  multiply_inv_ip_(a,eqns);
}
  
template<class T1,class T2> inline void multiply_inv_ip_(Matrix<T2>& a,BaseList<T1>& eqns)
{
  ScratchList<size_t> indx(a.rows());
  if (!LU_decomp(a,indx)) 
    throw SingularMatrix<T2>("multiply_inv: matrix singular?",a);
  lubksb_ip(eqns,a,indx);
}

void multiply_inv_ip2(cmatrix& a,BaseList<complex>& eqns)
{
  multiply_inv_ip_(a,eqns);
}

void multiply_inv_ip2(Matrix<double>& a,BaseList<complex>& eqns)
{
  multiply_inv_ip_(a,eqns);
}

void multiply_inv_ip2(Matrix<double>& a,BaseList<double>& eqns)
{
  multiply_inv_ip_(a,eqns);
}

Matrix<double> multiply_inv(const Matrix<double>& av,const Matrix<double>& eqns)
{
  Matrix<double> a(av,mxflag::normal);
  Matrix<double> res(eqns,mxflag::temporary);
  multiply_inv_ip_(a,res);
  return res;
}

cmatrix multiply_inv(const Matrix<double>& av,const cmatrix& eqns)
{
  Matrix<double> a(av,mxflag::normal);
  cmatrix res(eqns,mxflag::temporary);
  multiply_inv_ip_(a,res);
  return res;
}

cmatrix multiply_inv(const cmatrix& av,const cmatrix& eqns)
{ 
  cmatrix a(av,mxflag::normal);
  cmatrix res(eqns,mxflag::temporary);
  multiply_inv_ip_(a,res);
  return res;
}

List<double> multiply_inv(const Matrix<double>& av,const BaseList<double>& eqns)
{
  Matrix<double> a(av,mxflag::normal);
  List<double> res(eqns,mxflag::temporary);
  multiply_inv_ip_(a,res);
  return res;
}
 
List<complex> multiply_inv(const Matrix<double>& av,const BaseList<complex>& eqns)
{
  Matrix<double> a(av,mxflag::normal);
  List<complex> res(eqns,mxflag::temporary);
  multiply_inv_ip_(a,res);
  return res;
}
 
List<complex> multiply_inv(const cmatrix& av,const BaseList<complex>& eqns)
{
  cmatrix a(av,mxflag::normal);
  List<complex> res(eqns,mxflag::temporary);
  multiply_inv_ip_(a,res);
  return res;
}
 
void multiply_inv_ip(Matrix<double>& a,Matrix<double>& eqns)
{
  multiply_inv_ip_(a,eqns);
}
  
void multiply_inv_ip(Matrix<double>& a,cmatrix& eqns)
{
  multiply_inv_ip_(a,eqns);
}
  
void multiply_inv_ip(cmatrix& a,cmatrix& eqns)
{
  multiply_inv_ip_(a,eqns);
}
    
} //namespace libcmatrix
