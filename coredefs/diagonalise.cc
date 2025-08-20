/* Diagonalisation routines */
/* It is important to note that these only apply to normal matrices
   i.e. unitary or hermitian (orthogonal or symmetric for real matrices) matrices only */

/* Matrix diagonalisation routines from GAMMA */

#include "config.h"
#include "cmatrix_external.h"
#include <cstdlib>
#include <cmath>
#include "List.h"
#include "cmatrix.h"
#include "rmatrix.h"
#include "ScratchList.h"

#define LCM_LAPACK_CONVERGENCE_CHECK 1

#ifdef LCM_LAPACK_CONVERGENCE_CHECK
#include <float.h>
#endif

namespace libcmatrix {

  eigensystem_controller::eigensystem_controller()
    : tolerance(0.0),
      effectivezero(0.0),
      throwexception(false),
      verbose(0),
      chebyshev_iterations(LCM_CHEBYSHEV_ITERATIONS),
      useexternal(false),
      nonorthogonal_warning("eigensystem: eigenvectors are not properly orthogonal",&lcm_base_warning)
  {}

  eigensystem_controller cmatrix_eigensystem_controller; //adjust behaviour of general complex diagonalise

  Warning<> invalidmaxiter_warning("Failed to parse LCM_MAXITERATIONS_PER_EIGENVALUE as floating point number >0",&lcm_base_warning);

#ifdef LCM_LAPACK_CONVERGENCE_CHECK

  /* simpler version of dlamch for the case of IEEE754-compliant FPU module by Piotr Luszczek S.
     taken from http://www.mail-archive.com/numpy-discussion@lists.sourceforge.net/msg02448.html */
  
#ifndef DBL_DIGITS
#define DBL_DIGITS 53
#endif
  
  double dlamch(char cm)
  {
    static const double eps = DBL_EPSILON;
    static double smallv=0.0;
    
    switch (tolower(cm)) {
    case 'b':
      return FLT_RADIX;
    case 'e':
      return eps;
    case 'l':
      return DBL_MAX_EXP;
    case 'm':
      return DBL_MIN_EXP;	
    case 'n':
      return DBL_DIGITS;
    case 'o':
      return DBL_MAX;
    case 'p':
      return eps * FLT_RADIX;
    case 'r':
      return (FLT_ROUNDS < 2);
    case 's': {
      /* Use SMALL plus a bit, to avoid the possibility of rounding causing overflow
	 when computing  1/sfmin. */
      if (smallv==0.0) {
	double sfmin = DBL_MIN;	
	smallv = 2. / DBL_MAX;
	if (smallv <= sfmin)
	  smallv = sfmin * (1 + eps);
      }
      return smallv;
    }
    case 'u':
      return DBL_MIN;
    }
    throw InvalidParameter("dlamch: unknown argument");
  }   

  static const char CONVTYPE[]="LAPACK";
#else
  static const char CONVTYPE[]="NR";
#endif

  inline bool isconverged(double e, double v1, double v2)
  {
#ifdef LCM_LAPACK_CONVERGENCE_CHECK
    static const double safmin=dlamch('s');
    static const double eps=dlamch('e');
    static const double eps2=eps*eps;
    const double abse=fabs(e);
    return ( (abse*abse)<=(eps2*fabs(v1*v2)+safmin) );
#else
    const double dd=fabs(v1)+fabs(v2);
    double volatile temp=fabs(e)+dd;
    return (temp==dd);
#endif
  }

  namespace {
    const char LCM_ES[]="eigensystem";

    size_t get_max_iter(size_t n)
    {
      static double scale_fac=0.0;
      if (scale_fac==0.0) {
	const char* envval=getenv("LCM_MAXITERATIONS_PER_EIGENVALUE");
	if (envval && *envval) {
	  if ((sscanf(envval,"%lf",&scale_fac)!=1) || (scale_fac<=0.0))
	    invalidmaxiter_warning.raise();
	}
	if (scale_fac==0.0)
	  scale_fac=10;	
      }
      return size_t(0.5+ceil(scale_fac*n));
    }

    template<class T> void dump_matrix(const Matrix<T>& a, const char* name, std::ostream& ostr =std::cerr)
    {
      spy(ostr,a);
      ostr.flush();
      try {
	write_matrix(name,a); //trap errors
      } catch (std::exception& exc) {
	std::cerr << "dump_matrix failed: " << exc.what() << std::endl;
      }
    }
    
    inline double AbsNorm(const complex& z)
    { return fabs(real(z))+fabs(imag(z)); }
    
  }

/*
     this subroutine is a translation of a complex analogue of
     the algol procedure orthes, num. math. 12, 349-368(1968)
     by martin and wilkinson.
     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

     given a complex general matrix, this subroutine
     reduces a submatrix situated in rows and columns
     low through igh to upper hessenberg form by
     unitary similarity transformations.

     on input
        low and igh are integers determined by the balancing
          subroutine  cbal.  if  cbal  has not been used,
          set low=1, igh=n.

        a contains the complex input matrix.

     on output

        a of the hessenberg matrix.  information
          about the unitary transformations used in the reduction
          is stored in the remaining triangles under the
          hessenberg matrix.

        ort contains further information about the
          transformations.  only elements low through igh are used.

*/

static void corth(cmatrix& a, complex* ort, int low, int igh)
{
  int i,j,m;
  double f,g,h,sc;
  complex z;
  size_t n=a.rows();


  for (m=low+1; m<igh; m++) {
    h=0;
    ort[m]=0;
    sc=0;
    for (i=m; i<=igh; i++) 
      sc+=AbsNorm( a(i,m-1) );
    
    if (sc>0.0) {
      for (i=igh;i>=m;i--) {
	ort[i]= a(i,m-1) /sc;
	h+=norm(ort[i]);
      }
      g=std::sqrt(h);
      f=std::sqrt(norm(ort[m]));
      if (f==0.0) {
	ort[m]=g;
	real( a(m,m-1), sc);
      }
      else {
	h+=f*g;
	g/=f;
	ort[m]*=(1.0+g);
      }

      for (j=m;j<n;j++) {
	z=0.0;
	for (i=igh;i>=m;--i) 
	  mla_conj(z,a(i,j),ort[i]);
	z/=h;
	for (i=m;i<=igh;i++) 
	  a(i,j)-=z*ort[i];
      }

      for (i=0;i<=igh;i++) {
	z=0.0;
	for (j=igh;j>=m;--j) 
	  mla(z,ort[j],a(i,j));
	z/=h;
	for (j=m;j<=igh;j++) 
	  a(i,j)-= conj_multiply(ort[j],z);
      }

      ort[m]*=sc;
      a(m,m-1)*= -g;
    }
  }
}

/*
     this subroutine is a translation of a unitary analogue of the
     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
     and wilkinson.
     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
     the unitary analogue substitutes the qr algorithm of francis
     (comp. jour. 4, 332-345(1962)) for the lr algorithm.

     this subroutine finds the eigenvalues and eigenvectors
     of a complex upper hessenberg matrix by the qr
     method.  the eigenvectors of a complex general matrix
     can also be found if  corth  has been used to reduce
     this general matrix to hessenberg form.

     on input

        low and igh are integers determined by the balancing
          subroutine  cbal.  if  cbal  has not been used,
          set low=1, igh=n.

        ort contains information about the unitary trans-
          formations used in the reduction by  corth, if performed.
          only elements low through igh are used.  if the eigenvectors
          of the hessenberg matrix are desired, set ortr(j) and
          orti(j) to 0.0d0 for these elements.

        h contains the complex upper hessenberg matrix.
          their lower triangles below the subdiagonal contain further
          information about the transformations which were used in the
          reduction by  corth, if performed.  if the eigenvectors of
          the hessenberg matrix are desired, these elements may be
          arbitrary.

     on output

        ort, and the upper hessenberg portions of h
          have been destroyed.

        w contains the eigenvalues.  if an error
          exit is made, the eigenvalues should be correct
          for indices ierr+1,...,n.

        z comtains the eigenvectors.  the eigenvectors
          are unnormalized.  if an error exit is made, none of
          the eigenvectors has been found.

        ierr is set to
          zero       for normal return,
          j          if the limit of 30*n iterations is exhausted
                     while the j-th eigenvalue is being sought.

*/

/* pH: Some fluff removed, but code needs origin to be properly shifted to
C-style 0. Very tricky without messing up somewhere ... */

static int comqr3(cmatrix& a, int low, int igh,
		     complex* ort, BaseList<complex> w, cmatrix& z,
		     int flag)
{
  size_t n=a.rows();
  int i,j,k,l=0,m,en,ii,jj,ll,nn,its,iend,enm1;
  double u,s,ss,nor;
  complex sc,x,y,t,zz;

  const double machep=1.0e-13;

  z.identity(n);

  /* accumulate transformations */
  iend=igh-low-1;
  if (iend>=0) {
      for (ii=1;ii<=iend;ii++) {
	  i=igh-ii;
	  if (ort[i-1]!=0.0) {
	      if (a(i-1,i-2)!=0.0) {
		nor=real(conj_multiply(a(i-1,i-2),ort[i-1]));
		//nor=real(a(i-1,i-2))*real(ort[i-1])+imag(a(i-1,i-2))*imag(ort[i-1]);
		  for (k=i+1;k<=igh;k++) 
		    ort[k-1]=a(k-1,i-2);
		  for (j=i;j<=igh;j++) {
		      sc=0.0;
		      for (k=i;k<=igh;k++) 
			//mla_conj(sc,z(k-1,j-1),ort[k-1]);
				 sc+=conj_multiply(ort[k-1],z(k-1,j-1));
		      sc/=nor;
		      for (k=i;k<=igh;k++) 
			z(k-1,j-1)+=sc*ort[k-1];
		    }
		}
	    }
	}
      /* real subdiagonal elements */
      l=low+1;
      for (i=l;i<=igh;i++) {
	ll= (i+1<igh) ? i+1 : igh;
	if (imag(a(i-1,i-2))!=0.0) {
	  nor=std::sqrt(norm(a(i-1,i-2)));
	  y=a(i-1,i-2)/nor;
	  a(i-1,i-2)=nor;
	  for (j=i;j<=n;j++) 
	    a(i-1,j-1)=conj_multiply(y,a(i-1,j-1));
	  for (j=1;j<=ll;j++) 
	    a(j-1,i-1)*=y;
	  for (j=low;j<=igh;j++) 
	    z(j-1,i-1)*=y;
	}
      }
    }

  /* isolated roots */
  for (i=1;i<=n;i++)
    if ((i<low)||(i>igh))
      w(i-1)=a(i-1,i-1);

  en=igh;
  t=0.0;

  /* iterate for eigenvalues */
  while (en>=low) {
      its=0;
      enm1=en-1;
      while (1) {
	  for (ll=low;ll<=en;ll++) {
	      l=en+low-ll;
	      if (l==low)
		break;
	      s=AbsNorm(a(l-2,l-2))+AbsNorm(a(l-1,l-1));
	      if (s<1.0)
		s=1.0;
	      if (fabs(real(a(l-1,l-2)))<=machep*s)
		break;
	  }
	  if (l==en)
	    break;
	  if (its==30)
	    throw Failed(LCM_ES);
	  /* form shift */
	  if ((its==10)||(its==20)) 
	    sc=fabs(real(a(en-1,enm1-1))) 
	      + fabs(real(a(enm1-1,en-3)));
	  else {
	      sc=a(en-1,en-1);
	      x=a(enm1-1,en-1)*real(a(en-1,enm1-1));
	      if (x!=0.0) {
		  y=(a(enm1-1,enm1-1)-sc)/2.0;
		  zz = sqrt(y*y+x);
		  if (real(y)*real(zz)+imag(y)*imag(zz)<0.0) 
		    zz = -zz;
		  zz = x/(y+zz);
		  sc-=zz;
	      }
	  }
	  for (i=low;i<=en;i++) 
	    a(i-1,i-1)-=sc;
	  t+=sc;
	  its+=1;
	  /* QR decomposition */
	  for (i=l+1;i<=en;i++) {
	      s=real(a(i-1,i-2));
	      nor=std::sqrt( norm(a(i-2,i-2))+s*s );
	      x=a(i-2,i-2)/nor;
	      w(i-2)=x;
	      a(i-2,i-2)=nor;
	      s /= nor;
	      a(i-1,i-2)=complex(0,s);
	      for (j=i;j<=n;j++) {
		y=a(i-2,j-1);
		zz=a(i-1,j-1);
		a(i-2,j-1)= conj_multiply(x,y)+s*zz;
		a(i-1,j-1)= x*zz - s*y;
	      }
	  }
	  //s=imag(a(en-1,en-1));
	  //isreal=(fabs(imag(a(en-1,en-1)))<machep);
	  const bool isreal=(imag(a(en-1,en-1))==0.0);
	  if (!isreal) {
	    nor=std::sqrt(norm(a(en-1,en-1)));
	    sc=conj(a(en-1,en-1)/nor);
	    a(en-1,en-1)=nor;
	    if (en!=n) 
	      for (j=en+1;j<=n;j++)
		a(en-1,j-1)*=sc;
	  }
	  /* calculate RQ */
	  for (j=l+1;j<=en;j++) {
	    x=w(j-2);
	    ss=imag(a(j-1,j-2));
	    const complex cx=conj(x);
	    for (i=1;i<=j;i++) {
	      y=real(a(i-1,j-2));
	      zz=a((i-1),j-1);
	      if (i!=j) {
		imag(y, imag(a(i-1,j-2)));
		imag(a(i-1,j-2), real(x)*imag(y)+imag(x)*real(y)+ss*imag(zz));
	      }
	      real( a(i-1,j-2), real(x)*real(y)-imag(x)*imag(y)+ss*real(zz));
	      a(i-1,j-1)=cx*zz-ss*y;
	    }
	    for (i=low;i<=igh;i++) {
	      y=z(i-1,j-2);
	      zz=z(i-1,j-1);
	      z(i-1,j-2)=x*y+ss*zz;
	      z(i-1,j-1)=cx*zz-ss*y;
	    }
	  }
	  sc=conj(sc);
	  if (!isreal) {
	    for (i=1;i<=en;i++) 
	      a(i-1,en-1)*=sc;
	    for (i=low;i<=igh;i++) 
	      z(i-1,en-1)*=sc;
	  }
      }
      /* a root found */
      a(en-1,en-1)+=t;
      w(en-1)=a(en-1,en-1);
      en=enm1;
  }

  if (flag) {
      /* all roots found, backsubstitution of eigenvectors */
      nor=0.0;
      for (i=1;i<=n;i++) 
	for (j=i;j<=n;j++) 
	  nor+=AbsNorm(a(i-1,j-1));
      if ((n==1)||(nor==0.0))
	return 0;
      for (nn=2;nn<=n;nn++) {
	en=n+2-nn;
	x=w(en-1);
	a(en-1,en-1)=1;
	enm1=en-1;
	for (ii=1;ii<=enm1;ii++) {
	  i=en-ii;
	  zz=a(i-1,en-1);
	  if (i!=enm1) {
	    for (j=i+1;j<=enm1;j++) 
	      zz+= a(i-1,j-1)*a(j-1,en-1);
	  }
	  y=x-w(i-1);
	  if ((fabs(real(y))<machep)&&(fabs(imag(y))<machep)) 
	    real( y, machep*nor);
	  a(i-1,en-1)=zz/y;
	}
      }
  }
  /* eigenvectors of isolated root */
  for (i=1;i<=n-1;i++) {
    if ((i<low)||(i>igh)) 
      for (j=i+1;j<=n;j++) 
	z((i-1),j-1)=a((i-1),j-1);
  }
  /* multiply for eigenbase */
  for (jj=low;jj<=n-1;jj++) {
    j=n+low-jj;
    if (j-1<igh) m=j-1;
    else m=igh;
    for (i=low;i<=igh;i++) {
      zz=z(i-1,j-1);
      for (k=low;k<=m;k++) 
	zz+=z(i-1,k-1)*a(k-1,j-1);
      z(i-1,j-1)=zz;
    }
  }
  if (flag) {
    /* normalize vectors */
    for (i=0;i<n;i++) {
      u=0.0;
      for (j=0;j<n;j++) 
	u+=AbsNorm(z(j,i));
      s=0.0;
      for (j=0;j<n;j++) 
	s+=norm(z(j,i)/u);
      s=u*std::sqrt(s);
      for (j=0;j<n;j++) 
	z(j,i)/=s;
    }
  }
  return 0;
}


/********************************************************
*							*
*	The routine diag combines the routines cred,	*
*	rred and tqli ( using the function sign )	*
*	to diagonalize a complex Hermitian matrix.	*
*       It uses the routines corth and comqr3 to        *
*       diagonalise general complex matrices            *
*							*
*	During the course of the transformations	*
*	exclusively unitary ( orthogonal ) transfor-	*
*	mations are used. Errors come mainly from	*
*	rounding. The resulting eigenbase deviates from	*
*	unitarity only negligibly.			*
*							*
********************************************************/

  void eigensystem(cmatrix& d,List<complex>& vals,const cmatrix& a, cmatrix* tmp)
{
  vals.create(a.rows());
  eigensystem(d,static_cast< BaseList<complex>& >(vals),a,tmp);
}
  
  void cleanzero_ip(BaseList<complex> vec, double tol)
  {
    if (tol==0.0)
      return;
    const complex zero(0.0);
    const double tol2(tol*tol);
    for (size_t i=vec.size();i--;) {
      complex& z(vec(i));
      if (norm(z)<tol2)
	z=zero;
    }
  }

  void cleanzero_ip(BaseList<double> vec, double tol)
  {
    if (tol==0.0)
      return;
    for (size_t i=vec.size();i--;) {
      double& v(vec(i));
      if (fabs(v)<tol)
	v=0.0;
    }
  }

namespace {

  void checkout(const cmatrix& d, const BaseList<complex>& vals, const cmatrix& av)
  {
    if (cmatrix_eigensystem_controller.tolerance==0.0)
      return;
    cmatrix checkm;
    conj_transpose_multiply(checkm,d,d);
    if (hasoffdiagonal(checkm,cmatrix_eigensystem_controller.tolerance)) {
      if (cmatrix_eigensystem_controller.throwexception) {
	if (cmatrix_eigensystem_controller.verbose) {
	  dump_matrix(av,"eigensystemInput");
	  dump_matrix(checkm,"eigensystemOrthTest");
	}
	cmatrix_eigensystem_controller.nonorthogonal_warning.raiseas(BaseWarning::RaiseException);
      }
      else {
	cmatrix_eigensystem_controller.nonorthogonal_warning.raise();
	if (cmatrix_eigensystem_controller.verbose) {
	  std::cerr << "Input matrix\n" << av;
	  std::cerr << "V'*V:\n" << checkm;
	  std::cerr << "Eigenvalues: " << vals << '\n';
	}
      }
    }
    else {
      if (cmatrix_eigensystem_controller.verbose>1)
	std::cout << "Passed eigenvector orthogonality check\n";
    }
  }

  size_t checkin(const cmatrix& av, const BaseList<complex>& vals)
  {    
    if (!issquare(av))
      throw NotSquare(LCM_ES);
    const size_t n=av.rows();
    if (vals.length()!=n)
      throw Mismatch(LCM_ES);
    return n;
  }
}

  void internal_eigensystem(cmatrix& d,BaseList<complex> vals,const cmatrix& av, cmatrix* wspace)
{
  const size_t n=checkin(av,vals);

  cmatrix a_tmp;
  cmatrix* ap(wspace ? wspace : &a_tmp);
  *ap=av; //make copy (aliasing of av and wspace will be detected)

  //optionally cleanup not-quite-zero entries
  cleanzero_ip(ap->row(),cmatrix_eigensystem_controller.effectivezero);

  ScratchList<complex> tmp(n,complex(0.0));
  corth(*ap,tmp.vector(),0,n-1);	// To upper Hessenberg form
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "eigensystem: Upper Hessenberg form:\n" << (*ap) << std::endl;
  comqr3(*ap,1, n, tmp.vector(),vals,d,1); //!< To diagonal form
  checkout(d,vals,av);
}

#ifdef LCM_F2C_SETUP

static char errm[256];

void lapack_eigensystem(cmatrix& V, BaseList<complex> eigs, const cmatrix& a0, cmatrix* wspace)
{
  integer n=checkin(a0,eigs);
  integer info=0;

  cmatrix a_tmp;
  cmatrix& a(wspace ? *wspace : a_tmp);
  transpose(a,a0);
  V.create(n,n);
  
#ifdef USE_SUNPERFACML
  zgeev('N','V',n,reinterpret_cast<complex*>(a.vector()),n,reinterpret_cast<complex*>(eigs.vector()),(complex*)NULL,1,reinterpret_cast<complex*>(V.vector()),n,&info);
#else
  complex worktmp=-1000.0;
  ScratchList<double> rwork(2*n);
  integer one=1;
  integer minusone=-1;
  static integer lastn=-1;
  integer lwork=-1;

  if (lastn!=n) {
    zgeev_("N","V",&n,reinterpret_cast<complex*>(a.vector()),&n,reinterpret_cast<complex*>(eigs.vector()),(complex*)NULL,&one,reinterpret_cast<complex*>(V.vector()),&n,(complex *)(&worktmp),&minusone,rwork.vector(),&info);

    lwork=int(0.5+libcmatrix::real(worktmp));
    lastn=n;
  }
  ScratchList<complex> work(lwork);

  ::zgeev_("N","V",&n,reinterpret_cast<complex*>(a.vector()),&n,reinterpret_cast<complex*>(eigs.vector()),(complex*)NULL,&one,reinterpret_cast<complex*>(V.vector()),&n,work.vector(),&lwork,lapack_pass(rwork),&info);
#endif

  if (info) {
    if (info>0)
      sprintf(errm,"zgeev: failed to find %i eigenvalues",info);
    else
      sprintf(errm,"zgeev: error in parameter %i",-info);
    throw Failed(errm);
  }

  V.transpose();
  checkout(V,eigs,a0);
}
#endif

void eigensystem(cmatrix& d,BaseList<complex> vals,const cmatrix& av, cmatrix* wspace)
{
#ifdef HAVE_LIBLAPACK
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "Complex general matrix diagonalisation using: " << (cmatrix_eigensystem_controller.useexternal ? "LAPACK" : "internal") << " routine\n";
  if (cmatrix_eigensystem_controller.useexternal)
    lapack_eigensystem(d,vals,av,wspace);
  else
#endif
    internal_eigensystem(d,vals,av,wspace);
}

  bool eigensystem(cmatrix& D, BaseList<complex> cycle_eigs, const cmatrix& U, double ptol, cmatrix* tmp)
{
  if (!issquare(U))
    throw NotSquare(LCM_ES);
  const size_t n=U.rows();
  if (cycle_eigs.length()!=n)
    throw Mismatch(LCM_ES);
  if (n==1) {
    D.create(1,1,complex(1.0));
    cycle_eigs(0)=U(0,0);
    return false;
  }

  if (!hasoffdiagonal(U,ptol)) {
    diag(cycle_eigs,U);
    D.identity(U.rows());
    return true;
  }
  eigensystem(D,cycle_eigs,U,tmp);
  return false;
}

/* Diagonalisation of Hermitian matrices */

/********************************************************
*							*
*	The routine cred reduces a complex Hermitian	*
*	matrix to complex Hermitian tridiagonal form	*
*	using the Housholder algorithm.			*
*							*
********************************************************/

static void cred(cmatrix& a, cmatrix& z)
{
  int i,j,l;
  int n=a.rows();
  double s;			// norm of the vector a(l+1,l) to a(n-1,l)
  double sw;			// norm of the vextor a(l+2,l) to a(n-1,l)
  double t;
  
  ScratchList<complex> w(n);

  /* create unit matrix */
  z.identity(n);

  /* Housholder reduction */
  for (l=0;l<n-2;l++) {

    double u=0.0;			// column scaling factor
    for (i=l+2;i<n;i++) 
      u+=AbsNorm(a(i,l));

    if (u>0.0) {
      const complex &all=a(l+1,l);
      u+=AbsNorm(all);
   
      sw=0.0;
      for (i=l+2;i<n;i++) 
	sw+=norm(a(i,l)/u);
      s=sw+norm(all/u);
      s=u*std::sqrt(s);
      if (all!=0.0) {
	t=std::sqrt(norm(all));
	w(l+1)=all*(1.0+s/t);
      }
      else 
	w(l+1)= -s;
      s=sw+norm(w(l+1)/u);
      s=std::sqrt(2/s)/u;
      
      for (i=l+2;i<n;i++) 
	w(i)=a(i,l)*s;
      w(l+1)*=s;

      complex ss;
      /* transform matrix from the left */
      for (i=l;i<n;i++) {
	ss=0.0;
	for (j=l+1;j<n;j++)
	  mla_conj(ss,a(j,i),w(j));
	ss=-ss;
	  
	for (j=l+1;j<n;j++)
	  mla(a(j,i),ss,w(j));
	  //  a(j,i)-=ss*w(j);
      }
      /* transform matrix from the right */
      for (i=l;i<n;i++) {
	BaseList<complex> ar=a.row(i);
	ss=0.0;
	for (j=l+1;j<n;j++)
	  mla(ss,w(j),ar(j));
	ss=-ss;
	for (j=l+1;j<n;j++)
	  mla_conj(ar(j),ss,w(j));
	  //	  a(i,j)-=conjmul(w(j),ss);
      }
      
      /* accumulate transformations */
      for (i=0;i<n;i++) {
	BaseList<complex> zr=z.row(i);
	ss=0.0;
	for (j=l+1;j<n;j++)
	  mla(ss,w(j),zr(j));
	ss=-ss;
	for (j=l+1;j<n;j++)
	  mla_conj(zr(j),ss,w(j));
      }
    }
  }
}


/********************************************************
*							*
*	The routine rred reduces a complex Hermitian	*
*	tridiagonal matrix to real symmetric		*
*	tridiagonal form using a diagonal unitary	*
*	transformation.					*
*							*
********************************************************/

static void rred(cmatrix& a, cmatrix& eigvectors)
{
  int i,l;
  size_t n=a.rows();
  double s;
  complex ss;

  for (l=1;l<n;l++) {
    s=std::sqrt(norm(a(l,l-1)));
    if (s!=0) {		// if there is something to be done
      ss=conj(a(l,l-1))/s;
      a(l,l-1)=a(l-1,l)=s;
      if (l<n-1) {		// transform
        a(l,l+1)*=ss;
        a(l+1,l)*=conj(ss);
      }
      for (i=0;i<n;i++)	// accumulate transform
        eigvectors(i,l)*=conj(ss);
    }
  }
}

/********************************************************
*							*
*	The function sign is the C++ implementation	*
*	of the standard FORTRAN library function.	*
*							*
********************************************************/
inline double sign(double a,double b) {return (b<0.0)?-fabs(a):fabs(a);}

/********************************************************
*							*
*	The routine tqli uses the QL algorithm for the	*
*	input real symmetric tridiagonal matrix.	*
*							*
********************************************************/

static void tqli(const cmatrix& a, cmatrix& z, double* vallist)
{
  int i,k,l,m,iter;
  size_t n=a.rows();
  double s,r,p,g,f,c,b;

  double x;
  complex xx,yy;

  const size_t maxiter=get_max_iter(n);

  cmatrix zback;
  if (cmatrix_eigensystem_controller.verbose)
    zback=z;

  ScratchList<double> E(n);
  double* e=E.vector();
 
  for (i=0;i<n-1;i++)	// copy to workspace
    e[i]=real(a(i,i+1));

  for (i=n;i--;)
    vallist[i]=real(a(i,i));

  e[n-1]=0.0;
  iter=0;
  for (l=0;l<n;l++) {		// loop for eigenvalues
    do {

      // This loop looks crazy...
      m=l;
      while (1) {		// look for an eigenvalue
	if (m==n-1)
	  break;	// all eigenvalues found
	if (isconverged(e[m],vallist[m],vallist[m+1]))
	    break;
	m++;
      }
      if (m!=l) {
	if (iter>maxiter) {
	  char scratch[256];
	  if (cmatrix_eigensystem_controller.verbose)
	    dump_matrix(zback,"tqliDump");
	  snprintf(scratch,sizeof(scratch),"tqli (complex) failed after %" LCM_PRI_SIZE_T_MODIFIER "u iterations",maxiter);
	  throw Failed(scratch);
	// this happens almost surely only on erroneous input
	}
	
	iter++;
	// form shift
	g=(vallist[l+1]-vallist[l])/(2.0*e[l]);
	r=std::sqrt(g*g+1.0);
	x=e[l]/(g+sign(r,g));
	g=vallist[m]-vallist[l]+x;
	s=1.0;
	c=1.0;
	p=0.0;
	// QL step
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  if (fabs(f)>=fabs(g)) {
	    c=g/f;
	    r=std::sqrt(c*c+1.0);
	    e[i+1]=f*r;
	    s=1.0/r;
	    c*=s;
	  }
	  else {
	    s=f/g;
	    r=std::sqrt(s*s+1.0);
	    e[i+1]=g*r;
	    c=1.0/r;
	    s*=c;
	  }
	  g=vallist[i+1]-p;
	  r=(vallist[i]-g)*s+2.0*c*b;
	  p=s*r;
	  vallist[i+1]=g+p;
	  g=c*r-b;
	  // accumulate transforms
	  for (k=n;k--;) {
	    complex& zi =z(k,i);
	    complex& zi1=z(k,i+1);
	    xx=zi; xx*=s;
	    yy=zi1; yy*=s;
	    //	    xx=s*zi;
	    zi*=c; zi-=yy;
		    //	    zi=c*zi-s*zi1;
	    zi1*=c; zi1+=xx;
		    // 	    zi1=xx+c*zi1;
	  }
	}
	vallist[l]-=p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m!=l);		// all eigenvalues found
  }
    if (cmatrix_eigensystem_controller.verbose)
      printf("tqli (complex) complete in %i iterations (convergence check: %s)\n",iter,CONVTYPE);
}

/* Explicit diagonalisation for 2x2 hermitian matrices */

static void hermitian_eigensystem2(cmatrix& V,double eigs[2],const cmatrix& M)
{
  V.create(2,2);

  const double a=real(M(0,0));
  const double b=real(M(1,1));

  double sum=(a+b)/2;
  double z=(a-b)/2;

  const complex &c=M(0,1);

  complex phase(1.0);
  double r;

  if (imag(c)) {
    r=std::sqrt(norm(c));
    phase=c/r;
  }
  else
    r=real(c);
    
  double theta=M_PI/4;
  double lambda;

  if (z==0)
    lambda=r;
  else {
    theta=atan2(r,z)/2;
    lambda=z*std::cos(2*theta)+r*std::sin(2*theta);
  }

  double ctheta=std::cos(theta);
  double stheta=std::sin(theta);
  
  V(1,0)=stheta;
  V(1,1)=ctheta;
  V(0,0)=ctheta*phase;
  V(0,1)=-stheta*phase;

  eigs[0]=sum+lambda;
  eigs[1]=sum-lambda;
}

const char HE[]="hermitian_eigensystem";

  void hermitian_eigensystem(cmatrix& d,List<double>& vals,const cmatrix& a, cmatrix* tmp)
{
  if (!a)
    throw Undefined(HE);
  vals.create(a.rows());
  hermitian_eigensystem(d,static_cast< BaseList<double> >(vals),a,tmp);
}
  
  void hermitian_eigensystem(cmatrix& d,BaseList<double> vals,const cmatrix& av, cmatrix* tmp)
{
  if (!issquare(av))
    throw NotSquare(HE);
  const size_t n=av.rows();
  if (n!=vals.length())
    throw Mismatch(HE);

  if (n==2) {
    hermitian_eigensystem2(d,vals.vector(),av);
    return;
  }
  
  cmatrix a_tmp;
  cmatrix* ap(tmp ? tmp : &a_tmp);
  *ap=av; //make copy (aliasing of av and tmp will be detected)

  d.create(n,n);

  cred(*ap,d);		// To Hermitian tridiagonal form
  rred(*ap,d);		// To real symm. tridiagonal form
  tqli(*ap,d,vals.vector());	// To diagonal form
}

/* Diagonalisation of real symmetric matrices, adapted from Numerical Recipes */

//static const size_t MAX_ITER=30;

inline double norm(double a) { return a*a; }

static void tred2(rmatrix& A,BaseList<double>& D, BaseList<double>& E)
{
  if (!issquare(A))
    throw NotSquare("tred2");
  const size_t n=A.rows();
  if (D.length()!=n || E.length()!=n)
    throw Mismatch("tred2");
  if (n<2)
    throw InvalidParameter("tred2");

  int i,j,k;

  for (i=n-1;i>=1;i--) {
    const int l=i-1;
    double h=0.0;

    if (l>0) {
      double scale=0.0;
      //      for (k=i;k--;) scale+=fabs(A(i,k));
      for (k=l;k--;)
	scale+=fabs(A(i,k));

      if (scale==0.0)
	E(i)=A(i,l);
      else {
	scale+=fabs(A(i,l));
	for (k=0;k<i;k++) {
	  A(i,k)/=scale;
	  h+=norm(A(i,k));
	}
	double f=A(i,l);
	double g=-sign(std::sqrt(h),f);
	E(i)=scale*g;
	h-=f*g;
	A(i,l)=f-g;
	f=0.0;
	for (j=0;j<i;j++) {
	  A(j,i)=A(i,j)/h;
	  g=0.0;
	  for (k=0;k<=j;k++)
	    g+=A(j,k)*A(i,k);
	  if (l>j) {
	    for (k=j+1;k<i;k++)
	      g+=A(k,j)*A(i,k);
	  }
	  E(j)=g/h;
	  f+=E(j)*A(i,j);
	}
	const double hh=f/(h+h);
	for (j=0;j<i;j++) {
	  f=A(i,j);
	  g=E(j)-hh*f;
	  E(j)=g;
	  for (k=0;k<=j;k++)
	    A(j,k)-=f*E(k)+g*A(i,k);
	}
      }
    }
    else
      E(i)=A(i,l);
    D(i)=h;
  }
  D(0)=0.0;
  E(0)=0.0;

  for (i=0;i<n;i++) {
    if (D(i)!=0.0) {
      for (j=0;j<i;j++) {
	double g=0.0;
	for (k=0;k<i;k++)
	  g+=A(k,j)*A(i,k);
	for (k=0;k<i;k++)
	  A(k,j)-=g*A(k,i);
      }
    }
    D(i)=A(i,i);
    A(i,i)=1.0;
    for (j=0;j<i;j++)
      A(i,j)=A(j,i)=0.0;
  }
}

static void tqli(BaseList<double>& D, BaseList<double>& E, rmatrix& Z)
{
  const size_t n=D.length();
  if (!issquare(Z))
    throw NotSquare("tqli");
  if (Z.rows()!=n || E.length()!=n)
    throw Mismatch("tqli");

  if (n<2)
    return;

  rmatrix Zback;
  if (cmatrix_eigensystem_controller.verbose)
    Zback=Z;

  const size_t maxiter=get_max_iter(n);

  int i,k,l,m;
  int iter=0;
  
  for (i=1;i<n;i++)
    E(i-1)=E(i);
  E(n-1)=0.0;

  for (l=0;l<n;l++) {
    for (;;iter++) {
      for (m=l;m<n-1;m++) {
	if (isconverged(E(m),D(m),D(m+1)))
	  break;
	    //	double DD=fabs(D(m))+fabs(D(m+1));
	    //	double volatile temp=fabs(E(m))+DD;
	    //if (temp==DD)
	    //  break;
      }
      if (m==l)
	break;

      if (iter>maxiter) {
	char scratch[256];
	if (cmatrix_eigensystem_controller.verbose)
	  dump_matrix(Zback,"tqliDump");
	snprintf(scratch,sizeof(scratch),"tqli (real) failed after %" LCM_PRI_SIZE_T_MODIFIER "u iterations",maxiter);
	throw Failed(scratch);
	// this happens almost surely only on erroneous input
      }
      
      double g=(D(l+1)-D(l))/(2*E(l));
      double r=std::sqrt(g*g+1.0);
      g=D(m)-D(l)+E(l)/(g+sign(r,g));
      double s=1.0;
      double c=1.0;
      double p=0.0;
      
      for (i=m-1;i>=l;i--) {
	double f=s*E(i);
	double b=c*E(i);
	if (fabs(f)>=fabs(g)) {
	  c=g/f;
	  r=std::sqrt(c*c+1.0);
	  E(i+1)=f*r;
	  s=1.0/r;
	  c*=s;
	}
	else {
	  s=f/g;
	  r=std::sqrt(s*s+1.0);
	  E(i+1)=g*r;
	  c=1.0/r;
	  s*=c;
	}
	g=D(i+1)-p;
	r=(D(i)-g)*s+2*c*b;
	p=s*r;
	D(i+1)=g+p;
	g=c*r-b;
	for (k=0;k<n;k++) {
	  double* const ZMK=Z.vector(k);
	  const double x=ZMK[i+1]; //Z(k,i+1);
	  const double y=ZMK[i]; //Z(k,i);
	  ZMK[i+1]=s*y+c*x;
	  ZMK[i]=c*y-s*x;
	}
      }
      D(l)-=p;
      E(l)=g;
      E(m)=0.0;
    }
  }
  if (cmatrix_eigensystem_controller.verbose)
    printf("tqli (complex) complete in %i iterations (convergence check: %s)\n",iter,CONVTYPE);  
}

/* Explicit diagonalisation for 2x2 real symmetric matrices */

static void hermitian_eigensystem2(Matrix<double>& V,double eigs[2],const Matrix<double>& M)
{
  V.create(2,2);
  //  const double sum=(M[0][0]+M[1][1])/2;
  const double h=(M(0,0)-M(1,1));
  const double c=M(0,1);
    
  double t;

  if (h+100.0*fabs(c)==h) {// off-diagonal is zero to machine precn 
    // special case for h=0
    if (h==0.0) {
      eigs[0]=M(0,0);
      eigs[1]=M(1,1);
      V(1,0)=V(0,1)=0.0;
      V(1,1)=V(0,0)=1.0;
      return;
    }
    t=c/h;
  }
  else {
    const double theta=0.5*h/c;
    t=1/(fabs(theta)+std::sqrt(1.0+theta*theta));
    if (theta<0.0) t=-t;
  }
  const double ctheta=1/std::sqrt(1.0+t*t);
  const double stheta=t*ctheta;
  
  V(1,0)=-(V(0,1)=ctheta);
  V(1,1)=V(0,0)=stheta;

  const double diff=t*c;

  eigs[0]=M(1,1)-diff;
  eigs[1]=M(0,0)+diff;
}

void hermitian_eigensystem(rmatrix& V,BaseList<double> evals,const rmatrix& a)
{
  if (!issquare(a))
    throw NotSquare(HE);
  const size_t n=a.rows();
  if (evals.length()!=n)
    throw Mismatch(HE);

  switch (n) {
  case 1:
    evals(0)=a(0,0);
    V.create(n,n); V=1.0;
    return;
  case 2:
    hermitian_eigensystem2(V,evals.vector(),a);
    return;
  }
  
  ScratchList<double> E(n);

  V=a;

  //optionally cleanup not-quite-zero entries
  cleanzero_ip(V.row(),cmatrix_eigensystem_controller.effectivezero);

  tred2(V,evals,E);

  cleanzero_ip(V.row(),cmatrix_eigensystem_controller.effectivezero);
  tqli(evals,E,V);
}

void hermitian_eigensystem(rmatrix& V,List<double>& E,const rmatrix& A)
{
  E.create(A.rows());
  hermitian_eigensystem(V,static_cast< BaseList<double> >(E),A);
}

namespace {
  inline double labs(double x) { return fabs(x); }
  inline double labs(const complex& z) { return std::sqrt(norm(z)); }
  inline double lreal(double x) { return x; }
  inline double lreal(const complex& z) { return z.real(); }
}

//! return eigenvalue range for matrix with real eigenvalues using Gershgorin's circle theorem
template<class T> std::pair<double,double> gershgorin_minmax_(const Matrix<T>& H) 
{
  double min=1e30;
  double max=-1e30;
  if (!issquare(H))
    throw NotSquare("max_eigenvalue");
  const size_t n=H.rows();
  for (size_t i=0;i<n;i++) {
    double rad=0.0;
    for (size_t j=n;j--;) {
      if (i!=j)
	rad+=labs(H(i,j));
    }
    const double centre=lreal(H(i,i));
    if ((i==0) || ((centre+rad)>max))
      max=centre+rad;
    if ((i==0) || ((centre-rad)<min))
      min=centre-rad;
  }
  return std::pair<double,double>(min,max);
}

std::pair<double,double> hermitian_eigenvalue_range(const cmatrix& a)
{ return gershgorin_minmax_(a); }

std::pair<double,double> eigenvalue_range(const rmatrix& a)
{ return gershgorin_minmax_(a); }

}//namespace libcmatrix
