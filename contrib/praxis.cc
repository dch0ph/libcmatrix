//#include <limits>
#include "rmatrix.h"
#include "optim.h"

/*   Comment from praxis.f

     praxis returns the minimum of the function f(x,n) of n variables
     using the principal axis method.  the gradient of the function is
     not required.

     for a description of the algorithm, see chapter seven of
     "Algorithms for finding zeros and extrema of functions without
     calculating derivatives" by Richard P Brent.

     the parameters are:
     t0       is a tolerance.  praxis attempts to return praxis=f(x)
              such that if x0 is the true local minimum near x, then
              norm(x-x0) < t0 + squareroot(machep)*norm(x).
     machep   is the machine precision, the smallest number such that
              1 + machep > 1.  machep should be 16.**-13 (about
              2.22d-16) for real*8 arithmetic on the ibm 360.
     h0       is the maximum step size.  h0 should be set to about the
              maximum distance from the initial guess to the minimum.
              (if h0 is set too large or too small, the initial rate of
              convergence may be slow.)
     n        (at least two) is the number of variables upon which
              the function depends.
     prin     controls the printing of intermediate results.
              if prin=0, nothing is printed.
              if prin=1, f is printed after every n+1 or n+2 linear
              minimizations.  final x is printed, but intermediate x is
              printed only if n is at most 4.
              if prin=2, the scale factors and the principal values of
              the approximating quadratic form are also printed.
              if prin=3, x is also printed after every few linear
              minimizations.
              if prin=4, the principal vectors of the approximating
              quadratic form are also printed.
     x        is an array containing on entry a guess of the point of
              minimum, on return the estimated point of minimum.
     f(x,n)   is the function to be minimized.  f should be a real*8
              function declared external in the calling program.
     fmin     is an estimate of the minimum, used only in printing
              intermediate results.
     The approximating quadratic form is
              q(x') = f(x,n) + (1/2) * (x'-x)-transpose * a * (x'-x)
     where x is the best estimate of the minimum and a is
              inverse(v-transpose) * d * inverse(v)
     (v(*,*) is the matrix of search directions; d(*) is the array
     of second differences).  If f has continuous second derivatives
     near x0, a will tend to the hessian of f at x0 as x approaches x0.

     It is assumed that on floating-point underflow the result is set
     to zero.
     The user should observe the comment on heuristic numbers after
     the initialization of machine dependent numbers.
*/

// machine dependent numbers:

namespace libcmatrix {

  using ::std::sqrt;

  struct praxis_status {
    size_t nf,nl;
    double fx;
    const double ldt;
    double dmin;
    double qd0,qd1,qf1;
    double qa,qb,qc;

    praxis_status(double fxv,double ldtv,double dminv) : nf(0), nl(0), fx(fxv), ldt(ldtv), dmin(dminv), qd0(0.0), qf1(fxv) {};
  };
  
//static const double machep=numeric_limits<double>::epsilon();
static const double machep=1e-16;
static const double small=machep*machep;
static const double vsmall=small*small;
static const double large=1.0/small;
static const double vlarge=1.0/vsmall;
static const double M2=sqrt(machep);
static const double M4=sqrt(M2);

/*   heuristic numbers:
     If the axes may be badly scaled (which is to be avoided if
     possible), then set scbd=10.  otherwise set scbd=1.
     If the problem is known to be ill-conditioned, set illc=true.
     otherwise set illc=false.
     ktm is the number of iterations without improvement before the
     algorithm terminates.  ktm=4 is very cautious; usually ktm=1
     is satisfactory. */

// praxis_print prints information about the progress of the iteration.

void praxis_print(const BaseList<double>& x,int verbose,double fmin,const praxis_status& status)
{
  ::std::cout << "Number of linear searches: " << status.nl << '\n';
  ::std::cout << "Number of function evaluations: " << status.nf << '\n';
  ::std::cout << "Smallest value found: " << status.fx << '\n';

  ::std::cout << "log(f(x)-" << fmin << ") ";
  if (status.fx<=fmin)
    ::std::cout << "is undefined\n";
  else
    ::std::cout << "= " << log10(status.fx-fmin) << '\n';

  if (x.length()<5 && verbose>3) ::std::cout << "X: " << x << '\n';
}

  /*  FLIN is the function of one variable to be minimized by MINNY.

    In fact, what is happening is that the scalar function F(X),
    where X is an N dimensional vector, is being minimized along a 
    fixed line.

    Modified: 15 March 2000 */

double flin(int j, double l,const BaseMinFunction& func, const BaseList<double>& x, const rmatrix& V, const BaseList<double>& q0, const BaseList<double>& q1, praxis_status& status)
{
  const size_t N=x.length();
  ScratchList<double> t(N);

  if (j>=0) {
    // The search is linear.
    t=x;
    mla(t,l,V.row(j));
  }
  else {
    // The search is along a parabolic space curve.
    const double& qd0= status.qd0;
    const double& qd1= status.qd1;

    const double qa = ( l*(l-qd1) ) / ( qd0*(qd0+qd1) );
    const double qb = ( (l+qd0) * (qd1-l) ) / (qd0*qd1);
    const double qc = ( l*(l+qd0) ) / ( qd1*(qd0+qd1) );

    for (int i=N;i--;)
      t(i) = ( qa*q0(i) + qb*x(i) ) + qc*q1(i);
  }

  // The function evaluation counter NF is incremented.
  status.nf++;
  return func(t);
}

  /* an improved version of minfit (see Golub and Reinsch, 1969)
     restricted to m=n,p=0.
     Rhe singular values of the array AB are returned in Q and AB is
     overwritten with the orthogonal matrix V such that U.diag(Q) = AB.V,
     where U is another orthogonal matrix. */

void minfit(double tol,rmatrix& AB,BaseList<double> Q)
{
  if (!issquare(AB)) throw NotSquare("minfit");
  const int N=AB.cols();
  if (Q.length()!=N) throw Mismatch("minfit");

  if (N==1) {
    Q(0)=AB(0,0);
    AB(0,0)=1.0;
    return;
  }

  ScratchList<double> E(N);

  //    Householder's reduction to bidiagonal form

  double eps=machep;
  double g=0.0;
  double x=0.0;
  double y=1e30;
  int i;

  for (i=0;i<N;i++) {
    const int l=i+1;

    E(i)=g;

    double s=0.0;
    for (int j=i;j<N;j++) s+=norm(AB(i,j));

    g=0.0;
    if (s>=tol) {
      const double ABii=AB(i,i);
      g= (ABii>=0.0) ? -sqrt(s) : sqrt(s);
      const double h=g*ABii-s;
      AB(i,i)=ABii-g;
      for (int j=l;j<N;j++) {
	double f=0.0;
	for (int k=i;k<N;k++) f+=AB(i,k)*AB(j,k);
	f/=h;
	for (int k=i;k<N;k++) AB(j,k)+=f*AB(i,k);
      }
    }
    Q(i)=g;
    if (i!=N-1) {
      s=0.0;
      for (int j=l;j<N;j++) s+=norm(AB(j,i));
      g=0.0;
      if (s>=tol) {
	const double f=AB(i+1,i);
	g= (f>=0.0) ? -sqrt(s) : sqrt(s);
	const double h=f*g-s;
	AB(i,i+1)=f-g;
	for (int j=l;j<N;j++) E(j)=AB(j,i)/h;
	for (int j=l;j<N;j++) {
	  s=0.0;
	  for (int k=l;k<N;k++) s+=AB(k,j)*AB(k,i);
	  for (int k=l;k<N;k++) AB(k,j)+=s*E(k);
	}
      }
    }
    y=fabs(Q(i))+fabs(E(i));
  }
  if (y>x) x=y;
  
  // Accumulation of right-hand transformations...
  AB(N-1,N-1)=1.0;
  g=E(N-1);
  int l=N-1;
  for (int ii=1;i<N;i++) {
    i=N-ii-1;
    if (g!=0.0) {
      const double h=AB(i+1,i)*g;
      for (int j=l;j<N;j++) AB(i,j)=AB(j,i)/h;
      for (int j=l;j<N;j++) {
	double s=0.0;
	for (int k=l;k<N;k++) s+=AB(k,i)*AB(j,k);
      }
    }
    for (int j=l;j<N;j++) AB(i,j)=AB(j,i)=0.0;
    AB(i,i)=1.0;
    g=E(i);
    l=i;
  }

  // Diagonalization of the bidiagonal form...
  eps*=x;
  for (int kk=0;kk<N;kk++) {
    const int k = N-kk-1;
    int kt=0;
    double z;
    for (;;) {
      kt++;
      if (kt>30) {
	E(k)=0.0;
	throw Failed("minfit: QR failed to converge");
      }

      double c=0.0;
      double s=1.0;
      double h;

      for (int ll2=0;ll2<k;ll2++) {
	l = k -ll2-1;
	if (fabs(E(l))<=eps) goto line120;
	if (l && (fabs(Q(l-1))<=eps)) break;
      }
      // Cancellation of E(l) if l>1...
      for (i=l;i<=k;i++) {
	const double f=s*E(i);
	E(i)*=c;
	if (fabs(f)<=eps) break;
	g=Q(i);
	// Q(I) = H = DSQRT(G*G + F*F)...
	if (fabs(f)<fabs(g))
	  h = fabs(g)*sqrt(1.0+norm(f/g));
	else
	  h = f ? fabs(f)*sqrt(1.0+norm(g/f)) : 0.0;
	
	Q(i)=h;
	if (h==0.0) g=h=1.0;
	c=g/h;
	s=-f/h;
      }
    line120:
      // Test for convergence...
      z=Q(k);
      if (l==k) break;
      
      // Shift from bottom 2*2 minor...
      x=Q(l);
      y=Q(k-1);
      g=E(k-1);
      h=E(k);
      double f=((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y);
      g=sqrt(f*f+1.0);
      const double temp = (f>=0.0) ? f+g : f-g;
      f=((x-z)*(x+z)+h*(y/temp-h))/x;
      
      // Next QR transformation.
      c = 1.0;
      s = 1.0;
      
      for (int i=l+1;i<k;i++) {
	g = E(i);
	y = Q(i);
	h = s * g;
	g *=c;
	
	if (fabs(f) < fabs(h))
	  z = fabs(h)*sqrt(1.0 + norm(f/h));
	else
	  z = (f==0.0) ? 0.0 : fabs(f)*sqrt(1.0 + norm(h/f));
	
	E(i-1) = z;
	if (z==0.0) f=z=1.0;
	
	c=f/z;
	s=h/z;
	f = x * c + g * s;
	g = -x * s + g * c;
	h = y * s;
	y *= c;
	
	for (int j=0;j<N;j++) {
	  x = AB(i-1,j);
	  z = AB(i,j);
	  AB(i-1,j) = x * c + z * s;
	  AB(i,j) = -x * s + z * c;
	}
	
	if ( fabs(f) < fabs(h) )
	  z = fabs(h) * sqrt( 1.0 + norm(f/h) );
	else 
	  z = (f==0.0) ? 0.0 : fabs(f) * sqrt( 1.0 + norm(h/f) );
	
	Q(i-1) = z;
	
	if (z==0.0) f=z=1.0;
	
	c = f / z;
	s = h / z;
	f = c * g + s * y;
	x = - s * g + c * y;
      }
      E(l)=0.0;
      E(k)=f;
      Q(k)=x;
    }

    if (z<0.0) {
      Q(k)=-z;
      minus_ip(AB.row(k));
    }
  }
}

/*  MINNY minimizes a scalar function of N variables along a line.

    MINNY minimizes F from x in the direction v(*,j) unless
    j is less than 1, when a quadratic search is made in the plane
    defined by q0,q1,x.

    d2 is either zero or an approximation to half f".
    on entry, x1 is an estimate of the distance from x to the minimum
    along v(*,j) (or, if j = 0, a curve).  on return, x1 is the distance
    found.
    if fk = true, then f1 is flin(x1), otherwise x1 and f1 are ignored
    on entry unless final fx is greater than f1.
    nits controls the number of times an attempt will be made to halve
    the interval. */

void minny ( int j, size_t nits, double d2, double x1, double f1, bool fk, const BaseMinFunction& func, const BaseList<double>& x, double t, double h, const rmatrix& V, const BaseList<double>& q0, const BaseList<double>& q1, praxis_status& status)
{
  double& fx = status.fx;
  size_t& nf = status.nf;
  const double& ldt=status.ldt;

  if (!issquare(V)) throw NotSquare("minny");
  const int N=V.rows();

  double sf1=f1;
  double sx1 = x1;
  int k = 0;
  double xm = 0.0;
  double fm = fx;
  double f0 = fx;
  bool dz = ( d2 < machep );

  // Find the step size.

  double s=sqrt(norm(x));
  double temp= dz ? status.dmin : d2;
  double t2 = M4*sqrt(fabs(fx)/temp + s*ldt) + M2*ldt;
  s = M4*s + t;
  if (dz && (t2>s)) t2=s;

  if (t2<small) t2=small;
  if (t2>0.01*h) t2=0.01*h;

  if (fk && (f1<=fm)) {
    xm=x1;
    fm=f1;
  }

  if (!fk || (fabs(x1)<t2)) {
    double x1= (x1>=0.0) ? t2 : -t2;
    f1 = flin ( j, x1, func, x, V, q0, q1, status );
  }

  if (f1<=fm) {
    xm=x1;
    fm=f1;
  }

  // Evaluate FLIN at another point and estimate the second derivative.
 line4:

  if (dz) {
    const double x2= (f0>=f1) ? 2.0*x1 : -x1;
    const double f2 = flin ( j, x2, func, x, V, q0, q1, status );

    if (f2<=fm) {
      xm = x2;
      fm = f2;
    }

    d2 = ( x2*(f1-f0) - x1*(f2-f0) ) / ( (x1*x2) * (x1-x2) );
  }

  // Estimate the first derivative at 0.

  const double d1 = (f1-f0) / x1- x1*d2;
  dz=true;

  // Predict the minimum.

  double x2;
  if (d2<=small)
    x2=(d1>=0.0) ? -h : h;
  else
    x2 = ( -0.5*d1 ) / d2;
  
  if (fabs(x2)>h)
    x2 = (x2<=0.0) ? -h : h;

  // Evaluate F at the predicted minimum.

  double f2;
  for (;;) {
    f2 = flin (j, x2, func, x, V, q0, q1, status );

    if ( k>=nits || (f2<=f0) ) break;

    // No success, so try again.
    k++;
    if ( (f0<f1) && (x1*x2>0.0) ) goto line4;
    x2*=0.5;
  }

  // Increment the one-dimensional search counter.
  status.nl++;

  if (f2>fm)
    x2=xm;
  else
    fm=f2;
  
  // Get a new estimate of the second derivative.

  if ( fabs( x2*(x2-x1) ) > small )
    d2 = ( x2*(f1-f0) - x1*(fm-f0) ) / ( (x1*x2) * (x1-x2) );
  else {
    if (k>0) d2=0.0;
  }

  if (d2<small) d2=small;

  x1 = x2;
  fx = fm;

  if ( sf1 < fx ) {
    fx = sf1;
    x1 = sx1;
  }

  // Update X for linear but not parabolic search.
  if (j>=0) mla(x,x1,V.row(j));
}

/* QUAD seeks to minimize the scalar function F along a particular curve.

    The minimizer to be sought is required to lie on a curve defined
    by Q0, Q1 and X.*/

void quad(const BaseMinFunction& func, BaseList<double> x, double t, double h, const rmatrix& V, BaseList<double> q0, BaseList<double> q1, praxis_status& status )
{
  double& qf1=status.qf1;
  double& qd0=status.qd0;
  double& qd1=status.qd1;
  double& fx=status.fx;
  double& qa=status.qa;
  double& qb=status.qb;
  double& qc=status.qc;
  
  SWAP(fx,qf1);

  qd1 = 0.0;
  const size_t N=x.length();
  for (int i=N;i--;) {
    const double s = x(i);
    const double l = q1(i);
    x(i) = l;
    q1(i) = s;
    qd1 += norm(s-l);
  }
  qd1 = sqrt(qd1);

  double l = qd1;
  double s = 0.0;

  if ( (qd0 <= 0.0) || (qd1 <= 0.0) || (status.nl < 3*N*N) ) {
    fx = qf1;
    qa = 0.0;
    qb = qa;
    qc = 1.0;
  }
  else {
    minny(-1, 2, s, l, qf1, true, func, x, t, h, V, q0, q1, status );
      
    qa = ( l*(l-qd1) ) / ( qd0*(qd0+qd1) );
    qb = ( (l+qd0) * (qd1-l) ) / (qd0*qd1);
    qc = ( l*(l+qd0) ) / ( qd1*(qd0+qd1) );
  }

  qd0 = qd1;

  for (int i=0;i<N;i++) {
    double s = q0(i);
    q0(i) = x(i);
    x(i) = ( qa*s + qb*x(i) ) + qc*q1(i);
  }
}

double praxis(BaseList<double> x,const BaseMinFunction& func,const BaseList<size_t>& actord,double h0,int verbose, double tol,int kterm,double scbd,bool illcond)
{
  if (kterm<1) throw InvalidParameter("Praxis: kterm must be >=1");

  const double ldfac= illcond ? 0.1 : 0.01;
  int kt=0;
  int nl=0;
  int nf=0;
  double t=small+fabs(tol);
  double t2=t;

  double h=100.0*t;
  if (h<h0) h=h0;

  praxis_status status(func(x),h,small);

  //The first set of search directions. V is the identity matrix.....

  const size_t N=x.length();
  rmatrix V;
  V.identity(N);

  ScratchList<double> D(N);
  D(0)=0.0;
  
  ScratchList<double> q0(x);
  ScratchList<double> q1(x);

  if (verbose) praxis_print(x,verbose,fmin);

  // main loop starts here.....
  for (;;) {
    double sf=D(0);
    D(0)=0.0;
    double s=0.0;
    // minimize along the first direction V(*,1).
    
    minny(1,2,D(0),S,fx,false,func,x,t,h,status);
    if (S<=0.0) {
      for (int i=N;i--;) V(0,i)=-V(1,i);
    }
    if (sf<=0.9*D(0) || 0.9*sf >=D(0)) {
      for (int i=1;i<N;i++) D(i)=0.0;
    }
    
    // inner loop starts here
    for (int k=1;k<N;k++) {
      
      y=x;
      sf=fx;
      if (kt>0) illcond=true;
      
    line80:
      kl=k;
      df=0.0;
      
      /* A random step follows (to avoid resolution valleys).
	 praxis assumes that random returns a random number uniformly
	 distributed in (0,1) */
      
      if (illcond) {
	for (int i=0;i<N;i++) {
	  
	  s=(0.1*ldt+t2*pow(10.0,kt))*(random(1.0)-0.5);
	  z(i)=s;
	  mla(x,s,V.row(i));
	}
	fx=func(x);
	nf++;
      }
      
      // Minimize along the "non-conjugate" directions V(*,K),...,V(*,N)
      
      for (int k2=k;k<N;k++) {
	sl=fx;
	s=0.0;
	minny(k2,2,D(k2),S,fx,false,func,x,t,hv,q0,q1,status);
	s = illcond ? D(k2)*square(s+z(k2)) : sl-fx; 
	
	if (df<=s) {
	  df=s;
	  kl=k2;
	}
      }
      
      if (!illcond && (df<fabs(100*machep*fx))) {
	/*.....If there was not much improvement on the first try, set
	  illc=true and start the inner loop again..... */	
	illcond=true;
	goto line80;
      }
      if ((k==2) && (verbose>1)) std::cout << "The second difference array: " << D << '\n';
      
      // Minimize along the "conjugate" directions V(*,1),...,V(*,K-1)
      
      for (int k2=0;k2<k;k2++)
	minny(k2,2,D(k2),0.0,fx,false,func,x,t,h,status);
      
      double f1=fx;
      double fx=sf;
      double lds=0.0;
      for (int i=0;i<N;i++) {
	double& xi=x(i);
	double& yi=y(i);
	double sl=xi-yi;
	xi=yi;
	yi=sl;
	lds+=sl*sl;
      }
      lds=sqrt(lds);
      
      if (lds>small) {
	
	/* Discard direction V(*,kl).
	   If no random step was taken, V(*,kl) is the "non-conjugate"
	   direction along which the greatest improvement was made..... */
	
	klmk=kl-k;
	for (int ii=0;i<klmk;i++) {
	  i=kl-ii;
	  V.row(i+1)=V.row(i);
	  D(i+1)=D(i);
	}
	
	D(k)=0;
	V.row(k)=y/lds;
	
	// minimize along the new "conjugate" direction v(*,k), which is
	// the normalized vector:  (new x) - (0ld x).....
	
	minny(k,4,D(k),lds,f1,true,func,x,t,hv,q0,q1,status);
	
	if (lds<=0.0) {
	  lds=-lds;
	  minus_ip(V.row(k));
	}
      }
      ldt*=ldfac;
      if (ldt<lds) ldt=lds;
      if (verbose) praxis_print(x,verbose,fmin);
      
      t2=M2*sqrt(norm(x))+t;
      
      // See whether the length of the step taken since starting the
      // inner loop exceeds half the tolerance.....
      
      if (ldt>0.5*t2) kt=-1;
      kt++;
      if (kt>kterm) { 
	if (verbose) std::cout << "X: " << x << '\n';
	return fx;
      }
    } // Inner loop ends
    
    // Try quadratic extrapolation in case we are in a curved valley.
    
    quad(func,x,t,h,v,q0,q1);
    double dn=0.0;
    for (int i=0;i<N;i++) {
      D(i)=1.0/sqrt(D(i));
      if (dn < D(i)) dn=D(i);
    }
    if (verbose>3) std::cout << "The new direction vectors\n" << V << '\n';
    
    for (int j=0;j<N;j++) V.row(j)*=D(j)/dn;
    
    // Scale the axes to try to reduce the condition number.....
    
    if (scalebd > 1.0) {
      S=vlarge;
      for (int i=0;i<N;i++) {
	sl=0.0;
	for (int j=0;j<N;j++) sl+=V(j,i)*V(j,i);
	z(i)=sqrt(sl);
	if (z(i)<M4) z(i)=M4;
	if (s>z(i)) s=z(i);
      }
      for (int i=0;i<N;i++) {
	sl=s/z(i);
	z(i)=1.0/sl;
	if (z(i)>scalebd) {
	  sl=1.0/scalebd;
	  z(i)=scalebd;
	}
	for (int j=0;j<N;j++) V(j,i)*=sl;
      }
    }
    
    // Calculate a new set of orthogonal directions before repeating the main loop.
    // first transpose v for minfit:
    
    V.transpose();
    
    /* Call minfit to find the singular value decomposition of v.
       This gives the principal values and principal directions of the
       approximating quadratic form without squaring the condition
       number */
    
    minfit(V,vsmall,D);
    
    // Unscale the axes
    
    if (scalebd>1.0) {
      for (int i=0;i<N;i++) {
	const double scale=z(i);
	for (int j=0;j<N;j++) V(j,i)*=scale;
      }
      for (int i=0;i<N;i++) {
	const double scale=sqrt(norm(V.row(i)));
	D(i)*=scale;
	V.row(i)*=1.0/scale;
      }
    }
    
    for (int i=0;i<N;i++) {
      const double dni=dn*D(i);
      if (dni>large)
	D(i)=vsmall;
      else
	D(i) = (dni<small) ? vlarge : 1.0/(dni*dni);
    }
    
    // Sort the eigenvalues and eigenvectors.....
    ScratchList<size_t> indx=sort_index(D);
    D=D(indx);
    V=V(slice(),indx);

    double dmin=D(N-1);
    if (dmin<small) dmin=small;
    status.dmin=dmin;
    illcond=(M2*D(0)>dmin);

    if (verbose>1) {
      if (scalebd>1.0) std::cout << "The scale factors: " << z << '\n';
      std::cout << "Principal values of the quadratic form: " << D << '\n';
      if (verbose>3) std::cout << "The principal axes\n" << V << '\n';
    }
  }
}

} //namespace libcmatrix
