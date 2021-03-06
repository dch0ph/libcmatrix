/* non-linear fitting adapted from Numerical Recipes */

#include <cstdlib>
#include "optim.h"
#include "List.h"
#include "rmatrix.h"

#define ASYMMETRIC_BYDEFAULT false
const double difference_ratio=1e-2;

//DOTHREAD stuff is broken
//#if defined(LCM_REENTRANT) && defined(HAVE_LIBPTHREAD)
//#define DOTHREAD 1
//#include "cmatrix_threads.h"
//#include "ttyio.h"
//#endif

namespace libcmatrix {

 //  class LevMarqCalc {
//   public:
//     LevMarqCalc(LevMarqFit& fitfuncv)
//       : fitfunc_(fitfuncv) {}

//     void operator()(size_t start, size_t end, size_t nthr)
//     {
//       if (start+1!=end) throw InternalError("calc_func");
//       fitfunc_.run(nthr);
//     }
//   private:
//     LevMarqFit& fitfunc_;
//   };


static void calc_grad_nothread(rmatrix &grad,const BaseFitFunction& fobj,const BaseList<double> &p,const BaseList<size_t> &actord, const BaseList<double> &errs,BaseList<double> ptry,BaseList<double> ytrial,int verbose =0,bool asym =false)
{
  if (p.length()!=errs.length())
    throw Mismatch("calc_grad: parameters/errors");

  const size_t actpars=actord.length();

  grad.create(actpars,ytrial.length());

  if (asym)
    fobj(ytrial,p);

  for (size_t j=actpars;j--;) {
    BaseList<double> gradr=grad.row(j);

    const size_t which=actord(j);
    ptry=p;
    const double dv=difference_ratio*errs(which);
    if (asym) { // only calculate finite difference to one side; quicker but (even) less rigorous
      ptry(which)=p(which)+dv;
      fobj(gradr,ptry);
      gradr-=ytrial;
      gradr/=dv;
    }
    else {
      ptry(which)=p(which)-dv;
      fobj(ytrial,ptry);
      ptry(which)=p(which)+dv;
      fobj(gradr,ptry);
      gradr-=ytrial;
      gradr/=2*dv;
    }
    if (verbose>2)
      std::cout << "Gradient " << which << ": " << gradr << '\n';
  }
}

void LevMarqFit::calc_grad(const BaseList<double>& p_)
{
// #ifdef DOTHREAD
//   if (tcon.get_max_threads()) {

//     const size_t actpars=variable_parameters();
//     const size_t npars=total_parameters();

//     if (p.length()!=npars) throw Mismatch("calc_grad");
//     grad.create(actpars,datapoints());

//     if (asym) {
//       paras.create(actpars+1,npars);
//       stores.create(actpars+1);
//       for (size_t k=0;k<=actpars;k++) {
// 	BaseList<double> cpars=paras.row(k);
// 	cpars=p_;
// 	if (k==actpars)
// 	  stores(k).create(ytrial);
// 	else {
// 	  stores(k).create(grad.row(k));
// 	  const size_t which=actord(k);
// 	  const double dv=1e-4*errs(which);
// 	  cpars(which)+=dv;
// 	}
//       }
      
//       LevMarqCalc threadfunc(*this);
//       tcon.run(threadfunc,actpars+1,1);

//       for (size_t j=actpars;j--;) {
// 	BaseList<double> gradr=grad.row(j);
// 	const double dv=1e-4*errs(actord(j));
// 	gradr-=ytrial;
// 	gradr/=dv;
//       }
//     }
//     else {
//       paras.create(2*actpars,npars);
//       stores.create(2*actpars);

//       size_t k=0;
//       size_t k2=0;
//       for (;k<actpars;k++,k2+=2) {
// 	const size_t which=actord(k);
// 	const double dv=1e-4*errs(which);
// 	BaseList<double> cparsp=paras.row(k2);
// 	cparsp=p_; cparsp(which)-=dv;
// 	stores(k2).create(ytmp.row(k));
// 	BaseList<double> cparsm=paras.row(k2+1);
// 	cparsm=p_; cparsm(which)+=dv;
// 	stores(k2+1).create(grad.row(k));
//       }

//       LevMarqCalc threadfunc(*this);
//       tcon.run(threadfunc,2*actpars,1);

//       for (size_t j=actpars;j--;) {
// 	BaseList<double> gradr=grad.row(j);	
// 	const double dv=1e-4*errs(actord(j));
// 	gradr-=ytmp.row(j);
// 	gradr/=2*dv;
//       }
//     }
//   }
//   else
// #endif
    calc_grad_nothread(grad,fobj,p_,actord,errs,pscr,ytrial,verbose_,asym);
}

double LevMarqFit::mqrcof(rmatrix& alpha_,BaseList<double> beta_,const BaseList<double> &p_)
{
  size_t j,k,n;

  const size_t actpars=variable_parameters();
  const size_t npts=datapoints();
  
  calc_grad(p_);
  if (!asym)
    fobj(ytrial,p_);
  havevalidspec=true;
  subtract(dy,data,ytrial);

  alpha_.create(actpars,actpars);

  for (j=0;j<actpars;j++) {
    const BaseList<double> gradj=grad.row(j);
    for (k=j;k<actpars;k++) {
      const BaseList<double> gradk=grad.row(k);
      double sum=0.0;
      double sbeta=0.0;
      for (n=npts;n--;) {
	const double sgrad=gradj(n)/w2(n);
	sum+=gradk(n)*sgrad;
	sbeta+=dy(n)*sgrad;
      }
      alpha_(j,k)=sum;
      beta_(j)=sbeta;
    }
  }
  for (j=1;j<actpars;j++)
    for (k=0;k<j;k++)
      alpha_(j,k)=alpha_(k,j);

  double chisqr=0.0;
  for (n=npts;n--;)
    chisqr+=dy(n)*dy(n)/w2(n);
  
  return chisqr/npts; //!< adjusted to return normalised chi^2  Nov 06
}

void hessian(rmatrix &a,const BaseFitFunction& fobj, const BaseList<double> &p,const BaseList<size_t> &actord, const BaseList<double> &errs,const BaseList<double> &weights, int verbose)
{
  size_t j,k,n;

  const size_t actpars=actord.length();
  const size_t npts=weights.length();

  ScratchList<double> ptry(p.length());
  ScratchList<double> ytrial(npts);

  rmatrix grad;
  calc_grad_nothread(grad,fobj,p,actord,errs,ptry,ytrial,verbose,ASYMMETRIC_BYDEFAULT);

  BaseList<double>& w2(ytrial);
  w2=weights; w2*=weights;

  a.create(actpars,actpars);

  for (j=0;j<actpars;j++) {
    const BaseList<double> gradj=grad.row(j);
    for (k=j;k<actpars;k++) {
      const BaseList<double> gradk=grad.row(k);
      double sum=0.0;
      for (n=npts;n--;)
	sum+=gradk(n)*gradj(n)/w2(n);
      a(j,k)=sum;
    }
  }
  for (j=1;j<actpars;j++)
    for (k=0;k<j;k++)
      a(j,k)=a(k,j);
}

void covariance(rmatrix &a,const BaseFitFunction& fobj, const BaseList<double> &p,const BaseList<size_t> &actord, const BaseList<double> &errs,const BaseList<double> &weights, int)
{
  rmatrix tmp;
  hessian(tmp,fobj,p,actord,errs,weights);
  try {
    a=inv(tmp);
  } catch (MatrixException& exc) {
    std::cerr << "covariance: failed to invert Hessian matrix - function independent of one or more parameters?\n";
    bool donesuspect=false;
    const int indexbase=cmatrix_ostream_controller(std::cerr).indexbase;
    for (size_t i=0;i<tmp.rows();i++) {
      if (tmp(i,i)==0.0) {
	if (!donesuspect) {
	  std::cerr << "Suspect *variable* parameter(s):";
	  donesuspect=true;
	}
	std::cerr << ' ' << (i+indexbase);
      }
    }
    if (donesuspect)
      std::cerr << std::endl;	
    throw;
  }
}

double LevMarqFit::initialise()
{
  const size_t npars=p.length();
  const size_t actpars=actord.length();
  const size_t npts=data.length();

  if (actpars>npars || actpars<1)
    throw InvalidParameter("fitdata");
  if (npts<actpars)
    throw Failed("fitdata");
  if (npars!=errs.length())
    throw Mismatch("fitdata: number of parameters doesn't match number of error steps");
  for (size_t i=w2.size();i--;) {
    if (w2(i)==0.0)
      throw InvalidParameter("fitdata: noise vector contains zero or negative elements");
  }

  asym=false;
  lambda=0.001;
  verbose_=1;
  ptry.create(npars);
  pscr.create(npars);
  datastore.create(2*npts);
  ytrial.create(npts,datastore.vector());
  dy.create(npts,datastore.vector()+npts);
  ytmp.create(actpars,npts);
  beta.create(actpars);
  havevalidspec=false;

  for (size_t i=0;i<actpars;i++) {
    const double step=errs(actord(i));
    if (!step)
      throw InvalidParameter("LevMarqFit: error must not be zero");
    if (verbose_>1)
      std::cerr << "Delta for parameter " << i << ": " << (asym ? "+" : "+/-") << (step*difference_ratio) << '\n';    
  }
// #ifdef DOTHREAD
//   const char *envp=getenv("OMP_NUM_THREADS");
//   if (envp) threads(parseint(envp));
// #endif

  return (ochisqr=mqrcof(alpha,beta,p));
}

void LevMarqFit::threads(int n)
{
#ifdef DOTHREAD
  tcon.create(n);
#else
  ::std::cerr << "WARNING: calc to LevMarqFit::threads ignored" << ::std::endl;
#endif
}

LevMarqFit::LevMarqFit(double &rochisqr,const BaseList<double> &p_,const BaseFitFunction& fobj_,const BaseList<double> &data_,const BaseList<size_t> &actord_, const BaseList<double> &errs_, const BaseList<double> &weights_) : fobj(fobj_), data(data_), actord(actord_), errs(errs_), p(p_)
{
  if (data_.length()!=weights_.length())
    throw Mismatch("fitdata");
  w2=weights_; w2*=weights_;

  rochisqr=initialise();
}

LevMarqFit::LevMarqFit(double &rochisqr,const BaseList<double> &p_,const BaseFitFunction& fobj_,const BaseList<double> &data_,const BaseList<size_t> &actord_, const BaseList<double> &errs_, double weight) : fobj(fobj_), data(data_), actord(actord_), errs(errs_), p(p_)
{
  w2.create(data.length());
  w2=weight*weight;

  rochisqr=initialise();
}


double LevMarqFit::next()
{
  move();
  if (verbose_) {
    std::cout << last_params() << " \t";
    if (verbose_>1)
      std::cout << current_lambda() << " \t";
    std::flush(std::cout);
  }
  return evaluate();
}

void LevMarqFit::move()
{
  size_t j;

  const size_t actpars=variable_parameters();

  tmp=alpha;
  for (j=actpars;j--;)
    tmp(j,j)=alpha(j,j)*(1.0+lambda);
  //std::cout << "Alpha\n" << tmp << std::endl;
  da=beta;
  //std::cout << "Beta: " << beta << std::endl;

  //  try {
    multiply_inv_ip2(tmp,da);
//   }
//   catch (SingularMatrix& exc) {
//     //matrix was singular
//     std::cerr << "Inversion of Hessian failed: contents in HessianDump\n";
//     try {
//       spy(std::cerr,tmp);
//       std::cerr << tmp << std::endl;
//       write_matrix("HessianDump",tmp); //trap errors
//     } catch (...) {}
//     throw;
//   }
  
  //std::cout << "change: " << da << std::endl;
  
  ptry=p;
  for (j=actpars;j--;) {
    const size_t which=actord(j);
    ptry(which)=p(which)+da(j);
  }
}

double LevMarqFit::evaluate()
{
  const double chisqr=mqrcof(tmp,da,ptry);
  
  if (chisqr<ochisqr) {
    lambda*=0.1;
    alpha=tmp;
    beta=da;
    p=ptry;
    ochisqr=chisqr;
    if (asym)
      havevalidspec=true;
  }
  else {
    lambda*=10.0;
    havevalidspec=false; //invalidate ytrial
  }
  return chisqr;
}

BaseList<double> LevMarqFit::current_trial()
{
  if (!havevalidspec) { 
    fobj(ytrial,p);
    havevalidspec=true;
  } 
  return ytrial;
}

double LevMarqFit::final(rmatrix& covar,BaseList<double> p_)
{
  p_=p;
  covar=inv(alpha);
  return ochisqr;
}

BaseWarning optimisation_base_warning(BaseWarning::Silent); //!< by default warnings are silent

Warning<> fitsmallimprovement_warning("Small (<5%) improvement in chi^2 observed. Fitting against relatively large non-fitted background without increasing convergence tolerance or masking out non-fitted regions will introduce inaccuracies.",&optimisation_base_warning);

static optimisation_info dofitdata(rmatrix& covar,BaseList<double> p,LevMarqFit& fitobj,double ochisqr,int verbose,double chistop,int maxsteps)
{
  optimisation_info info;
  info.converged=false;
  fitobj.verbose(verbose);
  fitobj.asymmetric(ASYMMETRIC_BYDEFAULT);

//   if (verbose>1) {
//     std::cout << "Asymmetric gradient calculation: " << (ASYMMETRIC_BYDEFAULT ? "yes\n" : "no\n");
//     //    if (fitobj.threads()) std::cout << "Threads used: " << fitobj.threads() << "\n";
//   }

  if (maxsteps<0 || (chistop<=0) || (chistop>=1.0))
    throw InvalidParameter("fitdata");

  size_t steps=0;
  const double firstchisqr=ochisqr;

  if (verbose) {
    std::cout << "Initial (normalised) chi^2: " << ochisqr << " (calculated on " << fitobj.datapoints() << " data points)" << std::endl;
    if (maxsteps>1) {
      if (verbose>1)
	std::cout << "Stop parameter (fractional change in chi^2): " << chistop << std::endl;

      std::cout << "Step  Parameters... ";
      if (verbose>1)
	std::cout << "lambda ";
      std::cout << "chi^2";
      if (verbose>1)
	std::cout << "stop ";
      std::cout << '\n';    
    }
  }

  while (steps<maxsteps) {
    steps++;
    info.iterations=steps;
    if (verbose)
      ::std::cout << steps << ": ";
    const double chisqr=fitobj.next();
    const double stopcrit=(ochisqr-chisqr)/ochisqr;
    if (verbose) {
      std::cout << chisqr;
      if (verbose>1)
	std::cout << " (" << stopcrit << ')';
      std::cout << std::endl;
    }
    
    if (chisqr>ochisqr)
      continue; //chi^2 increasing! - carry on...
    
    if (stopcrit<chistop) {
      info.optimum=fitobj.final(covar,p);
      info.converged=true;
      const double fracdrop=(firstchisqr-chisqr)/firstchisqr;
      if ((fracdrop<0.05) && (chisqr>2))
	fitsmallimprovement_warning.raise();

      return info;
    }

    ochisqr=chisqr;
  }
  info.optimum=ochisqr;
  p=fitobj.current_params();
  info.converged=false;
  return info;
}

optimisation_info fitdata(rmatrix &covar,BaseList<double> p, const BaseFitFunction& fobj, const BaseList<double> &data,const BaseList<size_t> &actord, const BaseList<double> &errs, const BaseList<double> &weights,int verbose,double chistop,int maxsteps)
{
  double ochisqr;
  LevMarqFit fitobj(ochisqr,p,fobj,data,actord,errs,weights);
  return dofitdata(covar,p,fitobj,ochisqr,verbose,chistop,maxsteps);
}

optimisation_info fitdata(rmatrix &covar,BaseList<double> p, const BaseFitFunction& fobj, const BaseList<double> &data,const BaseList<size_t> &actord, const BaseList<double> &errs, double weight,int verbose,double chistop,int maxsteps)
{
  double ochisqr;
  LevMarqFit fitobj(ochisqr,p,fobj,data,actord,errs,weight);
  return dofitdata(covar,p,fitobj,ochisqr,verbose,chistop,maxsteps);
}

}//namespace libcmatrix
