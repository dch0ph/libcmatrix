/* Multi-dimensional minimisation using Simplex method
   (from Numerical Recipes p.293) */

#include <cmath>
#include "List.h"
#include "ScratchList.h"
#include "optim.h"

namespace libcmatrix {

static const float alpha=1.0;
static const float beta=0.5;
static const float cgamma=2.0;

void init_amoeba(List<double> &y, const BaseMinFunction& fobj,const rmatrix &p)
{
  y.create(p.rows());
  for (int i=y.length();i--;)
    y(i)=fobj(p.row(i));
}

void init_amoeba(List<double> &y, rmatrix &p, const BaseMinFunction& fobj, const BaseList<double> &p0,const BaseList<size_t> &actord,const BaseList<double> &errs)
{
  const int npars=p0.length();
  const int actpars=actord.length();
  if (actpars>npars || actpars<1)
    throw InvalidParameter("init_amoeba");

  if (errs.length()!=npars)
    throw Mismatch("init_amoeba");
  p.create(actpars+1,npars);
  
  p.row(0)=p0;
  for (int i=0;i<actpars;i++) {
    const int which=actord(i);
    p.row(i+1)=p0;
    p(i+1,which)+=errs(which);
  }
  init_amoeba(y,fobj,p);
}

optimisation_info simplex_fitdata(rmatrix &covar,BaseList<double> p, const BaseFitFunction& fobj,const BaseList<double> &data,const BaseList<size_t> &actord,const BaseList<double> &errs,double weight,int verbose,double chistop,int maxsteps)
{
  List<double> weights(data.length(),weight);
  return simplex_fitdata(covar,p,fobj,data,actord,errs,weights,verbose,chistop,maxsteps);
}

// struct chi2block {
//   const BaseFitFunction& fobj;
//   const BaseList<double> &weights;
//   List<double> tdata;
//   const BaseList<double> &fdata;

//   chi2block(const BaseFitFunction&,const BaseList<double> &data,const BaseList<double> &_weights);
// };

chi2func::chi2func(const BaseFitFunction& fobj_,const BaseList<double> &_data,const BaseList<double> &_weights) : fobj(fobj_), weights(_weights), fdata(_data)
{
  const int dpts=_weights.length();
  if (dpts!=fdata.length())
    throw Mismatch("chi2func");
  tdata.create(dpts);
}

double chi2func::operator()(const BaseList<double> &pars) const
{
  fobj(tdata,pars);

  double chi2=0.0;
  for (int i=tdata.length();i--;) {
    const double diff=(tdata(i)-fdata(i))/weights(i);
    chi2+=diff*diff;
  }
  return chi2;
}
  
//! returns iterations used which will be <max_iter if converged
int amoeba(const BaseMinFunction& fobj,rmatrix& p, BaseList<double> y, const BaseList<size_t>& actord, double tol,int verbose,int max_iter)
{
  if (max_iter<1)
    throw InvalidParameter("amoeba: iterations must be >0");

  const int actpars=actord.length();
  const int verts=p.rows();

  const int totpars=p.cols();
  if (verts!=actpars+1 || verts!=y.length())
    throw Mismatch("amoeba");

  int i,j;

  if ((tol<=0.0) || (tol>=1.0))
    throw InvalidParameter("amoeba: stop criterion");

  if (verbose>1)
    std::cout << "Stopping criterion: " << tol << std::endl;

  ScratchList<double> pall(3*totpars);
  BaseList<double> pbar(totpars,pall.vector());
  pbar=p.row(0);
  BaseList<double> pr(totpars,pall.vector()+totpars);
  pr=p.row(0);
  BaseList<double> prr(totpars,pall.vector()+2*totpars);
  prr=p.row(0);

  for (int iter=0;;) {
    int ilo=1;
    int ihi,inhi;
    
    if (y(0)>y(1)) {
      ihi=1;
      inhi=2;
    }
    else {
      ihi=2;
      inhi=1;
    }
    for (i=0;i<verts;i++) {
      if (y(i)<y(ilo))
	ilo=i;
      if (y(i)>y(ihi)) {
	inhi=ihi;
	ihi=i;
      }
      else {
	if (y(i)>y(inhi)) {
	  if (i!=ihi)
	    inhi=i;
	}
      }
    }
    
    const double rtol=2*fabs(y(ihi)-y(ilo))/(fabs(y(ihi))+fabs(y(ilo)));
    if (rtol<tol)
      return iter;
    iter++;
    if (iter>=max_iter)
      return max_iter;
    
    pbar=0.0;
    
    for (i=verts;i--;) {
      if (i!=ihi)
	pbar+=p.row(i);
    }
    pbar/=actpars;

    for (j=actpars;j--;) {
      const int which=actord(j);
      pr(which)=(1+alpha)*pbar(which)-alpha*p(ihi,which);
    }

    if (verbose)
      std::cout << iter << ": " << pr << " \t" << std::flush;

    const double ypr=fobj(pr);

    if (verbose)
      std::cout << ypr << std::endl;

    if (ypr<=y(ilo)) {
      for (j=actpars;j--;) {
	const int which=actord(j);
	prr(which)=cgamma*pr(which)+(1.0-cgamma)*pbar(which);
      }
      const double yprr=fobj(prr);
      if (yprr<y(ilo)) {
	p.row(ihi)=prr;
	y(ihi)=yprr;
      }
      else {
	p.row(ihi)=pr;
	y(ihi)=ypr;
      }
    }
    else {
      if (ypr>=y(inhi)) {
	if (ypr<y(ihi)) {
	  p.row(ihi)=pr;
	  y(ihi)=ypr;
	}
	for (j=actpars;j--;) {
	  const int which=actord(j);
	  prr(which)=beta*p(ihi,which)+(1.0-beta)*pbar(which);
	}
	const double yprr=fobj(prr);
	if (yprr<y(ihi)) {
	  p.row(ihi)=prr;
	  y(ihi)=yprr;
	}
	else {
	  for (i=verts;i--;) {
	    if (i!=ilo) {
	      for (j=actpars;j--;) {
		const int which=actord(j);
		pr(which)=0.5*(p(i,which)+p(ilo,which));
	      }
	      p.row(i)=pr;
	      y(i)=fobj(pr);
	    }
	  }
	}
      }
      else {
	p.row(ihi)=pr;
	y(ihi)=ypr;
      }
    }
  }
}

template<class T> int min_index(const BaseList<T>& y)
{
  int i=y.length()-1;
  if (i<0)
    throw Failed("min_index");

  int besti=i--;

  for (;i--;) {
    if (y(i)<y(besti))
      besti=i;
  }
  return besti;
}

optimisation_info simplex(BaseList<double> p0,const BaseMinFunction& fobj,const BaseList<size_t> &actord,const BaseList<double> &errs,int verbose,double stop_val,int maxsteps)
{
  List<double> y;
  rmatrix p;
  
  init_amoeba(y,p,fobj,p0,actord,errs);
  optimisation_info info;
  info.iterations=amoeba(fobj,p,y,actord,stop_val,verbose,maxsteps);
  info.converged=(info.iterations<maxsteps);  
  const int besti=min_index(y);
  p0=p.row(besti);
  info.optimum=y(besti);
  return info;
}

Warning<> optim_covariance_warning("calculation of covariance matrix failed (badly defined problem?)",&lcm_base_warning);

optimisation_info simplex_fitdata(rmatrix& covar, BaseList<double> p0, const BaseFitFunction& fobj, const BaseList<double>& data, const BaseList<size_t>& actord, const BaseList<double>& errs, const BaseList<double>& weights, int verbose, double chistop, int maxsteps)
{
  chi2func chi2(fobj,data,weights);

  const optimisation_info info(simplex(p0,chi2,actord,errs,verbose,chistop,maxsteps));

  //if hessian c
  try {
    hessian(covar,fobj,p0,actord,errs,weights);
    covar=inv(covar);
  } catch (...) {
    optim_covariance_warning.raise();
    covar.clear();
  }

  return info;
}

}//namespace libcmatrix
