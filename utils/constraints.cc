#include "optim.h"

namespace libcmatrix {

namespace {
  void to_external(BaseList<double> pvals, const BaseList<double>& p, BaseList<Parameter>& pars, int verbose)
  {
    if (p.size()!=pvals.size())
      throw Mismatch("to_external");
    for (size_t i=p.size();i--;) {
      Parameter& par(pars(i));
      if (!par.isfixed())
	par.set_internal(p(i));    
      pvals(i)=par.get();
    }
    if (verbose>1)
      std::cout << "Mapped internal parameter values: " << p << " to external values: " << pvals << '\n';
  }
}

  Parameter::Parameter(const char* namev, double valv, double errv)
    : name_(namev), isfixed_(false)
  { setraw(valv); error(errv); }

  Parameter::Parameter(double valv, double errv)
    : isfixed_(false)
  { setraw(valv); error(errv); }

  void Parameter::error(double errv)
  {
    if (errv<=0.0)
      throw InvalidParameter("Parameter::error must be >0");
    error_=errv;
  }

  void Parameter::constrain(const BaseBoundFunction& func) 
  { 
    if (!func.isvalid(externalvalue_))
      throw Failed("Parameter::constrain: current parameter setting outside new constraint");
    constrainfuncp_.reset(func.clone());
    internalvalue_=func.external_to_internal(externalvalue_);
  } //!< apply constraint function (copied)

  void Parameter::unconstrain()
  {
    constrainfuncp_.clear();
    internalvalue_=externalvalue_;
  }

  void Parameter::setraw(double v) 
  {
    if (!isvalid(v))
      throw InvalidParameter("Parameter::set: attempt to set parameter outside constraint");      
    externalvalue_=v;
    internalvalue_=external_to_internal(v);
  }
  
  void Parameter::set(double v) 
  {
    if (externalvalue_!=v)
      setraw(v);
  }

  void Parameter::set_internal(double v) 
  {
    if (internalvalue_!=v) {
      internalvalue_=v;
      externalvalue_=internal_to_external(v);
    }
  }

  class ConstrainedFitFunction : public BaseFitFunction
  {
  public:
    ConstrainedFitFunction(const BaseFitFunction&, BaseList<Parameter>, int =0);
    void operator()(BaseList<double>, const BaseList<double>&) const;
  private:
    mutable ScratchList<double> pvals;
    const BaseFitFunction& rawfunc;
    mutable BaseList<Parameter> pars;
    int verbose;
  };
  
  ConstrainedFitFunction::ConstrainedFitFunction(const BaseFitFunction& rawfuncv, BaseList<Parameter> parsv, int verbosev)
    : pvals(parsv.size()),
      rawfunc(rawfuncv), pars(parsv),
      verbose(verbosev) {}

  void ConstrainedFitFunction::operator()(BaseList<double> dest, const BaseList<double>& p) const
  {
    to_external(pvals,p,pars,verbose);
    rawfunc(dest,pvals);
  }
  
  class ConstrainedMinFunction : public BaseMinFunction
  {
  public:
    ConstrainedMinFunction(const BaseMinFunction&, BaseList<Parameter>, int =0);
    double operator()(const BaseList<double>&) const;
  private:
    mutable ScratchList<double> pvals;
    const BaseMinFunction& rawfunc;
    mutable BaseList<Parameter> pars;
    int verbose;
  };
  
  ConstrainedMinFunction::ConstrainedMinFunction(const BaseMinFunction& rawfuncv, BaseList<Parameter> parsv, int verbosev)
  : pvals(parsv.size()),
    rawfunc(rawfuncv), pars(parsv),
    verbose(verbosev) {}

double ConstrainedMinFunction::operator()(const BaseList<double>& p) const
{
  to_external(pvals,p,pars,verbose);
  return rawfunc(pvals);
}

std::ostream& operator<< (std::ostream& ostr, const SimpleBoundsState& a)
{
  bool haveconstraint=false;
  if (a.islowerset()) {
    ostr << "Minimum: " << a.lower();
    haveconstraint=true;
  }
  if (a.isupperset()) {
    if (haveconstraint)
      ostr << "  ";
    ostr << "Maximum: " << a.upper();
    haveconstraint=true;
  }
  if (!haveconstraint)
    ostr << "Unconstrained";
  return ostr;
}

double SimpleBoundsState::lower() const 
{
  if (!islowerset())
    throw Failed("SimpleBoundsState: lower bound unset");
  return min_;
}

double SimpleBoundsState::upper() const 
{
  if (!isupperset())
    throw Failed("SimpleBoundsState: upper bound unset");
  return max_;
}

bool SimpleBoundsState::lower(double val)
{
  bool needupdate=false;
  if (havestate_!=MIN) {
    havestate_=MIN;
    needupdate=true;
    min_=val;
  }
  else 
    update(min_,val,needupdate);

  if (needupdate)
    doupdate();
  return needupdate;
}

  bool SimpleBoundsState::unconstrain()
  {
    if (havestate_!=NONE) {
      havestate_=NONE;
      return true;
    }
    return false;
  }

bool SimpleBoundsState::upper(double val)
{
  bool needupdate=false;
  if (havestate_!=MAX) {
    havestate_=MAX;
    max_=val;
    needupdate=true;
  }
  else
    update(max_,val,needupdate);

  if (needupdate)
    doupdate();
  return needupdate;
}

bool SimpleBoundsState::lowerupper(double minval, double maxval)
{
  if (maxval<=minval)
    throw InvalidParameter("SimpleBoundsState: maximum <= minimum!");
  bool needupdate=false;
  if (havestate_!=BOTH) {
    max_=maxval;
    min_=minval;
    havestate_=BOTH;
    needupdate=true;
  }
  else {
    update(max_,maxval,needupdate);
    update(min_,minval,needupdate);
  }
  if (needupdate)
    doupdate();
  return needupdate;
}

void SimpleBoundsState::doupdate()
{
  switch (havestate_) {
  case NONE:
    funcp.clear();
    break;
  case MIN:
    funcp.reset(new LowerBoundFunction(min_));
    break;
  case MAX:
    funcp.reset(new UpperBoundFunction(max_));
    break;
  case BOTH:
    funcp.reset(new LowerUpperBoundFunction(min_,max_));
    break;
  }
}

struct optim_interface {
  List<size_t> actord;
  ScratchList<double> p;
  ScratchList<double> errs;
  int verbose;
  bool ignoreconstraints;
  optim_interface(const BaseList<Parameter>&, int, bool =false);
  void copyback(BaseList<Parameter>) const; //!< copy parameter values back
};

optim_interface::optim_interface(const BaseList<Parameter>& pars, int verbosev, bool ignorev)
  : actord(pars.size()), p(pars.size()), errs(pars.size()), 
    verbose(verbosev),
    ignoreconstraints(ignorev) 
{
  actord.create(0U);
  for (size_t i=0;i<pars.size();i++) {
    const Parameter& par(pars(i));
    if (!par.isfixed())
      actord.push_back(i);
    p(i)=ignoreconstraints ? par.get() : par.get_internal();
    errs(i)=par.error();
  }
}

void optim_interface::copyback(BaseList<Parameter> dest) const
{
  if (dest.size()!=p.size())
    throw Mismatch("optim_interface::copyback");
  if (ignoreconstraints) {
    for (size_t i=p.size();i--;)
      dest(i).set(p(i));
  }
  else {
    for (size_t i=p.size();i--;)
      dest(i).set_internal(p(i));
  }
}

void hessian(rmatrix & d,const BaseFitFunction& func, const BaseList<Parameter>& p, const BaseList<double>& weights, int verbose, bool ignoreconstraints)
{
  optim_interface pinfo(p,verbose,ignoreconstraints);
  hessian(d,func,pinfo.p,pinfo.actord,pinfo.errs,weights,verbose);
}

void covariance(rmatrix & d,const BaseFitFunction& func, const BaseList<Parameter> &p,const BaseList<double> &weights, int verbose, bool ignoreconstraints)
{
  optim_interface pinfo(p,verbose,ignoreconstraints);
  if (ignoreconstraints)
    covariance(d,func,pinfo.p,pinfo.actord,pinfo.errs,weights,verbose);
  else {
    ConstrainedFitFunction tfunc(func,p,verbose);
    covariance(d,tfunc,pinfo.p,pinfo.actord,pinfo.errs,weights,verbose);
  }
}

optimisation_info fitdata(rmatrix &covar, BaseList<Parameter> p, const BaseFitFunction& func, const BaseList<double> &data,const BaseList<double> &weights,int verbose,double chistop,int maxsteps)
{
  optim_interface pinfo(p,verbose);
  ConstrainedFitFunction tfunc(func,p,verbose);
  const optimisation_info info(fitdata(covar,pinfo.p,tfunc,data,pinfo.actord,pinfo.errs,weights,verbose,chistop,maxsteps));
  pinfo.copyback(p);
  return info;
}

optimisation_info fitdata(rmatrix &covar, BaseList<Parameter> p, const BaseFitFunction& func, const BaseList<double> &data, double weight,int verbose,double chistop,int maxsteps)
{
  optim_interface pinfo(p,verbose);
  ConstrainedFitFunction tfunc(func,p,verbose);
  const optimisation_info info(fitdata(covar,pinfo.p,tfunc,data,pinfo.actord,pinfo.errs,weight,verbose,chistop,maxsteps));
  pinfo.copyback(p);
  return info;
}

optimisation_info simplex(BaseList<Parameter> p,const BaseMinFunction& func, int verbose,double stop_val,int maxstep)
{
  optim_interface pinfo(p,verbose);
  ConstrainedMinFunction tfunc(func,p,verbose);
  const optimisation_info info(simplex(pinfo.p,tfunc,pinfo.actord,pinfo.errs,verbose,stop_val,maxstep));
  pinfo.copyback(p);
  return info;
}

  optimisation_info simplex_fitdata(rmatrix& covar, BaseList<Parameter> p,const BaseFitFunction& func, const BaseList<double>& data, const BaseList<double>& weights, int verbose, double chistop, int maxsteps)
{
  optim_interface pinfo(p,verbose);
  ConstrainedFitFunction tfunc(func,p,verbose);
  const optimisation_info info(simplex_fitdata(covar,pinfo.p,tfunc,data,pinfo.actord,pinfo.errs,weights,verbose,chistop,maxsteps));
  pinfo.copyback(p);
  return info;
}

optimisation_info simplex_fitdata(rmatrix& covar, BaseList<Parameter> p,const BaseFitFunction& func, const BaseList<double>& data, double weight,int verbose, double chistop, int maxsteps)
{
  optim_interface pinfo(p,verbose);
  ConstrainedFitFunction tfunc(func,p,verbose);
  const optimisation_info info(simplex_fitdata(covar,pinfo.p,tfunc,data,pinfo.actord,pinfo.errs,weight,verbose,chistop,maxsteps));
  pinfo.copyback(p);
  return info;
}

double LowerBoundFunction::external_to_internal(double v) const { 
  const double x=v-a_+1.0;
  return std::sqrt(x*x-1.0);
}
  double LowerBoundFunction::internal_to_external(double v) const { return a_-1.0+std::sqrt(v*v+1.0); }

double UpperBoundFunction::external_to_internal(double v) const { 
  const double x=a_-v+1.0;
  return std::sqrt(x*x-1.0);
}
  double UpperBoundFunction::internal_to_external(double v) const { return a_+1.0-std::sqrt(v*v+1.0); }

LowerUpperBoundFunction::LowerUpperBoundFunction(double lowerv,double upperv)
  : lower_(lowerv), upper_(upperv)
{
  if (upper_<=lower_)
    throw InvalidParameter("LowerUpperBoundFunction");
  mid_=0.5*(lower_+upper_);
  range_=0.5*(upper_-lower_);
}

double LowerUpperBoundFunction::external_to_internal(double v) const { return asin((v-mid_)/range_); }

  double LowerUpperBoundFunction::internal_to_external(double v) const { return mid_+range_*std::sin(v); }

} //namespace libcmatrix
