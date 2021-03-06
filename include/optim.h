#ifndef LCM_OPTIM_H_
#define LCM_OPTIM_H_

#include "rmatrix.h"
#include "List.h"
#include "ScratchList.h"
#ifdef HAVE_LIBMINUIT
#if MINUITVER==2
#define LCM_MINUITNAMESPACE ROOT::Minuit2
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameters.h>
#else
#define LCM_MINUITNAMESPACE
#include <Minuit/FCNBase.h>
#include <Minuit/MnUserParameters.h>
#endif
#include <vector>
#endif
#include "smartptr.h"
#include "Warnings.h"

namespace libcmatrix {

  extern BaseWarning optimisation_base_warning;
  extern Warning<> fitsmallimprovement_warning;

  typedef void (*P_FITFUNC)(BaseList<double>, const BaseList<double>&, void*);
  typedef double (*P_MINFUNC)(const BaseList<double> &,void*);
 
  struct optimisation_info {
    double optimum;
    size_t iterations;
    bool converged;
    operator double() const { return optimum; }
  };

 struct BaseFitFunction {
   virtual ~BaseFitFunction() {}
   virtual void operator()(BaseList<double>, const BaseList<double>&) const =0;
 };

 struct BaseMinFunction {
   virtual ~BaseMinFunction() {}   
   virtual double operator()(const BaseList<double>&) const =0;
 };

class FitFunction : public BaseFitFunction {
 public:
  explicit FitFunction(P_FITFUNC func_, void* pblock_ =NULL) : func(func_), pblock(pblock_) {
    if (!func)
      throw InvalidParameter("FitFunction(): NULL function");
  }

   void operator()(BaseList<double> dest, const BaseList<double>& pars) const {
     (*func)(dest,pars,pblock);
   }
 private:
   P_FITFUNC func;
   void* pblock;
 };

 class MinFunction : public BaseMinFunction {
 public:
   explicit MinFunction(P_MINFUNC func_, void* pblock_ =NULL) : func(func_), pblock(pblock_) {
    if (!func)
      throw InvalidParameter("FitFunction(): NULL function");
   }

   double operator()(const BaseList<double>& pars) const 
   { return (*func)(pars,pblock); }

 private:
   P_MINFUNC func;
   void* pblock;
 };

class chi2func : public BaseMinFunction {
public:
  chi2func(const BaseFitFunction&, const BaseList<double>& data, const BaseList<double>& weights_);
  ~chi2func() {}

  double operator()(const BaseList<double>&) const;

private:
  const BaseFitFunction& fobj;
  const BaseList<double>& weights;
  mutable List<double> tdata; //chi2func object can't be shared between threads
  const BaseList<double>& fdata;
};

 class FitBase {
 public:
   virtual ~FitBase() {};
   virtual double next() =0;
   virtual double final(rmatrix&, List<double>&) =0;
 };

 // NB Object is not re-entrant; do not share between threads
  class LevMarqFit : public FitBase {
   const BaseFitFunction& fobj;
   BaseList<double> data;
   ScratchList<size_t> actord;
   BaseList<double> errs;
   ScratchList<double> p;

   List<double> datastore;
   BaseList<double> ytrial;
   BaseList<double> dy;
   rmatrix ytmp;
   ScratchList<double> ptry;
   ScratchList<double> pscr;
   ScratchList<double> beta;
   ScratchList<double> da;

   List< BaseList<double> > stores;
   rmatrix paras;
   bool asym;
   double lambda;
   ScratchList<double> w2;
   int verbose_;
   rmatrix alpha,tmp,grad;
   double ochisqr;
   bool havevalidspec; //indicates whether ytrial does contain func(p)
/* #ifdef HAVE_LIBPTHREAD */
/*    thread_controller tcon; */
/* #endif */

   void calc_grad(const BaseList<double>&);
   double initialise();
   //   double initialise_(); //debug
   double mqrcof(rmatrix&,BaseList<double>, const BaseList<double>&);
   void threads(int);
   
 public:
   LevMarqFit(double &,const BaseList<double> &p,const BaseFitFunction&, const BaseList<double> &data,const BaseList<size_t> &actord, const BaseList<double> &errs, const BaseList<double> &weights);
   LevMarqFit(double &,const BaseList<double> &p,const BaseFitFunction&, const BaseList<double> &data,const BaseList<size_t> &actord, const BaseList<double> &errs, double weight);
   //   LevMarqFit(double &,const BaseList<double> &p,const BaseFitFunction&, const BaseList<double> &data,const BaseList<size_t> &actord, const BaseList<double> &errs, double weight,bool); //debug
   void move();
   double evaluate();
   double next();
   double current_lambda() const { return lambda; }
   BaseList<double> current_trial();
   BaseList<double> current_params() const { return p; }
   BaseList<double> last_params() const { return ptry; }
   void verbose(int v_) { verbose_=v_; }
   void asymmetric(bool asym_) { asym=asym_; }
   double final(rmatrix&,BaseList<double>);
   
/*    void run(size_t nthr) { */
/*       BaseList<double> stored(stores(nthr)); */
/*       const BaseList<double> pars(paras.row(nthr)); */
/*       std::cout << "Thread " << nthr << " evaluating par set " << nthr << ": " << pars << std::endl; */
/*       fobj(stored,pars); */
/*    } */

   size_t total_parameters() { return errs.length(); }
   size_t variable_parameters() { return actord.length(); }
   size_t datapoints() { return data.length(); }
   double final(rmatrix& a,List<double>& p_) { p_.create(total_parameters()); return final(a,static_cast<BaseList<double>& >(p_)); }
/*    size_t threads() const { */
/* #ifdef HAVE_LIBPTHREAD */
/*      return tcon.get_max_threads(); */
/* #else */
/*      return 0; */
/* #endif */
/*    } */

 };

  #define LCM_DEFAULT_FIT_ITERATIONS 1000
  #define LCM_DEFAULT_SIMPLEX_ITERATIONS 50

void hessian(rmatrix &,const BaseFitFunction&,const BaseList<double> &p,const BaseList<size_t> &actord, const BaseList<double> &errs,const BaseList<double> &weights, int verbose =0);

void covariance(rmatrix &,const BaseFitFunction&,const BaseList<double> &p,const BaseList<size_t> &actord, const BaseList<double> &errs,const BaseList<double> &weights, int verbose =0);

optimisation_info fitdata(rmatrix &covar, BaseList<double> p, const BaseFitFunction&, const BaseList<double> &data,const BaseList<size_t> &actord, const BaseList<double> &errs, const BaseList<double> &weights,int verbose,double chistop,int maxsteps =LCM_DEFAULT_FIT_ITERATIONS);
  
  optimisation_info fitdata(rmatrix &covar, BaseList<double> p, const BaseFitFunction&, const BaseList<double> &data,const BaseList<size_t> &actord, const BaseList<double> &errs, double weight,int verbose,double chistop,int maxsteps =LCM_DEFAULT_FIT_ITERATIONS);

  optimisation_info simplex(BaseList<double> p0,const BaseMinFunction&, const BaseList<size_t>& actord, const BaseList<double>& errs,int verbose,double stop_val,int maxstep);

  optimisation_info simplex_fitdata(rmatrix&, BaseList<double> p,const BaseFitFunction&, const BaseList<double>& data,const BaseList<size_t>&,  const BaseList<double> &errs, const BaseList<double> &weights,int verbose, double chistop, int maxsteps =LCM_DEFAULT_SIMPLEX_ITERATIONS);

  optimisation_info simplex_fitdata(rmatrix &,BaseList<double> p,const BaseFitFunction&, const BaseList<double> &data,const BaseList<size_t> &, const BaseList<double> &errs, double weight,int verbose, double chistop, int maxsteps =LCM_DEFAULT_SIMPLEX_ITERATIONS);
  
  extern Warning<> optim_covariance_warning;

//New Minuit interface
#ifdef HAVE_LIBMINUIT
  template<class F> class MinuitAdaptor : public LCM_MINUITNAMESPACE::FCNBase {
  public:
    explicit MinuitAdaptor(const F& func) : func_(func) {}
#if MINUITVER==2
    double Up() const { return 1.0; }
#else
    double up() const { return 1.0; }
#endif
    double operator()(const std::vector<double>&) const;
  private:
    F func_;
  };
  
  template<class F> double MinuitAdaptor<F>::operator()(const std::vector<double>& in) const 
  {
    const double* asvec(&(in.front()));
    const size_t n=in.size();
    if (&(in.back())-asvec==n-1) //check that storage is contiguous
      return func_(BaseList<double>(n,const_cast<double*>(asvec)));
    //slow fall-back
    List<double> aslist;
    aslist.create(n,in.begin());
    return func_(aslist);
  }
#endif

  //Extended variable interface

  //! (abstract) base class from bound functions
  struct BaseBoundFunction {
    virtual ~BaseBoundFunction() {}
    virtual double external_to_internal(double) const =0;
    virtual double internal_to_external(double) const =0;
    virtual bool isvalid(double) const =0;
    virtual BaseBoundFunction* clone() const =0;
  };

  //Bound functions matching MINUIT usage

  //! constrain lower bound
  class LowerBoundFunction : public BaseBoundFunction {
  public:
    LowerBoundFunction(double av) : a_(av) {}
    double external_to_internal(double) const;
    double internal_to_external(double) const;
    BaseBoundFunction* clone() const { return new LowerBoundFunction(*this); }
    bool isvalid(double v) const { return (v>=a_); }
  private:
    double a_;
  };

  //! constrain upper bound
  class UpperBoundFunction : public BaseBoundFunction {
  public:
    UpperBoundFunction(double av) : a_(av) {}
    double external_to_internal(double) const;
    double internal_to_external(double) const;
    BaseBoundFunction* clone() const { return new UpperBoundFunction(*this); }
    bool isvalid(double v) const { return (v<=a_); }
  private:
    double a_;
  };

  //! constrain parameter range
  class LowerUpperBoundFunction : public BaseBoundFunction {
  public:
    LowerUpperBoundFunction(double lowerv,double upperv);
    double external_to_internal(double) const;
    double internal_to_external(double) const;
    BaseBoundFunction* clone() const { return new LowerUpperBoundFunction(*this); }
    bool isvalid(double v) const { return (v<=upper_) && (v>=lower_); }
  private:
    double lower_,upper_,mid_,range_;
  };
  
  class Parameter {
  public:
    Parameter(const char* name, double valv, double errv);
    Parameter(double valv, double errv);

    const std::string& name() const { return name_; }
    std::string& name() { return name_; }
    bool isnamed() const { return (*(name_.c_str())!='\0'); } //!< return \c true if name has been set

    void fix() { isfixed_=true; } //!< fix parameter value
    void release() { isfixed_=false; } //!< release parameter
    bool isfixed() const { return isfixed_; } //!< return \c true if value fixed

    double operator()() const { return externalvalue_; }
    double get() const { return externalvalue_; } //!< return current (external) value
    double get_internal() const { return internalvalue_; } //!< return current translated value
    void set(double);
    void set_internal(double);

    void error(double);
    double error() const { return error_; }

    void constrain(const BaseBoundFunction&); //!< use constraint function (copied)
    void unconstrain(); //!< remove any constraint function
    bool isconstrained() const { return !!constrainfuncp_; } //!< \c true if value is constrained
    double internal_to_external(double v) const { return !constrainfuncp_ ? v : constrainfuncp_->internal_to_external(v); }
    double external_to_internal(double v) const { return !constrainfuncp_ ? v : constrainfuncp_->external_to_internal(v); }
    bool isvalid(double v) const { return !constrainfuncp_ ? true : constrainfuncp_->isvalid(v); }

  private:
    std::string name_;
    double externalvalue_;
    double internalvalue_;
    double error_;
    bool isfixed_;
    smartptr<BaseBoundFunction> constrainfuncp_;
    void setraw(double);
  };

  class SimpleBoundsState {
  public:
    SimpleBoundsState() : havestate_(NONE) {}
    SimpleBoundsState(double minv, double maxv) : havestate_(NONE) { lowerupper(minv,maxv); }
    
    bool unconstrain(); //!< release 
    bool lower(double); //!< set lower bound (only), \c true if updated
    bool upper(double); //!< set upper bound (only)
    bool lowerupper(double,double); //!< set double bound
    bool islowerset() const { return (havestate_ & MIN); }
    bool isupperset() const { return (havestate_ & MAX); }
    bool isconstrained() const { return (havestate_==NONE); }
    double lower() const;
    double upper() const;
    const BaseBoundFunction* function() const { return funcp.get(); }
    enum state_t { NONE=0, MIN=1, MAX=2, BOTH=3 };
#ifdef HAVE_LIBMINUIT
    void apply(LCM_MINUITNAMESPACE::MnUserParameters&, size_t) const; //!< apply constraint to MINUIT parameter
#endif
    
  private:
    state_t havestate_;
    double min_,max_;
    smartptr<BaseBoundFunction> funcp;
    
    template<class T> static void update(T& dest, const T val, bool& needup) {
      if (dest!=val) {
	dest=val;
	needup=true;
      }
    }
    void doupdate(); //!< update stored function pointer
  };

  std::ostream& operator<< (std::ostream&, const SimpleBoundsState&);

#ifdef HAVE_LIBMINUIT
  // rather large for inline, but means library doesn't need recompilation for MINUIT support
  inline void SimpleBoundsState::apply(LCM_MINUITNAMESPACE::MnUserParameters& pars, size_t i) const {
    switch (havestate_) {
    case NONE:
#if MINUITVER==2
      pars.RemoveLimits(i);
      break;
    case MIN:
      pars.SetLowerLimit(i,min_);
      break;
    case MAX:
      pars.SetUpperLimit(i,max_);
      break;
    case BOTH:
      pars.SetLimits(i,min_,max_);
#else
      pars.removeLimits(i);
      break;
    case MIN:
      pars.setLowerLimit(i,min_);
      break;
    case MAX:
      pars.setUpperLimit(i,max_);
      break;
    case BOTH:
      pars.setLimits(i,min_,max_);
#endif
      break;
    }
  }
#endif

  //Extended interface versions

  void hessian(rmatrix& d,const BaseFitFunction& func, const BaseList<Parameter>& p, const BaseList<double>& weights, int verbose =0, bool ignoreconstraints =false);
  void covariance(rmatrix& d,const BaseFitFunction& func, const BaseList<Parameter> &p,const BaseList<double> &weights, int verbose =0, bool ignoreconstraints =false);
    
  optimisation_info fitdata(rmatrix &covar, BaseList<Parameter> p, const BaseFitFunction& func, const BaseList<double> &data,const BaseList<double> &weights,int verbose,double chistop,int maxsteps =LCM_DEFAULT_FIT_ITERATIONS);

  optimisation_info fitdata(rmatrix &covar, BaseList<Parameter> p, const BaseFitFunction& func, const BaseList<double> &data, double weight,int verbose,double chistop,int maxsteps =LCM_DEFAULT_FIT_ITERATIONS);

  optimisation_info simplex(BaseList<Parameter> p0,const BaseMinFunction& func, int verbose,double stop_val,int maxstep);

  optimisation_info simplex_fitdata(rmatrix& covar, BaseList<Parameter> p,const BaseFitFunction& func, const BaseList<double>& data, const BaseList<double>& weights, int verbose, double chistop, int maxsteps =LCM_DEFAULT_SIMPLEX_ITERATIONS);

  optimisation_info simplex_fitdata(rmatrix& covar, BaseList<Parameter> p,const BaseFitFunction& func, const BaseList<double>& data, double weight,int verbose, double chistop, int maxsteps =LCM_DEFAULT_SIMPLEX_ITERATIONS);

} //namespace libcmatrix

#endif
