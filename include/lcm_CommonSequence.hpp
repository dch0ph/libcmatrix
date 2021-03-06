#ifndef lcm_CommonSequence_hpp_
#define lcm_CommonSequence_hpp_

/*! \file
  \brief  Objects common to sequence propagator generators */

#include <list>
#include <map>
#include "lcm_FunctionObject.h"
#include "BlockedMatrix.h"
#include "HamiltonianStore.h"
#include "UnionHolder.h"
#include "BaseMetaPropagation.h"
#include "simple_counter.h"
#include "Warnings.h"

namespace libcmatrix {

 class Validator {
  public:
    Validator() : value(0) {}
    typedef long type;
    void next() { value++; }
    type operator()() const { return value; }
   static const long invalid=-1; //!< "guaranteed" to be invalid
  public:
    type value;
  };

  template<class T> struct PulseGenerator_traits {};

  struct transient_state {
    BlockedMatrix<complex> U;
    double tip,phase;

    bool empty() const { return !U; }
    void clear() { U.clear(); }
  };

  class PulseGeneratorBase {
  public:
   typedef BlockedMatrix<complex> operator_type;
   typedef ListList<complex> diag_operator_type;

    template<class Sys> explicit PulseGeneratorBase(const Sys&, nuclei_spec =NULL_NUCLEUS, int =0, simple_counter* =NULL);
    template<class Sys> explicit PulseGeneratorBase(const Sys&, const operator_type&, nuclei_spec =NULL_NUCLEUS, int =0, simple_counter* =NULL);
    virtual ~PulseGeneratorBase() {}

    virtual void reset() {}

    virtual PulseGeneratorBase* clone() const =0;

    typedef RealComplexHolder< BlockedMatrix<double>, BlockedMatrix<complex> > RCblockedmatrix;

    virtual void transient_amplitude(double) =0;
    bool hastransients() const { return (transamp_!=0.0); }
    virtual double transient_duration() const { return 0.0; }
    bool hassystem() const { return hassys_; }
   const operator_type& Hsystem() const { return Hsys_; }
    void Hsystem(const operator_type& Hsysv); //NB invalidates any RFEvents using this pgen

    //non-zero if unprocessed transient
    virtual bool isfinished() const { return true; }
    virtual void clear_transient() {}
    virtual void clear_transient(operator_type&, int =-1) {}
    virtual void store_transient(transient_state&) const {}
    virtual void restore_transient(const transient_state&) { throw Failed("PulseGenerator: can't set transients for this type"); }

    void operator()(operator_type&, double,phase_spec) const;

    virtual void operator()(operator_type&,const operator_type& H,double time,double offset,double angle,phase_spec, int =false) const =0;

    //no RF; default is pulse of zero amplitude
    virtual void operator()(operator_type& U, const operator_type& H, double dt) const
    { operator()(U,H,dt,0.0,0.0,0.0); } 

    void operator()(cmatrix& U,double time,double offset,double angle,phase_spec pspec, int =0) const;

    virtual void operator()(operator_type& U,double time,double offset,double angle,phase_spec, int =0) const =0;

    void operator()(cmatrix& U, double tip, phase_spec pspec) const;

   cmatrix operator()(double tip,phase_spec pspec) const {
     cmatrix U(size_,size_,mxflag::temporary); 
     (*this)(U,tip,pspec);
     return U;
   }

    void operator()(cmatrix& U,const cmatrix& H,double time,double offset,double angle,phase_spec pspec, int =0) const;

    cmatrix operator()(const cmatrix& H,double time,double offset,double angle,phase_spec pspec, int flags =0) const {
     cmatrix U(size_,size_,mxflag::temporary);
     (*this)(U,H,time,offset,angle,pspec,flags);
     return U;
   }

    cmatrix operator()(double time,double offset,double angle,phase_spec pspec, int flags) const {
     cmatrix U(size_,size_,mxflag::temporary);
     (*this)(U,time,offset,angle,pspec,flags);
     return U;
   }

    virtual void H_and_Us(operator_type&, operator_type&, diag_operator_type&, double,double,double,double, int) const =0;
    // { throw Failed("PulseGenerator::H_and_Us invalid with time-dependent transients"); }

    void apply_flags(operator_type& U,double time,double offset,double angle,double phase, int flags, int which =-1) const;
    void apply_flags(operator_type& U,double vRF,double offset,double phase, int flags, int which =-1) const;

    static Warning<> partial_incompatible_warning;

   void correct_partial(operator_type& U,double time,double offres, int which =-1) const;

    void correct_partial(cmatrix& U, double time, double offres) const;
    simple_counter* cache_counter() { return cachep_; }
    const simple_counter* cache_counter() const { return cachep_; }

    int verbose() const { return verbose_; }
    const ListList<double>& Fz() const { return Fz_; }
    const RCblockedmatrix& Fx() const { return Fx_; }
    Validator::type validation_key() const { return validator_(); }    

    void rotatez_ip_cached(operator_type&, double) const; //!< apply z-rotation of \a p using (possibly) cached values

    const SpinOpGeneratorBase& generator() const { return opgen_; }
    const matrix_partition_set& partitioning() const { return partitions_; }

  protected:
    operator_type Hsys_;
    bool hassys_;
    size_t nuc_;
    const SpinOpGeneratorBase& opgen_;
    int verbose_;
    simple_counter* cachep_;
    ListList<double> Fz_;
    RCblockedmatrix Fx_;
    size_t mzblocks_,eigblocks_,totblocks_;
    size_t size_;
    Validator validator_;
    matrix_partition_set partitions_;
    double transamp_;

    mutable Mutex<ThreadingActive> cache_lock_;    
    //caches eigenbasis (mutex protected)
    mutable ListList<double> eigs;
    mutable RCblockedmatrix V;

    void ensure_eigenbasis() const;

    void rawpulsegen(operator_type&,const operator_type& H,double time,double offset,double angle,double, int =0) const;
    void H_and_U(operator_type&, diag_operator_type&, double,double,double,double, int) const;
    void make_transient(BlockedMatrix<complex>&,double time,double angle, double phase) const;
    void make_transient(BlockedMatrix<complex>&,double vRF, double phase) const;

    //Used if single-threaded
    mutable ListList<complex> Utmp_;
    mutable ListList<complex> ceigs;

  private:
    template<class Sys> void common_init(const Sys&, nuclei_spec);
    void common_init_();
    //    void raw_correct_partial(operator_type& U,double phasev, int which =-1) const;
    void apply_raw_transient(operator_type& U, const operator_type& Utrans, int which =-1) const;
  };

  class ZshiftCache {
  public:
    ZshiftCache(const ListList<double>& Fzv, int verbosev =0)
      : Fz_(Fzv), verbose_(verbosev) {}
    typedef std::map<double,ListList<complex> > store_t;
    const ListList<complex>& operator()(double) const;
    void operator()(BlockedMatrix<complex>& U, double p) const { rotatez_ip(U,(*this)(p)); }
  private:
    mutable Mutex<ThreadingActive> lock_;    
    mutable store_t cache_;
    const ListList<double> Fz_;
    int verbose_;
  };
    
  class PulseGenerator : public PulseGeneratorBase {
 public:
    template<class Sys> explicit PulseGenerator(const Sys& sys, nuclei_spec whichn =NULL_NUCLEUS, int verbosev =0, simple_counter* cachepv =NULL) : PulseGeneratorBase(sys,whichn,verbosev,cachepv) {}
    template<class Sys> explicit PulseGenerator(const Sys& sys, const operator_type& H, nuclei_spec whichn =NULL_NUCLEUS, int verbosev =0, simple_counter* cachepv =NULL) : PulseGeneratorBase(sys,H,whichn,verbosev,cachepv) {}
    
    PulseGeneratorBase* clone() const { return new PulseGenerator(*this); }

    using PulseGeneratorBase::operator();

    void operator()(operator_type& U, const operator_type& H, double dt) const { propagator(U,H,dt); }

    void operator()(operator_type&,const operator_type& H,double time,double offset,double angle,phase_spec,int =0) const;
    void operator()(operator_type& U,double time,double offset,double angle,phase_spec,int =0) const;

    //    void operator()(operator_type&, double tip,phase_spec pspec) const;
						
   //Not strictly necessary, but avoids copying cached values
    PulseGenerator(const PulseGenerator& a) : PulseGeneratorBase(static_cast<const PulseGeneratorBase&>(a)) {}

    void transient_amplitude(double);

 private:

    void H_and_Us(operator_type&, operator_type&, diag_operator_type&, double,double,double,double,int =0) const;

    mutable double lastangle_;
    mutable double lastoffset_;
    mutable double lastdur_;
    mutable double lasttrans_;
    mutable operator_type Ustash_;
  };

  template<> struct PulseGenerator_traits<PulseGenerator> {
    static const bool allowideal=true;
    static const bool timedependent=false;
  };

  class RFEvent {
  public:
    RFEvent(double t, double nomdt, PulseGeneratorBase* pgenp =NULL)
      : pgenp_(pgenp), duration_(t), nomduration_(nomdt), isvalid_(!pgenp) {
      if ((t<0) || (nomduration_<0))
	throw InvalidParameter("RFEvent: duration cannot be negative");
    }
    explicit RFEvent(PulseGeneratorBase* pgenp =NULL)
      : pgenp_(pgenp), duration_(0.0), nomduration_(0.0), isvalid_(!pgenp) {}

    virtual ~RFEvent() {}

    typedef PulseGeneratorBase::operator_type operator_type;
    typedef PulseGeneratorBase::diag_operator_type diag_operator_type;

    virtual void add_Hrf(operator_type&) const {}

    double duration() const { return duration_; }
    virtual void duration(double) { throw Failed("duration: can't set duration for this RFEvent"); }
    void nominal_duration(double); //NB changes *nominal* duration
    virtual bool isnull() const =0;

    PulseGeneratorBase& pulse_generator() { 
      if (pgenp_)
	return *pgenp_;
      throw Failed("RFEvent: no pulse generator"); }

    const PulseGeneratorBase& pulse_generator() const {
      if (pgenp_)
	return *pgenp_;
      throw Failed("RFEvent: no pulse generator"); }
    bool haspulse_generator() const { return (pgenp_!=NULL); }
    
    virtual RFEvent* clone() const =0;

    virtual void transforms() =0;
    bool hassystem() const { return pulse_generator().hassystem(); }
    virtual bool needs_correct() const { return false; }
    virtual void correct_partial(operator_type&, double,int) const
    { throw InternalError("correct_partial shouldn't be called for this RFEvent"); }
    virtual void apply_event(operator_type&, operator_type&, int) const =0;
    virtual void apply_start(operator_type&, int) const
    { throw InternalError("apply_start shouldn't be called for this RFEvent"); }
    virtual void apply_end(operator_type&, operator_type&, int) const
    { throw InternalError("apply_end shouldn't be called for this RFEvent"); }

    virtual double offset() const { return 0.0; }
    virtual void offset(double)
    { throw Failed("offset cannot be set for this RFEvent"); }

    virtual void print(std::ostream&) const =0;
    virtual void synchronise(double &,double &,double, char, double, const UnaryFunction<double,double>&);

    void ensurevalid() const;
    void flagdirty() { isvalid_=false; }
    virtual bool isvalid() const { return isvalid_; }
    
  protected:
    PulseGeneratorBase* pgenp_;
    double duration_;
    double nomduration_;
  private:
    bool isvalid_;
    Validator::type pgenkey_;
  };

  inline std::ostream& operator<< (std::ostream& ostr,const RFEvent& a) {
    a.print(ostr);
    return ostr;
  }

  class SoftPulseBase;
  class HardPulse;

  class RFPulseEvent : public RFEvent {
protected:

    RFPulseEvent(PulseGeneratorBase& pgen, double t,double nomdt,double rfv,double phasev, int flagsv =0)
      : RFEvent(t,nomdt,&pgen),
	flags_(flagsv),
	rf_(rfv), phase_(phasev) {}

    int flags_;
  double rf_;
  double phase_;

    struct propagator_set {      
      propagator_set() : cachep(NULL) {}
      ~propagator_set() { clear(); }
      void create(const SoftPulseBase&, double, double);
      void create(const HardPulse&);
      void check_cache(const RFPulseEvent&);
      bool empty() const { return Utotal.empty() && Hrf.empty(); }
      void clear();

      operator_type Ubegin;
      operator_type Utotal;
      operator_type Hrf;
      diag_operator_type Uend;

      simple_counter* cachep;
      bool destroy_when_finished;

      void apply_start(operator_type&, int which) const;
      void apply_event(operator_type&, operator_type&, int, bool) const;
      void apply_end(operator_type&, operator_type&, int, bool) const;
      void add_Hrf(operator_type& H) const {      
	if (!!Hrf) //!< Hrf may be empty if rf_ is zero
	  H+=Hrf;
      }
      void print(std::ostream&) const;
      size_t used() const;

      friend class SoftPulseBase;
      friend class HardPulse;
    };

      propagator_set Uset;

public:
    enum { COHERENT=1, TRANSIENTS=2 };

  double phase() const { return phase_; }
    void phase(phase_spec pspec);
        
    double rf() const { return rf_; }
    void rf(double);

    bool isvalid() const { return RFEvent::isvalid() && (!Uset.empty() || isnull()); }

    double tip() const { return 2.0*M_PI*nomduration_*rf_; }
    bool isnull() const { return (nomduration_*rf_==0.0); }
    int flags() const { return flags_; }
    
    void apply_start(operator_type& U, int which) const;
    void apply_event(operator_type& U, operator_type& H, int which ) const;
    void apply_end(operator_type& U, operator_type& H, int which) const;
    void add_Hrf(operator_type& H) const {      
      ensurevalid();
      Uset.add_Hrf(H);
    }
  };

  class RFChangeOffset : public RFEvent {
  public: 
    RFChangeOffset(double offsetv, PulseGeneratorBase* pgenp =NULL)
      : RFEvent(pgenp), offset_(offsetv) {}

    RFEvent* clone() const { return new RFChangeOffset(*this); }
    void print(std::ostream&) const;
    
    void transforms() { multiply(Hrf,offset_,pulse_generator().Fz()); }
    void apply_event(operator_type&, operator_type& H, int) const 
    { H+=Hrf; }

    double offset() const { return offset_; }
    void offset(double);

    bool isnull() const { return (offset_==0.0); }
  private:
    double offset_;
    ListList<double> Hrf;
  };


class HardPulse : public RFPulseEvent {
public:
  HardPulse(PulseGeneratorBase& pgen, double t,double rfv, phase_spec pspec)
    : RFPulseEvent(pgen,0.0,t,rfv,pspec()) {}

  void print(std::ostream&) const;

  RFEvent* clone() const { return new HardPulse(*this); }

 protected:

private:
  void transforms() {  Uset.create(*this); }
};

class SoftPulseBase : public RFPulseEvent {
public:  
  double offset() const { return offres_; }
  void offset(double);

protected:
  SoftPulseBase(PulseGeneratorBase& pgen, double t,double rfv,phase_spec pspec,double offresv, int flagsv)
    : RFPulseEvent(pgen,t,t,rfv,pspec(),flagsv),
      offres_(offresv) {}

  void duration(double);    
  bool needs_correct() const { return (offres_!=0.0) && ((flags_ & COHERENT)==0); }
  void correct_partial(operator_type&, double,int) const;

  double offres_;
};

class SoftPulse : public SoftPulseBase {
public:  
  SoftPulse(PulseGeneratorBase& pgen,double t,double rfv,phase_spec pspec,double offresv =0, int flagsv =0)
    : SoftPulseBase(pgen,t,rfv,pspec(),offresv,flagsv) {}

  void print(std::ostream &) const;

  RFEvent* clone() const { return new SoftPulse(*this); }

private:
  void transforms();
};

class CWPulse : public SoftPulseBase {
public:
  CWPulse(PulseGeneratorBase& pgen,double t,double rfv,phase_spec pspec, double offresv =0.0, int flagsv =0)
    : SoftPulseBase(pgen,t,rfv,pspec(),offresv,flagsv), actduration_(t) {}

  void synchronise(double& ,double& ,double, char,double, const UnaryFunction<double,double>&); //overwrite default sync method
  void print(std::ostream &) const;

  RFEvent* clone() const { return new CWPulse(*this); }

 private:
  double actduration_;
  void transforms();
};

  struct TimedEvent {
    TimedEvent(RFEvent* eventp_,double timev,char syncv ='|')
      : eventp(eventp_), start_(0.0), end_(-1.0)  {
      time(timev);
      if (!eventp_)
	throw InvalidParameter("TimedEvent: null RFEvent pointer");
      if (syncv!='|' && eventp->duration())
	// std::cerr << "TimedEvent warning: synchronisation guide is ignored for finite duration events\n";
	syncv='|'; //distinguish sync of sequence building from final application
      
      sync_=syncv;
    }
    
    double time() const { return time_; }
    double duration() const { return end_-start_; }
    void time(double timev) {
      if (timev<0.0)
	throw InvalidParameter("TimedEvent::time: cannot be negative");
      time_=timev;
    }
    double start() const { return start_; }
    double end() const { return end_; }
    bool inevent(double timev) const {
      return ((timev>start_) && (timev<end_));
    }
    RFEvent* event() { return eventp; }
    const RFEvent* event() const { return eventp; }
    
    RFEvent* eventp;
    double time_;
    char sync_;
    
    double start_;
    double end_;
    
    void dump() const;
  };

  std::ostream& operator << (std::ostream&, const TimedEvent&);

  class RoundFunc : public UnaryFunction<double,double> {
  public:
    RoundFunc(double resv =0.0, double basev =0.0, double tolv =1e-9) { resolution(resv,basev,tolv); }
    void resolution(double resv, double basev, double tolv) {      
      if (resv<0.0)
	throw InvalidParameter("RoundFunc::resolution");
      res_=resv;
      base_=basev;
      tol_=tolv;
    }
    double operator()(double t) const {
      return res_ ? base_+res_*ceil((t-tol_-base_)/res_) : t;
    }
  private:
    double res_,base_,tol_;
  };

  class Spectrometer {
  public:
    Spectrometer()
      : tres_(0.0) {}
    double time_resolution() const { return tres_; }
    void time_resolution(double tresv) {
      if (tresv<0.0)
	throw InvalidParameter("Spectrometer::time_resolution");
      tres_=tresv;
    }
  private:
    double tres_;
  };
  

  class Sequence : public ::std::list<TimedEvent> {
  public:
    Sequence(double periodv, const Spectrometer& spectro =Spectrometer()) { 
      init(spectro);
      period(periodv);
    }
    
   Sequence(const Spectrometer& spectro =Spectrometer())
       : period_(0.0), scale_(1.0)
   { init(spectro);
     reset(); applyscale(); }

    TimedEvent& push_back(RFEvent* ev, double start,char sync ='|');
    void push_back(const BaseList<RFEvent*>&,double start,char ='|');
    void push_back(const Sequence& a) { 
      insert(end(),a.begin(),a.end());
      period(period_); //!< not sure why this is here!
    }
   void time_resolution(double res, double tol, double base =0.0) { 
     if (!empty())
       throw Failed("time_resolution: can't be changed part way through sequence");
     roundfunc_.resolution(res,base,tol);
   }
    double push_delay(double,double);
   void period(double);
   double period() const { return period_; }
   double duration() const;

   void autosync(double,double);
   void autosync(double tol) { autosync(tol,duration()-tol); }

   PulseGeneratorBase* pulse_generator() const { return pgenp_; }

    static BaseWarning sync_warning; //!< base for sync warnings
    static Warning<> truncated_transient_warning;
    static Warning<> boundary_event_warning;

 private:
   double period_;
   double scale_;
   double lastend_;
   bool havepgen_;
   PulseGeneratorBase* pgenp_;
   RoundFunc roundfunc_;

   void init(const Spectrometer&);
   void reset() { havepgen_=false; pgenp_=NULL; lastend_=0.0; }
   void check_flush(double);
   void syncadd(TimedEvent&);
   void applyscale();
   void dump() const;
 };

  std::ostream& operator<< (std::ostream&, const Sequence&);

  // implementation details

  template<class Sys> PulseGeneratorBase::PulseGeneratorBase(const Sys& a,const operator_type& Hsysv, nuclei_spec whichn, int verbosev, simple_counter* cachepv)
    : Hsys_(Hsysv), hassys_(true), opgen_(a.generator()), verbose_(verbosev), cachep_(cachepv)
    {
      common_init(a,whichn);
    }

   template<class Sys> PulseGeneratorBase::PulseGeneratorBase(const Sys& a,nuclei_spec whichn, int verbosev, simple_counter* cachepv) 
     : hassys_(false), opgen_(a.generator()), verbose_(verbosev), cachep_(cachepv)
    {
      common_init(a,whichn);
    }

  template<class Sys> void PulseGeneratorBase::common_init(const Sys& a,nuclei_spec whichn)
    {
      nuc_=whichn();
      if (nuc_==NULL_NUCLEUS)
	nuc_=spinsystem_glue_<Sys>::spinsystem(a).homonucleus();
      Fz_=a.diag_Fz(nuc_);
      if (Fz_.empty())
	throw Failed("PulseGenerator: nucleus not present in spin system!");
      block_pattern blkspec;
      if (spinsystem_glue_<Sys>::ispurereal) {
	BlockedMatrix<complex> tmp;
	a.mla(tmp,blkspec,1.0,operator_spec(whichn,'x'));
	real(Fx_.set_real(),tmp);
      }
      else
	a.mla(Fx_.set_complex(),blkspec,1.0,operator_spec(whichn,'x'));
 
      if (!blkspec.isdiagonal())
	throw Failed("PulseGenerator: RF channel cannot have mz blocking"); 
      common_init_();
   }
 
}

#endif
