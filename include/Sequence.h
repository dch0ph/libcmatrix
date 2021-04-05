#ifndef _Sequence_h_
#define _Sequence_h_

/*! \file
  \brief  General sequence RFEvents and propagator generators */

#include "NMR.h"
#include "UnionHolder.h"
#include "lcm_CommonSequence.hpp"
#include "MetaPropagation.h"

namespace libcmatrix {

  extern BaseWarning lcm_sequence_warning; //!< base warning for Sequence related events

//   class PulseGenerator_SimpleTransients : public PulseGeneratorBase {
//  public:
//     template<class Sys> explicit PulseGenerator_SimpleTransients(const Sys& sys, nuclei_spec whichn =NULL_NUCLEUS, int verbosev =0)
//       : PulseGeneratorBase(sys,whichn,verbosev) { transient_amplitude(0.0); }

//     bool hastransients() const { return true; }
//     void reset() { trans_state_.clear(); }
//     bool isfinished() const { return trans_state_.empty(); }
//     void clear_transient();
//     void clear_transient(operator_type&, int =-1);
//     void store_transient(transient_state& store) const { store=trans_state_; }
//     void restore_transient(const transient_state& store) { trans_state_=store; }

//     PulseGeneratorBase* clone() const { return new PulseGenerator_SimpleTransients(*this); }

//     void operator()(operator_type&,const operator_type& H,double time,double offset,double angle,phase_spec) const;

//     void operator()(operator_type&,double,double,double,phase_spec) const {
//       throw InternalError("PulseGenerator_SimpleTransients: Hamiltonian unspecified");
//     }

//     void H_and_Us(operator_type&, operator_type&, diag_operator_type&, double time,double offres,double angle,double phase,bool) const;

//  private:
//     mutable transient_state trans_state_;
//     mutable BlockedMatrix<complex> Unewtrans_,UBtmp_;

//     void make_transient(BlockedMatrix<complex>&,double time,double offset, double angle, double phase) const;
//     void prop(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, int =-1) const;   
//  };

//   template<> struct PulseGenerator_traits<PulseGenerator_SimpleTransients> {
//     static const bool allowideal=false;
//     static const bool timedependent=false;
//   };

  template<class Sys> class BaseSequencePropagator;

  struct sequence_state {
    enum { FINISHED, INACTIVE, ACTIVE, TRANSIENT };

    sequence_state() : pgenp_(NULL), state(INACTIVE) {}
    
    PulseGeneratorBase* pgenp_;
    int state;
    double cycletime;
    double evend;
    Sequence::iterator iter;
    Sequence::iterator iterend;

    void terminate() {
      state = isfinished() ? FINISHED : TRANSIENT;
    }
    //used to check for residual transients
    bool isfinished() const { return pgenp_ ? pgenp_->isfinished() : true; }
    void reset(Sequence&,double);
    void clear_transient(BlockedMatrix<complex>&, int);
    void print(std::ostream&) const;
  };

  struct BaseSequencePropagator_  : public BaseMetaPropagator 
  {

    template<class TypeH> BaseSequencePropagator_(const TypeH& Hsysv,double intdtv, double originv, const BaseList<Sequence>& seqsv, const BaseList<double>& periodsv,double tolv, int verbosev, int flagsv)
     : BaseMetaPropagator(Hsysv,verbosev),
       intdt_(intdtv), origin_(originv),
       seqs_(seqsv),
       periods_(periodsv),
       uniqueperiod_(-1.0),
       synchint_(0.0)
    { create(tolv); }

    template<class TypeH> BaseSequencePropagator_(const TypeH& Hsysv, double intdtv, double originv, const BaseList<Sequence>& seqsv, double periodv,double tolv, int verbosev, int flagsv)
      : BaseMetaPropagator(Hsysv,verbosev,flagsv),
	intdt_(intdtv), origin_(originv),
	seqs_(seqsv),
	uniqueperiod_(periodv),
	synchint_(0.0)
    { create(tolv); }

    double origin() const { return origin_; }
    void origin(double originv) { origin_=originv; }
    double period(size_t chan) const { return periods_.empty() ? uniqueperiod_ : periods_(chan); }
    double period() const { return periods_.empty() ? uniqueperiod_ : 0.0; }
    bool hasuniqueperiod() const { return periods_.empty(); }
    void synchronisation_hint(double);
    double synchronisation_hint() const { return synchint_; }

    void propagators(cmatrix&, size_t, double, double,size_t,size_t) const {
      throw Failed("Can't evaluate diagonal propagators for RF sequence");
    }

    typedef cmatrix operator_type;
    
    void operator()(BlockedMatrix<complex>& U, double t1, double t2) const {
      calc_U(U,t1,t2,-1);
    }

    void operator()(cmatrix& U,double t1,double t2, size_t mzeig, size_t eig) const;
    
    cmatrix operator() (double t1,double t2, size_t mzblk =0, size_t eigblk =0) {
      cmatrix U(mxflag::temporary);
      operator()(U,t1,t2,mzblk,eigblk);
      return U;
    }
    
    const double intdt_;
    double origin_;
    BaseList<Sequence> seqs_;
    ScratchList<double> periods_;
    double uniqueperiod_;
    double synchint_;

    size_t channels;
    int offreslev;
    double lasttime_,tol;
    bool started;
    bool allowquick_;

    List<sequence_state> status;
    BlockedMatrix<complex> Hcur;
    mutable BlockedMatrix<complex> UBtmp_seq_;
    mutable CommonHolder<cmatrix,BlockedMatrix<complex> > UtmpSeqs_;

    void advance(size_t);
    void process(size_t);
    void reset(size_t,double);
    void reset(double);
    void create(double);
    bool isfinished() const;

    void print(std::ostream&) const;
    void calc_U(BlockedMatrix<complex>&, double, double,int which) const;
    void prop_U_(BlockedMatrix<complex>&, double, double,int which) const;
    void delay_propagate(BlockedMatrix<complex>& U, double dt,int which =-1);

    virtual void genpropagator(BlockedMatrix<complex>& U, double t1, double t2, int which) const =0;
    virtual bool Hisconstant() const =0;
    virtual double Hperiod() const =0;
  };

  template<class Sys> struct BaseSequencePropagator : public BaseSequencePropagator_ {
   typedef typename Sys::result_type baseH_type;

    BaseSequencePropagator(const Sys& Hsysv, double intdtv, double originv, const BaseList<Sequence>& seqsv, const BaseList<double>& periodsv,double tolv, int verbosev, int flagsv)
      : BaseSequencePropagator_(Hsysv,intdtv,originv,seqsv,periodsv,tolv,verbosev,flagsv),
	Hsys_(Hsysv) { create(); }

    BaseSequencePropagator(const Sys& Hsysv, double intdtv, double originv, const BaseList<Sequence>& seqsv, double periodv,double tolv, int verbosev, int flagsv)
      : BaseSequencePropagator_(Hsysv,intdtv,originv,seqsv,periodv,tolv,verbosev,flagsv),
	Hsys_(Hsysv) { create(); }

    BaseSequencePropagator(const Sys& Hsysv, double intdtv, double originv, const Sequence& seqv, double periodv,double tolv, int verbosev, int flagsv)
      : BaseSequencePropagator_(Hsysv,intdtv,originv,BaseList<Sequence>(1,const_cast<Sequence*>(&seqv)),periodv,tolv,verbosev,flagsv),
	Hsys_(Hsysv) { create(); }

       BaseMetaPropagator* clone() const { return new BaseSequencePropagator<Sys>(*this); }    

    void genpropagator(BlockedMatrix<complex>& U, double t1, double t2, int which) const { 
      if (Hcur.empty()) {
	if (!Hsys_)
	  U.clear();
	else
	  common_propagator(U,Hsys_,t1,t2,intdt_,which);
      }
      else {
	if (!Hsys_)
	  ::libcmatrix::propagator(U,Hcur,t2-t1);
	else
	  Hadded_propagator(U,Hsys_,Hcur,t1,t2,intdt_,which);
      }
    }

    bool Hisconstant() const { return Ham_traits<Sys>::isconstant; }
    double Hperiod() const { return Hperiod_(Bool2Type<Ham_traits<Sys>::isconstant>()); }

 private:
    void create() {
      if (!Ham_traits<Sys>::isconstant && (intdt_==0.0))
	throw InvalidParameter("SequencePropagator: integration timestep unspecified for time-dependent Hamiltonian");
      allowquick_=(channels==1) && Ham_traits<Sys>::isconstant;
    }
    double Hperiod_(Bool2Type<true>) const { return 0.0; }
    double Hperiod_(Bool2Type<false>) const { return Hsys_.period(); }

  protected:
    const Sys& Hsys_;
  };

  template<class T> std::ostream& operator<< (std::ostream& ostr, const BaseSequencePropagator<T>& a) {
    a.print(ostr);
    return ostr;
  }

  class SequencePropagator : public MinimalPropagator {
 public:

   bool operator!() const { return false; } //always defined
   
   template<class TypeH> SequencePropagator(const TypeH& Hsysv,const Sequence& seq,double originv, double periodv =0,double tolv =1e-9,int verbosev =0, int flagsv =0)
     : MinimalPropagator(periodv),
       objp_(new BaseSequencePropagator<TypeH>(Hsysv,0.0,originv,seq,periodv,tolv,verbosev,flagsv)) {}

   template<class TypeH> SequencePropagator(const TypeH& Hsysv,const BaseList<Sequence>& seqsv, double originv, const BaseList<double>& periodsv,double tolv =1e-9, int verbosev =0, int flagsv =0)
     : MinimalPropagator(0.0),
     objp_(new BaseSequencePropagator<TypeH>(Hsysv,0.0,originv,seqsv,periodsv,tolv,verbosev,flagsv)) {}

   template<class TypeH> SequencePropagator(const TypeH& Hsysv,const BaseList<Sequence>& seqsv, double originv, double periodv,double tolv =1e-9, int verbosev =0, int flagsv =0)
     : MinimalPropagator(periodv), 
     objp_(new BaseSequencePropagator<TypeH>(Hsysv,0.0,originv,seqsv,periodv,tolv,verbosev,flagsv)) {}

   template<class TypeH> SequencePropagator(const TypeH& Hsysv,
					    double intdtv,const Sequence& seq,double originv, double periodv,
					    double tolv =1e-9,int verbosev =0, int flagsv =0)
     : MinimalPropagator(periodv), 
     objp_(new BaseSequencePropagator<TypeH>(Hsysv,intdtv,originv,seq,periodv,tolv,verbosev,flagsv)) {}

   template<class TypeH> SequencePropagator(const TypeH& Hsysv,
					    double intdtv,const BaseList<Sequence>& seqs,double originv, 
					    const BaseList<double>& periodsv,double tolv =1e-9,int verbosev =0, int flagsv =0)
     : MinimalPropagator(0.0),
     objp_(new BaseSequencePropagator<TypeH>(Hsysv,intdtv,originv,seqs,periodsv,tolv,verbosev,flagsv)) {}

   template<class TypeH> SequencePropagator(const TypeH& Hsysv,
					    double intdtv,const BaseList<Sequence>& seqs,double originv, 
					    double periodv,double tolv =1e-9,int verbosev =0, int flagsv =0)
     : MinimalPropagator(periodv),
     objp_(new BaseSequencePropagator<TypeH>(Hsysv,intdtv,originv,seqs,periodv,tolv,verbosev,flagsv)) {}

   const SpinOpGeneratorBase& generator() const { return objp_->generator(); }

    void propagators(BaseList<cmatrix> Us, double t1, double t2, size_t mzblk, size_t eigblk) const {
      objp_->propagators(Us,t1,t2,mzblk,eigblk);
    }

    void propagators(BaseList<cmatrix>, double, double) const { 
      throw Failed("Can't return propagator for Blocked Hamiltonian in simple matrix"); 
    } 
    
  void operator()(BlockedMatrix<complex>& U, double t1, double t2) const {
    objp_->operator()(U,t1,t2);
  }

  void operator()(cmatrix& U, double t1, double t2, size_t mzblk, size_t eigblk) const {
    objp_->operator()(U,t1,t2,mzblk,eigblk);
  }

  void operator()(cmatrix&, double, double) const {
    throw Failed("Can't return propagator for Blocked Hamiltonian in simple matrix");
  }

    void partitioning(const matrix_partition_set& partv) { objp_->partitioning(partv); }

    static Warning<> hintnoperiod_warning; //!< synchronisation hint present, but no unique period
    static Warning<> syncfailed_warning; //!< synchronisation hint didn't work!

    void synchronisation_hint(double v) { dynamic_cast<BaseSequencePropagator_*>(objp_.get())->synchronisation_hint(v); }
    double synchronisation_hint() const { return dynamic_cast<BaseSequencePropagator_*>(objp_.get())->synchronisation_hint(); }

  private:
   smartptr<BaseMetaPropagator> objp_;
 };
  
  template<> struct Propagator_traits<SequencePropagator> {
    static const bool allowsetH=false;
    static const bool hassynchronisation=true;
  };

} //namespace libcmatrix

#endif

