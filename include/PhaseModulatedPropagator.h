#ifndef LCM_PhaseModulatedPropagator_h_
#define LCM_PhaseModulatedPropagator_h_

/*! \file
 \brief  Header file for phase modulated RF sequences
*/

#include "simple_counter.h"
#include "lcm_CommonSequence.hpp"
#include "MetaPropagation.h"
#include "Warnings.h"

namespace libcmatrix {

  class PhaseModulation;

#define LCM_DEFAULT_PHASEMODULATION_LIMIT 2

  //! structure storing step of phase modulation
  struct PM_phase_step {
    size_t whichstate;
    double phase; //!< RF phase
    const ListList<complex>* zshiftp; //!< pointer to z-rotation

    PM_phase_step() : whichstate(-1), zshiftp(NULL) {}
    PM_phase_step(size_t whichstatev, double phasev, const ListList<complex>* zshiftpv =NULL)
      : whichstate(whichstatev), phase(phasev), zshiftp(zshiftpv) {}
  };
  
  class PhaseModulation 
  {
    struct PM_aux {
      double phase;
      int flags;
      const PulseGeneratorBase* pgenp;
    };

  public:
    PhaseModulation(const BaseList<Sequence>& seqsv, double periodv, double tolv, int verbosev =0, size_t maxnv =LCM_DEFAULT_PHASEMODULATION_LIMIT)
      : seqs(seqsv), period_(periodv), tol_(tolv), verbose_(verbosev), maxn_(maxnv)
    { validate(); }

    PhaseModulation(const Sequence& seqv, double periodv, double tolv, int verbosev =0, size_t maxnv =LCM_DEFAULT_PHASEMODULATION_LIMIT)
      : seqs(BaseList<Sequence>(1,const_cast<Sequence*>(&seqv))),
	period_(periodv), tol_(tolv), verbose_(verbosev), maxn_(maxnv)
    { validate(); }
    
    struct PM_channel {
      PM_channel() : pgenp(NULL), zcachep(NULL), flags(-1) {}
      bool operator!() const { return (flags<0); }
      const PulseGeneratorBase* pgenp;  
      smartptr<ZshiftCache,false> zcachep;
      int flags;
      void create(const PulseGeneratorBase&, int, int);
    };

    typedef std::pair<double,double> PM_state;
    friend struct PM_phase_step;
    friend struct BasePMSPropagator_;
    friend std::ostream& operator<< (std::ostream&, const PhaseModulation&);
    
    void validate(); //!< check that ::Sequence is simple phase modulation
    void build(int);
    size_t multiplier() const;
    bool operator!() const { return (mult_==0); } //!< \a true if ::build has not been called
    double period() const { return period_; } //!< return overall period of PhaseModulation
    
    const BaseList<PM_phase_step> operator()(size_t n) const { return steps.row(n); }
    
    const size_t levels() const { return states.rows(); }
    const BaseList<PM_state> level(size_t n) const { return states.row(n); }
    const BaseList<PM_channel>& pulse_channels() const { return channels_; }
    size_t channels() const { return seqs.size(); }
    double timebase() const { return tauchar; }
    size_t timesteps() const { return nchar; }
    double tick() const;
    size_t roundtime(double, double) const;
    size_t roundtime(double dur) const { return roundtime(dur,tauchar); }
    bool isCW() const { return (tauchar==0.0); }
    size_t size() const { return steps.rows(); }
    
    struct iterator;

    static BaseWarning base_warning;
    static Warning<> notransients_warning;
    static Warning<> partitioning_change_warning;
    static Warning<> badsync_warning;
    static Warning<> cache_warning;
    static Warning<> offset_change_warning;

  private:
    double tauchar;
    bool allsame;
    size_t nchar;
    BaseList<Sequence> seqs;
    double period_;
    double tol_;
    int verbose_;
    size_t mult_;
    size_t maxn_;

    Matrix<PM_state> states;
    Matrix<PM_phase_step> steps;
    List<PM_channel> channels_;
  };

  std::ostream& operator<< (std::ostream&, const PhaseModulation&);

  struct BasePMSPropagator_ : public BaseMetaPropagator
  {
    template<class TypeH> BasePMSPropagator_(const TypeH& Hsysv, const matrix_partition_set* partpv, double intdtv, double originv, const PhaseModulation& pmseqv, double toffsetv, simple_counter* counterp, double tolv, int verbosev, int flagsv, size_t explicitnv, double tcommonv)
      : BaseMetaPropagator(Hsysv,verbosev,flagsv),
	intdt_(intdtv), origin_(originv),
	pmseq_(pmseqv),
	counterp_(counterp),
	tcommon_(tcommonv)
    { create(partpv,toffsetv,tolv,explicitnv); }
    
    ~BasePMSPropagator_() { clear(); }

    double origin() const { return origin_; }
    void origin(double originv) { origin_=originv; }
    double period() const { return pmseq_.period(); }
    
    void propagators(cmatrix&, size_t, double, double,size_t,size_t) const {
      throw Failed("Can't evaluate diagonal propagators for RF sequence");
    }
    
    typedef cmatrix operator_type;
    typedef RealComplexHolder< BlockedMatrix<double>, BlockedMatrix<complex> > RCblockedmatrix;
    
    void operator()(BlockedMatrix<complex>& U, double t1, double t2) const {
      calc_U(U,t1,t2,-1);
    }
    
    void operator()(cmatrix& U,double t1,double t2, size_t mzeig, size_t eig) const {
      BlockedMatrix<complex> lUtmp;
      const size_t which=indexer_(mzeig,eig);
      calc_U(lUtmp,t1,t2,which);
      lUtmp.swap(U);      
    }
    
    cmatrix operator() (double t1,double t2, size_t mzblk =0, size_t eigblk =0) const {
      cmatrix U(mxflag::temporary);
      operator()(U,t1,t2,mzblk,eigblk);
      return U;
    }
    
    //! clear propagator caches and subtract from cache count
    void clear();
    template<class T> void clearcaches(List< List< BlockedMatrix<T> > >&);
    void doclear(size_t);
    void doclaim(size_t) const;

    void clear_temporaries() const {
      clear_temporary(Utmp);
      clear_temporary(Utmp2);
    }

    void clear_temporary(BlockedMatrix<complex>&) const;

    template<class T> bool createcache(List< BlockedMatrix<T> >& cache, size_t& n, size_t needed);
    
    void calc_U(BlockedMatrix<complex>&, double t1, double t2, int which) const;
    void calc_U_H_(BlockedMatrix<complex>&, double t1, double t2, int which) const;

    virtual void explicit_propagator(BlockedMatrix<complex>& U, double t1, double t2, size_t ntick, int which) const =0;
    virtual void explicit_propagate_H(BlockedMatrix<complex>& U, double t1, double t2, size_t ntick, int which) const =0;
    virtual void Hsys_offset(double) =0;
        
    void applyphaseshifts(BlockedMatrix<complex>&, const BaseList<PM_phase_step>&, int) const;

    double intdt_;
    bool isstretched_;
    double origin_;
    const PhaseModulation& pmseq_;
    double dt_;
    double toffsetcache_; //!< time offset for which propagators have been calculated
    double toffsetcur_; //!< time offset for which propagators are being requested
    mutable simple_counter* counterp_;
    double tcommon_;
    double Hsys_period_;
    double tol_;    
    
    mutable BlockedMatrix<complex> Utmp,Utmp2;
    mutable ListList<complex> ceigs;
    size_t n_;
    mutable List< List<BlockedMatrix<complex> > > Ucaches;
    List< List<BlockedMatrix<complex> > > ccaches;
    List< List<BlockedMatrix<double> > > rcaches;    
    List< List< ListList<double> > > eigcaches;
    const matrix_partition_set* cached_partitionp_;
    bool cacheok_;

    bool cacheok() const { return cacheok_; }

    void create(const matrix_partition_set*, double gammav, double tolv, size_t explicitnv);
    bool empty() const { return Ucaches.empty(); }
    void Hsys_offset_(double);

    void calctime(size_t& ticklow_rf, double& offset, double t) const;
    void calctimes(size_t& ticklow_sys, size_t& ticklow_rf, double& offset, double t) const;

    static void rotatez_ip(BlockedMatrix<complex>&, const ListList<complex>&, int);
    static bool ismultiple(double,double,double);
    double ticktime() const { return pmseq_.isCW() ? dt_ : pmseq_.tick(); }
  };

template<class Sys> struct BasePMSPropagator : public BasePMSPropagator_ {
  typedef typename Sys::result_type baseH_type;
  
  BasePMSPropagator(const Sys& Hsysv, const matrix_partition_set* partpv, double intdtv, double originv, const PhaseModulation& pmseqv, double toffsetv, simple_counter* counterpv, double tolv, int verbosev, int flagsv, size_t explicitnv, double tcommonv)
    : BasePMSPropagator_(Hsysv,partpv,intdtv,originv,pmseqv,toffsetv,counterpv,tolv,verbosev,flagsv,explicitnv,tcommonv),
      Hsys_(Hsysv)
  { makecaches(); }

  void Hsys_offset(double toffsetv) {
    Hsys_offset_(toffsetv);
    if (empty())
      makecaches();
  }

  BaseMetaPropagator* clone() const { return new BasePMSPropagator<Sys>(*this); }      
  
private:
  bool makecache(List< BlockedMatrix<complex> >&, const BaseList<PhaseModulation::PM_state>&);//, const BaseList<PhaseModulation::PM_channel>& channels);
  bool makeHcache(List<baseH_type>&, List< ListList<double> >&, const BaseList<PhaseModulation::PM_state>&);
  
  List< List< BlockedMatrix<double> > >& getcaches(Type2Type< BlockedMatrix<double> >) { return rcaches; }
  List< List< BlockedMatrix<complex> > >& getcaches(Type2Type< BlockedMatrix<complex> >) { return ccaches; }
  const List< List< BlockedMatrix<double> > >& getcaches(Type2Type< BlockedMatrix<double> >) const { return rcaches; }
  const List< List< BlockedMatrix<complex> > >& getcaches(Type2Type< BlockedMatrix<complex> >) const { return ccaches; }
  List< List<baseH_type> >& getcaches() { return getcaches(Type2Type<baseH_type>()); }
  const List< List<baseH_type> >& getcaches() const { return getcaches(Type2Type<baseH_type>()); }

  void makecaches();
    
  static void makeH(baseH_type& Hrf, const BaseList<PhaseModulation::PM_state>& states, const BaseList<PhaseModulation::PM_channel>& channels);

  void explicit_propagator(BlockedMatrix<complex>& U, double t1, double t2, size_t ntick, int which) const;
  void explicit_propagate_H(BlockedMatrix<complex>& U, double t1, double t2, size_t ntick, int which) const; 

  protected:
    const Sys& Hsys_;
  };

  class PhaseModulatedPropagator {
  public:
    
    bool operator!() const { return false; } //always defined
    
    template<class TypeH> PhaseModulatedPropagator(const TypeH& Hsysv, const matrix_partition_set* partpv, double intdtv, double originv, const PhaseModulation& seq,double toffsetv =0.0, simple_counter* counterpv =NULL, double tolv =1e-9,int verbosev =0, int flagsv =0, size_t explicitnv =0, double tcommonv =0.0)
      : objp_(new BasePMSPropagator<TypeH>(Hsysv,partpv,intdtv,originv,seq,toffsetv,counterpv,tolv,verbosev,flagsv,explicitnv,tcommonv)) {}

    const SpinOpGeneratorBase& generator() const { return objp_->generator(); }

    void propagators(BaseList<cmatrix> Us, double t1, double t2, size_t mzblk, size_t eigblk) const {
      objp_->propagators(Us,t1,t2,mzblk,eigblk);
    }
    
    void propagators(BaseList<cmatrix> Us, double t1, double t2) const { 
      objp_->propagators(Us,t1,t2,0,0); 
    } 
    
    void operator()(BlockedMatrix<complex>& U, double t1, double t2) const {
      objp_->operator()(U,t1,t2);
    }
    
    void operator()(cmatrix& U, double t1, double t2, size_t mzblk, size_t eigblk) const {
      objp_->operator()(U,t1,t2,mzblk,eigblk);
    }
    
    void operator()(cmatrix& U, double t1, double t2) const {
      objp_->operator()(U,t1,t2,0,0);
    }

    void Hsys_offset(double t) { dynamic_cast<BasePMSPropagator_*>(objp_.get())->Hsys_offset(t); }

    void partitioning(const matrix_partition_set& partv) { objp_->partitioning(partv); }
    
    bool cacheok() const { return dynamic_cast<BasePMSPropagator_*>(objp_.get())->cacheok(); }

    double period() const { return dynamic_cast<BasePMSPropagator_*>(objp_.get())->period(); }

  private:
    smartptr<BaseMetaPropagator> objp_;
  };

} //namespace
#endif
