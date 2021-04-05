#include "PhaseModulatedPropagator.h"
#include "timer.h"
#include <sstream>

namespace libcmatrix {

#include "lcm_accumulate.hpp"
  
  namespace {
    const double rad_to_deg=180.0/M_PI;

    template<typename T> int findrow(const Matrix<T>& a, const BaseList<T>& b)
    {
      for (size_t i=a.rows();i--;) {
	if (a.row(i)==b)
	  return i;
      }
      return -1;
    }
    
    std::ostream& prettyprint_interval(double t1, double t2, std::ostream& ostr =std::cout)
    {
      prettyprint_time(t1,ostr) << " to ";
      return prettyprint_time(t2,ostr);
    }

  }
    
  struct PhaseModulation::iterator {
    struct seq_state {
      Sequence::const_iterator curp_;
      Sequence::const_iterator end_;
      bool on;
      bool isfinished() const { return (curp_==end_); }
      void create(const Sequence&);
    };
    
    iterator(const PhaseModulation&);
    bool next(List<PM_state>&, List<PM_aux>&, double& dur);
    
    const BaseList<Sequence>& seqs_;
    const double tol_;
    ScratchList<seq_state> states_;
    double delay_;
    double lasttime_,period_;
  };
  
  std::ostream& operator<< (std::ostream& ostr, const PhaseModulation& a)
  {
    for (size_t j=0;j<a.states.rows();j++) {
      ostr << j << ": ";
      for (size_t i=0;i<a.states.cols();i++) {
	const PhaseModulation::PM_state& cstate(a.states(j,i));
	ostr << "Amplitude=" << (cstate.first/1e3) << " kHz";
	if (cstate.second)
	  ostr << " (offset=" << (cstate.second/1e3) << " kHz)";
	ostr << "  ";
      }
      ostr << '\n';
    }
    ostr << "CW: " << (a.isCW() ? "Yes\n" : "No\n");
    if (!a.isCW()) {
      ostr << "Number of steps: " << a.nchar << '\n';
      ostr << "Characteristic time: ";
      prettyprint_time(a.tauchar,ostr) << '\n';
      ostr << "All time steps identical: " << (a.allsame ? "Yes\n" : "No\n");
    }
    if (a.mult_) {
      ostr << "Phase modulation pattern {<phase on channel> ... [<state>] }...: ";
      for (size_t i=0;i<a.size();i++) {
	const BaseList<PM_phase_step> csteps(a(i));
	for (size_t j=0;j<csteps.size();j++)
	  ostr << (rad_to_deg*csteps(j).phase) << ' ';
	ostr << '[' << csteps.front().whichstate << "] ";
      }
      ostr << '\n';
      if (a.tauchar) {
	ostr << "Step time: ";
	prettyprint_time(a.tick()) << '\n';
      }
    }
    else
      ostr << "Step time: unset\n";
    return ostr;
  }
  
  double PhaseModulation::tick() const
  {
    if (mult_==0)
      throw Failed("PhaseModulation::tick: tick interval not set");
    return tauchar/mult_;
  }
  
  size_t PhaseModulation::roundtime(double x, double y) const
  {
    const size_t n=int(0.5+(x/y));
    return (fabs(n*y-x)>tol_) ? 0 : n;
  }
  
  PhaseModulation::iterator::iterator(const PhaseModulation& pmod)
    : seqs_(pmod.seqs),
      tol_(pmod.tol_),
      states_(seqs_.size()),
      lasttime_(0.0), period_(pmod.period())
  {
    for (size_t i=seqs_.size();i--;)
      states_(i).create(seqs_(i));

    assert(period_>0.0);
    assert(tol_>0.0);
  }
  
//   void PM_phase_step::print(std::ostream& ostr, const PhaseModulation& pm) const 
//   {
//     ostr << (rad_to_deg*phase);
//     if (pm.amp2>=0.0)
//       ostr << '(' << (isamp1 ? '1' : '2') << ')';
//   }
  
  void PhaseModulation::iterator::seq_state::create(const Sequence& seq)
  {
    curp_=seq.begin();
    end_=seq.end();
    on=false;
  }

  bool PhaseModulation::iterator::next(List<PM_state>& fstates, List<PM_aux>& aux, double& dur) 
  {
    fstates.create(seqs_.size());
    aux.create(seqs_.size());
      
    bool havetime=false;
    double nexttime=-1e6;
    
    for (size_t k=0;k<states_.size();k++) {
      seq_state& cstate(states_(k));
      if (cstate.isfinished()) {
	if (((period_-lasttime_)>tol_) && !havetime) { //!< delay at end
	  nexttime=period_;
	  havetime=true;
	}
      }
      else {
	const TimedEvent& tev=*(cstate.curp_);
	if (cstate.on || (fabs(lasttime_-tev.start())<tol_)) { //!< if at start of pulse, switch to on state
	  cstate.on=true;
	  if (!havetime || (tev.end()<nexttime)) {
	    nexttime=tev.end();
	    havetime=true;
	  }
	}
	else {	    
	  if (!havetime || (tev.start()<nexttime)) {
	    havetime=true;
	    nexttime=tev.start();
	  }
	}
      }
    }
    if (!havetime) //!< failed to find anything
      return false;
    
    dur=nexttime-lasttime_;
    if (dur<tol_)
      throw Failed("PhaseModulation: negative / zero duration encountered");
    
    for (size_t k=0;k<states_.size();k++) {
      seq_state& cstate(states_(k));
      PM_state& fstate(fstates(k));
      if (cstate.isfinished()) {
	fstate=std::pair<double,double>(0.0,0.0);
	continue;
      }
      const TimedEvent& tev=*(cstate.curp_);
      if (cstate.on) {
	if (fabs(nexttime-tev.end())<tol_) {
	  cstate.on=false;
	  ++(cstate.curp_);
	}
      }
      else {
	if (fabs(lasttime_-tev.start())>tol_) {
	  assert(tev.start()>lasttime_);
	  fstate=std::pair<double,double>(0.0,0.0);
	  continue;
	}
	cstate.on=true;
      }
      const RFPulseEvent* evp=dynamic_cast<const RFPulseEvent*>(tev.event());
      if (!evp)
	throw Failed("PhaseModulation: sequence elements can only be simple pulses");
      PM_aux& caux(aux(k));
      caux.pgenp=&(evp->pulse_generator());
      if (evp->needs_correct())
	throw InvalidParameter("PhaseModulation: pulses must be on-resonance (or phase coherent)");
      caux.flags=evp->flags();
      const double amp=evp->rf();
      if (amp<0.0)
	throw InvalidParameter("PhaseModulation: amplitude cannot be negative (change phase instead)");
      fstate=PM_state(amp,evp->offset());
      caux.phase=evp->phase();
    }
    lasttime_=nexttime;
    return true;
  }
  
  size_t PhaseModulation::multiplier() const 
  { 
    //if (steps.empty())
    //  throw Undefined("PhaseModulation::multiplier"); // taken out 18/11/09 - steps may be empty if CW?  Don't *need* to throw exception if not built
    return mult_;
  }

  void PhaseModulation::build(int mult)
  {
    if (mult<=0)
      throw InvalidParameter("PhaseModulation::build: multiple must be >=1");
    if (mult==mult_)
      return;
    
    mult_=mult;
    const size_t nsteps=mult*nchar;
    iterator iter(*this);  
    List<PM_state> fstate;
    List<PM_aux> aux;
    double dur;
    steps.create(nsteps,channels());

    size_t count=0;    
    while (iter.next(fstate,aux,dur)) {

      int whichstate=findrow(states,fstate);
      if (whichstate<0) 
	throw InternalError("PhaseModulation: sequence amplitudes/offsets have changed!");
      
      size_t n=1;
      if (!allsame) {
	n=roundtime(dur);
	if (n==0)
	  throw InternalError("PhaseModulation: sequence timing has changed!");
      }

      n*=mult;
      for (size_t k=channels();k--;) {
	const double amp=fstate(k).first;
	const double phase=aux(k).phase;
	const ListList<complex>* zshiftp=NULL;
	if (amp && phase)
	  zshiftp=&((*(channels_(k).zcachep))(phase));
	const PM_phase_step curstep(whichstate,phase,zshiftp);
	for (size_t i=n;i--;)	  
	  steps(count+i,k)=curstep;
      }
      count+=n;
    }
    assert(count==nsteps);
  }

  BaseWarning PhaseModulation::base_warning(lcm_base_warning);

  Warning<> PhaseModulation::notransients_warning("PhaseModulation is incompatible with phase transients - results are likely to be erroneous",&base_warning);

  void PhaseModulation::PM_channel::create(const PulseGeneratorBase& pgen, int flagsv, int verbose)
  {
    flags=flagsv;    
    if (pgen.hastransients())
      notransients_warning.raise();
    pgenp=&pgen;
    zcachep.reset(new ZshiftCache(pgen.Fz(),verbose));
  }

  void PhaseModulation::validate()
  {
    if (tol_<=0.0)
      throw InvalidParameter("PhaseModulation: tolerance must be >0");
    if (period_<=tol_)
      throw InvalidParameter("PhaseModulation: period must be >tolerance");
    
    tauchar=0.0;
    mult_=0;
    allsame=true;
    
    iterator iter(*this);
    List<PM_state> state;
    List<PM_aux> aux;
    channels_.create(seqs.size());
    states.create(maxn_,channels());
    states.resize(0,channels());
    double dur;
    bool lisCW=true;
    nchar=0;
    while (iter.next(state,aux,dur)) {
      nchar++;
      if (!tauchar)
	tauchar=dur;
      else {    
	lisCW=false;
	if (tauchar!=dur) {
	  allsame=false;
	  if (dur<tauchar)
	    tauchar=dur;
	}
      }      
      if (findrow(states,state)<0) {
	if (states.rows()==maxn_)
	  throw Failed("PhaseModulation: maximum of RF states exceeded");
	states.push_back_row(state);
      }
      for (size_t k=state.size();k--;) {
	if (state(k).first) {
	  PM_channel& curchan(channels_(k));
	  if (!curchan)
	    curchan.create(*(aux(k).pgenp),aux(k).flags,(verbose_>1) ? 1 : 0);
	  else {
	    if (curchan.flags!=aux(k).flags)
	      throw Failed("PhaseModulation: pulses must be of the same type");
	    if (curchan.pgenp!=aux(k).pgenp)
	      throw Failed("PhaseModulation: sequence elements must share the same pulse generator");
	  }
	}
      }
    }
    if (!lisCW)
      assert(states.rows()>0);
    
    if (!allsame) {
      //! Not all elements same length - need to check compatibility
      //! This is unlikely so don't bother that we are reworking through the sequence again
      iterator iter(*this);
      nchar=0;
      while (iter.next(state,aux,dur)) {
	const size_t n=roundtime(dur);
	if (n==0)
	  throw Failed("PhaseModulation: sequence durations are not all multiples of smallest element");
	nchar+=n;
      }
    }

    const double ttime=nchar*tauchar;
    if (lisCW)
      tauchar=0.0; 

    if ((fabs(ttime-period_)>tol_) && !lisCW)
      throw Failed("PhaseModulation: period does not match sequence duration");    
  }
 
  namespace {
    //! sensible version of % that always returns +ve
    size_t remainder(int n, size_t m)
    {
      n = n % m;
      return (n<0) ? n+m : n;
    }
  }

    void BasePMSPropagator_::doclear(size_t cleared)
    {
      assert(counterp_!=NULL);
      (*counterp_)-=cleared;
      if (verbose_>1)
	std::cout << "BasePMSPropagator: releasing " << cleared << " bytes\n";
    }

  void BasePMSPropagator_::doclaim(size_t claim) const
  {
    if (counterp_) {
      (*counterp_)+=claim;
      if (verbose_>1)
	std::cout << "BasePMSPropagator: claiming " << claim << " bytes\n";
    }
  }

  template<class T> void BasePMSPropagator_::clearcaches(List< List< BlockedMatrix<T> > >& cache)
  {
    if (cache.empty())
      return;
    if (counterp_) {
      for (size_t j=cache.size();j--;) {
	List< BlockedMatrix<T> >& curcache(cache(j));
	if (curcache.size())
	  doclear(curcache.size()*simple_counter::needed(curcache.front().row()));
      }
    }
    cache.clear(); //prevents being cleared twice
  }

  void BasePMSPropagator_::clear()
  {
    if (isstretched_) {
      clearcaches(ccaches);
      clearcaches(rcaches);
      if (counterp_ && !Ucaches.empty()) {
	for (size_t j=Ucaches.size();j--;) {
	  List< BlockedMatrix<complex> >& curcache(Ucaches(j));
	  size_t toclear=0;
	  for (size_t k=curcache.size();k--;)
	    toclear+=simple_counter_needed(curcache(k));
	  doclear(toclear);
	}
      }
    }
    else
      clearcaches(Ucaches);
  }

  void BasePMSPropagator_::calctime(size_t& ticklow_rf, double& offset, double t) const
  {
    const double ticktim = ticktime();
    const double seqt=t-origin_;
    const int i_ticklow_rf = int((seqt+tol_)/ticktim);
    offset = seqt-i_ticklow_rf*ticktim;
    ticklow_rf= remainder(i_ticklow_rf,pmseq_.size());
 }

  void BasePMSPropagator_::calctimes(size_t& ticklow_sys, size_t& ticklow_rf, double& offset, double t) const
  {
    calctime(ticklow_rf,offset,t);
    if (isstretched_)
      throw InternalError("calctimes");

    const double ticktim = ticktime();    
    const double syst=t+toffsetcur_-toffsetcache_-origin_;
    const int i_ticklow_sys = int((syst+tol_)/ticktim);
    const double offset2=syst-i_ticklow_sys*ticktim;
    if (fabs(offset-offset2)>tol_)
      throw Failed("calctimes: RF and rotation out of sync");
    ticklow_sys = remainder(i_ticklow_sys,n_);
  }

  //! clear temporary if cache is full
  void BasePMSPropagator_::clear_temporary(BlockedMatrix<complex>& a) const
  {
    if (!a.empty() && counterp_ && (counterp_->slots(simple_counter::needed(a.row()))==0))
      a.clear();
  }

  void BasePMSPropagator_::rotatez_ip(BlockedMatrix<complex>& U, const ListList<complex>& zfacs, int which)
  {
    if (which<0)
      ::libcmatrix::rotatez_ip(U,zfacs);
    else {
      if (U.size()==1)
	::libcmatrix::rotatez_ip(U.front(),zfacs(which));
      else
	::libcmatrix::rotatez_ip(U(which),zfacs(which));
    }
  }
  
void BasePMSPropagator_::calc_U_H_(BlockedMatrix<complex>& U, double t1, double t2, int which) const
  {
    const double origt1=t1;
    const double ticktim=ticktime();
    const size_t ticks_rf=pmseq_.size();
    assert(ticks_rf!=0);

    size_t tick_rf;
    double offset1;

    calctime(tick_rf,offset1,t1);

    U.clear(); 

    if (verbose_)
      std::cout << "Starting at RF tick " << tick_rf << "  Need initial fill-in: " << ((offset1>tol_) ? "Yes" : "No") << '\n';

    if (offset1>tol_) {//!< need initial fill in?
      const double endt=t1+offset1;
      explicit_propagate_H(U,t1,endt,tick_rf,which);
      if ((++tick_rf)==ticks_rf)
	tick_rf=0;
      t1=endt;
    }    
    for (;;) {
      double step=t2-t1;
      if (step<tol_)
	break;
      if (step>ticktim-tol_)
	step=ticktim;
      explicit_propagate_H(U,t1,t1+step,tick_rf,which);
      if ((++tick_rf)==ticks_rf)
	tick_rf=0;
      t1+=step;
    }
    
    if (verbose_) {
      std::cout << "Calculated overall propagator from ";
      prettyprint_interval(origt1,t2) << '\n';
      if (verbose_>1)
	std::cout << U;      
    }
    clear_temporaries();
  }

  Warning<> PhaseModulation::partitioning_change_warning("PhaseModulation: partitioning has changed between caching and point of use",&base_warning);
  Warning<> PhaseModulation::badsync_warning("PhaseModulation: start / end times not synchronised to cached times",&base_warning);    
  

  struct SpecialHelper {
    size_t nspecial;
    const BlockedMatrix<complex>* Uspecialp;
    BlockedMatrix<complex> Utmp;
    bool verbose;

    SpecialHelper(bool verbosev =false) : verbose(verbosev) { reset(); }

    void reset();
    void flush(BlockedMatrix<complex>&, int, bool&, BlockedMatrix<complex>&);
    void add(const BlockedMatrix<complex>&);
    static void doflush(BlockedMatrix<complex>&, int, bool&, const BlockedMatrix<complex>&, BlockedMatrix<complex>&);
  };

  void SpecialHelper::reset() {
    nspecial=0;
    Uspecialp=NULL;
  }

  void SpecialHelper::doflush(BlockedMatrix<complex>& U, int which, bool& isfirst, const BlockedMatrix<complex>& Uadd, BlockedMatrix<complex>& Utmp2)
  {
    if (&Uadd==&Utmp2)
      throw ArgumentClash("SpecialHelper");
    if (isfirst) {
      BaseMetaPropagator::set(U,Uadd,which);
      isfirst=false;
    }
    else
      accumulate_(U,Uadd,which,&Utmp2);
  }

  void SpecialHelper::flush(BlockedMatrix<complex>& U, int which, bool& isfirst, BlockedMatrix<complex>& Utmp2)
  {
    if (Uspecialp) {
      switch (nspecial) {
      case 0:
	throw InternalError("SpecialHelper:3");
      case 1:
	doflush(U,which,isfirst,*Uspecialp,Utmp2);
	break;
      default:
	pow(Utmp,*Uspecialp,nspecial);
	doflush(U,which,isfirst,Utmp,Utmp2);
      }
      if (verbose)
	std::cout << '^' << nspecial;
    }
    reset();
  }
    
  void SpecialHelper::add(const BlockedMatrix<complex>& U)
  {
    nspecial++;
    if (Uspecialp) {
      if (Uspecialp!=&U)
	throw InternalError("SpecialHelper");
    }
    else {
      if (nspecial!=1)
	throw InternalError("SpecialHelper:2");
      Uspecialp=&U;
    }
  }
      
void BasePMSPropagator_::calc_U(BlockedMatrix<complex>& U, double t1, double t2, int which) const
  {
   if (t2<t1)
      throw InvalidParameter("BasePMSPropagator: end before start!");

    if (partitioning()!=cached_partitionp_)
      PhaseModulation::partitioning_change_warning.raise();

    if (isstretched_) {
      calc_U_H_(U,t1,t2,which);
      return;
    }
    const bool allowcombine=(flags_ & combinepropagators);
    const double origt1=t1;
    const double ticktim=ticktime();
    const size_t ticks_rf=pmseq_.size();
    const size_t ticks_sys=n_;
    assert(ticks_sys!=0);
    assert(ticks_rf!=0);
    //! start and end times within sequence
    size_t tick_sys,tick_rf;
    size_t tick2_sys,tick2_rf;
    double offset1,offset2;

    calctimes(tick_sys,tick_rf,offset1,t1);
    calctimes(tick2_sys,tick2_rf,offset2,t2);

    if (verbose_)
      std::cout << "Starting at system tick " << tick_sys << " RF tick " << tick_rf << "  Need initial fill-in: " << ((offset1>tol_) ? "Yes" : "No") << '\n';

    bool roundwarn=false;
    bool isfirst=true;
    if (offset1>tol_) {//!< need initial fill in?
      roundwarn=true;
      double endt=t1+ticktim-offset1;
      if (endt>t2-tol_) {
	endt=t2;
	offset2=0.0;
      }
      explicit_propagator(U,t1,endt,tick_rf,which);
      if ((++tick_rf)==ticks_rf)
	tick_rf=0;
      if ((++tick_sys)==ticks_sys)
	tick_sys=0;
      isfirst=false;
      t1=endt;
    }

    const double endt=t2-offset2; //end of tick
    const double evolt=endt-t1;
    if (fabs(evolt)>tol_) { //!< non-zero tick period?
      size_t nticks=pmseq_.roundtime(evolt,ticktim);
      if (nticks==0) {
	if (verbose_)
	  std::cerr << "Evolution time of " << (evolt*1e6) << " us does not divide into tick time of " << (ticktim*1e6) << std::endl;
	throw InternalError("PhaseModulatedPropagator: evolution time is not a multiple of tick period");
      }
      if (verbose_>1)
	std::cout << "Applying cached propagators:";
      SpecialHelper shelper(verbose_>1);
      int lasttick=-1;
      for (size_t tickc=0;tickc<nticks;tickc++) {
	const BaseList<PM_phase_step> psteps(pmseq_(tick_rf));
	const size_t whichstate=psteps.front().whichstate;
	const BaseList< BlockedMatrix<complex> > cache(Ucaches(whichstate));
	if (verbose_>1) {
	  std::cout << ' ' << tick_sys;
	  if (Ucaches.size()>1)
	    std::cout << '(' << (whichstate+1) << ')';
	  //	  std::cout << (rad_to_deg*pstep.phase) << ')';
	}
	if (tick_sys>=cache.size()) {
	  if (verbose_>1)
	    std::cout << 'X'; //!< flag not cached
	  const double lt1=t1+tickc*ticktim;
	  shelper.flush(U,which,isfirst,Utmp2);
	  BlockedMatrix<complex>* Udestp=isfirst ? &U : &Utmp;
	  explicit_propagator(*Udestp,lt1,lt1+ticktim,tick_rf,which);
	  if (isfirst)
	    isfirst=false;
	  else
	    U&=Utmp;
	}
	else {
	  const BlockedMatrix<complex>* Usourcep=&(cache(tick_sys));
	  bool allowspecial=allowcombine && ((lasttick<0) || (tick_sys==lasttick));
	  if (verbose_>1)
	    std::cout << '[';
	  for (size_t j=0;j<psteps.size();j++) {
	    const PM_phase_step& pstep(psteps(j));
	    if (verbose_>1)
	      std::cout << (j ? "," : "") << (rad_to_deg*pstep.phase);
	    const ListList<complex>* zshiftp=pstep.zshiftp;
	    if (zshiftp) {
	      allowspecial=false;
	      if (Usourcep!=&Utmp) {
		Utmp=cache(tick_sys);
		Usourcep=&Utmp;
	      }
	      rotatez_ip(Utmp,*zshiftp,which); //!< a tad wasteful if more than one phase shift (unlikely)
	    }
	  }
	  if (verbose_>1)
	    std::cout << ']';
 	  if (allowspecial)
 	    shelper.add(*Usourcep);
 	  else {
 	    shelper.flush(U,which,isfirst,Utmp2);
 	    shelper.doflush(U,which,isfirst,*Usourcep,Utmp2);
// 	  if (isfirst) {
// 	    BaseMetaPropagator::set(U,*Usourcep,which);
// 	    isfirst=false;
// 	  }
// 	  else
// 	    accumulate_(U,*Usourcep,which,&Utmp2);
	  }	  
	}
#ifndef NDEBUG
	if (verbose_>1) 
	  std::cout << '\n' << U;
#endif
	lasttick=tick_sys;

	if ((++tick_rf)==ticks_rf)
	  tick_rf=0;
	if ((++tick_sys)==ticks_sys)
	  tick_sys=0;
      }
      shelper.flush(U,which,isfirst,Utmp2);
      if (verbose_>1)
	std::cout << std::endl;
    }

    if ((tick_sys!=tick2_sys) || (tick_rf!=tick2_rf)) {
      std::cerr << "Finishing at system tick " << tick_sys << " RF tick " << tick_rf << "  Expecting system tick " << tick2_sys << " RF tick " << tick2_rf << "  Need final fill-in: " << ((offset2>tol_) ? "Yes" : "No") << '\n';
      throw InternalError("Final tick state does not match expected");
    }

    if (verbose_)
      std::cout << "Finishing at system tick " << tick_sys << " RF tick " << tick_rf << "  Need final fill-in: " << ((offset2>tol_) ? "Yes" : "No") << '\n';

    if (offset2>tol_) {//!< need final fill in?
      roundwarn=true;
      BlockedMatrix<complex>* Udestp=isfirst ? &U : &Utmp;

      explicit_propagator(*Udestp,endt,t2,tick2_rf,which);
      if (isfirst)
	isfirst=false;
      else
	U&=Utmp;      
    }

    if (isfirst) //!< fallback in case no propagator created
      U.clear(); 

    if (roundwarn)
      PhaseModulation::badsync_warning.raise();

    if (verbose_>1) {
      std::cout << "Propagator from ";
      prettyprint_interval(origt1,t2) << '\n' << U;
    }
    clear_temporaries();
  }
  

// template<class T> std::ostream& operator<< (std::ostream& ostr, const BasePMSPropagator<T>& a) {
//   a.print(ostr);
//   return ostr;
// }
  
  Warning<> PhaseModulation::cache_warning("PhaseModulation: cache too small to accommodate request",&base_warning);

  template<class T> bool BasePMSPropagator_::createcache(List< BlockedMatrix<T> >& cache, size_t& n, size_t needed)
  {
    n=n_;
    bool ok=true;
    if (counterp_) {
      const size_t avail=counterp_->slots(needed);
      if (avail<n_) {
	ok=false;
	n=avail;
	if (verbose_) {
	  char buf[256];
	  snprintf(buf,sizeof(buf),".  Filling %lu of %lu items of %g K.  Would need to increase cache by %g M",(unsigned long)n,(unsigned long)n_,needed/1024.0,(n_-n)*needed/(1024.0*1024.0));
	  PhaseModulation::cache_warning.raise(buf);
	}
      }
    }
    if (n) {
      cache.create(n); 
      doclaim(n*needed); //!< NB. assumes that creation of propagators will be successful
    }
    return ok;
  }

    void BasePMSPropagator_::create(const matrix_partition_set* partpv, double toffsetv, double tolv, size_t explicitnv)
    {
      if (tolv<0.0 || (tolv>0.5))
	throw InvalidParameter("PhaseModulatedPropagator: ridiculous tolerance value");
      isstretched_=(explicitnv>0);
      n_=explicitnv;
      tol_=tolv;
      if (intdt_<=0.0)
	throw InvalidParameter("PhaseModulatedPropagator: integration timestep invalid or unspecified");   
            
      if (!pmseq_)
	throw Failed("PhaseModulatedPropagator: modulation timestep has not been set");
      Hsys_offset_(toffsetv);
      if (partpv)
	partitioning(*partpv);
      //     if (Ham_traits<Sys>::isconstant && verbose)
      //       std::cerr << "Warning: PhaseModulatedPropagator applied to time-independent Hamiltonian.  Why not use SequencePropagator?\n";
    }

  bool BasePMSPropagator_::ismultiple(double x, double y, double tol)
  {
    return (fabs(x-y*floor(0.5+x/y))<tol);
  }

  Warning<> PhaseModulation::offset_change_warning("PhaseModulation: change of time offset",&base_warning);
    
  void BasePMSPropagator_::Hsys_offset_(double toffsetv)
  {
    if (empty()) {
      toffsetcur_=toffsetcache_=toffsetv;
      return;
    }
    const double shift=toffsetv-toffsetcache_;
    //! use a higher tolerance here to avoid Failed exceptions later on from rounding errors
    if (!ismultiple(shift,dt_,0.5*tol_)) { //! ??? Previously ticktime()
      if (verbose_ && PhaseModulation::offset_change_warning.enabled()) {
	std::ostringstream str(std::ostringstream::out);
	str << " (";
	prettyprint_time(shift,str) << ") is not multiple of cached time interval (";
	prettyprint_time(dt_,str) << ").  Cached propagators will be discarded.";
	PhaseModulation::offset_change_warning.raise(str.str().c_str());
      }
      clear();
      toffsetcache_=toffsetv;
    }    
    toffsetcur_=toffsetv;
  }

  void mlaHrf(BlockedMatrix<complex>& Hrf, double amp, const PulseGeneratorBase& pgen)
  {
    if (pgen.Fx().iscomplex())
      mla(Hrf,amp,pgen.Fx().get_complex());
    else
      mla(Hrf,amp,pgen.Fx().get_real());
  }

  void mlaHrf(BlockedMatrix<double>& Hrf, double amp, const PulseGeneratorBase& pgen)
  {
    if (pgen.Fx().iscomplex())
      throw Failed("mlaHrf: complex RF Hamiltonian can't be used with real system Hamiltonian");
    mla(Hrf,amp,pgen.Fx().get_real());
  }

  template<class Sys> void BasePMSPropagator<Sys>::makeH(baseH_type& Hrf, const BaseList<PhaseModulation::PM_state>& states, const BaseList<PhaseModulation::PM_channel>& channels)
  {
    assert(channels.size()==states.size());
    for (size_t k=states.size();k--;) {
      const double amp=states(k).first;
      const PulseGeneratorBase* pgenp(channels(k).pgenp);
      if (amp) {
	if (!pgenp)
	  throw InternalError("make_H: missing PulseGenerator");
	mlaHrf(Hrf,amp,*pgenp);
	const double offset=states(k).second;
	if (offset)
	  mla(Hrf,offset,pgenp->Fz());
      }
    }      
  }

  template<class Sys> bool BasePMSPropagator<Sys>::makecache(List< BlockedMatrix<complex> >& cache, const BaseList<PhaseModulation::PM_state>& states) //, const BaseList<PhaseModulation::PM_channel>& channels)
  {
    const BaseList<PhaseModulation::PM_channel>& channels(pmseq_.pulse_channels());

    baseH_type Hrf;
    makeH(Hrf,states,channels);
    size_t n=n_;

    bool ok=true;
    for (size_t i=0;i<n;i++) {
      const double t1=i*dt_ + origin_; //+toffsetcache_;  Modified 20/1/07 - Hamiltonian should be *already adjusted* for gamma shift, toffsetcache_ records this shift
      if (!Hrf)
	common_propagator(Utmp,Hsys_,t1,t1+dt_,intdt_,-1);
      else
	Hadded_propagator(Utmp,Hsys_,Hrf,t1,t1+dt_,intdt_,-1);
      if (!Utmp)
	throw Failed("PhaseModulation: null propagator detected.  Integration time step too small?");
      if (cache.empty()) {
	ok&=createcache(cache,n,simple_counter::needed(Utmp.row()));
	if (n==0)
	  break;
      }
      cache(i).swap(Utmp);
    }
    return ok;
  }

  template<class Sys> bool BasePMSPropagator<Sys>::makeHcache(List<baseH_type>& cache, List< ListList<double> >& eigcache, const BaseList<PhaseModulation::PM_state>& states)
  {
    const BaseList<PhaseModulation::PM_channel>& channels(pmseq_.pulse_channels());

    baseH_type Hrf;
    makeH(Hrf,states,channels);
    size_t n=n_;

    bool ok=true;
    for (size_t i=0;i<n;i++) {
      const double tmid=(i+0.5)*dt_ + origin_; //+toffsetcache_;  Modified 20/1/07 - Hamiltonian should be *already adjusted* for gamma shift, toffsetcache_ records this shift
      baseH_type H(Hsys_(tmid));
      if (!!Hrf)
	H+=Hrf;
      if (cache.empty()) {
	ok&=createcache(cache,n,simple_counter::needed(H.row()));
	if (n==0)
	  break;
	eigcache.create(n);
      }
      hermitian_eigensystem(cache(i),eigcache(i),H);
    }
    return ok;
  }

  template<class Sys> void BasePMSPropagator<Sys>::makecaches()
  {
    LCM_STATIC_CHECK( !Ham_traits<Sys>::isconstant, PhaseModulatedPropagator_restricted_to_time_dependent_Hamiltonians );
    Hsys_period_=Hsys_.period();
    if (Hsys_period_<=0.0)
      throw InvalidParameter("PhaseModulatedPropagator: period of system Hamiltonian");
    
    if (!isstretched_) {
      if (pmseq_.isCW())
	n_=pmseq_.multiplier();
      else
	n_=pmseq_.roundtime(Hsys_period_,pmseq_.tick());
      if (n_==0) {
	char buffer[256];
	LCM_SNPRINTF(buffer,sizeof(buffer),"PhaseModulatedPropagator: base time interval (%g us) does not divide into Hamiltonian period (%g us)",1e6*pmseq_.tick(),1e6*Hsys_period_);
	throw Failed(buffer);
      }
    }
    dt_=Hsys_period_/n_;
    const size_t levels=pmseq_.levels();
    if (isstretched_) {
      getcaches().create(levels);
      eigcaches.create(levels);
    }
    Ucaches.create(levels);
				     
    cacheok_=true;
    toffsetcache_=toffsetcur_;
    if (verbose_) {
      std::cout << "PhaseModulation: creating " << (isstretched_ ? "eigenbasis" : "propagator") << " caches using " << n_ << " steps over a period of ";
      prettyprint_time(Hsys_period_) << " with time offset of ";
      prettyprint_time(toffsetcache_) << " and time origin at t=";
      prettyprint_time(origin_) << '\n';
    }

    for (size_t j=levels;j--;) {
      if (isstretched_) {
	List< List<baseH_type> >& curcaches(getcaches());
	cacheok_&=makeHcache(curcaches(j),eigcaches(j),pmseq_.level(j));
      }
      else
	cacheok_&=makecache(Ucaches(j),pmseq_.level(j));//,pmseq_.pulse_channels());
    }

    cached_partitionp_=partitioning();
  }
    
//   template<class Sys> template<class HType> void BasePMSPropagator<Sys>::makerawU(BlockedMatrix<complex>& U, double t1, double t2, double amp, const RCblockedmatrix& Fx, int which, Type2Type< BlockedMatrix<double> >) const
//   {
//   }

//   template<class Sys> void BasePMSPropagator<Sys>::makerawU(BlockedMatrix<complex>& U, double t1, double t2, double amp, const RCblockedmatrix& Fx, int which, Type2Type< BlockedMatrix<complex> >) const
//   {
//     BlockedMatrix<complex> Hrf(Fx.get_complex());
//     if (which<0)      
//       Hrf*=amp;
//     else
//       Hrf(which)*=amp;
//     Hadded_propagator(U,Hsys_,Hrf,t1,t2,intdt_,which);
//   }
  
  template<class Sys> void BasePMSPropagator<Sys>::explicit_propagator(BlockedMatrix<complex>& U, double t1, double t2, size_t ntick, int which) const
  {
    if (verbose_) {
      std::cout << "Explicitly calculating propagator from ";
      prettyprint_interval(t1,t2);
      if (verbose_>1)
	std::cout << " (" << which << ')';
      std::cout << '\n';
    }
    baseH_type Hrf;

    const BaseList<PM_phase_step> psteps(pmseq_(ntick));
    const size_t level=psteps.front().whichstate;
    const BaseList<PhaseModulation::PM_channel>& chans(pmseq_.pulse_channels());
    const BaseList<PhaseModulation::PM_state>& leveldesc(pmseq_.level(level));
    makeH(Hrf,leveldesc,chans);
    if (!Hrf)
      common_propagator(U,Hsys_,t1,t2,intdt_,which);
    else {
      Hadded_propagator(U,Hsys_,Hrf,t1,t2,intdt_,which);
      applyphaseshifts(U,psteps,which);
    }
    if (verbose_>1)
      std::cout << U;
  }

  void BasePMSPropagator_::applyphaseshifts(BlockedMatrix<complex>& U, const BaseList<PM_phase_step>& psteps, int which) const
    {
      const size_t level=psteps.front().whichstate;
      const BaseList<PhaseModulation::PM_state>& leveldesc(pmseq_.level(level));
      const BaseList<PhaseModulation::PM_channel>& chans(pmseq_.pulse_channels());

      if (verbose_>1)
	std::cout << '(';
      for (size_t j=0;j<psteps.size();j++) {
	const PulseGeneratorBase* pgenp(chans(j).pgenp);
	if (!pgenp)
	  continue;
	const PM_phase_step& pstep(psteps(j));
	const PhaseModulation::PM_state& clevel(leveldesc(j));
	pgenp->apply_flags(U,clevel.first,clevel.second,pstep.phase,chans(j).flags,which);
	const ListList<complex>* zshiftp=pstep.zshiftp;
	if (zshiftp)
	  rotatez_ip(U,*zshiftp,which);
	if (verbose_>1) {
	  if (j)
	    std::cout << ':';
	  std::cout << (pstep.phase*rad_to_deg);
	}
      }   
      if (verbose_>1)
	std::cout << ')';
    }   

  template<typename Sys> void BasePMSPropagator<Sys>::explicit_propagate_H(BlockedMatrix<complex>& U, double t1, double t2, size_t ntick, int which) const
  {
    const BaseList<PM_phase_step> psteps(pmseq_(ntick));
    const size_t level=psteps.front().whichstate;
    const List<baseH_type>& curcache(getcaches()(level));
    List< BlockedMatrix<complex> >& curUcache(Ucaches(level));
//     const bool iscached=(ntick<curcache.size());
//     if (!iscached) {
//       explicit_propagator(Utmp,t1,t2,ntick,which);
//       accumulate_(U,Utmp,which,&Utmp2);
//       return;
//     }

    if (verbose_) {
      std::cout << "Assembling propagator from ";
      prettyprint_interval(t1,t2);
      if (verbose_>1)
	std::cout << " (" << which << ')';
      std::cout << '\n';
    }

    const List< ListList<double> >& reigs(eigcaches(level));
    for (;;) {
      double step=t2-t1;
      if (step<tol_)
	break;
      if (step>intdt_-tol_)
	step=intdt_;
      const double t=t1+step/2.0;
      const size_t i = remainder(int(t/dt_),n_);
      if (verbose_>1) {
	std::cout << "Time mid point: ";
	prettyprint_time(t) << "  Propagator: " << i;
      }
      if (i>=curcache.size()) {
	explicit_propagator(Utmp,t1,t2,ntick,which);
	if (verbose_>1)
	  std::cout << " [explicit]\n";
	if (!U)
	  U.swap(Utmp);
	else {
	  multiply(Utmp2,Utmp,U);
	  U.swap(Utmp2);
	}
      }
      else {
	const bool trycache=(fabs(tcommon_-step)<tol_); //!< "zero" step already caught, so can't match if tcommon is zero
	if (trycache && !curUcache.empty() && !!curUcache(i)) {
	  Utmp=curUcache(i);
	  if (verbose_>1)
	    std::cout << " [from cache]";
	}
	else {
	  ::libcmatrix::propagator(ceigs,reigs(i),step);
	  unitary_simtrans(Utmp,ceigs,curcache(i));
	  if (trycache && (!counterp_ || !counterp_->full())) {
	    if (curUcache.empty())
	      curUcache.create(n_);
	    curUcache(i)=Utmp;
	    if (verbose_>1)
	      std::cout << " [cached]";
	    doclaim(simple_counter_needed(Utmp));
	  }
	}
	if (verbose_>1)
	  std::cout << '\n';
	applyphaseshifts(Utmp,psteps,which);
	accumulate_(U,Utmp,which,&Utmp2);
      }
      t1+=step;
    }
    if (verbose_>1)
      std::cout << U;
  }

#define LCM_PMSPROPAGATOR_INSTANTIATE(X)\
  template void BasePMSPropagator< X >::makecaches();\
    template void BasePMSPropagator< X >::explicit_propagate_H(BlockedMatrix<complex>& U, double t1, double t2, size_t, int) const;\
    template void BasePMSPropagator< X >::explicit_propagator(BlockedMatrix<complex>& U, double t1, double t2, size_t, int) const;

  LCM_PMSPROPAGATOR_INSTANTIATE( BlockedSpinningHamiltonian<double> )
    LCM_PMSPROPAGATOR_INSTANTIATE( BlockedSpinningHamiltonian<complex> )

}  
