#undef LCM_SUPPRESS_VIEWS
#include "Sequence.h"
#include "lcm_basethreads.h"
#include <cmath>
#include <algorithm>
#include "PartitionedMatrix.h"
#include "timer.h"

namespace libcmatrix {

#include "lcm_accumulate.hpp"

  static const double rad_to_deg=180.0/M_PI;
  static const double deg_to_rad=M_PI/180.0;
  static const double microtol=5e-11; //tolerance for deciding whether two times are equal *through rounding error only*
  static const double synctol=1e-6; //tolerance for deciding whether ratio is close enough to integer

  void PulseGeneratorBase::common_init_()
  {
    if (verbose_) {
      if (verbose_>1) {
	std::cout << "Pulse Fx:\n" << Fx_;
	std::cout << "Pulse Fz: " << Fz_ << '\n';
      }    
      else
	std::cout << "Initialising pulse generator with " << (Fx_.iscomplex() ? "complex" : "real") << " Fx operator\n";
    }
    mzblocks_=opgen_.mzblocks();
    eigblocks_=opgen_.eigblocks();
    totblocks_=mzblocks_*eigblocks_;
    size_ = (totblocks_==1) ? opgen_.size(0,0) : 1; //default must be valid (caught later)
    transamp_=0.0;

    partitions_.clear();
    const ListList<size_t>& diagstr(opgen_.partitioned_diagonal_structure());
    if (diagstr.empty())
      return;

    const block_pattern blkspec(opgen_,operator_spec(nuclei_spec(nuc_),'+'),false); //!< get block pattern for +1 coherence using total blocking
     Matrix<bool> present;
     blkspec.makepresent(present);
    
    for (size_t majorblk=0;majorblk<diagstr.size();majorblk++) {
      const BaseList<size_t> fullsizes(diagstr(majorblk));
      partitions_.push_back(matrix_partition(present,fullsizes));
    }
  }

//   void PulseGenerator_SimpleTransients::prop(BlockedMatrix<complex>& U, const BlockedMatrix<complex>& Umult, int which) const
//   {
//     accumulate_(U,Umult,which,&UBtmp_);
//   }
	
//   void PulseGenerator_SimpleTransients::clear_transient()
//   {
//     if (verbose_)
//       std::cout << "Transient flushed\n";
//     reset();
//   }

//   void PulseGenerator_SimpleTransients::clear_transient(operator_type& U, int which)
//   {
//     if (!!(trans_state_.U)) {
//       prop(U,trans_state_.U,which);
//       if (verbose_>1) {
// 	std::cout << "Post-transient:\n" << trans_state_.U;
// 	std::cout << "Accumulated propagator\n" << U;
//       }
//     }
//     else {
//       if (verbose_>1)
// 	std::cout << "No transient active\n";
//     }
//     clear_transient();
//   }

  void PulseGeneratorBase::operator()(cmatrix& U,double time,double offset,double angle,phase_spec pspec, int flags) const
  {
    if (totblocks_>1)
      throw Failed("PulseGenerator(): only valid if unblocked");
    U.create(size_,size_);
    operator_type Utmp(U,mxflag::nondynamic);
    (*this)(Utmp,time,offset,angle,pspec,flags);
  }

  void PulseGeneratorBase::operator()(cmatrix& U, double tip, phase_spec pspec) const 
  {
    if (totblocks_>1)
      throw Failed("PulseGenerator(): only valid if unblocked");
    U.create(size_,size_);
    operator_type Utmp(U,mxflag::nondynamic);
    (*this)(Utmp,tip,pspec);
  }

  const ListList<complex>& ZshiftCache::operator()(double phase) const
  {
   lock_.acquire(); //!< double-locking unnecessary - not a disaster if entry is repeated
    store_t::iterator curp(cache_.find(phase));
    if (curp==cache_.end()) {
      curp=(cache_.insert(store_t::value_type(phase,ListList<complex>()))).first;
      ListList<complex>& dest(curp->second);
      rotatezfacs(dest,Fz_,phase);
      if (verbose_)
	std::cout << "ZshiftCache: inserting factors for phase=" << (rad_to_deg*phase) << ": " << dest << '\n';
    }
    lock_.release();
    return curp->second;
  }

  void PulseGeneratorBase::operator()(cmatrix& U,const cmatrix& H,double time,double offset,double angle,phase_spec pspec, int flags) const 
  {
    if (totblocks_>1)
      throw Failed("PulseGenerator(): only valid if unblocked");
    U.create(size_,size_);
    operator_type Utmp(U,mxflag::nondynamic);
    const operator_type Htmp(H,mxflag::nondynamic);
    (*this)(Utmp,Htmp,time,offset,angle,pspec,flags);
  }

  void PulseGeneratorBase::correct_partial(cmatrix& U, double time, double offres) const {
    if (totblocks_>1)
      throw Failed("PulseGenerator(): only valid if unblocked");
    U.create(size_,size_);
    operator_type Utmp(U,mxflag::nondynamic);
    correct_partial(Utmp,time,offres,0);
  }

  void PulseGeneratorBase::Hsystem(const operator_type& Hsysv)
  {
    if (!isfinished())
      throw Failed("Hsystem: invalid while PulseGenerator active");
    Hsys_=Hsysv;
    hassys_=true;
    validator_.next(); //bump check value
  }

  void RFEvent::ensurevalid() const
  {
    if (!pgenp_ || (isvalid() && (pgenp_->validation_key()==pgenkey_)))
      return;

    RFEvent* nonconst_this=const_cast<RFEvent*>(this);      
    nonconst_this->transforms();    
    nonconst_this->pgenkey_=pgenp_->validation_key();
    nonconst_this->isvalid_=true;
  }

//   void RFPulseEvent::propagate(operator_type& U,const operator_type& Hcur, int which) const
//   {
//     ensurevalid();
//     if (Ustart.empty()) {
//       if (isnull()) //check for null pulse
// 	return;
//       propagator(Ustart,Hcur,duration_);
//     }
//     apply_start(U,which);
//   }
  
  void RFPulseEvent::phase(phase_spec pspec)
  {
    const double phasev=pspec();
    if (phasev!=phase_) {
      phase_=phasev;
      flagdirty();
    }
  }

  void RFPulseEvent::rf(double rfv)
  {
    if (rf_!=rfv) {
      rf_=rfv;
      flagdirty();
    }
  }

  void SoftPulseBase::offset(double offsetv)
  { 
    if (offsetv!=offres_) { 
      offres_=offsetv;
      flagdirty();
    }
  }

  void RFChangeOffset::offset(double offsetv)
  { 
    if (offsetv!=offset_) { 
      offset_=offsetv;
      flagdirty();
    }
  }

  void SoftPulseBase::duration(double dt)
  {
    if (dt<0.0)
      throw InvalidParameter("SoftPulseBase: duration cannot be negative");
    duration_=nomduration_=dt;
    flagdirty();
  }

  void RFEvent::nominal_duration(double t)
  { //NB changes *nominal* duration
    if (t==nomduration_)
      return;
    if (t<0.0)
      throw InvalidParameter("nominal duration cannot be <0");
    nomduration_=t;
    flagdirty();
  }

  void RFPulseEvent::propagator_set::apply_start(operator_type& U, int which) const
  {
    //    if (!Ubegin)
    //  return;
    if (!!Ubegin)
      accumulate_(U,Ubegin,which);
//     if (which>=0)
//       accumulate_(U,Ubegin(which));
//     else
//       U&=Ubegin;
  }

  void RFPulseEvent::apply_start(operator_type& U, int which) const {
    ensurevalid();
    Uset.apply_start(U,which); 
  }
  void RFPulseEvent::apply_event(operator_type& U, operator_type& H, int which ) const { 
    ensurevalid();
    if (Uset.empty()) {
      if (isnull())
	return;
      throw InternalError("apply_event: U undefined for pulse");
    }
    Uset.apply_event(U,H,which,flags_ & RFPulseEvent::TRANSIENTS);
  }
  void RFPulseEvent::apply_end(operator_type& U, operator_type& H, int which) const { 
    Uset.apply_end(U,H,which,flags_ & RFPulseEvent::TRANSIENTS);
  }

  void RFPulseEvent::propagator_set::apply_event(operator_type& U, operator_type&, int which, bool hastrans) const
  {
    accumulate_(U,Utotal,which);
//     if (which>=0)
//       accumulate_(U,Ustart(which));
//     else
//       U&=Ustart;
    if (destroy_when_finished)
      const_cast<RFPulseEvent::propagator_set*>(this)->clear();
  }

  void RFPulseEvent::propagator_set::apply_end(operator_type& U, operator_type& H, int which, bool hastrans) const
  {
    if (!!Hrf)
      H-=Hrf;
    if (Uend.empty())
      return;
 //    if (which>=0)
//       accumulate_(U,Uend(which));
//     else
//       U&=Uend;
    accumulate_(U,Uend,which);
    if (!!Ubegin) {
      if (!hastrans || U.empty())
	throw InternalError("apply_end");
      operator_type Utmp;
      if (which>=0)
	conj_transpose_multiply(Utmp,Ubegin(size_t(which)),U);
      else
	conj_transpose_multiply(Utmp,Ubegin,U);
      U.swap(Utmp);
    }
    if (destroy_when_finished)
      const_cast<RFPulseEvent::propagator_set*>(this)->clear();
  }

//   void RFEvent::dump() const {
//     print(std::cout);
//   }

 std::ostream& operator<< (std::ostream& ostr, const Sequence& a)
 {
   const Sequence::const_iterator aend=a.end();
   Sequence::const_iterator astart=a.begin();
   while (astart!=aend) {
     ostr << (*astart);
     ++astart;
   }
   return ostr;
 }

  void Sequence::dump() const {
    std::cout << (*this);
  }

void PulseGeneratorBase::ensure_eigenbasis() const
{
  if (!V) {
    cache_lock_.acquire();
    if (!V) {
      PulseGeneratorBase& nonconst=const_cast<PulseGeneratorBase&>(*this);
      if (nonconst.Fx_.iscomplex())
	hermitian_eigensystem(nonconst.V.set_complex(),nonconst.eigs,nonconst.Fx_.get_complex());
      else
	hermitian_eigensystem(nonconst.V.set_real(),nonconst.eigs,nonconst.Fx_.get_real());
      if (verbose_>1)
	std::cout << "Diagonalised RF Hamiltonian\n" << Fx_ << '\n';
    }
    cache_lock_.release();
  }
}
  
  //ideal pulse
//   void PulseGenerator::operator() (operator_type& U,double angle,phase_spec pspec) const
//   {
//     rawpulsegen(U,angle,pspec());
//   }
  
  void PulseGeneratorBase::operator()(operator_type& U, double angle, phase_spec pspec) const
  {
    ensure_eigenbasis();
    
    const BaseList<double> eigsv=eigs.row();
#ifdef LCM_REENTRANT
    ListList<complex> useceigs;
#else
    ListList<complex>& useceigs(ceigs);
#endif
    if (useceigs.empty())
      useceigs.duplicate_structure(eigs);
    
    BaseList<complex> ceigsv(useceigs.row());
    for (size_t n=eigsv.size();n--;) 
      ceigsv(n)=expi(-angle*eigsv(n));
    //std::cout << "Eigs: " << ceigs << '\n';
    
    if (V.iscomplex())
      unitary_simtrans(U,useceigs,V.get_complex());
    else
      unitary_simtrans(U,useceigs,V.get_real());
    // std::cout << "Before ph\n" << U;
    const double phase(pspec());
    if (phase)
      rotatez_ip(U,Fz_,phase);
    
    if (verbose_)
      std::cout << "Pulse tip: " << (angle*rad_to_deg) << "  Phase: " << (phase*rad_to_deg) << '\n' << U;
  }

//   void PulseGenerator_SimpleTransients::operator()(operator_type& U, const operator_type& H, double time,double offset,double angle,phase_spec pspec) const
//   {
//     const double phase=pspec();
//     if (angle) {
//       rawpulsegen(U,H,time,offset,angle,0.0);
//       if (phase)
// 	rotatez_ip(U,Fz_,phase);
//     }
//     else
//       U.clear();

//     PulseGenerator_SimpleTransients* mutable_this(const_cast<PulseGenerator_SimpleTransients* >(this));
//     mutable_this->make_transient(U,time,offset,angle,phase);
//     //transient left
//   }

  void PulseGenerator::operator()(operator_type& U, double time,double offset,double angle,phase_spec pspec, int flags) const
{
  if (time==0.0) {
    U.clear();
    return;
  }
  const double trans= (flags & RFPulseEvent::TRANSIENTS) ? transamp_ : 0.0;
  if (!Ustash_ || (time!=lastdur_) || (offset!=lastoffset_) || (angle!=lastangle_) || (trans!=lasttrans_)) {
    cache_lock_.acquire();
    if (!Ustash_ || (time!=lastdur_) || (offset!=lastoffset_) || (angle!=lastangle_)) {
      rawpulsegen(Ustash_,Hsys_,time,offset,angle,0.0,flags);
      lastdur_=time;
      lastoffset_=offset;
      lastangle_=angle;
      lasttrans_=trans;
    }
    cache_lock_.release();
  }
  else {
    if (verbose_>1)
      std::cout << "Reusing last pulse propagator\n";
  }
  const double phase=pspec();
  U=Ustash_;
  if (phase)
    rotatez_ip(U,Fz_,phase);
}

  void PulseGenerator::transient_amplitude(double transampv)
  {
    //    if (!isfinished())
    //  throw Failed("transient_amplitude: can't modify while transients active");
    if (transamp_!=transampv) {
      transamp_=transampv; 
      validator_.next(); //stored propagators will be invalid
      Ustash_.clear(); //!< clear any stashed propagator
    }
  }

  void PulseGenerator::H_and_Us(operator_type& Ustart, operator_type& Hrf, diag_operator_type& Uend, double time,double offres,double angle,double phase, int flags) const 
{
  if (time<=0.0)
    throw InvalidParameter("H_and_Us: duration cannot be <=0");
  H_and_U(Hrf,Uend,time,offres,angle,phase,flags & RFPulseEvent::COHERENT);
  if (flags & RFPulseEvent::TRANSIENTS)
    make_transient(Ustart,time,angle,phase);
  else
    Ustart.clear();
}

// void PulseGenerator_SimpleTransients::H_and_Us(operator_type& Ustart, operator_type& Hrf, diag_operator_type& Uend,double time,double offres,double angle,double phase, bool coherent) const
// {
//   H_and_U(Hrf,Uend,time,offres,angle,phase,coherent);
//   Ustart.clear();
//   make_transient(Ustart,time,offres,angle,phase);
// }


//   void PulseGenerator_SimpleTransients::make_transient(BlockedMatrix<complex>& U,double dur,double offset, double angle, double phase) const
//   {
//     const double vRF=angle/(2000.0*M_PI*dur); //in kHz
//     const double transtip=transamp_*vRF*deg_to_rad;
    
//     static bool donewarning=false;
//     if (!donewarning && (transtip!=0.0) && (offset!=0.0)) {
//       std::cerr << "PulseGenerator_SimpleTransients: incompatible with (variable) transmitter offsets\n";
//       donewarning=true;
//     }
    
//     if (!isfinished() && (trans_state_.phase==phase) && (transtip==trans_state_.tip)) {
//       if (verbose_)
// 	std::cout << "Merging transient with previous RF\n";
//       return;
//     }

//     trans_state_.phase=phase;
//     trans_state_.tip=transtip;
    
//     if (!isfinished()) {
//       if (verbose_) {
// 	std::cout << "Accumulating post-transient\n";
// 	if (verbose_>1)
// 	  std::cout << trans_state_.U; 
//       }
//       prop(U,trans_state_.U);
//     }

//     if (transtip) {
//       PulseGeneratorBase::operator()(Unewtrans_,transtip,phase+M_PI/2);
//       prop(U,Unewtrans_);
//       if (verbose_) {
// 	std::cout << "Applying pre-transient with tip=" << (transtip*rad_to_deg) << " deg, phase=" << ((phase+M_PI/2)*rad_to_deg) << '\n';
// 	if (verbose_>1) {
// 	  std::cout << Unewtrans_;
// 	  std::cout << "Accumulated propagator:\n" << U;
// 	}
//       }
//       conj_transpose(trans_state_.U,Unewtrans_);
//     }
//     else
//       trans_state_.clear();
//   }

void PulseGeneratorBase::make_transient(BlockedMatrix<complex>& U,double dur, double angle, double phase) const
{
  const double vRF=angle/(2*M_PI*dur);
  make_transient(U,vRF,phase);
}

void PulseGeneratorBase::make_transient(BlockedMatrix<complex>& U,double vRF,double phase) const
{
  const double transtip=transamp_*vRF*0.001*deg_to_rad;
  if (verbose_>1)
    std::cout << "Creating transient of phase " << (phase*rad_to_deg) << " and tip angle " << (transtip*rad_to_deg) << '\n';
  if (transtip)
    PulseGeneratorBase::operator()(U,transtip,phase+M_PI/2);
  else
    U.clear();
}

  void PulseGenerator::operator()(operator_type& U,const operator_type& H,double time,double offset,double angle,phase_spec pspec, int flags) const
{
  if (hassystem() && !!H)
    throw Failed("PulseGenerator: ambiguous system Hamiltonian");
  rawpulsegen(U,H,time,offset,angle,pspec(),flags);
}

  void PulseGeneratorBase::rawpulsegen(operator_type& U,const operator_type& H,double time,double offset,double angle,double phase, int flags) const
{
  if (time<=0)
    throw InvalidParameter("PulseGenerator()");

  const double amp=angle/(2*M_PI*time);

  operator_type& Heff(U); //can use U to store Heff
  if (Fx_.iscomplex())
    multiply(Heff,amp,Fx_.get_complex());
  else
    multiply(Heff,amp,Fx_.get_real());

  if (phase)
    rotatez_ip(Heff,Fz_,phase);

  if (!!H)
    add_ip(Heff,H);

  if (offset)
    mla(Heff,-offset,Fz_);

  propagator(U,Heff,time);

  apply_flags(U,time,offset,angle,phase,flags);
}

void PulseGeneratorBase::apply_raw_transient(operator_type& U, const operator_type& Utrans, int which) const
{
  if (!Utrans)
    return;
  if (which<0) {
    operator_type Utmp;
    multiply(Utmp,U,Utrans);
    conj_transpose_multiply(U,Utrans,Utmp);
  }
  else {
    cmatrix Utmp;
    const cmatrix& Utranssel(Utrans(size_t(which)));;
    multiply(Utmp,U.front(),Utranssel);
    conj_transpose_multiply(U.front(),Utranssel,Utmp);
  }
}

  void PulseGeneratorBase::apply_flags(operator_type& U,double time,double offset,double angle,double phase, int flags, int which) const
{
  if (offset && !(flags & RFPulseEvent::COHERENT)) {
#ifdef LCM_REENTRANT
    ListList<complex> Utmpuse_;
#else
    ListList<complex>& Utmpuse_(Utmp_);
#endif
    propagator(Utmpuse_,Fz_,time*offset);
    
    if (which<0)
      U*=Utmpuse_;
    else
      (U.front())*=Utmpuse_(which);
  }
  if (flags & RFPulseEvent::TRANSIENTS) {
    operator_type Utrans;
    make_transient(Utrans,time,angle,phase);	
    apply_raw_transient(U,Utrans,which);
  }
}

  void PulseGeneratorBase::apply_flags(operator_type& U,double vRF,double offset,double phase, int flags, int which) const
{
  if (offset && !(flags & RFPulseEvent::COHERENT))
    throw Failed("apply_flags: can't perform off-resonance correction");

  if (flags & RFPulseEvent::TRANSIENTS) {
    operator_type Utrans;
    make_transient(Utrans,vRF,phase);	
    apply_raw_transient(U,Utrans,which);
  }
}

void PulseGeneratorBase::H_and_U(operator_type& Hrf, diag_operator_type& Uend,double time,double offres,double angle,double phase, int flags) const
{
  if ((time==0.0) || (angle==0.0)) {
    Hrf.clear();
    Uend.clear();
    return;
  }
  const double amp=angle/(2*M_PI*time);
  if (Fx_.iscomplex())
    multiply(Hrf,amp,Fx_.get_complex());
  else
    multiply(Hrf,amp,Fx_.get_real());
  rotatez_ip(Hrf,Fz_,phase);
  bool needcorr=false;
  if (offres) {
    mla(Hrf,offres,Fz_);
    needcorr=!(flags & RFPulseEvent::COHERENT);
  }
  if (needcorr)
    propagator(Uend,Fz_,time*offres);
  else
    Uend.clear();

  if (verbose_) {
    std::cout << "Created Hrf for pulse duration: ";
    prettyprint_time(time) << "  Angle: " << (angle*180.0/M_PI) << "  Phase: " << (phase*180.0/M_PI) << '\n' << Hrf;
    if (needcorr) 
      std::cout << "Off-resonance correction: " << Uend << '\n';
  }
}

inline double frac(double x) { return x-floor(x); }

Warning<> PulseGeneratorBase::partial_incompatible_warning("correct_partial not applicable with defined system Hamiltonian",&lcm_base_warning);

  void PulseGeneratorBase::correct_partial(operator_type& U,double toffres,double offres,int which) const
{
  if (hassystem())
    partial_incompatible_warning.raise();

  const double phaseev=frac(offres*toffres);
  if (verbose_) {
    std::cout << "Correcting for frame shift of " << 360.0*phaseev << " degrees from ";
    prettyprint_time(toffres,std::cout) << " spent " << (offres/1e3) << " kHz off-resonance\n";
  }
  //  raw_correct_partial(U,phasev,which);
  //void PulseGeneratorBase::raw_correct_partial(operator_type& U,double phasev,int which) const

#ifdef LCM_REENTRANT
  ListList<complex> Utmpuse_;
#else
  ListList<complex>& Utmpuse_(Utmp_);
#endif
  propagator(Utmpuse_,Fz_,phaseev);
  if (which<0)
    U&=Utmpuse_;
  else
    accumulate_(U,Utmpuse_(which));
}
  
void CWPulse::transforms()
{
  const double angle=2*M_PI*rf()*actduration_;
  Uset.create(*this,actduration_,angle);
}

void SoftPulse::transforms()
{
  Uset.create(*this,duration_,tip());
//   if (!pgenp_->hassystem())
//     pgenp_->H_and_Us(Ustart,Hrf,Uend,duration_,offres_,tip(),phase_,phase_coherent_,hastransients_);
//   else {
//     Uend.clear();
//     Hrf.clear();
//     (*pgenp_)(Ustart,duration_,offres_,tip(),phase_,hastransients_);
//   }
}

  void SoftPulseBase::correct_partial(operator_type& U,double toffres,int which) const
{
  if (flags_ & RFPulseEvent::COHERENT)
    throw InternalError("correct_partial shouldn't be called for phase coherent events");
  pgenp_->correct_partial(U,toffres,offres_,which);
}

void HardPulse::print(std::ostream& ostr) const
{
  ostr << "Tip: " << (tip()/deg_to_rad) << ", phase: " << (phase()/deg_to_rad) << '\n';
  if (pgenp_->verbose()>1)
    Uset.print(ostr);
}

  void RFChangeOffset::print(std::ostream& ostr) const
  {
    ostr << "Delta offset: " << (offset_/1e3) << " kHz";
    if (!Hrf.empty() && (pgenp_->verbose()>1))
      ostr << "  Delta Hrf:\n" << Hrf;
    else
      ostr << '\n';
  }
   
//Default sychronisation method
void RFEvent::synchronise(double& start, double& end, double time,char sync_, double, const UnaryFunction<double,double>& syncf)
{
  if (duration_) {
    switch (sync_) {
    case '|':
      start=syncf(time-duration_/2);
      end=syncf(time+duration_/2);
      break;
    case '+':
      start=syncf(time);
      end=syncf(time+duration_);
      break;
    case '-':
      end=syncf(time);
      start=syncf(time-duration_);
      break;
    default:
      throw Failed("synchronise: unknown sync type");
    }
    if (start>=end)
      throw Failed("synchronise: finite duration event rounded to a duration <0!");
    nominal_duration(end-start);
  }
  else 
    start=end=syncf(time);
}

  void RFPulseEvent::propagator_set::print(std::ostream& ostr) const
  {
    if (!Utotal.empty())
      ostr << "U overall:\n" << Utotal;
    if (!Ubegin.empty())
      ostr << "U initial:\n" << Ubegin;
    if (!Hrf.empty())
      ostr << "Hrf\n" << Hrf << '\n';
    if (!Uend.empty())
      ostr << "U finish (offset correction): " << Uend << '\n';
  }

void SoftPulse::print(std::ostream& ostr) const
{
  ostr << "Tip: " <<(tip()/deg_to_rad) << "  duration: ";
  prettyprint_time(duration_,ostr);
  if (nomduration_!=duration_) {
    ostr << "[nominal: ";
    prettyprint_time(nomduration_,ostr) << ']';
  }
  ostr << "  phase: " << (phase()/deg_to_rad);
  if (offset()) {
    ostr << "  off-resonance: " << (offset()*1e-3) << " kHz";
    if (flags_ & COHERENT)
      ostr << " (coherent)";
  }
  if (flags_ & TRANSIENTS)
    ostr << "  +transients";
  ostr << '\n';
  if (pulse_generator().verbose()>1)
    Uset.print(ostr);//, pulse_generator().hassystem());
}

void CWPulse::synchronise(double& start,double& end,double time,char, double scale, const UnaryFunction<double,double>& syncf)
{
  start=syncf(time);
  end=syncf(time+duration_*scale);
  if (start>end)
    throw Failed("synchronise: finite duration event rounded to a duration <0!");
  const double newdur=end-start;
  if (newdur!=actduration_) {
    actduration_=newdur;
    flagdirty();
  }
}

void CWPulse::print(std::ostream& ostr) const
{
  ostr << "vRF: " << (rf()/1e3) << " kHz  duration: ";
  prettyprint_time(duration_,ostr);
  if (fabs(duration_-actduration_)>microtol) {
    ostr << " [actual: ";
    prettyprint_time(actduration_,ostr) << ']';
  }
  ostr << "  phase: " << (phase()/deg_to_rad);
  if (offres_) {
    ostr << "  off-resonance: " << (offres_*1e-3) << " kHz";
    if (flags_ & COHERENT)
      ostr << " (coherent)";
  }
  if (flags_ & TRANSIENTS)
    ostr << "  +transients";
  ostr << '\n';
}

void TimedEvent::dump() const
{
  std::cout << (*this);
}

std::ostream& operator<< (std::ostream& ostr,const TimedEvent& a)
{
  if (a.start_<=a.end_) {
    if (a.start_==a.end_) {
      ostr << "At: ";
      prettyprint_time(a.start_,ostr);
    }
    else {
      ostr << "Start: ";
      prettyprint_time(a.start_,ostr) << "  End: ";
      prettyprint_time(a.end_,ostr);
    }
  }
  else {
    ostr << "Time: ";
    prettyprint_time(a.time_,ostr);
  }
  ostr << ' ' << a.sync_ << "  ";
  if (!a.eventp)
    return ostr << "<Empty>\n";
  else
    return ostr << *(a.eventp);
}
  
BaseWarning Sequence::sync_warning(lcm_base_warning);

void Sequence::autosync(double mintime,double maxtime)
{
  {
    static Warning<> lsync_warning("start synchronised event at end of sequence",&sync_warning);
    const reverse_iterator send(rend());
    reverse_iterator start(rbegin());
    while ((send!=start) && (start->end()>maxtime)) {
      switch (start->sync_) {
      case '|':
	start->sync_='-';
	break;
      case '+':
	if (start->start()>=mintime)
	  lsync_warning.raise();
	break;
      }
      ++start;
    }
  }
  { //do + second (for zero duration sequences)
    static Warning<> lsync_warning("end synchronised event at start of sequence",&sync_warning);
    const iterator send(end());
    iterator start(begin());
    while ((start!=send) && (start->start()<mintime)) {
      switch (start->sync_) {
      case '|':
	start->sync_='+';
	break;
      case '-':
	if (start->end()<=maxtime)
	  lsync_warning.raise();
	break;
      }
      ++start;
    }
  }
}
 
void Sequence::init(const Spectrometer& spectro)
{
  roundfunc_.resolution(spectro.time_resolution(),0.0,1e-9);
}

void Sequence::push_back(const BaseList<RFEvent*>& evlist,double timev,char sync)
{
  const size_t n=evlist.length();
  for (size_t i=0;i<n;i++) {
    RFEvent* const evp=evlist(i);
    push_back(evp,timev,sync);
    timev+=evp->duration();
  }
}

double Sequence::duration() const
{
  Sequence* mutable_this=const_cast<Sequence*>(this);
  mutable_this->applyscale(); //need to refresh in case duration of events has changed
  return lastend_;
}

TimedEvent& Sequence::push_back(RFEvent* evp, double timev,char sync) 
{
  ::std::list<TimedEvent>::push_back(TimedEvent(evp,timev,sync));
  TimedEvent& ref(back());
  syncadd(ref);
  return ref;
}

Warning<> Sequence::truncated_transient_warning("finite duration transient being truncated by overly short delay",&lcm_base_warning);

double Sequence::push_delay(double start, double dur)
{
  if (dur<0.0)
    throw InvalidParameter("Sequence::push_delay");
  const double actstart(roundfunc_(start));
  if (pgenp_) {
    check_flush(actstart);
    if (pgenp_->transient_duration()>dur)
      truncated_transient_warning.raise();
    pgenp_->clear_transient();
  }    
  lastend_=roundfunc_(start+dur);
  return lastend_-actstart;
}

void Sequence::check_flush(double start)
{
  if (fabs(start-lastend_)>1e-8)
    throw Failed("Sequence: elements must be contiguous if transients are active (or missing sequence reset?)");
}

void Sequence::syncadd(TimedEvent& ref)
{
  RFEvent* eventp=ref.eventp;
  eventp->synchronise(ref.start_,ref.end_,ref.time_,ref.sync_,scale_,roundfunc_);

  if (!havepgen_) {
    if (eventp->haspulse_generator()) {
      pgenp_=&(eventp->pulse_generator());
      if (!(pgenp_->hastransients()))
	pgenp_=NULL;
      else
	pgenp_->reset();
    }
    else
      pgenp_=NULL;
    havepgen_=true;
  }
  else {
    if (pgenp_) {
      if (pgenp_!=&(eventp->pulse_generator()))
	throw Failed("Sequence:: can't switch PulseGenerators if transients active");
      check_flush(ref.start_);
    }
  }
//   if (pgenp_)
//     eventp->transforms(); //always recompute
//   else
//     eventp->ensurevalid();
  
  //  ref.synchronise(scale_);
  if (ref.start_<lastend_)
    throw Failed("Sequence: event overlaps with previous");
  lastend_=ref.end_;
}

void Sequence::period(double periodv)
{
  if (periodv<=0.0)
    throw Failed("synchronise: bad time scaling factor");
  period_=periodv;
  applyscale();
}

void Sequence::applyscale()
{
  scale_ = period_ ? period_ : 1.0;

  const iterator iend=end();
  iterator istart=begin();

  reset();

  while (istart!=iend) {
    syncadd(*istart);
    ++istart;
  }
}

  void RFPulseEvent::propagator_set::create(const HardPulse& ev)
  {
    clear();
    ev.pulse_generator()(Utotal,ev.tip(),ev.phase_);
    check_cache(ev);
  }

  void RFPulseEvent::propagator_set::create(const SoftPulseBase& ev, double actdur, double angle)
  {
    if (actdur==0.0) {
      clear();
      return;
    }
    const PulseGeneratorBase& pgen(ev.pulse_generator());
    if (pgen.hassystem()) 
      pgen(Utotal,actdur,ev.offset(),angle,ev.phase(),ev.flags());    
    else
      Utotal.clear();
    pgen.H_and_Us(Ubegin,Hrf,Uend,actdur,ev.offset(),angle,ev.phase(),ev.flags());
    check_cache(ev);
  }
  
  void RFPulseEvent::propagator_set::check_cache(const RFPulseEvent& ev)
  {
    cachep=const_cast<simple_counter*>(ev.pulse_generator().cache_counter());
    if (cachep) {
      (*cachep)+=used();
      destroy_when_finished=cachep->full();
    }
    else
      destroy_when_finished=false;
  }

  size_t RFPulseEvent::propagator_set::used() const
  {
    //    return simple_counter::needed(Ubegin)+simple_counter::needed(Utotal)+simple_counter::needed(Hrf);
    return simple_counter_needed(Ubegin)+simple_counter_needed(Utotal)+simple_counter_needed(Hrf);
  }

  void RFPulseEvent::propagator_set::clear()
  {
    if (cachep)
      (*cachep)-=used();
    Ubegin.clear();
    Utotal.clear();
    Hrf.clear();
  }
      
  class isafter : public std::unary_function<TimedEvent,bool> {
    const double t_;
  public:
    isafter(double tv) : t_(tv) {};
    bool operator() (const TimedEvent& x) const {
      return (x.start_>=t_);
    }
  };

  void BaseSequencePropagator_::synchronisation_hint(double val)
  {
    if (val<0.0)
      throw InvalidParameter("Sequence: synchronisation hint must be >=0");
    synchint_=val;
  }

void BaseSequencePropagator_::reset(size_t i,double start)
{
  const double shifted=start-origin_;
  const bool checkadvance=(shifted>0);
  
  Sequence& cseq(seqs_(i));
  PulseGeneratorBase* pgenp=cseq.pulse_generator();
  if (pgenp)
    pgenp->reset(); //ensure pulse generator is clean
  sequence_state& cstate(status(i));
  const double per=period(i);
  const double cycletime=per && checkadvance ? origin_+floor(shifted/per)*per : origin_;
  cstate.reset(cseq,cycletime);
  
  if (checkadvance) {
    const double startmod=start-cstate.cycletime;
    cstate.iter=::std::find_if(cseq.begin(),cstate.iterend,isafter(startmod));
    if (cstate.iter!=cseq.begin()) {
      //      cstate.state=sequence_state::FINISHED;
      --(cstate.iter);
      if (cstate.iter->inevent(startmod))
	process(i);
    }
  }
}

  bool BaseSequencePropagator_::isfinished() const
  {
    for (size_t i=status.size();i--;) {
      if (!status(i).isfinished())
	return false;
    }
    return true;
  }

void BaseSequencePropagator_::create(double tolv)
{
  if (tolv<0.0 || (tolv>0.5))
    throw InvalidParameter("SequencePropagator: ridiculous tolerance value");
  tol=tolv;
  channels=seqs_.length();
  if (!hasuniqueperiod() && (channels!=periods_.length()))
    throw Mismatch("SequencePropagator()");
  if (channels==0)
    throw Failed("SequencePropagator: no active channels!");
  status.create(channels);
  for (size_t i=channels;i--;) {
    if (period(i)<0.0)
      throw InvalidParameter("SequencePropagator: period cannot be negative");
  }
  started=false;
}

  void sequence_state::reset(Sequence& cseq, double cycletimev)
  {
    pgenp_=cseq.pulse_generator();
    cycletime=cycletimev;
    if (pgenp_)
      pgenp_->reset();
    iterend=cseq.end();
    iter=cseq.begin();
    state = (iter==iterend) ? FINISHED : INACTIVE;
  }

  void sequence_state::clear_transient(BlockedMatrix<complex>& U, int which)
  {
    if (pgenp_->transient_duration())
      throw Failed("Can't handle finite duration transients yet");
    pgenp_->clear_transient(U,which);
#ifndef NDEBUG
    if (!pgenp_->isfinished())
      throw InternalError("Failed to clear transient");
#endif
    state= (iter==iterend) ? FINISHED : INACTIVE;
  }

void sequence_state::print(std::ostream& ostr) const
{
  switch (state) {
  case FINISHED: ostr << "Terminated"; break;
  case ACTIVE: ostr << "On"; break;
  case INACTIVE: ostr << "Resting"; break;
  case TRANSIENT: ostr << "Transient"; break;
  default: throw InternalError("Unknown channel state");
  }
  ostr << "  time base: ";
  prettyprint_time(cycletime,ostr) << '\n';
}

void BaseSequencePropagator_::print(std::ostream& ostr) const
{
  ostr << "Time origin: ";
  prettyprint_time(origin_,ostr) << '\n';
  ostr << "Synchronisation hint: ";
  if (synchint_)
    ostr << "none\n";
  else
    ostr << (synchint_*1e6) << " us\n";

  if (started) {
    ostr << "Last time: ";
    prettyprint_time(lasttime_,ostr) << '\n';
    ostr << "Off resonance pulses: " << offreslev << '\n';
    for (size_t i=0;i<channels;i++) {
      ostr << i << ": ";
      status(i).print(ostr);
    }
  }
  else
    ostr << "Not active\n";
}

void BaseSequencePropagator_::process(size_t which)
{
  sequence_state& cstatus(status(which));
  const RFEvent& ev=*(cstatus.iter->eventp);
  if (ev.duration()==0) 
    return;
  cstatus.state=sequence_state::ACTIVE;
  if (ev.needs_correct()) 
    offreslev++;
  ev.add_Hrf(Hcur);
//   if (ev.hassystem()) {
//     if (allowquick_)
//       Hcur=ev.H();
//     else
//       throw Failed("SequencePropagator::process is confused");
//   }
//   else {
//     if (!!ev.H())
//       add_ip(Hcur,ev.H());
//   }
}

void BaseSequencePropagator_::advance(size_t which)
{
  sequence_state& cstate(status(which));
  if (++cstate.iter==cstate.iterend) {
    if (period(which)) {
      cstate.iter=seqs_(which).begin(); //reset iterator for periodic sequence
      if (verbose_)
	std::cout << "Channel " << which << ": loop\n";
      cstate.cycletime+=period(which);
    }
    else {
      cstate.terminate();
      if (verbose_)
	std::cout << "Channel " << which << ": finish\n";
    }
  }
}

void BaseSequencePropagator_::reset(double start)
{
  Hcur.clear();
  offreslev=0;
  if (origin_>start+tol) {
    std::cout << "Warning: SequencePropagator reset to ";
    prettyprint_time(start) << " which is before sequence start (";
    prettyprint_time(origin_) << ")\n";
  }
  for (size_t i=channels;i--;)
    reset(i,start);
  lasttime_=start;
  started=true;
}

namespace {
  int check_sync(double nfloat,double tol =synctol)
  {
    if (nfloat<=0.0)
      throw InvalidParameter("check_sync");
    int n=int(nfloat+0.5);
    return (fabs(nfloat-n)>tol) ? 0 : n;
  }
}

  BaseWarning lcm_sequence_warning(lcm_base_warning);
  Warning<> SequencePropagator::hintnoperiod_warning("Synchronisation hint provided, but no unique sequence period",&lcm_sequence_warning);
  Warning<> SequencePropagator::syncfailed_warning("Synchronisation hint was not multiple of RF and Hamiltonian periods!",&lcm_sequence_warning);

/* Strictly speaking function is non-const, but state changes are not "permanent"
   It is definitely non-reentrant however! */
void BaseSequencePropagator_::calc_U(BlockedMatrix<complex>& U,double t1,double t2, int which) const
{  
  if (!started || (fabs(t1-lasttime_)>tol)) {
    BaseSequencePropagator_* mutable_this=const_cast<BaseSequencePropagator_*>(this);
    mutable_this->reset(t1);
  }
  
  //! can we accumulate multiple propagators?
  if (!(flags_ & BaseMetaPropagator::nosynchronisation)) {
    double per=0.0;
    if (hasuniqueperiod()) {
      const double Hper=Hperiod();
      if (synchint_) {
	if (((Hper==0.0) || check_sync(synchint_/Hper)) && check_sync(synchint_/period()))
	  per=synchint_;
	else
	  SequencePropagator::syncfailed_warning.raise();
      }
      else {
	per=period();
	if (Hper) {
	  int n;
	  if (Hper>per) {
	    n=check_sync(Hper/per); //!< corrected bad problem here 17/8/
	    per=Hper;
	  }
	  else
	    n=check_sync(per/Hper);
	  if (n==0)
	    per=0.0; //!< no sync;
	}
      }
    }
    else {
      if (synchint_)
	SequencePropagator::hintnoperiod_warning.raise();
    }
    if (verbose_) {
      std::cout << "Characteristic period: ";
      if (per)
	prettyprint_time(per) << '\n';
      else
	std::cout << "none\n";
    }
    if (per && (t2-t1+tol>2.0*per)) {
      //    if (Hisconstant() && hasuniqueperiod() && per && (t2-t1+tol>2.0*per)) {    
      const int reptimes=int((t2-t1+tol)/per);
      BlockedMatrix<complex> Utmp2; //!< must use different temporary
      prop_U_(Utmp2,t1,t1+per,which);
      pow(U,Utmp2,reptimes);
      if (verbose_) {
	std::cout << "Accumulating propagator " << reptimes << " times\n";
	if (verbose_>1)
	  std::cout << U;
      }
      t1+=reptimes*per;
      if (fabs(t2-t1)<tol) //NB Check goes here NOT in prop_U_ otherwise zero duration sequences get screwed
	return;
    }
    else
      U.clear();
  }
  else {
    U.clear();
    if (verbose_)
      std::cout << "Characteristic period: <ignored> (synchronisation disabled)\n";
  }
  prop_U_(U,t1,t2,which);
}

Warning<> Sequence::boundary_event_warning("Sequence: event exactly on boundary!",&lcm_base_warning);

void BaseSequencePropagator_::prop_U_(BlockedMatrix<complex>& U, double t1, double t2, int which) const
{
  if (t2<t1)
    throw InvalidParameter("SequencePropagator: start time is after end time!");

  BaseSequencePropagator_* mutable_this=const_cast<BaseSequencePropagator_*>(this);
  // static bool donewarn=false;

  for (;;) {
    double nextev=t2;
    int whichseq=-1;
    double evtime;
    
    for (size_t i=channels;i--;) {
      sequence_state& cstatus(mutable_this->status(i));
      int& cstate=cstatus.state;

      if (cstate==sequence_state::FINISHED)
	continue;

      if (cstatus.iter==cstatus.iterend) { //terminated except for transients?
	//if we haven't found a new event, check for uncleared transients
	if (cstatus.evend<t2-tol) {
	  if (verbose_)
	    std::cout << "Clearing transient on channel " << i << '\n';
	  cstatus.clear_transient(U,which);
	  if (cstate==sequence_state::FINISHED)
	    continue;
	}
	else
	  continue; //ignore transient at end of evaluation period (belongs in next)
      }

      TimedEvent& ev=*(cstatus.iter);

      if ((cstate==sequence_state::TRANSIENT) && (cstatus.cycletime+ev.start_>lasttime_+tol) && (cstatus.evend<t2-tol)) { //next event is in future
	if (verbose_)
	  std::cout << "Clearing transient on channel " << i << '\n';
	cstatus.clear_transient(U,which);
      }
	
      if (ev.start_!=ev.end_) { //finite duration events
	if (cstate==sequence_state::ACTIVE) {
	  evtime=cstatus.cycletime+ev.end_;
	  if (evtime-tol<nextev) {
	    whichseq=i;
	    nextev=evtime;
	  }
	}
	else {
	  evtime=cstatus.cycletime+ev.start_;
	  if (evtime+tol<nextev) {
	    whichseq=i;
	    nextev=evtime;
	  }
	}
      }
      else {
	evtime=cstatus.cycletime + ev.start_;
	double evtol;
	switch (ev.sync_) {
	case '+':
	  evtol=-tol;
	  break;
	case '-':
	  evtol=tol;
	  break;
	default:
	  evtol=tol/2.0;
	}
	if (evtime<nextev+evtol) {
	  whichseq=i;
	  nextev=evtime;
	  if ((ev.sync_!='-') && (evtime==t2) && (t1!=t2))
	    Sequence::boundary_event_warning.raise();
	}
      }
    }
    if (whichseq<0)
      break; //no event found

    sequence_state& cstatus(mutable_this->status(whichseq));
    int& cstate=cstatus.state;
    const TimedEvent& tev(*(cstatus.iter));
    RFEvent* const evp(cstatus.iter->eventp);
    
    if (fabs(nextev-lasttime_)>tol)
      mutable_this->delay_propagate(U,nextev,which);
    
    if (cstate==sequence_state::ACTIVE) { //finish event
      if (verbose_) {
	std::cout << "Channel " << whichseq << ": finish event at t=";
	prettyprint_time(nextev) << '\n';
      }
      evp->apply_end(U,mutable_this->Hcur,which);
      if ((channels==1) && !!Hcur && (fabs(real(Hcur.row().front()))<1e-6))
	mutable_this->Hcur.clear(); //if only channel and no offset can simply delete Hrf
      if (evp->needs_correct())
	mutable_this->offreslev--;
      cstatus.state = cstatus.isfinished() ? sequence_state::INACTIVE : sequence_state::TRANSIENT;
      mutable_this->advance(whichseq);
      mutable_this->lasttime_=nextev;
    }
    else { //start event
      if (evtime<lasttime_-tol)
	throw Failed("Event out of order!");

      if (evp->duration()) {
	if (verbose_)
	cstate=sequence_state::INACTIVE; //clear any (about to be squashed) transient
	cstatus.evend=cstatus.cycletime+tev.end_;
	if (allowquick_ && evp->pulse_generator().hassystem() && (cstatus.evend-tol<t2)) {//Can we use stored Us?
	  evp->apply_event(U,mutable_this->Hcur,which);
	  mutable_this->advance(whichseq);
	  mutable_this->lasttime_=cstatus.evend;
	  if (verbose_) {
	    std::cout << "Channel " << whichseq << ": applying finite event at t=";
	    prettyprint_time(nextev) << '\n';
	  }
	  if (!cstatus.isfinished())
	    cstate=sequence_state::TRANSIENT;
	}
	else {
	  if (verbose_) {
	    std::cout << "Channel " << whichseq << ": starting finite event at t=";
	    prettyprint_time(nextev) << '\n';
	  }
	  evp->apply_start(U,which);// check for initial transient tilt
	  mutable_this->process(whichseq);
	}
      }
      else { //zero duration
	if (verbose_) {
	  std::cout << "Channel " << whichseq << ": applying zero duration event at t=";
	  prettyprint_time(nextev) << '\n';
	}
	evp->apply_event(U,mutable_this->Hcur,which);
	mutable_this->advance(whichseq);
      }
    }
    if (verbose_>1) {
      std::cout << "Accumulated propagator at t=";
      prettyprint_time(lasttime_) << '\n' << U << '\n';
    }
  }//(infinite loop)

  if (fabs(lasttime_-t2)>tol) {
    mutable_this->delay_propagate(U,t2,which);
    //    if (verbose_>1) 
    //  std::cout << "Accumulated propagator at t=" << (t2*1e6) << " us\n" << U << '\n';
  }

  //correct if leaving while off-resonance
  int todo=offreslev;
  size_t i=0;
  while (todo) {
    sequence_state& cstatus(mutable_this->status(i));
    if (cstatus.state==sequence_state::ACTIVE) {
      const TimedEvent& tev=*(cstatus.iter);
      const RFEvent& ev=*(tev.eventp);
      if (ev.needs_correct()) {
	const double evstart=tev.start_+cstatus.cycletime;
	const double start=(evstart<=t1) ? t1 : evstart; //how long have we spent here?
	const double toffres=t2-start;
	if (toffres<0)
	  throw InternalError("This shouldn't happen!");
	if (toffres>tol) {
	  ev.correct_partial(U,toffres,which);
	  if (verbose_) {
	    std::cout << "After off-resonance correction on channel " << i << "\n";
	    if (verbose_>1) 
	      std::cout << U;
	  }
	}
	todo--;
      }
    }
    i++;
  }
  mutable_this->lasttime_=t2;
}

void BaseSequencePropagator_::delay_propagate(BlockedMatrix<complex>& U,double end, int which)
{
  if (fabs(end-lasttime_)>tol) {
    if (end<lasttime_)
      throw InternalError("end time before start time!");
    if (verbose_) {
      std::cout << "Propagating from ";
      prettyprint_time(lasttime_) << " to ";
      prettyprint_time(end) << '\n';
    }
    BlockedMatrix<complex>& Utmp(UtmpSeqs_(Type2Type<BlockedMatrix<complex> >()));
    genpropagator(Utmp,lasttime_,end,which);
    if (!!Utmp)
      U&=Utmp;
    if (verbose_>1) {
      std::cout << "Accumulated propagator at t=";
      prettyprint_time(end) << '\n' << U << '\n';
    }
  }
  lasttime_=end;
}

void BaseSequencePropagator_::operator()(cmatrix& U,double t1,double t2, size_t mzeig, size_t eig) const 
{
  BlockedMatrix<complex> Utmp;
  const size_t which=indexer_(mzeig,eig);
  calc_U(Utmp,t1,t2,which);
  Utmp.swap(U);      
}

#define LCM_BASESEQUENCE_INSTANTIATE(X) \
template void BaseSequencePropagator< X >::create();\
template void BaseSequencePropagator< X >::genpropagator(BlockedMatrix<complex>& U, double t1, double t2, int which) const;

LCM_BASESEQUENCE_INSTANTIATE( BlockedStaticHamiltonian<double> )
LCM_BASESEQUENCE_INSTANTIATE( BlockedStaticHamiltonian<complex> )
LCM_BASESEQUENCE_INSTANTIATE( BlockedSpinningHamiltonian<double> )
LCM_BASESEQUENCE_INSTANTIATE( BlockedSpinningHamiltonian<complex> )
LCM_BASESEQUENCE_INSTANTIATE( BlockedDiagonalSpinningHamiltonian )
LCM_BASESEQUENCE_INSTANTIATE( BlockedDiagonalStaticHamiltonian )

} //namespace libcmatrix
