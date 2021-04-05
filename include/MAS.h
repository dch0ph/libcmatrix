#ifndef _MAS_H_
#define _MAS_H_

#include "space_T.h"
#include "NMR.h"
#include "MultiMatrix.h"
#include "List.h"
#include <cmath>
#include <stdlib.h>
#include "smartptr.h"
#include "PartitionedMatrix.h"

namespace libcmatrix {

extern const double MAGIC_ANGLE;
  
class RotorInfo {
  double _rotor_angle;
  double _rotor_gamma;
  bool ismagic;

  Tensor<double> _dvalues;
  cmatrix d2values_;

protected:
  void angle(double);
public:
  explicit RotorInfo(int rank_ =2,double =MAGIC_ANGLE,double =0.0);
  RotorInfo(const RotorInfo&, int);

  int rank() const { return _dvalues.rank(); }
  double angle() const { return _rotor_angle; }
  double orientation() const { return _rotor_gamma; }
  void orientation(double gamma_) { _rotor_gamma=gamma_; d2values_.clear(); }
  const Tensor<double>& dvalues() const { return _dvalues; }
  const cmatrix& d2values() const;
};

 extern const RotorInfo MASRotorInfo;

 extern double lcm_MAS_timingtol; //!< rounding error within which rotor period is assumed to be complete (typically 1 ns)
 extern Warning<> rotordrift_warning; //!< warning if timing has drifted significantly (compared to integration period)

 class IntervalSamplerBase {
 public:
   virtual ~IntervalSamplerBase() {}
   IntervalSamplerBase(bool randomisev =false) : randomise_(randomisev) {}
   virtual IntervalSamplerBase* clone() const =0;
   virtual double operator() (double t1, double t2) =0;
   virtual bool isuniform() const =0;
   typedef randomseed_t synchronise_t;
   synchronise_t get_synchronisation() const;
   virtual void set_synchronisation(synchronise_t) =0;
   void reset();
 private:
   bool randomise_;
 };
 
 class IntervalSampler : public IntervalSamplerBase {
 public:
   IntervalSampler(double jitterv =0.0, double cortimv =0.0, bool randomisev =false);
   IntervalSamplerBase* clone() const { return new IntervalSampler(*this); } //!< Note not exact clone as new object will have reset generators
   bool isuniform() const { return (jitter_==0.0); }
   double operator() (double,double);
   double jitter() const { return jitter_; }
   void jitter(double, double =0.0);
   void set_synchronisation(synchronise_t);
   //   static Warning<> excessivejitter_warning; //!< jitter is large compared to interval
   static Warning<> invalidreset_warning; //!< attempt to reset to different start time

 private:
   GaussianGenerator<> ggen_;
   double jitter_;
   double lastmidt_; //!< last value of midt (if jitter active) - bogus -ve value if unset
   double lastjitter_; //!< last jitter value
   double lastresettime_; //!< store when initially reset
   bool needsreset_; //!< flag that generator needs resetting
 };

extern IntervalSampler lcm_defsampler; //!< default sampler - can be shared as carries no state

 struct RotorSpec {
   smartptr<const RotorInfo,false> infop_;
   smartptr<const IntervalSamplerBase> samplerp_;
   double rotor_speed_;
   double rotor_phase_;
   double rotor_period_;
   double omegar_;

   explicit RotorSpec(const RotorInfo& infov) 
     : infop_(&infov,mxflag::nondynamic), 
       rotor_speed_(0.0), rotor_phase_(0.0), rotor_period_(0.0), omegar_(0.0) {}

   RotorSpec(const RotorInfo& infov, double rotor_speedv, double rotor_phasev)
     :  infop_(&infov,mxflag::nondynamic), 
	rotor_phase_(rotor_phasev) {
     rotor_speed(rotor_speedv);
   }

   RotorSpec(const RotorInfo& infov, const IntervalSamplerBase& samplerv) 
     : infop_(&infov,mxflag::nondynamic), 
       samplerp_(&samplerv,mxflag::nondynamic),
       rotor_speed_(0.0), rotor_phase_(0.0), rotor_period_(0.0), omegar_(0.0) {}

   RotorSpec(const RotorInfo& infov, const IntervalSamplerBase& samplerv, double rotor_speedv, double rotor_phasev)
     :  infop_(&infov,mxflag::nondynamic), 
	samplerp_(&samplerv,mxflag::nondynamic),
	rotor_phase_(rotor_phasev) {
     rotor_speed(rotor_speedv);
   }

   RotorSpec(const RotorSpec& spec, int rankv);

   const RotorInfo& info() const { return *infop_; }
   int rank() const { return infop_->rank(); }
   double rotor_speed() const { return rotor_speed_; }
   void rotor_speed(double);
   double period() const { return rotor_period_; }
   const IntervalSamplerBase* samplerp() const { return samplerp_.get(); }
   IntervalSamplerBase* samplerp() { return const_cast<IntervalSamplerBase*>(samplerp_.get()); }

   double rotor_phase() const { return rotor_phase_; }
   void rotor_phase(double rotor_phasev) {
     rotor_phase_=rotor_phasev;
   }
   bool isuniform() const { return (!samplerp_ || (samplerp_->isuniform())); }

   double phase(double t) const {
     if (omegar_==0.0)
       throw Failed("RotorSpec: rotor speed unset");
     return omegar_*t+rotor_phase_;
   }

   bool operator!() const { return (omegar_==0.0); }

   //!< doesn't check rotor instability
   bool operator==(const RotorSpec& a) const {
     return (infop_==a.infop_) && (rotor_speed_==a.rotor_speed_) && (rotor_phase_==a.rotor_phase_);
   }

   bool operator!=(const RotorSpec& a) const {
     return (rotor_speed_!=a.rotor_speed_) || (rotor_phase_!=a.rotor_phase_) || (infop_!=a.infop_);
   }

 };

 class SpinningHamiltonian;
 class RealSpinningHamiltonian;
 class DiagonalSpinningHamiltonian;

 class DynamicPhase : public RotorSpec {
public:
  DynamicPhase(double speed,double rotor_phase,int nobsv =0, const RotorInfo& =MASRotorInfo);
  DynamicPhase(double speed,double rotor_phase,int nobsv,const space_T& A,const RotorInfo& =MASRotorInfo);
   DynamicPhase(const RealSpinningHamiltonian&, size_t, int nobsv =0);
   DynamicPhase(const SpinningHamiltonian&, size_t, int nobsv =0);
   DynamicPhase(const DiagonalSpinningHamiltonian&, size_t, int nobsv =0);

  double component0() const { return B0; }
  complex component(int m) const;

  void rotor_phase(double lphase) { RotorSpec::rotor_phase(lphase); update(); }
  void rotor_speed(double speed) { RotorSpec::rotor_speed(speed); update(); }

   double rotor_phase() const { return rotor_phase_; }
   double rotor_speed() const { return rotor_speed_; }

  void observations(int);
  size_t observations() const { return nobs; }

  double instant_phase(double t) const;
  void tensor(const space_T&);
  
  double isotropic(int n) const {
    if (nobs)
      return B0*(n*phase_step);
    else
      throw Failed("DynamicPhase: no set number of observations");
  }
  double anisotropic(int n) const {
    if (nobs)
      return phases_table(n % phases_table.length());
    else
      throw Failed("DynamicPhase: no set number of observations"); 
  }
  double operator() (int n) const {
    return isotropic(n)+anisotropic(n);
  }

  double isotropic(double t1,double t2) const {
    return B0*MINUS_TWO_PI*(t2-t1);
  }
  double anisotropic(double,double) const; 
  double operator() (double t1,double t2) const {
    return isotropic(t1,t2)+anisotropic(t1,t2);
  }

  void propagator(BaseList<complex>, const BaseList<double>&, double iphase) const;

  List<complex> propagator(const BaseList<double>& eigs, double iphase) const {
    List<complex> U(eigs.length(),mxflag::temporary);
    propagator(U,eigs,iphase);
    return U;
  }

  void propagators(cmatrix&, const BaseList<double>&) const;

  cmatrix propagators(const BaseList<double>& eigs) const {
    cmatrix U(mxflag::temporary);
    propagators(U,eigs);
    return U;
  }

  void propagators(BaseList<complex>) const;

  List<complex> propagators() const {
    if (!nobs)
      throw Failed("DynamicPhase: no set number of observations"); 
    List<complex> U(nobs,mxflag::temporary);
    propagators(U);
    return U;
  }


  friend std::ostream& operator<< (std::ostream&,const DynamicPhase&);

 private:
  size_t nobs;
  List<complex> Bvalues;

  double phase_step;
  cmatrix factor_table;
  List<double> phases_table;
  double B0;

  void update();
  void update_tensor();
  void make_factortable();
  void dump() const;
};
   
 template<class T> struct SpinningHamiltonianBase : public UnaryFunction<T,double>, public RotorSpec {
    explicit SpinningHamiltonianBase(const RotorSpec& info_) 
      : RotorSpec(info_) {} 
    explicit SpinningHamiltonianBase(const RotorInfo& info_)
      : RotorSpec(info_) {}

    explicit SpinningHamiltonianBase(double rotor_speedv, double rotor_phasev =0.0, const RotorInfo& info_ =MASRotorInfo) 
      : RotorSpec(info_,rotor_speedv,rotor_phasev) {} 

   SpinningHamiltonianBase(double rotor_speedv, double rotor_phasev, const RotorInfo& info_, const IntervalSamplerBase& samplerv) 
     : RotorSpec(info_,samplerv,rotor_speedv,rotor_phasev) {} 

   SpinningHamiltonianBase(const RotorSpec& spec, int rankv)
     : RotorSpec(spec,rankv) {}

 };

 class SpinningHamiltonian : public SpinningHamiltonianBase<cmatrix> {
public:
  explicit SpinningHamiltonian(const RotorInfo& info_ =MASRotorInfo)
    : SpinningHamiltonianBase<cmatrix>(info_) { clear(); }

   explicit SpinningHamiltonian(double rotor_speedv, double rotor_phasev =0.0, const RotorInfo& info_ =MASRotorInfo)
     : SpinningHamiltonianBase<cmatrix>(rotor_speedv,rotor_phasev, info_) { clear(); }

     SpinningHamiltonian(double rotor_speedv, double rotor_phasev, const RotorInfo& info_, const IntervalSamplerBase& samplerv)
       : SpinningHamiltonianBase<cmatrix>(rotor_speedv,rotor_phasev, info_, samplerv) { clear(); }

   explicit SpinningHamiltonian(const RotorSpec& spec_)
     : SpinningHamiltonianBase<cmatrix>(spec_) { clear(); }

   SpinningHamiltonian(const SpinningHamiltonian&,const BaseList<size_t>&);

  friend class RealSpinningHamiltonian;
  friend class DiagonalSpinningHamiltonian;
  friend class DynamicPhase;
  
  void tensor(const Tensor<cmatrix>&);
  void tensor(const space_T&);
  void clear();

  void add(const space_T& A, const cmatrix& H) { _add(A,H,H.rows()); }
  void add(const space_T& A, const rmatrix& H) { _add(A,H,H.rows()); }
  void add(const space_T&, const BaseList<cmatrix>&);
  void add(const space_T&, const BaseList<rmatrix>&);
  void add(const space_T&, const MultiMatrix<double,3>&);
  void add(const space_T& A, const BaseList<double>& H) { _add(A,H,H.size()); }
  template<class T> void add(double coup, const T& H) { mla(component(0),coup,H); }

  template<class T2> SpinningHamiltonian& operator+= (const Matrix<T2>& H) { component(0)+=H; return *this; }
  template<class T2> SpinningHamiltonian& operator-= (const Matrix<T2>& H) { component(0)-=H; return *this; }
  SpinningHamiltonian& operator+= (const BaseList<double>& H) { component(0)+=H; return *this; }
  SpinningHamiltonian& operator-= (const BaseList<double>& H) { component(0)-=H; return *this; }
  SpinningHamiltonian& operator+= (const SpinningHamiltonian&); 
  SpinningHamiltonian& operator-= (const SpinningHamiltonian&);
  SpinningHamiltonian& operator+= (const DiagonalSpinningHamiltonian&);
  SpinningHamiltonian& operator-= (const DiagonalSpinningHamiltonian&);

  size_t size() const { return tensor_values.empty() ? 0 : tensor_values.front().rows(); }
  bool operator!() const { return tensor_values.empty(); }

  void operator() (cmatrix&,double) const;
  cmatrix operator() (double t) const {
    cmatrix tmp(mxflag::temporary);
    (*this)(tmp,t);
    return tmp;
  }

  inline cmatrix& component(int);
  inline const cmatrix& component(int) const;

   size_t usage() const;

   int rank() const { return tensor_values.empty() ? 0 : (tensor_values.size()-1)/2; } //!< reflects actual rank (overrides RotorSpec)

 private:
  List<cmatrix> tensor_values; //Use List<cmatrix> rather than MultiMatrix<complex,3> as this allows H(0) to be undefined for MAS

  template<class T> void _add(const space_T&, const T&, size_t);
  template<class M> void add_ns(const space_T&, const M&);
  bool ensure_tensor(size_t);
  void squeeze(double& H0,BaseList<complex> Hm,size_t which,double tol) const;
  void dump() const;
};

std::ostream& operator << (std::ostream&,const SpinningHamiltonian &);

 class RealSpinningHamiltonian : public SpinningHamiltonianBase<rmatrix> {
public:
  explicit RealSpinningHamiltonian(const RotorInfo& info_ =MASRotorInfo)
    : SpinningHamiltonianBase<rmatrix>(info_) { clear(); }

    explicit RealSpinningHamiltonian(const RotorSpec& spec_)
      : SpinningHamiltonianBase<rmatrix>(spec_) { clear(); }

    explicit RealSpinningHamiltonian(double rotor_speedv, double rotor_phasev =0.0, const RotorInfo& info_ =MASRotorInfo)
      : SpinningHamiltonianBase<rmatrix>(rotor_speedv,rotor_phasev,info_) { clear(); }

      RealSpinningHamiltonian(double rotor_speedv, double rotor_phasev, const RotorInfo& info_, const IntervalSamplerBase& samplerv)
	: SpinningHamiltonianBase<rmatrix>(rotor_speedv,rotor_phasev,info_,samplerv) { clear(); }

  RealSpinningHamiltonian(const RealSpinningHamiltonian&,const BaseList<size_t>&);

  friend class SpinningHamiltonian;
  friend class DiagonalSpinningHamiltonian;
  friend class DynamicPhase;
  
  void clear();

  void add(const space_T& A, const rmatrix &H) { _add(A,H,H.rows()); }
  void add(const space_T&, const BaseList<rmatrix>&);
  void add(const space_T& A, const BaseList<double>& H) { _add(A,H,H.size()); }
  void add(double coup, const BaseList<double>& H) { mla(component0(),coup,H); }
  void add(double coup, const rmatrix& H) { mla(component0(),coup,H); }

  void tensor(const space_T&);
  void tensor(const Tensor<cmatrix>&);

  RealSpinningHamiltonian& operator+= (const rmatrix& H) { H0+=H; return *this; }
  RealSpinningHamiltonian& operator-= (const rmatrix& H) { H0-=H; return *this; }
  RealSpinningHamiltonian& operator+= (const BaseList<double>& H) { H0+=H; return *this; }
  RealSpinningHamiltonian& operator-= (const BaseList<double>& H) { H0-=H; return *this; }
  RealSpinningHamiltonian& operator+= (const RealSpinningHamiltonian&);
  RealSpinningHamiltonian& operator-= (const RealSpinningHamiltonian&);
  RealSpinningHamiltonian& operator+= (const DiagonalSpinningHamiltonian&);
  RealSpinningHamiltonian& operator-= (const DiagonalSpinningHamiltonian&);

  size_t size() const { return tensor_values.dimension(1);}
  bool operator!() const { return H0.empty() && tensor_values.empty(); }
  int rank() const { return tensor_values.empty() ? 0 : tensor_values.dimension(0U); } //!< reflects actual rank (overrides RotorSpec)

  template<class T> void operator() (Matrix<T>&, double) const;
  rmatrix operator() (double t) const {
    rmatrix tmp(mxflag::temporary); 
    RealSpinningHamiltonian::operator()(tmp,t);
    return tmp;
  }

  rmatrix& component0() { return H0; }
  const rmatrix& component0() const { return H0; }
  inline cmatrix component(int m); //NB these return non-dynamic cmatrices
  inline const cmatrix component(int m) const;

  size_t usage() const { return 2*tensor_values.size()+H0.size(); }
 private:
  rmatrix H0;
  MultiMatrix<complex,3> tensor_values;

  template<class T> void _add(const space_T&, const T&, size_t);
  bool ensure_tensor(size_t);
  void dump() const;
};

std::ostream& operator << (std::ostream&,const RealSpinningHamiltonian &);

  class DiagonalSpinningHamiltonian : public SpinningHamiltonianBase< List<double> > {
public:
  friend class SpinningHamiltonian;
  friend class RealSpinningHamiltonian;
  friend class DynamicPhase;

  explicit DiagonalSpinningHamiltonian(const RotorInfo& info_ =MASRotorInfo)
    : SpinningHamiltonianBase< List<double> >(info_) { clear(); }

    explicit DiagonalSpinningHamiltonian(const RotorSpec&);

    explicit DiagonalSpinningHamiltonian(double rotor_speedv, double rotor_phasev =0.0, const RotorInfo& info_ =MASRotorInfo)
      : SpinningHamiltonianBase< List<double> >(rotor_speedv,rotor_phasev, info_) { clear(); }

     DiagonalSpinningHamiltonian(double rotor_speedv, double rotor_phasev, const RotorInfo& info_, const IntervalSamplerBase& samplerv)
       : SpinningHamiltonianBase< List<double> >(rotor_speedv,rotor_phasev, info_, samplerv)// { clear(); }
       { throw Failed("Can't currently create diagonal spinning Hamiltonian with timing jitter"); }

      DiagonalSpinningHamiltonian(const SpinningHamiltonian& Ham, const BaseList<double>& Zeeman, const BaseList<size_t>& order, double degfac =1e-8)
	: SpinningHamiltonianBase< List<double> >(Ham,4) { interactions(Ham,Zeeman,order,degfac); } //2nd order

  DiagonalSpinningHamiltonian(const SpinningHamiltonian&, size_t);
  DiagonalSpinningHamiltonian(const RealSpinningHamiltonian&, size_t);
  DiagonalSpinningHamiltonian(const DiagonalSpinningHamiltonian&, size_t);
  DiagonalSpinningHamiltonian(const DiagonalSpinningHamiltonian&, const BaseList<size_t>&);
  
  void clear();
  void interactions(const SpinningHamiltonian& Ham, const BaseList<double>& Zeeman, const BaseList<size_t>&, double);
  void interactions(const RealSpinningHamiltonian& Ham, const BaseList<double>& Zeeman, const BaseList<size_t>&, double);
  void add(const space_T&, const BaseList<double>&);
  void add(double coup, const BaseList<double>& H) { mla(component0(),coup,H); }
  void addQ2(double, const space_T&, const BaseList<double>&, const BaseList<double>&, bool =false);

  DiagonalSpinningHamiltonian& operator+= (const BaseList<double>& H) { H0+=H; return *this; }
  DiagonalSpinningHamiltonian& operator-= (const BaseList<double>& H) { H0-=H; return *this; }
  DiagonalSpinningHamiltonian& operator+= (const DiagonalSpinningHamiltonian&);
  DiagonalSpinningHamiltonian& operator-= (const DiagonalSpinningHamiltonian&);

  int rank() const { return tensor_values.rows(); } //!< reflects actual rank (overrides RotorSpec)
  size_t size() const;
  bool operator!() const { return H0.empty() && tensor_values.empty(); }

  void operator()(BaseList<double>,double) const;
  void operator()(List<double>& dest,double t) const {
    dest.create(size());
    operator()( static_cast< BaseList<double>& >(dest),t);
  }
  List<double> operator() (double t) const {
    List<double> tmp(mxflag::temporary);
    operator()(tmp,t);
    return tmp;
  }

  List<double>& component0() { return H0; }
  const List<double>& component0() const { return H0; }
  inline BaseList<complex> component(int m);
  inline const BaseList<complex> component(int m) const;

  size_t usage() const { return 2*tensor_values.size()+H0.size(); }

 private:
  List<double> H0;
  cmatrix tensor_values;

  void dump() const;
  template<class HType> void interactions_(const cmatrix[5], const HType& Ham, const BaseList<double>& Zeeman, const BaseList<size_t>&, double);
};

std::ostream& operator << (std::ostream&,const DiagonalSpinningHamiltonian&);

 void unitary_simtrans(SpinningHamiltonian&, const SpinningHamiltonian&, const cmatrix&);
 void unitary_simtrans(SpinningHamiltonian&, const RealSpinningHamiltonian&, const cmatrix&);
 void unitary_simtrans(SpinningHamiltonian&, const DiagonalSpinningHamiltonian&, const cmatrix&);

 cmatrix gammareduce(const BaseList<cmatrix>&, size_t);
 void gammareduce(BaseList<cmatrix>, const BaseList<cmatrix>&, size_t);

  //Not thread-safe : do not share between threads
 template<class M =cmatrix> class HomogeneousPropagator : public PropGen_t {
 public:
   enum mode_t { Diagonalisation, Chebyshev, LiouvilleSpace };

  template<class HType> HomogeneousPropagator(const HType& Ham_,double intdtv,const BaseList<double>& Zeemanv, int verbosev =0)
    : Hamp(&Ham_), size_(Ham_.size()),
      period(Ham_.period()), intdt(intdtv), Zeeman(Zeemanv), mode_(Diagonalisation),
      verbose_(verbosev) {}
 
   template<class HType> HomogeneousPropagator(const HType& Ham_,double intdtv, mode_t modev =Diagonalisation, const matrix_partition* partpv =NULL, int verbosev =0)
     : Hamp(&Ham_), period(Ham_.period()), intdt(intdtv), mode_(modev), partp(partpv), verbose_(verbosev) {}

   PropGen_t* clone() const { return new HomogeneousPropagator(*this); }
   size_t size() const { return size_; }
 void verbose(int verbosev) { verbose_=verbosev; }

  void operator() (cmatrix&, double,double) const;
  cmatrix operator() (double t1,double t2) const {
    cmatrix U(mxflag::temporary);
    operator()(U,t1,t2);
    return U;
  };


 private:
  const SpinningHamiltonianBase<M>* Hamp;
 mutable smartptr<IntervalSamplerBase> Samplerp; //!< cached local copy of sampler function
 size_t size_;
 double period;
 mutable M H;
  mutable cmatrix Utmp;
  const double intdt;
  const ScratchList<double> Zeeman;

   mode_t mode_;
   const matrix_partition* partp;
   int verbose_;

   void prop_U(cmatrix&, double,double) const;
 };
   
template<class M> class HomogeneousPropagator_second : public PropGen_t
{
public:
  HomogeneousPropagator_second(const SpinningHamiltonianBase<M>& Hamv, double intdtv, const ListList<size_t>& blkstrv, const BaseList<double>& Hzeemantv, int verbosev =0)
    : Hamp(&Hamv),
      intdt(intdtv), blkstr(blkstrv), Hzeemant(Hzeemantv), verbose_(verbosev) {}

  void operator()(cmatrix&, double, double) const;

  cmatrix operator()(double t1, double t2) const {
    cmatrix U;
    (*this)(U,t1,t2);
    return U;
  }
  PropGen_t* clone() const { return new HomogeneousPropagator_second(*this); }
  size_t size() const { return Hzeemant.size(); }

private:
  const SpinningHamiltonianBase<M>* Hamp;
  mutable smartptr<IntervalSamplerBase> Samplerp; //!< cached local copy of sampler function
  const double intdt;
  const ListList<size_t>& blkstr;
  const BaseList<double>& Hzeemant;
  int verbose_;

  mutable cmatrix Utmp;
  mutable M H,V;
  mutable List<double> eigs;
  mutable List<complex> ceigs;

  cmatrix* scratch_space() const { return NULL; }
};

 class DiagonalInhomogeneousPropagator : public DiagPropGen_t
 {
   const BaseList<DynamicPhase> dphases;
   const rmatrix Hs;
   smartptr<const DiagonalSpinningHamiltonian,false> Hspinp;
   const ScratchList<double> Hcon;
   int verbose_;

 public:
  DiagonalInhomogeneousPropagator(const DynamicPhase& dphasev,const BaseList<double>& Hv,const BaseList<double>& Hconv =BaseList<double>(), int verbosev =0)
    : dphases(1,const_cast<DynamicPhase*>(&dphasev)), Hs(1,Hv.length(),Hv.vector()), Hcon(Hconv), verbose_(verbosev) {}

  DiagonalInhomogeneousPropagator(const BaseList<DynamicPhase>& dphasesv,const rmatrix& Hsv, const BaseList<double>& Hconv =BaseList<double>(), int verbosev =0)
    : dphases(dphasesv), Hs(Hsv), Hcon(Hconv), verbose_(verbosev) {}

  DiagonalInhomogeneousPropagator(const DynamicPhase& dphasev,const BaseList<double>& Hv, int verbosev)
    : dphases(1,const_cast<DynamicPhase*>(&dphasev)), Hs(1,Hv.length(),Hv.vector()), verbose_(verbosev) {}

  DiagonalInhomogeneousPropagator(const BaseList<DynamicPhase>& dphasesv,const rmatrix& Hsv, int verbosev)
    : dphases(dphasesv), Hs(Hsv), verbose_(verbosev) {}

    DiagonalInhomogeneousPropagator(const DiagonalSpinningHamiltonian& Hspinv, int verbosev =0)
      : Hspinp(&Hspinv,mxflag::nondynamic), verbose_(verbosev) {}

   DiagPropGen_t* clone() const { return new DiagonalInhomogeneousPropagator(*this); }

   size_t size() const;
   size_t interactions() const { return dphases.length(); }

   void operator()(List<complex>& U, double t1, double t2) const;
   void operator()(BaseList<complex>, double,double) const;   
   List<complex> operator() (double t1,double t2) const;
 };
   
 class InhomogeneousPropagator : public PropGen_t {
  DiagonalInhomogeneousPropagator dobj;
  const cmatrix VC;
  const rmatrix VR;

 public:
  InhomogeneousPropagator(const DynamicPhase& dphasev,const BaseList<double>& Hv,const BaseList<double>& Hconv =BaseList<double>(),const rmatrix& VRv =rmatrix(), int verbosev =0)
    : dobj(dphasev,Hv,Hconv,verbosev), VR(VRv) {}

  InhomogeneousPropagator(const DynamicPhase& dphasev,const BaseList<double>& Hv,const BaseList<double>& Hconv,const cmatrix& VCv, int verbosev =0)
    : dobj(dphasev,Hv,Hconv,verbosev), VC(VCv) {}

    InhomogeneousPropagator(const BaseList<DynamicPhase>& dphasesv,const rmatrix& Hsv, const BaseList<double>& Hconv =BaseList<double>(), const rmatrix& VRv =rmatrix(), int verbosev =0)
    : dobj(dphasesv,Hsv,Hconv,verbosev), VR(VRv) {}

  InhomogeneousPropagator(const BaseList<DynamicPhase>& dphasesv,const rmatrix& Hsv, const BaseList<double>& Hconv,const cmatrix& VCv, int verbosev =0)
    : dobj(dphasesv,Hsv,Hconv,verbosev), VC(VCv) {}

    InhomogeneousPropagator(const DiagonalSpinningHamiltonian& Hspinv, const rmatrix& VRv =rmatrix(), int verbosev =0)
      : dobj(Hspinv,verbosev), VR(VRv) {}

    InhomogeneousPropagator(const DiagonalSpinningHamiltonian& Hspinv, const cmatrix& VCv, int verbosev =0)
      : dobj(Hspinv,verbosev), VC(VCv) {}

  InhomogeneousPropagator(const DynamicPhase& dphasev, const BaseList<double>& Hv, int verbosev)
    : dobj(dphasev,Hv,verbosev) {}

   InhomogeneousPropagator(const BaseList<DynamicPhase>& dphasesv,const rmatrix& Hsv, int verbosev)
     : dobj(dphasesv,Hsv,verbosev) {}

   InhomogeneousPropagator(const DiagonalSpinningHamiltonian& Hspinv, int verbosev)
      : dobj(Hspinv,verbosev) {}

  bool havetransform() const { return !(VC.empty() && VR.empty()); }
  size_t interactions() const { return dobj.interactions(); }
  
  void operator() (cmatrix&, double,double) const;
  cmatrix operator() (double t1,double t2) const {
    cmatrix U(mxflag::temporary);
    operator()(U,t1,t2);
    return U;
  };
};
   
  template<class T1,class T2> inline void mlad(Matrix<T1>& a, T1 b,const T2& H)
    {
      mla(a,b,H);
    }

  template<class T1> inline void mlad(Matrix<T1>& a, T1 b, const double& H)
    {
      b*=H;
      for (size_t m=a.rows();m--;)
	a(m,m)+=b;
    }

  template<class T> bool hasodd(const Tensor<T>& A)
  {
    int l=A.max_rank();
    if (l<0) throw Undefined("hasodd");
    if (l & 1) return true;
    for (l--;l>=1;l-=2) {
      if (A.have_rank(l))
	return true;
    }
    return false;
  }

  extern const char NOODD[];

const cmatrix RealSpinningHamiltonian::component(int m) const
{ 
#ifndef NDEBUG
  if ((m<1) || ::std::abs(m)>rank())
    throw BadIndex("component",m,rank());
#endif
  return tensor_values(size_t(m-1));
}
 
cmatrix RealSpinningHamiltonian::component(int m)
{
#ifndef NDEBUG
  if ((m<1) || ::std::abs(m)>rank())
    throw BadIndex("component",m,rank());
#endif
  return tensor_values(size_t(m-1));
}
const BaseList<complex> DiagonalSpinningHamiltonian::component(int m) const
{ 
#ifndef NDEBUG
  if ((m<1) || ::std::abs(m)>rank())
    throw BadIndex("component",m,rank());
#endif
  return tensor_values.row(size_t(m-1));
}
 
BaseList<complex> DiagonalSpinningHamiltonian::component(int m)
{
#ifndef NDEBUG
  if ((m<1) || ::std::abs(m)>rank())
    throw BadIndex("component",m,rank());
#endif
  return tensor_values.row(size_t(m-1));
}
const cmatrix& SpinningHamiltonian::component(int m) const
{ 
#ifndef NDEBUG
  if (::std::abs(m)>rank())
    throw BadIndex("component",m,rank());
#endif
  return tensor_values(size_t(m+rank()));
}
 
cmatrix& SpinningHamiltonian::component(int m)
{
#ifndef NDEBUG
  if (::std::abs(m)>rank())
    throw BadIndex("component",m,rank());
#endif
  return tensor_values(size_t(m+rank()));
}

} //namespace libcmatrix
#endif
