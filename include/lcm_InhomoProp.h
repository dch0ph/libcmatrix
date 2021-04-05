//Not to be included directly
#include "Warnings.h"

namespace libcmatrix {

class DynamicPhase;

//Base object for both time- and frequency- domain MAS inhomogeneous propagation

  class MASInhomoObj : public SwapStore<InhomoStash>, public PeriodicObj {  
public:
  friend struct InhomoHelper_;
  
    static Warning<> ignoring_period_warning;

  bool arebothreal() const { return (!row().iscomplex() && !col().iscomplex()); }
/*   bool operator!() const { */
/*     return SwapStore<InhomoStash>::operator!() || !(row().valid) || !(col().valid()); */
/*   } */

  bool explicit_Us() const {
    return row().explicit_Us() || col().explicit_Us();
  }
  
  void set_Hs(const rmatrix&);
  void set_Hs(char,const rmatrix&);
  template<class M> void set_H(char sel,const M& H) {
    setRC(sel).set_H(H);
  }
  template<class M> void set_H(const M& H) {
    setdiagonal().set_H(H);
  }

  //set propagators explicitly
  void set_Us(const cmatrix&, double =0.0);
  void set_Us(char,const cmatrix&, double =0.0);
  
  size_t interactions() const { 
    const size_t nr=row().interactions();
    if (nr!=col().interactions())
      throw Mismatch("interactions");
    return nr;
  }

  friend std::ostream& operator<< (std::ostream&,const MASInhomoObj&);

  void reset() { updatephases(); }
    void sideband(int);

  void set_phases(const BaseList<double>& phases,double =0.0);
  void set_phases(const rmatrix& phases,double =0.0);
  void set_phases(const BaseList<DynamicPhase>&);
  void set_phases(const DynamicPhase& ph) { set_phases(BaseList<DynamicPhase>(1,const_cast<DynamicPhase*>(&ph))); }

 protected:
  MASInhomoObj(size_t nobsv,bool tdomainv,int verbosev)
    : PeriodicObj(nobsv,verbosev),
      tdomain(tdomainv), restrictsideband(false) {};

  mutable List<complex> sideband_;
    bool tdomain;
    
    bool restrictsideband;
    int valsideband;
    size_t valrestrict;

  rmatrix phase;
    virtual void set_steps(int);

  void updatephases() const;
  void invalidatephases();
  complex calcps(size_t,size_t) const; //not strictly const!
  complex getp(size_t,const complex&) const;
};

class BaseMASInhomoSpectrum : public MASInhomoObj {
public:
  double frequency(size_t lr,size_t ls) const;
  double frequency(size_t lr,size_t ls,size_t n) const {
    if (n>=this->observations()) 
      throw BadIndex("frequency");
    return frequency(lr,ls)+n*this->speed();
  }
 friend std::ostream& operator<< (std::ostream&,const BaseMASInhomoSpectrum&);  

protected:
  BaseMASInhomoSpectrum(int nobsv, int verbosev)
    : MASInhomoObj(nobsv,false,verbosev), scale_factor(0.0), finished(true) {}

  InhomoIter iter_;
  double scale_factor;
  bool usefft,finished;
  double viso,diff,nextfreq;
  size_t ncount;
  bool havenext;

  List<complex> ft_facs;
  List<complex> scr;

  void lock() { finished=true; }
  void reset_();
  void set_steps(int);
  void calcamps_(size_t,size_t);
  void invalidatephases() { row().invalidatephases(); col().invalidatephases(); }
};

 template<typename T> class BaseSingleInhomo;
 template<typename T> std::ostream& operator<< (std::ostream&,const BaseSingleInhomo<T>&);

template<typename T> class BaseSingleInhomo : public BaseMASInhomoSpectrum, public SpectrumIterator<complex> {
 public:
   typedef T transition_type;
   
  bool operator() (complex& amp,double& freq);

  void add(BaseList<complex>&,size_t,size_t);
  void add(List<complex>&, size_t,size_t);
  //  double zero_frequency() const { return isdiagonal() ? real(A(0U,0U)) : 0.0; }

  void reset();  
  friend std::ostream& operator<< <>(std::ostream&,const BaseSingleInhomo &);

  const Matrix<T>& amplitudes() const {
    verifyclean();
    return A;
  }

 protected:
  BaseSingleInhomo(int verbosev)
    : BaseMASInhomoSpectrum(0,verbosev) {}

  Matrix<T> A;
  complex ampfac,nextampfac,nextamp;

 private:

  void calcamps(size_t lr,size_t ls) { calcamps_(lr,ls); }
};

//   template<class Obj> class InhomogeneousAccumulator : public SpectrumIterator<typename Obj::amplitude_type> {
//   public:
//     typedef typename Obj::amplitude_type amplitude_type;
//     InhomogeneousAccumulator(Obj& objv) : obj_(objv), finished_(true) {}

//     void add();
//     typedef MultiMatrix<amplitude_type,3> store_type;
//     bool operator()(T& amp, double& f);   
//   private:
//     Obj& obj_;
//     bool finished_;
//     size_t r_,c_,n_;
//     size_t curr_,curc_,curn_;
//     store_type amps;
//     void advance();
//   };

//   template<class Obj> bool operator()(amplitude_type& amp, double& f)
//   {
//     while (!finished_) {
//       amp=amps(curr_,curc_,curn_);
//       if (amp!=amplitude_type(0.0)) {
// 	f=obj_.frequency(curr_,curc_);
// 	advance();
// 	return true;
//       }
//     } 
//     return false;
//   }

//   template<class Obj> void advance() {
//     if (++curn_==n_) {
//       curn_=0;
//       if (++curc_==c_) {
// 	curc_=0;
// 	if (++curr_==r_)
// 	  finished_=true;
//       }
//     }
//   }

//   template<class Obj> void
//   InhomogeneousAccumulator<Obj>::add() {
//     if (amps.empty()) {
//       r_=obj_.rows();
//       c_=obj_.cols();
//       n_=obj_.observations();
//       amps.create(r_,c_,n_,amplitude_type(0.0));
//       curr_=curc_=curn_=0;
//       finished_=false;
//     }
//     else {
//       if (finished_)
// 	throw Failed("InhomogeneousAccumulator: can't add new transitions after iterator used");
//       if ((r_!=obj_.rows()) || (c_!=obj_.cols()))
// 	throw Failed("InhomogeneousAccumulator: object has changed shape!");
//     }
//       for (size_t r=r_;r--;)
// 	for (size_t c=r;c--;)
// 	  obj_.add(amps(r,c),r,c);
//   }      
    
class InhomogeneousSpectrum : public BaseSingleInhomo<complex> {
 public:
  InhomogeneousSpectrum(int verbosev =0)
    : BaseSingleInhomo<complex>(verbosev) {}

    void observe(const cmatrix& sigma0,const cmatrix& det);
    void observe(const rmatrix& sigma0,const rmatrix& det);
    void observe(const BaseList<double>& sigma0,const BaseList<double>& det);
};

template<> struct SignalGenerator_traits<InhomogeneousSpectrum> {
  static const size_t input_dimensionality=2;
  static const bool timedomain=false;
  static const bool allowH=false; //Both H and U false - flags not normal behaviour!
  static const bool allowU=false;
  static const bool allowreal=true;
  static const bool allowED=false;
  static const bool allownonED=true;
  static const bool allowdiagonal=true;
  static const bool gamma=false;
};

class InhomogeneousFID : public MASInhomoObj, public InhomoTDBase_ {
 public:
  InhomogeneousFID(int verbosev =0)
    : MASInhomoObj(0,true,verbosev) {}

    friend std::ostream& operator<< (std::ostream&, const InhomogeneousFID&);
   
 protected:
  cmatrix A;
  void observe();
  void observe(const rmatrix& sigma0det);
  void observe(const cmatrix& sigma0det);
  void observe(const BaseList<double>& sigma0det);
  void observe(const rmatrix& sigma0, const rmatrix& detect);
  void observe(const cmatrix& sigma0, const cmatrix& detect);
  void observe(const BaseList<double>& sigma0, const BaseList<double>& detect);
  
  void add_FID_hermitian_(BaseList<double> FID,double =1.0) const;
  void add_FID_hermitian_(BaseList<complex> FID,complex =complex(1.0)) const;
  void add_FID_(BaseList<complex> FID,complex scale =complex(1.0)) const;
};

template<> struct SignalGenerator_traits<InhomogeneousFID> {
  static const size_t input_dimensionality=2;
  static const bool timedomain=true;
  static const bool allowH=false; //!< Both H and U false - flags not normal behaviour!
  static const bool allowU=false;
  static const bool usegenerator=false; //!< distinguish from AsychronousFID
  static const bool allowreal=true;
  static const bool allowED=true;
  static const bool allowdiagonal=true;
  static const bool gamma=false;
};
 
class InhomogeneousSpectrumED : public BaseSingleInhomo<double> {
 public:
  InhomogeneousSpectrumED(int verbosev =0)
    : BaseSingleInhomo<double>(verbosev) {}

    void observe(const cmatrix& sigma0det);
    void observe(const rmatrix& sigma0det);
    void observe(const BaseList<double>& sigma0det);
    void observe();
};

template<> struct SignalGenerator_traits<InhomogeneousSpectrumED> {
  static const size_t input_dimensionality=2;
  static const bool timedomain=false;
  static const bool allowH=false; //Both H and U false - flags not normal behaviour!
  static const bool allowU=false;
  static const bool allowreal=true;
  static const bool allowED=true;
  static const bool allownonED=false;
  static const bool allowdiagonal=true;
  static const bool gamma=false;
};

 template<typename T> class BaseGammaInhomo;
 template<typename T> std::ostream& operator<< (std::ostream&,const BaseGammaInhomo<T>&);

template<typename T> class BaseGammaInhomo : public BaseMASInhomoSpectrum, public SpectrumIterator<T>
  {
public:
  void add(BaseList<T>&,size_t,size_t);
  void add(List<T>&, size_t,size_t);
  bool operator() (T& amp,double& freq);
    //    double zero_frequency() const { return isdiagonal() ? real(A(0U,0U)) : 0.0; }

  friend std::ostream& operator<< <>(std::ostream&,const BaseGammaInhomo&);
  void reset();

  typedef T transition_type;

  const Matrix<T>& amplitudes() const {
    verifyclean();
    return A;
  }

protected:

   BaseGammaInhomo(int nobsv, int verbosev)
     : BaseMASInhomoSpectrum(nobsv,verbosev), amps(nobsv) {
     if (nobsv<=0)
       throw Failed("Attempt to initial Gamma object with observations <=0!");
   }

  Matrix<T> A;

private:
  T ampfac,nextampfac,nextamp;
  ScratchList<double> amps;

  void calcamps(size_t,size_t);
};

class GammaInhomogeneousSpectrumED : public BaseGammaInhomo<double> {
 public:
  GammaInhomogeneousSpectrumED(int nobsv, int verbosev =0)
    : BaseGammaInhomo<double>(nobsv,verbosev) {}

    void observe(const cmatrix& sigma0det);
    void observe(const rmatrix& sigma0det);
    void observe(const BaseList<double>& sigma0det);
    void observe();
};

template<> struct SignalGenerator_traits<GammaInhomogeneousSpectrumED> {
  static const size_t input_dimensionality=2;
  static const bool timedomain=false;
  static const bool allowH=false; //Both H and U false - flags not normal behaviour!
  static const bool allowU=false;
  static const bool allowreal=true;
  static const bool allowED=true;
  static const bool allownonED=false;
  static const bool allowdiagonal=true;
  static const bool gamma=true;
};
 
class GammaInhomogeneousSpectrum : public BaseGammaInhomo<complex> {
 public:
  GammaInhomogeneousSpectrum(int nobsv, int verbosev =0)
    : BaseGammaInhomo<complex>(nobsv,verbosev) {}

    void observe(const cmatrix& sigma0, const cmatrix& det);
    void observe(const rmatrix& sigma0, const rmatrix& det);
    void observe(const BaseList<double>& sigma0,const BaseList<double>& det);
};

template<> struct SignalGenerator_traits<GammaInhomogeneousSpectrum> {
  static const size_t input_dimensionality=2;
  static const bool timedomain=false;
  static const bool allowH=false; //Both H and U false - flags not normal behaviour!
  static const bool allowU=false;
  static const bool allowreal=true;
  static const bool allowED=false;
  static const bool allownonED=true;
  static const bool allowdiagonal=true;
  static const bool gamma=true;
};
 
class GammaInhomogeneousFID : public MASInhomoObj, public InhomoTDBase_ {

  cmatrix A;
  void observe();
  void observe(const rmatrix& sigma0det);
  void observe(const cmatrix& sigma0det);
  void observe(const BaseList<double>& sigma0det);
  void observe(const rmatrix& sigma0, const rmatrix& detect);
  void observe(const cmatrix& sigma0, const cmatrix& detect);
  void observe(const BaseList<double>& sigma0, const BaseList<double>& detect);

 public:
  GammaInhomogeneousFID(int nobsv, int verbosev =0)
    : MASInhomoObj(nobsv,true,verbosev) {
     if (nobsv<=0)
       throw Failed("Attempt to initial GammaInhomogeneousFID with observations <=0!");
  }
    friend std::ostream& operator<< (std::ostream&, const GammaInhomogeneousFID&);
protected:
  void add_FID_(BaseList<complex>, complex =complex(1.0)) const;
  void add_FID_hermitian_(BaseList<double>, double =1.0) const;
  void add_FID_hermitian_(BaseList<complex>, complex =complex(1.0)) const;
};

template<> struct SignalGenerator_traits<GammaInhomogeneousFID> {
  static const size_t input_dimensionality=2;
  static const bool timedomain=true;
  static const bool allowH=false; //!< Both H and U false - flags not normal behaviour!
  static const bool allowU=false;
  static const bool usegenerator=false; //!< distinguish from AsychronousFID
  static const bool allowreal=true;
  static const bool allowED=true;
  static const bool allowdiagonal=true;
  static const bool gamma=true;
};
 
} //namespace libcmatrix
