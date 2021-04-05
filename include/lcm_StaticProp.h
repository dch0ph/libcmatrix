//Not to be included directly

#include "ScratchList.h"
#include "lcm_Propagation-common.h"

// Static (block) propagation
namespace libcmatrix {

struct StashU : public cmatrix {
  bool iselement() const { return (rows()==1); }
  size_t size() const { return rows(); }

  void operator= (const cmatrix& Uv) {
    if (!issquare(Uv)) throw NotSquare("invalid propagator");
    cmatrix::operator=(Uv);
  }
  void operator= (const complex& Uv) { create(1,1,Uv); }
};

struct StashH : public RCmatrix {
  StashH() : dtstore(0.0) {} //!< flag dt not set

  ScratchList<double> eff_freqs;
  ScratchList<complex> evol_facs;
  double dtstore;

  bool operator!() const { return evol_facs.empty(); }
  friend std::ostream& operator << (std::ostream&,const StashH&);
  void finish_set_(double);

  size_t size() const { 
    if (!*this)
      throw Failed("Hamiltonian not set");
    return evol_facs.size();
  }
  bool iselement() const { 
    if (!*this)
      throw Failed("Hamiltonian not set");
    return RCmatrix::operator!();
  }

  void set_H(const cmatrix&,const BaseList<double>&,double);
  void set_H(const rmatrix&,const BaseList<double>&,double);
  void set_H(const cmatrix&,double);
  void set_H(const rmatrix&,double);  
  void set_H(const BaseList<double>&,double);
  void set_H(double,double);
  void set_U(const cmatrix&);
};

  class StaticFID_H : public InhomoTDBase_, public SwapStore<StashH> {
public:
    explicit StaticFID_H(int verbosev =0) 
      : larmor_f(0.0), min_f(0.0), max_f(0.0), verbose_(verbosev) {}

    bool arebothreal() const { return (!row().iscomplex() && !col().iscomplex()); }

  friend struct InhomoHelper_;

  void larmor(double larmor_);
  void larmor(double larmor_,double,double);

    template<class M> void set_H(char sel,const M& H,double dt) { setRC(sel).set_H(H,dt); }
    template<class M> void set_H(const M& H,double dt) { setdiagonal().set_H(H,dt); }
    template<class M> void set_H(char sel,const M& V,const BaseList<double>& eigs, double dt) { setRC(sel).set_H(V,eigs,dt); }
    void set_U(char sel, const cmatrix& U) { setRC(sel).set_U(U); }
    void set_U(const cmatrix& U) { setdiagonal().set_U(U); }
    void set_H(const cmatrix& V,const BaseList<double>& eigs, double dt) { setdiagonal().set_H(V,eigs,dt); }
    void set_H(const rmatrix& V,const BaseList<double>& eigs, double dt) { setdiagonal().set_H(V,eigs,dt); }

private:
  cmatrix A;
  static void reset() {} //empty reset callback
  void observe();
  void observe(const rmatrix& sigma0det);
  void observe(const cmatrix& sigma0det);
  void observe(const BaseList<double>& sigma0det);
  void observe(const rmatrix& sigma0, const rmatrix& detect);
  void observe(const cmatrix& sigma0, const cmatrix& detect);
  void observe(const BaseList<double>& sigma0, const BaseList<double>& detect);

  double get_dwell() const;

  void add_FID_(BaseList<complex>, complex) const;
    //virtual so can't be templated
    inline void add_FID_hermitian_(BaseList<double> FID_, double scale) const
    { add_FID_hermitian__(FID_,scale); }
    inline void add_FID_hermitian_(BaseList<complex> FID_, complex scale) const
    { add_FID_hermitian__(FID_,scale); }

    template<class T> void add_FID_hermitian__(BaseList<T>, T scale) const;

  double larmor_f,min_f,max_f;
  int verbose_;
};

template<> struct SignalGenerator_traits<StaticFID_H> {
  static const size_t input_dimensionality=2;
  static const bool timedomain=true;
  static const bool allowH=true;
  static const bool allowU=true;
  static const bool allowreal=true;
  static const bool allowED=true;
  static const bool allowdiagonal=true;
  static const bool gamma=false; //Gamma averaging not relevant, but indicates that no initialisation argument is needed
};

 class StaticBase_ : public SwapStore<InhomoStash> {
 protected:
  int verbose_;
  InhomoIter iter_;
  BaseList<double> eff_freqsR;
  BaseList<double> eff_freqsC;
  double larmor_f,min_f,max_f;

  template<class T> bool get_trans(T&, double&, const Matrix<T>&);
 public:
  StaticBase_(int verbosev =0)
    : verbose_(verbosev), larmor_f(0.0), min_f(0.0), max_f(0.0) {}

    friend struct InhomoHelper_;
    
    void larmor(double);
    void larmor(double,double,double);

    bool arebothreal() const { return (!row().iscomplex() && !col().iscomplex()); }
    
    template<class M> void set_U(char sel,const M& U, double period) {
      setRC(sel).set_U(U,period);
    }
    template<class M> void set_H(char sel,const M& H) {
      setRC(sel).set_H(H);
    }
   template<typename T> void set_H(char sel,const Matrix<T>& V, const BaseList<double>& eigs) {
     setRC(sel).set_H(V,eigs);
   }
   template<class M> void set_H(const M& H) {
     setdiagonal().set_H(H);
   }
   template<typename T> void set_H(const Matrix<T>& V, const BaseList<double>& eigs) {
     setdiagonal().set_H(V,eigs);
   }
   void set_U(const cmatrix& U,double period) {
     setdiagonal().set_U(U,period);
   }
};
 

class StaticSpectrum : public StaticBase_, public SpectrumIterator<complex> {
  cmatrix AC;
public:
  explicit StaticSpectrum(int verbosev =0)
    : StaticBase_(verbosev) {}

  void observe(const cmatrix& sigma0, const cmatrix& det);
  void observe(const rmatrix& sigma0, const rmatrix& det);
  void observe(const BaseList<double>& sigma0,const BaseList<double>& det);
  
  bool operator() (complex& amp,double& freq);
  void print(std::ostream& =std::cout) const; 
  void reset();
};

template<> struct SignalGenerator_traits<StaticSpectrum> {
  static const bool timedomain=false;
  static const bool allowreal=false;
  static const bool allowdiagonal=false;
  static const size_t input_dimensionality=2;
  static const bool allowED=false;
  static const bool allownonED=true;
  static const bool allowH=true;
  static const bool allowU=true;
  static const bool gamma=false; //Gamma averaging not relevant, but indicates that no initialisation argument is needed
};

 class StaticFID_U : public SwapStore<StashU> {
public:
  explicit StaticFID_U(int =0) {}

  template<class M> void set_U(char sel,const M& U) { setRC(sel)=U; }
  void set_U(const cmatrix& U) { setdiagonal()=U; }

   List<complex> FID(int npoints) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d,complex(1.0)); return d; }
   List<complex> FID(int npoints,const cmatrix& sigma0,const cmatrix& detect) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d,complex(1.0),sigma0,detect); return d; }
   List<double> FID_hermitian(int npoints,const cmatrix& sigma0,const cmatrix& detect) { List<double> d(npoints,0.0,mxflag::temporary); this->add_FID_hermitian(d,1.0,sigma0,detect); return d; }  

   void add_FID(BaseList<complex>, complex scale);
   void add_FID(BaseList<complex>, complex scale, const cmatrix&, const cmatrix&);
   template<class T> void add_FID_hermitian(BaseList<T>, T, const cmatrix& ,const cmatrix&);

private:
  cmatrix sigma,tmp;
};

template<> struct SignalGenerator_traits<StaticFID_U> {
  static const size_t input_dimensionality=2;
  static const bool timedomain=true;
  static const bool allowH=false;
  static const bool allowU=true;
  static const bool allowdiagonal=false;
  static const bool allowreal=false;
  static const bool allowED=false;
  static const bool gamma=false; //Gamma averaging not relevant, but indicates that no initialisation argument is needed
};

inline std::ostream& operator<< (std::ostream& ostr,const StaticSpectrum &a) { a.print(ostr); return ostr; }

// special case for matching detection operator and initial density matrix

class StaticSpectrumED : public StaticBase_, public SpectrumIterator<double> {
  rmatrix AR;
public:
  explicit StaticSpectrumED(int verbosev =0)
    : StaticBase_(verbosev) {}

  void observe(const cmatrix& sigma0det);
  void observe(const rmatrix& sigma0det);
  void observe(const BaseList<double> &sigma0det);
  void observe();
  void reset();
  bool operator() (double& amp,double& freq);

  void print(std::ostream& =std::cout) const; 
};

template<> struct SignalGenerator_traits<StaticSpectrumED> {
  static const bool timedomain=false;
  static const bool allowreal=false;
  static const bool allowdiagonal=true;
  static const size_t input_dimensionality=2;
  static const bool allowED=true;
  static const bool allownonED=false;
  static const bool allowH=true;
  static const bool allowU=true;
  static const bool gamma=false; //Gamma averaging not relevant, but indicates that no initialisation argument is needed
};

inline std::ostream& operator<< (std::ostream& ostr,const StaticSpectrumED &a) { a.print(ostr); return ostr; }

} //namespace libcmatrix

