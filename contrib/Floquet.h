/* Utility functions for dealing with Floquet simulations */
#undef LCM_SUPPRESS_VIEWS
#include "MAS.h"
#include "Histogram.h"
#include "BasePropagation.h"
#include "lcm_Propagation-common.h"

namespace libcmatrix {

  cmatrix floquet_hamiltonian(const SpinningHamiltonian&,int,double);
  cmatrix floquet_hamiltonian(const RealSpinningHamiltonian&,int,double);
  cmatrix floquet_hamiltonian(const BaseList<complex>& Hels,int N,double v);
  rmatrix floquet_hamiltonian(const BaseList<double>& Hels,int N,double v);

  cmatrix floquet_hamiltonian(const BaseList<cmatrix>&,int,double,const BaseList<cmatrix>&,int,double);
  rmatrix floquet_hamiltonian(const BaseList<rmatrix>&,int,double,const BaseList<rmatrix>&,int,double);
  
  cmatrix floquet_operator(const cmatrix&,int);
  rmatrix floquet_operator(const rmatrix&,int);

  double get_frequency(const cmatrix&,int);
  double get_frequency(const rmatrix&,int);
  void set_frequency(cmatrix&, int, double);
  void set_frequency(rmatrix&, int, double);

  size_t hilbert_order(size_t fdim,int order);

  typedef range floquet_block_t;

  inline range floquet_block(int N,int nh,int r) {
    if (r>N || (r<-N)) throw BadIndex("floquet_block");
    const int start=(N+r)*nh;
    return range(start,start+nh-1);
  }

  struct FloquetTDStash {
    RCmatrix V;
    List<double> eigs;
    List<complex> U;

    void set_H(const rmatrix& FHam,double dt) {
      hermitian_eigensystem(V.set_real(),eigs,FHam);
      propagator(U,eigs,dt);
    }
    void set_H(const cmatrix& FHam,double dt) {
      hermitian_eigensystem(V.set_complex(),eigs,FHam);
      propagator(U,eigs,dt);
    }
    size_t size() const { return eigs.length(); }
    void clear() { V.clear(); eigs.clear(); }    
    bool operator!() const { return eigs.empty(); }
  };

  struct FloquetFDStash {
    RCmatrix V;
    List<double> eigs;
    void operator= (const rmatrix& FHam) { hermitian_eigensystem(V.set_real(),eigs,FHam); }
    void operator= (const cmatrix& FHam) { hermitian_eigensystem(V.set_complex(),eigs,FHam); }
    size_t size() const { return eigs.length(); }
    void clear() { V.clear(); eigs.clear(); }
    bool operator!() const { return eigs.empty(); }
  };

  inline std::ostream& operator<< (std::ostream& ostr, const FloquetFDStash& a) {
    ostr << "Eigenvalues: " << a.eigs << '\n';
    ostr << "Transform matrix\n" << a.V << '\n';
    return ostr;
  }

  class BaseFloquetFID : public SwapStore<FloquetTDStash> {
   protected:
    int N;
    int verbose_;

    BaseFloquetFID(bool isgammav, int verbosev)
      : N(0), verbose_(verbosev), isgamma(isgammav) {}

  public:
    void set_H(char sel,const cmatrix& H,int Nv,double dt) { setRC(sel).set_H(H,dt); freqdt=dt*get_frequency(H,Nv); N=Nv; }
    void set_H(char sel,const rmatrix& H,int Nv,double dt) { setRC(sel).set_H(H,dt); freqdt=dt*get_frequency(H,Nv); N=Nv; }
    void set_H(const cmatrix& H,int Nv,double dt) { setdiagonal().set_H(H,dt); freqdt=dt*get_frequency(H,Nv); N=Nv; }
    void set_H(const rmatrix& H,int Nv,double dt) { setdiagonal().set_H(H,dt); freqdt=dt*get_frequency(H,Nv); N=Nv; }

    List<complex> FID(int npoints) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d); return d; }
    List<complex> FID(int npoints,const cmatrix& sigma0, const cmatrix& detect) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d,1.0,sigma0,detect); return d; }
    List<complex> FID(int npoints,const rmatrix& sigma0, const rmatrix& detect) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d,1.0,sigma0,detect); return d; }
     
    void add_FID(BaseList<complex>, double scale, const cmatrix&, const cmatrix&);
    void add_FID(BaseList<complex>, double scale, const rmatrix&, const rmatrix&);
    void add_FID(BaseList<complex>, double scale =1.0);

  private:
   template<class T> void add_FID_(BaseList<complex>, double, const Matrix<T>&, const Matrix<T>&);
    template<class TH,class TM> void add_FID__(BaseList<complex>& FID,double scale,const Matrix<TH>& VR,const Matrix<TH>& VC,const Matrix<TM>& sigma0,const Matrix<TM>& detect);
    template<class T> void BaseFloquetFID::add_FID__(BaseList<complex>& FID,double scale,const Matrix<T>& VR,const Matrix<T>& VC);

    double freqdt;
    const bool isgamma;
  };

  struct FloquetFID : public BaseFloquetFID {
    FloquetFID(int verbosev =0) : BaseFloquetFID(false,verbosev) {}
  };

  template<> struct SignalGenerator_traits<FloquetFID> {
    static const bool timedomain=true;
    static const bool allowdiagonal=false;
    static const bool allowreal=true;
    static const bool allowED=false;
    //can't really define input_dimensionality, allowH!
  };

  class GammaFloquetFID : public BaseFloquetFID { 
  public:
    GammaFloquetFID(int verbosev =0)
      : BaseFloquetFID(true,verbosev) {}

    List<double> FID_hermitian(int npoints,const cmatrix& sigma0, const cmatrix& detect) { List<double> d(npoints,0.0,mxflag::temporary); this->add_FID_hermitian(d,1.0,sigma0,detect); return d; }
    List<double> FID_hermitian(int npoints,const rmatrix& sigma0, const rmatrix& detect) { List<double> d(npoints,0.0,mxflag::temporary); this->add_FID_hermitian(d,1.0,sigma0,detect); return d; }
    void add_FID_hermitian(BaseList<double>, double scale, const cmatrix&, const cmatrix&);
    void add_FID_hermitian(BaseList<double>, double scale, const rmatrix&, const rmatrix&);
    void add_FID_hermitian(BaseList<complex>, complex scale, const cmatrix&, const cmatrix&);
    void add_FID_hermitian(BaseList<complex>, complex scale, const rmatrix&, const rmatrix&);

  private:
    template<class Td,class T> void add_FID_hermitian_(BaseList<Td>, Td, const Matrix<T>&, const Matrix<T>&);
    template<class Td,class TH,class TM> void add_FID_hermitian__(BaseList<Td>&,Td,const Matrix<TH>&,const Matrix<TM>&,const Matrix<TM>&);
  };

  template<> struct SignalGenerator_traits<GammaFloquetFID> {
    static const bool timedomain=true;
    static const bool allowdiagonal=false;
    static const bool allowreal=true;
    static const bool allowED=false;
  };

  class BaseFloquetSpectrum : public SpectrumIterator<complex>, public SwapStore<FloquetFDStash> {
 private:
  InhomoIter iter_;
  int verbose_;
  size_t rs,cs;
  int N,Nval;
  double f,freq;
  const bool isgamma;
  bool finished;
  RCmatrix sigma0t,detect,detectt,tmp,tmp2;
  BaseList<double> reigs,ceigs;

  void update_detect();
  template<class T> void observe_(const Matrix<T>&,const Matrix<T>&);
  template<class T> void set_frequency(const Matrix<T>& H, int order) { 
    const double freqv=get_frequency(H,order);
    if (freqv==0.0) throw InvalidParameter("FloquetSpectrum: frequency cannot be zero");
    freq=freqv;
    N=order;
  }

public:
  BaseFloquetSpectrum(bool isgammav,int verbosev)
    : verbose_(verbosev), freq(0),isgamma(isgammav), finished(true) {}

  void observe();
  void observe(const cmatrix& sigma0, const cmatrix& det) { observe_(sigma0,det); }
  void observe(const rmatrix& sigma0, const rmatrix& det) { observe_(sigma0,det); }
  
  void set_H(const cmatrix&, int order);
  void set_H(const rmatrix&, int order);
  void set_H(char,const cmatrix&, int order);
  void set_H(char,const rmatrix&, int order);

  bool operator() (complex& amp,double& freq);
  friend std::ostream& operator<< (std::ostream&,const BaseFloquetSpectrum&);
  void reset();
};

  struct FloquetSpectrum : public BaseFloquetSpectrum {
    FloquetSpectrum(int verbosev =0) : BaseFloquetSpectrum(false,verbosev) {}
  };
  struct GammaFloquetSpectrum : public BaseFloquetSpectrum {
    GammaFloquetSpectrum(int verbosev =0) : BaseFloquetSpectrum(true,verbosev) {}
  };

template<> struct SignalGenerator_traits<BaseFloquetSpectrum> {
  static const bool timedomain=false;
  static const bool allowreal=true;
  static const bool allowdiagonal=false;
  static const bool allowED=false;
  static const bool allownonED=true;
};

  std::ostream& operator<< (std::ostream&,const BaseFloquetSpectrum&);

} //namespace libcmatrix


