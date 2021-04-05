//Not to be included directly

#include "lcm_MultiHolder.h"
#include "lcm_Propagation-common.h"

namespace libcmatrix {

struct HomoStashTD {
  MultiHolder<complex> Us;

  HomoStashTD() {}

  bool operator!() const { return !Us; }
  bool iselement() const { return (Us.dimensions()==1); }
  size_t size() const { return (Us.dimensions()==3) ? Us.multimatrix3().dimension(1) : 1; }
  void clear() { if (Us.dimensions()>1) Us.clear(); }

  void operator= (const BaseList<cmatrix>& Usv) { Us=Usv; }
  void operator= (const BaseList<complex>& Usv) { Us=Usv; }

  friend std::ostream& operator<< (std::ostream& ostr, const HomoStashTD& a) { return ostr << a.Us; }
};

struct HomoStash : public HomoStashTD {
  List<double> eff_freq;
  HomoStash() : eff_freq(SCRATCH_SIZE) {}

  void operator= (const BaseList<cmatrix>& Us);
  void operator= (const BaseList<complex>& Us);

  friend std::ostream& operator << (std::ostream&, const HomoStash&);
};

struct HomoStashGammaTD : public HomoStashTD {
  List<complex> cycle_eigs;
  cmatrix D;
  
  HomoStashGammaTD() : cycle_eigs(SCRATCH_SIZE) {}

  void operator= (const BaseList<cmatrix>&);
  void operator= (const BaseList<complex>&);

  friend std::ostream& operator << (std::ostream&, const HomoStashGammaTD&);

  void clear() { HomoStashTD::clear(); D.clear(); }
};

 template<class T> class BasePeriodicSpectrum;
 template<class T> std::ostream& operator<< (std::ostream&, const BasePeriodicSpectrum<T>&);

  template<class T> class BasePeriodicSpectrum : public SpectrumIterator<T>, public SwapStore<HomoStash>, public PeriodicObj {
protected:
   explicit BasePeriodicSpectrum(size_t nobsv, bool usegammastepsv, int verbosev)
     : PeriodicObj(nobsv,verbosev), usegammasteps_(usegammastepsv), eigrestrict(false) {}
    //     row().showwarning=col().showwarning=(verbosev>0);
  
    bool usegammasteps_;
  
   List<T> sideband_;
   List<complex> gft;
   MultiHolder<complex> sigma_table,detect_table;
   
   virtual bool calcamps(size_t,size_t) =0;
   void set_steps(int);

   size_t dimR,dimC;
   size_t r,s,ncount;
   double diff,fscale;
   bool usefft,finished;
   List<complex> ft_facs;
   List<complex> scr;
   bool eigrestrict;
   int valsideband;
   size_t valrestrict;
   
   bool advance();
   void lock() { finished=true; }
   
   size_t whichrs() const { return (dimR==1) ? s : r; }

 public:

   template<class U> void add(BaseList<U> amps,U scale,size_t lr,size_t ls);
   template<class U> void add(List<U>& amps,U scale,size_t lr,size_t ls);

  void set_Us(const BaseList<cmatrix>&, double);
  void set_Us(const BaseList<complex>&, double);
  void set_Us(char sel,const BaseList<cmatrix>&, double);
  void set_Us(char sel,const BaseList<complex>&, double);

  double frequency(size_t lr,size_t ls) const
  { return fscale*(col().eff_freq(ls)-row().eff_freq(lr)); }

  double frequency(size_t lr,size_t ls,size_t n) const {
    if (n>=observations())
      throw BadIndex("frequency");
    return frequency(lr,ls)+n*speed();
  }
  void sideband(int);
   void reset();

  bool operator() (T& amp,double& freq);

  friend std::ostream& operator<< <>(std::ostream&, const BasePeriodicSpectrum&);
 };
 
 class PeriodicSpectrum : public BasePeriodicSpectrum<complex>
{
   bool calcamps(size_t,size_t);

public:
  explicit PeriodicSpectrum(int verbosev =0)
    : BasePeriodicSpectrum<complex>(0,false,verbosev) {}

  void observe(const cmatrix&, const cmatrix&);
  void observe(const rmatrix&, const rmatrix&);
  void observe(const BaseList<double>&, const BaseList<double>&);

  void observe(const cmatrix&);
  void observe(const rmatrix&);
  void observe(const BaseList<double>&);
  void observe();
};

template<> struct SignalGenerator_traits<PeriodicSpectrum> {
  static const bool timedomain=false;
  static const bool allowreal=true;
  static const bool allowdiagonal=true;
  static const size_t input_dimensionality=3;
  static const bool allowH=false;
  static const bool allowU=true;
  static const bool allowED=true;
  static const bool allownonED=true;
  static const bool gamma=false;
};

class GammaPeriodicSpectrum : public BasePeriodicSpectrum<complex> {
  bool calcamps(size_t,size_t);

public:
  explicit GammaPeriodicSpectrum(int nobsv, int verbosev =0)
    : BasePeriodicSpectrum<complex>(nobsv,false,verbosev) {}

  void observe(const cmatrix&, const cmatrix&);
  void observe(const rmatrix&, const rmatrix&);
  void observe(const BaseList<double>&, const BaseList<double>&);
};

template<> struct SignalGenerator_traits<GammaPeriodicSpectrum> {
  static const bool timedomain=false;
  static const bool allowreal=true;
  static const bool allowdiagonal=true;
  static const size_t input_dimensionality=3;
  static const bool allowED=false;
  static const bool allownonED=true;
  static const bool allowH=false;
  static const bool allowU=true;
  static const bool gamma=true;
};

 class GammaPeriodicSpectrumED : public BasePeriodicSpectrum<double>
{
  bool calcamps(size_t,size_t);

public:
  explicit GammaPeriodicSpectrumED(int nobsv, int verbosev =0)
    : BasePeriodicSpectrum<double>(nobsv,true,verbosev) {}

  void observe(const cmatrix&);
  void observe(const rmatrix&);
  void observe(const BaseList<double>&); 
  void observe();
};

template<> struct SignalGenerator_traits<GammaPeriodicSpectrumED> {
  static const bool timedomain=false;
  static const bool allowreal=true;
  static const bool allowdiagonal=true;
  static const size_t input_dimensionality=3;
  static const bool allowED=true;
  static const bool allownonED=false;
  static const bool allowH=false;
  static const bool allowU=true;
  static const bool gamma=true;
};

  class PeriodicFID : public SwapStore<HomoStashTD>, public PeriodicObj {
 public:
  explicit PeriodicFID(int verbosev =0)
    : PeriodicObj(0,verbosev) {}

    List<complex> FID(int npoints,const cmatrix& sigma0,const cmatrix& detect) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d,complex(1.0),sigma0,detect); return d; }
   List<complex> FID(int npoints) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d); return d; }
    List<double> FID(int npoints,const BaseList<double>& sigma0,const BaseList<double>& detect) { List<double> d(npoints,0.0,mxflag::temporary); this->add_FID(d,1.0,sigma0,detect); return d; }
    List<double> FID_hermitian(int npoints,const cmatrix& sigma0,const cmatrix& detect) { List<double> d(npoints,0.0,mxflag::temporary); this->add_FID_hermitian(d,1.0,sigma0,detect); return d; }  

  void set_Us(const BaseList<cmatrix>&);
  void set_Us(const BaseList<complex>&);
  void set_Us(char sel,const BaseList<cmatrix>&);
  void set_Us(char sel,const BaseList<complex>&);

    void add_FID(BaseList<complex>, complex =complex(1.0)) const;
  void add_FID(BaseList<complex>, complex, const cmatrix&, const cmatrix&) const;
  void add_FID(BaseList<double> FIDv, double scale, const BaseList<double>& sigma0, const BaseList<double>& detect) const { add_FID_(FIDv,scale,sigma0,detect); }
  void add_FID(BaseList<complex> FIDv, complex scale, const BaseList<double>& sigma0, const BaseList<double>& detect) const { add_FID_(FIDv,scale,sigma0,detect); }
  template<class T> void add_FID_hermitian(BaseList<T>, T scale, const cmatrix& sigma0, const cmatrix& detect) const;

  friend std::ostream& operator << (std::ostream&, const PeriodicFID&);

 private:
   //  PeriodicObj periodicobj_;
    //  int verbose_;
  mutable MultiHolder<complex> detect_tab;
  mutable cmatrix tmp;
  template<class T> void add_FID_(BaseList<T>, T, const BaseList<double>&, const BaseList<double>&) const;
};

template<> struct SignalGenerator_traits<PeriodicFID> {
  static const bool timedomain=true;
  static const bool allowdiagonal=true;
  static const bool allowreal=false;
  static const bool allowED=false;
  static const size_t input_dimensionality=3;
  static const bool allowH=false;
  static const bool allowU=true;
  static const bool gamma=false;
};

  class GammaPeriodicFID : public SwapStore<HomoStashGammaTD>, public PeriodicObj {
 public:
  explicit GammaPeriodicFID(int nobsv, int verbosev =0);

  void set_Us(const BaseList<cmatrix>&);
  void set_Us(const BaseList<complex>&);
  void set_Us(char sel,const BaseList<cmatrix>&);
  void set_Us(char sel,const BaseList<complex>&);

    int reduction_factor() const { return reduction_factor_; }
    void reduction_factor(int);

  List<complex> FID(int npoints,const cmatrix& sigma0,const cmatrix& detect) {
    List<complex> d(npoints,complex(0.0),mxflag::temporary);
    this->add_FID(d,complex(1.0),sigma0,detect);
    return d;
  }
  List<complex> FID(int npoints,const cmatrix& sigma0det) { 
    List<complex> d(npoints,complex(0.0),mxflag::temporary);
    this->add_FID(d,complex(1.0),sigma0det);
    return d;
  }
  List<double> FID(int npoints,const BaseList<double>& sigma0det) {
    List<double> d(npoints,0.0,mxflag::temporary); 
    this->add_FID(d,1.0,sigma0det);
    return d;
  }
  List<complex> FID(int npoints) {
    List<complex> d(npoints,complex(0.0),mxflag::temporary);
    this->add_FID(d);
    return d;
  }
  List<double> FID(int npoints,const BaseList<double>& sigma0,const BaseList<double>& detect) {
    List<double> d(npoints,0.0,mxflag::temporary);
    this->add_FID(d,1.0,sigma0,detect);
    return d;
  }
  List<double> FID_hermitian(int npoints,const cmatrix& sigma0,const cmatrix& detect) {
    List<double> d(npoints,0.0,mxflag::temporary);
    this->add_FID_hermitian(d,1.0,sigma0,detect); 
    return d;
  }  
  List<double> FID_hermitian(int npoints,const cmatrix& a) {
    List<double> d(npoints,0.0,mxflag::temporary); 
    this->add_FID_hermitian(d,1.0,a);
    return d;
  }  

    void add_FID(BaseList<complex>, complex =complex(1.0)) const;
  void add_FID(BaseList<complex>, complex, const cmatrix&) const;
  void add_FID(BaseList<complex>, complex, const cmatrix&, const cmatrix&) const;
  template<class T> void add_FID_hermitian(BaseList<T>, T scale, const cmatrix& sigma0det) const;
  template<class T> void add_FID_hermitian(BaseList<T>, T scale, const cmatrix& sigma0, const cmatrix& detect) const;
  void add_FID(BaseList<double> FIDv, double scale, const BaseList<double>& sigma0, const BaseList<double>& detect) const { add_FID_(FIDv,scale,sigma0,detect); }
  void add_FID(BaseList<complex> FIDv, complex scale, const BaseList<double>& sigma0, const BaseList<double>& detect) const { add_FID_(FIDv,scale,sigma0,detect); }
  void add_FID(BaseList<double> FIDv, double scale, const BaseList<double>& sigma0det) const { add_FID_(FIDv,scale,sigma0det); }
  void add_FID(BaseList<complex> FIDv, complex scale, const BaseList<double>& sigma0det) const { add_FID_(FIDv,scale,sigma0det); }

  friend std::ostream& operator << (std::ostream&,const GammaPeriodicFID&);

 private:
    int reduction_factor_; //!< time sampling factor (default 1)
  mutable MultiHolder<complex> sigma0_tab,detect_tab,S;
  mutable cmatrix R,tmp;

    complex getR(size_t,size_t) const; //!< return propagation factor including time scaling
  void make_R() const; //not really a const function!
  template<class T,class S> void add_FID_(BaseList<T>, S) const;
  template<class T> void add_FID_(BaseList<T>, T scale, const BaseList<double>&, const BaseList<double>&) const;
  template<class T> void add_FID_(BaseList<T>, T scale, const BaseList<double>&) const;

    template<class T> void doadd_hermitian(BaseList<T> FIDv,T scale, bool forceherm =false) const;
    void doadd(BaseList<complex> FIDv,complex scale) const;
    void doadd_compressed(BaseList<complex> FIDv,complex scale) const;

};

template<> struct SignalGenerator_traits<GammaPeriodicFID> {
  static const bool timedomain=true;
  static const bool allowdiagonal=true;
  static const bool allowreal=false;
  static const bool allowED=true;
  static const size_t input_dimensionality=3;
  static const bool allowH=false;
  static const bool allowU=true;
  static const bool gamma=true;
};

} //namespace libcmatrix
