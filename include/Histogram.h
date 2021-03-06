#ifndef _Histogram_h_
#define _Histogram_h_

//Utility classes for accumulation of histograms
//Initialised with BaseList to which data is added using the add() member function

#include "BaseList.h"

namespace libcmatrix {

  inline double mod(double x,double y)
    {
      if (x<0) {
	double xn=floor((-x)/y)*y;
	return y+x+xn;
      }
      else {
	double xn=floor(x/y)*y;
	return x-xn;
      }
    }

  //! structure describing spectral histogram
  struct HistogramSpec {
    HistogramSpec()
      : nbins(0), rangef(0.0), minf(0.0) {}

    HistogramSpec(size_t nbinsv, double rangev) //!< default min frequency is -half bin width i.e. zero frequency is in middle of first bin
      : nbins(nbinsv), rangef(rangev), minf(-0.5*(rangev/nbins)) { validate(); }

    HistogramSpec(size_t nbinsv, double rangev, double minv)
      : nbins(nbinsv), rangef(rangev), minf(minv) { validate(); }
    
    size_t nbins; //!< number of bins
    double rangef; //!< spectral width
    double minf; //!< minimum frequency

    double midpoint(size_t bin) const {
      if (bin>=nbins)
	throw BadIndex("Histogram: bin index out of range",bin,nbins);
      return minf+((bin+0.5)/nbins)*rangef;
    }

    void validate() {
      if (rangef<=0.0 || nbins<1)
	throw InvalidParameter("HistogramSpec");
    }

    bool operator== (const HistogramSpec& a) const { return (nbins==a.nbins) && (rangef==a.rangef) && (minf==a.minf); }
    bool operator!= (const HistogramSpec& a) const { return (nbins!=a.nbins) || (rangef!=a.rangef) || (minf!=a.minf); }
  };

  template<typename T> class BaseHistogram : protected HistogramSpec {
  public:
    using HistogramSpec::midpoint;

     BaseHistogram(BaseList<T>& datav, double rangev)
       : HistogramSpec(datav.size(),rangev),
 	data(datav) { create(); }

    BaseHistogram(BaseList<T>& datav, double rangev, double minv)
      : HistogramSpec(datav.size(),rangev,minv),
	data(datav) { create(); }

    virtual ~BaseHistogram() {}

    virtual void add(const T&,double) =0;
    void addlost(const T& v) { lostv++; lostsumv+=v; } //!< functions added 1/4/15
    void resetlost() { lostv=0U; lostsumv=T(0); }
    BaseList<T> row() const { return data; }

    const HistogramSpec& specification() const { return *this; }

    double range() const { return rangef; }
    double minimum() const { return minf; }
    double maximum() const { return minf+rangef; }

    size_t lost() const { return lostv; }
    T lostsum() const { return lostsumv; }

  protected:
    
    BaseList<T> data;
    double scalef,deltaf;
    size_t lostv;
    T lostsumv;
    
  private:
    BaseHistogram<T>(const BaseHistogram<T>&);
    BaseHistogram<T>& operator=(const BaseHistogram<T>&); //!< make this private to prevent copies

    void create() {
      scalef=nbins/rangef;
      deltaf=rangef/nbins;
      resetlost();
      //      lostv=0U;
      //lostsumv=T(0);
    }
  };

#define LCM_HISTOGRAM_INTERFACE(T)\
    using BaseHistogram<T>::scalef;\
    using BaseHistogram<T>::rangef;\
    using BaseHistogram<T>::deltaf;\
    using BaseHistogram<T>::minf;\
    using BaseHistogram<T>::nbins;\
    using BaseHistogram<T>::data;\
    using BaseHistogram<T>::lostv;\
    using BaseHistogram<T>::lostsumv

  template<typename T> class Histogram : public BaseHistogram<T> {
  public:
    LCM_HISTOGRAM_INTERFACE(T);

    Histogram(BaseList<T> datav,double maxfv)
      : BaseHistogram<T>(datav,maxfv) {}

    Histogram(BaseList<T> datav,double rangev, double minfv)
      : BaseHistogram<T>(datav,rangev,minfv) {}
    
    void add(const T& a,double f) {
      const int nbin=int(scalef*(f-minf));
      if (nbin>=nbins || nbin<0)
	this->addlost(a);
      else
	data(nbin)+=a;
    }
  };
  
  template<typename T> class FoldingHistogram : public BaseHistogram<T> {
  public:
    LCM_HISTOGRAM_INTERFACE(T);

    FoldingHistogram(BaseList<T> datav,double maxfv)
      : BaseHistogram<T>(datav,maxfv) {}

    FoldingHistogram(BaseList<T> datav,double rangev, double minfv)
      : BaseHistogram<T>(datav,rangev,minfv) {}
    
    inline double fold(double f) { return mod(f-minf,rangef); }

    void add(const T& a,double f) {
      size_t nbin=int(scalef*fold(f));
      data((nbin==nbins) ? 0U : nbin)+=a;
    }
  };

  template<typename T> class InterpHistogram : public BaseHistogram<T> {
  public:
    InterpHistogram(BaseList<T>& datav,double maxfv)
      : BaseHistogram<T>(datav,maxfv) {}

    InterpHistogram(BaseList<T>& datav,double rangev, double minfv)
      : BaseHistogram<T>(datav,rangev,minfv) {}
    
    LCM_HISTOGRAM_INTERFACE(T);

    void add(const T& a,double f);
  };
  
  template<typename T> class FoldingInterpHistogram : public BaseHistogram<T> {
  public:
    FoldingInterpHistogram(BaseList<T> datav,double maxfv)
      : BaseHistogram<T>(datav,maxfv) { }

    FoldingInterpHistogram(BaseList<T> datav,double rangev,double minfv)
      : BaseHistogram<T>(datav,rangev,minfv) { }
    
    LCM_HISTOGRAM_INTERFACE(T);

    void add(const T& a,double F);
  };

//implementation details

  template<class T> void 
    InterpHistogram<T>::add(const T& a,double f) {
    f-=minf;
    int binlow=int(scalef*f-0.5);
    if (binlow>=nbins) {
      this->addlost(a);
      return;
    }
    const double frac=(f-binlow*deltaf)/deltaf-0.5; //DODGY
    const T bino=frac*a;
    if (binlow<0)
      lostsumv+=a-bino;
    else
      data(binlow)+=a-bino;
    if (++binlow==nbins)
      lostsumv+=bino;
    else
      data(binlow)+=bino;
  }
  
  template<typename T> void
    FoldingInterpHistogram<T>::add(const T& a,double f)
    { 
      f=mod(f-minf,rangef);
      int binlow=int(scalef*f-0.5);
      const double frac=(f-binlow*deltaf)/deltaf-0.5;
      const T bino=frac*a;
      if (binlow<0)
	data(nbins-1)+=a-bino;
      else
	data(binlow)+=a-bino;
      if (++binlow==nbins)
	data.front()+=bino;
      else
	data(binlow)+=bino;
    }
  
} //namespace libcmatrix

#endif
