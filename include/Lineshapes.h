#ifndef LCM_Lineshapes_h_
#define LCM_Lineshapes_h_

/*! \file
  \brief Functionality for direct addition of lineshapes to frequency histograms

*/

#include "List.h"
#include "Histogram.h"
#include <map>

namespace libcmatrix {

  //! simple struct defining parameters for Lineshape function

  struct LineshapeSpec {
    enum { nofold=1, //!< lineshape is not folded in histogram
	   verbose=2, //!< enable verbose output (otherwise silent)
	   variable=4 //!< allow output size to vary
    };
    
    LineshapeSpec() : 
      cutoff(0.0),
      resolution_steps(10U), //!< default is subdivision into 10 steps 
      cache_maximum(1U), //!< default cache single lineshape
      flags(0) {}
    
    double cutoff; //!< lineshape foot cutoff as fraction of max (0 no cutoff)
    //double frequency_resolution;
    size_t resolution_steps; //!< subdivision of histogram bins
    size_t cache_maximum; //!< max number of lineshapes to cache
    int flags; //!< flags

    void validate() const; //!< throw ::InvalidParamater exception if any values are out of range
  };
  
  //! (abstract) base class for lineshape generators
  class BaseLineshape 
  {
  public:
    explicit BaseLineshape(const LineshapeSpec&); //!< initialise from general lineshape parameters
    virtual ~BaseLineshape() {}

    virtual BaseLineshape* clone() const =0;
    
    virtual int operator()(List<double>& dest, const HistogramSpec& spec, double freq =0.0) const =0;
    
    virtual void print(std::ostream& ostr =std::cout) const =0; //!< dump state
    
  protected:
    LineshapeSpec spec_; //!< copy of general parameters
    typedef int hash_t; //!< hash type for identifying distinct lineshapes
    hash_t hash(double offset, int) const; //!< create hash from (fractional) offset and bin offset
    void reversehash(double& offset, int&, hash_t) const; //!< undo hash
    //    static void copyroll(BaseList<double>& dest, const BaseList<double>& source, size_t width, int roll); //!< roll lineshape centred on middle bin and with specified half-width 
    virtual void reset() =0;  //!< clear caches etc.
    void flags(int);
    void reset_counts(); //!< clear call counts etc.
    
    mutable size_t callcount_; //!< count of function calls;
    mutable size_t hitcount_; //!< cache hit count
  };
  
  //! class for generating Lorentzian and/or Gaussians (NOT thread-safe)

  class LorentzianGaussian : public BaseLineshape
  {
  public:
    static const double PureLorentzian;
    static const double PureGaussian;
    
    LorentzianGaussian(double lw, double frac =PureLorentzian, const LineshapeSpec& =LineshapeSpec()); //!< initialise with linewidth \a lw and G/L fraction \a frac
    
    BaseLineshape* clone() const { return new LorentzianGaussian(*this); }
    
    void lw(double); //!< set linewidth
    void fraction(double); //!< set G/L fraction
    
    int operator()(List<double>&, const HistogramSpec& spec, double freq =0.0) const; //!< lineshape generator
    
    void reset(); //!< clear cache etc.
    //    bool operator!() const { return cachet_.empty(); }
    void print(std::ostream& ostr =std::cout) const; //!< output state
    
  private:
    //! private structure for cache entry
    struct cache_t {
      List<double> gauss; //!< unit gaussian lineshape
      List<double> lorentz; //!< unit lorentzian lineshape
      void reset(); //!< invalidate cached lineshapes
      void swap(cache_t&); //!< swap contents
      void make(const LorentzianGaussian&, hash_t); //!< force make 
      void ensure(const LorentzianGaussian&, hash_t); //!< update invalid entries
      void operator()(List<double>& dest, double frac) const; //!< create lineshape
    };
    
    double lw_; //!< linewidth
    double frac_; //!< G/L fraction
    double a_,a2_,ls_; //!< cached values for lorentzval
    double gs_,ga_; //!< cached values for gaussval

    double lorentzval(double f) const; //!< return Lorentzian at offset frequency \a f
    double lorentz0() const { return ls_/a2_; } //!< maximum of Lorentzian
    double gaussval(double f) const;
    double gauss0() const { return gs_; } //!< maximum of Gaussian

    typedef std::map<BaseLineshape::hash_t,cache_t> cachemap_t; //!< type of lineshape cache
    mutable cachemap_t cachemap_; //!< lineshape cache
    mutable cache_t tmp_; //!< temporary for unit lineshapes
    mutable List<double> cachet_; //!< temporary for total lineshape
    mutable HistogramSpec speccache_; //!< specifications of cached lineshapes
    mutable size_t Lwidth_; //!< half-width of Lorentzian (0 if not used)
    mutable size_t Gwidth_; //!< half-width of gaussian (0 if not used)
    
    size_t gethalfwidth() const; //!< return half width of current lineshape (0 if unknown)
    size_t getmid() const { return speccache_.nbins/2; } //!< return central bin of histogram
    void getbin(size_t&, double&, hash_t) const; //!< get actual bin position, frequency offset
    void setlw_(double); //!< \internal set linewidth
    void setfraction_(double); //!< \internal set G/L fraction
    void makel_(List<double>&, size_t, double) const; //!< make Lorentzian of unit intensity
    void makeg_(List<double>&, size_t, double) const; //!< make Gaussian of unit intensity
    typedef double (LorentzianGaussian::*LineshapeFunc_t)(double) const;
    size_t make_(List<double>& dest,size_t mid, double offset, const char* name, LineshapeFunc_t func,double zeroval) const; //!< raw lineshape creation
  };

  //! class for accumulating histogram with distinct linewidth
  template<typename T> class LineshapeHistogram : public BaseHistogram<T>
  {
  public:
    LineshapeHistogram(BaseList<T> datav, BaseLineshape& lshapev, double maxv)
      : BaseHistogram<T>(datav,maxv), lineshape_(lshapev) { reset(); }
    
    LineshapeHistogram(BaseList<T> datav, BaseLineshape& lshapev, double rangev, double minfv)
      : BaseHistogram<T>(datav,rangev, minfv), lineshape_(lshapev) { reset(); }
    
    using BaseHistogram<T>::data;

    void add(const T& a, double f) {
      const int offset=lineshape_(tmp_,this->specification(),f);
      if (offset>=0) {
	assert(offset+tmp_.size()<=data.size());
	BaseList<T> dest(data(range(offset,offset+tmp_.size()-1)));
	mla(dest,a,tmp_);
      }
      else
	this->addlost(a);
    }
    
    // size_t lost() const { return lostv; }
    //T lostsum() const { return lostsumv; }

  private:
    void reset() { BaseHistogram<T>::resetlost(); }

    BaseLineshape& lineshape_; //!< non-const ref not essential but prevents inappropriate copy and reduces risk of referencing a temporary object
    List<double> tmp_;
    //    size_t lostv;
    //T lostsumv;
  };

  void cubic_interpolate(BaseList<double> dest, double newstart, double newdx, const BaseList<double>& source, double oldstart, double olddx, double padding =0.0);

  List<double> cubic_interpolate(size_t newnp, double newstart, double newdx, const BaseList<double>& source, double oldstart, double olddx, double padding =0.0);

  void cubic_interpolate(BaseList<double>, const BaseList<double>&);
  List<double> cubic_interpolate(size_t, const BaseList<double>&);

  void cubic_interpolate(BaseList<complex> dest, double newstart, double newdx, const BaseList<complex>& source, double oldstart, double olddx, complex padding =complex(0.0));

  List<complex> cubic_interpolate(size_t newnp, double newstart, double newdx, const BaseList<complex>& source, double oldstart, double olddx, complex padding =complex(0.0));

  void cubic_interpolate(BaseList<complex>, const BaseList<complex>&);
  List<complex> cubic_interpolate(size_t, const BaseList<complex>&);

  extern int verbose_interpolate; //!< verbosity for interpolation

} //namespace libcmatrix

#endif
