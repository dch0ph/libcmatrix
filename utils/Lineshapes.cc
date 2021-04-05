/*! \file
  \brief Functionality for direct addition of lineshapes to frequency histograms
*/

#undef LCM_SUPPRESS_VIEWS
#include "Lineshapes.h"
#include <cmath>

namespace libcmatrix {

  int verbose_interpolate=0;

  BaseLineshape::hash_t BaseLineshape::hash(double off, int binkey) const
  {
    if (binkey<0) 
      return -hash(off,-binkey);
    else
      return int(0.5+spec_.resolution_steps*(binkey+off));
  }

  void BaseLineshape::reversehash(double& offset, int& binkey, hash_t hash) const
  {
    if (hash<0) {
      reversehash(offset,binkey,-hash);
      binkey=-binkey;
    }
    else {
      binkey=hash/int(spec_.resolution_steps);
      offset=double(hash-spec_.resolution_steps*binkey)/spec_.resolution_steps;
    }
  }

  const double LorentzianGaussian::PureLorentzian=0.0;
  const double LorentzianGaussian::PureGaussian=1.0;

  void LorentzianGaussian::cache_t::reset()
  {
    gauss.create(0U);
    lorentz.create(0U);
  }

 void LorentzianGaussian::cache_t::make(const LorentzianGaussian& obj, hash_t key)
 {
   reset(); //!< invalidate cached values
   ensure(obj,key);
 }

  void LorentzianGaussian::cache_t::swap(cache_t& a)
  {
    gauss.swap(a.gauss);
    lorentz.swap(a.lorentz);
  }
  
  double LorentzianGaussian::lorentzval(double f) const {
    return ls_/(a2_+f*f); 
  }

  double LorentzianGaussian::gaussval(double f) const {
    return gs_*std::exp(-ga_*f*f); 
  }

void LineshapeSpec::validate() const
{
  if ((cutoff<0.0) || (cutoff>=1.0))
    throw InvalidParameter("LineshapeSpec: cutoff must be between 0 and 1");
}

  BaseLineshape::BaseLineshape(const LineshapeSpec& spec)
    : spec_(spec)
{
  spec.validate(); //!< need to verify here since LineshapeSpec is a simple struct
  reset_counts();
}

  void LorentzianGaussian::cache_t::ensure(const LorentzianGaussian& obj, hash_t key)
  {
    double offset;
    size_t bin;
    obj.getbin(bin,offset,key);
    if ((obj.frac_!=PureLorentzian) && gauss.empty())
      obj.makeg_(gauss,bin,offset);
    if ((obj.frac_!=PureGaussian) && lorentz.empty())
      obj.makel_(lorentz,bin,offset);
  }

void LorentzianGaussian::getbin(size_t& umid, double& offset, hash_t key) const
{
  int mid;
  reversehash(offset,mid,key);
  mid+=getmid();
  const size_t nbins=speccache_.nbins;
  assert((mid>=0) && (mid<nbins));
  umid=mid;
  offset*=speccache_.rangef/nbins;
}

  void LorentzianGaussian::cache_t::operator()(List<double>& dest, double frac) const
  {
    if (frac==PureGaussian)
      dest=gauss;
    else {
      dest=lorentz;
      if (frac!=PureLorentzian) {
	dest*=frac;
	mla(dest,1.0-frac,gauss);
      }
    }
  }

#define CALL_MEMBER_FN(obj,fptr) ((obj).*(fptr))

  size_t LorentzianGaussian::make_(List<double>& dest,size_t mid, double offset, const char* name, LineshapeFunc_t func,double zeroval) const
  {
    const double cut=zeroval*spec_.cutoff;
    const size_t nbins=speccache_.nbins;
    const double deltaf=speccache_.rangef/nbins;
    dest.create(nbins);
    if (name && (spec_.flags & LineshapeSpec::verbose))
      std::cout << "Creating " << name << " with width " << lw_ << " centered on bin " << mid << " with frequency offset " << offset << '\n';

    if (spec_.cutoff)
      dest=0.0; //!< only zero if no cutoff
    int width=0;

    if (offset==0.0) {//!< symmetrical
      dest(size_t(mid))=zeroval;
      for (width=1;;width++) {
	offset+=deltaf;
	const double v=CALL_MEMBER_FN(*this,func)(offset);
	if (v<cut)
	  break;
	const size_t upperbin=mid+width;
	bool failed=true;
	if (upperbin<nbins) {
	  dest(upperbin)=v;
	  failed=false;
	}
	const int lowerbin=mid-width;
	if (lowerbin>=0) {
	  dest(size_t(lowerbin))=v;
	  failed=false;
	}
	if (failed) //!< break out if finished
	  break; 
      }
    }
    else { //!< asymmetrical
      double loffset=-offset;
      double roffset=deltaf-offset;
      for (width=0;;width++) {       //! could be a bit smarter and not calculate one side when out of range, but this wouldn't gain much
	const double lv=CALL_MEMBER_FN(*this,func)(loffset);
	const double rv=CALL_MEMBER_FN(*this,func)(roffset);
	if ((lv<cut) && (rv<cut))
	  break;
	const size_t upperbin=mid+1+width;
	bool failed=true;
	if (upperbin<nbins) {
	  dest(upperbin)=rv;
	  failed=false;
	}
	const int lowerbin=mid-width;
	if (lowerbin>=0) {
	  dest(size_t(lowerbin))=lv;
	  failed=false;
	}
	if (failed) //!< break out if finished
	  break; 
	loffset-=deltaf;
	roffset+=deltaf;
      }
    }
    return width;
  }

  void LorentzianGaussian::makel_(List<double>& dest, size_t mid, double offset) const
  {
    const size_t width=make_(dest,mid,offset,"Lorentzian",&LorentzianGaussian::lorentzval,lorentz0());
    if (width>Lwidth_)
      Lwidth_=width;
  }

  void LorentzianGaussian::makeg_(List<double>& dest, size_t mid, double offset) const
  {
    const size_t width=make_(dest,mid,offset,"Gaussian",&LorentzianGaussian::gaussval,gauss0());
    if (width>Gwidth_)
      Gwidth_=width;
  }

  size_t LorentzianGaussian::gethalfwidth() const
  {
    if (frac_==PureLorentzian)
      return Lwidth_;
    if (frac_==PureGaussian)
      return Gwidth_;
    if ((Lwidth_==0U) || (Gwidth_==0U))
      return 0U;
    return std::max(Lwidth_,Gwidth_);
  }
  
  void LorentzianGaussian::reset()
  {
    BaseLineshape::reset_counts();
    //! clear cache *contents* but not memory itself
    const cachemap_t::iterator end(cachemap_.end());
    cachemap_t::iterator iter(cachemap_.begin());
    while (iter!=end) {
      iter->second.reset();
      ++iter;
    }
    cachet_.create(0U);
    Lwidth_=Gwidth_=0U;
  }
  
  void BaseLineshape::flags(int flagsv)
  { 
    int diff=spec_.flags ^ flagsv;
    static const int trivialmask=LineshapeSpec::verbose;
    if (diff) {
      spec_.flags=flagsv;      
      if (diff & ~trivialmask)
	reset();
    }
  }

void BaseLineshape::reset_counts()
{
  callcount_=hitcount_=0;
}
  
  void LorentzianGaussian::print(std::ostream& ostr) const
  {
    ostr << "Linewidth: " << lw_ << " Hz   G/L fraction: " << (100.0*frac_) << "%\n";
    ostr << "Lorentzian half-width: " << Lwidth_ << " bins   Gaussian half-width: " << Gwidth_ << " bins\n";
    ostr << "Cached lineshapes: " << cachemap_.size() << '\n';
    ostr << "Call count: " << callcount_ << "  Cache usage: " << hitcount_ << " (" << ((100.0*hitcount_)/callcount_) << "%)\n";
  }

  int LorentzianGaussian::operator()(List<double>& dest, const HistogramSpec& histspec, double freq) const
  {
    if (lw_==0.0)
      throw Failed("LorentzianGaussian: lineshape requested before linewidth set");
    //    if (histspec.nbins!=dest.size())
    //  throw Mismatch("LorentzianGaussian::make");

    if (histspec!=speccache_) {
      LorentzianGaussian* nonconst_this=const_cast< LorentzianGaussian*>(this);
      nonconst_this->reset(); //!< clear cache if histogram specification has changed
      speccache_=histspec;
    }
         
    const double ifreq=freq;
    
    const bool nofolding=spec_.flags & LineshapeSpec::nofold;
    const bool isverbose=spec_.flags & LineshapeSpec::verbose;
    const bool allowvar=spec_.flags & LineshapeSpec::variable;
    const size_t nbins=histspec.nbins;
    if (!nofolding)
      freq=mod(freq-histspec.minf,histspec.rangef);
    
    int bin;
    double offset=0.0;
    const double deltaf=histspec.rangef/nbins;
    const double scalef=1.0/deltaf;
	
    if (spec_.resolution_steps<2)
      bin = nofolding ? int(scalef*(freq-histspec.minf)) : int(scalef*freq);
    else {
      if (nofolding) {
	const double relfreq=freq-histspec.minf;
	if (relfreq<0) //!< always round down
	  bin=-int(scalef*-relfreq)-1;
	else
	  bin=int(scalef*relfreq);
	offset=(freq-histspec.minf)/deltaf-bin;
      }
      else {
	bin=int(scalef*freq);
	offset=freq/deltaf-bin;
      }
      //!< catch rounding problems
      if (offset>=1.0) {
	offset-=1.0;
	bin++;
      }
    }
    const int mid=getmid();
    if (isverbose)
      std::cout << "Initial frequency " << ifreq << " reduced to bin " << bin << ", fractional offset " << offset << '\n';
    assert((offset>=0.0) && (offset<=1.0));
    if ((bin<0) || (bin>=nbins)) {
      if (isverbose)
	std::cerr << "Warning: frequency outside histogram range (ignored)\n";
      if (allowvar)
	dest.create(0U);
      else
	dest.create(nbins,0.0);
      return -1;
    }
    int binkey=0;
    int roll=bin-mid;
    int halfwidth=gethalfwidth();
    if (nofolding && ((halfwidth==0U) || (spec_.cutoff==0) || (std::abs(binkey)+halfwidth>=mid))) {
      roll=0;
      binkey=bin-mid;
    }
    
    const hash_t hashkey(hash(offset,binkey));
    if (isverbose) {
      std::cout << "Post-folding binkey " << binkey << ", roll " << roll << '\n';
      std::cout << "Hashed to " << hashkey << '\n';
      double noff;
      int nbin;
      reversehash(noff,nbin,hashkey);
      const double offshift=spec_.resolution_steps*fabs(noff-offset);
      if ((offshift>1.0) || (nbin!=binkey)) {
	std::cerr << "Reverse hash failed: binkey " << nbin << ", offset " << offset << std::endl;
	throw InternalError("LorentzianGaussian");
      }
    }
    callcount_++;
    cachemap_t::iterator iter(cachemap_.find(hashkey));
    cache_t* sourcep=&tmp_;
    if (iter==cachemap_.end()) {
      tmp_.make(*this,hashkey);
      if (cachemap_.size()<spec_.cache_maximum) {
	cache_t& curentry(cachemap_[hashkey]);
	curentry.swap(tmp_);
	sourcep=&curentry;
	if (isverbose)
	  std::cout << "Creating cached lineshape for hash: " << hashkey << '\n';
      }
      else {
	if (isverbose)
	  std::cout << "Cache full\n";
      }
    }
    else {
      sourcep=&(iter->second);
      sourcep->ensure(*this,hashkey);
      if (isverbose)
	std::cout << "Using cached lineshape for hash: " << hashkey << '\n';
      hitcount_++;
    }

    (*sourcep)(cachet_,frac_); //!< create lineshape
    
    halfwidth=gethalfwidth(); //!< re-read as it might have been updated
    if (halfwidth==0)
      throw InternalError("LorentzianGaussian");

    const int origin=mid+roll;
    const int left=origin-halfwidth;
    const int right=origin+halfwidth;
    const bool doesfit=(right<nbins) && (left>=0);
    const BaseList<double> sourceshape(cachet_(range(mid-halfwidth,mid+halfwidth)));

    if (!allowvar || !doesfit) {
      if (roll==0)
	dest=cachet_;
      else {
	dest.create(nbins);
	if (isverbose) 
	  std::cout << "Making lineshape with halfwidth " << halfwidth << " bins and roll factor " << roll << '\n';
	if (doesfit) {//!< shape fits entirely 
	  dest=0.0;
	  BaseList<double> destsub(dest(range(left,right)));
	  destsub=sourceshape;
	}
	else {
	  if (roll<0)
	    roll+=nbins;
	  if ((roll<0) || (roll>=nbins))
	    throw InternalError("LorentzianGaussian(2)");
	  BaseList<double> dest1(dest(range(roll,nbins-1)));
	  dest1=cachet_(range(0U,nbins-roll-1));
	  BaseList<double> dest2(dest(range(0U,roll-1)));
	  dest2=cachet_(range(nbins-roll,nbins-1));
	}
      }
      return 0;
    }
    //variable width
    if (isverbose) 
      std::cout << "Making (restricted) lineshape with halfwidth " << halfwidth << " bins\n";
    dest=sourceshape;
    return left;
  }  
  
  LorentzianGaussian::LorentzianGaussian(double lw, double frac, const LineshapeSpec& specv)
    : BaseLineshape(specv), Lwidth_(0U), Gwidth_(0U)
  {
    setlw_(lw);
    setfraction_(frac);
  }
  
  void LorentzianGaussian::lw(double nlw)
  {
    if (nlw!=lw_) {
      setlw_(nlw);
      reset();
    }
  }
  
  void LorentzianGaussian::fraction(double nfrac)
  {
    if (nfrac!=frac_) {
      setfraction_(nfrac);
      cachet_.create(0U); //!< only invalidate total cache
    }
  }
  
  void LorentzianGaussian::setlw_(double nlw)
  {
    if (nlw<=0.0)
      throw InvalidParameter("LorentzianGaussian: linewidth cannot be <=0.0");
    lw_=nlw;
    a_=lw_/2.0;
    a2_=a_*a_;
    ls_=a_/M_PI;
    static const double ln2=std::log(2.0);
    ga_=ln2/a2_;
    gs_=std::sqrt(ga_/M_PI);    
  }
  
  void LorentzianGaussian::setfraction_(double nfrac)
  {
    if ((nfrac<0.0) || (nfrac>1.0))
      throw InvalidParameter("LorentzianGaussian: L/G fraction must be between 0 and 1");
    frac_=nfrac;
  }

} //namespace libcmatrix
