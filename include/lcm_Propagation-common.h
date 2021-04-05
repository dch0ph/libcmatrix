//Not to be included directly
#ifndef Propagation_common_h
#define Propagation_common_h

#include "BasePropagation.h"
#include "cmatrix.h"
#include "rmatrix.h"
#include "UnionHolder.h"
#include "Warnings.h"

namespace libcmatrix {

 const double offdiagonal_tolerance=1e-6;

  extern Warning<InvalidParameter> propagation_closetodiagonal_warning;

  typedef RealComplexHolder<Matrix<double> ,Matrix<complex> > RCmatrix;

/*   template<class M> class RotatingFrame; */
/*   template<class M> std::ostream& operator<< (std::ostream&, const RotatingFrame<M>&); */

/*   template<class M> class RotatingFrame : public SpectrumIterator<typename M::amplitude_type> { */
/*     M& obj_; */
/*     const double rate_; */
/*     const bool quadrature_; */
/*     double minrate,maxrate; */
/*   public: */
/*     RotatingFrame(M& objv,double ratev, bool quadraturev =true) */
/*       : obj_(objv), rate_(ratev), quadrature_(quadraturev) { */
/*       if (ratev>0) { */
/* 	minrate=-ratev/2; maxrate=ratev/2; */
/*       } */
/*       else { */
/* 	maxrate=-ratev/2; minrate=ratev/2; */
/*       } */
/*     } */

/*     bool operator() (typename M::amplitude_type& amp, double& freq) { */
/*       if (!rate_) */
/* 	return obj_(amp,freq); */
/*       double tfreq; */
/*       while (obj_(amp,tfreq)) { */
/* 	freq=tfreq-rate_; */
/* 	if (freq>minrate && freq<maxrate) return true; */
/* 	if (!quadrature_) { */
/* 	  freq=tfreq+rate_; */
/* 	  if (freq>minrate && freq<maxrate) return true; */
/* 	} */
/*       } */
/*       return false; */
/*     } */

/*     friend std::ostream& operator<< <>(std::ostream&,const RotatingFrame<M>&); */
/*   }; */

/*     template<class M> std::ostream& operator<< (std::ostream& ostr,const RotatingFrame<M>& obj) */
/*     { */
/*       ostr << obj.obj_; */
/*       if (obj.rate_) */
/* 	return ostr << "\nWith frame rotating at " << (obj.rate_/1e6) << " MHz\n"; */
/*       return ostr << "\nWith null rotating frame correction\n"; */
/*     } */

//very primitive!
 struct InhomoIter {

   InhomoIter() { lock(); }

   void reset(size_t dimRv,size_t dimCv,bool offdiagv);
   bool advance();
   bool isfinished() { return finished_; }
   void lock() { finished_=true; }

   friend std::ostream& operator<< (std::ostream&, const InhomoIter&);

   size_t dimR,dimC;
   size_t r,s;
   bool finished_;
   bool offdiagonly;
 };

struct InhomoStash : public RCmatrix {

  rmatrix eff_freqs;
  cmatrix phasestore;
  bool isvalid;

  friend std::ostream& operator << (std::ostream&, const InhomoStash&);

  InhomoStash() { isvalid=false; }

  size_t interactions() const { return eff_freqs.rows(); }
  bool explicit_Us() const //true if propagators have been set explicitly
  { return eff_freqs.empty() && !phasestore.empty(); }
  size_t size() const {
    if (eff_freqs.empty()) {
      if (phasestore.empty())
	throw Failed("Hamiltonian not set");
      return phasestore.rows();
    }
    return eff_freqs.cols();
  }
  bool operator!() const
  { return eff_freqs.empty() && phasestore.empty(); } //this overwrites corresponding operator in RCmatrix
  void clear() { RCmatrix::clear(); eff_freqs.clear(); }
  //Either eff_freqs or phasestore will be set if object initialised
  //clear must also be overridden!

  bool iselement() const { 
    if (!(*this)) throw Failed("InhomoStash::iselement: Hamiltonian not set");
    return (type()==NONE);
  }
  void updatephases(const rmatrix&,bool);
  void invalidatephases() { isvalid=false; }
  //  void write(FILE*, const char*) const;
  //int read(FILE*);

  void set_H(const cmatrix&);
  void set_H(const rmatrix&);
  void set_H(const BaseList<double>&);
  void set_H(const cmatrix&,const BaseList<double>&);
  void set_H(const rmatrix&,const BaseList<double>&);
  void set_H(double);

  void set_U(const cmatrix& U, double period);
  void set_U(const complex& U, double period);
  void set_U(const BaseList<complex>& U, double period);

  void set_Us(const cmatrix&, double period);
};

struct InhomoHelper_ {
  template<class M> static void observe(M&, cmatrix&, const cmatrix&, const cmatrix&);
  template<class M> static void observe(M&, cmatrix&, const rmatrix&, const rmatrix&);
  template<class M> static void observe(M&, cmatrix&, const BaseList<double>&, const BaseList<double>&);
  template<class M> static void observe(M&, cmatrix&, const cmatrix&);
  template<class M> static void observe(M&, rmatrix&, const cmatrix&);
  template<class M> static void observe(M&, rmatrix&, const rmatrix&);
  template<class M> static void observe(M&, cmatrix&, const rmatrix&);
  template<class M> static void observe(M&, rmatrix&, const BaseList<double>&);
  template<class M> static void observe(M&, cmatrix&, const BaseList<double>&);
  template<class M> static void observe(M&, rmatrix&);
  template<class M> static void observe(M&, cmatrix&);

  template<class M> static void process(M&, rmatrix&, const cmatrix&);
  template<class M> static void process(M&, cmatrix&, const cmatrix&, const cmatrix&);
  template<class M> static void process(M&, cmatrix&, const rmatrix&, const rmatrix&);
  template<class M> static void process(M&, cmatrix&, const cmatrix&);
  template<class M> static void process(M&, cmatrix&, const rmatrix&);
  template<class M> static void process(M&, rmatrix&);
  template<class M> static void process(M&, cmatrix&);
};
  
  class InhomoTDBase_ {
    public:

   template<class T> List<complex> FID(int npoints,const Matrix<T>& sigma0,const Matrix<T>& detect) {
     List<complex> d(npoints,complex(0.0),mxflag::temporary);
     this->add_FID(d,1.0,sigma0,detect);
     return d;
   }

   template<class T> List<complex> FID(int npoints,const Matrix<T>& sigma0det) {
     List<complex> d(npoints,complex(0.0),mxflag::temporary);
     this->add_FID(d,1.0,sigma0det);
     return d;
   }

   template<class T> List<double> FID_hermitian(int npoints,const Matrix<T>& sigma0,const Matrix<T>& detect) {
     List<double> d(npoints,0.0,mxflag::temporary);
     this->add_FID_hermitian(d,1.0,sigma0,detect);
     return d;
   }  

   template<class T> List<double> FID_hermitian(int npoints,const Matrix<T>& a) { 
     List<double> d(npoints,0.0,mxflag::temporary);
     this->add_FID_hermitian(d,1.0,a);
     return d;
   }  

   List<double> FID(int npoints,const BaseList<double>& sigma0,const BaseList<double>& detect) {
     List<double> d(npoints,0.0,mxflag::temporary);
     this->add_FID(d,1.0,sigma0,detect);
     return d;
   }

   List<double> FID(int npoints,const BaseList<double>& sigma0det) {
     List<double> d(npoints,0.0,mxflag::temporary);
     this->add_FID(d,1.0,sigma0det);
     return d;
   }

   List<complex> FID(int npoints) {
     List<complex> d(npoints,complex(0.0),mxflag::temporary);
     this->add_FID(d); return d;
   }
     
   void add_FID(BaseList<double> FID_, double scale, const BaseList<double>& sigma0det) { 
     observe(sigma0det);
     add_FID_hermitian_(FID_,scale);
   }
   
   void add_FID(BaseList<double> FID_, double scale, const BaseList<double>& sigma0, const BaseList<double>& detect) {
     observe(sigma0,detect);
     add_FID_hermitian_(FID_,scale);
   }

   void add_FID(BaseList<complex> FID_, complex scale, const BaseList<double>& sigma0det) {
     observe(sigma0det);
     add_FID_(FID_,scale); 
   }

   void add_FID(BaseList<complex> FID_, complex scale, const BaseList<double>& sigma0, const BaseList<double>& detect) {
     observe(sigma0,detect);
     add_FID_(FID_,scale);
   }

   template<class T> void add_FID(BaseList<complex> FID_, complex scale, const Matrix<T>& sigma0det) {
     observe(sigma0det);
     add_FID_(FID_,scale);
   }

   template<class T> void add_FID(BaseList<complex> FID_, complex scale, const Matrix<T>& sigma0, const Matrix<T>& detect) {
     observe(sigma0,detect);
     add_FID_(FID_,scale);
   }

    void add_FID(BaseList<complex> FID_, complex scale =complex(1.0)) {
     observe();
     add_FID_(FID_,scale);
   }

   template<class T,class M> void add_FID_hermitian(BaseList<M> FID_, M scale, const Matrix<T>& sigma0det) { 
     observe(sigma0det);
     add_FID_hermitian_(FID_,scale);
   }

   template<class T,class M> void add_FID_hermitian(BaseList<M> FID_, M scale, const Matrix<T>& sigma0, const Matrix<T>& detect) {
     observe(sigma0,detect); 
     add_FID_hermitian_(FID_,scale);
   }

  protected:
    virtual void observe() =0;
    virtual void observe(const rmatrix&) =0;
    virtual void observe(const cmatrix&) =0;
    virtual void observe(const BaseList<double>&) =0;
    virtual void observe(const rmatrix&, const rmatrix&) =0;
    virtual void observe(const cmatrix&, const cmatrix&) =0;
    virtual void observe(const BaseList<double>&, const BaseList<double>&) =0;
    
    virtual void add_FID_(BaseList<complex>, complex) const =0;
    virtual void add_FID_hermitian_(BaseList<double>, double) const =0;
    virtual void add_FID_hermitian_(BaseList<complex>, complex) const =0;
  };

 class PeriodicObj {
 public:
   PeriodicObj(size_t nobsv =0, int verbosev =0)
     : nobs(nobsv), verbose_(verbosev), isgamma(nobsv!=0), period_(0.0), speed_(0), gamma_steps(0) {}

  double period() const { return period_; }
  double speed() const {
    if (!period_)
      throw Failed("PeriodicObj: period unset");
    return speed_; }

  size_t observations() const { return nobs; }
  size_t gammasteps() const { return gamma_steps; }

  friend std::ostream& operator<< (std::ostream&, const PeriodicObj&);

 protected:
   bool set_steps(int);

  void period(double periodv) { 
    if (periodv<=0.0)
      throw InvalidParameter("PeriodicObj: period cannot be <0");
    period_=periodv;
    speed_=1.0/period_;
  }
  
  size_t nobs;
   size_t verbose_;
  bool isgamma;
  double period_;
  double speed_;
  size_t gamma_steps;
  //  size_t intsteps;
 };

} //namespace libcmatrix
#endif
