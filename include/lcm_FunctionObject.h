#ifndef _FunctionObject_h_
#define _FunctionObject_h_

//#include "ScratchList.h"
#include <ostream>

namespace libcmatrix {

  template<class T,class T2> struct UnaryFunction : public ::std::unary_function<T2,T> {
    virtual ~UnaryFunction() {};
    
    virtual void operator() (T& d,T2 x) const { d=(*this)(x); }
    virtual T operator() (T2 x) const =0;
    //    virtual UnaryFunction* clone() const =0;
  };

/*  template<class T,class T2> class ConstUnaryFunction : public UnaryFunction<T,T2> { */
/*    const T val; */
/*  public: */
/*    ConstUnaryFunction(const T& val_) : val(val_) {}; */
/*    T operator()(T2) const { return val; } */
/*  }; */

//   template<class T> struct NullUnaryFunction : public UnaryFunction<T,T> {
//     T operator()(T x) { return x; }
//     //    UnaryFunction* clone() const { return new NullUnaryFunction<T>(); }
//   };
    
template<class T,class T2,class T3> class BinaryFunction : public ::std::binary_function<T2,T3,T> {
public:
  virtual ~BinaryFunction() {};
  
  virtual void operator() (T& d,T2 x1,T3 x2) const { d=(*this)(x1,x2); }
  virtual T operator() (T2,T3) const =0;
  //  virtual BinaryFunction* clone() const =0;
};

/*  template<class T,class T2,class T3> class ConstBinaryFunction : public BinaryFunction<T,T2,T3> { */
/*    T val; */
/*  public: */
/*    ConstBinaryFunction(const T& val_) : val(val_) {}; */
/*    T operator()(T2,T3) const { return val; } */
/*  }; */

 //Fourier series...
/* class RealFourierFunc : public UnaryFunction<double,double> { */
/*   const ScratchList<complex> coeffs; */
/*   double mod_freq; */

/* public: */
/*   RealFourierFunc(const BaseList<complex> &coeffs_,double tcycle_) : coeffs(coeffs_) { */
/*     if (tcycle_<=0.0) throw InvalidParameter("RealFourierFunc: cycle time must be >0"); */
/*     mod_freq=1.0/tcycle_; */
/*   } */
  
/*   double operator()(double t) const { */
/*     const complex fac=expi(-2*M_PI*mod_freq*t); */
/*     complex f=fac; */
/*     double ph=real(coeffs(0)); */

/*     const size_t n=coeffs.length(); */
/*     for (size_t i=1;i<n;i++) { */
/*       real_mla(ph,coeffs(i),f); */
/*       f*=fac; */
/*     } */
/*     return ph; */
/*   } */
/* }; */


//Regularly spaced table lookup (Failed thrown if outside range 0 <= t <= tcycle) 
/*  class ContinuousIndex : public UnaryFunction<size_t,double> { */
/*    double mfac,tol,tcycle; */
/*    size_t nbins; */
   
/*  public: */
/*    ContinuousIndex(int nbins_,double tcycle_) { */
/*       if (tcycle_<=0.0) throw InvalidParameter("ContinuousIndex: cycle time must be >0"); */
/*       tcycle=tcycle_; */
/*       if (nbins_<1) throw InvalidParameter("ContinuousIndex: list empty"); */
/*       nbins=nbins_; */
/*       mfac=nbins/tcycle_; */
/*       tol=1e-8*tcycle_; */
/*    } */
/*    ContinuousIndex() { nbins=0; mfac=0.0; } */
/*    size_t operator()(double t) const { */
/*      int bin=int(t*mfac); */
/*      if (bin<0) { */
/*        if (t<-tol) throw Failed("TableFunc: out of range"); */
/*        return 0; */
/*      } */
/*      else { */
/*        if (bin>nbins-1) { */
/* 	 if (t>tcycle+tol) throw Failed("TableFunc: out of range"); */
/* 	 return nbins-1; */
/*        } */
/*        return bin; */
/*      } */
/*    } */
/*  }; */

/*  template<class L,class Index> class Table : public UnaryFunction<typename L::value_type,typename Index::argument_type> { */
/*    Index which; */
/*    const L& vals; */
/*  public: */
/*    Table(const L& vals_,Index which_) : vals(vals_), which(which_) {} */
/*    typename L::value_type operator()(typename Index::argument_type x) const { return vals(which(x)); } */
/*  }; */

/*  template<class T> struct ContinuousTable : public Table<T,ContinuousIndex> { */
/*    ContinuousTable(const T& vals_,double tcycle_) : Table<T,ContinuousIndex>(vals_,ContinuousIndex(vals_.length(),tcycle_)) {} */
/*  }; */

/*  //Cyclic table lookup */
/*  class ContinuousCyclicIndex : public UnaryFunction<size_t,double> { */
/*    double sfac; */
/*    size_t nbins; */

/*  public: */
/*    ContinuousCyclicIndex(int nbins_,double tcycle_) { */
/*       if (tcycle_<=0.0) throw InvalidParameter("ContinuousCyclicIndex: cycle time must be >0"); */
/*       if (nbins_<1) throw InvalidParameter("ContinuousCyclicIndex: list empty"); */
/*       nbins=nbins_; */
/*       sfac=1.0/tcycle_; */
/*    } */
/*    ContinuousCyclicIndex() { nbins=0; sfac=0.0; } */
/*    size_t operator()(double t) const { */
/*      t*=sfac; */
/*      int bin=int(nbins*(t-floor(t))); */
/*      return (bin==nbins) ? 0 : bin; */
/*    } */
/*  }; */

/*  template<class T> struct ContinuousCyclicTable : public Table<T,ContinuousCyclicIndex> { */
/*    ContinuousCyclicTable(const T& vals_,double tcycle_) : Table<T,ContinuousCyclicIndex>(vals_,ContinuousCyclicIndex(vals_.length(),tcycle_)) {} */
/*  }; */

/*  template<class T> class CyclicTable : public UnaryFunction<typename T::value_type,size_t> { */
/*    const size_t n; */
/*    const T& vals; */
/*  public: */
/*    CyclicTable(const T& vals_) : n(vals_.length()), vals(vals_) {} */
/*    T operator()(size_t m) const { return vals(m % n); } */
/*  }; */

 //Table lookup with irregular spacings; only efficient for sequential access
/* class IrregularIndex : public UnaryFunction<size_t,double> { */

/*   mutable size_t lastindex; */
/*   double tol; */
/*   ScratchList<double> starts; */
/* protected: */
/*   double tcycle; */

/* public: */
/*   IrregularIndex(const BaseList<double>& durs_) : starts(durs_.length()+1) { */
/*     const size_t n=durs_.length(); */
/*     if (n==0) throw Failed("IrregularIndex"); */
/*     starts(0)=0.0; */
/*     for (size_t i=0;i<n;i++) { */
/*       if (durs_(i)<=0.0) throw Failed("IrregularIndex: duration is <=0.0"); */
/*       starts(i+1)=starts(i)+durs_(i); */
/*     } */
/*     tcycle=starts(n); */
/*     tol=1e-8*tcycle; */
/*     lastindex=0; */
/*   } */

/*   size_t operator()(double x) const { */
/*     if (x>=starts(lastindex+1)) { */
/*       if (x>tcycle) { */
/* 	if (x>tcycle+tol) throw Failed("IrregularIndex"); */
/* 	return (lastindex=starts.length()-2); */
/*       } */
/*       while (x>=starts(lastindex+1)) lastindex++; */
/*       return lastindex; */
/*     } */
/*     if (x<=0.0) { */
/*       if (x<-tol) throw Failed("IrregularIndex"); */
/*       return (lastindex=0); */
/*     }  */
/*     while (x<starts(lastindex)) lastindex--; */
/*     return lastindex; */
/*   } */
/* }; */

/* template<class T> class IrregularTable : public UnaryFunction<T,double> { */
/*   IrregularIndex which; */
/*   const BaseList<T>& vals; */

/* public: */
/*   IrregularTable(const BaseList<T>& vals_, const BaseList<double>& durs_) : vals(vals_), which(durs_) { */
/*     if (durs_.length()!=vals_.length()) throw Mismatch("IrregularTable"); */
/*   } */
/*   T operator()(double t) const { return vals(which(t)); } */
/* }; */
  
//Cyclic table lookup with irregular spacings
/* class IrregularCyclicIndex : public IrregularIndex { */
/*   const double sfac; */
/* public: */
/*   IrregularCyclicIndex(const BaseList<double>& durs_) : IrregularIndex(durs_), sfac(1.0/tcycle) {} */

/*   size_t operator()(double x) const { */
/*     if (x>=tcycle || x<=0.0) { */
/*       double t=x*sfac; */
/*       x=tcycle*(t-floor(t)); */
/*     } */
/*     return IrregularIndex::operator()(x); */
/*   } */
/* }; */

/*  template<class T> struct IrregularCyclicTable : public Table<T,IrregularCyclicIndex> { */
/*    IrregularCyclicTable(T& vals_, const BaseList<double>& durs_) : Table<T,IrregularCyclicIndex>(vals_,IrregularCyclicIndex(durs_)) { */
/*      if (durs_.length()!=vals_.length()) throw Mismatch("IrregularCyclicFunc"); */
/*    } */
/*  }; */
 
} //namespace libcmatrix

#endif
