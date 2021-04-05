#ifndef lcm_AsynchronousPropagation_h_
#define lcm_AsynchronousPropagation_h_

#include "BasePropagation.h"
#include "MultiMatrix.h"
#include "cmatrix.h"
//! for definition of PropGen_t
#include "NMR.h" 
#include "smartptr.h"

namespace libcmatrix {

//   struct StashPoint : public MultiMatrix<complex,3> {
//     size_t size() const {
//       if (empty())
// 	throw Failed("StashPoint: propagators not set");
//       return dimension(1);
//     }
//   };
  struct StashPoint {
    void set(const PropGen_t& adaptorv, double inittv, double dtv) { 
      adaptorp.reset(adaptorv.clone());
      initt_=inittv;
      dt_=dtv;
    }
    size_t size() const { return adaptorp->size(); }
    double starttime(size_t n) const { return initt_+n*dt_; }

    void operator()(cmatrix& U, size_t n) const
    {
      const double t1=starttime(n);
      (*adaptorp)(U,t1,t1+dt_); //!< smartptr will throw exception if unset
    }

    void clear() { adaptorp.clear(); }
    smartptr<const PropGen_t> adaptorp;
    double initt_,dt_;
  };

  class AsynchronousFID : public SwapStore<StashPoint> {
 public:
    AsynchronousFID(double inittv, double dtv, int verbosev =0)
      : initt_(inittv), dt_(dtv), verbose_(verbosev) {}

    void set_Ugenerator(char sel,const PropGen_t& Ugen) {
      setRC(sel).set(Ugen,initt_,dt_);
    }
    void set_Ugenerator(const PropGen_t& Ugen) {
      setdiagonal().set(Ugen,initt_,dt_);
    }

//   void set_Us(char sel,const BaseList<cmatrix>& Us) {
//     MultiMatrix<complex,3>& where=setRC(sel);
//     where=Us;
//   }
//   void set_Us(const BaseList<cmatrix>& Us) {
//     MultiMatrix<complex,3>& where=setdiagonal();
//     where=Us;
//   }
//   size_t nincs() {
//     const size_t nincs_=row().dimension(0);
//     if (nincs_!=col().dimension(0))
//       throw Mismatch("AsynchronousFID propagator lists");
//     return nincs_;
//   }
    List<complex> FID(int npoints) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d,complex(1.0)); return d; }
    List<complex> FID(int npoints,const cmatrix& sigma0,const cmatrix& detect) { List<complex> d(npoints,complex(0.0),mxflag::temporary); this->add_FID(d,complex(1.0),sigma0,detect); return d; }
   List<double> FID_hermitian(int npoints,const cmatrix& sigma0,const cmatrix& detect) { List<double> d(npoints,0.0,mxflag::temporary); this->add_FID_hermitian(d,1.0,sigma0,detect); return d; }  

  void add_FID(BaseList<complex>, complex scale);
  void add_FID(BaseList<complex>, complex scale, const cmatrix&, const cmatrix&);
  template<class T> void add_FID_hermitian(BaseList<T> FID, T, const cmatrix&, const cmatrix&);

 private:
    double initt_,dt_;
    int verbose_;
    void propagateRC(size_t n);
    double starttime(size_t n) const { return initt_+n*dt_; }
    cmatrix sigma,tmp1,tmp2;

};

template<> struct SignalGenerator_traits<AsynchronousFID> {
  static const bool timedomain=true;
  static const bool allowreal=false;
  static const bool allowED=false;
  static const bool allowdiagonal=false;
  static const size_t input_dimensionality=3;
  static const bool allowH=false;
  static const bool allowU=false;
  static const bool usegenerator=true;
  static const bool gamma=false;
};

} //namespace libcmatrix

#endif
