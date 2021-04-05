#ifndef _spin_system_h_
#define _spin_system_h_

/* General spin system, creating full Hilbert space matrices by direct product */

#include "basespin_system.h"
#include "cmatrix.h"

namespace libcmatrix {
    
class spin_system : public basespin_system {
public:
  spin_system(int n, const char *label ="1H") : basespin_system(n,spin(label)) {}
  spin_system(int n, const spin& spinv) : basespin_system(n,spinv) {}
  spin_system(const basespin_system&, size_t);

  void mla_I(cmatrix &,double scale,const BaseList<char>&) const;
  void mla_I(cmatrix &,double scale,size_t,const BaseList<char>&) const;
  void mla_I(cmatrix &,double scale,size_t n,char op) const;
  void mla_I(cmatrix &,double scale,size_t,char,size_t,char) const;

  void mla_ST(cmatrix &,double scale,size_t n,size_t,size_t) const;
  void mla_ST(List<double>&,double scale,size_t n,size_t) const;

  void get_states(BaseList<size_t>&, state_t) const;
  state_t from_states(const BaseList<size_t>&) const;

  basespin_system* clone() const { return new spin_system(*this); }
  spin_system* clone(size_t n) const { return new spin_system(*this,n); }
  
  bool isspinhalfonly() const;

 private:
  void I(cmatrix &,size_t,char,double =1.0) const;
  void I(cmatrix &,size_t,size_t,size_t,double =1.0) const;
  void I(cmatrix &,size_t,char,size_t,char,double =1.0) const;

  void ST(cmatrix &,size_t,size_t,size_t,double =1.0) const; //!< single-transition operator generator
  void ST(List<double>&, size_t,size_t,double =1.0) const;

  state_t rawpermute(state_t, const Permutation&) const;

  void rawmla_Iz(BaseList<float_t> &,double scale,size_t) const;
  void rawmla_Iz(BaseList<float_t> &,double scale,size_t,size_t) const;

  void dump() const;
};

 template<> struct spin_system_traits<spin_system> {
   static const bool spinhalfonly=false;
 };

} //namespace libcmatrix

#endif

