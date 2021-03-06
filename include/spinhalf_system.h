#ifndef _spinhalf_system_h_
#define _spinhalf_system_h_

#include "cmatrix.h"
#include "List.h"
#include "basespin_system.h"
#include <iostream>

namespace libcmatrix {

class spinhalf_system : public basespin_system {
public:
  //explicit spinhalf_system(int n) : basespin_system(n,0.5,true) { docreate(n); }
  spinhalf_system(int n,const char *def ="1H") : basespin_system(n,spin(def)) { docreate(n); }
  //spinhalf_system(int n,const BaseList<state_t> &l) : basespin_system(n,0.5,true) { docreate(n,l); }
  //spinhalf_system(int n,const BaseList<state_t> &rl,const BaseList<state_t> &cl) : basespin_system(n,0.5,true) { docreate(n,rl,cl); }
  spinhalf_system(int n,const char* def, const BaseList<state_t>& l) : basespin_system(n,spin(def)) { docreate(n,l); }
  spinhalf_system(int n,const char* def, const BaseList<state_t>& rl, const BaseList<state_t>& cl) : basespin_system(n,spin(def)) { docreate(n,rl,cl); }
  spinhalf_system(int n,const char* def, float_t);
  spinhalf_system(int n,const char* def, float_t, float_t);
  spinhalf_system(const basespin_system&, size_t);

  basespin_system* clone() const { return new spinhalf_system(*this); }
  spinhalf_system* clone(size_t n) const { return new spinhalf_system(*this,n); }
  //  basespin_system* make(int n, const char* def) const { return new spinhalf_system(n,def); }

  size_t rows() const { return _rows; }
  size_t cols() const { return cstates.length(); }

  state_t mask(size_t n) const { if (n>=_nspins) throw BadIndex("mask"); else return state_t(1)<<(_nspins-n-1); }

  friend std::ostream& operator << (std::ostream &,const spinhalf_system &);

  BaseList<state_t> brastates() const { return isdiagonal() ? cstates : rstates; }
  BaseList<state_t> ketstates() const { return cstates; }

  List<state_t> Lstates() const;
  List<state_t> mzstates(const char *,float_t) const;

  bool iscomplete() const { return iscomplete_; }
  bool isdiagonal() const { return (rstates.length()==0); }

  void expandsuperop(Matrix<float_t>&, size_t, const Matrix<float_t>&) const;
  void expandsuperop(cmatrix&, size_t, const cmatrix&) const;

  void mla_I(cmatrix &,double,size_t,char) const;
  void mla_I(cmatrix &,double,size_t,char,size_t,char) const;
  void mla_I(cmatrix &,double,const BaseList<char> &) const;
  void mla_I(cmatrix &,double,size_t,const BaseList<char>&) const { throw Failed("not for spin-1/2 only"); }

  void permutation_vectorH(BaseList<state_t>&, const Permutation&) const;

  void nucleus(size_t,const spin&);
  bool isspinhalfonly() const { return true; }

  static state_t apply_permutation(state_t orig, const Permutation& permvec);

 private:

  List<state_t> cstates;
  List<state_t> rstates;
  List<int> reverse;

  bool iscomplete_;
  size_t _nspins;
  size_t _rows;

  void mla_Ip(cmatrix &,double,state_t) const;
  void mla_Im(cmatrix &,double,state_t) const;
  void mla_Iy(cmatrix &,double,state_t) const;
  void mla_Ix(cmatrix &,double,state_t) const;
  void rawmla_Iz(cmatrix &,double,state_t) const;

  void ensurematrix(cmatrix &) const;
  void createstates(List<state_t> &,const BaseList<size_t> &);

  void ensurereverse();

  void docreate(size_t);
  void docreate(size_t,const BaseList<state_t> &);
  void docreate(size_t,const BaseList<state_t> &,const BaseList<state_t> &);

  void print(std::ostream& = std::cout) const;
    
  state_t rawpermute(state_t orig, const Permutation& perm) const {
    if (perm.size()!=nspins())
      throw Mismatch("rawpermute");
    return apply_permutation(orig, perm);
  }

  void rawmla_Iz(BaseList<float_t> &,double,size_t) const;
  void rawmla_Iz(BaseList<float_t> &,double,size_t,size_t) const;

  void dump() const;
};

  //return with bit corresponding to spin n (in system of nspins) set to bit (should be 0/1)
inline state_t
tostate(size_t bit,size_t nspins,size_t n) {
  return state_t(bit) << (nspins-n-1);
}

  //return bit mask corresponding to spin n (in system of nspins) - equivalent to tostate(1,nspins,n)
inline state_t
maskelement(size_t nspins,size_t n) {
  return state_t(1) << (nspins-n-1);
}

//! count non-zero bits in state i.e. number of beta states
 size_t countnonzero(state_t);
 
  //return Inz for given state (system of nspins) 
inline double
Izelement(int nspins,state_t state,size_t n) {
  return ( state & maskelement(nspins,n) ) ? -0.5 : 0.5;
}

  //return Inz.Imz for given state (system of nspins)
inline double
Izelement(int nspins,state_t state,size_t n,size_t m) {
  return ( (!(state & maskelement(nspins,n))) ^ (!(state & maskelement(nspins,m))) ) ? -0.25 : 0.25;
}

  /* 'raw' functions can be used if maskelement is precalculated i.e. Izelement_raw(state & mark) */

  //translate between spin state and Iz (zero -> +0.5, non-zero -> -0.5)
inline double
Izelement_raw(state_t selbit) {
  return (selbit ? -0.5 : 0.5);
}

  //translate between spin states and I1z.I2z (+1/4 if states same, -1/4 if different)
inline double
Izelement_raw(state_t selbit1,state_t selbit2) {
  return ( (!selbit1) ^ (!selbit2) ) ? -0.25 : 0.25;
}

double setval(char);
double notsetval(char);

complex Ielement(int,state_t,state_t,size_t,char);
complex Ielement(int,state_t,state_t,size_t,char,size_t,char);
complex Ielement(int,state_t,state_t,const BaseList<char> &);

bool isdiagop(char);

void printbin(std::ostream &,state_t,int);
void printbin(std::ostream &,const BaseList<state_t> &,int);

List<state_t> mzstates(size_t nspins,float_t mz);
inline state_t Lstate(size_t nspins,state_t bra,state_t ket) { return (bra<<nspins) | ket; }

 template<> struct spin_system_traits<spinhalf_system> {
   static const bool spinhalfonly=true;
 };

} //namespace libcmatrix

#endif
