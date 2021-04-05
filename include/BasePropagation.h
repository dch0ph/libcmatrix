#ifndef BasePropagation_h_
#define BasePropagation_h_

#include "BaseList.h"

namespace libcmatrix {

/* Base class for iterating over spectrum */
  template<class T> class SpectrumIterator {
  public:
    SpectrumIterator() {}
    virtual ~SpectrumIterator() {};
    
    typedef T amplitude_type;
    virtual bool operator () (T& amp,double& freq) =0;
  };

void add_FID(BaseList<complex> FID,complex,double f);
void add_FID(BaseList<double> FID,complex,double f);

inline bool isrow(char rc)
{
  switch (rc) {
  case 'r': case 'R': return true;
  case 'c': case 'C': return false;
  }
  throw InvalidParameter("Must be one of r,R,c,R");
}

/* Generic holder for row/column "objects" */

template<typename T> class SwapStore {
 private:
  T stash1;
  T stash2;
  bool rpoint1,cpoint1;
  bool dirty_;

 public:
  SwapStore()
    : rpoint1(true), cpoint1(false), dirty_(false) {}

  const T& row() const { return rpoint1 ? stash1 : stash2; }
  const T& col() const { return cpoint1 ? stash1 : stash2; }

  size_t rows() const { return row().size(); }
  size_t cols() const { return col().size(); }
  bool iselement() const { return (row().iselement() || col().iselement()); }
  bool operator! () const { return (!(row()) || !(col())); }

  void swap();
  void shuffle(char);

  bool isdiagonal() const { return (rpoint1==cpoint1); }

  const T& operator()(char rc) const {
    verifyRC();
    return (isrow(rc) ? row() : col());
  }

  void verifyRC() const { if (rpoint1==cpoint1) 
      throw Failed("Object being used for diagonal blocks"); 
  }

  void clear() { stash1.clear(); stash2.clear(); dirty_=true; }
  void clear(char which) { (*this)(which).clear(); dirty_=true; }
  void copy(char which) {
    if (isdiagonal())
      throw Failed("SwapStore::copy");
    switch (which) {
    case 'R': 
      row()=col();
      break;
    case 'C':
      col()=row();
      break;
    default:
      throw InvalidParameter("SwapStore::copy");
    }
  }

 protected:
  T& row() { return rpoint1 ? stash1 : stash2; }
  T& col() { return cpoint1 ? stash1 : stash2; }

  T& operator()(char rc) {
    verifyRC();
    return (isrow(rc) ? row() : col());
  }

  void verifyclean() const { if (dirty_) 
      throw Failed("Object has been altered during operation!");
  }

  void setclean() { dirty_=false; }
  bool isdirty() const { return dirty_; }

  T& setdiagonal() { 
    rpoint1=cpoint1=true;
    stash2.clear(); //remove unused
    dirty_=true;
    return stash1;
  }

  T& setRC(char which) {
    const bool isr=isrow(which);
    if (rpoint1==cpoint1) {
      rpoint1=true; //set to off-diagonal
      cpoint1=false;
      if (isr) //clear stash that is *not* about to be addressed
	stash2.clear();
      else
	stash1.clear();
    }
    dirty_=true;
    return (isr ? row() : col());
  }
};

 template<class T> std::ostream& operator<< (std::ostream& ostr, const SwapStore<T>& a) {
   if (a.isdiagonal())
     return ostr << "Diagonal:\n" << a.row();
   else
     return ostr << "Rows:\n" << a.row() << "\nColumns:\n" << a.col(); 
 }

 //Implementation details
template<typename T> void SwapStore<T>::swap()
{
  verifyRC();
  rpoint1=!rpoint1;
  cpoint1=!cpoint1;
  dirty_=true;
  //  reset();
}

template<typename T> void SwapStore<T>::shuffle(char which)
{
  if (isrow(which)) {
    if (rpoint1==cpoint1)
      rpoint1=!rpoint1; // switch row away from current data
    else {
      rpoint1=!rpoint1;
      cpoint1=!cpoint1; // swap rows and columns
    }
    row().clear();
  }
  else {
    if (rpoint1==cpoint1)
      cpoint1=!cpoint1; // switch col away from current data
    else {
      rpoint1=!rpoint1;
      cpoint1=!cpoint1; // swap rows and columns
    }
    col().clear();
  }
  dirty_=true;
} 

template<class T> struct SignalGenerator_traits {};

 extern const char LCM_INCOMPAT[];

} //namespace libcmatrix
#endif
