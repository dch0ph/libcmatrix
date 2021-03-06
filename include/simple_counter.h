#ifndef LCM_simple_counter_h_
#define LCM_simple_counter_h_

/*! \file
 \brief  Object keeping track of bytes used
*/
#include "basedefs.h"
#include "Warnings.h"

namespace libcmatrix {

  class simple_counter {
  public:
    simple_counter(long limitv, long usedv =0);

    simple_counter& operator+=(long);
    simple_counter& operator-=(long valv) { return operator+=(-valv); }
    
    size_t slots(size_t n) const { 
      return (used_>=limit_) ? 0 : (limit_-used_)/n;
    }
    
    void limit(long);
    long limit() { return limit_; }
    long operator()() const { return used_; }
    bool full() const { return (used_>=limit_); }

    template<class T> static long needed(const T& a) {
      return a.size()*sizeof(LCM_VAL(T));
    }
    static Warning<InternalError> negative_warning;

  private:
    long limit_;
    long used_;
  };

  //! overrideable template
  template<class T> size_t simple_counter_needed(const T& a) { return simple_counter::needed(a); }
}

#endif  
  
