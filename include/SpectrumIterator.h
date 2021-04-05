#ifndef _SpectrumIterator_h_
#define _SpectrumIterator_h_

/* Base class for iterating over spectrum */

namespace libcmatrix {
  template<class T> class SpectrumIterator {
  public:
    SpectrumIterator() {}
    virtual ~SpectrumIterator() {};
    
    typedef T amplitude_type;
    virtual bool operator () (T& amp,double& freq) =0;
  };
}

#endif
