#ifndef _Euler_h_
#define _Euler_h_

#include <iostream>

namespace libcmatrix {

  struct Euler {
    double alpha;
    double beta;
    double gamma;
    Euler() : alpha(0.0), beta(0.0), gamma(0.0) {}
    Euler(double a,double b,double g) : alpha(a), beta(b), gamma(g) {}
    
    friend std::ostream &operator << (std::ostream &,const Euler &);
    
    Euler& operator*= (double x) { alpha*=x; beta*=x; gamma*=x; return *this; }
    Euler& operator/= (double x) { alpha/=x; beta/=x; gamma/=x; return *this; }

    friend bool operator== (const Euler& x, const Euler& y)
    { return (x.beta==y.beta) && (x.alpha==y.alpha) && (x.gamma==y.gamma); }

    friend bool operator!= (const Euler& x, const Euler& y)
    { return  (x.beta!=y.beta) || (x.alpha!=y.alpha) || (x.gamma!=y.gamma); }

    friend bool operator! (const Euler& x)
    { return (x.beta!=0.0) || (x.alpha!=0.0) || (x.gamma!=0.0); }

    Euler operator- () const;
  };


  class Euler_controller {
  public:
    Euler_controller();
    int verbose;
    double trig_tolerance; //!< tolerance for trig functions to be within +/- 1 
    double gimbal_tolerance; //!< tolerance for sin beta =0 (gimbal lock)
    double verify_tolerance; //!< tolerance for matrix norm when comparing recreated DCM (0 for no check)
    double asymmetry_tolerance; //!< force asymmetry to be zero if below this tolerance (not strictly Euler angles)
    size_t eigenvalue_flip_pattern; //!< eigenvector flip pattern (debug only)
    double filter_trigvalue(double) const; //!< constrain input to lie within -1 to 1
  };
  
  extern Euler_controller cmatrix_euler_controller; //!< default values of parameters
  
  template<typename T> class Matrix;
  Euler DCM_to_Euler(const Matrix<double>&, const Euler_controller& =cmatrix_euler_controller); //!< extract Euler angles from direction cosine matrix
  
}

#endif
