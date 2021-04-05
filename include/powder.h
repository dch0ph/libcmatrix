#ifndef _Powder_h_
#define _Powder_h_

/* Powder averaging methods */

#include "Matrix.h"
#include "Euler.h"
#include "smartptr.h"
#include "Warnings.h"

namespace libcmatrix {

  enum range_t { sphere, hemisphere, octant };

class PowderMethod {
public:
  PowderMethod() : angles_(0) { max=current=0; }
  PowderMethod(size_t _max, size_t anglesv)
    : angles_(anglesv)
  { create(_max); }

  enum include_t { middle, start, end, both };

  virtual ~PowderMethod() {}
  virtual PowderMethod* clone() const =0;

  static double alphamax(range_t);
  static double betamax(range_t);
  
  //  static range_t range_type(const char*);
  size_t angles() const { return angles_; }

  static void include_parameters(double&, double&, double, size_t, include_t);

  size_t orientations() const { return max; }
  void orientation(Euler& dest, double&, size_t) const;

  virtual void index_to_orientation(Euler &,double &,size_t) const =0;
  
  bool next(Euler &,double &);
  void reset() { current=0; }

private:
  size_t max,current;

protected:
  size_t angles_;
  void create(size_t);
  double normal;

  static void alpha_parameters(double&, double&, size_t, range_t, include_t);
};

// Dummy method with only one orientation!
class PowderSingle : public PowderMethod {
public:
  PowderSingle(const Euler& anglev =Euler(0,0,0))
    : PowderMethod(1,0), angle(anglev) {}

  void index_to_orientation(Euler& dest,double& scale, size_t) const
    { scale=1.0; dest=angle; }

  PowderMethod* clone() const { return new PowderSingle(*this); }

private:
  Euler angle;
};

// Equal steps over alpha and beta
class PlanarGrid : public PowderMethod {
public:
  PlanarGrid(size_t asteps, size_t bsteps, range_t =sphere, include_t =middle);
  void index_to_orientation(Euler&, double&, size_t) const;
  PowderMethod* clone() const { return new PlanarGrid(*this); }

private:
  size_t asteps;
  size_t bsteps;
  double astep,bstep;
  bool issingle;
  double aoffset,boffset;

  double getbeta(size_t whichb) const;
};

int fibonacci(int);

// ZCW sampling on alpha, beta grid
class PlanarZCW : public PowderMethod {
public:
  PlanarZCW(size_t, range_t = sphere, include_t =middle);
  void index_to_orientation(Euler &,double &,size_t) const;
  PowderMethod* clone() const { return new PlanarZCW(*this); }

  static size_t N_to_orientations(size_t);
  static int orientations_to_N(size_t);

private:
  int g1,g2;
  double ascale,bscale,aoffset,boffset;

  double getbeta(size_t) const;
};

//! ZCW sampling on alpha, beta, gamma grid
class ZCW3 : public PowderMethod {
public:
  ZCW3(size_t, range_t = sphere, include_t =middle, double =2*M_PI);
  void index_to_orientation(Euler &,double &,size_t) const;
  PowderMethod* clone() const { return new ZCW3(*this); }

  static size_t N_to_orientations(size_t);
  static int orientations_to_N(size_t);

private:
  int mod1,mod2,mod3;
  double ascale,bscale,gscale,aoffset,boffset;

  double getbeta(size_t) const;
};

//! ZCW sampling with alpha, beta taken over sphere
class SphericalZCW : public PowderMethod {
public:
  SphericalZCW(size_t, range_t =sphere, include_t =middle);
  void index_to_orientation(Euler &,double &,size_t) const;
  PowderMethod* clone() const { return new SphericalZCW(*this); }

  static size_t N_to_orientations(size_t n) { return PlanarZCW::N_to_orientations(n); }
  static int orientations_to_N(size_t n) { return PlanarZCW::orientations_to_N(n); }

private:
  int g1,g2;
  double ascale,aoffset,bscale,boffset;
};

class ExplicitSampling : public PowderMethod {
public:
  ExplicitSampling(const Matrix<double>&);
  void index_to_orientation(Euler&, double&, size_t) const;
  PowderMethod* clone() const { return new ExplicitSampling(*this); }

private:
  const Matrix<double>& angweights;
};

// Regular sample of alpha, beta treated as solid angles  
class SphericalGrid : public PowderMethod {
public:
  SphericalGrid(size_t asteps,size_t bsteps, range_t =sphere, include_t =middle);
  void index_to_orientation(Euler&, double&, size_t) const;
  PowderMethod* clone() const { return new SphericalGrid(*this); }

private:
  size_t asteps;
  size_t bsteps;
  double astep,bstep,aoffset,boffset;
};

class WithGamma : public PowderMethod {
public:
  WithGamma(const PowderMethod&, size_t);
  void index_to_orientation(Euler&, double&, size_t) const;
  PowderMethod* clone() const { return new WithGamma(*this); }
  static Warning<> singlecombine_warning;
private:
  const smartptr<PowderMethod> methp;
  size_t gsteps;
  double scale;
};

} //namespace libcmatrix

#endif 
