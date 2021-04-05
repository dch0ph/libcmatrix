/* Statistics stuff */

#include <cstdlib>
#include <ctime>
#include "basedefs.h"
#include "BaseList.h"
#include "cmatrix_utils.h"
#include <iostream>
#include "lcm_basethreads.h"

namespace libcmatrix {

  static const int lcm_default_seed=1;

// static int seedval=-1;
  
// int get_seed_default()
// {
//   return -1;
// }

  int RandomGenerator::get_seed_random()
  {
    static int calls=0;
    calls++; //!< fix problem if get_seed_random called too frequently! 24/1/14
    return int((calls+time(NULL)) &0x7fff);
  }

  int RandomGenerator::get_seed_default()
  {
    return lcm_default_seed;
  }

/* Routines taken from Numerical Recipes - new interface to gauss() */

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

  template<class Gen> void GaussianGenerator<Gen>::set_correlationtime(double cortimv, double chartimv)
  {
    if (cortimv<0.0)
      throw InvalidParameter("GaussianGenerator: correlation time cannot be <0");
    cortim_=cortimv;
    if (cortim_==0.0) {
      lme2_=-2.0;
      weightfac_=0.0;
    }
    else
      set_chartime(chartimv);
    reset();
  }

  /* Original source of this code (CGAUSS) to generate coloured Gaussian noise is unclear 
     Comments refer to article by Fox et al., Physical Review A vol.38(1988) page 5938. */

  template<class Gen> void GaussianGenerator<Gen>::set_chartime(double chartimv)
  {
    if (chartimv<0.0)
      throw InvalidParameter("GaussianGenerator: characteristic time cannot be <0");
    chartim_=chartimv;
    if (chartimv!=0.0) {
      if (cortim_==0.0) 
	throw InvalidParameter("GaussianGenerator: cannot set characteristic time when generating white noise");
      weightfac_=std::exp(-chartim_/cortim_);
      //      lme2_=(1.0-weightfac_*weightfac_)*(-2.0/cortim_);
      lme2_=(1.0-weightfac_*weightfac_)*-4.0; //!< adjusted expression keeps amplitude fixed
#ifndef NDEBUG
      std::cout << "GaussianGenerator: setting prev-timepoint weight factor to " << weightfac_ << "  amp scaling factor to " << std::sqrt(-lme2_) << '\n';
#endif
    }
    iset=false; //!< force re-start of Guassian generator
  }

  template<class Gen> void GaussianGenerator<Gen>::reset()
  {
    iset=false;
    prev_=0.0;
  }

  template<class Gen> double GaussianGenerator<Gen>::operator() ()
  {
    double retval=gset;
    double v1=0.0;
    if  (iset)
      iset=false;
    else {
      double r,v2;
      do {
	v1=2.0*rgen_()-1.0;
	v2=2.0*rgen_()-1.0;
	r=v1*v1+v2*v2;
      } while (r >= 1.0);
      const double fac=std::sqrt(lme2_*std::log(r)/r);
      gset=v1*fac;
      iset=true;
      retval=v2*fac;
    }
    if (cortim_!=0.0) {
      if (chartim_==0.0)
	throw Failed("GaussianGenerator: characteristic time unset");
      if (iset && (prev_==0.0)) { //!< if no previous time-point, get one from cached value
	prev_=v1; //!< unscaled value better reflects amplitude
	iset=false;
      }
      retval+=prev_*weightfac_;
      prev_=retval;
    }
    return retval;
  }

  template<class Gen> double GaussianGenerator<Gen>::operator() (double chartimv)
  {
    if ((cortim_!=0.0) && (std::fabs(chartimv-chartim_)>1e-12))
      set_chartime(chartimv);
    return (*this)();
  }

  template class GaussianGenerator<RandomGenerator>;

  static RandomGenerator rgen;
  static GaussianGenerator<> ggen; //!< contains its own RandomGenerator but only initialised as required
  
  // // Assume that seedval can be set atomically
  void set_seed(int val)
  {
    rgen.set_seed(-std::abs(val));
  }
  
//   int get_seed_random()
//   {
//     return int(time(NULL) &0x7fff);
//   }

  void set_seed()
  {
    rgen.set_seed(RandomGenerator::get_seed_random());
  }

  RandomGenerator::RandomGenerator() : seedval(-lcm_default_seed) {}
  RandomGenerator::RandomGenerator(randomseed_t seedval) { set_seed(seedval); }
  
  void RandomGenerator::set_seed(randomseed_t val)
  {
    if (val<=0)
      throw InvalidParameter("RandomGenerator: seed value cannot be <=0");
    seedval=-val;
  }
  
  double RandomGenerator::operator() ()
  {
    int j;
  
    if (seedval<0) {
      ix1=(IC1-seedval) % M1;
      ix1=(IA1*ix1+IC1) % M1;
      ix2=ix1 % M2;
      ix1=(IA1*ix1+IC1) % M1;
      ix3=ix1 % M3;
      for (j=1;j<=97;j++) {
	ix1=(IA1*ix1+IC1) % M1;
	ix2=(IA2*ix2+IC2) % M2;
	r[j]=(ix1+ix2*RM2)*RM1;
      }
      seedval=1; //!< flag initialised by setting seed to +ve value
    }
    ix1=(IA1*ix1+IC1) % M1;
    ix2=(IA2*ix2+IC2) % M2;
    ix3=(IA3*ix3+IC3) % M3;
    j=1 + ((97*ix3)/M3);
    if (j > 97 || j < 1)
      throw InternalError("nr_random");
    const double temp=r[j];
    r[j]=(ix1+ix2*RM2)*RM1;
    return temp;
  }

double random(double sigma, bool unsafe)
{
  static Mutex<ThreadingActive> lock;
  if (!unsafe)
    lock.acquire();
  const double res=sigma*rgen();
  if (!unsafe)
    lock.release();
  return res;
}

double gauss(double sigma, bool unsafe)
{
  static Mutex<ThreadingActive> lock;
  if (!unsafe)
    lock.acquire();

  const double res=ggen();

  if (!unsafe)
    lock.release();

  return res;
}

// double nr_random(int *idum, bool unsafe)
// {
//   static Mutex<ThreadingActive> lock;

//   if (!unsafe)
//     lock.acquire();
//   const double res=nr_random_unsafe(idum);
//   if (!unsafe)
//     lock.release();
//   return res;
// }

}//namespace libcmatrix
