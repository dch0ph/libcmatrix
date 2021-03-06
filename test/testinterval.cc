
#include "Warnings.h"
#include "cmatrix_utils.h"
#include "Histogram.h"
#include "MAS.h"
#include "ttyio.h"

using namespace libcmatrix;

//bool debug=false;

void test_interval(const char* name, IntervalSamplerBase& sampler, double period, double maxdt, size_t repeats)
{
  // if (debug)
  //   std::cout << name << ":\nRead out before set: " << sampler() << '\n';

  // sampler.set(period,maxdt,tol);

  // std::cout << "Intervals: ";

  //  IntervalSamplerBase::synchronise_t syncval=sampler.get_synchronisation();

  for (size_t i=0;i<repeats;i++) {
    // sampler.set_synchronisation(syncval);

    std::cout << "\nRepeat " << (i+1) << '\n';
    double t=0.0;
    while (t<period) {
      double dt=period-t;
      if (dt>maxdt)
	dt=maxdt;
      const double expmidt=t+0.5*dt;
      const double midt=sampler(t,t+dt);
      std::cout << (t*1e6) << "  " << (1e6*(t+dt)) << "  " << (midt*1e6) << "  " << ((midt-expmidt)*1e9) << " ns\n";
      t+=dt;
    }
  }

//   double tot=0.0;
//   for (;;) {
//     const double dt=sampler();
//     if (dt==0)
//       break;
//     tot+=dt;
//     std::cout << (dt*1e6) << " ";
//   }
//   std::cout << "\nTotal: " << (tot*1e6) << " us" << std::endl;

  // if (debug)
  //   std::cout << "Read out after exhausted: " << sampler() << '\n';
}

double getrotorperiod(IntervalSamplerBase& sampler, double period, double maxdt)
{
  const size_t nsteps=size_t(0.5+period/maxdt);
  const double dt=period/nsteps;
  double t=0.0;
  double totalvelocity=0.0;
  for (size_t i=0;i<nsteps;i++) {    
    const double expmidt=t+0.5*dt;
    const double midt=sampler(t,t+dt); //!< note that time out is a proxy for the rotor *phase* rather than time per se
    const double velocityscalefactor=1.0+(midt-expmidt)/expmidt;
    totalvelocity+=velocityscalefactor;
    t+=dt;
  }
  const double meanvelocityscale=totalvelocity/nsteps;
  // std::cout << meanvelocityscale << '\n';
  return period/meanvelocityscale;
}

int main(int argc,const char *argv[])
{
  int count=1;
  
  const double rotor_speed=50e3;
  const double rotor_period=1.0/rotor_speed;

  std::cout << "Rotor period: " << (1e6*rotor_period) << " us\n";

  const double maxdt=getfloat(argc,argv,count,"Maxdt (us)? ",1.0)*1e-6;
  // const double tol=1e-9;

  const double jitter=getfloat(argc,argv,count,"Jitter (ns)? ",0.0)*1e-9;
  double cortime=0.0;
  if (jitter)
    cortime=getfloat(argc,argv,count,"Correlation time (us) [0 if none]? ",0.0)*1e-6;
  bool randomise=false;

  if (jitter)
    randomise=getlogical(argc,argv,count,"Randomise? ");

  IntervalSampler regsampler(jitter,cortime,randomise);
  test_interval("Sampling",regsampler,rotor_period,maxdt,2);

  const int repeats=getint(argc,argv,count,"Repeats for rotor period distribution? ",10000);
  const size_t nbins=20;
  const double period_range=0.1*rotor_period;

  List<size_t> histstore(int(nbins),size_t(0));
  Histogram<size_t> periodhist(histstore,period_range,rotor_period-period_range/2.0);
  for (size_t i=0;i<repeats;i++) {
    const double period=getrotorperiod(regsampler,rotor_period,maxdt);
    periodhist.add(1,period);
  }
  for (size_t i=0;i<nbins;i++)
    std::cout << (1e6*periodhist.midpoint(i)) << " us: " << histstore(i) << '\n';

  std::cout << "\nSamples lost: " << periodhist.lostsum() << '\n';
  //  IntervalSamplerGaussian gsampler(jitter);
  //test_interval("Gaussian jitter sampling",gsampler,rotor_period,maxdt,tol);

  return 0;
}

  //  static Warning<> usedagain_warning;

  //protected:
  //double int_,maxdt_,tol_;
  //double curt_;
  //  bool donezero;

  //  void set_base(double,double, double);
  //  virtual void set_implementation() {}
  //  virtual double get_implementation() =0;

// void IntervalSamplerBase::set(double intv, double maxdtv, double tolv)
// {
//   int_=intv;
//   maxdt_=maxdtv;
//   tol_=tolv;
//   if ((int_ < 0.0) || (maxdt_<=tolv) || (tolv<=0.0))
//     throw InvalidParameter("IntervalSampler: invalid time interval / maximum sampling interval");
//   curt_=0.0;
//   donezero=false;
//   set_implementation();
// }

// Warning<> IntervalSamplerBase::usedagain_warning("IntervalSampler used incorrectly - called before set or re-used after exhausted",&lcm_base_warning);

// double IntervalSamplerBase::operator() () 
// {
//   if (donezero) {
//     usedagain_warning.raise();
//     return 0.0;
//   }
//   const double retval=get_implementation();
//   if (retval==0.0)
//     donezero=true;
//   return retval;   
// }

//   static const double jitterfac=1.5;

//   class IntervalSamplerGaussian : public IntervalSamplerBase {
//   public:
//     IntervalSamplerGaussian(double jitterv =0) { jitter(jitterv); }
//     IntervalSamplerBase* clone() { return new IntervalSamplerGaussian(*this); }
//     void set_implementation();
//     double get_implementation();
//     double jitter() const { return jitter_; }
//     void jitter(double);

//     static Warning<> excessivejitter_warning; //!< jitter is large compared to interval

//   private:
//     double jitter_;
//     bool usejitter_;
//     double usemaxdt_;
//   };

//   void IntervalSamplerGaussian::set_implementation()
//   {
//     usejitter_=(jitterfac*(jitter_+tol_)<int_);
//     if (usejitter_)
//       usemaxdt_=maxdt_-jitterfac*jitter_;
//     else {
//       if (jitter_)
// 	excessivejitter_warning.raise();
//       usemaxdt_=maxdt_;
//     }
//   }

//   double IntervalSamplerGaussian::get_implementation()
//   {
//     double dt=int_-curt_;
//     if (dt>maxdt_+tol_) {
//       dt=usemaxdt_;
//       if (usejitter_)
// 	dt+=gauss()*jitter_;
//       curt_+=dt;
//     }
//     else
//       curt_=int_;
//     return dt;
//   }

//     const double dt=int_-curt_;
//   if (dt>maxdt_+tol_) {
//     curt_+=maxdt_;
//     return maxdt_;
//   }
//   curt_=int_;
//   return dt;
