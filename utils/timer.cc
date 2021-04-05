
#include "basedefs.h"
#include "timer.h"
//#ifdef HAVE_SYS_TIMES_H
//#include <sys/times.h>
//#endif
#include <cmath>

namespace libcmatrix {
  
#ifndef LCM_GETTIME_CPUCLOCK
#ifndef LCM_GETTIME_WALLCLOCK
#define LCM_NO_GETTIME
#endif
#endif

#ifndef LCM_NO_GETTIME
  namespace {
    double getresolution(clockid_t clk_id) {
      timespec spec;
      if (clock_getres(clk_id,&spec))
	throw Failed("Failed to get timer resolution");
      return spec.tv_sec+spec.tv_nsec*1e-9;
    }
  }
#endif

#ifdef LCM_GETTIME_CPUCLOCK
  static const double cpu_timer_resolution_(getresolution(LCM_GETTIME_CPUCLOCK));
#else
  //#ifdef HAVE_SYSCONF
//static const size_t clocks_per_sec(sysconf(_SC_CLK_TCK));
//#else
#ifdef CLOCKS_PER_SEC
  static const size_t clocks_per_sec(CLOCKS_PER_SEC);
#else
#ifdef CLK_TCK
  static const size_t clocks_per_sec(CLK_TCK);
#else
#error "Failed to define system clock resolution"
#endif
#endif
  static const double cpu_timer_resolution_(1.0/double(clocks_per_sec));
#endif 

  template<> double timer<CPUTimer>::resolution() { return libcmatrix::cpu_timer_resolution_; }

#ifdef LCM_GETTIME_WALLCLOCK
  static const double wall_timer_resolution_(getresolution(LCM_GETTIME_WALLCLOCK));
#else
  static const double wall_timer_resolution_=0.0;
#endif 

  template<> double timer<WallTimer>::resolution() { return wall_timer_resolution_; }
  
  template<> double timer<CPUTimer>::toseconds(const time_type& t)
  {
#ifdef LCM_GETTIME_CPUCLOCK
    return t.tv_sec+1e-9*t.tv_nsec;
#else
    return t*cpu_timer_resolution_;
#endif
  }

  template<> void timer<CPUTimer>::gettime(time_type& dest) {
    //struct tms buffer;      
    //times(&buffer);
    //return buffer.tms_utime;
#ifdef LCM_GETTIME_CPUCLOCK
    clock_gettime(LCM_GETTIME_CPUCLOCK,&dest);
#else
    dest=std::clock();
#endif
  }

  template<> void timer<WallTimer>::gettime(time_type& dest) {
#ifdef LCM_GETTIME_WALLCLOCK
    clock_gettime(LCM_GETTIME_WALLCLOCK,&dest);
#else
#ifdef HAVE_SYS_TIMES_H
    gettimeofday(&dest,NULL);
#else
    dest=time(0);
#endif
#endif
  }

  template<> double timer<WallTimer>::toseconds(const time_type& tp)
  {
#ifdef LCM_GETTIME_WALLCLOCK
    return double(tp.tv_sec)+(tp.tv_nsec/1e9);
#else
#ifdef HAVE_SYS_TIMES_H
    return double(tp.tv_sec)+(tp.tv_usec/1e6);
#else
    return double(tp);
#endif
#endif
  }

  template<timer_type Type> void accumulating_timer<Type>::reset() { 
    timer_zero(accumulated_);
    active_=false;
    entries_=0;
  }

  template<timer_type Type> void accumulating_timer<Type>::enter()
  {
    if (active_)
      throw Failed("accumulating_timer::enter: timer re-entered");
    timer<Type>::gettime(timepoint_);
    active_=true;
    entries_++;
  }

  template<timer_type Type> accumulating_timer<Type>::guard::guard(accumulating_timer<Type>& stopwatch, bool allowrec) 
    : stopwatchp_(NULL) {
    if (!allowrec || !stopwatch.isactive()) {
      stopwatch.enter();
      stopwatchp_=&stopwatch;
    }
  }

  template<timer_type Type> accumulating_timer<Type>::guard::~guard()
  {
    if (stopwatchp_)
      stopwatchp_->leave();
  }
  
  template<timer_type Type> void accumulating_timer<Type>::leave()
  { 
    if (!active_)
      throw Failed("accumulating_timer::leave: timer not active");
    time_type tmp;
    timer<Type>::gettime(tmp);
    timer_subtract(tmp,timepoint_);
    timer_add(accumulated_,tmp);
    active_=false;
  }

  template class timer<CPUTimer>;
  template class timer<WallTimer>;
  template class accumulating_timer<CPUTimer>;
  template class accumulating_timer<WallTimer>;

  std::ostream& ostream_controller::printtime(std::ostream& ostr, double eta, size_t switch_level) const
  {
    static const size_t levels=6;
    static const char* units[levels]={"h","m","s","ms","us","ns"};
    static const double multiplier[levels]={3600.0,60.0,1.0,1e-3,1e-6,1e-9};
    static const int sub[levels]={2,2,2,3,3,3};
    static const size_t default_level=2; //!< default if t is too small be represented 
      
    if (switch_level>=levels)
      throw InvalidParameter("printtime (switch_level)");

    if (eta<0) { //!< negative time?
      ostr << '-';
      return printtime(ostr,-eta,switch_level);
    }
    
    size_t level=default_level; //!< default if t=0
    if (eta) {      
      int digits= (timeprecision<0) ? ostr.precision() : timeprecision;

      bool donext=false;      
      for (level=0;level<levels;level++) {
	const double curmult=multiplier[level];
	if (donext || (eta>=curmult)) {
	  if (donext)
	    ostr << ' '; //!< space between time elements
	  if (digits<=sub[level])
	    return ostr << floor(0.5+eta/curmult) << ' ' << units[level];
	  if (!donext && (level>=switch_level))
	    break;
	  ostr << floor(eta/curmult) << ' ' << units[level];	
	  digits-=sub[level];
	  eta=fmodf(eta,curmult);
	  donext=true;
	}
      }
      //! fall through
      if (level==levels)
	level--;
//       else { //!< if using single unit, use default level (s?) i.e. avoid 1.23 h
// 	if (level<default_level)
// 	  level=default_level;
//       }
    }
    return ostr << (eta/multiplier[level]) << ' ' << units[level];    
  }

}
//#endif

