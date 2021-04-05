#ifndef _timer_h_
#define _timer_h_

#include "basedefs.h"
#include <sys/types.h>
#include <ctime>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TIMES_H
#include <sys/time.h>
#endif

#ifdef CLOCK_PROCESS_CPUTIME_ID
#define LCM_GETTIME_CPUCLOCK CLOCK_PROCESS_CPUTIME_ID
#endif

#ifdef CLOCK_MONOTONIC
#define LCM_GETTIME_WALLCLOCK CLOCK_MONOTONIC
#else
#ifdef CLOCK_REALTIME
#define LCM_GETTIME_WALLCLOCK CLOCK_REALTIME
#endif
#endif

#ifndef LCM_DEFAULT_TIMER_TYPE
#define LCM_DEFAULT_TIMER_TYPE CPUTimer
#endif

namespace libcmatrix {

  //timers can only be defined in Unix-like environments with sys/time.h
  
  enum timer_type { CPUTimer, WallTimer };
  
  template<timer_type> struct timer_traits {};
  
  template<typename T> void timer_subtract(T& a, const T& b) { a-=b; }
  template<typename T> void timer_add(T& a, const T& b) { a+=b; }
  template<typename T> void timer_zero(T& a) { a=T(0); }
  template<typename T1, typename T2> inline void timer_subtract_(T1& asubsec, T2& asec, const T1 bsubsec, const T2 bsec, const int maxvalr) {
    const T1 maxval(maxvalr); 
    asubsec-=bsubsec;
    asec-=bsec;
    if (asubsec<0) {
      asubsec+=maxval;
      asec-=1;
    }
  }
  template<typename T1, typename T2> inline void timer_add_(T1& asubsec, T2& asec, const T1 bsubsec, const T2 bsec, const int maxvalr) {
    const T1 maxval(maxvalr); 
    asubsec+=bsubsec;
    asec+=bsec;
    if (asubsec>=maxval) {
      asubsec-=maxval;
      asec+=1;
    }
  }

#ifdef CLOCK_REALTIME
#define LCM_NANOSECONDS 1000000000
  inline void timer_subtract(timespec& a, const timespec& b) {
    timer_subtract_(a.tv_nsec,a.tv_sec,b.tv_nsec,b.tv_sec,LCM_NANOSECONDS);
  }
  inline void timer_add(timespec& a, const timespec& b) {
    timer_add_(a.tv_nsec,a.tv_sec,b.tv_nsec,b.tv_sec,LCM_NANOSECONDS);
  }
  inline void timer_zero(timespec& a) {
    a.tv_nsec=0;
    a.tv_sec=0;
  }
#endif
#ifdef HAVE_SYS_TIMES_H
#define LCM_MICROSECONDS 1000000
  inline void timer_subtract(timeval& a, const timeval& b) {
    timer_subtract_(a.tv_usec,a.tv_sec,b.tv_usec,b.tv_sec,LCM_MICROSECONDS);
  }
  inline void timer_add(timeval& a, const timeval& b) {
    timer_add_(a.tv_usec,a.tv_sec,b.tv_usec,b.tv_sec,LCM_MICROSECONDS);
  }
  inline void timer_zero(timeval& a) {
    a.tv_usec=0;
    a.tv_sec=0;
  }
#endif
  
  template<> struct timer_traits<CPUTimer> {
#ifdef LCM_GETTIME_CPUCLOCK
    typedef timespec time_type;
#else
    typedef clock_t time_type;
#endif
  };
  template<> struct timer_traits<WallTimer> {
#ifdef LCM_GETTIME_WALLCLOCK
    typedef timespec time_type;
#else
#ifdef HAVE_SYS_TIMES_H
    typedef timeval time_type;
#else
    typedef time_t time_type;
#endif
#endif
  };

  template<timer_type Type =LCM_DEFAULT_TIMER_TYPE> class timer {
  public:
    typedef typename timer_traits<Type>::time_type time_type;

    timer() { reset(); }
    void reset() { gettime(time_store); }
    static double resolution(); //!< allowed to return zero if resolution unknown
    static double toseconds(const time_type&);
    double operator()() const {
      gettime(time_tmp);
      timer_subtract(time_tmp,time_store);
      return toseconds(time_tmp); }
    static void gettime(time_type&);    
    
  private:
    time_type time_store;
    mutable time_type time_tmp;
  };
  
  template<timer_type Type =LCM_DEFAULT_TIMER_TYPE> class accumulating_timer {
  public:
    typedef typename timer_traits<Type>::time_type time_type;
    accumulating_timer() { reset(); }
    double operator()() const { return timer<Type>::toseconds(accumulated_); }
    size_t entries() const { return entries_; }
    double resolution() const { return timer<Type>::resolution(); }

    bool isactive() const { return active_; }
    void reset();
    void enter();
    void leave();

    class guard {
    public:
      guard(accumulating_timer<Type>&, bool allowrec =false);
      ~guard();
    private:
      accumulating_timer<Type>* stopwatchp_;
    };

  private:  
    time_type accumulated_,timepoint_;
    size_t entries_;
    bool active_;
  };  
  //#endif
  
  typedef timer<WallTimer> wtimer; //!< typedef for programs using wtimer
//   class wtimer {
//     time_t time_store;
//   public:
//     wtimer() { reset(); }
//     void reset() { time_store=time(0); }
//     double operator()() const { return (time(0)-time_store); }
//   };

  //! pretty print time \a t to stream \a ostr
  inline std::ostream& prettyprint_time(double t, std::ostream& ostr =std::cout, size_t switch_level =2) {
    cmatrix_ostream_controller(ostr).printtime(ostr,t,switch_level);
    return ostr;
  }

} //namespace libcmatrix
#endif

