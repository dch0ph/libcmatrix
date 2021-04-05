#ifndef lcm_basethreads_h_
#define lcm_basethreads_h_

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif

//Required for Failed
#include "basedefs.h"
#include "Warnings.h"
#include "timer.h"

namespace libcmatrix {

#ifdef LCM_REENTRANT
  const bool ThreadingActive=true;
#else
  const bool ThreadingActive=false;
#endif

  template<bool =true> class Mutex {
  public:
    void acquire() {}
    void release() {}
  };

  template<> class Mutex<true> {
  public:
    static Warning<> mutexnotfunctional_warning;

#ifdef HAVE_PTHREAD_H
    Mutex() {
      pthread_mutex_init(&lock_,NULL);
    }
    ~Mutex() {
      pthread_mutex_destroy(&lock_);
    }
    void acquire() {
      pthread_mutex_lock(&lock_);
    }
    void release() {
      pthread_mutex_lock(&lock_);
    }
    pthread_mutex_t& operator()() { return lock_; }
    const pthread_mutex_t& operator()() const { return lock_; }

  private:
    pthread_mutex_t lock_;
#else
    Mutex() { mutexnotfunctional_warning.raise(); }
    void acquire() {}
    void release() {}
#endif
  };

  //automatically acquires lock and releases when destroyed (exception safe)
  template<bool Active =true> class MutexLock {
  public:
    MutexLock(Mutex<Active>& lockv)
      : lock_(lockv) { lock_.acquire(); }
      ~MutexLock() { lock_.release(); }
  private:
      Mutex<Active>& lock_;
  };

  //empty class when not locking
  template<> class MutexLock<false> {
  public:
    MutexLock(Mutex<false>&) {}
  };

extern "C" {
  typedef void (*P_DESTFUNC)(void*);
  void* worker(void *);
}

#ifdef LCM_REENTRANT
template<class T> class key_holder {
  pthread_key_t key;
public:
  key_holder(P_DESTFUNC kill_func =0) {
    if (int errcode=pthread_key_create(&key,kill_func))
      throw Failed(strerror(errcode));
  }
  key_holder(P_DESTFUNC kill_func,T* ivalue) {
    if (int errcode=pthread_key_create(&key,kill_func))
      throw Failed(strerror(errcode));
    set(ivalue);
  }
  void set(T* value) { 
    if (int errcode=pthread_setspecific(key,value))
      throw Failed(strerror(errcode));
  }
  T* get() const
    { return (T*)pthread_getspecific(key); }

  ~key_holder() { pthread_key_delete(key); }
  // NB cannot call destructors!  Should only be called once 
};
#else
template<class T> class key_holder {
  T* value;
  P_DESTFUNC kill_funcp;
public:
  explicit key_holder(P_DESTFUNC kill_func =0,T* ivalue =NULL)
    : value(ivalue),kill_funcp(kill_func) {}
  void set(T* ivalue) { value=ivalue; }
  T* get() const { return value; }
  ~key_holder() { 
    if (kill_funcp && value)
      (*kill_funcp)((void*)value);
  }
  // NB cannot call destructors!  Should only be called once 
};
#endif

typedef void (*P_THREADFUNC)(size_t,size_t,size_t);

 struct BaseThreadFunction {
   virtual ~BaseThreadFunction() {}
   virtual void operator()(size_t,size_t,size_t) const =0;
 };


class ThreadFunction : public BaseThreadFunction {
 public:
  explicit ThreadFunction(P_THREADFUNC func_) : func(func_) {
    if (!func)
      throw InvalidParameter("ThreadFunction(): NULL function");
  }

  void operator()(size_t st, size_t en, size_t nthr) const
  { (*func)(st,en,nthr); }

 private:
   P_THREADFUNC func;
 };

#ifndef LCM_PARALLEL_DEFAULT_VERBOSE
#ifdef NDEBUG
#define LCM_PARALLEL_DEFAULT_VERBOSE 0
#else
#define LCM_PARALLEL_DEFAULT_VERBOSE 1
#endif
#endif
  
  class base_parallel_controller {
  public:
    base_parallel_controller(int verbosev =LCM_PARALLEL_DEFAULT_VERBOSE)
      : id_(0), verbose_(verbosev), log(NULL) {}

    ~base_parallel_controller();

    int get_thread_num() const { return id_; }
    bool ammaster() const { return (id_==0); }
    int get_processors() const { return nprocs_; } //!< get number of processors
    int get_workers() const { return nworkers_; } //!< get number of actual workers
    void verbose(int verbosev) { verbose_=verbosev; }
    static Warning<> logopenfailed_warning; //!< failed to open log file
    
  protected:
    size_t nprocs_;
    size_t nworkers_;
    int id_;
    int verbose_;

    void error(const char* mess) const {
      writelog(mess);
      throw InternalError(mess);
    }
    void openlog(const char*);
    void writelog(const char*) const;
    void closelog();
    double time() const { return timer_(); }

  private:
    FILE* log;    
    timer<WallTimer> timer_;
  };


} //namespace libcmatrix
    
#endif
