#ifndef lcm_fork_controller_h_
#define lcm_fork_controller_h_

/*! \file
 \brief  Parallisation using fork and shared memory
*/

//Simplifies matters!
#ifdef LCM_REENTRANT
#error "Fork_controller can't be combined with multi-threading."
#endif

#include <cassert>
#include <sys/types.h>
#include "lcm_basethreads.h"
#include "BaseList.h"

#ifndef HAVE_PTHREAD_H
#error Fork_controller needs functional pthread
#endif

namespace libcmatrix {

  class shared_memory_obj {
  public:
    shared_memory_obj() : rawptr_(NULL), alignedptr_(NULL), alloc_(0), shmid_(-1) {}
    ~shared_memory_obj() { zap(); }
    void* attach(int,size_t,size_t);
    void* operator()() { return alignedptr_; }
    void* ensure(size_t, size_t);
    int shmid() const { return shmid_; }
    void kill() { zap(); }
  private:
    void* rawptr_;
    void* alignedptr_;
    size_t alloc_;
    int shmid_;
    bool amowner_;
    void zap();
    void* rawattach(int,bool,size_t,size_t);

    static Warning<> shmat_warning; //!< failed to attach shared memory
    static Warning<> shmget_warning; //!< failed to claim shared memory
  };
    
  class Fork_controller : public base_parallel_controller {
  public:
    Fork_controller(size_t, int =LCM_PARALLEL_DEFAULT_VERBOSE);
    ~Fork_controller();
    
    int get_thread_num() const { return id_; }
    bool ammaster() const { return (get_thread_num()==0); }
    int get_processors() const { return nprocs_; }

    void start();
    bool join();
    void join_kill(); //!< join outstanding threads and cleanup
    bool isactive() const;
    void run(const BaseThreadFunction&, size_t, size_t);
    template<typename T> bool verify(const BaseList<T>&);

    template<typename T> void sum_join(BaseList<T> dest, const BaseList<T>& source);
    template<typename T> void sum_join(T& dest, const T& source) {
      duplicate_structure(dest,source);
      sum_join(dest.row(),source.row());
    }

    static Warning<> destroywhileactive_warning; //!< Fork_controller destroyed while threads active

    template<typename T> void broadcast(BaseList<T>) { throw InternalError("broadcast not implemented for Fork_controller"); }

    static const long pollns=500; //!< nanoseconds to wait during synchronisation polling 
    size_t checkjobstates(size_t) const;

private:
    shared_memory_obj shdctrl;
    shared_memory_obj shddata;
    BaseList<pid_t> pids;
    size_t jobstate;
    BaseList<size_t> jobstates;
    size_t databytes_; //!< record of data claim
    size_t alignment_; 

    void* ensuredata(size_t, size_t);
    static void poll();
};

  // implementation details below here

  template<typename T> bool Fork_controller::verify(const BaseList<T>& source)
    {
      const size_t n(source.size());
      return (ensuredata(n*(nprocs_-1)*sizeof(T),memory_traits<T>::alignment)!=NULL);
    }

  template<typename T> void Fork_controller::sum_join(BaseList<T> dest, const BaseList<T>& source) {
    const size_t n(source.size());
    if (n!=source.size())
      throw Mismatch("Fork_controller::sum",n,source.size());
    T* const datap=reinterpret_cast<T*>(ensuredata(n*(nprocs_-1)*sizeof(T),memory_traits<T>::alignment));
    if (datap==NULL)
      throw Failed("sum_join: failed to allocate shared memory buffer");
    jobstate++;
    if (!ammaster()) {
      BaseList<T> tmpd(n,datap+(id_-1)*n);
      jobstates(id_-1)=jobstate;
      tmpd=source;
    }
    if (!join())
      throw Failed("sum_join: at least one process terminated abnormally - calculation invalid");
    if (checkjobstates(jobstate))
      throw Failed("sum_join: not all process provided data");
    dest=source; 
    for (size_t i=nprocs_-1;i--;)
      dest+=BaseList<T>(n,datap+i*n); 
  }
  
} //namespace libcmatrix

#endif

