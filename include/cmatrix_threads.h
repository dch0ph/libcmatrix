#ifndef cmatrix_threads_h_
#define cmatrix_threads_h_

#include "config.h"

#if !defined(HAVE_LIBPTHREAD)
#error "This won't work! You don't have a functional pthread or didn't enable threads!"
#endif

#include <cstring>
#include "lcm_basethreads.h"
#include <pthread.h>
#include "List.h"

namespace libcmatrix {

 class thread_controller;

struct thread_block {
  thread_controller* tconp;
  size_t nthread;
  pthread_t tid;
};

class thread_controller {
  pthread_mutex_t lock;
  pthread_cond_t start_cond,done_cond;
  size_t chunk;
  size_t current;
  size_t todo;
  size_t notdone;
  size_t workers;
  const BaseThreadFunction* fobj;
  List<thread_block> tblocks;

  bool in_parallel_;

  thread_controller(const thread_controller&);
  void operator= (const thread_controller&); // can't copy objects
	      
 public:
  explicit thread_controller(int n)
    : workers(0), fobj(NULL)
    { create(n); }

  thread_controller() 
    : workers(0), in_parallel_(false) {}

  ~thread_controller();

  void create(int n);

  void run(const BaseThreadFunction&,size_t,size_t);
  void start(const BaseThreadFunction&,size_t,size_t);
  void wait();

  bool in_parallel() const
  { return in_parallel_; }

  size_t get_num_threads() const
  { return (in_parallel() ? workers : 0); }

  size_t get_max_threads() const
  { return workers; }

  int get_thread_num() const;

  friend void *worker(void *);
};

} //namespace libcmatrix
      
#endif 
