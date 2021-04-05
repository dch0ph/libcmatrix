#include "cmatrix.h"
#include "rmatrix.h"
#include <cstdio>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <cstdlib>
#include <iostream>
#include <cstring>

#ifdef LCM_REENTRANT
#define ISREENTRANT
#endif

#include "cmatrix_threads.h"

namespace libcmatrix {

#ifdef ISREENTRANT

extern "C" void *worker(void *tblockp)
{
  thread_block &tblock=*((thread_block *)tblockp);
  thread_controller& tcon=*(tblock.tconp);
  const size_t nthread=tblock.nthread;  

  size_t start_pow,end_pow,this_chunk;

  if (!tcon.fobj) throw Failed("thread_controller: function not set");

// #ifdef MAKE_CSTACK
//   cmatrix_new_thread();
// #endif
// #ifdef MAKE_RSTACK
//   rmatrix_new_thread();
// #endif

  for (;;) {
    pthread_mutex_lock(&tcon.lock);
    
    while (tcon.current==tcon.todo)
      pthread_cond_wait(&tcon.start_cond,&tcon.lock);
    
    start_pow=tcon.current;
    this_chunk = (start_pow+tcon.chunk>tcon.todo) ? tcon.todo-start_pow : tcon.chunk;
    
    end_pow=start_pow+this_chunk;
  
    tcon.current=end_pow;
    
    pthread_mutex_unlock(&tcon.lock);
    
    (*(tcon.fobj))(start_pow,end_pow,nthread);

    pthread_mutex_lock(&tcon.lock);
    tcon.notdone-=end_pow-start_pow;
    if (tcon.notdone==0) pthread_cond_signal(&tcon.done_cond);
    pthread_mutex_unlock(&tcon.lock);
  }
  return NULL; //doesn't actually return, but keeps compiler happy
}

  //a bit stupid
  int thread_controller::get_thread_num() const 
  {
    const pthread_t self=pthread_self();
    for (size_t i=workers;i--;) {
      if (pthread_equal(self,tblocks(i).tid)) return i;
    }
    return -1;
  }

  void thread_controller::create(int trywork)
{
  if (workers) throw Failed("thread_controller::create: threads already created");
  if (trywork<1) throw InvalidParameter("thread_controller");

  in_parallel_=false;
  
  pthread_cond_init(&start_cond,NULL);
  pthread_cond_init(&done_cond,NULL);
  pthread_mutex_init(&lock,NULL);

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);

  size_t ccon=pthread_getconcurrency();

  tblocks.create(trywork);
  for (size_t i=0;i<trywork;i++) {
    thread_block& tblock=tblocks(i);
    tblock.tconp=this;
    tblock.nthread=i;
    if (int errcode=pthread_create(&(tblock.tid),NULL,worker,&tblock)) 
      throw Failed(strerror(errcode));
    workers++;
  }

  pthread_setconcurrency(workers+ccon);
}

void thread_controller::run(const BaseThreadFunction& fobjv,size_t orients,size_t rchunk)
{
  start(fobjv,orients,rchunk);
  wait();
}

void thread_controller::start(const BaseThreadFunction& fobjv,size_t orients,size_t rchunk)
{
  if (in_parallel_)
    throw Failed("start: already active");

  fobj=&fobjv;

  pthread_mutex_lock(&lock);
  
  in_parallel_=true;
  current=0;
  todo=notdone=orients;
  chunk=rchunk;
  
  pthread_cond_broadcast(&start_cond);
}

void thread_controller::wait()
{
  if (!in_parallel_) return;

  while (notdone) pthread_cond_wait(&done_cond,&lock);
  in_parallel_=false;
  pthread_mutex_unlock(&lock);
}

thread_controller::~thread_controller()
{
  if (!workers) return;
  for (size_t i=0;i<workers;i++) 
    pthread_cancel(tblocks(i).tid);
  pthread_cond_destroy(&start_cond);
  pthread_cond_destroy(&done_cond);
  pthread_mutex_destroy(&lock);
}

#else

thread_controller::~thread_controller() {}
void thread_controller::wait() {}
void thread_controller::run(const BaseThreadFunction& fobj_,size_t orients,size_t) { fobj_(0,orients-1,0); }

void thread_controller::create(int)
{
  in_parallel_=false;
  workers=1;
}
#endif

}
// int get_workers()
// {
//   const char *envp=getenv("PARALLEL");
//   if (envp) {
//     int threads;
//     if (sscanf(envp,"%i",&threads)!=1)
//       throw Failed("PARALLEL not a valid number");
//     return threads;
//   }
// #ifdef _SC_NPROCESSORS_ONLN
//   return sysconf(_SC_NPROCESSORS_ONLN);
// #else
//   return 0;
// #endif
// }
