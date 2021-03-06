#include "Fork_controller.h"
#include <sys/shm.h>
#include <sys/wait.h>
#include <time.h>
#include <cstdio>

namespace libcmatrix {

  Fork_controller::Fork_controller(size_t np, int verbosev)
    : base_parallel_controller(verbosev),
    jobstate(0),
      databytes_(0)
  {
    if (np<2)
      throw InvalidParameter("Fork_controller: must be initialised with >1 processor!");
    nprocs_=np;
    nworkers_=np;
    const size_t slaves=np-1;
    const size_t ctrlbytes=sizeof(int)+slaves*(sizeof(pid_t)+sizeof(size_t));
    void* ctrlp=shdctrl.ensure(ctrlbytes,sizeof(size_t));
    int* shdkey=reinterpret_cast<int*>(ctrlp);
    *shdkey=-1;
    //create vectors in shared control segment
    char* bytep=reinterpret_cast<char*>(ctrlp);
    pids.create(slaves,reinterpret_cast<pid_t*>(sizeof(int)+bytep));
    jobstates.create(slaves,reinterpret_cast<size_t*>(sizeof(int)+bytep+slaves*sizeof(pid_t)));
    if (((char*)(jobstates.vector())+slaves)>((char*)ctrlp+ctrlbytes))
      throw InternalError("Fork_controller: memory usage overrun");
    pids=pid_t(0); //flag not set
  }

  bool Fork_controller::isactive() const { return pids( id_ ? (id_-1) : size_t(0)); }

  size_t Fork_controller::checkjobstates(size_t check) const
  {
    size_t failed=0;
    for (size_t i=nprocs_-1;i--;) {
      if (jobstates(i)!=check)
	failed++;
    }
    return failed;
  }

  void Fork_controller::start()
  {
    if (isactive())
      throw Failed("Fork_controller::start: fork already done");
    //    jobstates=0;
    const size_t slaves=nprocs_-1;
    size_t pid=0;
    for (size_t i=slaves;i--;) {
      pid=fork();
      if (pid==-1)
	throw Failed("Fork_controller failed to create slaves");
      if (pid)
	pids(i)=pid;
      else 
       	break; // must be child, so break out
    }
    if (pid==0) { //must  be slave
      pid=getpid();
      bool allset;
      do {
	allset=true;
	size_t i=slaves;
	for (;i--;) {
	  if (pids(i)==pid) {
	    id_=i+1;
	    break;
	  }
	  else {
	    if (pids(i)==pid_t(0))
	      allset=false; //!< flag not all PIDS stored
	  }
	}
      } while (!allset);
      if (id_==0)
	throw Failed("Fork_controller: failed to identify myself");
    }

#ifndef NDEBUG
    char scratch[256];
    sprintf(scratch,"Fork_controllerLOG%i",id_);
    openlog(scratch);
#endif
  }

  bool Fork_controller::join()
  {
    if (!isactive())
      throw Failed("Fork_controller:join called when no processes active");

    writelog("joining");
    closelog();

    if (ammaster()) {
      bool ok=true;
      int status=-1;
      for (size_t i=nprocs_-1;i--;) {
	if (waitpid(pids(i),&status, 0)==-1)
	  perror("waitpid");
	if (verbose_) {
	  char scratch[256];
	  snprintf(scratch,sizeof(scratch),"%" LCM_PRI_SIZE_T_MODIFIER "u exit status: %i",i,status);
	  writelog(scratch);
	}
	if (!WIFEXITED(status) || WEXITSTATUS(status))
	  ok=false;
      }
      pids=0; //!< flag not active
      return ok;
    }
    else {
      exit(0);
    }
  }

  void Fork_controller::join_kill()
  {
    if (isactive())
      join();
    if (ammaster()) { //!< should be redundant
      shddata.kill(); //!< explicitly release shared memory
      shdctrl.kill();
    }
  }

  Warning<> Fork_controller::destroywhileactive_warning("Fork_controller destroyed while still active - processes not cleaned up",&lcm_base_warning);

  Fork_controller::~Fork_controller()
  {
    if (isactive())
      destroywhileactive_warning.raisesafe(); //!< Needs to be exception-safe in destructor, and too late to join
    //      join(); //!< close threads
  }

void Fork_controller::run(const BaseThreadFunction& fobj,size_t totsteps,size_t)
{
  if (!isactive())
    throw Failed("Fork_controller::run: called before controller start called");

  size_t chunksize=totsteps/nprocs_;
  const size_t start=id_*chunksize;
  const size_t end=(id_==nprocs_-1) ? totsteps : (id_+1)*chunksize;

  fobj(start,end,id_); //calculate orientations
  writelog("Finished run");
}

  void shared_memory_obj::zap()
  {
    if (rawptr_) {
      shmdt(rawptr_);
      if (amowner_) {//!< only delete segment if created here
	if (shmctl(shmid_,IPC_RMID,NULL)) { //!< errors are silent by default
#ifndef NDEBUG
	  fprintf(stderr,"Error deleting shared memory segment (ID: %i)\n",shmid_);
	  perror("shmctl");  
#endif
	}
      }
      rawptr_=NULL;
      alloc_=0;
    }
  }

  void* shared_memory_obj::attach(int shmidv,size_t n,size_t alignment)
  {
    if (rawptr_)
      zap();
    return rawattach(shmidv,false,n,alignment);
  }

  Warning<> shared_memory_obj::shmat_warning("Attaching shared memory to process failed",&lcm_base_warning);

  void* shared_memory_obj::rawattach(int shmidv, bool amowner, size_t n, size_t alignment)
  {
    if (rawptr_)
      throw Failed("shared_memory_obj::attach: can't resize");
    void* rawptr=shmat(shmidv,NULL,0);
    if (rawptr==(void*)-1) {
#ifndef NDEBUG
      perror("shmmat");
#endif
      shmat_warning.raise();
      return NULL;
    }
    shmid_=shmidv;    
    amowner_=amowner;
    const ptrdiff_t alignmask= alignment ? (alignment-1) : 0;
    rawptr_=rawptr; //!< Exception-safety: only assign after check
    alignedptr_=alignment ? reinterpret_cast<void*>( reinterpret_cast<ptrdiff_t>(rawptr_) & ~alignmask ) : rawptr_;
    if (alignedptr_==NULL)
      throw InternalError("shared_memory_obj");
    alloc_=n;
    return alignedptr_;
  }

  Warning<> shared_memory_obj::shmget_warning("Request for shared memory failed (try reducing number of data size / processors used).",&lcm_base_warning);
  
  void* shared_memory_obj::ensure(size_t n, size_t alignment)
  {
    ptrdiff_t alignmask=0;
    if (alignment) {
      alignmask=alignment-1;
      if (alignment & alignmask)
	throw InvalidParameter("shared_memory_obj: alignment must be power of 2");
      n+=alignmask;
    }
    if (n<=alloc_)
      return alignedptr_;
    //    if (rawptr_ && wasattached_)
    if (rawptr_)
      throw Failed("Can't resize shared memory block");
    //    zap();
#ifndef NDEBUG
    std::cout << "Creating new shared memory segment\n";
#endif
    const int shmid = shmget(IPC_PRIVATE, n, IPC_CREAT | 0777);
    if (shmid<0) {
      char buf[256];
#ifndef NDEBUG
      perror("shmget");
#endif
      snprintf(buf,sizeof(buf)," Claim: %" LCM_PRI_SIZE_T_MODIFIER "u bytes",n);
      shmget_warning.raise(buf);
      return NULL;
    }
#ifndef NDEBUG
    fprintf(stderr,"Created shared memory segment (ID: %i)\n",shmid);
#endif
    return rawattach(shmid,true,n,alignment);
  }

  void Fork_controller::poll()
  {
    static timespec tspec;
    tspec.tv_sec=0;
    tspec.tv_nsec=pollns;

    nanosleep(&tspec,NULL);
  }

  void* Fork_controller::ensuredata(size_t totbytes, size_t alignment)
  {
    if (totbytes==0)
       throw InvalidParameter("Fork_controller::ensuredata");    

    if ((shddata()!=NULL) && (totbytes==databytes_) && (alignment_==alignment))
      return shddata();

    volatile int* shdkeyp=reinterpret_cast<volatile int*>(shdctrl());
    void* shddatap=NULL;

    if (ammaster()) {
       shddatap=shddata.ensure(totbytes,alignment);  //!< creates shared memory segment (will fail if request exceeds currently allocated)
       *shdkeyp=shddata.shmid();
    }
     else {
       while (*shdkeyp==-1) //ensure shared memory segment has been created
	 poll();
       shddatap=shddata.attach(*shdkeyp,totbytes,alignment); //!< may return NULL if shmget failed
    }
    if (shddatap) {
      databytes_=totbytes;
      alignment_=alignment;
    }
    else
      databytes_=0;

    return shddatap;
  }

}//namespace libcmatrix
