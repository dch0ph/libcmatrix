#include "cmatrix_MPI.h"
#include "List.h"

namespace libcmatrix {

  static const double mintime=30; //!< warn if job completes in less than X seconds

  //const MPI_Datatype MPI_type<double>::type;
  //  const MPI_Datatype MPI_type<libcmatrix::complex>::type;
  //  const MPI_Datatype MPI_type<int>::type;

  MPI_Datatype MPI_type<double>::type=MPI_DOUBLE;
  MPI_Datatype MPI_type<libcmatrix::complex>::type=MPI_DOUBLE_COMPLEX;
  MPI_Datatype MPI_type<int>::type=MPI_INT;

  MPI_controller::MPI_controller(int& argc, char**& argv, int verbosev)
    : base_parallel_controller(verbosev),
    jobcount(0)
  {
    MPI::Init(argc,argv);
    id_=MPI::COMM_WORLD.Get_rank();
    nprocs_=MPI::COMM_WORLD.Get_size();
    nworkers_=nprocs_-1;
#ifndef NDEBUG
    sprintf(scratch,"MPI_controllerLOG%i",id_);
    openlog(scratch);
#endif
    if (verbose_) {
      char scratch[MPI_MAX_PROCESSOR_NAME+64];
      char name[MPI_MAX_PROCESSOR_NAME];
      int len;
      MPI_Get_processor_name(name,&len);
      snprintf(scratch,sizeof(scratch),"Created process %d of %" LCM_PRI_SIZE_T_MODIFIER "u on %s\n",id_,nprocs_,name);
      writelog(scratch);
    }
  }

  MPI_controller::~MPI_controller()
  {
    MPI::Finalize();
  } //clean up

  //  Warning<> MPI_controller::called_by_slave_warning("sum called by slave",&lcm_base_warning);

  bool MPI_controller::sendchunk(int dest)
  {
    if (this->notdone==0)
      return false;

    const size_t this_chunk = (current+chunk>todo) ? todo-current : chunk;
    buffer[0]=current;
    buffer[1]=current+this_chunk;
    if (verbose_) {
      sprintf(scratch,"Sending steps %u to %u to worker %i",buffer[0],buffer[1],dest);
      writelog(scratch);
    }
    //    const int retcode=MPI::COMM_WORLD.Send(buffer,2,MPI::UNSIGNED,dest,jobcount);
    const int retcode=MPI_Send(buffer,2,MPI::UNSIGNED,dest,jobcount,MPI_COMM_WORLD);
    if (retcode!=MPI_SUCCESS) { //!< probably not useful - doesn't mean message received
      sprintf(scratch,"Failed to send work to worker %i",dest);
      writelog(scratch);
      return false;
    }
      
    notdone-=this_chunk;
    current+=this_chunk;
    return true;
  }

  Warning<> MPI_controller::notallstarted_warning("Not all MPI workers successfully initiated - workers either failed to start or job was over too quickly",&lcm_base_warning);
  Warning<> MPI_controller::tooshort_warning("Calculation will be much more efficient run conventionally using -disable:parallel. ",&lcm_base_warning);

void MPI_controller::run(const BaseThreadFunction& fobj,size_t orients,size_t rchunk)
{
  jobcount++; //inc job counter
  
  if (ammaster()) {
    current=0;
    todo=notdone=orients;
    chunk=rchunk;
    num_workers=0;
    size_t count_workers=0;

    // for (int i=1;i<=nprocs_-1;i++) {
    //   if (this->notdone) {
    // 	if (sendchunk(i))
    // 	  num_workers++;
    // 	else
    // 	  break;
    //   }
    // }
    List<bool> started(nprocs_,false); //!< could be nprocs_-1, but no point in unnecessarily risking fence-post error
    const double time_start=time();
    
    while ( (this->notdone) || (num_workers) || ((time()-time_start)<mintime)) {
      writelog("Waiting for request");
      MPI::COMM_WORLD.Recv(buffer,0,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
      const int dest=status.Get_source();
      const bool newworker=(started(dest)==false);
      if (verbose_) {
	snprintf(scratch,sizeof(scratch),"Received request from %i (%s worker)",dest, (newworker ? "new" : "existing"));
	writelog(scratch);
      }
      if (newworker) {
	started(dest)=true;
	count_workers++;
	num_workers++;
      }
      if (status.Get_tag()!=jobcount)
	error("MPI_controller: out of sequence packet (master)");
      if (this->notdone)
	sendchunk(dest);
      else {
	num_workers--;
	if (verbose_) {
	  sprintf(scratch,"Sending job done to %i",dest);
	  writelog(scratch);
	}
	MPI::COMM_WORLD.Send(buffer,0,MPI::UNSIGNED,dest,0);
      }
    }
    const double timetaken=time()-time_start;
    if (timetaken<mintime) {
      sprintf(scratch,"Job finished in only %g seconds!",timetaken);
      tooshort_warning.raise(scratch);
    }
    if (count_workers!=nprocs_-1) {
      sprintf(scratch," (%" LCM_PRI_SIZE_T_MODIFIER "u out of %" LCM_PRI_SIZE_T_MODIFIER "u)",count_workers,nprocs_-1);
      notallstarted_warning.raise(scratch);
    }
  }
  else { //worker?
    for (;;) {
      writelog("Asking for instructions");
      MPI::COMM_WORLD.Send(buffer,0,MPI::UNSIGNED,master,jobcount);
      writelog("Waiting for instructions");
      MPI::COMM_WORLD.Recv(buffer,2,MPI::UNSIGNED,master,MPI::ANY_TAG,status);
      const int tag=status.Get_tag();
      if (tag==0)
	break;
      if (verbose_) {
	sprintf(scratch,"Worker %i running %u to %u" ,get_thread_num(),buffer[0],buffer[1]);
	writelog(scratch);
      }
      if (tag!=jobcount)
	error("MPI_controller: out of sequence packet (worker)");
      fobj(buffer[0],buffer[1],get_thread_num());
    }
  }
  writelog("Finished run");
}

}//namespace libcmatrix
