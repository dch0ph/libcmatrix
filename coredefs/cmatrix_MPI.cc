#include "cmatrix_MPI.h"

namespace libcmatrix {

  const MPI_Datatype MPI_type<double>::type;
  const MPI_Datatype MPI_type<libcmatrix::complex>::type;
  const MPI_Datatype MPI_type<int>::type;

  MPI_controller::MPI_controller(int& argc, char**& argv)
  {
    MPI::Init(argc,argv);
    id_=MPI::COMM_WORLD.Get_rank();
    nprocs_=MPI::COMM_WORLD.Get_size();
    jobcount=0;
#ifndef NDEBUG
    sprintf(scratch,"MPI_controllerLOG%i",id_);
    log=fopen(scratch,"w");
#else
    log=NULL;
#endif
  }

  MPI_controller::~MPI_controller()
  {
    MPI::Finalize(); //resource leak if this throws
    if (log)
      fclose(log);
  } //clean up

  Warning<> MPI_controller::called_by_slave_warning("sum called by slave",&lcm_base_warning);

  bool MPI_controller::sendchunk(int dest)
  {
    if (this->notdone==0)
      return false;

    const size_t this_chunk = (current+chunk>todo) ? todo-current : chunk;
    buffer[0]=current;
    buffer[1]=current+this_chunk;
#ifndef NDEBUG
    sprintf(scratch,"Sending %u to %u to %i",buffer[0],buffer[1],dest);
    writelog(scratch);
#endif
    MPI::COMM_WORLD.Send(buffer,2,MPI::UNSIGNED,dest,jobcount);
    notdone-=this_chunk;
    current+=this_chunk;
    return true;
  }

void MPI_controller::run(const BaseThreadFunction& fobj,size_t orients,size_t rchunk)
{
  jobcount++; //inc job counter
  
  if (ammaster()) {
    current=0;
    todo=notdone=orients;
    chunk=rchunk;
    num_workers=0;

    for (int i=1;i<=nprocs_-1;i++) {
      if (sendchunk(i))
	num_workers++;
      else
	break;
    }

    while (num_workers) {
      writelog("Waiting for result packet");
      MPI::COMM_WORLD.Recv(buffer,0,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
      const int dest=status.Get_source();
#ifndef NDEBUG
      sprintf(scratch,"Received packet from %i",dest);
      writelog(scratch);
#endif
      
      if (status.Get_tag()!=jobcount)
	error("MPI_controller: out of sequence packet (master)");
      if (!sendchunk(dest)) {
	num_workers--;
#ifndef NDEBUG
	sprintf(scratch,"Sending job done to %i",dest);
	writelog(scratch);
#endif
	MPI::COMM_WORLD.Send(buffer,0,MPI::UNSIGNED,dest,0);
      }
    }
    
  }
  else { //worker?
    for (;;) {
      writelog("Waiting for instructions");
      MPI::COMM_WORLD.Recv(buffer,2,MPI::UNSIGNED,master,MPI::ANY_TAG,status);
      const int tag=status.Get_tag();
      if (tag==0)
	break;
#ifndef NDEBUG
      sprintf(scratch,"Worker %i running %u to %u" ,get_thread_num(),buffer[0],buffer[1]);
      writelog(scratch);
#endif
      if (tag!=jobcount)
	error("MPI_controller: out of sequence packet (worker)");
      fobj(buffer[0],buffer[1],get_thread_num());
      MPI::COMM_WORLD.Send(buffer,0,MPI::UNSIGNED,master,jobcount);
    }
  }
  writelog("Finished run");
}

}//namespace libcmatrix
