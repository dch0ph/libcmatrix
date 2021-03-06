#ifndef lcm_cmatrix_MPI_h_
#define lcm_cmatrix_MPI_h_

//Simplifies matters!
#ifdef LCM_REENTRANT
#error "libcmatrix can't combine MPI with multi-threading."
#endif

#include "lcm_basethreads.h"
#include "Warnings.h"
#include "mpi.h"

namespace libcmatrix {

  // Type can't be created non-const in some MPI implementations (OPENMPI)
  template<typename T> struct MPI_type;
  template<> struct MPI_type<libcmatrix::complex> {
    static MPI_Datatype type;
  };
  template<> struct MPI_type<double> {
    static MPI_Datatype type;
  };
  template<> struct MPI_type<int> {
    static MPI_Datatype type;
  };
 
  class MPI_controller : public base_parallel_controller {
public:
    MPI_controller(int& argc, char**& argv, int verbose =LCM_PARALLEL_DEFAULT_VERBOSE);
    ~MPI_controller();
  
  void run(const BaseThreadFunction&, size_t, size_t);

  template<typename T> void sum(BaseList<T> dest, const BaseList<T>& source) const;
  template<typename T> void sum(T& dest, const T& source) const {
    duplicate_structure(dest,source);
    sum(dest.row(),source.row());
  }

    template<typename T> void broadcast(BaseList<T>) const;
    void sync() const { MPI::COMM_WORLD.Barrier(); }

    //  static Warning<> called_by_slave_warning;
    
    static Warning<> notallstarted_warning; //!< warn if not all MPI processes started successfully
    static Warning<> tooshort_warning; //!< warn if job completes too quickly

private:
  static const int master=0;
    int jobcount;

  unsigned int buffer[2];
  MPI::Status status;
  char scratch[256];

  //only used by master process
  size_t current;
  size_t todo;
  size_t notdone;
  size_t chunk;
  int num_workers;

  bool sendchunk(int dest); 
};

  // implementation details below here

  template<typename T> void MPI_controller::sum(BaseList<T> dest, const BaseList<T>& source) const {
    //if (ammaster()) {
      //      called_by_slave_warning.raise();
    const size_t n(dest.size());
    if (n==0)
	throw Undefined("MPI_controller::sum");
    if (n!=source.size())
	throw Mismatch("MPI_controller::sum");
    MPI::COMM_WORLD.Reduce(source.vector(),dest.vector(),n,MPI_type<T>::type,MPI::SUM,0);
  }

  template<typename T> void MPI_controller::broadcast(BaseList<T> buf) const 
  {
    MPI::COMM_WORLD.Bcast(buf.vector(),buf.size(),MPI_type<T>::type,0);
  }

} //namespace libcmatrix

#endif
