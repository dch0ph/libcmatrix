/* test MPI programming with libcmatrix */

#include "timer.h"
#include "List.h"
#include <iostream>
#include "cmatrix_MPI.h"

libcmatrix::List<double> sumvec(1,0.0);

using namespace std;
using namespace libcmatrix;

void thread_func(size_t start,size_t end,size_t nthread)
{
  cout << "(" << nthread << "): ";
  //  counts(nthread)++;
  for (size_t j=start;j<end;j++) {
    cout << j << " ";
    sumvec.front()+=j;
    sleep(1);
  }
}

int main(int argc, char* argv[])
{
  try {

    libcmatrix::MPI_controller tcon(argc,argv,1);

    timer<CPUTimer> stopwatch;
    timer<WallTimer> wstopwatch;
    //cout << "Started" << endl;

    const size_t n=32;

    tcon.run(ThreadFunction(thread_func),n,5);
    cout << endl;
    
    libcmatrix::List<double> destvec;
    tcon.sum(destvec,sumvec);
    if (tcon.ammaster()) {
      cout << "Sum: " << destvec.front() << "   (Expect " << ((n*(n-1)/2)) << ")\n";
      //cout << "CPU time: " << stopwatch() << endl;
      cout << "Wall-clock time: " << wstopwatch() << " s" << endl;
    //    cout << "Times each processor called: " << counts << endl;
    }
  }
  catch (MatrixException& exc) {
    cerr << exc << endl;
  }

  return 0;
}
