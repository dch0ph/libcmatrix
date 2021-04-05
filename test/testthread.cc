#include <unistd.h>
#include "cmatrix_threads.h"
#include "timer.h"
#include "Matrix.h"

using namespace std;
using namespace libcmatrix;

void thread_func(size_t start,size_t end,size_t nthread)
{
  cout << "(" << nthread << "): ";
  counts(nthread)++;
  for (size_t j=start;j<end;j++) {
    cout << j << " ";
    sleep(1);
  }
}

const int NTHREADS=4; // use 4 threads

List<int> counts(NTHREADS,0);

void thread_func2(size_t,size_t,size_t nthread)
{
  cout << "(" << nthread << "): ";
  counts(nthread)++;
  volatile double b=2.0;
  double a;
  for (size_t n=5;n--;) {
    Matrix<double> tmp(5,5); //trigger creation 
  }
  matrix_traits<double>::allocator::print(std::cout);
  for (size_t n=5000000;n--;) a=b;
  //sleep(1);
  //cout << "thread_func2: ";
}

int main()
{
  try {

    thread_controller tcon(NTHREADS);    
    
    timer<CPUTimer> stopwatch;
    timer<WallTimer> wstopwatch;
    tcon.run(ThreadFunction(thread_func),40,5);
    cout << endl;

    cout << "CPU time: " << stopwatch() << endl;
    cout << "Wall-clock time: " << wstopwatch() << endl;
    cout << "Times each thread called: " << counts << endl;
  }
  catch (MatrixException &exc) {
    cerr << exc << endl;
  }
  return 0;
}
