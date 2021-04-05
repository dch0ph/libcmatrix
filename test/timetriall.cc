/* Time a load of operations */

#include "cmatrix_complex.h"
#include "List.h"
#include "timer.h"
#include "ttyio.h"

using namespace libcmatrix;

timer<> stopwatch;

void reset()
{
  stopwatch.reset();
}

int UNITST;

void print_usage(int ops)
{
  const double t=stopwatch();
  std::cout << (t*1e6/UNITST) << " us   " << (ops/(stopwatch()*1e6)) << " Megaops/s" << std::endl;
}

void dop(const List<double> &l) { std::cout << "List<double>: " << l << std::endl; }

#define TYPE complex

int main(int argc,const char *argv[])
{
  int i;
  int count=1;
  const double Mflops=getfloat(argc,argv,count,"Megaops (0 for 1 iter)? ",100);
  const int len=getint(argc,argv,count,"List length? ",32);
  UNITST=Mflops ? int(Mflops*1e6/len) : 1;
  const bool verbose=(len<64);
  List<TYPE> A(len);
  List<TYPE> B(len);
  List<TYPE> C(len);

  std::cout << "Simple type: " << (type_traits<TYPE>::trivialconstructor ? "Yes" : "No") << std::endl;

  for (i=len;i--;) {
    B(i)=i;
    C(i)=i;
  }

  const int ops=UNITST*len;
  std::cout << UNITST << " repeats\n";

//   std::cout << "\nInit A(B): ";
//   reset();
//   for (i=UNITST;i--;) List<TYPE> D(B);
//   print_usage(ops);

//   std::cout << "\nCopy A=B: ";
//   reset();
//   for (i=UNITST;i--;) A=B;
//   print_usage(ops);
//   if (verbose) std::cout << A << std::endl;

//   std::cout << "\nAddition (A=B+C): ";
//   reset();
//   for (i=UNITST;i--;) A=B+C;
//   print_usage(ops);
//   if (verbose) std::cout << A << std::endl;

//   std::cout << "\nCopy add (A=B; A+=C): ";
//   reset();
//   for (i=UNITST;i--;) {
//     A=B;
//     A+=C;
//   }
//   print_usage(ops);
//   if (verbose) std::cout << A << std::endl;
  
//   std::cout << "\nCopy subtract (A=B; A-=C): ";
//   reset();
//   for (i=UNITST;i--;) {
//     A=B;
//     A-=C;
//   }
//   print_usage(ops);
//   if (verbose) std::cout << A << std::endl;
  
  std::cout << "\nSupplied-destination add(A,B,C): ";
  reset();
  for (i=UNITST;i--;) add(A,B,C);
  print_usage(ops);
  if (verbose) std::cout << A << std::endl;

  const TYPE c(3.0);

  std::cout << "\nAdd constant (A+=c): ";
  A=B;
  reset();
  for (i=UNITST;i--;) A+=c;
  print_usage(ops);
  if (verbose) std::cout << A << std::endl;
  
//   std::cout << "\nSubtract constant (A-=c): ";
//   A=B;
//   reset();
//   for (i=UNITST;i--;) A-=c;
//   print_usage(ops);
//   if (verbose) std::cout << A << std::endl;

//   std::cout << "\nAdd constant (A=B+c): ";
//   reset();
//   for (i=UNITST;i--;) A=B+c;
//   print_usage(ops);
//   if (verbose) std::cout << A << std::endl;
  
  A=TYPE(0);
  std::cout << "\nmla: ";
  reset();
  for (i=UNITST;i--;) mla(A,B,C);
  print_usage(ops);
  if (verbose) std::cout << A << std::endl;

  return 0;
}
