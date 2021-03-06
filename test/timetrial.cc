/* Time a load of operations */

#include <cstring>
#include <cstdlib>
#include "NMR.h"
#include "timer.h"
#include "ttyio.h"

using namespace libcmatrix;

timer<> stopwatch;

void reset()
{
  //  reset_stats();
  stopwatch.reset();
}

void print_usage(int times,int simple =0)
{
  std::cout << "Time taken: " << (1e6*stopwatch()/times) << " us\n";
  if (!simple)
    matrix_traits<complex>::allocator::print(std::cout);
}

int main(int argc,const char* argv[])
{
#ifdef NDEBUG
  try {
#endif
    //    cmatrix_resize_stack(0); //this will disable cmatrix stacking

    int count=1;

  spin_system ax(3);

  cmatrix H01=10.0*spin_dipolar(ax,0,1)+15.0*spin_dipolar(ax,1,2);
  cmatrix H02=100.0*spin_dipolar(ax,0,2);

  const slice cut1(0,4);
  const slice cut2(4,4);
    
  const double Mflops=getfloat(argc,argv,count,"Mflops (0 for 1 iter)? ",100);

  const int n=H01.rows();
  const int UNIT=Mflops ? int(Mflops*1e6/(2*n*n)) : 1;
  const int UNITST=Mflops ? int (Mflops*1e6/(6*n*n*n)) : 1;
  const int COMPLEX_ADDS=UNIT*n*n;

  int i;

//   reset();
//   for (i=UNIT;i--;) {
//     cmatrix tmp=H01+H02;
//   }
//   std::cout << "Copy: "; print_usage();

//   reset();
//   for (i=UNIT;i--;) {
//     cmatrix tmp(H01+H02);
//   }
//   std::cout << "Init: "; print_usage();

  std::cout << "Complex number operations: " << COMPLEX_ADDS << "\n\n";

  //complex a(0.0,0.0);
  //complex b(2.0,3.0);
  //complex c(4.0,3.0);

  //These tests are misleading since maths may be optimised away!

//   reset();
//   for (i=COMPLEX_ADDS;i--;) a=b+c;
//   print_usage(UNIT,1);
//   std::cout << "Complex add a=b+c: " << a << std::endl << std::endl;

//   reset();
//   for (i=COMPLEX_ADDS;i--;) {
//     a=b;
//     a+=c;
//   }
//   print_usage(UNIT,1);
//   std::cout << "Copy add (a=b; a+=c): " << a << std::endl << std::endl;

//   reset();
//   for (i=COMPLEX_ADDS;i--;) a=b*c;
//   print_usage(UNIT,1);
//   std::cout << "Multiply a=b*c: " << a << std::endl << std::endl;

//   std::cout << "expi: " << UNIT << std::endl;
//   reset();
//   for (i=UNIT;i--;) expi(1.0);
//   print_usage(UNIT,1);
//   std::cout << "expi(b): " << expi(1.0) << std::endl;
  
//   std::cout << "\nmla:\n";
//   a=0.0;
//   reset();
//   for (i=COMPLEX_ADDS;i--;) mla(a,b,c);
//   print_usage(UNIT,1);
//   std::cout << "mla a+=b*c: " << a << std::endl << std::endl;

  std::cout.precision(2);

  cmatrix A;
  List<double> vec(n,2.0);

  std::cout << "\nCopy (A=B):" << std::endl;
  reset();
  for (i=UNIT;i--;)
    A=H01;
  print_usage(UNIT);
  std::cout << A << std::endl;

  std::cout << "Matrix addition: " << UNIT << std::endl;
  reset();
  for (i=UNIT;i--;) A=H01+H02;
  print_usage(UNIT);
  std::cout << A << std::endl;

  std::cout << "\nCopy and add (A=B; A+=C):" << std::endl;
  reset();
  for (i=UNIT;i--;) {
    A=H01;
    A+=H02;
  }
  print_usage(UNIT);
  std::cout << A << std::endl;

//   std::cout << "\nA+=C:" << std::endl;
//   reset();
//   for (i=UNIT;i--;)
//     _Apply_ip2<2,2>::add(A,H02);
//   print_usage(UNIT);

  std::cout << "\nSupplied-dest addition (add(A,B,C)):" << std::endl;
  reset();
  for (i=UNIT;i--;) add(A,H01,H02);
  print_usage(UNIT);
  std::cout << A << std::endl;

  const rmatrix rB=real(H01);
  rmatrix rA;
  std::cout << "\nMixed arithmetic I: A=real(B)+C:" << std::endl;
  reset();
  for (i=UNIT;i--;) A=rB+H02;
  print_usage(UNIT);
  std::cout << A << std::endl;

  std::cout << "\nMixed arithmetic II: add(A,real(B),C):" << std::endl;
  reset();
  for (i=UNIT;i--;) add(A,rB,H02);
  print_usage(UNIT);
  std::cout << A << std::endl;
  
//   std::cout << "\nFunctional Supplied-dest addition:" << std::endl;
//   reset();
//   for (i=UNIT;i--;) addf(A,H01,H02);
//   print_usage(UNIT);
//   std::cout << A << std::endl;

  std::cout << "\nMultiplication (A=B*C):" << std::endl;
  reset();
  for (i=UNITST;i--;) A=H01*H02;
  print_usage(UNITST);
  std::cout << A << std::endl;

  std::cout << "\nSupplied dest multiplication (multiply(A,B,C)):" << std::endl;
  reset();
  for (i=UNITST;i--;) multiply(A,H01,H02);
  print_usage(UNITST);
  std::cout << A << std::endl;

  std::cout << "\nmultiply and add (A+=2.0*B): " << std::endl;
  A.kill();
  reset();
  for (i=UNIT;i--;) A+=2.0*H01;
  print_usage(UNIT);
  std::cout << A << std::endl;

  std::cout << "\nMLA (mla(A,2.0,B)): " << std::endl;
  A.kill();
  reset();
  for (i=UNIT;i--;) mla(A,2.0,H01);
  print_usage(UNIT);
  std::cout << A << std::endl;

  std::cout << "\nMLA (mla(A,complex(2.0),B): " << std::endl;
  A.kill();
  const complex c2(2.0,0.0);
  reset();
  for (i=UNIT;i--;) mla(A,c2,H01);
  print_usage(UNIT);
  std::cout << A << std::endl;

  List<cmatrix> HB01(2);
  HB01(0)=H01(cut1,cut1); HB01(1)=H01(cut2,cut2);

  std::cout << "\n'Blocked' H01:\n" << HB01 << std::endl;

  List<cmatrix> HB02(2);
  HB02(0)=H02(cut1,cut1); HB02(1)=H02(cut2,cut2);

  List<cmatrix> T01(2);

//   std::cout << "\nBlock addition (A=B+C): " << std::endl;
//   reset();
//   for (i=UNIT;i--;)
//     T01=HB01+HB02;
//   print_usage(UNIT);
//   std::cout << T01 << std::endl;
  
  std::cout << "\nblock copy and add addition (A=B; A+=C):" << std::endl;
  reset();
  for (i=UNIT;i--;) {
    T01=HB01;
    T01+=HB02;
  }
  print_usage(UNIT);
  std::cout << T01 << std::endl;

  std::cout << "\nexplicit block supplied-dest. add (add(A,B,C)):" << std::endl;
  reset();
  for (i=UNIT;i--;) { add(T01(0),HB01(0),HB02(0)); add(T01(1),HB01(1),HB02(1)); }
  print_usage(UNIT);
  std::cout << T01 << std::endl;

//   std::cout << "\nblock supplied-dest. subtract (subtract(A,B,C)):" << std::endl;
//   reset();
//   for (i=UNIT;i--;) subtract(T01,HB01,HB02);
//   print_usage(UNIT);
//   std::cout << T01 << std::endl;

//   std::cout << "\nblock functional supplied-dest. add:" << std::endl;
//   reset();
//   for (i=UNIT;i--;) addf(T01,HB01,HB02);
//   print_usage(UNIT);
//   std::cout << T01 << std::endl;

  std::cout << "\nSimilarity transform: " << UNITST << std::endl;
  std::cout << "Before:\n" << H01 << std::endl;

  cmatrix B=propagator(H02,1e-3);
  std::cout << B << std::endl;
  List<cmatrix> BB(2);
  BB(0)=propagator(HB02(0),1e-3);
  BB(1)=propagator(HB02(1),1e-3);

  reset();
  for (i=UNITST;i--;) unitary_simtrans(A,H01,B);
  print_usage(UNITST);
  std::cout << A << std::endl;

  std::cout << "\nSimilarity itransform: " << UNITST << std::endl;
  reset();
  for (i=UNITST;i--;) unitary_isimtrans(A,H01,B);
  print_usage(UNITST);

  std::cout << A << std::endl;

  std::cout << "\nblocked similarity itransform: " << UNITST << std::endl;
  reset();
  for (i=UNITST;i--;) {
    unitary_isimtrans(T01(0),HB01(0),BB(0));
    unitary_isimtrans(T01(1),HB01(1),BB(1));
  }
  print_usage(UNITST);

  std::cout << T01 << std::endl;

#ifdef NDEBUG
  }
  catch (MatrixException &exc) {
    std::cerr << exc << std::endl;
  }
#endif
  return 0;
}
