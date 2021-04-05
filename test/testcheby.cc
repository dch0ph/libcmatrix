/* Test of Chebyshev propagation */

#include "NMR.h"
#include "BlockedMatrix.h"
#include "PartitionedMatrix.h"
#include "timer.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

#define TYPE complex

int main(int argc, const char* argv[])
{
  // Matrix<TYPE> As(3,3,ScratchList<TYPE>(TYPE(20.0), TYPE(5.0), TYPE(0.0),
  //				   TYPE(5.0), TYPE(30.0), TYPE(10.0),
  //					   TYPE(0.0), TYPE(10.0), TYPE(40.0)));
  Matrix<TYPE> As(4,4,TYPE(0.0));
  As(0,1)=As(1,0)=As(2,3)=As(3,2)=27250;
  As(0,2)=As(2,0)=As(1,3)=As(3,1)=31250;
					
  //   As(0,0)=20.0;
//   As(1,1)=30.0;
//   As(2,2)=40.0;
//   // As(3,3)=2.0;
//   As(0,1)=As(1,0)=5.0;
//   As(1,2)=As(2,1)=10.0;

  int count=1;
  //const size_t nrep=1;
  const size_t nrep=getint(argc,argv,count,"multiply factor? ",2);
  const size_t ntimes=getint(argc,argv,count,"repeat times? ",1);
  const double dt=getfloat(argc,argv,count,"dt? ",1e-3);

    Matrix<TYPE> A(As);
//   const size_t Arows=As.rows();
//   const size_t nrows=Arows*nrep;
//   rmatrix A(nrows,nrows,0.0);
//   const rmatrix tAs(transpose(As));
//   for (size_t i=nrep-1;i--;) {
//     const range rsel(Arows*i,Arows*i+Arows-1);
//     const range csel(Arows*(i+1),Arows*(i+2)-1);
//     A(rsel,csel)=As;
//     A(csel,rsel)=tAs;
//   }
  if (A.rows()<10)
    cout << "A\n" << A << '\n'; 
  else
    spy(A);

  //  const List<size_t> div(nrep,Arows,mxflag::normal);
  //const List<size_t> div(ExplicitList<2,size_t>(2U,1U));
  //cout << "Selection: " << div << '\n';
  //matrix_partition part(div,div,A);
  //cout << "Partition:\n" << part << '\n';
  //  PartitionedMatrix<TYPE> Ablocked(div,div,A);
  // cout << "A blocked\n" << Ablocked << '\n';

//   SubMatrix<double> Asub(A,range(0,1),range(1,2));
//   cout << "Asub\n" << Asub << '\n'; 
//   Asub+=5;

//   cout << "A\n" << A << '\n';


  cmatrix Ustart(identity(A.rows()));
  cmatrix U,Utmp;

  timer<> stopwatch;

  cout.precision(8);
  cmatrix_eigensystem_controller.tolerance=1e-6;
  cmatrix_eigensystem_controller.verbose=(A.rows()<20) ? 2 : 1;

  Matrix<TYPE> V;
  List<double> eigs;
  hermitian_eigensystem(V,eigs,A);
  cout << "Max eig: " << max(eigs) << "  Min eig: " << min(eigs) << '\n';

   for (size_t n=ntimes;n--;) {
     Utmp=Ustart;
     propagator(Utmp,A,dt);
     multiply(U,Utmp,Ustart);
   }
   cout << "U (conventional)\n";
   if (nrep<3)
     cout << Utmp;
   if (ntimes>1)
     cout << "Time per step: " << (stopwatch()*1e6/ntimes) << " us\n";
  
  stopwatch.reset();
  //for (size_t N=1;N<=10;N++) {
  size_t N=7;
  cmatrix_eigensystem_controller.chebyshev_iterations=N;
    for (size_t n=ntimes;n--;) {
      //      U=Ustart;
      chebyshev_propagator(U,A,dt,NULL);//&part);
      //chebyshev_propagate(U,Ablocked,dt,N);
    }
    cout << "U (Chebyshev " << N << ")\n";
    if (nrep<3)
      cout << U;
    if (ntimes>1)
      cout << "Time per step: " << (stopwatch()*1e6/ntimes) << " us\n";
    cout << '\n';
    //}
  
  return 0;
}
