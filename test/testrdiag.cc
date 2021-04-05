#include <stdlib.h>
#include "cmatrix_utils.h"
//#include "diagonal.h"
#include "NMR.h"
#include "ScratchList.h"
//#include "hacked_sunperf.h"

using namespace libcmatrix;
using namespace std;

extern "C" {
  extern void dsyevx_(char,char,char,const int *,double *,const int *, const double *, const double *, const int *, const int *, const double*, const int *, double *,double *, const int *,double *,const int *,int *,int *,int *);
		 }

void lapack_eigensystem(rmatrix &V,List<double> &evals,rmatrix a)
{
  if (!issquare(a))
    throw NotSquare("lapack_eigensystem");
  const int n=a.rows();
  
  evals.create(n);
  V.create(n,n);

  int nfound,info,ifail;

  const int lwork=8*n;
  ScratchList<double,SCRATCH_SIZE> work(lwork);
  ScratchList<int,SCRATCH_SIZE> iwork(5*n);

  //  dsyevx('V','A','U',n,a.vector(),n,0.0,0.0,0,0,0.0,&nfound,evals.vector(),V.vector(),n,&ifail,&info);
  
  double vl=0.0;
  double vu=0.0;
  int il=0;
  int iu=0;
  const double abstol=0.0;

  dsyevx_('V','A','U',&n,a.vector(),&n,&vl,&vu,&il,&iu,&abstol,&nfound,evals.vector(),V.vector(),&n,work.vector(),&lwork,iwork.vector(),&ifail,&info);

  if (info) {
    cerr << "dsyevx failed: " << info << endl;
    exit(1);
  }
  cerr << "Optimal work: " << work(0) << endl;
}

int main()
{

  //spin_system ax(2,"2H");
  //List<double> fred=diag_spin_quadrupolar(ax,1);
  //cout << fred << endl;
  
  int i;
  const int n=4;
  const int reptimes=1;

  rmatrix a(n,n);
  
  for (i=0;i<n;i++) {
    for (int j=0;j<=i;j++)
      a(i,j)=(a(j,i)=random(1.0));
  }

  //cmatrix ca(a);
  //cmatrix VC;

  cout.precision(3);
  cout << "A:\n" << a << endl;

  rmatrix V;
  List<double> evals;

  for (i=reptimes;i--;)
    lapack_eigensystem(V,evals,a);
    //propagator(VC,a,1e-6);
    //hermitian_eigensystem(VC,evals,ca);

  cout << "Eigenvalues: " << evals << endl;
  //cout << VC << endl;
  
  //V.transpose();
  cout << "Eigenvectors:\n" << V << endl;

  cout << "Unitary check:\n" << transpose(V)*V << endl;

  rmatrix d;
  unitary_simtrans(d,evals,V);
  //  cout << "Check:\n" << V*evals*transpose(V) << endl;
  cout << "Check:\n" << d << endl;
}
