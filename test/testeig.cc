#ifdef LCM_USE_EXTERNAL
#include "config.h"
#include "cmatrix_external.h"
#endif
#include "cmatrix.h"
#include "cmatrix_utils.h"
#include "matlabio.h"
#include "timer.h"
#include "ScratchList.h"
#include "ttyio.h"

using namespace std;
using namespace libcmatrix;

double defzero=1e-10;

//!< hidden declarations
namespace libcmatrix {
  void lapack_eigensystem(cmatrix& V, BaseList<complex> eigs, const cmatrix& a0, cmatrix* wspace =NULL);
  void internal_eigensystem(cmatrix& V, BaseList<complex> eigs, const cmatrix& a0, cmatrix* wspace =NULL);
}

template<class T> void make_block(Matrix<T>& a,double zeroel)
{
  const size_t n(a.rows());
  Matrix<T> tmp(2*n,2*n,0.0);
  if (zeroel) {
    zeroel*=2.0;
    BaseList<T> tmpr(tmp.row());
    for (size_t i=tmp.size();i--;) 
      tmpr(i)=zeroel*(gauss()-0.5);
  }
  const int base=2*n-1;
  for (size_t r=n;r--;) 
    for (size_t c=n;c--;)
      tmp(r,c)=tmp(base-r,base-c)=a(r,c);
  a.swap(tmp);
}

const int bigsize=4;  //must be even for testing block diagonal

inline rmatrix conj_transpose(const rmatrix &a) { return transpose(a); }
//inline void conj_transpose_multiply(rmatrix& a, const rmatrix& b, const rmatrix& c) { transpose_multiply(a,b,c); }

timer<> stopwatch;

template<class T,class T2> void print_tests(T2 &D,const BaseList<T> &eig,const T2 &start,size_t times)
{
  if (times>1)
    cout << "Time per matrix (us): " << stopwatch()*1e6/times << endl;
  cout << "Initial matrix:\n" << start << endl;

  cout << "Eigenvalues: " << eig << endl;    
  cout << "Eigenvectors\n" << D << endl;
  cout << "Unitarity of eigenvectors\n" << (conj_transpose(D)*D) << endl;

  T2 tmp;
  unitary_simtrans(tmp,eig,D);
  cout << "Regenerated matrix\n" << tmp << endl;
  cout << "Difference from original matrix\n" << (tmp-start) << endl;

  cout << endl;
}  

template<class T,class T2> void print_tests_gen(T2 &D,const BaseList<T> &eig,const T2 &start,size_t times)
{
  if (times>1)
    cout << "Time per matrix (us): " << stopwatch()*1e6/times << endl;
  cout << "Initial matrix:\n" << start << endl;

  cout << "Eigenvalues: " << eig << endl;    
  cout << "Eigenvectors (non-unitary)\n" << D << endl;

  T2 tmp;
  simtrans(tmp,eig,D);
  cout << "Regenerated matrix\n" << tmp << endl;
  cout << "Difference from original matrix\n" << (tmp-start) << endl;

  cout << endl;
}  

int main(int argc,const char *argv[])
{
  enum { SYMMETRIC, HERMITIAN, UNITARY, CSYMMETRIC, GENERAL};

  try {    
    cmatrix D;
    rmatrix RD;
    List<complex> ceig;
    List<double> reig;
    double tolerance=0.0;

    int count=1;
    const size_t times=getint(argc,argv,count,"Repeat times? ",1);
    const bool blockd=getlogical(argc,argv,count,"Block diagonal? ",false);
    if (blockd) {
      defzero=getfloat(argc,argv,count,"Zero size? ",0.0);
      tolerance=getfloat(argc,argv,count,"Detect tolerance? ",1e-5);
      cmatrix_eigensystem_controller.effectivezero=getfloat(argc,argv,count,"Effective zero? ",10*defzero);
    }
    cmatrix_eigensystem_controller.verbose=1;

    cout.precision(4);

    cout << "Test of matrix diagonalisation\n" << endl;
    
    //set_seed();
      
    for (int type=SYMMETRIC;type<=GENERAL;type++) {

      switch (type) {
      case SYMMETRIC:
	cout << "real symmetric"; break;
 //      case REAL:
//   	cout << "real general"; break;
      case HERMITIAN:
	cout << "hermitian"; break;
      case CSYMMETRIC:
	cout << "complex symmetric"; break;
      case UNITARY:
	cout << "unitary"; break;
      case GENERAL:
	cout << "general complex"; break;
      }
      cout << " matrices\n";

      const bool isreal((type==SYMMETRIC));

      for (int s=0;s<2;s++) {
	const int size=(s && !blockd) ? bigsize : 2;
	
	ceig.create(size);
	reig.create(size);

	cmatrix start(size,size);
	rmatrix rstart(size,size);
	
	for (int i=0;i<size;i++) {
	  switch (type) {
	  case SYMMETRIC:
	    rstart(i,i)=gauss();
	    for (int j=i+1;j<size;j++)
	      rstart(j,i)=rstart(i,j)=gauss();
	    break;
//   	  case REAL:
//   	    for (int j=0;j<size;j++)
//   	      rstart(i,j)=gauss();
//   	    break;
	  case HERMITIAN: case UNITARY:
	    start(i,i)=gauss();
	    for (int j=i+1;j<size;j++)
	      start(j,i)=conj(start(i,j)=complex(0.0,gauss()));//gauss()));
	    break;
	  case CSYMMETRIC:
	    start(i,i)=gauss();
	    for (int j=i+1;j<size;j++)
	      start(j,i)=start(i,j)=complex(gauss(),gauss());
	    break;
	  case GENERAL:
	    for (int j=0;j<size;j++)
	      start(i,j)=complex(gauss(),gauss());
	    break;
	  default:
	    throw Failed("Unknown matrix type");
	  }
	}

	if (type==UNITARY)
	  start=hermitian_exp(start,complex(0,1)); // exp(iA);

	if (blockd & s) {
	  if (isreal)
	    make_block(rstart,defzero);
	  else
	    make_block(start,defzero);
	}

	cout << "Created\n";

	cmatrix_eigensystem_controller.tolerance = (type==UNITARY) ? tolerance : 0.0;

	stopwatch.reset();

	switch (type) {
	case SYMMETRIC:
	  for (size_t loop=times;loop--;)
	    hermitian_eigensystem(RD,reig,rstart);
	  print_tests(RD,reig,rstart,times);
	  break;
	  //  	case REAL:
//  	  eigensystem(RD,reig,rstart);
//  	  print_tests_gen(RD,reig,rstart);
//  	  break;
	case HERMITIAN:
	  for (size_t loop=times;loop--;)
	    hermitian_eigensystem(D,reig,start);
	  print_tests(D,reig,start,times);
	  break;
	case UNITARY:
	  for (size_t loop=times;loop--;)
	    internal_eigensystem(D,ceig,start);
	  print_tests(D,ceig,start,times);
	  break;
	case GENERAL: case CSYMMETRIC:
	  for (size_t loop=times;loop--;)
	    internal_eigensystem(D,ceig,start);
	  print_tests_gen(D,ceig,start,times);
	  break;
	  //default:
	  //throw Failed("Unknown matrix type");
	}

#ifdef LCM_USE_EXTERNAL
	stopwatch.reset();
	cout << "LAPACK equivalent\n";

	switch (type) {
 	case SYMMETRIC: case HERMITIAN:
	  cout << "(Not implemented)\n\n";
	  break;
//  	  for (size_t loop=times;loop--;)
//  	    lapack_symmetric_eigensystem(RD,reig,rstart);
//  	  print_tests(RD,reig,rstart,times);
//  	  break;
// 	case REAL:
//  	  for (size_t loop=times;loop--;)
//  	    lapack_eigensystem(RD,reig,rstart);
//  	  print_tests(RD,reig,rstart,times);
//  	  break;
// 	case HERMITIAN:
// 	  for (size_t loop=times;loop--;)
// 	    lapack_hermitian_eigensystem(D,reig,start);
// 	  print_tests(D,reig,start,times);
// 	  break;
 	case UNITARY:
 	  for (size_t loop=times;loop--;)
 	    lapack_eigensystem(D,ceig,start);
 	  print_tests(D,ceig,start,times);
 	  break;
 	case GENERAL: case CSYMMETRIC:
 	  for (size_t loop=times;loop--;)
 	    lapack_eigensystem(D,ceig,start);
 	  print_tests_gen(D,ceig,start,times);
 	  break;
	}
#endif
      }
    }
  }
  catch (MatrixException &exc) {
    cerr << exc << endl;
  }

  return 0;
}

// void lapack_hermitian_eigensystem(rmatrix& V, List<double>& eigs, const rmatrix& a)
// {
//   if (!issquare(a)) throw NotSquare("lapack_hermitian_eigensystem");

//   int n=a.rows();
//   eigs.create(n);
//   int info=0;

//   V=a;

//   double worktmp=-1000.0;
//   int minusone=-1;

//   //get optimal workspace
//   dsyev_("V",CMATRIX_UPLO,&n,V.vector(),&n,eigs.vector(),&worktmp,&minusone,&info);

//   int lwork=int(0.5+worktmp);
//   ScratchList<double> work(lwork);

//   dsyev_("V",CMATRIX_UPLO,&n,V.vector(),&n,eigs.vector(),work.vector(),&lwork,&info);
  
//   if (info==0) return;
//   if (info>0) throw Failed("dsyev: failed to converge");
//   sprintf(errm,"dsyev: error in parameter %i",-info);
//   throw Failed(errm);
// }

// void lapack_hermitian_eigensystem(cmatrix& V, List<double>& eigs, const cmatrix& a)
// {
//   if (!issquare(a)) throw NotSquare("lapack_hermitian_eigensystem");

//   int n=a.rows();
//   eigs.create(n);
//   int info=0;

//   V=a;

//   int lwork=2*n;
//   ScratchList<complex> work(lwork);
//   ScratchList<double> rwork(3*n);

//   zheev_("V",CMATRIX_UPLO,&n,reinterpret_cast<complex*>(V.vector()),&n,eigs.vector(),work.vector(),&lwork,rwork.vector(),&info);
//   V.conj();
//   if (info==0) return;
//   if (info>0) throw Failed("zheev: failed to converge");
//   sprintf(errm,"zheev: error in parameter %i",-info);
//   throw Failed(errm);
// }
