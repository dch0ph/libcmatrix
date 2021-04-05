/* Real matrices 
   Note that rmatrix is equivalent to Matrix<float_t> */

#ifndef _rmatrix_h_
#define _rmatrix_h_

#include "Matrix.h"
#include <cstdio>

namespace libcmatrix {

typedef Matrix<float_t> rmatrix;

rmatrix identity(int);

rmatrix inv(const rmatrix&);
void inv(rmatrix &,const rmatrix&);

rmatrix multiply_inv(const rmatrix&,const rmatrix&);
List<float_t> multiply_inv(const rmatrix&,const BaseList<float_t>&);
void multiply_inv_ip2(rmatrix&, rmatrix&); // NB *Both* arguments are modified
void multiply_inv_ip2(rmatrix&, BaseList<float_t>&);

void unitary_simtrans(rmatrix &,const rmatrix &,const rmatrix &,rmatrix* tmp =NULL);
void unitary_isimtrans(rmatrix &,const rmatrix &,const rmatrix &,rmatrix* tmp =NULL);
void unitary_simtrans(rmatrix &,const rmatrix &,const rmatrix &,const rmatrix &,rmatrix* tmp =NULL);
void unitary_isimtrans(rmatrix &,const rmatrix &,const rmatrix &,const rmatrix &,rmatrix* tmp =NULL);
void unitary_isimtrans(rmatrix&, const BaseList<float_t>&, const rmatrix&,rmatrix* tmp =NULL);
void unitary_simtrans(rmatrix&, const BaseList<float_t>&, const rmatrix&,rmatrix* tmp =NULL);

template<> inline void Matrix<double>::unitary_simtrans(const rmatrix& b,rmatrix* tmp) { ::libcmatrix::unitary_simtrans(*this,*this,b,tmp); }
template<> inline void Matrix<double>::unitary_isimtrans(const rmatrix& b,rmatrix* tmp) { ::libcmatrix::unitary_isimtrans(*this,*this,b,tmp); }


void hermitian_eigensystem(rmatrix&, BaseList<float_t>, const rmatrix&);
void hermitian_eigensystem(rmatrix&, List<float_t>&, const rmatrix&);
  
  rmatrix kronecker(const rmatrix &,const rmatrix &);
  void kronecker(rmatrix&, const rmatrix &,int);
  void kronecker(rmatrix&, int,const rmatrix &);
  void kronecker(rmatrix&, const rmatrix&, const rmatrix&);
  void kronecker_transpose(rmatrix&, int, const rmatrix&);
  void kronecker_transpose(rmatrix&, const rmatrix&, const rmatrix&);

bool hasoffdiagonal(const rmatrix&, double =0.0);

 double stddev(const BaseList<double>&); //!< doesn't really belong here!

void read_matrix(rmatrix&, const char*);
void read_matrix(rmatrix&, FILE*);
void read_vector(List<float_t>&, const char*);
void read_vector(List<float_t>&, FILE*);
void read_vector_ascii(List<float_t>&, FILE*);

void write_vector(FILE *,const BaseList<float_t> &,int =mxflag::doublep);
void write_vector(const char *,const BaseList<float_t> &,int =mxflag::doublep);
void write_matrix(FILE *,const rmatrix &,const char * =NULL,int =(mxflag::doublep | mxflag::block));
void write_matrix(const char *,const rmatrix &,const char * =NULL,int =(mxflag::doublep | mxflag::block));

//inline double norm(const rmatrix& a) { return norm(a.row()); }
 inline rmatrix enorm(const rmatrix& a) { rmatrix d(mxflag::temporary); enorm(d,a); return d; }

  //void spy(std::ostream&,const rmatrix &,double =1e-10);
  //inline void spy(const rmatrix &a,double tol =1e-10) { spy(std::cout,a,tol); }

void rmatrix_new_thread();

} // namespace libcmatrix
#endif

