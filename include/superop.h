// Random functions for super operators

#include "cmatrix.h"
#include "rmatrix.h"
#include "PartitionedMatrix.h"

namespace libcmatrix {

  void superop_multiply(cmatrix &,const cmatrix &,const cmatrix &);
  void superop_multiply(rmatrix &,const rmatrix &,const rmatrix &);
  cmatrix superop_multiply(const cmatrix &,const cmatrix &);
  rmatrix superop_multiply(const rmatrix &,const rmatrix &);
  
  void commutator(cmatrix&, const cmatrix&);
  cmatrix commutator(const cmatrix &);
  cmatrix double_commutator(const cmatrix &);
  void commutator(rmatrix&, const rmatrix &);
  rmatrix commutator(const rmatrix &);
  rmatrix double_commutator(const rmatrix &);

  void commutator(cmatrix&, const cmatrix&, const BaseList<size_t>&, const BaseList<size_t>&);
  cmatrix commutator(const cmatrix&, const BaseList<size_t>&, const BaseList<size_t>&);
  void commutator(rmatrix&, const rmatrix&, const BaseList<size_t>&, const BaseList<size_t>&);
  rmatrix commutator(const rmatrix&, const BaseList<size_t>&, const BaseList<size_t>&);

  List<double> commutator(const BaseList<double>&);
  List<complex> commutator(const BaseList<complex>&);
  void commutator(BaseList<complex> d, const BaseList<complex>& a);
  void commutator(BaseList<double> d, const BaseList<double>& a);
    
  extern bool lcm_pade_use; //!< use Pade approximants for Liouville space propagation

  cmatrix propagatorL(const cmatrix&, double, const matrix_partition* =NULL); //!< Liouville space propagator
  void propagatorL(cmatrix&, const cmatrix&, double, const matrix_partition* =NULL); //!< Liouville space propagator
  void pade_propagatorL(cmatrix&, const cmatrix&, double, const matrix_partition* =NULL);   //! use Pade approximants to determine Liouville propagator

  extern Warning<> partitioning_ignored_warning;
}
