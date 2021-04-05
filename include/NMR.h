#ifndef _NMR_h_
#define _NMR_h_

#include "rmatrix.h"
#include "space_T.h"
#include "List.h"
#include "cmatrix_utils.h"
#include "spin_system.h"
#include "lcm_FunctionObject.h"
#include "ListList.h"
#include "Warnings.h"

namespace libcmatrix {

  // simple object for specifying phases either as absolute number or as "x", "y", "-x" etc.
  class phase_spec {
  public:
    phase_spec() : value(0.0) {}
    phase_spec(double _value) : value(_value) {}
    phase_spec(char);
    phase_spec(const char *);
    double operator()() const { return value; }
    
  private:
    double value;
  };

const double mu0=M_PI*4e-7;
const double hbar=1.05457266e-34;

double dipolar_coupling(double gamma1,double gamma2,double r);
 inline double dipolar_coupling(nuclei_spec nuc1, nuclei_spec nuc2, double r) 
   { return dipolar_coupling(gamma(nuc1),gamma(nuc2),r); }

 //! return r corresponding to dipolar coupling (added 2017-11-07)
 double dipolar_coupling_to_r(double d, double gamma1,double gamma2);
 inline double dipolar_coupling_to_r(double d, nuclei_spec nuc1, nuclei_spec nuc2) 
   { return dipolar_coupling_to_r(d, gamma(nuc1),gamma(nuc2)); }


cmatrix spin_dipolar(const basespin_system &,size_t,size_t);
List<double> diag_spin_dipolar(const basespin_system&, size_t,size_t);
List<double> diag_spin_truncated_dipolar(const basespin_system&, size_t,size_t);

 cmatrix spin_J(const basespin_system &,size_t,size_t);
 List<double> diag_spin_weakJ(const basespin_system&, size_t,size_t);

cmatrix spin_CS(const basespin_system &,size_t);
inline List<double> diag_spin_CS(const basespin_system& sys,size_t n) { return diag_Iz(sys,n); }

  List<double> diag_spin_quadrupolar(const basespin_system &,size_t); //!< Quadrupolar Hamiltonian in terms of quadrupole *splitting*
  List<double> diag_spin_quadrupolar_Cq(const basespin_system &,size_t, size_t order =0); //!< Quadrupolar Hamiltonian in terms of Cq/chi
cmatrix spin_quadrupolar(const basespin_system &,size_t);

 enum ns_flag { NS_NONE, NS_FIRST, NS_BOTH };

 template<typename T,size_t N> class MultiMatrix;

  void spin_A2_ns(MultiMatrix<double,3>&, const basespin_system&, size_t,size_t,ns_flag);
inline void spin_quadrupolar_ns(MultiMatrix<double,3>& d, const basespin_system& sys, size_t i)
{ spin_A2_ns(d,sys,i,i,NS_BOTH); }
  
  extern Warning<InvalidParameter> NMR_asymmetry_warning;
  extern Warning<> propagation_emptyhamiltonian_warning;
  extern Warning<InvalidParameter> NMR_axis_ordering_warning;

  enum ordering_convention_t { convention_Haeberlen,
			       convention_NQR,
			       convention_IUPACIncreasing, 
			       convention_IUPACDecreasing };

  inline bool isorderingIUPAC(ordering_convention_t conv) {  return  (conv==convention_IUPACIncreasing) || (conv==convention_IUPACDecreasing); }

  void cartesian_to_eigensystem_symmetric(double& XX, double& YY, double& ZZ, Matrix<double>& V, const Matrix<double>& A, ordering_convention_t =convention_Haeberlen, const Euler_controller& =cmatrix_euler_controller);
  void cartesian_to_PAS_symmetric(double& iso11, double& aniso22, double& asym33, Euler& PAS, const Matrix<double>& A, ordering_convention_t =convention_Haeberlen, const Euler_controller& =cmatrix_euler_controller, bool doscale =false); //!< last parameter refers to scaling associated with odd libcmatrix spherical tensor convention (only used if reqd)
  inline void active_rotate_DCM_tensor(Matrix<double>& dest, const Matrix<double>& R, const Matrix<double>& A) { unitary_simtrans(dest,A,R); }
  
 space_T spatial_tensor(double iso11, double delta22, double eta33, ordering_convention_t =convention_Haeberlen);
  space_T spatial_tensor(double delta,double eta=0, ordering_convention_t =convention_Haeberlen);
//  space_T spatial_tensor_ns(double iso,double delta,double eta);
//  space_T spatial_tensor_ns(double delta,double eta=0);

  bool anisotropy_to_cartesian(double&, double&, double&, double,double,double, ordering_convention_t =convention_Haeberlen);  //!< returns \c false if asymmetry out of range
  inline void anisotropy_to_cartesian(double& xx,double& yy,double &zz, double aniso,double asym, ordering_convention_t convention =convention_Haeberlen) { anisotropy_to_cartesian(xx,yy,zz,0.0,aniso,asym,convention); }

  bool cartesian_to_anisotropy(double& iso, double& aniso, double& asym, double xx, double yy, double zz, ordering_convention_t =convention_Haeberlen);
  

 cmatrix spin_rf(const basespin_system&,nuclei_spec,phase_spec, double =1.0);

 int validateL(const space_T&,size_t =2);
  int getL(size_t);
 //void add_Hamiltonian_ns(cmatrix&,const space_T&,const BaseList<cmatrix>&);

 cmatrix Upulse(const basespin_system &,nuclei_spec,double,phase_spec);
 inline cmatrix Upulse(const basespin_system &sys,double angle,phase_spec phase) { return Upulse(sys,NULL_NUCLEUS,angle,phase); }
 void Upulse(cmatrix &,const rmatrix& Fx,const BaseList<double>& Fz,double,phase_spec);
 cmatrix Upulse(const rmatrix& Fx,const BaseList<double>& Fz,double angle,phase_spec);

/* cmatrix Usoftpulse(const basespin_system &,const cmatrix &H,nuclei_spec,double time,double freq,double angle,phase_spec); */
/* inline cmatrix Usoftpulse(const basespin_system &sys,const cmatrix &H,double time,double freq,double angle,phase_spec phase) { return Usoftpulse(sys,H,NULL_NUCLEUS,time,freq,angle,phase); } */
/* void Usoftpulse(cmatrix&, const cmatrix& Fx, const BaseList<double>& Fz,const cmatrix& H,double time,double freq,double angle,phase_spec); */

cmatrix rotatez(const cmatrix&, const BaseList<double>& Fz, double angle);
void rotatez_ip(cmatrix&, const BaseList<double>& Fz, double angle);
 inline void rotatez_ip(cmatrix& U, const BaseList<complex>& zrot) { U.unitary_simtrans(zrot); }

void rotatezfacs(BaseList<complex>, const BaseList<double>& Fz, double angle);
void rotatezfacs(cmatrix&, const BaseList<double>& Fz, const BaseList<double>& angles);
inline void rotatezfacs(List<complex>& zrot, const BaseList<double>& Fz,double angle)
  { zrot.create(Fz.length()); rotatezfacs( static_cast< BaseList<complex>& >(zrot), Fz,angle); }
  
const double MINUS_TWO_PI=-2*M_PI;
 const double TWO_PI=2*M_PI;

inline complex propagator(double a,double t) { return expi(MINUS_TWO_PI*t*a); }
void propagator(cmatrix&, const cmatrix&, double);
cmatrix propagator(const cmatrix&, double);

 void propagator_ns(cmatrix&, const cmatrix&, double,const BaseList<double>&);
 void propagator_ns(cmatrix&, const rmatrix&, double,const BaseList<double>&);

  void propagator_ns_Hrf(cmatrix&, const cmatrix&, const cmatrix&, double,const BaseList<double>&);
  void propagator_ns_Hrf(cmatrix&, const rmatrix&, const rmatrix&, double,const BaseList<double>&);
  void propagator_ns_Hrf(cmatrix&, const rmatrix&, const cmatrix&, double,const BaseList<double>&);

  void propagator_second(cmatrix&, const cmatrix&, double, const ListList<size_t>&, const BaseList<double>&);
  void propagator_second(cmatrix&, const rmatrix&, double, const ListList<size_t>&, const BaseList<double>&);

  void propagator_second_Hrf(cmatrix&, const cmatrix&, const cmatrix&, double, const ListList<size_t>&, const BaseList<double>&);
  void propagator_second_Hrf(cmatrix&, const rmatrix&, const rmatrix&, double, const ListList<size_t>&, const BaseList<double>&);
  void propagator_second_Hrf(cmatrix&, const rmatrix&, const cmatrix&, double, const ListList<size_t>&, const BaseList<double>&);
   
 void hermitian_eigensystem_ns(rmatrix& d,List<double>&, const rmatrix& a,const BaseList<double>& Hzeeman);
 void hermitian_eigensystem_ns(cmatrix& d,List<double>&, const cmatrix& a,const BaseList<double>& Hzeeman);
 void hermitian_eigenvalues_second(List<double>&, const cmatrix&, const BaseList<double>& Hzeeman, double =1e-8);
 void hermitian_eigenvalues_second(BaseList<double>, const cmatrix&, const BaseList<double>& Hzeeman, double =1e-8);
 void hermitian_eigenvalues_second(List<double>&, const rmatrix&, const BaseList<double>& Hzeeman, double =1e-8);
 void hermitian_eigenvalues_second(BaseList<double>, const rmatrix&, const BaseList<double>& Hzeeman, double =1e-8);
  void hermitian_eigensystem_blocked(cmatrix& vectors, List<double>& eigs, const cmatrix& H, const ListList<size_t>&); //!< \note \a vectors is empty if no diagonalisation was required
 void hermtian_eigensystem_blocked(rmatrix& vectors, List<double>& eigs, const rmatrix& H, const ListList<size_t>&);
 void hermitian_eigensystem_second(cmatrix&, List<double>&, const cmatrix&, const ListList<size_t>&, const BaseList<double>&);
 void hermitian_eigensystem_second(rmatrix&, List<double>&, const rmatrix&, const ListList<size_t>&, const BaseList<double>&);

void propagator(cmatrix&, const rmatrix&, double);
cmatrix propagator(const rmatrix&, double);

double diag_propagator(const complex& ,double);
void diag_propagator(BaseList<double>,const BaseList<complex>&, double);
void diag_propagator(List<double>&, const BaseList<complex>&, double);
void diag_propagator(cmatrix&, List<double>&, const cmatrix&, double, double tol =0.0);
void diag_propagator(cmatrix&, BaseList<double>, const cmatrix&, double, double tol =0.0);

void propagator(BaseList<complex>, const BaseList<double>&, double);
void propagator(List<complex>&, const BaseList<double>&, double);
//List<complex> propagator(const BaseList<double>&, double);

void evolution_matrix(cmatrix&, const BaseList<double>&, double);
void evolution_matrix(cmatrix&, const BaseList<complex>&);
void evolution_matrix(cmatrix&, const BaseList<complex>&, const BaseList<complex>&);
void evolution_matrix(BaseList<complex>&, const complex&, const BaseList<complex>&);
void evolution_matrix(BaseList<complex>&, const BaseList<complex> &,const complex&);

  struct PropGen_t : public BinaryFunction<cmatrix,double,double> 
  {
    virtual PropGen_t* clone() const =0;
    virtual size_t size() const =0;
  };
  struct DiagPropGen_t : public BinaryFunction<List<complex>,double,double>
  {
    virtual DiagPropGen_t* clone() const =0;
    virtual size_t size() const =0;
  };
  // typedef BinaryFunction<cmatrix,double,double> PropGen_t;
  // typedef BinaryFunction<List<complex>,double,double> DiagPropGen_t;

void propagators(BaseList<cmatrix>, const PropGen_t&, double t1,double t2);
void propagators(cmatrix&, size_t,const DiagPropGen_t&, double t1,double t2);

 template<class M> class StaticPropagator : public PropGen_t {
  const M Ham_;
  const List<double> Hzeeman_;
 public:
  StaticPropagator(const M& Hamv) : Ham_(Hamv) {}
    StaticPropagator(const M& Hamv, const BaseList<double>& Hzeemanv)
      : Ham_(Hamv), Hzeeman_(Hzeemanv) {}

   PropGen_t* clone() const { return new StaticPropagator(*this); }
   size_t size() const { return Ham_.rows(); }
   
  void operator() (cmatrix& U,double t1,double t2) const
  { 
    if (Hzeeman_.empty())
      propagator(U,Ham_,t2-t1);
    else
      propagator_ns(U,Ham_,t2-t1);
  }

  cmatrix operator() (double t1,double t2) const {
    cmatrix U(mxflag::temporary);
    (*this)(U,t1,t2);
    return U;
  }
};

  char operator_intersection(char,char); //return true if 'detect' and 'sigma0' operators intersect

  void coherencematrix_ip(Matrix<bool>&,const BaseList<double>&, const BaseList<int>&, int =0);
  inline void coherencematrix_ip(Matrix<bool>& mask, const BaseList<double>& zop, int coher, int maxcoher =0) { coherencematrix_ip(mask,zop,ExplicitList<1,int>(coher),maxcoher); }

 Matrix<bool> coherencematrix(const BaseList<double>&, const BaseList<int>&); 
 Matrix<bool> coherencematrix(const basespin_system&, const BaseList<int>&);
 inline Matrix<bool> coherencematrix(const basespin_system& sys, int which) { return coherencematrix(sys,BaseList<int>(1,&which)); }
  int coherencematrix_range(const BaseList<double>&, const BaseList<int>& cohers);

 Matrix<bool> coherencematrix(const basespin_system&, const BaseList<size_t>&, const ListList<int>&);
 Matrix<bool> coherencematrix(const basespin_system&, const BaseList<size_t>&, const BaseList<int>&);
 Matrix<bool> coherencematrix(const basespin_system&, nuclei_spec, const BaseList<int>&);
 inline Matrix<bool> coherencematrix(const basespin_system& sys, nuclei_spec spec, int which) 
   { return coherencematrix(sys,spec,BaseList<int>(1,&which)); }

} //namespace libcmatrix
#endif

