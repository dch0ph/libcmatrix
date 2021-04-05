// routines for non-secular / second-order Hamiltonians
#undef LCM_SUPPRESS_VIEWS
#include "NMR.h"
#include "MAS.h"
#include "MultiMatrix.h"
#include "tensorop.h"

namespace libcmatrix {

/*WARNING: Tricky to check if Hzeeman 'matches' supplied Hamiltonian (a)
  If they don't, the subtraction of Larmor frequencies will return rubbish */

void subtract_zeeman(BaseList<double> eigs,const BaseList<double> Hzeeman)
{
  const size_t n=Hzeeman.size();
  if (eigs.size()!=n)
    throw Mismatch("subtract_zeeman");

  //Subtract off Zeeman eigenstates - don't know ordering!
#ifndef NDEBUG
  std::cout << "Original eigenvalues: " << eigs << '\n';
  std::cout << "Zeeman eigenvalues: " << Hzeeman << '\n';
#endif
  List<size_t> which(n);
  which=slice(0,n);
  size_t curindex=0;
  double best_diff;
  double last_diff=-1.0;
  bool warn=false;
  while (which.size()>1) {
    const double val=eigs(curindex);
    size_t best_index=which.size()-1;
    best_diff=fabs(Hzeeman(which(best_index))-val);
    for (size_t index=best_index;index--;) {
      const double diff=fabs(Hzeeman(which(index))-val);
      if (diff<best_diff) {
	last_diff=best_diff;
	best_diff=diff;
	best_index=index;
      }
    }
    if ((last_diff>0.0) && (best_diff*2>last_diff))
      warn=true;
    eigs(curindex++)-=Hzeeman(which(best_index));
    which(best_index)=which.back(); 
    which.pop_back();
  }
  eigs(curindex)-=Hzeeman(which.front());
#ifndef NDEBUG
  std::cout << "RF eigenvalues: " << eigs << '\n';
#endif
  if (warn)
    std::cerr << "Warning propagator_ns: assignment of eigenvalues may be ambiguous\n";
}

void hermitian_eigensystem_ns(cmatrix& vectors,List<double>& eigs, const cmatrix& a,const BaseList<double>& Hzeeman)
{
  if (!issquare(a))
    throw NotSquare("hermitian_eigensystem_ns");
  if (a.rows()!=Hzeeman.size())
    throw Mismatch("hermitian_eigensystem_ns");
  cmatrix tmp(a);
  tmp+=Hzeeman;
  //Diagonalise Hamiltonian in non-rotating frame
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "H in lab frame:\n" << tmp;
  hermitian_eigensystem(vectors,eigs,tmp);
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "V:\n" << vectors;
  subtract_zeeman(eigs,Hzeeman);
}

void hermitian_eigensystem_ns(rmatrix& vectors,List<double>& eigs, const rmatrix& a,const BaseList<double>& Hzeeman)
{
  if (!issquare(a))
    throw NotSquare("hermitian_eigensystem_ns");
  if (a.rows()!=Hzeeman.size())
    throw Mismatch("hermitian_eigensystem_ns");
  rmatrix tmp(a);
  tmp+=Hzeeman;
  //Diagonalise Hamiltonian in non-rotating frame
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "H in lab frame:\n" << tmp;
  hermitian_eigensystem(vectors,eigs,tmp);
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "V:\n" << vectors;
  subtract_zeeman(eigs,Hzeeman);
}

template<typename T> void propagator_ns_(cmatrix& d,const Matrix<T>& a,double t,const BaseList<double>& Hzeeman)
{
  Matrix<T> vectors;
  List<double> eigs;
  hermitian_eigensystem_ns(vectors,eigs,a,Hzeeman);
  const size_t n=Hzeeman.size();
  //Diagonalise Hamiltonian in non-rotating frame
  ScratchList<complex> ceigs(n);
  const double scale=-2*M_PI*t;
  for (size_t i=n;i--;)
    ceigs(i)=expi(scale*eigs(i));
  unitary_simtrans(d,ceigs,vectors);
}

void propagator_ns(cmatrix& d,const rmatrix& a,double t,const BaseList<double>& Hzeeman)
{
  propagator_ns_(d,a,t,Hzeeman);
}

void propagator_ns(cmatrix& d,const cmatrix& a,double t,const BaseList<double>& Hzeeman)
{
  propagator_ns_(d,a,t,Hzeeman);
}

template<typename T1,typename T2> void propagator_ns_Hrf_(cmatrix& U, const Matrix<T1>& Hsys, const Matrix<T2>& Hrf, double t,const BaseList<double>& Hzeeman)
{
  Matrix<T1> vectors;
  Matrix<T2> HrfT;
  Matrix<complex> Utmp;
  List<double> eigs;
  hermitian_eigensystem_ns(vectors,eigs,Hsys,Hzeeman);
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "Result of hermitian_eigensystem_ns:\nEigenvectors:\n" << vectors << "Eigenvalues: " << eigs << '\n';
  unitary_isimtrans(HrfT,Hrf,vectors);
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "Hrf in eigenbasis\n" << HrfT << '\n';
  HrfT+=eigs;
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "H_total\n" << HrfT << '\n';
  propagator(Utmp,HrfT,t);
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "H_total (in eigenbasis):\n" << HrfT << "U (in eigenbasis):\n" << Utmp;
  unitary_simtrans(U,Utmp,vectors);
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "U from propagator_Hrf_ns:\n" << U << '\n';
}

void propagator_ns_Hrf(cmatrix& d,const rmatrix& Hsys, const rmatrix& Hrf, double t,const BaseList<double>& Hzeeman)
{
  propagator_ns_Hrf_(d,Hsys,Hrf,t,Hzeeman);
}

void propagator_ns_Hrf(cmatrix& d,const cmatrix& Hsys, const cmatrix& Hrf, double t,const BaseList<double>& Hzeeman)
{
  propagator_ns_Hrf_(d,Hsys,Hrf,t,Hzeeman);
}

void propagator_ns_Hrf(cmatrix& d,const rmatrix& Hsys, const cmatrix& Hrf, double t,const BaseList<double>& Hzeeman)
{
  propagator_ns_Hrf_(d,Hsys,Hrf,t,Hzeeman);
}

  namespace {
    inline double lreal(double x) { return x; }
    inline double lreal(const complex& z) { return z.real(); }
  }

// NB returns empty matrix if no diagonalisation was required
template<class T> void hermitian_eigensystem_blocked_(Matrix<T>& V, List<double>& eigs, const Matrix<T>& H, const ListList<size_t>& blkstr)
{
  if (!issquare(H))
    throw NotSquare("hermitian_eigensystem");
  const int n=H.rows();
  if ((n!=blkstr.row().size()))
    throw Mismatch("hermitian_eigensystem");

  eigs.create(n);
  //  V.create(n,n,T(0.0));
  V.clear();

  size_t ptr=0;
  Matrix<T> tmp;
  Matrix<T> vecs;
  for (size_t blk=0;blk<blkstr.size();blk++) {
    const BaseList<size_t> curblk(blkstr(blk));
    const size_t curn(curblk.size());
    if (curn==1) {
      const size_t ind(curblk.front());
      eigs(ptr)=lreal(H(ind,ind));
      if (!!V)
	V(ptr,curblk.front())=1.0;
    }
    else {
      tmp=H(curblk,curblk);
      if (hasoffdiagonal(tmp)) {
	BaseList<double> ceigs(curn,eigs.vector()+ptr);
	hermitian_eigensystem(vecs,ceigs,tmp);
	if (!V)
	  V.create(n,n,T(0.0));
	V(range(ptr,ptr+curn-1),curblk)=vecs;
      }
      else {
	for (size_t i=tmp.rows();i--;) {
	  eigs(ptr+i)=lreal(tmp(i,i)); // similarly don't need to create V explicitly
	  if (!!V)
	    V(ptr+i,curblk(i))=1.0;
	}
      }
    }
    ptr+=curn;
  }  
}

void hermitian_eigensystem_blocked(cmatrix& vectors, List<double>& eigs, const cmatrix& H, const ListList<size_t>& blkstr)
{
  hermitian_eigensystem_blocked_(vectors,eigs,H,blkstr);
}

void hermitian_eigensystem_blocked(rmatrix& vectors, List<double>& eigs, const rmatrix& H, const ListList<size_t>& blkstr)
{
  hermitian_eigensystem_blocked_(vectors,eigs,H,blkstr);
}

void hermitian_eigenvalues_second(List<double>& eigs, const cmatrix& H, const BaseList<double>& Hzeeman, double degfac)
{
  eigs.create(H.rows());
  hermitian_eigenvalues_second(static_cast< BaseList<double> >(eigs),H,Hzeeman,degfac);
}

void hermitian_eigenvalues_second(List<double>& eigs, const rmatrix& H, const BaseList<double>& Hzeeman, double degfac)
{
  eigs.create(H.rows());
  hermitian_eigenvalues_second(static_cast< BaseList<double> >(eigs),H,Hzeeman,degfac);
}

void hermitian_eigenvalues_second(BaseList<double> eigs, const cmatrix& H, const BaseList<double>& Hzeeman, double degfac)
{
  if (degfac<0.0)
    throw InvalidParameter("hermitian_eigenvalues_second: degeneracy criterion cannot be <=0");
  const double degfac2=degfac*degfac;
  const bool ignore=(degfac==0.0);

  if (!issquare(H))
    throw NotSquare("hermitian_eigenvalues_second");
  const int n=H.rows();
  if (eigs.size()!=n)
    throw Mismatch("hermitian_eigenvalues_second");
  for (size_t i=0;i<n;i++) {
    double val(real(H(i,i)));
    for (size_t j=0;j<n;j++) {
      if (i==j)
	continue;
      const double diff=Hzeeman(i)-Hzeeman(j);
      const double cor=norm(H(i,j));
      if (cor>degfac2) {
	if (fabs(diff)<=degfac) {
	  if (!ignore)
	    throw Failed("hermitian_eigenvalues_second: degenerate eigenvalues");
	}
	else
	  val+=cor/diff;
      }
    }
    eigs(i)=val;
  }
}

//inline double sqr(double x) { return x*x; }

void hermitian_eigenvalues_second(BaseList<double> eigs, const rmatrix& H, const BaseList<double>& Hzeeman, const double degfac)
{
  if (degfac<0.0)
    throw InvalidParameter("hermitian_eigenvalues_second: degeneracy criterion cannot be <=0");
  const double degfac2=degfac*degfac;
  const bool ignore=(degfac==0.0);

  if (!issquare(H))
    throw NotSquare("hermitian_eigenvalues_second");
  const int n=H.rows();
  if (eigs.size()!=n)
    throw Mismatch("hermitian_eigenvalues_second");
  for (size_t i=0;i<n;i++) {
    double val(H(i,i));
    for (size_t j=0;j<n;j++) {
      if (i==j)
	continue;
      const double diff=Hzeeman(i)-Hzeeman(j);
      const double cor=norm(H(i,j));
      if (cor>degfac2) {
	if (fabs(diff)<=degfac) {
	  if (!ignore)
	    throw Failed("hermitian_eigenvalues_second: degenerate eigenvalues");
	}
	else
	  val+=cor/diff;
      }
    }
    eigs(i)=val;
  }
}

template<class M> void hermitian_eigensystem_second_(M& V, List<double>& seigs, const M& H, const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
#ifndef NDEBUG
  std::cout << "Hsys\n" << H;
  std::cout << "Block structure: " << blkstr << '\n';
#endif
  hermitian_eigensystem_blocked(V,seigs,H,blkstr);
  M Ht;
  const M& Huse= !V ? H : Ht;
  if (!!V)
    unitary_isimtrans(Ht,H,V);
#ifndef NDEBUG
  std::cout << "Hsys(0) eigenvalues: " << seigs << '\n';
  std::cout << "Hsys(1) in eigenbasis of *Zeeman* Hamiltonian\n" << Huse;
#endif
  hermitian_eigenvalues_second(seigs,Huse,Hzeemant,0.0); //ignore degeneracy
}

void hermitian_eigensystem_second(cmatrix& V, List<double>& seigs, const cmatrix& H, const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  hermitian_eigensystem_second_(V,seigs,H,blkstr,Hzeemant);
}

void hermitian_eigensystem_second(rmatrix& V, List<double>& seigs, const rmatrix& H, const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  hermitian_eigensystem_second_(V,seigs,H,blkstr,Hzeemant);
}

template<> cmatrix* HomogeneousPropagator_second<cmatrix>::scratch_space() const { return &H; }

#include "lcm_accumulate.hpp"

template<class M> void HomogeneousPropagator_second<M>::operator()(cmatrix& U,double t1,double t2) const
{
  if (t1>t2)
    throw InvalidParameter("HomogeneousPropagator_second: disordered time values");
  bool haveflag=false;
  cmatrix Utmptmp;

  IntervalSamplerBase* usesamplerp=&lcm_defsampler;
  if (Hamp->samplerp() ) {
    if (!Samplerp)
      Samplerp.reset( (Hamp->samplerp() )->clone());       
    usesamplerp=Samplerp.get();
  }
  
  const double maxdttol=intdt+lcm_MAS_timingtol;

  for (;t1<t2;) {
    double nextt1=t2;
    if (nextt1-t1>maxdttol)
      nextt1=t1+intdt;

    const double midt=(*usesamplerp)(t1,nextt1);

    (*Hamp)(H,midt);
    hermitian_eigensystem_second(V,eigs,H,blkstr,Hzeemant);
    propagator(ceigs,eigs,nextt1-t1);
    if (cmatrix_eigensystem_controller.verbose>1)
      std::cout << "Propagator (in eigenbasis) at 'midpoint' of " << (1e6*t1) << " and " << (1e6*nextt1) << ", t=" << (midt*1e6) << " us: " << ceigs << '\n';
  if (!V)
      accumulate_(U,ceigs,Utmptmp,haveflag);
    else {
      if (cmatrix_eigensystem_controller.verbose>1)
	std::cout << "Eigenbasis\n" << V;
      isunitary(V,"V");	
      unitary_simtrans(Utmp,ceigs,V,scratch_space());
      accumulate_(U,Utmp,Utmptmp,haveflag);      
    }
    if (cmatrix_eigensystem_controller.verbose>1)
      std::cout << "Accumulated U from timestep at 'midpoint' of " << (1e6*t1) << " and " << (1e6*nextt1) << ", t=" << (midt*1e6) << " us\n" << U;
    isunitary(U);
    t1=nextt1;
  }
  assert(haveflag);
}

template<class T> void propagator_second_(cmatrix& U, const Matrix<T>& H,double dt,const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  Matrix<T> V;
  const size_t n(Hzeemant.size());
  List<double> eigs;
  hermitian_eigensystem_second(V,eigs,H,blkstr,Hzeemant);
  ScratchList<complex> Udiag(n);
  propagator(Udiag,eigs,dt);  
  if (!V)
    full(U,Udiag);
  else
    unitary_simtrans(U,Udiag,V);
}

template<class T1,class T2> void propagator_second_Hrf_(cmatrix& U, const Matrix<T1>& Hsys, const Matrix<T2>& Hrf, double dt,const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  Matrix<T1> V;
  Matrix<T2> HrfT;
  Matrix<complex> Utmp;
  List<double> eigs;
  hermitian_eigensystem_second(V,eigs,Hsys,blkstr,Hzeemant);
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "From hermitian_eigensystem_second:\nEigenvectors\n" << V << "Eigenvalues: " << eigs << '\n';
  if (!V)
    HrfT=Hrf;
  else
    unitary_isimtrans(HrfT,Hrf,V);
  HrfT+=eigs;
  cmatrix& Uuse= !V ? U : Utmp;
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "Hrf + Hsys in eigenbasis\n" << HrfT << '\n';
  propagator(Uuse,HrfT,dt); // calc propagator in Zeeman eigenbasis
  isunitary(Uuse);
  if (!!V)
    unitary_simtrans(U,Uuse,V); // transform to starting frame
  if (cmatrix_eigensystem_controller.verbose>1)
    std::cout << "Output from propagator_second_rf\n" << U << '\n';
}

void propagator_second(cmatrix& U, const rmatrix& H,double dt,const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  propagator_second_(U,H,dt,blkstr,Hzeemant);
}

void propagator_second(cmatrix& U, const cmatrix& H,double dt,const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  propagator_second_(U,H,dt,blkstr,Hzeemant);
}

void propagator_second_Hrf(cmatrix& U, const rmatrix& Hsys, const rmatrix& Hrf, double dt,const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  propagator_second_Hrf_(U,Hsys,Hrf,dt,blkstr,Hzeemant);
}

void propagator_second_Hrf(cmatrix& U, const cmatrix& Hsys, const cmatrix& Hrf, double dt,const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  propagator_second_Hrf_(U,Hsys,Hrf,dt,blkstr,Hzeemant);
}

void propagator_second_Hrf(cmatrix& U, const rmatrix& Hsys, const cmatrix& Hrf, double dt,const ListList<size_t>& blkstr, const BaseList<double>& Hzeemant)
{
  propagator_second_Hrf_(U,Hsys,Hrf,dt,blkstr,Hzeemant);
}

template class HomogeneousPropagator_second<rmatrix>;
template class HomogeneousPropagator_second<cmatrix>;

// static const double two3rd=2.0/3.0;
// static const double sqrt3rd=std::sqrt(1.0/3);

// space_T spatial_tensor_ns(double iso,double d,double asym)
// {
//   double xx,yy,zz;
//   anisotropy_to_cartesian(xx,yy,zz,-sqrt3rd*iso,d*two3rd,asym);
//   return A2(xx,yy,zz);
// }

// space_T spatial_tensor_ns(double d,double asym)
// {
//   double xx,yy,zz;
//   anisotropy_to_cartesian(xx,yy,zz,0.0,d*two3rd,asym);
//   return A2(xx,yy,zz,2);
// }

//static const double sqrttwo3rd=::std::sqrt(2.0/3.0);
static const double sqrtsix=::std::sqrt(6.0);

void spin_A2_ns(MultiMatrix<double,3>& d, const basespin_system& sys,size_t i, size_t j,ns_flag nstype)
{
  static const double T1factor=::std::sqrt(3.0/2.0);
  switch (nstype) {
  case NS_FIRST: {
    d.create(3,sys.rows(),sys.cols());
    // these are scaled by sqrt(3/2) from the normal tensor ops
    rmatrix d0(d(0U),mxflag::nondynamic);
    real(d0,I(sys,i,'-',j,'z')); d0*=T1factor; //d0*=0.5;
    rmatrix d1(d(1U),mxflag::nondynamic);
    //! Isn't this assuming truncation across heteronuclei ?
    real(d1,I(sys,i,'z',j,'z')); d1*=2.0; //d1*=sqrttwo3rd;
    rmatrix d2(d(2U),mxflag::nondynamic);
    real(d2,I(sys,i,'+',j,'z')); d2*=-T1factor; //d2*=-0.5;
  }
    //    d*=sqrttwo3rd;
    break;
  case NS_BOTH:
    d.create(5,sys.rows(),sys.cols());
    for (int m=-2;m<=2;m++) {
      rmatrix dm(d(size_t(m+2)));
      real(dm,T2(sys,i,j,2,m));
    }
    d*=sqrtsix; //!< 'compensate' for 1/sqrt(6) scaling of A2
    break;
  default:
    throw InvalidParameter("spin_A2_ns: nonsecular flag");
  }

}

} //namespace libcmatrix
