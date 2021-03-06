#include "AsynchronousPropagation.h"
#include "HomoProp_shared.h"
#include "timer.h"

namespace libcmatrix {

// inline void doscale(cmatrix& sigma,complex scale)
// {
//   if (imag(scale))
//     sigma*=scale;
//   else {
//     double rscale=real(scale);
//     if (rscale!=1.0)
//       sigma*=rscale;
//   }
// }

  // static const char ASYNC_MISMATCH[]="AsynchronousFID: FID length vs propagator list";

  void AsynchronousFID::propagateRC(size_t n)
  {
    if (verbose_) {
      const double t1=starttime(n);
      std::cout << "Asynchronous propagation from ";
      prettyprint_time(t1) << " to ";
      prettyprint_time(t1+dt_) << '\n';
    }
    row()(tmp1,n);
    if (!sigma)
      tmp2.swap(tmp1);
    else
      multiply(tmp2,tmp1,sigma);    
    if (!isdiagonal())
      col()(tmp1,n);
    multiply_conj_transpose(sigma,tmp2,tmp1);
  }

void AsynchronousFID::add_FID(BaseList<complex> FIDv,complex scale,const cmatrix& sigma0,const cmatrix& det)
{
  if (!aretranspose(sigma0,det))
    throw Mismatch("add_FID: sigma0 and detect operator incompatible");

  //  const size_t n=FIDv.size();
//   if (1+nincs()!=n)
//     throw Mismatch(ASYNC_MISMATCH);
  
  multiply(sigma,scale,sigma0);
  complex val=trace_multiply(det,sigma);
  if (verbose_) {
    std::cout << "Density matrix at time step 0\n" << sigma << '\n';
    std::cout << "Adding " << val << " to FID\n";
  }
  FIDv.front()+=trace_multiply(det,sigma);
  for (size_t i=1;i<FIDv.size();i++) {
    propagateRC(i-1);
    val=trace_multiply(det,sigma);
    if (verbose_) {
      std::cout << "Density matrix at time step " << i << '\n' << sigma << '\n';
      std::cout << "Adding " << val << " to FID\n";
    }
    FIDv(i)+=val;
  }
}

void AsynchronousFID::add_FID(BaseList<complex> FIDv, complex scale)
{
  if (isdiagonal())
    throw Failed("add_FID: only valid for off-diagonal propagation");

  //const size_t n=FIDv.size();
 //  if (1+nincs()!=n)
//     throw Mismatch(ASYNC_MISMATCH);

  mla(FIDv.front(),scale,float_t(row().size()));
  sigma.clear();
  for (size_t i=1;i<FIDv.size();i++) {
    propagateRC(i-1);
    mla(FIDv(i),scale,trace(sigma));
  }
  
//   mla(FIDv.front(),scale,VR(0).rows());

//   for (size_t i=1;i<n;i++) {
//     multiply_conj_transpose(sigma,VR(i-1),VC(i-1));
//     mla(FIDv(i),scale,trace(sigma));
//   }
}

template<class T> void AsynchronousFID::add_FID_hermitian(BaseList<T> FIDv, T scale,const cmatrix& sigma0,const cmatrix& det)
{
  if (!isdiagonal())
    throw Failed("add_FID_hermitian: only valid for diagonal blocks");

  // const size_t n=FIDv.size();
//   if (1+nincs()!=n)
//     throw Mismatch(ASYNC_MISMATCH);

  // const StashPoint& Ugen(row());
  //  const MultiMatrix<complex,3>& U=row();

  multiply(sigma,scale,sigma0);
  FIDv.front()+=hermitian_trace_multiply(det,sigma);
  for (size_t i=1;i<FIDv.size();i++) {
    //    row()(tmp1,i-1);
    //multiply(tmp2,tmp1,sigma);    
    //multiply_conj_transpose(sigma,tmp2,tmp1);
    propagateRC(i-1);
    FIDv(i)+=hermitian_trace_multiply(det,sigma);
  }
}

template void AsynchronousFID::add_FID_hermitian(BaseList<double>, double, const cmatrix&, const cmatrix&);
template void AsynchronousFID::add_FID_hermitian(BaseList<complex>, complex, const cmatrix&, const cmatrix&);

}//namespace libcmatrix
