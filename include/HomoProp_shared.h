#include "lcm_MultiHolder.h"
#include "BasePropagation.h"
#include "cmatrix.h"

namespace libcmatrix {

  template<class T> void transpose_emultiply(Matrix<T>& dest,const Matrix<T>& a,const Matrix<T>& b)
    {
      if (issame(&dest,&a)) throw ArgumentClash("transpose_emultiply");
      if (!b) throw Undefined("emultiply");
      const int rs=b.rows();
      const int cs=b.cols();
      dest.create(rs,cs);
      if ( (a.rows()!=cs) || (a.cols()!=rs) )
	throw Mismatch("emultiply");
      
      int point=0;
      for (int r=0;r<rs;r++) {
	for (int c=0;c<cs;c++) {
	  dest(point)=b(point)*a(c,r);
	  point++;
	}
      }
    }
  
  void conj_mla(cmatrix&, const cmatrix&, const cmatrix&);
  
  inline void unitary_isimtrans_identity(cmatrix& d,const cmatrix& VL,const cmatrix& VR)
    { conj_transpose_multiply(d,VL,VR); }
  
  template<class T> void squeeze(Matrix<T>& A)
    {
      A(0,0)=trace(A);
      const int dim=A.rows();
      for (int i=1;i<dim;i++) A(i,i)=0.0;
    }

  void unitary_isimtrans_identity(cmatrix&, const BaseList<complex>& VL, const BaseList<complex>& VR);

  inline void unitary_isimtrans_identity(cmatrix &d,const BaseList<complex> &VL,const cmatrix &VR)
    { conj_multiply(d,VL,VR); }

  inline complex unitary_isimtrans_identity(const complex& VL, const complex& VR) { return conj_multiply(VL,VR); }

  inline void unitary_isimtrans_identity(cmatrix &d,const cmatrix &VL,const BaseList<complex> &VR)
    { conj_transpose_multiply(d,VL,VR); }

  inline void unitary_isimtrans_identity(Matrix<double> &d,const Matrix<double> &VR,const Matrix<double> &VC)
    { transpose_multiply(d,VR,VC); }

//need to store in full matrix
template<typename T> void unitary_isimtransLR(cmatrix& d, const BaseList<complex>& VL, const BaseList<T>& a, const BaseList<complex>& VR)
{
  size_t n=a.length();
  if ((n!=VL.length()) || (n!=VR.length())) throw Mismatch("unitary_isimtrans");
  d.create(n,n);
  d=complex(0.0);
  for (;n--;) d(n,n)=a(n)*conj_multiply(VL(n),VR(n));
}

 void unitary_simtrans_gen_ip(cmatrix& ,const cmatrix& VR,const cmatrix& VC,cmatrix&);

template<typename T,typename F> void transform_mh(MultiHolder<complex>& table,const F& func,const MultiHolder<complex>& UTRs,const MultiHolder<complex>& UTCs,const Matrix<T>& a)
{ 
  const size_t dimR=UTRs.dimensions();
  const size_t dimC=UTCs.dimensions();

  if (dimR>1 && dimC>1) {
    table.dimensions(3);
    MultiMatrix<complex,3>& dest=table.multimatrix3();
    func(dest,UTRs.multimatrix3(),UTCs.multimatrix3(),a);
    return;
  }
  //vector
  table.dimensions(2);
  cmatrix& dest=table.matrix();

  if (dimR==3) 
    func(dest,UTRs.multimatrix3(),UTCs.list(),a.row());
  else {
    if (dimC==1)
      func(dest,UTRs.list(),UTCs.list(),a(0,0));
    else
      func(dest,UTRs.list(),UTCs.multimatrix3(),a.row());
  }
}

template<class F> void transform_mh(MultiHolder<complex>& htable,const F& func,const MultiHolder<complex>& UTRs,const MultiHolder<complex>& UTCs)
{
  if (UTRs.dimensions()!=UTCs.dimensions())
    throw Failed("do_transform_mh: rows and columns are incompatible");
  
  switch (UTRs.dimensions()) {
  case 3: 
    htable.dimensions(3);
    func(htable.multimatrix3(),UTRs.multimatrix3(),UTCs.multimatrix3());
    break;
  case 1:
    htable.dimensions(2);
    func(htable.matrix(),UTRs.list(),UTCs.list());
    break;
  default:
    throw Failed("do_transform_mh: unhandled combination");
  }
}

//Only diagonal blocks
template<typename F> void transform_mh(MultiHolder<complex>& table,const F& func,const MultiHolder<complex>& UTs,const BaseList<double>& a)
{
  if (UTs.dimensions()==3) {
    table.dimensions(3);
    func(table.multimatrix3(),UTs.multimatrix3(),UTs.multimatrix3(),a);
  }
  else {
    if (a.size()!=1) throw Mismatch("transform_mh");
    table.create(UTs.list().size());
    table=complex(a(0));
  }
  throw InternalError("transform_mh");
}

} //namespace libcmatrix
