#include "Floquet.h"

namespace libcmatrix {

  namespace {
    template<class T> double real_(const T& v) { return real(v); }
    template<> double real_(const double& v) { return v; }
    template<class T> double norm_(const T& v) { return norm(v); }
    template<> double norm_(const double& v) { return v*v; }
  }

  template<typename T> Matrix<T> floquet_operator_(const Matrix<T>& a,int N) 
  {
    if (N<0)
      throw InvalidParameter("floquet_operator: order cannot be negative");
    if (!issquare(a))
      throw NotSquare("floquet_operator");
    if (N==0)
      return a;
    const int n=a.rows();
    const int orders=2*N+1;

    Matrix<T> d(n*orders,n*orders,T(0.0),mxflag::temporary);
    range sel(0,n-1);
    for (int i=orders;i--;) { //reproduce input matrix along diagonal
      d(sel,sel)=a;
      sel+=n;
    }
    return d;
  }

  cmatrix floquet_operator(const cmatrix& a,int N) { return floquet_operator_(a,N); }
  rmatrix floquet_operator(const rmatrix& a,int N) { return floquet_operator_(a,N); }

  static int process_order(int N,int length) {
    if (N<0 || (length==0))
      throw InvalidParameter("floquet_hamiltonian: order cannot be negative");
    if ((length-1) % 2)
      throw InvalidParameter("floquet_hamiltonian: input vector must have length 2n+1");

    const int inorder=(length-1)/2;
    if (N<inorder)
      throw InvalidParameter("floquet_hamiltonian: Floquet order less than Hamiltonian order!");
    return inorder;
  }

//   template<typename T> Matrix<T> floquet_hamiltonian_(const BaseList<Matrix<T> >& as,int N,double v) 
//   {
//     const int inorder=process_order(N,as.length());
//     if (N==0) return as(inorder);

//     const int n=as(0).rows(); //non-square or mismatched matrices will result in Mismatch exceptions later
//     const int orders=2*N+1;

//     Matrix<T> d(n*orders,n*orders,T(0.0),mxflag::temporary);
//     for (int rank=0;rank<=inorder;rank++) {
//       for (int r=rank;r<orders;r++) {
// 	const range selr(r*n,(r+1)*n-1);
// 	if (rank==0) {
// 	  d(selr,selr)=as(inorder);
// 	  const double freq=v*(r-inorder);
// 	  for (int j=n;j--;) {
// 	    const size_t ind=selr(j);
// 	    d(ind,ind)+=freq; //add multiple of cycle frequency
// 	  }
// 	}
// 	else {
// 	  const int c=r-rank;
// 	  const range selc(c*n,(c+1)*n-1);
// 	  d(selr,selc)=as(inorder-rank);
// 	  d(selc,selr)=as(inorder+rank);
// 	}
//       }
//     }
//     return d;
//   }

//   cmatrix floquet_hamiltonian(const BaseList<cmatrix>& a,int N,double v) { return floquet_hamiltonian_(a,N,v); }
//   rmatrix floquet_hamiltonian(const BaseList<rmatrix>& a,int N,double v) { return floquet_hamiltonian_(a,N,v); }

  cmatrix floquet_hamiltonian(const SpinningHamiltonian& H,int N,double v)
  {
    const int inorder=H.rank();
    if (N<inorder) 
      throw InvalidParameter("floquet_hamiltonian: Floquet order less than Hamiltonian order!");

    const cmatrix& H0=H.component(0);
    const int n=H0.rows(); //non-square or mismatched matrices will result in Mismatch exceptions later
    const int orders=2*N+1;

    cmatrix d(n*orders,n*orders,complex(0.0),mxflag::temporary);
    for (int rank=0;rank<=inorder;rank++) {
      for (int r=rank;r<orders;r++) {
	const range selr(r*n,(r+1)*n-1);
	if (rank==0) {
	  d(selr,selr)=H0;
	  const double freq=v*(r-N);
	  for (size_t j=n;j--;) {
	    const size_t ind=selr(j);
	    d(ind,ind)+=freq; //add multiple of cycle frequency
	  }
	}
	else {
	  const int c=r-rank;
	  const range selc(c*n,(c+1)*n-1);
	  d(selr,selc)=H.component(-rank);
	  d(selc,selr)=H.component(rank);
	}
      }
    }
    return d;
  }

  template<class T> Matrix<T> floquet_hamiltonian_(const BaseList<T>& Hels,int N,double v)
  {
    const int inorder=(Hels.length()-1)/2;
    if (inorder<0 || (2*inorder+1!=Hels.length()))
      throw InvalidParameter("floquet_hamiltonian: invalid Hamiltonian specification");
    if (N<inorder)
      throw InvalidParameter("floquet_hamiltonian: Floquet order less than Hamiltonian order!");
    
    const size_t orders=2*N+1;
    
    Matrix<T> d(orders,orders,T(0.0),mxflag::temporary);

    const T H0=Hels(inorder);
    for (int m=-N;m<=N;m++) d(m+N,m+N)=H0+m*v;

    for (int rank=1;rank<=inorder;rank++) {
      const T Hvalp=Hels(rank+inorder);
      const T Hvalm=Hels(-rank+inorder);
      for (size_t r=rank;r<orders;r++) {
	const size_t c=r-rank;
	d(r,c)=Hvalm;
	d(c,r)=Hvalp;
      }
    }
    return d;
  }

  cmatrix floquet_hamiltonian(const BaseList<complex>& Hels,int N,double v)
  {
    return floquet_hamiltonian_(Hels,N,v);
  }

  rmatrix floquet_hamiltonian(const BaseList<double>& Hels,int N,double v)
  {
    return floquet_hamiltonian_(Hels,N,v);
  }

  cmatrix floquet_hamiltonian(const RealSpinningHamiltonian& H,int N,double v)
  {
    const int inorder=H.rank();
    if (N<inorder) 
      throw InvalidParameter("floquet_hamiltonian: Floquet order less than Hamiltonian order!");

    const rmatrix& H0=H.component0();
    const int n=H0.rows(); //non-square or mismatched matrices will result in Mismatch exceptions later
    const int orders=2*N+1;

    cmatrix d(n*orders,n*orders,complex(0.0),mxflag::temporary);
    for (int rank=0;rank<=inorder;rank++) {
      for (int r=rank;r<orders;r++) {
	const range selr(r*n,(r+1)*n-1);
	if (rank==0) {
	  d(selr,selr)=H0;
	  const double freq=v*(r-N);
	  for (size_t j=n;j--;) {
	    const size_t ind=selr(j);
	    d(ind,ind)+=freq; //add multiple of cycle frequency
	  }
	}
	else {
	  const int c=r-rank;
	  const range selc(c*n,(c+1)*n-1);
	  const cmatrix Hc=H.component(rank);
	  conj_ip(d(selr,selc)=Hc);
	  d(selc,selr)=Hc;
	}
      }
    }
    return d;
  }

  template<typename T> Matrix<T> floquet_hamiltonian_(const BaseList<Matrix<T> >& as1,int N1,double v1,const BaseList<Matrix<T> >& as2,int N2,double v2) 
  {
    const int inorder1=process_order(N1,as1.length());
    const int inorder2=process_order(N2,as2.length());

    Matrix<T> H0(as1(inorder1));
    H0+=as2(inorder2);

    const int n=H0.rows();
    //non-square or mismatched matrices will result in Mismatch exceptions later

    const int orders1=2*N1+1;
    const int orders2=2*N2+1;

    const int totsize=n*n*orders1*orders2;

    Matrix<T> d(totsize,totsize,T(0.0),mxflag::temporary);

    for (int rank1=0;rank1<=inorder1;rank1++) {
      for (int r1=rank1;r1<orders1;r1++) {
	const int c1=r1-rank1;
	for (int r2=0;r2<orders2;r2++) {
	  const int totr=r1*orders2+r2;
	  const range selr(totr*n,(totr+1)*n-1);

	  if (rank1==0) {
	    d(selr,selr)=H0;
	    const double freq=(r1-N1)*v1;
	    for (size_t j=n;j--;) {
	      const size_t ind=selr(j);
	      d(ind,ind)+=freq; //add multiple of cycle frequency
	    }
	  }
	  else {
	    const int totc=c1*orders2+r2;
	    const range selc(totc*n,(totc+1)*n-1);

	    d(selr,selc)=as1(inorder1-rank1);
	    d(selc,selr)=as1(inorder1+rank1);
	  }
	}
      }
    }
    for (int r1=0;r1<orders1;r1++) {
      for (int rank2=0;rank2<=inorder2;rank2++) {
	for (int r2=rank2;r2<orders2;r2++) {
	    
	  const int totr=r1*orders2+r2;
	  const range selr(totr*n,(totr+1)*n-1);

	  if (rank2==0) {
	    const double freq=(r2-N2)*v2;
	    for (size_t j=n;j--;) {
	      const size_t ind=selr(j);
	      d(ind,ind)+=freq; //add multiple of cycle frequency
	    }
	  }
	  else {
	    const int c2=r2-rank2;
	    const int totc=r1*orders2+c2;
	    const range selc(totc*n,(totc+1)*n-1);

	    d(selr,selc)=as2(inorder2-rank2);
	    d(selc,selr)=as2(inorder2+rank2);
	  }
	}
      }
    }
    return d;
  }

  int floquet_check(int Fdim,int nh) {
    const int forderm1=Fdim/nh-1;
    if ( (Fdim % nh) || (forderm1 % 2))
      throw Mismatch("Floquet matrix is not 2N+1 * Hilbert space size");
    return forderm1/2;
  }

 cmatrix floquet_hamiltonian(const BaseList<cmatrix>& as1,int N1,double v1,const BaseList<cmatrix>& as2,int N2,double v2) { return floquet_hamiltonian_(as1,N1,v1,as2,N2,v2); }

  template<class T1,class T3> void conj_transpose_multiply_(Matrix<T1>& d,const rmatrix& a,const Matrix<T3>& b) { transpose_multiply(d,a,b); }
  template<class T1,class T3> void conj_transpose_multiply_(Matrix<T1>& d,const cmatrix& a,const Matrix<T3>& b) { conj_transpose_multiply(d,a,b); }

  template<class Ta,class Tb,class Tc> void transform_MN(Matrix<Ta>& out, const Matrix<Tb>& in,const Matrix<Tc>& VR,int m,const Matrix<Tc>& VC,int n,Matrix<Ta>& tmp)
  {
    const size_t c=in.cols();
    const size_t r=in.rows();
    const range blockm(m*c,c*(m+1)-1);
    const range blockn(n*r,r*(n+1)-1);
    conj_transpose_multiply_(tmp,VR.rows(blockn),in);
    multiply(out,tmp,VC.rows(blockm));
  }    

  template<class Td,class T> void transform_MN(Matrix<Td>& out, int N,const Matrix<T>& VR,int m,const Matrix<T>& VC,int n)
  {
    const int dim=VR.rows()/(2*N+1);
    const range blockm(m*dim,dim*(m+1)-1);
    const range blockn(n*dim,dim*(n+1)-1);
    conj_transpose_multiply_(out,VR.rows(blockn),VC.rows(blockm));
  }    

  template<class Ta,class Tb> void make_detectm(Matrix<Ta>& detectt,const Matrix<Tb>& detect,int m,int N)
  {
    multiply(detectt,2*N-abs(m)+1,detect);
  }

  template<class Ta,class Tb,class Tc> void make_detectm(Matrix<Ta>& detectt,const Matrix<Tb>& detect,const Matrix<Tc>& VR,const Matrix<Tc>& VC,int m,int N,Matrix<Ta>& tmp,Matrix<Ta>& tmp2,int verbose)
  {
    //construct transformed detection matrix
    int stoff,endoff;
    if (m<0) {
      stoff=-m;
      endoff=2*N;
    }
    else {
      stoff=0;
      endoff=2*N-m;
    }
    bool isfirst=true;
    for (int offs=stoff;offs<=endoff;offs++) {
      if (isfirst) {
	transform_MN(detectt,detect,VR,offs,VC,offs+m,tmp);
	isfirst=false;
      }
      else {
	transform_MN(tmp2,detect,VR,offs,VC,offs+m,tmp);
	detectt+=tmp2;
      }
    }
    if (verbose)
      std::cout << "Detect(m=" << m <<") in eigenbasis:\n" << detectt << std::endl;
  }

  template<class T> void make_detectm(Matrix<T>& detectt,int rs,const Matrix<T>& VR,const Matrix<T>& VC,int m,int N,Matrix<T>& tmp,int verbose)
  {
    //construct transformed detection matrix
    int stoff,endoff;
    if (m<0) {
      stoff=-m;
      endoff=2*N;
    }
    else {
      stoff=0;
      endoff=2*N-m;
    }
    bool isfirst=true;
    for (int offs=stoff;offs<=endoff;offs++) {
      if (isfirst) {
	transform_MN(detectt,rs,VR,offs,VC,offs+m);
	isfirst=false;
      }
      else {
	transform_MN(tmp,rs,VR,offs,VC,offs+m);
	detectt+=tmp;
      }
    }
    if (verbose)
      std::cout << "Detect(m=" << m <<") in eigenbasis:\n" << detectt << std::endl;
  }

  void make_detectm(RCmatrix& detectt,const rmatrix& detect,const RCmatrix& VR,const RCmatrix& VC,int n,int N,RCmatrix& tmp,RCmatrix& tmp2,int verbose)
  {
    switch (VR.type()) {
    case RCmatrix::NONE:
      detectt.set_real()=detect;
      break;
    case RCmatrix::REAL:
      make_detectm(detectt.set_real(),detect,VR.get_real(),VC.get_real(),n,N,tmp.set_real(),tmp2.set_real(),verbose); 
      break;
    case RCmatrix::COMPLEX:
      make_detectm(detectt.set_complex(),detect,VR.get_complex(),VC.get_complex(),n,N,tmp.set_complex(),tmp2.set_complex(),verbose);
      break;
    }
  }

  void make_detectm(RCmatrix& detectt,const cmatrix& detect,const RCmatrix& VR,const RCmatrix& VC,int n,int N,RCmatrix& tmp,RCmatrix& tmp2,int verbose)
  {
    cmatrix& d=detectt.set_complex();
    cmatrix& tmpc=tmp.set_complex();
    cmatrix& tmp2c=tmp2.set_complex();
    if (VR.type()!=VC.type())
      throw Mismatch("FloquetSpectrum: mixed real/complex Hamiltonians");
    switch (VR.type()) {
    case RCmatrix::NONE:
      d=detect; break;
    case RCmatrix::REAL:
      make_detectm(d,detect,VR.get_real(),VC.get_real(),n,N,tmpc,tmp2c,verbose);
      break;
    case RCmatrix::COMPLEX:
      make_detectm(d,detect,VR.get_complex(),VC.get_complex(),n,N,tmpc,tmp2c,verbose);
      break;
    }
  }

  void make_detectm(RCmatrix& detectt,const RCmatrix& detect,int n,int N)
  {
    switch (detect.type()) {
    case RCmatrix::REAL:
      make_detectm(detectt.set_real(),detect.get_real(),n,N);
      break;
    case RCmatrix::COMPLEX:
      make_detectm(detectt.set_complex(),detect.get_complex(),n,N);
      break;
    default:
      throw Failed("FloquetSpectrum: detection matrix unset");
    }
  }

  void make_detectm(RCmatrix& detectt,int rs,const RCmatrix& VR,const RCmatrix& VC,int n,int N,RCmatrix& tmp,int verbose)
  {
    switch (VR.type()) {
    case RCmatrix::NONE:
      throw Failed("make_detectm");
    case RCmatrix::REAL:
      make_detectm(detectt.set_real(),rs,VR.get_real(),VC.get_real(),n,N,tmp.set_real(),verbose);
      break;
    case RCmatrix::COMPLEX:
      make_detectm(detectt.set_complex(),rs,VR.get_complex(),VC.get_complex(),n,N,tmp.set_complex(),verbose);
      break;
    }
  }

  void BaseFloquetSpectrum::update_detect()
  {
    iter_.reset(reigs.length(),ceigs.length(),false);
    if (isgamma && (Nval!=0))
      throw Failed("This shouldn't happen!");
    const RCmatrix& VR=row().V;
    const RCmatrix& VC=col().V;
    
    if (!detect) {
      switch (VR.type()) {
      case RCmatrix::REAL:
	make_detectm(detectt.set_real(),rs,VR.get_real(),VC.get_real(),Nval,N,tmp.set_real(),verbose_);
	break;
      case RCmatrix::COMPLEX:
	make_detectm(detectt.set_complex(),rs,VR.get_complex(),VC.get_complex(),Nval,N,tmp.set_complex(),verbose_);
	break;
      default:
	throw Failed("FloquetSpectrum: Hamiltonians not set");
      }
    }
    else {
      switch (detect.type()) {
      case RCmatrix::REAL:
	make_detectm(detectt,detect.get_real(),VR,VC,Nval,N,tmp,tmp2,verbose_);
	break;
      case RCmatrix::COMPLEX:
	make_detectm(detectt,detect.get_complex(),VR,VC,Nval,N,tmp,tmp2,verbose_);
	break;
      default:
	throw Failed("FloquetSpectrum: observation matrix not set");
      }
    }
    f=Nval*freq;
  }

  template <class T> void doaddFloquetFID(BaseList<complex>& FIDv,double scale,const BaseList<complex>& UR,const BaseList<complex>& UC,const Matrix<T>& sigma0t,const Matrix<T>& detectt,double mfreqdt,int verbose)
  {
    const size_t npts=FIDv.length();
    const size_t nFr=UR.length();
    const size_t nFc=UC.length();
    double tamp=0.0;

    const complex ph=expi(2.0*M_PI*mfreqdt);
    for (size_t i=0;i<nFr;i++) {
      for (size_t j=0;j<nFc;j++) {
	complex amp=detectt(i,j)*sigma0t(j,i);
	if (norm(amp)>=1e-12) {
	  amp*=scale;
	  const complex v=conj_multiply(UR(i),UC(j))*ph;
	  if (verbose>1) {
	    std::cout << i << "," << j << " adding " << amp << " with phase factor " << v << '\n';
	    tamp+=norm(amp);
	  }
	  for (size_t k=0;k<npts;k++) {
	    FIDv(k)+=amp;
	    amp*=v;
	  }
	}
      }
    }
    if (verbose>1)
      std::cout << "Total amp: " << tamp << '\n';
  }

  template <class Td,class T> void doaddFloquetFID_hermitian(BaseList<Td>& FIDv,Td scale,const BaseList<complex>& U,const Matrix<T>& sigma0t,const Matrix<T>& detectt,int verbose)
  {
    const size_t npts=FIDv.length();
    const size_t nF=U.length();

    Td zerofreq(0.0);
    size_t i;
    for (i=nF;i--;) zerofreq+=real_(detectt(i,i)*sigma0t(i,i));
    zerofreq*=scale;
    //double tamp=fabs(zerofreq);

    scale*=2.0;
    if (zerofreq!=0.0)
      FIDv+=zerofreq;

    for (i=1;i<nF;i++) {
      for (size_t j=i;j--;) {
	Td amp(real_(detectt(i,j)*sigma0t(j,i)));
	if (norm_(amp)>=1e-12) {
	  amp*=scale;
	  complex camp(amp);
	  const complex v=conj_multiply(U(i),U(j));
	  if (verbose>1) {
	    std::cout << i << "," << j << " adding " << amp << " with phase factor " << v << '\n';
	    //	    tamp+=fabs(amp);
	  }
	  for (size_t k=0;k<npts;k++) {
	    FIDv(k)+=real(camp);
	    camp*=v;
	  }
	}
      }
    }
    //    if (verbose>1)
    //  std::cout << "Total amp: " << tamp << '\n';
  }

  //Based on FID_vega of GAMMA
  template<class TH,class TM> void BaseFloquetFID::add_FID__(BaseList<complex>& FIDv,double scale,const Matrix<TH>& VR,const Matrix<TH>& VC,const Matrix<TM>& sigma0,const Matrix<TM>& detect)
  {
    if (!arematching(sigma0,detect))
      throw Mismatch("Floquet::observe");
    
    if ( (floquet_check(VR.rows(),sigma0.rows())!=N) || (floquet_check(VC.cols(),detect.rows())!=N))
      throw Mismatch("add_FloquetFID");
    
    typedef typename promote_trait<TH,TM>::value_type TR;
    Matrix<TR> sigma0t,tmp,tmp2,detectt;
    
    transform_MN(sigma0t,sigma0,VR,N,VC,N,tmp);

    if (verbose_>1) {
      std::cout << "Transform matrix (row):\n" << VR << std::endl;
      std::cout << "sigma0 in eigenbasis:\n" << sigma0t << std::endl;
    }
    if (freqdt==0.0)
      freqdt=1.0/(2*N+1);
    
    const BaseList<complex>& UR=row().U;
    const BaseList<complex>& UC=col().U;

    if (isgamma) {
      make_detectm(detectt,detect,VR,VC,0,N,tmp,tmp2,verbose_);
      doaddFloquetFID(FIDv,scale,UR,UC,sigma0t,detectt,0,verbose_);
    }
    else {
      for (int m=-N;m<=N;m++) {
	make_detectm(detectt,detect,VR,VC,m,N,tmp,tmp2,verbose_);
	doaddFloquetFID(FIDv,scale,UR,UC,sigma0t,detectt,m*freqdt,verbose_);
      }
    }
  }
      
  template<class Td,class TH,class TM> void GammaFloquetFID::add_FID_hermitian__(BaseList<Td>& FIDv,Td scale,const Matrix<TH>& V,const Matrix<TM>& sigma0,const Matrix<TM>& detect)
  {
    if (!isdiagonal()) 
      throw Failed("add_FID_hermitian: only valid for diagonal blocks");
    if (!arematching(sigma0,detect))
      throw Mismatch("Floquet::observe");
    if (floquet_check(V.rows(),sigma0.rows())!=N)
      throw Mismatch("Floquet::observe");

    typedef typename promote_trait<TH,TM>::value_type TR;
    Matrix<TR> sigma0t,tmp,detectt,tmp2;
    
    transform_MN(sigma0t,sigma0,V,N,V,N,tmp);
    make_detectm(detectt,detect,V,V,0,N,tmp,tmp2,verbose_);
    doaddFloquetFID_hermitian(FIDv,scale,row().U,sigma0t,detectt,verbose_);
  }
      
  template<class T> void BaseFloquetFID::add_FID__(BaseList<complex>& FIDv,double scale,const Matrix<T>& VR,const Matrix<T>& VC)
  {
    Matrix<T> sigma0t,detectt,tmp;
    transform_MN(sigma0t,N,VR,N,VC,N);

    if (verbose_) {
      std::cout << "Transform matrix (row):\n" << VR << std::endl;
      std::cout << "sigma0 in eigenbasis:\n" << sigma0t << std::endl;
    }
    if (freqdt==0.0)
      freqdt=1.0/(2*N+1);
    const BaseList<complex>& UR=row().U;
    const BaseList<complex>& UC=col().U;
    const size_t dim=UR.length();
    if (dim!=UC.length())
      throw Mismatch("FloquetFID::add_FID");

    if (isgamma) {
      make_detectm(detectt,dim,VR,VC,0,N,tmp,verbose_);
      doaddFloquetFID(FIDv,scale,UR,UC,sigma0t,detectt,0,verbose_);
    }
    else {
      for (int m=-N;m<=N;m++) {
	make_detectm(detectt,dim,VR,VC,m,N,tmp,verbose_);
	doaddFloquetFID(FIDv,scale,UR,UC,sigma0t,detectt,m*freqdt,verbose_);
      }
    }
  }
      
template<class T> void BaseFloquetFID::add_FID_(BaseList<complex> FIDv,double scale,const Matrix<T>& sigma0, const Matrix<T>& detect)
{
  const RCmatrix& VR=row().V;
  const RCmatrix& VC=col().V;
  switch (VR.type()) {
  case RCmatrix::REAL: add_FID__(FIDv,scale,VR.get_real(),VC.get_real(),sigma0,detect); break;
  case RCmatrix::COMPLEX: add_FID__(FIDv,scale,VR.get_complex(),VC.get_complex(),sigma0,detect); break;
  default: throw Failed("FloquetFID::add_FID_");
  }
}

  template<class Td, class T> void GammaFloquetFID::add_FID_hermitian_(BaseList<Td> FIDv, Td scale,const Matrix<T>& sigma0, const Matrix<T>& detect)
{
  if (!isdiagonal()) 
    throw Failed("add_FID_hermitian: only valid for diagonal blocks");
  const RCmatrix& V=row().V;
  switch (V.type()) {
  case RCmatrix::REAL: add_FID_hermitian__(FIDv,scale,V.get_real(),sigma0,detect); break;
  case RCmatrix::COMPLEX: add_FID_hermitian__(FIDv,scale,V.get_complex(),sigma0,detect); break;
  default:
    throw Failed("GammaFloquetFID::add_FID_hermitian_");
  }
}

void BaseFloquetFID::add_FID(BaseList<complex> FIDv,double scale)
{
  if (!isdiagonal()) 
    throw Failed("add_FID: not valid for diagonal blocks");
  const RCmatrix& VR=row().V;
  const RCmatrix& VC=col().V;
  switch (VR.type()) {
  case RCmatrix::REAL: add_FID__(FIDv,scale,VR.get_real(),VC.get_real()); break;
  case RCmatrix::COMPLEX: add_FID__(FIDv,scale,VR.get_complex(),VC.get_complex()); break;
  default: throw Failed("FloquetFID::add_FID_");
  }
}

//   template<class TH,class TM> void add_FloquetFID_(BaseList<complex>& FID, double scale, const Matrix<TH>& FHam, const Matrix<TM>& sigma0, const Matrix<TM>& detect,double freq,double dt,bool isgamma)
//   {
//     Matrix<TH> V;
//     const size_t n=FHam.rows();
//     ScratchList<double> eigs(n);
//     hermitian_eigensystem(V,eigs,FHam);
//     ScratchList<complex> U(n);
//     propagator(U,eigs,dt);

// #ifndef NDEBUG
//     std::cout << "Floquet eigenvalues: " << eigs << '\n' << "Floquet phase factors: " << U << '\n';
// #endif

//     add_FloquetFID__(FID,scale,V,V,U,U,sigma0,detect,freq*dt,isgamma);
//    }

//   void add_FloquetFID(BaseList<complex>& FID, double scale, const cmatrix& FHam, const cmatrix& sigma0, const cmatrix& detect,double freq,double dt)
//   { add_FloquetFID_(FID,scale,FHam,sigma0,detect,freq,dt,false); }

//   void add_FloquetFID(BaseList<complex>& FID, double scale, const rmatrix& FHam, const cmatrix& sigma0, const cmatrix& detect,double freq,double dt)
//   { add_FloquetFID_(FID,scale,FHam,sigma0,detect,freq,dt,false); }
  
//   void add_FloquetFID(BaseList<complex>& FID, double scale, const cmatrix& FHam, const rmatrix& sigma0, const rmatrix& detect,double freq,double dt)
//   { add_FloquetFID_(FID,scale,FHam,sigma0,detect,freq,dt,false); }

//   void add_FloquetFID(BaseList<complex>& FID, double scale, const rmatrix& FHam, const rmatrix& sigma0, const rmatrix& detect,double freq,double dt)
//   { add_FloquetFID_(FID,scale,FHam,sigma0,detect,freq,dt,false); }

//   template<class T> void add_FloquetFID_(BaseList<complex>& FID, double scale, const cmatrix& FU, const Matrix<T>& sigma0, const Matrix<T>& detect,double freqdt,bool isgamma)
//   {
//     cmatrix V;
//     ScratchList<complex> U(FU.rows());
//     eigensystem(V,U,FU);
// #ifndef NDEBUG
//     std::cout << "Floquet phase factors: " << U << '\n';
// #endif
//     add_FloquetFID__(FID,scale,V,V,U,U,sigma0,detect,freqdt,isgamma);
//    }

//   void add_FloquetFID(BaseList<complex>& FID, double scale, const cmatrix& FU, const cmatrix& sigma0, const cmatrix& detect,double freqdt)
//   { add_FloquetFID_(FID,scale,FU,sigma0,detect,freqdt,false); }

//   void add_FloquetFID(BaseList<complex>& FID, double scale, const cmatrix& FU, const rmatrix& sigma0, const rmatrix& detect,double freqdt)
//   { add_FloquetFID_(FID,scale,FU,sigma0,detect,freqdt,false); }

  void transform_MN(RCmatrix& d,const cmatrix& in,const RCmatrix& VR,int m,const RCmatrix& VC,int n,RCmatrix& tmp)
  {
    cmatrix& dc=d.set_complex();
    cmatrix& tmpc=tmp.set_complex();
    if (VR.type()!=VC.type())
      throw Mismatch("row/column have different types");
    switch (VR.type()) {
    case RCmatrix::NONE:
      dc=in; break;
    case RCmatrix::REAL:
      transform_MN(dc,in,VR.get_real(),m,VC.get_real(),n,tmpc); break;
    case RCmatrix::COMPLEX:
      transform_MN(dc,in,VR.get_complex(),m,VC.get_complex(),n,tmpc); break;
    }
  }

  void transform_MN(RCmatrix& d,int dim,const RCmatrix& VR,int m,const RCmatrix& VC,int n)
  {
    if (VR.type()!=VC.type())
      throw Mismatch("row/column have different types");
    switch (VR.type()) {
    case RCmatrix::NONE:
      throw Failed("This combination shouldn't happen!");
    case RCmatrix::REAL:
      transform_MN(d.set_real(),dim,VR.get_real(),m,VC.get_real(),n); break;
    case RCmatrix::COMPLEX:
      transform_MN(d.set_complex(),dim,VR.get_complex(),m,VC.get_complex(),n); break;
    }
  }

  void transform_MN(RCmatrix& d,const rmatrix& in,const RCmatrix& VR,int m,const RCmatrix& VC,int n,RCmatrix& tmp)
  {
    if (VR.type()!=VC.type())
      throw Mismatch("row/column have different types");
    switch (VR.type()) {
    case RCmatrix::NONE:
      d.set_real()=in; break;
    case RCmatrix::REAL:
      transform_MN(d.set_real(),in,VR.get_real(),m,VC.get_real(),n,tmp.set_real()); break;
    case RCmatrix::COMPLEX:
      transform_MN(d.set_complex(),in,VR.get_complex(),m,VC.get_complex(),n,tmp.set_complex()); break;
    }
  }

  template<class T> void BaseFloquetSpectrum::observe_(const Matrix<T>& sigma0,const Matrix<T>& detectv)
  {
    if (!N)
      throw Failed("FloquetSpectrum: Hamiltonian not set");
    if (!arematching(sigma0,detectv))
      throw Mismatch("FloquetSpectrum::observe");
    rs=sigma0.rows(); cs=sigma0.cols();
    if ((floquet_check(rows(),rs)!=N) || (floquet_check(cols(),cs)!=N))
      throw Mismatch("Floquet dimensions don't match");
    transform_MN(sigma0t,sigma0,row().V,N,col().V,N,tmp);
    detect=detectv;
    reset();
  }

template void BaseFloquetSpectrum::observe_(const cmatrix&, const cmatrix&);
template void BaseFloquetSpectrum::observe_(const rmatrix&, const rmatrix&);

  void BaseFloquetSpectrum::observe()
  {
    if (!N)
      throw Failed("FloquetSpectrum: Hamiltonian not set");
    if (isdiagonal())
      throw Failed("Not valid for diagonal blocks");
    if (cols()!=rows())
      throw Mismatch("Dimensions don't match");
    rs=cs=rows();
    transform_MN(sigma0t,N,row().V,N,col().V,N);
    detect.clear();
    reset();
  }

  bool BaseFloquetSpectrum::operator() (complex& amp,double& lfreq)
  {
    if (finished)
      return false;

    double v;

    const bool iscomplex=detectt.iscomplex();
    for (;;) {
      bool isOK;
      if (iscomplex) {
	amp=(detectt.get_complex())(iter_.r,iter_.s)*(sigma0t.get_complex())(iter_.s,iter_.r);
	isOK=(norm(amp)>1e-12);
      }
      else {
	amp=(v=(detectt.get_real())(iter_.r,iter_.s)*(sigma0t.get_real())(iter_.s,iter_.r));
	isOK=(fabs(v)>1e-6);
      }
      if (isOK)
	lfreq=ceigs(iter_.s)-reigs(iter_.r)-f;

      if (!iter_.advance()) {
        if ((Nval==N) || isgamma) {
	  finished=true;
	  return isOK;
	}
	Nval++;
	update_detect();
      }
      else {
	if (isOK)
	  return true;
      }
    }    
  }
    
void BaseFloquetSpectrum::reset()
{ 
  reigs.create(row().eigs); ceigs.create(col().eigs);
  Nval=isgamma ? 0 : -N;
  update_detect();
}

 //  template<class TH,class TM> void add_FloquetSpectrum__(BaseHistogram<complex>& Spec,double scale,const Matrix<TH>& V, const BaseList<double>& eigs,const Matrix<TM>& sigma0,const Matrix<TM>& detect,double freq)
//   {
//     if (!issquare(sigma0) || !issquare(detect)) throw NotSquare("add_FloquetFID");

//     const int nh=sigma0.rows();
//     if (nh!=detect.rows()) throw Mismatch("add_FloquetFID");
//     const int nF=V.rows();
//     const int N=floquet_check(nF,nh);

//     typedef typename promote_trait<TH,TM>::value_type TR;
//     Matrix<TR> tmp,tmp2,sigma0t,detectt;

//     transform_MN(sigma0t,sigma0,V,N,V,N,tmp);

// #ifndef NDEBUG
//     std::cout << "sigma0 in eigenbasis:\n" << sigma0t << std::endl;
// #endif

//     for (int m=-N;m<=N;m++) {
//       make_detectm(detectt,detect,V,V,m,N,tmp,tmp2);
//       const double f=m*freq;
//       for (int i=0;i<nF;i++) {
// 	for (int j=0;j<nF;j++) {
// 	  complex amp=detectt(i,j)*sigma0t(j,i);
// 	  if (norm(amp)>=1e-12) {
// 	    amp*=scale;
// 	    const double tfreq=eigs(j)-eigs(i)-f;
// #ifndef NDEBUG
// 	    std::cout << i << "," << j << " adding " << amp << " of " << tfreq << '\n';
// #endif
// 	    Spec.add(amp,tfreq);
// 	  }
// 	}
//       } 
//     }
//   }

//   template<class TH,class TM> void add_FloquetSpectrum_(BaseHistogram<complex>& Spec,double scale,const Matrix<TH>& FHam, const Matrix<TM>& sigma0, const Matrix<TM>& detect,double freq)
//   {
//     Matrix<TH> V;
//     ScratchList<double> eigs(FHam.rows());
//     hermitian_eigensystem(V,eigs,FHam);

//     add_FloquetSpectrum__(Spec,scale,V,eigs,sigma0,detect,freq);
//    }

//   void add_FloquetSpectrum(BaseHistogram<complex>& Spec,double scale,const cmatrix& FHam, const cmatrix& sigma0, const cmatrix& detect,double freq)
//   { add_FloquetSpectrum_(Spec,scale,FHam,sigma0,detect,freq); }

//   //  void add_FloquetSpectrum(BaseHistogram<complex>& Spec,double scale,const rmatrix& FHam, const cmatrix& sigma0, const cmatrix& detect,double freq)
//   //  { add_FloquetSpectrum_(Spec,scale,FHam,sigma0,detect,freq); }

//   void add_FloquetSpectrum(BaseHistogram<complex>& Spec,double scale,const cmatrix& FHam, const rmatrix& sigma0, const rmatrix& detect,double freq)
//   { add_FloquetSpectrum_(Spec,scale,FHam,sigma0,detect,freq); }

//   template<class T> void add_FloquetSpectrum_(BaseHistogram<complex>& Spec,double scale,const cmatrix& FU, const Matrix<T>& sigma0, const Matrix<T>& detect,double freq,double dt)
//   {
//     cmatrix V;
//     ScratchList<complex> U(FU.rows());
//     eigensystem(V,U,FU);
//     ScratchList<double> eigs(U.length());
//     diag_propagator(eigs,U,dt);

//     add_FloquetSpectrum__(Spec,scale,V,eigs,sigma0,detect,freq);
//   }

//  void add_FloquetSpectrum(BaseHistogram<complex>& Spec,double scale,const cmatrix& FU, const cmatrix& sigma0, const cmatrix& detect,double freq,double dt)
//   { add_FloquetSpectrum_(Spec,scale,FU,sigma0,detect,freq,dt); }

//  void add_FloquetSpectrum(BaseHistogram<complex>& Spec,double scale,const cmatrix& FU, const rmatrix& sigma0, const rmatrix& detect,double freq,double dt)
//   { add_FloquetSpectrum_(Spec,scale,FU,sigma0,detect,freq,dt); }

void BaseFloquetSpectrum::set_H(const cmatrix& H, int order)
{
  get_frequency(H,order);
  setdiagonal()=H;
  N=order;
}

void BaseFloquetSpectrum::set_H(char sel,const cmatrix& H, int order)
{
  get_frequency(H,order);
  setRC(sel)=H;
  N=order;
}

void BaseFloquetSpectrum::set_H(const rmatrix& H, int order)
{
  get_frequency(H,order);
  setdiagonal()=H;
  N=order;
}

void BaseFloquetSpectrum::set_H(char sel,const rmatrix& H, int order)
{
  get_frequency(H,order);
  setRC(sel)=H;
  N=order;
}

size_t hilbert_order(size_t totdim,int order)
{
  const size_t fsize=2*order+1;
  if ((order<1) || (totdim % fsize))
    throw Failed("hilbert_order: invalid Floquet matrix/order");
  return totdim/fsize;
}

template<class T> inline double get_frequency_(const Matrix<T>& H,int order)
{
  if (!issquare(H)) 
    throw NotSquare("get_frequency");
  const size_t dim=hilbert_order(H.rows(),order);
  const size_t base0=order*dim;
  return real_(H(base0,base0)-H(base0-dim,base0-dim));
}

template<class T> inline void set_frequency_(Matrix<T>& H,int N,double freq)
{
  if (!issquare(H)) 
    throw NotSquare("set_frequency");
  const size_t n=hilbert_order(H.rows(),N);
  IndirectList<T,slice> Hdiag(H.diag());
  const size_t base0=N*n;
  for (size_t j=n;j--;) {
    const size_t basej=base0+j;
    const double Hval=real_(Hdiag(basej));
    for (int r=1;r<=N;r++) {
      const double freqv=freq*r;
      Hdiag(basej+r*n)=Hval+freqv;
      Hdiag(basej-r*n)=Hval-freqv;
    }
  }
}

double get_frequency(const rmatrix& H,int order) { return get_frequency_(H,order); }
double get_frequency(const cmatrix& H,int order) { return get_frequency_(H,order); }

void set_frequency(rmatrix& H,int order,double freq) { return set_frequency_(H,order,freq); }
void set_frequency(cmatrix& H,int order,double freq) { return set_frequency_(H,order,freq); }

void BaseFloquetFID::add_FID(BaseList<complex> FIDv,double scale,const rmatrix& sigma0, const rmatrix& detect)
{
  add_FID_(FIDv,scale,sigma0,detect);
}

void BaseFloquetFID::add_FID(BaseList<complex> FIDv,double scale,const cmatrix& sigma0, const cmatrix& detect)
{
  add_FID_(FIDv,scale,sigma0,detect);
}

void GammaFloquetFID::add_FID_hermitian(BaseList<double> FIDv,double scale,const rmatrix& sigma0, const rmatrix& detect)
{
  add_FID_hermitian_(FIDv,scale,sigma0,detect);
}

void GammaFloquetFID::add_FID_hermitian(BaseList<double> FIDv,double scale,const cmatrix& sigma0, const cmatrix& detect)
{
  add_FID_hermitian_(FIDv,scale,sigma0,detect);
}

void GammaFloquetFID::add_FID_hermitian(BaseList<complex> FIDv,complex scale,const rmatrix& sigma0, const rmatrix& detect)
{
  add_FID_hermitian_(FIDv,scale,sigma0,detect);
}

void GammaFloquetFID::add_FID_hermitian(BaseList<complex> FIDv,complex scale,const cmatrix& sigma0, const cmatrix& detect)
{
  add_FID_hermitian_(FIDv,scale,sigma0,detect);
}

std::ostream& operator<< (std::ostream& ostr, const BaseFloquetSpectrum& a)
{
  ostr << static_cast<const SwapStore<FloquetFDStash>&>(a);
  if (!a.finished) {
    ostr << "Left eigenvalues: " << a.reigs << '\n';
    ostr << "Right eigenvalues: " << a.ceigs << '\n';
  }
  return ostr << a.iter_;
}

}//namespace libcmatrix
