//#define LCM_USE_SSECOMPLEX
#include <iostream>
#include <complex>
#include "cmatrix_utils.h"
#include "cmatrix_complex.h"

//#define CTYPE std::complex<double>
//#define CTYPE libcmatrix::complex
//#define CTYPE SSEcomplex

// struct padded {
//   char x;
//   const libcmatrix::complex z;
//   padded(char xv, const libcmatrix::complex& zv) : x(xv), z(zv) {}
// };

// libcmatrix::complex paddeddouble_(const padded z)
// {
//   return (z.z)*2.0;
// }

// libcmatrix::complex paddeddouble(const libcmatrix::complex &z)
// {
//   padded padz('X',z);
//   return paddeddouble_(padz);
// }

using namespace libcmatrix;

complex adddouble(char,const complex a)
{
  return a;
}

namespace libcmatrix {

template<> struct doesconj_ip< std::complex<double> >
 { void operator()( std::complex<double>& v) const { v= std::complex<double>(v.real(),-v.imag()); } };

}

template<typename CTYPE> void testextra(CTYPE&, const CTYPE&, double, Bool2Type<false>)  {}
template<typename CTYPE> void testextra(CTYPE& a, const CTYPE& b, double S3, Bool2Type<true>) 
{
  //  const CTYPE da=paddeddouble(a); //!< try to force unaligned stack
  //std::cout << "2*a: " << da << '\n';
  const CTYPE added=adddouble('C',a);
  std::cout << "Add double: " << added << '\n';

  //CTYPE as[4];
  //BaseList<CTYPE> alist(3,as);
  // BaseList<CTYPE> unalist(3,(CTYPE*)(2+(char*)as));
  double bs[4];
  BaseList<double> blist(3,bs);
  // BaseList<double> unblist(3,(double*)(16+(char*)bs));
  //std::cout << "unblist: " << (unblist.vector()) << '\n';
  //alist=CTYPE(2.0,1.0);
  blist=double(3.0);
  List<CTYPE> listr=a*blist;
  std::cout << "listr: " << listr.vector() << '\n';
  std::cout << "alist * unblist: " << listr << '\n';

  std::cout << "real(a*b): " << real_multiply(a,b) << "  " << real(a*b) << '\n';
  std::cout << "a*conj(b): " << multiply_conj(a,b) << '\n';
  std::cout << "S*conj(a): " << multiply_conj(S3,a) << '\n';
  std::cout << "real(conj(a)*b): " << real_conj_multiply(a,b) << "  " << real(conj_multiply(a,b)) << '\n';

  mla(a,b,b);
  std::cout << "a+=b*b: " << a << '\n';
  mla(a,S3,b);
  std::cout << "a+=S3*b: " << a << '\n';
  
  real(a,5.0);
  imag(a,3.0);
  std::cout << "a.real(5.0), a.imag(3.0): " << a << '\n';

  std::cout << "a==a: " << (a==a) << '\n';
  std::cout << "a==b: " << (a==b) << '\n';
  std::cout << "a!=a: " << (a!=a) << '\n';
  std::cout << "a!=b: " << (a!=b) << '\n';

  static const CTYPE complex_zero(0.0);
  const CTYPE pureimag(0.0,1.0);
  std::cout << "a==0: " << (a==complex_zero) << '\n';
  std::cout << "a!=0: " << (a!=complex_zero) << '\n';
  std::cout << "pureimag==0: " << (pureimag==complex_zero) << '\n';
  std::cout << "pureimag!=0: " << (pureimag!=complex_zero) << '\n';
  
  CTYPE d(0.0);
  mla(d,a,b);
  std::cout << "mla(d,a,b): " << d << '\n';
}

template<typename CTYPE, bool Extra> void testcomplex(const char* name)
{
  std::cout << "\nTesting " << name << "\n";
  CTYPE a(2.0,3.0);
  CTYPE b(3.0,-1.0);
  CTYPE a1(1.0,2.0);
  CTYPE b1(4.0);
  CTYPE a2(3.0,4.0);
  CTYPE b2(-1.0);
  CTYPE a3(-2.0,1.0);
  double rawds[2];
  const double S3=2.0*random(1.0);
  double* Sunalignedp=(double*)(2+(char*)rawds);
  *Sunalignedp=S3;
  std::cout << "Aligned double: " << (&S3) << "  Unaligned: " << (Sunalignedp) << "\n";
  
  std::cout << "a: " << a << '\n';
  std::cout << "b: " << b << '\n';
  //  a=conj(a);
  doesconj_ip<CTYPE> doesconj;
  doesconj(a);
  std::cout << "conj(a): " << a << '\n';
  doesconj(a);
  //  a=conj(a);
  a=-a;
  std::cout << "-a: " << a << '\n';
  a=-a;  
  std::cout << "a+b: " << (a+b) << '\n';
  std::cout << "a-b: " << (a-b) << '\n';
  std::cout << "a*b: " << (a*b) << '\n';
  std::cout << "a/b: " << (a/b) << '\n';
  std::cout << "norm(a): " << norm(a) << '\n';
 
  std::cout << "S (random start): " << S3 << '\n';
  const CTYPE asum(a+*Sunalignedp);
  //const CTYPE asum(a+S3);
  std::cout << "a+Su: " << asum << "  S+a: " << (S3+a) << '\n';
  std::cout << "a-Su: " << (a-*Sunalignedp) << "  S-a: " << (S3-a) << '\n';
  std::cout << "a*Su: " << (a*(*Sunalignedp)) << "  S*a: " << (S3*a) << '\n';
  std::cout << "a/Su: " << (a/(*Sunalignedp)) << "  S/a: " << (S3/a) << '\n';

  testextra(a,b,S3,Bool2Type<Extra>());
//   libcmatrix::timer<libcmatrix::WallTimer> stopwatch;

//   for (unsigned long long i=3000000000ull;i--;) {
//     a+=b;
//     a1+=b1;
//     a2+=b2;
//     a3+=S3;
//   }
//   std::cout << "Time taken: " << stopwatch() << " s\n";
//   std::cout << a << ' ' << a1 << ' ' << a2 << ' ' << a3 << std::endl;
}

int main()
{	
#ifdef LCM_USE_SSECOMPLEX
	std::cout << "Using SSE complex: Yes\n";
#else
	std::cout << "Using SSE complex: No\n";
#endif
  testcomplex<libcmatrix::complex, true>("libcmatrix::complex");
  testcomplex<std::complex<double>, false >("std::complex");
}
