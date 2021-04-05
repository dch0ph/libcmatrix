#include "cmatrix.h"
#include "rmatrix.h"
#include "cmatrix_utils.h"
#include "timer.h"
#include "ttyio.h"
#include "cmatrix_external.h"

//using namespace std;
using namespace libcmatrix;

#define libcmatrix_mvmul multiply

bool showres=false;
bool issym=false;
int nmflops=0;

const char CTMstr[]="CTM";

#ifdef USE_ATLAS

static const complex one(1.0);
static const complex zero(0.0);

void lapack_hermitian_multiply(RawMatrix<complex>& d, const RawMatrix<complex>& a, const RawMatrix<complex>& b, bool acc)
{
  size_t ar,ac,br,bc;
  getrowcols(ar,ac,a,false);
  getrowcols(br,bc,b,false);

  if (ac!=br)
    throw Mismatch("lapack_hermitian_multiply",ar,ac,br,bc);

// #ifdef USE_SUNPERFACML
//   zgemm(Atrans ? 'Y' : 'N',
// 	Btrans ? 'Y' : 'N',
// 	ar,bc,ac,
// 	lapack_pass(one),
// 	lapack_pass(a),a.step(),
// 	lapack_pass(b),b.step(),
// 	lapack_pass(acc ? one : zero),
// 	lapack_pass(d),d.step());
// #else
//#ifdef USE_ATLAS
  cblas_zhemm(LCM_ORDER,
	      CblasLeft,CblasUpper,
	      ar,ac,
	lapack_pass(one),
	      lapack_pass(a),a.step(),
	      lapack_pass(b),b.step(),
	      lapack_pass(acc ? one : zero),
	      lapack_pass(d),d.step());
}
  //#endif


void lapack_hermitian_multiply(RawMatrix<double>& d, const RawMatrix<double>& a, const RawMatrix<double>& b, bool acc)
{
  size_t ar,ac,br,bc;
  getrowcols(ar,ac,a,false);
  getrowcols(br,bc,b,false);
   
  if (ac!=br)
    throw Mismatch("lapack_hermitian_multiply",ar,ac,br,bc);

  //#ifdef USE_ATLAS
  cblas_dsymm(LCM_ORDER,
	      CblasLeft,CblasUpper,
	ar,bc,
	1.0,
	      lapack_pass(a),a.step(),
	      lapack_pass(b),b.step(),
	      acc ? 1.0 : 0.0,
	      lapack_pass(d),d.step());
}

void lapack_hermitian_square(RawMatrix<double>& d, const RawMatrix<double>& a, bool acc)
{
  //size_t ar,ac;
  //getrowcols(ar,ac,a,false);
   
  if (a.rows()!=a.cols())
    throw NotSquare("lapack_hermitian_square");

  //#ifdef USE_ATLAS
  cblas_dsyrk(LCM_ORDER,
	      CblasUpper,
	      CblasNoTrans,
	      a.rows(),a.cols(),
	1.0,
	      lapack_pass(a),a.step(),
	      acc ? 1.0 : 0.0,
	      lapack_pass(d),d.step());
}

void lapack_hermitian_square(RawMatrix<complex>& d, const RawMatrix<complex>& a, bool acc)
{
  //size_t ar,ac;
  //getrowcols(ar,ac,a,false);
   
  if (a.rows()!=a.cols())
    throw NotSquare("lapack_hermitian_square");

  //#ifdef USE_ATLAS
  cblas_zherk(LCM_ORDER,
	      CblasUpper,
	      CblasNoTrans,
	      a.rows(),a.cols(),
	      1.0,
	      lapack_pass(a),a.step(),
	      acc ? 1.0 : 0.0,
	      lapack_pass(d),d.step());
}
#endif

// template<class T> void fill_baseline(Matrix<T>& d)
// {
//   const int cs=d.cols();
//   for (int r=d.rows();r--;) {
//     for (int c=cs;c--;) d(r,c)=0.0;
//   }
// }

// template<class T> void fill(Matrix<T>& d)
// {
//   const int cs=d.cols();
//   const int rs=d.rows();
//   for (int r=rs;r--;) {
//     for (int c=cs;c--;) d(index(rs,r,c))=0.0;
//   }
// }

void conj_transpose_multiply_baseline(cmatrix &to,const cmatrix &T,const cmatrix &a)
{
  if (!T)
    throw Undefined(CTMstr);
  if (issame(to,T) || issame(to,a))
    throw ArgumentClash(CTMstr);

  const size_t n=a.rows();
  if (n!=T.rows()) 
    throw Mismatch(CTMstr);

  const size_t ca=a.cols();
  const size_t ct=T.cols();
  to.create(ct,ca);
  complex *destp=to.vector();

  for (size_t i=0;i<ct;i++) {
    for (size_t k=0;k<ca;k++) {
      complex sum(0.0,0.0);
      for (size_t j=n;j--;) 
	mla_conj(sum,a(j,k),T(j,i));
      *destp++=sum;
    }
  }
}

extern const char MCTstr[]="MCT";

template<class Type> void multiply_conj_transpose_baseline(cmatrix &to,const cmatrix &a,const Matrix<Type> &T)
{
  if ( &to==&a || (void *)&to==(void *)&T)
    throw ArgumentClash(MCTstr);
  const Type *sTp=T.vector();
  const complex *ap=a.vector();

  const size_t n=a.cols();
  if (n!=T.cols())
    throw Mismatch(MCTstr);

  const size_t ra=a.rows();
  const size_t rt=T.rows();
  to.create(ra,rt);

  complex *destp=to.vector();

  for (size_t i=0;i<ra;i++) {
    const Type *Tp=sTp;
    for (size_t k=0;k<rt;k++) {
      complex sum(0.0,0.0);
      for (size_t j=n;j--;) 
	mla_conj(sum,ap[j],Tp[j]);
      *destp++=sum;
      Tp+=n;
    }
    ap+=n;
  }
}

#ifdef USE_ATLAS

void lapack_copy(BaseList<complex> y, const BaseList<complex>& x)
{
  const size_t n=x.length();
  if (n!=y.length())
    throw Mismatch("lapack_copy");
  cblas_zcopy(n,x.vector(),1,y.vector(),1);
}

void lapack_copy(BaseList<double> y, const BaseList<double>& x)
{
  const size_t n=x.length();
  if (n!=y.length())
    throw Mismatch("lapack_copy");
  cblas_dcopy(n,x.vector(),1,y.vector(),1);
}

void lapack_mla(BaseList<complex> y,const complex& alpha,const BaseList<complex>& x)
{
  const size_t n=x.length();
  if (n!=y.length())
    throw Mismatch("lapack_mla");
  cblas_zaxpy(n,lapack_pass(alpha),x.vector(),1,y.vector(),1);
}

void lapack_mla(BaseList<double> y,double alpha,const BaseList<double>& x)
{
  const size_t n=x.length();
  if (n!=y.length())
    throw Mismatch("lapack_mla");
  cblas_daxpy(n,alpha,x.vector(),1,y.vector(),1);
}

void lapack_mvmul(List<complex>& d, const Matrix<complex>& a, const BaseList<complex>& b,bool Atrans =false)
{
  size_t ar,ac;
  if (Atrans) {
    ar=a.cols();
    ac=a.rows();
  }
  else {
    ar=a.rows();
    ac=a.cols();
  }
  size_t blen=b.length();

  if (ac!=blen)
    throw Mismatch("lapack_mvmul");
  d.create(ar);

  static complex one(1.0);
  static complex zero(0.0);

  cblas_zgemv(LCM_ORDER,
	      Atrans ? CblasConjTrans : CblasNoTrans,
	ar,ac,
	lapack_pass(one),
	      lapack_pass(a),a.cols(),
	lapack_pass(b),1,
	lapack_pass(zero),
	lapack_pass(d),1);
}

#endif

template<typename T> inline void multiply_baseline(Matrix<T>& d, const Matrix<T>& a,const Matrix<T>& b)
{
  RawMatrix<T> rawd(d);
  multiply_naive_direct_(rawd,RawMatrix<T>(a),RawMatrix<T>(b),false);
}

template<> inline void multiply_baseline(Matrix<complex>& d, const Matrix<complex>& a,const Matrix<complex>& b)
{
  RawMatrix<complex> rawd(d);
  multiply_naive_ref_(rawd,RawMatrix<complex>(a),RawMatrix<complex>(b),false);
}

template<class Iter_d, class Iter_a, class Iter_b,class Function> inline void apply2_baseline(Iter_d dstart, Iter_d dend, Function mop,Iter_a astart, Iter_b bstart)
{
  if (dstart==dend)
    return;
  mop(*dstart,*astart,*bstart);
  //prefer to use prefix ++ since generally faster for non-trivial iterators
  while ((++dstart)!=dend)
    mop(*dstart,*(++astart),*(++bstart));
}

template<class Iter_a, class Iter_b,class Function> inline void applyip_baseline(Function mop,Iter_a astart, Iter_a aend,Iter_b bstart)
{
  if (astart==aend)
    return;
  mop(*astart,*bstart);
  while ((++astart)!=aend)
    mop(*astart,*(++bstart));
}

template<typename T1,typename T2,typename T3,typename Function> inline void apply2_baseline(Matrix<T1>& d,Function mop,const Matrix<T2>& a, const Matrix<T3>& b)
{
  if (!arematching(a,b))
    throw Mismatch("apply2");
  d.create(a.rows(),a.cols());
  apply2_baseline(d.begin(),d.end(),mop,a.begin(),b.begin());
}

template<typename T1,typename T2,typename Function> inline void applyip_baseline(Function mop,Matrix<T1>& d, const Matrix<T2>& b)
{
  if (!d) {
    d=b;
    return;
  }
  if (!arematching(d,b)) 
    throw Mismatch("applyip");
  applyip_baseline(mop,d.begin(),d.end(),b.begin());
}

template<typename T> inline void add_baseline(Matrix<T>& d,const Matrix<T>& a, const Matrix<T>& b)
{
  apply2_baseline(d,doesadd_sd<T,T,T>(),a,b);
}

template<typename T> inline void addip_baseline(Matrix<T>& d, const Matrix<T>& b)
{
  applyip_baseline(doesadd_ip<T,T>(),d,b);
}

void fillrandom(BaseList<double> A)
{
  for (size_t i=A.length();i--;)
    A(i)=random(1.0);
}

void fillrandom(Matrix<double>& A)
{
  if (issym) {
    const size_t n=A.rows();
    for (size_t i=n;i--;) {
      A(i,i)=random(1.0);
      for (size_t j=i;j--;) 
	A(j,i)=A(i,j)=random(1.0);
    }
  }
  else
    fillrandom(A.row());
}

void fillrandom(BaseList<complex> A)
{
  for (size_t i=A.length();i--;)
    A(i)=complex(random(1.0),random(1.0));
}

void fillrandom(Matrix<complex>& A)
{
  if (issym) {
    const size_t n=A.rows();
    for (size_t i=n;i--;) {
      A(i,i)=random(1.0);
      for (size_t j=i;j--;) 
	A(j,i)=conj(A(i,j)=complex(random(1.0),random(1.0)));
    }
  }
  else
    fillrandom(A.row());
}

#define MAXOPNAME 16

enum op_t { ADD, MULTIPLY, CTM, MCT, MVMUL, MLA, COPY, XOVER };
enum mat_t { T_REAL, T_SYMMETRIC, T_COMPLEX, T_HERMITIAN };
const char* opnames[]={"","MM","TM","MT","MV","MLA", ""};

inline bool complexity(op_t op, mat_t mtype)
{
  return  (op==CTM || op==MCT) ? true : (mtype>=T_COMPLEX);
}

void getopname(char* buf,op_t op, mat_t mtype)
{
  snprintf(buf,MAXOPNAME,"%c%s",(complexity(op,mtype) ? 'Z' : 'D'),opnames[op]);
}

#define MAXMETHODS 4

struct optimise_t {
  op_t op;
  mat_t mtype;
  char opname[MAXOPNAME];
  bool iscomplex;
  bool issquare;
  bool compbaseline;
  bool self;
  char methods[MAXMETHODS+1];

  optimise_t(op_t opv, mat_t typev);
};

optimise_t::optimise_t(op_t opv, mat_t typev =T_COMPLEX)
  : op(opv), mtype(typev) 
{
  iscomplex=complexity(opv,typev);
  getopname(opname,opv,typev);
  issquare=true;
  compbaseline=false;
  self=false;
  strcpy(methods,"IE");
}

std::pair<int,int> doopt(const optimise_t& opdesc, size_t start, size_t end, size_t factor, bool isexp)
{
  const op_t op=opdesc.op;
  //  const bool isbinary=(op!=TRANSPOSE && op!=MVMUL);
  const bool isbinary=(op!=MVMUL);
  //  const bool istransp=(op==TRANSPOSE || op==MULTIPLY);
  const bool istransp=(op==MULTIPLY);

  std::cout << "n x m\ttimes\tflops\tMflps/s\t" << (opdesc.compbaseline ? "Check\n" : "\n");

  size_t xover=0;
  size_t xover_val=0;
  for (size_t n=start;n<=end;) {

    const size_t m=(opdesc.issquare ? 1 : 2)*n;
  
    Matrix<double> Ar(n,m);
    Matrix<complex> Ac(n,m);    
    Matrix<double> Br;
    Matrix<complex> Bc;
    Matrix<double> Cr(n,n);
    Matrix<complex> Cc(n,n);
    Matrix<double> Dr(n,n);
    Matrix<complex> Dc(n,n);
    List<double> Vr(m);
    List<complex> Vc(m);
    List<double> Wr;
    List<complex> Wc;

    if (istransp) {
      Br.create(m,n);
      Bc.create(m,n);
    }
    else {
      Br.create(n,m);
      Bc.create(n,m);
    }

    if (opdesc.iscomplex) {
      fillrandom(Ac);
      if (showres)
	std::cout << "A\n" << Ac << std::endl;
      if (isbinary) {
	if (opdesc.self)
	  Bc=Ac;
	else {
	  fillrandom(Bc);
	  if (showres)
	    std::cout << "B\n" << Bc << std::endl;
	}
      }
      else
	fillrandom(Vc);
    }
    else {
      fillrandom(Ar);
      if (showres)
	std::cout << "A\n" << Ar << std::endl;
      if (isbinary) {
	if (opdesc.self)
	  Br=Ar;
	else {
	  fillrandom(Br);
	  if (showres)
	    std::cout << "B\n" << Br << std::endl;
	}
      }
      else
	fillrandom(Vr);
    }

    if (opdesc.compbaseline) {
      switch (op) {
 //      case TRANSPOSE: 
// 	if (iscomplex)
// 	  _transpose_naive(Dc,Ac);
// 	else
// 	  _transpose_naive(Dr,Ar);
// 	break;
      case ADD:
	if (opdesc.iscomplex)
	  add_baseline(Dc,Ac,Bc);
	else
	  add_baseline(Dr,Ar,Br);
	break;
      case MULTIPLY:
	if (opdesc.iscomplex)
	  multiply_baseline(Dc,Ac,Bc);
	else
	  multiply_baseline(Dr,Ar,Br);
	break;
      case CTM:
	conj_transpose_multiply_baseline(Dc,Ac,Bc);
	break;
      case MCT:
	multiply_conj_transpose_baseline(Dc,Ac,Bc);
	break;
//       case INDEX:
// 	if (iscomplex)
// 	  fill_baseline(Dc);
// 	else
// 	  fill_baseline(Dr);
// 	break;
      default: throw InternalError("Invalid operation");
      }
    }    

    double fops;

    switch (op) {
      //case TRANSPOSE:
    case ADD: case COPY: fops=(opdesc.iscomplex ? 2 : 1)*m*double(n); break;
    case MVMUL: fops=(opdesc.iscomplex ? 8 : 2)*m*double(n); break;
      //NB mla uses real scalar, so fops=4 rather than 8
    case MLA: fops=(opdesc.iscomplex ? 4 : 2)*m*double(n); break;
    case MULTIPLY: case MCT: case CTM: fops=(opdesc.iscomplex ? 8 : 2)*m*n*double(n); break;
    default: throw InternalError("Invalid operation");
    }

    size_t times=int(1000000*(double(nmflops)/fops));
    if (times<1)
      times=1;

    std::cout << n << "x" << m << "\t" << times << "\t" << fops; std::cout.flush();

    if (opdesc.iscomplex)
      Cc=complex(0.0);
    else
      Cr=0.0;

#ifdef LCM_USE_EXTERNAL
    bool useint;
#endif
    double ires=0;
    double eres=0;
    const size_t nmethods=strlen(opdesc.methods);
    
    for (size_t ntry=0;ntry<nmethods;ntry++) {
      const char meth=opdesc.methods[ntry];

      bool foundop=true;
      timer<CPUTimer> stopwatch;
      switch (op) {
//       case TRANSPOSE:
// 	switch (meth) {
// 	case 'B':
// 	  if (opdesc.iscomplex)
// 	    for (size_t loop=times;loop--;)
// 	      _transpose_naive(Cc,Ac);
// 	  else
// 	    for (size_t loop=times;loop--;)
// 	      _transpose_naive(Cr,Ar);
// 	  break;
// 	case 'U': case 'I':
// 	  if (opdesc.iscomplex)
// 	    for (size_t loop=times;loop--;)
// 	      transpose(Cc,Ac);
// 	  else
// 	    for (size_t loop=times;loop--;) 
// 	      transpose(Cr,Ar);
// 	  break;
// 	default:
// 	  foundop=false;
// 	}
// 	break;
	
      case ADD:
	switch (meth) {
	case 'B':
	  if (opdesc.iscomplex)
	    for (size_t loop=times;loop--;)
	      addip_baseline(Cc,Ac); //add_baseline(Cc,Ac,Bc)
	  else
	    for (size_t loop=times;loop--;)
	      addip_baseline(Cr,Ar);//add_baseline(Cr,Ar,Br);
	  break;
	case 'U': case 'I':
	  if (opdesc.iscomplex)
	    for (size_t loop=times;loop--;)
	      Cc+=Ac;//add(Cc,Ac,Bc);
	  else
	    for (size_t loop=times;loop--;)
	      Cr+=Ar;//add(Cr,Ar,Br);
	  break;
#ifdef LCM_USE_EXTERNAL
	case 'E':
	  if (opdesc.iscomplex) {
	    //	    useint=(Ac.size()<LCM_INTERNAL_ZMLA);
	    for (size_t loop=times;loop--;)
	      lapack_mla(Cc.row(),complex(1.0),Ac.row());
	  }
	  else {
	    //	    useint=(Ar.size()<LCM_INTERNAL_DMLA);
	    for (size_t loop=times;loop--;)
	      lapack_mla(Cr.row(),1.0,Ar.row());
	  }
	  break;
#endif
	default:
	  foundop=false;
	}
	break;
      
      case COPY:
	switch (meth) {
	case 'U': case 'I':
	  if (opdesc.iscomplex)
	    for (size_t loop=times;loop--;)
	      Cc=Ac;
	  else
	    for (size_t loop=times;loop--;)
	      Cr=Ar;
	  break;
#ifdef LCM_USE_EXTERNAL
	case 'E':
	  if (opdesc.iscomplex) {
	    //	    useint=(Ac.size()<LCM_INTERNAL_ZMLA);
	    for (size_t loop=times;loop--;)
	      lapack_copy(Cc.row(),Ac.row());
	  }
	  else {
	    //	    useint=(Ar.size()<LCM_INTERNAL_DMLA);
	    for (size_t loop=times;loop--;)
	      lapack_copy(Cr.row(),Ar.row());
	  }
	  break;
#endif
	default:
	  foundop=false;
	}
	break;
      
      case MULTIPLY:
	switch (meth) {
	case 'B':
	  if (opdesc.iscomplex)
	    for (size_t loop=times;loop--;)
	      multiply_baseline(Cc,Ac,Bc);
	  else
	    for (size_t loop=times;loop--;) 
	      multiply_baseline(Cr,Ar,Br);
	  break;
	case 'I':
	  if (opdesc.iscomplex)
	    for (size_t loop=times;loop--;) {
	      RawMatrix<complex> d(Cc);
	      cmatrix_multiply(d,Ac,Bc,false);
	    }
	  else
	    for (size_t loop=times;loop--;) 
	      cmatrix_multiply(Cr,Ar,Br,false);
	  break;
	case 'U':
	  if (opdesc.iscomplex)
	    for (size_t loop=times;loop--;)
	      multiply(Cc,Ac,Bc);
	  else
	    for (size_t loop=times;loop--;) 
	      multiply(Cr,Ar,Br);
	  break;
#ifdef LCM_USE_EXTERNAL
	case 'E':
	  if (opdesc.iscomplex) {
	    xover_val=LCM_INTERNAL_ZMM;
	    useint=issimple(Ac,Bc,xover_val);
	    if (issym) {
	      for (size_t loop=times;loop--;) {
		RawMatrix<complex> Craw(Cc);
		if (opdesc.self)
		  lapack_hermitian_square(Craw,RawMatrix<complex>(Ac),false);
		else
		  lapack_hermitian_multiply(Craw,RawMatrix<complex>(Ac),RawMatrix<complex>(Bc),false);
	      }
	    }
	    else {
	      for (size_t loop=times;loop--;) {
		RawMatrix<complex> Craw(Cc);
		lapack_multiply(Craw,RawMatrix<complex>(Ac),RawMatrix<complex>(Bc),false);
	      }
	    }
	  }
	  else {
	    xover_val=LCM_INTERNAL_DMM;
	    useint=issimple(Ar,Br,xover_val);
	    if (issym) {
	      for (size_t loop=times;loop--;) {
		RawMatrix<double> Craw(Cr);
		if (opdesc.self)
		  lapack_hermitian_square(Craw,RawMatrix<double>(Ar),false);
		else
		  lapack_hermitian_multiply(Craw,RawMatrix<double>(Ar),RawMatrix<double>(Br),false);
	      }
	    }
	    else {
	      for (size_t loop=times;loop--;) {
		RawMatrix<double> Craw(Cr);
		lapack_multiply(Craw,RawMatrix<double>(Ar),RawMatrix<double>(Br),false);
	      }
	    }
	  }
	  break;
#endif
	default:
	  foundop=false;
	}
	break;

    case MVMUL:
      switch (meth) {
      case 'I': case 'U':
	if (opdesc.iscomplex)
	  for (size_t loop=times;loop--;)
	    multiply(Wc,Ac,Vc);
	else
	  for (size_t loop=times;loop--;) { 
	    Wr.create(Ar.rows()); 
	    cmatrix_multiply(Wr,Ar,Vr);
	  }
	break;
#ifdef USE_ATLAS
      case 'E':
	if (opdesc.iscomplex) {
	  xover_val=LCM_INTERNAL_ZMV;
	  useint=issimple(Ac,Bc,xover_val);
	  for (size_t loop=times;loop--;)
	    lapack_mvmul(Wc,Ac,Vc);
	}
	else {
	  xover_val=LCM_INTERNAL_DMV;
	  useint=issimple(Ar,Br,xover_val);
	  for (size_t loop=times;loop--;)
	    lapack_multiply(Wr,Ar,Vr);
	}
	break;
#endif
      default:
	foundop=false;
      }
      break;
      
      case CTM:
	switch (meth) {
	case 'B':
	  for (size_t loop=times;loop--;)
	    conj_transpose_multiply_baseline(Cc,Ac,Bc);
	  break;
	case 'I':
	  for (size_t loop=times;loop--;)
	    cmatrix_CTM(Cc,Ac,Bc);
	  break;
	case 'U':
	  for (size_t loop=times;loop--;)
	    conj_transpose_multiply(Cc,Ac,Bc);
	  break;
#ifdef LCM_USE_EXTERNAL
	case 'E':
	  xover_val=LCM_INTERNAL_ZTM;
	  useint=issimple(Ac,Bc,xover_val);
	  for (size_t loop=times;loop--;)
	    lapack_CTM(Cc,Ac,Bc);
	  break;
#endif
	default:
	  foundop=false;
	}
	break;

      case MCT:
	switch (meth) {
	case 'B':
	  for (size_t loop=times;loop--;) 
	    multiply_conj_transpose_baseline(Cc,Ac,Bc);
	  break;
	case 'I':
	  for (size_t loop=times;loop--;) 
	    cmatrix_MCT(Cc,Ac,Bc);
	  break;
	case 'U':
	  for (size_t loop=times;loop--;) 
	    multiply_conj_transpose(Cc,Ac,Bc);
	  break;
#ifdef LCM_USE_EXTERNAL
	case 'E':
	  xover_val=LCM_INTERNAL_ZMT;
	  useint=issimple(Ac,Bc,xover_val);
	  for (size_t loop=times;loop--;) 
	    lapack_MCT(Cc,Ac,Bc);
	  break;
#endif
	default:
	  foundop=false;
	}
	break;

//       case INDEX:
// 	switch (meth) {
// 	case 'B':
// 	  if (opdesc.iscomplex)
// 	    for (size_t loop=times;loop--;) fill_baseline(Dc);
// 	  else
// 	    for (size_t loop=times;loop--;) fill_baseline(Dr);
// 	  break;
// 	case 'U': case 'I':
// 	  if (opdesc.iscomplex)
// 	    for (size_t loop=times;loop--;) fill(Dc);
// 	  else
// 	    for (size_t loop=times;loop--;) fill(Dr);
// 	  break;
// 	default:
// 	  foundop=false;
// 	}
// 	break;

      case MLA:
	switch (meth) {
	case 'I': case 'B': case 'U':
	  if (opdesc.iscomplex) {
	    BaseList<complex> Acr=Ac.row(); 
	    for (size_t loop=times;loop--;)
	      mla(Acr,2.0,Bc.row());
	  }
	  else {
	    BaseList<double> Arr=Ar.row();
	    for (size_t loop=times;loop--;)
	      mla(Arr,2.0,Br.row());
	  }
	  break;
#ifdef USE_ATLAS
	case 'E':
	  if (opdesc.iscomplex) {
	    const complex two(2.0);
	    xover_val=LCM_INTERNAL_ZMLA;
	    useint=(Ac.size()<xover_val);
	    for (size_t loop=times;loop--;)
	      lapack_mla(Ac.row(),two,Bc.row());
	  }
	  else {
	    xover_val=LCM_INTERNAL_DMLA;
	    useint=(Ar.size()<xover_val);
	    for (size_t loop=times;loop--;)
	      lapack_mla(Ar.row(),2.0,Br.row());
	  }
	  break;
#endif
	default:
	  foundop=false;
	}
	break;
      }

      if (foundop) {
	double t=stopwatch()/times;

	if (op==MLA) {
	  if (opdesc.iscomplex)
	    Cc=Ac;
	  else
	    Cr=Ar;
	}
	
	const double res=fops/(t*1e6);
	std::cout << "\t" << res;
	if (opdesc.compbaseline)
	  std::cout << "\t" << (opdesc.iscomplex ? norm(Dc-Cc) : norm(Dr-Cr));
	
	if (showres) {
	  if (opdesc.iscomplex)
	    std::cout << "\nOutput\n" << Cc << std::endl;
	  else
	    std::cout << "\nOutput\n" << Cr << std::endl;
	}
	else
	  std::cout.flush();

	switch (meth) {
	case 'I': ires=res; break;
	case 'E': eres=res; break;
	}
      }
      else {
	std::cout << "\tN/A";
	std::cout.flush();
      }
    }
#ifdef LCM_USE_EXTERNAL
    if (eres)
      std::cout << (useint ? "\town" : "\tExternal");
#endif
    std::cout << '\n';

    if (eres>ires && ires && (xover==0))
      xover=n;

    if (isexp)
      n*=factor;
    else
      n+=factor;
  }
  if (xover==start)
    xover=-1;

  return std::pair<int,int>(xover,xover_val);
}

void find_xover(FILE* fp, optimise_t opdesc)
{
  const int huge=20000;
  int xover=0;
  for (;;) {
    std::pair<int,int> res=doopt(opdesc,2,256,2,true);
    xover=res.first;
    if (xover<=0) {
      std::cerr << "Failed to find crossover for " << opdesc.opname << ". Using dummy value\n";
      xover=(xover==0) ? huge : 0;
      break;
    }
    else {
      int step=(xover>=32) ? 2 : 1;
      res=doopt(opdesc,xover/2,xover,step,false);
      xover=res.first;
      if (xover<=0)
	std::cout << "Inconsistent crossover for " << opdesc.opname << ". Repeating...\n";
      else
	break;
    }
  }
  if ((opdesc.op=='L') && (xover<huge))
    xover*=xover; //!< MLA uses number of matrix elements
  
  std::cout << "Crossover for " << opdesc.opname << ": " << xover << "\n\n";

  fprintf(fp,"#define LCM_INTERNAL_%s %i\n",opdesc.opname,xover);
}


int main(int argc,const char *argv[])
{
  try {

  int count=1;
  nmflops=getint(argc,argv,count,"Mflops? ",500);

  std::cout << "A - Addition\nM - multiply\nC - conj_transpose_multiply\nX - multiply_conj_transpose\nV - matrix/vector multiply\nL - mla\nY - copy\nO - optimise crossover\n";
  const op_t op=(op_t)getoption(argc,argv,count,"Operation? ","AMCXVLYO");

  if (op==XOVER) {
#ifndef LCM_USE_EXTERNAL
    std::cerr << "Can't determine cross-over as no external library\n";
    return 1;
#endif
    char buffer[240];
    char fname[256];
    getstring(argc,argv,count,"System type for file? ",buffer,sizeof(buffer),"default");
    sprintf(fname,"%s_xover.h",buffer);
    FILE* fp=fopen(fname,"wa");
    fprintf(fp,"/* Internal / external cross-overs created by testops %i */\n",nmflops);
    find_xover(fp,optimise_t(MULTIPLY, T_REAL));
    find_xover(fp,optimise_t(MULTIPLY, T_COMPLEX));
    find_xover(fp,optimise_t(CTM));
    find_xover(fp,optimise_t(MCT));
    find_xover(fp,optimise_t(MLA,T_REAL));
    find_xover(fp,optimise_t(MLA,T_COMPLEX));
    find_xover(fp,optimise_t(MVMUL,T_REAL));
    find_xover(fp,optimise_t(MVMUL,T_COMPLEX));
    fclose(fp);
    return 0;
  }

  std::cout << "R - general real\nS - real symmetric\nC - general complex\nH - Hermitian\n";
  const mat_t mtype=(mat_t)getoption(argc,argv,count,"Matrix type? ","RSCH");
  //  const bool iscomplex = complexity(op,mtype);
  issym= (mtype==T_SYMMETRIC) || (mtype==T_HERMITIAN);
  optimise_t opdesc(op,mtype);
  opdesc.self=issym ? getlogical(argc,argv,count,"A=B? ",false) : false;

#ifdef LCM_USE_EXTERNAL
  const size_t methods=4;
  std::cout << "I - Internal method\nE - External (LAPACK/BLAS) method\nU - method Used\nB - Baseline method\n";
#else
  const size_t methods=2;
  std::cout << "U - method Used\nB - Baseline method\n";
#endif
  
  getstring(argc,argv,count,"Methods to time? ",opdesc.methods,methods+1,"U");
  opdesc.compbaseline = getlogical(argc,argv,count,"Compare result with baseline algorithm? ",false);
  opdesc.issquare = issym ? true : getlogical(argc,argv,count,"Square matrices? ",true);  
  const bool isexp = getlogical(argc,argv,count,"Multiplicative step factor? ");
  const size_t minstep=isexp ? 2 : 1;
  const size_t start= getint(argc,argv,count,"Start? ",2);
  const size_t end= getint(argc,argv,count,"End? ",512);
  const size_t factor = getint(argc,argv,count,"Step factor? ",minstep);
  if (factor<minstep) {
    std::cerr << "Step factor cannot be less than " << minstep << '\n';
    exit(1);
  }
  std::pair<int,int> res=doopt(opdesc,start,end,factor,isexp);
  int xover=res.first;

  if (xover) {
    if (xover>start) {
      if (op=='L')
	xover*=xover; //!< MLA uses number of matrix elements
      std::cout << "From these results\n#define LCM_INTERNAL_" << opdesc.opname << " " << xover << "\nin config.h\n";
    }
  }
  if (res.second)
    std::cout << "Current value of LCM_INTERNAL_"  << opdesc.opname << " is " << res.second << '\n';

  } catch (MatrixException& exc) {
    std::cerr << exc << std::endl;
    return 1;
  }

  return 0;
}
