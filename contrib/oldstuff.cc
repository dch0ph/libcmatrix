/* sort is based on piksrt from Numerical Recipies section 8.1 to sort an array
   by straight insertion. */

template <class T> void sort(T arr[],int n)
{
  int i,j;
  T a;
  
  for (j=1;j<n;j++) {
    a=arr[j];
    i=j-1;
    while (i >= 0 && arr[i] > a) {
      arr[i+1]=arr[i];
      i--;
    }
    arr[i+1]=a;
  }
}


//Heap sort adapted from Numerical Recipes			  
template<class M,class F> List<size_t> sort_index(const M& A,const F& compare)
{
  typedef LCM_VAL(M) T;
  const size_t n=A.size();
  List<size_t> indx(n,mxflag::temporary);

  int i,j;
  for (j=n;j--;) indx(j)=j;
  int l=n/2+1;
  int ir=n;
  const T* Q;
  size_t indxt;

  for (;;) {
    if (l>1) {
      l--;
      indxt=indx(l-1);
      Q=&(A(indxt));
    }
    else {
      indxt=indx(ir-1);
      Q=&(A(indxt));
      indx(ir-1)=indx(0);
      ir--;
      if (ir==1) {
	indx(0)=indxt;
	return indx;
      }
    }
    i=l;
    j=2*l;
    while (j<=ir) {
      if (j<ir) {
	if ( compare( A(indx(j-1)),A(indx(j)))) j++;
      }
      if ( compare(*Q,A(indx(j-1)))) {
	indx(i-1)=indx(j-1);
	i=j;
	j+=j;
      }
      else
	j=ir+1;
    }
    indx(i-1)=indxt;
  }
}

template<class M> List<size_t> sort_index(const M& A) {
  return sort_index(A,std::less<LCM_VAL(M)>()); }

//old MINUIT interface

// extern "C" {
//   typedef void (*P_MINUITLINKFUNC)(int *,double *,double *,double *,int *,const ::libcmatrix::MinuitFunction&);

//   void f_init();
//   void f_exit();
//   void f_setarg();
//   void f_setsig();
//   void minuit_( P_MINUITLINKFUNC,void *);
//   void minuitinitio_();
// }




// struct minuit {
//   enum states {
//     notused=0,
//     initial=1,
//     main=2,
//     final=3
//   };
// };

//  typedef double (*P_MINUITFUNC)(const BaseList<double> &,int);

//  extern "C" {
//   void minuitlink(int *,double *,double *,double *,int *,const BaseMinFunction&);
//  }


//  class MinuitFunction : public BaseMinFunction {
//  public:
//    explicit MinuitFunction(P_MINUITFUNC func_) : func(func_) {
//     if (!func) throw InvalidParameter("MinuitFunction(): NULL function");
//    }
//    //   double operator()(const BaseList<double>& pars, int state) const { return (*func)(pars,state); }
//    double operator()(const BaseList<double>& par) const { return (*func)(pars); }

//  private:
//    P_MINUITFUNC func;
//  };

//  void minuitinit(const char * =NULL);
//  void minuitexit();
//  void minuitcall(const MinuitFunction&);
