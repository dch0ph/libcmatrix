#ifndef _baselist_h_
#define _baselist_h_

#include "basedefs.h"
#include "cmatrix_complex.h"
//full iostream will need to be included at some point
#include <iosfwd>
#include <cstdarg>

namespace libcmatrix {

  class range {
  public:
    range() 
      : start_(0), size_(0) {}

    range(int first,int last)
      : start_(first), size_(last-first+1) {
      if (size_<0)
	throw InvalidParameter("range: negative length!");
    }    

    range(int sizev)
      : start_(0), size_(sizev) {
      if (size_<0)
	throw InvalidParameter("range: negative length!");
    }    

    bool empty() const { return (size_==0); }
    size_t start() const { return start_; }
    size_t size() const { return size_; }
    size_t max() const { return start_+size_-1; }

    size_t operator()(size_t ind) const {
#if (BOUNDS_CHECK)
      if (ind>=size_)
	throw BadIndex("range()",ind,size_);
#endif
      return start_+ind;
    }

  private:
    size_t start_;
    int size_;
  };

  inline size_t max(const range& sl) { return sl.max(); }

  template<> struct type_traits<range> {
    static const bool trivialconstructor=true;
    static const size_t dimensionality=1;
    static const int rank=0;
    typedef size_t value_type;
  };

#ifndef LCM_SUPPRESS_VIEWS
  class slice;

  template<typename T,typename Map> class IndirectList;

  template<typename T, size_t Rank,typename Map> struct LCM_List_Select_ {};
  template<typename T,typename Map> struct LCM_List_Select_<T,0,Map> {
    typedef T& return_type;
    static T& func(BaseList<T>& a,size_t n) {
      return a.operator()(n); }
  };
  template<typename T,typename Map> struct LCM_List_Select_<T,1,Map> {
    typedef IndirectList<T,Map> return_type;
    static return_type func(BaseList<T>& a,const Map& map) {
      return return_type(a,map); }
  };
  template<typename T> struct LCM_List_Select_<T,1,range> {
    typedef BaseList<T> return_type;
    static return_type func(BaseList<T>& a,const range& sel) {
      return return_type(sel.size(),a.vector()+sel.start()); }
  };

  template<typename T,size_t Rank,typename Map> struct const_LCM_List_Select_ {};
  template<typename T,typename Map> struct const_LCM_List_Select_<T,0,Map> {
    typedef const T& return_type;
    static const T& func(const BaseList<T>& a,size_t n) {
      return a.operator()(n); }
  };
  template<typename T,typename Map> struct const_LCM_List_Select_<T,1,Map> {
    typedef const IndirectList<T,Map> return_type;
    static return_type func(const BaseList<T>& a,const Map& map) {
      return return_type(a,map); }
  };
  template<typename T> struct const_LCM_List_Select_<T,1,range> {
    typedef const BaseList<T> return_type;
    static return_type func(const BaseList<T>& a,const range& sel) {
      return return_type(sel.size(),a.vector()+sel.start()); }
  };

#endif

/*   template<class T,class T2> struct assign__ { */
/*     static inline BaseList<T>& func(BaseList<T>& d,const T2& a) { */
/*       assign(d,a); */
/*       return d; } */
/*   }; */
/*   template<class T> struct assign__<T,const T*> { */
/*     static inline BaseList<T>& func(BaseList<T>& d,const T* a) { */
/*       return (d= (const_cast<T*>(a))); } //catch case of assign from const T* */
/*   }; */

  template<typename T> class BaseList : public ::std::unary_function<size_t,T> {
public:
  BaseList() { clear(); }

  BaseList(int n,T* addr) { 
    create(n,addr); }

  // definitions for Standard C++ Library containers
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;  
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef int difference_type;
  typedef size_t size_type;

  template<class T2> inline BaseList<T>& operator=(const T2& a) {
    assign(*this,a);
    return *this;
  }
    //    return assign__<T,T2>::func(*this,a); }

  inline BaseList<T>& operator=(T* addr) {
    memop_t::assign(datast,datast+nitems,addr);
    return *this; }

/*   inline const BaseList<T>& operator=(const T* addr) { */
/*     memop_t::assign(datast,datast+nitems,addr); */
/*     return *this; } */

  BaseList<T>& operator=(const BaseList<T>&);

  void create(int n,T* addr) { 
    if (n<0)
      throw InvalidParameter("BaseList<T>::create"); 
    nitems=n; datast=addr;
  }

  void create(const BaseList& a) { 
    nitems=a.nitems; datast=a.datast;
  }

  void create(size_t n) {
    if (n!=nitems)
      throw Mismatch("Cannot resize BaseList<T>");
  }

  BaseList<T> row() { return *this; }
  const BaseList<T>& row() const { return *this; }

/*   void replace(size_t n,const T& v) { */
/*     ::libcmatrix::construct(operator()(n),v); */
/*   } */

//     void clear() {
//       throw Failed("Can't explicitly clear BaseList");
//     }

#ifndef LCM_SUPPRESS_VIEWS
    
    template<class Map> typename LCM_List_Select_<T,LCM_DIM(Map),Map>::return_type operator() (const Map& map) {
      return LCM_List_Select_<T,LCM_DIM(Map),Map>::func(*this,map); }

    template<class Map> typename const_LCM_List_Select_<T,LCM_DIM(Map),Map>::return_type operator() (const Map& map) const {
    return const_LCM_List_Select_<T,LCM_DIM(Map),Map>::func(*this,map); }
#endif

  inline T& operator()(size_t);
  inline const T& operator()(size_t) const;

  inline T& at(size_t);
  inline const T& at(size_t) const;
  
  size_t length() const { return nitems; }
  T* vector() const { return datast; }
  
  T* data() { return datast; }
  T* data() const { return datast; }
  
  T& operator[](size_t n) { return datast[n]; }
  const T& operator[](size_t n) const { return datast[n]; }

  //  bool operator!() const { return (nitems==0); }
  bool empty() const { return (nitems==0); }
  size_t size() const { return nitems; }
  void swap(BaseList<T>& a) {
    ::std::swap(nitems,a.nitems);
    ::std::swap(datast,a.datast);
  }

    size_t dimension(size_t n) const {
      if (n!=0)
	throw BadIndex("BaseList::dimension");
      return nitems; }
  
  const_iterator begin() const { return datast; }
  const_iterator cbegin() const { return datast; }
  iterator begin() { return datast; }
  const_iterator end() const { return datast+nitems; }
  const_iterator cend() const { return datast+nitems; }
  iterator end() { return datast+nitems; }

  void print() const
    { std::cout << *this << std::endl; }

  BaseList<T> truncate(size_t new_items) const
    { return BaseList<T>(new_items,datast); }

  template<class T2> BaseList<T>& operator+= (const T2& a) { return add_ip(*this,a); }
  template<class T2> BaseList<T>& operator-= (const T2& a) { return subtract_ip(*this,a); }
  template<class T2> BaseList<T>& operator*= (const T2& a) { return multiply_ip(*this,a); }
  template<class T2> BaseList<T>& operator&= (const T2& a) { return premultiply_ip(*this,a); }
  template<class T2> BaseList<T>& operator/= (const T2& a) { return divide_ip(*this,a); }

  inline T& back();
  inline const T& back() const;
  inline T& front() { 
    if (datast) return *datast;
    throw Undefined("BaseList<T>::front");
  }
  inline const T& front() const { 
    if (datast) return *datast;
    throw Undefined("BaseList<T>::front");
  }

  inline void clear() {
    nitems=0; datast=0;
  }

protected:
  typedef memop<T,type_traits<T>::trivialconstructor> memop_t;

  //! fudge to get around unnecessary warnings
  BaseList(const BaseList<T>& a, mxflag::tempflag flag) {
    if (flag==mxflag::temporary) {
      create(a);
      const_cast<BaseList<T>& >(a).clear();
    }
  }

  size_t nitems;
  T* datast;

  void dump() const { ::libcmatrix::print(*this,std::cout); }
};

 template<typename T> struct type_traits< BaseList<T> > {
   static const bool trivialconstructor=false; //although ctor is trivial, it is not equivalent to assign!
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef T value_type;
 };

template<size_t n,class T> class ExplicitList : public BaseList<T> {
  char stat[n*sizeof(T)];
  typedef memop<T,type_traits<T>::trivialconstructor> memop_t;
public:
// changed from const T2& v0 (undefined behaviour) for va_args with references
  template<class T2> ExplicitList(T2 v0, ... )
    : BaseList<T>(n,(T*)(stat)) {
    T* tstat=BaseList<T>::vector();
    new (tstat++) T(v0);
    va_list va;
    va_start(va,v0);
    for (int i=1;i<n;i++)
      new(tstat++) T(va_arg(va,T2));
    va_end(va);
  }
  ExplicitList(const ExplicitList<n,T>& a)
    : BaseList<T>(n,(T*)(stat)) { 
    memop_t::construct(BaseList<T>::datast,BaseList<T>::datast+n,const_cast<T*>(a.datast));
  }
    ~ExplicitList() {
      memop_t::destroy(BaseList<T>::vector(),BaseList<T>::vector()+n);
    }
 };

 template<size_t N,typename T> struct type_traits< ExplicitList<N,T> > {
   static const bool trivialconstructor= type_traits<T>::trivialconstructor;
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef T value_type;
 };

template<typename T,typename T2> struct doespremultiply_ip< BaseList<T>,T2> {
  void operator()(BaseList<T>& d, const T2& a) const { premultiply_ip(d,a); }
};

template<typename T1,typename T2,typename T3> void divide(BaseList<T1>&, const BaseList<T2>&, const BaseList<T3>&);

 template<class T,class T2> inline bool operator== (const BaseList<T>& a, const T2& b) {
   return areequal(a,b); }

 template<class T,class T2> inline bool operator!= (const BaseList<T>& a, const T2& b) {
   return arenotequal(a,b); }

  bool areequal(const BaseList<complex>&, const BaseList<complex>&, double tol);

/*  template<typename T1,typename T2> LCM_NEWTYPE(T1,T2,true) dot(const T1& a,const T2& b) */
/*    { */
/*      ENSURE_N(LCM_DIM(T1),1); */
/*      if (!arematching(a,b)) throw Mismatch("dot"); */
/*      if (!a) throw Undefined("dot"); */
     
/*      typename T1::const_iterator astart=a.begin(); */
/*      typename T1::const_iterator aend=a.end(); */
/*      typename T2::const_iterator bstart=b.begin(); */
     
/*      LCM_NEWTYPE(T1,T2,true) res=(*astart)*(*bstart); */
/*      while (++astart!=aend) mla(res,*astart,*++bstart); */
/*      return res; */
/*    } */
     
 template<typename T> std::ostream &operator << (std::ostream& ostr, const BaseList<T>& a) {
   print(a,ostr);
   return ostr;
 }

bool isvalid_indexlist(const BaseList<size_t> &,size_t);

template<typename T1,typename T2,typename T3> struct doesmla< BaseList<T1>, T2, BaseList<T3> > : public binary_function_ip<BaseList<T1>,T2,BaseList<T3> > {
  inline void operator()(BaseList<T1>& d, const T2& a, const BaseList<T3>& b) const { 
    mla(d,a,b); }
};

void real_mla(BaseList<double>, complex, const BaseList<complex>&);
void real_mla(BaseList<complex>, complex, const BaseList<complex>&);
 
// Implementation details below

 template<class T> T& BaseList<T>::at(size_t i) { 
    if (i>=nitems)
      throw BadIndex("BaseList<T>::at",i,nitems);
    else
      return datast[i]; }

 template<class T> const T& BaseList<T>::at(size_t i) const {
    if (i>=nitems)
      throw BadIndex("BaseList<T>::at",i,nitems);
    else
      return datast[i]; }

#if (BOUNDS_CHECK)
 template<class T> T& BaseList<T>::operator()(size_t i) { 
    if (i>=nitems)
      throw BadIndex("BaseList<T>()",i,nitems);
    else
      return datast[i]; }

 template<class T> const T& BaseList<T>::operator()(size_t i) const {
    if (i>=nitems)
      throw BadIndex("BaseList<T>()",i,nitems);
    else
      return datast[i]; }

#else
 template<class T> T& BaseList<T>::operator()(size_t i) {
   return datast[i]; }

 template<class T> const T& BaseList<T>::operator()(size_t i) const {
    return datast[i];}
#endif

template<typename T> T& BaseList<T>::back() {
  if (nitems)
    return datast[nitems-1];
  throw Failed("BaseList<T>::back");
}

template<typename T> const T& BaseList<T>::back() const {
  if (nitems)
    return datast[nitems-1];
  throw Failed("BaseList<T>::back");
}

template<class T> BaseList<T>& BaseList<T>::operator= (const BaseList<T>& a)
{
  if (this!=&a) {
    if (a.length()!=nitems)
      throw Mismatch("BaseList::=",a.length(),nitems);
    //    T* avector=const_cast<T*>(a.vector());
    memop_t::assign(datast,datast+nitems,a.vector());
  }
  return *this;
}

template<typename T1,typename T2,typename T3> void divide(BaseList<T1> d,const BaseList<T2>& a, const BaseList<T3>& b)
{
   apply_sd(d,doesdivide_sd<T1,T2,T3>(),a,b);
}

template<class M> void sort_ip(M& a) { ::std::sort(a.begin(),a.end()); }
template<class M,class F> void sort_ip(M& a,const F& compare) { ::std::sort(a.begin(),a.end(),compare); }
template<class M> void stable_sort_ip(M& a) { ::std::stable_sort(a.begin(),a.end()); }
template<class M,class F> void stable_sort_ip(M& a,const F& compare) { ::std::stable_sort(a.begin(),a.end(),compare); }

template<class T> inline void negate_ip(BaseList<T> a) { apply_ip(doesnegate_ip<LCM_VAL(T)>(),a); }

template<class T1,class T2> LCM_NEWTYPE(T1,T2,true) trace_multiply(const BaseList<T1>& b,const BaseList<T2>& c)
{
  size_t n=b.length();
  if (n!=c.length())
    throw Mismatch("trace_multiply");
  LCM_NEWTYPE(T1,T2,true) sum(0);
  for (;n--;)
    mla(sum,b(n),c(n));
  return sum;
}

template<class T> complex dot(const BaseList<complex>& b,const BaseList<T>& c)
{
  size_t n=b.length();
  if (n!=c.length())
    throw Mismatch("dot");
  complex sum(0.0);
  for (;n--;)
    conj_mla(sum,b(n),c(n));
  return sum;
}

template<class T> bool areequal(const BaseList<T>& a, const BaseList<T>& b, double tol) 
{
  size_t n=a.size();
  if (n!=b.size())
    return false;
  for (;n--;) {
    if (!areequal(a(n),b(n),tol))
      return false;
  }
  return true;
}

//complex trace_multiply_conj(const BaseList<complex>&, const BaseList<complex>&);

#ifdef LCM_USE_EXTERNTEMPLATE
#ifdef LCM_LEVEL1_INSTANT
#define LCM_L1_EXTERN
#else
#define LCM_L1_EXTERN extern
#endif
LCM_L1_EXTERN template class BaseList<double>;
LCM_L1_EXTERN template class BaseList<complex>;
LCM_L1_EXTERN template class BaseList<size_t>;
#endif

} // namespace libcmatrix

#ifndef LCM_SUPPRESS_VIEWS
#include "IndirectList.h"
#endif

#endif
