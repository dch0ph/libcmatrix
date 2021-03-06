#ifndef List_h_
#define List_h_

#include <cassert>
#include "BaseList.h"
#include "lcm_basethreads.h"

namespace libcmatrix {

//   template<class T> struct list_traits {
//     typedef DefaultAllocator<> allocator;
//   };

template<typename T> class List : public BaseList<T> {
public:
  List(mxflag::tempflag =mxflag::normal);
  explicit List(int);
  
  using BaseList<T>::nitems;
  using BaseList<T>::datast;

  //! this confuses VC++ which doesn't implement ANSI standard lookup ...
  template<class T2> explicit List(const T2& a, mxflag::tempflag tflag =mxflag::normal)
  { construct_(Int2Type<LCM_DIM(T2)>(),a,tflag); }

  List(const List<T>&, mxflag::tempflag =mxflag::normal);
  List(const BaseList<T>&, mxflag::tempflag =mxflag::normal);

  List(int,const T&, mxflag::tempflag =mxflag::normal);
  template<class T2> List(const T2&, const T2&);
  List(int,T*, mxflag::tempflag =mxflag::normal);

  ~List();

  bool istemporary() const
  { return (store_type==mxflag::temporary); }
  
  bool isdynamic() const
  { return (store_type!=mxflag::nondynamic); }
  
  bool hasclaim() const 
  { return (this->datast!=NULL) && (store_type!=mxflag::nondynamic); }

  void ensureempty() { 
    if (this->size())
      *this=T(0);
  }

  void kill();

  void create(int);
  void create(int,const T&);
  template<class F> void create(int n, const F& iter) {
    BaseList<T>::create(n,n ? alloc_t::get_iter(n,n,iter) : NULL);
    allocitems=n;
  }

  template<class T2> void set_dimensions(const T2& a) {
    LCM_STATIC_CHECK( LCM_DIM(T2)==1 , set_dimensions_from_non_1D_object );
    List<T>::create(a.size());
  }  
    
  template<class T2> List<T>& operator=(const T2& a) {
    assign(*this,a);
    return *this;
  }

  inline List<T>& operator=(const BaseList<T>&);
  inline List<T>& operator=(const List<T>& a) {
    if (this!=&a)
      operator=(static_cast< const BaseList<T>& >(a));
    return *this;
  }

  List<T>& operator=(T* a) { BaseList<T>::operator=(a); return *this; }

  template<typename T2> void push_back(const T2& val) { push_back_(val,Int2Type<LCM_DIM(T2)>()); }
  void push_back(size_t n, const T&);
  inline void push_back(const T& val =T()) { push_back_(val,Int2Type<0>()); }
  inline void pop_back();

#ifdef LCM_USE_CXX11
  template<typename... Args_> void emplace_back(Args_&&...);
#endif

  //back comes from BaseList

  template<typename T2> List<T>& operator+= (const T2& a) { return add_ip(*this,a); }
  template<typename T2> List<T>& operator-= (const T2& a) { return subtract_ip(*this,a); }
  template<typename T2> List<T>& operator*= (const T2& a) { return multiply_ip(*this,a); }
  template<typename T2> List<T>& operator&= (const T2& a) { return premultiply_ip(*this,a); }
  template<typename T2> List<T>& operator/= (const T2& a) { return divide_ip(*this,a); }

  BaseList<T> truncate(size_t new_items) const {
    if (new_items>allocitems)
      throw Failed("truncate: Can't truncate beyond allocated list size");
    return BaseList<T>::truncate(new_items);
  }
  void clear() { kill(); }
  void resize(int);
  size_t capacity() const { return allocitems; }
  void reserve(size_t);
  void swap(List<T>&);
  T* erase(T*); //delete pointed to element.  Note that element ordering is *not* preserved!

 private:
  size_t allocitems;
  mxflag::tempflag store_type;

  typedef MemoryFiller<T,type_traits<T>::trivialconstructor, DefaultAllocator< memory_traits<T>::alignment > > alloc_t;

  void zap() { 
    allocitems=(this->nitems)=0;
    this->datast=NULL;
  }

  inline void push_back_(const T&, Int2Type<0>);
  template<typename T2> void push_back_(const T2&, Int2Type<1>);

  void replacevector(size_t n, size_t allocitems_,T* newd) {
    if (this->datast) {
      assert(store_type!=mxflag::nondynamic); //this should be caught earlier
	//	alloc_t::release(newd,allocitems_,n);
	//throw Failed("Cannot resize nondynamic list");
      //      }
      alloc_t::release(this->datast,allocitems,this->nitems);
    }
    BaseList<T>::create(n,newd);
    allocitems=allocitems_;
  }

  template<bool HaveMemory> void set(mxflag::tempflag tflag) {
    if (!HaveMemory && (tflag==mxflag::nondynamic))
      throw Failed("DynamicList<T>(): cannot create empty nondynamic object"); 
#ifdef LCM_ALLOW_TEMPORARY
    store_type=tflag;
#else
    store_type=(tflag==mxflag::temporary) ? mxflag::normal : tflag;
#endif
  }

  void construct_(Int2Type<0>, int n, mxflag::tempflag flag) {  
    set<false>(flag);
    BaseList<T>::create(n,n ? alloc_t::get(n,n) : NULL); 
    allocitems=n;
  }
  
  template<class T2> void construct_(Int2Type<1>, const T2& a, mxflag::tempflag flag) {
    set<false>(flag);
    create(a.size(),a.begin());
   }
    
  void rawcreate(const BaseList<T>&, mxflag::tempflag);
};

 template<typename T> inline std::ostream& operator<< (std::ostream& ostr, const List<T>& a) {
    if (a.istemporary()) ostr << "[T] ";
    return ostr << static_cast<const BaseList<T>& >(a);
  }

 template<typename T> struct type_traits< List<T> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef T value_type;
 };

//Rules for new object operations

template<class T,size_t N1,size_t N2> struct new_trait_ {
  typedef void value_type; //shouldn't really be needed!
};
template<class T> struct new_trait_<T,0,0> { typedef T value_type; };
template<class T> struct new_trait_<T,1,1> { typedef List<T> value_type; };
template<class T> struct new_trait_<T,1,0> { typedef List<T> value_type; };
template<class T> struct new_trait_<T,0,1> { typedef List<T> value_type; };
//Can't define matrix OP vector result

template<class T1,class T2,bool autopro> struct new_trait {
  typedef typename new_trait_< LCM_NEWTYPE(LCM_VAL(T1),LCM_VAL(T2),autopro),LCM_DIM(T1),LCM_DIM(T2)>::value_type value_type;
};
#define LCM_NEWOBJECT(T1,T2) typename new_trait<T1,T2,true>::value_type

template<class T> List<T> operator-(const BaseList<T>&);

template<class T> inline List<T> operator+(const T& a,const BaseList<T>& b)
{ return b+a; }

template<class T> inline List<T> operator-(const T& a,const BaseList<T>& b) {
  List<T> d(mxflag::temporary);
  subtract(d,a,b);
  return d;
}

//NB behaviour if ptr is not a valid iterator is undefined
 template<class T> T* List<T>::erase(T* ptr)
   {
     const size_t ind=ptr-BaseList<T>::vector();
     if (ind>=this->size())
       throw BadIndex("erase");
     *ptr=this->back();
     pop_back();
     return ptr;
   }

  // class complex;

 inline List<double> real(const BaseList<complex>& a) {
   List<double> d(mxflag::temporary);
   real(d,a);
   return d; }

 inline List<double> imag(const BaseList<complex>& a) {
   List<double> d(mxflag::temporary);
   imag(d,a);
   return d; }

 inline List<double> enorm(const BaseList<complex>& a) {
   List<double> d(mxflag::temporary);
   enorm(d,a);
   return d; }

 inline List<double> enorm(const BaseList<double>& a) {
   List<double> d(mxflag::temporary);
   enorm(d,a);
   return d; }

 template<class T> inline List<T> conj(const BaseList<T>& a) {
   List<T> d(mxflag::temporary);
   conj(d,a);
   return d; }

template<typename T,typename T2> struct doespremultiply_ip< List<T>,T2> {
  void operator()(List<T>& d, const T2& a) const
  { premultiply_ip(d,a); }
};

  //implementation details

  template<class T> List<T>::List(mxflag::tempflag flag)
    : BaseList<T>(), allocitems(0)
  { set<false>(flag); }

  template<class T> template<class T2> List<T>::List(const T2& start_, const T2& end_) 
    : BaseList<T>(), allocitems(0), store_type(mxflag::normal) 
  {
    T2 start(start_);
    while (start!=end_) {
      push_back(*start);
      ++start;
    }
  }

  template<class T> List<T>::List(int n)
  { construct_(Int2Type<0>(),n,mxflag::normal); }

 template<class T> List<T>::List(int n,const T& val, mxflag::tempflag tflag)
   {
     set<false>(tflag);
     BaseList<T>::create(n,n ? alloc_t::get(n,n,val) : NULL);
     allocitems=n;
   }

  template<class T> void List<T>::rawcreate(const BaseList<T>& a, mxflag::tempflag tflag)
  {
    const size_t n=a.size();
    set<false>(tflag);
    if (isdynamic())
      BaseList<T>::create(n,n ? alloc_t::get(n,n,a.vector()) : NULL);
    else
      BaseList<T>::create(n,a.vector());
    allocitems=n;
  }

   template<class T> List<T>::List(int n, T* valp, mxflag::tempflag tflag)
     {
       rawcreate(BaseList<T>(n,valp),tflag);
     }

  template<class T> List<T>::List(const BaseList<T>& a, mxflag::tempflag tflag)
  {
    rawcreate(a,tflag);
  }

     //NB We have allocated *allocitems* entries
     template<class T> List<T>::~List() { 
       if (hasclaim())
	 alloc_t::release(this->datast,allocitems,this->nitems); //only nitems items have been created
     }
     
     template<class T> void List<T>::kill() { 
       if (store_type==mxflag::nondynamic) 
	 throw Failed("Can't explicitly kill non-dynamic objects"); 
       if (this->datast) {
	 alloc_t::release(this->datast,allocitems,this->nitems);
	 zap();
       }
     }

  template<typename T> List<T>& List<T>::operator=(const BaseList<T>& a)
  {
    if (a.empty())
      clear();
    else {
      const size_t n(a.size());
      if (n>allocitems) {
	if (store_type==mxflag::nondynamic)
	  throw Failed("Cannot resize nondynamic List");
	replacevector(n,n,alloc_t::get(n,n,a.vector()));
	allocitems=n;
      }
      else { //NB Not exception safe
	if (n<this->nitems)
	  alloc_t::memop_t::destroy(this->datast+n,this->datast+this->nitems); //kill off excess
	else
	  alloc_t::memop_t::construct(this->datast+this->nitems,this->datast+n,a.vector()+this->nitems);
	this->nitems=n;
	alloc_t::memop_t::assign(this->datast,this->datast+n,a.vector());
      }
    }
    return *this;
  }
    
 //NB memory only claimed for *active* elements of a i.e. may not be exact copy!
  template<class T> List<T>::List(const List<T>& a, mxflag::tempflag tflag)
    :  BaseList<T>()
  {
       if ((a.store_type==mxflag::nondynamic) || (tflag==mxflag::nondynamic))
	 throw Failed("Can't create create dynamic object from non-dynamic or vice versa; use reference!");
       set<true>(tflag);
       if (a.datast==NULL) {
	 allocitems=0;
	 return;
       }
#ifdef LCM_ENABLE_TEMPORARY
       if (a.store_type==mxflag::temporary) {
	 BaseList<T>::create(a.nitems,a.datast);
	 allocitems=a.allocitems;
	 (const_cast< List<T>& >(a)).zap();
       }
       else {
#endif
	 if (store_type==mxflag::nondynamic) {
	   BaseList<T>::create(a.nitems,a.datast);
	   allocitems=a.allocitems;
	 }
	 else {
	   const size_t n=a.nitems;
	   BaseList<T>::create(n,alloc_t::get(n,n,a.datast));
	   allocitems=n;
	 }
#ifdef LCM_ENABLE_TEMPORARY
       }
#endif
     }
 
template<typename T> void List<T>::reserve(size_t n) //change allocitems while preserving contents
{
  if ((store_type==mxflag::nondynamic) && (n!=nitems))
    throw Failed("Cannot resize nondynamic list");
  if (n>allocitems)
    replacevector(nitems,n,nitems ? alloc_t::get(n,nitems,datast) : alloc_t::get_uninit(n));
}

template<typename T> void List<T>::create(int n) 
{
  if ((store_type==mxflag::nondynamic) && (n!=this->nitems))
    throw Failed("Cannot resize nondynamic list");
  if (n>allocitems)
    replacevector(n,n,alloc_t::get(n,n));
  else {
    if (n<0)
      throw InvalidParameter("List size cannot be negative");
  }
  this->nitems=n;
}
   
 template<typename T> void List<T>::resize(int new_items) {
   if (new_items==this->nitems)
     return;
   if (new_items<0)
     throw InvalidParameter("List size cannot be negative");
   if (store_type==mxflag::nondynamic)
     throw Failed("non-dynamic objects can't be resized");
   if (new_items>this->nitems) {
     if (new_items>allocitems)
       reserve(new_items);
     alloc_t::memop_t::construct(this->datast+this->nitems,this->datast+new_items);
   }
   else
     alloc_t::memop_t::destroy(this->datast+new_items,this->datast+this->nitems);
   this->nitems=new_items;
 }
   
template<typename T> void List<T>::create(int n,const T& a)
{ 
  if (n!=this->nitems) {
    if (store_type==mxflag::nondynamic)
      throw Failed("Can't resize dynamic object");

    if (n>allocitems) {
      replacevector(n,n,alloc_t::get(n,n,a));
      return;
    }
    alloc_t::memop_t::destroy(this->datast+n,this->datast+allocitems);
    this->nitems=n;
  }
  alloc_t::memop_t::assign(this->datast,this->datast+n,a);
}

template<typename T> void List<T>::push_back_(const T& a, Int2Type<0>)
  {
    if (nitems==allocitems)
      //extension will fail if memory is static
      reserve((allocitems>8) ? size_t(LCM_EXTEND_STACK*allocitems) : LCM_DEFAULT_STACK);
    alloc_t::memop_t::construct(datast[nitems],a);
    nitems++;
  }

#ifdef LCM_USE_CXX11
 template<typename T> template<typename... Args_> void List<T>::emplace_back(Args_&&... args_)
   {
    if (nitems==allocitems)
      //extension will fail if memory is static
      reserve((allocitems>8) ? size_t(LCM_EXTEND_STACK*allocitems) : LCM_DEFAULT_STACK);
    alloc_t::memop_t::construct(datast[nitems], std::forward<Args_>(args_)...);
    nitems++;
  }

#endif

template<typename T> template<typename T2> void List<T>::push_back_(const T2& a, Int2Type<1>)
{
  const size_t n(a.size());
  const size_t newsize=nitems+n;
  if (newsize>allocitems) {
    const size_t extendto= (allocitems>8) ? size_t(LCM_EXTEND_STACK*allocitems) : LCM_DEFAULT_STACK;
    reserve(newsize<extendto ? extendto : newsize);
  }
  alloc_t::memop_t::construct_iter(datast+nitems,datast+newsize,a.begin());
  nitems+=n;
 }

  template<typename T> void List<T>::push_back(size_t n, const T& v)
  {
    const size_t newsize=nitems+n;
    if (newsize>allocitems) {
      const size_t extendto= (allocitems>8) ? size_t(LCM_EXTEND_STACK*allocitems) : LCM_DEFAULT_STACK;
      reserve(newsize<extendto ? extendto : newsize);
    }
    alloc_t::memop_t::construct(datast+nitems,datast+newsize,v);
    nitems+=n;
  }

template<typename T> void List<T>::pop_back()
  { if (this->nitems==0)
      throw Failed("List<T>::pop_back");
    else {
      this->nitems--;
      alloc_t::memop_t::destroy((this->datast)[this->nitems]);
    }
  }

 template<typename T> void List<T>::swap(List<T>& a)
   {
//      if ( (store_type==mxflag::nondynamic) ^ (a.store_type==mxflag::nondynamic))
//        throw Failed("swap: cannot mix dynamic and nondynamic objects");
     if (store_type!=a.store_type)
       throw Failed("swap: cannot swap objects with different storage types");
     BaseList<T>::swap(a);
     ::std::swap(allocitems,a.allocitems);
   }
 
template<typename Ta,typename Tb,typename Tc> inline void divide(List<Ta>& a,const BaseList<Tb>& b,const BaseList<Tc>& c) {
  a.create(b.size());
  divide( static_cast< BaseList<Ta> >(a),b,c);
}

template<class T> inline List<T> operator- (const BaseList<T>& a)
{
  List<T> d(mxflag::temporary);
  ::libcmatrix::negate(d,a);
  return d;
}

template<class T1,class T2> inline typename new_trait<BaseList<T1>,T2,true>::value_type
operator+ (const BaseList<T1>& a,const T2& b)
 { 
   typename new_trait<BaseList<T1>,T2,true>::value_type d(mxflag::temporary);
   add(d,a,b);
   return d;
 }

template<class T1,class T2> inline typename new_trait<BaseList<T1>,T2,true>::value_type
operator- (const BaseList<T1>& a,const T2& b)
 { 
   typename new_trait<BaseList<T1>,T2,true>::value_type d(mxflag::temporary);
   subtract(d,a,b);
   return d;
 }

template<class T1,class T2> inline typename new_trait<BaseList<T1>,BaseList<T2>,true>::value_type
operator* (const BaseList<T1>& a,const BaseList<T2>& b)
{
  typename new_trait<BaseList<T1>,BaseList<T2>,true>::value_type d(mxflag::temporary);
  emultiply(d,a,b);
  return d;
}

template<class T> inline List<T> operator* (const T& a,const BaseList<T>& b)
{
  List<T> d(mxflag::temporary);
  emultiply(d,a,b);
  return d;
}

//Handle these special cases since general templates above assume BaseList is 1st argument
inline List<complex> operator* (double a,const BaseList<complex>& b)
{
  List<complex> d(mxflag::temporary);
  emultiply(d,a,b);
  return d;
}

inline List<complex> operator* (const complex& a,const BaseList<double>& b)
{
  List<complex> d(mxflag::temporary);
  emultiply(d,a,b);
  return d;
}

template<class T1,class T2> inline typename new_trait<BaseList<T1>,T2,true>::value_type
operator/ (const BaseList<T1>& a,const T2& b)
{
  typename new_trait<BaseList<T1>,T2,true>::value_type d(mxflag::temporary);
  divide(d,a,b);
  return d;
}

template<class M> inline List< typename M::value_type > sort(const M& A) {
  List< typename M::value_type > B(A,mxflag::temporary); 
  sort_ip(B);
  return B; }

template<class M> inline List< typename M::value_type > stable_sort(const M& A) {
  List< typename M::value_type > B(A,mxflag::temporary); 
  stable_sort_ip(B);
  return B; }

template<class M,class F> inline List< typename M::value_type > sort(const M& A,const F& compare) {
  List< typename M::value_type > B(A,mxflag::temporary);
  sort_ip(B,compare);
  return B; }

// idea adapted from http://www.ddj.com/cpp/184401318

template<class random_iterator> class IndexedComparison : public std::binary_function<size_t,size_t,bool>
{
public:
  IndexedComparison(random_iterator begin)
    : begin_(begin) {}
  
  bool operator() (size_t i, size_t j) const
  { return *(begin_+i) < *(begin_+j); }
   
private:
  random_iterator const begin_;
};

template<class T> void indexed_sort(BaseList<size_t> ind, const T& data) 
{
  for (size_t i=ind.size();i--;)
    ind(i)=i;
  std::sort(ind.begin(),ind.end(),
	    IndexedComparison<typename T::const_iterator>(data.begin()));
}

template<class T> List<size_t> indexed_sort(const T& data) 
{
  List<size_t> ind(data.size(),mxflag::temporary);
  for (size_t i=ind.size();i--;)
    ind(i)=i;
  std::sort(ind.begin(),ind.end(),
	    IndexedComparison<typename T::const_iterator>(data.begin()));
  return ind;
}

#ifdef LCM_USE_EXTERNTEMPLATE
LCM_L1_EXTERN template class List<double>;
LCM_L1_EXTERN template class List<complex>;
LCM_L1_EXTERN template class List<size_t>;
#endif

} //namespace libcmatrix

#endif //List_h_

