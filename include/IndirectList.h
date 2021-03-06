#ifndef _IndirectList_h_
#define _IndirectList_h_

#include <cstdlib>

namespace libcmatrix {

  template<typename T,typename Map =BaseList<size_t> > class IndirectList {
    
  public:
    IndirectList(const BaseList<T>& datast, const Map& index, int rfac)
      : datast_(datast), index_(index), length_(index.size()) { 
      if (length_==0)
	throw InvalidParameter("IndirectList cannot be null");
      if (rfac<1)
	throw InvalidParameter("IndirectList scaling cannot be 0 or negative\n");
      if (rfac!=1)
	index_*=rfac;
    }

    IndirectList(const BaseList<T>& datast, const Map& index)
      : datast_(datast), index_(index), length_(index.size()) { 
      if (length_==0)
	throw InvalidParameter("IndirectList cannot be null");
    }
    
    bool operator!() const { return false; }
    void create(size_t r) { if (r!=length_)
	throw Mismatch("IndirectList::create"); }
    T* vector() const { return datast_.vector(); }
    void clear() { throw Failed("IndirectList can't be explicitly cleared"); }

    template<class T2> LCM_INLINE IndirectList<T,Map>& operator= (const T2&);
    template<class T2> LCM_INLINE IndirectList<T,Map>& operator+= (const T2&);
    template<class T2> LCM_INLINE IndirectList<T,Map>& operator-= (const T2&);
    template<class T2> LCM_INLINE IndirectList<T,Map>& operator*= (const T2&);
    template<class T2> LCM_INLINE IndirectList<T,Map>& operator/= (const T2&);
    
   class base_iterator_ {
   public:
     base_iterator_(const base_iterator_& x)
       : datast_(x.datast_), cindex_(x.cindex_), index_(x.index_) {}
     base_iterator_(const BaseList<T>& datast,const Map& index,size_t cindex)
       : datast_(datast), cindex_(cindex), index_(index) {}
     bool operator== (const base_iterator_& x) const { return (x.cindex_==cindex_); }
     bool operator!= (const base_iterator_& x) const { return (x.cindex_!=cindex_); }
     bool operator> (const base_iterator_& x) const { return (cindex_>x.cindex_); }
     bool operator< (const base_iterator_& x) const { return (cindex_<x.cindex_); }
     bool operator>= (const base_iterator_& x) const { return (cindex_>=x.cindex_); }
     bool operator<= (const base_iterator_& x) const { return (cindex_<=x.cindex_); }
     ptrdiff_t operator- (const base_iterator_& x) const { return (cindex_-x.cindex_); }

   protected:
     BaseList<T> datast_;
     size_t cindex_;
     const Map& index_; //presence of references means operator= is invalid, copy construction OK (which is reasonable)
   };
   friend class base_iterator_;

   struct iterator : public base_iterator_, public ::std::iterator< ::std::random_access_iterator_tag,T> {
     using base_iterator_::datast_;
     using base_iterator_::index_;
     iterator(const BaseList<T>& datast,const Map& index,size_t cindex) : base_iterator_(datast,index,cindex) {}
     T& operator*() { return datast_((this->index_)(this->cindex_)); }
     iterator& operator++() { this->cindex_++; return *this; }
     iterator operator++(int) { this->cindex_++; return iterator(datast_,index_,this->cindex_-1); }
     iterator& operator--() { this->cindex_--; return *this; }
     iterator operator--(int) { this->cindex_--; return iterator(datast_,index_,this->cindex_+1); }
     iterator& operator+=(ptrdiff_t i) { this->cindex_+=i; return *this; }
     iterator& operator-=(ptrdiff_t i) { this->cindex_-=i; return *this; }
     iterator operator+ (ptrdiff_t i) const { return iterator(datast_,index_,this->cindex_+i); }
     iterator operator- (ptrdiff_t i) const { return iterator(datast_,index_,this->cindex_-i); }
     ptrdiff_t operator- (const iterator& x) const { return (int)(this->cindex_)-(int)x.cindex_; }
     T& operator[] (ptrdiff_t i) { return datast_(index_(this->cindex_+i)); }
   };

   struct const_iterator : public base_iterator_, public ::std::iterator< ::std::random_access_iterator_tag,T,ptrdiff_t,const T*,const T&> {
     using base_iterator_::datast_;
     using base_iterator_::index_;
     const_iterator(const BaseList<T>& datast,const Map& index,size_t cindex) : base_iterator_(datast,index,cindex) {}
     const T& operator*() const { return datast_((this->index_)(this->cindex_)); }
     const_iterator& operator++() { this->cindex_++; return *this; }
     const_iterator operator++(int) { this->cindex_++; return const_iterator(datast_,index_,this->cindex_-1); }
     const_iterator& operator--() { this->cindex_--; return *this; }
     const_iterator operator--(int) { this->cindex_--; return const_iterator(datast_,index_,this->cindex_+1); }
     const_iterator& operator+=(ptrdiff_t i) { this->cindex_+=i; return *this; }
     const_iterator& operator-=(ptrdiff_t i) { this->cindex_-=i; return *this; }
     const_iterator operator+ (ptrdiff_t i) const { return const_iterator(datast_,index_,this->cindex_+i); }
     const_iterator operator- (ptrdiff_t i) const { return const_iterator(datast_,index_,this->cindex_-i); }
     ptrdiff_t operator- (const const_iterator& x) const { return (int)(this->cindex_)-(int)x.cindex_; }
     const T& operator[] (ptrdiff_t i) const { return datast_(index_(this->cindex_+i)); }
   };

   iterator begin() { return iterator(datast_,index_,0); }
   iterator end() { return iterator(datast_,index_,length_); }
   const_iterator begin() const { return const_iterator(datast_,index_,0); }
   const_iterator end() const { return const_iterator(datast_,index_,length_); }
   
   friend class BaseList<T>;

  static bool empty() { return false; } //Not allowed empty IndirectList!
   size_t size() const { return length_; }
   size_t dimension(size_t n) const {
#ifndef NDEBUG
     if (n!=0)
       throw BadIndex("IndirectList::dimension");
#endif
     return size(); }

  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;  

   T& operator() (size_t r) {
#if (BOUNDS_CHECK)
     if (r>=length_)
       throw BadIndex("IndirectList()",r,length_);
#endif
     return datast_(index_(r));
   }
   const T& operator() (size_t r) const {
#if (BOUNDS_CHECK)
     if (r>=length_)
       throw BadIndex("IndirectList()",r,length_);
#endif
     return datast_(index_(r));
   }
 
  private:
  IndirectList<T,Map>& operator= (const IndirectList<T,Map>&);
  /* IndirectList's can't be assigned 
     Copy construction OK, if slow (default constructor works) */

  BaseList<T> datast_;
  Map index_;
  const size_t length_;
  };


 template<typename T,typename Map> std::ostream& operator<< (std::ostream& ostr, const IndirectList<T,Map>& slar)
   {
     print(slar,ostr);
     return ostr;
   }
 
 template<typename T,typename Map> struct type_traits< IndirectList<T,Map> > {
   static const int rank=0;
   static const bool trivialconstructor=false;
   static const size_t dimensionality=1;
   typedef T value_type;
 };
    
 //implementation details below here
  
template<class T,class Map> template<class T2> IndirectList<T,Map>& IndirectList<T,Map>::operator+= (const T2& a) { return add_ip(*this,a); } 
template<class T,class Map> template<class T2> IndirectList<T,Map>& IndirectList<T,Map>::operator-= (const T2& a) { return subtract_ip(*this,a); } 
template<class T,class Map> template<class T2> IndirectList<T,Map>& IndirectList<T,Map>::operator*= (const T2& a) {
  apply_ip2(doesmultiply_ip<T,T2>(),*this,a);
  return *this;
} 
template<class T,class Map> template<class T2> IndirectList<T,Map>& IndirectList<T,Map>::operator/= (const T2& a) {
  apply_ip2(doesdivide_ip<T,T2>(),*this,a);
  return *this;
} 
 template<class T,class Map> template<typename T2> IndirectList<T,Map>& IndirectList<T,Map>::operator= (const T2& a) {
   assign(*this,a);
   return *this;
 }

  class slice : public std::unary_function<size_t,size_t> {
  public:
    slice() : stride_(1) { start_=size_=max_=0; }
    slice(int startv, int sizev, int stridev =1) { create(startv,sizev,stridev); }
    slice(const range& a) //!< allow implicit conversion
    { create(a.start(),a.size(),1); }

    bool empty() const { return (size_==0); }
    size_t start() const { return start_; }
    size_t size() const { return size_; }
    int stride() const { return stride_; }
    size_t max() const { return max_; }

    static const size_t* vector() { return NULL; } //!< dummy for issame check

    size_t operator()(size_t ind) const {
#if (BOUNDS_CHECK)
      if (ind>=size_)
	throw BadIndex("slice()",ind,size_);
#endif
      return start_+ind*stride_;
    }

    slice& operator*= (int);
    slice& operator+= (int);

    class const_iterator :  public ::std::iterator< ::std::random_access_iterator_tag,size_t,ptrdiff_t,const size_t*,const size_t&> {
    public:
     const_iterator(size_t curv,size_t strv) : cur(curv),str(strv) {}

     const_iterator& operator++() { cur+=str; return *this; }
     const_iterator operator++(int) { cur+=str; return const_iterator(cur-str,str); }
     const_iterator& operator--() { cur-=str; return *this; }
     const_iterator operator--(int) { cur-=str; return const_iterator(cur+str,str); }
     const_iterator& operator+=(ptrdiff_t i) { cur+=i*str; return *this; }
     const_iterator& operator-=(ptrdiff_t i) { cur-=i*str; return *this; }
     const_iterator operator+ (ptrdiff_t i) const { return const_iterator(cur+i*str,str); }
     const_iterator operator- (ptrdiff_t i) const { return const_iterator(cur-i*str,str); }

     size_t operator*() const { return cur; }
     bool operator== (const const_iterator& x) const { return (x.cur==cur); }
     bool operator!= (const const_iterator& x) const { return (x.cur!=cur); }
     bool operator> (const const_iterator& x) const { return (str>0 ? x.cur>cur : cur>x.cur); }
     bool operator< (const const_iterator& x) const { return (str>0 ? x.cur<cur : cur>x.cur); }
     bool operator>= (const const_iterator& x) const { return (str>0 ? x.cur>=cur : cur>=x.cur); }
     bool operator<= (const const_iterator& x) const { return (str>0 ? x.cur<=cur : cur<=x.cur); }
     ptrdiff_t operator- (const const_iterator& x) const { return (int(cur)-int(x.cur))/str; }
     size_t operator[] (ptrdiff_t i) const { return cur+i*str; }
    private:
     size_t cur;
     int str;
   };
    typedef const_iterator iterator;

    const_iterator begin() const { return const_iterator(start_,stride_); }
    const_iterator begin() { return const_iterator(start_,stride_); }
    const_iterator end() const { return const_iterator(max_+stride_,stride_); }
    const_iterator end() { return const_iterator(max_+stride_,stride_); }

  private:
    size_t start_,size_,max_;
    int stride_; //must be signed, is never zero

    void create(int,int,int);
 };

template<> struct type_traits<slice> {
  static const bool trivialconstructor=true;
  static const size_t dimensionality=1;
  static const int rank=0;
  typedef size_t value_type;
};

  inline std::ostream& operator<< (std::ostream& ostr, const slice& sl) {
    return ostr << "(size:" << sl.size() << "  start:" << sl.start() << "  stride: " << sl.stride() << ')';
  }

  inline size_t max(const slice& sl) { return sl.max(); }


  inline slice operator* (const range& a,int v) { return slice(a.start()*v,a.size(),v); }
  inline slice operator* (int v,const range& a) { return slice(a.start()*v,a.size(),v); }
  inline slice operator* (const slice& a,int v) { return slice(a.start()*v,a.size(),a.stride()*v); }
  inline slice operator* (int v,const slice& a) { return slice(a.start()*v,a.size(),a.stride()*v); }

}//namespace libcmatrix

#endif //#ifdef _IndirectList_h_
