#ifndef _ScratchList_h_
#define _ScratchList_h_

//#include "List.h"
//#include "CachedAllocator.h"

namespace libcmatrix {

//  template<class T> struct scratchlist_traits {
//     typedef DefaultAllocator allocator;
//   };

/* template<> struct scratchlist_traits<COMPLEX_T> {  */
/*   typedef CachedAllocator<ThreadingActive,COMPLEX_T> allocator;  */
/* };  */
/* template<> struct scratchlist_traits<float_t> {  */
/*   typedef CachedAllocator<ThreadingActive,float_t> allocator;  */
/* };  */

template<typename T,int size_ =SCRATCH_SIZE> class ScratchList : public BaseList<T> {
public:
typedef MemoryFiller<T,type_traits<T>::trivialconstructor, DefaultAllocator< memory_traits<T>::alignment > > alloc_t;
  typedef typename alloc_t::memop_t memop_t;

  ScratchList() : BaseList<T>() {}
  ScratchList(int n,const T& v) {
    docreate_uninit(n);
    memop_t::construct(this->datast,this->datast+n,v);
  }
  ~ScratchList() { clear(); }

  void clear() {
    if ((this->nitems)>size_)
      alloc_t::release(this->vector(),this->nitems,this->nitems);
    this->nitems=0;
  }

  template<class T2> explicit ScratchList(const T2& a) {
    make_(Int2Type<LCM_DIM(T2)>(),a);
  }
  ScratchList(const ScratchList<T,size_>& a) {  //Mustn't use default constructor!
    docreate_iter(a.size(),a.vector());
  }

  template<class T2> ScratchList(const T2& v0,const T2& v1) {
    docreate_uninit(2);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
  }
  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2) {
    docreate_uninit(3);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
  }
  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3) {
    docreate_uninit(4);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
  }
  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3,const T2& v4) {
    docreate_uninit(5);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
    memop_t::construct((this->datast)[4],v4);
  }
  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3,const T2& v4,const T2& v5) {
    docreate_uninit(6);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
    memop_t::construct((this->datast)[4],v4);
    memop_t::construct((this->datast)[5],v5);
  }
  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3,const T2& v4,const T2& v5,const T2& v6) {
    docreate_uninit(7);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
    memop_t::construct((this->datast)[4],v4);
    memop_t::construct((this->datast)[5],v5);
    memop_t::construct((this->datast)[6],v6);
  }
  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3,const T2& v4,const T2& v5,const T2& v6,const T2& v7) {
    docreate_uninit(8);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
    memop_t::construct((this->datast)[4],v4);
    memop_t::construct((this->datast)[5],v5);
    memop_t::construct((this->datast)[6],v6);
    memop_t::construct((this->datast)[7],v7);
  }
  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3,const T2& v4,const T2& v5,const T2& v6,const T2& v7,const T2& v8) {
    docreate_uninit(9);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
    memop_t::construct((this->datast)[4],v4);
    memop_t::construct((this->datast)[5],v5);
    memop_t::construct((this->datast)[6],v6);
    memop_t::construct((this->datast)[7],v7);
    memop_t::construct((this->datast)[8],v8);
  }
  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3,const T2& v4,const T2& v5,const T2& v6,const T2& v7,const T2& v8,const T2& v9) {
    docreate_uninit(10);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
    memop_t::construct((this->datast)[4],v4);
    memop_t::construct((this->datast)[5],v5);
    memop_t::construct((this->datast)[6],v6);
    memop_t::construct((this->datast)[7],v7);
    memop_t::construct((this->datast)[8],v8);
    memop_t::construct((this->datast)[9],v9);
  }

  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3,const T2& v4,const T2& v5,const T2& v6,const T2& v7,const T2& v8,const T2& v9, const T2& v10) {
    docreate_uninit(11);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
    memop_t::construct((this->datast)[4],v4);
    memop_t::construct((this->datast)[5],v5);
    memop_t::construct((this->datast)[6],v6);
    memop_t::construct((this->datast)[7],v7);
    memop_t::construct((this->datast)[8],v8);
    memop_t::construct((this->datast)[9],v9);
    memop_t::construct((this->datast)[10],v10);
  }

  template<class T2> ScratchList(const T2& v0,const T2& v1,const T2& v2,const T2& v3,const T2& v4,const T2& v5,const T2& v6,const T2& v7,const T2& v8,const T2& v9, const T2& v10, const T2& v11) {
    docreate_uninit(12);
    memop_t::construct((this->datast)[0],v0);
    memop_t::construct((this->datast)[1],v1);
    memop_t::construct((this->datast)[2],v2);
    memop_t::construct((this->datast)[3],v3);
    memop_t::construct((this->datast)[4],v4);
    memop_t::construct((this->datast)[5],v5);
    memop_t::construct((this->datast)[6],v6);
    memop_t::construct((this->datast)[7],v7);
    memop_t::construct((this->datast)[8],v8);
    memop_t::construct((this->datast)[9],v9);
    memop_t::construct((this->datast)[10],v10);
    memop_t::construct((this->datast)[11],v11);
  }


  ScratchList<T,size_>& operator= (const ScratchList<T,size_>& a) { // Not exception-safe
    create(a.size());
    BaseList<T>::operator=(a);
    return *this;
  }

  template<class T2> ScratchList<T,size_>& operator= (const T2& a) {
    create(a.size());
    BaseList<T>::operator=(a);
    return *this;
  }
  ScratchList<T,size_>& operator= (const T& a) {
    BaseList<T>::operator=(a);
    return *this;
  }
  void create(int n) {
    if (n!=this->nitems) {
      clear();
      docreate(n);
    }
  }

private:
#if (SCRATCH_SIZE!=0)

#ifdef LCM_NEED_128BIT_ALIGN
char unaligned_stat[size_*sizeof(T)+15];
#else
long stat[(size_*sizeof(T)+sizeof(long)-1)/sizeof(long)];
#endif

  void docreate_uninit(int n) {
    T* p= (n>size_)
      ? alloc_t::get_uninit(n) 
#ifdef LCM_NEED_128BIT_ALIGN
      : reinterpret_cast<T*>( (size_t(unaligned_stat)+15) & ~size_t(15));
#else
    : reinterpret_cast<T*>(stat);
#endif
    BaseList<T>::create(n,p);
  }
#else
  void docreate_uninit(int n) {
    BaseList<T>::create(n,alloc_t::get_uninit(n));
  }
#endif

  void docreate(int n) {
    docreate_uninit(n);
    memop_t::construct(this->datast,this->datast+n);
  }
  template<class Iter> void docreate_iter(int n,Iter v) {
    if (n) {
      docreate_uninit(n);
      memop_t::construct_iter(this->datast,this->datast+n,v);
    }
    else
      BaseList<T>::create(0U,(T*)NULL);
  }

  void make_(Int2Type<0>, int n) {
    docreate(n);
  }
  template<class T2> void make_(Int2Type<1>, const T2& v) {
    docreate_iter(v.size(),v.begin());
  }
};

 template<typename T, int N> struct type_traits< ScratchList<T,N> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef T value_type;
 };

} //namespace libcmatrix 

#endif

