#ifndef DynamicList_h_
#define DynamicList_h_

#include <cstdio>
#include "BaseList.h"
#include "lcm_basethreads.h"
#include "Warnings.h"

namespace libcmatrix {
  
  extern Warning<> dynamiclist_alias_warning;

  template<class T, class AllocatorPolicy = DefaultAllocator< memory_traits<T>::alignment > >
    class DynamicList : public BaseList<T> {
    public:

     inline bool istemporary() const
     { return (store_type==mxflag::temporary); }
     
     inline bool isdynamic() const
     { return (store_type!=mxflag::nondynamic); }
     
     inline bool hasclaim() const
    { return (this->datast!=NULL) && (store_type!=mxflag::nondynamic); }
     
     ~DynamicList() {
       if (hasclaim())
	 alloc_t::release(this->datast,this->nitems,this->nitems);
     }
     
     void clear() {
       if (!(this->datast))
	 return;
       if (store_type==mxflag::nondynamic) 
	 throw Failed("Can't explicitly kill non-dynamic objects"); 
       alloc_t::release(this->datast,this->nitems,this->nitems);
       zap();
     }

    explicit DynamicList(mxflag::tempflag tflag =mxflag::normal)
       : BaseList<T>()
     {
       set(true,tflag); //'have memory' to avoid nondynamic tripping out
     }
     
     template<class T2> DynamicList(const T2& a, mxflag::tempflag tflag)
    { construct_(Int2Type<LCM_DIM(T2)>(),a,tflag); }

    DynamicList(int n, const T& v, mxflag::tempflag tflag) {
       set(false,tflag);
       BaseList<T>::create(n,alloc_t::get(n,n,v));
     }
        
    template<class T2> DynamicList(int n, const T2& v, mxflag::tempflag tflag) {
       set(false,tflag);
       BaseList<T>::create(n,alloc_t::get(n,n,T(v)));
     }
     
     DynamicList(int n, T* vp, mxflag::tempflag tflag) {
       set(true,tflag);
       docreate(n,vp);
     }
     
     DynamicList(int n, const T* vp, mxflag::tempflag tflag) {
       set(true,tflag);
       docreate(n,const_cast<T*>(vp));
     }
     
     //NB temporary_ attribute is NOT exchanged!
     void swap(DynamicList& a) {
       if ( (store_type==mxflag::nondynamic) ^ (a.store_type==mxflag::nondynamic)) {
	 if (store_type!=mxflag::nondynamic)
	   throw Failed("swap: cannot swap non-dynamic object into dynamic object");
	 *this=a;
	 a.clear(); //not a strict swap; flag this by clearing a (which is probably just a temporary)
       }
       else
	 BaseList<T>::swap(a);
     }

    //zero initially for safety
    DynamicList(const DynamicList& a, mxflag::tempflag tflag)
    :  BaseList<T>(0,NULL) {
//        if ((a.store_type==mxflag::nondynamic) ^ (tflag==mxflag::nondynamic))
// 	 throw Failed("Can't create create dynamic object from non-dynamic or vice versa; use reference!");
       set(true,tflag);
       if (a.datast==NULL) return;
#ifdef LCM_ENABLE_TEMPORARY
       if (a.store_type==mxflag::temporary) {
	 BaseList<T>::create(a.nitems,a.datast);
	 (const_cast< DynamicList<T>& >(a)).zap();
       }
       else {
#endif
	 docreate(a.nitems,a.datast);
#ifdef LCM_ENABLE_TEMPORARY
       }
#endif
     }
     
    //zero initially for safety
    DynamicList(const DynamicList& a)
    :  BaseList<T>(0,NULL), store_type(a.store_type) {
#ifndef NDEBUG
      if ((store_type==mxflag::nondynamic) && (a.datast!=NULL))
	dynamiclist_alias_warning.raise();
#endif      
#ifdef LCM_ENABLE_TEMPORARY
      if (store_type==mxflag::temporary) {
	BaseList<T>::create(a.nitems,a.datast);
	(const_cast< DynamicList<T>& >(a)).zap();
	store_type=mxflag::normal; //type is now normal
      }
      else {
#endif
	if (a.datast)
	  docreate(a.nitems,a.datast);
#ifdef LCM_ENABLE_TEMPORARY
      }
#endif
    }
     
     DynamicList& operator= (const DynamicList& a) {
       if (this!=&a) {
	 if (a.empty())
	   clear();
	 else
	   create(a.nitems,a.datast);
       }
       return *this;
     }
    template<class T2> DynamicList& operator= (const T2& a) {
      assign(*this,a);
      return *this;
    }
    
    bool operator!() const { return this->empty(); }
     
     bool create(int n) {
       if (this->nitems==n)
	 return false;
       replacevector(n,alloc_t::get(n,n));
       return true;
     }

   void create(int n, const T& v) {
     if ((store_type==mxflag::nondynamic) && (n!=this->nitems))
       throw Mismatch("DynamicList::create");
     if (n==this->nitems)
       alloc_t::memop_t::assign(this->datast,this->datast+n,v);
     else
       replacevector(n,alloc_t::get(n,n,v));
   }

    void create(size_t n, T* newd) {
      if ((store_type==mxflag::nondynamic) && (n!=this->nitems))
	throw Mismatch("DynamicList::create");
      if (n==this->nitems)
	alloc_t::memop_t::assign(this->datast,this->datast+n,newd);
      else {
	T* newv=alloc_t::get(n,n,newd);
	//	std::cout << "created: " << newv << std::endl;
	replacevector(n,newv);
      }
    } 
    
    template<class Iter> void create_iter(size_t n, const Iter& iter) {
      if ((store_type==mxflag::nondynamic) && (n!=this->nitems))
	throw Mismatch("DynamicList::create_iter");
      if (n==this->nitems)
	alloc_t::memop_t::assign(this->datast,this->datast+n,iter);
      else
	replacevector(n,alloc_t::get_iter(n,n,iter));
    }

      void reserve(size_t); //change allocitems while preserving contents

   private:
     mxflag::tempflag store_type;

    typedef MemoryFiller<T,type_traits<T>::trivialconstructor,AllocatorPolicy> alloc_t;

    void zap() { 
      this->nitems=0;
      this->datast=NULL;
    }

    void replacevector(size_t n, T* newd) { //only valid for pointers created by alloc_t::get
     if (this->datast) {
       if (store_type==mxflag::nondynamic) {
	 if (newd)
	   alloc_t::release(newd,n,n);
	 throw Failed("DynamicList<T>: cannot resize nondynamic object");
       }
       alloc_t::release(this->datast,this->nitems,this->nitems);
     }
     BaseList<T>::create(n,newd);
   }

     void docreate(int n, T* vp)
     {
       BaseList<T>::create(n,(store_type==mxflag::nondynamic) ? vp : alloc_t::get(n,n,vp));
     }
     
    void construct_(Int2Type<0>, int n, mxflag::tempflag flag) {  
      set(false,flag);
      BaseList<T>::create(n,alloc_t::get(n,n)); 
    }

    template<class T2> void construct_(Int2Type<1>, const T2& a, mxflag::tempflag flag) {
      set(false,flag);
      const size_t n=a.size();
      BaseList<T>::create(n,n ? alloc_t::get_iter(n,n,a.begin()) : NULL);
    }

     void set(bool havememory, mxflag::tempflag tflag) {
       if (!havememory && (tflag==mxflag::nondynamic))
	 throw Failed("DynamicList<T>(): cannot create empty nondynamic object"); 
#ifdef LCM_ALLOW_TEMPORARY
       store_type=tflag;
#else
       store_type=(tflag==mxflag::temporary) ? mxflag::normal : tflag;
#endif
     }
 };

  template<typename T,typename Alloc> void DynamicList<T,Alloc>::reserve(size_t n) //change allocitems while preserving contents
  {
    if ((store_type==mxflag::nondynamic) && (n!=this->nitems))
      throw Failed("Cannot resize nondynamic object");
    if (n>this->nitems) {
      if (this->nitems) {
	const size_t olditems=this->nitems;
	replacevector(n,alloc_t::get(n,olditems,this->datast));
	alloc_t::memop_t::construct(this->datast+olditems,this->datast+n);
      }
      else
	replacevector(n,alloc_t::get(n,n));
    }
  }
  
  template<typename T,typename Alloc> struct type_traits< DynamicList<T,Alloc> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef T value_type;
 };

  //IO-related exceptions
  struct NoComplexRead : public Failed {
    NoComplexRead(const char* errmv =NULL) : Failed(errmv ? errmv : "Can't read complex matrix into non-complex") {}
  };
  struct FileCorrupt : public Failed {
    FileCorrupt(const char* errmv =NULL) : Failed(errmv ? errmv : "Cannot parse file: corrupt?") {}
  };
  struct Unsupported : public Failed {
    Unsupported(const char* errmv) : Failed(errmv) {}
  };
  struct OpenFailed : public Failed {
    OpenFailed(const char* errmv) : Failed(errmv) {}
  };

  struct FileHandle {
    FileHandle(const char* name, const char* open)
      : fp(fopen(name,open)) {
      if (!fp)
	throw OpenFailed(name);
    }
    ~FileHandle() { fclose(fp); }
    FILE* operator()() const { return fp; }
    FILE* fp;
  };

#ifdef LCM_USE_EXTERNTEMPLATE
LCM_L1_EXTERN template class DynamicList<double>;
LCM_L1_EXTERN template class DynamicList<complex>;
LCM_L1_EXTERN template class DynamicList<size_t>;
#endif

} //namespace libcmatrix

#endif
