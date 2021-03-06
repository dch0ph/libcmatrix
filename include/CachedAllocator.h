/* Caches addresses of blocks of fixed size
   Borrows from SmallObjAllocator (Modern C++ design)
   but uses new for actual memory allocation */

#ifndef _CachedAllocator_h_
#define _CachedAllocator_h_

#include <cassert>
#include <vector>
#include <iosfwd>
#include <algorithm>
#include "lcm_basethreads.h"
#include "List.h"
#include "Warnings.h"

namespace libcmatrix {

  //NB check mutex
  template<class T, class Arena> class SingletonHolder {
  public:
    template<bool MustExist> static T& Instance(Bool2Type<MustExist>) {
      static Mutex<ThreadingActive> lock;
      if (!pInstance_) {
	if (MustExist) {
	  std::cerr << "Dead SingletonHolder detected" << std::endl;
	  throw InternalError("This shouldn't happen!:2");
	}
	lock.acquire();
	if (!pInstance_) {
	  if (destroyed_) {
	    std::cerr << "Attempt to use destroyed SingletonHolder" << std::endl;
	    throw InternalError("This shouldn't happen!");
	  }
	  Create();
	}
	lock.release();
      }
      return *pInstance_;
    }
    
  private:
    SingletonHolder();
    virtual ~SingletonHolder() {
      pInstance_=NULL;
      destroyed_=true;
    }
    static void Create() {
      static T theInstance;
      pInstance_=&theInstance;
    }
    typedef T InstanceType;
    static InstanceType* pInstance_;
    static bool destroyed_;
  };

  template<class T, class A> bool SingletonHolder<T,A>::destroyed_=false;
  template<class T, class A> typename SingletonHolder<T,A>::InstanceType* SingletonHolder<T,A>::pInstance_=NULL;

    class FixedAllocator {
    private:
      const size_t size_;
      //use basic allocator, otherwise we'll get in a mess...
      List<void*> stack_;
      size_t failedcount;
      size_t successcount;

      //Make these private - we don't really want to be copying all those pointers...
      FixedAllocator& operator=(const FixedAllocator& a);
      FixedAllocator(const FixedAllocator&); 
    public:
      FixedAllocator(size_t sizev, size_t totalv);
      ~FixedAllocator();

      size_t BlockSize() const { return size_; }
      size_t BlockCount() const { return stack_.size(); }
      void* allocate();
      void deallocate(void*);
      void print(std::ostream&) const;
    };

namespace { // anonymous 
  struct CompareFixedAllocatorSize
    : std::binary_function<const FixedAllocator* const, size_t, bool> {
    bool operator()(const FixedAllocator* const x, size_t numBytes) const {
      return x->BlockSize() < numBytes;
    }
    //! keeps debugging code of Visual C++ compiler happy
#ifdef _DEBUG
    bool operator()(size_t numBytes,const FixedAllocator* const x) const {
      return !((*this)(x,numBytes));
    }
    bool operator()(const FixedAllocator* const, const FixedAllocator* const) const {
      return false;
    }    
#endif
  };
} 

  class CachedAllocator_ {
    size_t max_;
    size_t total_;

    typedef std::vector<FixedAllocator*> pool_t;
    pool_t pool_;
    FixedAllocator* pLastAlloc_;
    FixedAllocator* pLastDealloc_;
    bool destroyed_;
    size_t oversize_claimed_words_, oversize_released_words_;

    CachedAllocator_(const CachedAllocator_&);
    CachedAllocator_& operator=(const CachedAllocator_&);

    static Warning<> post_destruction_warning;
        
  public:
    CachedAllocator_(int maxkv =LCM_CACHE_MAXK, int totalkv =LCM_CACHE_TOTALK);
    ~CachedAllocator_();
    void setcachelimits(int maxkv, int totalkv);
    void* allocate(size_t);
    void deallocate(void*, size_t);
    void print(std::ostream&) const;
    size_t allocated_words() const; //!< return net number of words currently allocated
  };

  //Arena allows multiple CachedAllocators to be created
  template<bool,class Arena =NullType> class CachedAllocator {
    public:
    typedef SingletonHolder<CachedAllocator_,Arena> holder_t;

    static void* allocate(size_t length)
    { return holder_t::Instance(Bool2Type<false>() ).allocate(length); }

    static void deallocate(void* p,size_t length)
    { holder_t::Instance( Bool2Type<true>() ).deallocate(p,length); }

    static void print(std::ostream& ostr) {
      holder_t::Instance(Bool2Type<false>() ).print(ostr);
    }

    static size_t allocated_words() {
      return holder_t::Instance(Bool2Type<false>() ).allocated_words();
    }

  };

  template<class Arena> class CachedAllocator<true,Arena> {
  public:
  typedef SingletonHolder<key_holder<CachedAllocator_>,Arena> holder_t;

  static void* allocate(size_t length)
  {
    key_holder<CachedAllocator_>& holder(holder_t::Instance(Bool2Type<false>() ));
    CachedAllocator_* allocp(holder.get()); //should now be thread-specific i.e. locking not reqd
    if (!allocp) { //stack doesn't exist?
      allocp=new CachedAllocator_();
      holder.set(allocp);
    }
    return allocp->allocate(length);
  }

  static void deallocate(void* p, size_t length) {
    CachedAllocator_* allocp(holder_t::Instance( Bool2Type<true>() ).get());
    if (allocp)
      allocp->deallocate(p,length);
    else
      CachedAllocator_::post_destruction_warning.raise();
  }
  
  static size_t allocated_words() {
    const CachedAllocator_* allocp(holder_t::Instance(Bool2Type<false>() ).get());
    return allocp ? allocp->allocated_words() : 0;
  }

  static void print(std::ostream& ostr) {
    CachedAllocator_* allocp(holder_t::Instance(Bool2Type<false>() ).get());
    if (allocp) {
      ostr << "Allocator address: " << allocp << '\n';
      allocp->print(ostr);
    }
    else
      ostr << "CachedAllocator_ not created";
  }
  };

}//namespace libcmatrix

#endif // _CachedAllocator_h_
