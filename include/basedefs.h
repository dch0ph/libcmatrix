// This file should be included by all libcmatrix functions
#ifndef _basedefs_h_
#define _basedefs_h_

#include "config.h"
#include <iostream>
#include <iomanip>
#include <new>
#include <stdexcept>
#include <cmath>
#include <cstring> // for memcpy
#include <string>
#include <cstddef> // for size_t
#include <functional>
#include <iterator>
#include <algorithm>

// #ifdef LCM_DEBUG_ALLOCATOR
// #include <malloc.h>
// #ifdef DMALLOC
// #include <dmalloc.h>
// #endif
// #ifdef DBALLOC
// #include <dballoc.h>
// #endif
// #endif

//if LCM_INLINE not already defined, define on basis of LCM_FORCE_INLINE
#ifndef LCM_INLINE
#ifdef LCM_FORCE_INLINE
#define LCM_INLINE inline
#else
#define LCM_INLINE
#endif
#endif

#define LCM_RESTRICT
//#define LCM_RESTRICT __restrict__

#define LCM_ENUM_CAST(X) int(X)
#define LCM_ITER_CAT(X) typename ::std::iterator_traits< X >::iterator_category

#ifdef LCM_USE_NULLPTR
#define LCM_NULL nullptr
#else
#define LCM_NULL NULL
#endif

namespace libcmatrix {

  extern const char cmatrix_abi_version[]; //!< string (X.Y.Z) identifying library ABI (X major incompatability, Y minor breaking changes, Z fixes only) N.B. first defined for 3.1.0

// Exception mix-in class;
  class MatrixException {
  public:
    explicit MatrixException(const char* type_)
    : type(type_) { hook(); }
    MatrixException(const char* type_, const char* errm_)
      : type(type_), errm(errm_) { hook(); }
    MatrixException(const char* type_, const std::string& errm_)
      : type(type_), errm(errm_) { hook(); }

    friend std::ostream& operator << (std::ostream&, const MatrixException&);
    const char* what() const { return errm.c_str(); }
 private:
   const char* type;
    void hook() const; //!< hook for debugging
  protected:
    std::string errm;
  };

 std::ostream& operator << (std::ostream&, const MatrixException&);

// Exceptions shared by many components
  struct Failed : public MatrixException, public std::runtime_error { 
    virtual ~Failed() throw() {};
    explicit Failed(const char* errm_)
      : MatrixException("Failed",errm_),
	runtime_error(errm_) {}
  };

  struct InvalidParameter : public MatrixException, public std::invalid_argument {
    virtual ~InvalidParameter() throw() {};
    explicit InvalidParameter(const char* errm_)
      : MatrixException("Invalid parameter",errm_),
	std::invalid_argument(errm_) {}
  };

  struct InternalError : public MatrixException, public std::logic_error {
    virtual ~InternalError() throw() {};
    explicit InternalError(const char* errm_)
      : MatrixException("Internal error",errm_),
	std::logic_error(errm_) {}
  };

 struct BadIndex : public MatrixException, public std::out_of_range {
   virtual ~BadIndex() throw() {};
   explicit BadIndex(const char* errm_)
     : MatrixException("Bad index",errm_),
       std::out_of_range(errm_) {}
   BadIndex(const char* errm_, int arg_, size_t max_);
   BadIndex(const char* errm_, int arg1_, size_t max1_, int arg2_, size_t max2_);
   BadIndex(const char* errm_, int arg1_, size_t max1_, int arg2_, size_t max2_, int arg3_, size_t max3_);
  };

  struct Mismatch : public MatrixException, public std::domain_error {
    virtual ~Mismatch() throw() {};
    explicit Mismatch(const char* errm_) 
      : MatrixException("Mismatch",errm_),
	std::domain_error(errm_) {}
    Mismatch(const char* errm_, size_t v1, size_t v2);
    Mismatch(const char* errm_, size_t r1, size_t c1, size_t r2, size_t c2);
  };

  struct Undefined : public MatrixException, public std::domain_error {
    virtual ~Undefined() throw() {};
    explicit Undefined(const char* errm_)
      : MatrixException("Undefined",errm_),
	std::domain_error(errm_) {}
  };

  struct ArgumentClash : public MatrixException, public std::invalid_argument {
    virtual ~ArgumentClash() throw() {};
    explicit ArgumentClash(const char* errm_)
      : MatrixException("Input used as output",errm_),
	std::invalid_argument(errm_) {}
  };

  struct NotSquare : public MatrixException, public std::domain_error {
    virtual ~NotSquare() throw() {};
    explicit NotSquare(const char* errm_)
      : MatrixException("Matrix is not square/undefined",errm_),
	std::domain_error(errm_) {}
  };

//common forward declarations.  Put in one place in case defns change

#ifdef LCM_USE_STDCOMPLEX
  typedef ::std::complex<double> complex;
#else
  class complex;
#endif

template<typename T> class Matrix;
template<typename T> class BaseList;
 template<typename T> class List;

//Core of type/math engine
template<typename T> struct type_traits {
  static const bool trivialconstructor=false; 
  static const size_t dimensionality=0;
  static const size_t rank=0;
  typedef T value_type;
};

  // used to control memory allocation
  template<typename T> struct memory_traits {
    static const size_t alignment=0;
  };
    
 template<bool trivial> struct isundefined_ {
   template<class T> inline static bool isundefined(const T& a) {
     return a.empty(); }
 };
 template<> struct isundefined_<true> { 
   template<class T> inline static bool isundefined(const T&) {
     return false; }
 };

template<class T> inline bool isundefined(const T& a) {
  return isundefined_< type_traits<T>::trivialconstructor>::isundefined(a); }

#define LCM_VAL(T) typename type_traits< T >::value_type
#define LCM_DIM(T) type_traits< T >::dimensionality

 template<typename T,size_t N> struct index_dim_ { static const size_t dim=0; };
 template<> struct index_dim_<size_t,1> { static const size_t dim=1; };

#define LCM_INDEX_DIM(T) index_dim_<LCM_VAL(T),LCM_DIM(T)>::dim

#define FUNDAMENTAL_(T,R) template<> struct type_traits< T > {\
 static const bool trivialconstructor=true; \
 static const size_t dimensionality=0;      \
 static const size_t rank=R; \
 typedef T value_type; };

   FUNDAMENTAL_(bool,1) 
   FUNDAMENTAL_(char,1)
   FUNDAMENTAL_(unsigned char,1)
   FUNDAMENTAL_(short,10)
   FUNDAMENTAL_(int,100)
   FUNDAMENTAL_(unsigned int,200)
   FUNDAMENTAL_(long,300)
   FUNDAMENTAL_(unsigned long,400)
   FUNDAMENTAL_(float,500)
   FUNDAMENTAL_(double,600)
   FUNDAMENTAL_(long double,700)

     template<typename T> struct isinteger {
       static const bool result=(type_traits<T>::rank>=100) && (type_traits<T>::rank<=400);
     };

/* Adapted from promote.h of Blitz++ */

template<class T,bool autopro> struct autopromote_trait { typedef T value_type; };

#define LCM_DECLARE_AUTOPROMOTE(T1,T2)     \
    template<> struct autopromote_trait<T1,true> { typedef T2 value_type; };

// These are the odd cases where small integer types are
// automatically promoted to int or unsigned int for arithmetic.
LCM_DECLARE_AUTOPROMOTE(bool, int)
LCM_DECLARE_AUTOPROMOTE(char, int)
LCM_DECLARE_AUTOPROMOTE(unsigned char, int)
LCM_DECLARE_AUTOPROMOTE(short int, int)
LCM_DECLARE_AUTOPROMOTE(short unsigned int, unsigned int)

template<class T1, class T2,int promoteToT1> struct _bz_promote2 {};
template<class T1, class T2> struct _bz_promote2<T1,T2,1> { typedef T1 value_type; };
template<class T1, class T2> struct _bz_promote2<T1,T2,0> { typedef T2 value_type; };

//promote_trait is only valid for base types - it does not recurse through container types
 template<class T1_orig, class T2_orig,bool autopro =true> struct promote_trait {
    // Handle promotion of small integers to int/unsigned int
    typedef typename autopromote_trait<T1_orig,autopro>::value_type T1;
    typedef typename autopromote_trait<T2_orig,autopro>::value_type T2;

    enum {
    // True if T1 is higher ranked
      T1HigherRank =
        (LCM_ENUM_CAST(type_traits<T1>::rank) >
	 LCM_ENUM_CAST(type_traits<T2>::rank)) ? 1 : 0,

    // True if we know ranks for both T1 and T2
      knowRanks =
        LCM_ENUM_CAST(type_traits<T1>::rank) && LCM_ENUM_CAST(type_traits<T2>::rank)
    };

    enum {
      isfirstbetter = (LCM_ENUM_CAST(knowRanks) ? LCM_ENUM_CAST(T1HigherRank) : -100)
    };
    typedef typename _bz_promote2<T1,T2,isfirstbetter>::value_type value_type; 
};
 template<class T,bool autopro> struct promote_trait<T,T,autopro> {
    typedef typename autopromote_trait<T,autopro>::value_type value_type;
 };

#define LCM_NEWTYPE(T1,T2,autopro) typename promote_trait< T1,T2,autopro >::value_type
 
 //Pinched from Loki / Modern C++ Design
 template<bool> struct CompileTimeError;
 template<> struct CompileTimeError<true> {};
#define LCM_STATIC_CHECK(expr, msg) \
   { libcmatrix::CompileTimeError<((expr) != 0)> Error_##msg; (void)Error_##msg; }
#define LCM_STATIC_ERROR(msg) \
 { libcmatrix::CompileTimeError<false> Error_##msg; (void)Error_##msg; }

template<typename T> 
  struct Type2Type {
    typedef T OrginalType;
  };

template<int v>
  struct Int2Type {
    enum { value=v };
  };

template<bool v>
  struct Bool2Type {
    enum { value=0 };
  };
template<>
  struct Bool2Type<true> {
    enum { value=1 };
  };

 struct NullType {
   bool operator!() const { return true; }
 };

//const stripper
 template<class T> struct Unconst {
   typedef T Result;
 };
 template<class T> struct Unconst<const T*> {
   typedef T* Result;
 };
 template<class T> struct Unconst<const T> {
   typedef T Result;
 };

 template<size_t N1,size_t N2> struct AreMatching_ {}; //Missing func will cause compilation to fail

  template<> struct AreMatching_<0,0> {};
  // template<class T1,class T2> inline static bool func(const T1&,const T2&) {} };
//    return true; } };

 template<> struct AreMatching_<1,1> { template<class T1,class T2> inline static bool func(const T1& a,const T2& b) {
   return (a.size()==b.size()); } };

 template<> struct AreMatching_<2,2> { template<class T1,class T2> inline static bool func(const T1& a,const T2& b) {
   return (a.rows()==b.rows()) && (a.cols()==b.cols()); } };

 template<size_t N> struct AreMatching_<N,N> { template<class T1,class T2> inline static bool func(const T1& a,const T2& b) {
   for (size_t i=N;i--;) {
     if (a.dimension(i)!=b.dimension(i))
       return false;
   }
   return true; }
 };

 template<typename T1,typename T2> inline bool arematching(const T1& a,const T2& b) {
   return AreMatching_<LCM_DIM(T1),LCM_DIM(T2)>::func(a,b); }

  template<size_t N1,size_t N2,typename M1, typename M2, bool> struct issame_ { 
  template<class T1,class T2> static bool issame(const T1&,const T2&) {
    return false; }
};
  template<size_t N, typename M> struct issame_<N,N,M,M,false> {
  template<class T1,class T2> static bool issame(const T1& a, const T2& b) {
    // return false; }
    return !(a.empty()) && (a.vector()==b.vector()); }
  template<class T> static bool issame(const T& a,const T& b) {
    return !(a.empty()) && (a.vector()==b.vector()); }
};
  template<typename M> struct issame_<0,0,M,M,false> { 
   template<class T1,class T2> static bool issame(const T1&,const T2&) {
     return false; }
 };

/* template<class M> struct issame_<0,M,false,0,M,true> {  */
/*   template<class T1,class T2> static bool issame(const T1&,const T2&) { */
/*     return false; } */
/* }; */

/* template<class M> struct issame_<0,M,true,0,M,false> {  */
/*   template<class T1,class T2> static bool issame(const T1&,const T2&) { */
/*     return false; } */
/* }; */
/* template<class M> struct issame_<0,M,true,0,M,true> {  */
/*   template<class T1,class T2> static bool issame(const T1&,const T2&) { */
/*     return false; } */
/* }; */

template<class T1,class T2> inline bool issame(const T1& a,const T2& b) {
  //  return issame_<LCM_DIM(T1),LCM_VAL(T1),type_traits<T1>::trivialconstructor,LCM_DIM(T2),LCM_VAL(T2),type_traits<T2>::trivialconstructor>::issame(a,b); }
  return issame_<LCM_DIM(T1),LCM_DIM(T2), LCM_VAL(T1), LCM_VAL(T2), type_traits<T1>::trivialconstructor && type_traits<T2>::trivialconstructor>::issame(a,b);
}

 template<size_t N1,size_t N2> struct dup_struct_ {};
 template<size_t N> struct dup_struct_<N,N> {
   template<typename T1,typename T2> inline static void func(T1& d, const T2& a) {
     d.set_dimensions(a);
   }
 };
 template<> struct dup_struct_<2,2> {
   template<typename T1,typename T2> inline static void func(T1& d, const T2& a) {
     d.create(a.rows(),a.cols());
   }
 };
 template<> struct dup_struct_<1,1> {
   template<typename T1,typename T2> inline static void func(T1& d, const T2& a) {
     if (issame(d,a))
       throw ArgumentClash("duplicate_structure<1,1>");
//      if (isundefined(a))
//        throw Undefined("duplicate_structure<1,1>");
     d.create(a.size());
   }
 };
 template<> struct dup_struct_<0,0> {
   template<typename T1,typename T2> inline static void func(T1&, const T2&) {}
 }; 

 template<typename T1,typename T2> inline void duplicate_structure(T1& d, const T2& a) {
   dup_struct_<LCM_DIM(T1),LCM_DIM(T2)>::func(d,a); }

  /* Default allocator (does nothing clever) */
  template<size_t Align =0> struct DefaultAllocator;

#ifdef LCM_NEED_128BIT_ALIGN
#include <mm_malloc.h>

  template<size_t Align> struct DefaultAllocator {
    static void* allocate(size_t length) {
      void* p=_mm_malloc(length,Align);
      if (!p)
	throw std::bad_alloc();
      return p;
    }
    static void deallocate(void* p,size_t) { _mm_free(p); }    
    static void print(std::ostream&) {}
    static size_t allocated_words() { return size_t(0); } //!< don't know
  };
#endif
#ifdef LCM_DEBUG_ALLOCATOR
  template<> struct DefaultAllocator<0> {
    static void* allocate(size_t length) {
      void* p=malloc(length);
      if (!p)
	throw std::bad_alloc();
      return p;
    }
    static void deallocate(void* p,size_t) { free(p); }    
    static void print(std::ostream&) {}
    static size_t allocated_words() { return size_t(0); } //!< don't know
  };
#else
  template<> struct DefaultAllocator<0> {
    static void* allocate(size_t length) { return operator new(length); }
    static void deallocate(void* p,size_t) { operator delete(p); }    
    static void print(std::ostream&) {}
    static size_t allocated_words() { return size_t(0); } //!< don't know
  };
#endif

  inline void lcm_verify_getargs_(int n, int ninit) {
    if ((n<ninit) || (n<0))
      throw InvalidParameter("Bad vector initialise: claim too small or <0");
  }

  template<class T, bool =false> struct memop {

      static void construct(T* s, T* e) {
	T* b=s;
	while (b!=e) {
	  try {
	    new(b) T();
	    b++;
	  } catch(...) {
	    destroy(s,b);
	    throw;
	  }
	}
      }

      static void construct(T* s, T* e, const T& v) {
	T* b=s;
	while (b!=e) {
	  try {
	    new(b) T(v);
	    b++;
	  } catch(...) {
	    destroy(s,b);
	    throw;
	  }
	}
      }

      static void assign(T* s, T* e, const T v) {
	T* b=s;
	while (b!=e) {
	  try {
	    *b=v;
	    b++;
	  } catch(...) {
	    destroy(s,b);
	    throw;
	  }
	}
      }

     static void construct(T* s, T* e, T* vp) {
	T* b=s;
	while (b!=e) {
	  try {
	    new(b) T(*vp++);
	    b++;
	  } catch(...) {
	    destroy(s,b);
	    throw;
	  }
	}
      }

     template<typename Iter> static void construct_iter(T* s, T* e, Iter iter) {
       T* b=s;
       while (b!=e) {
	 try {
	   new(b) T(*iter);
	   ++iter;
	   b++;	 
	 } catch(...) {
	   destroy(s,b);
	   throw;
	 }
       }
     }
 
    inline static void construct(T& d,const T& v) { 
      new(&d) T(v); 
    } 

    template<class Iter> static void assign(T* b,T* e,Iter s) {
      T* start=b;
      while (b!=e) {
	try {
	  *b=*s;
	  ++s;
	  b++;
	} catch(...) {
	  destroy(start,b);
	  throw;
	}
      }
    }

    static void destroy(T* b,T* e) {
      while (b!=e) {
	b->~T();
	++b;
      }
    }

    static void destroy(T& v) { v.~T(); }
  };

  template<class T, bool fastinit, class AllocatorClass =DefaultAllocator<> > 
    class MemoryFiller {
   
  private:
  typedef AllocatorClass alloc_t;
    
  public:
      typedef memop<T,type_traits<T>::trivialconstructor> memop_t;

    static T* get_uninit(size_t n) { 
      //claimed_+=((n*sizeof(T))>>2);
      return n ? static_cast<T*>(alloc_t::allocate(n*sizeof(T))) : (T*)NULL;
    } 
      static void release_uninit(T* p, size_t n) {
	//	released_+=((n*sizeof(T))>>2);
	alloc_t::deallocate(p,n*sizeof(T));
      }

    static void release(T* p,size_t n,size_t used) { 
      memop_t::destroy(p,p+used);
      release_uninit(p,n);
    }

    static T* get(int n, int ninit) {
      lcm_verify_getargs_(n,ninit);
      T* datast=get_uninit(n);
      try {
	memop_t::construct(datast,datast+ninit);
      } catch(...) {
	release_uninit(datast,n);
	  throw;
      }
      return datast;
    }

    static T* get(int n, int ninit, const T& v) {
      lcm_verify_getargs_(n,ninit);
      T* datast=get_uninit(n);
      try {
	memop_t::construct(datast,datast+ninit,v);
      } catch(...) {
	release_uninit(datast,n);
	throw;
      }
      return datast;
    }

    static T* get(int n, int ninit, T* vp) {
      lcm_verify_getargs_(n,ninit);
      T* datast=get_uninit(n);
      try {
	memop_t::construct(datast,datast+ninit,vp);
      } catch(...) {
	release_uninit(datast,n);
	throw;
      }
      return datast;
    }
    
   template<typename Iter> static T* get_iter(int n, int ninit, const Iter& s) {
      lcm_verify_getargs_(n,ninit);
      T* datast=get_uninit(n);
      try {
	memop_t::construct_iter(datast,datast+ninit,s);
      } catch(...) {
	release_uninit(datast,n);
	throw;
      }
      return datast;
    }

      static void print(std::ostream& ostr) {
	//ostr << "Words claimed (general): " << claimed_ << "   released: " << released_ << '\n';
      }      
    };
 
#ifndef LCM_NO_EXPLOIT_TRIVIAL
  //specialised for trivial constructors

  template<class T> struct memop<T,true> {

    static void construct(T*,T*) {}

    static void construct(T* b,T* e,const T v) {
      while (b!=e)
	*b++=v;
    }
    static void assign(T* b,T* e,const T v) {
      while (b!=e)
	*b++=v;
    }
    static void construct(T& d,const T& v) { 
      d=v; 
    } 
    static void construct(T* b, T* e, T* s) { 
      memcpy(b,s,(e-b)*sizeof(T)); 
    } 
    static void assign(T* b, T* e, const T* s) {
      memcpy(b,s,(e-b)*sizeof(T));
    }
    static void assign(T* b, T* e, T* s) {
      memcpy(b,s,(e-b)*sizeof(T));
    }
    template<typename Iter> inline static void construct_iter(T* b,T* e,Iter s) { 
      for (;b!=e;++s) 
 	*b++=*s; 
    } 
    
    template<typename Iter> inline static void assign(T* b,T* e,Iter s) { 
      for (;b!=e;++s) 
 	*b++=*s; 
    } 
    static void destroy(T*,T*) {}
    static void destroy(T&) {}
  };

  template<class T, class AllocatorClass > 
    class MemoryFiller<T,true,AllocatorClass> {

  private:
    typedef AllocatorClass alloc_t;
      //      static size_t claimed_,released_;
      
  public:
    typedef memop<T,true> memop_t;

      //      MemoryFiller() { claimed_=released_=0; }

    static T* get_uninit(size_t n) { 
      //      claimed_+=((n*sizeof(T))>>2);
      return n ? static_cast<T*>(alloc_t::allocate(n*sizeof(T))) : (T*)NULL; 
    } 
      static void release_uninit(T* p, size_t n) {
	//	released_+=((n*sizeof(T))>>2);
	alloc_t::deallocate(p,n*sizeof(T));
      }
      
    static void release(T* p,size_t n, size_t) { 
      release_uninit(p,n);
    }

    static T* get(int n, int ninit) {
      lcm_verify_getargs_(n,ninit);
      return get_uninit(n);
    }
    static T* get(int n, int ninit, const T& v) {
      lcm_verify_getargs_(n,ninit);
      T* datast=get_uninit(n);
      memop_t::construct(datast,datast+ninit,v);
      return datast;
    }
    static T* get(int n, int ninit, T* vp) {
      lcm_verify_getargs_(n,ninit);
      T* datast=get_uninit(n);
      memop_t::assign(datast,datast+ninit,vp);
      return datast;
    }
    template<typename Iter> static T* get_iter(int n, int ninit, Iter s) {
      lcm_verify_getargs_(n,ninit);
      T* datast=get_uninit(n);
      T* b=datast;
      T* e=datast+ninit;
      try {
	while (b!=e) { //Although operations on T will not throw, iterator is not guaranteed!
	  *b++=*s;
	  ++s;
	}
      } catch(...) {
	release(datast,0,n);
	throw;
      }
      return datast;
    }

      static void print(std::ostream& ostr) {
	//	ostr << "Words claimed (trivial): " << claimed_ << "  released: " << released_ << "  diff: " << (claimed_-released_) << '\n';
      }
  };
#endif

// Some functions take flags of the form mxflag::allownull etc.

namespace mxflag {
  enum oflags { 
    binary= 0x01, // Generate binary data files
    doublep=0x02, // Output double precision
    norowsep=0x08, // Suppress insertion of newlines between rows
    block = 0x10, // Complete row per output line
    rmat =0x20 //write in RMAT format (internal)
  };

  enum ioflags {
    allownull=0x01,
    confirmoverwrite=0x02,
    openbinary=0x04,
    defaultstdout=0x08
  };

  enum tensflag {
    NOTUSED=-1, //!< avoids warning with use of temporary
    maximum=0,
    all=1,
    none=2
  };

  enum tempflag {
    normal=0,
    temporary=-1,
    nondynamic=-2
  };
} //namespace mxflag


//Similar to STL functionals but in-place

 template<class Tout,class Tin> struct unary_function_ip {
   typedef Tout result_type;
   typedef Tin argument_type;
 };

 template<class Tout,class Tin1,class Tin2> struct binary_function_ip {
   typedef Tin1 first_argument_type;
   typedef Tin2 second_argument_type;
   typedef Tout result_type;
 };

template<class T> struct LCM_iterator_cat_ { typedef ::std::bidirectional_iterator_tag iterator_category; };
template<> struct LCM_iterator_cat_< ::std::forward_iterator_tag> { typedef ::std::forward_iterator_tag iterator_category; };

 template<typename T1,typename T2> struct doesnegate : public ::std::unary_function<T2,T1> {
   T1 operator()(const T2& a) const {
     return -a; }
 };
 template<typename T> struct doesnegate_ip {
   typedef T argument_type;
   void operator()(T& a) const {
     a=-a;
   }
 };

template<typename T1,typename T2> struct doesassign : public ::std::unary_function<T2,T1> {
  T1 operator()(const T2& a) const {
    return T1(a); }
};

  template<typename T> struct checksnonzero : public ::std::unary_function<T,bool> {
    bool operator()(const T& a) const { return (a!=T(0)); }
  };

  template<typename T> struct checksnonzero_tolerance : public ::std::unary_function<T,bool> {
    checksnonzero_tolerance(T tolv) : tol_(tolv) {
      if (tol_<0.0)
	throw InvalidParameter("checksnonzero_tolerance");
    }
    bool operator()(T a) const { return (fabs(a)>tol_); }
    T tol_;
  };

//use standard STL functionals for logical ops (logical_and etc.)
template<typename T1,typename T2,typename T3> struct doesand : public ::std::binary_function<T2,T3,T1> {
  inline T1 operator()(const T2& a, const T3& b) const {
    return a & b; }
};

template<typename T1,typename T2,typename T3> struct doesor : public ::std::binary_function<T2,T3,T1> {
  inline T1 operator()(const T2& a, const T3& b) const {
    return a | b; }
};


template<typename T1,typename T2> struct doespremultiply_ip : public unary_function_ip<T1,T2> {
  inline void operator()(T1& d, const T2& a) const {
    d=a*d;
  }
};

#define LCM_MAKE_IP(NAME,OP)\
template<typename T1,typename T2> struct does##NAME##_ip : public unary_function_ip<T1,T2> {\
  inline void operator()(T1& d, const T2& a) const { d OP##= a; }\
};

  LCM_MAKE_IP(divide,/)
    LCM_MAKE_IP(add,+)
    LCM_MAKE_IP(subtract,-)
    LCM_MAKE_IP(multiply,*)
    LCM_MAKE_IP(and,&)
    LCM_MAKE_IP(or,|)
    LCM_MAKE_IP(xor,^)

//N.B. special definition of *= bool: *= false replaces value with '0'
template<typename T> struct doesmultiply_ip<T,bool> : public unary_function_ip<T,bool> {
  inline void operator()(T& d, bool a) const {
    if (!a)
      d=T(0);
  }
};

template<typename T1,typename T2,typename T3> struct doesadd_sd : public binary_function_ip<T1,T2,T3> {
  inline void operator()(T1& d, const T2& a, const T3& b) const {
    d=a+b;
  }
};

template<typename T1,typename T2,typename T3> struct doessubtract_sd : public binary_function_ip<T1,T2,T3> {
  inline void operator()(T1& d,const T2& a, const T3& b) const {
    d=a-b;
  }
};

template<typename T1,typename T2,typename T3> struct doesmultiply_sd : public binary_function_ip<T1,T2,T3> {
  inline void operator()(T1& d,const T2& a, const T3& b) const {
    d=a*b;
  }
};

  template<typename T1,typename T2,typename T3> struct doesmultiply_sd< Matrix<T1>, Matrix<T2>,Matrix<T3> > : public binary_function_ip< Matrix<T1>,Matrix<T2>,Matrix<T3> > {
    inline void operator()(Matrix<T1>& d,const Matrix<T2>& a, const Matrix<T3>& b) const {
      multiply(d,a,b);
    }
  };

//Not important, but necessary to avoid ambiguities
template<> struct doesmultiply_sd<bool,bool,bool> : public binary_function_ip<bool,bool,bool> {
  inline void operator()(bool& d,bool a, bool b) const {
    d=a & b;
  }
};
template<typename T> struct doesmultiply_sd<T,T,bool> : public binary_function_ip<T,T,bool> {
  static const T zero=0;
  inline void operator()(T& d,const T& a, bool b) const {
    d=b ? a : zero;
  }
};
template<typename T> struct doesmultiply_sd<T,bool,T> : public binary_function_ip<T,bool,T> {
  static const T zero=0;
  inline void operator()(T& d,bool b,const T& a) const {
    d=b ? a : zero;
  }
};

template<typename T1,typename T2,typename T3> struct doesdivide_sd : public binary_function_ip<T1,T2,T3> {
  inline void operator()(T1& d,const T2& a, const T3& b) const {
    d=a/b;
  }
};

template<typename T1,typename T2,typename T3> struct doesmla : public binary_function_ip<T1,T2,T3> {
  inline void operator()(T1& d, const T2& a, const T3& b) const {
    d+=a*b;
  }
};

template<class Tagd,class Taga,class Tagb> struct iters3_ {
  template<class Iter_d,class Iter_a,class Iter_b,class Function> static inline void apply2(Iter_d dstart, const Iter_d& dend, Function mop,Iter_a astart, Iter_b bstart)
  {
    if (dstart==dend)
      return;
    *dstart=mop(*astart,*bstart);
    //prefer to use prefix ++ since generally faster for non-trivial iterators
    while ((++dstart)!=dend)
      *dstart=mop(*(++astart),*(++bstart));
  }
  template<class Iter_d,class Iter_a,class Iter_b,class Function> static inline void apply_ip2(Function mop,Iter_d dstart, const Iter_d& dend, Iter_a astart, Iter_b bstart)
  {
    if (dstart==dend)
      return;
    mop(*dstart,*astart,*bstart);
    //prefer to use prefix ++ since generally faster for non-trivial iterators
    while ((++dstart)!=dend)
      mop(*dstart,*(++astart),*(++bstart));
  }
};

template<class Tag1,class Tag2> struct iters2_ {
  template<class Iter_a, class Iter_b,class Function> static inline void applyip2(Function mop,Iter_a astart, const Iter_a& aend,Iter_b bstart)
  {
    if (astart==aend)
      return;
    mop(*astart,*bstart);
    while ((++astart)!=aend)
      mop(*astart,*(++bstart));
  }

  template<class Iter_d, class T, class Iter_b, class Function> static inline void apply2nd(Iter_d dstart, const Iter_d& dend, Function mop,const T& a, Iter_b bstart)
  {
    if (dstart==dend)
      return;
    *dstart=mop(a,*bstart);
    while ((++dstart)!=dend)
      *dstart=mop(a,*(++bstart));
  }

  template<class Iter_d, class Iter_a, class Tb, class Function> static inline void apply1st(Iter_d dstart, const Iter_d& dend, Function mop,Iter_a astart,Tb b)
  {
    if (dstart==dend)
      return;
    *dstart=mop(*astart,b);
    while ((++dstart)!=dend)
      *dstart=mop(*(++astart),b);
  }

  template<class Iter_d, class T, class Iter_b, class Function> static inline void apply_ip1st(Function mop,Iter_d dstart, const Iter_d& dend, Iter_b bstart, T c)
  {
    if (dstart==dend)
      return;
    mop(*dstart,*bstart,c);
    while ((++dstart)!=dend)
      mop(*dstart,*(++bstart),c);
  }

  template<class Iter_d, class T, class Iter_b, class Function> static inline void apply_ip2nd(Function mop,Iter_d dstart, const Iter_d& dend, const T& c, Iter_b bstart)
  {
    if (dstart==dend)
      return;
    mop(*dstart,c,*bstart);
    while ((++dstart)!=dend) 
      mop(*dstart,c,*(++bstart));
  }

  template<class Iter_d, class Iter_b, class Function> static inline void apply(Iter_d dstart, const Iter_d& dend, Function mop, Iter_b bstart)
  {
    if (dstart==dend)
      return;
    *dstart=mop(*bstart);
    while ((++dstart)!=dend)
      *dstart=mop(*(++bstart));
  }
};

template<class Tag> struct iters1_ {
  template<class Iter,class Function> inline static void applyip(Function mop,Iter astart, const Iter& aend)
  {
    if (astart==aend)
      return;
    mop(*astart);
    while ((++astart)!=aend)
      mop(*astart);
  }
  template<class Iter,class T> inline static void assign(Iter astart, const Iter& aend,const T v)
  {
    if (astart==aend)
      return;
    *astart=v;
    while ((++astart)!=aend)
      *astart=v;
  }
  template<class Iter, class T, class Function> inline static void apply_ip1st(Function mop,Iter astart,const Iter& aend,const T& b)
  {
    if (astart==aend)
      return;
    mop(*astart,b);
    while ((++astart)!=aend)
      mop(*astart,b);
  }
};

#ifndef LCM_DISABLE_UNROLL

template<> struct iters3_< ::std::random_access_iterator_tag,::std::random_access_iterator_tag,::std::random_access_iterator_tag> {
  template<class Iter_d, class Iter_a, class Iter_b,class Function> static inline void apply2(Iter_d dstart, const Iter_d& dend, Function mop, Iter_a astart, Iter_b bstart)
  {
    const ptrdiff_t n=dend-dstart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      *dstart=mop(*astart,*bstart); 
      dstart[1]=mop(astart[1],bstart[1]); 
      dstart[2]=mop(astart[2],bstart[2]); 
      dstart+=3;
      astart+=3;
      bstart+=3;
      break;
    case 2:
      *dstart=mop(*astart,*bstart); 
      dstart[1]=mop(astart[1],bstart[1]); 
      dstart+=2;
      astart+=2;
      bstart+=2;
      break;
    case 1:
      *dstart++=mop(*astart++,*bstart++); 
      break;
    }
    while (dstart!=dend) {
      *dstart=mop(*astart,*bstart); 
      dstart[1]=mop(astart[1],bstart[1]); 
      dstart[2]=mop(astart[2],bstart[2]); 
      dstart[3]=mop(astart[3],bstart[3]); 
      dstart+=4;
      astart+=4;
      bstart+=4;
    }
  }

  template<class Iter_d, class Iter_a, class Iter_b,class Function> static inline void apply_ip2(Function mop,Iter_d LCM_RESTRICT dstart, const Iter_d& LCM_RESTRICT dend,Iter_a LCM_RESTRICT astart, Iter_b LCM_RESTRICT bstart)
  {
    const ptrdiff_t n=dend-dstart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      mop(*dstart,*astart,*bstart); 
      mop(dstart[1],astart[1],bstart[1]); 
      mop(dstart[2],astart[2],bstart[2]); 
      dstart+=3;
      astart+=3;
      bstart+=3;
      break;
    case 2:
      mop(*dstart,*astart,*bstart); 
      mop(dstart[1],astart[1],bstart[1]); 
      dstart+=2;
      astart+=2;
      bstart+=2;
      break;
    case 1:
      mop(*dstart++,*astart++,*bstart++); 
      break;
    }
    while (dstart!=dend) {
      mop(*dstart,*astart,*bstart);
      mop(dstart[1],astart[1],bstart[1]);
      mop(dstart[2],astart[2],bstart[2]);
      mop(dstart[3],astart[3],bstart[3]);
      dstart+=4;
      astart+=4;
      bstart+=4;
    }
  }
};

template<> struct iters2_< ::std::random_access_iterator_tag,::std::random_access_iterator_tag> {
  template<class Iter_a, class Iter_b,class Function> static inline void applyip2(Function mop,Iter_a astart, const Iter_a LCM_RESTRICT & aend,Iter_b LCM_RESTRICT bstart)
  {
    const ptrdiff_t n=aend-astart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      mop(*astart,*bstart); 
      mop(astart[1],bstart[1]); 
      mop(astart[2],bstart[2]); 
      astart+=3;
      bstart+=3;
      break;
    case 2:
      mop(*astart,*bstart); 
      mop(astart[1],bstart[1]); 
      astart+=2;
      bstart+=2;
      break;
    case 1:
      mop(*astart,*bstart); 
      ++astart;
      ++bstart;
      break;
    }
    while (astart!=aend) {
      mop(*astart,*bstart);
      mop(astart[1],bstart[1]);
      mop(astart[2],bstart[2]);
      mop(astart[3],bstart[3]);
      astart+=4;
      bstart+=4;
    }
  }

  template<class Iter_a, class Iter_b,class T,class Function> static inline void apply_ip1st(Function mop,Iter_a astart, const Iter_a& aend,Iter_b bstart, T c)
  {
    const ptrdiff_t n=aend-astart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      mop(*astart,*bstart,c); 
      mop(astart[1],bstart[1],c); 
      mop(astart[2],bstart[2],c); 
      astart+=3;
      bstart+=3;
      break;
    case 2:
      mop(*astart,*bstart,c); 
      mop(astart[1],bstart[1],c); 
      astart+=2;
      bstart+=2;
      break;
    case 1:
      mop(*astart,*bstart,c); 
      ++astart;
      ++bstart;
      break;
    }
    while (astart!=aend) {
      mop(*astart,*bstart,c);
      mop(astart[1],bstart[1],c);
      mop(astart[2],bstart[2],c);
      mop(astart[3],bstart[3],c);
      astart+=4;
      bstart+=4;
    }
  }

  template<class Iter_a, class Iter_b,class T,class Function> static inline void apply_ip2nd(Function mop,Iter_a astart, const Iter_a& aend, const T& c, Iter_b bstart)
  {
    const ptrdiff_t n=aend-astart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      mop(*astart,c,*bstart); 
      mop(astart[1],c,bstart[1]); 
      mop(astart[2],c,bstart[2]); 
      astart+=3;
      bstart+=3;
      break;
    case 2:
      mop(*astart,c,*bstart); 
      mop(astart[1],c,bstart[1]); 
      astart+=2;
      bstart+=2;
      break;
    case 1:
      mop(*astart,c,*bstart); 
      ++astart;
      ++bstart;
      break;
    }
    while (astart!=aend) {
      mop(*astart,c,*bstart);
      mop(astart[1],c,bstart[1]);
      mop(astart[2],c,bstart[2]);
      mop(astart[3],c,bstart[3]);
      astart+=4;
      bstart+=4;
    }
  }

  template<class Iter_d, class Iter_b, class Function> static inline void apply(Iter_d dstart, const Iter_d& dend, Function mop, Iter_b bstart)
  {
    const ptrdiff_t n=dend-dstart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      *dstart=mop(*bstart); 
      dstart[1]=mop(bstart[1]); 
      dstart[2]=mop(bstart[2]); 
      dstart+=3;
      bstart+=3;
      break;
    case 2:
      *dstart=mop(*bstart); 
      dstart[1]=mop(bstart[1]); 
      dstart+=2;
      bstart+=2;
      break;
    case 1:
      *dstart=mop(*bstart); 
      ++dstart;
      ++bstart;
      break;
    }
    while (dstart!=dend) {
      *dstart=mop(*bstart); 
      dstart[1]=mop(bstart[1]); 
      dstart[2]=mop(bstart[2]); 
      dstart[3]=mop(bstart[3]); 
      dstart+=4;
      bstart+=4;
    }
  }

  template<class Iter_d, class T, class Iter_b, class Function> static inline void apply2nd(Iter_d dstart, const Iter_d& dend, Function mop,const T& a, Iter_b bstart)
  {
    const ptrdiff_t n=dend-dstart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      *dstart=mop(a,*bstart); 
      dstart[1]=mop(a,bstart[1]); 
      dstart[2]=mop(a,bstart[2]); 
      dstart+=3;
      bstart+=3;
      break;
    case 2:
      *dstart=mop(a,*bstart); 
      dstart[1]=mop(a,bstart[1]); 
      dstart+=2;
      bstart+=2;
      break;
    case 1:
      *dstart=mop(a,*bstart); 
      ++dstart;
      ++bstart;
      break;
    }
    while (dstart!=dend) {
      *dstart=mop(a,*bstart); 
      dstart[1]=mop(a,bstart[1]); 
      dstart[2]=mop(a,bstart[2]); 
      dstart[3]=mop(a,bstart[3]); 
      dstart+=4;
      bstart+=4;
    }
  }

  template<class Iter_d, class Iter_b, class T,class Function> static inline void apply1st(Iter_d dstart, const Iter_d& dend, Function mop,Iter_b bstart, T a)
  {
    const ptrdiff_t n=dend-dstart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      *dstart=mop(*bstart,a); 
      dstart[1]=mop(bstart[1],a); 
      dstart[2]=mop(bstart[2],a); 
      dstart+=3;
      bstart+=3;
      break;
    case 2:
      *dstart=mop(*bstart,a); 
      dstart[1]=mop(bstart[1],a); 
      dstart+=2;
      bstart+=2;
      break;
    case 1:
      *dstart=mop(*bstart,a); 
      ++dstart;
      ++bstart;
      break;
    }
    while (dstart!=dend) {
      *dstart=mop(*bstart,a); 
      dstart[1]=mop(bstart[1],a); 
      dstart[2]=mop(bstart[2],a); 
      dstart[3]=mop(bstart[3],a); 
      dstart+=4;
      bstart+=4;
    }
  }

};

template<> struct iters1_< ::std::random_access_iterator_tag> {
  template<class Iter,class Function> inline static void applyip(Function mop,Iter astart, const Iter& aend)
  {
    const ptrdiff_t n=aend-astart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      mop(*astart); 
      mop(astart[1]); 
      mop(astart[2]); 
      astart+=3;
      break;
    case 2:
      mop(*astart); 
      mop(astart[1]); 
      astart+=2;
      break;
    case 1:
      mop(*astart); 
      ++astart;
      break;
    }
    while (astart!=aend) {
      mop(*astart);
      mop(astart[1]);
      mop(astart[2]);
      mop(astart[3]);
      astart+=4;
    }
  }
  template<class Iter,class T,class Function> inline static void apply_ip1st(Function mop,Iter astart, const Iter& aend,T b)
  {
    const ptrdiff_t n=aend-astart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      mop(*astart,b); 
      mop(astart[1],b); 
      mop(astart[2],b); 
      astart+=3;
      break;
    case 2:
      mop(*astart,b); 
      mop(astart[1],b); 
      astart+=2;
      break;
    case 1:
      mop(*astart,b); 
      ++astart;
      break;
    }
    while (astart!=aend) {
      mop(*astart,b);
      mop(astart[1],b);
      mop(astart[2],b);
      mop(astart[3],b);
      astart+=4;
    }
  }
  template<class Iter,class T> inline static void assign(Iter astart, const Iter& aend,const T v)
  {
    const ptrdiff_t n=aend-astart; //must assume that ptrdiff_t is a simple type
    switch (n & 3) {
    case 3:
      *astart=v;
      astart[1]=v;
      astart[2]=v;
      astart+=3;
      break;
    case 2:
      *astart=v;
      astart[1]=v;
      astart+=2;
      break;
    case 1:
      *astart=v;
      ++astart;
      break;
    }
    while (astart!=aend) {
      (*astart)=v;
      astart[1]=v;
      astart[2]=v;
      astart[3]=v;
      astart+=4;
    }
  }
};
#endif

//rules for in-place operations

//In-place, unary
template<size_t dim> struct Apply_ip {
  template<class Function,class T> inline static void func(Function func,T& a) { 
    iters1_<LCM_ITER_CAT(typename T::iterator)>::applyip(func,a.begin(),a.end());
  }
};
template<> struct Apply_ip<0> {
  template<class Function,class T> inline static void func(Function func,T& a) {
    func(a); }
};

template<class T,class Function> inline void apply_ip(Function mop,T& a)
{ Apply_ip<LCM_DIM(T)>::func(mop,a); }


//In-place, binary

template<class Function,class T1,class T2> T1& diagop_(Function mop,T1& a,const T2& b)
{
  if (!issquare(a))
    throw NotSquare("diagop");
  const size_t n=a.rows();
  if (n!=b.size())
    throw Mismatch("diagop");
  typename T2::const_iterator start=b.begin();
  for (size_t i=0;i<n;i++)
    mop(a(i,i),*start++); //postfix ++ may be inefficient, but these are not time-critical ops
  return a;
}


template<size_t dim1, size_t dim2> struct Apply_ip2 {};
template<size_t N> struct Apply_ip2<N,N> {
  template<class T1,class T2,class Function> inline static void func(Function func,T1& a, const T2& b) {
    if (!arematching(a,b))
      throw Mismatch("Apply_ip2");
    iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator)>::applyip2(func,a.begin(),a.end(),b.begin());
  }
  template<class T1,class T2,class Function> inline static T1& func_allownull(Function func,T1& a, const T2& b) {
    if (a.empty())
      a=b;
    else {
      if (!arematching(a,b))
	throw Mismatch("func_allownull<N,N>");
      iters2_<LCM_ITER_CAT(typename T1::iterator),
	LCM_ITER_CAT(typename T2::const_iterator)>::applyip2(func,a.begin(),a.end(),b.begin());
    }
    return a;
  }
  template<class T1,class T2> inline static T1& subtract(T1& a, const T2& b) {
    if (a.empty()) {
      duplicate_structure(a,b);
      iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator)>::apply(a.begin(),a.end(),doesnegate<LCM_VAL(T1),LCM_VAL(T2)>(),b.begin());
    }
    else {
      if (!arematching(a,b))
	throw Mismatch("subtract_ip<N,N>");
      iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator)>::applyip2(doessubtract_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a.begin(),a.end(),b.begin());
    }
    return a;
  }
};
template<size_t N> struct Apply_ip2<N,0> {
  template<class T1,class T2,class Function> LCM_INLINE static void func(Function func,T1& a, const T2& b) {
    iters1_<LCM_ITER_CAT(typename T1::iterator)>::apply_ip1st(func,a.begin(),a.end(),b);
  }
  template<class T1,class T2,class Function> LCM_INLINE static T1& func_allownull(Function func,T1& a, const T2& b) {
    iters1_<LCM_ITER_CAT(typename T1::iterator)>::apply_ip1st(func,a.begin(),a.end(),b);
    return a;
  }
  template<class T1,class T2> LCM_INLINE static T1& subtract(T1& a, const T2& b) {
    iters1_<LCM_ITER_CAT(typename T1::iterator)>::apply_ip1st(doessubtract_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a.begin(),a.end(),b);
    return a;
  }
};

template<> struct Apply_ip2<2,1> {
  template<class T1,class T2,class Function> inline static T1& func_allownull(Function func,T1& a, const T2& b) {
    if (a.empty()) {
      const size_t n=b.size();
      a.create(n,n,(LCM_VAL(T1))0);
      typename T2::const_iterator start=b.begin();
      for (size_t i=0;i<n;i++)
	a(i,i)=*start++;
      return a;
    }
    return diagop_(func,a,b);
  }

  template<class T1,class T2> LCM_INLINE static T1& subtract(T1& a, const T2& b) {
    if (a.empty()) {
      const size_t n=b.size();
      a.create(n,n,(LCM_VAL(T1))0);
      typename T2::const_iterator start=b.begin();
      for (size_t i=0;i<n;i++) 
	a(i,i)=-(*start++);
      return a;
    }
    else
      return diagop_(doessubtract_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b);
  }
};

template<> struct Apply_ip2<0,0> {
  template<class T1,class T2,class Function> inline static void func(Function func,T1& a, const T2& b) {
    func(a,b); }

  template<class T1,class T2,class Function> inline static T1& func_allownull(Function func,T1& a, const T2& b)
  { return func(a,b); }

  template<class T1,class T2> inline static T1& subtract(T1& a, const T2& b)
  { return doessubtract_ip<T1,T2>::operator()(a,b); }
};

template<class T1,class T2,class Function> inline T1& apply_ip2(Function mop,T1& a,const T2& b)
{ Apply_ip2<LCM_DIM(T1),LCM_DIM(T2)>::func(mop,a,b); return a; }

template<size_t, size_t, size_t> struct Apply_ip3 {};
template<size_t N> struct Apply_ip3<N,N,N> {
  template<class T1,class T2,class T3,class Function> inline static void func(Function func,T1& a, const T2& b,const T3& c) {
    if (!arematching(a,b) || !arematching(a,c))
      throw Mismatch("Apply_ip3");
    iters3_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator),LCM_ITER_CAT(typename T3::const_iterator)>::
      apply_ip2(func,a.begin(),a.end(),b.begin(),c.begin());
  }
};

template<size_t N> struct Apply_ip3<N,N,0> {
  template<class T1,class T2,class T3,class Function> LCM_INLINE static void func(Function func,T1& a, const T2& b,const T3& c) {
    if (!arematching(a,b))
      throw Mismatch("Apply_ip3");
    iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator)>::
      apply_ip1st(func,a.begin(),a.end(),b.begin(),c);
  }
};

template<size_t N> struct Apply_ip3<N,0,N> {
  template<class T1,class T2,class T3,class Function> LCM_INLINE static void func(Function func,T1& a, const T2& b,const T3& c) {
    if (!arematching(a,c))
      throw Mismatch("Apply_ip3");
    iters2_< LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator) >::
      apply_ip2nd(func,a.begin(),a.end(),b,c.begin());
  }
};

template<> struct Apply_ip3<0,0,0> {
  template<class T1,class T2,class T3,class Function> inline static void func(Function func,T1& a, const T2& b,const T3& c)
  { func(a,b,c); }
};

template<class T1,class T2,class T3,class Function> inline void apply_ip3(Function mop,T1& a,const T2& b,const T3& c)
{ Apply_ip3<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(mop,a,b,c); }

template<class T1,class T2> inline T1& add_ip(T1& a,const T2& b)
{ return Apply_ip2<LCM_DIM(T1),LCM_DIM(T2)>::func_allownull(doesadd_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b); }

template<class T1,class T2> inline T1& subtract_ip(T1& a,const T2& b)
{ return Apply_ip2<LCM_DIM(T1),LCM_DIM(T2)>::subtract(a,b); }

//rules for in-place multiply
template<size_t dim1, size_t dim2> struct _Multiply_ip {
  template<class T1,class T2> inline static T1& post(T1& a, const T2& b) {
    Apply_ip2<dim1,dim2>::func(doesmultiply_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b);
    return a;
  }
  template<class T1,class T2> inline static T1& pre(T1& a, const T2& b) {
    Apply_ip2<dim1,dim2>::func(doespremultiply_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b);
    return a;
  }
};

template<> struct _Multiply_ip<2,2> {
  template<class T1,class T2> LCM_INLINE static T1& post(T1& a, const T2& b) {
    if (a.empty())
      return (a=b);
    T1 tmp;
    multiply(tmp,a,b);
    a.swap(tmp);
    return a;
  }

  template<class T1,class T2> LCM_INLINE static T1& pre(T1& a, const T2& b) {
    if (a.empty())
      return (a=b);
    T1 tmp;
    multiply(tmp,b,a);
    a.swap(tmp);
    return a;
  }

};

template<> struct _Multiply_ip<2,1> {
  template<class T1,class T2> inline static T1& post(T1& a, const T2& b)
  {
    if (a.empty())
      full(a,b);
    else {
      const size_t n=b.size();
      if (n!=a.cols())
	throw Mismatch("multiply_ip<2,1>");
      
      typename T1::iterator dcur=a.begin();
      
      for (size_t i=a.rows();i--;) {
	for (size_t j=0;j<n;j++) {
	  (*dcur)*=b(j);
	  ++dcur;
	}
      }
    }
    return a;
  }
  template<class T1,class T2> static T1& pre(T1& a, const T2& b)
  {
    if (a.empty())
      full(a,b);
    else {
      const size_t m=b.size();
      if (m!=a.rows())
	throw Mismatch("multiply_ip<2,1>");
      const size_t n=a.cols();
      
      typename T1::iterator dcur=a.begin();
      
      for (size_t i=0;i<m;i++) {
	const LCM_VAL(T2)& bi=b(i);
	for (size_t j=n;j--;) {
	  (*dcur)*=bi;
	  ++dcur;
	}
      }
    }
    return a;
  }
};


template<> struct _Multiply_ip<1,1> {
  template<class T1,class T2> LCM_INLINE static void copy_if_match(T1& a, const T2& b, Bool2Type<true>) {
    assign(a,b);
  }
  template<class T1,class T2> LCM_INLINE static void copy_if_match(T1& a, const T2&, Bool2Type<false>) {
    throw Undefined("pre-multiply");
  }

  template<class T1,class T2> LCM_INLINE static T1& post(T1& a, const T2& b) {
    if (a.size())
      return apply_ip2(doesmultiply_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b);
    copy_if_match(a,b,Bool2Type< LCM_DIM( LCM_VAL(T1) ) == LCM_DIM( LCM_VAL(T2) ) >());
    return a;
  }
  template<class T1,class T2> LCM_INLINE static T1& pre(T1& a, const T2& b) {
    if (a.size())
      return apply_ip2(doespremultiply_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b);
    copy_if_match(a,b,Bool2Type< LCM_DIM( LCM_VAL(T1) ) == LCM_DIM( LCM_VAL(T2) ) >());
    return a;
  }
};

template<class T1,class T2> inline T1& multiply_ip(T1& a,const T2& b)
{ return _Multiply_ip<LCM_DIM(T1),LCM_DIM(T2)>::post(a,b); }

template<class T1,class T2> inline T1& premultiply_ip(T1& a,const T2& b)
{ return _Multiply_ip<LCM_DIM(T1),LCM_DIM(T2)>::pre(a,b); }

#define LCM_DECLARE_IP(OP)\
template<class T1,class T2> inline T1& OP##_ip(T1& a, const T2& b)\
{ return apply_ip2(does##OP##_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b); }

#define LCM_DECLARE_IP_ALLOWNULL(OP)\
template<class T1,class T2> inline T1& OP##_ip(T1& a, const T2& b)\
{ return Apply_ip2<LCM_DIM(T1),LCM_DIM(T2)>::func_allownull(does##OP##_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b); }

LCM_DECLARE_IP(divide)
  LCM_DECLARE_IP_ALLOWNULL(or)
  LCM_DECLARE_IP(and)
  LCM_DECLARE_IP_ALLOWNULL(xor)

// template<class T1,class T2> inline T1& divide_ip(T1& a, const T2& b)
// { return apply_ip2(doesdivide_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b); }

template<class T1,class T2> inline void emultiply_ip(T1& a, const T2& b)
{ apply_ip2(doesmultiply_ip<LCM_VAL(T1),LCM_VAL(T2)>(),a,b); }

template<size_t dim1, size_t dim2, size_t dim3> struct Apply2_ {};

template<size_t N> struct Apply2_<N,N,N> {
  template<class T1,class T2,class T3,class Function> LCM_INLINE static void func(T1& d, Function mop, const T2& a, const T3& b) { 
    if (a.empty())
      throw Undefined("Apply2<N,N,N>");
    if (!arematching(a,b))
      throw Mismatch("Apply2<N,N,N>");
    duplicate_structure(d,a);
    iters3_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator),LCM_ITER_CAT(typename T3::const_iterator)>::apply2(d.begin(),d.end(),mop,a.begin(),b.begin());
  }
};
template<size_t N> struct Apply2_<N,N,0> {
  template<class T1,class T2,class T3,class Function> LCM_INLINE static void func(T1& d, Function mop, const T2& a, const T3& b) { 
    duplicate_structure(d,a);
    iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator)>::
      apply1st(d.begin(),d.end(),mop,a.begin(),b);
  }
};
template<size_t N> struct Apply2_<N,0,N> {
  template<class T1,class T2,class T3,class Function> LCM_INLINE static void func(T1& d, Function mop, const T2& a, const T3& b) { 
    duplicate_structure(d,b);
    iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T3::const_iterator)>::
      apply2nd(d.begin(),d.end(),mop,a,b.begin());
  }
};
template<> struct Apply2_<0,0,0> {
  template<class T1,class T2,class T3,class Function> inline static void func(T1& d, Function mop, const T2& a, const T3& b)
  { mop(d,a,b); }
};

template<class T1,class T2,class T3,class Function> inline void apply(T1& d,Function mop,const T2& a,const T3& b)
{ Apply2_<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(d,mop,a,b); }

template<size_t dim1,size_t dim2> struct Apply_ {};

template<size_t N> struct Apply_<N,N> {
  template<class T1,class T2,class Function> inline static void func(T1& a, Function func, const T2& b) {
//     if (b.empty())
//       a.clear();
//     else {
      duplicate_structure(a,b);
      iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator)>::apply(a.begin(),a.end(),func,b.begin());
      //    }
  }
  template<class T1,class T2> inline static void assign(T1& a,const T2& b) {
//     if (b.empty())
//       a.clear();
//     else {
    if (issame(a,b))
      throw ArgumentClash("assign");
    duplicate_structure(a,b);
    iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator)>::apply(a.begin(),a.end(),doesassign<LCM_VAL(T1),LCM_VAL(T2)>(),b.begin());
      //    }
  }    
};
template<size_t N> struct Apply_<N,0> {
  template<class T1,class T2> LCM_INLINE static void assign(T1& a, const T2& b) {
    if (a.empty())
      throw Undefined("assign<N,0>");
    iters1_<LCM_ITER_CAT(typename T1::iterator)>::assign(a.begin(),a.end(),b);
  }
};
template<> struct Apply_<0,0> {
  template<class T1,class T2,class Function> inline static void func(T1& d, Function mop, const T2& a) {
    d=mop(a);
  }

  template<class T1,class T2> inline static void assign(T1& a,const T2& b) {
    a=b;
  }
};

template<class T1,class T2,class Function> inline void apply(T1& d,Function mop,const T2& a)
{ Apply_<LCM_DIM(T1),LCM_DIM(T2)>::func(d,mop,a); }

template<typename T1,typename T2> T1& assign(T1& a,const T2& b) {
  Apply_<LCM_DIM(T1),LCM_DIM(T2)>::assign(a,b);
  return a;
}
template<class T1,class T2> inline void negate(T1& d, const T2& a) {
  Apply_<LCM_DIM(T1),LCM_DIM(T2)>::func(d,doesnegate<LCM_VAL(T1),LCM_VAL(T2)>(),a);
}

template<class T> inline void negate_ip(T& a) {
  Apply_ip<LCM_DIM(T)>::func(doesnegate_ip<LCM_VAL(T)>(),a); 
 }

template<size_t,size_t,size_t> struct Multiply_ {};

template<size_t N> struct Multiply_<N,0,N> {
  template<class T1,class T2,class T3> LCM_INLINE static void func(T1& d, const T2& a, const T3& b) {
    d=b;
    premultiply_ip(d,a);
  }
};

template<size_t N> struct Multiply_<N,N,0> {
  template<class T1,class T2,class T3> LCM_INLINE static void func(T1& d, const T2& a, const T3& b) {
    d=a;
    multiply_ip(d,b);
  }
};

template<> struct Multiply_<2,2,2> {
  template<class T1,class T2,class T3> static void func(T1& dest,const T2& a,const T3& b)
  {
    if (issame(&dest,&a) || issame(&dest,&b))
      throw ArgumentClash("Multiply<2,2,2>");
    
    const size_t ca=a.cols();
    const size_t cb=b.cols();
    const size_t ra=a.rows();
    
    if (ca!=b.rows())
      throw Mismatch("multiply");
    
    dest.create(ra,cb);
    
    typedef LCM_VAL(T1) dest_type;
    //typedef LCM_VAL(T2) a_type;
    //typedef LCM_VAL(T3) b_type;
        
    dest_type sum;
    for (size_t i=0;i<ra;i++) {
      for (size_t j=0;j<cb;j++) {
	sum=a(i,0)*b(0,j);
	for (size_t k=1;k<ca;k++)
	  mla(sum,a(i,k),b(k,j));
	dest(i,j)=sum;
      }
    }
  }
};

template<> struct Multiply_<2,1,1> {
  template<class T1,class T2,class T3> static void func(T1& d, const T2& a, const T3& b) {
    const size_t na=a.size();
    const size_t nb=b.size();
    d.create(na,nb);
    
    typedef LCM_VAL(T2) a_type;
    
    for (size_t ia=na;ia--;) {
      const a_type& aval=a(ia);
      for (size_t ib=nb;ib--;)
	d(ia,ib)=aval*b(ib);
    }
  }
};

template<> struct Multiply_<2,1,2> {
  template<class T1,class T2,class T3> static void func(T1& dest,const T2& b,const T3& a)
  {
    if (issame(&dest,&a))
      _Multiply_ip<2,1>::pre(dest,b);
    else {
      const size_t m=b.size();
      const size_t n=a.cols();
      if (m!=a.rows())
	throw Mismatch("Multiply<2,1,2>");
      dest.create(m,n);
      
      typedef LCM_VAL(T2) b_type;
      typename T1::iterator dcur=dest.begin();
      typename T3::const_iterator acur=a.begin();
      
      for (size_t i=0;i<m;i++) {
	const b_type& bi=b(i);
	for (size_t j=n;j--;) {
	  *dcur=bi*(*acur);
	  ++dcur;
	  ++acur;
	}
      }
    }
  }
};

template<> struct Multiply_<2,2,1> {
  template<class T1,class T2,class T3> static void func(T1& dest,const T2& a,const T3& b)
    {
      if (issame(&dest,&a))
	_Multiply_ip<2,1>::post(dest,b);
      else {
	const size_t n=b.size();
	if (n!=a.cols()) 
	  throw Mismatch("multiply<2,2,1>");
	const size_t m=a.rows();
	
	dest.create(m,n);
	
	typename T1::iterator dcur=dest.begin();
	typename T2::const_iterator acur=a.begin();
	
	for (size_t i=m;i--;) {
	  for (size_t j=0;j<n;j++) {
	    *dcur=(*acur)*b(j);
	    ++dcur;
	    ++acur;
	  }
	}
      }
    }
};

  template<> struct Multiply_<1,1,2> {
    template<class Td,class Tb,class Ta> static void func(Td& dest,const Tb& b,const Ta& a)
    {
      if (a.empty()) 
	throw Undefined("multiply<1,1,2>");
      if (issame(&dest,&b))
	throw ArgumentClash("multiply<1,1,2>");

      const size_t rs=a.rows();
      if (rs!=b.size())
	throw Mismatch("multiply<1,1,2>");
      size_t n=a.cols();
      dest.create(n);
      
      const typename Tb::value_type& b0=b(0);

      for (;n--;) {
	typename Td::value_type& d=dest(n);
	d=b0*a(0,n);
	for (size_t j=1;j<rs;j++)
	  mla(d,b(j),a(j,n));
      }
    }
  };

  template<> struct Multiply_<1,2,1> {
    template<class Td,class Ta,class Tb> static void func(Td& dest,const Ta& a,const Tb& b)
    {
      if (a.empty())
	throw Undefined("multiply<1,2,1>");
      if (issame(&dest,&b))
	throw ArgumentClash("multiply<1,2,1>");
      
      const size_t cs=a.cols();
      size_t n=a.rows();
      if (cs!=b.size())
	throw Mismatch("multiply<1,2,1>");
      dest.create(n);
      
      const typename Tb::value_type& b0=b(0);
      for (;n--;) {
	typename Td::value_type& d=dest(n);
	d=a(n,0)*b0;
	for (size_t m=1;m<cs;m++) 
	  mla(d,a(n,m),b(m));
      }
    }
  };

template<> struct Multiply_<1,1,1> {
  template<class T1,class T2,class T3> LCM_INLINE static void func(T1& d, const T2& a, const T3& b) {
    if (!arematching(a,b))
      throw Mismatch("Multiply<1,1,1>");
    duplicate_structure(d,a);
    iters3_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator),LCM_ITER_CAT(typename T3::const_iterator)>::apply_ip2(doesmultiply_sd<LCM_VAL(T1),LCM_VAL(T2),LCM_VAL(T3)>(),d.begin(),d.end(),a.begin(),b.begin());
    //    d=a;
    //multiply_ip(d,b);
  }
};

template<> struct Multiply_<0,0,0> {
  template<class T1,class T2,class T3> inline static void func(T1& d, const T2& a, const T3& b) {
    doesmultiply_sd<T1,T2,T3>()(d,a,b);
  }
};

 template<class T1,class T2,class T3> inline void multiply(T1& d,const T2& a,const T3& b) {
   Multiply_<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(d,a,b);
 }

template<size_t,size_t,size_t> struct Mla_ {
  template<class T1,class T2,class T3> static void func(T1&, const T2&, const T3&) {
    LCM_STATIC_ERROR( mla_incompatible_dimensionalities );
  }
};

template<size_t N> struct Mla_<N,0,N> {
  template<class T1,class T2,class T3> LCM_INLINE static void func(T1& d, const T2& a, const T3& b) {
    if (b.empty())
      throw Undefined("Mla<N,0,N>");
    if (d.empty())
      multiply(d,a,b);
    else {
      if (!arematching(d,b))
	throw Mismatch("Mla<N,0,N>");

      iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T3::const_iterator)>::
	apply_ip2nd(doesmla<LCM_VAL(T1),T2,LCM_VAL(T3)>(),d.begin(),d.end(),a,b.begin());
    }
  }
};

template<> struct Mla_<2,0,1> {
  template<class T1,class T2,class T3> static void func(T1& d, const T2& a, const T3& b) {
    size_t n=b.size();
    typename T3::const_iterator bstart=b.begin();
    if (d.empty()) {
      if (n==0)
	return;
      const LCM_VAL(T1) dest_zero(0);
      d.create(n,n,dest_zero);
      doesmultiply_sd<LCM_VAL(T1),T2,LCM_VAL(T3)> fobj; 
      for (size_t i=0;i<n;i++)
	fobj(d(i,i),a,*bstart++);
    }
    else {
      if (!issquare(d))
	throw NotSquare("mla");
      if (b.size()!=d.rows())
	throw Mismatch("Mla<2,0,1>");
      doesmla< LCM_VAL(T1), T2, LCM_VAL(T3)> fobj;
      for (size_t i=0;i<n;i++)
	fobj(d(i,i),a,*bstart++);
    }
  }
};

template<> struct Mla_<1,1,1> {
  template<class T1,class T2,class T3> LCM_INLINE static void func(T1& d, const T2& a, const T3& b) {
    if (d.size()==0)
      Multiply_<1,1,1>::func(d,a,b);
    else {
      if (!arematching(a,b) || !arematching(a,d))
	throw Mismatch("Mla<1,1,1>");
      iters3_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator),LCM_ITER_CAT(typename T3::const_iterator)>::apply_ip2(doesmla<LCM_VAL(T2),LCM_VAL(T3),LCM_VAL(T1)>(),d.begin(),d.end(),a.begin(),b.begin());
    }
  }
};

template<> struct Mla_<0,0,0> {
  template<class T1,class T2,class T3> inline static void func(T1& d, const T2& a, const T3& b) {
    doesmla<T1,T2,T3>()(d,a,b);
  }
};

template<class T1,class T2,class T3> inline void mla(T1& d,const T2& a,const T3& b) {
  Mla_<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(d,a,b);
}

template<class T,size_t N> struct Print_ {
  static void print(const T& a,std::ostream& ostr) {
    ostr << '[';
    for (size_t i=0;i<N-1;i++) 
      ostr << a.dimension(i) << 'x';
    ostr << a.dimension(N-1) << ']';
  }
};
template<class T> struct Print_<T,0> { 
  static void print(const T& a,std::ostream& ostr) { 
    ostr << a;
  }
};
template<class T> struct Print_<T,1> {
  static void print(const T& a, std::ostream& ostr) {
    const size_t n=a.size();
    if (n==0)
      ostr << "<empty>";
    else {
      ostr << "(";
      typename T::const_iterator iter=a.begin();
      for (size_t i=0;i<n;i++)
	ostr << (*iter++) << (i==n-1 ? ')' : ',');
    }
  }
};

  //! preferences for stream output
  struct ostream_controller {
    ostream_controller();
    int matrixprecision; //!< precision for matrix output (-1 use default, 0 spy)
    enum complexview_t { pair, //!< (real,imag)
			 withi, //!< R+Ii
			 compact //!< R+Ii omitting zero components
    };
    complexview_t complexview; //!< format for complex output    
    int timeprecision; //!< precision for time output
    double spytolerance;
    int indexbase;

    std::ostream& print(std::ostream&, const complex&, size_t =0) const;
    std::ostream& printtime(std::ostream&, double, size_t) const;
    bool complexfixedwidth() const { return (complexview!=compact); }
  };

  struct setmatrixprecision {
    explicit setmatrixprecision(int precv) : prec_(precv) {}
    int prec_;
  };

  std::ostream& operator<< (std::ostream&, const setmatrixprecision&);

  struct settimeprecision {
    explicit settimeprecision(int precv) : prec_(precv) {}
    int prec_;
  };

  std::ostream& operator<< (std::ostream&, const settimeprecision&);

void cmatrix_ostream_controller(std::ostream&, ostream_controller&);
ostream_controller& cmatrix_ostream_controller(std::ostream& =std::cout);

void spy(std::ostream&, const Matrix<double>&, double =1e-10);
void spy(std::ostream&, const Matrix<complex>&, double =1e-10);
void spy(std::ostream&, const Matrix<bool>&);


template<class T> bool tryspy_(const T&, std::ostream&, double) { return false; }
template<> bool tryspy_(const Matrix<double>&, std::ostream&, double tol);
template<> bool tryspy_(const Matrix<bool>&, std::ostream&, double);


template<class T> struct Print_<T,2> {

  static void print(const T& a,std::ostream& ostr) {
    if (a.empty()) {
      ostr << "<empty> [" << a.rows() << " x " << a.cols() << "]\n";
      return;
    }
    ostream_controller& ctrl(cmatrix_ostream_controller(ostr));    
    if ((ctrl.matrixprecision==0) && tryspy_(a,ostr,ctrl.spytolerance))
      return;

    // Float output is barely readable without fixed format
    //    const LCM_IOS::fmtflags oldflags=ostr.setf(LCM_IOS::floatfield,LCM_IOS::fixed);
    int width=(ctrl.matrixprecision>=0) ? ctrl.matrixprecision : ostr.precision(); //!< default to stream precision if not explicitly set
    int pad=ostr.width()-width;
    if (pad<1)
      pad=1;

    typename T::const_iterator iter=a.begin();
    for (size_t i=0;i<a.rows();i++) {
      for (size_t j=0;j<a.cols();j++) {
	ostr << std::setw(width) << (*iter);
	++iter;
	for (size_t k=pad;k--;)
	  ostr << ' ';
      }
      ostr << '\n';
    }
    //    ostr.flags(oldflags);
  }
};
template<class T> struct Print_<T,3> {
  static void print(const T& a,std::ostream& ostr) {
    if (a.empty()) {
      ostr << "<empty>\n";
      return;
    }
    for (size_t i=0;i<a.dimension(0);i++)
      ostr << "slice " << i << '\n' << a(i);
  }
};

template<class T> std::ostream& print(const T& a, std::ostream& ostr =std::cout) { 
  Print_<T,LCM_DIM(T)>::print(a,ostr);
  return ostr; }

//Here we insist output is already of correct size
template<size_t dim1, size_t dim2, size_t dim3> struct Apply_sd_ {};

template<size_t N> struct Apply_sd_<N,N,N> {
  template<class T1,class T2,class T3,class Function> LCM_INLINE static void func(Function mop, T1& d, const T2& a, const T3& b) {
    if (!arematching(a,b)) 
      throw Mismatch("apply_sd<N,N,N>");
    duplicate_structure(d,a);
    iters3_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator),LCM_ITER_CAT(typename T3::const_iterator)>::apply_ip2(mop,d.begin(),d.end(),a.begin(),b.begin());
  }
};
template<size_t N> struct Apply_sd_<N,N,0> {
  template<class T1,class T2,class T3,class Function> LCM_INLINE static void func(Function mop,T1& d, const T2& a, const T3& b) {
    duplicate_structure(d,a);
    iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T2::const_iterator)>::
      apply_ip1st(mop,d.begin(),d.end(),a.begin(),b);
  }
};
template<size_t N> struct Apply_sd_<N,0,N> {
  template<class T1,class T2,class T3,class Function> LCM_INLINE static void func(Function mop, T1& d, const T2& a, const T3& b) {
    duplicate_structure(d,b);
    iters2_<LCM_ITER_CAT(typename T1::iterator),LCM_ITER_CAT(typename T3::const_iterator)>::
      apply_ip2nd(mop,d.begin(),d.end(),a,b.begin());
  }
};
template<> struct Apply_sd_<0,0,0> {
  template<class T1,class T2,class T3,class Function> inline static void func(Function mop,T1& d, const T2& a, const T3& b) { 
    mop(d,a,b);
  }
};

template<class Function,class T1,class T2,class T3> void apply_sd(T1& d,Function mop,const T2& a, const T3& b) {
  Apply_sd_<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(mop,d,a,b);
}

template<class T1,class T2,class T3> inline void add(T1& d,const T2& a,const T3& b) {
  Apply_sd_<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(doesadd_sd<LCM_VAL(T1),LCM_VAL(T2),LCM_VAL(T3)>(), d, a, b);
}

template<class T1,class T2,class T3> inline void subtract(T1& d,const T2& a,const T3& b) {
  Apply_sd_<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(doessubtract_sd<LCM_VAL(T1),LCM_VAL(T2),LCM_VAL(T3)>(), d, a, b);
}

template<typename T1,typename T2,typename T3> void emultiply(T1& d,const T2& a, const T3& b) {
  Apply_sd_<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(doesmultiply_sd<LCM_VAL(T1),LCM_VAL(T2),LCM_VAL(T3)>(), d, a, b); 
}

 template<typename T1,typename T2,typename T3> void divide(T1& d,const T2& a, const T3& b) {
   Apply_sd_<LCM_DIM(T1),LCM_DIM(T2),LCM_DIM(T3)>::func(doesdivide_sd<LCM_VAL(T1),LCM_VAL(T2),LCM_VAL(T3)>(), d, a, b); 
 }

template<class T1,class T2> bool areequal(const T1& a,const T2& b)
{
  if (!arematching(a,b))
    return false;
  typename T1::const_iterator astart=a.begin();
  typename T2::const_iterator bstart=b.begin();
  size_t n=a.size();
  for (;n--;) {
    if (*astart!=*bstart) 
      return false;
    ++astart; 
    ++bstart;
  }
  return true;
}

bool areequal(double,double,double);

template<class T1,class T2> bool areequal(const T1& a,const T2& b, double tol)
{
  if ((tol>=1.0) || (tol<=0.0))
    throw InvalidParameter("areequal");
  if (!arematching(a,b))
    return false;
  typename T1::const_iterator astart=a.begin();
  typename T2::const_iterator bstart=b.begin();
  size_t n=a.size();
  for (;n--;) {
    if (!areequal(*astart,*bstart,tol))
      return false;
    ++astart; 
    ++bstart;
  }
  return true;
}

template<class T1,class T2> bool arenotequal(const T1& a,const T2& b)
{
  if (!arematching(a,b))
    return true;
  typename T1::const_iterator astart=a.begin();
  typename T2::const_iterator bstart=b.begin();
  size_t n=a.size();
  for (;n--;) {
    if (*astart!=*bstart) 
      return true;
    ++astart;
    ++bstart;
  }
  return false;
}

template<class T> LCM_VAL(T) sum(const T& a)
{
  typename T::const_iterator start=a.begin();
  const typename T::const_iterator end=a.end();
  if (start==end)
    throw Failed("sum: container is empty!");
  LCM_VAL(T) s=*start;
  while (++start!=end)
    s+=*start;
  return s;
}

template<class T,class F> typename F::result_type sum(const T& a,F func)
{
  typename T::const_iterator start=a.begin();
  const typename T::const_iterator end=a.end();
  if (start==end)
    throw Failed("sum: container is empty!");

  typename F::result_type s=func(*start);
  while (++start!=end)
    s+=func(*start);
  return s;
}

template<typename T> inline bool ismember(const T& a, const typename T::value_type& v)
{
  const typename T::const_iterator end(a.end());
  return (std::find(a.begin(),end,v)!=end);
}

template<class T,class F> LCM_VAL(T) max(const T& a,F func)
{
  typename T::const_iterator start=a.begin();
  const typename T::const_iterator end=a.end();
  if (start==end)
    throw Failed("max: container is empty!");

  LCM_VAL(T) maxv=*start;
  while (++start!=end) {
    if (func(maxv,*start))
      maxv=*start;
  }
  return maxv;
}

template<class T,class F> LCM_VAL(T) min(const T& a,F func)
{
  typename T::const_iterator start=a.begin();
  const typename T::const_iterator end=a.end();
  if (start==end)
    throw Failed("min: container is empty!");

  LCM_VAL(T) minv=*start;
  while (++start!=end) {
    if (func(*start,minv))
      minv=*start;
  }
  return minv;
}

template<class T> LCM_VAL(T) max(const T& a) {
  return max(a,std::less<LCM_VAL(T)>()); }

template<class T> LCM_VAL(T) min(const T& a) {
  return min(a,std::less<LCM_VAL(T)>()); }

template<typename T> struct doessum : public ::std::unary_function<T,LCM_VAL(T)> {
  inline LCM_VAL(T) operator()(const T& a) const { 
    return sum(a); }
};
template<typename T,typename F =std::less<LCM_VAL(T)> > struct doesmin : public ::std::unary_function<T,LCM_VAL(T)> {
  inline LCM_VAL(T) operator()(const T& a) const {
    return min(a,F()); }
};

 template<size_t N> struct Sum_ {
   template<class Func,class T> static inline typename Func::result_type func(T& a,Func fobj) { 
     return sum(a,fobj); }
 };
 template<> struct Sum_<0> {
   template<class Func,class T> static inline typename Func::result_type func(T& a,Func fobj) { 
     return fobj(a); }
 };

//interface to std::count
template<typename T1,typename T2> inline size_t count(const T1& a,const T2& v)
{ return ::std::count(a.begin(),a.end(),v); }

template<typename T> bool isconstant(const T& a)
{
  typename T::const_iterator start=a.begin();
  const typename T::const_iterator end=a.end();
  
  if (start==end)
    throw Undefined("isconstant");
  typename T::const_reference a0=*start;

  while (++start!=end) {
    if (*start!=a0)
      return false;
  }
  return true;
}

} //namespace libcmatrix

#endif
