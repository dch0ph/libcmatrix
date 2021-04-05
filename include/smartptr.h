#ifndef smartptr_h_
#define smartptr_h_

/* Simple smart pointer for polymorphic types (must have clone() function)
   If flag is mxflag::normal (default), smartptr will take over responsibility for memory allocation NB object must have been created via new!
   If nondynamic, pointer is freely copied i.e. object is assumed permanent and shareable */

#include "basedefs.h"
 
namespace libcmatrix {
      
  template<class T, bool isPolyMorphic =true> class smartptr {
  public:
    typedef T value;

    smartptr() : ptr_(NULL), nondynamic_(true) {}
    
    explicit smartptr(T* ptrv, mxflag::tempflag tflag =mxflag::normal) 
      : ptr_(ptrv), nondynamic_(tflag==mxflag::nondynamic) {}
    
    smartptr(const smartptr& a)
  : ptr_(create(a)),
      nondynamic_(a.nondynamic_) {}
    
    ~smartptr() 
      { if (!nondynamic_)
	  delete ptr_;
      }
    
    smartptr& operator= (const smartptr& a)
    {
      if (this!=&a) {
	T* newptr=create(a);
	if (!nondynamic_)
	  delete ptr_;
	ptr_=newptr;
	nondynamic_=a.nondynamic_;
      }
      return *this;
    }

    T* get() const { 
      return ptr_;
    }
    
    //NB no 'release' as there is no means to distinguish dynamic / nondynamic

  void reset(T* ptrv, mxflag::tempflag tflag =mxflag::normal) 
  {
    if (!nondynamic_)
      delete ptr_;
    ptr_=ptrv;
    nondynamic_=(tflag==mxflag::nondynamic);
  }

    bool operator!() const { return !ptr_; }
    T& operator* () const {
      validate();
      return *ptr_;
    }
    T* operator-> () const {
      validate();
      return ptr_;
    }

    void clear() {
      if (!nondynamic_)
	delete ptr_;
      ptr_=NULL;
    }

  bool operator== (const smartptr<T,isPolyMorphic>& a) const
  { return (ptr_==a.ptr_); }

  bool operator!= (const smartptr<T,isPolyMorphic>& a) const
  { return (ptr_!=a.ptr_); }

  private:
    T* ptr_;
    bool nondynamic_;

  void validate() const { 
    if (!ptr_)
      throw Failed("smartptr: attempt to dereference null pointer");
  }

    template<class Tin> static T* cast_(Tin* a) { return dynamic_cast<T*>(a); }
    static T* cast_(T* a) { return a; }

  static T* create(const smartptr<T,false>& a) {
    return ((a.nondynamic_ || (a.ptr_==NULL)) ? a.ptr_ : new T(*(a.ptr_)));
  }
  static T* create(const smartptr<T,true>& a) {
    return ((a.nondynamic_ || (a.ptr_==NULL)) ? a.ptr_ : cast_((a.ptr_)->clone()));
  }

  };

}
#endif
