#ifndef _MatrixStack_h_
#define _MatrixStack_h_

#include <cmath>
#include "List.h"
#ifdef _REENTRANT
#include "cmatrix_threads.h"
extern "C" {
  typedef void (*P_DELFUNC)(void*);
}
#endif

namespace libcmatrix {

template<typename T> class MatrixStack {
    List< Stack<T*> > stacks;
    const bool tidy_on_exit;
    size_t max_buf;
    
    void clear();
    void docreate(size_t);
  public:
    MatrixStack(size_t size, bool tidy_ =true) : tidy_on_exit(tidy_) { docreate(size); }
    ~MatrixStack() { if (tidy_on_exit) clear(); }
    
    void resize(size_t newsize) { clear(); docreate(newsize); }
    
    Stack<T*>* operator()(size_t i) { return (i<max_buf) ? &(stacks(i)) : NULL; }
  };

template<typename T> void MatrixStack<T>::docreate(size_t max_items)
{
  max_buf=size_t(std::sqrt(double(max_items)));
  
  stacks.create(max_buf,Stack<T*>(0));
  for (size_t i=1;i<max_buf;i++) {
    const size_t number=max_items/(i*i);
    stacks(i).reserve(number);
  }
}

template<typename T> void MatrixStack<T>::clear()
{
  for (size_t i=1;i<max_buf;i++) {
    Stack<T*>& stack=stacks(i);
    const size_t length=i*i;
    while (!stack.empty()) {
      T* data=stack.back();
      rawalloc::destroy(data,data+length);
      rawalloc::release_memory(data);
      stack.pop_back();
    }
  }
}

#ifdef _REENTRANT

  template<typename T> class ThreadedMatrixStack {
  public:
    ThreadedMatrixStack(int stack_size_, P_DELFUNC pfunc) : stack_size(stack_size_), stack_holder(pfunc,new MatrixStack<T>(stack_size_)) {}

    Stack<T*>* operator()(size_t r) {
      MatrixStack<T>* stackp=stack_holder.get();
      if (!stackp) { //stack doesn't exist?
	stackp=new MatrixStack<T>(stack_size);
	stack_holder.set(stackp);
      }
      return (*stackp)(r);
    }
    void resize(size_t newsize) {
      MatrixStack<T>* stackp=stack_holder.get();
      if (stackp) stackp->resize(newsize);
      stack_size=newsize;
    }

  private:
    size_t stack_size;
    key_holder< MatrixStack<T> > stack_holder;
  };

#endif

}//namespace libcmatrix

#endif // _MatrixStack_h_
