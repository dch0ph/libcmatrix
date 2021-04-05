#ifndef _MultiHolder_h_
#define _MultiHolder_h_

/* Simple type holding data of different dimensionality
   in single object */

#include "List.h"
#include "Matrix.h"
#include "MultiMatrix.h"
#include "UnionHolder.h"
#include <iostream>

namespace libcmatrix {

#define LCM_MULTIHOLDER_DIMS 4

#define LCM_HOLDER(T) UnionHolder<5,T,List<T>,Matrix<T>,MultiMatrix<T,3>,MultiMatrix<T,4> >

 template<class T> class MultiHolder : public LCM_HOLDER(T) {
  public:
   MultiHolder() : LCM_HOLDER(T)() {}

     using LCM_HOLDER(T)::get;
     using LCM_HOLDER(T)::set;
     using LCM_HOLDER(T)::type;
     using LCM_HOLDER(T)::operator();

     template<class T2> explicit MultiHolder(const T2& a) : LCM_HOLDER(T)(a) {}

       explicit MultiHolder(const BaseList< Matrix<T> >& a)
	 : LCM_HOLDER(T)(Int2Type<4>()) { get(Int2Type<4>())=a; }

       explicit MultiHolder(const BaseList<T>& a)
	 : LCM_HOLDER(T)(Int2Type<2>()) { get(Int2Type<2>())=a; }
	 
       template<class T2> MultiHolder& operator= (const T2& a)
	 { LCM_HOLDER(T)::operator=(a); return *this; }

       MultiHolder& operator= (const MultiHolder<T>& a)
	 { LCM_HOLDER(T)::operator=( static_cast< LCM_HOLDER(T) >(a)); return *this; }

       MultiHolder& operator= (const BaseList<T>& a)
	 { LCM_HOLDER(T)::set(Int2Type<2>())=a; return *this; }

       MultiHolder& operator= (const BaseList< Matrix<T> >& a) { 
	 LCM_HOLDER(T)::set(Int2Type<4>())=a;
	 return *this;
       }

       int dimensions() const { return int(type())-1; }
       void dimensions(int ndim) { type(ndim+1); }
       //  size_t dimension(size_t) const;
       //N.B. MultiHolder object may be initialised but contain uninitialised data!
       
       T& create() { return set(Int2Type<1>()); }
       List<T>& create(int r)
	 { List<T>& tmp(set(Int2Type<2>())); tmp.create(r); return tmp; }
       Matrix<T>& create(int r,int s)
	 { Matrix<T>& tmp(set(Int2Type<3>())); tmp.create(r,s); return tmp; }
       MultiMatrix<T,3>& create(int r,int s,int t)
	 { MultiMatrix<T,3>& tmp(set(Int2Type<4>())); tmp.create(r,s,t); return tmp; }
       MultiMatrix<T,4>& create(int r,int s,int t,int u)
	 { MultiMatrix<T,4>& tmp(set(Int2Type<5>())); tmp.create(r,s,t,u); return tmp; }
       
       typedef T* iterator;
       typedef const T* const_iterator;
       iterator begin();
       const_iterator begin() const;
       iterator end();
       const_iterator end() const;
       
       typedef T value_type;
       typedef T& reference;
       typedef const T& const_reference;  
       typedef int difference_type;
       typedef size_t size_type;
       
       size_t size() const;
       bool empty() const { return (size()==0); }
              
       T& scalar() { return operator()(Int2Type<1>()); }
       const T& scalar() const { return operator()(Int2Type<1>()); }
       
       List<T>& list() { return operator()(Int2Type<2>()); }
       const List<T>& list() const { return operator()(Int2Type<2>()); }
       
       Matrix<T>& matrix() { return operator()(Int2Type<3>()); }
       const Matrix<T>& matrix() const { return operator()(Int2Type<3>()); }
       
       MultiMatrix<T,3>& multimatrix3() { return operator()(Int2Type<4>()); }
       const MultiMatrix<T,3>& multimatrix3() const { return operator()(Int2Type<4>()); }
       
       MultiMatrix<T,4>& multimatrix4() { return operator()(Int2Type<5>()); }
       const MultiMatrix<T,4>& multimatrix4() const { return operator()(Int2Type<5>()); }
  };

// Implementation details below
 
template<class T> size_t MultiHolder<T>::size() const
{
  switch (type()) {
  case 1: return 1;
  case 2: return get(Int2Type<2>()).size();
  case 3: return get(Int2Type<3>()).size();
  case 4: return get(Int2Type<4>()).size();
  case 5: return get(Int2Type<5>()).size();
  default: throw Undefined("MultiHolder<T>::size()");
  }
}
 
 template<class T> T* MultiHolder<T>::begin()
   {
     switch (type()) {
     case 1: return &(get(Int2Type<1>()));
     case 2: return get(Int2Type<2>()).begin();
     case 3: return get(Int2Type<3>()).begin();
     case 4: return get(Int2Type<4>()).begin();
     case 5: return get(Int2Type<5>()).begin();
     default: throw Undefined("MultiHolder<T>::begin()");
     }
   }

 template<class T> const T* MultiHolder<T>::begin() const
   {
     switch (type()) {
     case 1: return &(get(Int2Type<1>()));
     case 2: return get(Int2Type<2>()).begin();
     case 3: return get(Int2Type<3>()).begin();
     case 4: return get(Int2Type<4>()).begin();
     case 5: return get(Int2Type<5>()).begin();
     default: throw Undefined("MultiHolder<T>::begin()");
     }
   }

 template<class T> T* MultiHolder<T>::end()
   {
     switch (type()) {
     case 1: return 1+&(get(Int2Type<1>()));
     case 2: return get(Int2Type<2>()).end();
     case 3: return get(Int2Type<3>()).end();
     case 4: return get(Int2Type<4>()).end();
     case 5: return get(Int2Type<5>()).end();
     default: throw Undefined("MultiHolder<T>::end()");
     }
   }
 
 template<class T> const T* MultiHolder<T>::end() const
   {
     switch (type()) {
     case 1: return 1+&(get(Int2Type<1>()));
     case 2: return get(Int2Type<2>()).end();
     case 3: return get(Int2Type<3>()).end();
     case 4: return get(Int2Type<4>()).end();
     case 5: return get(Int2Type<5>()).end();
     default: throw Undefined("MultiHolder<T>::end()");
     }
   }
      
template<class T1,class T2> bool arematching(const MultiHolder<T1>& a,const MultiHolder<T2>& b)
{
  const size_t n=a.dimensions();
  if (n!=b.dimensions()) return false;
  switch (n) {
  case 0: return true;
  case 1: return arematching(a.list(),b.list());
  case 2: return arematching(a.matrix(),b.matrix());
  case 3: return arematching(a.multimatrix3(),b.multimatrix3());
  case 4: return arematching(a.multimatrix4(),b.multimatrix4());
  }
  throw Undefined("arematching");
}

template<class T> void spy(std::ostream& ostr,const MultiHolder<T>& a,double tol =1e-10)
{
  switch (a.dimensions()) {
  case 0:
    spy(ostr,a.scalar(),tol);
    break;
  case 1:
    spy(ostr,a.list(),tol);
    break;
  case 2:
    spy(ostr,a.matrix(),tol);
    break;
  case 3:
    spy(ostr,a.listmatrix(),tol);
    break;
  case 4:
    spy(ostr,a.matrixmatrix(),tol);
    break;
  }
  ostr << "<undefined>" << std::endl;
}

template<typename T> struct type_traits<MultiHolder<T> > {
  static const bool known=true;
  static const bool trivialconstructor=false; 
  //attempt to check dimensionality or rank creates compile time error
  typedef T value_type;
};

} //namespace libcmatrix

#endif
