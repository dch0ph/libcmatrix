#ifndef _PermutedMultiMatrix_h_
#define _PermutedMultiMatrix_h_

#include "MultiMatrix.h"

namespace libcmatrix {

  //can't be created for N<2
 template<typename T,size_t N> class PermutedMultiMatrix : private BaseList<T> {
 public:

   PermutedMultiMatrix(const MultiMatrix<T,N>&, const BaseList<size_t>& =BaseList<size_t>());
   template<typename T2> PermutedMultiMatrix& operator= (const T2& a) { assign(*this,a); return *this; }

   typedef T value_type;
   typedef T& reference;
   typedef const T& const_reference;  
   typedef size_t size_type;

   size_t size() const { return BaseList<T>::length(); }
   bool operator!() const { return false; }
   static bool empty() { return false; }
   size_t rows() const { LCM_STATIC_CHECK(N==2, Non2D_object_passed_to_rows); return dim[1]; }
   size_t cols() const { LCM_STATIC_CHECK(N==2, Non2D_object_passed_to_cols); return dim[0]; }
   
  size_t dimension(size_t n) const {
#ifndef NDEBUG
    if (n>=N) throw BadIndex("MultiMatrix::dimension");
#endif
    return dim[N-n-1]; }

  size_t index(size_t r,size_t s) const { 
    LCM_STATIC_CHECK(N==2,Not2D);
#ifndef NDEBUG
    if (r>=dim[1] || s>=dim[0]) throw BadIndex("MultiMatrix::check_bounds");
#endif
    return r*mults[1]+s*mults[0]; }

  size_t index(size_t r,size_t s,size_t t) const { 
    LCM_STATIC_CHECK(N==3,Not3D);
#ifndef NDEBUG
    if (r>=dim[2] || s>=dim[1] || t>=dim[0]) throw BadIndex("MultiMatrix::check_bounds");
#endif
    return r*mults[2]+s*mults[1]+t*mults[0]; }

  size_t index(size_t r,size_t s,size_t t,size_t u) const {
    LCM_STATIC_CHECK(N==4,Not4D);
#ifndef NDEBUG
    if (r>=dim[3] || s>=dim[2] || t>=dim[1] || u>=dim[0]) throw BadIndex("MultiMatrix::check_bounds");
#endif
    return r*mults[3]+s*mults[2]+t*mults[1]+u*mults[0]; }

  T& operator()(size_t r,size_t s) { return BaseList<T>::operator()(index(r,s)); }
  T& operator()(size_t r,size_t s,size_t t) { return BaseList<T>::operator()(index(r,s,t)); }
  T& operator()(size_t r,size_t s,size_t t,size_t u) { return BaseList<T>::operator()(index(r,s,t,u)); }
  const T& operator()(size_t r,size_t s) const { return BaseList<T>::operator()(index(r,s)); }
  const T& operator()(size_t r,size_t s,size_t t) const { return BaseList<T>::operator()(index(r,s,t)); }
  const T& operator()(size_t r,size_t s,size_t t,size_t u) const { return BaseList<T>::operator()(index(r,s,t,u)); }

  friend class _base_iterator;
  class _base_iterator {
  protected:
    T* const data_;
    const size_t* const dim;
    const size_t* const mults;
    T* cpos;
    int cindex[N];

    void up() { 
      if (++(cindex[0])!=dim[0]) return;      
      int i=0;
      while (i<N-1) {
	cindex[i++]=0;
	if (++(cindex[i])!=dim[i]) { reset(); return; }
      }
    }
    void down() {
      if (--(cindex[0])>=0) return;
      int i=0;
      while (i<N-1) {
	cindex[i]=dim[i]-1;
	if (--(cindex[++i])>=0) { reset(); return; }
      }
      //iterator now invalid
    }

    void reset() { 
      cpos=data_;
      for (int i=1;i<N;i++) cpos+=mults[i]*cindex[i];
    }

   public:
     _base_iterator(const PermutedMultiMatrix<T,N>& a,bool isatendv)
       : data_(const_cast<T*>(a.vector())), dim(a.dim), mults(a.mults) { 
       if (isatendv) {
	 for (int i=N;i--;) cindex[i]=dim[i]-1;
	 cindex[0]++;
	 cpos=data_+a.length();
       }
       else {
	 for (int i=N;i--;) cindex[i]=0;
	 cpos=data_;
       }
     }
    
     bool operator!= (const _base_iterator& x) const { return (cindex[0]!=x.cindex[0]) || (cpos!=x.cpos); }
     bool operator== (const _base_iterator& x) const { return (cindex[0]==x.cindex[0]) && (cpos==x.cpos); }
   };

   struct iterator : public _base_iterator, public ::std::iterator< ::std::bidirectional_iterator_tag,T> {
     iterator(PermutedMultiMatrix<T,N>& a,bool isatendv) : _base_iterator(a,isatendv) {}
    
     T& operator*() const { return cpos[cindex[0]*mults[0]]; }
     //Don't use postfix forms if possible!
     iterator& operator++() { up(); return *this; }
     iterator operator++(int) { iterator tmp(*this); up(); return tmp; }
     iterator& operator--() { down(); return *this; }
     iterator operator--(int) { iterator tmp(*this); down(); return tmp; }
  };

   struct const_iterator : public _base_iterator, ::std::iterator< ::std::bidirectional_iterator_tag,T,ptrdiff_t,const T*,const T&> {
     const_iterator(const PermutedMultiMatrix<T,N>& a,bool isatendv) : _base_iterator(a,isatendv) {}
    
     const T& operator*() const {
       //       std::cout << "Base: " << (cpos-data_) << "  Index[0]: " << cindex[0] << "  offset: " << cindex[0]*mults[0]; 
       return cpos[cindex[0]*mults[0]]; }
     //Don't use postfix forms if possible!
     const_iterator& operator++() { up(); return *this; }
     const_iterator operator++(int) { const_iterator tmp(*this); up(); return tmp; }
     const_iterator& operator--() { down(); return *this; }
     const_iterator operator--(int) { const_iterator tmp(*this); down(); return tmp; }
  };

   iterator begin() { return iterator(*this,false); }
   const_iterator begin() const { return const_iterator(*this,false); }
   iterator end() { return iterator(*this,true); }
   const_iterator end() const { return const_iterator(*this,true); }

 private:
   PermutedMultiMatrix();

   //N.B. These are stored in reverse order to MultiMatrix!
   size_t dim[N];
   size_t mults[N];
 };

 template<typename T,size_t N> struct type_traits< PermutedMultiMatrix<T,N> > {
  static const bool known=true;
  static const bool trivialconstructor=false; 
  static const size_t dimensionality=N;
  static const size_t rank=0;
  typedef T value_type;
 };
 
 template<typename T,size_t N> std::ostream& operator<< (std::ostream& ostr, const PermutedMultiMatrix<T,N>& a)
   {
     Print_< PermutedMultiMatrix<T,N> ,N>::print(a,ostr);
     return ostr;
   }

  template<typename T,size_t N> PermutedMultiMatrix<T,N>::PermutedMultiMatrix(const MultiMatrix<T,N>& a,const BaseList<size_t>& order) : BaseList<T>(a.row()) {
    LCM_STATIC_CHECK(N>1, PermutedMultiMatrix_applied_to_lessthan_2D);
     int i;
     if (order.length()) {
       if (order.length()!=N) throw Mismatch("MultiMatrix::iterator");
       for (i=N;i--;) dim[0]=0;
       for (i=N;i--;) { 
	 const int which=order(i);
	 const int nwhich=N-which-1;
	 if (which>=N || dim[nwhich]) throw InvalidParameter("MultiMatrix::iterator: bad permutation vector");
	 dim[nwhich]=a.dim[which];
	 mults[nwhich]=a.mults[which];
       }
     }
     else {
       for (i=N;i--;) {
	 dim[i]=a.dim[i];
	 mults[i]=a.mults[i];
       }
     }
   }

}// namespace libcmatrix

#endif
