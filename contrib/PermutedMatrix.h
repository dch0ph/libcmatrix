#ifndef _PermutedMatrix_h_
#define _PermutedMatrix_h_

#include "Matrix.h"

namespace libcmatrix {

 template<typename T> class PermutedMatrix {
 public:
   PermutedMatrix(const Matrix<T>& a) : data_(const_cast<T*>(a.vector())), r_(a.rows()), c_(a.cols()) {}
   PermutedMatrix(T* datav,int rv,int cv) : data_(datav), r_(rv), c_(cv) { if ( (rv<1) || (cv<1)) throw InvalidParameter("PermutedMatrix: can't have zero dimension"); }

   typedef T value_type;
   typedef T& reference;
   typedef const T& const_reference;  
   typedef size_t size_type;

   static bool empty() { return false; }
   size_t rows() const { return c_; }
   size_t cols() const { return r_; }
   size_t size() const { return r_*c_; }
   bool operator!() const { return false; }

   void create(size_t r,size_t s) { if (r!=c_ || s!=r_) throw Mismatch("PermutedMatrix::create"); }

   size_t dimension(size_t n) const {
     switch (n) {
     case 0: return c_;
     case 1: return r_;
     default:
       throw BadIndex("PermutedMatrix::dimension");
     }
   }
   
   template<typename T2> PermutedMatrix& operator= (const T2& a) { assign(*this,a); return *this; }

   class _base_iterator {
   protected:
     T* const data_;
     const size_t rows_,cols_;

     T* ccol;
     int r_;
     size_t c_;
    
     void up() { if (++r_==rows_) reset(0,c_+1); }
     void down() { if (--r_<0) reset(rows_-1,c_-1); }
     void reset(int r,int c) { r_=r; c_=c; ccol=data_+c_; }

   public:
     _base_iterator(const PermutedMatrix<T>& a,int r,int c) : data_(a.data_), rows_(a.r_), cols_(a.c_) { reset(r,c); }
     bool operator== (const _base_iterator& x) const { return (x.c_==c_) && (x.r_==r_); }
     bool operator!= (const _base_iterator& x) const { return (x.c_!=c_) || (x.r_!=r_); }
   };
   friend class _base_iterator;

   struct iterator : public _base_iterator, public ::std::iterator< ::std::bidirectional_iterator_tag,T> {
     iterator(const PermutedMatrix<T>& a,size_t r,size_t c) : _base_iterator(a,r,c) {}
     T& operator*() const { return ccol[r_*cols_]; }
     iterator& operator++() { up(); return *this; }
     iterator operator++(int) { iterator tmp=*this; up(); return tmp; }
     iterator& operator--() { down(); return *this; }
     iterator operator--(int) { iterator tmp=*this; down(); return tmp; }
   };

   struct const_iterator : public _base_iterator, public ::std::iterator< ::std::bidirectional_iterator_tag,T,ptrdiff_t,const T*,const T&> {
     const_iterator(const PermutedMatrix<T>& a,size_t r,size_t c) : _base_iterator(a,r,c) {}
     const T& operator*() const { return ccol[r_*cols_]; }
     const_iterator& operator++() { up(); return *this; }
     const_iterator operator++(int) { const_iterator tmp=*this; up(); return tmp; }
     const_iterator& operator--() { down(); return *this; }
     const_iterator operator--(int) { const_iterator tmp=*this; down(); return tmp; }
   };

   iterator begin() { return iterator(*this,0,0); }
   iterator end() { return iterator(*this,r_,0); }
   const_iterator begin() const { return const_iterator(*this,0,0); }
   const_iterator end() const { return const_iterator(*this,r_,0); }

   typedef T value_type;
   typedef T& reference;
   typedef const T& const_reference;  

   T& operator() (size_t r,size_t s) {
#if (CHECK_BOUNDS)
     if (r>=c_ || s>=c_) throw BadIndex("PermutedMatrix()");
 #endif
     return data_[s*c_+r];
   }
   const T& operator() (size_t r,size_t s) const {
#if (CHECK_BOUNDS)
     if (r>=c_ || s>=r_) throw BadIndex("PermutedMatrix()");
 #endif
     return data_[s*c_+r];
   }
 
 private:
   PermutedMatrix();

   T* const data_;
   const size_t r_,c_;
 };

 template<typename T> struct type_traits< PermutedMatrix<T> > {
   static const bool known=true;
   static const int rank=0;
   static const bool isfundamental=false;
   static const size_t dimensionality=2;
   typedef T value_type;
 };

 template<typename T> std::ostream& operator<< (std::ostream& ostr, const PermutedMatrix<T>& a)
   {
     print(a,ostr);
     return ostr;
   }

}// namespace libcmatrix

#endif
