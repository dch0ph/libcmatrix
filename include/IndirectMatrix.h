#ifndef _IndirectMatrix_h_
#define _IndirectMatrix_h_

namespace libcmatrix {

 template<typename T, typename Map1 =BaseList<size_t>, typename Map2 =BaseList<size_t> > class IndirectMatrix {
 public:
   IndirectMatrix(T*, int, const Map1& map1,const Map2& map2);

   size_t rows() const { return rlength_; }
   size_t cols() const { return clength_; }
   size_t size() const { return rlength_*clength_; }
   static bool empty() { return false; }

   void create(size_t r,size_t s) { 
     if (r!=rlength_ || s!=clength_)
       throw Mismatch("IndirectMatrix::create");
   }

   size_t dimension(size_t n) const {
     switch (n) {
     case 0: return rows();
     case 1: return cols();
     default:
       throw BadIndex("IndirectMatrix::dimension");
     }
   }
   
   template<typename T2> IndirectMatrix& operator= (const T2& a) {
     assign(*this,a);
     return *this;
   }
//    IndirectMatrix& operator= (const IndirectMatrix<T,Map1,Map2>& a) {
//      assign(*this,a);
//      return *this;
//    }
   template<class T2> IndirectMatrix& operator+= (const T2& a) { return add_ip(*this,a); }
   template<class T2> IndirectMatrix& operator-= (const T2& a) { return subtract_ip(*this,a); }
   template<class T2> IndirectMatrix& operator*= (const T2& a) { return multiply_ip(*this,a); }
   template<class T2> IndirectMatrix& operator/= (const T2& a) { 
     apply_ip(doesdivide_ip<T,LCM_VAL(T2)>(),*this,a);
     return *this;
   }

//    T* vector() { return data_; }
    const T* vector() const { return data_; }
//    T* vector(size_t r) { return data_+c_*rindex_(r); }
//    const T* vector(size_t r) const { return data_+c_*rindex_(r); }

   static void clear() { throw Failed("IndirectMatrix can't be explicitly cleared"); }

   class base_iterator_ {
   protected:
     T* data_;
     T* crow;
     size_t r_;
     int c_;
   
     const size_t cols_;
     const Map1& rindex_;
     const size_t rlength_;
     const Map2& cindex_;
     const size_t clength_;

     void up() {
       if (++c_==clength_) 
	 reset(r_+1,0);
     }
     void down() {
       if (--c_<0)
	 reset(r_-1,clength_-1);
     }
     void reset(size_t r,int c) { r_=r; c_=c;
				  if (r<rlength_)
				    crow=data_+cols_*rindex_(r_);
     } // Need to catch iterator falling out of range

   public:
     base_iterator_(const IndirectMatrix<T,Map1,Map2>& a,size_t r,size_t c)
       : data_(a.data_), cols_(a.c_), rindex_(a.rindex_), rlength_(a.rlength_), cindex_(a.cindex_), clength_(a.clength_) { reset(r,c); }
     bool operator== (const base_iterator_& x) const { return (x.c_==c_) && (x.r_==r_); }
     bool operator!= (const base_iterator_& x) const { return (x.c_!=c_) || (x.r_!=r_); }
   };
   friend class base_iterator_;

   struct iterator : public base_iterator_, ::std::iterator< ::std::bidirectional_iterator_tag,T> {
     using base_iterator_::cindex_;
     using base_iterator_::crow;
     using base_iterator_::c_;
     iterator(IndirectMatrix<T,Map1,Map2>& a,size_t r,size_t c) : base_iterator_(a,r,c) {}
     T& operator*() const { return crow[cindex_(c_)]; }
     iterator& operator++() { this->up(); return *this; }
     iterator operator++(int) { iterator tmp(*this); this->up(); return tmp; }
     iterator& operator--() { this->down(); return *this; }
     iterator operator--(int) { iterator tmp(*this); this->down(); return tmp; }
   };

   struct const_iterator : public base_iterator_, ::std::iterator< ::std::bidirectional_iterator_tag,T,ptrdiff_t,const T*,const T&> {
     using base_iterator_::cindex_;
     using base_iterator_::crow;
     using base_iterator_::c_;
     const_iterator(const IndirectMatrix<T,Map1,Map2>& a,size_t r,size_t c) : base_iterator_(a,r,c) {}
     const T& operator*() const { return crow[cindex_(c_)]; }
     const_iterator& operator++() { this->up(); return *this; }
     const_iterator operator++(int) { const_iterator tmp(*this); this->up(); return tmp; }
     const_iterator& operator--() { this->down(); return *this; }
     const_iterator operator--(int) { const_iterator tmp(*this); this->down(); return tmp; }
   };

   iterator begin() { return iterator(*this,0,0); }
   iterator end() { return iterator(*this,rlength_,0); }
   const_iterator begin() const { return const_iterator(*this,0,0); }
 const_iterator end() const { return const_iterator(*this,rlength_,0); }
 
 typedef T value_type;
 typedef T& reference;
 typedef const T& const_reference;  
 
   T& operator() (size_t r,size_t s) {
#if (BOUNDS_CHECK)
     if (r>=rlength_ || s>=clength_)
       throw BadIndex("IndirectMatrix()",r,rlength_,s,clength_);
 #endif
     return data_[c_*rindex_(r)+cindex_(s)];
   }
   const T& operator() (size_t r,size_t s) const {
#if (BOUNDS_CHECK)
     if (r>=rlength_ || s>=clength_)
       throw BadIndex("IndirectMatrix()",r,rlength_,s,clength_);
 #endif
     return data_[c_*rindex_(r)+cindex_(s)];
   }
 
 private:
   IndirectMatrix();

   T* data_;
   const size_t c_;
   const Map1 rindex_;
   const size_t rlength_;
   const Map2 cindex_;
   const size_t clength_;
 };

 template<typename T,typename Map1,typename Map2> struct type_traits< IndirectMatrix<T,Map1,Map2> > {
   static const int rank=0;
   static const bool trivialconstructor=false;
   static const size_t dimensionality=2;
   typedef T value_type;
 };

 template<typename T,typename Map1,typename Map2> std::ostream& operator<< (std::ostream& ostr, const IndirectMatrix<T,Map1,Map2>& a)
   {
     print(a,ostr);
     return ostr;
   }

 template<typename T,typename Map1,typename Map2> IndirectMatrix<T,Map1,Map2>::IndirectMatrix(T* data, int cv,const Map1& rindex, const Map2& cindex)
     : data_(data), c_(cv), rindex_(rindex), rlength_(rindex.size()), cindex_(cindex), clength_(cindex.size()) 
     {
       if (size()==0)
	 throw InvalidParameter("IndirectMatrix cannot be null");
     }

}// namespace libcmatrix

#endif
