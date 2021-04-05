#ifndef Matrix_h_
#define Matrix_h_

#include "List.h"
#include "DynamicList.h"
#if (MAKE_CSTACK>0) || (MAKE_RSTACK>0)
#include "CachedAllocator.h"
#endif

namespace libcmatrix {

  template<class T> class Matrix;
  template<class T> class RawMatrix;
  template<class T> class BlockedMatrix;
template<class T,size_t N> class MultiMatrix;
template<typename T,typename Map1,typename Map2> class IndirectMatrix; 

 template<typename T> struct type_traits< Matrix<T> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=2;
   static const size_t rank=0;
   typedef T value_type;
 };

  template<typename T> struct matrix_traits {
    typedef DefaultAllocator< memory_traits<T>::alignment > allocator;
  };

//This is all we require to make cmatrix used a CachedAllocator - pretty funky...
//Type is only used to create separate arenas for complex and real
#if (MAKE_CSTACK>0)
template<> struct matrix_traits<COMPLEX_T> { 
  typedef CachedAllocator<ThreadingActive,COMPLEX_T> allocator; 
}; 
#endif
#if (MAKE_RSTACK>0)
template<> struct matrix_traits<float_t> { 
  typedef CachedAllocator<ThreadingActive,float_t> allocator; 
}; 
#endif

template<class T> bool isdefined(const Matrix<T>&);

 struct MatrixHelper_ {
   template<typename Map> static inline void check_bounds1(const Map& map,size_t maxv);
   static inline void check_bounds0(size_t i,size_t maxv);
   template<class Map> static inline const Map& expand(const Map& map, size_t maxv) {
     check_bounds1(map,maxv);
     return map;
   }
#ifndef LCM_SUPPRESS_VIEWS
   static inline void check_boundssl(const slice&,size_t);

   //create new range object
   static inline const range expand(const range& map, size_t maxv) {
     if (map.size()) {
       check_bounds1(map,maxv);
       return map;
     }
     else
       return range(0,maxv-1);
   }
#endif
   
 };

#ifndef LCM_SUPPRESS_VIEWS
  template<typename T,size_t Rank1,typename Map1,size_t Rank2,typename Map2> struct LCM_Matrix_Select {};
  template<typename T,typename Map1,typename Map2> struct LCM_Matrix_Select<T,0,Map1,1,Map2> {
    typedef IndirectList<T,Map2> return_type;
    static return_type func(Matrix<T>& a,size_t r,const Map2& map) { 
      MatrixHelper_::check_bounds0(r,a.rows());
      //      MatrixHelper_::check_bounds1(map,a.cols());
      const size_t offset=r*a.cols();
      return return_type(BaseList<T>(a.size()-offset,a.vector()+offset,map));
    }
  };
  template<typename T,typename Map1> struct LCM_Matrix_Select<T,0,Map1,1,range> {
    typedef BaseList<T> return_type;
    static return_type func(Matrix<T>& a,size_t r,const range& map) { 
      MatrixHelper_::check_bounds0(r,a.rows());
      if (map.size()) {
	MatrixHelper_::check_boundssl(map,a.cols());
	return BaseList<T>(map.size(),a.vector()+r*a.cols()+map.start());
      }
      else
	return BaseList<T>(a.cols(),a.vector()+r*a.cols());
    }
  };
  template<typename T,typename Map1,typename Map2> struct LCM_Matrix_Select<T,1,Map1,1,Map2> {
    typedef IndirectMatrix<T,Map1,Map2> return_type;
    static return_type func(Matrix<T>& a,const Map1& map1,const Map2& map2) { 
      return return_type(a.vector(),a.cols(),MatrixHelper_::expand(map1,a.rows()),MatrixHelper_::expand(map2,a.cols()));    }
  };
    
  template<typename T,size_t Rank1,typename Map1,size_t Rank2,typename Map2> struct LCM_Matrix_const_Select {};
  template<typename T,typename Map1,typename Map2> struct LCM_Matrix_const_Select<T,0,Map1,0,Map2> {
    typedef const T& return_type;
    static const T& func(const Matrix<T>& a,size_t r,size_t s) { return a.operator()(r,s); }
  };
  template<typename T,typename Map1,typename Map2> struct LCM_Matrix_const_Select<T,1,Map1,0,Map2> {
    typedef const IndirectList<T,Map1> return_type;
    static return_type func(const Matrix<T>& a,const Map1& map,size_t s) {
      //      MatrixHelper_::check_bounds1(map,a.rows());
      MatrixHelper_::check_bounds0(s,a.cols());
      return return_type(BaseList<T>(a.size()-s,a.vector()+s),map,a.cols()); 
    }
  };
  template<typename T,typename Map> struct LCM_Matrix_const_Select<T,1,range,0,Map> {
    typedef const IndirectList<T,slice> return_type;
    static return_type func(const Matrix<T>& a,const range& map,size_t s) {
      MatrixHelper_::check_bounds0(s,a.cols());
      if (map.size()) {
	//	MatrixHelper_::check_bounds1(map,a.rows());
	return return_type(BaseList<T>(a.size()-s,a.vector()+s),map,a.cols()); 
      }
      else
	return return_type(BaseList<T>(a.size()-s,a.vector()+s),slice(0,a.rows(),a.cols()));
    }
  };
  template<typename T,typename Map1,typename Map2> struct LCM_Matrix_const_Select<T,0,Map1,1,Map2> {
    typedef const IndirectList<T,Map2> return_type;
    static return_type func(const Matrix<T>& a,size_t r,const Map2& map) { 
      MatrixHelper_::check_bounds0(r,a.rows());
      //      MatrixHelper_::check_bounds1(map,a.cols());
      const size_t offset=r*a.cols();
      return return_type(BaseList<T>(a.size()-offset,a.vector()+offset),map);
    }
  };
  template<typename T,typename Map1> struct LCM_Matrix_const_Select<T,0,Map1,1,range> {
    typedef const BaseList<T> return_type;
    static return_type func(const Matrix<T>& a,size_t r,const range& map) { 
      MatrixHelper_::check_bounds0(r,a.rows());
      if (map.size()) {
	MatrixHelper_::check_boundssl(map,a.cols());
	return BaseList<T>(map.size(),const_cast<T*>(a.vector())+r*a.cols()+map.start());
      }
      else
	return BaseList<T>(a.cols(),const_cast<T*>(a.vector())+r*a.cols());
    }
  };
  template<typename T,typename Map1,typename Map2> struct LCM_Matrix_const_Select<T,1,Map1,1,Map2> {
    typedef const IndirectMatrix<T,Map1,Map2> return_type;
    static return_type func(const Matrix<T>& a,const Map1& map1,const Map2& map2) { 
      return return_type(const_cast<T*>(a.vector()),a.cols(),MatrixHelper_::expand(map1,a.rows()),MatrixHelper_::expand(map2,a.cols()));
    }
  };
    
  template<typename T,typename Map1,typename Map2> struct LCM_Matrix_Select<T,0,Map1,0,Map2> {
     typedef T& return_type;
     static T& func(Matrix<T>& a,size_t r,size_t s) { return a(r,s); }
   };
  template<typename T,typename Map1,typename Map2> struct LCM_Matrix_Select<T,1,Map1,0,Map2> {
     typedef IndirectList<T,Map1> return_type;
     static return_type func(Matrix<T>& a,const Map1& map,size_t s) {
       //       MatrixHelper_::check_bounds1(map,a.rows());
       MatrixHelper_::check_bounds0(s,a.cols());
       return return_type(BaseList<T>(a.size()-s,a.vector()+s),map,a.cols()); 
     }
   };
  template<typename T,typename Map> struct LCM_Matrix_Select<T,1,range,0,Map> {
     typedef IndirectList<T,slice> return_type;
     static return_type func(Matrix<T>& a,const range& map,size_t s) {
       MatrixHelper_::check_bounds0(s,a.cols());
       if (map.size()) {
	 MatrixHelper_::check_bounds1(map,a.rows());
	 return return_type(BaseList<T>(a.size()-s,a.vector()+s),map,a.cols()); 
       }
       else
	 return return_type(BaseList<T>(a.size()-s,a.vector()+s),slice(0,a.rows(),a.cols()));
     }
   }; 
  template<typename T> struct LCM_Matrix_Select<T,1,range,1,range> {
    typedef RawMatrix<T> return_type;
    static return_type func(Matrix<T>& a,const range& rmap, const range& cmap) {
      return RawMatrix<T>(a,rmap,cmap);
    }
  }; 
  template<typename T> struct LCM_Matrix_const_Select<T,1,range,1,range> {
    typedef const RawMatrix<T> return_type;
    static return_type func(const Matrix<T>& a,const range& rmap, const range& cmap) {
      return RawMatrix<T>(a,rmap,cmap);
    }
  }; 
#endif

 template<typename T> class Matrix {
public:
  explicit Matrix(mxflag::tempflag tflag =mxflag::normal)
    : store_(tflag), r_(0), c_(0) {}

  Matrix(int rows_,int cols_,mxflag::tempflag tflag =mxflag::normal) 
    : store_(rows_*cols_,tflag), r_(rows_), c_(cols_) {}

  Matrix(int rows_,int cols_,const T& def, mxflag::tempflag tflag =mxflag::normal)
    : store_(rows_*cols_,def,tflag), r_(rows_), c_(cols_) {}

  Matrix(int rows_,int cols_,const T* ptr, mxflag::tempflag tflag =mxflag::normal)
    : store_(rows_*cols_,ptr,tflag), r_(rows_), c_(cols_) {}

  Matrix(int rows_,int cols_, const BaseList<T>& a, mxflag::tempflag tflag =mxflag::normal)
    : store_(a.size(),a.vector(),tflag), r_(rows_), c_(cols_) {
    if (r_*c_!=store_.size())
      throw Mismatch("Matrix(): vector initialisation doesn't match matrix size", store_.size(),r_*c_);
  }

  inline Matrix(const Matrix&);
  inline Matrix(const Matrix&, mxflag::tempflag);

  template<class T2> explicit inline Matrix(const T2& a, mxflag::tempflag tflag =mxflag::normal) : store_(tflag)
    { LCM_STATIC_CHECK( LCM_DIM(T2)==2, Matrix_construct_from_non2D_object );
      store_.create_iter(a.rows()*a.cols(),a.begin());
      r_=a.rows();
      c_=a.cols();
    }

     inline bool istemporary() const
   { return store_.istemporary(); }
     
     inline bool isdynamic() const
   { return store_.isdynamic(); }
     
  inline Matrix<T>& operator=(const Matrix<T>&);
  inline Matrix<T>& operator=(const T& v);

  template<class T2> inline Matrix<T>& operator=(const T2& a) {
    assign(*this,a);
    return *this; }

  //void move(size_t r,size_t c,List<T>&); //"transfer" list into matrix

  void clear() {
    store_.clear();
    r_=c_=0;
  }
  //void ensureempty() { store_.ensureempty(); }

  inline void kill() { clear(); }

  // functions that "replace" an existing matrix
  inline void create(int rows,int cols);
  template<class T2> inline void create(int rows,int cols,T2*);
  inline void create(int rows,int cols,const T&);

  template<class T2> void set_dimensions(const T2& a) {
    //    typedef typename matching_N<2,LCM_DIM(T2)>::type check_type;
    LCM_STATIC_CHECK( LCM_DIM(T2)==2, set_dimensions_from_non2D_object );
    create(a.rows(),a.cols());
  }
      
  size_t rows() const { return r_; }
  size_t cols() const { return c_; }
   size_t size() const { return store_.size(); }

  bool operator!() const { return (store_.vector()==0); }

  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;  
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef int difference_type;
  typedef size_t size_type;

  T* vector() { 
    if (store_.vector())
      return store_.vector();
    throw Undefined("Matrix<T>::vector");
  }
  const T* vector() const {
    if (store_.vector()) 
      return store_.vector();
    throw Undefined("Matrix<T>::vector");
  }

  T* vector(size_t n) {
#if (BOUNDS_CHECK) 
    if (n>=r_)
      throw BadIndex("Matrix<T>::vector",n,r_);
#endif
    return store_.vector()+n*c_;
  }

  const T* vector(size_t n) const { 
#if (BOUNDS_CHECK)
    if (n>=r_)
      throw BadIndex("Matrix<T>::vector",n,r_);
#endif
    return const_cast<const T*>(vector()+n*c_);
  }

  T* begin() { return store_.vector(); }
  const T* begin() const { return store_.vector(); }
   T* end() { return store_.vector()+store_.size(); }
   const T* end() const { return store_.vector()+store_.size(); }

  bool empty() const { return (store_.vector()==0); }

   void resize(int,int);
   template<class T2> void push_back_row(const T2&);

  size_t dimension(size_t n) const {
    switch (n) {
    case 0: return r_;
    case 1: return c_;
    default: throw BadIndex("Matrix::dimension");
    }
  }

  void swap(Matrix<T>&);
   void swaporcopy(Matrix<T>&); //!< as ::swap but will copy if objects are nondynamic

  BaseList<T> row() const {
    if (empty())
      throw Undefined("Matrix<T>::row"); 
    return store_;
    //    return BaseList<T>(store_.size(),store_.vector());
  }

#ifndef LCM_SUPPRESS_VIEWS
  Matrix<T> rows(const range& ra) {
    if (ra.max()>=r_)
      throw BadIndex("rows",ra.max(),r_);
    return Matrix<T>(ra.size(),c_,store_.vector()+ra.start()*c_,mxflag::nondynamic);
  }
  const Matrix<T> rows(const range& ra) const {
    if (ra.max()>=r_)
      throw BadIndex("rows",ra.max(),r_);
    return Matrix<T>(ra.size(),c_,store_.vector()+ra.start()*c_,mxflag::nondynamic);
  }

  inline IndirectList<T,slice> diag(int =0);
  inline const IndirectList<T,slice> diag(int =0) const;

  template<typename Map1,typename Map2> typename LCM_Matrix_Select<T,LCM_DIM(Map1),Map1,LCM_DIM(Map2),Map2>::return_type operator() (const Map1& map1,const Map2& map2) 
  { return LCM_Matrix_Select<T,LCM_DIM(Map1),Map1,LCM_DIM(Map2),Map2>::func(*this,map1,map2); }

  template<typename Map1,typename Map2> typename LCM_Matrix_const_Select<T,LCM_DIM(Map1),Map1,LCM_DIM(Map2),Map2>::return_type operator() (const Map1& map1,const Map2& map2) const
  { return LCM_Matrix_const_Select<T,LCM_DIM(Map1),Map1,LCM_DIM(Map2),Map2>::func(*this,map1,map2); }

#endif
  //LCM_SUPPRESS_VIEWS

#if (BOUNDS_CHECK)
   T& operator()(size_t row_,size_t col_) {
     if (row_>=r_ || col_>=c_)
       throw BadIndex("Matrix<T>()",row_,r_,col_,c_); 
     return store_[row_*c_+col_];
   }
   const T& operator()(size_t row_,size_t col_) const {
     if (row_>=r_ || col_>=c_)
       throw BadIndex("Matrix<T>()",row_,r_,col_,c_);
     return store_[row_*c_+col_];
   }
   BaseList<T> row(size_t index_) const {
     if (index_>=r_)
       throw BadIndex("Matrix<T>::row",index_,r_);
     return BaseList<T>(c_,store_.vector()+index_*c_);
   }     
#else
   BaseList<T> row(size_t index_) const {
     return BaseList<T>(c_,store_.vector()+index_*c_);
   }
   T& operator()(size_t r,size_t c) {
     return store_[r*c_+c];
   }
   const T& operator()(size_t r,size_t c) const {
     return store_[r*c_+c];
   }
#endif

  void allrows(List<BaseList<T> >&);
  inline List<BaseList<T> > allrows()
    { List<BaseList<T> > res(mxflag::temporary); 
      allrows(res);
      return res; }

  friend bool isdefined <> (const Matrix<T>&);
  friend class BlockedMatrix<T>;

  template<class T2> LCM_INLINE Matrix<T>& operator+= (const T2& b) { return add_ip(*this,b); }
  template<class T2> LCM_INLINE Matrix<T>& operator-= (const T2& b) { return subtract_ip(*this,b); }
  template<class T2> LCM_INLINE Matrix<T>& operator*= (const T2& b) { return multiply_ip(*this,b); }
  template<class T2> LCM_INLINE Matrix<T>& operator&= (const T2& b) { return premultiply_ip(*this,b); }
  //only divide operation is divide by constant
  inline Matrix<T>& operator/= (const T&);
  template<class T2> LCM_INLINE void emultiply(T2& b) { emultiply_ip(*this,b); }
  
  static size_t total_allocated_words() { return matrix_traits<T>::allocator::allocated_words(); }

  void identity(int);
  void identity();
  void transpose();

  void conj(); //no definition provided.  Link will fail if not defined
  void conj_transpose();
  void simtrans(const Matrix<T>&);
  void isimtrans(const Matrix<T>&);
  void unitary_simtrans(const Matrix<T>&, Matrix<T>* =NULL);
  void unitary_isimtrans(const Matrix<T>&, Matrix<T>* =NULL);
  void unitary_simtrans(const BaseList<T>&);
  void unitary_isimtrans(const BaseList<T>&);

   class base_permuted_iterator_ {
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
     base_permuted_iterator_(const Matrix<T>& a,int r,int c) : data_(const_cast<T*>(a.vector())), rows_(a.rows()), cols_(a.cols()) { reset(r,c); }
     bool operator== (const base_permuted_iterator_& x) const { return (x.c_==c_) && (x.r_==r_); }
     bool operator!= (const base_permuted_iterator_& x) const { return (x.c_!=c_) || (x.r_!=r_); }
   };

   struct permuted_iterator : public base_permuted_iterator_, public ::std::iterator< ::std::bidirectional_iterator_tag,T> {
     using base_permuted_iterator_::r_;
     using base_permuted_iterator_::cols_;
     using base_permuted_iterator_::ccol;
     permuted_iterator(Matrix<T>& a,size_t r,size_t c) : base_permuted_iterator_(a,r,c) {}
     T& operator*() const { return ccol[r_*cols_]; }
     permuted_iterator& operator++() { this->up(); return *this; }
     permuted_iterator operator++(int) { permuted_iterator tmp(*this); this->up(); return tmp; }
     permuted_iterator& operator--() { this->down(); return *this; }
     permuted_iterator operator--(int) { permuted_iterator tmp(*this); this->down(); return tmp; }
   };

   struct const_permuted_iterator : public base_permuted_iterator_, public ::std::iterator< ::std::bidirectional_iterator_tag,T,ptrdiff_t,const T*,const T&> {
     using base_permuted_iterator_::r_;
     using base_permuted_iterator_::cols_;
     using base_permuted_iterator_::ccol;
     const_permuted_iterator(const Matrix<T>& a,size_t r,size_t c) : base_permuted_iterator_(a,r,c) {}
     const T& operator*() const { return ccol[r_*cols_]; }
     const_permuted_iterator& operator++() { this->up(); return *this; }
     const_permuted_iterator operator++(int) { const_permuted_iterator tmp(*this); this->up(); return tmp; }
     const_permuted_iterator& operator--() { this->down(); return *this; }
     const_permuted_iterator operator--(int) { const_permuted_iterator tmp(*this); this->down(); return tmp; }
   };

   permuted_iterator permuted_begin() { return permuted_iterator(*this,0,0); }
   permuted_iterator permuted_end() { return permuted_iterator(*this,r_,0); }
   const_permuted_iterator permuted_begin() const { return const_permuted_iterator(*this,0,0); }
   const_permuted_iterator permuted_end() const { return const_permuted_iterator(*this,r_,0); }

private:

  void dump() const { ::std::cout << *this; }

   DynamicList<T,typename matrix_traits<T>::allocator> store_;
   size_t r_, c_;
};

 template<typename T> Matrix<T>::Matrix(const Matrix<T>& a) : store_(a.store_)
{
  r_=a.r_;
  c_=a.c_;
 }

   template<typename T> Matrix<T>::Matrix(const Matrix<T>& a, mxflag::tempflag tflag) : store_(a.store_,tflag)
{
  r_=a.r_;
  c_=a.c_;
 }

   template<typename T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& a)
     {
       if (this != &a) {
	 store_=a.store_;
	 r_=a.r_;
	 c_=a.c_;
       }
       return *this;
     }
   
//new object rules
template<class T> struct new_trait_<T,2,2> { typedef Matrix<T> value_type; };
template<class T> struct new_trait_<T,2,0> { typedef Matrix<T> value_type; };
template<class T> struct new_trait_<T,0,2> { typedef Matrix<T> value_type; };

 template<class T> std::ostream& operator<< (std::ostream& ostr, const Matrix<T>& a) {
#ifndef NDEBUG
   if (a.istemporary())
     ostr << "[temporary]\n";
   else {
     if (!a.isdynamic()) 
       ostr << "[nondynamic]\n";
   }
#endif
   return print(a,ostr);
 }

  class range;
 
  //! minimal matrix definition for interfacing with BLAS level routines
  template<class T> class RawMatrix {
  public:
    //! create reference from full matrix type
    //* NB allow implicit conversion from full matrix type
    RawMatrix(const Matrix<T>& a)  
      : r_(a.rows()), c_(a.cols()), data_(const_cast<T*>(a.vector()))
    { step_=c_; }

    RawMatrix(const Matrix<T>&, const range&); //!< select diagonal submatrix
    RawMatrix(const Matrix<T>& a, const range&,const range&);
    RawMatrix(size_t rv, size_t cv, T* datav, size_t stepv =0);

    void create(size_t r, size_t c) {
      if ((r!=r_) || (c!=c_))
	throw Mismatch("RawMatrix::create",r,c,r_,c_);
    }
    static bool empty() { return false; }
    size_t rows() const { return r_; }
    size_t cols() const { return c_; }
    size_t step() const { return step_; }
    template<class T2> RawMatrix<T>& operator= (const T2& v);
    // bool iscontiguous() const { return (step_==c_); }

    const T* vector() const { return data_; }
    T* vector() { return data_; }

   T& operator()(size_t row_,size_t col_) {
#if (BOUNDS_CHECK)
     if (row_>=r_ || col_>=c_)
       throw BadIndex("Matrix<T>()",row_,r_,col_,c_); 
#endif
     return data_[row_*step_+col_];
   }
   const T& operator()(size_t row_,size_t col_) const { 
#if (BOUNDS_CHECK)
     if (row_>=r_ || col_>=c_)
       throw BadIndex("Matrix<T>()",row_,r_,col_,c_);
#endif
     return data_[row_*step_+col_];
   }

   class base_iterator_ {
   protected:
     RawMatrix<T>& data_;
     T* crow;
     size_t r_;
     int c_;
   
     void up() {
       if (++c_==data_.c_) 
	 reset(r_+1,0);
     }
     void down() {
       if (--c_<0)
	 reset(r_-1,data_.c_-1);
     }
     void reset(size_t r,int c) {
       r_=r; c_=c;
       crow=data_.data_+r_*data_.step_;
     } 

   public:
     base_iterator_(RawMatrix<T>& a,size_t r,size_t c)
       : data_(a) { reset(r,c); }
     bool operator== (const base_iterator_& x) const { return (x.c_==c_) && (x.r_==r_); }
     bool operator!= (const base_iterator_& x) const { return (x.c_!=c_) || (x.r_!=r_); }
   };

   friend class base_iterator_;

   struct iterator : public base_iterator_, ::std::iterator< ::std::bidirectional_iterator_tag,T> {
     using base_iterator_::crow;
     using base_iterator_::c_;
     iterator(RawMatrix<T>& a,size_t r,size_t c) : base_iterator_(a,r,c) {}
     T& operator*() const { return crow[c_]; }
     iterator& operator++() { this->up(); return *this; }
     iterator operator++(int) { iterator tmp(*this); this->up(); return tmp; }
     iterator& operator--() { this->down(); return *this; }
     iterator operator--(int) { iterator tmp(*this); this->down(); return tmp; }
   };

   struct const_iterator : public base_iterator_, ::std::iterator< ::std::bidirectional_iterator_tag,T,ptrdiff_t,const T*,const T&> {
     using base_iterator_::crow;
     using base_iterator_::c_;
     const_iterator(const RawMatrix<T>& a,size_t r,size_t c) : base_iterator_(const_cast< RawMatrix<T>& >(a),r,c) {}
     const T& operator*() const { return crow[c_]; }
     const_iterator& operator++() { this->up(); return *this; }
     const_iterator operator++(int) { const_iterator tmp(*this); this->up(); return tmp; }
     const_iterator& operator--() { this->down(); return *this; }
     const_iterator operator--(int) { const_iterator tmp(*this); this->down(); return tmp; }
   };

    iterator begin() { return iterator(*this,0,0); }
    iterator end() { return iterator(*this,r_,0); }
    const_iterator begin() const { return const_iterator(*this,0,0); }
    const_iterator end() const { return const_iterator(*this,r_,0); }
    
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;  

  private:
    size_t r_,c_,step_;
    T* data_;

    const BaseList<T> row() const { return BaseList<T>(r_*c_,data_); }
    BaseList<T> row() { return BaseList<T>(r_*c_,data_); }
  };

 template<typename T> struct type_traits< RawMatrix<T> > {
   static const int rank=0;
   static const bool trivialconstructor=false;
   static const size_t dimensionality=2;
   typedef T value_type;
 };

  template<class T> template<class T2> RawMatrix<T>& RawMatrix<T>::operator= (const T2& v)
  {
    assign(*this,v);
//     if (!iscontiguous())
//       throw Failed("RawMatrix<T>::= not supported on non-contiguous matrices");
//     BaseList<T> asrow(row());
//     asrow=v;
    return *this;
  }

  template<class T> RawMatrix<T>::RawMatrix(size_t rv, size_t cv, T* datav, size_t stepv)
    : r_(rv), c_(cv), step_(stepv), data_(datav)
  {
    if (step_) {
      if (step_<c_)
	throw InvalidParameter("RawMatrix: data step is less than number of columns!");
    }
    else
      step_=c_;
  }
  
template<class T> inline bool issquare(const T& a)
{
  LCM_STATIC_CHECK( LCM_DIM(T)==2, issquare_applied_to_non2D_object );
    if (!a)
      throw Undefined("issquare");
    return (a.rows()==a.cols()); 
  }

 template<class T1,class T2> inline bool issimple(const T1& a, const T2&,size_t len) { return (a.cols()<len); }
// template<class T1,class T2> inline bool issimple(const T1& a, const BaseList<T2>&,size_t len) { return (a.cols()<len); }

template<typename T1,typename T2,typename T3> void multiply_transpose(Matrix<T1>&, const BaseList<T2>&, const Matrix<T3>&);
template<typename T1,typename T2,typename T3> void transpose_multiply(Matrix<T1>&, const Matrix<T2>&, const BaseList<T3>&);
template<typename T1,typename T2,typename T3> void transpose_multiply(Matrix<T1>&, const Matrix<T2>&, const Matrix<T3>&);
template<typename T1,typename T2,typename T3> void transpose_multiply(BaseList<T1>&, const Matrix<T2>&, const BaseList<T3>&);

  template<class T> inline bool isdefined(const Matrix<T>& a) { return (!!a); }

 void print_cmatrix(std::ostream& ostr, const Matrix<complex>&);

 template<> struct Print_<Matrix<complex>,2> {
   static inline void print(const Matrix<complex>& a, std::ostream& ostr)
   { print_cmatrix(ostr,a); }
 };

template <> inline Matrix<float_t>& Matrix<float_t>::operator/= (const float_t& a) { return ((*this)*=1.0/a); }
 template<> void Matrix<complex>::conj();
 template<> void Matrix<complex>::conj_transpose();
 template<> void Matrix<complex>::simtrans(const Matrix<complex>&);
 template<> void Matrix<complex>::isimtrans(const Matrix<complex>&);
 template<> void Matrix<complex>::unitary_simtrans(const Matrix<complex>&, Matrix<complex>*);
 template<> void Matrix<complex>::unitary_isimtrans(const Matrix<complex>&, Matrix<complex>*);
 template<> void Matrix<float_t>::unitary_simtrans(const Matrix<float_t>&, Matrix<double>*);
 template<> void Matrix<float_t>::unitary_isimtrans(const Matrix<float_t>&, Matrix<double>*);
 template<> void Matrix<complex>::unitary_simtrans(const BaseList<complex>&);
 template<> void Matrix<complex>::unitary_isimtrans(const BaseList<complex>&);

 template<class T,class T2> inline bool operator== (const Matrix<T>& a, const T2& b) { return areequal(a,b); }
 template<class T,class T2> inline bool operator!= (const Matrix<T>& a, const T2& b) { return arenotequal(a,b); }

 //template<class T> inline T sum(const Matrix<T>& a) { return sum(a.row()); }
template<class T> T trace(const Matrix<T>&);
template<class T1,class T2> void diag(T1&, const T2&, int =0); 
template<class T1,class T2> void full(Matrix<T1>&, const T2&);
template<class T> inline List< LCM_VAL(T) > diag(const T& a) { List< LCM_VAL(T) > d(a.rows(),mxflag::temporary); diag(d,a); return d;}
template<class T> inline Matrix< LCM_VAL(T) > full(const T& a) { Matrix<LCM_VAL(T) > d(mxflag::temporary); full(d,a); return d; }

template<class T> Matrix<T> transpose(const Matrix<T>&);
template<class T> Matrix<T> operator- (const Matrix<T>&);
template<class T> Matrix<T> operator/ (const Matrix<T>&, const T&);

template<class T> void row(List< BaseList<T> >&ind,const BaseList< Matrix<T> >&a,size_t row)
{
  size_t rs=a.size();
  ind.create(rs);
  for (;rs--;)
    ind(rs)=a(rs).row(row);
}

template<class T1,class T2> bool aretranspose(const T1& a,const T2& b) { 
  return (!a.empty()) && (a.rows()==b.cols()) && (a.cols()==b.rows());
}
  
// Implementation details

//! return first bad index found, 0 if none (assumes map is non-empty)
template<typename Map> size_t find_bad_index(const Map& map,size_t maxv)
{
  for (size_t i=map.size();i--;) {
    if (map(i)>=maxv)
      return map(i);
  }
  return 0;
}

#if (BOUNDS_CHECK)
 template<typename Map> void MatrixHelper_::check_bounds1(const Map& map,size_t maxv)
   {
     const size_t badi=find_bad_index(map,maxv);
     if (badi)
       throw BadIndex("Matrix<T>::check_bounds",badi,maxv);
   }
 void MatrixHelper_::check_bounds0(size_t i,size_t maxv) {
   if (i>=maxv)
     throw BadIndex("Matrix<T>::check_bounds",i,maxv);
 }
#ifndef LCM_SUPPRESS_VIEWS
 void MatrixHelper_::check_boundssl(const slice& r,size_t maxv) { 
   if (r.max()>=maxv) 
     throw BadIndex("Matrix<T>::check_bounds",r.max(),maxv);
 }
#endif
#else
 template<typename Map> void MatrixHelper_::check_bounds1(const Map&,size_t) {}
 void MatrixHelper_::check_bounds0(size_t,size_t) {}
#ifndef LCM_SUPPRESS_VIEWS
 void MatrixHelper_::check_boundssl(const slice&,size_t) {}
#endif
#endif

 template<class T> Matrix<T>& Matrix<T>::operator=(const T& v) {
   T* d=store_.vector(); 
   for (size_t i=size();i--;)
     *d++=v; 
   return *this;
 }
 

template<class T1,class T2> void transpose(T1& d,const T2& a)  
{
  LCM_STATIC_CHECK( (LCM_DIM(T1)==2) && (LCM_DIM(T2)==2) , Transpose_arguments_must_be_2D );
  if (!d.empty() && issame(d,a))
    throw ArgumentClash("transpose");

  size_t ar=a.rows();
  size_t ac=a.cols();

  d.create(ac,ar);
  for (;ar--;)
    for (size_t c=ac;c--;)
      d(c,ar)=a(ar,c);
}

template<class T1,class T2> void conj_transpose_(T1& d,const T2& a, Type2Type<double>)
{
  transpose(d,a);
}

template<class T1,class T2> void conj_transpose_(T1& d,const T2& a, Type2Type<complex>)  
{
  LCM_STATIC_CHECK( (LCM_DIM(T1)==2) && (LCM_DIM(T2)==2) , Conj_transpose_arguments_must_be_2D );
  if (!d.empty() && issame(d,a))
    throw ArgumentClash("conj_transpose");

  size_t ar=a.rows();
  size_t ac=a.cols();

  d.create(ac,ar);
  for (;ar--;)
    for (size_t c=ac;c--;)
      d(c,ar)=conj(a(ar,c));
}

template<class T1,class T2> void conj_transpose(T1& d,const T2& a)  
{
  conj_transpose_(d,a,Type2Type< LCM_VAL(T2) >());
}

// template<class T> void transpose(Matrix<T>& d,const Matrix<T>& a)
// {
//   _transpose_blocked(d,a);
// }

template<class T> Matrix<T> transpose(const Matrix<T>& a)
{
  Matrix<T> d(mxflag::temporary);
  transpose(d,a);
  return d;
}

  template<class T1,class T2,class T3> inline void multiply_naive_ref_(RawMatrix<T1>& dest,const RawMatrix<T2>& a,const RawMatrix<T3>& b, bool acc)
{
  const size_t ca=a.cols();
  const size_t cb=b.cols();
  const size_t ra=a.rows();
  const size_t idc=dest.step();
  const size_t iac=a.step();
  const size_t ibc=b.step();

  if (ca!=b.rows())
    throw Mismatch("multiply",ra,ca,b.rows(),cb);

  T1 *destp=dest.vector();
  const T2 *nad=a.vector();
  const T3 *bd=b.vector();

  for (size_t i=0;i<ra;i++) {
    for (size_t j=0;j<cb;j++) {
      T1& sum(destp[j]);
      if (!acc)
	sum=T1(0);
      const T3 *bp=bd+j;
      for (size_t k=0;k<ca;k++,bp+=ibc)
	mla(sum,nad[k],*bp);
    }
    nad+=iac;
    destp+=idc;
  }
}

  template<class T1,class T2,class T3> inline void multiply_naive_direct_(RawMatrix<T1>& dest,const RawMatrix<T2>& a,const RawMatrix<T3>& b, bool acc)
{
  const size_t ca=a.cols();
  const size_t cb=b.cols();
  const size_t ra=a.rows();
  const size_t idc=dest.step();
  const size_t iac=a.step();
  const size_t ibc=b.step();

  if (ca!=b.rows())
    throw Mismatch("multiply",ra,ca,b.rows(),cb);

  T1 *destp=dest.vector();
  const T2 *nad=a.vector();
  const T3 *bd=b.vector();

  for (size_t i=0;i<ra;i++) {
    for (size_t j=0;j<cb;j++) {
      T1 sum(acc ? destp[j] : T1(0));
      const T3 *bp=bd+j;
      for (size_t k=0;k<ca;k++,bp+=ibc)
	mla(sum,nad[k],*bp);
      destp[j]=sum;
    }
    nad+=iac;
    destp+=idc;
  }
}

  template<class T1,class T2,class T3> inline void multiply_ref_(RawMatrix<T1>& dest,const RawMatrix<T2>& a,const RawMatrix<T3>& b, bool acc)
{
  const size_t ca=a.cols();
  const size_t cb=b.cols();
  const size_t ra=a.rows();

  if (ca!=b.rows())
    throw Mismatch("multiply",ra,ca,b.rows(),cb);

  if (!acc)
    dest=T1(0);

  doesmla<T1,T2,T3> mlaobj;

  const size_t idc=dest.step();
  const size_t iac=a.step();
  const size_t ibc=b.step();
  const T3* const bp=b.vector();
  T1* dp=dest.vector();
  const T2* ap=a.vector();

  size_t ib=0;
  while (ib<ra) {
    size_t fini=ib+LCM_BLOCK_FACTOR;
    if (fini>ra)
      fini=ra;
    
    size_t kb=0;
    while (kb<ca) {
      size_t fink=kb+LCM_BLOCK_FACTOR;
      if (fink>ca)
	fink=ca;

      const size_t duffk=(fink-kb) & 3;
      const size_t thisfink=fink-duffk;
      
      size_t jb=0;
      while (jb<cb) {
	size_t finj=jb+LCM_BLOCK_FACTOR;
	if (finj>cb)
	  finj=cb;

	const T3 *bpstart=bp+kb*ibc;
	
	for (size_t i=ib;i<fini;i++) {
	  T1* dptri=dp+i*idc;
	  const T2* aptri=ap+i*iac;
	  for (size_t j=jb;j<finj;j++) {
	    T1& sum=dptri[j];
	    const T3* lbp=bpstart+j;
	    
	    size_t k=kb;
	    for (;k<thisfink;) {
	      mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	      mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	      mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	      mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	    }
	    switch (duffk) {
	    case 3: mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	    case 2: mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	    case 1: mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	    }
	  }
	}
	jb=finj;
      }
      kb=fink;
    }
    ib=fini;
  }
}

  template<class T1,class T2,class T3> inline void multiply_direct_(RawMatrix<T1>& dest,const RawMatrix<T2>& a,const RawMatrix<T3>& b, bool acc)
{
  const size_t ca=a.cols();
  const size_t cb=b.cols();
  const size_t ra=a.rows();

  if (ca!=b.rows())
    throw Mismatch("multiply",a.rows(),a.cols(),b.rows(),b.cols());
  if ((dest.rows()!=ra) || (dest.cols()!=cb))
    throw Mismatch("multiply: destination has wrong size");

  if (!acc)
    dest=T1(0);

  const size_t idc=dest.step();
  const size_t iac=a.step();
  const size_t ibc=b.step();
  const T3* const bp=b.vector();
  T1* dp=dest.vector();
  const T2* ap=a.vector();

  doesmla<T1,T2,T3> mlaobj;

  size_t ib=0;
  while (ib<ra) {
    size_t fini=ib+LCM_BLOCK_FACTOR;
    if (fini>ra)
      fini=ra;
    
    size_t kb=0;
    while (kb<ca) {
      size_t fink=kb+LCM_BLOCK_FACTOR;
      if (fink>ca)
	fink=ca;

      const size_t duffk=(fink-kb) & 3;
      const size_t thisfink=fink-duffk;

      size_t jb=0;
      while (jb<cb) {
	size_t finj=jb+LCM_BLOCK_FACTOR;
	if (finj>cb)
	  finj=cb;

	const T3* LCM_RESTRICT bpstart=bp+kb*ibc;

	for (size_t i=ib;i<fini;i++) {
	  T1* LCM_RESTRICT dptri=dp+i*idc;
	  const T2* LCM_RESTRICT aptri=ap+i*iac;
	  for (size_t j=jb;j<finj;j++) {
	    T1 sum=dptri[j];
	    T1 sum2=0.0;
	    const T3* LCM_RESTRICT lbp=bpstart+j;
	    
	    size_t k=kb;
	    for (;k<thisfink;) {
	      mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	      mlaobj(sum2,aptri[k],*lbp); k++; lbp+=ibc;
	      mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	      mlaobj(sum2,aptri[k],*lbp); k++; lbp+=ibc;
	    }
	    switch (duffk) {
	    case 3: mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	    case 2: mlaobj(sum2,aptri[k],*lbp); k++; lbp+=ibc;
	    case 1: mlaobj(sum,aptri[k],*lbp); k++; lbp+=ibc;
	    }
	    dptri[j]=sum+sum2;
	  }
	}
	jb=finj;
      }
      kb=fink;
    }
    ib=fini;
  }
}

//! post-checking entry point
template<class T1,class T2,class T3> void multiply(RawMatrix<T1>& d, const RawMatrix<T2>& a,const RawMatrix<T3>& b, bool acc =false)
{
  if (issimple(a,b,LCM_NAIVE_BELOW))
    multiply_naive_direct_(d,a,b,acc);
  else
    multiply_direct_(d,a,b,acc);
}

// template<class T1,class T2,class T3> void multiply(RawMatrix<T1>& d, const RawMatrix<T2>& a,const RawMatrix<T3>& b, bool acc =false)
// {
//   if (issame(d,a) || issame(d,b))
//     throw ArgumentClash("multiply");
//   if ((d.rows()!=a.rows()) || (d.cols()!=b.cols()))
//     throw Mismatch("multiply: source and destination",a.rows(),b.cols(),d.rows(),d.cols());
//   multiply_(d,a,b,acc);
// }

  template<class T1,class T2,class T3> void multiply(Matrix<T1>& d, const Matrix<T2>& a,const Matrix<T3>& b)
  {
    if (!!d && (issame(d,a) || issame(d,b)))
      throw ArgumentClash("multiply");

    d.create(a.rows(),b.cols());
    RawMatrix<T1> draw(d);
    multiply(draw,RawMatrix<T2>(a),RawMatrix<T3>(b),false);
  }

//#if (LCM_COMPLEX_BYREF!=0)
  template<class T2,class T3> void multiply(RawMatrix<complex>& d, const RawMatrix<T2>& a,const RawMatrix<T3>& b, bool acc)
{
  if (issimple(a,b,LCM_NAIVE_BELOW))
    multiply_naive_ref_(d,a,b,acc);
  else
    multiply_ref_(d,a,b,acc);
}

//#endif

void cmatrix_multiply(RawMatrix<complex>&, const RawMatrix<complex>&, const RawMatrix<complex>&, bool);
void cmatrix_multiply(Matrix<complex>&, const RawMatrix<complex>&, const RawMatrix<complex>&, bool);

void cmatrix_multiply(RawMatrix<double>&, const RawMatrix<double>&, const RawMatrix<double>&,bool);
void cmatrix_multiply(Matrix<double>& d, const RawMatrix<double>& a, const RawMatrix<double>& b, bool acc);

#ifdef LCM_USE_EXTERNAL
void lapack_multiply(RawMatrix<complex>&, const RawMatrix<complex>&, const RawMatrix<complex>&, bool, bool, bool);
inline void lapack_multiply(RawMatrix<complex>& d, const RawMatrix<complex>& a, const RawMatrix<complex>& b, bool acc) { lapack_multiply(d,a,b,acc,false,false); }

void lapack_multiply(RawMatrix<double>&, const RawMatrix<double>&, const RawMatrix<double>&, bool, bool, bool);
inline void lapack_multiply(RawMatrix<double>& d, const RawMatrix<double>& a, const RawMatrix<double>& b, bool acc) { lapack_multiply(d,a,b,acc,false,false); }

  template<> inline void multiply(RawMatrix<complex>& d, const RawMatrix<complex>& a, const RawMatrix<complex>& b, bool acc)
{
  if (issimple(a,b,LCM_INTERNAL_ZMM))
    cmatrix_multiply(d,a,b,acc);
  else
    lapack_multiply(d,a,b,acc);
}

  template<> inline void multiply(RawMatrix<double>& d, const RawMatrix<double>& a, const RawMatrix<double>& b, bool acc)
{
  if (issimple(a,b,LCM_INTERNAL_DMM))
    cmatrix_multiply(d,a,b,acc);
  else
    lapack_multiply(d,a,b,acc);
}
#else
  template<> inline void multiply(RawMatrix<double>& d, const RawMatrix<double>& a, const RawMatrix<double>& b, bool acc) {  cmatrix_multiply(d,a,b,acc); }

  template<> inline void multiply(RawMatrix<complex>& d, const RawMatrix<complex>& a, const RawMatrix<complex>& b, bool acc) { cmatrix_multiply(d,a,b,acc); }

#endif

template<class T> T trace(const Matrix<T>& a)
{
  if (!issquare(a))
    throw NotSquare("trace");

  size_t i=a.rows()-1;  
  T sum=a(i,i);
  for (;i--;)
    sum+=a(i,i);
  return sum;
}

template<class T1,class T2> inline typename new_trait<Matrix<T1>,T2,true>::value_type
operator+ (const Matrix<T1>& a,const T2& b)
 { 
   typename new_trait<Matrix<T1>,T2,true>::value_type d(mxflag::temporary);
   add(d,a,b);
   return d;
 }

template<class T1,class T2> inline typename new_trait<Matrix<T1>,T2,true>::value_type
operator- (const Matrix<T1>& a,const T2& b)
 { 
   typename new_trait<Matrix<T1>,T2,true>::value_type d(mxflag::temporary);
   subtract(d,a,b);
   return d;
 }

template<class T1,class T2> inline typename new_trait<Matrix<T1>,T2,false>::value_type
operator| (const Matrix<T1>& a,const T2& b)
 { 
   typedef typename new_trait<Matrix<T1>,T2,false>::value_type return_type;
   return_type d(mxflag::temporary);
   apply(d,doesor<LCM_VAL(return_type),T1,LCM_VAL(T2)>(),a,b);
   return d;
 }

template<class T1,class T2> inline typename new_trait<Matrix<T1>,T2,false>::value_type
operator& (const Matrix<T1>& a,const T2& b)
 { 
   typedef typename new_trait<Matrix<T1>,T2,false>::value_type return_type;
   return_type d(mxflag::temporary);
   apply(d,doesand<LCM_VAL(return_type),T1,LCM_VAL(T2)>(),a,b);
   return d;
 }

template<class T> Matrix<T> operator/ (const Matrix<T>& a,const T& b)
{
  Matrix<T> d(a,mxflag::temporary);
  return d/=b;
}

template<class T> Matrix<T> operator+ (const T& b,const Matrix<T>& a)
{
  Matrix<T> d(a,mxflag::temporary);
  return d+=b;
}

template<class T> Matrix<T> operator- (const T& b,const Matrix<T>& a)
{
  Matrix<T> d(a,mxflag::temporary);
  return d-=b;
}

template<class T> Matrix<T> operator* (const Matrix<T>& a,const T& b)
{
  Matrix<T> d(b,mxflag::temporary);
  return d*=a;
}

 template<typename T1,typename T2> struct doespremultiply_ip<Matrix<T1>,Matrix<T2> > : public unary_function_ip<Matrix<T1>,Matrix<T2> > {
  inline void operator()(Matrix<T1>& d, const Matrix<T2>& a) const {
    Matrix<T1> tmp;
    multiply(tmp,a,d);
    d.swap(tmp);
  }
};

 template<typename T1,typename T2> struct doespremultiply_ip<Matrix<T1>,BaseList<T2> > : public unary_function_ip<Matrix<T1>,BaseList<T2> > {
  inline void operator()(Matrix<T1>& d, const BaseList<T2>& a) const {
    multiply(d,a,d);
  }
};

 template<typename T1,typename T2> struct doesmultiply_ip<Matrix<T1>,Matrix<T2> > : public unary_function_ip<Matrix<T1>,Matrix<T2> > {
  inline void operator()(Matrix<T1>& d, const Matrix<T2>& a) const {
    Matrix<T1> tmp;
    multiply(tmp,d,a);
    d.swap(tmp);
  }
};

 template<typename T1,typename T2> struct doesmultiply_ip<Matrix<T1>,BaseList<T2> > : public unary_function_ip<Matrix<T1>,Matrix<T2> > {
  inline void operator()(Matrix<T1>& d, const BaseList<T2>& a) const {
    multiply(d,d,a);
  }
};

/* N.B. Mixed operations assume smaller object second for add/subtract but first for multiplication i.e. 2*A but A+2.
   We define both <op>(T2 b,Matrix<T1> a) and <op>(Matrix<T1> a,T2 b) since this would be ambiguous */

template<class T1,class T2> inline LCM_NEWOBJECT(Matrix<T1>,T2) operator* (const T2& b,const Matrix<T1>& a)
 { 
   LCM_NEWOBJECT(Matrix<T1>,T2) d(mxflag::temporary);
   multiply(d,b,a);
   return d;
 }

template<class T1,class T2> Matrix<LCM_NEWTYPE(T1,T2,true)> emultiply(const Matrix<T1>& a,const Matrix<T2>& b)
{
  Matrix<LCM_NEWTYPE(T1,T2,true)> d(mxflag::temporary);
  emultiply(d,a,b);
  return d;
}

// unitary minus
template<class T> Matrix<T> operator- (const Matrix<T>& a)
{
  Matrix<T> d(mxflag::temporary);
  negate(d,a);
  return d;
}

template<class Td,class Ta,class Tb> void cmatrix_multiply(BaseList<Td> dest,const Matrix<Ta>& a,const BaseList<Tb>& b)
{
  if (!a)
    throw Undefined("multiply");
  const size_t cs=a.cols();
  const size_t len=dest.size();
  if ( (cs!=b.size()) || (len!=a.rows()))
    throw Mismatch("multiply",a.rows(),a.cols(),b.size(),b.size());
  if (issame(dest,b))
    throw ArgumentClash("multiply");

  Td* destp=dest.vector();
  const Ta* amr=a.vector();
  const Tb* source=b.vector();
  
  for (size_t n=len;n--;) {
    Td& d=*destp++;
    size_t m=cs-1;
    d=source[m]*amr[m];
    for (;m--;)
      mla(d,source[m],amr[m]);
    amr+=cs;
  }
}

template<class Td,class Ta,class Tb> void cmatrix_multiply(BaseList<Td> dest,const BaseList<Ta>& a,const Matrix<Tb>& b)
{
  if (!b)
    throw Undefined("multiply");
  const size_t rs=b.rows();
  const size_t len=dest.size();
  if ( (rs!=a.size()) || (len!=b.cols()))
    throw Mismatch("multiply",a.size(),a.size(),b.rows(),b.cols());
  if (issame(dest,b))
    throw ArgumentClash("multiply");

  for (size_t n=len;n--;) {
    Td& d=dest(n);
    size_t m=rs-1;
    multiply(d,a(m),b(m,n));
    for (;m--;)
      mla(d,a(m),b(m,n));
  }
}

template<class Td,class Ta,class Tb> inline void multiply(BaseList<Td> dest,const Matrix<Ta>& a,const BaseList<Tb>& b) { cmatrix_multiply(dest,a,b); }

template<class Td,class Ta,class Tb> inline void multiply(BaseList<Td> dest,const BaseList<Ta>& a,const Matrix<Tb>& b) { cmatrix_multiply(dest,a,b); }

#ifdef LCM_USE_EXTERNAL
//dgemv specialisation
void lapack_multiply(BaseList<double>&, const Matrix<double>&, const BaseList<double>&, bool =false);
void lapack_multiply(BaseList<double>&, const BaseList<double>&, const Matrix<double>&, bool =false);
//zgemv specialisation
void lapack_multiply(BaseList<complex>&, const Matrix<complex>&, const BaseList<complex>&, bool =false);
void lapack_multiply(BaseList<complex>&, const BaseList<complex>&, const Matrix<complex>&, bool =false);

template<> inline void multiply(BaseList<double>& d, const Matrix<double>& a, const BaseList<double>& b)
{
  if (issimple(a,b,LCM_INTERNAL_DMV))
    cmatrix_multiply(d,a,b);
  else
    lapack_multiply(d,a,b);
}

template<> inline void multiply(BaseList<double>& d, const BaseList<double>& b, const Matrix<double>& a)
{
  //assume same cross-over for MV and VM
  if (issimple(a,b,LCM_INTERNAL_DMV))
    cmatrix_multiply(d,b,a);
  else
    lapack_multiply(d,a,b,true);
}

template<> inline void multiply(BaseList<complex>& d, const Matrix<complex>& a, const BaseList<complex>& b)
{
  if (issimple(a,b,LCM_INTERNAL_ZMV))
    cmatrix_multiply(d,a,b);
  else
    lapack_multiply(d,a,b);
}

template<> inline void multiply(BaseList<complex>& d, const BaseList<complex>& b, const Matrix<complex>& a)
{
  if (issimple(a,b,LCM_INTERNAL_ZMV))
    cmatrix_multiply(d,b,a);
  else
    lapack_multiply(d,a,b,true);
}
#endif

template<class Td,class Ta,class Tb> void multiply(List<Td>& dest,const BaseList<Tb>& b,const Matrix<Ta>& a)
{
  if (!a)
    throw Undefined("multiply");
  dest.create(a.cols());
  multiply( static_cast< BaseList<Td>& >(dest),b,a);
}

extern const char TM[];
extern const char MT[];

template<typename T1,typename T2,typename T3> void multiply_transpose(Matrix<T1>& to,const Matrix<T2>& a,const Matrix<T3>& T)
{
  if (!!to && (issame(to,a) || issame(to,T)))
    throw ArgumentClash(MT);
  const T3 *sTp=T.vector();
  const T2 *ap=a.vector();

  const size_t n=a.cols();
  if (n!=T.cols())
    throw Mismatch(MT,a.rows(),a.cols(),T.rows(),T.cols());

  const size_t ra=a.rows();
  const size_t rt=T.rows();
  to.create(ra,rt);

  T1 *destp=to.vector();

  for (size_t i=0;i<ra;i++) {
    const T3 *Tp=sTp;
    for (size_t k=0;k<rt;k++) {
      T1& sum=*destp++;
      size_t j=n-1;
      sum=ap[j]*Tp[j];
      for (;j--;) mla(sum,ap[j],Tp[j]);
      Tp+=rt;
    }
    ap+=n;
  }
}

template<typename T1, typename T2,typename T3> void multiply_transpose(Matrix<T1>& to,const BaseList<T2>& D,const Matrix<T3>& V)
{
  if (!V)
    throw Undefined(MT);
  if (!!to && issame(to,V))
    throw ArgumentClash(MT);

  const size_t n=D.size();
  if (V.cols()!=n)
    throw Mismatch(MT,D.size(),D.size(),V.rows(),V.cols());
  const size_t m=V.rows();

  to.create(n,m);
  T1 *destp=to.vector();

  for (size_t i=0;i<n;i++) {
    const T2& v=D(i);
    for (size_t k=0;k<m;k++) *destp++=V(k,i)*v;
  }
}

template<typename T1,typename T2,typename T3> void transpose_multiply(Matrix<T1>& to,const Matrix<T2>& T,const Matrix<T3>& a)
{
  if (!T)
    throw Undefined(TM);
  if (!!to && (issame(to,T) || issame(to,a)))
    throw ArgumentClash(TM);

  const size_t n=a.rows();
  if (n!=T.rows())
    throw Mismatch(TM,T.rows(),T.cols(),a.rows(),a.cols());

  const size_t ca=a.cols();
  const size_t ct=T.cols();
  to.create(ct,ca);
  T1 *destp=to.vector();

  for (size_t i=0;i<ct;i++) {
    for (size_t k=0;k<ca;k++) {
      T1 sum(0);
      for (size_t j=n;j--;) mla(sum,a(j,k),T(j,i));
      *destp++=sum;
    }
  }
}

template<typename T1,typename T2,typename T3> void transpose_multiply(Matrix<T1>& to,const Matrix<T2>& V,const BaseList<T3>& D)
{
  if (!V)
    throw Undefined(TM);
  if (!!to && issame(to,V))
    throw ArgumentClash(TM);

  const size_t n=D.size();
  if (V.rows()!=n)
    throw Mismatch(TM,V.rows(),V.cols(),D.size(),D.size());
  const size_t m=V.cols();

  to.create(m,n);

  for (size_t k=n;k--;) {
    const T3& v=D(k);
    const T2 *Vr=V.vector(k);
    for (size_t i=m;i--;) to(i,k)=Vr[i]*v;
  }
}

template<typename T1,typename T2,typename T3> void transpose_multiply(BaseList<T1>& to,const Matrix<T2>& T, const BaseList<T3>& a)
{
  if (!T)
    throw Undefined(TM);
  if (!!to && issame(to,a))
    throw ArgumentClash(TM);
  
  const size_t n=a.size();
  const size_t ct=T.cols();
  if (n!=T.rows() || ct!=to.size())
    throw Mismatch(TM,T.rows(),T.cols(),a.size(),a.size());
  T1 *destp=to.vector();
  
  for (size_t i=0;i<ct;i++) {
    T1 sum(0);
    for (size_t j=n;j--;) mla(sum,T(j,i),a(j));
    *destp++=sum;
  }
}

template<class T1,class T2> void diag(T1& d,const T2& a, int which)
{
  LCM_STATIC_CHECK( LCM_DIM(T2)==2, diag_destination_is_not_1D );
  if (!issquare(a))
    throw NotSquare("diag");
  size_t n=a.rows();
  if (::std::abs(which)>=n)
    throw BadIndex("diag");
  if (which<0) {
    n+=which;
    d.create(n);
    for (;n--;)
      d(n)=a(n-which,n);
  }
  else {
    n-=which;
    d.create(n);
    for (;n--;)
      d(n)=a(n,n+which);
  }
}

template<class T1,class T2> void full(Matrix<T1>& a,const T2& v)
{
  LCM_STATIC_CHECK( LCM_DIM(T2)==1, full_called_with_non_1D_argument );
  const size_t n=v.size();
  a.create(n,n,T1(0));
  for (size_t i=0;i<n;i++)
    a(i,i)=v(i);
}

 template<class T> class WithinTol {};

   template<> class WithinTol<double> : public ::std::unary_function<double,bool> {
   public:
     WithinTol(double tol_) : tol(tol_) {
       if (tol<0.0) throw InvalidParameter("WithinTol: tolerance can't be negative");
     }
     bool operator()(double x) const { return (fabs(x)<tol); }
   private:
     double tol;
   };

   template<> class WithinTol<complex> : public ::std::unary_function<complex,bool> {
   public:
     WithinTol(double tol_) : toltol(tol_) {
       if (tol_<0.0) throw InvalidParameter("WithinTol: tolerance can't be negative");
     }
     bool operator()(const complex& x) const { return (norm(x)<toltol); }
   private:
     double toltol;
   };
				     

template<class T> bool issymmetric(const Matrix<T>& a) {
  if (!issquare(a)) throw NotSquare("issymmetric");
  const size_t dim=a.rows();
  for (size_t r=1;r<dim;r++) {
    for (size_t c=r;c--;) {
      if (a(r,c)!=a(c,r)) return false;
    }
  }
  return true;
}

template<class T> bool issymmetric(const Matrix<T>& a, double tol) {
  if (!issquare(a)) throw NotSquare("issymmetric");
  const size_t dim=a.rows();
  WithinTol<T> compfunc(tol);
  for (size_t r=1;r<dim;r++) {
    for (size_t c=r;c--;) {
      if (!compfunc(a(r,c)-a(c,r))) return false;
    }
  }
  return true;
}

template<class T1,class T2> LCM_NEWTYPE(T1,T2,true) trace_multiply(const Matrix<T1>& b,const Matrix<T2>& c)
{
  const size_t cs=c.cols();
  const size_t rs=c.rows();

  if (b.cols()!=rs || b.rows()!=cs)
    throw Mismatch("trace_multiply",b.rows(),b.cols(),c.rows(),c.cols());

  LCM_NEWTYPE(T1,T2,true) sum(0);

  for (size_t i=cs;i--;) {
    const BaseList<T1> bi=b.row(i);
    for (size_t j=rs;j--;)
      mla(sum,bi(j),c(j,i));
  }
  return sum;
}

template<class T1,class T2> Matrix<T1> operator+ (const BaseList<T2>& b, const Matrix<T1>& a)
{
  Matrix<T1> d(a,mxflag::temporary);
  return d+=b;
}

template<class T1,class T2> Matrix<T1> operator- (const BaseList<T2>& b, const Matrix<T1>& a)
{
  Matrix<T1> d(mxflag::temporary);
  negate(d,a);
  return d+=b;
}

  //implementation details

  template<typename T> void Matrix<T>::create(int rows_,int cols_)
  {
    store_.create(rows_*cols_);
    r_=rows_; c_=cols_;
  }
  
  template<typename T> void Matrix<T>::resize(int rows_,int cols_)
  {
    const int new_items=rows_*cols_;
    const size_t old_items=store_.size();
    if (new_items!=old_items) {
      if (new_items<0)
	throw InvalidParameter("Matrix size cannot be negative");
      store_.reserve(new_items); //!< Note that contents are initialised
    }
    r_=rows_; c_=cols_;
  }

  
  template<typename T> template<typename T2> void Matrix<T>::push_back_row(const T2& a)
  {
    LCM_STATIC_CHECK( LCM_DIM(T2)==1, push_back_row_argument_must_be_one_dimensional );
    if (a.size()!=cols())
      throw Mismatch("push_back_row");
    resize(rows()+1,cols());
    BaseList<T> drow(row(rows()-1));
    drow=a;
  }

  template<typename T> void Matrix<T>::create(int rows_,int cols_,const T& v)
  {
    store_.create(rows_*cols_,v);
    r_=rows_; c_=cols_;
  }
  
  template<typename T> template<class T2> void Matrix<T>::create(int rows_,int cols_, T2* vp)
  {
    typedef typename Unconst<T2*>::Result unconst_t; 
    store_.create(rows_*cols_,const_cast<unconst_t>(vp));
    r_=rows_; c_=cols_;
  }
  
template<class T> void Matrix<T>::allrows(List<BaseList<T> >& res)
{
  res.create(r_);
  for (size_t r=r_;r--;)
    res(r).create(row(r));
}

template<class T> void Matrix<T>::swap(Matrix<T>& a)
{ 
  store_.swap(a.store_);
  ::std::swap(r_,a.r_);
  ::std::swap(c_,a.c_);
}

//template<class T> Matrix<T>& Matrix<T>::operator/= (const T& b)
//{
//  apply_ip(doesdivide_ip<T,T>(),*this,b);
//  return *this;
//}

template<class T> void Matrix<T>::transpose()
{
  if ((rows()==1) || (cols()==1)) { //special case of vector
    ::std::swap(r_,c_);
    return;
  }
  Matrix<T>& data=*this;
  if (!issquare(data)) {
    Matrix<T> tmp;
    ::libcmatrix::transpose(tmp,*this);
    swap(tmp);
    return;
  }
  T tmp;
  for (size_t r=rows()-1;r>0;r--) {
    T *vec=data.vector(r);
    for (size_t c=r;c--;) {
      tmp=vec[c];
      vec[c]=data(c,r);
      data(c,r)=tmp;
    }
  }
}

template<class T> void Matrix<T>::identity(int n)
{
  create(n,n);
  identity();
}

template<class T> void Matrix<T>::identity()
{
  *this=T(0);
  static const T one=T(1); // create object (may be complex)
  for (size_t i=rows();i--;) (*this)(i,i)=one;
}

template<class T> class SingularMatrix : public Failed {
public:
  explicit SingularMatrix(const char* errmv, const Matrix<T>& v) 
    : Failed(errmv),
      v_(v) {}
  virtual ~SingularMatrix() throw() {};
  const Matrix<T>& matrix() const { return v_; }
private:
  Matrix<T> v_;
};

#ifndef LCM_SUPPRESS_VIEWS
template<class T> const IndirectList<T,slice> Matrix<T>::diag(int which) const
{
  if (!issquare(*this)) throw NotSquare("diag");
  const size_t n=rows();
  if (::std::abs(which)>=n)
    throw BadIndex("diag");
  if (which<0)
    return IndirectList<T,slice>(row(),slice(-which*n,n+which,n+1));
  else
    return IndirectList<T,slice>(row(),slice(which,n-which,n+1));
}

template<class T> IndirectList<T,slice> Matrix<T>::diag(int which)
{
  if (!issquare(*this)) throw NotSquare("diag");
  const size_t n=rows();
  if (::std::abs(which)>=n)
    throw BadIndex("diag");
  if (which<0)
    return IndirectList<T,slice>(row(),slice(-which*n,n+which,n+1));
  else
    return IndirectList<T,slice>(row(),slice(which,n-which,n+1));
}
#endif

} //namespace libcmatrix

#ifndef LCM_SUPPRESS_VIEWS

#include "IndirectMatrix.h"

namespace libcmatrix {
  template<class T> RawMatrix<T>::RawMatrix(const Matrix<T>& a, const range& sel)
    : step_(a.cols()) 
  {
    r_=c_=sel.size();
    const size_t off=sel.start();
    data_=const_cast<T*>(a.vector(off))+off;
  }
  
  template<class T> RawMatrix<T>::RawMatrix(const Matrix<T>& a, const range& rsel, const range& csel)
    : r_(rsel.size()), c_(csel.size()), step_(a.cols())
  {
    data_=const_cast<T*>(a.vector(rsel.start()))+csel.start();
  }

#ifdef LCM_USE_EXTERNTEMPLATE
#ifdef LCM_LEVEL2_INSTANT_REAL
#define LCM_L2_EXTERN_REAL
#else
#define LCM_L2_EXTERN_REAL extern
#endif
#ifdef LCM_LEVEL2_INSTANT_COMPLEX
#define LCM_L2_EXTERN_COMPLEX
#else
#define LCM_L2_EXTERN_COMPLEX extern
#endif
LCM_L2_EXTERN_REAL template class RawMatrix<double>;
LCM_L2_EXTERN_COMPLEX template class RawMatrix<complex>;
LCM_L2_EXTERN_REAL template class Matrix<double>;
LCM_L2_EXTERN_COMPLEX template class Matrix<complex>;
#endif

}

#endif

#endif //Matrix_h_

// template<class T> inline void _transpose_blocked(Matrix<T>& d,const Matrix<T>& a)
// {
//   if (!!d && issame(d,a))
//     throw ArgumentClash("transpose");

//   size_t ar=a.rows();
//   size_t ac=a.cols();
  
//   d.create(ac,ar);
  
//   size_t r=0;
//   while (r<ar) {
//     size_t finr=r+4;
//     if (finr>ar) {
//       for (size_t c=ac;c--;) {
// 	T* dptrc=d.vector(c);
//  	for (size_t sr=r;sr<ar;sr++)
//  	  dptrc[sr]=a(sr,c);
//       }
//       return;
//     }
//     else {
//       for (size_t c=ac;c--;) {
// 	T* dptrc=d.vector(c);
// 	for (size_t sr=r;sr<finr;sr+=4) {
// 	  dptrc[sr]=a(sr,c);
// 	  dptrc[sr+1]=a(sr+1,c);
// 	  dptrc[sr+2]=a(sr+2,c);
// 	  dptrc[sr+3]=a(sr+3,c);
// 	}
//       }
//     }
//     r=finr;
//   }
// }
