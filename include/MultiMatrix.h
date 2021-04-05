#ifndef _MultiMatrix_h_
#define _MultiMatrix_h_

#include "ListList.h"
#include "Matrix.h"
#include "IndirectList.h"

//Max dimensions for MultiMatrix
#define LCM_MULTIMATRIX_DIMS 4

namespace libcmatrix {

  template<class T,size_t N> class PermutedMultiMatrix;
  template<size_t,size_t> struct Make_ {};
  template<class T,size_t spare,size_t M> struct selectN {};
  template<class T,size_t spare,size_t M> struct const_selectN {};

  struct create_ {
    template<typename T,size_t N> static inline void func(MultiMatrix<T,N>& a,const BaseList<size_t>& dims)
    { a.dimensions(dims); }

    template<typename T> static inline void func(MultiMatrix<T,0>& a,const T& v)
    { a.fill_create(v); }

    template<typename T> static inline void func(MultiMatrix<T,0>& a, T* v)
    { a.pointer_create(v); }

    template<typename T> static inline void func(MultiMatrix<T,1>& a,int n)
    { a.raw_create(n); }

    template<typename T> static inline void func(MultiMatrix<T,1>& a,int n,const T& v)
    { a.fill_create(n,v); }

    template<typename T> static inline void func(MultiMatrix<T,1>& a,int n,T* v)
    { a.pointer_create(n,v); }

    template<typename T> static inline void func(MultiMatrix<T,2>& a,int r,int s)
    { a.raw_create(r,s); }

    template<typename T> static inline void func(MultiMatrix<T,2>& a,int r,int s,const T& v)
    { a.fill_create(r,s,v); }

    template<typename T> static inline void func(MultiMatrix<T,2>& a,int r,int s, T* v)
    { a.pointer_create(r,s,v); }

    template<typename T> static inline void func(MultiMatrix<T,3>& a,int r,int s,int t) { a.raw_create(r,s,t); }

    template<typename T> static inline void func(MultiMatrix<T,3>& a,int r,int s,int t,const T& v)
    { a.fill_create(r,s,t,v); }

    template<typename T> static inline void func(MultiMatrix<T,3>& a,int r,int s,int t,T* v)
    { a.pointer_create(r,s,t,v); }

    template<typename T> static inline void func(MultiMatrix<T,4>& a,int r,int s,int t,int u)
    { a.raw_create(r,s,t,u); }

    template<typename T> static inline void func(MultiMatrix<T,4>& a,int r,int s,int t,int u,const T& v)
    { a.fill_create(r,s,t,u,v); }

    template<typename T> static inline void func(MultiMatrix<T,4>& a,int r,int s,int t,int u,T* v)
    { a.pointer_create(r,s,t,u,v); }
  };
  
#ifndef LCM_SUPPRESS_VIEWS
  template<class T,size_t N> class IndirectMultiMatrix;
  template<class T,size_t N> std::ostream& operator<< (std::ostream&, const IndirectMultiMatrix<T,N>&);

  struct Missing_;
  struct MultiSummary;
  
  template<class T1, class T2 =Missing_, class T3 =Missing_, class T4 =Missing_> struct MultiIndex {
    enum { rank = LCM_INDEX_DIM(T1) + LCM_INDEX_DIM(T2) + LCM_INDEX_DIM(T3) + LCM_INDEX_DIM(T4) };
  };
#endif

#if (BOUNDS_CHECK)
#define LCM_CheckBoundsDefault CheckBounds
#else
#define LCM_CheckBoundsDefault CheckNone
#endif

  template<size_t N> struct CheckNone {
    inline static void check(bool, const char*) {}
    inline static void check_bounds_open(const size_t [], size_t) {}
    inline static void check_bounds_open(const size_t [], size_t, size_t) {}
    inline static void check_bounds_open(const size_t [], size_t, size_t, size_t) {}
    inline static void check_bounds_open(const size_t [], size_t, size_t, size_t, size_t) {}
    inline static void check_bounds(const size_t [], size_t) {}
    inline static void check_bounds(const size_t [], size_t, size_t) {}
    inline static void check_bounds(const size_t [], size_t, size_t, size_t) {}
    inline static void check_bounds(const size_t [], size_t, size_t, size_t, size_t) {} 
  };

  template<size_t N> struct CheckBounds {
    inline static void check(bool cond, const char* errm =NULL) {
      if (!cond)
	throw BadIndex(errm);
    }

    inline static void check_bounds_open(const size_t dim[], size_t r) {
      LCM_STATIC_CHECK( N>=1, Dimensions_dont_match );
      if (r>=dim[0])
	throw BadIndex("check_bounds",r,dim[0]);
    }

    inline static void check_bounds_open(const size_t dim[],size_t r,size_t s) {
    LCM_STATIC_CHECK( N>=2, Dimensions_dont_match );
    if (r>=dim[0] || s>=dim[1])
      throw BadIndex("check_bounds",r,dim[0],s,dim[1]);
    }

    inline static void check_bounds_open(const size_t dim[], size_t r,size_t s,size_t t) {
      LCM_STATIC_CHECK( N>=3, Dimensions_dont_match );
      if (r>=dim[0] || s>=dim[1] || t>=dim[2])
	throw BadIndex("check_bounds",r,dim[0],s,dim[1],t,dim[2]);
    }

    inline static void check_bounds_open(const size_t dim[], size_t r,size_t s,size_t t,size_t u) {
      LCM_STATIC_CHECK( N>=4, Dimensions_dont_match );
      if (r>=dim[0] || s>=dim[1] || t>=dim[2] || u>=dim[3])
	throw BadIndex("check_bounds");
    }

    inline static void check_bounds(const size_t dim[],size_t r) {
      LCM_STATIC_CHECK( N==1, Dimensions_dont_match );
      if (r>=dim[0])
	throw BadIndex("check_bounds",r,dim[0]);
    }

    inline static void check_bounds(const size_t dim[],size_t r,size_t s) {
      LCM_STATIC_CHECK( N==2, Dimensions_dont_match );
      if (r>=dim[0] || s>=dim[1])
	throw BadIndex("check_bounds",r,dim[0],s,dim[1]);
    }

    inline static void check_bounds(const size_t dim[],size_t r,size_t s,size_t t) {
      LCM_STATIC_CHECK( N==3, Dimensions_dont_match );
      if (r>=dim[0] || s>=dim[1] || t>=dim[2])
	throw BadIndex("check_bounds",r,dim[0],s,dim[1],t,dim[2]);
    }

    inline static void check_bounds(const size_t dim[],size_t r,size_t s,size_t t,size_t u) {
      LCM_STATIC_CHECK( N==4, Dimensions_dont_match );
      if (r>=dim[0] || s>=dim[1] || t>=dim[2] || u>=dim[3])
	throw BadIndex("check_bounds");
    }

  };

  template<size_t N, template<size_t> class CheckClass =LCM_CheckBoundsDefault>
  class Indexer : public CheckClass<N> {
  public:
    size_t multiplier(size_t n) const {
      this->check(n<N,"Indexer::multiplier");
      return mults[n]; }

   template<int M> size_t multiplier(Int2Type<M>) const {
      LCM_STATIC_CHECK( M<N, Indexer_multipler_bad_index );
      return mults[M];
    }

    size_t dimension(size_t n) const {
      this->check(n<N,"Indexer::dimension");
      return dim[n]; }

    template<int M> size_t dimension(Int2Type<M>) const {
      LCM_STATIC_CHECK( M<N, Indexer_dimension_bad_index );
      return dim[M];
    }
    
    size_t size() const 
    { return mults[0]*dim[0]; }

    const BaseList<size_t> dimensions() const {
      return BaseList<size_t>(N,const_cast<size_t*>(dim));
    }
    
    void resetzero() { 
      for (size_t i=N;i--;)
	dim[i]=0;
      mults[N-1]=1;
    } // Contents not set
    
    template<class T2> inline void copydims(const T2& a)
    {
      LCM_STATIC_CHECK( N==LCM_DIM(T2), copy_Dimensions_dont_match );
      size_t sofar=1;
      for (size_t i=N;i--;) {
	mults[i]=sofar;
	sofar*=(dim[i]=a.dimension(i));
      }      
    }

  Indexer() { resetzero(); }
  Indexer(size_t r) { (void)setdims(r); }
  Indexer(size_t r, size_t s) { (void)setdims(r,s); }
  Indexer(size_t r, size_t s, size_t t) { (void)setdims(r,s,t); }
  Indexer(size_t r, size_t s, size_t t, size_t u) { (void)setdims(r,s,t,u); }

    Indexer(const BaseList<size_t>& inds) { (void)setdims(inds); }

    size_t setdims(const BaseList<size_t>& a) {
      if (a.size()>N)
	throw InvalidParameter("Indexer(): too many dimensions");
      size_t sofar=1;
      const size_t offset=N-a.size();
      for (size_t i=a.size();i--;) {
	mults[i+offset]=sofar;
	sofar*=(dim[i+offset]=a(i));
      }      
      for (size_t i=offset;i--;) {
	mults[i]=sofar;
	dim[i]=1;
      }
      return sofar;
    }
	 
    size_t setdims(size_t r) {
      LCM_STATIC_CHECK( N==1, Dimensions_dont_match );
      dim[0]=r;
      mults[0]=1;
      return r;
    }
    size_t setdims(size_t r,size_t s) {
      LCM_STATIC_CHECK( N==2, Dimensions_dont_match );
      if (r<1)
	throw InvalidParameter("Indexer: dimensions must be >0");
      dim[0]=r; mults[0]=dim[1]=s; mults[1]=1;
      return r*s;
    }
    size_t setdims(size_t r,size_t s,size_t t) {
      LCM_STATIC_CHECK( N==3 , Dimensions_dont_match );
      if (r<=0 || s<=0)
	throw InvalidParameter("Indexer: dimensions must be >0");
      dim[0]=r; dim[1]=s;
      mults[1]=dim[2]=t; mults[0]=s*t; mults[2]=1;
      return r*mults[0];
    }
    
    size_t setdims(size_t r,size_t s,size_t t,size_t u) {
      LCM_STATIC_CHECK( N==4, Dimensions_dont_match );
      if (r<=0 || s<=0 || t<=0)
	throw InvalidParameter("MultiMatrix<T,N>: dimensions must be >0");
      dim[0]=r; dim[1]=s; dim[2]=t;
      mults[2]=dim[3]=u; mults[0]=s*(mults[1]=t*u); mults[3]=1;
      return r*mults[0]; 
    }
    
    template< template<size_t> class CheckBounds_> void swap(Indexer<N,CheckBounds_>& a) {
      for (size_t i=N;i--;) {
	::std::swap(mults[i],a.mults[i]);
	::std::swap(dim[i],a.dim[i]);
      }
    }

    size_t operator()(size_t r) const
    { this->check_bounds_open(dim,r);
      return r*mults[0]; }

    size_t operator()(size_t r,size_t s) const
    { this->check_bounds_open(dim,r,s);
      return r*mults[0]+s*mults[1]; }

    size_t operator()(size_t r,size_t s,size_t t) const
    { this->check_bounds_open(dim,r,s,t); 
      return r*mults[0]+s*mults[1]+t*mults[2]; }

    size_t operator()(size_t r,size_t s,size_t t,size_t u) const
    { this->check_bounds_open(dim,r,s,t,u);
      return r*mults[0]+s*mults[1]+t*mults[2]+u*mults[3]; }
    
    size_t index(size_t r) const
    { this->check_bounds(dim,r);
      return r; }

    size_t index(size_t r,size_t s) const
    { this->check_bounds(dim,r,s);
      return r*mults[0]+s; }

    size_t index(size_t r,size_t s,size_t t) const
    { this->check_bounds(dim,r,s,t);
      return r*mults[0]+s*mults[1]+t; }

    size_t index(size_t r,size_t s,size_t t,size_t u) const
    { this->check_bounds(dim,r,s,t,u);
      return r*mults[0]+s*mults[1]+t*mults[2]+u; }
    
    void reverse(size_t& r, size_t ind) const {
      LCM_STATIC_CHECK(N==1, Indexer_non1D_object);
      this->check(ind<dim[0],"Indexer::reverse");
      r=ind;
    }

    void reverse(size_t& r, size_t& s, size_t ind) const {
      LCM_STATIC_CHECK(N==2, Indexer_non2D_object);
      this->check(ind<size(),"Indexer::reverse");
      r=ind / mults[0];
      s=ind-r*mults[0];
    }

    void reverse(size_t& r, size_t& s, size_t& t, size_t ind) const {
      LCM_STATIC_CHECK(N==3, Indexer_non3D_object);
      this->check(ind<size(),"Indexer::reverse");
      r=ind / mults[0];
      ind-=r*mults[0];
      s=ind / mults[1];
      t=ind-s*mults[1];
    }

    void reverse(size_t& r, size_t& s, size_t& t, size_t& u, size_t ind) const {
      LCM_STATIC_CHECK(N==4, Indexer_non4D_object);
      this->check(ind<size(),"Indexer::reverse");
      r=ind / mults[0];
      ind-=r*mults[0];
      s=ind / mults[1];
      ind-=s*mults[1];
      t=ind / mults[2];
      u=ind-t*mults[2];
    }

    class permuted_iterator;
        
  permuted_iterator permuted_begin(const BaseList<size_t>& order =BaseList<size_t>()) const
  { return permuted_iterator(*this,order); }

  permuted_iterator permuted_end(const BaseList<size_t>& =BaseList<size_t>()) const
  { return permuted_iterator(*this); }

  bool operator==(const Indexer<N>& a) const {
    for (size_t i=N;i--;) {
      if (dim[i]!=a.dim[i])
	return false;
    }
    return true;
  }

  bool operator!=(const Indexer<N>& a) const {
    return !(*this==a);
  }

    void print(std::ostream& ostr) const { ostr << "Dimensions: " << BaseList<size_t>(N,(size_t*)dim) << "  Multipliers: " << BaseList<size_t>(N,(size_t*)mults); }

  private:
  size_t dim[N ? N : 1];
  size_t mults[N ? N : 1];
  };

  template<size_t N> std::ostream& operator<< (std::ostream& ostr, const Indexer<N>& a) {
    a.print(ostr);
    return ostr; }

  template<class T,size_t N> class PermutedIterator;
  template<class T,size_t N> class MultiMatrix;
  template<class T,size_t N> std::ostream& operator<< (std::ostream&, const MultiMatrix<T,N>&);
  
template<class T,size_t N> class MultiMatrix {
public:
  MultiMatrix(mxflag::tempflag tflag =mxflag::normal)
    : store_(tflag) {
    indexer_.resetzero();
  }
  template<typename T2> MultiMatrix(const T2& a,mxflag::tempflag tflag =mxflag::normal)
    : store_(tflag) { 
    Make_<N,LCM_DIM(T2)>::create(*this,a);
  }
  MultiMatrix(const MultiMatrix<T,N>& a,mxflag::tempflag tflag =mxflag::normal) 
    : store_(a.row(),tflag) {
    copydims(a);
  }
  template<class T2> MultiMatrix(int rs,const T2& v,mxflag::tempflag tflag =mxflag::normal)
    : store_(tflag) {
    create(rs,v);
  }
  template<class T2> MultiMatrix(int rs,int cs,const T2& v,mxflag::tempflag tflag =mxflag::normal)
    : store_(tflag) {
    create(rs,cs,v);
  }
  template<class T2> MultiMatrix(int rs,int cs,int ts,const T2& v,mxflag::tempflag tflag =mxflag::normal)
    : store_(tflag) {
    create(rs,cs,ts,v);
  }
  template<class T2> MultiMatrix(int rs,int cs,int ts,int us,const T2& v,mxflag::tempflag tflag =mxflag::normal)
    : store_(tflag) {
    create(rs,cs,ts,us,v);
  }
  template<typename T2> MultiMatrix<T,N>& operator= (const T2& a) {
    Make_<N,LCM_DIM(T2)>::assign(*this,a); return *this;
  }

  MultiMatrix& operator= (const MultiMatrix<T,N>& a) {
    if (this!=&a) {
      store_=a.row();
      copydims(a);
    }
    return *this;
  }

  typedef DynamicList<T> StorageType;

  friend class PermutedMultiMatrix<T,N>;
  template<class,size_t,size_t> friend struct selectN;
  template<class,size_t,size_t> friend struct const_selectN;

  // definitions for Standard C++ Library containers
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;  
  typedef int difference_type;
  typedef size_t size_type;

  size_t size() const { return store_.size(); }
  void clear() { store_.clear(); indexer_.resetzero(); }
  bool empty() const { return store_.empty(); }
  bool operator!() const { return (store_.size()==0); }
  //  void ensureempty() { store_.ensureempty(); }

  size_t rows() const { 
    LCM_STATIC_CHECK( N>1, rows_applied_to_less_than_1d_object );
    return size()/indexer_.dimension(Int2Type<N-1>());
  }

  size_t cols() const { 
    LCM_STATIC_CHECK( N==2, cols_applied_to_non2D_object );
    return indexer_.dimension(Int2Type<1>());
  }

  size_t dimension(size_t n) const {
    return indexer_.dimension(n);
  }

  const BaseList<size_t> dimensions() const {
    return indexer_.dimensions();
  }

  template<class T2> void set_dimensions(const T2& a) {
    indexer_.copydims(a);
    if (N)
      store_.create(a.size());
  }

  void dimensions(const BaseList<size_t>& dims) { 
    store_.create(indexer_.setdims(dims));
  }
   
  template<class T2> inline void create(const T2& v) { create_::func(*this,v); }
  template<class T2> inline void create(int r,const T2& v) { create_::func(*this,r,v); }
  template<class T2> inline void create(int r,int s,const T2& v) { create_::func(*this,r,s,v); }
  template<class T2> inline void create(int r,int s,int t,const T2& v) { create_::func(*this,r,s,t,v); }
  template<class T2> inline void create(int r,int s,int t,int u,const T2& v) { create_::func(*this,r,s,t,u,v); }

  void raw_create(int n) { store_.create(indexer_.setdims(n)); }
  void fill_create(const T& v) { store_.create(1,v); }
  void fill_create(int n,const T& v) { store_.create(indexer_.setdims(n),v); }
  void raw_create(int r,int s) { store_.create(indexer_.setdims(r,s)); }
  void fill_create(int r,int s,const T& v) { store_.create(indexer_.setdims(r,s),v); }
  void raw_create(int r,int s,int t) { store_.create(indexer_.setdims(r,s,t)); }
  void fill_create(int r,int s,int t,const T& v) { store_.create(indexer_.setdims(r,s,t),v); }
  void raw_create(int r,int s,int t,int u) { store_.create(indexer_.setdims(r,s,t,u)); }
  void fill_create(int r,int s,int t,int u,const T& v) { store_.create(indexer_.setdims(r,s,t,u),v); }

  void pointer_create(T* v) { store_.create(1,v); }
  void pointer_create(int n, T* v) { store_.create(indexer_.setdims(n),v); }
  void pointer_create(int r,int s, T* v) { store_.create(indexer_.setdims(r,s),v); }
  void pointer_create(int r,int s,int t, T* v) { store_.create(indexer_.setdims(r,s,t),v); }
  void pointer_create(int r,int s,int t,int u, T* v) { store_.create(indexer_.setdims(r,s,t,u),v); }

  typedef typename StorageType::iterator iterator;
  typedef typename StorageType::const_iterator const_iterator;

  typedef PermutedIterator<T,N> permuted_iterator;
  typedef PermutedIterator<const T,N> const_permuted_iterator;

  iterator begin() { return store_.begin(); }
  const_iterator begin() const { return store_.begin(); }
  iterator end() { return store_.end(); }
  const_iterator end() const { return store_.end(); }

  permuted_iterator permuted_begin(const BaseList<size_t>& order =BaseList<size_t>())
  { return PermutedIterator<T,N>(store_.vector(),indexer_,order); }

  permuted_iterator permuted_end(const BaseList<size_t>& =BaseList<size_t>())
  { return PermutedIterator<T,N>(indexer_); }

  const_permuted_iterator permuted_begin(const BaseList<size_t>& order =BaseList<size_t>()) const
  { return PermutedIterator<const T,N>(const_cast<const T*>(store_.vector()),indexer_,order); }

  const_permuted_iterator permuted_end(const BaseList<size_t>& =BaseList<size_t>()) const
  { return PermutedIterator<const T,N>(indexer_); }
  
  bool istemporary() const { return store_.istemporary(); }
  bool isdynamic() const { return store_.isdynamic(); }
  friend std::ostream& operator<< <>(std::ostream&,const MultiMatrix<T,N>&);

  void swap(MultiMatrix<T,N>& a) {
    store_.swap(a.store_); //swap contents;
    indexer_.swap(a.indexer_);
  }

  T* vector() const { return store_.vector(); }

  T& select() {
    return *(store_.vector());
  }
  T& select(size_t r) {
    return store_(r);
  }
  T& select(size_t r,size_t s) {
    return store_(indexer_(r,s));
  }
  T& select(size_t r,size_t s,size_t t) {
    return store_(indexer_(r,s,t));
  }
  T& select(size_t r,size_t s,size_t t,size_t u) {
    return store_(indexer_(r,s,t,u));
  }
  const T& select() const {
    return *(store_.vector());
  }

  const T& select(size_t r) const {
    return store_(r);
  }

  const T& select(size_t r,size_t s) const {
    return store_(indexer_(r,s));
  }

  const T& select(size_t r,size_t s,size_t t) const {
    return store_(indexer_(r,s,t));
  }

  const T& select(size_t r,size_t s,size_t t,size_t u) const {
    return store_(indexer_(r,s,t,u));
  }

#ifndef LCM_SUPPRESS_VIEWS
  template<class Map> typename selectN< T,N-1,MultiIndex<Map>::rank >::return_type
    operator()(const Map& map) {
    return selectN< T,N-1,MultiIndex<Map>::rank >::func(*this,map);
  }
  
  template<class Map> typename const_selectN< T,N-1,MultiIndex<Map>::rank >::return_type
    operator()(const Map& map) const {
    return const_selectN< T,N-1,MultiIndex<Map>::rank >::func(*this,map);
  }
  
  template<class Map1,class Map2> typename selectN< T,N-2,MultiIndex<Map1,Map2>::rank >::return_type
    operator()(const Map1& map1,const Map2& map2) {
    return selectN< T,N-2,MultiIndex<Map1,Map2>::rank >::func(*this,map1,map2);
  }

  template<class Map1,class Map2> typename const_selectN< T,N-2,MultiIndex<Map1,Map2>::rank >::return_type
    operator()(const Map1& map1,const Map2& map2) const {
    return const_selectN< T,N-2,MultiIndex<Map1,Map2>::rank >::func(*this,map1,map2);
  }

  template<class Map1,class Map2,class Map3> typename selectN< T,N-3,MultiIndex<Map1,Map2,Map3>::rank >::return_type
    operator()(const Map1& map1,const Map2& map2,const Map3& map3) {
    return selectN<T,N-3,MultiIndex<Map1,Map2,Map3>::rank >::func(*this,map1,map2,map3);
  }

  template<class Map1,class Map2,class Map3> typename const_selectN< T,N-3,MultiIndex<Map1,Map2,Map3>::rank >::return_type
    operator()(const Map1& map1,const Map2& map2,const Map3& map3) const {
    return const_selectN<T,N-3,MultiIndex<Map1,Map2,Map3>::rank >::func(*this,map1,map2,map3);
  }
#else
  typename selectN<T,N-1,0>::return_type operator()(size_t r) {
    return selectN<T,N-1,0>::func(*this,r);
  }
  
  typename const_selectN<T,N-1,0>::return_type operator()(size_t r) const {
    return const_selectN<T,N-1,0>::func(*this,r);
  }
  
  typename selectN<T,N-2,0>::return_type operator()(size_t r, size_t s) {
    return selectN<T,N-2,0>::func(*this,r,s);
  }

  typename const_selectN<T,N-2,0>::return_type operator()(size_t r, size_t s) const {
    return const_selectN<T,N-2,0>::func(*this,r,s);
  }

  typename selectN<T,N-3,0>::return_type operator()(size_t r, size_t s, size_t t) {
    return selectN<T,N-3,0>::func(*this,r,s,t);
  }

  typename const_selectN<T,N-3,0>::return_type operator()(size_t r, size_t s, size_t t) const {
    return const_selectN<T,N-3,0>::func(*this,r,s,t);
  }
#endif
  //LCM_SUPPRESS_VIEWS

  T& operator() (size_t r,size_t s,size_t t,size_t u) {
    return store_(indexer_(r,s,t,u));
  }
  
  const T& operator() (size_t r,size_t s,size_t t,size_t u) const {
    return store_(indexer_(r,s,t,u));
  }

  typename selectN<T,N,0>::return_type operator()() {
    return selectN<T,N,0>::func(*this);
  }

  typename const_selectN<T,N,0>::return_type operator()() const {
    return const_selectN<T,N,0>::func(*this);
  }

  typename selectN<T,N-1,0>::return_type back() {
    return selectN<T,N-1,0>::back(*this);
  }

  typename selectN<T,N-1,0>::return_type front() {
    return selectN<T,N-1,0>::front(*this);
  }

  typename const_selectN<T,N-1,0>::return_type back() const {
    return const_selectN<T,N-1,0>::back(*this);
  }

  typename const_selectN<T,N-1,0>::return_type front() const {
    return const_selectN<T,N-1,0>::front(*this);
  }

  BaseList<T>& row() {
    return store_;
  }

  const BaseList<T>& row() const {
    return store_;
  }

  const BaseList<T> row(size_t r) const {
    const size_t lastdim=indexer_.dimension(Int2Type<N-1>());
    const size_t offset=r*lastdim;
    if (offset>=store_.size())
      throw BadIndex("MultiMatrix::row",r,store_.size()/lastdim);
    return BaseList<T>(lastdim,store_.vector()+offset);
  }

  BaseList<T> row(size_t r) { //const fudge to avoid code duplication
    return const_cast<const MultiMatrix<T,N>* >(this)->row(r);
  }

  template<class T2> MultiMatrix<T,N>& operator+= (const T2& b) { return add_ip(*this,b); }
  template<class T2> MultiMatrix<T,N>& operator-= (const T2& b) { return subtract_ip(*this,b); }
  template<class T2> MultiMatrix<T,N>& operator*= (const T2& b) { return multiply_ip(*this,b); }
  template<class T2> MultiMatrix<T,N>& operator/= (const T2& b) { return divide_ip(*this,b); }
   
 private:
  StorageType store_;
  Indexer<N,LCM_CheckBoundsDefault> indexer_;

  template<class T2> inline void copydims(const T2& a)
    {
    if (issame(vector(),a.vector()))
      throw ArgumentClash("MultiMatrix::copydims"); //should never be constructing 
    indexer_.copydims(a);
    }
};


  template<class T,size_t N,class T2> inline bool operator== (const MultiMatrix<T,N>& a, const T2& b) {
    return areequal(a,b);
  }
  template<class T,size_t N,class T2> inline bool operator!= (const MultiMatrix<T,N>& a, const T2& b) {
    return arenotequal(a,b);
  }

template<class T> void spy(std::ostream& ostr,const MultiMatrix<T,3>& a,double tol =1e-10)
{
  if (!a) {
    ostr << "<undefined>" << std::endl;
    return;
  }
  if (a.dimension(3)<a.dimension(2))
    spyrow(ostr,a,tol);
  else {
    for (size_t i=0;i<a.dimension(0);i++)
      ostr << "matrix " << i << '\n' << a(i) << std::endl;
  }
}

template<class T> void spyrow(std::ostream& ostr,const MultiMatrix<T,3>& a,double tol)
{
  const size_t rows=a.dimension(1);
  for (size_t r=0;r<rows;r++) {
    for (size_t c=0;c<a.dimension(0);c++) {
      spy(ostr,a(c,r),tol);
      ostr << ' ';
    }
    if (r==rows-1)
      ostr << std::endl;
    else
      ostr << '\n';
  }
}

template<class T> void spy(std::ostream& ostr,const MultiMatrix<T,4>& a,double tol =1e-10)
{
  if (!a) {
    ostr << "<undefined>" << std::endl;
    return;
  }
  for (size_t r=0;r<a.dimension(0);r++) {
    spyrow(ostr,a(r),tol);
    ostr << std::endl;
  }
}

// Implementation details below

  template<size_t N> struct Make_<N,0> { 
    template<typename T1,typename T2> static inline void assign(MultiMatrix<T1,N>& a,const T2& v) { a.row()=v; }
  };
  template<size_t N> struct Make_<N,N> {
    template<typename T1,typename T2> static inline void assign(MultiMatrix<T1,N>& a,const T2& v) { ::libcmatrix::assign(a,v); }
    template<typename T1,typename T2> static inline void create(MultiMatrix<T1,N>& a,const T2& v) {
      if (a.isdynamic())
	a.store_.create(v.size(),v.begin());
      else
	throw Failed("Can't construct non-dynamic MultiMatrix from this object");
      a.set_dimensions(v);
    }
    template<typename T> static inline void create(MultiMatrix<T,N>& a,const Matrix<T>& v) {
      a.pointer_create(v.rows(),v.cols(),v.vector());
    }
    template<typename T> static inline void create(MultiMatrix<T,N>& a,const BaseList<T>& v) {
      a.pointer_create(v.size(),v.vector());
    }
  };
  //special case for copying List<Matrix> into MultiMatrix<3>
  template<> struct Make_<3,1> {
    template<typename T1> static inline void assign(MultiMatrix<T1,3>& a,const BaseList< Matrix<T1> >& v) {
      size_t n=v.size();
      const Matrix<T1>& v0(v(0));
      if (!v0)
	throw Undefined("MultiMatrix<T> create from list of Matrix<T>");
      a.raw_create(n,v0.rows(),v0.cols());
      for (;n--;)
	a(n)=v(n);
    }
  };
  template<> struct Make_<1,0> {
    template<typename T1> static inline void create(MultiMatrix<T1,1>& a,int n) { a.raw_create(n); }
  };
  template<> struct Make_<0,0> { 
    template<typename T1,typename T2> static inline void assign(MultiMatrix<T1,0>& a,const T2& v) { a.select()=v; }
    template<typename T1,typename T2> static inline void create(MultiMatrix<T1,0>& a,const T2& v) { a.fill_create(v); }
  };

struct Missing_ {
    typedef size_t iterator;
    typedef size_t const_iterator;
    static size_t begin() { return 0; }
    static size_t end() { return 0; }
    static size_t size() { return 0; }
};

 template<> struct type_traits<Missing_> {
   static const bool dimensionality=0;
   typedef void* value_type; //needed to keep LCM_INDEX_DIM happy
 };

 template<class T> struct selectN<T,0,0> { 
   typedef T& return_type;

   static T& func(MultiMatrix<T,0>& a) {
     return a.select(); 
   }
   static T& func(MultiMatrix<T,1>& a,size_t r) {
     return a.select(r); 
   }
   static T& func(MultiMatrix<T,2>& a,size_t r,size_t s) {
     return a.select(r,s);
   }
   static T& func(MultiMatrix<T,3>& a,size_t r,size_t s,size_t t) {
     return a.select(r,s,t);
   }
   static T& func(MultiMatrix<T,4>& a,size_t r,size_t s,size_t t,size_t u) {
     return a.select(r,s,t,u);
   }
   static T& back(MultiMatrix<T,1>& a) {
     return a.select(a.dimension(0)-1);
   }
   static T& front(MultiMatrix<T,1>& a) {
     return a.row().front();
   }
 };

template<class T> struct selectN<T,1,0> { 
   typedef BaseList<T> return_type;

  static BaseList<T> func(MultiMatrix<T,1>& a) {
    return BaseList<T>(a.dimension(0),a.vector());
  }
  static BaseList<T> func(MultiMatrix<T,2>& a,size_t r) {
    return BaseList<T>(a.dimension(1),a.vector()+a.indexer_(r));
  }
  static BaseList<T> func(MultiMatrix<T,3>& a,size_t r,size_t s) {
    return BaseList<T>(a.dimension(2),a.vector()+a.indexer_(r,s));
  }
  static BaseList<T> func(MultiMatrix<T,4>& a,size_t r,size_t s,size_t t) {
    return BaseList<T>(a.dimension(3),a.vector()+a.indexer_(r,s,t));
  }
  static BaseList<T> back(MultiMatrix<T,2>& a) {
    return BaseList<T>(a.dimension(1),a.vector()+a.size()-a.dimension(1));
  }
  static BaseList<T>  front(MultiMatrix<T,2>& a) {
    return BaseList<T>(a.dimension(1),a.vector());
  }
 };

 template<class T> struct selectN<T,2,0> { 
   typedef Matrix<T> return_type;

   static Matrix<T> func(MultiMatrix<T,2>& a) {
     return Matrix<T>(a.dimension(0),a.dimension(1),a.vector(),mxflag::nondynamic);
   }
   static Matrix<T> func(MultiMatrix<T,3>& a,size_t r) {
     return Matrix<T>(a.dimension(1),a.dimension(2),a.vector()+a.indexer_(r),mxflag::nondynamic);
   }
   static Matrix<T> func(MultiMatrix<T,4>& a,size_t r,size_t s) {
     return Matrix<T>(a.dimension(2),a.dimension(3),a.vector()+a.indexer_(r,s),mxflag::nondynamic);
   }
   static Matrix<T> back(MultiMatrix<T,3>& a) {
     return Matrix<T>(a.dimension(1),a.dimension(2),a.vector()+a.size()-a.multiplier(0),mxflag::nondynamic);
   }
   static Matrix<T> front(MultiMatrix<T,3>& a) {
     return Matrix<T>(a.dimension(1),a.dimension(2),a.vector(),mxflag::nondynamic);
   }
 };

 template<class T> struct selectN<T,3,0> { 
   typedef MultiMatrix<T,3> return_type;

   static MultiMatrix<T,3> func(const MultiMatrix<T,3>& a) {
     return MultiMatrix<T,3>(a,mxflag::nondynamic);
   }
   static MultiMatrix<T,3> func(const MultiMatrix<T,4>& a,size_t r) {
     return MultiMatrix<T,3>(a.dimension(1),a.dimension(2),a.dimension(3),a.vector()+a.indexer_(r),mxflag::nondynamic);
   }
 };

 template<class T> struct selectN<T,4,0> { 
   typedef MultiMatrix<T,4> return_type;
   static MultiMatrix<T,4> func(const MultiMatrix<T,4>& a) {
     return MultiMatrix<T,4>(a,mxflag::nondynamic);
   }
 };

 template<class T> struct const_selectN<T,0,0> { 
   typedef const T& return_type;
   
   static const T& func(const MultiMatrix<T,1>& a,size_t r) {
     return a.select(r);
   }
   static const T& func(const MultiMatrix<T,2>& a,size_t r,size_t s) {
     return a.select(r,s);
   }
   static const T& func(const MultiMatrix<T,3>& a,size_t r,size_t s,size_t t) {
     return a.select(r,s,t);
   }
   static const T& func(const MultiMatrix<T,4>& a,size_t r,size_t s,size_t t,size_t u) {
     return a.select(r,s,t,u);
   }
   static const T& back(const MultiMatrix<T,1>& a) {
     return a.select(a.dimension(0)-1);
   }
   static const T& front(const MultiMatrix<T,1>& a) {
     return a.row().front();
   }
 };
  
 template<class T> struct const_selectN<T,1,0> { 
   typedef const BaseList<T> return_type;

   static const BaseList<T> func(const MultiMatrix<T,1>& a) {
     return BaseList<T>(a.dimension(0),a.vector());
   }
   static const BaseList<T> func(const MultiMatrix<T,2>& a,size_t r) {
     return BaseList<T>(a.dimension(1),a.vector()+a.indexer_(r));
   }
   static const BaseList<T> func(const MultiMatrix<T,3>& a,size_t r,size_t s) {
     return BaseList<T>(a.dimension(2),a.vector()+a.indexer_(r,s));
   }
   static const BaseList<T> func(const MultiMatrix<T,4>& a,size_t r,size_t s,size_t t) {
     return BaseList<T>(a.dimension(3),a.vector()+a.indexer_(r,s,t));
   }
   static const BaseList<T> back(const MultiMatrix<T,2>& a) {
     return BaseList<T>(a.dimension(1),a.vector()+a.size()-a.dimension(1));
   }
   static const BaseList<T> front(const MultiMatrix<T,2>& a) {
     return BaseList<T>(a.dimension(1),a.vector());
   }
 };

 template<class T> struct const_selectN<T,2,0> { 
   typedef const Matrix<T> return_type;

   static const Matrix<T> func(const MultiMatrix<T,2>& a) {
     return Matrix<T>(a.dimension(0),a.dimension(1),a.vector(),mxflag::nondynamic);
   }
   static const Matrix<T> func(const MultiMatrix<T,3>& a,size_t r) {
     return Matrix<T>(a.dimension(1),a.dimension(2),a.vector()+a.indexer_(r),mxflag::nondynamic);
   }
   static const Matrix<T> func(const MultiMatrix<T,4>& a,size_t r,size_t s) {
     return Matrix<T>(a.dimension(2),a.dimension(3),a.vector()+a.indexer_(r,s),mxflag::nondynamic);
   }
   static const Matrix<T> back(const MultiMatrix<T,3>& a) {
     return Matrix<T>(a.dimension(1),a.dimension(2),a.vector()+a.size()-a.indexer_.multiplier(0),mxflag::nondynamic);
   }
   static const Matrix<T> front(const MultiMatrix<T,3>& a) {
     return Matrix<T>(a.dimension(1),a.dimension(2),a.vector(),mxflag::nondynamic);
   }
 };

 template<class T> struct const_selectN<T,3,0> { 
   typedef const MultiMatrix<T,3> return_type;

   static const MultiMatrix<T,3> func(const MultiMatrix<T,3>& a) {
     return MultiMatrix<T,3>(a,mxflag::nondynamic);
   }
   static const MultiMatrix<T,3> func(const MultiMatrix<T,4>& a,size_t r) {
     return MultiMatrix<T,3>(a.dimension(1),a.dimension(2),a.dimension(3),a.vector()+a.indexer_(r),mxflag::nondynamic);
   }
   static const Matrix<T> back(const MultiMatrix<T,4>& a) {
     return MultiMatrix<T,3>(a.dimension(1),a.dimension(2),a.dimension(3),a.vector()+a.size()-a.multiplier(0),mxflag::nondynamic);
   }
   static const Matrix<T> front(const MultiMatrix<T,4>& a) {
     return MultiMatrix<T,3>(a.dimension(1),a.dimension(2),a.dimension(3),a.vector(),mxflag::nondynamic);
   }
 };

 template<class T> struct const_selectN<T,4,0> { 
   typedef const MultiMatrix<T,4> return_type;
   static const MultiMatrix<T,4> func(const MultiMatrix<T,4>& a) {
     return MultiMatrix<T,4>(a,mxflag::nondynamic);
   }
 };

template<class T,size_t N> std::ostream& operator<< (std::ostream& ostr,const MultiMatrix<T,N>& a)
{
  if (a.istemporary()) ostr << "[T]" << (N>1 ? '\n' : ' ');
  return print(a,ostr);
}

 template<typename T,size_t N> struct type_traits< MultiMatrix<T,N> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=N;
   static const size_t rank=0;
   typedef T value_type;
 };

template<class T,size_t N> struct new_trait_<T,0,N> { typedef MultiMatrix<T,N> value_type; };
template<class T,size_t N> struct new_trait_<T,N,0> { typedef MultiMatrix<T,N> value_type; };
template<class T,size_t N> struct new_trait_<T,N,N> { typedef MultiMatrix<T,N> value_type; };

 template<class T1,size_t N,class T2> inline typename new_trait<MultiMatrix<T1,N>,T2,true>::value_type
   operator+ (const MultiMatrix<T1,N>& a,const T2& b) 
   {
     typename new_trait<MultiMatrix<T1,N>,T2,true>::value_type d(mxflag::temporary);
     add(d,a,b);
     return d;
   }

#ifndef LCM_SUPPRESS_VIEWS

 struct MultiSummary {
  ScratchList< List<size_t> > lists_;
  size_t offset;
  size_t userank;
  
  MultiSummary(size_t N) : lists_(N), offset(0), userank(0) {}

  template<class Map> void addmap(size_t dimmult, const Map& map)
    {  multiply(lists_(userank++),dimmult,map); }

   void addmap(size_t dimmult, const slice& map)
   {  lists_(userank++)=dimmult*map; }

   void addmap(size_t dimmult, const range& map)
   {  addmap(dimmult,slice(map)); }

   template<class Map> void add(MultiSummary& obj,size_t dimsize, size_t dimmult,const Map& map, Int2Type<1>)
   {
#if (BOUNDS_CHECK)
     for (size_t i=map.size();i--;) {
       if (map(i)>=dimsize)
	 throw BadIndex("MultiMatrix()",map(i),dimsize);
     }
#endif
     obj.addmap(dimmult,map);
   }

   void add(MultiSummary& obj,size_t dimsize, size_t dimmult,const range& map, Int2Type<1>)
   {
     if (map.size()) {
#if (BOUNDS_CHECK)
       if (map.max()>=dimsize)
	 throw BadIndex("range",map.max(),dimsize);
#endif
       obj.addmap(dimmult,map);
     }
     else
       obj.addmap(dimmult,range(0,dimsize-1));
   }

#if (BOUNDS_CHECK)
   void add(MultiSummary& obj, size_t dimsize, size_t dimmult, size_t i, Int2Type<0>) {
     if (i>=dimsize)
       throw BadIndex("MultiMatrix()",i,dimsize);
     obj.addindex(dimmult,i);
   }
#else
   void add(MultiSummary& obj, size_t, size_t dimmult, size_t i, Int2Type<0>) {
     obj.addindex(dimmult,i);
   }
#endif

   void addindex(size_t dimmult, size_t i) { offset+=dimmult*i; }
};
  
 template<class T,size_t M> struct selectN<T,0,M> {
    typedef IndirectMultiMatrix<T,M> return_type;

    template<class Map> static return_type func(MultiMatrix<T,1>& a,const Map& map) {
      typedef IndirectMultiMatrix< T, MultiIndex<Map>::rank > return_type;
      MultiSummary summary(MultiIndex<Map>::rank);
      summary.add(a.dimension(0), a.multiplier(0),map,Int2Type<LCM_INDEX_DIM(Map)>() );
      return return_type(a,summary);
    }
    template<class Map1,class Map2> static return_type func(MultiMatrix<T,2>& a,const Map1& map1,const Map2& map2) {
      typedef IndirectMultiMatrix< T, MultiIndex<Map1,Map2>::rank > return_type;
      MultiSummary summary(MultiIndex<Map1,Map2>::rank);
      summary.add(a.dimension(0), a.multiplier(0),map1,Int2Type<LCM_INDEX_DIM(Map1)>());
      summary.add(a.dimension(1), a.multiplier(1),map2,Int2Type<LCM_INDEX_DIM(Map2)>());
      return return_type(a,summary);
    }
    template<class Map1,class Map2,class Map3> static return_type func(MultiMatrix<T,3>& a,const Map1& map1,const Map2& map2,const Map3& map3) {
      typedef IndirectMultiMatrix< T, MultiIndex<Map1,Map2,Map3>::rank > return_type;
      MultiSummary summary(MultiIndex<Map1,Map2,Map3>::rank);
      summary.add(a.dimension(0), a.multiplier(0),map1, Int2Type<LCM_INDEX_DIM(Map1)>());
      summary.add(a.dimension(1), a.multiplier(1),map2, Int2Type<LCM_INDEX_DIM(Map2)>());
      summary.add(a.dimension(2), a.multiplier(2),map3, Int2Type<LCM_INDEX_DIM(Map3)>());
      return return_type(a,summary);
    }
  };

 template<class T,size_t M> struct const_selectN<T,0,M> {
    typedef const IndirectMultiMatrix<T,M> return_type;

    template<class Map> static return_type func(const MultiMatrix<T,1>& a,const Map& map) {
      typedef const IndirectMultiMatrix<T,MultiIndex<Map>::rank> return_type;
      MultiSummary summary(MultiIndex<Map>::rank);
      summary.add(a.dimension(0), a.multiplier(0),map,Int2Type<LCM_INDEX_DIM(Map)>() );
      return return_type(a,summary);
    }
    template<class Map1,class Map2> static return_type func(const MultiMatrix<T,2>& a,const Map1& map1,const Map2& map2) {
      typedef const IndirectMultiMatrix<T,MultiIndex<Map1,Map2>::rank> return_type;
      MultiSummary summary(MultiIndex<Map1,Map2>::rank);
      summary.add(a.dimension(0), a.multiplier(0),map1,Int2Type<LCM_INDEX_DIM(Map1)>());
      summary.add(a.dimension(1), a.multiplier(1),map2,Int2Type<LCM_INDEX_DIM(Map2)>());
      return return_type(a,summary);
    }
    template<class Map1,class Map2,class Map3> static return_type func(const MultiMatrix<T,3>& a,const Map1& map1,const Map2& map2,const Map3& map3) {
      typedef const IndirectMultiMatrix<T,MultiIndex<Map1,Map2,Map3>::rank> return_type;
      MultiSummary summary(MultiIndex<Map1,Map2,Map3>::rank);
      summary.add(a.dimension(0), a.multiplier(0),map1, Int2Type<LCM_INDEX_DIM(Map1)>());
      summary.add(a.dimension(1), a.multiplier(1),map2, Int2Type<LCM_INDEX_DIM(Map2)>());
      summary.add(a.dimension(2), a.multiplier(2),map3, Int2Type<LCM_INDEX_DIM(Map3)>());
      return return_type(a,summary);
    }
  };

 template<class T,size_t N> class IndirectMultiMatrix {
   BaseList<T> data_;
   const ListList<size_t> lists_;
   size_t size_;
  
 public:
   template<size_t M> IndirectMultiMatrix(const MultiMatrix<T,M>& a,const MultiSummary& summary) :
     data_(a.size()-summary.offset,a.vector()+summary.offset), 
     lists_(summary.lists_)
   {
     size_=1;
     for (size_t n=lists_.size();n--;) 
       size_*=lists_.size(n);
   }

   typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;  

  friend std::ostream& operator<< <>(std::ostream&,const IndirectMultiMatrix<T,N>&);

  bool operator! () const { return false; }
  size_t dimension(size_t n) const { return lists_.size(n); }
  size_t size() const { return size_; }
  size_t rows() const {
    LCM_STATIC_CHECK( N==2, rows_applied_to_non2D_object );
    return lists_.size(0);
  }
  size_t cols() const { 
    LCM_STATIC_CHECK( N==2, rows_applied_to_non2D_object );
    return lists_.size(1);
  }

  template<typename T2> IndirectMultiMatrix<T,N>& operator= (const T2& a) {
    return assign(*this,a);
  }
  template<class T2> IndirectMultiMatrix<T,N>& operator+= (const T2& b) {
    return add_ip(*this,b);
  }
  template<class T2> IndirectMultiMatrix<T,N>& operator-= (const T2& b) {
    return subtract_ip(*this,b);
  }
  template<class T2> IndirectMultiMatrix<T,N>& operator*= (const T2& b) {
    return multiply_ip(*this,b);
  }
  template<class T2> IndirectMultiMatrix<T,N>& operator/= (const T2& b) {
    apply_ip(doesdivide_ip<T,LCM_VAL(T2)>(),*this,b);
    return *this; }

   class iterator : public ::std::iterator< ::std::forward_iterator_tag,T> {
     DirectSum_iterator<size_t> iter_;
     T* const data_;
     
   public:
     iterator(T* datav,const ListList<size_t>& indexv,bool isatendv) : data_(datav), iter_(indexv,isatendv) {}
     
     T& operator*() const { return data_[*iter_]; }
     iterator& operator++() { iter_++; return *this; }
     //Don't use postfix form if possible - horribly slow!
     iterator& operator++(int) { iterator tmp(*this); iter_++; return tmp; }
     
     bool operator== (const iterator& x) const { return (iter_==x.iter_); }
     bool operator!= (const iterator& x) const { return (iter_!=x.iter_); }
   };

   class const_iterator : public ::std::iterator< ::std::forward_iterator_tag,T,ptrdiff_t,const T*,const T& > {
     DirectSum_iterator<size_t> iter_;
     const T* const data_;
     
   public:
     const_iterator(const T* datav,const ListList<size_t>& indexv,bool isatendv) : data_(datav), iter_(indexv,isatendv) {}
     
     const T& operator*() const { return data_[*iter_]; }
     const_iterator& operator++() { iter_++; return *this; }
     //Don't use postfix form if possible - horribly slow!
     const_iterator operator++(int) { const_iterator tmp(*this); iter_++; return tmp; }
     
     bool operator== (const const_iterator& x) const { return (iter_==x.iter_); }
     bool operator!= (const const_iterator& x) const { return (iter_!=x.iter_); }
   };

  iterator begin() {
    return iterator(data_.vector(),lists_,false);
  }

  iterator end() {
    return iterator(data_.vector(),lists_,true);
  }

   const_iterator begin() const {
     return const_iterator(data_.vector(),lists_,false);
   }

   const_iterator end() const {
     return const_iterator(data_.vector(),lists_,true);
   }

   T& operator()(size_t r) {
     LCM_STATIC_CHECK( N==1, MultiMatrix_Dimensions_dont_match );
     return data_(lists_(0U,r));
   }
   T& operator()(size_t r,size_t s) {
     LCM_STATIC_CHECK( N==2, MultiMatrix_Dimensions_dont_match );
     return data_(lists_(0U,r)+lists_(1U,s));
   }
   T& operator()(size_t r,size_t s,size_t t) {
     LCM_STATIC_CHECK( N==3, MultiMatrix_Dimensions_dont_match );
     return data_(lists_(0U,r)+lists_(1U,s)+lists_(2U,t)); 
   }
   T& operator()(size_t r,size_t s,size_t t,size_t u) {
     LCM_STATIC_CHECK( N==4, MultiMatrix_Dimensions_dont_match );
     return data_(lists_(0U,r)+lists_(1U,s)+lists_(2U,t)+lists_(3U,u));
   }

   const T& operator()(size_t r) const {
     LCM_STATIC_CHECK( N==1, MultiMatrix_Dimensions_dont_match );
     return data_(lists_(0U,r)); 
   }
   const T& operator()(size_t r,size_t s) const { 
     LCM_STATIC_CHECK( N==2, MultiMatrix_Dimensions_dont_match );
     return data_(lists_(0U,r)+lists_(1U,s));
   }
   const T& operator()(size_t r,size_t s,size_t t) const { 
     LCM_STATIC_CHECK( N==3, MultiMatrix_Dimensions_dont_match );
     return data_(lists_(0U,r)+lists_(1U,s)+lists_(2U,t));
   }
   const T& operator()(size_t r,size_t s,size_t t,size_t u) const { 
     LCM_STATIC_CHECK( N==4, MultiMatrix_Dimensions_dont_match );
     return data_(lists_(0U,r)+lists_(1U,s)+lists_(2U,t)+lists_(3U,u));
   }
 };

  template<class T,size_t N> std::ostream& operator<< (std::ostream& ostr,const IndirectMultiMatrix<T,N>& a)
  {
#ifndef NDEBUG
    ostr << "Index: " << a.lists_ << '\n';
#endif
    print(a,ostr);
    return ostr;
  }

 template<typename T,size_t N> struct type_traits< IndirectMultiMatrix<T,N> > {
   static const bool trivialconstructor=false; 
   static const size_t dimensionality=N;
   static const size_t rank=0;
   typedef T value_type;
 };


#endif
 //LCM_SUPPRESS_VIEWS

//  template<size_t N, template<size_t> class CheckClass =LCM_CheckBoundsDefault>
template<size_t N, template<size_t> class CheckClass> //!< default removed 27/1/16
  class Indexer<N,CheckClass>::permuted_iterator : public ::std::iterator< ::std::bidirectional_iterator_tag,size_t> 
  {
    private:
    size_t dim[N];
    size_t mults[N];
    int cpos;
    int cindex[N];
    
    void recalc() { 
      cpos=0;
      for (int i=N;i--;)
	cpos+=mults[i]*cindex[i];
    }

   public:

    permuted_iterator(const Indexer<N>& a, const BaseList<size_t>& order)
    { 
      LCM_STATIC_CHECK(N>1, PermutedIterator_applied_to_lessthan_2D);
      size_t i;
      if (order.length()) {
	if (order.length()!=N)
	  throw Mismatch("Indexer::permuted_iterator");
	for (i=N;i--;)
	  cindex[i]=1; //mark not set
	for (i=N;i--;) {
	  const int which=order(N-i-1);
	  if (which>=N || (cindex[which]==0))
	    throw InvalidParameter("Indexer::permuted_iterator: bad permutation vector");
	  dim[i]=a.dim[which];
	  mults[i]=a.mults[which];
	  cindex[which]=0;
	}
      }
      else { //this actually corresponds to reversing order!
	for (i=N;i--;) {
	  dim[i]=a.dim[i];
	  mults[i]=a.mults[i];
	  cindex[i]=0;
	}
      }
      cpos=0;
    }
    //end iterator - don't need to fill in rest
    permuted_iterator(const Indexer<N>& a) { cpos=a.size(); }

    size_t index(size_t i) const {
      if (i>=N)
	throw BadIndex("permuted_iterator::index");
      return cindex[i];
    }
    size_t multiplier(size_t i) const {
      if (i>=N)
	throw BadIndex("permuted_iterator::multiplier");
      return mults[i];
    }

    void reset() {
      for (int i=N;i--;)
	cindex[i]=0;
      cpos=0;
    }

    bool operator!= (const permuted_iterator& x) const { return (cpos!=x.cpos); }
    bool operator== (const permuted_iterator& x) const { return (cpos==x.cpos); }

    size_t operator*() const { return cpos; }

     //Don't use postfix forms if possible!
    permuted_iterator operator++(int) { permuted_iterator tmp(*this); ++(*this); return tmp; }
    permuted_iterator operator--(int) { permuted_iterator tmp(*this); --(*this); return tmp; }

    permuted_iterator& operator++() {
      if (++(cindex[0])!=dim[0]) {
	cpos+=mults[0];
	return *this;
      }
      int i=0;
      while (i<N-1) {
	cindex[i++]=0;
	if (++(cindex[i])!=dim[i])
	  break;
      }
      recalc();
      return *this;
    }

    permuted_iterator& operator--() {
      if (--(cindex[0])>=0) {
	cpos-=mults[0];
	return *this;
      }
      int i=0;
      while (i<N-1) {
	cindex[i]=dim[i]-1;
	if (--(cindex[++i])>=0) {
	  recalc();
	  return *this;
	}
      }
      throw Failed("Indexer::permuted_iterator underflow");
    }
  };

  class MultiIterator {
  public:
    MultiIterator(const BaseList<size_t>&);
    size_t size() const { return dims_.size(); }
    size_t total() const { return total_; }
    bool next(BaseList<size_t>);
    bool next(List<size_t>& dest) {
      if (dest.empty())
	dest.create(size());
      return next( static_cast< BaseList<size_t> >(dest));
    }
    void reverse(BaseList<size_t>, size_t) const;
    void reverse(List<size_t>& dest, size_t index) const {
      dest.create(size());
      reverse( static_cast< BaseList<size_t> >(dest), index);
    }
    void reset();
  private:
    ScratchList<size_t,10> dims_;
    ScratchList<size_t,10> mults_;
    ScratchList<size_t,10> cur_;
    bool finished_;
    size_t total_;
  };
     
 template<class T,size_t N> class PermutedIterator : public ::std::iterator< ::std::bidirectional_iterator_tag,T> {
 public:
   PermutedIterator(T* sourcev, const Indexer<N>& indexerv, const BaseList<size_t>& order)
     : source_(sourcev), iter_(indexerv.permuted_begin(order)) {}

     PermutedIterator(const Indexer<N>& indexerv) 
       : iter_(indexerv.permuted_end()) {}

     //NB we don't bother to compare source - comparing iterators derived from different objects is undefined
     bool operator!= (const PermutedIterator& x) const { return (iter_==x.iter_); }
     bool operator== (const PermutedIterator& x) const { return (iter_!=x.iter_); }
     
     PermutedIterator& operator++() { ++iter_; return *this; }
     PermutedIterator operator++(int) { PermutedIterator<T,N> tmp(*this); ++iter_; return tmp; }
     PermutedIterator& operator--() { --iter_; return *this; }
     PermutedIterator operator--(int) { PermutedIterator<T,N> tmp(*this); --iter_; return tmp; }

     T& operator*() const { return source_[*iter_]; }

 private:
     T* source_;
     typename Indexer<N>::permuted_iterator iter_;
 };

//   template<size_t N> template<size_t M> Indexer<N>::Indexer(const Indexer<M>& a)
//   {
//     if (M<=N) {
//       size_t sofar=1;
//       size_t d=N;
//       for (size_t s=M;s--;) {
// 	mults[--d]=sofar;
// 	sofar*=(dim[d]=a.dimension(s));
//       }
//       for (;d--;) {
// 	mults[d]=sofar;
// 	dim[d]=1;
//       }
//     }
//     else {
//       for (size_t s=M-1;s>=N;s--) {
// 	if (a.dimension(s)>1)
// 	  throw Failed("Indexer<N>: can't shrink");
//       }
//       setdims(BaseList<size_t>(N,a.dims));
//     }
//   }

    
} //namespace libcmatrix

#endif
