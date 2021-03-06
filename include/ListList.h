#ifndef _ListList_h_
#define _ListList_h_

#include "List.h"
#include "ScratchList.h"

namespace libcmatrix {

template<class T> class ListList;
template<class T> std::ostream& operator << (std::ostream&, const ListList<T>&);

template<class T> class ListList {
public:
  explicit ListList(mxflag::tempflag tflag =mxflag::normal);
  ListList(int maxblocks_,int maxitems_);
  ListList(int maxblocks_,int maxitems_,const T&);
  ListList(const BaseList<size_t>& sizes, const T&);
  template<class T2> ListList(const BaseList<size_t>& sizes,const BaseList<T2>&);
  ListList(const ListList<size_t>&, const BaseList<T>&);

  ListList(const ListList<T>& a)
    : store_(a.store_)
    , str(a.str)
    {}

  template<class T2> explicit ListList(const ListList<T2>& a)
      : store_(a.row())
      , str(a.structure()) {}

  template<class T2> explicit ListList(const BaseList<T2>& a);
  explicit ListList(const BaseList<size_t>& sizes) {
    create(sizes);
  }

  void swap(ListList<T>& a) {
    store_.swap(a.store_); //do more 'risky' swap first
    str.swap(a.str);
  }
    
  void create(int n) {
    if (n!=size())
      throw Mismatch("ListList::create: blocks",n,size());
  }
  void create(int maxblocks_,int maxitems_);

  void create(const BaseList<size_t>& sizes) {
    createstr(sizes);
    store_.create(str.back()); }

  void create(const BaseList<size_t>& sizes, const T& v) {
    createstr(sizes);
    store_.create(str.back(),v); }

  template<typename T2> inline void duplicate_structure(const ListList<T2>& a) {
    str=a.structure();
    store_.create(a.items());
  }

  BaseList<T> front() {
    return operator()(size_t(0));
  }
  const BaseList<T> front() const {
    return operator()(size_t(0));
  }

  BaseList<T> back() { 
    return operator()(size()-1); }

  const BaseList<T> back() const {
    return operator()(size()-1); }

  void pop_back() { 
    if (str.size()==1)
      throw Failed("ListList<T>::pop_back");
    else {
      str.pop_back();
      store_.resize(str.back());
    }
  }

  template<class T2> void push_back_(Int2Type<0>, const T2& v)
  {
    store_.push_back(v);
    if (str.empty())
      str.push_back(1);
    else
      str.back()++;
  }

  template<class T2> void push_back_(Int2Type<1>, const T2& a) { 
    const size_t n=a.size();
    str.push_back(str.back()+n);
    store_.reserve(items()+n);//make sure there's enough room for a
    const typename T2::const_iterator lend(a.end());
    typename T2::const_iterator start(a.begin());    
    while (start!=lend) {
      store_.push_back(*start);
      ++start;
    }
  }
  
  void reserve(size_t blksv, size_t itemvs) { 
    str.reserve(blksv);
    store_.reserve(itemvs);
  }

  void push_back(size_t n, const T& v) {
    str.push_back(str.back()+n);
    store_.push_back(n,v);
  }
    
  template<class T2> void push_back(const T2& v)
  { push_back_(Int2Type< LCM_DIM(T2) >(),v); }

  void push_back() //start new block
  {
    str.push_back(str.empty() ? 0 : str.back());
  }
  
  typedef BaseList<T> value_type;
  typedef BaseList<T> reference;
  typedef const BaseList<T> const_reference;  
  typedef int difference_type;
  typedef size_t size_type;

  ListList<T>& operator=(const T& val) {
    store_=val;
    return *this; }

  template<class T2> ListList<T>& operator=(const ListList<T2>& a) { 
    this->str=a.structure(); 
    store_=a.row();
    return *this; }

  template<class T2> ListList<T>& operator= (const BaseList<T2>&);

  size_t size(size_t i) const { 
    if (i>=size())
      throw BadIndex("ListList<T>::size",i,size()); 
    return str(i+1)-str(i); }

  size_t offset(size_t i) const {
#if (BOUNDS_CHECK)
    if (i>=size())
      throw BadIndex("ListList<T>::offset",i,size()); 
#endif
    return str(i); }

  size_t offset(size_t i, size_t j) const {
#if (BOUNDS_CHECK)
    if (j>=size(i))
      throw BadIndex("ListList<T>::offset",j,size(i)); 
#endif
    return str(i)+j; }

  size_t size() const {
    return str.size()-1; } //needed so that ListList looks like 1D object

  size_t rows() const {
    return str.size()-1; }

  void clear() {
    str.create(1); //reset structure
    str=0;
    store_.clear();
  }

  bool empty() const { 
    return (size()==0); }
  
  BaseList<T> operator()(size_t i) {
    return BaseList<T>(size(i),store_.vector()+str(i)); }

  const BaseList<T> operator()(size_t i) const {
    return BaseList<T>(size(i),store_.vector()+str(i)); }

  BaseList<T> row(size_t i) {
    return BaseList<T>(size(i),store_.vector()+str(i)); }

  const BaseList<T> row(size_t i) const {
    return BaseList<T>(size(i),store_.vector()+str(i)); }

  ListList<T> operator()(const BaseList<size_t>&) const;
  T& operator()(size_t,size_t);
  const T& operator()(size_t,size_t) const;
  T& operator()(const std::pair<size_t,size_t>& a) { return (*this)(a.first,a.second); }
  const T& operator()(const std::pair<size_t,size_t>& a) const { return (*this)(a.first,a.second); }

  BaseList<T> row() { return store_; }
  const BaseList<T> row() const { return store_; }

  void allrows(List<BaseList<T> >&);
  inline List<BaseList<T> > allrows() { 
    List<BaseList<T> > res(mxflag::temporary);
    allrows(res);
    return res; }

  void allrows(List<const BaseList<T> >&) const;

  //! problematic without move semantics due to destruction of consts
  //  inline List<const BaseList<T> > allrows() const {
  //  List<const BaseList<T> > res(mxflag::temporary); 
  //  allrows(res);
  //  return res; }

  size_t items() const { 
    return store_.size(); }

  friend std::ostream& operator<< <> (std::ostream&, const ListList<T>&);

  template<class T2> ListList<T>& operator+= (const T2& b) {
    return add_ip(*this,b); }
  template<class T2> ListList<T>& operator-= (const T2& b) {
    return subtract_ip(*this,b); }
  template<class T2> ListList<T>& operator*= (const T2& b) {
    return multiply_ip(*this,b); }
  template<class T2> ListList<T>& operator/= (const T2& b) {
    apply_ip(doesdivide_ip<value_type,LCM_VAL(T2)>(),*this,b);
    return *this; }
  template<class T2> ListList<T>& operator&= (const T2& b) {
    return premultiply_ip(*this,b); }

  //Specialisation for constants (can treated as flat vector)
  ListList<T>& operator+= (const T& b) {
    store_.operator+=(b);
    return *this; }

  ListList<T>& operator-= (const T& b) {
    store_.operator-=(b);
    return *this; }

  ListList<T>& operator*= (const T& b) {
    store_.operator*=(b);
    return *this; }

  ListList<T>& operator/= (const T& b) {
    store_.operator/=(b);
    return *this; }

  class iterator;
  class const_iterator;

   iterator begin() { 
     return iterator(*this,0); }

   iterator end() {
     return iterator(*this,size()); }

   const_iterator begin() const { 
     return const_iterator(*this,0); }

   const_iterator end() const {
     return const_iterator(*this,size()); }

   const List<size_t>& structure() const { 
     return str; }

  bool operator== (const ListList<T>& b) const {
    return (store_==b.store_) && (str==b.str);
  }

  bool operator!= (const ListList<T>& b) const {
    return (store_!=b.store_) || (str!=b.str);
  }

 private:
  void createstr(const BaseList<size_t>&);
  void dump() const { std::cout << *this; }

  List<T> store_; //List (unlike DynamicList) allows ListList to be used like a stack
  List<size_t> str;
};

 template<typename T> struct type_traits< ListList<T> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef BaseList<T> value_type;
 };


  template<class T> class ListList<T>::iterator
  : public ::std::iterator< ::std::bidirectional_iterator_tag,BaseList<T>,ptrdiff_t,BaseList<T>*,BaseList<T>& > 
  {
   public:
     iterator(ListList<T>& obj_,size_t ptr) : 
       pobj(&obj_), cptr(ptr) {}

     BaseList<T>& operator*() { 
       tmp_.create((*pobj)(cptr));
       return tmp_; }

     bool operator== (const iterator& x) const {
       return (x.cptr==cptr); }
     bool operator!= (const iterator& x) const { 
       return (x.cptr!=cptr); }

     iterator& operator++() { 
       cptr++; 
       return *this; }

     iterator operator++(int) {
       iterator tmp(*this);
       cptr++;
       return tmp; }

     iterator& operator--() { 
       cptr--;
       return *this; }

     iterator operator--(int) {
       iterator tmp(*this);
       cptr--;
       return tmp; }

   private:
     //iterator() {}
     ListList<T>* const pobj;
     BaseList<T> tmp_;
     size_t cptr;
   };

  template<class T> class ListList<T>::const_iterator
  : public ::std::iterator< ::std::random_access_iterator_tag,BaseList<T>,ptrdiff_t,const BaseList<T>*,const BaseList<T> >
  {
   public:
     const_iterator(const ListList<T>& obj_,size_t ptr) :
       pobj(&obj_), cptr(ptr) {}

     const BaseList<T> operator*() const {
       return (*pobj)(cptr); }
     bool operator== (const const_iterator& x) const { 
       return (x.cptr==cptr); }
     bool operator!= (const const_iterator& x) const { 
       return (x.cptr!=cptr); }
     bool operator> (const const_iterator& x) const {
       return (x.cptr>cptr); }
     bool operator< (const const_iterator& x) const {
       return (x.cptr<cptr); }
     bool operator>= (const const_iterator& x) const {
       return (x.cptr>=cptr); }
     bool operator<= (const const_iterator& x) const { 
       return (x.cptr<=cptr); }

     const_iterator& operator++() { 
       cptr++;
       return *this; }

     const_iterator operator++(int) {
       const_iterator tmp(*this);
       cptr++;
       return tmp; }

     const_iterator& operator--() { 
       cptr--;
       return *this; }

     const_iterator operator--(int) {
       const_iterator tmp(*this);
       cptr--; 
       return tmp; }

     const_iterator& operator+=(ptrdiff_t i) {
       cptr+=i;
       return *this; }

     const_iterator& operator-=(ptrdiff_t i) { 
       cptr-=i;
       return *this; }

     const_iterator operator+ (ptrdiff_t i) const {
       const_iterator tmp(*this);
       return tmp+=i; }

     const_iterator operator- (ptrdiff_t i) const { 
       const_iterator tmp(*this);
       return tmp-=i; }

     ptrdiff_t operator- (const_iterator& x) const {
       return cptr-x.cptr; }

     const BaseList<T> operator[] (ptrdiff_t i) const { 
       return (*pobj)(cptr+i); }

   private:
     //     const_iterator() {}
     const ListList<T>* const pobj;
     size_t cptr;
   };

 template<class T,class T2> void duplicate_structure(ListList<T>& a, const T2& b) 
 { a.duplicate_structure(b); }

 template<class T> class DirectSum_iterator : public ::std::iterator< ::std::forward_iterator_tag,T> {
   const ListList<T>& a;
   ScratchList<size_t> cur;
   T* lastrow;
   size_t& lastindex;
   size_t lastsize;
   T cursum;
   bool isatend;

   void resetcursum() {
     cursum=0;
     for (size_t i=cur.size()-1;i--;)
       cursum+=a(i,cur(i));
   }

   DirectSum_iterator& operator=(const DirectSum_iterator&)   
     { throw InternalError("DirectSum_iterator: can't do this!"); }

 public:
   DirectSum_iterator(const ListList<T>& av,bool isatendv =false) :
     a(av),
     cur(a.size(),0),
     lastrow(a.back().vector()), 
     lastindex(cur.back()),
     lastsize(a.size(a.size()-1)),
     cursum(0),
     isatend(isatendv) {}

   DirectSum_iterator(const DirectSum_iterator& b) : 
     a(b.a),
     cur(b.cur),
     lastrow(b.lastrow),
     lastindex(cur.back()),
     lastsize(b.lastsize),
     cursum(b.cursum), 
     isatend(b.isatend) {}
   
   T operator*() const { 
     return cursum+lastrow[lastindex]; }

   DirectSum_iterator& operator++() { 
     if (++lastindex!=lastsize) 
       return *this;
     lastindex=0;
     for (size_t i=cur.size()-1;i--;) {
       if (++cur(i)!=a.size(i)) {
	 resetcursum();
	 return *this;
       }
       cur(i)=0;
     }
     isatend=true;
   }
   //Don't use postfix form if possible - horribly slow!
   DirectSum_iterator operator++(int) { 
     DirectSum_iterator tmp(*this);
     ++(*this); 
     return tmp; }

   bool operator== (const DirectSum_iterator& x) const { 
     if (isatend!=x.isatend)
       return false; //quick check for end
     return isatend ? true : (cur==x.cur); 
   }
   bool operator!= (const DirectSum_iterator& x) const { 
     if (isatend!=x.isatend) 
       return true;
     return isatend ? false : (cur!=x.cur); 
   }
 };

 template<class T> class Matrix;

  ListList<size_t> find_blocks(const Matrix<complex> &,double tol =0.0);
  template<typename T,class F> ListList<size_t> find_blocks(const BaseList<T>&, const F&);
  template<typename T> ListList<size_t> find_blocks(const BaseList<T>& a) { return find_blocks(a,std::equal_to<T>()); }
  ListList<size_t> find_blocks(const BaseList<double>&, double tol);
  
 template<class T,class T2> ListList<T> operator+ (const ListList<T>& a, const T2& b)
   {
     return ListList<T>(a)+=b;
   }

template<class T,class T2> ListList<T> operator- (const ListList<T>& a, const T2& b)
   {
     return ListList<T>(a)-=b;
   }
            
// Implementation details below

 template<class T> ListList<T>::ListList(mxflag::tempflag tflag) :
   store_(tflag),
   str(1,size_t(0),tflag) {}

   template<class T> ListList<T>::ListList(int maxblocks_,int maxitems_) :
     store_(maxitems_),
     str(maxblocks_+1)
       {
	 store_.resize(0);
	 str.create(1,size_t(0));
       } //allocate space, but leave empty

     template<class T> void ListList<T>::create(int maxblocks_,int maxitems_)
       {
	 store_.reserve(maxitems_);
	 store_.resize(0);
	 str.reserve(maxblocks_+1);
	 str.create(1,size_t(0));
       }

  template<class T> ListList<T>::ListList(int maxblocks_,int maxitems_,const T& v) : 
    store_(maxitems_,v),
    str(maxblocks_+1)
  {
    store_.create(0,v);
    str.create(1,size_t(0));
  }

  template<class T> ListList<T>::ListList(const BaseList<size_t>& sizes, const T& v) {
    createstr(sizes);
    store_.create(str.back(),v);
  }

  template<class T> template<class T2>
  ListList<T>::ListList(const BaseList<size_t>& sizes,const BaseList<T2>& a) 
    : store_(a)
  {
    createstr(sizes);
    if (str.back()!=a.size())
      throw Mismatch("ListList<T>(): items",str.back(),a.size());
  }

 template<class T> void ListList<T>::createstr(const BaseList<size_t>& sizes)
   {
     const size_t blks=sizes.size();
     str.create(blks+1);
     
     size_t point=0;
     size_t blk;
     for (blk=0;blk<blks;blk++) {
       str(blk)=point;
       point+=sizes(blk);
     }
     str(blk)=point;
   }

 template<class T> template<class T2> ListList<T>::ListList(const BaseList<T2>& a) :
   str(1,size_t(0),mxflag::normal)
{
  LCM_STATIC_CHECK( LCM_DIM(T2)==1, ListList_source_not_list_of_1D_vectors);
  for (size_t i=0;i<a.size();i++)
    push_back(a(i));
}

 template<class T> template<class T2> ListList<T>& ListList<T>::operator= (const BaseList<T2>& a)
{
  LCM_STATIC_CHECK( LCM_DIM(T2)==1, ListList_source_not_list_of_1D_vectors);
  str.create(1,size_t(0));
  create(0);
  for (size_t i=0;i<a.size();i++)
    push_back(a(i));
  return *this;
}

template<class T> ListList<T>::ListList(const ListList<size_t>& list, const BaseList<T>& z)
{
  create(list.structure());
  BaseList<T> dasvec=row();
  const BaseList<size_t> iasvec=list.row();
  for (size_t i=iasvec.size();i--;) 
    dasvec(i)=z(iasvec(i));
}

template<class T> ListList<T> ListList<T>::operator()(const BaseList<size_t>& which) const
{
  const size_t blks=which.size();
  List<size_t> sizes(blks);
  size_t i;
  for (i=blks;i--;) 
    sizes(i)=size(which(i));
  
  ListList<T> d;
  d.create(sizes);
  for (i=blks;i--;)
    d(i)=(*this)(which(i));
  return d;
}

#if (BOUNDS_CHECK)
template<class T> T& ListList<T>::operator()(size_t i,size_t j) {
  if (j<size(i))
    return store_(str(i)+j);
  else
    throw BadIndex("ListList::()",j,size(i)); }

template<class T> const T& ListList<T>::operator()(size_t i,size_t j) const {
  if (j<size(i))
    return store_(str(i)+j);
  else
    throw BadIndex("ListList::()",j,size(i)); }

#else
template<class T> T& ListList<T>::operator()(size_t i,size_t j) { 
  return store_(str(i)+j); }

template<class T> const T& ListList<T>::operator()(size_t i,size_t j) const { 
  return store_(str(i)+j); }
#endif

template<class T> std::ostream& operator<< (std::ostream &ostr,const ListList<T> &list)
{
  const size_t n=list.size();
  if (n==0)
    return ostr << "[Empty]";
  for (size_t i=0;i<n;i++) {
    ostr << list(i);
    if (i!=n-1) ostr << ',';
  }
  return ostr;
}

template<class T> void ListList<T>::allrows(List<BaseList<T> >& res)
{
  size_t n=size();
  res.create(n);
  for (;n--;) 
    res(n).create(operator()(n));
}

template<class T> void ListList<T>::allrows(List<const BaseList<T> >& res) const
{
  size_t n=size();
  List< BaseList<T> >& ncres(reinterpret_cast< List< BaseList<T> >& >(res));
  ncres.create(n);
  for (;n--;) {
    const BaseList<T> crow(res(n));
    BaseList<T>& trow(const_cast< BaseList<T>& >(crow));
    trow.create(operator()(n));
  }
}

  template<typename T,class F> ListList<size_t> find_blocks(const BaseList<T>& a,const F& matchobj)
{
  const size_t n=a.length();
  if (n==0)
    throw Undefined("find_blocks");

  ScratchList<size_t> tmpsize(n);
  ScratchList<size_t> ind(n);
  ScratchList<bool> used(n,false);

  size_t i;
  size_t index=0;
  size_t blockcount=0;
  size_t searchstart=0;

  for (;;) {
    for (i=searchstart;i<n && used(i);i++);
    if (i==n)
      break;
    searchstart=i+1;
    const T& value(a(i));
    const size_t indexstore=index;

    used(i)=true;
    ind(index++)=i;

    while (++i<n) {
      if (matchobj(a(i),value)) {
	used(i)=true;
	ind(index++)=i;
      }
    }
    tmpsize(blockcount)=index-indexstore;
    blockcount++;
  }
  assert(index==n);

  return ListList<size_t>(tmpsize.truncate(blockcount),ind);
}

#ifdef LCM_USE_EXTERNTEMPLATE
#ifdef LCM_LEVEL3_INSTANT
#define LCM_L3_EXTERN
#else
#define LCM_L3_EXTERN extern
#endif
LCM_L3_EXTERN template class ListList<double>;
LCM_L3_EXTERN template class ListList<complex>;
LCM_L3_EXTERN template class ListList<size_t>;
#endif

} //namespace libcmatrix

#endif
