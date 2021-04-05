#ifndef _Tensor_h_
#define _Tensor_h_

#include "DynamicList.h"
#include "Warnings.h"

namespace libcmatrix {

  struct BadRank : public MatrixException, public std::invalid_argument {
    virtual ~BadRank() throw() {};
    explicit BadRank(const char* errm_)
      : MatrixException("Invalid/missing tensor rank",errm_),
	std::invalid_argument(errm_) {}
    BadRank(const char* errm_, int rank_, int maxrank_);
  };

  extern Warning<> tensor_ignore_warning;

  //class complex;
template<class T> class Tensor;
template<class T> std::ostream& operator << (std::ostream &,const Tensor<T> &);

template<class T> class Tensor {
public:

  Tensor() : _rank(-1) {}
    explicit Tensor(int rank, mxflag::tensflag flag =mxflag::maximum);
    Tensor(int rank,const T&, mxflag::tensflag flag =mxflag::maximum);
  Tensor(const Tensor& a);
  Tensor& operator=(const Tensor&);
  ~Tensor();
  void clear();

  void swap(Tensor<T>&);
  inline void create(int,mxflag::tensflag =mxflag::none);

  int rank() const { return _rank; }
  int max_rank() const;

  template<class T2> Tensor& operator= (const T2&);

  Tensor<T>& operator*= (double);
  Tensor<T>& operator/= (double v) { return((*this)*=1.0/v); }

  Tensor<T>& operator+= (const Tensor &);
  Tensor<T>& operator-= (const Tensor &);
  bool operator! () const { return (rank()<0); }

  const BaseList<T> operator()(int lrank) const { 
    if (have_rank(lrank))
      return BaseList<T>(2*lrank+1,_entries[lrank]-lrank);
    throw BadRank("Tensor<T>()",lrank,_rank);
  } 
  BaseList<T> operator()(int lrank) { 
    if (have_rank(lrank))
      return BaseList<T>(2*lrank+1,_entries[lrank]-lrank);
    throw BadRank("Tensor<T>()",lrank,_rank);
  } 

  const T* vector(int lrank) const { 
    if (have_rank(lrank))
      return _entries[lrank];
    throw BadRank("Tensor<T>::vector",lrank,_rank);
  } 
  T* vector(int lrank) { 
    if (have_rank(lrank))
      return _entries[lrank];
    throw BadRank("Tensor<T>::vector",lrank,_rank);
  } 

#if (BOUNDS_CHECK)
  const T& operator()(int lrank,int order_) const {
    check_rankorder(lrank,order_);
    if (!_entries[lrank])
      throw BadRank("Tensor<T>()",lrank,_rank);
    return _entries[lrank][order_];
  }
  T& operator()(int lrank,int order_) {
    check_rankorder(lrank,order_);
    if (!_entries[lrank])
      throw BadRank("Tensor<T>()",lrank,_rank);
    return _entries[lrank][order_];
  }
#else
  const T& operator()(size_t lrank,int order_) const { return _entries[lrank][order_]; }
  T& operator()(size_t lrank,int order_) { return _entries[lrank][order_]; }
#endif

  bool have_rank(int) const;
  void ensure_rank(int);
  void ensure_rank(int, const T&);
  void clear(int);
  void check_rankorder(int,int) const;

  friend std::ostream& operator << <>(std::ostream&, const Tensor<T>&);
  static int max_rank(const Tensor<T>&, const Tensor<T>&);

  void negate();

 private:
  DynamicList<T> store_;
  T** _entries;
  int _rank;

  void print(std::ostream& =std::cout) const;
  void printsimple(std::ostream& =std::cout) const;

  void make_ranks(mxflag::tensflag);
  void create_tensor(int);
  void create_structure(int);
};

template<> std::ostream& operator << (std::ostream &ostr,const Tensor<complex> &a);
template<> std::ostream& operator << (std::ostream &ostr,const Tensor<double> &a);

template<class T> std::ostream& operator << (std::ostream &ostr,const Tensor<T> &a) { a.print(ostr); return ostr; }

// Only scaling by double is permitted
template<class T> Tensor<T> operator* (const Tensor<T> &,double);
template<class T> Tensor<T> operator/ (const Tensor<T> &a,double v) { return a*(1.0/v); }
template<class T> Tensor<T> operator* (double b,const Tensor<T> &a) { return a*b; }

template<class T> Tensor<T> operator- (const Tensor<T>&);
template<class T> Tensor<T> operator+ (const Tensor<T>&, const Tensor<T>&);
template<class T> Tensor<T> operator- (const Tensor<T>&, const Tensor<T>&);

// Return smallest of two ranks compatible with active ranks
 template<class T> int Tensor<T>::max_rank(const Tensor<T>& a,const Tensor<T>& b)
{
  const int rank_a=a.rank();
  const int rank_b=b.rank();
  const int min_rank= (rank_a>rank_b) ? rank_b : rank_a;  

  return (a.max_rank()>min_rank) ? rank_a : rank_b;
}

  template<typename T> bool operator==(const Tensor<T>&, const Tensor<T>&);
  template<typename T> bool operator!=(const Tensor<T>& a, const Tensor<T>& b) { return !(a==b); }
  
// implementation details only below

  template<class T> bool Tensor<T>::have_rank(int lrank) const {
    if (lrank<0)
      throw BadRank("have_rank",lrank,_rank);
    return (lrank<=_rank) && (_entries[lrank]!=NULL);
  }
   
  template<class T> void Tensor<T>::ensure_rank(int lrank) { 
    if ((lrank>=0) && (lrank<=_rank))
      _entries[lrank]=store_.vector()+lrank*(lrank+1);
    else 
      throw BadRank("ensure_rank",lrank,_rank);
  }

  template<class T> void Tensor<T>::ensure_rank(int lrank, const T& def) { 
    if ((lrank>=0) && (lrank<=_rank)) {
      T* vec=(_entries[lrank]=store_.vector()+lrank*(lrank+1));
      for (int i=-lrank;i<=lrank;i++)
	vec[i]=def;
    }
    else
      throw BadRank("ensure_rank",lrank,_rank);
  }

  template<class T> void Tensor<T>::clear(int lrank) {
    if ((lrank>=0) && (lrank<=_rank))
      _entries[lrank]=NULL;
    else 
      throw BadRank("clear",lrank,_rank); }

 template<class T> void Tensor<T>::create_tensor(int maxr)
{
  if (maxr<0)
    throw BadRank("create_tensor");
  store_.create((maxr+1)*(maxr+1));
  create_structure(maxr);
}

  template<class T> void Tensor<T>::create_structure(int newrank)
  {
    _entries= new T* [newrank+1];
    _rank=newrank; //only set if claim worked
    for (int i=0;i<=_rank;)
      _entries[i++]=NULL;
  }

template<class T> Tensor<T> operator* (const Tensor<T>& a,double v)
{
  Tensor<T> b(a);
  b*=v;
  return b;
}

template<class T> Tensor<T> operator- (const Tensor<T>& a)
{
  Tensor<T> d(a);
  d.negate();
  return d;
}

template<class T> void Tensor<T>::negate()
{
  for (int i=rank();i>=0;i--) {
    if (have_rank(i))
      negate_ip((*this)(i));
  }
}

 template<typename T> struct doesnegate_ip< Tensor<T> > {
   typedef Tensor<T> argument_type;
   void operator()(Tensor<T>& a) const { a.negate(); }
 };

template<class T> Tensor<T> operator+ (const Tensor<T>& a,const Tensor<T>& b)
{
  const int rank=Tensor<T>::max_rank(a,b);

  Tensor<T> c(rank,mxflag::none);

  for (int l=rank;l>=0;l--) {
    if (a.have_rank(l)) {
      c.ensure_rank(l);
      if (b.have_rank(l)) {
	BaseList<T> dest(c(l));
	add(dest,a(l),b(l));
      }
      else
	c(l)=a(l);
    }
    else {
      if (b.have_rank(l)) {
	c.ensure_rank(l);
	c(l)=b(l);
      }
    }
  }
  return c;
}

template<class T> inline void Tensor<T>::check_rankorder(int lrank,int order_) const
{
  if (lrank<0 || lrank>_rank)
    throw BadRank("check_rankorder",lrank,_rank);
  if (order_<-lrank || order_>lrank)
    throw BadIndex("check_rankorder",order_,lrank);
}

 template<class T> inline void Tensor<T>::create(int lrank,mxflag::tensflag flag)
{
  if (_rank!=lrank) {
    Tensor<T> tmp(lrank,flag); //this is exception safe
    swap(tmp);
  }
}

  template<typename T> bool operator== (const Tensor<T>& a, const Tensor<T>& b)
  {
    int l=a.max_rank();
    if (l!=b.max_rank())
      return false;
    for (;l>=0;l--) {
      const bool isactive(a.have_rank(l));
      if (isactive ^ b.have_rank(l))
	return false;
      if (isactive && (a(l)!=b(l)))
	return false;
    }
    return true;
  }

template<class T> Tensor<T> operator- (const Tensor<T>& a, const Tensor<T>& b)
{
  const int ran=Tensor<T>::max_rank(a,b);

  Tensor<T> c(ran,mxflag::none);

  for (int l=ran;l>=0;l--) {
    if (a.have_rank(l)) {
      c.ensure_rank(l);
      if (b.have_rank(l)) {
	BaseList<T> dest(c(l));
	subtract(dest,a(l),b(l));
      }
      else
	c(l)=a(l);
    }
    else {
      if (b.have_rank(l)) {
	c.ensure_rank(l);
	BaseList<T> dest(c(l));
	negate(dest,b(l));
      }
    }
  }
  return c;
}

template<class T> Tensor<T>::Tensor(int maxr,mxflag::tensflag flag)
{
  create_tensor(maxr);
  make_ranks(flag);
}
 
  template<class T> Tensor<T>::Tensor(int maxr,const T& v, mxflag::tensflag flag)
  {
    create_tensor(maxr);
    if (flag==mxflag::none)
      tensor_ignore_warning.raise();
    else {
      make_ranks(flag);
      (*this)=v;
    }
  }
 
template<class T> void Tensor<T>::make_ranks(mxflag::tensflag flag)
{
  switch (flag) {
  case mxflag::none: case mxflag::temporary:
    break;
  case mxflag::maximum:
    ensure_rank(_rank);
    break;
  case mxflag::all:
    for (int i=0;i<=_rank;i++)
      ensure_rank(i);
    break;
  default:
    throw InvalidParameter("Tensor<T>::docreate");
  }
}

  template<class T> Tensor<T>::Tensor(const Tensor<T>& a) : store_(a.store_) {
    create_structure(a._rank);
    for (int i=0;i<=_rank;i++) {
      if (a._entries[i])
	ensure_rank(i);
    }
  }

  template<class T> void Tensor<T>::swap(Tensor<T>& a) {
    store_.swap(a.store_);
    ::std::swap(_entries,a._entries);
    ::std::swap(_rank,a._rank);
  }

template<class T> template<class T2> Tensor<T>& Tensor<T>::operator= (const T2& val)
{
  if (_rank<0)
    throw Undefined("Tensor<T>::=");
  for (int i=0;i<=_rank;i++) {
    if (have_rank(i))
      (*this)(i)=val;
  }
  return *this;
}
 
template<class T> int Tensor<T>::max_rank() const
{
  int maxran=rank();
  while (maxran>=0 && !_entries[maxran]) maxran--;
  return maxran;
}

template<class T> Tensor<T>& Tensor<T>::operator=(const Tensor<T>& a)
{
  if (a.rank()<0) {
    clear();
    return *this;
  }
  if (this !=&a) {
    if (_rank < a.max_rank()) {
      clear();
      create_tensor(a._rank);
    }
    const int userank=std::max(a._rank,_rank); //!< modified 13/5/07 to use *maximum*
    // (a._rank>_rank) ? _rank : a._rank;

    if (userank>=0) {
      for (int i=0;i<=userank;i++) {
	if (a.have_rank(i)) {
	  ensure_rank(i);
	  for (int j=-i;j<=i;j++)
	    _entries[i][j]=a._entries[i][j];
	}
	else
	  clear(i);
      }
    }
  }
  return *this;
}

 template<class T> Tensor<T>::~Tensor()
   {
     if (_rank>=0)
       delete[] _entries;
   }

template<class T> void Tensor<T>::clear()
{
  store_.clear();
  this->~Tensor();
  _rank=-1;
}
    
template<class T> void Tensor<T>::printsimple(std::ostream &ostr) const
{
  if (rank()<0) {
    ostr << "<undefined>" << std::endl;
    return;
  }

  const Tensor<T>& a=*this;
  for (int i=0;i<=a.rank();i++) {
    if (!a.have_rank(i)) continue;

    ostr << "Rank " << i << ": ";
    for (int j=-i;j<=i;j++)
      ostr << a(i,j) << ' ';
    ostr << std::endl;
  }
}
    
template<class T> void Tensor<T>::print(std::ostream &ostr) const
{
  if (rank()<0) {
    ostr << "<undefined>" << std::endl;
    return;
  }

  const Tensor<T> &a=*this;
  for (int i=0;i<=a.rank();i++) {
    if (!a.have_rank(i)) continue;

    for (int j=-i;j<=i;j++)
      ostr << "Component " << i << "," << j << '\n' << a(i,j) << '\n';
    ostr << std::endl;
  }
}
    
template<class T> Tensor<T>& Tensor<T>::operator*= (double v)
{
  Tensor<T> &a=*this;

  for (int l=rank();l>=0;l--) {
    if (a.have_rank(l)) {
      for (int m=-l;m<=l;m++) a(l,m)*=v;
    }
  }
  return a;
}

template<class T1,class T2,class T3> void mla(Tensor<T1>& d, const T2& v, const Tensor<T1>& a)
{
  if (a.max_rank()>d.rank()) {
    Tensor<T1> tmp(a);
    tmp*=v;
    tmp+=d;
    d.swap(tmp);
    return;
  }
  for (int l=d.rank();l>=0;l--) {
    if (!a.have_rank(l))
      continue;
    if (d.have_rank(l))
      mla(d(l),v,a(l));
    else {
      d.ensure_rank(l);
      multiply(d(l),v,a(l));
    }
  }
}

template<class T> Tensor<T>& Tensor<T>::operator+= (const Tensor<T>& a)
{
  // If rank of a is bigger then need to expand original tensor
  if (a.max_rank()>rank()) {
    Tensor<T> tmp(a);
    tmp+=*this;
    swap(tmp);
    return *this;
  }

  for (int l=rank();l>=0;l--) {
    if (!a.have_rank(l))
      continue;
    if (have_rank(l))
      (*this)(l)+=a(l);
    else {
      ensure_rank(l);
      (*this)(l)=a(l);
    }
  }
  return *this;
}

template<class T> Tensor<T>& Tensor<T>::operator-= (const Tensor<T>& a)
{
  if (a.max_rank()>rank()) {
    Tensor<T> tmp(*this-a);
    swap(tmp);
    return *this;
  }

  for (int l=rank();l>=0;l--) {
    if (!a.have_rank(l)) continue;
    if (have_rank(l))
      (*this)(l)-=a(l);
    else {
      ensure_rank(l);
      negate((*this)(l),a(l));
    }
  }
  return *this;
}

} //namespace libcmatrix

#endif
