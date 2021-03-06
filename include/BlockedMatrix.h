#ifndef _BlockedMatrix_h_
#define _BlockedMatrix_h_

#include "Matrix.h"
#include "ListList.h"
#include "simple_counter.h"

namespace libcmatrix {

  template<class T> class BlockedMatrix;

  template<size_t> struct lcm_post_multiply_ {
    template<typename T,typename T2> static void apply(BlockedMatrix<T>& a, const T2& b);
  };
  template<> struct lcm_post_multiply_<0> {
    template<typename T,typename T2> static void apply(BlockedMatrix<T>& a, const T2& b) { a.row()*=b; }
  };

  template<class T> class BlockedMatrix : public DynamicList< Matrix<T> > {
  public:
    BlockedMatrix(mxflag::tempflag tflag =mxflag::normal) : store_(tflag), active_(0) {} 

      BlockedMatrix(const Matrix<T>&, mxflag::tempflag tflag =mxflag::normal);
    BlockedMatrix(const BaseList<size_t>&);
    BlockedMatrix(const BaseList<size_t>&, const BaseList<size_t>&);
    BlockedMatrix(const BaseList<size_t>&, const T&);
    BlockedMatrix(const BaseList<size_t>&, const BaseList<size_t>&, const T&);
    template<class T2> BlockedMatrix(const BaseList< Matrix<T2> >& a) { *this=a; }
    template<class T2> BlockedMatrix(const BlockedMatrix<T2>&);

    BlockedMatrix(const BlockedMatrix&); //Need special forms
    BlockedMatrix<T>& operator= (const BlockedMatrix<T>&);
    template<class T2> BlockedMatrix<T>& operator= (const BlockedMatrix<T2>&);
    template<class T2> BlockedMatrix<T>& operator= (const BaseList< Matrix<T2> >&);
    BlockedMatrix<T>& operator= (const T& v) {
      store_=v;
      return *this;
    }

    typedef Matrix<T> value_type;

    void clear() {
      store_.clear();
      active_=0;
      DynamicList< Matrix<T> >::clear();
    }

    template<class T2> void duplicate_structure(const T2& source, bool transpose =false)
    {
      store_.create(getsize(source));
      buildmatrices(source,transpose);
    }

    template<class T2> void duplicate_structure(const ListList<T2>&);
    template<class T2> void duplicate_structure(const BaseList< Matrix<T2> >&);
    template<class T2> void duplicate_structure(const Matrix<T2>& a) {
      create(ExplicitList<1,size_t>(a.rows()),ExplicitList<1,size_t>(a.cols()));
    }

    template<class T2> BlockedMatrix& operator+=(const T2& a);
    template<class T2> BlockedMatrix& operator+=(const BlockedMatrix<T2>& a);
    template<class T2> BlockedMatrix& operator-=(const BlockedMatrix<T2>& a);
    template<class T2> BlockedMatrix& operator*=(const T2& a) { lcm_post_multiply_<LCM_DIM(T2)>::apply(*this,a); return *this; }
    template<class T2> BlockedMatrix& operator&=(const T2& a);
    template<class T2> void emultiply(const T2& a) { emultiply_ip(*this,a); }

    //overwrites create present in List
    void create(const BaseList<size_t>&);
    void create(const BaseList<size_t>&, const BaseList<size_t>&); 
    void create(const BaseList<size_t>&, const T&);
    void create(const BaseList<size_t>&, const BaseList<size_t>&, const T&); 

    //Not strictly necessary but enables uniformity with more complex blocked types
    size_t rows(size_t i) const
    { return (*this)(i).rows(); }

    size_t cols(size_t i) const
    { return (*this)(i).cols(); }

    size_t active() const { return active_; }

    BaseList<T> row() const {
      if (this->empty()) 
	throw Undefined("BlockedMatrix<T>::row"); 
      return store_;
    }

    void swap(BlockedMatrix& a) {
      DynamicList< Matrix<T> >::swap(a);
      store_.swap(a.store_); //mess if exception thrown here...
      ::std::swap(active_,a.active_);
    }
    void swap(Matrix<T>&);

    bool ismatching(const BaseList<size_t>&) const;
    bool ismatching(const BaseList<size_t>&, const BaseList<size_t>&) const;
    
    void unitary_simtrans(const BlockedMatrix<complex>&)
    { LCM_STATIC_ERROR( unitary_simtrans_undefined_for_this_BlockedMatrix ); }
    void unitary_isimtrans(const BlockedMatrix<complex>&)
    { LCM_STATIC_ERROR( unitary_simtrans_undefined_for_this_BlockedMatrix ); }
    void unitary_simtrans(const ListList<complex>&)
    { LCM_STATIC_ERROR( unitary_simtrans_undefined_for_this_BlockedMatrix ); }
    void unitary_isimtrans(const ListList<complex>&)
    { LCM_STATIC_ERROR( unitary_simtrans_undefined_for_this_BlockedMatrix ); }

    void identity();

    template<class T2> void set(const BlockedMatrix<T2>& a) {
      duplicate_structure(a);
      assign(store_,a.row());
    }

    template<class T2> void set(const ListList<T2>& a) { full(*this,a); }

    template<class T2> void set(const Matrix<T2>& a) {
      duplicate_structure(a);
      assign(store_,a.row());
    }

    //Use same allocator as Matrix<T>
    typedef  DynamicList<T,typename matrix_traits<T>::allocator> storage_t;
    const storage_t& raw_storage() const { return store_; }

    std::ostream& dump_content(std::ostream&) const; //!< dump raw content

  private:
    storage_t store_;
    size_t active_;

    void dump() const { ::std::cout << *this; }

    template<class T2> static size_t getsize(const BaseList< Matrix<T2> >&);
    template<class T2> static size_t getsize(const BlockedMatrix<T2>& a) { return a.row().size(); }

    template<class T2> void buildmatrices(const BlockedMatrix<T2>&, bool =false);
    void rawcreate(const BaseList<size_t>&, const BaseList<size_t>&);
    void count_active();
  };

 template<typename T> struct type_traits< BlockedMatrix<T> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef Matrix<T> value_type;
 };

  template<typename T1,typename T2> void duplicate_structure(BlockedMatrix<T1>& d, const T2& a) { d.duplicate_structure(a); }
  

 template<class T1,class T2> void full(BlockedMatrix<T1>& dest, const ListList<T2>& source)
 {
   dest.duplicate_structure(source);
   for (size_t i=source.size();i--;) {
     if (source.size(i))
       full(dest(i),source(i));
   }
 }

 template<typename T> size_t simple_counter_needed(const BlockedMatrix<T>& a) { return simple_counter::needed(a.raw_storage()); }

// template<class T1,class T2> void full(BlockedMatrix<T1>& d, const BaseList<T2>& a, const ListList<size_t>& blocking)
// {
//   d.duplicate_structure(blocking);
//   d=T1(0);
//   for (size_t blk=d.size();blk--;) {
//     cmatrix dcur(d(blk));
//     const BaseList<size_t> curind(blocking(blk));
//     for (size_t k=curind.size();k--;)
//       dcur(k,k)=a(curind(k));
//   }
// }

//   template<class T1,class T2> void add_ip(BlockedMatrix<T1>& d, const BaseList<T2>& a, const ListList<size_t>& blocking)
//   {
//     if (!d) {
//       full(d,a,blocking);
//       return;
//     }
//     size_t blk=d.size();
//     if (blk!=blocking.size())
//       throw Mismatch("BlockedMatrix::add_ip");
//     for (;blk--;) {
//       cmatrix dcur(d(blk));
//       const BaseList<size_t> curind(blocking(blk));
//       for (size_t k=curind.size();k--;)
// 	dcur(k,k)+=a(curind(k));
//     }
//   }
    
 inline double norm(const BlockedMatrix<complex>& a) { return norm(a.row()); }

  //Specialised 'NMR' functions    
 template<class T> void propagator(BlockedMatrix<complex>&, const BlockedMatrix<T>&, double dt);
  template<class T> void propagator_ns(BlockedMatrix<complex>&, const BlockedMatrix<T>&, double dt, const ListList<double>&);
  template<class T> void propagator_second(BlockedMatrix<complex>&, const BlockedMatrix<T>&, double dt, const BaseList< ListList<size_t> >&, const ListList<double>&);
  template<class T1,class T2> void propagator_ns_Hrf(BlockedMatrix<complex>&, const BlockedMatrix<T1>&, const BlockedMatrix<T2>&, double dt, const ListList<double>&);
  template<class T1,class T2> void propagator_second_Hrf(BlockedMatrix<complex>&, const BlockedMatrix<T1>&, const BlockedMatrix<T2>&, double dt, const BaseList< ListList<size_t> >&, const ListList<double>&);
 void propagator(ListList<complex>&, const ListList<double>&, double dt);
 template<class T> void unitary_simtrans(BlockedMatrix<complex>&, const BlockedMatrix<T>&, const BlockedMatrix<complex>&);
 template<class T> void unitary_isimtrans(BlockedMatrix<complex>&, const BlockedMatrix<T>&, const BlockedMatrix<complex>&);
 
 template<class T> void transpose(BlockedMatrix<T>&, const BlockedMatrix<T>&);
 template<class T> void conj_transpose(BlockedMatrix<T>&, const BlockedMatrix<T>&);
 template<class T> T trace(const BlockedMatrix<T>&);
  //double hermitian_trace(const BlockedMatrix<complex>&);
  // double hermitian_trace_multiply(const BlockedMatrix<complex>&, const BlockedMatrix<complex>&);
  // template<class T> T trace_multiply(const BlockedMatrix<T>&, const BlockedMatrix<T>&);
  //complex trace_multiply_conj_transpose(const BlockedMatrix<complex>&, const BlockedMatrix<complex>&);
 
//  inline complex trace_conj_transpose_multiply(const BlockedMatrix<complex>& a, const BlockedMatrix<complex>& b)
//    { return conj(trace_multiply_conj_transpose(a,b)); }

//   template<class T1,class T2> inline LCM_NEWOBJECT(BlockedMatrix<T1>,T2) operator* (const T2& b,const BlockedMatrix<T1>& a)
//  { 
//    LCM_NEWOBJECT(BlockedMatrix<T1>,T2) d;
//    multiply(d,b,a);
//    return d;
//  }

 template<class T1,class T2,class T3> void mla(BlockedMatrix<T1>& d, const T2& a, const ListList<T3>& b)
   {
     size_t n=b.size();
     if (!d) {
       ListList<T1> tmp(b);
       tmp*=a;
       full(d,tmp);
     }
     else {
       if (d.size()!=n)
	 throw Mismatch("mla");
       for (;n--;)
	 mla(d(n),a,b(n));
     }
   }

  //implementation details

 template<class T> template<class T2> BlockedMatrix<T>& BlockedMatrix<T>::operator+=(const T2& a) {
   if (!*this)
     set(a);
   else {
     if (this->size()!=a.size())
       throw Mismatch("BlockedMatrix+=");
     if (active()!=this->size())
       throw Failed("BlockedMatrix<T>+=: invalid on defective BlockedMatrix");
     for (size_t i=this->size();i--;)
       (*this)(i)+=a(i);
   }
   return *this;
 }

 template<class T> template<class T2> BlockedMatrix<T>& BlockedMatrix<T>::operator+=(const BlockedMatrix<T2>& a) {
   if (!*this)
     set(a);
   else {
     if (!arematching(*this,a))
       throw Mismatch("BlockedMatrix+=");
     store_+=a.row();
   }
   return *this;
 }

 template<class T> template<class T2> BlockedMatrix<T>& BlockedMatrix<T>::operator&=(const T2& a) {
   if (!*this)
     set(a);
   else {
     if (this->size()!=a.size())
       throw Mismatch("BlockedMatrix&=");
     for (size_t i=this->size();i--;) {
       Matrix<T>& curm((*this)(i));
       if (!curm) {
	 if (!a(i).empty())
	   throw Mismatch("BlockedMatrix&=");
       }
       else
	 curm&=a(i);
     }
   }
   return *this;
 }

  template<size_t N> template<typename T,typename T2> void lcm_post_multiply_<N>::apply(BlockedMatrix<T>& a, const T2& b)
  {
    if (!a)
      a.set(b);
    else {
      if (a.size()!=b.size())
	throw Mismatch("BlockedMatrix*=");
      for (size_t i=a.size();i--;) {
	Matrix<T>& curm(a(i));
	if (!curm) {
	  if (!b(i).empty())
	    throw Mismatch("BlockedMatrix*=");
	}
	else
	  curm*=b(i);
      }
    }
  }

 template<> void BlockedMatrix<complex>::unitary_simtrans(const BlockedMatrix<complex>&);
 template<> void BlockedMatrix<complex>::unitary_isimtrans(const BlockedMatrix<complex>&);
 template<> void BlockedMatrix<complex>::unitary_simtrans(const ListList<complex>&);
 template<> void BlockedMatrix<complex>::unitary_isimtrans(const ListList<complex>&);

 template<class T> void transpose(BlockedMatrix<T>& a, const BlockedMatrix<T>& b)
   {
     a.duplicate_structure(b,true);
     for (size_t i=b.size();i--;) {
       if (!!b(i))
	 transpose(a(i),b(i));
     }
   }

 template<class T> void conj_transpose(BlockedMatrix<T>& a, const BlockedMatrix<T>& b)
   {
     a.duplicate_structure(b,true);
     for (size_t i=b.size();i--;) {
       if (!!b(i))
	 conj_transpose(a(i),b(i));
     }
   }

//utility functions for transforming propagators etc.
BlockedMatrix<complex> rotatez(const BlockedMatrix<complex>&, const ListList<double>& Fz, double angle);
void rotatez_ip(BlockedMatrix<complex>&, const ListList<double>& Fz, double angle);
 inline void rotatez_ip(BlockedMatrix<complex>& U, const ListList<complex>& zrot) { U.unitary_simtrans(zrot); }
void rotatezfacs(ListList<complex>&, const ListList<double>& Fz, double angle);

 //Not a full structure check: ensures that data sizes are compatible and no. of active blocks is same
 template<class T1,class T2> bool arematching(const BlockedMatrix<T1>& a, const BlockedMatrix<T2>& b)
 {
   return ((a.row().size()==b.row().size()) && (a.active()==b.active()));
 }

 template<class T> template<class T2> BlockedMatrix<T>& BlockedMatrix<T>::operator-=(const BlockedMatrix<T2>& a) {
   if (!*this) {
     duplicate_structure(a);
     ::libcmatrix::negate(store_,a.row());
   }
   else {
     if (!arematching(*this,a))
       throw Mismatch("BlockedMatrix-=");
     store_-=a.row();
   }
   return *this;
 }

 template<class T> BlockedMatrix<T>::BlockedMatrix(const Matrix<T>& a, mxflag::tempflag tflag)
   : store_(a.size(),a.vector(),tflag)
   //NB create via raw pointer - otherwise data always copied
   {
     DynamicList< Matrix<T> >::create(1,Matrix<T>(a.rows(),a.cols(),store_.vector(),mxflag::nondynamic));
     count_active();
   }

  template<class T> BlockedMatrix<T>::BlockedMatrix(const BaseList<size_t>& str) : DynamicList< Matrix<T> >()
    { create(str); }

  template<class T> BlockedMatrix<T>::BlockedMatrix(const BaseList<size_t>& rstr, const BaseList<size_t>& cstr) : DynamicList< Matrix<T> >()
    { create(rstr,cstr); }

    template<class T> BlockedMatrix<T>::BlockedMatrix(const BaseList<size_t>& str, const T& val) : DynamicList< Matrix<T> >()
      { create(str,val); }

  template<class T> BlockedMatrix<T>::BlockedMatrix(const BaseList<size_t>& rstr, const BaseList<size_t>& cstr, const T& val) : DynamicList< Matrix<T> >()
  { create(rstr,cstr,val); }
  
  template<class T> BlockedMatrix<T>::BlockedMatrix(const BlockedMatrix<T>& a)
    : DynamicList< Matrix<T> >(), store_(a.store_)
  { buildmatrices(a); }

  template<class T> template<class T2> BlockedMatrix<T>::BlockedMatrix(const BlockedMatrix<T2>& a)
    : DynamicList< Matrix<T> >(), store_(a.raw_storage(),mxflag::normal)
  { buildmatrices(a); }

    template<class T> template<class T2> void BlockedMatrix<T>::duplicate_structure(const ListList<T2>& a)
      {
	size_t n=a.size();
	ScratchList<size_t> sizes(n);
	for (;n--;)
	  sizes(n)=a.size(n);
	create(sizes);
      }

    template<class T> template<class T2> void BlockedMatrix<T>::duplicate_structure(const BaseList< Matrix<T2> >& a)
      {
	size_t n=a.size();
	ScratchList<size_t> rsizes(n);
	ScratchList<size_t> csizes(n);
	for (;n--;) {
	  rsizes(n)=a(n).rows();
	  csizes(n)=a(n).cols();
	}
	create(rsizes,csizes);
      }

   template<class T> template<class T2> size_t BlockedMatrix<T>::getsize(const BaseList< Matrix<T2> >& source) 
     {
       size_t tot=0;
       for (size_t k=source.size();k--;)
	 tot+=source(k).size();
       return tot;
     }

   template<class T> BlockedMatrix<T>& BlockedMatrix<T>::operator= (const BlockedMatrix<T>& a)
  {
    if (this!=&a) {
      store_=a.store_;
      buildmatrices(a);
    }
    return *this;
  }

   template<class T> template<class T2> BlockedMatrix<T>& BlockedMatrix<T>::operator= (const BlockedMatrix<T2>& a)
  {
    store_=a.row();
    buildmatrices(a);
    return *this;
  }

   inline void real(BlockedMatrix<double>& d, const BlockedMatrix<complex>& a) {
     d.duplicate_structure(a);
     BaseList<double> drow(d.row());
     real(drow,a.row());
   }

   inline BlockedMatrix<double> real(const BlockedMatrix<complex>& a) {
     BlockedMatrix<double> d;
     real(d,a);
     return d;
   }

   template<class T> template<class T2> void BlockedMatrix<T>::buildmatrices(const BlockedMatrix<T2>& a, bool transpose) 
  {
    T* addr=store_.vector();
    active_=a.active();
    const size_t n=a.size();
    if (this->size()!=n)
      DynamicList< Matrix<T> >::create(n,Matrix<T>(mxflag::nondynamic));
    size_t count=0;
    for (size_t i=0;i<n;i++) {
      size_t r,c;
      if (transpose) {
	r=a(i).cols();
	c=a(i).rows();
      }
      else {
	r=a(i).rows();
	c=a(i).cols();
      }
      const size_t items=r*c;
      Matrix<T> tmp(r,c,items ? addr : NULL,mxflag::nondynamic);
      (*this)(i).swap(tmp);
      if (items) {
	addr+=items;
	count++;
      }
    }
    assert(count==active());
  }

  template<class T> bool BlockedMatrix<T>::ismatching(const BaseList<size_t>& str) const
  {
    if (str.size()!=this->size())
      return false;
    for (size_t j=this->size();j--;) {
      const size_t n=str(j);
      if ((rows(j)!=n) || (cols(j)!=n))
	return false;
    }
    return true;
  }
      
  template<class T> bool BlockedMatrix<T>::ismatching(const BaseList<size_t>& rstr, const BaseList<size_t>& cstr) const
  {
    const size_t n=rstr.size();
    if (n!=cstr.size())
      throw Mismatch("BlockedMatrix<T>::ismatching");

    if (n!=this->size())
      return false;
    for (size_t j=n;j--;) {
      if ((rows(j)!=rstr(j)) || (cols(j)!=cstr(j)))
	return false;
    }
    return true;
  }

  template<class T> void BlockedMatrix<T>::create(const BaseList<size_t>& str, const T& val) 
    { 
      if (!ismatching(str))
	rawcreate(str,str);
      store_=val;
    }

  template<class T> void BlockedMatrix<T>::create(const BaseList<size_t>& rstr, const BaseList<size_t>& cstr, const T& val) 
    { 
      if (!ismatching(rstr,cstr))
	rawcreate(rstr,cstr);
      store_=val;
    }

  template<class T> void BlockedMatrix<T>::create(const BaseList<size_t>& str) 
    { 
      if (!ismatching(str))
	rawcreate(str,str);
    }

  template<class T> void BlockedMatrix<T>::create(const BaseList<size_t>& rstr, const BaseList<size_t>& cstr) 
    { 
      if (!ismatching(rstr,cstr))
	rawcreate(rstr,cstr);
    }

  template<class T> void BlockedMatrix<T>::rawcreate(const BaseList<size_t>& rstr, const BaseList<size_t>& cstr)
    {
      const size_t n=rstr.size();      
      if (n!=cstr.size()) 
       throw Mismatch("BlockedMatrix<T>::create"); 

      if (n!=this->size()) {
	DynamicList< Matrix<T> >::clear();
	DynamicList< Matrix<T> >::create(n,Matrix<T>(mxflag::nondynamic));
      }
      
      size_t i;
      size_t totalitems=0;
      for (i=n;i--;)
	totalitems+=rstr(i)*cstr(i);
      store_.create(totalitems);
      T* addr=store_.vector();
      active_=0;
      for (i=0;i<n;i++) {
	const size_t items(rstr(i)*cstr(i));
	Matrix<T> tmp(rstr(i),cstr(i),items ? addr : NULL,mxflag::nondynamic);
	(*this)(i).swap(tmp);
	if (items) {
	  active_++;
	  addr+=items;
	}
      }
    }
     
  template<class T> void BlockedMatrix<T>::swap(Matrix<T>& a) {
    if (this->size()>1)
      throw Failed("BlockedMatrix::Matrix swap invalid");
    store_.swap(a.store_);
    Matrix<T> tmp(a.rows(),a.cols(),store_.vector(),mxflag::nondynamic);
    if (this->empty()) {
      DynamicList< Matrix<T> >::create(1,tmp);
      a.clear();
    }
    else {
      std::swap(a.r_,this->front().r_);
      std::swap(a.c_,this->front().c_); //highly dodgy
      this->front().swap(tmp);
    }
    count_active();
  }

  template<class T> void BlockedMatrix<T>::count_active() {
    active_=0;
    for (size_t i=this->size();i--;) {
      if (!!(*this)(i))
	active_++;
    }
  }

  template<class T> std::ostream& operator<< (std::ostream& ostr, const BaseList< Matrix<T> >& a) {
    if (a.empty())
      return ostr << "<empty>\n";
    for (size_t i=0;i<a.size();i++)
      ostr << a(i) << '\n';
    return ostr;
  }

  template<class T> void spy(std::ostream& ostr, const BlockedMatrix<T>& a)
  {
    bool done=false;
    for (size_t i=0;i<a.size();i++) {
      if (!!a(i)) {
	if (a.size()>1)
	  ostr << "Block " << i << '\n';
	spy(ostr,a(i));
	done=true;
      }
    }
    if (!done)
      ostr << "<empty>\n";
  }

  template<class T> std::ostream& BlockedMatrix<T>::dump_content(std::ostream& ostr) const
    {
      return ostr << static_cast< const BaseList< Matrix<T> >& >(*this);
    }

  template<class T> std::ostream& operator<< (std::ostream& ostr, const BlockedMatrix<T>& a) 
    {
      if (a.empty())
	return ostr << "<undefined>\n";
      
      ostr << "Active blocks: " << a.active() << '/' << a.size() << '\n';
      return ostr << static_cast< const BaseList< Matrix<T> >& >(a);
      //      return a.dump_content(ostr);
    }

  template<class T> template<class T2> BlockedMatrix<T>& BlockedMatrix<T>::operator= (const BaseList< Matrix<T2> >& a)
  {
    duplicate_structure(a);
    for (size_t i=a.size();i--;) {
      Matrix<T> curm((*this)(i),mxflag::nondynamic);
      curm=a(i);
    }
    return *this;
  }

  template<class T> T trace(const BlockedMatrix<T>& a)
    {
      if (!a)
	throw Undefined("trace");
      T sum(0);
      for (size_t i=a.size();i--;) {
	if (!!a(i))
	  sum+=trace(a(i));
      }
      return sum;
    }

  template<typename T> void BlockedMatrix<T>::identity()
  {
    for (size_t i=this->size();i--;)
      (*this)(i).identity();
  }

  template<typename T> void pow(BlockedMatrix<T>&, const BlockedMatrix<T>&, int); //!< only instantiated for complex?

  template<class T> void hermitian_eigensystem(BlockedMatrix<T>&, ListList<double>&, const BlockedMatrix<T>&);

  template<class Td,class T,class T2> void unitary_simtrans(BlockedMatrix<Td>& d, const T2& a, const BlockedMatrix<T>& V)
   {
     if (!arematching(a,V))
       throw Mismatch("unitary_simtrans");
     d.duplicate_structure(V);
     for (size_t k=d.size();k--;) {
       Matrix<Td> dk(d(k),mxflag::nondynamic);
       if (!dk.empty())
	 unitary_simtrans(d(k),a(k),V(k));
     }
   }

  template<class T> void conj_transpose_multiply(BlockedMatrix<complex>& d, const BlockedMatrix<complex>& a, const BlockedMatrix<T>& b)
  {
    if (!arematching(a,b))
      throw Mismatch("conj_transpose_multiply");
    d.duplicate_structure(b);
    for (size_t k=d.size();k--;) {
      Matrix<complex> dk(d(k),mxflag::nondynamic);
      if (!dk.empty())
	conj_transpose_multiply(dk,a(k),b(k));
    }
  }

  template<class T> void multiply_conj_transpose(BlockedMatrix<complex>& d, const BlockedMatrix<complex>& a, const BlockedMatrix<T>& b)
  {
    if (!arematching(a,b))
      throw Mismatch("multiply_conj_transpose");
    d.duplicate_structure(b);
    for (size_t k=d.size();k--;) {
      Matrix<complex> dk(d(k),mxflag::nondynamic);
      if (!dk.empty())
	multiply_conj_transpose(dk,a(k),b(k));
    }
  }

#ifdef LCM_USE_EXTERNTEMPLATE
  //!< don't currently instantiate double version as not all ops defined
  //LCM_L3_EXTERN template class BlockedMatrix<double>;
LCM_L3_EXTERN template class BlockedMatrix<complex>;
#endif

} //namespace libcmatrix

#endif
