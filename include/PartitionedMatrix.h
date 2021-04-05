#ifndef lcm_PartitionedMatrix_h_
#define lcm_PartitionedMatrix_h_

/*! \file
  \brief Partitioning of matrices into non-zero blocks
  
  (Not very complete; sufficient for Chebyshev propagation) */

#include "cmatrix.h"
#include "ScratchList.h"
#include "IndirectList.h"
#include "BlockedMatrix.h"
#include "Warnings.h"

namespace libcmatrix {
  
  class diagonal_partition {
  public:
    diagonal_partition() {}
    explicit diagonal_partition(const BaseList<size_t>&);
    double density() const; //!< return non-zero fraction (based on *blocks* not *elements*)
    size_t row_size() const { return offs_.size()-1; }
    size_t rows() const { return offs_.back(); }
    size_t rows(size_t i) const {
#if (BOUNDS_CHECK)
      if (i>=row_size())
	throw BadIndex("diagonal_partition::rows",i,row_size());
#endif
      return offs_(i+1)-offs_(i); 
    }

    range index_range(size_t i) const {
#if (BOUNDS_CHECK)
      if (i>=row_size())
	throw BadIndex("diagonal_partition::index_range",i,row_size());
#endif
      return range(offs_(i),offs_(i+1)-1);
    }

    static void makepartition(ScratchList<size_t>& offs, const BaseList<size_t>& sizes);
    template<class T1,class T2,class T3> void multiply(Matrix<T1>&, const Matrix<T2>&, const Matrix<T3>&) const;

    template<class T> static bool isnull(const T& a, double tol =0.0) {
      return tol
	? isabsent(a,checksnonzero_tolerance< LCM_VAL(T) >(tol))
	: isabsent(a,checksnonzero< LCM_VAL(T) >()); 
    }

    template<class T,class F> static bool isabsent(const T& a, const F& func) {
      const typename T::const_iterator end(a.end());
      return (std::find_if(a.begin(),end,func)==end);
    }
 
    template<class T> bool validate(const Matrix<T>& a, double tol) const; //!< return \c true if partitioning is valid for matrix \a a within tolerance \a tol

    friend class matrix_partition;
    friend std::ostream& operator<< (std::ostream& ostr, const diagonal_partition& a);
  private:
    ScratchList<size_t> offs_;

    template<class T, class F> bool validate_(const Matrix<T>& a, const F&) const;
  };

  class matrix_partition;

  class matrix_partition_base {
  public:
    matrix_partition_base() {}
    matrix_partition_base(const BaseList<size_t>& sizes) { create(sizes); }
    matrix_partition_base(const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes) { create(rsizes,csizes); }

    bool operator!() const { return roffs_.empty(); } //!< undefined only if partition structure not created
    size_t rows() const { return roffs_.back(); }
    size_t rows(size_t i) const { return roffs_(i+1)-roffs_(i); }
    size_t cols() const { return coffs_.back(); }
    size_t cols(size_t i) const { return coffs_(i+1)-coffs_(i); }
    size_t row_size() const { return roffs_.size()-1; }
    size_t col_size() const { return roffs_.size()-1; }
//     size_t totalblocks() const {
//       if (!isdiagonal_)
// 	throw Failed("matrix_partition::totalblocks: only valid for diagonal matrices");
//       return totalblocks_;
//     }
//     size_t totalrows() const { 
//       if (!isdiagonal_)
// 	throw Failed("matrix_partition::totalrows: only valid for diagonal matrices");
//       return totalrows_;
//     }

    bool isdiagonal() const { return isdiagonal_; }
    //    bool isantidiagonal() const { return (totalblocks_!=0); }

    range row_range(size_t i) const {
#if (BOUNDS_CHECK)
      if (i>=row_size())
	throw BadIndex("row_range",i,row_size());
#endif
      return range(roffs_(i),roffs_(i+1)-1);
    }
    range column_range(size_t i) const {
#if (BOUNDS_CHECK)
      if (i>=col_size())
	throw BadIndex("column_range",i,col_size());
#endif
      return range(coffs_(i),coffs_(i+1)-1);
    }

    void create(const BaseList<size_t>&);
    void create(const BaseList<size_t>&, const BaseList<size_t>&);

    friend std::ostream& operator<< (std::ostream&, const matrix_partition&);

  private:
    ScratchList<size_t> roffs_;
    ScratchList<size_t> coffs_;
    bool isdiagonal_;

    template<class T> RawMatrix<T> extract(const Matrix<T>& a, size_t r, size_t c) const {
      return RawMatrix<T>(a,row_range(r),column_range(c));
    }
  };

  //! just check matching number of blocks - mismatch of block sizes will be caught further down
  inline bool arematching(const matrix_partition_base& a, const matrix_partition_base& b)
  {
    return (a.row_size()==b.row_size()) && (a.col_size()==b.col_size());// && (a.isantidiagonal()==b.isantidiagonal());
  }


  class matrix_partition : public matrix_partition_base {
  public:
    matrix_partition() : density_(-1.0) {}
    matrix_partition(const BaseList<size_t>& sizes) : matrix_partition_base(sizes) {}
    matrix_partition(const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes) : matrix_partition_base(rsizes,csizes) {}
    template<class T> matrix_partition(const BaseList<size_t>&, const BaseList<size_t>&, const Matrix<T>&);
    matrix_partition(const Matrix<bool>& present, const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes);
    matrix_partition(const Matrix<bool>& present, const BaseList<size_t>& sizes);//, size_t totalblocksv =0);

    void adddiagonal();
    template<class T> bool validate(const Matrix<T>& a, double tol) const; //!< return \c true if partitioning is valid for matrix \a a within tolerance \a tol
    matrix_partition& operator|=(const matrix_partition&);

    template<class T> void findblocks(const Matrix<T>&);
    void create(const BaseList<size_t>&, const BaseList<size_t>&);
    void create(const BaseList<size_t>&);//, size_t =0);
    matrix_partition& operator= (const matrix_partition_base&);

    friend std::ostream& operator<< (std::ostream& ostr, const matrix_partition& a);
    
    //    template<class T1,class T2,class T3> void partitioned_multiply(Matrix<T1>&, const Matrix<T2>&, const Matrix<T3>&) const;

    const Matrix<bool>& present() const { return present_; }
    bool present(size_t j, size_t k) const { return present_(j,k); }
    bool& present(size_t j, size_t k) { return present_(j,k); }

    double density() const; //!< return non-zero fraction (based on *blocks* not *elements*)
    void swappresent(matrix_partition&);
    void zero() { present_=false; density_=0.0; }
    void invalidate() const { density_=-1.0; }
    
  private:
    Matrix<bool> present_;
    mutable double density_;
    
    template<class T> static T front(const Matrix<T>& a) { return a(0U,0U); }
    template<class T> static T front(const RawMatrix<T>& a) { return *(a.vector()); }

    //    template<class T1,class T2,class T3> void partitioned_multiply_(Matrix<T1>& d, const T2& a, const Matrix<T3>& b) const;
    template<class T, class F> bool validate_(const Matrix<T>& a, const F&) const;
      
    void create();
    void cleanempty();
  };

  inline bool arematching(const matrix_partition& a, const matrix_partition& b)
  {
    return (a.row_size()==b.row_size()) && (a.col_size()==b.col_size());
  }

  std::ostream& operator<< (std::ostream&, const diagonal_partition&);
  std::ostream& operator<< (std::ostream&, const matrix_partition&);
  
  struct diagonal_partition_set;

  struct matrix_partition_set;
  //! needs to be defined for |=
  template<> struct type_traits<matrix_partition_set> {
    static const bool trivialconstructor=false;
    static const size_t dimensionality=1;
    static const size_t rank=0;
    typedef matrix_partition value_type;
  };

  struct matrix_partition_set : public List<matrix_partition> {
    matrix_partition_set() {}
    matrix_partition_set& operator|= (const matrix_partition_set& a) { return or_ip(*this,a); }
    double density() const;
    void adddiagonal();
  };

  std::ostream& operator<< (std::ostream& ostr, const matrix_partition_set&);

  struct diagonal_partition_set : public List<diagonal_partition> {
    template<class T1,class T2,class T3> void multiply(BlockedMatrix<T1>&, const BlockedMatrix<T2>&, const BlockedMatrix<T3>&) const;
    double density() const;
  };

//! type to store generalised blocked matrix
/** Full matrix is stored to facilitate linear operations
    present/absent matrix keeps multiplication efficient */

  template<class T> class PartitionedMatrix : protected Matrix<T> {
  public:
    PartitionedMatrix() {}
    PartitionedMatrix(const Matrix<T>&, const matrix_partition&);
    PartitionedMatrix(const Matrix<T>&, const matrix_partition&, mxflag::tempflag);

//   PartitionedMatrix(const Matrix<bool>& present, const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes, const T& v) { 
//     create(present,rsizes,csizes); 
//     row()=v;
//   }
       
    const matrix_partition& partitioning() const { return partition_; }
    matrix_partition& partitioning() { return partition_; }
    using Matrix<T>::row;
    using Matrix<T>::rows;
    using Matrix<T>::cols;
    using Matrix<T>::empty;
    using Matrix<T>::operator!;

    //    size_t rows() const { return partition_.rows(); }
    //size_t cols() const { return partition_.rows(); }
    //size_t size() const { return store_.size(); }

//     size_t dimension(size_t n) const {
//       switch (n) {
//       case 0: return rows();
//       case 1: return cols();
//       default: throw BadIndex("PartitionedMatrix::dimension");
//       }
//     }

//    void clear() { store_.clear(); }
    
    PartitionedMatrix<T>& operator-= (const PartitionedMatrix<T>& a) {
      Matrix<T>::operator-= (static_cast<const Matrix<T>& >(a));
      partition_|=a.partition_; //!< N.B. don't try to spot cancellation
      return *this;
    }

    PartitionedMatrix<T>& operator+= (const PartitionedMatrix<T>& a) {
      Matrix<T>::operator+=(static_cast<const Matrix<T>& >(a));
      partition_|=a.partition_;
      return *this;
    }

    Matrix<T>& asmatrix() { return static_cast< Matrix<T>& >(*this); }
    const Matrix<T>& asmatrix() const { return static_cast<const Matrix<T>& >(*this); }

    void create(const matrix_partition_base& a) {
      partition_=a;
      Matrix<T>::create(a.rows(),a.cols());
    }
    
    template<class T2> PartitionedMatrix& operator*=(const T2& v) {
      LCM_STATIC_CHECK( LCM_DIM(T2)==0, PartitionedMatrix_in_place_multiplication_must_be_scalar ); 
      Matrix<T>::operator*= (v);
      return *this;
    }
         
    void print(std::ostream&) const;
    
    void zero() { 
      Matrix<T>::operator=(T(0));
      partition_.zero();
    }
       
    void swap(PartitionedMatrix<T>& a) {
      if (arematching(partition_,a.partitioning()))
	partition_.swappresent(a.partition_);
      else
	std::swap(partition_,a.partitioning());
      //	throw Mismatch("swap (PartitionedMatrix): incompatible partitioning");
      Matrix<T>::swap(a);
    }

    RawMatrix<T> block(size_t r, size_t c) { return RawMatrix<T>(*this,partition_.row_range(r),partition_.column_range(c)); }
    const RawMatrix<T> block(size_t r, size_t c) const { return RawMatrix<T>(*this,partition_.row_range(r),partition_.column_range(c)); }
    
//   void create(const Matrix<bool>& present, const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes);

  private:
    matrix_partition partition_;
  };

  extern Warning<> chebyshev_convergence_warning;
  extern Warning<> chebyshev_nonunitary_warning;

  void chebyshev_propagator(Matrix<complex>&, const Matrix<double>&, double, const matrix_partition* =NULL);
  void chebyshev_propagator(Matrix<complex>&, const Matrix<complex>&, double, const matrix_partition* =NULL);
  
  void propagator(Matrix<complex>&, const Matrix<double>&, double, const diagonal_partition*);
  void propagator(Matrix<complex>&, const Matrix<complex>&, double, const diagonal_partition*);

  void propagator(BlockedMatrix<complex>&, const BlockedMatrix<double>&, double, const diagonal_partition_set&);
  void propagator(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, double, const diagonal_partition_set&);

  /* implementation details below */

  template<class T1,class T2,class T3> void diagonal_partition::multiply(Matrix<T1>& d, const Matrix<T2>& a, const Matrix<T3>& b) const
  {
    if (!!d && (issame(d,a) || issame(d,b)))
      throw ArgumentClash("matrix_partition::multiply");
    if (!issquare(a) || !issquare(b))
      throw NotSquare("matrix_partition::multiply");
    if ((rows()!=b.rows()) || (a.cols()!=b.rows()))
      throw Mismatch("matrix_partition::multiply");
    d.create(rows(),b.cols(),T1(0));

    for (size_t r=rows();r--;) {
      const range sel(index_range(r));
      if (sel.size()) {
	RawMatrix<T1> dsub(d,sel);
	::libcmatrix::multiply(dsub,RawMatrix<T2>(a,sel),RawMatrix<T3>(b,sel));
      }
    }    
  }


//   template<class T1,class M2,class T3> void matrix_partition::partitioned_multiply_(Matrix<T1>& d, const M2& a, const Matrix<T3>& b) const
//   {
//     typedef typename M2::value_type T2;
//     if (!!d && issame(d,b))
//       throw ArgumentClash("partitioned_multiply");
//     if (cols()!=b.rows())
//       throw Mismatch("partitioned_multiply");
//     d.create(rows(),b.cols(),T1(0));
    
//     for (size_t r=present_.rows();r--;) {
//       const range drsel(row_range(r));
//       if (drsel.size()==0)
// 	continue;
//       for (size_t c=present_.cols();c--;) {
// 	if (!present_(r,c))
// 	  continue;
// 	const range brsel(row_range(c));

//       if ((drsel.size()==1) && (brsel.size()==1)) { //!< specialisation for single-element
// 	BaseList<T1> drow(d.row(drsel.start()));
// 	const T2 val(matrix_partition::front(extract(a,r,c)));
// 	mla(drow,val,b.row(brsel.start()));
//       }
//       else {
// 	const RawMatrix<T2> suba(extract(a,r,c));
// 	for (size_t k=present_.cols();k--;) {
// 	  const range csel(column_range(k));
// 	  RawMatrix<T1> subd(d,drsel,csel);
// 	  const RawMatrix<T3> subb(b,brsel,csel);  
// 	  ::libcmatrix::multiply(subd,suba,subb,true);
// 	}
//       }
//     }
//   }      
// }

//   template<class T1,class T2,class T3> void matrix_partition::partitioned_multiply(Matrix<T1>& d, const Matrix<T2>& a, const Matrix<T3>& b) const
//   {
//     if (!!d && issame(d,a))
//       throw ArgumentClash("partitioned_multiply");
//    partitioned_multiply_(d,a,b);
//   }


  template<class T> std::ostream& operator<< (std::ostream& ostr, const PartitionedMatrix<T>& a) {
    a.print(ostr);
    return ostr;
  }
  
  template<class T> void PartitionedMatrix<T>::print(std::ostream& ostr) const
  {
    if (empty()) {
      ostr << "<empty>\n";
      return;
    }
    const Matrix<bool>& present(partition_.present());
    spy(ostr,present);
    for (size_t r=0;r<present.rows();r++) {
      for (size_t c=0;c<present.cols();c++) {
	if (present(r,c)) {
	  ostr << "Block " << r << ',' << c << ":\n";
	  ::libcmatrix::print(block(r,c),ostr) << '\n';
	}
     }
   }
 }


  template<class T> matrix_partition::matrix_partition(const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes, const Matrix<T>& a)
  {
    create(rsizes,csizes);
    findblocks(a);
  }

//   namespace {
//     void intercept() 
//     {
//       std::cout << "intercept\n";
//     }
//   }

  template<class T, class F> bool matrix_partition::validate_(const Matrix<T>& a, const F& func) const
  {
    if ((a.rows()!=rows()) || (a.cols()!=cols()))
      throw Mismatch("matrix_partition::validate");
    for (size_t r=present_.rows();r--;) {
      if (rows(r)==0) //!< trap empty row
	continue;
      for (size_t c=present_.cols();c--;) {
	if (!present_(r,c) && cols(c)) {
	  const RawMatrix<T> asub(a,row_range(r),column_range(c));
	  if (!diagonal_partition::isabsent(asub,func))
	    return false;
	}
      }
    }
    return true;	
 }

  template<class T, class F> bool diagonal_partition::validate_(const Matrix<T>& a, const F& func) const
  {
    if (!issquare(a))
      throw NotSquare("diagonal_partition::validate");    
    if ((a.rows()!=rows()))
      throw Mismatch("diagonal_partition::validate");
    for (size_t r=rows();r--;) {
      if (!isabsent(RawMatrix<T>(a,index_range(r)),func))
	return false;	  
    }
    return true;	
 }

  template<class T> bool matrix_partition::validate(const Matrix<T>& a, double tol) const
  {
    return tol
      ? validate_(a,checksnonzero_tolerance<T>(tol))
      : validate_(a,checksnonzero<T>()); 
  }

  template<class T> bool diagonal_partition::validate(const Matrix<T>& a, double tol) const
  {
    return tol
      ? validate_(a,checksnonzero_tolerance<T>(tol))
      : validate_(a,checksnonzero<T>()); 
  }

  template<class T> void matrix_partition::findblocks(const Matrix<T>& a)
  {
   if (!a)
     throw Undefined("matrix_partition::find_blocks: input matrix is empty");
   
   const size_t rblks=present_.rows();
   const size_t cblks=present_.cols();

   for (size_t r=0;r<rblks;r++) {
     const range rsel(row_range(r));
     if (rsel.size()==0) //!< empty row range
       continue;
     for (size_t c=0;c<cblks;c++) {    
       const range csel(column_range(c));
       if (csel.size() && !diagonal_partition::isnull(a(rsel,csel)))
	 present_(r,c)=true;
     }
   }   
  }
  
  template<class T1,class T2,class T3> void diagonal_partition_set::multiply(BlockedMatrix<T1>& d, const BlockedMatrix<T2>& a, const BlockedMatrix<T3>& b) const
  {
    size_t n=size();
    if ((n!=a.size()) || (n!=b.size()))
      throw Mismatch("diagonal_partition_set::multiply");
    d.duplicate_structure(a);
    for (;n--;) {
      Matrix<complex> dn(d(n),mxflag::nondynamic);
      (*this)(n).multiply(dn,a(n),b(n));
    }
  }

  template<class T> PartitionedMatrix<T>::PartitionedMatrix(const Matrix<T>& a, const matrix_partition& part, mxflag::tempflag flag)
    : Matrix<T>(a,flag),
      partition_(part)
  {
    if ((a.rows()!=partition_.rows()) || (a.cols()!=partition_.cols()))
      throw Mismatch("PartitionedMatrix",a.rows(),a.cols(),part.rows(),part.cols());
  }

  template<class T> PartitionedMatrix<T>::PartitionedMatrix(const Matrix<T>& a, const matrix_partition& part)
    : partition_(part)
  {
    if ((a.rows()>partition_.rows()) && (a.cols()>partition_.cols()))
      Matrix<T>::operator= (RawMatrix<T>(a,range(0,part.rows()-1),range(0,part.cols()-1)));
    else {
      if ((a.rows()!=partition_.rows()) || (a.cols()!=partition_.cols()))
	throw Mismatch("PartitionedMatrix",a.rows(),a.cols(),part.rows(),part.cols());
      Matrix<T>::operator= (a);
    }
  }


  template<typename T> struct type_traits< PartitionedMatrix<T> > {
    static const int rank=0;
    static const bool trivialconstructor=false;
    static const size_t dimensionality=2;
    typedef T value_type;
  };

  template<class T1,class T2,class T3> void raw_multiply(PartitionedMatrix<T1>& d, const PartitionedMatrix<T2>& a, const PartitionedMatrix<T3>& b, bool acc, bool isherm) 
  {
    const matrix_partition& apart(a.partitioning());    
    const matrix_partition& bpart(b.partitioning());    
    //! this will catch mismatching antidiagonal
    if (!arematching(apart,bpart))
      throw Mismatch("multiply (PartitionedMatrix)");
    if (!d) {
      d.create(a.partitioning());
      d.asmatrix()=T1(0);
    }
    else {
      if (!arematching(d.partitioning(),bpart))
	throw Mismatch("multiply (PartitionedMatrix)");
      if (!acc)
	d.zero();
    }
    if (a.cols()!=b.rows())
      throw Mismatch("multiply (PartitionedMatrix)",a.cols(),b.rows());
    if (isherm && !issquare(a))
      throw NotSquare("multiply (PartitionedMatrix)");
//     if (d.isantidiagonal() && !isherm)
//       throw Failed("multiply (PartitionedMatrix): antisymmetric input must also be hermitian");
    
    matrix_partition& dpart(d.partitioning());
    //    const size_t totalblocks=dpart.isantidiagonal() ? dpart.totalblocks() : 0;
    //    const size_t maxrows= dpart.isantidiagonal() ? dpart.totalblocks() : dpart.row_size();
  
    for (size_t j=dpart.row_size();j--;) {
      const range jrows(bpart.row_range(j));
      if (jrows.size()==0)
	continue;
      const size_t end = isherm ? j+1 : dpart.col_size();
//       if (totalblocks) {
// 	size_t newend=totalblocks-end+1;
// 	if (newend<end)
// 	  newend=end;
//       }
      for (size_t k=end;k--;) {
	const range kcols(bpart.column_range(k));
	if (kcols.size()==0)
	  continue;
	bool found=false;
	RawMatrix<T1> djk(d.asmatrix(),jrows,kcols);
	for (size_t m=apart.col_size();m--;) {
	  if (apart.present(j,m) && bpart.present(m,k)) {
	    const range mcols(bpart.column_range(m));
	    const range mrows(bpart.row_range(m));	      
	    multiply(djk,RawMatrix<T2>(a.asmatrix(),jrows,mcols),RawMatrix<T3>(b.asmatrix(),mrows,kcols),true);
	    found=true;
	  }
	}
	if (found) {
	  dpart.present(j,k)=true;
	  if (isherm && (j!=k)) {
	    RawMatrix<T1> dkj(d.asmatrix(),kcols,jrows);
	    conj_transpose(dkj,djk);
	    dpart.present(k,j)=true;
	  }	  
	}
      }
    }
    d.partitioning().invalidate();
  }

  template<class T1,class T2,class T3> void multiply(PartitionedMatrix<T1>& d, const PartitionedMatrix<T2>& a, const PartitionedMatrix<T3>& b) { raw_multiply(d,a,b,false,false); }

  template<class T1,class T2,class T3> inline void hermitian_multiply(T1& d, const T2& a, const T3& b, bool =false) { multiply(d,a,b); } //!< by default, don't do anything different

  //  template<class T1,class T2,class T3> inline void hermitian_multiply(PartitionedMatrix<T1>& d, const PartitionedMatrix<T2>& a, PartitionedMatrix<T3>& b) { raw_multiply(d,a,b,false,true); }

} //namespace libcmatrix

#endif

//   List< Matrix<T> > tmp(1); //create one matrix for workspace
//   const size_t rblks=partition_.present_.rows();
//   const size_t cblks=partition_.present_.cols();

//   //  size_t active=0;
//   for (size_t r=0;r<rblks;r++) {
//     const range rsel(partition_.row_range(r));
//     if (rsel.size()==0) //!< empty row range
//       continue;
//     for (size_t c=0;c<cblks;c++) {    
//       const range csel(partition_.column_range(c));
//       if (csel.size()==0) //!< empty column range
// 	continue;
//       Matrix<T>& tmpm(tmp.back());
//       tmpm=a(rsel,csel); //!< extract sub matrix
//       if (!matrix_partition::isnull(tmpm.row())) {
// 	partition_.present_(r,c)=tmp.size()-1; //!< set offset
// 	//	active++;
// 	tmp.push_back(Matrix<T>()); //!< create new workspace
//       }
//     }
//   }
//   tmp.pop_back(); //!< remove workspace
//   if (tmp.size()) 
//     store_=tmp; //!< copy data
//   else
//     clear();
// }
	
// template<class T> void PartitionedMatrix<T>::create(const Matrix<bool>& present, const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes)
// {
//   const size_t active=partition_.create(present,rsizes,csizes);
//   ScratchList<size_t> rs(active);
//   ScratchList<size_t> cs(active);
//   const size_t rblks=rsizes.size();
//   const size_t cblks=csizes.size();
//   size_t count=0;
//   for (size_t r=0;r<rblks;r++) {
//     for (size_t c=0;c<cblks;c++) {
//       int offset_=partition_.present_(r,c);
//       if (offset_>=0) {
// 	rs(count)=rsizes(r);
// 	cs(count)=csizes(c);
// 	count++;
//       }
//     }
//   }
//   assert(active==count);
//   store_.create(rs,cs);
// }
