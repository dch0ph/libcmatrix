#ifndef _BlockedMatrixTensor_h_
#define _BlockedMatrixTensor_h_

#include "Matrix.h"
#include "Tensor.h"

namespace libcmatrix {

  template<class T> class MatrixTensor {
  public:
    MatrixTensor<T>() : size_(0) {}
    MatrixTensor<T>(size_t sizev) 
    { create(sizev); }

    bool operator!() const { return (size_==0); }

    MatrixTensor& operator= (const T& a) {
      switch (size_) {
      case 0:
	throw Undefined("MatrixTensor=");
      case 1:
	element_=a;
	break;
      default:
	tensor_=a;
	break;
      }
      return *this;
    }

    void create(size_t sizev, int rank, mxflag::tensflag flag =mxflag::maximum) {
      size_=sizev;
      switch (sizev) {
      case 0:
	tensor_.clear();
	break;
      case 1:
	tensor_.clear();
	element_.create(rank,flag);
	break;
      default:
	tensor_.create(rank,flag);
	for (int r=0;r<=rank;r++) {
	  if (tensor_.have_rank(r))
	    create_rank(r);
	}
      }
    }
    
    void clear() {
      tensor_.clear();
      size_=0;
    }

    size_t rows() const { return size_; }
    size_t cols() const { return size_; }

    bool iselement() const
    { return (size_==1); }

    bool ismatrix() const
    { return (size_>1); }

    Tensor<T>& element() {
      if (!iselement())
	throw Failed("MatrixTensor::element: not an element");
      return element_;
    }

    const Tensor<T>& element() const {
      if (!iselement())
	throw Failed("MatrixTensor::element: not an element");
      return element_;
    }

    T& element(int l, int m) {
      if (!iselement())
	throw Failed("MatrixTensor::element: not an element");
      return element_(l,m);
    }

    const T& element(int l, int m) const {
      if (!iselement())
	throw Failed("MatrixTensor::element: not an element");
      return element_(l,m);
    }

    Tensor<Matrix<T> >& matrix() {
      if (!ismatrix())
	throw Failed("MatrixTensor::matrix: not an matrix");
      return tensor_;
    }

    const Tensor<Matrix<T> >& matrix() const {
      if (!ismatrix())
	throw Failed("MatrixTensor::matrix: not an matrix");
      return tensor_;
    }

    Matrix<T>& matrix(int l, int m) {
      if (!ismatrix())
	throw Failed("MatrixTensor::matrix: not an matrix");
      return tensor_(l,m);
    }

    const Matrix<T>& matrix(int l, int m) const {
      if (!ismatrix())
	throw Failed("MatrixTensor::matrix: not an matrix");
      return tensor_(l,m);
    }

    bool have_rank(int l) const { //slight cheat - relies on undefined tensor not throwing on have_rank
      return ismatrix() ? tensor_.have_rank(l) : element_.have_rank(l);
    }

    void ensure_rank(int l) {
      switch (size_) {
      case 0:
	break;
      case 1:
	element_.ensure_rank(l);
	break;
      default:
	tensor_.ensure_rank(l);
	create_rank(l);
      }
    }

    void ensure_rank(int l, const T& def) {
      switch (size_) {
      case 0:
	break;
      case 1:
	element_.ensure_rank(l,def);
	break;
      default:
	tensor_.ensure_rank(l);
	create_rank(l);
	tensor_(l)=def;
      }
    }

  private:
    Tensor<Matrix<T> > tensor_;
    Tensor<T> element_;
    size_t size_;

    void create_rank(int r) {
      for (int m=-r;m<=r;m++)
	tensor_(r,m).create(size_,size_);
    }
  };

    template<class T> struct BlockedMatrixTensor : public DynamicList<MatrixTensor<T> > {
      BlockedMatrixTensor(mxflag::tempflag tflag =mxflag::normal)
	: DynamicList<MatrixTensor<T> >(tflag) {}
      BlockedMatrixTensor(const BaseList<size_t>& str, int rank =2, mxflag::tensflag flag =mxflag::maximum) 
	: DynamicList<MatrixTensor<T> >() { create(str, rank, flag); }
      BlockedMatrixTensor(const BaseList<size_t>& str, int rank, const T& val, mxflag::tensflag flag =mxflag::maximum)
	: DynamicList<MatrixTensor<T> >() { create(str,rank,val,flag); }

      //overwrites create present in List
      void create(const BaseList<size_t>&, int =2, mxflag::tensflag =mxflag::maximum);
      void create(const BaseList<size_t>&, int rank, const T&, mxflag::tensflag =mxflag::maximum);

      //      void ensure(const BaseList<size_t>&);

      void ensure_rank(int l) {
	for (size_t i=this->size();i--;)
	  (*this)(i).ensure_rank(l);
      }

      void ensure_rank(int l, const T& def) {
	for (size_t i=this->size();i--;)
	  (*this)(i).ensure_rank(l,def);
      }

      BlockedMatrixTensor<T>& operator= (const T& val) {
	for (size_t i=this->size();i--;) {
	  if (!!(*this)(i))
	    (*this)(i)=val; //need to avoid meaningless assignment to missing block
	}
	return *this;
      }

      //      bool operator!() const { return DynamicList<MatrixTensor<T> >::operator!(); }
      bool ismatching(const BaseList<size_t>&) const;      

    };

  template<class T> std::ostream& operator<< (std::ostream& ostr, const MatrixTensor<T>& a) {
    switch (a.rows()) {
    case 0:
      ostr << "<null>\n";
      break;
    case 1:
      ostr << a.element() << '\n';
      break;
    default:
      ostr << a.matrix();
    }
    return ostr;
  }

 template<typename T> struct type_traits< BlockedMatrixTensor<T> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef MatrixTensor<T> value_type;
 };

  //Implementation details

  template<class T> bool BlockedMatrixTensor<T>::ismatching(const BaseList<size_t>& str) const
  {
    for (size_t i=str.size();i--;) {
      if (str(i)!=(*this)(i).rows())
	return false;
    }
    return true;
  }
      
  template<class T> void BlockedMatrixTensor<T>::create(const BaseList<size_t>& str, int rank, const T& val, mxflag::tensflag flag) 
  { 
    create(str,rank,flag);
    *this=val;
  }

  template<class T> void BlockedMatrixTensor<T>::create(const BaseList<size_t>& str, int rank, mxflag::tensflag flag) 
     { 
       const size_t n=str.size();      
       DynamicList<MatrixTensor<T> >::create(n);
       for (size_t i=n;i--;) //empty unused matrices
	 (*this)(i).create(str(i),rank,flag);
     }
   
  template<class T> std::ostream& operator<< (std::ostream& ostr, const BlockedMatrixTensor<T>& a) {
    for (size_t i=0;i<a.size();i++) {
      switch (a(i).rows()) {
      case 0:
	ostr << "<null>\n";
	break;
      case 1:
	ostr << a(i).element() << '\n';
	break;
      default:
	ostr << a(i).matrix() << '\n';
      }
    }
    return ostr;
  }
    
} //namespace libcmatrix

#endif
