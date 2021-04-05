/*! \file
  \brief Declares templated functions for writing to Matlab V5 files */

#include "matlabio.h"
#include "lcm_MultiHolder.h"
#include "lcm_unary_composition.h"

namespace libcmatrix {
  extern const char MATLAB_BADWRITE[];

  void matlab_pad(FILE*, int);

  template<class Iter> void WriteMATLAB5Iter(FILE* fp,Iter start,size_t nitems)
  {
    typedef typename Unconst< typename ::std::iterator_traits<Iter>::value_type >::Result T;
    const UINT32_t tagbytes=WriteMATLAB5Tag_raw(fp,matlab_traits<T>::dtype,nitems*sizeof(T),true);
    for (size_t i=0;i<nitems;i++) {
      const T x=*start;
      if (fwrite(&x,sizeof(T),1,fp)!=1)
	throw Failed(MATLAB_BADWRITE);
      ++start;
    }
    matlab_pad(fp,tagbytes+nitems*sizeof(T));
  }
  
  template<class T> void WriteMATLAB5(FILE* fp,const Matrix<T>& data,const char* name)
  {
    const long initpos=WriteMATLAB5MatrixHeader(fp,matlab_traits<T>::dclass,false,ScratchList<size_t>(data.rows(),data.cols()),name);
    if (!data)
      WriteMATLAB5Iter(fp,data.row().begin(),data.size());
    else
      WriteMATLAB5Iter(fp,data.permuted_begin(),data.size());
    TidyMATLAB5Header(fp,initpos);
  }

  template<> void WriteMATLAB5(FILE*, const Matrix<complex>&, const char*);
  
  template<class T> void WriteMATLAB5cells_(FILE* fp, const T& data, const char* name)
  {
    size_t dim(data.size());
    if (dim==0)
      throw Undefined("WriteMATLAB5: List<>");
    //Matlab doesn't like a single dimension...
    const long initpos=WriteMATLAB5MatrixHeader(fp,mxCELL_CLASS,false,ScratchList<size_t>(size_t(1),dim),name);
    for (size_t i=0;i<dim;i++)
      WriteMATLAB5(fp,data(i),"");
    TidyMATLAB5Header(fp,initpos);
  }

  template<class T> void WriteMATLAB5(FILE* fp, const BaseList< Matrix<T> >& data, const char* name)
  {
    WriteMATLAB5cells_(fp,data,name);
  }
  
  template<class T> void WriteMATLAB5(FILE* fp, const ListList<T>& data, const char* name)
  {
    WriteMATLAB5cells_(fp,data,name);
  }
    
  template<class T,size_t N> void WriteMATLAB5(FILE* fp,const MultiMatrix<T,N>& data,const char* name)
  {
    const long initpos=WriteMATLAB5MatrixHeader(fp,matlab_traits<T>::dclass,false,data.dimensions(),name);
    WriteMATLAB5Iter(fp,data.permuted_begin(),data.size());
    TidyMATLAB5Header(fp,initpos);
  }
  
  template<size_t N> void WriteMATLAB5(FILE* fp,const MultiMatrix<complex,N>& data,const char* name)
  {
    const long initpos=WriteMATLAB5MatrixHeader(fp,mxDOUBLE_CLASS,true,data.dimensions(),name);
    const typename MultiMatrix<complex,N>::const_permuted_iterator start(data.permuted_begin());
    WriteMATLAB5Iter(fp,compose_unary_iterator(doesreal<>(),start),data.size());
    WriteMATLAB5Iter(fp,compose_unary_iterator(doesimag<>(),start),data.size());
    TidyMATLAB5Header(fp,initpos);
  }

  template<class T> void WriteMATLAB5(FILE* fp,const MultiHolder<T>& a,const char* name)
  {
    switch (a.dimensions()) {
    case 1: WriteMATLAB5(fp,a.list(),name); return;
    case 2: WriteMATLAB5(fp,a.matrix(),name); return;
    case 3: WriteMATLAB5(fp,a.multimatrix3(),name); return;
    case 4: WriteMATLAB5(fp,a.multimatrix4(),name); return;
    }
    throw Failed("WriteMATLAB: MultiHolder dimension not supported");
  }
    
}
