/*! \file
  \brief Declares templated functions for reading Matlab V5 files */

#include "matlabio.h"
#include "lcm_MultiHolder.h"

namespace libcmatrix {

  template<size_t N> void matlab_controller::ReadArray5(MultiMatrix<double,N>& dest,bool iscomplex)
  {
    if (iscomplex)
      throw NoComplexRead();
    read_array_matlab5(dest.permuted_begin(),dest.size());
    flush();
  }

  template<class T> void matlab_controller::ReadMATLAB5(Matrix<T>& dest)
  {
    bool iscomplex;
    matlab_t type;
    ScratchList<size_t> dimvals(2);
    ReadMATLAB5Header(dimvals,type,iscomplex);
    if (type!=ARRAY)
      throw Failed("ReadMatlab5: expecting matrix");
    dest.create(dimvals(0),dimvals(1));
    ReadArray5(dest,iscomplex);
  }

  template<class T> void matlab_controller::ReadMATLAB5(List<T>& dest)
  {
    bool iscomplex;
    matlab_t type;
    ScratchList<size_t> dimvals(1);
    ReadMATLAB5Header(dimvals,type,iscomplex);
    if (type!=ARRAY)
      throw Failed("ReadMatlab5: expecting matrix");
    dest.create(dimvals.front());
    ReadArray5(dest,iscomplex);
  }

  template<class T,size_t N> void matlab_controller::ReadMATLAB5(MultiMatrix<T,N>& dest)
{
  bool iscomplex;
  matlab_t type;
  ScratchList<size_t> dimvals(N);
  ReadMATLAB5Header(dimvals,type,iscomplex);
  if (type!=ARRAY)
    throw Failed("ReadMatlab5: expecting matrix");
  dest.dimensions(dimvals);
  ReadArray5(dest,iscomplex);
}

  template<class T> void matlab_controller::ReadMATLAB5(List< Matrix<T> >& dest)
{
  const long initpos=ftell();
  bool iscomplex;
  matlab_t type;
  ScratchList<size_t> dimvals(1);
  ReadMATLAB5Header(dimvals,type,iscomplex);  
  const long end=endpos;
  int dimension;
  switch (type) {
  case ARRAY:
    dimension=1;
    fseek(initpos);
    break;
  case CELL:
    dimension=dimvals.front();
    break;
  default:
    throw Failed("ReadMATLAB5: can't read (char) array into List< Matrix<T> >");
  }
  dest.create(dimension);
  for (size_t i=0;i<dimension;i++) {
    Matrix<T>& curm(dest(i));
    ReadMATLAB5(curm);
  }
  endpos=end; //restore original end check
  flush();
}

  template<class T> void matlab_controller::ReadMATLAB5(MultiHolder<T>& dest)
{
  List<size_t> dims;
  bool iscomplex;
  matlab_t type;
  ReadMATLAB5Header(dims,type,iscomplex);
  if (type!=ARRAY)
    throw Failed("ReadMatlab5: expecting matrix");
  switch (dims.size()) {
    //case 0: return ReadArray5(dest.create(),iscomplex,isbigend,fp);
  case 1: ReadArray5(dest.create(dims(0)),iscomplex); break;
  case 2: ReadArray5(dest.create(dims(0),dims(1)),iscomplex); break;
  case 3: ReadArray5(dest.create(dims(0),dims(1),dims(2)),iscomplex); break;
  case 4: ReadArray5(dest.create(dims(0),dims(1),dims(2),dims(3)),iscomplex); break;
  }
  throw Failed("ReadMATLAB: Invalid dimensionality for MultiHolder<T>");
}
  
}
