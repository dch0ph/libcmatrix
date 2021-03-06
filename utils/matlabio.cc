#include "cmatrix_utils.h"
#include "lcm_writematlab.cc"

namespace libcmatrix {

  struct doessingle : public ::std::unary_function<double,float> {
   float operator()(double v) const { return v; }
 };

  static const bool ambigend(ambigendian());

  void matlab_controller::ensurewrite()
  {
    ensureopen();
    if (isread)
      throw Failed("matlab_controller: file is open for reading");
  }

  void matlab_controller::close()
  {
    if (fp) {
      fclose(fp);
      fp=NULL;
    }
  }

  void matlab_controller::ensureread() {
    ensureopen();
    if (!isread)
      throw Failed("matlab_controller: file is open for writing");
  }

FILE* open_matlab(const char *fname,const char *mode)
{
  if (fname==NULL)
    throw InvalidParameter("open_matlab: null filename");
  ScratchList<char> tmp_(5+strlen(fname));
  char* tmp(tmp_.vector());
  sprintf(tmp,"%s.mat",fname);
  return file_open(tmp,mode);
}

  //  const int32_t matlab_controller_maxstructlength

const char MATLAB_BADWRITE[]="WriteMATLAB: write failed";

static void write_header(FILE *fp,int rows,int cols,const char *name,bool iscomplex)
{
  Fmatrix mat;

  mat.type=ambigend ? 1000 : 0000;
  mat.mrows=rows;
  mat.ncols=cols;
  mat.imagf=iscomplex ? 1 : 0;
  mat.namlen=1+strlen(name);

  if (fwrite(&mat,sizeof(Fmatrix),1,fp)!=1) 
    throw Failed(MATLAB_BADWRITE);
  fwrite(name,sizeof(char),mat.namlen,fp);
}

//! 9/2/2016 unnecessary specialisations removed 

// const char MATLAB_NOTBOOL4[]="bool write to Matlab 4 format not supported - use V5 format";

// template<> void WriteMATLAB4(FILE*, const BaseList<bool>&, const char*)
// {
//   throw Failed(MATLAB_NOTBOOL4);
// }
// template<> void WriteMATLAB4(FILE*, const Matrix<bool>&, const char*)
// {
//   throw Failed(MATLAB_NOTBOOL4);
// }

template<> void WriteMATLAB4(FILE* fp,const Matrix<complex>& a,const char* name)
{
  if (!a) {
    write_header(fp,0,0,name,false);
    return;
  }
  
  write_header(fp,a.rows(),a.cols(),name,true);

  size_t count=0;
  size_t i,j;
  double tmp;
  
  for (i=0;i<a.cols();i++) {
    for (j=0;j<a.rows();j++) {
      tmp=real(a(j,i));
      count+=fwrite(&tmp,sizeof(double),1,fp);
    }
  }
  for (i=0;i<a.cols();i++) {
    for (j=0;j<a.rows();j++) {
      tmp=imag(a(j,i));
      count+=fwrite(&tmp,sizeof(double),1,fp);
    }
  }
  if (count!=2*a.size()) 
    throw Failed(MATLAB_BADWRITE);
}

template<> void WriteMATLAB4(FILE* fp,const BaseList<complex>& a,const char* name)
{
  const size_t points=a.size();
  if (points==0) {
    write_header(fp,0,0,name,false);
    return;
  }
  
  write_header(fp,1,points,name,true);

  size_t count=0;
  size_t i;
  double tmp;
  
  for (i=0;i<points;i++) {
    tmp=real(a(i));
    count+=fwrite(&tmp,sizeof(double),1,fp);
  }
  for (i=0;i<points;i++) {
    tmp=imag(a(i));
    count+=fwrite(&tmp,sizeof(double),1,fp);
  }
  if (count!=2*a.size()) 
    throw Failed(MATLAB_BADWRITE);
}

template<> void WriteMATLAB4(FILE* fp,const BaseList<double>& a,const char* name)
{
  const size_t points=a.size();
  if (points==0) {
    write_header(fp,0,0,name,false);
    return;
  }
  write_header(fp,1,points,name,false);

  if (fwrite(a.vector(),sizeof(double),points,fp)!=points)
    throw Failed(MATLAB_BADWRITE);
}

template<> void WriteMATLAB4(FILE* fp,const Matrix<double>& a,const char* name)
{
  if (!a) {
    write_header(fp,0,0,name,false);
    return;
  }

  write_header(fp,a.rows(),a.cols(),name,false);

  size_t count=0;
  
  for (size_t i=0;i<a.cols();i++) {
    for (size_t j=0;j<a.rows();j++)
      count+=fwrite(&(a(j,i)),sizeof(double),1,fp);
  }
  
  if (count!=a.size()) 
    throw Failed(MATLAB_BADWRITE);
}

void matlab_pad(FILE* fp, int bytes)
{
  const char dummy[]="\0\0\0\0\0\0\0\0";
  const int left=bytes & 7;
  if (left) 
    fwrite(dummy,1,8-left,fp);
}

UINT32_t WriteMATLAB5Tag_raw(FILE* fp,UINT32_t type,UINT32_t bytes, bool trycompress)
{
  if (trycompress && (bytes>0) && (bytes<5)) {
    int16_t bytes16[2];
    if (ambigend) {
      bytes16[0]=bytes;
      bytes16[1]=type;
    }
    else {
      bytes16[1]=bytes;
      bytes16[0]=type;
    }
    fwrite(bytes16,2,2,fp);
    return 4U;
  }
  //    WriteMATLAB5Tag_raw(fp,type,bytes);
  fwrite(&type,4,1,fp);
  fwrite(&bytes,4,1,fp);
  return 8U;
}

  void WriteMATLAB5Tag(FILE* fp,UINT32_t type,UINT32_t bytes,const void* data, bool trycompress)
{
  const UINT32_t totbytes=bytes+WriteMATLAB5Tag_raw(fp,type,bytes,trycompress);
  if (bytes)
    fwrite(data,1,bytes,fp);
  matlab_pad(fp,totbytes);
}

static const UINT32_t complex_mask=8;

void WriteMATLAB5Header(FILE* fp,const char* title)
{
  //Write 128 byte header
  const int headsize=strlen(title);
  if (headsize>=124)
    throw InvalidParameter("WriteMATLAB5Header: header comment too long");
  fputs(title,fp);
  const char empty[128]={' '};
  fwrite(empty,1,124-headsize,fp);
  UINT16_t versend[2];
  versend[0]=0x0100;
  versend[1]=('M' << 8) | 'I';
  fwrite(&versend,2,2,fp);
}

namespace {
  void seek_(FILE* fp, long pos)
  {
    if (fseek(fp,pos,SEEK_SET)) {
      perror("TidyMATLAB5Header seek set");
      fflush(stderr);
      throw Failed("TidyMATLABHeader: seek set failed");
    }
  }

  void writedims_(FILE* fp, const BaseList<size_t>& dims)
  {
    const size_t ndims=dims.size();
    const ScratchList<int32_t> dimsnew(dims); //copy into new array in case size_t is not 32 bit
    WriteMATLAB5Tag(fp,miINT32,4*ndims,dimsnew.vector());
  }

}
  
  long matlab_controller::ftell() const 
  {
    ensureopen();
    return ::std::ftell(fp);
  }
  
void matlab_controller::fseek(long pos) const
{
  ensureopen();
  seek_(fp,pos);
}

  long WriteMATLAB5MatrixHeader(FILE* fp,UINT32_t sourcetype,bool iscomplex,const BaseList<size_t>& dims,const char *name)
  { 
    long initpos=ftell(fp); //remember start
    WriteMATLAB5Tag(fp,miMATRIX,0,NULL);
    
    UINT32_t flagdata[2];
    flagdata[0]=(iscomplex ? complex_mask<<8 : 0) | sourcetype;
    flagdata[1]=0;
    WriteMATLAB5Tag(fp,miUINT32,8,&flagdata);

    writedims_(fp,dims);
    
    WriteMATLAB5Tag(fp,miINT8,strlen(name),name,true);  //write array name
    
    return initpos;
  }

//   void matlab_controller::fend() const
//   {
//     if (::std::fseek(fp,0,SEEK_END))
//       throw Failed("matlab_controller: seek end failed");
//   }

void TidyMATLAB5Header(FILE* fp,long initpos)
{
  const long curpos=ftell(fp);
  const int written=curpos-initpos-8;
  seek_(fp,initpos+4);
  fwrite(&written,4,1,fp); //set number of bytes written
  if (fseek(fp,curpos,SEEK_SET))
    throw Failed("matlab_controller: fp restore failed");
}

// void TidyMATLAB5Header(FILE* fp, long initpos, const BaseList<size_t>& dims)
// {
//   const int written=ftell(fp)-initpos-8;
//   seek_(fp,initpos+4);
//   fwrite(&written,4,1,fp); //set number of bytes written
//   seek_(fp,initpos+24); //!< we know that header has not been compressed - we wrote it!
//   writedims_(fp,dims);
//   if (fseek(fp,0,SEEK_END))
//     throw Failed("TidyMATLABHeader: seek end failed");  
// }

template<class Iter> void WriteMATLAB5Iter_single(FILE* fp,Iter start,size_t nitems)
{
  WriteMATLAB5Iter(fp,compose_unary_iterator(doessingle(),start),nitems);
}

  template<> void WriteMATLAB5_single(FILE* fp,const BaseList<double>& data,const char* name)
  {
    const List<float> singledata(data);
    WriteMATLAB5(fp,singledata,name);
  }

template<> void WriteMATLAB5(FILE* fp,const BaseList<complex>& data,const char* name)
{
  const long initpos=WriteMATLAB5MatrixHeader(fp,matlab_traits<double>::dclass,true,ScratchList<size_t>(size_t(1),data.size()),name);
  const BaseList<complex>::const_iterator start(data.begin());
  WriteMATLAB5Iter(fp,compose_unary_iterator(doesreal<>(),start),data.size());
  WriteMATLAB5Iter(fp,compose_unary_iterator(doesimag<>(),start),data.size());
  TidyMATLAB5Header(fp,initpos);
}

template<> void WriteMATLAB5_single(FILE* fp,const BaseList<complex>& data,const char* name)
{
  const long initpos=WriteMATLAB5MatrixHeader(fp,matlab_traits<float>::dclass,true,ScratchList<size_t>(size_t(1),data.size()),name);
  const BaseList<complex>::const_iterator start(data.begin());
  WriteMATLAB5Iter_single(fp,compose_unary_iterator(doesreal<>(),start),data.size());
  WriteMATLAB5Iter_single(fp,compose_unary_iterator(doesimag<>(),start),data.size());
  TidyMATLAB5Header(fp,initpos);
}

template<> void WriteMATLAB5_single(FILE* fp,const Matrix<double>& data,const char* name)
{
  Matrix<float> datasingle(data);
  WriteMATLAB5(fp,datasingle,name);
}

  template<class T> void WriteMATLAB5_single_cells_(FILE* fp, const T& data, const char* name)
  {
    size_t dim(data.size());
    if (dim==0)
      throw Undefined("WriteMATLAB5_single: List<>");
    //Matlab doesn't like a single dimension...
    const long initpos=WriteMATLAB5MatrixHeader(fp,mxCELL_CLASS,false,ScratchList<size_t>(size_t(1),dim),name);
    for (size_t i=0;i<dim;i++)
      WriteMATLAB5_single(fp,data(i),"");
    TidyMATLAB5Header(fp,initpos);
  }
  

 template<> struct doesreal<float,complex> : public  ::std::unary_function<complex,float> {
   float operator()(const complex& z) const { return real(z); }
 };

 template<> struct doesimag<float,complex> : public  ::std::unary_function<complex,float> {
   float operator()(const complex& z) const { return imag(z); }
 };

template<> void WriteMATLAB5(FILE* fp,const Matrix<complex>& data,const char* name)
{
  const long initpos=WriteMATLAB5MatrixHeader(fp,mxDOUBLE_CLASS,true,ScratchList<size_t>(data.rows(),data.cols()),name);
  const size_t items=data.size();
  if (items) {
    const Matrix<complex>::const_permuted_iterator start=data.permuted_begin();
    WriteMATLAB5Iter(fp,compose_unary_iterator(doesreal<>(),start),items);
    WriteMATLAB5Iter(fp,compose_unary_iterator(doesimag<>(),start),items);
  }
  else {
    WriteMATLAB5Iter(fp,(const double*)NULL,items);
    WriteMATLAB5Iter(fp,(const double*)NULL,items);

  }
  TidyMATLAB5Header(fp,initpos);
}

template<> void WriteMATLAB5_single(FILE* fp,const Matrix<complex>& data,const char* name)
{
  const long initpos=WriteMATLAB5MatrixHeader(fp,mxSINGLE_CLASS,true,ScratchList<size_t>(data.rows(),data.cols()),name);
  const size_t items=data.size();
  if (items) {
    const Matrix<complex>::const_permuted_iterator start=data.permuted_begin();
    WriteMATLAB5Iter(fp,compose_unary_iterator(doesreal<float,complex>(),start),items);
    WriteMATLAB5Iter(fp,compose_unary_iterator(doesimag<float,complex>(),start),items);
  }
  else {
    WriteMATLAB5Iter(fp,(const double*)NULL,items);
    WriteMATLAB5Iter(fp,(const double*)NULL,items);

  }
  TidyMATLAB5Header(fp,initpos);
}

  template<class T> void WriteMATLAB5_single(FILE* fp, const BaseList< Matrix<T> >& data, const char* name)
  {
    WriteMATLAB5_single_cells_(fp,data,name);
  }
  
  template<class T> void WriteMATLAB5_single(FILE* fp, const ListList<T>& data, const char* name)
  {
    WriteMATLAB5_single_cells_(fp,data,name);
  }

template<size_t N> void WriteMATLAB5_single(FILE* fp,const MultiMatrix<double,N>& data,const char* name)
{
  const long initpos=WriteMATLAB5MatrixHeader(fp,matlab_traits<float>::dclass,false,data.dimensions(),name);
  WriteMATLAB5Iter_single(fp,data.permuted_begin(),data.size());
  TidyMATLAB5Header(fp,initpos);
}

template<size_t N> void WriteMATLAB5_single(FILE* fp,const MultiMatrix<complex,N>& data,const char* name)
{
  const long initpos=WriteMATLAB5MatrixHeader(fp,mxSINGLE_CLASS,true,data.dimensions(),name);
  const typename MultiMatrix<complex,N>::const_permuted_iterator start(data.permuted_begin());
  WriteMATLAB5Iter_single(fp,compose_unary_iterator(doesreal<>(),start),data.size());
  WriteMATLAB5Iter_single(fp,compose_unary_iterator(doesimag<>(),start),data.size());
  TidyMATLAB5Header(fp,initpos);
}

  template<class T> void WriteMATLAB5_single(FILE* fp,const MultiHolder<T>& a,const char* name)
  {
    switch (a.dimensions()) {
    case 1: WriteMATLAB5(fp,a.list(),name); return;
    case 2: WriteMATLAB5(fp,a.matrix(),name); return;
    case 3: WriteMATLAB5(fp,a.multimatrix3(),name); return;
    case 4: WriteMATLAB5(fp,a.multimatrix4(),name); return;
    }
    throw Failed("WriteMATLAB: MultiHolder dimension not supported");
  }

#define LCM_INSTANTIATE_WRITEMULTI(T)\
  template void WriteMATLAB5(FILE*, const MultiMatrix<T,3>&, const char*); \
  template void WriteMATLAB5(FILE*, const MultiMatrix<T,4>&, const char*);

#define LCM_INSTANTIATE_WRITEMULTI_SINGLE(T)\
  template void WriteMATLAB5_single(FILE*, const MultiMatrix<T,3>&, const char*); \
  template void WriteMATLAB5_single(FILE*, const MultiMatrix<T,4>&, const char*);

#define LCM_INSTANTIATE_WRITE(T)\
  template void WriteMATLAB5(FILE*, const MultiHolder<T>&, const char*);\
  template void WriteMATLAB5(FILE*, const Matrix<T>&, const char*);\
  template void WriteMATLAB5(FILE*, const BaseList<T>&, const char*);\
  template void WriteMATLAB5(FILE*, const ListList<T>&, const char*);\
  template void WriteMATLAB5(FILE*, const BaseList< Matrix<T> >&, const char*);\
  LCM_INSTANTIATE_WRITEMULTI(T)

//  template void WriteMATLAB5_single(FILE*, const Matrix<T>&, const char*);
//  template void WriteMATLAB5_single(FILE*, const BaseList<T>&, const char*);

#define LCM_INSTANTIATE_WRITE_SINGLE(T)\
  template void WriteMATLAB5_single(FILE*, const MultiHolder<T>&, const char*);\
  template void WriteMATLAB5_single(FILE*, const ListList<T>&, const char*);\
  template void WriteMATLAB5_single(FILE*, const BaseList< Matrix<T> >&, const char*);\
  LCM_INSTANTIATE_WRITEMULTI_SINGLE(T)

LCM_INSTANTIATE_WRITE(bool)
  LCM_INSTANTIATE_WRITE(double)
  template void WriteMATLAB5(FILE*, const MultiHolder<complex>&, const char*);
template void WriteMATLAB5(FILE*, const ListList<complex>&, const char*);
template void WriteMATLAB5(FILE*, const ListList<size_t>&, const char*);
template void WriteMATLAB5(FILE*, const BaseList< Matrix<complex> >&, const char*);
LCM_INSTANTIATE_WRITEMULTI(complex)
//LCM_INSTANTIATE_WRITE(complex)
LCM_INSTANTIATE_WRITE_SINGLE(double)
LCM_INSTANTIATE_WRITE_SINGLE(complex)

//template void WriteMATLAB5(FILE*, const Matrix<bool>&, const char*);
//template void WriteMATLAB5(FILE*, const BaseList< Matrix<bool> >&, const char*);

//template void WriteMATLAB4(FILE*, const Matrix<complex>&, const char*);
//template void WriteMATLAB4(FILE*, const Matrix<double>&, const char*);
//template void WriteMATLAB4(FILE*, const Matrix<bool>&, const char*);
//template void WriteMATLAB4(FILE*, const BaseList<complex>&, const char*);
//template void WriteMATLAB4(FILE*, const BaseList<double>&, const char*);

matlab_controller::matlab_controller(const char* fname, int version_, int flags_, const char* header)
  : isread(false), singleprec(flags_ & singleprecision)
  {
    bool notnew=false;
    if (flags_ & append) {
      FILE* lfp=open_matlab(fname,"rb");
      if (lfp) {
	const int actversion=get_version(lfp);
	if (version_) {
	  if (actversion!=version) {
	    fclose(lfp);
	    throw Failed("matlab_controller: version mismatch on append");
	  }
	}
	else
	  version_=actversion;
	notnew=true;
	fclose(lfp);
      }
    }

    const char* opentype = (flags_ & append) ? "ab" : "wb";
    fp=open_matlab(fname,opentype);
    if (!fp)
      throw Failed("matlab_controller: failed to open file");

    if ((version_!=4) && (version_!=5))
      throw InvalidParameter("matlab_controller: version must be 4 or 5");

    version=version_;

    if ((version==5) && !notnew)
      WriteMATLAB5Header(fp,header ? header : "MATLAB 5.0 MAT-file.  Written by libcmatrix");
    else {
      if (header)
	ignoring_comment_warning.raise();
    }
  }

Warning<> matlab_controller::ignoring_comment_warning("matlab_controller: header comment ignored (V4 or appending existing file",&lcm_io_warning);

template<> void WriteMATLAB5Data(FILE* fp,const BaseList<char>& data)
{
  const size_t nitems=data.size();
  const size_t towrite=nitems*sizeof(UINT16_t);
  const UINT32_t tagbytes=WriteMATLAB5Tag_raw(fp,matlab_traits<char>::dtype,towrite,true);
  size_t written=0;
  for (size_t i=0;i<nitems;i++) {
    const UINT16_t ui16(data(i));
    written+=fwrite(&ui16,sizeof(UINT16_t),1,fp);
  }
  if (written!=nitems)
    throw Failed(MATLAB_BADWRITE);
  matlab_pad(fp,tagbytes+towrite);
}

  template<> void WriteMATLAB5(FILE* fp,const BaseList<char>& data,const char* name)
  {
    const long initpos=WriteMATLAB5MatrixHeader(fp,matlab_traits<char>::dclass,false,ScratchList<size_t>(size_t(1),data.size()),name);
    WriteMATLAB5Data(fp,data);
    TidyMATLAB5Header(fp,initpos);
  }

  template<> void WriteMATLAB4(FILE*, const BaseList<char>&, const char*)
  {
    throw Failed("Write char array into V4 not supported");
  }

template<> void WriteMATLAB5(FILE* fp, const BaseList<const char*>& data, const char* name)
{
  WriteMATLAB5cells_(fp,data,name);
}

} //namespace libcmatrix
