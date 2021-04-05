#undef LCM_SUPPRESS_VIEWS
#include "matlabio.h"
#include "lcm_MultiHolder.h"
#include "cmatrix_utils.h"

namespace libcmatrix {

  struct Matlab5Header {
    char text[124];
    UINT16_t version;
    char M;
    char I;
  };
  
  namespace {
    const bool ambigend(ambigendian());

    struct doessingle : public ::std::unary_function<double,float> {
      float operator()(double v) const { return v; }
    };
    
    int pad(int nbytes)
    {
      const int left= nbytes & 7;
      return left ? nbytes+8-left : nbytes;
    }
    
    char errorbuf[256]; //!< scratch space for errors

    struct Matlab4Tag {
      Matlab4Tag(FILE*);
      Fmatrix data;
      bool isbigend;
      bool issingle;
    }; 

  }

  void pad(FILE* fp, int nbytes)
  {
    const int left=nbytes & 7;
    if (left)
      ::std::fseek(fp,8-left,SEEK_CUR);
  }

  template<int N> void Matlab5Tag::endianswap(bool isbigend)
  {
    if (isbigend!=ambigend) {
      char* ptr(data.get());
      ::libcmatrix::endianswap(ptr);
    }
  }
  
  template void Matlab5Tag::endianswap<2>(bool);
  template void Matlab5Tag::endianswap<4>(bool);
  template void Matlab5Tag::endianswap<8>(bool);

  matlab_controller::matlab_controller(const char* fname)
    : padmode(LCM_DEFAULT_MATLAB_PADMODE), isread(true), endpos(0)
  {
    fp=open_matlab(fname,"rb");
    if (!fp)
      throw Failed("matlab_controller: failed to open file");
    version=get_version(fp);
    if (version==5) {
      Matlab5Header hdr;
      if (fread(&hdr,sizeof(Matlab5Header),1,fp)!=1)
	throw Failed("matlab_controller: failed to read MATLAB header");
      isbigend=(hdr.M=='M');
    }
  }

int matlab_controller::get_version(FILE *fp)
{
  int version=5;
  file_position_lock markobj(fp); //will reset file pointer on exit
  for (int i=4;i--;) { // if any of the first 4 bytes are zero, file is V4
    if (getc(fp)==0) {
      version=4;
      break;
    }
  }
  return version;
}

  void ReadMATLAB4Header(size_t& rows,size_t& cols,bool& iscomplex,bool& needswap,matlab_element_t& dtype,FILE* fp)
{
  Matlab4Tag tag(fp);

  cols=tag.data.ncols;
  rows=tag.data.mrows;
  iscomplex=tag.data.imagf;
  needswap=(tag.isbigend != ambigend);
  dtype=(tag.issingle ? miSINGLE : miDOUBLE);
}
    
void Matlab5Tag::read_all(FILE* fp,bool isbigend)
{
  const size_t oldbytes=nbytes;
  read_header(fp,isbigend);

  int nread;
  //  if ((nbytes==0) && !iscompressed)
  //  std::cerr << "Warning: zero size but uncompressed Matlab tag at file offset " << ftell(fp) << '\n';
  if (!iscompressed) {
    if (nbytes<=8) {
      data.reset(chars,mxflag::nondynamic);
      if (nbytes==0)
	return;
    }
    else {
      if (nbytes>oldbytes) //don't claim memory if have enough
	data.reset(new char[nbytes]);
    }
    nread=fread(data.get(),nbytes,1,fp);
    pad(fp,nbytes);
  }
  else {
    data.reset(chars+4,mxflag::nondynamic);
    nread=fread(data.get(),4,1,fp);
  }
  if (nread!=1)
    throw FileCorrupt("Matlab5Tag::read_all");
}

inline int getint(char* chars,bool needswap)
{
  int32_t& asint=*reinterpret_cast<int32_t*>(chars);
  if (needswap) 
    endianswap(asint);
  return asint;
}

inline UINT16_t getshort(char* chars,bool needswap)
{
  UINT16_t& asshort=*reinterpret_cast<UINT16_t*>(chars);
  if (needswap) 
    endianswap(asshort);
  return asshort;
}

void Matlab5Tag::read_header(FILE* fp,bool isbigend)
{
  if (fread(chars,4,1,fp)!=1)
    throw FileCorrupt();

  const int whichbytes= isbigend ? 0 : 2;
  iscompressed= (chars[whichbytes] || chars[whichbytes+1]);
  const bool needswap=(isbigend!=ambigend);

  if (iscompressed) {
    if (isbigend) {
      nbytes=getshort(chars,needswap);
      dtype=static_cast<matlab_element_t>(getshort(chars+2,needswap));
    }
    else {
      nbytes=getshort(chars+2,needswap);
      dtype=static_cast<matlab_element_t>(getshort(chars,needswap));
    }
    return;
  }
  // read rest of tag
  if (fread(chars+4,4,1,fp)!=1)
    throw FileCorrupt(); 
  dtype=static_cast<matlab_element_t>(getint(chars,needswap));
  nbytes=getint(chars+4,needswap);
}
  
static const UINT32_t complex_mask=8;

template<class Iter,typename T> static int raw_read_(Iter,int,FILE*, bool,int)
{ throw Failed("Matlab::raw_read_: can't read into this data type"); }

struct base_reader_ {
  FILE* fp;
  bool needswap;
  base_reader_(FILE* fp_, bool needswap_) : fp(fp_), needswap(needswap_) {}
};

template<typename SourceType> struct reader_ : public base_reader_ {
  reader_(FILE* fp_,bool needswap_) : base_reader_(fp_,needswap_) {}

  template<typename T> void get(T& dest) {
    const int nread=fread(&tmp,sizeof(tmp),1,fp);
    assert(nread==1);
    dest=tmp;
  }
 
  SourceType tmp;
};
  
template<> struct reader_<float> : public base_reader_ {
  reader_<float>(FILE* fp_,bool needswap_) : base_reader_(fp_,needswap_) {}

  template<typename T> void get(T& dest) {
    get(tmp);
    dest=tmp;
  }
  
  void get(float& x) {
    const int nread=fread(&x,sizeof(x),1,fp);
    assert(nread==1);
    if (needswap)
      endianswap(x);
  }

  void get(bool& dest) {
    get(tmp);
    dest = tmp ? true : false;
  }

  float tmp;
};

template<> struct reader_<double> : public base_reader_ {
  reader_<double>(FILE* fp_,bool needswap_) : base_reader_(fp_,needswap_) {}

  template<typename T> void get(T& dest) {
    get(tmp);
    dest=tmp;
  }
  
  void get(double& x) {
    const int nread=fread(&x,sizeof(x),1,fp);
    assert(nread==1);
    if (needswap)
      endianswap(x);
  }

  void get(bool& dest) {
    get(tmp);
    dest = tmp ? true : false;
  }

  double tmp;
};
  
  template<class Iter,class T> void raw_read_(Iter dest, int n, FILE* fp, bool needswap, Type2Type<T>, bool iscompressed)
{
  reader_<T> reader(fp,needswap);
  for (size_t i=n;i--;++dest)
    reader.get(*dest);
  int nbytes=n*sizeof(T);
  if (iscompressed)
    nbytes+=4;
  pad(fp,nbytes);
}

  template<class Iter> void raw_read(const Iter& dest,int n,FILE* fp,bool needswap, matlab_element_t dtype, bool iscompressed)
{
  switch (dtype) {
  case miSINGLE:
    raw_read_(dest,n,fp,needswap,Type2Type<float>(),iscompressed);    
    break;

  case miDOUBLE:
    raw_read_(dest,n,fp,needswap,Type2Type<double>(),iscompressed);
    break;
    
  case miUINT8:
    raw_read_(dest,n,fp,needswap,Type2Type<UINT8_t>(),iscompressed);
    break;
    
  case miUINT16:
    raw_read_(dest,n,fp,needswap,Type2Type<UINT16_t>(),iscompressed);
    break;
    
  case miINT16:
    raw_read_(dest,n,fp,needswap,Type2Type<int16_t>(),iscompressed);
    break;
    
  case miINT32:
    raw_read_(dest,n,fp,needswap,Type2Type<int32_t>(),iscompressed);
    break;

  case miUINT32:
    raw_read_(dest,n,fp,needswap,Type2Type<UINT32_t>(),iscompressed);
    break;

  case miCOMPRESSED:
    throw Failed("Matlab read: compressed data not supported");

  case miINT64: case miUINT64:
    throw Failed("Matlab read: found 64 bit integers in numeric array (non-standard!)");

  case miUTF8: case miUTF16: case miUTF32:
    throw Failed("Matlab read: found UTF data in numeric array (non-standard!)");

  default:
    LCM_SNPRINTF(errorbuf,sizeof(errorbuf),"Matlab read: unsupported data type: %i",dtype);
    throw Failed(errorbuf);
  }
}

template<class Iter> void matlab_controller::read_array_matlab5(const Iter& start, size_t nitems)
{
  Matlab5Tag datatag;
  datatag.read_header(fp,isbigend);
  raw_read(start,nitems,fp,isbigend!=ambigend,datatag.dtype,datatag.iscompressed);
}

//postfix iterators but these should be pointers anyway
template<class Iterout,class Iterin> void zipcomplex(Iterout dest,int items,Iterin real,Iterin imag)
{
  for (;items--;)
    *dest++=complex(*real++,*imag++);
}

  template<> void matlab_controller::ReadArray5(Matrix<complex>& dest,bool iscomplex)
{
  const size_t nr=dest.rows();
  const size_t nc=dest.cols();

  rmatrix reals(nr,nc);
  if (!reals)
    read_array_matlab5((double*)NULL,reals.size());
  else    
    read_array_matlab5(reals.permuted_begin(),reals.size());
  if (iscomplex) {
    rmatrix imags(nr,nc);
    if (!imags)
      read_array_matlab5((double*)NULL,imags.size());
    else
      read_array_matlab5(imags.permuted_begin(),imags.size());
    zipcomplex(dest.begin(),dest.size(),reals.begin(),imags.begin());
  }
  else
    dest=reals;

  flush();
}
 
  template<size_t N> void matlab_controller::ReadArray5(MultiMatrix<complex,N>& dest,bool iscomplex)
{
  MultiMatrix<double,N> reals;
  reals.set_dimensions(dest);
  read_array_matlab5(reals.permuted_begin(),reals.size());
  if (iscomplex) {
    MultiMatrix<double,N> imags;
    imags.set_dimensions(dest);
    read_array_matlab5(imags.permuted_begin(),imags.size());
    zipcomplex(dest.begin(),dest.size(),reals.begin(),imags.begin());
  }
  else
    dest=reals;

  flush();
}
 
  template<> void matlab_controller::ReadArray5(List<complex>& dest,bool iscomplex)
{
  const size_t n=dest.size();
  ScratchList<double> reals(n);
  read_array_matlab5(reals.begin(),n);
  if (iscomplex) {
    ScratchList<double> imags(n);
    read_array_matlab5(imags.begin(),n);
    zipcomplex(dest.begin(),n,reals.begin(),imags.begin());
  }
  else
    dest=reals;

  flush();
}
 
  template<class T> void matlab_controller::ReadArray5(Matrix<T>& dest,bool iscomplex)
{
  if (iscomplex)
    throw NoComplexRead();
  if (!dest)
    read_array_matlab5((T*)NULL,dest.size());
  else
    read_array_matlab5(dest.permuted_begin(),dest.size());
  flush();
}

  template<> void matlab_controller::ReadArray5(List<double>& dest,bool iscomplex)
{
  if (iscomplex)
    throw NoComplexRead();
  read_array_matlab5(dest.begin(),dest.size());
  flush();
}
 
  template<size_t N> void matlab_controller::ReadArray5(MultiMatrix<double,N>& dest,bool iscomplex)
{
  if (iscomplex)
    throw NoComplexRead();
  read_array_matlab5(dest.permuted_begin(),dest.size());
  flush();
}
  
//   void Matlab5Tag::skip(FILE* fp, bool isbigend)
//   {
//     read_header(fp,isbigend);
//     const int nskip=iscompressed ? 4 : pad(nbytes);
//     if (fseek(fp,nskip,SEEK_CUR))
//       throw FileCorrupt("Matlab5Tag::skip");
//   }

bool matlab_controller::next()
{
  if (version!=5)
    throw Failed("matlab_controller::next: only valid for V5 file");
  ensureread();
  file_position_lock lock(fp);
  const long nextpos=read_size(false);
  if (nextpos==0) {
    lock.unlock(); //!< doesn't really matter whether fp is reset or not
    close();
    return false;
  }
  fseek(nextpos);
  lock.unlock();
  return true;
}

Warning<> matlab_controller::not_a_matrix_warning("ReadMATLAB5: Object is not a matrix",&lcm_io_warning);
Warning<> matlab_controller::padding_warning("ReadMATLAB5: Number of dimensions in read object less than expected",&lcm_io_warning);
Warning<> matlab_controller::compressed_warning("ReadMATLAB5: compressed data encountered (unsupported)",&lcm_io_warning);

long matlab_controller::read_size(bool checkvalid) 
{
  try {
    Matlab5Tag tag;
    tag.read_header(fp,isbigend);
    
    if (checkvalid) {
      switch (tag.dtype) {
      case miCOMPRESSED:
	compressed_warning.raise();
	return 0;
	
      case miMATRIX:
	break;

      default: {
	char buf[64];
	snprintf(buf,sizeof(buf)," (type=%i)",(int)(tag.dtype));
	not_a_matrix_warning.raise(buf);
	return 0;
      }
      }
    }
    if (tag.iscompressed)
      return ftell()+4;
    
    const size_t nbytes=(tag.nbytes+7) & ~7; //pad out to 8 bytes
    return ftell()+nbytes;
  } catch (...) {
    return 0;
  }
}

void matlab_controller::flush()
{
  static Warning<> underrun_warning("matlab_controller warning: read has underrun by more than 7 bytes",&lcm_io_warning);

  const size_t curpos=ftell();
  if (curpos>endpos) {
    LCM_SNPRINTF(errorbuf,sizeof(errorbuf),"matlab_controller: read overrun by %lu bytes",curpos-endpos);
    throw Failed(errorbuf);
  }
  if (endpos-curpos>7)
    underrun_warning.raise();
  if (::std::fseek(fp,endpos,SEEK_SET))
    throw Failed("flush: fseek");
}

void matlab_controller::read_matrix_header(header_info& info)
{
  ensureread();
  Matlab5Tag AF(fp,isbigend); // array flags

  info.dataclass=AF[isbigend ? 3 : 0];
  switch (info.dataclass) {
  case mxCELL_CLASS:
    info.type=CELL;
    break;
  case mxCHAR_CLASS:
    info.type=CHAR;
    break;
  case mxSTRUCT_CLASS:
    info.type=STRUCT;
    break;
  case mxOBJECT_CLASS: case mxSPARSE_CLASS:
    LCM_SNPRINTF(errorbuf,sizeof(errorbuf),"ReadMATLAB5: Unexpected/unhandled class (e.g. object, sparse): %i",info.dataclass);
    throw Failed(errorbuf);
  default:
    info.type=ARRAY;
  }

  const int Flags=int(AF[isbigend ? 2 : 1]);
  //  std::cout << "data[1]: " << AF.data[1]  << "  data[2]: " << AF.data[2]  << "  Flags: " << Flags << "\n";
  info.iscomplex= Flags & complex_mask;
  
  Dims.read_all(fp,isbigend);

  //  size_t ND=Dims.nbytes/4;
  
//   size_t* dimvals=reinterpret_cast<size_t*>(Dims.data.get());BAD
//   if (isbigend != ambigend) {
//     for (size_t i=ND;i--;)
//       endianswap(dimvals[i]);
//   }
  //info.dims.create(ND,dimvals);
  info.dims=Dims.row<int32_t>(isbigend);
}

bool matlab_controller::find(header_info& info, const char* name)
{
  file_position_lock lock(fp); //!< restore original file pos unless match found
  while (peek(info)) {
    if (strcmp(name,info.name)==0) {
      lock.unlock();
      return true;
    }
    if (!next())
      return false;
  }
  return false;
}

bool matlab_controller::peek(header_info& info)
{
  if (version!=5)
    throw Failed("matlab_controller::peek: only valid for V5 file");

  file_position_lock markobj(fp); //will reset file pointer at exit

  if (read_size()==0)
    return false;

  read_matrix_header(info);
  arrayname.read_all(fp,isbigend);
  const size_t n=arrayname.nbytes;
  lastname.create(n+1);
  memcpy(lastname.vector(),arrayname.data.get(),n); //not necessarily null terminated!
  lastname(n)='\0';
  info.name=lastname.vector();

  if ((info.type==CELL) && !!info) {
    header_info tmp;
    if (read_size()==0) 
      throw Failed("matlab_controller::peek");
    read_matrix_header(tmp);
    info.iscomplex=tmp.iscomplex;
  }  
    
  return true;
}

void matlab_controller::ReadMATLAB5Header(BaseList<size_t> dims,matlab_t& type, bool& iscomplex)
{
  endpos=read_size();
  if (endpos==0)
    throw Failed("ReadMATLAB5Header: end of file");

  header_info info;
  read_matrix_header(info);
  iscomplex=info.iscomplex;
  type=info.type;
  size_t ND=info.dims.size();
  const size_t destdims=dims.size();
  size_t* dimvals=info.dims.vector();

  if (ND<destdims) {
    size_t start,end;
    switch (padmode) {
    case PADNONE:
      padding_warning.raiseas(BaseWarning::RaiseException);
    case PADROW:
      start=destdims-ND;
      end=destdims-1;
      break;
    case PADCOLUMN:
      start=0;
      end=ND-1;
      break;
    }
    padding_warning.raise();
    dims=1;
    BaseList<size_t> sel(dims(range(start,end)));
    sel=dimvals;
  }
  else {
    //allow initial single element dimensions to be ignored
    while ( (ND>dims.size()) && (*dimvals==1)) {
      dimvals++;
      ND--;
    }
    
    if (ND!=dims.size()) {
      LCM_SNPRINTF(errorbuf,sizeof(errorbuf),"ReadMATLABHeader: expecting %" LCM_PRI_SIZE_T_MODIFIER "u dimensions, got %" LCM_PRI_SIZE_T_MODIFIER "u",dims.size(),ND);
      throw Failed(errorbuf);
    }

    dims=dimvals;
  }

  if (!next())
    throw Failed("ReadMATLAB5Header: failed to find next tag");
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

Warning<> matlab_controller::wide_characters_warning("ReadMATLAB5: wide characters detected",&lcm_io_warning);

void matlab_controller::ReadMATLAB5(List<char>& dest)
{
  bool iscomplex;
  matlab_t type;
  ScratchList<size_t> dimvals(1);
  ReadMATLAB5Header(dimvals,type,iscomplex);
  if (type!=CHAR)
    throw Failed("ReadMatlab5: expecting character string");
  const size_t nitems(dimvals.front());
  dest.create(nitems+1); //reserve space for extra char (e.g. null terminator)
  dest.pop_back();
  Matlab5Tag datatag;
  datatag.read_header(fp,isbigend);
  switch (datatag.dtype) {
  case miUINT16:
    break;
  case miUTF8: case miUTF16: case miUTF32:
    throw Failed("ReadMATLAB5: UTF character strings are not supported");
  default:
    throw FileCorrupt("ReadMATLAB5: expected UINT16 for char array");
  }
  const bool needswap(isbigend!=ambigend);
  UINT16_t ui16;
  bool warn=false;
  for (size_t i=0;i<nitems;i++) {
    int nread=fread(&ui16,sizeof(UINT16_t),1,fp);
    assert(nread==1);
    if (needswap) 
      endianswap(ui16);
    if (ui16>255) {
      warn=true;
      ui16&=~0xFF;
    }
    dest(i)=ui16;
  }
  if (warn)
    wide_characters_warning.raise();

  flush();
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
   
Matlab4Tag::Matlab4Tag(FILE* fp)
{
  if (fread(&data,sizeof(data),1,fp)!=1)
    throw FileCorrupt();
  
  int type=data.type;
  const int floattype=(type-type%1000);

  switch (floattype) {
  case 0: isbigend=false; break;
  case 1000: isbigend=true; break;
  default:    
    LCM_SNPRINTF(errorbuf,sizeof(errorbuf),"Matlab4Tag: unsupported machine type (%i)",floattype/1000);
    throw Unsupported(errorbuf); //!< unsupported float type
  }

  type-=floattype;

  const int otype=(type-type%100);
  type-=otype;

  const int ptype=(type-type%10);
  type-=ptype;

  switch (ptype) {
  case 0: issingle=false; break;
  case 10: issingle=true; break;
  default:
    LCM_SNPRINTF(errorbuf,sizeof(errorbuf),"Matlab4Tag: unsupported data type (%i)",ptype/10);
    throw Unsupported(errorbuf);
  }

  if (type!=0) {
    LCM_SNPRINTF(errorbuf,sizeof(errorbuf),"Matlab4Tag: unsupported object type (%i)",type);
    throw Unsupported(errorbuf);
  }

  if (isbigend!=::libcmatrix::ambigend) {
    ::libcmatrix::endianswap(data.type);
    ::libcmatrix::endianswap(data.mrows);
    ::libcmatrix::endianswap(data.ncols);
    ::libcmatrix::endianswap(data.imagf);
    ::libcmatrix::endianswap(data.namlen);
  }
  ::std::fseek(fp,data.namlen,SEEK_CUR);
}

void ReadData4(cmatrix& dest,FILE* fp,bool iscomplex,bool needswap,matlab_element_t dtype)
{
  const size_t rows=dest.rows();
  const size_t cols=dest.cols();
  const size_t items=rows*cols;

  rmatrix reals(rows,cols);
  raw_read(reals.permuted_begin(),items,fp,needswap,dtype,false);

  if (iscomplex) {
    rmatrix imags(rows,cols);
    raw_read(imags.permuted_begin(),items,fp,needswap,dtype,false);
    zipcomplex(dest.begin(),items,reals.begin(),imags.begin());
  }
  else 
    dest=reals;
}
  
void ReadData4(BaseList<complex>& dest,FILE* fp,bool iscomplex,bool needswap,matlab_element_t dtype)
{
  const int points=dest.size();

  ScratchList<double> store(iscomplex ? 2*points : points);
  BaseList<double> reals(points,store.vector());
  raw_read(reals.begin(),points,fp,needswap,dtype,false);

  if (iscomplex) {
    BaseList<double> imags(points,store.vector()+points);
    raw_read(imags.begin(),points,fp,needswap,dtype,false);

    for (size_t i=points;i--;)
      dest(i)=complex(reals(i),imags(i));
  }
  else 
    dest=reals;
}
  
void ReadData4(rmatrix& dest,FILE* fp,bool iscomplex,bool needswap,matlab_element_t dtype)
{
  if (iscomplex)
    throw NoComplexRead();
  raw_read(dest.permuted_begin(),dest.size(),fp,needswap,dtype,false);
}
  
void ReadData4(BaseList<double>& dest,FILE* fp,bool iscomplex,bool needswap,matlab_element_t dtype)
{
  if (iscomplex)
    throw NoComplexRead();
  raw_read(dest.begin(),dest.size(),fp,needswap,dtype,false);
}
  
template<class T> void ReadMATLAB4(Matrix<T>& dest, FILE* fp)
{
  size_t rows,cols;
  bool needswap,iscomplex;
  matlab_element_t dtype;
  ReadMATLAB4Header(rows,cols,iscomplex,needswap,dtype,fp);

  if ((rows==0) || (cols==0)) {
    dest.kill();
    return;
  }
  dest.create(rows,cols);
  ReadData4(dest,fp,iscomplex,needswap,dtype);
}

template<> void ReadMATLAB4(List<char>&, FILE*) {
  throw Failed("Read char array from V4 not supported");
}

template<class T> void ReadMATLAB4(List<T>& dest,FILE* fp)
{
  size_t rows,cols;
  bool needswap,iscomplex;
  matlab_element_t dtype;
  ReadMATLAB4Header(rows,cols,iscomplex,needswap,dtype,fp);

  if ((rows==0) || (cols==0)) {
    dest.kill();
    return;
  }
  if ((rows!=1) && (cols!=1))
    throw Failed("Data is not one-dimensional");

  dest.create(rows*cols);
  ReadData4(dest,fp,iscomplex,needswap,dtype);
}

//Instantiate forms required
template void matlab_controller::ReadMATLAB5(MultiHolder<complex>&);
template void matlab_controller::ReadMATLAB5(MultiHolder<double>&);
template void matlab_controller::ReadMATLAB5(MultiMatrix<double,3>&);
template void matlab_controller::ReadMATLAB5(MultiMatrix<complex,3>&);
template void matlab_controller::ReadMATLAB5(MultiMatrix<double,4>&);
template void matlab_controller::ReadMATLAB5(MultiMatrix<complex,4>&);
template void matlab_controller::ReadMATLAB5(Matrix<complex>&);
template void matlab_controller::ReadMATLAB5(Matrix<double>&);
template void matlab_controller::ReadMATLAB5(Matrix<bool>&);
template void matlab_controller::ReadMATLAB5(List<complex>&);
template void matlab_controller::ReadMATLAB5(List<double>&);
template void matlab_controller::ReadMATLAB5(List< Matrix<complex> >&);
template void matlab_controller::ReadMATLAB5(List< Matrix<double> >&);
template void matlab_controller::ReadMATLAB5(List< Matrix<bool> >&);

template void ReadMATLAB4(Matrix<complex>&, FILE*);
template void ReadMATLAB4(Matrix<double>&, FILE*);
template void ReadMATLAB4(List<complex>&, FILE*);
template void ReadMATLAB4(List<double>&, FILE*);

} //namespace libcmatrix
