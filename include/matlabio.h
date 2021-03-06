#ifndef LCM_MATLABIO_H
#define LCM_MATLABIO_H

#include <cstdio>
#include "cmatrix.h"
#include "rmatrix.h"
#include "ttyio.h"
#include "lcm_binaryutils.h"
#include "ListList.h"
#include "ScratchList.h"
#include "smartptr.h"

#define LCM_DEFAULT_MATLAB_PADMODE PADNONE

namespace libcmatrix {
  
  template<class T> class MultiHolder;
  template<class T,size_t N> class MultiMatrix;

  template<class T> void WriteMATLAB4(FILE*, const Matrix<T>&, const char*) {
    throw Failed("Write of this Matrix<T> into V4 not supported");
  }
  template<> void WriteMATLAB4(FILE*, const Matrix<complex>&, const char*);
  template<> void WriteMATLAB4(FILE*, const Matrix<double>&, const char*);
  //  template<> void WriteMATLAB4(FILE*, const Matrix<bool>&, const char*);  //!< 9/2/2016 unnecessary specialisation removed

  template<class T> void WriteMATLAB5(FILE*, const Matrix<T>&, const char*);
  template<class T> void WriteMATLAB5_single(FILE*, const Matrix<T>&, const char*);
  template<class T> void WriteMATLAB4(FILE*, const BaseList<T>&, const char*) {
    throw Failed("Write of this List<T> into V4 not supported");
  }
  template<class T> void WriteMATLAB5(FILE*, const ListList<T>&, const char*);
  template<class T> void WriteMATLAB5(FILE*, const BaseList<T>&, const char*);
  template<class T> void WriteMATLAB5_single(FILE*, const BaseList<T>&, const char*);
  template<class T> void WriteMATLAB5_single(FILE*, const ListList<T>&, const char*);
  template<class T> void WriteMATLAB5_single(FILE*, const BaseList< Matrix<T> >&, const char*);
  template<> void WriteMATLAB4(FILE*, const BaseList<char>&, const char*);
  template<> void WriteMATLAB5(FILE*, const BaseList<char>&, const char*);
  inline void WriteMATLAB4(FILE* fp, const char* data, const char* name) { WriteMATLAB4(fp,BaseList<char>(strlen(data)+1,const_cast<char*>(data)),name); }
  inline void WriteMATLAB5(FILE* fp, const char* data, const char* name) { WriteMATLAB5(fp,BaseList<char>(strlen(data)+1,const_cast<char*>(data)),name); }
  template<> void WriteMATLAB4(FILE*, const BaseList<complex>&, const char*);
  template<> void WriteMATLAB4(FILE*, const BaseList<double>&, const char*); //!< added 2/2/16
  template<> void WriteMATLAB5(FILE*, const BaseList<complex>&, const char*);
  template<> void WriteMATLAB5(FILE*, const BaseList<const char*>&, const char*); //!< add 3/8/16
  template<class T> void WriteMATLAB5(FILE*, const MultiHolder<T>&, const char*);
  template<class T,size_t N> void WriteMATLAB4(FILE* ,const MultiMatrix<T,N>&,const char*) { throw Failed("Can't write MultiMatrix to V4 format"); }
  template<class T,size_t N> void WriteMATLAB5(FILE* ,const MultiMatrix<T,N>&,const char*);
  template<size_t N> void WriteMATLAB5(FILE* ,const MultiMatrix<complex,N>&,const char*);
  template<size_t N> void WriteMATLAB5_single(FILE*,const MultiMatrix<double,N>&,const char*);
  template<size_t N> void WriteMATLAB5_single(FILE*,const MultiMatrix<complex,N>&,const char*);
  template<class T> void WriteMATLAB4(FILE*, const ListList<T>&, const char*) { throw Failed("Can't write ListList<> to V4 format"); }
  template<class T> void WriteMATLAB4(FILE*, const BaseList< Matrix<T> >&, const char*) { throw Failed("Can't write List<Matrix> to V4 format"); }

  long WriteMATLAB5MatrixHeader(FILE*, UINT32_t,bool,const BaseList<size_t>&, const char*);
  UINT32_t WriteMATLAB5Tag_raw(FILE* fp,UINT32_t type,UINT32_t bytes, bool =false);

// Matlab 5 element types
  enum matlab_element_t {
    miINT8=1,
    miUINT8=2,
    miINT16=3,
    miUINT16=4,
    miINT32=5,
    miUINT32=6,
    miSINGLE=7,
    miDOUBLE=9,
    miINT64=12,
    miUINT64=13,
    miMATRIX=14,
    miCOMPRESSED=15,
    miUTF8=16,
    miUTF16=17,
    miUTF32=18
  };
  
  //Array (class) types - note these take different values!
  enum matlab_class_t {
    mxCELL_CLASS=1,
    mxSTRUCT_CLASS=2,
    mxOBJECT_CLASS=3,
    mxCHAR_CLASS=4,
    mxSPARSE_CLASS=5,
    mxDOUBLE_CLASS=6,
    mxSINGLE_CLASS=7,
    mxINT8_CLASS=8,
    mxUINT8_CLASS=9,
	mxINT16_CLASS=10,
	mxUINT16_CLASS=11,
    mxINT32_CLASS=12,
    mxUINT32_CLASS=13,
    mxINT64_CLASS=14,
    mxUINT64_CLASS=15
  };

  template<class T> struct matlab_allowsingle {
    static const bool result=false;
  };
  template<> struct matlab_allowsingle<double> {
    static const bool result=true;
  };
  template<> struct matlab_allowsingle<complex> {
    static const bool result=true;
  };
  template<typename T> struct matlab_allowsingle< Matrix<T> > { //necessary for BlockedMatrix<T>
    static const bool result=matlab_allowsingle<T>::result;
  };
    					
  template<class T> struct matlab_default {
    static const int version = (LCM_DIM(T)>2) ? 5 : 4;
  };
  template<class T> struct matlab_default< BaseList< Matrix<T> > >  {
    static const int version = 5;
  };
  template<class T> struct matlab_default< List< Matrix<T> > >  {
    static const int version = 5;
  };
  template<class T> struct matlab_traits {};
  template<> struct matlab_traits<double> {
    static const matlab_class_t dclass=mxDOUBLE_CLASS;
    static const matlab_element_t dtype=miDOUBLE;
  };
  template<> struct matlab_traits<float> {
    static const matlab_class_t dclass=mxSINGLE_CLASS;
    static const matlab_element_t dtype=miSINGLE;
  };
  template<> struct matlab_traits<bool> {
    static const matlab_class_t dclass=mxUINT8_CLASS;
    static const matlab_element_t dtype=miUINT8;
  };
  template<> struct matlab_traits<char> {
    static const matlab_class_t dclass=mxCHAR_CLASS;
    static const matlab_element_t dtype=miUINT16;
  };
  template<> struct matlab_traits<int8_t> {
    static const matlab_class_t dclass=mxINT8_CLASS;
    static const matlab_element_t dtype=miINT8;
  };
  template<> struct matlab_traits<int16_t> {
    static const matlab_class_t dclass=mxINT16_CLASS;
    static const matlab_element_t dtype=miINT16;
  };
  template<> struct matlab_traits<int32_t> {
    static const matlab_class_t dclass=mxINT32_CLASS;
    static const matlab_element_t dtype=miINT32;
  };
  template<> struct matlab_traits<UINT8_t> {
    static const matlab_class_t dclass=mxUINT8_CLASS;
    static const matlab_element_t dtype=miUINT8;
  };
  template<> struct matlab_traits<UINT16_t> {
    static const matlab_class_t dclass=mxUINT16_CLASS;
    static const matlab_element_t dtype=miUINT16;
  };
  template<> struct matlab_traits<UINT32_t> {
    static const matlab_class_t dclass=mxUINT32_CLASS;
    static const matlab_element_t dtype=miUINT32;
  };

#ifndef LCM_DISALLOW_64BIT_TYPES
  template<> struct matlab_traits<int64_t> {
    static const matlab_class_t dclass=mxINT64_CLASS;
    static const matlab_element_t dtype=miINT64;
  };
  template<> struct matlab_traits<UINT64_t> {
    static const matlab_class_t dclass=mxUINT64_CLASS;
    static const matlab_element_t dtype=miUINT64;
  };
#endif
// workaround for MacOS using a distinct type for size_t
// This is a bodge since it assumes that size_t is 8 bytes
#ifdef LCM_MACOS_SIZET_BODGE
  template<> struct matlab_traits<size_t> {
    static const matlab_class_t dclass=mxUINT64_CLASS;
    static const matlab_element_t dtype=miUINT64;
  };
#endif

  template<class T> void ReadMATLAB4(Matrix<T>&, FILE*);
  template<class T> void ReadMATLAB4(List<T>&, FILE*);
  template<> void ReadMATLAB4(List<char>&, FILE*);
  template<class T,size_t N> void ReadMATLAB4(MultiMatrix<T,N>&,FILE*) { throw Failed("Can't read Matlab V4 into MultiMatrix<T,N>"); }
  template<class T> void ReadMATLAB4(MultiHolder<T>&, FILE*) { throw Failed("Can't read Matlab V4 into MultiHolder<T>"); }
  template<class T> void ReadMATLAB4(List< Matrix<T> >&, FILE*) { throw Failed("Can't read Matlab V4 into List<Matrix>"); }

  struct Matlab5Tag {
    char chars[8];
    int nbytes;
    matlab_element_t dtype;
    bool iscompressed;
    smartptr<char,false> data;

    Matlab5Tag() : nbytes(0) {} //ensure data length is valid
    Matlab5Tag(FILE* fp, bool isbigend) : nbytes(0) { read_all(fp,isbigend); }
    
    void read_header(FILE*, bool);
    void read_all(FILE*, bool);
    void skip(FILE*, bool);
    char operator[] (size_t i) const { return (data.get())[i]; }    

    template<typename T> const BaseList<T> row() const {
      if ((matlab_traits<T>::dtype!=dtype) || (nbytes % sizeof(T)))
	throw Failed("Matlab5Tag: data not compatible with read type");
      return BaseList<T>((nbytes/sizeof(T)),(T*)(data.get()));
    }
    template<int N> void endianswap(bool);

    template<typename T> const BaseList<T> row(bool isbigend) {
      //      typedef typename matlab_read_traits<T>::value_type value_type;
      endianswap< sizeof(T) >(isbigend);
      return row<T>();
    }

  };

  class matlab_controller;
  const int32_t matlab_controller_maxstructlength=32;

  struct Matlab5WriteStoreBase {
    Matlab5WriteStoreBase(const char*);
    virtual ~Matlab5WriteStoreBase() {}
    virtual void operator()(matlab_controller&) const =0;
    virtual Matlab5WriteStoreBase* clone() const =0;
    
    char name[matlab_controller_maxstructlength];
  };
  

  class matlab_controller {
  public:
    matlab_controller(const char* fname,int version, int flags =0, const char* header =NULL); //!< write constructor
    matlab_controller(const char* fname);//!< read constructor
    ~matlab_controller() { close(); }

    enum matlab_t { ARRAY, CELL, CHAR, STRUCT };
    enum { singleprecision=1, append=2 };
    enum padmode_t { PADNONE, PADROW, PADCOLUMN }; //!< how to handle cases where no. of dimensions < destination

    //    static const int32_t maxstructlength=32;

    struct header_info {
      const char* name;
      ScratchList<size_t,4> dims;
      bool iscomplex;
      matlab_t type;
      int dataclass;
      
      bool operator! () const { return (dims.empty() || (dims.front()==0)); }    
    };

    void ensurewrite();
    void ensureread();

    void write(const char* str, const char* name) {
      //cast is harmless
      write(BaseList<char>(strlen(str)+1,const_cast<char*>(str)),name);
    }

    template<class T> void write(const T&, const char* name);
    template<class T> void read(T&);

    bool next();
    bool peek(header_info&);
    bool find(header_info&, const char*);
    void close();
    bool isopen() const { return (fp!=NULL); }
    void ensureopen() const {
      if (!isopen())
	throw Failed("matlab_controller:: file has been closed");
    }
    FILE* file_pointer() const {
      ensureopen();
      return fp; 
    }
    long ftell() const;
    bool isreading() const { 
      ensureopen();
      return isread; }
    void fseek(long) const;

    static Warning<> wide_characters_warning;
    static Warning<> ignoring_comment_warning;
    static Warning<> underrun_warning;
    static Warning<> not_a_matrix_warning;
    static Warning<> padding_warning;
    static Warning<> compressed_warning;

    padmode_t padmode;

    class composite {
    public:
      composite(matlab_controller&, const char*, matlab_t =STRUCT, bool =false);
      composite(matlab_controller&, const char*, matlab_t, int, bool =false);
      composite(matlab_controller&, const char*, matlab_t, const BaseList<size_t>&, bool =false);

      composite(composite&, const char*, matlab_t =STRUCT, bool =false);
      composite(composite&, const char*, matlab_t, int, bool =false);
      composite(composite&, const char*, matlab_t, const BaseList<size_t>&, bool =false);

      composite(matlab_controller&); //!< read
      composite(composite&); //!< read
      
      ~composite() { close(); }

      static Warning<> ignoring_item_name_warning;
      static Warning<> write_overrun_warning;
      static Warning<> intruding_write_warning;
      static Warning<> name_length_warning;
      static Warning<> write_mismatch_warning;
      
      matlab_controller& controller() { return ctrl_; }
      const matlab_controller& controller() const { return ctrl_; }
      
      template<typename T> void write(T, const char* name ="");
      bool peek(header_info&);
      bool find(header_info&, const char*);
      
      template<typename T> void read(T&);

      bool next();
      bool ok_to_close() const; //!< return false if object cannot be successfully closed
      void close();      
      void ensurewrite() const;
      void ensureread() const;
      void ensureactive() const;
      void ensureknown() const;
      
    private:
      void checkupdate() { checkpos_=ctrl_.ftell(); }
      void checkposition() const;
      void validate(matlab_t);
      void write_header(matlab_t, const BaseList<size_t>&, const char*);
      void checkwritestart(const char*) const;
      void increment();
      void read_init();
      
      matlab_controller& ctrl_;
      bool disablecheck_;
      matlab_class_t class_;
      size_t max_;
      size_t count_;
      bool active_;
      List< smartptr<Matlab5WriteStoreBase> > writestack_;
      long initpos_,checkpos_,namespos_;
      size_t namelen_;
      char namebuf[matlab_controller_maxstructlength];
    }; // declaration of matlab_controller::composite

  private:
    FILE* fp;
    int version;
    const bool isread;
    bool isbigend;
    bool singleprec;

    long endpos;
    Matlab5Tag Dims; //used to (temporarily) cache read results
    Matlab5Tag arrayname;
    List<char> lastname;

    static int get_version(FILE*);

    //disallow copy
    matlab_controller(const matlab_controller&);
    matlab_controller& operator= (const matlab_controller&); 

    template<class T> void WriteMATLAB5_(const T& a, const char* name, Bool2Type<false>) {
      WriteMATLAB5(fp,a,name);
    }
    template<class T> void WriteMATLAB5_(const T& a, const char* name, Bool2Type<true>) {
      if (singleprec)
	WriteMATLAB5_single(fp,a,name);
      else
	WriteMATLAB5(fp,a,name);
    }    
    void ReadMATLAB5Header(BaseList<size_t> dims, matlab_t&, bool& iscomplex);
    void read_matrix_header(header_info&);
    long read_size(bool =true);
    void flush();
    template<class Iter> void read_array_matlab5(const Iter& start, size_t nitems);
    template<class T> void ReadArray5(Matrix<T>& dest,bool iscomplex);
    template<class T> void ReadArray5(List<T>& dest,bool iscomplex);
    template<size_t N> void ReadArray5(MultiMatrix<complex,N>& dest,bool iscomplex);
    template<size_t N> void ReadArray5(MultiMatrix<double,N>& dest,bool iscomplex);

    template<class T> void ReadMATLAB5(List<T>&);
    template<class T> void ReadMATLAB5(Matrix<T>&);
    template<class T,size_t N> void ReadMATLAB5(MultiMatrix<T,N>&);
    template<class T> void ReadMATLAB5(MultiHolder<T>&);
    template<class T> void ReadMATLAB5(List< Matrix<T> >&);
    void ReadMATLAB5(List<char>&);
  };

  FILE* open_matlab(const char* fname,const char* mode);

  void TidyMATLAB5Header(FILE*, long);
  void WriteMATLAB5Tag(FILE*, UINT32_t,UINT32_t bytes,const void*, bool trycompress =false);

  template<typename T> void WriteMATLAB5Data(FILE* fp,const BaseList<T>& data) {
    WriteMATLAB5Tag(fp,matlab_traits<T>::dtype,data.size()*sizeof(T),data.vector(),true);
  }

  template<> void WriteMATLAB5Data(FILE* fp,const BaseList<char>&);

  template<typename T> void WriteMATLAB5(FILE* fp,const BaseList<T>& data,const char* name)
  {
    const long initpos=WriteMATLAB5MatrixHeader(fp,matlab_traits<T>::dclass,false,ScratchList<size_t>(size_t(1),data.size()),name);
    WriteMATLAB5Data(fp,data);
    TidyMATLAB5Header(fp,initpos);
  }

  template<typename T> void WriteMATLAB5(FILE* fp,const BaseList< Matrix<T> >& data,const char* name);

  template<class M> void WriteMATLAB(const char* fname,const M& a,const char* name =NULL, int version =0)
  {
    const char* nameptr= name ? name : getbasename(fname);
    if (!version)
      version=matlab_default<M>::version;
    matlab_controller ctrl(fname,version);
    ctrl.write(a,nameptr);
  }
  
  template<class M> void ReadMATLAB(M& dest,const char* name)
  {
    matlab_controller ctrl(name);
    ctrl.read(dest);
  }

  void ReadMATLABHeader(int& rows,int& cols,bool& iscomplex,FILE *);
    
  template<typename T> struct matlab_composite_traits {
    typedef T write_store_type;
  };
  template<> struct matlab_composite_traits<const char*> {
    typedef std::string write_store_type;
  };

  template<typename T> struct Matlab5WriteStore : public Matlab5WriteStoreBase {
    template<typename Ta> Matlab5WriteStore(const Ta& datav, const char* namev)
      : Matlab5WriteStoreBase(namev),
	data(datav) {}
    
    void operator()(matlab_controller& ctrl) const { ctrl.write(data,""); } //!< don't use name - stored elsewhere in field
    Matlab5WriteStoreBase* clone() const { return new Matlab5WriteStore(*this); }
    
    T data;
  };

  template<> struct Matlab5WriteStore<std::string> : public Matlab5WriteStoreBase {
    template<typename Ta> Matlab5WriteStore(const Ta& datav, const char* namev)
      : Matlab5WriteStoreBase(namev),
	data(datav) {}
    
    void operator()(matlab_controller& ctrl) const { ctrl.write(data.c_str(),""); } //!< don't use name - stored elsewhere in field
    Matlab5WriteStoreBase* clone() const { return new Matlab5WriteStore(*this); }
    
    std::string data;
  };

  template<typename T> void matlab_controller::composite::read(T& a)
  {
    ensureread();
    if (count_==max_)
      throw Failed("matlab_controller::composite::read: object is exhausted");
    ctrl_.read(a);
    count_++;
  }

  template<typename T> void matlab_controller::composite::write(T a, const char* name)
  {
    checkwritestart(name);
    if (max_ || (class_==mxCELL_CLASS)) {      
      ctrl_.write(a,(class_==mxSTRUCT_CLASS) ? "" : name);
      count_++;
      checkupdate();
    }
    else {
      writestack_.push_back(smartptr<Matlab5WriteStoreBase>());
      writestack_.back().reset(new Matlab5WriteStore< typename matlab_composite_traits<T>::write_store_type >(a,name));
    }
  }
  
  template<class T> void matlab_controller::write(const T& a, const char* name) 
  {
    if (version==4)
      WriteMATLAB4(fp,a,name);
    else
      WriteMATLAB5_(a,name,Bool2Type<matlab_allowsingle< LCM_VAL(T) >::result>());
  }
  
  template<class T> void matlab_controller::read(T& a)
  {
    ensureread();
    const long initpos=ftell();
    try {
      if (version==4)
	::libcmatrix::ReadMATLAB4(a,fp);
      else
	ReadMATLAB5(a);
    } catch(...) {
      ::std::fseek(fp,initpos,SEEK_SET); //if read fails, reset to start of header
      throw;
    }
  }

  struct Fmatrix {
    int32_t type;
    int32_t mrows;
    int32_t ncols;
    int32_t imagf;
    int32_t namlen;
  };
  
} //namespace libcmatrix

#endif
