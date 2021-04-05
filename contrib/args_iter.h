#ifndef lcm_args_iter_h_
#define lcm_args_iter_h_

// alternate interface

#ifdef HAVE_SSTREAM
#include <sstream>
#define LCM_STRTYPE std::istringstream
#else
#include <strstream>
#define LCM_STRTYPE std::istrstream
#endif

namespace libcmatrix {

 class args_iter {
   size_t cpos,epos;
   const char **argv;
   std::istream& istr;
   std::ostream& ostr;
   bool haveostr;
   LCM_STRTYPE* cstrp;

   char tmp[LINE_MAXLEN];
   char fnametmp[LINE_MAXLEN];

   std::istream& update();
   void attach(const char *);
   void init();
   bool doprompt() const { return (haveostr && (cpos>=epos)); }
   void next() { if (cpos<epos) cpos++; }

 public:
   args_iter(size_t start,size_t end,const char *argv[]);
   args_iter(size_t start,size_t end,const char *argv[],std::istream&);
   args_iter(size_t start,size_t end,const char *argv[],std::istream&,std::ostream&);
   args_iter();
   args_iter(std::istream&);
   args_iter(std::istream&, std::ostream&);
   ~args_iter() { if (cstrp) delete cstrp; }

   template<class T> std::istream& operator()(T&,const char *prompt);
   template<class T,class F> std::istream& operator()(T&,const char *prompt,const F& pred);
   template<class T> std::istream& operator()(T&,const char *prompt,const T& def);
   template<class T,class F> std::istream& operator()(T&,const char *prompt,const F& pred, const T& def);

   int count() const { return cpos; }
   void count(int npos) {
     if (npos>epos)
       throw InvalidParameter("args_iter::count");
     cpos=npos;
     update();
   }

   std::istream& operator()() {
     if (!cstrp)
       throw Failed("No input");
     return *cstrp;
   }
   
   std::istream& operator()(char*, size_t,size_t,const char *prompt,const char * =NULL);
   std::istream& operator()(char* dest,size_t max,const char *prompt,const char *def =NULL) { return operator()(dest,0,max,prompt,def); }
   std::istream& operator()(size_t&, const char *prompt,const char *);
   std::istream& operator()(size_t&, const char *prompt,const char *,size_t);
   std::istream& operator()(FILE* &, const char *prompt,const char *mode, char*, size_t, int =mxflag::confirmoverwrite);
   inline std::istream& operator()(FILE* &fp, const char *prompt,const char *mode, int flags =mxflag::confirmoverwrite) { return operator()(fp,prompt,mode,fnametmp,LINE_MAXLEN,flags); }

   template<class T> T get(const char *prompt) {
     T res;
     (*this)(res,prompt);
     return res;
   }
   template<class T,class F> T get(const char *prompt,const F& preddef) {
     T res;
     (*this)(res,prompt,preddef);
     return res;
   }
   template<class T,class F> T get(const char *prompt,const F& pred,const T& def) {
     T res;
     (*this)(res,prompt,pred,def);
     return res;
   }
   size_t get(const char* prompt,const char* opts) {
     size_t res;
     (*this)(res,prompt,opts);
     return res;
   }
   size_t get(const char* prompt,const char* opts, size_t def) {
     size_t res;
     (*this)(res,prompt,opts,def);
     return res;
   }
 };

 template<class T> class Within;
 template<class T> std::ostream& operator<< (std::ostream&, const Within<T>&);

template<class T> class Within : public ::std::unary_function<T,bool> {
  public:
  Within(const T& minv, const T& maxv)
    : min_(minv), max_(maxv) {
    if (maxv<=minv)
      throw InvalidParameter("Within");
  }
  bool operator()(const T& v) const {
    return (v>=min_) && (v<=max_); }

  friend std::ostream& operator<< <>(std::ostream&, const Within<T>&);
 private:
  const T min_,max_;
};
 
 template<class T> void writestream(std::ostream& ostr,const T& v) { ostr << v; }
 template<class T> bool readstream(T& v,std::istream& in) { return (in >> v); }

 template<> inline void writestream(std::ostream& ostr,const bool& v) { ostr << (v ? 'y' : 'n'); }
 template<> bool readstream(bool&, std::istream&);

 //"Compatability" functions to match traditional getXXX functions
 inline bool getlogical(args_iter& iter, const char* prompt) { bool res; iter(res,prompt); return res; }
 inline bool getlogical(args_iter& iter, const char* prompt, bool def) { bool res; iter(res,prompt,def); return res; }

 inline int getint(args_iter& iter, const char* prompt) { int res; iter(res,prompt); return res; }
 inline int getint(args_iter& iter, const char* prompt, int def) { int res; iter(res,prompt,def); return res; }
 inline int getint(args_iter& iter, const char* prompt, int min, int max) { int res; iter(res,prompt,Within<int>(min,max)); return res; }
 inline int getint(args_iter& iter, const char* prompt, int min, int max, int def) { return iter.get(prompt,Within<int>(min,max),def); }
 
 inline float getfloat(args_iter& iter, const char* prompt) { float res; iter(res,prompt); return res; }
 inline float getfloat(args_iter& iter, const char* prompt, float def) { float res; iter(res,prompt,def); return res; }
 inline float getfloat(args_iter& iter, const char* prompt, float min, float max) { float res; iter(res,prompt,Within<float>(min,max)); return res; }
 inline float getfloat(args_iter& iter, const char* prompt, float min, float max, float def) { return iter.get(prompt,Within<float>(min,max),def); }
 
 inline char* getstring(args_iter& iter, const char* prompt, char* dest, size_t max, const char* def =NULL)
  { iter(dest,max,prompt,def); return dest; }

 inline size_t getoption(args_iter& iter, const char* prompt, const char* opts) { return iter.get(prompt,opts); }
 inline size_t getoption(args_iter& iter, const char* prompt, const char* opts, size_t def) { return iter.get(prompt,opts,def); }

 // Implementation details below

template<class T> std::istream& args_iter::operator()(T& dest,const char *prompt)
{
  if (!doprompt()) {
    std::istream& cistr=update();
    if (readstream(dest,cistr)) {
      next();
      if (haveostr) {
	ostr << prompt;
	writestream(ostr,dest);
	ostr << ::std::endl;
      }
      return cistr;
    }
    ::std::cerr << prompt << " Failed to parse input" << ::std::endl;
    throw Failed("Failed to parse input");
  }
  for (;;) {
    ostr << prompt;
    std::istream& cistr=update();
    if (readstream(dest,cistr)) return cistr;
  }
}

template<class T> struct output_traits {
  static const bool isprintable=type_traits<T>::known;
};

//Not all predicates (e.g. from std library) have an output method - hence this hack
template<bool> struct print_ifpossible_ {
  template<class F> static inline void print(std::ostream& ostr,const F& obj) { ostr << obj; }
};
template<> struct print_ifpossible_<false> {
  template<class F> static inline void print(std::ostream& ostr,const F&) {}
};

 template<class F> void print_ifpossible(std::ostream& ostr,const F& pred) { print_ifpossible_<output_traits<F>::isprintable>::print(ostr,pred); }

template<class T,class F> std::istream& args_iter::operator()(T& dest,const char *prompt,const F& pred)
{
  if (doprompt()) {
    for (;;) {
      std::istream& cistr=operator()(dest,prompt);
      if (pred(dest))
	return cistr;
      ostr << "Failed input criterion: ";
      print_ifpossible(ostr,pred);
      ostr << ::std::endl;
    }
  }
  std::istream& cistr=operator()(dest,prompt);
  if (!pred(dest)) {
    ::std::cerr << prompt << "Value (" << dest << ") failed input criterion ";
    print_ifpossible(::std::cerr,pred);
    ::std::cerr << ::std::endl;
    throw Failed("Parameter out of range");
  }
  return cistr;
}

template<class T,class F> std::istream& args_iter::operator()(T& dest,const char *prompt,const F& pred,const T& def)
{
  if (!pred(def))
    throw InvalidParameter("default fails input criterion!");

  if (doprompt()) {
    for (;;) {
      std::istream& cistr=operator()(dest,prompt,def);
      if (pred(dest))
	return cistr;
      ostr << "Failed input criterion: ";
      print_ifpossible(ostr,pred);
      ostr << ::std::endl;
    }
  }
  std::istream& cistr=operator()(dest,prompt,def);
  if (!pred(dest)) {
    ::std::cerr << prompt << "Value (" << dest << ") failed input criterion ";
    print_ifpossible(::std::cerr,pred);
    ::std::cerr << ::std::endl;
    throw Failed("Parameter out of range");
  }
  return cistr;
}

template<class T> std::istream& args_iter::operator()(T& dest,const char *prompt,const T& def)
{
  if (!doprompt()) {
    std::istream& cistr=update();

    char c=cistr.get();
    const bool isdef=((c=='-') && (cistr.peek()<0));

    if (isdef)
      dest=def;
    else {
      cistr.putback(c);
      if (!readstream(dest,cistr)) {
	::std::cerr << prompt << " Failed to parse input" << ::std::endl;
	throw Failed("Failed to parse input");
      }
    }
    next();
    if (haveostr) {
      ostr << prompt;
      writestream(ostr,dest);
      ostr << ::std::endl;
    }
    return cistr;
  }
  for (;;) {
    ostr << prompt << "[";
    writestream(ostr,def);
    ostr << "] ";
    std::istream& cistr=update();
    if (cistr.peek()<0) {
      dest=def;
      return cistr;
    }
    if (readstream(dest,cistr)) return cistr;
  }
}

template<class T> std::ostream& operator<< (std::ostream& ostr, const Within<T>& a)
{ return ostr << a.min_ << " <= x <= " << a.max_; }

template<class T> struct output_traits< Within<T> > { static const bool isprintable=true; };

template<class T> struct GreaterThanEqual : public ::std::unary_function<T,bool> {
  GreaterThanEqual(const T& minv) : min_(minv) {}

  bool operator()(const T& v) const {
    return (v>=min_); }

  const T min_;
};

template<class T> std::ostream& operator<< (std::ostream& ostr, const GreaterThanEqual<T>& a) {
  return ostr << "x >= " << a.min_; }

template<class T> struct output_traits< GreaterThanEqual<T> > { static const bool isprintable=true; };

} //namespace libcmatrix

#endif
