#ifndef _basespin_system_h_
#define _basespin_system_h_

#include <iostream>
#include "List.h"
#include "cmatrix.h"

namespace libcmatrix {

  //! sanity check on number of spin system states to prevent things blowing up later
#define NMRSIM_MAX_STATES 10000

const size_t NULL_NUCLEUS=0;
extern const size_t MAX_NUCLEUS;
  extern const size_t H_NUCLEUS;
 extern const size_t DEFAULT_NUCLEUS;

 //! default nucleus table
#ifndef LCM_DEFAULT_NUCLEUS_PROPERTIES
#define LCM_DEFAULT_NUCLEUS_PROPERTIES "HarrisIUPAC"
#endif

 void set_nucleus_properties(const char*); //!< change which table of nuclear spin properties to use

typedef unsigned long state_t;

size_t degeneracy(float);
  size_t maximumnuc();
size_t labeltonuc(const char *); //!< convert nucleus name into internal id
const char* nuctolabel(size_t nuc);  //!< return isotope name given internal id \a nuc
  void define_nucleus(const char* label, size_t deg, double gamma); //!< define new nucleus type


  struct nuclei_spec {
    size_t value;

    nuclei_spec(size_t);
    nuclei_spec(const char* label_) : value(labeltonuc(label_)) {}
    size_t operator()() const { return value; }
  };

  inline bool operator== (const nuclei_spec& spec1, const nuclei_spec& spec2) 
    { return spec1.value==spec2.value; }

  inline bool operator!= (const nuclei_spec& spec1, const nuclei_spec& spec2) 
    { return spec1.value!=spec2.value; }

  inline std::ostream& operator<< (std::ostream& ostr, const nuclei_spec& a)
    { return ostr << nuctolabel(a.value); }


 double gamma(const nuclei_spec&);

  template<class T> class replicator {
  public:
    typedef T value_type;
        
    replicator(const BaseList<T>& objv, size_t nv)
      : obj_(objv), n_(nv), size_(objv.size()*nv) {
      if (nv<1)
	throw InvalidParameter("replicator: zero count");
    }
    
    class iterator 
      : public ::std::iterator< ::std::forward_iterator_tag,T,ptrdiff_t,const T*,const T&> {
    public:
      iterator(Bool2Type<true>, const BaseList<T>& objv,size_t nv) 
	: start_(objv.vector()), stop_(start_+objv.size()),
	  curp_(start_), curmaj_(nv) {}

      iterator(Bool2Type<false>, const BaseList<T>& objv,size_t) 
	: curp_(objv.vector()), curmaj_(-1) {}

      const T& operator*() { return *curp_; }
      iterator& operator++() { up(); return *this; }
      iterator operator++(int) { iterator tmp(*this); up(); return tmp; }

      bool operator==(const iterator& x) { return (curp_==x.curp_) && (curmaj_==x.curmaj_); }
      bool operator!=(const iterator& x) { return (curp_!=x.curp_) || (curmaj_!=x.curmaj_); }
      
    private:
      const T* const start_;
      const T* const stop_;
      const T* curp_;
      int curmaj_;

      void up() {
	if (++curp_==stop_) {
	  curp_=start_;
	  curmaj_--;
	}
      }
    };

    typedef iterator const_iterator;

    iterator begin() const { return iterator(Bool2Type<true>(),obj_,n_); }
    iterator end() const { return iterator(Bool2Type<false>(),obj_,n_); }
    
    size_t size() const { return size_; }

  private:
    const BaseList<T>& obj_;
    const size_t n_;
    const size_t size_;
  };

 template<typename T> struct type_traits< replicator<T> > {
   static const bool trivialconstructor=false;
   static const size_t dimensionality=1;
   static const size_t rank=0;
   typedef T value_type;
 };

class spin {
public:
  explicit spin(size_t nuc =DEFAULT_NUCLEUS, bool restrict =false) { nucleus(nuc,restrict); }
  explicit spin(const char *lab, bool restrict =false) { isotope(lab,restrict); }
  
  double gamma() const;
  float qn() const;
  size_t deg() const;
  size_t effective_deg() const { return restrict_ ? 2 : deg(); }

  void isotope(const char *lab, bool restrict =false) {
    nucleus(labeltonuc(lab),restrict);
  }
  const char* isotope() const {
    return nuctolabel(nucleus_);  }

  size_t nucleus() const { return nucleus_; }
  void nucleus(size_t nuc, bool restrict =false);

  friend std::ostream& operator<< (std::ostream&, const spin&);
  bool isrestricted() const { return restrict_; }

 private:
  size_t nucleus_;
  bool restrict_;
};

 inline bool operator== (const spin& spin1, const spin& spin2)
   { return (spin1.nucleus()==spin2.nucleus()); }

 inline bool operator== (const spin& a, const char* label)
   { return (strcmp(a.isotope(),label)==0); }

 inline bool operator!= (const spin& spin1, const spin& spin2)
   { return (spin1.nucleus()!=spin2.nucleus()); }

 inline bool operator!= (const spin& a, const char* label)
   { return (strcmp(a.isotope(),label)!=0); }

  class Permutation : public DynamicList<size_t> {
  public:
    explicit Permutation(const BaseList<size_t>&);
    explicit Permutation(const char*);
    Permutation() : order_(1) {}

    size_t order() const;
    template<typename T> void apply(BaseList<T>, const BaseList<T>&) const;
    template<typename T> void apply(List<T>& d, const BaseList<T>& a) const
    { d.create(this->size()); apply(static_cast<BaseList<T> >(d),a); }
    template<typename T> void apply(Matrix<T>&, const Matrix<T>&) const;
    static bool isidentity(const BaseList<size_t>&);
    static bool isvalid(const BaseList<size_t>&);
  private:
    mutable size_t order_;
  };

  std::ostream& operator<< (std::ostream&, const Permutation&);

  class spinorder_spec {
  public:

    typedef std::pair<state_t,state_t> subspec_t;
    
    spinorder_spec(size_t nspinsv =0, bool iscoherv =false) { init(nspinsv, iscoherv); }
    
    spinorder_spec(size_t nspinsv, const BaseList<size_t>& whichv , const BaseList<size_t>& ns)
      { init(nspinsv), add(whichv,ns); }

    spinorder_spec(size_t nspinsv, const BaseList<size_t>& whichv , const BaseList<int>& ns)
      { init(nspinsv,true), add(whichv,ns); }
    
    spinorder_spec(size_t nspinsv, const BaseList<size_t>& ns) { create(nspinsv,ns); }
    spinorder_spec(size_t nspinsv, const BaseList<int>& ns) { create(nspinsv,ns); }
    
    spinorder_spec(size_t nspinsv, const BaseList<size_t>& whichv, size_t n)
      { init(nspinsv), add(whichv,BaseList<size_t>(1,&n)); }
    spinorder_spec(size_t nspinsv, const BaseList<size_t>& whichv, int n)
      { init(nspinsv,true), add(whichv,BaseList<int>(1,&n)); }
    
    spinorder_spec(size_t nspinsv, size_t m) { create(nspins,BaseList<size_t>(1,&m)); }
    spinorder_spec(size_t nspinsv, int m) { create(nspins,BaseList<int>(1,&m)); }

    void add(const BaseList<size_t>&, const BaseList<size_t>&);
    void add(const BaseList<size_t>&, const BaseList<int>&);

    bool overlappingranges() const; //!< return \c true if index ranges overlap (may be error)

    bool isactive(size_t) const; //!< return \c true if spin is part of order selection
    bool operator()(state_t, state_t) const; //!< return \c true if coherence between two states falls within order specification
    size_t size() const { return specs.size(); }
    const subspec_t& item(size_t n) const { return specs(n); }

    bool isspinorder() const { return !iscoher; }
    bool iscoherence(double, state_t) const; //!< is coherence within order list

    size_t bit_to_spin(size_t i) const { return nspins-int(i)-1; }
    void print(std::ostream&) const;
    bool allspinsused() const { return (totalspinmask==((1<<nspins)-1)); }

    static state_t compress_spinindices(size_t nspins, const BaseList<size_t>&, state_t spinmask =0);

  private:
    List<subspec_t> specs;
    size_t nspins;
    state_t totalspinmask;
    //    state_t totalordermask; //!< cumulative order/coherence mask makes little sense
    bool iscoher;

    void init(size_t nspinsv, bool =false);

    void create(size_t, const BaseList<size_t>&);
    void create(size_t, const BaseList<int>&);
    template<class T> void add_(size_t, state_t, const BaseList<T>& ordersv);

    void dumpbits(std::ostream& , state_t, size_t =0, bool reverse =false) const;
    void dumpcoherences(std::ostream& , state_t) const;
    state_t coherencemask(int) const;
    state_t compressorders(size_t, const BaseList<size_t>&) const;
    state_t compressorders(size_t, const BaseList<int>&) const;
    void registerorders(state_t, state_t);
 };


  inline std::ostream& operator << (std::ostream& ostr, const spinorder_spec& z) { z.print(ostr); return ostr; }

class basespin_system : protected List<spin> {
public:
  basespin_system(int nspins, const spin&);

  typedef spin value_type;
  typedef cmatrix operator_type;
  typedef List<double> diag_operator_type;
  
  virtual ~basespin_system() {};
  virtual basespin_system* clone() const =0;
  virtual basespin_system* clone(size_t) const =0;

  size_t dim_range(size_t,size_t) const;

  size_t length() const { 
    return List<spin>::length(); } //makes it look like a container!

  size_t nspins() const {
    return List<spin>::length();
  }
  size_t nspins(nuclei_spec) const;

  void isotope(size_t n, const spin& nuc) { 
    nucleus(n,nuc);
  }
//   void isotope(size_t n,nuclei_spec nuc) {
//     nucleus(n,spin(nuc()));
//   }
//   void isotope(size_t,size_t,nuclei_spec);
  void isotope(size_t,size_t,const spin&);

  size_t homonucleus() const {
    if (ishomonuclear())
      return List<spin>::front().nucleus();
    throw Failed("homonucleus: spin system is not homonuclear");
  }

  bool ishomonuclear() const {
    return isconstant(static_cast<const List<spin>& >(*this)); }

  virtual bool isspinhalfonly() const =0;

  const spin& operator()(size_t n) const { 
    if (n>=length())
      throw BadIndex("Bad spin index",n,length());
    return List<spin>::operator()(n); }

  spin& operator()(size_t n) { 
    if (n>=length()) 
      throw BadIndex("Bad spin index",n,length());
    return List<spin>::operator()(n); }
 
  size_t size() const;
  virtual size_t rows() const { 
    return size(); }
  virtual size_t cols() const {
    return size(); }
  virtual bool isdiagonal() const { return true; }

  void validate_diagonal() const {
    if (!isdiagonal()) 
      throw Failed("Operation not valid for diagonal blocks"); 
  }

  virtual void permutation_vectorH(BaseList<state_t>&, const Permutation&) const;
  void permutation_vectorL(BaseList<state_t>&, const Permutation&) const;
  Matrix<float_t> permutation_matrixH(const Permutation&) const;
  Matrix<float_t> permutation_matrixL(const Permutation&) const;

  state_t permute(state_t state,const Permutation& permvec) const { 
    validate_diagonal();
    return rawpermute(state,permvec); }

  void expandop_tmp(cmatrix&, size_t, cmatrix&) const;
  void expandop(cmatrix&, size_t, const cmatrix&) const;
  void expandop(Matrix<float_t>&, size_t, const Matrix<float_t>&) const;
  void mla_expandop(BaseList<double> &,size_t,const BaseList<double> &) const;

  virtual void mla_I(cmatrix &,double,const BaseList<char>&) const =0;
  virtual void mla_I(cmatrix &,double,size_t,char) const =0;
  void mla_I(cmatrix&, double, size_t, size_t, size_t) const;

  //! single transition operators
  virtual void mla_ST(cmatrix&, double, size_t, size_t, size_t) const
    { throw Failed("class does not implement single transition operators"); }
  virtual void mla_ST(List<double>&, double, size_t, size_t) const
    { throw Failed("class does not implement single transition operators"); }

  virtual void mla_I(cmatrix &,double,size_t,char,size_t,char) const =0;
  virtual void mla_I(cmatrix &,double,size_t,const BaseList<char>&) const =0;

  void mla_F(cmatrix&, double scale, size_t nuc, size_t r, size_t c) const; //!< add single transition operator \a r,\a c for nucleus \a nuc
  void mla_F(cmatrix &d,double scale,char op) const {
    mla_F(d,scale,NULL_NUCLEUS,op);
  }
  void mla_F(cmatrix &,double,size_t,char) const;
  void mla_Fz(BaseList<double> &,double,size_t =NULL_NUCLEUS) const;

  void mla_Iz(BaseList<double>& dest,double scale,size_t which) const {
    rawmla_Iz(dest,scale,which);
  }
  void mla_Iz(BaseList<double>& dest,double scale,size_t which1,size_t which2) const { 
    rawmla_Iz(dest,scale,which1,which2);
  }
  void mla_Iz(List<double>& dest,double scale,size_t which) const;
  void mla_Iz(List<double>& dest,double scale,size_t which1, size_t which2) const;
  void mla_Fz(List<double>& dest,double scale,size_t which =NULL_NUCLEUS) const;

  static char squeezeST(size_t,size_t); //!< convert spin-1/2 single transition operator to corresponding 'normal' operator (a,+,-,b)

protected:
  basespin_system(const basespin_system& sys, size_t n)
    : List<spin>(replicator<spin>(sys,n)) {}
  
 private:
  virtual state_t rawpermute(state_t,const Permutation&) const =0;

  virtual void rawmla_Iz(BaseList<double> &,double,size_t) const =0;
  virtual void rawmla_Iz(BaseList<double> &,double,size_t,size_t) const =0;

  virtual void nucleus(size_t n, const spin&); //virtual to allow restriction of allowed spin types
  };

 template<class T> struct spin_system_traits {};
 
std::ostream& operator << (std::ostream&, const basespin_system&);

const cmatrix& spinhalfop(char);
  void spinop(cmatrix &,int,char,double =1.0, bool restrict =false);
  void zspinop(List<double> &,int,double =1.0, bool restrict =false);

cmatrix I(const basespin_system &,size_t,char);
  cmatrix I(const basespin_system &,size_t,size_t,size_t); //!< single transition operator
cmatrix I(const basespin_system &,size_t,const char []);
cmatrix I(const basespin_system &,const BaseList<char> &);
cmatrix I(const basespin_system &,size_t,char,size_t,char);
 
List<double> diag_Iz(const basespin_system &,size_t);
List<double> diag_Iz(const basespin_system &,size_t,size_t);
List<double> diag_Fz(const basespin_system &,nuclei_spec =NULL_NUCLEUS);

cmatrix F(const basespin_system&, nuclei_spec, char);
  cmatrix F(const basespin_system&, nuclei_spec, size_t, size_t);
   
inline cmatrix F(const basespin_system& sys, char op) {
  return F(sys,NULL_NUCLEUS,op); }

  // This form creates ambiguities
//   inline cmatrix F(const basespin_system& sys, size_t r, size_t c) {
//     return F(sys,NULL_NUCLEUS,r,c); }

template<typename T> void Permutation::apply(BaseList<T> d, const BaseList<T>& a) const
{  
  if (a.empty())
    throw Undefined("Permutation::apply");
  size_t i=this->size();
  if (i) {
    if ((i!=d.size()) || (i!=a.size()))
      throw Mismatch("Permutation::apply");
    for (;i--;)
      d(i)=a((*this)(i));
  }
  else
    d=a;
}

template<typename T> void Permutation::apply(Matrix<T>& d, const Matrix<T>& a) const
{  
  if (!!d && issame(d,a))
    throw ArgumentClash("Permutation::apply");
  const size_t M=a.rows();
  const size_t N=a.cols()/M;
  const size_t total=M*N;
  if (total!=a.cols())
    throw InvalidParameter("Permutation::apply");
  d.create(M,total);
  for (size_t j=M;j--;) {
    const size_t destj=(*this)(j);
    for (size_t k=M;k--;) {
      size_t sourcek=k;
      size_t destk=(*this)(k);
      while (sourcek<total) {
	d(destj,destk)=a(j,sourcek);
	sourcek+=M;
	destk+=M;
      }
    }
  }
}

} //namespace libcmatrix

#endif

