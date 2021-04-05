#ifndef LCM_HamiltonianStore_h_
#define LCM_HamiltonianStore_h_

#include <map>
#include "NMR.h"
#include "MultiMatrix.h"
#include "smartptr.h"
#include "Warnings.h"

namespace libcmatrix {
  
  struct MetaFlags {
    static const int minflagvalue=8;
    enum { UseMzSymmetry=8,
	   UseEigSymmetry=16, 
	   DontCacheBinary=32,
	   UsePartitioning=64,
	   ClassicSecondOrder=128
    };
    static bool extract(int& flags, int flag, const char* warnmess =NULL);
    static Warning<> UseMzSymmetry_shift_warning;
    static Warning<> UseMzSymmetry_heteronuclear_warning;
    static Warning<> UseMzSymmetry_oddblocks_warning;
  };

  template<typename T> class HamiltonianStore;

  typedef size_t interaction_t;

  const interaction_t I_INVALID=0;
  const interaction_t I_CS=1; //must start at 1
  const interaction_t I_DIPOLE=2;
  const interaction_t I_J=3;
  const interaction_t I_QUAD=4;
  const interaction_t I_AUTO=I_QUAD;
  
  const std::string& interaction_name(interaction_t);
  interaction_t interaction_register(const std::string&);
  interaction_t interaction_find(const std::string&);

  template<typename T> struct couplingstore {
    Matrix<double> A0;
    Matrix<T> A2;
    bool operator!() const { return (!A0 && !A2); }

    bool isnonnull(size_t i,size_t k) const { return (!!A0 && A0(i,k)) || (!!A2 && !!A2(i,k)); }
    bool areequal(size_t i,size_t k, size_t i2,size_t k2, double tol =1e-9) const {
      return (!A0 || ::libcmatrix::areequal(A0(i,k),A0(i2,k2),tol)) && (!A2 || ::libcmatrix::areequal(A2(i,k),A2(i2,k2),tol));
    }
    bool isequal(const couplingstore<T>& a, double tol) const {
      return ::libcmatrix::areequal(A0,a.A0,tol) && ::libcmatrix::areequal(A2,a.A2,tol); 
    }
    void set(size_t, size_t, const T&);
    size_t rows() const { return !A2 ? A0.rows() : A2.rows(); }
    size_t cols() const { return !A2 ? A0.cols() : A2.cols(); }
    void apply(const couplingstore&, const Permutation&);
  };

  struct CrystalStructure_iterator {
    virtual ~CrystalStructure_iterator() {}
    virtual void reset(size_t,size_t) =0;
    virtual bool next(size_t&, size_t&) =0;
    virtual CrystalStructure_iterator* clone() const =0;
  };

  class CrystalStructure {
  public:
    explicit CrystalStructure(int =1, const Permutation& =Permutation());
    explicit CrystalStructure(const BaseList<size_t>&, const Permutation& =Permutation());
    size_t dimensions() const { return ndims; }
    size_t dimension(size_t n) const { return indexer_.dimension(n); }
    template<int N> size_t dimension(Int2Type<N> n) const { return indexer_.dimension(n); }
    size_t ncells() const { return indexer_.size(); }
    size_t operator()(size_t, size_t, size_t) const;
    size_t inverse(size_t, size_t, size_t) const;
    template<int N> size_t inverse(size_t cell, Int2Type<N> n) const {
      return cell ? indexer_.dimension(n)-cell : 0;
    }
    size_t operator()(const BaseList<size_t>&) const;
    CrystalStructure_iterator* create_iterator(size_t) const;
    const Indexer<3>& indexer() const { return indexer_; }
    bool haspermutation() const { return !perm_.empty(); }
    const Permutation& permutation() const { return perm_; }

  private:
    Indexer<3> indexer_;
    Permutation perm_;
    size_t ndims;
    bool higher;
    friend std::ostream& operator << (std::ostream&, const CrystalStructure&);
  };

  std::ostream& operator << (std::ostream&, const CrystalStructure&);

  class HamiltonianStructure {
  public:
    template<typename T> HamiltonianStructure(const basespin_system&, const HamiltonianStore<T>&, const BaseList<interaction_t>& =BaseList<interaction_t>(), int flagsv =0);
    template<typename T> HamiltonianStructure(const basespin_system&, const HamiltonianStore<T>&, const BaseList<nuclei_spec>&, const BaseList<interaction_t>& =BaseList<interaction_t>(), int flagsv =0);
    template<typename T> HamiltonianStructure(const basespin_system&, const HamiltonianStore<T>&, const char* label, const BaseList<interaction_t>& =BaseList<interaction_t>(), int flagsv =0);

    size_t nspins() const { return iscoupled_.cols(); }  
    size_t nspins_cell() const { return M_; } 
    size_t ncells() const { return iscoupled_.cols()/M_; }
    size_t spin_to_spin(size_t sk) const { return sk % M_; } 
    bool isdiagonal() const { return isdiagonal_; }
    size_t quadrupole_order() const { return quadrupole_order_; }
    const Matrix<bool> iscoupled() const { return iscoupled_; }
    bool iscoupled(size_t i,size_t j) const { 
      if (i>=j)
	throw InvalidParameter("iscoupled");
      return iscoupled_(i,j);

    }
    bool isweak(interaction_t id) const { return (std::find(weakints_.begin(),weakints_.end(),id)!=weakints_.end()); }
    bool iscomplex() const { return (quadrupole_order_!=1) && !(flags_ & MetaFlags::ClassicSecondOrder); }
    bool isclassicsecondorder() const { return (flags_ & MetaFlags::ClassicSecondOrder) && (quadrupole_order_!=1); }
    ns_flag nonsecular_type(size_t&, size_t&) const;
    bool nonsecular(size_t m) const { return nonsecular_(m); }
    const basespin_system& spinsystem() const { return *sysp_; } 
    void blockingnuclei(const BaseList<nuclei_spec>&);
    const BaseList<nuclei_spec>& blockingnuclei() const { return blockingnuclei_; }
    bool isblocked(nuclei_spec nuc) const { return (std::find(blockingnuclei_.begin(),blockingnuclei_.end(),nuc)!=blockingnuclei_.end()); }
    int flags() const { return flags_; }

  private:
    smartptr<const basespin_system> sysp_;

    List<nuclei_spec> blockingnuclei_;
    ScratchList<interaction_t> weakints_;
    int flags_;
    size_t M_;
    bool isdiagonal_;
    size_t quadrupole_order_;
    List<bool> nonsecular_;
    Matrix<bool> iscoupled_;
    
    template<class T> void create(const HamiltonianStore<T>&);
    template<class T> void process_binary(const couplingstore<T>&);
  };

  std::ostream& operator<< (std::ostream&, const HamiltonianStructure&);

  template<class T> struct spinsystem_glue_ {
    static const basespin_system& spinsystem(const T& a) { return a; }
    static size_t ncells(const T& a) { return 1; }
    static const bool ispurereal=true;
    
    template<class M> static bool arematching(const T& a, const HamiltonianStore<M>& store) {
      return (store.ncells()==1) && (a.nspins()==store.nspins());
    }   
  };
    
  template<class T> inline bool isnonnull(const T&) { throw InternalError("isnonnull"); };
  template<> inline bool isnonnull(const double& v) { return (v!=0.0); }
  template<> inline bool isnonnull(const space_T& v) { return !!v; }

  template<class T> inline bool isnonnull(const std::pair<double,const T& >& a) { return ((a.first) || isnonnull2(a.second)); }

  template<class T> struct mapcache_ {
    typedef std::map<interaction_t,T> map_t;

    mapcache_(map_t& mapv)
      : map_(mapv),lasttype_(I_INVALID), lastp_(NULL) {}

    const T* find(interaction_t inttype) const {
      if (inttype!=lasttype_) {
	const typename map_t::const_iterator iter(map_.find(inttype));
	if (iter==map_.end())
	  return NULL; //don't update
	lastp_=const_cast<T*>(&(iter->second));
	lasttype_=inttype;
      }
      return lastp_;
    }
    T& operator[](interaction_t inttype) {
      if (inttype!=lasttype_) {
	lastp_=&(map_[inttype]);
	lasttype_=inttype;
      }
      return *lastp_;
    }

    map_t& map_;
    mutable interaction_t lasttype_;
    mutable T* lastp_;
  };

  template<class T> class HamiltonianStore {
  public:
    HamiltonianStore(int Mv, int Nv=1);
    HamiltonianStore(const HamiltonianStore<space_T>&, const Euler&, bool =false);
    HamiltonianStore(const HamiltonianStore<T>&, const Permutation&);

    void clear();
    size_t nspins() const { return total_; }
    size_t nspins_cell() const { return M_; }
    size_t ncells() const { return N_; }

    const couplingstore<T>& get_couplings(interaction_t id) const { return get_(couplingmap,id); }
    const List<T>& get_shifts(interaction_t id =I_CS) const { return get_(shiftmap,id); }

    const List<T>& get_quadrupole() const { return qstore; }
    const List<int>& get_quadrupole_order() const { return qorder; }
  
    template<class OpGen> bool ismatching(const OpGen& opgen) const {
      return (opgen.nspins()==nspins()) && (opgen.ncells()==ncells());
    }

    bool verify(std::ostream& ostr, double tol =1e-9) const { return verify(ostr,CrystalStructure(ExplicitList<1,size_t>(ncells())),tol); }
    bool verify(std::ostream&, const CrystalStructure&, double) const;
    bool isequal(const HamiltonianStore<T>&, double) const;

    bool haslinear() const { return !(shiftmap.empty()); }

    void split_couplings(Matrix<double>&, Matrix<T>&) const;

    void set_coupling(interaction_t inttype,size_t ni, size_t nj, double iso, const T& v);

    static bool isnonnull2(double v) { return (v!=0.0); }
    static bool isnonnull2(const space_T& v) { return v.have_rank(2); }

    bool isset(interaction_t id, size_t ni) const { return isnonnull(get_shift(ni,id)); }
    bool isset(interaction_t id, size_t ni, size_t nj) const {
      return (get_isotropic(id,ni,ni) || isnonnull(get_anisotropic(id,ni,nj)));
    }
    double get_isotropic(interaction_t, size_t ni, size_t nj) const;
    const T& get_anisotropic(interaction_t, size_t ni, size_t nj) const;

    void set_coupling(interaction_t, size_t, size_t, const T&);
    //    { throw Failed("set_coupling: not valid for HamiltonianStore<double>"); }

    void set_quadrupole(size_t ni, const T& v, size_t orderv);
    size_t get_quadrupole_order(size_t ni) const;

    const T& get_quadrupole(size_t ni) const {
      return get_unary(I_QUAD,ni);
    }

    void set_shift(size_t ni, const T& v, interaction_t inttype =I_CS) { set_unary(inttype,ni,v); }

    const T& get_shift(size_t ni, interaction_t inttype =I_CS) const {
      return get_unary(inttype,ni);
    }

    void print(std::ostream&) const;

    void invert_linear();

    size_t cell_to_spin(size_t n, size_t m) const 
    { return n*M_+m; }

    typedef std::map<interaction_t,couplingstore<T> > couplingmap_type;
    typedef std::map<interaction_t,List<T> > shiftmap_type;

    template<class M> const M& get_(const std::map<interaction_t,M>& map,interaction_t id)  const {
      const typename std::map<interaction_t,M>::const_iterator iter(map.find(id));
      if (iter==map.end())
	throw Failed("HamiltonianStore: interaction not present");
      return iter->second;
    }

    typename shiftmap_type::const_iterator shifts_end() const { return shiftmap.end(); }
    typename shiftmap_type::const_iterator shifts_begin() const { return shiftmap.begin(); }

    typename couplingmap_type::const_iterator couplings_end() const { return couplingmap.end(); }
    typename couplingmap_type::const_iterator couplings_begin() const { return couplingmap.begin(); }

    typename shiftmap_type::iterator shifts_end() { return shiftmap.end(); }
    typename shiftmap_type::iterator shifts_begin() { return shiftmap.begin(); }

    typename couplingmap_type::iterator couplings_end() { return couplingmap.end(); }
    typename couplingmap_type::iterator couplings_begin() { return couplingmap.begin(); }

    bool empty() const { return (couplingmap.empty() && shiftmap.empty() && qstore.empty()); }

  private:
    size_t M_,N_,total_;

    couplingmap_type couplingmap;
    shiftmap_type shiftmap;

    List<T> qstore; //Q-pole is special - coupling Hamiltonian but single spin
    List<int> qorder;

    static const T null_;
    
    mapcache_<couplingstore<T> > couplingmap_cache_;
    mapcache_<List<T> > shiftmap_cache_;

    template<class T2> static void create_null(List<T2>& d, size_t ni) { d.create(ni,T2(0)); }
    static void create_null(couplingstore<double>& d, size_t ni, size_t nj) { d.A2.create(ni,nj,0.0); }
    static void create_null(couplingstore<space_T>& d, size_t ni, size_t nj) { d.A2.create(ni,nj); }

    void docreate(int Mv,int Nv);
    
    static bool areequal(const couplingstore<T>& a, const couplingstore<T>& b, double tol) { return a.isequal(b,tol); }
    static bool areequal(const BaseList<T>& a, const BaseList<T>& b, double tol) { return ::libcmatrix::areequal(a,b,tol); }
    template<typename Key,typename Obj> static bool areequal(const std::map<Key,Obj>&, const std::map<Key,Obj>&, double tol);
    template<typename Key,typename Obj> static void apply(std::map<Key,Obj>&, const std::map<Key,Obj>&, const Permutation&);
    static void apply(List<T>& d, const BaseList<T>& a, const Permutation& perm) { perm.apply(d,a); }
    static void apply(couplingstore<T>& d, const couplingstore<T>& a, const Permutation& perm) { d.apply(a,perm); }
    static bool verify_(const couplingstore<T>&, size_t j, size_t k, size_t swapj, size_t swapk,const char*, std::ostream&, double);
    bool verify_(std::ostream&, const couplingstore<T>&, const char*, const CrystalStructure&, double) const;
  
    static void rotate(couplingstore<double>&, const couplingstore<space_T>&, const cmatrix&);
    static void rotate(couplingstore<space_T>&, const couplingstore<space_T>&, const cmatrix&);
    static void rotatesum(List<double>&, const List<space_T>&, const cmatrix&);
    static void rotatesum(List<space_T>&, const List<space_T>&, const cmatrix&);
  
    static void create_null(List<space_T>& d, size_t ni) { d.create(ni); }
    static void create_null(Matrix<space_T>& d, size_t ni, size_t nj) { d.create(ni,nj); }
  
    static void print_(std::ostream& ostr, const char* name, const BaseList<double>& store);
    static void print_(std::ostream& ostr, const char* name, const BaseList<space_T>& store);
    static void print_(std::ostream& ostr, const char* name, const couplingstore<double>& store);
    static void print_(std::ostream& ostr, const char* name, const couplingstore<space_T>& store);
    void printqpole_(std::ostream& ostr, const BaseList<double>& store) const;
    void printqpole_(std::ostream& ostr, const BaseList<space_T>& store) const;

    void set_unary(interaction_t,size_t ni, const T& v);
    const T& get_unary(interaction_t,size_t ni) const;
  };

  template<class T> std::ostream& operator<< (std::ostream& ostr, const HamiltonianStore<T>& a) {
    a.print(ostr);
    return ostr;
  }

} //namespace libcmatrix

#endif
