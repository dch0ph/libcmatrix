#ifndef _crystal_h_
#define _crystal_h_

#include "spinhalf_system.h"
#include "ListList.h"
#include "MetaPropagation.h"
#include "BlockedMatrix.h"
#include "BlockedMatrixTensor.h"
#include "geometry.h"
#include <iostream>

namespace libcmatrix {

const int INVALID=-10000;

//! if defined, use new vector_to_Euler function
#define LCM_USE_VECTOR_TO_EULER 1

class CrystalGeometry {
public:
  CrystalGeometry(size_t, const vector3&, const BaseList<vector3>&);// const Permutation& =Permutation());
  CrystalGeometry(const BaseList<size_t>&, const BaseList<vector3>&, const BaseList<vector3>&);//, const Permutation& =Permutation());

  vector3 operator()(size_t) const;
  vector3 operator()(size_t,size_t) const;

  void reverse(BaseList<size_t>, size_t&, size_t) const;
  void reverse(List<size_t>& d, size_t& k, size_t sk) const
  { d.create(dimensions()); reverse(static_cast< BaseList<size_t> >(d),k,sk); }

  void reverse(size_t& na, size_t& nb, size_t& nc,size_t& m,size_t sk) const { indexer_.reverse(na,nb,nc,m,sk); }
  size_t dimension(size_t n) const { return indexer_.dimension(n); }
  template<int N> size_t dimension(Int2Type<N> n) const { return indexer_.dimension(n); }
  size_t dimensions() const { return cell_vectors_.size(); }

  size_t nspins() const { return indexer_.size(); }
  size_t nspins_cell() const { return unit_cell_.size(); }
  void print(std::ostream&) const;

  size_t spin_to_spin(size_t sk) const { return sk % unit_cell_.size(); }
  size_t cell_to_spin(size_t nc,size_t m) const { return indexer_.index(0,0,nc,m); }
  size_t cell_to_spin(size_t nb,size_t nc,size_t m) const { return indexer_.index(0,nb,nc,m); }
  size_t cell_to_spin(size_t na,size_t nb,size_t nc,size_t m) const { return indexer_.index(na,nb,nc,m); }
  //  void permute(size_t&, size_t&, size_t, size_t) const;
  
  double dipolar_coupling(const basespin_system& sys,size_t j,size_t sk) const;
  double dipolar_coupling_unscaled(const basespin_system& sys,size_t j,size_t sk) const;
  space_T dipolar_tensor(const basespin_system& sys,size_t j,size_t sk) const;
  void dipolar_parameters(double& d, Euler& PAS,const basespin_system& sys,size_t j,size_t sk) const;

private:
  size_t ndims_;
  Indexer<4> indexer_;
  Indexer<3> cell_to_cell_;
  List<vector3> unit_cell_;
  List<vector3> cell_vectors_;
  //  Permutation perm_;

  //MUST MATCH corresponding functions in CrystalStructure
  template<int N> size_t inverse(size_t cell, Int2Type<N> n) const {
    return cell ? indexer_.dimension(n)-cell : 0;
  }
  //size_t inverse(size_t, size_t, size_t) const;
  
  void verify();
  template<int N> void mla(vector3& d,int n,Int2Type<N>) const;
};

  inline std::ostream& operator<< (std::ostream& ostr, const CrystalGeometry& a) {
    a.print(ostr);
    return ostr;
  }

  typedef unsigned int magicsize_t;

  class CrystalSymmetriseBase {
  public:
    CrystalSymmetriseBase(size_t N, size_t M, size_t neigs,int flags);

    virtual ~CrystalSymmetriseBase() {}
    virtual CrystalSymmetriseBase* clone() const =0;

    const ListList<state_t>& linkedstates() const {
      return linkedstates_; }
    size_t permutation_blocks() const { return linkedstates_.size(); }
    size_t actual_eigenvalues() const { return useeigs; }
    virtual void print(std::ostream&) const;
    bool isreal(size_t k) const { return isreal_(k); }
    size_t eigweight(size_t k) const {
      return (useeigsym_ && !isreal_(k)) ? 2 : 1;
    }
    virtual bool haspermutation() const =0;
    size_t eigenvalues() const { return neigs_; }
    size_t ncells() const { return ncells_; }

    virtual void addtoevals(BaseList<complex>&,size_t blki,size_t blkj,size_t r,size_t c) const =0;
    //virtual double addtoeval0(size_t rs,size_t cs) const =0;
    virtual complex symmetrise(Matrix<double>&, magicsize_t rs,magicsize_t cs,int k) const =0;
    virtual complex symmetrise(Matrix<complex>&, magicsize_t rs,magicsize_t cs,int k) const =0;
    const List<double>& scale_factors() const { return scalefacs; }
    magicsize_t blockshape(size_t blk) const { return blockshapes_(blk); }
    size_t blocksize(size_t blk) const { return linkedstates_.size(blk); }
    const List<size_t>& state_to_block() const { return state_to_block_; }
    const List<size_t>& state_to_index() const { return state_to_index_; }

    const BaseList<bool> which_eigenvalues(magicsize_t i) const { return haseig.row(i); }
    static void make_eigfacs(BaseList<complex>);
    size_t maxstates() const { return maxstates_; }

  protected:
    size_t neigs_;
    size_t ncells_;
    size_t total_;
    bool useeigsym_;
    size_t useeigs;
    size_t maxstates_;
    Matrix<bool> haseig;
    List<double> scalefacs;
    List<size_t> state_to_block_;
    List<size_t> state_to_index_;
    ListList<state_t> linkedstates_;   
    List<magicsize_t> blockshapes_;
    List<bool> isreal_;
    void makeindices();
    template<class F> void getsymmetry(F&);
  };

  class CrystalOpGenerator : public SpinOpGeneratorBase {
public:
  //Constructors
  CrystalOpGenerator(const spinhalf_system&, const CrystalStructure&, int flags_ =0, int verbosev =1); 
  CrystalOpGenerator(const spinhalf_system&, const CrystalStructure&, const BaseList<nuclei_spec>&, int flags_ =0, int verbosev =1);
  CrystalOpGenerator(const spinhalf_system&, const CrystalStructure&, const char* label, int flags_ =0, int verbosev =1);
    CrystalOpGenerator(const HamiltonianStructure&, const CrystalStructure&, int verbosev =1);

  SpinOpGeneratorBase* clone() const { return new CrystalOpGenerator(*this); }

  typedef ListList<double> diagonal_type;
  typedef BlockedMatrix<complex> output_type;
  typedef Matrix<complex> suboutput_type;
  typedef complex base_type; //operators are always complex
    static const bool allowproductop=false;
    
  const spin& operator()(size_t n) const {
    return (*sysp_)(n);
  }
  
    size_t size(size_t mzblk, size_t eigblk) const {
      return eigsizes(mzblk,eigblk);
    }
    size_t eigweight(size_t k) const {
      return symmetriserp_->eigweight(k);
    }
    
  const basespin_system& spinsystem() const 
    { return *sysp_; }
  
  template<class T,class C> static void add_Hquadrupolar_ns(T&, const C&, size_t) { 
    throw InternalError("add_Hquadrupolar_ns invalid for CrystalOpGenerator");
  }

    size_t index(size_t mzeig,size_t eig) const { return eigstr_.index(mzeig,eig); }
  
  friend std::ostream& operator<< (std::ostream&, const CrystalOpGenerator&);

  friend class CrystalSystem_iterator;
  friend class CrystalSystem_diagiterator;

  typedef BlockedMatrix<complex> operator_type;
  typedef ListList<double> diag_operator_type;

    const BaseList<size_t>& diagonal_structure() const { return diagstr_; }

  const ListList<double>& diag_Iz(size_t m) const
    { return tzops_(m); }
    const ListList<double>& diag_Fz(nuclei_spec nuc) const 
    { return fzcache_(*this,nuc()); }
  void rawdiag_Fz(ListList<double>&, size_t) const;

  void mla_Iz(ListList<double>&, double scale, size_t sel) const;
  void mla_Fz(ListList<double>&, double scale, nuclei_spec) const;
    void mla_Iz(ListList<double>&, double scale, const operator_spec&) const;

  void mla(BlockedMatrix<complex>& dest, block_pattern& blkspec, double scale, const operator_spec&) const;
  
    void add(BlockedMatrix<complex>& dest, block_pattern& blkspec, const productoperator_spec& spec) const;

  template<typename T> void create(BlockedMatrix<T>& dest, const block_pattern&) const;

   template<class M> void create(M& dest) const
    { SpinOpGeneratorBase::create(dest); }

  void add_A2(BlockedMatrixTensor<complex>&, const Matrix<space_T>&) const;
  void add_A0(BlockedMatrixTensor<complex>&, const Matrix<double>&) const;
  template<class T> void add_A2(BlockedMatrix<complex>&, const Matrix<T>&) const;
  void add_A0(BlockedMatrix<complex>&, const Matrix<double>&) const;
    
    template<class T,class C> void add_ns(T&, const C&, size_t, size_t, ns_flag) const { throw Failed("CrystalOpGenerator: non-secular interactions not supported"); }
    
    template<class OutType, class StoreType> void add_Hquadrupolar(OutType&, const StoreType&, size_t) const
    { throw Failed("Qaudrupolar interactions not supported"); }
    
    void add_Hcs(BlockedMatrixTensor<complex>&, const space_T&, size_t) const;
    void add_Hcs(space_T&, const space_T&, size_t) const;

 //static versions
    void add_Hcs(cmatrix&,  double, size_t, int) const;
    void add_Hcs(BlockedMatrix<complex>&,  double, size_t) const;
    void add_Hcs(ListList<double>&, double, size_t) const;

  void add_Hcs(MatrixTensor<complex>&, const BaseList<space_T>&, int) const;
  void add_Hcs(BlockedMatrixTensor<complex>&, const BaseList<space_T>&) const;

    template<class T, typename Coup> void add_Hquadrupolar2(T&, const Coup&, size_t) const { throw Failed("Can't use second order quadrupoles with CrystalOpGenerator"); }

private:
    smartptr<const spinhalf_system,false> sysp_; //replicated (required to establish block structure)
    smartptr<CrystalSymmetriseBase> symmetriserp_; 
    ListList<size_t> blocklists_;
    List<size_t> diagstr_;
    List<size_t> mzsizes;
    Matrix<size_t> eigsizes;
    Matrix<int> eigptrs;
    block_pattern diag_blkspec_;

    mutable FzCache_ fzcache_;
    ListList<double> Hzeeman_dummy_;

    void make_tzop(ListList<double>& dest, size_t which) const { make_tzop(dest,BaseList<size_t>(1,&which)); }
    void make_tzop(ListList<double>&, const BaseList<size_t>&) const;
  DynamicList< ListList<double> > tzops_;

  //Return (previously calculated) mz structure for a given nucleus (sum Fz for nuc==NULL_NUCLEUS)
  //const BaseList< ListList<size_t> > mzblocks(size_t nuc) const;

  int mz2val(state_t) const;
  int mz2val(size_t nuc,state_t) const;

    void create(const BaseList<nuclei_spec>& =BaseList<nuclei_spec>());
    void makemzinds(Matrix<size_t>&, List<size_t>&, const BaseList<size_t>&);
    size_t defaultnucleus_;
};


const double DEADBEEF=1e30;

  template<> struct spinsystem_glue_<CrystalOpGenerator> {
    static const basespin_system& spinsystem(const CrystalOpGenerator& a)
    { return a.spinsystem(); }

    static size_t ncells(const CrystalOpGenerator& a)
    { return a.ncells(); }

    static const bool ispurereal=false;
    
    template<class M> static bool arematching(const CrystalOpGenerator& a, const HamiltonianStore<M>& store) {
      return (store.ncells()==a.ncells()) && (a.nspins()==(store.ncells()*store.nspins()));
    }
  };
  
inline int validmz(float_t mz) { 
  float_t fres=floor(2*mz+0.5);
  if (fabs(fres-2.0*mz)>1e-3)
    throw InvalidParameter("Bad mz quantum number"); 
  return int(fres); }

//Need to specialise as Hamiltonian is constructed differently
  template<> void BlockedSpinningHamiltonianImpl<complex,CrystalOpGenerator>::interactions(const HamiltonianStore<space_T>&);

  template<class F> void CrystalSymmetriseBase::getsymmetry(F& obj)
{
  List<state_t> state_list;
  state_list.reserve(maxstates_);
  ScratchList<bool> used(maxstates_,false);
  ScratchList<size_t> sizes(maxstates_);
  
  size_t blkcount=0;
  
  for (;;) {
    size_t start;
    for (start=0;start<maxstates_ && used(start);start++); //find 1st unused state
    if (start==maxstates_)
      break; //used them all?
    
    size_t curlength=state_list.size();
    blockshapes_.push_back(obj.push_states(state_list,used,start));
    const size_t madesize=state_list.size()-curlength;
    assert(madesize>0);
    sizes(blkcount++)=madesize;

    const state_t comple=maxstates_-start-1;
    if (!used(comple)) {
      curlength=state_list.size();
      blockshapes_.push_back(obj.push_states(state_list,used,comple));
      sizes(blkcount++)=state_list.size()-curlength;
    }
  }
  linkedstates_=ListList<state_t>(sizes.truncate(blkcount),state_list);
  //chops up linked state list by sizes of blocks
  makeindices();
}  

#ifdef LCM_ENABLE_GENERICPERIODIC

struct general_creator;

class GeneralSymmetrise : public CrystalSymmetriseBase {
public:
  GeneralSymmetrise(const CrystalStructure&, int M, int flags_, int verbose);

  complex symmetrise(Matrix<double>& store, magicsize_t si, magicsize_t sj,int k) const { return symmetrise_(store,si,sj,k); }
  complex symmetrise(Matrix<complex>& store, magicsize_t si,magicsize_t sj,int k) const { return symmetrise_(store,si,sj,k); }

  CrystalSymmetriseBase* clone() const { return new GeneralSymmetrise(*this); }
 
  void addtoevals(BaseList<complex>& evals, size_t blki, size_t blkj, size_t r, size_t c) const;
  bool haspermutation() const { return !perm.empty(); }

  friend struct general_creator;
private:
  template<typename T> complex symmetrise_(Matrix<T>&, magicsize_t, magicsize_t, int k) const;
  List<size_t> eigtranslate; //translate external eigenvalue into internal representation
  List<int> eigreverse;
  
  magicsize_t getshape(size_t n0,size_t n1,size_t n2) const 
  { return indexer(n0-1,n1-1,n2-1); }

  void addeigenvalue(size_t& index, magicsize_t, size_t inteig);

  Indexer<3> indexer;
  Permutation perm;
  int verbose;
  Matrix<int> eigindex;
  ListList<int> extraeigindex;
  List<cmatrix> transforms;
};

#endif

} //namespace libcmatrix

#endif
