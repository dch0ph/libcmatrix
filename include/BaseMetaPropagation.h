#ifndef BaseMetaPropagation_h_
#define BaseMetaPropagation_h_

#include "basespin_system.h"
#include "BasePropagation.h"
#include "MultiMatrix.h"
#include "BlockedMatrix.h"
#include "smartptr.h"
#include "MAS.h"
#include <map>
#include "HamiltonianStore.h"
#include "PartitionedMatrix.h"
//Base objects for meta-propagation (simple or periodic)

namespace libcmatrix {
    
  template<class T> struct Ham_traits;
  //properties specific to Hamiltonian objects
  
   template<class T> struct Propagator_traits {
     static const bool allowsetH=false; //default is not H (i.e. propagator)
     static const bool hassynchronisation=false; // has synchronisation hint (distinct from periodicity of H)
   };
  //properties specific to objects used as propagator/H sources


  struct FzCache_ {
    std::map<size_t,ListList<double> > cache_;

    template<class OpGen> const ListList<double>& operator()(const OpGen& opgen, size_t nuc) {
      ListList<double>& val(cache_[nuc]);
      if (val.empty()) {
	try {
	  opgen.rawdiag_Fz(val,nuc);
	  if (val.empty())
	    throw Failed("Fz: nucleus not present in spin system!");
	} catch(...) {
	  cache_.erase(nuc); //tidy up map if create failed
	  throw;
	}
      }
      return val;
    }
  };

  class SpinOpGeneratorBase;

  struct BaseStructure {
    size_t eigblocks() const { return indexer_.dimension(Int2Type<1>()); }
    size_t mzblocks() const { return mzblocks_; }
    size_t actual_mzblocks() const { return indexer_.dimension(Int2Type<0>()); }
    bool usemzsym() const { return usemzsym_; }
    size_t totalblocks() const { return indexer_.size(); }

    size_t mzblocks_;
    bool usemzsym_;
    Indexer<2> indexer_;

    BaseStructure() : mzblocks_(0) {}
    BaseStructure(size_t mzblocksv, size_t eigblocksv, bool usemzsymv =false) 
     { set(mzblocksv,eigblocksv,usemzsymv); }

    template<class T> explicit BaseStructure(const T&);
 
    bool operator== (const BaseStructure& a) const {
      return (indexer_==a.indexer_) && (usemzsym_==a.usemzsym_);
    }
    bool operator!= (const BaseStructure& a) const {
      return (indexer_!=a.indexer_) || (usemzsym_==a.usemzsym_);
    }

    void set(size_t mzblocksv, size_t eigblocksv, bool usemzsymv);

    const Indexer<2>& indexer() const { return indexer_; }
    size_t index(size_t mzeig, size_t eig) const
    { return indexer_(mzeig,eig); }
    
    void reverse(size_t& mzeig, size_t& eig, size_t n) const
    { indexer_.reverse(mzeig,eig,n); }
  };

  inline std::ostream& operator<< (std::ostream& ostr, const BaseStructure& a)
    {
      return ostr << "mz blocks: " << a.mzblocks() << "   eigenvalue blocks: " << a.eigblocks() << '\n';
    }

  struct operator_spec;

  struct blocking_information {
    blocking_information(const BaseList<size_t>& a, const BaseList<size_t>& index, bool usemz =false) { create(a,index,usemz); }
    blocking_information() : totallevels(0) {}
    
    ScratchList<size_t> levels;
    List<size_t> steps;
    bool usemz;
    size_t totallevels;

    void create(const BaseList<size_t>&, const BaseList<size_t>&, bool);
    size_t types() const { return steps.size(); } 
    size_t majorblocks(size_t) const;
    size_t mzblocks(size_t which, size_t step =0) const;
  };

  std::ostream& operator<< (std::ostream&, const blocking_information&);

  class SpinOpGeneratorBase {
  public:
    explicit SpinOpGeneratorBase(const CrystalStructure& cstructv, const HamiltonianStructure& Hstruct, int verbosev =0)
      : cstruct(cstructv), structurep_(&Hstruct), flags_(Hstruct.flags()), verbose_(verbosev) { 
      if (cstruct.ncells()!=Hstruct.ncells())
	throw Mismatch("SpinOpGeneratorBase");
      init_(Hstruct.nspins_cell());
    }
    explicit SpinOpGeneratorBase(const CrystalStructure& cstructv, size_t Mv, int flagsv, int verbosev =0)
      : cstruct(cstructv), structurep_(NULL), flags_(flagsv), verbose_(verbosev) { init_(Mv); }

    virtual ~SpinOpGeneratorBase() {}
    
    virtual SpinOpGeneratorBase* clone() const =0;
    const SpinOpGeneratorBase& generator() const { return *this; } //in case BaseStructure initialised directly from OpGen

    virtual size_t index(size_t,size_t) const =0;
    virtual size_t size(size_t, size_t) const =0;
    virtual size_t eigweight(size_t) const =0;

    size_t mzblocks() const { return eigstr_.mzblocks(); }
    size_t eigblocks() const { return eigstr_.eigblocks(); }
    size_t actual_mzblocks() const { return eigstr_.actual_mzblocks(); }
    bool usemzsym() const { return eigstr_.usemzsym(); }
    size_t totalblocks() const { return actual_mzblocks()*eigblocks(); }
    const BaseStructure& structure() const { return eigstr_; }
    
    const blocking_information& active_blocking_info() const { return active_blocking_info_; }
    const blocking_information& inactive_blocking_info() const { return inactive_blocking_info_; }

    const ListList<double>& Hzeeman() const { return Hzeeman_; }
    const BaseList<double> Hzeeman(size_t mzblk,size_t eigblk) const {
      if (Hzeeman_.empty())
	throw Undefined("Hzeeman");
      return Hzeeman_(index(mzblk,eigblk));
    }
    const BaseList< ListList<size_t> > Hzeeman_structure() const { return zeemanstr_; }
    
    const ListList<size_t>& Hzeeman_structure(size_t mzblk,size_t eigblk) const {
      if (zeemanstr_.empty())
	throw Undefined("Hzeeman");      
      return zeemanstr_(index(mzblk,eigblk));
    }

    virtual const basespin_system& spinsystem() const =0;

    template<class M> void create(M& dest) const { 
      dest.create(diagonal_structure());
    }

    int flags() const { return flags_; }
    int verbose() const { return verbose_; }

    size_t cell_to_spin(size_t n, size_t m) const 
    { return cell_to_spin_(n,m); }
    size_t spin_to_spin(size_t sk) const { return sk % M_; }
    size_t ncells() const { return N_; }
    size_t nspins_cell() const { return M_; }
    size_t nspins() const { return cell_to_spin_.size(); }

    size_t quadrupole_order() const { 
      return (structurep_) ? structurep_->quadrupole_order() : 1;
    }
    ns_flag nonsecular_type(size_t& i, size_t& j) const {
      return (structurep_) ? structurep_->nonsecular_type(i,j) : NS_NONE;
    }
    bool nonsecular(size_t i) const {
      return (structurep_) ? structurep_->nonsecular(i) : false; 
    }
    bool isweak(interaction_t id) const {
      return (structurep_) ? structurep_->isweak(id) : false; 
    }

    bool isblocked() const { return (mzblocking_.size()>1) || !(partitioned_diagstr_.empty()); }
    bool isblocked(size_t nuc) const
    { return ispresent_.empty() ? false : ispresent_(nuc); }

    size_t nuclei() const { return blocknuc_.size(); }
    size_t nuctoindex(size_t nuc) const;
    size_t nuctoindex(const operator_spec&) const;

    size_t index(const BaseList<size_t>&) const;

    virtual const BaseList<size_t>& diagonal_structure() const =0;

    const ListList<size_t>& partitioned_diagonal_structure() const { return partitioned_diagstr_; }
    const List< ListList<size_t> >& mzblocking() const { return mzblocking_; }
    const BaseList<size_t> blockindices(size_t mzblk) const { return mzblocking_(mzblk).row(); }

    size_t size() const { return size_; }
        
    virtual void print(std::ostream&) const;
    
    double proton_frequency() const { return pfreq_; }
    virtual void proton_frequency(double pfreqv) { pfreq_=pfreqv; }

    const diagonal_partition_set& partitioning() const { return partitions_; }
    const diagonal_partition* partition(size_t k) const { 
      return partitions_.empty() ? (const diagonal_partition*)NULL : &(partitions_(k)); }

    bool isclassicsecondorder() const { return (flags() & MetaFlags::ClassicSecondOrder) && (quadrupole_order()!=1); }
    size_t effective_quadrupole_order() const { return (flags() & MetaFlags::ClassicSecondOrder) ? 1 : quadrupole_order(); } //!< return 1 if classic second order to pass as a normal Hamiltonian
    bool issecular() const { return (effective_quadrupole_order()==1); } //!< \return \c true if secular approximation is being used

  protected:
    const CrystalStructure cstruct;
    const HamiltonianStructure* structurep_;
    double pfreq_;
    size_t N_,M_;
    Indexer<2> cell_to_spin_;
    
    const int flags_;
    int verbose_;

    DynamicList<size_t> blocknuc_;
    List<bool> ispresent_;

    List< ListList<size_t> > mzblocking_;
    size_t submzblocks_;
    ListList<size_t> partitioned_diagstr_;
    size_t size_;
    DynamicList<size_t> sizes_;

    blocking_information active_blocking_info_,inactive_blocking_info_;

    BaseStructure eigstr_;
    ListList<double> Hzeeman_;
    List< ListList<size_t> > zeemanstr_;
    List< ListList<double> > secondordershifts1_,secondordershifts2_;

    List<int> nuctoindex_;
    List<size_t> indextonuc;
    Matrix<size_t> allmzinds_;
    List<size_t> maxmzinds_;

    diagonal_partition_set partitions_;

    double larmor(size_t) const;
    void init_(size_t Mv);
    void init_spinsys();
    void init_blockstr(size_t neigs, const BaseList<nuclei_spec>& =BaseList<nuclei_spec>());
    void init_partitioning(); //!< construct diagonal_partition_set from previously initialised partitioned_diag_str
    virtual void makemzinds(Matrix<size_t>&, List<size_t>&, const BaseList<size_t>&) =0;
 };


  template<class T> BaseStructure::BaseStructure(const T& source)
  {
    const SpinOpGeneratorBase& opgen(source.generator());
    set(opgen.mzblocks(),opgen.eigblocks(),opgen.usemzsym());
  }

  //! struct defining simple operator term
  struct operator_spec {
    size_t nuc; //!< nucleus id (\c NULL_NUCLEUS if not specified)
    size_t number; //!< spin index (-1 if not specified)
    char op; //!< operator or row index if single transition)
    char col; //!< single transition column ('x' if not specified)
    
    operator_spec(nuclei_spec whichn, char op_)
      : nuc(whichn()), number(-1), op(op_), col('x') {} //!< initialise F-type operator
    
    operator_spec(size_t number_, char op_)
      : nuc(NULL_NUCLEUS), number(number_), op(op_), col('x') {} //!< initialise I-type operator

    //! initialise I-type operator
    operator_spec(int number_, char op_)
      : nuc(NULL_NUCLEUS), number(number_), op(op_), col('x') {
      if (number_<0)
	throw InvalidParameter("operator_spec");
    }

    operator_spec(nuclei_spec whichn, size_t r, size_t c)
      : nuc(whichn()), number(-1) { initST(r,c); } //!< initialise F single transition operator

    operator_spec(size_t number_, size_t r, size_t c)
      : nuc(NULL_NUCLEUS), number(number_) { initST(r,c); } //!< initialise I single transition operator

    operator_spec(int number_, size_t r, size_t c)
      : nuc(NULL_NUCLEUS), number(number_) { 
      if (number_<0)
	throw InvalidParameter("operator_spec");
      initST(r,c);
    } //!< initialise I single transition operator

    //! return nucleus involved
    size_t nucleus(const basespin_system& sys) const {
      return (nuc==NULL_NUCLEUS) ? sys(number).nucleus() : nuc;
    }

    int coherence() const; //!< return coherence order
    bool issingletransition() const { return isvalidST(op); } //!< return \c true if holding a single transition operator
    bool issumoperator() const { return (nuc!=NULL_NUCLEUS); }
    bool isdiagonal() const { return (coherence()==0); }
    
    //! return \c true if parameter is valid for index (rather than operator symbol)
    /** \note We are making the (safe) assumption that all operator symbols have character codes of >31 */
    static bool isvalidST(char n) { return (n<32); }
    void initST(size_t r, size_t c); //!< \internal validate and initialise single transition operator
  };

  inline bool operator== (const operator_spec& a, const operator_spec& b)
  {
    return (a.nuc==b.nuc) && (a.number==b.number) && (a.op==b.op) && (a.col==b.col);
  }

  inline bool operator!= (const operator_spec& a, const operator_spec& b)
  {
    return (a.nuc!=b.nuc) || (a.number!=b.number) || (a.op!=b.op) || (a.col!=b.col);
  }

  inline bool arematching(const operator_spec& a, const operator_spec& b)
  {
    return (a.nuc==b.nuc) && (a.number==b.number) && operator_intersection(a.op,b.op) && (a.col==b.col);
  }
  
  operator_spec operator& (const operator_spec&, const operator_spec&);

  std::ostream& operator<< (std::ostream&, const operator_spec&);


  struct productoperator_spec {

    productoperator_spec() {}
    productoperator_spec(const operator_spec&, const basespin_system&);
    productoperator_spec(const operator_spec&); 
        
    bool operator!() const { return scales.empty(); }
    bool issimple() const { return (specs.items()==1); }
    const BaseList<operator_spec> front() const { return specs.front(); }
    
    size_t nucleus() const {
      if (scales.empty())
	throw Undefined("productoperator_spec");
      return nuc;
    }
    
    void push_back(complex, const BaseList<operator_spec>&, const basespin_system&);
    void push_back(complex, const operator_spec&);
    void push_back(complex scale, const operator_spec& spec, const basespin_system& sys)
    { push_back(scale,BaseList<operator_spec>(1,const_cast<operator_spec*>(&spec)),sys); }
    
    size_t size() const { return scales.size(); }
    const operator_spec& unique() const {
      if (specs.items()!=1)
	throw Failed("productoperator_spec::unique: not simple operator");
      return specs(0U,0U);
    }
    bool isdiagonal() const;
            
    productoperator_spec& operator&= (const productoperator_spec&);
    bool ismatching(const productoperator_spec&) const;
    void print(std::ostream&) const;
    void clear() { scales.clear(); specs.clear(); }

    void swap(productoperator_spec& a) {
      scales.swap(a.scales);
      specs.swap(a.specs);
      std::swap(nuc,a.nuc);
    }

    List<complex> scales;
    ListList<operator_spec> specs;
    size_t nuc;
  };
  
  bool operator== (const productoperator_spec&, const productoperator_spec&);
  bool operator!= (const productoperator_spec&, const productoperator_spec&);

  inline bool arematching(const productoperator_spec& a, const productoperator_spec& b)
  { return a.ismatching(b); }
  
  inline std::ostream& operator<< (std::ostream& ostr, const productoperator_spec& a)
  { a.print(ostr); return ostr; }
    
 
  struct block_pattern {
    block_pattern() { clear(); }
    block_pattern(const blocking_information& info) { create(info); }
    block_pattern(const blocking_information& infov, size_t whichv, int coherv, bool lishermv) { create(infov,whichv,coherv,lishermv); }
    block_pattern(const SpinOpGeneratorBase& opgen, bool active =true) { create(active ? opgen.active_blocking_info() : opgen.inactive_blocking_info()); }
    block_pattern(const SpinOpGeneratorBase& opgen, const operator_spec& opspec, bool active =true) { create(opgen,opspec,active); }
  block_pattern(const SpinOpGeneratorBase&, const productoperator_spec&);

    void clear();

    int coherence() const { return coher; }
    bool isdiagonal() const { return (coher==0); }
    bool isupper() const { return (coher>=0); }
    bool ishermitian() const { return isherm; }

    void transpose() { 
      if (!isherm) 
 	coher=-coher;
    }

    void makepresent(Matrix<bool>&) const;

    bool usemzsym() const { return usemz; }
    bool operator!() const { return (blocks==0); }
    block_pattern& operator&= (const block_pattern&);

    bool ismatching(const block_pattern& a) const {
      return (indexer==a.indexer) && (blocks==a.blocks); 
    }

    bool operator== (const block_pattern& b) const {
      return (coher==b.coher) && ismatching(b) && (isherm==b.isherm);
    }

    bool operator!= (const block_pattern& b) const {
      return (coher!=b.coher) || !ismatching(b) || (isherm!=b.isherm);
    }

    struct iterator;

    void create(const SpinOpGeneratorBase&, const operator_spec&, bool =true);
    //    void create();
    void create(const blocking_information&);
    void create(const blocking_information& infov, size_t which, int coher, bool isherm);

    enum { MZBLOCK=0, REPEAT, MAJOR };
    
    int coher;
    size_t blocks;
    bool isherm;
    bool usemz;
    size_t step;
    size_t whichstep;
    size_t majorstep;
    size_t maxind;
    Indexer<3> indexer;
  };

  std::ostream& operator<< (std::ostream&, const block_pattern&);

  struct block_pattern::iterator {
    iterator(const block_pattern&);

    void reset();
    bool next(size_t&, size_t&, size_t&); //!< changed Mar 13 to include stepping of mzeig
    bool next(size_t&, size_t&, size_t&, bool&);
    bool next(size_t&, size_t&, size_t&, double&);
    size_t getindex() const;

    const block_pattern& source;
    bool isbelow;
    Indexer<3>::permuted_iterator iter;
    const Indexer<3>::permuted_iterator end;
    int last;
    size_t lasttick;
    size_t mzeigSD;
  };

  // Set Ham and Operator traits for unblocked types

  template<> struct Ham_traits<RealSpinningHamiltonian> {
    typedef space_T coupling_type;
    static const bool allowRF=true;
    static const bool isconstant=false;
    static const bool isdiagonal=false;
    static const bool ishomogeneous=true;
  };

  template<> struct Ham_traits< List<RealSpinningHamiltonian> > {
    typedef space_T coupling_type;
    static const bool allowRF=true;
    static const bool isconstant=false;
    static const bool isdiagonal=false;
    static const bool ishomogeneous=true;
  };

  template<> struct Ham_traits<DiagonalSpinningHamiltonian> {
    typedef space_T coupling_type;
    static const bool allowRF=false;
    static const bool isconstant=false;
    static const bool isdiagonal=true;
    static const bool ishomogeneous=false;
  };

  template<> struct Ham_traits< List<DiagonalSpinningHamiltonian> > {
    typedef space_T coupling_type;
    static const bool allowRF=false;
    static const bool isconstant=false;
    static const bool isdiagonal=true;
    static const bool ishomogeneous=false;
  };

  template<> struct Ham_traits<SpinningHamiltonian> {
    typedef space_T coupling_type;
    static const bool allowRF=true;
    static const bool isconstant=false;
    static const bool isdiagonal=false;
    static const bool ishomogeneous=true;
  };

  template<> struct Ham_traits< List<SpinningHamiltonian> > {
    typedef space_T coupling_type;
    static const bool allowRF=true;
    static const bool isconstant=false;
    static const bool isdiagonal=false;
    static const bool ishomogeneous=true;
  };

   template<> struct Ham_traits< BlockedMatrix<complex> > { 
     typedef double coupling_type; 
     static const bool allowRF=true;
     static const bool isconstant=true;
     static const bool isdiagonal=false;
     static const bool ishomogeneous=false;
   }; 

   template<> struct Ham_traits<cmatrix> { 
     typedef double coupling_type; 
     static const bool allowRF=true;
     static const bool isconstant=true;
     static const bool isdiagonal=false;
     static const bool ishomogeneous=false;
   }; 

   template<> struct Ham_traits< BlockedMatrix<double> > { 
     typedef double coupling_type; 
     static const bool allowRF=false;
     static const bool isconstant=true;
     static const bool isdiagonal=false;
     static const bool ishomogeneous=false;
   }; 
 
  //Need this "basic" ones to deduce nature of propagator
   template<> struct Ham_traits<rmatrix> { 
     typedef double coupling_type; 
     static const bool allowRF=false;
     static const bool isconstant=true;
     static const bool isdiagonal=false;
     static const bool ishomogeneous=false;
   }; 
 
   template<> struct Ham_traits< BaseList<double> > { 
     typedef double coupling_type; 
     static const bool allowRF=false;
     static const bool isconstant=true;
     static const bool isdiagonal=true;
     static const bool ishomogeneous=false;
   }; 
 
   template<> struct Ham_traits< ListList<double> > { 
     typedef double coupling_type; 
     static const bool allowRF=false;
     static const bool isconstant=true;
     static const bool isdiagonal=true;
     static const bool ishomogeneous=false;
   }; 

#define LCM_INSTANTIATE_SPINOPGEN_CREATE(X)\
template void X::create(BlockedMatrix<double>&, const block_pattern&) const;\
template void X::create(BlockedMatrix<complex>&, const block_pattern&) const;\
template void X::create(BlockedMatrix<bool>&, const block_pattern&) const;
  
} //namespace libcmatrix

#endif

//   class block_pattern {
//   public:
//     block_pattern() { clear(); }
//     block_pattern(const SpinOpGeneratorBase& opgen, bool active =true) { create(opgen,active); }
//     block_pattern(const SpinOpGeneratorBase& opgen, const operator_spec& opspec)
//     { create(opgen,opspec); }
//     block_pattern(const SpinOpGeneratorBase&, const productoperator_spec&);
//     block_pattern(const SpinOpGeneratorBase& opgen, size_t nuc, int lcoher, bool lisherm, bool active =true) { create(opgen,nuc,lcoher,lisherm,active); }

//     bool isdiagonal() const { return (coher==0); }
//     bool isupper() const { return (coher>=0); }

//     void clear() {
//       coher=0;
//       mzstep=0;
//       blocks=0;
//     }

//     void transpose() { 
//       if (!isherm) 
//  	coher=-coher;
//     }

//     bool operator!() const { return (mzstep==0); }
//     block_pattern& operator&= (const block_pattern&);

//     bool ismatching(const block_pattern& a) const {
//       return (std::abs(coher)==std::abs(a.coher)) && (mzstep==a.mzstep) && (blocks==a.blocks);
//     }

//     void makepresent(Matrix<bool>&) const;

//   struct iterator {
//     iterator() : finished(true) {}
//     iterator(const block_pattern&);
//     bool next(size_t&, size_t&);
//     bool next(size_t&, size_t&, bool&);
//     void reset_block();

//     size_t crind,ccind;
//     size_t superbase;
//     int offset;
//     size_t majorincr;
//     size_t mzstep;
//     size_t smallblocks;
//     size_t smallblockcount;
//     size_t majorcount;
//     size_t curadd;
//     size_t blocks;
//     size_t maxind;
//     bool finished;    
//   };

//     int coher;
//     size_t blocks;
//     bool isherm;
//     bool hasmiddle;
//     size_t mzstep;

//     bool operator== (const block_pattern& b) const {
//       return (coher==b.coher) && (mzstep==b.mzstep) && (isherm==b.isherm);
//     }

//     bool operator!= (const block_pattern& b) const {
//       return (coher!=b.coher) || (mzstep!=b.mzstep) || (isherm!=b.isherm);
//     }

//     size_t mzlevels;
//   private:
//     size_t totlevels;
//     size_t maxind;

//     void create(const SpinOpGeneratorBase&, bool =true);
//     void create(const SpinOpGeneratorBase&, const operator_spec&);
//     void create(const SpinOpGeneratorBase&, size_t nuc, int coher, bool isherm, bool =true);
//   };
   
//   // block_pattern operator& (const block_pattern&, const block_pattern& b);
