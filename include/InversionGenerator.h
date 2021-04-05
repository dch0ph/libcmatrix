#ifndef LCM_InversionGenerator_h_
#define LCM_InversionGenerator_h_

#include "HamiltonianStore.h"
#include "lcm_CommonSequence.hpp"

namespace libcmatrix {

typedef std::pair<size_t,size_t> blkind_t;

template<typename T, typename IndexT =size_t> struct SparseEntry {
  IndexT index;
  T value;
  SparseEntry(const IndexT& rv, const T& vv)
    : index(rv), value(vv) {}
  //template<typename T2> SparseEntry& operator*= (const T2& v) { value*=v; return *this; }
  const T& operator()() const { return value; }
  T& operator()() { return value; }
};

//template<typename T, typename IndexT> class UniqueRowMatrix;
  
template<typename T, typename IndexT =size_t> class UniqueColumnMatrix
  : public DynamicList< SparseEntry<T,IndexT> > {
public:
  bool operator! () const { return this->empty(); }
  //  void conj_transpose(UniqueRowMatrix&) const;
};

template<typename T, typename IndexT =size_t> class UniqueRowMatrix
  : public DynamicList< SparseEntry<T,IndexT> > {
public:
  bool operator! () const { return this->empty(); }
  //  void conj_transpose(UniqueColumnMatrix&) const;
};

// template<typename T, typename IndexT> void UniqueColumnMatrix<T,IndexT>::conj_transpose(UniqueRowMatrix<T,IndexT>& dest) const
// {
//   typedef typename UniqueColumnMatrix<T,IndexT>::value_type value_type;
//   dest.create(this->size(),value_type(-1,-1));
//   for (size_t i=this->size();i--;) {
//     const value_type& source((*this)(i));
//     dest(i)=value_type(source.index,conj(source.value));
//   }
// }

template<typename T, typename IndexT> std::ostream& operator<< (std::ostream& ostr, const SparseEntry<T,IndexT>& a)
{ return ostr << a.index << ": " << a.value; }

template<typename T, typename IndexT> std::ostream& operator<< (std::ostream& ostr, const UniqueColumnMatrix<T,IndexT>& a)
{ 
  for (size_t i=0;i<a.size();i++)
    ostr << 'R' << a(i) << ' ';
  return ostr;
}

// template<typename T, typename IndexT> std::ostream& operator<< (std::ostream& ostr, const UniqueRowMatrix<T,IndexT>& a)
// { 
//   for (size_t i=0;i<a.size();i++)
//     ostr << 'C' << a(i) << ' ';
//   return ostr;
// }

class BlockedOperator;

class InversionGenerator {
public:
  template<class OpGen> InversionGenerator(const OpGen&, const BaseList<nuclei_spec>&, nuclei_spec);
  template<class OpGen> InversionGenerator(const OpGen&, const HamiltonianStructure&, const CrystalStructure&, nuclei_spec);
  void operator()(BlockedOperator&, const BlockedOperator&, double) const;
  static Warning<> nucleusnotfound_warning;

private:
  const SpinOpGeneratorBase& opgen;
  smartptr<const SpinOpGeneratorBase> useopgenp;
  const block_pattern nucblockpattern;
  BaseList< ListList<size_t> > toblocking;
  size_t nstates;
  List<blkind_t> destindex;
  List<blkind_t> fromindex;
  typedef UniqueColumnMatrix<complex> sparse_t;
  sparse_t elements;
  mutable sparse_t elementsz;
  mutable smartptr<ZshiftCache,false> zshiftp;
  mutable ListList<complex> zfacs;  
  smartptr<const HamiltonianStructure,false> Hstructp;

  void initialise(const PulseGenerator&);
  static void blockinglist(List<nuclei_spec>&, const BaseList<nuclei_spec>&, nuclei_spec);
  static void makeindex(BaseList<blkind_t>, const BaseList< ListList<size_t> >&);
  void unitary_simtrans_ip(sparse_t&, const ListList<complex>&) const;
  friend std::ostream& operator<< (std::ostream&, const InversionGenerator&);
};

template<class OpGen> InversionGenerator::InversionGenerator(const OpGen& opgenv, const BaseList<nuclei_spec>& blocking, nuclei_spec nucspec)
  : opgen(opgenv), nucblockpattern(opgenv,operator_spec(nucspec,'+')),
      toblocking(opgenv.mzblocking())
{
  List<nuclei_spec> useblocking;
  blockinglist(useblocking,blocking,nucspec);
  const int verbose=opgen.verbose();
  const OpGen* tuseopgenp=new OpGen(opgen.spinsystem(),useblocking,verbose,opgen.flags());
  useopgenp.reset(tuseopgenp);
  
  PulseGenerator pgen(*tuseopgenp,nucspec,verbose);
  initialise(pgen);
}

template<class OpGen> InversionGenerator::InversionGenerator(const OpGen& opgenv, const HamiltonianStructure& Hstruct, const CrystalStructure& cstruct, nuclei_spec nucspec)
  : opgen(opgenv), nucblockpattern(opgenv,operator_spec(nucspec,'+')),
      toblocking(opgenv.mzblocking())
{
  List<nuclei_spec> useblocking;
  blockinglist(useblocking,Hstruct.blockingnuclei(),nucspec);
  const int verbose=opgen.verbose();
  HamiltonianStructure* Hp=new HamiltonianStructure(Hstruct);
  Hp->blockingnuclei(useblocking);
  Hstructp.reset(Hp);
  const OpGen* tuseopgenp=new OpGen(*Hp,cstruct,verbose);
  useopgenp.reset(tuseopgenp);
  
  PulseGenerator pgen(*tuseopgenp,nucspec,verbose);
  initialise(pgen);
}

} //namespace libcmatrix

#endif
