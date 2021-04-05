//Only necessary if defining new MetaPropagation methods

#ifndef lcm_MetaPropagation_h_
#define lcm_MetaPropagation_h_

#include "MetaPropagation.h"

namespace libcmatrix {
  
  const double lcm_degeneracy_factor=1e-8;

  template<class OutType,class InType> struct reducer_ {};
  template<> struct reducer_<double,double> {
    double operator()(double v) const { return v; }
    double cs(double v) const { return v; }
  };
  template<> struct reducer_<space_T,space_T> {
    const space_T& operator()(const space_T& v) const { return v; }
    const space_T& cs(const space_T& v) const { return v; }
  };
  template<> struct reducer_<double,space_T> {
    double operator()(const space_T& A) const {
      return A.have_rank(2) ? real(A(2,0)) : 0.0;
    }
    double cs(const space_T& A) const {
      double res=0.0;
      if (A.have_rank(0))
	res+=real(A(0,0));
      if (A.have_rank(2))
	res+=real(A(2,0));
      return res;
    }
  };

  template<class OutType,class OpGen> void add_ns_(OutType& dest, const OpGen& opgen, const space_T& A, size_t j, size_t k, ns_flag nstype)
  {
    opgen.add_ns(dest,A,j,k,nstype);
  }


  template<class OutType,class OpGen> void add_ns_(OutType&, const OpGen&, double, size_t, size_t, ns_flag)
  {
    throw Failed("Can't specify non-secular interaction with simple coupling");
  }

template<class OutType, class OpGen, class StoreType> void 
add_Hamiltonian(OutType& dest, const OpGen& opgen, interaction_t id, const couplingstore<StoreType>& coups)
{
  if (!coups)
    return;

  const size_t cells=opgen.ncells();
  const size_t M=coups.rows();
  if (M*cells!=coups.cols())
    throw Mismatch("add_Hamiltonian (binary)");

  const bool forceweak=opgen.isweak(id);

  const Matrix<double> storeA0(coups.A0);
  const bool haveiso=!!storeA0;
  const Matrix<StoreType> storeA2(coups.A2);
  reducer_<typename Ham_traits<OutType>::coupling_type,StoreType> reduce;

  for (size_t i=0;i<M;i++) {
    for (size_t j=0;j<M;j++) {
      for (size_t offset=0;offset<cells;offset++) {
	if ((i>=j) && (offset==0))
	  continue; //quick check for i>=j within cell
	const size_t sk=opgen.cell_to_spin(offset,j);
	if (haveiso) {
	  const double iso=storeA0(i,sk);
	  if (iso)
	    opgen.add_A0(dest,iso,i,sk,forceweak);
	}
	const StoreType& aniso(storeA2(i,sk));
	if (isnonnull(aniso)) {
	  size_t li(i);
	  size_t lk(sk);
	  const ns_flag nstype=opgen.nonsecular_type(li,lk);
	  if (nstype!=NS_NONE)
	    add_ns_(dest,opgen,aniso,li,lk,nstype);
	  else
	    opgen.add_A2(dest,reduce(aniso),li,lk,forceweak);
	}
      }
    }
  }
}

template<class OutType, class OpGen, class StoreType> void   
  add_Hamiltonian(OutType& dest, const OpGen& opgen, interaction_t id, const List<StoreType>& coups)
{
  if (coups.empty())
    return;

  reducer_< typename Ham_traits<OutType>::coupling_type, StoreType> reduce;

  for (size_t i=coups.size();i--;) {
    const StoreType& coup=coups(i);
    if (isnonnull(coup)) {
      if (id==I_QUAD) {
	if (opgen.nonsecular(i))
	  opgen.add_ns(dest,coup,i,i,NS_BOTH);
	else {
	  if (opgen.isclassicsecondorder())
	    opgen.add_Hquadrupolar2(dest,coup,i); 
	  else
	    opgen.add_Hquadrupolar(dest,reduce(coup),i);
	}
      }
      else
	opgen.add_Hcs(dest,reduce.cs(coup),i);
    }
  } 
}
 
template<class OutType, class OpGen, class StoreType> void add_coupling_Hamiltonians(OutType& H, const OpGen& opgen, const HamiltonianStore<StoreType>& Hstore) 
{
  typedef typename HamiltonianStore<StoreType>::couplingmap_type::const_iterator iter_type;
  const iter_type iterend(Hstore.couplings_end());
  iter_type iter(Hstore.couplings_begin());
  while (iter!=iterend) {
    add_Hamiltonian(H,opgen,iter->first,iter->second);
    ++iter;
  }
}

template<class OutType, class OpGen, class StoreType> void add_Hamiltonian(OutType& H, const OpGen& opgen, const HamiltonianStore<StoreType>& Hstore) 
{
  if (!Hstore.ismatching(opgen))
    throw Mismatch("add_Hamiltonian: Operator generator and HamiltonianStore don't match");

  add_coupling_Hamiltonians(H,opgen,Hstore);

  typedef typename HamiltonianStore<StoreType>::shiftmap_type::const_iterator iter_type;
  const iter_type iterend(Hstore.shifts_end());
  iter_type iter(Hstore.shifts_begin());
  while (iter!=iterend) {
    add_Hamiltonian(H,opgen,iter->first,iter->second);
    ++iter;
  }

  add_Hamiltonian(H,opgen,I_QUAD,Hstore.get_quadrupole());
}
  
template<class BaseType,class OpGen> void BlockedStaticHamiltonianImpl<BaseType,OpGen>::interactions(const HamiltonianStore<double>& store)
{
  this->reset();
  BlockedMatrix<BaseType>& H(*this);
  add_Hamiltonian(H,opgen_,store);
}

template<class BaseType,class OpGen> void BlockedStaticHamiltonianImpl<BaseType,OpGen>::interactions(const HamiltonianStore<space_T>& store)
{
  this->reset();
  BlockedMatrix<BaseType>& H(*this);
  add_Hamiltonian(H,opgen_,store);
}

template<class OpGen> void BlockedDiagonalStaticHamiltonianImpl<OpGen>::interactions(const HamiltonianStore<double>& store)
{
  reset();
  ListList<double>& H(*this);
  add_Hamiltonian(H,opgen_,store);
}

template<class OpGen> void BlockedDiagonalStaticHamiltonianImpl<OpGen>::interactions(const HamiltonianStore<space_T>& store)
{
  ListList<double>& H(*this);
  if ((opgen_.quadrupole_order()==2) && !(opgen_.flags() & MetaFlags::ClassicSecondOrder)) {
    if (opgen_.eigblocks()!=1)
      throw Failed("Non-secular Hamiltonians not compatible with eigenvalue structure");
    BlockedMatrix<complex> Htmp;
    add_Hamiltonian(Htmp,opgen_,store);
    opgen_.create(H);
    for (size_t n=H.size();n--;)
      hermitian_eigenvalues_second(H(n),Htmp(n),opgen_.Hzeeman(n,0U));
  }
  else {
    reset();
    add_Hamiltonian(H,opgen_,store);
  }
}

template<class BaseType> template<class OpGen> BlockedStaticHamiltonian<BaseType>::BlockedStaticHamiltonian(const OpGen& opgen, const HamiltonianStore<double>& store, const BlockedMatrix<BaseType>& Hbase)
  : BaseStructure(opgen), 
  pImpl_(new BlockedStaticHamiltonianImpl<BaseType,OpGen>(opgen,Hbase))
{
  interactions(store);
}

template<class BaseType> template<class OpGen> BlockedStaticHamiltonian<BaseType>::BlockedStaticHamiltonian(const OpGen& opgen, const BlockedMatrix<BaseType>& Hbase)
  : BaseStructure(opgen),
    pImpl_(new BlockedStaticHamiltonianImpl<BaseType,OpGen>(opgen,Hbase)) {}


template<class OpGen> void BlockedDiagonalSpinningHamiltonianImpl<OpGen>::interactions(const HamiltonianStore<space_T>& store)
{
  List<DiagonalSpinningHamiltonian>& H(*this);
  if ((opgen_.quadrupole_order()==2) && !(opgen_.flags() & MetaFlags::ClassicSecondOrder)) {
    if (opgen_.eigblocks()!=1)
      throw Failed("Non-secular Hamiltonians not compatible with eigenvalue structure");
    //FUTURE: If RealSpinningHamiltonian supports NS, then ...
    //    List<typename lcm_promote_spin_<typename OpGen::base_type>::spinning_type> Htmp;
    ensure();
    List<SpinningHamiltonian> Htmp(H.size(),SpinningHamiltonian(rotor_speed_,rotor_phase_,rinfo_));
    add_Hamiltonian(Htmp,opgen_,store);
    for (size_t n=H.size();n--;)
      H(n).interactions(Htmp(n),opgen_.Hzeeman(n,0U),opgen_.Hzeeman_structure(n,0U).row(),lcm_degeneracy_factor);
  }
  else {
    reset();
    add_Hamiltonian(H,opgen_,store);
  }
}

// template<class OpGen> void BlockedDiagonalSpinningHamiltonianImpl<OpGen>::interactions(const BlockedSpinningHamiltonian<complex>& Hsource)
// {
//   List<DiagonalSpinningHamiltonian>& H(*this);
//   ensure();
//   for (size_t n=H.size();n--;)
//     H(n).interactions(Hsource(n,0U),opgen_.Hzeeman(n,0U),lcm_degeneracy_factor);
// }

template<class OpGen> BlockedDiagonalStaticHamiltonian::BlockedDiagonalStaticHamiltonian(const OpGen& opgen, const HamiltonianStore<double>& store)
  : BaseStructure(opgen),
    pImpl_(new BlockedDiagonalStaticHamiltonianImpl<OpGen>(opgen))
{
  interactions(store);
}

template<class OpGen> BlockedDiagonalStaticHamiltonian::BlockedDiagonalStaticHamiltonian(const OpGen& opgen)
 : BaseStructure(opgen),
   pImpl_(new BlockedDiagonalStaticHamiltonianImpl<OpGen>(opgen)) {}

// template<class BaseType,class OpGen> void BlockedHamiltonian<BaseType,OpGen>::dump() const
//   {
//     std::cout << (*this);
//   }

template<class OpGen> BlockedDiagonalSpinningHamiltonian::BlockedDiagonalSpinningHamiltonian(const OpGen& opgen, const HamiltonianStore<space_T>& store, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo)
  : BaseStructure(opgen),
    pImpl_(new BlockedDiagonalSpinningHamiltonianImpl<OpGen>(opgen,rotor_speedv,rotor_phasev,rinfo))
{
  interactions(store);
}

template<class OpGen> BlockedDiagonalSpinningHamiltonian::BlockedDiagonalSpinningHamiltonian(const OpGen& opgen, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo)
  : BaseStructure(opgen),
    pImpl_(new BlockedDiagonalSpinningHamiltonianImpl<OpGen>(opgen,rotor_speedv,rotor_phasev,rinfo)) {}


//! WARNING - this code fails when instantiated
/* template<class OpGen> BlockedFilter::BlockedFilter(const OpGen& opgen,const BaseList<int>& cohers) */
/*   : BaseStructure(opgen) */
/* { */
/*   if (opgen.blockingnuc().size()) */
/*     throw Failed("Can't initialise BlockedFilter with blocked coherence"); */
/*   const ListList<double>& nucFz(opgen.diag_Fz(opgen.spinsystem().homonucleus())); */
/*   const int maxcoher(coherencematrix_range(nucFz.row(),cohers)); */
/*   const block_pattern diagpattern(opgen); */
/*   opgen.create(store_,diagpattern); */
/*   store_=true; */
/*   for (size_t n=nucFz.size();n--;) */
/*     coherencematrix_ip(store_(n),nucFz(n),cohers,maxcoher); */
/* } */

template<class BaseType> template<class OpGen> BlockedSpinningHamiltonian<BaseType>::BlockedSpinningHamiltonian(const OpGen& opgen, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const BlockedMatrix<BaseType>& Hbase)
  : BaseStructure(opgen),
    pImpl_(new BlockedSpinningHamiltonianImpl<BaseType,OpGen>(opgen,rotor_speedv,rotor_phasev,rinfo,Hbase)) {}

template<class BaseType> template<class OpGen> BlockedSpinningHamiltonian<BaseType>::BlockedSpinningHamiltonian(const OpGen& opgen, const HamiltonianStore<space_T>& Hstore, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const BlockedMatrix<BaseType>& Hbase)
  : BaseStructure(opgen),
    pImpl_(new BlockedSpinningHamiltonianImpl<BaseType,OpGen>(opgen,rotor_speedv,rotor_phasev,rinfo,Hbase))
{ pImpl_->interactions(Hstore); }

template<class BaseType> template<class OpGen> BlockedSpinningHamiltonian<BaseType>::BlockedSpinningHamiltonian(const OpGen& opgen, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const IntervalSamplerBase& samplerv, const BlockedMatrix<BaseType>& Hbase)
  : BaseStructure(opgen),
    pImpl_(new BlockedSpinningHamiltonianImpl<BaseType,OpGen>(opgen,rotor_speedv,rotor_phasev,rinfo,samplerv,Hbase)) {}

template<class BaseType> template<class OpGen> BlockedSpinningHamiltonian<BaseType>::BlockedSpinningHamiltonian(const OpGen& opgen, const HamiltonianStore<space_T>& Hstore, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const IntervalSamplerBase& samplerv, const BlockedMatrix<BaseType>& Hbase)
  : BaseStructure(opgen),
    pImpl_(new BlockedSpinningHamiltonianImpl<BaseType,OpGen>(opgen,rotor_speedv,rotor_phasev,rinfo,samplerv,Hbase))
{ pImpl_->interactions(Hstore); }

#define LCM_INSTANTIATE_TYPED_HAMILTONIANS(X,T)\
  template BlockedStaticHamiltonian<T>::BlockedStaticHamiltonian(const X&, const HamiltonianStore<double>&, const BlockedMatrix<T>&); \
  template BlockedStaticHamiltonian<T>::BlockedStaticHamiltonian(const X&, const BlockedMatrix<T>&); \
  template BlockedSpinningHamiltonian<T>::BlockedSpinningHamiltonian(const X&, const HamiltonianStore<space_T>&, double, double, const RotorInfo&, const BlockedMatrix<T>&); \
  template BlockedSpinningHamiltonian<T>::BlockedSpinningHamiltonian(const X&, double, double, const RotorInfo&, const BlockedMatrix<T>&); \
template BlockedSpinningHamiltonian<T>::BlockedSpinningHamiltonian(const X&, double, double, const RotorInfo&, const IntervalSamplerBase&, const BlockedMatrix<T>&);

#define LCM_INSTANTIATE_DIAGONAL_HAMILTONIANS(X)\
  template BlockedDiagonalStaticHamiltonian::BlockedDiagonalStaticHamiltonian(const X&, const HamiltonianStore<double>&);\
  template BlockedDiagonalStaticHamiltonian::BlockedDiagonalStaticHamiltonian(const X&);\
  template BlockedDiagonalSpinningHamiltonian::BlockedDiagonalSpinningHamiltonian(const X&, const HamiltonianStore<space_T>&, double, double, const RotorInfo&);\
  template BlockedDiagonalSpinningHamiltonian::BlockedDiagonalSpinningHamiltonian(const X&, double, double, const RotorInfo&);

// template<class T1,class T2> void SpinOpGenerator::full(BlockedMatrix<T1>& d, const BaseList<T2>& a) const
// {
//   if (mzblocking_.empty()) {
//     create(d);
//     cmatrix d0(d.front());
//     ::libcmatrix::full(d0,a);
//   }
//   else
//     ::libcmatrix::full(d,a,mzblocking_);
// }

  //Implementation details

//   template<class Source, class PropUse, class TypeSD> bool
//   RawMetaSpectrum<Source,PropUse,TypeSD>::next() {
//     if (mzeig_==0) {
//       finished_=true;
//       return false;
//     }
//     if ((mzeig_==0) || (mzeigSD_==0))
//       throw Mismatch("MetaSpectrum: run out of blocks");
//     mzeig_--;
//     mzeigSD_--;
//     if (dirn_!='M')
//       shuffle();
//     if (!set_HUs(mzeig_,dirn_=='M'))
//       throw InternalError("MetaSpectrum::next: this shouldn't happen!");
//     set_sigma0dets();
//     return true;
//   }

//   template<class Source, class PropUse, class TypeSD>
//     template<class PropBase> BaseMetaObj<Source,PropUse,TypeSD>::BaseMetaObj(const PropBase& objv, Source* pgenv, mxflag::tempflag tflag,int nobsv, double startv, double endv, const TypeSD& sigma0v, const TypeSD& detectv, int verbosev)
//     : 	Hpgenp_(pgenv,tflag),
// 	sigma0_(sigma0v), detect_(detectv),
//     verbose_(verbosev),
// 	start_(startv), end_(endv),
// 	objlist_(Ham_functions::eigblocks(*pgenv),objv,mxflag::normal) {
      
//     common_init(nobsv,endv-startv,true);
//   }

//   template<class Source, class PropUse, class TypeSD>
//     template<class PropBase> BaseMetaObj<Source,PropUse,TypeSD>::BaseMetaObj(const PropBase& objv, Source* Hv, mxflag::tempflag tflag,double periodv, const TypeSD& sigma0v, const TypeSD& detectv, int verbosev)
//     : Hpgenp_(Hv,tflag),
//       sigma0_(sigma0v), detect_(detectv),
//     verbose_(verbosev),
//       start_(0.0), end_(periodv),
//       objlist_(Ham_functions::eigblocks(*Hv),objv,mxflag::normal)
//   {  
//     common_init(1,periodv,false);
//   }

//   template<class Source, class PropUse, class TypeSD>
//   void RawMetaSpectrum<Source,PropUse,TypeSD>::reset()
//     {
//       if (!blkspec_) { //No overlap between sigm0 and detect!
// 	finished_=true;
// 	return;
//       }

//     mzeig_=mzblocks_-1;
//     mzeigSD_=mzblocksSD_-1;

//     finished_=false;
//     havenext_=false;
//     const bool dodiag=(dirn_=='M') || ismiddleblock(mzeigSD_);
//     bool isset=set_HUs(mzeig_,dodiag);
//     if (!isset && !dodiag) {
//       shuffle();
//       if (mzeig_==0)
// 	throw Mismatch("MetaSpectrum: run out of sigma0/detect blocks");
//       mzeig_--;
//       isset=set_HUs(mzeig_,false);
//     }
//     if (!isset)
//       throw InternalError("MetaProp::reset: Failed to initialise objects!");
//     set_sigma0dets();
//     reset_blocks(blocks_-1);
//   }

//   namespace {
//     complex conj_(const complex& a) { return conj(a); }
//     double conj_(double a) { return a; }
//   };

//   template<class Source, class PropUse, class TypeSD>
//     bool RawMetaSpectrum<Source,PropUse,TypeSD>::operator() (typename PropUse::amplitude_type& amp, double& freq) {
//     if (finished_)
//       return false;

//     if (havenext_) {
//       amp=nextamp_;
//       freq=nextfreq_;
//       havenext_=false;
//       return true;
//     }

//     while (!(*curobjp_)(amp,freq)) {
//       if (blockpos_==0) {
// 	if (!next()) {
// 	  finished_=true;
// 	  return false;
// 	}
// 	reset_blocks(blocks_-1);
//       }
//       else
// 	reset_blocks(blockpos_-1);
//     }
//     amp*=totweight_; //any weighting from mz symmetry
//     if (blkspec_.isherm) { //if auto-hermitian, line up flipped amp,freq pair
//       havenext_=true;
//       nextamp_=conj_(amp);
//       nextfreq_=-freq;
//     }
//     return true;
//   }    

//   template<class Source, class PropUse, class TypeSD>
//   void RawMetaSpectrum<Source,PropUse,TypeSD>::set_sigma0det(size_t blk, Bool2Type<true>)
//   {
//     if (!!detect_)
//       std::cerr << "Warning: detection operator is being ignored\n";
//     if (!sigma0_) {
//       if (Ham_functions::size(*Hpgenp_,mzeigSD_,blk)) //this eigenvalue doesn't exist?
// 	objlist_(blk).observe();
//     }
//     else {
//       const cmatrix& subsigma0(sigma0_(mzeigSD_,blk));
//       if (subsigma0.size())
// 	objlist_(blk).observe(subsigma0);
//     }
//   }

//   template<class Source, class PropUse, class TypeSD>
//   void RawMetaSpectrum<Source,PropUse,TypeSD>::set_sigma0det(size_t blk, Bool2Type<false>)
//   {
//     if (!detect_)
//       throw Failed("MetaSpectrum: detect operator must be supplied");
//     const cmatrix& subsigma0(sigma0_(mzeigSD_,blk));
//     if (subsigma0.size())
//       objlist_(blk).observe(subsigma0,detect_(mzeigSD_,blk));
//   }

//   template<class Source, class PropUse, class TypeSD> void
//   BaseMetaObj<Source,PropUse,TypeSD>::set_U(PropUse& obj, Int2Type<2>, Bool2Type<false>, Bool2Type<true>, bool dodiag) {
//     if (dodiag)
//       obj.set_U(passlist_.front(),period_);
//     else
//       obj.set_U(dirn_,passlist_.front(),period_);
//   }

//   template<class Source, class PropUse, class TypeSD> void
//   BaseMetaObj<Source,PropUse,TypeSD>::set_U(PropUse& obj, Int2Type<3>, Bool2Type<false>, Bool2Type<true>, bool dodiag) {
//     if (dodiag)
//       obj.set_Us(passlist_,period_);
//     else
//       obj.set_Us(dirn_,passlist_,period_);
//   }      

//   template<class Source, class PropUse, class TypeSD> void
//   BaseMetaObj<Source,PropUse,TypeSD>::set_U(PropUse& obj, Int2Type<2>, Bool2Type<true>, Bool2Type<true>, bool dodiag) {
//     if (dodiag)
//       obj.set_U(passlist_.front());
//     else
//       obj.set_U(dirn_,passlist_.front());
//   }

//   template<class Source, class PropUse, class TypeSD> void
//   BaseMetaObj<Source,PropUse,TypeSD>::set_U(PropUse& obj, Int2Type<3>, Bool2Type<true>, Bool2Type<true>, bool dodiag) {
//     if (dodiag)
//       obj.set_Us(passlist_);
//     else
//       obj.set_Us(dirn_,passlist_);
//   }      

//   template<class Source,class PropUse,class TypeSD>
//   template<bool isTD,class T> void BaseMetaObj<Source,PropUse,TypeSD>::set_H_2nd_(PropUse& obj, Bool2Type<isTD> istd, size_t mzeig, size_t blk,bool dodiag,Type2Type<BlockedHamiltonian<T> >) const
//   {
//     const Matrix<T>& Huse((*Hpgenp_)(mzeig,blk));
//     const BaseList<double> Hzeeman(Ham_functions::Hzeeman(*Hgenp_,mzblk,eigblk));
//     if (Hgenp_->isdiagonal()) {
//       hermitian_eigenvalues_second(seigs_,Huse,Hzeeman);
//       set_H_(obj,istd,seigs_,dodiag);
//     }
//     else {
//       Matrix<T> V;
//       hermitian_eigensystem_second(V,seigs_,Huse,blkstr_,Hzeeman);
//       set_H_(obj,istd,V,seigs_,dodiag);
//     }
//   }      

//   template<class Source, class PropUse, class TypeSD>
//     template<bool isTD,bool allowU> bool BaseMetaObj<Source,PropUse,TypeSD>::set_HU(PropUse& obj, Bool2Type<true>, Bool2Type<allowU>, Bool2Type<isTD> istd, size_t mzeig, size_t blk, bool dodiag) {
//     switch (Ham_functions::quadrupole_order(*Hpgenp_)) {
//     case 1: {
//       const typename Source::suboutput_type& Huse((*Hpgenp_)(mzeig,blk));
//       set_H_(obj,istd,Huse,dodiag);
//     }
//       break;
//     case 0: {
//       typename Source::suboutput_type Huse((*Hpgenp_)(mzeig,blk));
//       Huse+=Ham_functions::Hzeeman((*Hpgenp_),mzeig,blk);
//       set_H_(obj,istd,Huse,dodiag);
//       if (obslarmor_)
// 	obj.larmor(obslarmor_);
//     }
//       break;
//     case 2:
//       set_H_2nd_(obj,istd,mzeig,blk,dodiag,Type2Type<Source>());
//       break;
//     default:
//       throw Failed("BaseMetaObj: unhandled quadrupole order");
//     }
//     return true;
//   }

// /*   template<class Source, class PropUse, class TypeSD> */
// /*   template<bool allowU> bool BaseMetaObj<Source,PropUse,TypeSD>::set_HU(PropUse& obj, Bool2Type<true>, Bool2Type<allowU>, Bool2Type<false>, size_t mzeig, size_t blk, bool dodiag) { */
// /*     const typename Source::suboutput_type& Huse=(*Hpgenp_)(mzeig,blk); */
// /*     if (dodiag) */
// /*       obj.set_H(Huse); */
// /*     else */
// /*       obj.set_H(dirn_,Huse); */
// /*     return true; */
// /*   } */

//   template<class Source, class PropUse, class TypeSD> bool
//   BaseMetaObj<Source,PropUse,TypeSD>::set_HUs(size_t mzeig, bool dodiag) 
//   {
//     bool isinit=false;
//     //if (mzeig==blocks_)
//     //  throw Failed("set_HUs: Hamiltonian iterator overrun");
//     for (size_t blk=blocks_;blk--;) {
//       PropUse& obj_=objlist_(blk);
//       if (Ham_functions::size(*Hpgenp_,mzeig,blk)) { //non-zero block?
// 	if (!set_HU(obj_ ,Bool2Type<SignalGenerator_traits<PropUse>::allowH && Propagator_traits<Source>::allowsetH >(), Bool2Type<SignalGenerator_traits<PropUse>::allowU>(), Bool2Type<SignalGenerator_traits<PropUse>::timedomain>(),mzeig,blk, dodiag)) {
// 	  set_U(obj_, Int2Type<SignalGenerator_traits<PropUse>::input_dimensionality>(), Bool2Type<SignalGenerator_traits<PropUse>::timedomain>(), Bool2Type<SignalGenerator_traits<PropUse>::allowU>(),dodiag);
// 	}
// 	if (!!obj_) 
// 	  isinit=true; 
//       }
//       else {
// 	if (dodiag)
// 	  obj_.clear(dirn_);
// 	else
// 	  obj_.clear();
//       }
//     }
//     return isinit; //only return false if no objects initialised (some may be unused)
//   }

//   template<class Source, class PropUse, class TypeSD> void
//     RawMetaSpectrum<Source,PropUse,TypeSD>::set_sigma0dets()
//   {
//     for (size_t blk=blocks_;blk--;) {
//       set_sigma0det(blk,Bool2Type<SignalGenerator_traits<PropUse>::allowED>());
//     }
//   }

//   template<class Source, class PropUse, class TypeSD>
//   void
//     BaseMetaObj<Source,PropUse,TypeSD>::initSD()
//   {
//     if (!sigma0_)
//       throw Failed("Case not handled yet");
//     if (blocks_!=sigma0_.eigblocks())
//       throw Mismatch("MetaObj: Hamiltonian and sigma0 have different eigenvalue structure");

//     mzblocksSD_=Ham_functions::actual_mzblocks(sigma0_);

//     if (!detect_)
//       blkspec_=sigma0_.blockstructure();
//     else {
//       if (!arematching(sigma0_,detect_))
// 	throw Mismatch("MetaObj: sigma0 / detect don't match");	  

//       const block_pattern& sigma0str(sigma0_.blockstructure());
//       const block_pattern& detectstr(detect_.blockstructure());

//       if (sigma0str.isherm == detectstr.isherm) {
// 	if (sigma0str.isherm) {
// 	  if (sigma0str.whichdiag == detectstr.whichdiag)
// 	    blkspec_=sigma0str;
// 	}
// 	else {
// 	  if (sigma0str.whichdiag == -detectstr.whichdiag)
// 	    blkspec_=sigma0str;
// 	}
//       }
//       else {//one is hermitian, one not
// 	if (::std::abs(sigma0str.whichdiag)==::std::abs(detectstr.whichdiag))
// 	  blkspec_=block_pattern( sigma0str.isherm ? -detectstr.whichdiag : sigma0str.whichdiag, sigma0str.skip, false);
//       }

//       if (!blkspec_)
// 	std::cerr << "Warning: sigma0 and detection operator do not overlap:\n" << sigma0str << detectstr;
//     }

//     if (verbose)
//       ::std::cout << "Block structure: " << blkspec_;

//     if (blkspec_.whichdiag && usemzsym_ && ((Ham_functions::mzblocks(*Hpgenp_) & 1)==0))
//       midblock_=mzblocks_/2;
//     else
//       midblock_=mzblocks_+10000; //this will never match mz eig

//     if (blkspec_.whichdiag==0)
//       dirn_='M';
//     else
//       dirn_=(blkspec_.whichdiag>0) ? 'R' : 'C';
//   }

//   template<class Source, class PropUse, class TypeSD> void
//   BaseMetaObj<Source,PropUse,TypeSD>::common_init(int nobsv, double periodv, bool havevalidv)
//   {
//     blocks_=Ham_functions::eigblocks(*Hpgenp_);
//     mzblocks_=Ham_functions::actual_mzblocks(*Hpgenp_);
//     usemzsym_=Ham_functions::usemzsym(*Hpgenp_);

//     obslarmor_=0.0;

//     if (nobsv<=0 || periodv<0.0)
//       throw InvalidParameter("MetaObj initialise");
//     nobs_=nobsv;
//     static const bool isinhomo= !SignalGenerator_traits<PropUse>::allowH && !SignalGenerator_traits<PropUse>::allowU;
//     static const bool need_period=!isinhomo && ((!SignalGenerator_traits<PropUse>::allowH) ^ (SignalGenerator_traits<PropUse>::timedomain));
//     if ((periodv==0.0) && need_period)
//       throw InvalidParameter("MetaObj: period must be specified (and >0)");
  
//     havevalidpgen_=havevalidv;
//     period_=periodv;

//     initSD();
//     //Special case for inhomogeneous objects (both allowH and allowU false)
//     initobjects( Bool2Type<isinhomo>() );
//   }

} //namespace libcmatrix

#endif


