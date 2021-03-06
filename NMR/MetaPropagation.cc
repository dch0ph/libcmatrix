#undef LCM_SUPPRESS_VIEWS
#include "lcm_MetaPropagation.h"
#include "lcm_basethreads.h"
#include "BaseMetaPropagation.h"
#include <string>

#ifdef NDEBUG
#define CHECKCOND (verbose_>1)
#else
#define CHECKCOND true
#endif

namespace libcmatrix {
  
  namespace {

    interaction_t interaction_count=I_AUTO;
    List<std::string> interaction_name_(I_AUTO+1);
    
    struct Proxy_ {
      Proxy_() {
	interaction_name_(I_CS)="chemical shift";
	interaction_name_(I_DIPOLE)="dipole coupling";
	interaction_name_(I_J)="J coupling";
	interaction_name_(I_QUAD)="quadrupole coupling";
      }
    };

    Proxy_ proxy_; //declare built-in interactions
  }

  template<size_t N> struct CrystalStructure_iterator_impl : public CrystalStructure_iterator {
    CrystalStructure_iterator_impl(const CrystalStructure& a, size_t M)
      : nspins_cell_(M), finished_(true) { 
      if ((N!=a.dimensions()) && !((N==1) && (a.dimensions()==0)))
	throw Mismatch("CrystalStructure_iterator_impl");
      init_(a.indexer());
    }

    CrystalStructure_iterator* clone() const { return new CrystalStructure_iterator_impl<N>(*this); }

    bool next(size_t& i, size_t& j) {
      if (finished_)
	return false;
      i=curi; j=curj;
      curi+=nspins_cell_;
      if (curi>=indexer_.size())
	finished_=true;
      else
	nextj_();
      return true;
    }
    void nextj_();
    void reset_();
    void init_(const Indexer<3>&);

    void reset(size_t i, size_t j) {
      if (i>=nspins_cell_)
	throw InvalidParameter("CrystalStructure_iterator: not in first cell");
      curi=i;
      curj=j;
      for (size_t k=N;k--;)
	celli[k]=0;
      reset_();
      finished_=false;
    }

    template<int M> void addwrap(Int2Type<M> which) {
      cellj[M]=celli[M]+relpos[M];
      if (cellj[M]>=indexer_.dimension(which))
	cellj[M]-=indexer_.dimension(which);
    }

    Indexer<N+1> indexer_;      
    size_t nspins_cell_;
    bool finished_;
    size_t curi,curj;
    size_t celli[N];
    size_t cellj[N];
    size_t basej;
    size_t relpos[N];
  };

  template<> void CrystalStructure_iterator_impl<1>::init_(const Indexer<3>& a) {
    indexer_.setdims(a.dimension(Int2Type<2>()),nspins_cell_);
  }
  template<> void CrystalStructure_iterator_impl<1>::reset_() {}
  template<> void CrystalStructure_iterator_impl<1>::nextj_() {
    static const size_t nspins_total(indexer_.size());
    curj+=nspins_cell_;
    if (curj>=nspins_total)
      curj-=nspins_total;
  }

  template<> void CrystalStructure_iterator_impl<2>::init_(const Indexer<3>& a) {
    indexer_.setdims(a.dimension(Int2Type<1>()),a.dimension(Int2Type<2>()),nspins_cell_);
  }
  template<> void CrystalStructure_iterator_impl<2>::reset_() {
    indexer_.reverse(relpos[0],relpos[1],basej,curj);
  }
  template<> void CrystalStructure_iterator_impl<2>::nextj_() {
    if (++celli[1]==indexer_.dimension(Int2Type<1>())) {
      celli[1]=0;
      celli[0]++;
    }
    addwrap(Int2Type<0>());
    addwrap(Int2Type<1>());
    curj=indexer_(cellj[0],cellj[1],basej);
  }

  template<> void CrystalStructure_iterator_impl<3>::init_(const Indexer<3>& a) {
    indexer_.setdims(a.dimension(Int2Type<0>()),a.dimension(Int2Type<1>()),a.dimension(Int2Type<2>()),nspins_cell_);
  }
  template<> void CrystalStructure_iterator_impl<3>::reset_() {
    indexer_.reverse(relpos[0],relpos[1],relpos[2],basej,curj);
  }
  template<> void CrystalStructure_iterator_impl<3>::nextj_() {
    if (++celli[2]==indexer_.dimension(Int2Type<2>())) {
      celli[2]=0;
      if (++celli[1]==indexer_.dimension(Int2Type<1>())) {
	celli[1]=0;
	celli[0]++;
      }
    }
    addwrap(Int2Type<0>());
    addwrap(Int2Type<1>());
    addwrap(Int2Type<2>());
    curj=indexer_(cellj[0],cellj[1],cellj[2],basej);
  }

  CrystalStructure::CrystalStructure(const BaseList<size_t>& dims, const Permutation& permv)
    : perm_(permv), ndims(dims.size()), higher(ndims>1)
  {
    const size_t cellsx= ndims ? dims.back() : 1;
    const size_t cellsy= higher ? dims(ndims-2) : 1;
    const size_t cellsz= (ndims>2) ? dims.front() : 1;
    indexer_.setdims(cellsz,cellsy,cellsx);
  }

  CrystalStructure::CrystalStructure(int N, const Permutation& permv)
    : indexer_(1,1,N), perm_(permv), ndims(1), higher(false) {}

  CrystalStructure_iterator* CrystalStructure::create_iterator(size_t M) const
  {
    switch (ndims) {
    case 1: case 0:
      return new CrystalStructure_iterator_impl<1>(*this,M);
    case 2:
      return new CrystalStructure_iterator_impl<2>(*this,M);
    case 3:
      return new CrystalStructure_iterator_impl<3>(*this,M);
    }
    throw InternalError("CrystalStructure::create_iterator");
  }

  size_t CrystalStructure::operator()(const BaseList<size_t>& ns) const
  {
    if (ns.size()!=ndims)
      throw Mismatch("CrystalStructure()");
    switch (ndims) {
    case 0:
      return 0;
    case 1:
      return ns.front();
    case 2:
      return indexer_(0,ns.front(),ns[1U]);
    default:
      return indexer_(ns.front(),ns[1U],ns[2U]);
    }
  }
    
  size_t CrystalStructure::operator()(size_t cz, size_t cy, size_t cx) const
  {
#ifndef NDEBUG
    if ((cz>=dimension(0)) || (cy>=dimension(1)) || (cx>=dimension(2)))
      throw BadIndex("CrystalSystem()");
#endif
    return higher ? indexer_(cz,cy,cx) : cx; //short-cut if only 1 dimension
  }
          
  std::ostream& operator << (std::ostream& ostr, const CrystalStructure& a)
  {
    if (a.ndims>2)
      ostr << a.dimension(Int2Type<0>()) << 'x';
    if (a.higher)
      ostr << a.dimension(Int2Type<1>()) << 'x';
    ostr << a.dimension(Int2Type<2>());
    if (a.haspermutation())
      ostr << "  permutation: " << a.permutation();
    return ostr;
  }
      
  const std::string& interaction_name(interaction_t n)
  {
    if (n>=interaction_name_.size())
      throw InvalidParameter("interaction_name: unknown interaction");
    return interaction_name_(n);
  }

  interaction_t interaction_register(const std::string& namestr)
  {
    if (find(interaction_name_.begin(),interaction_name_.end(),namestr)!=interaction_name_.end())
      throw Failed("interaction_register: interaction already exists");
    interaction_name_.push_back(namestr);
    return ++interaction_count;
  }

  interaction_t interaction_find(const std::string& name)
  {
    const std::string* idp=find(interaction_name_.begin(),interaction_name_.end(),name);
    if (idp==interaction_name_.end())
      throw Failed("interaction_find: interaction unknown");
    return idp-interaction_name_.begin()+1;  //BIT DODGY - assumes interactions are indexed from 1
  }
  
  productoperator_spec::productoperator_spec(const operator_spec& spec_)
    : scales(1,complex(1.0)),
      specs(ExplicitList<1,size_t>(1U),spec_),
      nuc(spec_.nuc)
  {
    //if (nuc==NULL_NUCLEUS)
    //  throw Failed("productoperator_spec: unable to identify nucleus type");
  }

  void productoperator_spec::push_back(complex scale,const operator_spec& spec_)
  {    
    const size_t lnuc=spec_.nuc;
    //    if (lnuc==NULL_NUCLEUS)
    //  throw Failed("productoperator_spec: unable to identify nucleus type");
    if (nuc!=lnuc)
      nuc=NULL_NUCLEUS;
    
    scales.push_back(scale);
    specs.push_back(BaseList<operator_spec>(1U,const_cast<operator_spec*>(&spec_)));
  }
  
  productoperator_spec::productoperator_spec(const operator_spec& spec_, const basespin_system& sys_) 
    : scales(1,complex(1.0)),
      specs(ExplicitList<1,size_t>(1U),spec_),
      nuc(spec_.nucleus(sys_)) {}
  
  bool operator== (const productoperator_spec& a, const productoperator_spec& b) {
    return (a.specs==b.specs) && (a.scales==b.scales);
  }
  
  bool operator!= (const productoperator_spec& a, const productoperator_spec& b) {
    return (a.specs!=b.specs) || (a.scales!=b.scales);
  }
  
  void productoperator_spec::push_back(complex scale,const BaseList<operator_spec>& spec, const basespin_system& sys)
  {
    if (spec.empty())
      throw InvalidParameter("productoperator_spec: empty operator definition");

    if (spec.size()>1) {
      for (size_t i=spec.size();i--;) {
	if (spec(i).nuc!=NULL_NUCLEUS)
	  throw InvalidParameter("productoperator_spec: can't use sum (F) operators in non-trivial product operators");
      }
    }
    
    size_t lnuc(spec.front().nucleus(sys));
    for (size_t i=1;i<spec.size();i++) {
      if (spec(i).nucleus(sys)!=lnuc) {
	lnuc=NULL_NUCLEUS;
	break;
      }
    }
    nuc = (specs.empty() || (nuc==lnuc)) ? lnuc : NULL_NUCLEUS;
    
    scales.push_back(scale);
    specs.push_back(spec);
  }

//   namespace {
//     std::ostream& pretty_print(std::ostream& ostr, double val, bool needintro)
//     {
//       if (needintro && (val>=0.0))
// 	ostr << '+';
//       if (val!=1.0)
// 	ostr << val;
//       return ostr;
//     }
//   }

  namespace {
    const char PUREERR[]="productoperator_spec: scale factor must be pure real or pure imaginary";
  }

  void productoperator_spec::print(std::ostream& ostr) const
  {
    if (size()==0) {
      ostr << "<empty>";
      return;
    }
    const LCM_IOS::fmtflags oldflags=ostr.flags();
    ostream_controller& ctrl(cmatrix_ostream_controller(ostr));
    const ostream_controller::complexview_t oldview(ctrl.complexview);
    ctrl.complexview=ostream_controller::compact; //!< switch to R+Ii format
    for (size_t i=0;i<size();i++) {
      const complex cscale(scales(i));
      if (cscale==complex(1.0)) {
 	if (i)
 	  ostr << " + ";
      }
      else {
	ostr.setf(LCM_IOS::showpos); //!< force sign print
	ostr << cscale << ' ';
      }
	// 	if (imag(cscale)) {
// 	  if (real(cscale))
// 	    throw Failed(PUREERR);
// 	  pretty_print(ostr,imag(cscale),i!=0) << " i ";
// 	}
// 	else
// 	  pretty_print(ostr,real(cscale),i!=0) << ' ';
//       }
      const BaseList<operator_spec> speci(specs(i));
      for (size_t j=0;j<speci.size();j++)
	ostr << speci(j) << ' ';
    }
    ostr.flags(oldflags);
    ctrl.complexview=oldview;
    if (nuc!=NULL_NUCLEUS)
      ostr << '(' << nuctolabel(nuc) << ')';
    else
      ostr << "(mixed)";
  }
  

bool productoperator_spec::ismatching(const productoperator_spec& b) const
{
  const BaseList<operator_spec>& br(b.specs.row());
  const BaseList<operator_spec>& r(specs.row());

  return ((r.size()==1) && (br.size()==1) && arematching(r.front(),br.front()));
}

productoperator_spec& productoperator_spec::operator&= (const productoperator_spec& a)
{
  BaseList<operator_spec> thisr(specs.row());
  const BaseList<operator_spec>& ar(a.specs.row());
  if ((thisr.size()!=1) || (ar.size()!=1))
    throw Failed("productoperator_spec&=: only valid for simple operators");
  thisr.front() = thisr.front() & ar.front();
  //NB scale factor unchanged
  return *this;
}

  template<class T> void hermitian_eigensystem(BlockedMatrix<T>& dest, ListList<double>& eigs, const BlockedMatrix<T>& a)
  {
    if (!a)
      throw Undefined("hermitian_eigensystem");
    const size_t n=a.size();
    ScratchList<size_t> sizes(n);
    size_t k;
    for (k=n;k--;)
      sizes(k)=a(k).rows();
    dest.create(sizes);
    eigs.create(sizes);
    for (k=n;k--;)
      hermitian_eigensystem(dest(k),eigs(k),a(k));
  }

  template void hermitian_eigensystem(BlockedMatrix<double>&, ListList<double>&, const BlockedMatrix<double>&);
  template void hermitian_eigensystem(BlockedMatrix<complex>&, ListList<double>&, const BlockedMatrix<complex>&);

  void BaseStructure::set(size_t mzblocksv, size_t eigblocksv, bool usemzsymv) {
    usemzsym_=usemzsymv;
    mzblocks_=mzblocksv;
    const size_t actualmz= usemzsymv ? (mzblocksv+1)/2 : mzblocksv;
    if (actualmz==0)
      throw InternalError("BaseStructure::set");
    indexer_=Indexer<2>(actualmz,eigblocksv);
  }
  
std::ostream& operator<< (std::ostream& ostr, const HamiltonianStructure& a)
{
  ostr << "Diagonal: " << (a.isdiagonal() ? "Yes\n" : "No\n");
  ostr << "Quadrupole order: ";
  switch (a.quadrupole_order()) {
  case 0: 
    ostr << "exact\n"; break;
  case 2:
    ostr << "2nd (" << ((a.flags() & MetaFlags::ClassicSecondOrder) ? "classic" : "generalised") << ")\n";
    break;
  default:
    ostr << "1st\n"; break;
  }

  if (a.nspins()>1) {
    ostr << "Coupling matrix\n";
    spy(ostr,a.iscoupled());
  }
  return ostr;
}

template<typename T> void HamiltonianStructure::process_binary(const couplingstore<T>& couplings)
{
  if (!couplings)
    return;
  const size_t M=couplings.rows();
  const size_t spins=couplings.cols();
  const basespin_system& sys(*sysp_);

  for (size_t i=M;i--;) {
    for (size_t j=i+1;j<spins;j++) {
      if (couplings.isnonnull(i,j)) {
      	iscoupled_(i,j)=true;
	const size_t spinj= j % M;
	if (sys(i)==sys(spinj))
	  isdiagonal_=false;
      }
    }
  }
}

  template<typename T> HamiltonianStructure::HamiltonianStructure(const basespin_system& sys, const HamiltonianStore<T>& Hstore, const BaseList<interaction_t>& weakintsv, int flagsv)
    : sysp_(&sys,mxflag::nondynamic), weakints_(weakintsv), flags_(flagsv)
{
  create(Hstore);
}

template<typename T> HamiltonianStructure::HamiltonianStructure(const basespin_system& sys, const HamiltonianStore<T>& Hstore, const BaseList<nuclei_spec>& nuclistv, const BaseList<interaction_t>& weakintsv,int flagsv)
  : sysp_(&sys,mxflag::nondynamic), blockingnuclei_(nuclistv), weakints_(weakintsv), flags_(flagsv)
{
  create(Hstore);
}

template<typename T> HamiltonianStructure::HamiltonianStructure(const basespin_system& sys, const HamiltonianStore<T>& Hstore, const char *label, const BaseList<interaction_t>& weakintsv, int flagsv)
  : sysp_(&sys,mxflag::nondynamic), blockingnuclei_(1,nuclei_spec(label)), weakints_(weakintsv), flags_(flagsv)
{
  create(Hstore);
}

void HamiltonianStructure::blockingnuclei(const BaseList<nuclei_spec>& newblocking)
{
  if (!isdiagonal_)
    blockingnuclei_=newblocking;
}

Warning<> MetaFlags::UseMzSymmetry_shift_warning("UseMzSymmetry incompatible with prescence of shift interactions",&lcm_base_warning);
//Warning<> blocking_nonsecular_warning("mz blocking in combination with non-secular Hamiltonians is only partially supported and not evaluated",&lcm_base_warning);
Warning<> diagonalsigma0detectmissed_warning("special case of identity sigma0/detect detected but not exploited for frequency domain calculations",&lcm_base_warning);

template<typename T> void HamiltonianStructure::create(const HamiltonianStore<T>& Hstore)
{
  M_=Hstore.nspins_cell();
  if (sysp_->nspins()!=M_)
    throw Mismatch("HamiltonianStructure");
  
  iscoupled_.create(M_,Hstore.nspins(),false);

  isdiagonal_=true;
  quadrupole_order_=1;
  bool haveqpole=false;

  if (Hstore.haslinear() && (flags_ & MetaFlags::UseMzSymmetry)) {
    MetaFlags::UseMzSymmetry_shift_warning.raise();
    flags_&=~MetaFlags::UseMzSymmetry;
  }
  
  //process couplings to determine which nuclei are coupled
  typedef typename HamiltonianStore<T>::couplingmap_type::const_iterator iter_type;
  const iter_type iterend(Hstore.couplings_end());
  iter_type iter(Hstore.couplings_begin());
  while (iter!=iterend) {
    if (!isweak(iter->first))
      process_binary(iter->second);
    ++iter;
  }
 
  const BaseList<T> qcouplings(Hstore.get_quadrupole());
  nonsecular_.create(M_,false);

  const bool classic=(flags_ & MetaFlags::ClassicSecondOrder);

  if (!qcouplings.empty()) {
    std::set<size_t> nonseclist;

    for (size_t j=M_;j--;) {
      if (isnonnull(qcouplings(j))) {
	if ((*sysp_)(j).deg()<3)
	  throw Failed("HamiltonianStructure: quadrupole coupling specified for non-quadrupolar nucleus");
	const size_t order=Hstore.get_quadrupole_order(j);
	if ((order!=1) && !classic) {
	  nonsecular_(j)=true;
	  nonseclist.insert((*sysp_)(j).nucleus());
	}
	if (haveqpole) {
	  if (order!=1) {
	    if ((quadrupole_order_!=1) && (order!=quadrupole_order_))
	      throw Failed("HamiltonianStructure: mixture of non-secular quadrupole types (0 & 2)");
	    quadrupole_order_=order;
	  }
	}
	else {
	  quadrupole_order_=order;
	  haveqpole=true;
	}
      }
    }

    if (!(nonseclist.empty()) && !(blockingnuclei_.empty())) {
      for (size_t j=blockingnuclei_.size();j--;) {
	if (nonseclist.count(blockingnuclei_(j)()))
	  throw Failed("HamiltonianStructure: cannot block nucleus subject to non-secular interactions");
      }
      //      blocking_nonsecular_warning.raise();
    }
  }

  if (classic) {
    if ((quadrupole_order_<1) || (quadrupole_order_>2))
      throw InvalidParameter("HamiltonianStructure: classic second order method is only compatible with quadrupole orders of 1 or 2");
  }
  else {
    if (quadrupole_order_!=1)
      isdiagonal_=false; //!< hamiltonian not diagonal if non-secular
  }

  if (isdiagonal_)
    blockingnuclei_.clear(); //suppress blocking if diagonal

//   if ((quadrupole_order_!=1) && (pfreq_==0))
//     throw InvalidParameter("HamiltonianStore: proton frequency must be set");
}

#define LCM_INSTANTIATE_HAMS(X)\
  template HamiltonianStructure::HamiltonianStructure(const basespin_system&, const HamiltonianStore<X>&, const BaseList<interaction_t>&, int);\
  template HamiltonianStructure::HamiltonianStructure(const basespin_system&, const HamiltonianStore<X>&,const BaseList<nuclei_spec>&, const BaseList<interaction_t>&, int);\
  template HamiltonianStructure::HamiltonianStructure(const basespin_system&, const HamiltonianStore<X>&,const char* label, const BaseList<interaction_t>&, int);

  LCM_INSTANTIATE_HAMS(double)
  LCM_INSTANTIATE_HAMS(space_T)

    const complex& BlockedOperator::identityscale(size_t mzeig, size_t blk) const 
{
  static const complex zero(0.0);
  return identityscale_.empty() ? zero : identityscale_(index(mzeig,blk)); 
} 

double norm(const BlockedOperator& a)
{
  const block_pattern& blkstr(a.blockstructure());
  block_pattern::iterator mziter(blkstr);
  const size_t eigblks(a.eigblocks());

  double sum=0.0;
  size_t r,c,mzeigSD;
  double middlescale;
  const double basescale=(blkstr.isherm ? 2.0 : 1.0);
  while (mziter.next(r,c,mzeigSD,middlescale)) {
    const double scale= basescale*middlescale;
    for (size_t blk=eigblks;blk--;) {
      const cmatrix& ablk(a(mzeigSD,blk));
      sum+=scale*norm(ablk);
    }
  }
  return sum;
}
  
void BlockedOperator::swap(BlockedOperator& a) 
{
  if (!arematching(*this,a))
    throw Mismatch("BlockedOperator::swap: operator structures don't match");
  store_.swap(a.store_);
}

template<class T> void BlockedOperator::scale_(T scale, const BlockedFilter& filt)
{
  if (!filt.ismatching(*this))
    throw Mismatch("BlockedOperator::scale");
  invalidatecache();
  const BlockedMatrix<bool> asblocked(filt.row());
  for (size_t i=asblocked.size();i--;) {
    cmatrix curm(store_(i),mxflag::nondynamic);
    const Matrix<bool> curb(asblocked(i),mxflag::nondynamic);
    if (!arematching(curm,curb))
      throw Mismatch("BlockedOperator::scale",curm.rows(),curm.cols(),curb.rows(),curb.cols());    
    BaseList<complex> rowm(curm.row());
    const BaseList<bool> rowb(curb.row());
    for (size_t i=rowm.size();i--;) {
      if (rowb(i))
	rowm(i)*=scale;
    }
  }
}

void BlockedOperator::scale(double scale, const BlockedFilter& filt)
{
  scale_(scale,filt);
}

void BlockedOperator::scale(const complex& scale, const BlockedFilter& filt)
{
  scale_(scale,filt);
}

void BlockedOperator::apply(const BlockedFilter& filt)
{
  if (!filt.ismatching(*this))
    throw Mismatch("BlockedOperator::apply");
  invalidatecache();
  const BlockedMatrix<bool> asblocked(filt.row());
  for (size_t i=asblocked.size();i--;) {
    cmatrix curm(store_(i),mxflag::nondynamic);
    curm.emultiply(asblocked(i));
  }
}

BlockedOperator& BlockedOperator::operator+=(const BlockedOperator& v)
{ 
  if (!arematching(*this,v))
    throw Mismatch("BlockedOperator+=");
  store_+=v.row();
  invalidatecache();
  return *this;
}

BlockedOperator& BlockedOperator::operator-=(const BlockedOperator& v)
{ 
  if (!arematching(*this,v))
    throw Mismatch("BlockedOperator-=");
  store_-=v.row();
  invalidatecache();
  return *this;
}

template<typename T> void mla_(BlockedOperator& d, T a, const BlockedOperator& b)
{
  if (!d) {
    d=b;
    d*=a;
  }
  else {
    if (!arematching(d,b))
      throw Mismatch("BlockedOperator:mla");
    mla(d.row(),a,b.row());
    d.invalidatecache();
  }
}

void mla(BlockedOperator& d, double a, const BlockedOperator& b) { mla_(d,a,b); }
void mla(BlockedOperator& d, complex a, const BlockedOperator& b) { mla_(d,a,b); }

void BlockedOperator::makeidentitycache() const
{
  if (!(*this))
    return;
  const size_t n=store_.size();
  identityscale_.create(n);
  for (size_t i=n;i--;)
    identityscale_(i)=findscale_(store_(i));  
}

complex BlockedOperator::findscale_(const Matrix<complex>& a)
{
  static const complex zero(0.0);
  if (!a || !issquare(a))
    return zero;
  complex scalefac(zero);//!< intialisation should not be necessary but avoids warning
  const size_t n=a.rows();
  for (size_t r=0;r<n;r++) {
    if (r==0)
      scalefac=a(0U,0U);
    else {
      if (a(r,r)!=scalefac)
	return zero; //!< diagonal elements not identical
    }
    for (size_t c=r;c--;) {
      if ((a(r,c)!=zero) || (a(c,r)!=zero))
	return zero;
    }
  }
  return scalefac;      
}

namespace {
  template<class T> struct locconj {
    T operator()(const T& v) { return v; }
  };
  template<> struct locconj<complex> {
    complex operator()(const complex& v) { return conj(v); }
  };

  template<class T> void full_( Matrix<T>& d, const BlockedMatrix<T>& source, const BaseStructure& opstructure, const SpinOpGeneratorBase& opgen, const block_pattern& blkstr, const T& empty =T(0))
  {
    const size_t nstates=opgen.size();
    const BaseList< ListList<size_t> > blocking(opgen.mzblocking());
    const Indexer<2>& Hindexer(opgen.structure().indexer());
    locconj<T> conjobj;
    d.create(nstates,nstates,empty);
    size_t r,c,mzeigSD;
    block_pattern::iterator mziter(blkstr);
    
    const bool isherm=blkstr.isherm;
    while (mziter.next(r,c,mzeigSD)) {
      for (size_t eblk=opgen.eigblocks();eblk--;) {
	const size_t Hblkr=Hindexer(r,eblk);
	const size_t Hblkc=Hindexer(c,eblk);
	const BaseList<size_t> rindices(blocking(Hblkr).row());
	const BaseList<size_t> cindices(blocking(Hblkc).row());
	const size_t ind(opstructure.index(mzeigSD,eblk));
	const Matrix<T> ablk(source(ind)); //(*this)(mzeigSD,eblk));
	d(rindices,cindices)=ablk;
	if (isherm) {
	  for (size_t r=ablk.rows();r--;)
	    for (size_t c=ablk.cols();c--;)
	      d(cindices(c),rindices(r))=conjobj(ablk(r,c));
	}
      }
    }
  }
} // anonymous namespace  

void BlockedOperator::full(cmatrix& d) const
{
  if (!opgenp_)
    throw Undefined("BlockedOperator: full");
  full_(d, store_, *this, *opgenp_, blockstructure() );
}

void BlockedFilter::full( Matrix<bool>& d) const
{
  full_(d, store_, *this, opgen_, blockstructure() );
}

      //NB full compatibility not verified
  void BlockedOperator::unitary_simtrans(const BlockedMatrix<complex>& U)
  {
    if (!opgenp_)
      throw Undefined("BlockedOperator: unitary_simtrans");

    if (blkspec_.isdiagonal()) {
      store_.unitary_simtrans(U);
      return;
    }
    block_pattern::iterator mziter(blkspec_);
    const size_t eigblks(eigblocks());
    if (U.size() % eigblks)
      throw Failed("propagator blocks incompatible with eigenvalue structure");

    size_t r,c,mzeigSD;
    //    bool ismiddle;
    const Indexer<2>& Hindexer(opgenp_->structure().indexer());
    
    while (mziter.next(r,c,mzeigSD)) {
      for (size_t blk=eigblks;blk--;) {
	cmatrix& sigmablk((*this)(mzeigSD,blk));
	if (sigmablk.empty())
	  continue;
	const size_t rindex(Hindexer(r,blk));
	const size_t cindex(Hindexer(c,blk));
	::libcmatrix::unitary_simtransLR(sigmablk,U(rindex),sigmablk,U(cindex));
      }
    }
    if ((mzeigSD+1)!=actual_mzblocks()) //!< assumes incrementing
      throw Mismatch("unitary_simtrans: sigma0/detect blocks",mzeigSD+1,actual_mzblocks());
  }

 complex trace_multiply_conj_transpose(const cmatrix& a, const cmatrix& b)
 {
   if (a.cols()!=b.cols())
     throw Mismatch("trace_multiply_conj_transpose");
   return conj(dot(a.row(),b.row()));
 }
  
 complex trace_multiply_conj_transpose(const BlockedMatrix<complex>& a, const BlockedMatrix<complex>& b)
 {
   complex sum(0.0);
   if (!arematching(a,b))
     throw Mismatch("trace_multiply_conj_transpose");
   for (size_t i=a.size();i--;) {
     if (!!b(i))
       sum+=trace_multiply_conj_transpose(a(i),b(i));
   }
   return sum;
 }

complex trace_multiply(const BlockedOperator& a, const BlockedOperator& b)
{
  const block_pattern& blkstra(a.blockstructure());
  const block_pattern& blkstrb(b.blockstructure());

  block_pattern tblkstrb(blkstrb);
  tblkstrb.transpose();
  if (!blkstra.ismatching(tblkstrb))
    return complex(0.0); //no overlap between a and b

  block_pattern::iterator mziter(blkstra);
  const size_t eigblks(a.eigblocks());

  complex sum(0.0);
  size_t r,c,mzeigSD;
  double scale;
    
  while (mziter.next(r,c,mzeigSD,scale)) {
    for (size_t blk=eigblks;blk--;) {
      const cmatrix& ablk(a(mzeigSD,blk));
      if (!a)
	continue;
      const cmatrix& bblk(b(mzeigSD,blk));
      
      if (blkstrb.isherm) {
	if (blkstra.isherm)
	  sum+=scale*2.0*real(trace_multiply_conj_transpose(ablk,bblk));
	else {
	  if (blkstra.coher>=0)
	    sum+=scale*trace_multiply_conj_transpose(ablk,bblk);
	  else
	    sum+=scale*trace_multiply(ablk,bblk);
	}
      }
      else {
	if (blkstra.isherm && (blkstrb.coher>=0))
	  sum+=scale*conj(trace_multiply_conj_transpose(bblk,ablk));
	else
	  sum+=scale*trace_multiply(ablk,bblk);
      }
    }
  }
  
  if ((mzeigSD+1)!=a.actual_mzblocks())
    throw Mismatch("unitary_simtrans: sigma0/detect blocks",mzeigSD+1,a.actual_mzblocks());

  return sum;
}

 template<> void BlockedMatrix<complex>::unitary_simtrans(const ListList<complex>& U)
   {
     if (!arematching(*this,U))
       throw Mismatch("unitary_simtrans");
     for (size_t i=size();i--;) {
       cmatrix& curblk((*this)(i));
       if (!!curblk)
	 curblk.unitary_simtrans(U(i));
     }
   }

 template<> void BlockedMatrix<complex>::unitary_isimtrans(const ListList<complex>& U)
   {
     if (!arematching(*this,U))
       throw Mismatch("unitary_isimtrans");
     for (size_t i=size();i--;) {
       cmatrix& curblk((*this)(i));
       if (!!curblk)
	 curblk.unitary_isimtrans(U(i));
     }
   }

 template<> void BlockedMatrix<complex>::unitary_simtrans(const BlockedMatrix<complex>& U)
   {
     if (!arematching(*this,U))
       throw Mismatch("unitary_simtrans");
     for (size_t i=size();i--;) {
       cmatrix& curblk((*this)(i));
       if (!!curblk)
	 curblk.unitary_simtrans(U(i));
     }
   }

 template<> void BlockedMatrix<complex>::unitary_isimtrans(const BlockedMatrix<complex>& U)
   {
     if (!arematching(*this,U))
       throw Mismatch("unitary_isimtrans");
     for (size_t i=size();i--;) {
       cmatrix& curblk((*this)(i));
       if (!!curblk)
	 curblk.unitary_isimtrans(U(i));
     }
   }

  template<class T> void propagator(BlockedMatrix<complex>& U, const BlockedMatrix<T>& H, double dt)
  {
    U.duplicate_structure(H);
    for (size_t k=H.size();k--;) {
      if (!!H(k))
	propagator(U(k),H(k),dt);
    }
  }

template<typename T> void propagator_ns(BlockedMatrix<complex>& U, const BlockedMatrix<T>& H, double dt, const ListList<double>& Hzeeman)
  {
    U.duplicate_structure(H);
    for (size_t k=H.size();k--;) {
      if (!!H(k)) {
	cmatrix Uk(U(k),mxflag::nondynamic);
	propagator_ns(Uk,H(k),dt,Hzeeman(k));
      }
    }
  }

template<typename T1,typename T2> void propagator_ns_Hrf(BlockedMatrix<complex>& U, const BlockedMatrix<T1>& Hsys,  const BlockedMatrix<T2>& Hrf, double dt, const ListList<double>& Hzeeman)
{
  U.duplicate_structure(Hsys);
  for (size_t k=Hsys.size();k--;) {
    if (!!Hsys(k)) {
      cmatrix Uk(U(k),mxflag::nondynamic);
      propagator_ns_Hrf(Uk,Hsys(k),Hrf(k),dt,Hzeeman(k));
    }
  }
}

template<typename T> void propagator_second(BlockedMatrix<complex>& U, const BlockedMatrix<T>& H, double dt, const BaseList< ListList<size_t> >& Zstruct, const ListList<double>& Hzeeman)
  {
    U.duplicate_structure(H);
    for (size_t k=H.size();k--;) {
      if (!!H(k)) {
	cmatrix Uk(U(k),mxflag::nondynamic);
	propagator_second(Uk,H(k),dt,Zstruct(k),Hzeeman(k));
      }
    }
  }

template<typename T1, typename T2> void propagator_second_Hrf(BlockedMatrix<complex>& U, const BlockedMatrix<T1>& Hsys,  const BlockedMatrix<T2>& Hrf, double dt, const BaseList< ListList<size_t> >& Zstruct, const ListList<double>& Hzeeman)
{
  U.duplicate_structure(Hsys);
  for (size_t k=Hsys.size();k--;) {
    if (!!Hsys(k)) {
      cmatrix Uk(U(k),mxflag::nondynamic);
      propagator_second_Hrf(Uk,Hsys(k),Hrf(k),dt,Zstruct(k),Hzeeman(k));
    }
  }
}

void propagator(ListList<complex>& U, const ListList<double>& H, double dt)
{
  U.duplicate_structure(H);
  BaseList<complex> Urow(U.row());
  propagator(Urow,H.row(),dt);
}

  template void propagator(BlockedMatrix<complex>&, const BlockedMatrix<double>&, double);
  template void propagator(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, double);

template void propagator_ns(BlockedMatrix<complex>&, const BlockedMatrix<double>&, double, const ListList<double>&);
template void propagator_ns(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, double, const ListList<double>&);
template void propagator_ns_Hrf(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, const BlockedMatrix<complex>&, double, const ListList<double>&);
template void propagator_ns_Hrf(BlockedMatrix<complex>&, const BlockedMatrix<double>&, const BlockedMatrix<complex>&, double, const ListList<double>&);
template void propagator_ns_Hrf(BlockedMatrix<complex>&, const BlockedMatrix<double>&, const BlockedMatrix<double>&, double, const ListList<double>&);

template void propagator_second(BlockedMatrix<complex>&, const BlockedMatrix<double>&, double, const BaseList< ListList<size_t> >&, const ListList<double>&);
template void propagator_second(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, double, const BaseList< ListList<size_t> >&, const ListList<double>&);
template void propagator_second_Hrf(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, const BlockedMatrix<complex>&, double, const BaseList< ListList<size_t> >&, const ListList<double>&);
template void propagator_second_Hrf(BlockedMatrix<complex>&, const BlockedMatrix<double>&, const BlockedMatrix<complex>&, double, const BaseList< ListList<size_t> >&, const ListList<double>&);
template void propagator_second_Hrf(BlockedMatrix<complex>&, const BlockedMatrix<double>&, const BlockedMatrix<double>&, double, const BaseList< ListList<size_t> >&, const ListList<double>&);

  template<class T> void unitary_simtrans(BlockedMatrix<complex>& d, const BlockedMatrix<T>& a, const BlockedMatrix<complex>& U)
  {
    if (!arematching(a,U))
      throw Mismatch("unitary_simtrans");    
    duplicate_structure(d,a);
    for (size_t k=a.size();k--;) {
      if (!!a(k))
	unitary_simtrans(d(k),a(k),U(k));
    }
  }

  template<class T> void unitary_isimtrans(BlockedMatrix<complex>& d, const BlockedMatrix<T>& a, const BlockedMatrix<complex>& U)
  {
    if (!arematching(a,U))
      throw Mismatch("unitary_isimtrans");    
    duplicate_structure(d,a);
    for (size_t k=a.size();k--;) {
      if (!!a(k))
	unitary_isimtrans(d(k),a(k),U(k));
    }
  }

  template void unitary_simtrans(BlockedMatrix<complex>&, const BlockedMatrix<double>&, const BlockedMatrix<complex>&);
  template void unitary_simtrans(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, const BlockedMatrix<complex>&);

  template void unitary_isimtrans(BlockedMatrix<complex>&, const BlockedMatrix<double>&, const BlockedMatrix<complex>&);
  template void unitary_isimtrans(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, const BlockedMatrix<complex>&);

   std::ostream& operator<< (std::ostream& ostr, const block_pattern& a)
   {
     ostr << "Coherence: " << a.coher << "  Blocks: " << a.blocks;
     //<< "  Mz step: " << a.mzstep << '\n';
     ostr << "  Hermitian: " << (a.isherm ? "Yes" : "No");
     ostr << "  Using mz symmetry: " << (a.usemzsym() ? "Yes" : "No");
#ifndef NDEBUG
     ostr << "  Maximum mzeig: " << a.maxind;
     //<< "  Has middle block: " << (a.hasmiddle ? "Yes" : "No");
     ostr << "  Indexer: " << a.indexer;
#endif
     return ostr << '\n';
   }

bool block_pattern::iterator::next(size_t& rind, size_t& cind, size_t& mzeigSDv, double& scale)
{
  bool ismiddle;
  const bool nextres=next(rind,cind,mzeigSDv,ismiddle);
  scale=(source.usemzsym() && !ismiddle) ? 2.0 : 1.0;
  return nextres;
}

bool block_pattern::iterator::next(size_t& rind, size_t& cind, size_t& mzeigSDv, bool& ismiddle)
{
  if (!next(rind,cind,mzeigSDv)) 
    return false;
  
  ismiddle=false;

  if (rind>=source.maxind) {
    assert(cind<source.maxind);
    rind=cind;
    ismiddle=true;
  }
  else {
    if (cind>=source.maxind) {
      assert(rind<source.maxind);
      cind=rind;
      ismiddle=true;
    }
  }
  return true;
}

bool productoperator_spec::isdiagonal() const
{
  const BaseList<operator_spec> asrow(specs.row());
  for (size_t i=asrow.size();i--;) {
    if (!asrow(i).isdiagonal())
      return false;
  }
  return true;
}
    
// bool productoperator_spec::issimple() const
// {
//   for (size_t i=specs.size();i--;) {
//     if (specs(i).size()>1)
//       return 
//   return (specs.items()==1); }
// }

block_pattern::block_pattern(const SpinOpGeneratorBase& opgen, const productoperator_spec& prodopspec)
 {
   const size_t nuc(prodopspec.nucleus());
   if (nuc!=NULL_NUCLEUS) {
     // if (prodopspec.issimple()) {
     bool failed=false;
     for (size_t i=0;i<prodopspec.specs.size();i++) {
       const BaseList<operator_spec> curop(prodopspec.specs(i));
       if (curop.size()>1) {
	 failed=true;
	 break;
       }
       if (i) {
	 if (*this!=block_pattern(opgen,curop.front())) {
	   failed=true;
	   break;
	 }
       }
       else
	 create(opgen,curop.front());
     }
     if (!failed)       
       return;

     if (!opgen.isblocked(nuc)) { //not blocked?
       create(opgen.active_blocking_info());
       return;
     }
   }
   if (opgen.isblocked() && !prodopspec.isdiagonal())
     throw Failed("Can't determine block_pattern for general off-diagonal productoperator: disable blocking");
   create(opgen.active_blocking_info());
 }
 
void blocking_information::create(const BaseList<size_t>& levelsv, const BaseList<size_t>& nucinds, bool usemzv)
{
  levels=levelsv;
  usemz=usemzv;
  steps.create(levels.size());
  if (usemz && (levels.size()>1))
    throw Failed("blocking_information: mz symmetry incompatible with heteronuclear systems");
  totallevels=1;
  for (size_t i=0;i<levels.size();i++) {
    const size_t ind(nucinds(i));
    steps(ind)=totallevels;
    totallevels*=levels(ind);
  }
  if (totallevels==0)
    throw InternalError("blocking_information");
}

std::ostream& operator<< (std::ostream& ostr, const blocking_information& a)
{
  if (a.totallevels==0)
    return ostr << "<undefined>\n";
  ostr << "Levels: " << a.levels << "  Multipliers: " << a.steps;
  return ostr << "  Use mz symmetry: " << (a.usemz ? "Yes\n" : "No\n");
}

size_t blocking_information::majorblocks(size_t which) const { 
  return totallevels / (levels(which)*steps(which)); 
}

size_t blocking_information::mzblocks(size_t which, size_t step) const 
{
  if (levels(which)<=step)
    throw Failed("blocking_information: coherence can't be supported");      
  return levels(which)-int(step);
}

void block_pattern::create(const SpinOpGeneratorBase& opgen, const operator_spec& opspec, bool active)
{
  const bool ishermv=(opspec.op=='x') || (opspec.op=='y');
  const size_t nuc=opspec.nucleus(opgen.spinsystem());
  if (opgen.isblocked(nuc) || !active) {
    const blocking_information& blkinfo(active ? opgen.active_blocking_info() : opgen.inactive_blocking_info());  
    create(blkinfo,
	   opgen.nuctoindex(nuc),
	   opspec.coherence(),
	   ishermv);
  }
  else {
    coher=0;
    isherm=false;
    maxind=blocks=opgen.actual_mzblocks();
    indexer.setdims(1,1,blocks);
    step=whichstep=majorstep=0;
    usemz=false;
    //  hasmiddle=false;
  }
}

block_pattern::iterator::iterator(const block_pattern& a)
  : source(a),
    isbelow(a.coher<0),
    iter(a.indexer,ScratchList<size_t>(0,1,2)),
    end(a.indexer),
    last(-1),
    mzeigSD(0)
{}

void block_pattern::iterator::reset()
{
  iter.reset();
  last=-1;
  mzeigSD=0;
}

//diagonal
void block_pattern::create(const blocking_information& infov)
{
  create(infov,0,0,false);
}

void block_pattern::create(const blocking_information& infov, size_t which, int coherv, bool lisherm)
{
  coher=coherv;
  isherm=lisherm;
  step=std::abs(coher);
  const size_t majorblocks=infov.majorblocks(which);
  const size_t repeatblocks=infov.steps(which);
  size_t mzblocks=infov.mzblocks(which,step);
  whichstep=step*repeatblocks;
  const size_t clevels=infov.levels(which);
  majorstep=repeatblocks*clevels;
  usemz=infov.usemz;
  if (usemz) {
    maxind=(clevels+1)/2;
    //    hasmiddle= (mzblocks & 1);    
    mzblocks= (mzblocks+1)/2;
    blocks=mzblocks;
  }
  else {
    blocks=mzblocks*repeatblocks*majorblocks;
    maxind=infov.totallevels;
    //    hasmiddle=false;
  }
  indexer.setdims(majorblocks,repeatblocks,(step ? 1 : 0)+mzblocks);
//   std::cout << "Actual mz blocks: " << infov.minorblocks(which,step) << '\n';
//   std::cout << "Repeats: " << iter.dimension(REPEAT) << '\n';
//   std::cout << "Major blocks: " << iter.dimension(MAJOR) << '\n';
}

bool block_pattern::iterator::next(size_t& r, size_t& c, size_t& mzeigSDv) 
{
  if (iter==end)
    return false;    
  mzeigSDv=mzeigSD++;

  if (source.step==0) {
    c=r=*iter;
    ++iter;
    return true;
  }
  if (last<0) {
    last=getindex();
    lasttick=iter.index(MZBLOCK);
    ++iter;
  }
  for (;;) {
    if (iter==end)
      return false;
    const size_t cur=getindex();
    const size_t curtick=iter.index(MZBLOCK);
    ++iter;
    if (curtick>lasttick) { //!< can't have reset
      if (isbelow) {
	r=cur;
	c=last;
      }
      else {
	c=cur;
	r=last;
      }
      last=cur;
      return true;
    }
    last=cur;
    lasttick=curtick;
  }
}

void block_pattern::clear()
{
  coher=whichstep=majorstep=step=blocks=maxind=0;
}

size_t block_pattern::iterator::getindex() const 
{
  return iter.index(MAJOR)*source.majorstep+iter.index(REPEAT)+iter.index(MZBLOCK)*source.whichstep;
}
  
int operator_spec::coherence() const
{
  if (issingletransition())
    return int(col)-int(number);
    
  switch (op) {
  case '+': case 'x': case 'y': case 'c': //!< corrected 18/11/09 to make 'c' a +1 coherence (rather than 0)
    return 1;
  case '-':
    return -1;
  case 'z': case 'a': case 'b':
    return 0;
  }
  throw InternalError("operator_spec::coherence");
}

// void block_pattern::create(const SpinOpGeneratorBase& opgen,size_t nuc, int lcoher, bool ishermv, bool isactive)
// {
//   //  const bool isactive=(which<0);
//   if ((lcoher==0) || (isactive && !opgen.isblocked(nuc))) {
//     create(opgen);
//     return;
//   }
//   const size_t index=opgen.nuctoindex(nuc);

//   mzstep=isactive ? opgen.active_mzstep_(index) : opgen.total_mzstep_(index);
//   isherm=ishermv;
//   coher=lcoher;
//   mzlevels=opgen.maxmzinds_(index);
//   totlevels= isactive ? opgen.mzblocks() : opgen.submzblocks_; //opgen.mzblocking_(which).size();
//   const size_t multiplier=totlevels / mzlevels;
//   blocks= multiplier*(mzlevels-::std::abs(coher));
//   if (blocks<=0)
//     throw Failed("block_pattern: invalid coherence");
//   if (blocks % mzstep)
//     throw InternalError("block_pattern");
//   if (isactive && opgen.usemzsym()) { //won't be set for heteronuclear system
//     hasmiddle=(blocks & 1);
//     blocks=(blocks+1)/2;
//   }
//   else
//     hasmiddle=false;
//   maxind= isactive ? opgen.actual_mzblocks() : totlevels;
// }

 void block_pattern::makepresent(Matrix<bool>& present) const
 {
   present.create(maxind,maxind,false);
   iterator mziter(*this);
   size_t r,c,mzeig;
   while (mziter.next(r,c,mzeig))
     present(r,c)=present(c,r)=true;
 }

block_pattern& block_pattern::operator&= (const block_pattern& a)
{
  if (ismatching(a)) {
    if (isherm) {
      if (!a.isherm) {
	isherm=false;
	coher=-a.coher;
      }
      return *this;
    } 
    if (a.isherm || (coher==-a.coher))
      return *this;
  }
  clear();
  return *this;
}
  
size_t SpinOpGeneratorBase::nuctoindex(const operator_spec& opspec) const
{
  return nuctoindex(opspec.nucleus(spinsystem())); 
}

void SpinOpGeneratorBase::init_(size_t Mv)
{
  if (verbose_>=MetaFlags::minflagvalue)
    std::cerr << "Warning: SpinOpGenerator initialised with ridiculous verbose level (" << verbose_ << ").  Confused with flags?";
  N_=cstruct.ncells();
  M_=Mv;
  cell_to_spin_.setdims(N_,M_);
}

size_t SpinOpGeneratorBase::nuctoindex(size_t nuc) const
{
  if (nuc==NULL_NUCLEUS)
    throw InvalidParameter("nuctoindex: nucleus not specified");
  const int ind=nuctoindex_(nuc);
  if (ind<0)
    throw InvalidParameter("Nucleus type not present");
  return ind;
}

namespace {
  struct ToMzstate : public ::std::unary_function<size_t,double> {
    ToMzstate(double maxI_) : maxI(maxI_) {}
    size_t operator()(double v) const { return size_t(0.5+maxI-v); }
    double maxI;
  };
}

void SpinOpGenerator::makemzinds(Matrix<size_t>& allmzinds, List<size_t>& mzstates, const BaseList<size_t>& nucs)
{
  const basespin_system& sys(spinsystem());
  const size_t nnucs=nucs.size();
  allmzinds.create(nnucs,sys.size());
  mzstates.create(nnucs);
  for (size_t n=nnucs;n--;) {
    const size_t nuc(nucs(n));
    double maxI=0.0;
    for (size_t i=sys.nspins();i--;) {
      const spin& cspin(sys(i));
      if (cspin.nucleus()==nuc)
	maxI+=(cspin.isrestricted() ? 0.5 : cspin.qn());
    }
    BaseList<size_t> dest(allmzinds.row(n));
    apply(dest,ToMzstate(maxI),libcmatrix::diag_Fz(sys,nuc));
    mzstates(n)=size_t(2*maxI+1.5);
  }
}

namespace {
  template<class F,class T> void
  makeindex(List<size_t>& index, BaseList<int>& rindex,const F& a,size_t len,size_t (T::*func)() const)
  {
    const size_t max=rindex.length();
    rindex=-1;
    
    index.create(max);
    
    size_t types=0;
    
    for (size_t n=len;n--;) {
      const size_t val=(a(n).*func)();
      if (val>max)
	throw Failed("makeindex: index out of range");
      if (rindex(val)<0) {
	rindex(val)=types;
	index(types)=val;
	types++;
      }
    }
    index.resize(types);
  }
}

  /* given set of mz indices (and maximum) and selection of which values to consider (which)
     divide into sets of increasing mz
     algorithm is crude search, but this function not time critical */
namespace {
  void findmzblocks(List< ListList<size_t> >& dest, size_t& submzblocks, const BaseList<int>& mzinds, const BaseList<int>& activeinds)
  {  
    const size_t blks=mzinds.size();
    
    List<size_t> blksizes(blks); blksizes.create(0);
    List<size_t> aslist(blks); aslist.create(0);

    submzblocks=0;
    size_t ptr=0;
    size_t mzind=0;
    while (ptr<blks) {
      size_t cursize=0;
      for (size_t j=0;j<blks;j++) {
	if (mzinds(j)==mzind) {
	  aslist.push_back(j);
	  ptr++;
	  cursize++;
	}
      }
      blksizes.push_back(cursize);
      mzind++;
    }
#ifndef NDEBUG
    std::cout << "aslist: " << aslist << '\n';
    std::cout << "blksizes: " << blksizes << '\n';
#endif
    const size_t mzblks=blksizes.size();
    size_t lastsizebreak=0;
    size_t lastblkbreak=0;
    size_t off=0;
    for (size_t j=0;j<mzblks;j++) {      
      const size_t nextoff=off+blksizes(j);
      if ((j==mzblks-1) || (activeinds(aslist(nextoff))!=activeinds(aslist(off)))) {
	const range selblks(lastblkbreak,nextoff-1);
	const BaseList<size_t> subblks(aslist(selblks));
	const range selsizes(lastsizebreak,j);
	const BaseList<size_t> subsizes(blksizes(selsizes));
	dest.push_back(ListList<size_t>(subsizes,subblks));
	lastblkbreak=nextoff;
	lastsizebreak=j+1;
	if (submzblocks) {
	  if (submzblocks!=subsizes.size())
	    throw InternalError("findmzblocks");
	}
	else
	  submzblocks=subsizes.size();
      }
      off=nextoff;
    }    
    assert(off==blks);
  }

  std::ostream& dumpnucs(const BaseList<size_t>& a, std::ostream& ostr =std::cout)
  {
    for (size_t i=0;i<a.size();i++) {
      if (i)
	ostr << ' ';
      ostr << nuctolabel(a(i));
    }
    return ostr;
  }
}

Warning<> MetaFlags::UseMzSymmetry_heteronuclear_warning("UseMzSymmetry only applicable to homonuclear system",&lcm_base_warning);

void SpinOpGeneratorBase::init_blockstr(size_t neigs, const BaseList<nuclei_spec>& nucs)
{
  const basespin_system& sys(spinsystem());
  nuctoindex_.create(MAX_NUCLEUS);
  //create index and reverse index taking a nucleus id and converting to an index
  makeindex(indextonuc,nuctoindex_,sys,M_,&spin::nucleus);
  const size_t ntypes=indextonuc.size();
  if (verbose_>1) {
    std::cout << "Nuclei found: ";
    dumpnucs(indextonuc) << '\n';
  }

  bool lusemzsym=false;
  if (flags_ & MetaFlags::UseMzSymmetry) {
    if (ntypes>1) 
      MetaFlags::UseMzSymmetry_heteronuclear_warning.raise();
    else
      lusemzsym=true;
  }
  size_=sys.size();

  //Matrix<size_t> allmzinds;
  //List<size_t> maxmzinds;
  makemzinds(allmzinds_,maxmzinds_,indextonuc);
  if (verbose_>1) {
    std::cout << "mz indices matrix\n" << allmzinds_;
    std::cout << "mz states: " << maxmzinds_ << '\n';
  }

  const size_t nblocked(nucs.size());

  List<size_t> unblocked;

  if (nblocked) {
    ispresent_.create(MAX_NUCLEUS,false);  
    blocknuc_.create(nblocked);
    for (size_t k=nblocked;k--;) {
      const size_t nuc=nucs(k)();
      if (ispresent_(nuc))
	throw Failed("SpinOpGeneratorBase: Duplicate blocking nucleus");
      ispresent_(nuc)=true;
      blocknuc_(k)=nuc;
    }
    for (size_t k=ntypes;k--;) {
      const size_t nuc=indextonuc(k);
      if (!ispresent_(nuc))
	unblocked.push_back(nuc);
    }
  }
  else
    unblocked=indextonuc;
  if (verbose_) {
    std::cout << "Unblocked nuclei: ";
    dumpnucs(unblocked) << '\n';
  }
  
  List<int> total_mzind;
  List<int> active_mzind(allmzinds_.cols(),0);
  List<size_t> active_levels(ntypes,1); //!< set to 1 for non-blocked
  List<size_t> inactive_levels(ntypes,1);
  List<size_t> nucinds(ntypes);

  size_t total_multiplier=1;
  size_t active_multiplier=1;
  const size_t nunblocked=unblocked.size();
  for (size_t k=0;k<ntypes;k++) {
    const bool isblocked=(k>=nunblocked);
    const size_t nuc=isblocked ? blocknuc_(k-nunblocked) : unblocked(k);    
    const size_t ind(nuctoindex(nuc));
    nucinds(k)=ind;
    const BaseList<size_t>& curinds(allmzinds_.row(ind));
    mla(total_mzind,total_multiplier,curinds);
    const size_t mulfac=maxmzinds_(ind);
    total_multiplier*=mulfac;
    if (isblocked) {
      mla(active_mzind,active_multiplier,curinds);
      active_levels(ind)=mulfac;
      active_multiplier*=mulfac;
    }
    else
      inactive_levels(ind)=mulfac;
  }
  inactive_blocking_info_.create(inactive_levels,nucinds,false);//!< don't enable mz symmetry for inactive blocking
  active_blocking_info_.create(active_levels,nucinds,lusemzsym);
  if (verbose_>1) {
    std::cout << "total_mzind: " << total_mzind << '\n';
    std::cout << "active_mzind: " << active_mzind << '\n';
    std::cout << "maxmzind: " << maxmzinds_ << '\n';

    std::cout << "Active mz blocking: " << active_blocking_info_;
    std::cout << "Inactive mz blocking: " << inactive_blocking_info_;
  }
  findmzblocks(mzblocking_,submzblocks_,total_mzind,active_mzind);
  if (verbose_)
    std::cout << "Blocking Fz structure: " << mzblocking_ << '\n';
  
  eigstr_=BaseStructure(mzblocking_.size(),neigs,lusemzsym);
}

  SpinOpGenerator::SpinOpGenerator(const basespin_system& sysv, const CrystalStructure& cstructv, int flagsv, int verbosev)
    : SpinOpGeneratorBase(cstructv,sysv.nspins(),flagsv,verbosev),
      sysp_(N_==1 ? &sysv : sysv.clone(N_),N_==1 ? mxflag::nondynamic : mxflag::normal)
{
  common_init();
}

  SpinOpGenerator::SpinOpGenerator(const basespin_system& sysv, int flagsv, int verbosev)
    : SpinOpGeneratorBase(CrystalStructure(),sysv.nspins(),flagsv,verbosev),
      sysp_(&sysv,mxflag::nondynamic)
{
  common_init();
}
  
  SpinOpGenerator::SpinOpGenerator(const basespin_system& sysv, const CrystalStructure& cstructv, const BaseList<nuclei_spec>& blocknucsv, int flagsv, int verbosev)
    : SpinOpGeneratorBase(cstructv,sysv.nspins(),flagsv,verbosev),
      sysp_(N_==1 ? &sysv : sysv.clone(N_),N_==1 ? mxflag::nondynamic : mxflag::normal)   
{ 
  common_init(blocknucsv);
}

  SpinOpGenerator::SpinOpGenerator(const basespin_system& sysv, const BaseList<nuclei_spec>& blocknucsv, int flagsv, int verbosev)
    : SpinOpGeneratorBase(CrystalStructure(),sysv.nspins(),flagsv,verbosev),
      sysp_(&sysv,mxflag::nondynamic)
{ 
  common_init(blocknucsv);
}

SpinOpGenerator::SpinOpGenerator(const basespin_system& sysv, const CrystalStructure& cstructv, const char* label, int flagsv, int verbosev)
  : SpinOpGeneratorBase(cstructv,sysv.nspins(),flagsv,verbosev),
    sysp_(N_==1 ? &sysv : sysv.clone(N_),N_==1 ? mxflag::nondynamic : mxflag::normal)   
{ 
  common_init(ExplicitList<1,nuclei_spec>(label));
}

  SpinOpGenerator::SpinOpGenerator(const HamiltonianStructure& structurev, const CrystalStructure& cstructv, int verbosev)
    : SpinOpGeneratorBase(cstructv,structurev,verbosev),
    sysp_( (N_==1) ? &(structurev.spinsystem()) : structurev.spinsystem().clone(N_),(N_==1) ? mxflag::nondynamic : mxflag::normal)
  {
    common_init(structurev.blockingnuclei());
  }
   
  void SpinOpGenerator::proton_frequency(double pfreqv)
  {
    if (pfreqv==0.0)
      throw InvalidParameter("SpinOpGenerator::proton_frequency");
    if (pfreqv!=pfreq_) {
      pfreq_=pfreqv;
      update();
    }
  }

double SpinOpGeneratorBase::larmor(size_t j) const
{
  static const double gamma1H=gamma("1H");
  if (!pfreq_)
    throw Failed("proton frequency must be set for non-secular Hamiltonians");    
  const double v=pfreq_*gamma(spinsystem()(j).nucleus())/gamma1H;
  assert(v!=0.0);
  return v;
}

  void SpinOpGenerator::update()
  {
    if (effective_quadrupole_order()==1)
      return;

    for (size_t j=nspins_cell();j--;) {
      if (structurep_->nonsecular(j))
	mla_Iz(Hzeeman_,larmor(j),j);
    }
    if (verbose()>1)
      std::cout << "Hzeeman: " << Hzeeman_ << " Hz\n";
    
    if ((quadrupole_order()==2) && zeemanstr_.empty()) {
      List<double> tmp(Hzeeman_.items());
      size_t n=Hzeeman_.size();
      zeemanstr_.create(n);
      for (;n--;) {
	BaseList<double> lHzeeman(Hzeeman_(n));
	ListList<size_t>& zeemanstrn(zeemanstr_(n));
	zeemanstrn=find_blocks(lHzeeman,1e-5);
	tmp=lHzeeman(zeemanstrn.row());
	lHzeeman=tmp;
      }
      if (verbose()>1)
	std::cout << "Hzeeman (sorted): " << Hzeeman_ << " Hz\n";
    }
  }

void SpinOpGeneratorBase::init_partitioning()
{
  const ListList<size_t>& diagstr(partitioned_diagonal_structure());
  if (diagstr.empty())
    throw InternalError("init_partitioned: called before partitioned diagonal structure has been created");
  if (verbose_>1)
    std::cout << "Partitioned diagonal structure: " << diagstr << '\n';

  bool nontrivial=false;
  for (size_t i=0;i<diagstr.size();i++) {
    if (diagstr(i).size()>1)
      nontrivial=true;
    partitions_.push_back(diagonal_partition(diagstr(i)));
  }
  if (!nontrivial) {
    if (verbose_>1)
      std::cout << "init_partitioning: partitioning is trivial\n";
    partitions_.clear();
  }
}

Warning<> MetaFlags::UseMzSymmetry_oddblocks_warning("UseMzSymmetry  may not work correctly with odd number of blocks in the detect operator",&lcm_base_warning);

void SpinOpGenerator::common_init(const BaseList<nuclei_spec>& nucs)
{
  if (flags() & MetaFlags::UseMzSymmetry)
    MetaFlags::UseMzSymmetry_oddblocks_warning.raise();

  init_blockstr(1,nucs);
  pfreq_=0.0;

  size_t actblocks=mzblocking_.size();
  sizes_.create(actblocks);
  for (;actblocks--;)
    sizes_(actblocks)=mzblocking_(actblocks).items();

  assert(partitioned_diagstr_.empty());
  if (flags() & MetaFlags::UsePartitioning) {
    List<size_t> sizes;
    for (size_t i=0;i<mzblocking_.size();i++) {
      const ListList<size_t>& curblk(mzblocking_(i));
      const size_t n=curblk.size();
      sizes.create(n);
      for (size_t j=n;j--;)
	sizes(j)=curblk.size(j);    
      partitioned_diagstr_.push_back(sizes);
    }
    init_partitioning();
  }
  if (!isblocked()) {//!< reset mzblocking to natural order (29/11/09)
    BaseList<size_t> dest(mzblocking_.front().row());
    for (size_t i=dest.size();i--;)
      dest(i)=i;
  }

  cstructiterp.reset(cstruct.create_iterator(M_));

  if (isclassicsecondorder()) {
    //if (nspins()!=1)
    //  throw Failed("Spin system must involve single spin for classic 2nd order treatment");    
    //! build second order shifts
    size_t n=nspins_cell();
    secondordershifts1_.create(n);
    secondordershifts2_.create(n);    
    List<double> tmp;
    for (;n--;) {
      if ((*sysp_)(n).deg()>2) {
	tmp=diag_spin_quadrupolar_Cq(*sysp_,n,1);
	chop_mla(secondordershifts1_(n),1.0,tmp);
	tmp=diag_spin_quadrupolar_Cq(*sysp_,n,2);
	chop_mla(secondordershifts2_(n),1.0,tmp);
      }
    }
  }
}

template<typename M> void SpinOpGenerator::create(BlockedMatrix<M>& dest, const block_pattern& blkspec) const
{
  size_t r,c;
  if (!isblocked()) {
    if (!blkspec.isdiagonal())
      throw Failed("SpinOpGenerator: off-diagonal block in unblocked matrix");
    r=size(0U,0U);
    dest.create(BaseList<size_t>(1,&r));
  }
  else {
    if (blkspec.blocks==0)
      throw Undefined("SpinOpGenerator::create: block specification is empty");
    block_pattern::iterator iter(blkspec);
    ScratchList<size_t> rstr(blkspec.blocks);
    ScratchList<size_t> cstr(blkspec.blocks);
    size_t k;
    while (iter.next(r,c,k)) {
      rstr(k)=blockindices(r).size();
      cstr(k)=blockindices(c).size(); 
    } 
    dest.create(rstr,cstr);
  }
}

void SpinOpGenerator::rawdiag_Fz(ListList<double>& dest, size_t nuc) const
{
  List<double> tmpr;
  const basespin_system& sys(spinsystem());
  const size_t Ncell=nspins_cell();
  for (size_t i=Ncell;i--;) {
    if (sys(i).nucleus()==nuc) {
      for (size_t n=i;n<nspins();n+=Ncell)
	sysp_->mla_Iz(tmpr,1.0,n);
    }
  }
  chop_mla(dest,1.0,tmpr);
}

static const char NOTSIMPLE[]="SpinOpGenerator:can't store blocked Hamiltonian in simple matrix";

template<class T> class IsGreater
{
public:
  IsGreater(double tolv)
    : tol_(tolv) {}
  bool operator()(T v) { return (fabs(v)>tol_); }
private:
  double tol_;
};

template<> class IsGreater<complex>
{
public:
  IsGreater(double tolv)
    : tol_(tolv) {}
  bool operator()(const complex& v) { return (fabs(v.real())>tol_) || (fabs(v.imag())>tol_); }
private:
  double tol_;
};


template<class T> size_t countnz(const BaseList<T>& a, double tol =1e-8)
{
  IsGreater<T> checker(tol);

  size_t count=0;
  for (size_t i=a.size();i--;) {
    if (checker(a(i)))
      count++;
  }
  return count;
}

Warning<> losing_elements_warning("Lost elements on blocking",&lcm_base_warning);

template<class M1,class M2> void SpinOpGenerator::chop_mla(Matrix<M1>& dest, double coup, const Matrix<M2>& source) const
{
  if (coup!=1.0)
    ::libcmatrix::mla(dest,coup,source);
  else
    dest+=source;
}

  template<class M1,class M2> void SpinOpGenerator::chop_mla(BlockedMatrix<M1>& dest, double coup, const Matrix<M2>& source,const block_pattern& blkspec) const
  {
    const bool neednew=!dest;
    if (neednew) 
      create(dest,blkspec);

    if (!isblocked()) {
      if (neednew)
 	multiply(dest.front(),coup,source);
      else
 	::libcmatrix::mla(dest.front(),coup,source);
      return;
    }

    Matrix<M2> tmp;
    size_t r,c,k;
    block_pattern::iterator iter(blkspec);
    while (iter.next(r,c,k)) {
      if (verbose_>1)
	std::cout << "Creating block " << k << " (row mz eig: " << r << ", col mz eig: " << c << ") from indices R" << blockindices(r) << " C" << blockindices(c) << '\n';
      if (neednew) {
	dest(k)=source(blockindices(r),blockindices(c));
	if (coup!=1.0)
	  dest(k)*=coup;
      }
      else {
	tmp=source(blockindices(r),blockindices(c));
	::libcmatrix::mla(dest(k),coup,tmp);
      }
    }
    if (neednew && CHECKCOND) {
      const size_t tnz=countnz(dest.row());
      const size_t startnz=countnz(source.row());
      if (tnz!=startnz) {
	char error[256];
	sprintf(error," final: %" LCM_PRI_SIZE_T_MODIFIER "u vs initial: %" LCM_PRI_SIZE_T_MODIFIER "u [hook in chop_mla at line %i of %s]",tnz,startnz,__LINE__,__FILE__);
	losing_elements_warning.raise(error);
      }
    }
  }

static void cut(MultiMatrix<double,3>& tmp,const MultiMatrix<double,3>& source, const BaseList<size_t>& blkstr)
{
  const size_t depth(source.dimension(0U));
  const size_t blksize(blkstr.size());
  tmp.create(depth,blksize,blksize);
  for (size_t j=depth;j--;)
    tmp(j)=source(j)(blkstr,blkstr);
}

void SpinOpGenerator::chop_ns(List< MultiMatrix<double,3> >& dest, const MultiMatrix<double,3>& source) const
{
  if (source.empty())
    throw Undefined("chop_ns");
  if (!isblocked()) {
    dest.create(1);
    dest.front()=source;
    return;
  }
  create_ns(dest);  
  for (size_t k=actual_mzblocks();k--;)
    cut(dest(k),source,blockindices(k));

  if (CHECKCOND) {
    size_t tnz=0;
    for (size_t k=actual_mzblocks();k--;)
      tnz+=countnz(dest(k).row());
    const size_t startnz=countnz(source.row());
    if (tnz!=startnz) {
      char error[256];
      sprintf(error," final: %" LCM_PRI_SIZE_T_MODIFIER "u vs initial: %" LCM_PRI_SIZE_T_MODIFIER "u [hook in chop_ns at line %i of %s]",tnz,startnz,__LINE__,__FILE__);
      losing_elements_warning.raise(error);
    }
  }
}

template<class M> void SpinOpGenerator::chop_mla(Matrix<M>& dest, double scale, const List<double>& source) const
{
  if (actual_mzblocks()>1)
    throw Failed(NOTSIMPLE);
  if (scale!=1.0) {
    if (isblocked()) {
      List<double> tmp(source(blockindices(0U)));
      ::libcmatrix::mla(dest,scale,tmp);
    }
    else
      ::libcmatrix::mla(dest,scale,source);
  }
  else {
    if (isblocked())
      dest+=source(blockindices(0U));
    else
      dest+=source;
  }
}

  template<class M> void SpinOpGenerator::chop_mla(BlockedMatrix<M>& dest, double scale, const List<double>& source) const
  {
    const bool neednew=!dest;
    if (neednew) 
      create(dest);

    List<double> tmp(128);
    for (size_t k=actual_mzblocks();k--;) {
      if (!isblocked())
	tmp=source;
      else
	tmp=source(blockindices(k));

      if (scale!=1.0)
	tmp*=scale;

      if (neednew)	
	::libcmatrix::full(dest(k),tmp);
      else
	dest(k)+=tmp;
    }
  }

void SpinOpGeneratorBase::print(std::ostream& ostr) const
{
  ostr << "Spin system: " << spinsystem() << '\n';
  ostr << eigstr_;
  const size_t n(blocknuc_.size());
  if (n==0) {
    ostr << "No mz blocking\n";
    return;
  }
  ostr << "Block structure: " << mzblocking_ << '\n';
//   for (size_t nuci=0;nuci<n;nuci++) {
//     const size_t nuc(blocknuc_(nuci));
//     ostr << nuctolabel(nuc) << "  Active step: " << active_mzstep_(nuctoindex(nuc)) << "  Total step: " << total_mzstep_(nuctoindex(nuc)) << '\n';
//   }
  ostr << "Diagonal partioning: " << partitions_ << '\n';
  if (!Hzeeman_.empty())
    ostr << "Zeeman Hamiltonian: " << Hzeeman_ << " Hz\n";
  ostr << "Proton frequency: ";
  if (proton_frequency())
    ostr << (proton_frequency()*1e-6) << " MHz\n";
  else
    ostr << "unset\n";
  const size_t indexbase=cmatrix_ostream_controller(ostr).indexbase;
  for (size_t i=0;i<secondordershifts1_.size();i++) {
    if (!secondordershifts1_(i).empty()) {
      ostr << "Second order shift for spin " << (i+indexbase) << " n=1: " << secondordershifts1_(i) << '\n';
      ostr << "Second order shift for spin " << (i+indexbase) << " n=2: " << secondordershifts2_(i) << '\n';
    }
  }
}

template<class C,class M> void SpinOpGenerator::chop_mla(RealSpinningHamiltonian& dest, const C& coup, const M& source) const 
  {
    if (actual_mzblocks()>1)
      throw Failed(NOTSIMPLE);
    dest.add(coup,source);
  }

template<class C,class M> void SpinOpGenerator::chop_mla(SpinningHamiltonian& dest, const C& coup, const M& source) const 
  {
    if (actual_mzblocks()>1)
      throw Failed(NOTSIMPLE);
    dest.add(coup,source);
  }

template<class C,class M> void SpinOpGenerator::chop_mla(DiagonalSpinningHamiltonian& dest, const C& coup, const M& source) const 
  {
    if (actual_mzblocks()>1)
      throw Failed(NOTSIMPLE);
    dest.add(coup,source);
  }

  template<class M,class C> void SpinOpGenerator::chop_mla(List<M>& dest, const C& coup, const List<double>& source) const 
  {
    if (!isblocked()) {
      dest.front().add(coup,source);
      return;
    }
    List<double> tmp(128);
    for (size_t k=actual_mzblocks();k--;) {
      tmp=source(blockindices(k));
      dest(k).add(coup,tmp);
    }
  }

  template<class C> void SpinOpGenerator::chop_mla(List<double>& dest, const C& coup, const List<double>& source) const 
  {
    if (!isblocked())
      ::libcmatrix::mla(dest,coup,source);
    else
      throw Failed("SpinOpGenerator: no block structure permitted");
  }

  template<class M,class C> void SpinOpGenerator::chop_mla(List<M>& dest, const C& coup, const rmatrix& source) const 
  {
    if (!isblocked()) {
      dest.front().add(coup,source);
      return;
    }
    rmatrix tmp;
    for (size_t k=actual_mzblocks();k--;) {
      tmp=source(blockindices(k),blockindices(k));
      dest(k).add(coup,tmp);
    }
  }

void SpinOpGenerator::add_ns_(List<SpinningHamiltonian>& dest, const space_T& coup, const BaseList< MultiMatrix<double,3> >& source) const 
  {
    for (size_t k=actual_mzblocks();k--;)
      dest(k).add(coup,source(k));
  }

void SpinOpGenerator::add_ns_(SpinningHamiltonian& dest, const space_T& coup, const BaseList< MultiMatrix<double,3> >& source) const 
  {
    if (actual_mzblocks()>1)
      throw Failed(NOTSIMPLE);
    dest.add(coup,source.front());
  }

static const double sqrtfac=-std::sqrt(0.5);

void SpinOpGenerator::add_ns_(BlockedMatrix<complex>& dest, const space_T& coup, const BaseList< MultiMatrix<double,3> >& source) const 
  {
    if (dest.empty()) {
      create(dest);
      dest=complex(0.0);
    }
    const int L=getL(source.front().dimension(0U)); //validateL(coup,source.front().dimension(0U));
    //    if (L) {
    assert(L<=2);
    for (size_t k=actual_mzblocks();k--;) {
      const MultiMatrix<double,3>& sourcek(source(k));
      cmatrix& destk(dest(k));
      if (coup.have_rank(2)) {
	for (int m=-L;m<=L;m++)      
	  ::libcmatrix::mla(destk,coup(2,m),sourcek(m+L));
      }
      if (coup.have_rank(0))
	::libcmatrix::mla(destk,sqrtfac*real(coup(0,0)),sourcek(L));
    }
  }

void SpinOpGenerator::add_ns_(Matrix<complex>& dest, const space_T& coup, const BaseList< MultiMatrix<double,3> >& source) const 
  {
    if (actual_mzblocks()>1)
      throw Failed(NOTSIMPLE);
    const MultiMatrix<double,3>& source0(source.front());
    const int L=getL(source0.dimension(0U));//validateL(coup,source0.dimension(0U));
    //if (L) {
    assert(L<=2);
    if (coup.have_rank(2)) {
      for (int m=-L;m<=L;m++)      
	::libcmatrix::mla(dest,coup(2,m),source0(m+L));
    }
    if (coup.have_rank(0))
      ::libcmatrix::mla(dest,sqrtfac*real(coup(0,0)),source0(L));
    //    }
  }

  void SpinOpGenerator::chop_mla(ListList<double>& dest, double scale, const List<double>& source) const
  { 
    const bool neednew=dest.empty();
    if (neednew)
      create(dest);

    List<double> tmp(128); //pre-allocate space for variable sized temporary
    if (!isblocked()) {
      BaseList<double> destk(dest.front());
      if (neednew)
	multiply(destk,scale,source);
      else
	::libcmatrix::mla(destk,scale,source);
    }
    else {
      for (size_t k=actual_mzblocks();k--;) {
	BaseList<double> destk(dest(k));
	const BaseList<size_t> indicesk(blockindices(k));
	if (neednew) {
	  destk=source(indicesk);
	  if (scale!=1.0)
	    destk*=scale;
	}
	else {
	  tmp=source(indicesk);
	  ::libcmatrix::mla(destk,scale,tmp);
	}
      }
    }
  }

void SpinOpGenerator::add(BlockedMatrix<complex>& dest, block_pattern& blkspec, const productoperator_spec& prodspec) const
  {
    static const complex iconst(0.0,1.0);
    blkspec=block_pattern(*this,prodspec); //will fail if non-trivial & blocked
    cmatrix tmp;
    
    for (size_t i=prodspec.size();i--;) {
      const BaseList<operator_spec> ospec(prodspec.specs(i));

      const complex cscale(prodspec.scales(i));    
      bool isimag=false;
      cmatrix* wspacep=&tmp;
      double rscale;
      cmatrix tmpr;
      if (imag(cscale)) {
	if (real(cscale))
	  throw Failed(PUREERR);
	rscale=imag(cscale);
	isimag=true;
	wspacep=&tmpr;
      }
      else
	rscale=real(cscale);
      
      switch (ospec.size()) {
      case 1:
	mla_(*wspacep,rscale,ospec.front());
	if (isimag)
	  ::libcmatrix::mla(tmp,iconst,tmpr);
	break;
      case 2:
	mla_(*wspacep,rscale,ospec(0U),ospec(1U));
	if (isimag)
	  ::libcmatrix::mla(tmp,iconst,tmpr);
	break;
      default: {
	cmatrix tmp1,tmp2,tmp3;
	for (size_t j=ospec.size();j--;) {
	  if (!!tmp2)
	    tmp2=complex(0.0);
	  mla_(tmp2,1.0,ospec(j));
	  if (tmp1.empty())
	    tmp1.swap(tmp2);
	  else {
	    multiply(tmp3,tmp1,tmp2);
	    tmp1.swap(tmp3);
	  }
	}
	::libcmatrix::mla(tmp,cscale,tmp1);
      }
      }
    }
    chop_mla(dest,1.0,tmp,blkspec);
  }      

void SpinOpGenerator::mla(BlockedMatrix<complex>& dest, block_pattern& blkspec, double scale,const operator_spec& opspec) const
{
  if ((opspec.nuc==NULL_NUCLEUS) && (opspec.number>=nspins_cell()))
    throw InvalidParameter("SpinOpGenerator::mla restricted to unit cell for periodic problems");
  blkspec=block_pattern(*this,opspec);
  if (verbose_>1)
    std::cout << "Block structure: " << blkspec << '\n';
  cmatrix tmp;
  mla_(tmp,scale,opspec);
  chop_mla(dest,scale,tmp,blkspec);
}

void SpinOpGenerator::mla_(cmatrix& dest,double scale,const operator_spec& opspec) const
{
  if (opspec.nuc==NULL_NUCLEUS) {
    for (size_t n=opspec.number;n<nspins();n+=nspins_cell()) {
      if (opspec.issingletransition())
	sysp_->mla_I(dest,scale,n,opspec.op,opspec.col);
      else
	sysp_->mla_I(dest,scale,n,opspec.op);
    }
  }
  else {
    if (opspec.issingletransition())
      sysp_->mla_F(dest,scale,opspec.nuc,opspec.op,opspec.col);
    else
      sysp_->mla_F(dest,scale,opspec.nuc,opspec.op);  
  }
}
  
void SpinOpGenerator::mla_(cmatrix& dest,double scale,const operator_spec& opspec1, const operator_spec& opspec2) const
  {
    const size_t ncell(nspins_cell());
    if ((opspec1.number>=ncell) || (opspec2.number>=ncell))
      throw InvalidParameter("SpinOpGenerator::mla restricted to unit cell for periodic problems");
    assert((opspec1.nuc==NULL_NUCLEUS) && (opspec2.nuc==NULL_NUCLEUS));

    size_t n1=opspec1.number;
    size_t n2=opspec2.number;

    const int diff=n2-n1;
    const size_t totspins=nspins();
    for (size_t cellbase=0;cellbase<ncells();cellbase++) {
      const size_t finali=cell_to_spin(cellbase,n1);
      size_t finalj=finali+diff;
      if (finalj>=totspins)
	finalj-=totspins;
      if (finali<finalj)
	sysp_->mla_I(dest,scale,finali,opspec1.op,finalj,opspec2.op);
    }
  }
  
  void SpinOpGenerator::mla_Fz(ListList<double>& dest,double scale,nuclei_spec whichn) const
  {
    ::libcmatrix::mla(dest,scale,diag_Fz(whichn));
  }
  
  void SpinOpGenerator::mla_Iz(ListList<double>& dest,double scale,size_t n) const
  {
     List<double> tmp;
     sysp_->mla_Iz(tmp,1.0,n);
     chop_mla(dest,scale,tmp);
  }
  
  template<class T,class C> void
  SpinOpGenerator::addrawCS(T& dest, const C& coup, size_t j) const
  {
    if (N_==1)
      chop_mla(dest,coup,diag_spin_CS(*sysp_,j));
    else {
      List<double> Hraw;
      for (size_t n=N_;n--;)
	Hraw+=diag_spin_CS(*sysp_,cell_to_spin(n,j));
      chop_mla(dest,coup,Hraw);
    }
  }

  template<class T,class C> void
  SpinOpGenerator::addrawQ(T& dest, const C& coup, size_t j) const
  {
    if (N_==1)
      chop_mla(dest,coup,diag_spin_quadrupolar_Cq(*sysp_,j));
    else {
      List<double> Hraw;
      for (size_t n=N_;n--;)
	Hraw+=diag_spin_quadrupolar_Cq(*sysp_,cell_to_spin(n,j));
      chop_mla(dest,coup,Hraw);
    }
  }

void SpinOpGenerator::makeraw_ns(List< MultiMatrix<double,3> >& dest, size_t i, size_t j, ns_flag nstype) const 
{  
  MultiMatrix<double,3> Hraw,Htmp;
  if (i==j) {
    for (size_t n=0;n<N_;n++) {
      const size_t sk(cell_to_spin(n,j));
      if (n) {
	spin_A2_ns(Htmp,*sysp_,sk,sk,nstype);
	Hraw+=Htmp;
      }
      else
	spin_A2_ns(Hraw,*sysp_,sk,sk,nstype);
    }
  }
  else {
    
    const size_t cells(ncells());
    const size_t totspins(nspins());
    const int diff=j-int(i);
    for (size_t cellbase=0;cellbase<cells;cellbase++) {
      const int finali=cell_to_spin(cellbase,i);
      int finalj=finali+diff;
      if (finalj>=totspins)
	finalj-=totspins;
      else {
	if (finalj<0)
	  finalj+=totspins;
      }
      // if (finali<finalj) {
      spin_A2_ns(Htmp,*sysp_,finali,finalj,nstype);
      Hraw+=Htmp;
      //      }
    }
  }
  if (!Hraw)
    throw InternalError("makeraw_ns");
  chop_ns(dest,Hraw);
}

template<class T> void
SpinOpGenerator::add_ns_(T& dest, const space_T& coup, size_t j, size_t sk, ns_flag nstype, Bool2Type<false>) const
{
  if (flags_ & MetaFlags::DontCacheBinary) {
    List< MultiMatrix<double,3> > tmp;
    makeraw_ns(tmp,j,sk,nstype);
    add_ns_(dest,coup,tmp);
    return;
  }
  if ((nstype==NS_BOTH) && (j>sk))
    ::std::swap(j,sk); //couplings must be symmetrical; increase chances of cache hit

  MutexLock<ThreadingActive> autolock(lock);

  if (A2store.empty())	
    A2store.create(nspins_cell(),nspins()); 

  List< MultiMatrix<double,3> >& use(A2store(j,sk).set(Int2Type<BLOCKEDTENSOR>()));
  if (use.empty())
    makeraw_ns(use,j,sk,nstype);
  add_ns_(dest,coup,use);
}

template<class T,class M> void SpinOpGenerator::add_(Matrix<T>& dest, double coup, const M& a) const
{
  if (actual_mzblocks()>1)
    throw Failed(NOTSIMPLE);
  ::libcmatrix::mla(dest,coup,a.front());
}

  template<class T> void
  SpinOpGenerator::add_Hcs(T& dest, const typename Ham_traits<T>::coupling_type& coup, size_t j) const
  {
    MutexLock<ThreadingActive> autolock(lock);
    if (CSstore.empty())
	CSstore.create(nspins()/ncells());
    ListList<double>& use(CSstore(j));
    if (use.size()==0) 
      addrawCS(use,1.0,j);
    add_(dest,coup,use);
  }
  
template<class In,class Out> void add_real_(Out& out, const In& in, Out& tmp) { 
  real(tmp,in);
  out+=tmp;
}

template<class T> void add_real_(T& out, const T& in, T&) {
  out+=in;
}

  template<class T> struct lcm_real_type_ { typedef Matrix<double> type; };
  template<> struct lcm_real_type_<ListList<double> > { typedef List<double> type; };

  template<class T,class F,class C,class ResultType> void
  SpinOpGenerator::addraw(T& dest, F func, const C& coup, size_t i, size_t j,Type2Type<ResultType>) const
  { 
    MutexLock<ThreadingActive> locallock(iterlock);
    CrystalStructure_iterator& iter(*cstructiterp);
    iter.reset(i,j);
    typename lcm_real_type_<ResultType>::type Hraw,Htmpr;
    size_t li,lj;
    while (iter.next(li,lj)) {
      if (li<lj) {
	add_real_(Hraw,(*func)(*sysp_,li,lj),Htmpr);
	if (verbose()>1)
	  std::cout << "Adding coupling between spin " << li << " and " << lj << '\n';
      }
    }
    chop_mla(dest,coup,Hraw);
  }

template<class T,class F,class ResultType> void
SpinOpGenerator::add_A0(T& dest,F func, double coup, size_t j, size_t sk, Type2Type<ResultType> restype) const
{
  if (flags_ & MetaFlags::DontCacheBinary) {
    addraw(dest,func,coup,j,sk,restype);
    return;
  }
  if (j>sk)
    ::std::swap(j,sk); //couplings must be symmetrical; increase chances of cache hit
  
  MutexLock<ThreadingActive> autolock(lock);
  if (A0store.empty())
    A0store.create(nspins()/ncells(),nspins()); 
  ResultType& use(A0store(j,sk).set(restype));
  if (use.empty())
    addraw(use,func,1.0,j,sk,restype);
  add_(dest,coup,use);
}

  //nasty bit of trickery to make isotropic component disappear ...
  template<typename T> struct fudge_isotropic_;
  template<> struct fudge_isotropic_<space_T>  {
    fudge_isotropic_(const space_T& A) : ref_(const_cast<space_T&>(A))
      { ref_.clear(0); }
    const space_T& operator()() const { return ref_; }
    ~fudge_isotropic_() { ref_.ensure_rank(0); }
    space_T& ref_;
  };
  template<> struct fudge_isotropic_<double> { 
    fudge_isotropic_(double v) : v_(v) {};
    double operator()() const { return v_; }
    double v_;
  };

  template<class T,class F,class C,class ResultType> void
SpinOpGenerator::add_A2(T& dest,F func, const C& coupall, size_t j, size_t sk, Type2Type<ResultType> restype) const
{
  fudge_isotropic_<C> coup(coupall);
  if (verbose()>1)
    std::cout << "Adding base coupling between " << j << " and " << sk << ": " << coup() << '\n';
  if (flags_ & MetaFlags::DontCacheBinary) {
    addraw(dest,func,coup(),j,sk,restype);
    return;
  }
  if (j>sk)
    ::std::swap(j,sk); //couplings must be symmetrical; increase chances of cache hit
  
  MutexLock<ThreadingActive> autolock(lock);
  if (A2store.empty())
    A2store.create(nspins()/ncells(),nspins()); 
  ResultType& use(A2store(j,sk).set(restype));
  if (use.empty())
    addraw(use,func,1.0,j,sk,restype);  
  add_(dest,coup(),use);
}

template<class T> void 
SpinOpGenerator::add_Hquadrupolar2(T& dest, const space_T& V, size_t i) const
{
  if (!(V.have_rank(2)))
    return;
  
  if (secondordershifts1_.empty() || secondordershifts1_(i).empty())
    throw InternalError("add_Hquadrupolar2");

  reducer_< typename Ham_traits<T>::coupling_type, space_T> reduce;
  add_Hquadrupolar(dest,reduce(V),i);

  // Adds (1/2 v_0) [ V(2,1)*V(2,-1)* (4*I^2 - 8 I_z^2 -1)I_z + V(2,2)*V(2,-2)* (2I^2 - 2I_z^2 -1)I_z ]
  // Cq factor is included in V
  // Use symmetry between +/- tensor elements to simplify
  const double scalef=3.0/larmor(i);
  const double shift2=scalef*norm(V(2,2));
  const double shift1=-scalef*norm(V(2,1));
  if (verbose()>1)
    std::cout << "Adding second order corrections: " << shift1 << " Hz (n=1) and " << shift2 << " Hz (n=2)\n";
  add_(dest,shift1,secondordershifts1_(i));
  add_(dest,shift2,secondordershifts2_(i));
}

// void BlockedDiagonalSpinningHamiltonianBase::add0_(double coeff, const ListList<double>& A)
// {
//   RotorHolder_<DiagonalSpinningHamiltonian>& dest(*this);
//   for (size_t i=size();i--;)
//     mla(dest(i).component0(),coeff,A(i));
// }

// void BlockedDiagonalSpinningHamiltonianBase::add_(size_t n, const complex& coeff, const ListList<double>& A)
// {
//   RotorHolder_<DiagonalSpinningHamiltonian>& dest(*this);
//   for (size_t i=size();i--;) {
//     BaseList<complex> cur(dest(i).component(n));
//     mla(cur,coeff,A(i));
//   }
// }

void SpinOpGenerator::addQ2_(List<DiagonalSpinningHamiltonian>& dest, size_t n, double scalef, const space_T& V, const ListList<double>& ops) const
{
  tmp_.create(5);
  coeffs_.create(5);

  const cmatrix& d2s(dest.front().info().d2values());
  for (int m=5;m--;)
    tmp_(m)=d2s(m,n+2)*V(2,m-2);

  coeffs_=complex(0.0);
  for (int m=5;m--;) {
    for (int mp=5;mp--;) {
      int diffm=m-mp;
      if (diffm>=0)
	conj_mla(coeffs_(diffm),tmp_(mp),tmp_(m));
    }
  }
  coeffs_*=scalef;
  for (size_t m=5;m--;) {
    if (verbose()>1)
      std::cout << "Second order correction from n=" << n << " to rank " << m << ": " << coeffs_(m) << " Hz\n";
    for (size_t i=dest.size();i--;) {
      if (m) {
	BaseList<complex> cur(dest(i).component(m));
	::libcmatrix::mla(cur,coeffs_(m),ops(i));
      }
      else {
	BaseList<double> cur(dest(i).component0());
	::libcmatrix::mla(cur,real(coeffs_.front()),ops(i));
      }
    }
  }
}

void SpinOpGenerator::add_Hquadrupolar2(List<DiagonalSpinningHamiltonian>& dest, const space_T& V, size_t i) const
{
  if (secondordershifts1_.empty() || secondordershifts1_(i).empty())
    throw InternalError("add_Hquadrupolar2");

  add_Hquadrupolar(dest,V,i);
  
  const double scalef=3.0/larmor(i); //1.0/(larmor(i)*3.0);
  addQ2_(dest,1,-scalef,V,secondordershifts1_(i));
  addQ2_(dest,2,scalef,V,secondordershifts2_(i));
}

template<class T> void
SpinOpGenerator::add_Hquadrupolar(T& dest, const typename Ham_traits<T>::coupling_type& coup, size_t j) const
{
  MutexLock<ThreadingActive> autolock(lock);
  if (A2store.empty())
    A2store.create(nspins()/ncells(),nspins()); 
  ListList<double>& use(A2store(j,j).set(Type2Type< ListList<double> >()));
  if (use.empty())
    addrawQ(use,1.0,j);
  add_(dest,coup,use);
}

  void SpinOpGenerator::print(std::ostream& ostr) const
  {
    SpinOpGeneratorBase::print(ostr);
    spydefined(ostr,"Shifts: ",CSstore);
    spydefined(ostr,"Rank 2 spin ops:\n",A2store);
    spydefined(ostr,"Rank 0 spin ops:\n",A0store);
  }
  
  void SpinOpGenerator::spydefined(std::ostream& ostr, const char* name, const Matrix<cache_type>& a)
  {
    if (!a)
      return;
    ostr << name;
    for (size_t i=0;i<a.rows();i++) {
      for (size_t j=0;j<a.cols();j++)
	ostr << (!a(i,j) ? '.' : 'X');
      ostr << '\n';
    }
  }
  
  void SpinOpGenerator::spydefined(std::ostream& ostr, const char* name, const List< ListList<double> >& a)
  {
    if (a.empty())
      return;
    ostr << name;
    for (size_t i=0;i<a.size();i++)
      ostr << (a(i).empty() ? '.' : 'X');
    ostr << '\n';
  }

//Add 1 to indices
template<class T> bool HamiltonianStore<T>::verify_(const couplingstore<T>& coups, size_t j, size_t k, size_t swapj, size_t swapk,const char* name, std::ostream& ostr, double tol)
{
  if (coups.areequal(j,k,swapj,swapk,tol))
    return true;
  ostr << "Incompatible/missing coupling: " << name << ' ' << (j+1) << ' ' << (k+1) << " & " << name << ' ' << (swapj+1) << ' ' << (swapk+1) << '\n';
  return false;
}

size_t CrystalStructure::inverse(size_t cell0, size_t cell1, size_t cell2) const
{
  const size_t swapcell0(inverse(cell0,Int2Type<0>()));
  const size_t swapcell1(inverse(cell1,Int2Type<1>()));
  const size_t swapcell2(inverse(cell2,Int2Type<2>()));
  return (*this)(swapcell0,swapcell1,swapcell2);
}

template<class T> bool HamiltonianStore<T>::verify_(std::ostream& ostr, const couplingstore<T>& coups, const char* name, const CrystalStructure& cstruct, double tol) const
{
  bool ok=true;
  for (int cell0=cstruct.dimension(0)/2;cell0>=0;cell0--) {
    for (int cell1=cstruct.dimension(1)/2;cell1>=0;cell1--) {
      for (int cell2=cstruct.dimension(2)/2;cell2>=0;cell2--) {
	const size_t cell=cstruct(cell0,cell1,cell2);
	if (cell==0) //not couplings *within* unit cell
	  break;
	const size_t swapcell=cstruct.inverse(cell0,cell1,cell2);
	const bool isswap=(cell==swapcell);
	for (size_t j=0;j<M_;j++) {
	  for (size_t k=0;k<M_;k++) {
	    if (!isswap || (j>k)) {
	      if (!verify_(coups,
			   j,cell_to_spin(cell,k),
			   k,cell_to_spin(swapcell,j),
			   name,ostr,tol))
		ok=false;
	    }
	  }
	}
      }
    }
  }
  return ok;
}

template<class T> bool HamiltonianStore<T>::isequal(const HamiltonianStore<T>& a, double tol) const 
{
  if (!::libcmatrix::areequal(qstore,a.qstore,tol))
    return false;
  if (qorder!=a.qorder)
    return false;
  if (!areequal(couplingmap,a.couplingmap,tol))
    return false;
  if (!areequal(shiftmap,a.shiftmap,tol))
    return false;
  return true;
}

template<class T> template<typename Key,typename Obj> bool HamiltonianStore<T>::areequal(const std::map<Key,Obj>& a, const std::map<Key,Obj>& b, double tol) 
{
  if (a.size()!=b.size())
    return false;
  const typename std::map<Key,Obj>::const_iterator iterend(a.end());
  typename std::map<Key,Obj>::const_iterator iter(a.begin());
  typename std::map<Key,Obj>::const_iterator iter2(b.begin());
  while (iter!=iterend) {
    if (( (iter->first) !=(iter2->first)) || !areequal(iter->second,iter2->second,tol))
      return false;
    ++iter;
  }
  return true;
}

  template<class T> void couplingstore<T>::apply(const couplingstore& a, const Permutation& perm)
  {
    if (!!(*this))
      throw Failed("couplingstore::apply"); //dest must be empty
    if (!!a.A0)
      perm.apply(A0,a.A0);
    if (!!a.A2)
      perm.apply(A2,a.A2);
  }

template<class T> template<typename Key,typename Obj> void HamiltonianStore<T>::apply(std::map<Key,Obj>& dest, const std::map<Key,Obj>& source, const Permutation& permvec)
{
  const typename std::map<Key,Obj>::const_iterator iterend(source.end());
  typename std::map<Key,Obj>::const_iterator iter(source.begin());
  while (iter!=iterend) {
    apply(dest[iter->first],iter->second,permvec);
    ++iter;
  }
}

  template<class T> bool HamiltonianStore<T>::verify(std::ostream& ostr, const CrystalStructure& cstruct, double tol) const
  {
    if (N_!=cstruct.ncells())
      throw Mismatch("HamiltonianStore::verify");
    if (N_<2)
      return true;
    
    const typename couplingmap_type::const_iterator iterend(couplings_end());
    typename couplingmap_type::const_iterator iter(couplings_begin());
    bool ok=true;
    while (iter!=iterend) {
      if (!verify_(ostr,iter->second,interaction_name(iter->first).c_str(),cstruct,tol))
	ok=false;
      ++iter;
    }

    const Permutation& perm(cstruct.permutation());
    if (!perm.empty()) {
      HamiltonianStore<T> permuted(*this,perm);
      if (!isequal(permuted,tol)) {
	ostr << "Permuted Hamiltonian (below) is not same as original\n" << permuted;
	return false;
      }
    }
      
    return ok;
  }
  
template<> const double HamiltonianStore<double>::null_=0.0;
template<> const space_T HamiltonianStore<space_T>::null_=space_T();

  template<class T> void HamiltonianStore<T>::split_couplings(Matrix<double>& A0s, Matrix<T>& A2s) const
  {
    bool haveA0=false;

    A2s.create(M_,total_,null_);

    typedef typename couplingmap_type::const_iterator iter_type;
    const iter_type iterend(couplings_end());
    iter_type iter(couplings_begin());
    
    reducer_<T,T> reduce;
    
    while (iter!=iterend) {
      const couplingstore<T>& store(iter->second);
      const Matrix<double>& storeA0(store.A0);
      const Matrix<T>& storeA2(store.A2);
      for (size_t i=M_;i--;) {
	for (size_t j=i+1;j<total_;j++) {
	  if (!!storeA0 && storeA0(i,j)) {
	    if (!haveA0) {
	      A0s.create(M_,total_,0.0);
	      haveA0=true;
	    }
	    A0s(i,j)+=store.A0(i,j);
	  }
	  if (!!storeA2) {
	    const T& coup(storeA2(i,j));
	    if (isnonnull(coup))
	      A2s(i,j)+=reduce(coup);
	  }
	}
      }
      ++iter;
    }
    if (!haveA0)
      A0s.clear();
  }

  template<class T> void HamiltonianStore<T>::docreate(int Mv,int Nv)
  {
    if ((Mv<=0) || (Nv<=0))
      throw InvalidParameter("HamiltonianStore");
    M_=Mv;
    N_=Nv;
    total_=M_*N_;
  }    
  
  template<class T> void HamiltonianStore<T>::print(std::ostream& ostr) const {
    ostr << "Spins: ";
    if (ncells()>1)
      ostr << ncells() << " x ";
    ostr << nspins_cell() << '\n';
    printqpole_(ostr,qstore);
    {    //dump shifts
      const typename shiftmap_type::const_iterator iterend(shifts_end());
      typename shiftmap_type::const_iterator iter(shifts_begin());
      while (iter!=iterend) {
	print_(ostr,interaction_name(iter->first).c_str(),iter->second);
	++iter;
      }
    }
    { //dump couplings
      const typename couplingmap_type::const_iterator iterend(couplings_end());
      typename couplingmap_type::const_iterator iter(couplings_begin());
      while (iter!=iterend) {
	print_(ostr,interaction_name(iter->first).c_str(),iter->second);
	++iter;
      }
    }
  }

  template<class T> void HamiltonianStore<T>::clear() {
    couplingmap.clear();
    shiftmap.clear();
    qstore.clear();
    qorder.clear();
  }

  template<class T> void HamiltonianStore<T>::print_(std::ostream& ostr, const char* name, const BaseList<double>& store) {
    if (!store.empty())
      ostr << name << ": " << store << '\n';
  }

  template<class T> void HamiltonianStore<T>::printqpole_(std::ostream& ostr, const BaseList<double>& store) const
  {
    if (store.empty()) return;
    for (size_t i=0;i<store.size();i++) {
      if (store(i))
	ostr << "Quadrupole coupling on " << i << ":  order: " << get_quadrupole_order(i) << ": " << (store(i)*1e-6) << " MHz\n";
    }
  }

ns_flag HamiltonianStructure::nonsecular_type(size_t& i, size_t& j) const
{
  if (nonsecular_(i))
    return nonsecular_(spin_to_spin(j)) ? NS_BOTH : NS_FIRST;
  if (nonsecular_(spin_to_spin(j))) {
    ::std::swap(i,j);
    return NS_FIRST;
  }
  return NS_NONE;
}

template<class T> void HamiltonianStore<T>::printqpole_(std::ostream& ostr, const BaseList<space_T>& store) const
{
  if (store.empty()) return;
  for (size_t i=0;i<store.size();i++) {
    if (!!store(i))
      ostr << "Quadrupole coupling on " << i << ":  order: " << get_quadrupole_order(i) << "\n" << store(i);
  }
}
  
 template<class T> void HamiltonianStore<T>::print_(std::ostream& ostr, const char* name, const BaseList<space_T>& store) {
   for (size_t i=0;i<store.size();i++) {
     if (!!store(i))
       ostr << name << " " << i << ":\n" << store(i) << '\n';
   }
 }
 
 template<class T> void HamiltonianStore<T>::print_(std::ostream& ostr, const char* name, const couplingstore<double>& store) { 
   if (!!store.A0)
     ostr << name << " (isotropic):\n" << store.A0 << '\n';
   if (!!store.A2)
     ostr << name << " (anisotropic):\n" << store.A2 << '\n';
 } 

  template<class T> void HamiltonianStore<T>::print_(std::ostream& ostr, const char* name, const couplingstore<space_T>& store) {
   if (!!store.A0)
     ostr << name << " (isotropic):\n" << store.A0 << '\n';
   const Matrix<space_T>& aniso(store.A2);
   for (size_t i=0;i<aniso.rows();i++) {
     for (size_t j=i+1;j<aniso.cols();j++) {
       if (!!aniso(i,j))
	 ostr << name << " " << i << "," << j << ":\n" << aniso(i,j) << '\n';
     }
   }
 }

  template<class T> void HamiltonianStore<T>::set_quadrupole(size_t ni, const T& v, size_t orderv)
  {
    set_unary(I_QUAD,ni,v);
    if (qorder.empty())
      qorder.create(M_,-1);
    int& which(qorder(ni));
    if (which<0)
      which=orderv;
    else {
      if (which!=orderv)
	throw Failed("set_quadrupole: can't change quadrupole order");
    }
  }
  
  template<class T> size_t HamiltonianStore<T>::get_quadrupole_order(size_t ni) const {
    if (ni>=qorder.size())
      throw BadIndex("get_quadrupole_order");
    const int order(qorder(ni));
    if (order<0)
      throw Failed("get_quadrupole_order: order unset");
    return order;
  }

template<class T> void HamiltonianStore<T>::set_coupling(interaction_t inttype, size_t ni, size_t nj, double iso, const T& v)
{
  if ((ni>=M_) || (nj>=total_))
    throw BadIndex("HamiltonianStore");
  if (ni==nj)
    throw InvalidParameter("HamiltonianStore: can't set coupling between same spin!");
  couplingstore<T>& storeref(couplingmap_cache_[inttype]);
  if (!storeref)
    create_null(storeref,M_,total_);
  if (iso) {
    if (!(storeref.A0))
      storeref.A0.create(M_,total_,0.0);
    (storeref.A0)(ni,nj)=iso;
    if (nj<M_)//within unit cell?
      (storeref.A0)(nj,ni)=iso;
  }
  if (isnonnull2(v)) {
    (storeref.A2)(ni,nj)=v;
    if (nj<M_)//within unit cell?
      (storeref.A2)(nj,ni)=v;
  }
}

template<> void HamiltonianStore<space_T>::set_coupling(interaction_t id, size_t ni, size_t nj, const space_T& v)
{
  const double iso=v.have_rank(0) ? real(v(0,0)) : 0.0;
  fudge_isotropic_<space_T> v2(v);
  set_coupling(id,ni,nj,iso,v2());
}

  template<class T> const T& HamiltonianStore<T>::get_anisotropic(interaction_t inttype, size_t ni, size_t nj) const
  {
    if ((ni>=M_) || (nj>=total_))
      throw BadIndex("HamiltonianStore");
    const couplingstore<T>* storep(couplingmap_cache_.find(inttype));
    if (storep) {
      const Matrix<T>& A2(storep->A2);
      if (!A2.empty())
	return A2(ni,nj);
    }
    return null_;
  }

  template<class T> double HamiltonianStore<T>::get_isotropic(interaction_t inttype, size_t ni, size_t nj) const
  {
    if ((ni>=M_) || (nj>=total_))
      throw BadIndex("HamiltonianStore");
    const couplingstore<T>* storep(couplingmap_cache_.find(inttype));
    if (storep) {
      const Matrix<double>& A0(storep->A0);
      if (!A0.empty())
	return A0(ni,nj);
    }
    return 0.0;
  }

  template<class T> void HamiltonianStore<T>::set_unary(interaction_t inttype,size_t ni, const T& v)
  {
    List<T>& storeref( (inttype==I_QUAD) ? qstore : shiftmap_cache_[inttype]);
    if (storeref.empty())
      create_null(storeref,M_);
    if (ni>=M_)
      throw BadIndex("HamiltonianStore");
    storeref(ni)=v;
  }

  template<class T> const T& HamiltonianStore<T>::get_unary(interaction_t inttype, size_t ni) const {
    if (ni>=M_)
      throw BadIndex("HamiltonianStore");
    if (inttype==I_QUAD)
      return qstore.empty() ? null_ : qstore(ni);
    const List<T>* storep(shiftmap_cache_.find(inttype));
    return storep ? (*storep)(ni) : null_;
  }

  template<class T> void HamiltonianStore<T>::invert_linear()
  {
    const typename shiftmap_type::iterator iterend(shifts_end());
    typename shiftmap_type::iterator iter(shifts_begin());
    while (iter!=iterend) {
      negate_ip(iter->second);
      ++iter;
    }
  }

  template<class T> HamiltonianStore<T>::HamiltonianStore(const HamiltonianStore<T>& a, const Permutation& permvec)
    : couplingmap_cache_(couplingmap), shiftmap_cache_(shiftmap)
  {
    docreate(a.nspins_cell(),a.ncells());
    if (!a.qstore.empty()) {
      permvec.apply(qstore,a.qstore);
      permvec.apply(qorder,a.qorder);
    }
    this->apply(shiftmap,a.shiftmap,permvec);
    this->apply(couplingmap,a.couplingmap,permvec);
  }

template<class T> HamiltonianStore<T>::HamiltonianStore(int Mv, int Nv)
  : couplingmap_cache_(couplingmap), shiftmap_cache_(shiftmap)
{ docreate(Mv,Nv); }

template<class T> HamiltonianStore<T>::HamiltonianStore(const HamiltonianStore<space_T>& source, const Euler& powder,bool compress)
  : qorder(source.get_quadrupole_order()),
    couplingmap_cache_(couplingmap), shiftmap_cache_(shiftmap)
{
  docreate(source.nspins_cell(), source.ncells());
  
  cmatrix RotM;
  if (powder!=Euler(0,0,0))
    RotM=D(2,powder);

  {
    typedef typename HamiltonianStore<space_T>::shiftmap_type::const_iterator iter_type;
    const iter_type iterend(source.shifts_end());
    iter_type iter(source.shifts_begin());
    List<T>* destp = compress ? &(shiftmap[I_CS]) : NULL;
    while (iter!=iterend) {
      if (!compress)
	destp=&(shiftmap[iter->first]);
      rotatesum(*destp,iter->second,RotM);
      ++iter;
    }
  }
  {
    typedef typename HamiltonianStore<space_T>::couplingmap_type::const_iterator iter_type;
    const iter_type iterend(source.couplings_end());
    iter_type iter(source.couplings_begin());
    while (iter!=iterend) {
      rotate(couplingmap[iter->first],iter->second,RotM);
      ++iter;
    }
  }
  rotatesum(qstore,source.get_quadrupole(),RotM);
}

  template<class T> void HamiltonianStore<T>::rotate(couplingstore<double>& dest, const couplingstore<space_T>& source, const cmatrix& RotM)
{
  if (!source.A0)
    dest.A0.clear();
  else
    dest.A0=source.A0;
  if (!source.A2) {
    dest.A2.clear();
    return;
  }

  const Matrix<space_T>& sourceA2(source.A2);
  Matrix<double>& destA2(dest.A2);
  const size_t M=sourceA2.rows();
  destA2.create(M,sourceA2.cols(),0.0);
  for (size_t r=sourceA2.rows();r--;) {
    for (size_t c=sourceA2.cols();(--c)>r;) {
      const space_T& A(sourceA2(r,c));
      if (!!A) {
	const double v=!RotM ? real(A(2,0)) : real(::libcmatrix::rotate(A,0,RotM));
	destA2(r,c)=v;
	if (c<M)
	  destA2(c,r)=v;
      }
    }
  }
}

  template<class T> void HamiltonianStore<T>::rotatesum(List<double>& dest, const List<space_T>& source, const cmatrix& RotM)
{
  size_t n=source.size();
  if (dest.empty()) 
    dest.create(n,0.0);
  for (;n--;) {
    const space_T& A(source(n));
    if (!!A) 
      dest(n)+=!RotM ? sumzero(A) : sumzero_rotate(A,RotM);
  } 
}

template<class T> void HamiltonianStore<T>::rotate(couplingstore<space_T>& dest, const couplingstore<space_T>& source, const cmatrix& RotM)
{
  if (!source.A0)
    dest.A0.clear();
  else
    dest.A0=source.A0;

  const Matrix<space_T>& sourceA2(source.A2); 
  Matrix<space_T>& destA2(dest.A2);
  if (!sourceA2) {
    destA2.clear();
    return;
  }
  if (!RotM) {
    destA2=sourceA2;
    return;
  }

  const size_t M=sourceA2.rows();
  destA2.create(M,sourceA2.cols());
  for (size_t r=sourceA2.rows();r--;) {
    for (size_t c=sourceA2.cols();(--c)>r;) {
      const space_T& A=sourceA2(r,c);
      if (!!A) {
	::libcmatrix::rotate(destA2(r,c),A,RotM);
	if (c<M)
	  destA2(c,r)=destA2(r,c);
      }
    }
  }
}

  template<class T> void HamiltonianStore<T>::rotatesum(List<space_T>& dest, const List<space_T>& source, const cmatrix& RotM)
{
  if (!RotM)
    dest+=source;
  else {
    size_t n=source.size();
    const bool neednew=dest.empty();
    if (neednew)
      dest.create(n);
    space_T tmp;
    for (;n--;) {
      const space_T& A(source(n));
      if (!!A) {
	if (neednew)
	  ::libcmatrix::rotate(dest(n),A,RotM);
	else {
	  ::libcmatrix::rotate(tmp,A,RotM);
	  dest(n)+=tmp;
	}
      }
    }
  }
}

  template class HamiltonianStore<space_T>;
  template class HamiltonianStore<double>;
  
void BlockedOperator::print(std::ostream& ostr, bool includeblocking) const {
  if (!*this)
    ostr << "<empty>\n";
  else {
    if (includeblocking) {
      ostr << blkspec_;
      ostr << "Is matched to identity matrix: " << (identityscale_.empty() ? "No\n" : "Yes\n");
    }

    block_pattern::iterator mziter(blkspec_);
    size_t r,c,mzeigSD;
    const size_t eigblks=eigblocks();    
    while (mziter.next(r,c,mzeigSD)) {
      for (size_t blk=0;blk<eigblks;blk++) {
	if (includeblocking) {
	  ostr << "mz eigenvalue: " << mzeigSD;
	  if (eigblks>1)
	    ostr << "  index eigenvalue: " << blk;
	  ostr << '\n';
	}
	ostr << (*this)(mzeigSD,blk);
      }
    }
  }
}

void BlockedOperator::dump() const 
{
  print(std::cout);
}

template<class BaseType> void BlockedSpinningHamiltonianBase<BaseType>::print(std::ostream& ostr) const
{
  const List<spinning_type>& aslist(*this);
  for (size_t k=0;k<aslist.size();k++)
    ostr << "Block " << k << ":\n" << aslist(k);
  ostr << '\n';
}

void BlockedDiagonalSpinningHamiltonianBase::operator()(ListList<double>& dest, double t) const
{
  if (dest.empty())
    opgen_base_.create(dest);
  else {
    if (dest.size()!=size())
      throw Mismatch("BlockedDiagonalSpinningHamiltonian");
  }
  for (size_t i=size();i--;)
    List<DiagonalSpinningHamiltonian>::operator()(i)(dest(i),t);
}

void BlockedDiagonalSpinningHamiltonianBase::print(std::ostream& ostr) const
{
  const List<DiagonalSpinningHamiltonian>& aslist(*this);
  for (size_t k=0;k<aslist.size();k++)
    ostr << "Block " << k << ":\n" << aslist(k);
  ostr << '\n';
}

// BlockedOperator& BlockedOperator::operator= (const BlockedOperator& a)
// {
//   static_cast<BaseStructure&>(*this)=a;
//   blkspec_=a.blkspec_;
//   if (&opgen_ != &(a.opgen_))
//     throw Failed("BlockedOperator: can't copy between operators based on different generators"); 
//   store_=a.store_;  
//   return *this;
// }

void BlockedOperator::copy_props(const BlockedOperator& a, bool dotranspose) 
{
  static_cast<BaseStructure&>(*this)=a;
  blkspec_=a.blkspec_;
  opgenp_= a.opgenp_;
  //  Hindexer_=a.Hindexer_;
  if (dotranspose)
    blkspec_.transpose();
}

  void transpose(BlockedOperator& a, const BlockedOperator& b)
  {
    if (b.blockstructure().isherm) {
      a=b;
      conj_ip(a.row());
    }
    else {
      transpose(a.row(),b.row());
      a.copy_props(b,true);
    }
  }

  void conj_transpose(BlockedOperator& a, const BlockedOperator& b)
  {
    if (b.blockstructure().isherm)
      a=b;
    else {
      conj_transpose(a.row(),b.row());
      a.copy_props(b,true);
    }
  }

template<class BaseType> void BlockedSpinningHamiltonianBase<BaseType>::reset()
{
  ensure();
  List<spinning_type>& H(*this);
  for (size_t n=H.size();n--;) {
    spinning_type& Hspin(H(n));
    Hspin.clear();
    if (!!Hbase_)
      Hspin+=Hbase_(n);
  }
}
 
template<class BaseType,class OpGen> void BlockedSpinningHamiltonianImpl<BaseType,OpGen>::interactions(const HamiltonianStore<space_T>& store)
{
  this->reset();
  List<typename BlockedSpinningHamiltonian<BaseType>::spinning_type>& H(*this);
  add_Hamiltonian(H,opgen_,store);
}
	   
template struct BlockedSpinningHamiltonianBase<double>;
template struct BlockedSpinningHamiltonianBase<complex>;

LCM_INSTANTIATE_TYPED_HAMILTONIANS(SpinOpGenerator,double)
LCM_INSTANTIATE_TYPED_HAMILTONIANS(SpinOpGenerator,complex)
LCM_INSTANTIATE_DIAGONAL_HAMILTONIANS(SpinOpGenerator)
LCM_INSTANTIATE_SPINOPGEN_CREATE(SpinOpGenerator)

  template void BlockedMatrix<double>::dump() const;	
  template void BlockedMatrix<complex>::dump() const;	

bool BlockedFilter::isnonzero() const
{
  const BaseList<bool> asrow(store_.row());
  for (size_t i=asrow.size();i--;) {
    if (asrow(i))
      return true;
  }
  return false;
}

void BlockedFilter::coherence_filter_ip(const ListList<double>& Fz, const spinorder_spec& spec, const state_t cohers)
{
  size_t r,c,mzeigSD;
  const bool isherm=blkspec_.ishermitian();
  block_pattern::iterator mziter(blkspec_);
    
  while (mziter.next(r,c,mzeigSD)) {
    for (size_t blk=eigblocks();blk--;) {
      Matrix<bool> d((*this)(mzeigSD,blk));
      if (!d)
	continue;
      const size_t rind(opgen_.index(r,blk));
      const size_t cind(opgen_.index(c,blk));
      const BaseList<double> rFz(Fz(rind));
      const BaseList<double> cFz(Fz(cind));
      if (rFz.size()!=d.rows())
	throw Mismatch("coherence_filter_ip (rows)",rFz.size(),d.rows());
      if (cFz.size()!=d.cols())
	throw Mismatch("coherence_filter_ip (cols)",cFz.size(),d.cols());
      for (size_t i=d.rows();i--;) {
	for (size_t j=d.cols();j--;) {
	  const int coher= rFz(i)-cFz(j);
	  const bool isallowed=spec.iscoherence(coher,cohers);
	  if (isherm) {
	    const bool isotherallowed=spec.iscoherence(-coher,cohers);
	    if (isallowed == isotherallowed)
	      throw Failed("BlockedFilter: hermitian block structure is incompatible with non-hermitian filter matrix");
	  }
	  d(i,j) &= isallowed;
	}
      }
    }
  }
}

std::ostream& operator<< (std::ostream& ostr, const BlockedFilter& a)
{
  ostr << a.blockstructure();
  return a.row().dump_content(ostr);
}

BlockedFilter::BlockedFilter(const SpinOpGenerator& opgen, const block_pattern& blkpattern, const spinorder_spec& spec)
  : blkspec_(blkpattern), opgen_(opgen)
{
  opgen.create(store_,blkpattern);
  static_cast< BaseStructure& >(*this)= BaseStructure(blkspec_.blocks,opgen.eigblocks());

  if (spec.isspinorder()) {
    const basespin_system& sys(opgen.spinsystem());
    if (!(sys.isspinhalfonly()))
      throw Failed("spin order filter only implemented for spin-1/2 nuclei");
    if (eigblocks()!=1)
      throw Failed("spin order filter not implemented for periodic spin systems");

    size_t r,c,m;

  //! NB this may not be correct for hermitian block structures!

    block_pattern::iterator iter(blkpattern);
    while (iter.next(r,c,m)) {
      if (!store_(m))
	continue;
      
      const BaseList<size_t> rstates(opgen.blockindices(r));    
      const BaseList<size_t> cstates(opgen.blockindices(c)); 
      Matrix<bool> curm(store_(m),mxflag::nondynamic);
      
#ifndef NDEBUG
      std::cout << "Constructing spin- / coherence- order for bra states: " << rstates << " and ket states " << cstates << '\n';
#endif
      const size_t nr=rstates.size();
      const size_t nc=cstates.size();
      if ((nr!=curm.rows()) || (nc!=curm.cols()))
	throw Mismatch("BlockedFilter from spinorder_spec",nr,nc,curm.rows(),curm.cols());
      for (size_t i=nr;i--;) {
	const state_t statei(rstates(i));
	for (size_t j=nc;j--;)
	  curm(i,j)=spec(statei,cstates(j));
      }
    }
  }
  else {
    store_=true;
    ListList<double> tmpFz;
    opgen.create(tmpFz);
    
    for (size_t i=spec.size();i--;) {
      const spinorder_spec::subspec_t& curspec(spec.item(i));
      tmpFz=0.0;
      size_t curi=0;
      state_t which=curspec.first;
      while (which) {
	if (which & 1)
	  opgen.mla_Iz(tmpFz,1.0,spec.bit_to_spin(curi));
	which >>= 1;
	curi++;
      }
#ifndef NDEBUG
      std::cout << "Fz for compressed spin indices " << curspec.first << ": " << tmpFz << '\n';
#endif
      coherence_filter_ip(tmpFz,spec,curspec.second);
    }
  }
}
//   const size_t nidblks=nucids.size();
//   for (size_t i=nidblks;i--;) {
//     if (opgen.isblocked(nucids(i)))
//       throw Failed("Can't initialise BlockedFilter with blocked coherence");
//   }
//   if (nidblks!=cohers.size())
//     throw Mismatch("BlockedFilter: number of nuclei types vs. number of coherence lists",nidblks,cohers.size());

//   const block_pattern blkpattern(opgen); //always diagonal
//   opgen.create(store_,blkpattern);
//   store_=true;

//   for (size_t n=nidblks;n--;) {
//     const ListList<double>& nucFz(opgen.diag_Fz(nucids(n)));
//     const BaseList<int> cohersn(cohers(n));
//     const int maxcoher=coherencematrix_range(nucFz.row(),cohersn);
//     for (size_t m=nucFz.size();m--;)
//       coherencematrix_ip(store_(m),nucFz(m),cohersn,maxcoher);
//   }

//! instantiate BlockedFilter
//template BlockedFilter::BlockedFilter(const SpinOpGenerator&, const BaseList<int>&);
//template BlockedFilter::BlockedFilter(const SpinOpGenerator&, const block_pattern&, const spinorder_spec&);

template<class T> void BaseMetaPropagator::step_propagator(cmatrix& U, const Matrix<T>& H, size_t which, double dt) const
{
  if (flags_ & usechebyshev)
    chebyshev_propagator(U,H,dt,partition(which));
  else
    ::libcmatrix::propagator(U,H,dt);//,partition(which));
}

template<class T> void BaseMetaPropagator::step_propagator(BlockedMatrix<complex>& U, const BlockedMatrix<T>& H, int, double dt) const
{
  if (flags_ & usechebyshev) {
    U.duplicate_structure(H);
    for (size_t which=U.size();which--;) 
      chebyshev_propagator(U(which),H(which),dt,partition(which));
  }
  else
    ::libcmatrix::propagator(U,H,dt);
}
  
template void BaseMetaPropagator::step_propagator(cmatrix&, const Matrix<double>&, size_t, double) const;
template void BaseMetaPropagator::step_propagator(cmatrix&, const Matrix<complex>&, size_t, double) const;

template void BaseMetaPropagator::step_propagator(BlockedMatrix<complex>&, const BlockedMatrix<double>&, int, double) const;
template void BaseMetaPropagator::step_propagator(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, int, double) const;

  void BaseMetaPropagator::operator()(BlockedMatrix<complex>& U, double t1, double t2) const
  {
    ensure(U);
    for (size_t blk=blkstr_.size();blk--;) {
      size_t mzblk,eigblk;
      indexer_.reverse(mzblk,eigblk,blk);
      (*this)(U(blk),t1,t2,mzblk,eigblk);
    }
  }

  void BaseMetaPropagator::propagators(BaseList<cmatrix>& Us, double t1, double t2, size_t mzeig, size_t blk) const
  {
    const size_t n=Us.size();
    if (n==0)
      throw InvalidParameter("propagators: empty destination list");
    const double dt=(t2-t1)/n;
    if (dt<0.0)
      throw InvalidParameter("BaseMetaPropagator: end time < start time");
    if ((dt==0.0) && (n>1))
      throw InvalidParameter("BaseMetaPropagator: multiple observations incompatible with zero dwell time");
    if (verbose_) {
      std::cout << "MetaPropagator: creating " << n << " propagator(s) for block (" << mzeig << ',' << blk << ") over the interval ";
      prettyprint_time(t1) << " to ";
      prettyprint_time(t2) << '\n';
    }
    double t=t1;
    for (size_t j=0;j<n;j++) {
      const double nextt=(j==n-1) ? t2 : t1+(j+1)*dt;
      (*this)(j ? tmp : Us.front(),t,nextt,mzeig,blk); //intdt,Bool2Type<Ham_traits<HType>::isconstant>());
      if (j) {
		if (!tmp)
			Us(j)=Us(j-1);
		else {
			if (!Us(j-1))
				Us(j)=tmp;
			else
				multiply(Us(j),tmp,Us(j-1));
		}
      }
      t=nextt;
    }
	// difficult to work with undefined propagators. Expand to identity at this point.
	for (size_t j=n;j--;) {
		if (!Us(j))
			Us(j).identity(opgen_base_.size(mzeig, blk));
	}
  }
  
  void BaseMetaPropagator::ensure(BlockedMatrix<complex>& U) const
  {
    if (U.empty())
      U.create(blkstr_);
    else {
      if (U.size()!=blkstr_.size())
	throw Mismatch("MetaPropagator()");
    }
  }
  
  void BaseMetaPropagator::ensure(ListList<complex>& U) const
  {
    if (U.empty())
      U.create(blkstr_);
    else {
      if (U.size()!=blkstr_.size())
	throw Mismatch("MetaPropagator()");
    }
  }

void BaseMetaPropagator::ensure(BlockedMatrix<complex>& U, size_t mzblk, size_t eigblk) const
{
  const size_t n(opgen_base_.size(mzblk,eigblk));
  if ((U.size()!=1) || (U.front().rows()!=n))
    U.create(ExplicitList<1,size_t>(n));
}

void BaseMetaPropagator::try_compress(const BlockedSpinningHamiltonian<complex>& Hspin)
{
  const SpinOpGeneratorBase& opgen(generator());
  //  if (opgen.totalblocks()>1)
  //  throw Failed("MetaPropagator: blocked second-order Hamiltonians not supported");
  Hdiagps_.create(opgen.totalblocks());
  for (size_t index=0;index<opgen.totalblocks();index++) {
    size_t mzblk,eigblk;
    opgen.structure().reverse(mzblk,eigblk,index);
    DiagonalSpinningHamiltonian* Hdiagp=new DiagonalSpinningHamiltonian(Hspin(mzblk,eigblk),opgen.Hzeeman(mzblk,eigblk),opgen.Hzeeman_structure(mzblk,eigblk).row());
    Hdiagps_(index).reset(Hdiagp);
    if (opgen.verbose()>1)
      std::cout << "Diagonal Hamiltonian, block " << index << "\n" << (*Hdiagp) << '\n';
  }
}

#include "lcm_CommonNMR.cc"

BlockedMatrix<complex> gammareduce(const BaseList< BlockedMatrix<complex> >& source, size_t n)
{
  BlockedMatrix<complex> U;
  gammareduce_(U,source,source.size(),n);
  return U;
}

void gammareduce(BaseList< BlockedMatrix<complex> > dest, const BaseList< BlockedMatrix<complex> >& source, size_t n)
{
  gammareduce_(dest,source,n);
}

size_t SpinOpGenerator::size(size_t mzblk, size_t eigblk) const
{
  if (eigblk)
    throw InvalidParameter("SpinOpGenerator has no eigenvalue structure");
  if (isblocked())
    return blockindices(mzblk).size();
  if (mzblk)
    throw BadIndex("SpinOpGenerator::size");
  return SpinOpGeneratorBase::size();
}

operator_spec operator& (const operator_spec& a, const operator_spec& b)
{
  operator_spec res(a);
  res.op=operator_intersection(a.op,b.op);
  if ((a.nuc!=b.nuc) || (a.number!=b.number) || !(res.op) || (a.col!=b.col))
    throw Failed("operator_spec&: operators do not intersect");
  return res;
}

std::ostream& operator<< (std::ostream& ostr, const operator_spec& a)
{
  if (a.nuc==NULL_NUCLEUS) {
    ostr << "I" << (a.number+cmatrix_ostream_controller(ostr).indexbase);

    if (operator_spec::isvalidST(a.op))
      ostr << ':';
  }
  else
    ostr << nuctolabel(a.nuc) << ':';
  
  if (operator_spec::isvalidST(a.op))
    return ostr << (1+size_t(a.op)) << ',' << (1+size_t(a.col));
  else
    return ostr << a.op;

  return ostr;
}
  
void operator_spec::initST(size_t r, size_t c)
{
  if (!isvalidST(r) || !isvalidST(c))
    throw InvalidParameter("operator_spec: invalid single transition specifier");      
  op=char(r);
  col=char(c);
}

void rotatez_ip(BlockedMatrix<complex>& U,const ListList<double>& Fz,double angle)
{
  if (U.size()!=Fz.size())
    throw Mismatch("rotatez_ip");
  for (size_t i=Fz.size();i--;)
    rotatez_ip(U(i),Fz(i),angle);
}

void rotatezfacs(ListList<complex>& zrot, const ListList<double>& Fz, double angle)
{
  zrot.duplicate_structure(Fz);
  for (size_t i=Fz.size();i--;)
    rotatezfacs(zrot(i),Fz(i),angle);
}

BlockedMatrix<complex> rotatez(const BlockedMatrix<complex>& U, const ListList<double>& Fz, double angle)
{
  BlockedMatrix<complex> dU(U);
  if (angle)
    rotatez_ip(dU,Fz,angle);
  return dU;
}


void BaseMetaPropagator::set(BlockedMatrix<complex>& Uacc, const BlockedMatrix<complex>& U, int which)
{
  if (which<0)
    Uacc=U;
  else
    Uacc.set(U(which));
}

bool MetaFlags::extract(int& flags, int flag, const char* warnmess) 
{
  const bool isset=((flags & flag)!=0);
  flags&=~flag;
  if ((warnmess!=NULL) && (flags!=0))
    std::cerr << warnmess << ": " << flags << '\n';
  return isset;
}

} //namespace libcmatrix

// BaseInhomogeneousH::BaseInhomogeneousH(const SpinOpGenerator& opgen, const HamiltonianStore<space_T>& store, double rotor_speed, double rotor_phase, int nobs, const RotorInfo& rinfo)
// {
//   if (!spinsystem_glue_<SpinOpGenerator>::arematching(opgen,store))
//     throw Mismatch("BaseInhomogeneousH");
//   if (opgen.mzblocks()!=1)
//     throw Failed("BaseInhomogeneousH: not compatible with mz blocking");
//   if (opgen.quadrupole_order()!=1)
//     throw Failed("BaseInhomogeneousH: not compatible with non-secular Hamiltonian");
  
//   BaseInhomogeneousH dest;

//   const int ints=find(opgen,store);
//   if (ints<=0)
//     throw Failed("BaseInhomogeneousH: Hamiltonian is not (spin) inhomogeneous");

//   Hs.create(ints,opgen.size(0,0));
//   Hs=0.0;
//   dphases.create(ints,DynamicPhase(rotor_speed,rotor_phase,nobs,rinfo));

//   //add interactions
//   find(opgen,store);
// }

//   //count number of inhomogeneous interactions 
// int BaseInhomogeneousH::find(const SpinOpGenerator& opgen, const HamiltonianStore<space_T>& store) {
//   const size_t N=store.ncells();
//   int count=0;
//   if (!find(count,I_DIPOLE,opgen,store.get_dcouplings(),N))
//     return -1;
//   if (!find(count,I_J,opgen,store.get_jcouplings(),N))
//     return -1;
//   if (!find(count,I_SHIFT,opgen,store.get_shifts(),N))
//     return -1;
//   if (!find(count,I_QUAD,opgen,store.get_quadrupole(),N))
//     return -1;
//   return count;
// }

// bool BaseInhomogeneousH::find(int& count, inter_t id, const SpinOpGenerator& opgen, const BaseList<space_T>& coups, size_t cells)
// {
//   if (coups.empty())
//     return true;

//   const size_t M=coups.size();
//   Indexer<2> cell_to_spin(cells,M);

//   const basespin_system& sys=spinsystem_glue_<SpinOpGenerator>::spinsystem(opgen);
  
//   const bool doadd=!!Hs;
//   for (size_t j=coups.size();j--;) {
//     if (!coups(j))
//       continue;
    
//     if (doadd) {
//       dphases(count).tensor(coups(j));
//       BaseList<double> Hrow(Hs.row(count));
      
//       for (size_t cell=cells;cell--;) {
// 	const size_t sk=cell_to_spin(cell,j);
	
// 	switch (id) {
// 	case I_SHIFT:
// 	  Hrow+=diag_spin_CS(sys,sk);
// 	  break;
// 	case I_QUAD:
// 	  Hrow+=diag_spin_quadrupolar(sys,sk);
// 	  break;
// 	default:
// 	  throw Failed("InhomogeneousHamiltonian: unknown unary interaction");
// 	}
//       }
//     }
//     count++;
//   }
//   return true;
// }

// bool BaseInhomogeneousH::find(int& count, inter_t id, const SpinOpGenerator& opgen, const Matrix<space_T>& coups, size_t cells)
// {
//   if (coups.empty())
//     return true;
  
//   const bool doadd=!!Hs;
//   const size_t M=coups.rows();
//   if (M*cells!=coups.cols())
//     throw Mismatch("find_interactions (binary)");

//   const basespin_system& sys=opgen.spinsystem();
  
//   for (size_t i=0;i<M;i++) {
//     for (size_t j=0;j<M;j++) {
//       for (size_t offset=0;offset<cells;offset++) {
// 	if ((i>=j) && (offset==0))
// 	  continue; //quick check for i>=j within cell
// 	const size_t sk=opgen.cell_to_spin(offset,j);
// 	const space_T& coup=coups(i,sk);
// 	if (!coup)
// 	  continue;
// 	if (opgen(i)==opgen(sk))
// 	  return false; //found homonuclear interaction
	
// 	if (doadd) {
// 	  dphases(count).tensor(coup);
// 	  BaseList<double> Hrow(Hs.row(count));
	  
// 	  for (size_t cellbase=0;cellbase<cells;cellbase++) {
// 	    const size_t finali=opgen.cell_to_spin(cellbase,i);
// 	    const size_t finalj=opgen.cell_to_spin((cellbase+offset) % cells,j);
// 	    if (finali>=finalj)
// 	      continue;
// 	    switch (id) {
// 	    case I_DIPOLE:
// 	      Hrow+=diag_spin_dipolar(sys,finali,finalj);
// 	      break;
// 	    case I_J:
// 	      Hrow+=diag_spin_weakJ(sys,finali,finalj);
// 	      break;
// 	    default:
// 	      throw InternalError("find_interactions (binary)");
// 	    }
// 	  }
// 	}
// 	count++;
//       }
//     }
//   }
//   return true;
// }

//   InhomogeneousHamiltonian::InhomogeneousHamiltonian(const SpinOpGenerator& opgen, const HamiltonianStore<space_T>& store, double rotor_speed, double rotor_phase, const RotorInfo& rinfo)
//     : BaseStructure(opgen),
//       homoinfo(opgen,store,rotor_speed,rotor_phase,0,rinfo) {
    
//     const rmatrix& source=homoinfo.Hs;
//     const size_t rows=source.rows();
    
//     const BlockingNuclei& bstr(opgen.blockingnuc());
//     const ListList<size_t>& inds(bstr.mzblocking());
//     if (inds.empty()) {
//       create(ExplicitList<1,size_t>(rows),ExplicitList<1,size_t>(bstr.size()));
//       this->front()=source;
//     }
//     else {
//       const size_t blocks=bstr.getblocks();
//       const DynamicList<size_t> rstr(blocks,rows,mxflag::normal);
//       DynamicList<size_t> cstr(blocks,mxflag::normal);
//       size_t k;
//       for (k=0;k<blocks;k++)
// 	cstr(k)=inds(k).size();
      
//       create(rstr,cstr);
      
//       for (k=blocks;k--;)
// 	(*this)(k)=source(range(),inds(k));
//     }
//   }
  
// std::ostream& operator<< (std::ostream& ostr, const BaseInhomogeneousH& a)
// {
//   if (!a)
//     ostr << "<undefined>\n";
//   else {
//     ostr << "Hamiltonians:\n" << a.Hs;
//     ostr << "Dynamic phases:\n";
//     for (size_t j=a.dphases.size();j--;)
//       ostr << a.dphases(j);
//   }
//   return ostr;
// }

// block_pattern::iterator::iterator(const block_pattern& blkspec)
// {
//   blocks=blkspec.blocks;
//   if (blocks==0)
//     throw InvalidParameter("block_pattern::iterator()");
  
//   if (blkspec.isdiagonal()) {
//     mzstep=1;
//     offset=0;
//     superbase=0;
//     majorincr=blkspec.blocks;
//     majorcount=1;
//     smallblocks=blkspec.blocks;
//   }
//   else {
//     mzstep=blkspec.mzstep;
//     offset=blkspec.coher*mzstep;
//     majorincr=blkspec.mzlevels*mzstep;
//     majorcount=blkspec.totlevels / majorincr;
// #ifndef NDEBUG
//     if ((blkspec.totlevels % majorincr) || (blkspec.blocks % majorcount))
//       throw InternalError("block_pattern::iterator");
// #endif
//     smallblocks=blkspec.blocks/majorcount;
//   }
//   superbase=0;
//   curadd=0;
//   majorcount--;
//   finished=false;
//   maxind=blkspec.maxind;
//   reset_block();
// }

// void block_pattern::iterator::reset_block()
// {
//   if (curadd==mzstep) {
//     superbase+=majorincr;
//     if (majorcount==0) {
//       finished=true;
//       return;
//     }
//     majorcount--;
//     curadd=0;
//   }

//   const size_t base=superbase+curadd;
//   if (offset>0) {
//     crind=base;
//     ccind=base+offset;
//   }
//   else {
//     crind=base-offset;
//     ccind=base;
//   }
//   smallblockcount=smallblocks-1;
// }

// bool block_pattern::iterator::next(size_t& rind, size_t& cind)
// {
//   if (finished)
//     return false;

//   rind=crind;
//   cind=ccind;
//   //increment
//   if (smallblockcount==0) {
//     curadd++;
//     reset_block();
//   }
//   else {
//     smallblockcount--;
//     crind+=mzstep;
//     ccind+=mzstep;
//   }
//   return true;
// }

