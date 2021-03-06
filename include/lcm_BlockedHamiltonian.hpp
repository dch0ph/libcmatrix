//internal file, all BlockedHamiltonians

template<class BaseType> struct BlockedStaticHamiltonianBase : public BlockedMatrix<BaseType>  {
  BlockedStaticHamiltonianBase(const SpinOpGeneratorBase& opgen_basev,const BlockedMatrix<BaseType>& Hbasev) : opgen_base_(opgen_basev), Hbase_(Hbasev) {}
  virtual ~BlockedStaticHamiltonianBase() {}

  void print(std::ostream& ostr) const { ostr << static_cast< const BlockedMatrix<BaseType>& >(*this); }

  const Matrix<BaseType>& operator()(size_t mzeig, size_t blk) const { return BlockedMatrix<BaseType>::operator()(opgen_base_.index(mzeig,blk)); }

  virtual BlockedStaticHamiltonianBase<BaseType>* clone() const =0;

  virtual void interactions(const HamiltonianStore<double>&) =0;
  virtual void interactions(const HamiltonianStore<space_T>&) =0;  

  void reset() {
    BlockedMatrix<BaseType>& H(*this);
    if (!Hbase_) {
      if (!!H)
	H=BaseType(0.0);
    }
    else
      H=Hbase_;
  }

  const SpinOpGeneratorBase& generator() const { return opgen_base_; }

  const SpinOpGeneratorBase& opgen_base_;
  const BlockedMatrix<BaseType> Hbase_;
};

template<class BaseType,class OpGen> class BlockedStaticHamiltonianImpl : public BlockedStaticHamiltonianBase<BaseType>
{
public:
  BlockedStaticHamiltonianImpl(const OpGen& opgenv, const BlockedMatrix<BaseType>& Hbasev)
    : BlockedStaticHamiltonianBase<BaseType>(opgenv,Hbasev), opgen_(opgenv) {}

  BlockedStaticHamiltonianBase<BaseType>* clone() const { return new BlockedStaticHamiltonianImpl(*this); }

  void interactions(const HamiltonianStore<double>&);
  void interactions(const HamiltonianStore<space_T>&);

private:
  const OpGen& opgen_;
};

template<class BaseType> class BlockedStaticHamiltonian : public BaseStructure
{
public:

  template<class OpGen> BlockedStaticHamiltonian(const OpGen&, const HamiltonianStore<double>&, const BlockedMatrix<BaseType>& =BlockedMatrix<BaseType>());
    
  template<class OpGen> BlockedStaticHamiltonian(const OpGen&, const BlockedMatrix<BaseType>& =BlockedMatrix<BaseType>());


  const SpinOpGeneratorBase& generator() const { return pImpl_->generator(); }

  typedef BlockedMatrix<BaseType> result_type;
  typedef BlockedMatrix<BaseType> output_type;
  typedef Matrix<BaseType> suboutput_type;
  typedef Matrix<BaseType> subresult_type;

  bool operator!() const { return pImpl_->operator!(); }

  const BlockedMatrix<BaseType>& operator()() const { return *pImpl_; }
  const Matrix<BaseType>& operator()(size_t mzblk, size_t eigblk) const { return pImpl_->operator()(mzblk,eigblk); }
    
  void interactions(const HamiltonianStore<double>& Hstore) { pImpl_->interactions(Hstore); }
  void interactions(const HamiltonianStore<space_T>& Hstore) { pImpl_->interactions(Hstore); }
  void print(std::ostream& ostr) const {
    ostr << static_cast<const BaseStructure&>(*this);
    pImpl_->print(ostr); }
  
private:
  smartptr<BlockedStaticHamiltonianBase<BaseType> > pImpl_;
};

  
template<class BaseType> std::ostream& operator<< (std::ostream& ostr, const BlockedStaticHamiltonian<BaseType>& a)
{
  a.print(ostr);
  return ostr;
}

  template<class T> void propagator(BlockedMatrix<complex>& U, const BlockedStaticHamiltonian<T>& H, double dt)
  {
    propagator(U,H(),dt);
  }

  template<class T> void propagator_ns(BlockedMatrix<complex>& U, const BlockedStaticHamiltonian<T>& H, double dt,const ListList<double>& Hzeeman)
  {
    propagator_ns(U,H(),dt,Hzeeman);
  }

  template<class T> void propagator_second(BlockedMatrix<complex>& U, const BlockedStaticHamiltonian<T>& H, double dt,const BaseList< ListList<size_t> >& Zstruct, const ListList<double>& Hzeeman)
  {
    propagator_second(U,H(),dt,Zstruct,Hzeeman);
  }

  template<class BaseType> struct Ham_traits<BlockedStaticHamiltonian<BaseType> > {
    typedef double coupling_type;
    static const bool isconstant=true;
    static const bool allowRF=true;
    static const bool isdiagonal=false;
    static const bool ishomogeneous=false;
  };

  template<class BaseType> struct Propagator_traits<BlockedStaticHamiltonian<BaseType> > {
    static const bool allowsetH=true;
  };

  template<class OutType, class OpGen, class StoreType> void add_Hamiltonian(OutType&, const OpGen&, const HamiltonianStore<StoreType>&);
  

struct BlockedDiagonalStaticHamiltonianBase : public ListList<double>  {
  explicit BlockedDiagonalStaticHamiltonianBase(const SpinOpGeneratorBase& opgen_basev)
    : opgen_base_(opgen_basev) {}

  virtual ~BlockedDiagonalStaticHamiltonianBase() {}

  size_t usage() const { return row().size(); }
  BaseList<double> operator()(size_t mzeig, size_t blk) const { return ListList<double>::operator()(opgen_base_.index(mzeig,blk)); }
  
  virtual BlockedDiagonalStaticHamiltonianBase* clone() const =0;

  virtual void interactions(const HamiltonianStore<double>&) =0;
  virtual void interactions(const HamiltonianStore<space_T>&) =0;

//   template<class BaseType> void interactions(const BlockedStaticHamiltonianBase<BaseType>&) {
//     if (opgen_base_.quadrupole_order()!=2)
//       throw InternalError("BlockedDiagonalStaticHamiltonian");
//     throw Failed("Not imp");
//   }

  const SpinOpGeneratorBase& generator() const { return opgen_base_; }
  void print(std::ostream& ostr) const { ostr << static_cast<const ListList<double>& >(*this) << '\n'; }
  void reset() { if (!empty()) { ListList<double>::operator=(0.0); } }
    
  const SpinOpGeneratorBase& opgen_base_;
};

template<class OpGen> class BlockedDiagonalStaticHamiltonianImpl : public BlockedDiagonalStaticHamiltonianBase
{
public:
  explicit BlockedDiagonalStaticHamiltonianImpl(const OpGen& opgenv)
    : BlockedDiagonalStaticHamiltonianBase(opgenv), opgen_(opgenv) {}
    
  BlockedDiagonalStaticHamiltonianBase* clone() const { return new BlockedDiagonalStaticHamiltonianImpl(*this); }

  void interactions(const HamiltonianStore<double>&);
  void interactions(const HamiltonianStore<space_T>&);

private:
  const OpGen& opgen_;
};
  
class BlockedDiagonalStaticHamiltonian : public BaseStructure {
public:
  
  template<class OpGen> BlockedDiagonalStaticHamiltonian(const OpGen&, const HamiltonianStore<double>&);
  template<class OpGen> explicit BlockedDiagonalStaticHamiltonian(const OpGen&);

  typedef ListList<double> result_type;
  typedef BaseList<double> subresult_type;
  typedef ListList<double> output_type;
  typedef BaseList<double> suboutput_type;

  const BaseList<double> operator()(size_t mzeig,size_t eig) const { return pImpl_->operator()(mzeig,eig); }
  const ListList<double>& operator()() const { return *pImpl_; }

  bool operator!() const { return pImpl_->empty(); }

  const SpinOpGeneratorBase& generator() const { return pImpl_->generator(); }
  void interactions(const HamiltonianStore<double>& Hstore) { pImpl_->interactions(Hstore); }
  void interactions(const HamiltonianStore<space_T>& Hstore) { pImpl_->interactions(Hstore); }
  void print(std::ostream& ostr) const {
    ostr << static_cast<const BaseStructure&>(*this);
    pImpl_->print(ostr); }
  size_t usage() const { return pImpl_->usage(); }
  
private:
  smartptr<BlockedDiagonalStaticHamiltonianBase> pImpl_;
};

  inline std::ostream& operator<< (std::ostream& ostr, const BlockedDiagonalStaticHamiltonian& a)
    {
      a.print(ostr);
      return ostr;
    }

template<> struct Ham_traits<BlockedDiagonalStaticHamiltonian> { 
  typedef double coupling_type; 
  static const bool allowRF=false;
  static const bool isconstant=true;
  static const bool isdiagonal=true;
  static const bool ishomogeneous=false;
}; 

template<> struct Propagator_traits<BlockedDiagonalStaticHamiltonian> {
  static const bool allowsetH=true;
};


template<typename BaseType> struct lcm_promote_spin_ {
  typedef SpinningHamiltonian spinning_type;
};
template<> struct lcm_promote_spin_<double> {
  typedef RealSpinningHamiltonian spinning_type;
};

template<class StoreType> struct RotorHolder_ : public List<StoreType> {
  RotorHolder_(double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo)
    : rinfo_(rinfo), samplerp_(NULL) ,rotor_speed_(rotor_speedv), rotor_phase_(rotor_phasev) {}
  
  RotorHolder_(double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const IntervalSamplerBase& sampler)
    : rinfo_(rinfo), samplerp_(&sampler) ,rotor_speed_(rotor_speedv), rotor_phase_(rotor_phasev) {}
  
  typedef StoreType spinning_type;

  bool operator!() const { return this->empty(); }
  size_t usage() const { return this->empty() ? 0 : this->size()*(this->front().usage()); }

  void rotor_phase(double phase_) {
    rotor_phase_=phase_;
    for (size_t k=this->size();k--;)
      (*this)(k).rotor_phase(phase_);
  }

  void ensure(int n) {
    if (this->empty()) {
      if (samplerp_) {
	const IntervalSamplerBase::synchronise_t syncval=samplerp_->get_synchronisation();
	List<StoreType>::create(n,StoreType(rotor_speed_,rotor_phase_,rinfo_,*samplerp_));
	for (size_t k=this->size();k--;)
	  ((*this)(k).samplerp())->set_synchronisation(syncval); //!< resets random number generators to consistent seed value
      }
      else
	List<StoreType>::create(n,StoreType(rotor_speed_,rotor_phase_,rinfo_));
    }
  }
  void clear() {
    for (size_t k=this->size();k--;)
      (*this)(k).clear();
  }

  double rotor_phase() const { return rotor_phase_; }
  double phase(double t) const {
    if (this->empty())
      throw Undefined("RotorHolder: not initialised");
    return this->front().phase(t);
  }

  void rotor_speed(double speed_) { 
    rotor_speed_=speed_;
    for (size_t k=this->size();k--;)
      (*this)(k).rotor_speed(speed_);
  }

  double rotor_period() const { return fabs(1.0/rotor_speed_); }

  const IntervalSamplerBase* samplerp() const { return samplerp_; }

  const RotorInfo& rinfo_;
  const IntervalSamplerBase* samplerp_;
  double rotor_speed_,rotor_phase_;
};

template<class BaseType> struct BlockedSpinningHamiltonianBase : public RotorHolder_< typename lcm_promote_spin_<BaseType>::spinning_type > 
{
  typedef typename lcm_promote_spin_<BaseType>::spinning_type spinning_type;
  BlockedSpinningHamiltonianBase(const SpinOpGeneratorBase& opgen_basev, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const BlockedMatrix<BaseType>& Hbasev) : RotorHolder_< typename lcm_promote_spin_<BaseType>::spinning_type >(rotor_speedv,rotor_phasev,rinfo), opgen_base_(opgen_basev), Hbase_(Hbasev) {
//     if (opgen_basev.isclassicsecondorder())
//       throw Failed("BlockedSpinningHamiltonian: not compatible with classic 2nd order treatments");
  }
  BlockedSpinningHamiltonianBase(const SpinOpGeneratorBase& opgen_basev, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const IntervalSamplerBase& samplerv, const BlockedMatrix<BaseType>& Hbasev) : RotorHolder_< typename lcm_promote_spin_<BaseType>::spinning_type >(rotor_speedv,rotor_phasev,rinfo,samplerv), opgen_base_(opgen_basev), Hbase_(Hbasev) {
//     if (opgen_basev.isclassicsecondorder())
//       throw Failed("BlockedSpinningHamiltonian: not compatible with classic 2nd order treatments");
  }
  virtual ~BlockedSpinningHamiltonianBase() {}

  const SpinOpGeneratorBase& generator() const { return opgen_base_; }

  void operator()(BlockedMatrix<BaseType>&, double) const;
  void print(std::ostream&) const;

  const spinning_type& operator()(size_t mzeig, size_t blk) const
  { return RotorHolder_<spinning_type>::operator()(opgen_base_.index(mzeig,blk)); }

  void reset();
  void ensure() {
    RotorHolder_<spinning_type>::ensure(opgen_base_.actual_mzblocks()*opgen_base_.eigblocks());
  }

  virtual BlockedSpinningHamiltonianBase<BaseType>* clone() const =0;

  virtual void interactions(const HamiltonianStore<space_T>&) =0;

  const SpinOpGeneratorBase& opgen_base_;
  const BlockedMatrix<BaseType> Hbase_;
};

template<class BaseType,class OpGen> class BlockedSpinningHamiltonianImpl : public BlockedSpinningHamiltonianBase<BaseType>
{
public:

  BlockedSpinningHamiltonianImpl(const OpGen& opgenv, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const BlockedMatrix<BaseType>& Hbasev)
    : BlockedSpinningHamiltonianBase<BaseType>(opgenv,rotor_speedv,rotor_phasev,rinfo,Hbasev), opgen_(opgenv) {}

  BlockedSpinningHamiltonianImpl(const OpGen& opgenv, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const IntervalSamplerBase& samplerv, const BlockedMatrix<BaseType>& Hbasev)
    : BlockedSpinningHamiltonianBase<BaseType>(opgenv,rotor_speedv,rotor_phasev,rinfo,samplerv,Hbasev), opgen_(opgenv) {}

  BlockedSpinningHamiltonianBase<BaseType>* clone() const { return new BlockedSpinningHamiltonianImpl(*this); }

  void interactions(const HamiltonianStore<space_T>&);

private:
  const OpGen& opgen_;
};

template<class BaseType> class BlockedSpinningHamiltonian : public BaseStructure, public UnaryFunction< BlockedMatrix<BaseType> , double>
{
public:
  typedef typename BlockedSpinningHamiltonianBase<BaseType>::spinning_type spinning_type;

  template<class OpGen> BlockedSpinningHamiltonian(const OpGen&, const HamiltonianStore<space_T>&, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo =MASRotorInfo, const BlockedMatrix<BaseType>& =BlockedMatrix<BaseType>());
    
  template<class OpGen> BlockedSpinningHamiltonian(const OpGen&, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const BlockedMatrix<BaseType>& =BlockedMatrix<BaseType>());

  template<class OpGen> BlockedSpinningHamiltonian(const OpGen&, const HamiltonianStore<space_T>&, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const IntervalSamplerBase&, const BlockedMatrix<BaseType>& =BlockedMatrix<BaseType>());
    
  template<class OpGen> BlockedSpinningHamiltonian(const OpGen&, double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo, const IntervalSamplerBase&, const BlockedMatrix<BaseType>& =BlockedMatrix<BaseType>());

  const SpinOpGeneratorBase& generator() const { return pImpl_->generator(); }

  BlockedMatrix<BaseType> operator()(double t) const {
    BlockedMatrix<BaseType> d;
    operator()(d,t);
    return d;
  }
  void operator()(BlockedMatrix<BaseType>& d, double t) const { pImpl_->operator()(d,t); }

  typedef BlockedMatrix<BaseType> result_type;
  typedef Matrix<BaseType> subresult_type;
  typedef spinning_type suboutput_type;

  const spinning_type& operator()(size_t mzeig,size_t eig) const { return pImpl_->operator()(mzeig,eig); }

  bool operator!() const { return pImpl_->operator!(); }
  size_t usage() const { return pImpl_->usage(); }

  void interactions(const HamiltonianStore<space_T>& Hstore) { pImpl_->interactions(Hstore); }
  void print(std::ostream& ostr) const { 
    ostr << static_cast<const BaseStructure&>(*this);
    pImpl_->print(ostr);
  }

  const IntervalSamplerBase* samplerp() const { return pImpl_->samplerp(); }
  void rotor_phase(double phase_) { return pImpl_->rotor_phase(phase_); }
  double rotor_phase() const { return pImpl_->rotor_phase(); }
  double phase(double t) const { return pImpl_->phase(t); }
  void rotor_speed(double speed_) { pImpl_->rotor_speed(speed_); }
  double period() const { return pImpl_->rotor_period(); }
  
private:
  smartptr<BlockedSpinningHamiltonianBase<BaseType> > pImpl_;
};

  
template<class BaseType> std::ostream& operator<< (std::ostream& ostr, const BlockedSpinningHamiltonian<BaseType>& a)
{
  a.print(ostr);
  return ostr;
}

  template<class BaseType> struct Ham_traits<BlockedSpinningHamiltonian<BaseType> > {
    typedef space_T coupling_type;
    static const bool isconstant=false;
    static const bool allowRF=true;
    static const bool isdiagonal=false;
    static const bool ishomogeneous=true;
  };

  template<class BaseType> struct Propagator_traits<BlockedSpinningHamiltonian<BaseType> > {
    static const bool allowsetH=false;
  };

struct BlockedDiagonalSpinningHamiltonianBase : public RotorHolder_<DiagonalSpinningHamiltonian>
 {
  BlockedDiagonalSpinningHamiltonianBase(const SpinOpGeneratorBase& opgen_basev,double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo =MASRotorInfo)
    : RotorHolder_<DiagonalSpinningHamiltonian>(rotor_speedv,rotor_phasev,rinfo),
    opgen_base_(opgen_basev) {}

  virtual ~BlockedDiagonalSpinningHamiltonianBase() {}

  void operator()(ListList<double>&, double) const;

  const DiagonalSpinningHamiltonian& operator()(size_t mzeig, size_t blk) const
  { return RotorHolder_<DiagonalSpinningHamiltonian>::operator()(opgen_base_.index(mzeig,blk)); }

   virtual BlockedDiagonalSpinningHamiltonianBase* clone() const =0;

  virtual void interactions(const HamiltonianStore<space_T>&) =0;
   
  void print(std::ostream&) const;
  const SpinOpGeneratorBase& generator() const { return opgen_base_; }

  void reset() {
    if (empty())
      ensure();
    else 
      clear();
  }
  
  const SpinOpGeneratorBase& opgen_base_;

   void ensure() {
     RotorHolder_<DiagonalSpinningHamiltonian>::ensure(opgen_base_.actual_mzblocks()*opgen_base_.eigblocks());
   }

   //   void add0_(double, const ListList<double>&);
   //void add_(size_t, const complex&, const ListList<double>&);
};

template<class OpGen> class BlockedDiagonalSpinningHamiltonianImpl : public BlockedDiagonalSpinningHamiltonianBase
{
public:
  BlockedDiagonalSpinningHamiltonianImpl(const OpGen& opgenv,double rotor_speedv, double rotor_phasev, const RotorInfo& rinfo)
    : BlockedDiagonalSpinningHamiltonianBase(opgenv,rotor_speedv,rotor_phasev,rinfo), opgen_(opgenv) {}
    
  BlockedDiagonalSpinningHamiltonianBase* clone() const { return new BlockedDiagonalSpinningHamiltonianImpl(*this); }

  void interactions(const HamiltonianStore<space_T>&);

private:
  const OpGen& opgen_;
};

class BlockedDiagonalSpinningHamiltonian : public UnaryFunction< ListList<double>,double>, public BaseStructure {
public:
  
  template<class OpGen> BlockedDiagonalSpinningHamiltonian(const OpGen& opgen, const HamiltonianStore<space_T>& store,double rotor_speed_, double rotor_phase_, const RotorInfo& rinfo =MASRotorInfo);
    
  template<class OpGen> BlockedDiagonalSpinningHamiltonian(const OpGen& opgen, double rotor_speed_, double rotor_phase_, const RotorInfo& rinfo =MASRotorInfo);

  typedef DiagonalSpinningHamiltonian suboutput_type;
  const DiagonalSpinningHamiltonian& operator()(size_t mzeig,size_t eig) const { return pImpl_->operator()(mzeig,eig); }

  void operator()(ListList<double>& d, double t) const { return pImpl_->operator()(d,t); }
  ListList<double> operator()(double t) const { 
    ListList<double> d;
    operator()(d,t);
    return d;
  }

  bool operator!() const { return pImpl_->operator!(); }
  size_t usage() const { return pImpl_->usage(); }

  const SpinOpGeneratorBase& generator() const { return pImpl_->generator(); }

  void interactions(const HamiltonianStore<space_T>& Hstore) { pImpl_->interactions(Hstore); }
  
  void print(std::ostream& ostr) const {
    ostr << static_cast<const BaseStructure&>(*this);
    pImpl_->print(ostr); }
 
  const IntervalSamplerBase* samplerp() const { return pImpl_->samplerp(); }
  void rotor_phase(double phase_) { return pImpl_->rotor_phase(phase_); }
  double rotor_phase() const { return pImpl_->rotor_phase(); }
  double phase(double t) const { return pImpl_->phase(t); }
  void rotor_speed(double speed_) { pImpl_->rotor_speed(speed_); }
  double period() const { return pImpl_->rotor_period(); }
  const BlockedDiagonalSpinningHamiltonianBase& implementation() const { return *pImpl_; }
  BlockedDiagonalSpinningHamiltonianBase& implementation() { return *pImpl_; }
 
private:
  smartptr<BlockedDiagonalSpinningHamiltonianBase> pImpl_;
};

  inline std::ostream& operator<< (std::ostream& ostr, const BlockedDiagonalSpinningHamiltonian& a)
    {
      a.print(ostr);
      return ostr;
    }

template<> struct Ham_traits<BlockedDiagonalSpinningHamiltonian> { 
  typedef space_T coupling_type; 
  static const bool allowRF=false;
  static const bool isconstant=false;
  static const bool isdiagonal=true;
  static const bool ishomogeneous=false;
}; 

template<> struct Propagator_traits<BlockedDiagonalSpinningHamiltonian> {
  static const bool allowsetH=false;
};

template<typename BaseType> void BlockedSpinningHamiltonianBase<BaseType>::operator()(BlockedMatrix<BaseType>& dest, double t) const
{
  if (!*this)
    throw Undefined("BlockedSpinningHamiltonian()");
  if (dest.empty())
    opgen_base_.create(dest);
  else {
    if (dest.size()!=this->size())
      throw Mismatch("BlockedSpinningHamiltonian");
  }
  const List<spinning_type>& aslist(*this);
  for (size_t k=aslist.size();k--;)
    (aslist(k))(dest(k),t);
}
