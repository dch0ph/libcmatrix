#ifndef LCM_MetaPropagator_h_
#define LCM_MetaPropagator_h_

/*! \file
 \brief  Base objects for high level propagation
*/

class MinimalPropagator {
public:
  MinimalPropagator(double periodv) : period_(periodv) {}
  virtual ~MinimalPropagator() {}
  virtual void operator()(cmatrix&, double, double, size_t, size_t) const =0;
  double period() const { return period_; } //!< return time over which propagator repeats (0 if none)
  template<class T> static double period(const T& source) { return getperiod(source,Bool2Type<Ham_traits<T>::isconstant>()); }

private:
  double period_;
  template<class T> static double getperiod(const T& source, Bool2Type<false>) { return source.period(); }
  template<class T> static double getperiod(const T&, Bool2Type<true>) { return 0.0; }
};

extern Warning<> diagonalsigma0detectmissed_warning;

  struct BaseMetaPropagator : public BaseStructure {
    const SpinOpGeneratorBase& opgen_base_;
    const BaseList<size_t>& blkstr_;  
    mutable smartptr<IntervalSamplerBase> Samplerp; //!< cached local copy of sampler function
    int verbose_;
    int flags_;
    const matrix_partition_set* general_partp_;
    double period_;

    enum { usechebyshev=1, //!< use Chebyshev methods rather than diagonalisation
	   combinepropagators=2, //!< allow repeated propagators to be combined
	   nosynchronisation=4 //!< don't look for synchronisation conditions / use hints
    };

    //! temporary objects - Propagators are assumed not to be shared
    mutable cmatrix tmp;
    mutable List<complex> dU;
    mutable ListList<complex> blocked_dU;
    mutable CommonHolder<rmatrix,cmatrix,BlockedMatrix<double>,BlockedMatrix<complex> > Htmps_;
    mutable CommonHolder<cmatrix,BlockedMatrix<complex> > Utmps_;
    List< smartptr<DiagonalSpinningHamiltonian,false> > Hdiagps_;

    virtual void operator()(BlockedMatrix<complex>&, double, double) const; //!< default provided
    virtual void operator()(cmatrix&, double, double, size_t, size_t) const =0;
    virtual ~BaseMetaPropagator() {}
    virtual BaseMetaPropagator* clone() const =0;
    
    void partitioning(const matrix_partition_set& part) { 
      if (!part.empty())
	general_partp_=&part;
    }
    const matrix_partition_set* partitioning() const { return general_partp_; }
    const matrix_partition* partition(size_t k) const { 
      return general_partp_ ? &((*general_partp_)(k)) : (const matrix_partition*)NULL; }

    template<class TypeH> void try_compress(const TypeH&) {} //do nothing for most Hamiltonians
    void try_compress(const BlockedSpinningHamiltonian<complex>&);

    const SpinOpGeneratorBase& generator() const { return opgen_base_; }
 
    template<class T> BaseMetaPropagator(const T& source, int verbosev, int flagsv)
      : BaseStructure(source),
	opgen_base_(source.generator()),
	blkstr_(opgen_base_.diagonal_structure()),
	verbose_(verbosev),
	flags_(flagsv),
	//	diagonal_partp_(source.generator().partitioning()),
	general_partp_(NULL),
	period_(MinimalPropagator::period(source)) {}

    template<class T> static double getperiod(const T& source, Bool2Type<false>) { return source.period(); }
    template<class T> static double getperiod(const T&, Bool2Type<true>) { return 0.0; }

    void ensure(BlockedMatrix<complex>&) const;
    void ensure(ListList<complex>&) const;
    void ensure(BlockedMatrix<complex>&, size_t, size_t) const; 
    void identity(BlockedMatrix<complex>& U) const {
      opgen_base_.create(U);
      U.identity();
    }

    template<class HType> void propagator_quad_(BlockedMatrix<complex>& U, const HType& H, double t1, double t2, double intdt, Bool2Type<false>) const
    {
      switch (opgen_base_.quadrupole_order()) {
      case 0:
	propagator(U,H,-1,t1,t2,intdt,opgen_base_.Hzeeman());
	break;
      case 2:
	propagator(U,H,-1,t1,t2,intdt,opgen_base_.Hzeeman_structure(),opgen_base_.Hzeeman());
	break;
      default:
	throw InternalError("propagator_quad_");
      }
    }

    template<class HType> void propagator_quad_(BlockedMatrix<complex>&, const HType&, double, double, double, Bool2Type<true>) const
    {
      throw Failed("Second order/exact quadrupolar Hamiltonian requested from diagonal Hamiltonian");
    }

    template<class HType> void common_propagator(BlockedMatrix<complex>& U, const HType& H, double t1, double t2, double intdt, int which) const
    {
      if (which<0) {
	if (!H) {
	  this->identity(U);
	  return;
	}
	if (opgen_base_.effective_quadrupole_order()==1) 
	  propagator(U,H,-1,t1,t2,intdt);
	else
	  propagator_quad_(U,H,t1,t2,intdt,Bool2Type<Ham_traits<HType>::isdiagonal>());
      }
      else {
	size_t mzblk,eigblk;
	indexer_.reverse(mzblk,eigblk,which);
	common_propagator(U,H,mzblk,eigblk,t1,t2,intdt);
      }
    }

    template<class HType> void common_propagator(BlockedMatrix<complex>& U, const HType& H, size_t mzblk, size_t eigblk, double t1, double t2, double intdt) const
    {
      ensure(U,mzblk,eigblk);
      if (!U.empty())
	common_propagator_(U.front(),H,mzblk,eigblk,t1,t2,intdt,Bool2Type<Ham_traits<HType>::isdiagonal>());
    }

    template<class HType> void common_propagator(cmatrix& U, const HType& H, size_t mzblk, size_t eigblk, double t1, double t2, double intdt) const
    {
      if (opgen_base_.size(mzblk,eigblk))
	common_propagator_(U,H,mzblk,eigblk,t1,t2,intdt,Bool2Type<Ham_traits<HType>::isdiagonal>());
      else
	U.clear();
    }

    template<class UType,class HType> void propagator(UType& U, const HType& H,int which, double t1, double t2, double intdt) const
    {
      propagator_(U,H,which,t1,t2,intdt,Bool2Type<Ham_traits<HType>::isconstant>(),Bool2Type<Ham_traits<HType>::ishomogeneous>(),Bool2Type<Ham_traits<HType>::isdiagonal>());
    } 

    template<class UType,class HType,class DType> void propagator(UType& U, const HType& H, int, double t1, double t2, double intdt, const DType& Hzeeman) const
    {
      propagator_exact_(U,H,t1,t2,intdt,Hzeeman,Bool2Type<Ham_traits<HType>::isconstant>());
    } 

    template<class UType,class HType,class LType,class DType> void propagator(UType& U, const HType& H, int, double t1, double t2, double intdt, const LType& blkstr, const DType& Hzeeman) const 
    {
      propagator_second_(U,H,t1,t2,intdt,blkstr,Hzeeman,Bool2Type<Ham_traits<HType>::isconstant>());
    }

    static void set(BlockedMatrix<complex>&, const BlockedMatrix<complex>&, int =1);
  
    virtual void propagators(cmatrix& dUs, size_t n, double t1, double t2,size_t mzeig,size_t blk) const =0;
    void propagators(BaseList<cmatrix>& Us, double t1, double t2, size_t mzeig, size_t blk) const;

    //homogeneous case
    template<class HType,bool isDiag> void propagator_(cmatrix& U, const HType& H, size_t which, double t1, double t2, double intdt, Bool2Type<false>, Bool2Type<true>, Bool2Type<isDiag>) const
    { 
      typedef HomogeneousPropagator<typename HType::result_type> pgen_t;
      pgen_t Ugen(H,intdt,(flags_ & usechebyshev) ? pgen_t::Chebyshev : pgen_t::Diagonalisation,partition(which),verbose_);
      Ugen(U,t1,t2);
    } 

    template<class HType,bool isDiag> void propagator_(BlockedMatrix<complex>& U, const HType& H, int, double t1, double t2, double intdt, Bool2Type<false>, Bool2Type<true>, Bool2Type<isDiag>) const
    { 
      ensure(U);
      size_t mzblk,eigblk;
      typedef HomogeneousPropagator<typename HType::subresult_type> pgen_t;
      const typename pgen_t::mode_t mode=(flags_ & usechebyshev) ? pgen_t::Chebyshev : pgen_t::Diagonalisation;
      for (size_t i=blkstr_.size();i--;) {
	indexer_.reverse(mzblk,eigblk,i);
	pgen_t Ugen(H(mzblk,eigblk),intdt,mode,partition(i),verbose_);
	cmatrix Ui(U(i),mxflag::nondynamic);
	Ugen(Ui,t1,t2);
      }
    } 

    template<class HType> void propagator_(BlockedMatrix<complex>& U, const HType& H, int, double t1, double t2, double, Bool2Type<true>, Bool2Type<false>, Bool2Type<true>) const
    { 
      ::libcmatrix::propagator(blocked_dU,H(),t2-t1);
      full(U,blocked_dU);
    }

    template<class HType> void propagator_exact_(cmatrix& U, const HType& H, double t1, double t2, double intdt, const BaseList<double>& Hzeeman, Bool2Type<false>) const
    { 
      HomogeneousPropagator<typename HType::result_type> Ugen(H,intdt,Hzeeman,verbose_);
      Ugen(U,t1,t2);
    } 

    template<class HType> void propagator_exact_(BlockedMatrix<complex>& U, const HType& H, double t1, double t2, double intdt, const ListList<double>& Hzeeman, Bool2Type<false>) const
    { 
      ensure(U);
      size_t mzblk,eigblk;
      for (size_t i=blkstr_.size();i--;) {
	indexer_.reverse(mzblk,eigblk,i);
	HomogeneousPropagator<typename HType::subresult_type> Ugen(H(mzblk,eigblk),intdt,Hzeeman(i),verbose_);
	cmatrix Ui(U(i));
	Ugen(Ui,t1,t2);
      }
    } 

    template<class HType> void propagator_(BlockedMatrix<complex>& U, const HType& H, int, double t1, double t2, double, Bool2Type<false>, Bool2Type<false>, Bool2Type<true>) const
    { 
      ensure(blocked_dU);
      size_t mzblk,eigblk;
      for (size_t i=blkstr_.size();i--;) {
	indexer_.reverse(mzblk,eigblk,i);
	DiagonalInhomogeneousPropagator Ugen(H(mzblk,eigblk),verbose_);
	BaseList<complex> dUi(blocked_dU(i));
	Ugen(dUi,t1,t2);
	full(U(i),dUi);
      }
    } 

    template<class UType,class HType,class DType> void propagator_exact_(UType& U, const HType& H, double t1, double t2, double, const DType& Hzeeman, Bool2Type<true>) const
    { 
      propagator_ns(U,H,t2-t1,Hzeeman);
    } 

    template<typename T> void propagator_exact_(BlockedMatrix<complex>& U, const Matrix<T>& H, double t1, double t2, double, const ListList<double>& Hzeeman, Bool2Type<true>) const
    { 
      U.duplicate_structure(H);
      propagator_ns(U.front(),H,t2-t1,Hzeeman.row());
    } 

    template<class UType,class HType,class LType,class DType> void propagator_second_(UType& U, const HType& H, double t1, double t2, double, const LType& blkstr, const DType& Hzeeman, Bool2Type<true>) const
    { 
      propagator_second(U,H,t2-t1,blkstr,Hzeeman);
    }

    template<class HType,class LType,class DType> void propagator_second_(cmatrix& U, const HType& H, double t1, double t2, double intdt, const LType& blkstr, const DType& Hzeeman, Bool2Type<false>) const
    { 
      HomogeneousPropagator_second<typename HType::result_type> Ugen(H,intdt,blkstr,Hzeeman,verbose_);
      Ugen(U,t1,t2);
    }

    template<class HType,class LType,class DType> void propagator_second_(BlockedMatrix<complex>& U, const HType& H, double t1, double t2, double intdt, const LType& blkstr, const DType& Hzeeman, Bool2Type<false>) const
    { 
      ensure(U);
      size_t mzblk,eigblk;
      for (size_t i=blkstr_.size();i--;) {
	indexer_.reverse(mzblk,eigblk,i);
	HomogeneousPropagator_second<typename HType::subresult_type> Ugen(H(mzblk,eigblk),intdt,blkstr(i),Hzeeman(i),verbose_);
	cmatrix Ui(U(i));
	Ugen(Ui,t1,t2);
      }
    }

    template<class T> void step_propagator(cmatrix& U, const Matrix<T>& H, size_t, double dt) const;
    template<class T> void step_propagator(BlockedMatrix<complex>& U, const BlockedMatrix<T>& H, int, double dt) const;

    //inhomogeneous case
    template<class HType> void propagator_(cmatrix& U, const HType& H, int, double t1, double t2, double, Bool2Type<false>, Bool2Type<false>, Bool2Type<true>) const
    { 
      DiagonalInhomogeneousPropagator Ugen(H,verbose_);
      Ugen(dU,t1,t2);
      full(U,dU);
    }  
    
    template<class HType,bool isHomo> void propagator_(BlockedMatrix<complex>& U, const HType& H, int which, double t1, double t2, double, Bool2Type<true>, Bool2Type<isHomo>, Bool2Type<false>) const
    { 
      step_propagator(U,H(),which,t2-t1);
    }

    template<class HType,bool isHomo> void propagator_(cmatrix& U, const HType& H, size_t which, double t1, double t2, double, Bool2Type<true>, Bool2Type<isHomo>, Bool2Type<false>) const
    { 
      step_propagator(U,H,which,t2-t1);
    }

    template<typename T,bool isHomo> void propagator_(BlockedMatrix<complex>& U, const Matrix<T>& H, size_t which, double t1, double t2, double, Bool2Type<true>, Bool2Type<isHomo>, Bool2Type<false>) const
    { 
      U.duplicate_structure(H);
      step_propagator(U.front(),H,which,t2-t1);
    }

    template<class HType,bool isHomo> void propagator_(cmatrix& U, const HType& H, int, double t1, double t2, double, Bool2Type<true>, Bool2Type<isHomo>, Bool2Type<true>) const
    { 
      ::libcmatrix::propagator(dU,H,t2-t1);
      full(U,dU);
    }

    template<class HType> void steppropagator(cmatrix& U,const HType& H, size_t which, double dt) const
    {
      size_t mzblk,eigblk;
      indexer_.reverse(mzblk,eigblk,which);
      switch (opgen_base_.effective_quadrupole_order()) {
      case 0:
	::libcmatrix::propagator_ns(U,H(which),dt,opgen_base_.Hzeeman(mzblk,eigblk));
	break;
      case 1:
	::libcmatrix::propagator(U,H(which),dt);
	break;
      case 2:
	::libcmatrix::propagator_second(U,H(which),dt,opgen_base_.Hzeeman_structure(mzblk,eigblk),opgen_base_.Hzeeman(mzblk,eigblk));
	break;
      default:
	throw InternalError("steppropagator: unhandled quadrupole order");
      }
    }

    template<class HSysType, class HrfType> void steppropagator_added_ns(cmatrix& U,const HSysType& Hsys, const HrfType& Hrf, size_t which, double dt, Bool2Type<false>) const
    {
      size_t mzblk,eigblk;
      indexer_.reverse(mzblk,eigblk,which);
      switch (opgen_base_.effective_quadrupole_order()) {
      case 0:
	::libcmatrix::propagator_ns_Hrf(U,Hsys(which),Hrf(which),dt,opgen_base_.Hzeeman(mzblk,eigblk));
	break;
      case 2:
	::libcmatrix::propagator_second_Hrf(U,Hsys(which),Hrf(which),dt,opgen_base_.Hzeeman_structure(mzblk,eigblk),opgen_base_.Hzeeman(mzblk,eigblk));
	break;
      default:
	throw InternalError("steppropagator_added_ns: unhandled quadrupole order");
      }
    }

    template<class HSysType, class HrfType> void steppropagator_added_ns(cmatrix&, const HSysType&, const HrfType&, size_t, double, Bool2Type<true>) const { throw InternalError("Non-secular Hamilonian requested from diagonal Hamiltonian"); }

    template<class HType> void steppropagator(BlockedMatrix<complex>& U, const HType& H, double dt) const
    {
      if (opgen_base_.issecular())
	step_propagator(U,H,-1,dt);
      else {
	switch (opgen_base_.quadrupole_order()) {
	case 0:
	  ::libcmatrix::propagator_ns(U,H,dt,opgen_base_.Hzeeman());
	  break;
	case 2:
	  ::libcmatrix::propagator_second(U,H,dt,opgen_base_.Hzeeman_structure(),opgen_base_.Hzeeman());
	  break;
	default:
	  throw InternalError("steppropagator(2): unhandled quadrupole order");
	}
      }
    }

    template<class HsysType, class HrfType> void steppropagator_added_ns(BlockedMatrix<complex>& U, const HsysType& Hsys, const HrfType& Hrf, double dt, Bool2Type<false>) const
    {
      switch (opgen_base_.quadrupole_order()) {
      case 0:
	::libcmatrix::propagator_ns_Hrf(U,Hsys,Hrf,dt,opgen_base_.Hzeeman());
	break;
      case 2:
	::libcmatrix::propagator_second_Hrf(U,Hsys,Hrf,dt,opgen_base_.Hzeeman_structure(),opgen_base_.Hzeeman());
	break;
      default:
	throw InternalError("steppropagator_added_ns(2): unhandled quadrupole order");
      }
    }

    template<class HsysType, class HrfType> void steppropagator_added_ns(BlockedMatrix<complex>&, const HsysType&, const HrfType&, double, Bool2Type<true>) const { throw InternalError("Non-secular Hamilonian requested from diagonal Hamiltonian"); }

    template<class HType> void steppropagator(BlockedMatrix<complex>& U, const HType& H, int which, double dt) const
    {
      if (which<0) {
	steppropagator(U,H,dt);
	return;
      }
      const size_t uwhich(which);
      size_t mzblk,eigblk;
      indexer_.reverse(mzblk,eigblk,uwhich);
      ensure(U,mzblk,eigblk);
      steppropagator(U.front(),H,uwhich,dt);
    }

    template<class HsysType, class HrfType, bool DiagType> void steppropagator_added_ns(BlockedMatrix<complex>& U, const HsysType& Hsys, const HrfType& Hrf, int which, double dt, Bool2Type<DiagType> dtype) const
    {
      if (which<0) {
	steppropagator_added_ns(U,Hsys,Hrf,dt,dtype);
	return;
      }
      const size_t uwhich(which);
      size_t mzblk,eigblk;
      indexer_.reverse(mzblk,eigblk,uwhich);
      ensure(U,mzblk,eigblk);
      steppropagator_added_ns(U.front(),Hsys,Hrf,uwhich,dt,dtype);
    }

    template<class UType,class HType,class HaddT> void Hadded_propagator_(UType& U, const HType& Hsys, const HaddT& Hadd, size_t which, double t1, double t2, double,Bool2Type<true>) const
     {
       if (Hsys.empty()) {
	 step_propagator(U,Hadd,which,t2-t1);
	 return;
       }
       if (opgen_base_.issecular()) {
	 HaddT& Htmp(Htmps_(Type2Type<HaddT>()));
	 Htmp=Hadd;
	 // if (!Hsys.empty())
	 Htmp+=Hsys;
	 steppropagator(U,Htmp,which,t2-t1);
       }
       else
	 steppropagator_added_ns(U,Hsys,Hadd,which,t2-t1,Bool2Type<Ham_traits<HType>::isdiagonal>());
     }

    template<class HType, class HaddT> void Hadded_propagator_(BlockedMatrix<complex>& U, const HType& Hsys, const HaddT& Hadd, int which, double t1, double t2, double, Bool2Type<true>) const
      {
	if (!Hsys) {
	  step_propagator(U,Hadd,which,t2-t1);
	  return;
	}
	if (opgen_base_.issecular()) {
	  HaddT& Htmp(Htmps_(Type2Type<HaddT>()));
	  Htmp=Hadd;
	  // if (!!Hsys)
	  Htmp+=Hsys();
	  steppropagator(U,Htmp,which,t2-t1);
	}
	else
	  steppropagator_added_ns(U,Hsys(),Hadd,which,t2-t1,Bool2Type<Ham_traits<HType>::isdiagonal>());
      }
   
//     template<typename T, class HaddT> void Hadded_propagator_(BlockedMatrix<complex>& U, const Matrix<T>& Hsys, const HaddT& Hadd, size_t which, double t1, double t2, double, const ListList<double>& Hzeeman,Bool2Type<true>) const
//       {
// 	if (Hadd.size()>1)
// 	  throw Failed("BaseMetaPropagator: can't squeeze BlockedMatrix into Matrix");
// 	HaddT& Htmp(Htmps_(Type2Type<HaddT>()));
// 	Htmp=Hadd;
// 	if (!!Hsys)
// 	  Htmp+=Hsys();
// 	U.duplicate_structure(Htmp);
// 	if (Hzeeman.empty()) 
// 	  step_propagator(U.front(),Htmp.front(),which,t2-t1);
// 	else {
// 	  Htmp+=Hzeeman;
// 	  ::libcmatrix::propagator_ns(U.front(),Htmp.front(),t2-t1,Hzeeman.row());
// 	}
//       }

  template<class HType,class subHcType,class subHType> static void getadd_(subHcType& Htmp, const HType& H, double t,const subHcType& Hadd, subHType& Htmpr)
    {
      H(Htmpr,t);
      if (Htmpr.empty())
	Htmp=Hadd;
      else
	add(Htmp,Htmpr,Hadd);
    }

    template<class HType> static void getadd_(cmatrix& Htmp, const HType& H, double t,const cmatrix& Hadd, List<double>& Htmpr)
    {
      H(Htmpr,t);
      Htmp=Hadd;
      if (!(Htmpr.empty()))
	Htmp+=Htmpr;
    }

  template<class HType,class subHcType> static void getadd_(subHcType& Htmp, const HType& H,double t,const subHcType& Hadd, subHcType&)
    {
      H(Htmp,t);
      Htmp+=Hadd;
    }
  
    template<class UType,class HType,class HaddType> void Hadded_propagator_(UType& U, const HType& Hsys,const HaddType& Hadd, int which, double t1, double t2, double intdt, Bool2Type<false>) const
   {
     const double maxdttol=intdt+lcm_MAS_timingtol;

     HaddType& Htmp(Htmps_(Type2Type<HaddType>()));
     UType& Utmp(Utmps_(Type2Type<UType>()));
     if (&U==&Utmp)
       throw ArgumentClash("Hadded_propagator_");

     typename HType::result_type Htmpr;

     if (verbose_) {
       std::cout << "Evaluating overall propagator from ";
       prettyprint_time(t1) << " to ";
       prettyprint_time(t2) << '\n';
     }
      
     IntervalSamplerBase* usesamplerp=&lcm_defsampler;
     if (Hsys.samplerp() ) {
       if (!Samplerp)
	 Samplerp.reset( (Hsys.samplerp())->clone());       
       usesamplerp=Samplerp.get();
     }

      bool first=true;
      for (;;) {
	double dt=t2-t1;
	if (fabs(dt)<lcm_MAS_timingtol) 
	  break;
	if (dt<0.0)
	  throw InternalError("Negative delay!");
	if (dt>maxdttol) 
	  dt=intdt;
	const double midt=(*usesamplerp)(t1,t1+dt);
	UType& Udest(first ? U : Utmp);
	
	if (opgen_base_.issecular()) {
	  getadd_(Htmp,Hsys,midt,Hadd,Htmpr);
	  if (verbose_>2) {
	    std::cout << "H at 'midpoint' of " << (1e6*t1) << " and " << (1e6*(t1+dt)) << ", t=" << (1e6*midt) << " us:\n";
	    std::cout<< Htmp << '\n';
	  }
	  steppropagator(Udest,Htmp,which,dt);
	}
	else {
	  Hsys(Htmpr,midt);	  
	  steppropagator_added_ns(Udest,Htmpr,Hadd,which,dt,Bool2Type<Ham_traits<HType>::isdiagonal>());
	}
	if (verbose_>2) {
	  std::cout << "U(";
	  prettyprint_time(t1) << ',';
	  prettyprint_time(t1+dt) << ")\n" << Udest << '\n';
	}
	if (!first)
	  U&=Udest;
	else
	  first=false;
	t1+=dt;
      }
      assert(!first);
    }

    template<class HsysT, class HaddT> void Hadded_propagator(BlockedMatrix<complex>& U, const HsysT& Hsys, const HaddT& Hcur, double t1, double t2, double intdt, int which) const
    {
      if (verbose_>1) {
	std::cout << "Current RF at ";
	prettyprint_time(t1) << ":\n";
	if (which<0)
	  std::cout << Hcur;
	else
	  std::cout << Hcur(which);
      }

      if (which<0) {
// 	if (opgen_base_.effective_quadrupole_order()==2)
// 	  throw Failed("BaseMetaPropagator: unhandled quadrupole order");
	Hadded_propagator_(U,Hsys,Hcur,-1,t1,t2,intdt,Bool2Type<Ham_traits<HsysT>::isconstant>());
      }
      else {
	size_t mzblk,eigblk;
	indexer_.reverse(mzblk,eigblk,which);
	ensure(U,mzblk,eigblk);
	//const size_t uwhich(which);
	if (!Hsys)
	  step_propagator(U.front(),Hcur(size_t(which)),size_t(which),t2-t1);
	else {
	  Hadded_propagator_(U,Hsys,Hcur,which,t1,t2,intdt,Bool2Type<Ham_traits<HsysT>::isconstant>());
// 	  switch (opgen_base_.effective_quadrupole_order()) {
// 	  case 1:
// 	    Hadded_propagator_(U.front(),Hsys(mzblk,eigblk),Hcur(uwhich),which,t1,t2,intdt,Bool2Type<Ham_traits<HsysT>::isconstant>());
// 	    break;
// 	  case 0:
// 	    Hadded_propagator_(U.front(),Hsys(mzblk,eigblk),Hcur(uwhich),which,t1,t2,intdt,Bool2Type<Ham_traits<HsysT>::isconstant>());
// 	    break;
// 	  default:
// 	    throw Failed("BaseMetaPropagator: unhandled quadrupole order");
// 	  }
	}
      }
      if (verbose_>1) {
	std::cout << "Propagator from ";
	prettyprint_time(t1) << " to ";
	prettyprint_time(t2) << '\n' << U;      
      }
    }

    template<class HType> void common_propagator_(cmatrix& U, const HType& H, size_t mzblk, size_t eigblk, double t1, double t2, double intdt,Bool2Type<true>) const //diagonal, so ignore quad order
    {
      assert(opgen_base_.effective_quadrupole_order()==1);
      propagator(U,H(mzblk,eigblk),indexer_(mzblk,eigblk),t1,t2,intdt);
    }

    //not diagonal, hence consider quad order
    template<class HType> void common_propagator_(cmatrix& U, const HType& H, size_t mzblk, size_t eigblk, double t1, double t2, double intdt, Bool2Type<false>) const
    {
      const size_t which(indexer_(mzblk,eigblk));
      switch (opgen_base_.effective_quadrupole_order()) {
      case 1:
	propagator(U,H(mzblk,eigblk),which,t1,t2,intdt);
	break;
      case 0:
	propagator(U,H(mzblk,eigblk),which,t1,t2,intdt,opgen_base_.Hzeeman(mzblk,eigblk));
	break;
      case 2:
	if (Hdiagps_.empty()) 
	  BaseMetaPropagator::propagator(U,H(mzblk,eigblk),which,t1,t2,intdt,opgen_base_.Hzeeman_structure(mzblk,eigblk),opgen_base_.Hzeeman(mzblk,eigblk));
	else {
	  DiagonalSpinningHamiltonian& Hdiag( *(Hdiagps_(which) ));
	  DiagonalInhomogeneousPropagator Ugen(Hdiag,verbose_);
	  Ugen(dU,t1,t2);
	  full(U,dU);
	}
	break;
      default:
	throw Failed("BaseMetaPropagator: unhandled quadrupole order");
      }
    }

    }; 

  template<class TypeH> class RawMetaPropagator : public BaseMetaPropagator {
  public:
    typedef typename TypeH::result_type baseH_type;

    RawMetaPropagator(const TypeH&, double intdtv, int verbosev, int flagsv);

    BaseMetaPropagator* clone() const { return new RawMetaPropagator<TypeH>(*this); }

    void propagators_(cmatrix& dUs, size_t n, double t1, double t2,size_t mzeig,size_t blk,Bool2Type<true>) const
    {
      DiagonalInhomogeneousPropagator propgen(H_(mzeig,blk),verbose_);
      ::libcmatrix::propagators(dUs,n,propgen,t1,t2);
    }

    void propagators_(cmatrix&, size_t, double, double,size_t,size_t,Bool2Type<false>) const
    { throw Failed("MetaPropagators:: can't create diagonal propagators from non-diagonal Hamiltonian"); }

    void propagators(cmatrix& dUs, size_t n, double t1, double t2,size_t mzeig,size_t blk) const
    {
      propagators_(dUs,n,t1,t2,mzeig,blk,Bool2Type<Ham_traits<TypeH>::isdiagonal && !Ham_traits<TypeH>::isconstant>());
    }

    void propagators(BaseList<cmatrix> Us, double t1, double t2, size_t mzblk, size_t eigblk) const
    {
      if (opgen_base_.size(mzblk,eigblk)==0) {
	for (size_t n=Us.length();n--;)
	  Us(n).clear();
      }
      else
	BaseMetaPropagator::propagators(Us,H_(mzblk,eigblk),t1,t2,intdt_);
    }

    void operator()(cmatrix& U, double t1, double t2, size_t mzblk, size_t eigblk) const
    {
      BaseMetaPropagator::common_propagator(U,H_,mzblk,eigblk,t1,t2,intdt_);
    }

  private:
    const TypeH& H_;
    double intdt_;
  };


//! bind Propagator method to lower-level PropGen_t from NMR.h
class PropagatorAdaptor : public PropGen_t {
public:
  template<class T> PropagatorAdaptor(const T& metapgenv, size_t mzblkv, size_t eigblkv)
    : metapgen_(metapgenv), size_(metapgenv.generator().size(mzblkv,eigblkv)),
      mzblk_(mzblkv), eigblk_(eigblkv) {}

  void operator()(cmatrix& U, double t1, double t2) {
    metapgen_(U,t1,t2,mzblk_,eigblk_);
  }

  cmatrix operator() (double t1,double t2) const {
    cmatrix U(mxflag::temporary);
    metapgen_(U,t1,t2,mzblk_,eigblk_);
    return U;
  }
  
  PropGen_t* clone() const { return new PropagatorAdaptor(*this); }
  size_t size() const { return size_; }

private:
  const MinimalPropagator& metapgen_;
  size_t size_;
  size_t mzblk_,eigblk_;
};

template<class TypeH> RawMetaPropagator<TypeH>::RawMetaPropagator(const TypeH& Hv, double intdtv, int verbosev, int flagsv)
    : BaseMetaPropagator(Hv,verbosev,flagsv),
    H_(Hv), intdt_(intdtv)
    {
      if (intdt_==0.0) {
	if (Ham_traits<TypeH>::ishomogeneous)
	  throw InvalidParameter("MetaPropagator: integration timestep unspecified for homogeneous Hamiltonian");
	if (flags_ & usechebyshev)
	  throw InvalidParameter("MetaPropagator: can't use Chebyshev propagation without integration time step");
      }
      if (opgen_base_.effective_quadrupole_order()==2)
	try_compress(Hv);
    }

class MetaPropagator : public MinimalPropagator {
  public:
    typedef cmatrix suboutput_type;
    typedef BlockedMatrix<complex> output_type;
    
    const SpinOpGeneratorBase& generator() const { return objp_->generator(); }

    template<class HType> MetaPropagator(const HType& Hv, double intdtv =0, int verbosev =0, int flagsv =0)
      : MinimalPropagator(MinimalPropagator::period(Hv)),
	objp_(new RawMetaPropagator<HType>(Hv,intdtv,verbosev,flagsv)) {}

    void propagators(BaseList<cmatrix> Us, double t1, double t2, size_t mzblk, size_t eigblk) const {
      objp_->propagators(Us,t1,t2,mzblk,eigblk);
    }

    void propagators(cmatrix& dUs, size_t n, double t1, double t2,size_t mzblk,size_t blk) const {
      objp_->propagators(dUs,n,t1,t2,mzblk,blk);
    }

    void operator()(BlockedMatrix<complex>& U, double t1, double t2) const {
      objp_->operator()(U,t1,t2);
    }

    void operator()(cmatrix& U, double t1, double t2, size_t mzblk, size_t eigblk) const {
      objp_->operator()(U,t1,t2,mzblk,eigblk);
    }

    void partitioning(const matrix_partition_set& partv) { objp_->partitioning(partv); }

  private:
    smartptr<BaseMetaPropagator> objp_;
  };


  template<class Source, class PropUse, class TypeSD> class BaseMetaObj {
  public:
    //Propagator-based
    template<class PropBase> BaseMetaObj(const PropBase&, const Source& pgenv, int nobsv, double startv, double endv, const TypeSD& sigma0v, const TypeSD& detectv,int verbosev);

    //    template<class PropBase> BaseMetaObj(const PropBase&, const Source& Hv, double periodv, const TypeSD& sigma0v, const TypeSD& detectv, int verbosev);

  void eigenvalue(size_t seleig) {
    if (seleig>=opgen_base_.eigblocks())
      throw InvalidParameter("BaseMetaObj::eigenvalue out of range");
    eigoffset_=seleig;
    objlist_.create(1,objlist_.front());
    useblocks_=1;
  }

//     void reduction_factor(int rv)
//     { 
//       check_reduction_valid(Int2Type<SignalGenerator_traits<PropUse>::input_dimensionality>(), Bool2Type<SignalGenerator_traits<PropUse>::allowU>());
//       reduction_factor_=rv;
//       if (rv<1)
// 	throw InvalidParameter("BaseMetaObj::reduction_factor cannot be <1");
//     }

//     size_t reduction_factor() { return reduction_factor_; }

    size_t current_mzeig() const { return mzeigSD_; }
        
  protected:

    const Source* Hpgenp_; 
    const SpinOpGeneratorBase& opgen_base_;
    const TypeSD sigma0_;
    const TypeSD detect_;
    int verbose_;
    double start_,end_;
    mutable DynamicList<PropUse> objlist_;
    block_pattern blkspec_;
  size_t useblocks_,mzblocks_,mzblocksSD_;
  size_t eigoffset_;
    bool ismiddle_;
    bool havevalidpgen_;
    DynamicList<cmatrix> passlist_;
    cmatrix dpasslist_;
    double period_;
    double obslarmor_;
    char dirn_;
    size_t skip_;
    size_t mzeig_;
    size_t nobs_;
    bool usemzsym_;
    smartptr<block_pattern::iterator,false> mziterp_;
    size_t lastrind_,lastcind_;
    size_t mzeigSD_;

    mutable List<double> seigs_;

    bool isspecialcase(complex& scale, size_t mzeig, size_t blk) const {
      static const complex zero(0.0);
      const complex& sigma0scale(sigma0_.identityscale(mzeig,blk));
      const complex& detectscale(detect_.identityscale(mzeig,blk)); //!< fixed from sigma0 19/11/2013
      if ((sigma0scale==zero) || (detectscale==zero))
	return false;
      scale*=sigma0scale*detectscale;
      if (sigma0_.ishermitian() && detect_.ishermitian())
	scale*=2.0; //!< add 19/11/2013
      return true;
    }
	
    bool doublemzblock() const 
    { return (usemzsym_ && !ismiddle_); }

    void shuffle(char dirn) {
      switch (dirn) {
      case 'C':
	lastrind_=lastcind_;
	break;
      case 'R':
	lastcind_=lastrind_;
	break;
      default:
	throw InvalidParameter("MetaPropagator::shuffle");
      }
      for (size_t i=useblocks_;i--;)
	objlist_(i).shuffle(dirn);
    }

    void larmor(double larmorv) {
      obslarmor_=larmorv;
    }

    bool inciter(bool& makeherm);

    template<class HType> void set_H_(PropUse& obj,Bool2Type<true>,const HType& H,char which) const
      {
	if (which=='M')
	  obj.set_H(H,period_);
	else
	  obj.set_H(which,H,period_);
      }
    
    template<class VType> void set_H_(PropUse& obj,Bool2Type<true>,const VType& V,const BaseList<double>& eigs,char which) const
      {
	if (which=='M')
	  obj.set_H(V,eigs,period_);
	else
	  obj.set_H(which,V,eigs,period_);
      }
    
    template<class HType> void set_H_(PropUse& obj,Bool2Type<false>,const HType& H,char which) const
      {
	if (which=='M')
	  obj.set_H(H);
	else
	  obj.set_H(which,H);
      }
  
    template<class VType> void set_H_(PropUse& obj,Bool2Type<false>,const VType& V,const BaseList<double>& eigs,char which) const
      {
	if (which=='M')
	  obj.set_H(V,eigs);
	else
	  obj.set_H(which,V,eigs);
      }
    
    template<bool isTD,bool isDiag> void set_H_2nd_(PropUse& obj, Bool2Type<isTD> istd, size_t mzeig, size_t blkeig, char which, Bool2Type<isDiag>) const
    { // assume that diagonal Hamiltonians are 2nd order Hamiltonians that have been `reduced'
       set_H_(obj,istd,(*Hpgenp_)(mzeig,blkeig),which);
       //       throw Failed("BaseMetaObj: 2nd order not valid for this Hamiltonian");
     }

    template<bool isTD> void set_H_2nd_(PropUse&, Bool2Type<isTD>, size_t, size_t,char, Bool2Type<false>) const;

    void set_HUs(size_t,char);

    template<class T, bool isTD> static void set_U(PropUse&, T, Bool2Type<isTD>, Bool2Type<false>, char)
    { throw Failed("No set_U method!\n"); }

//     template<int N, bool allowU> static void check_reduction_valid(Int2Type<N>, Bool2Type<allowU>)
//     { throw Failed("Reduction factor not valid for this object"); }

//     //! only allowed to set reduction factor for propagation using single U
//     static void check_reduction_valid(Int2Type<2>, Bool2Type<true>) {}
    
    void set_U(PropUse&, Int2Type<2>, Bool2Type<false>, Bool2Type<true>, char);
    void set_U(PropUse&, Int2Type<2>, Bool2Type<true>, Bool2Type<true>, char);
    void set_U(PropUse&, Int2Type<3>, Bool2Type<false>, Bool2Type<true>, char);
    void set_U(PropUse&, Int2Type<3>, Bool2Type<true>, Bool2Type<true>, char);

    //! Special case taking propagator generators e.g. AsynchronousFID
    void set_HUspecial(PropUse& obj_, Bool2Type<true>, size_t mzeig, size_t blk, char which) {
      const PropagatorAdaptor adaptor(*Hpgenp_,mzeig,blk);
      if (which=='M')
	obj_.set_Ugenerator(adaptor);
      else
	obj_.set_Ugenerator(which,adaptor);
    }

    //! Special case of inhomogeneous (spinning) Hamiltonian - TD
  void set_HUspecial(PropUse& obj_, Bool2Type<false>, size_t mzeig, size_t blk, char which) {
    Hpgenp_->propagators(dpasslist_,nobs_,start_,end_,mzeig,blk);
    if (which=='M')
      obj_.set_Us(dpasslist_);
    else
      obj_.set_Us(which,dpasslist_);
  }

  bool set_HU(PropUse& obj_, Bool2Type<false>, Bool2Type<false>, Bool2Type<true>, size_t mzeig, size_t blk, char which) {
    set_HUspecial(obj_, Bool2Type< SignalGenerator_traits<PropUse>::usegenerator>(),mzeig,blk,which);
    return true;
  }
  //not TD
  bool set_HU(PropUse& obj_, Bool2Type<false>, Bool2Type<false>, Bool2Type<false>, size_t mzeig, size_t blk, char which) {
    Hpgenp_->propagators(dpasslist_,nobs_,start_,end_,mzeig,blk);
    if (which=='M')
      obj_.set_Us(dpasslist_,period_); //non-trivial overhead - would be better to pass DiagSpinHam directly
    else
      obj_.set_Us(which,dpasslist_,period_);
    return true;
  }

    //allowH, allowU
    template<bool isTD> bool set_HU(PropUse&, Bool2Type<false>, Bool2Type<true>, Bool2Type<isTD>, size_t mzeig, size_t blk, char) { 
      passlist_.create(nobs_);
      Hpgenp_->propagators(passlist_,start_,end_, mzeig, blk);
      if (verbose_>1)
	std::cout << "Propagators for mz=" << mzeig << "  eig=" << blk << '\n' << passlist_;
      return false;
    }

    template<bool isTD,bool allowU> bool set_HU(PropUse&, Bool2Type<true>, Bool2Type<allowU>, Bool2Type<isTD>, size_t, size_t, char);

  private:
    void common_init(int nobsv, double periodv, bool havevalidpgenv =false);
    void initSD();
  };


  template<class Source, class PropUse, class TypeSD> class RawMetaSpectrum
    : public BaseMetaObj<Source,PropUse,TypeSD>, public SpectrumIterator<typename PropUse::amplitude_type> {
  public:

    RawMetaSpectrum(const PropUse& objv, const Source& genv, int nobsv, double startv, double endv, const TypeSD& sigma0v, const TypeSD& detectv, int verbosev)
      : BaseMetaObj<Source,PropUse,TypeSD>(objv,genv,nobsv,startv,endv,sigma0v,detectv,verbosev) {
      reset();
    }


     using BaseMetaObj<Source,PropUse,TypeSD>::dirn_;
     using BaseMetaObj<Source,PropUse,TypeSD>::lastrind_;
     using BaseMetaObj<Source,PropUse,TypeSD>::lastcind_;
     using BaseMetaObj<Source,PropUse,TypeSD>::inciter;
     using BaseMetaObj<Source,PropUse,TypeSD>::mziterp_;
      using BaseMetaObj<Source,PropUse,TypeSD>::ismiddle_;
    using BaseMetaObj<Source,PropUse,TypeSD>::doublemzblock;
      using BaseMetaObj<Source,PropUse,TypeSD>::useblocks_;
      using BaseMetaObj<Source,PropUse,TypeSD>::mzblocks_;
      using BaseMetaObj<Source,PropUse,TypeSD>::mzblocksSD_;
    using BaseMetaObj<Source,PropUse,TypeSD>::eigoffset_;
      using BaseMetaObj<Source,PropUse,TypeSD>::set_HUs;
      using BaseMetaObj<Source,PropUse,TypeSD>::opgen_base_;
      using BaseMetaObj<Source,PropUse,TypeSD>::verbose_;
      using BaseMetaObj<Source,PropUse,TypeSD>::shuffle;
      using BaseMetaObj<Source,PropUse,TypeSD>::sigma0_;
      using BaseMetaObj<Source,PropUse,TypeSD>::detect_;
      using BaseMetaObj<Source,PropUse,TypeSD>::objlist_;
      using BaseMetaObj<Source,PropUse,TypeSD>::blkspec_;
    
    bool operator() (typename PropUse::amplitude_type&, double&);

     void eigenvalue(size_t val) {
       // if (!finished_)
       //throw Failed("MetaSpectrum::eigenvalue: can't change while active");
       BaseMetaObj<Source,PropUse,TypeSD>::eigenvalue(val);
       reset();
     }

  private:
    size_t totweight_;
    bool makeherm_;
    //    size_t mzeigSD_;
    PropUse* curobjp_;
    bool finished_;
    int blockpos_;
    typename PropUse::amplitude_type nextamp_;
    double nextfreq_;
    bool havenext_;

    RawMetaSpectrum(const RawMetaSpectrum&); //disable copy ops
    RawMetaSpectrum& operator= (const RawMetaSpectrum&);
     
    void set_sigma0dets(size_t);
    void set_sigma0det(size_t, size_t, Bool2Type<true>, Bool2Type<false>);    
    void set_sigma0det(size_t, size_t, Bool2Type<true>, Bool2Type<true>);    
    void set_sigma0det(size_t, size_t, Bool2Type<false>, Bool2Type<true>);
  
    void reset_blocks(size_t newblock) {
      blockpos_=newblock;
      curobjp_=&(this->objlist_(blockpos_));
      const size_t cureig(eigoffset_+blockpos_);
      totweight_=this->opgen_base_.eigweight(cureig)*((this->doublemzblock)() ? 2 : 1);
      if (this->verbose_>1)
	std::cout << "  Eigenvalue: " << cureig  << "  Weight: " << totweight_ << '\n';
    }
    
    bool next() {
      if (!inciter(makeherm_))
	return false;
      set_sigma0dets(BaseMetaObj<Source,PropUse,TypeSD>::current_mzeig());
      //      mzeigSD_++;
      return true;
    }
    void reset();
  };

  template<typename AmpType> class MetaSpectrum : public SpectrumIterator<AmpType> {
  public:
    template<class PropUse,class MetaProp,class TypeSD>
    MetaSpectrum(const PropUse& puse, const MetaProp& genv, int nobsv, double startv, double endv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(), int verbosev =0)
      : metaobjp_(new RawMetaSpectrum<MetaProp,PropUse,TypeSD>(puse,genv,nobsv,startv,endv,sigma0v,detectv,verbosev)) {}
        
    template<class PropUse,class TypeH, class TypeSD>
    MetaSpectrum(const PropUse& puse,const TypeH& Hv, double periodv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(), int verbosev =0)
      : metaobjp_(new RawMetaSpectrum<TypeH,PropUse,TypeSD>(puse,Hv,periodv,sigma0v,detectv,verbosev)) {}

    template<class PropUse,class TypeH,class TypeSD>
    MetaSpectrum(const PropUse& puse,const TypeH& Hv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(), int verbosev =0)
      : metaobjp_(new RawMetaSpectrum<TypeH,PropUse,TypeSD>(puse,Hv,0.0,sigma0v,detectv,verbosev)) {}
            
    bool operator()(AmpType& amp, double& freq) 
    { return (*metaobjp_.get())(amp,freq); }

    void larmor(double larmorv)
    { this->metaobjp_->larmor(larmorv); }

    void eigenvalue(size_t valuev)
    { this->metaobjp_->eigenvalue(valuev); }

    void sideband(size_t valuev)
    { this->metaobjp_->sideband(valuev); }

  private:
    smartptr<SpectrumIterator<AmpType>,false> metaobjp_;
  };

  struct BaseFID_ {
    virtual void add_FID(BaseList<complex>& FID, complex scale) =0;
    virtual void larmor(double) =0;
    virtual ~BaseFID_() {}
  };
  
    template<class Source, class PropUse, class TypeSD> class RawMetaFID : public BaseMetaObj<Source,PropUse,TypeSD>, public BaseFID_ {
  public:

      using BaseMetaObj<Source,PropUse,TypeSD>::dirn_;
      using BaseMetaObj<Source,PropUse,TypeSD>::mziterp_;
      using BaseMetaObj<Source,PropUse,TypeSD>::inciter;
      using BaseMetaObj<Source,PropUse,TypeSD>::ismiddle_;
      using BaseMetaObj<Source,PropUse,TypeSD>::doublemzblock;
      using BaseMetaObj<Source,PropUse,TypeSD>::useblocks_;
      using BaseMetaObj<Source,PropUse,TypeSD>::eigoffset_;
      using BaseMetaObj<Source,PropUse,TypeSD>::opgen_base_;
      using BaseMetaObj<Source,PropUse,TypeSD>::verbose_;
      using BaseMetaObj<Source,PropUse,TypeSD>::shuffle;
      using BaseMetaObj<Source,PropUse,TypeSD>::sigma0_;
      using BaseMetaObj<Source,PropUse,TypeSD>::detect_;
      using BaseMetaObj<Source,PropUse,TypeSD>::objlist_;
      using BaseMetaObj<Source,PropUse,TypeSD>::blkspec_;

    template<class PropBase> RawMetaFID(const PropBase& objv, const Source& genv, int nobsv, double startv, double endv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(),int verbosev =0)
      : BaseMetaObj<Source,PropUse,TypeSD>(objv,genv,nobsv,startv,endv,sigma0v,detectv,verbosev) {}
     
	template<class PropBase> RawMetaFID(const PropBase& objv, const Source& Uv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(), int verbosev =0)
	  : BaseMetaObj<Source,PropUse,TypeSD>(objv,Uv,1,0.0,0.0,sigma0v,detectv, verbosev) {}

	  void larmor(double larmorv) { BaseMetaObj<Source,PropUse,TypeSD>::larmor(larmorv); }

    void add_FID(BaseList<complex>& FID, complex scale) {
      mziterp_->reset();
      size_t mzeigSD=0;
      bool makeherm;

      while (inciter(makeherm)) {
	const bool doubleblock(doublemzblock());
	const complex weightedscale(doubleblock ? 2.0*scale : scale);
	const bool isherm=makeherm || blkspec_.isherm;
	if (this->verbose_>1)
	  std::cout << "Setting sigma0/detect for mz index: " << mzeigSD << "  Double intensity: " << (doubleblock ? "Yes" : "No") << "  Hermitian: " << (isherm ? "Yes\n" : "No\n");
	for (size_t blkc=useblocks_;blkc--;) {
	  const size_t blk=blkc+eigoffset_;
	  const complex allweight(weightedscale*opgen_base_.eigweight(blk));
	  if (verbose_>1)
	    std::cout << "Eigenvalue: " << blk << "  Weight: " << allweight << '\n';
	  add_FID_(FID,allweight,mzeigSD,blkc,Bool2Type<SignalGenerator_traits<PropUse>::allowED>(),isherm);
	}
	mzeigSD++;
      }
      if (mzeigSD!=this->mzblocksSD_)
 	throw Mismatch("MetaFID: sigma0/detect blocks",mzeigSD,this->mzblocksSD_);
    }

  private:

      void add_FID_(BaseList<complex>& FID, complex scale, size_t mzeig, size_t blkc, Bool2Type<true>, bool isherm) {
	const size_t blk(blkc+eigoffset_);
	const cmatrix& sigma0(sigma0_(mzeig,blk));
	if ((sigma0.size()!=1) && this->isspecialcase(scale,mzeig,blk)) {
	  if (verbose_>1)
	    std::cout << "identity sigma0/detect special: adding " << scale << '\n';
	  objlist_(blkc).add_FID(FID,scale);
	  return;
	}
	if (verbose_>1)
	  std::cout << "sigma0\n" << sigma0 << '\n';
	
	if (!sigma0) //check for empty matrix (no check for matching detect, assume caught elsewhere)
	  return; 
	if (!detect_) {
	  if (isherm)
	    objlist_(blkc).add_FID_hermitian(FID,scale,sigma0);
	  else
	    objlist_(blkc).add_FID(FID,scale,sigma0);
	}
	else {
	  if (isherm)
	    objlist_(blkc).add_FID_hermitian(FID,scale,sigma0,detect_(mzeig,blk));
	  else
	    objlist_(blkc).add_FID(FID,scale,sigma0,detect_(mzeig,blk));
	}
      }

      void add_FID_(BaseList<complex>& FID, complex scale, size_t mzeig, size_t blkc, Bool2Type<false>, bool isherm) {
	const size_t blk(blkc+eigoffset_);
	if (!detect_)
	  throw Undefined("add_MetaFID: detect operator not specified");
	const cmatrix& sigma0(sigma0_(mzeig,blk));
	if (!sigma0) {//check for empty matrix (no check for matching detect, assume caught elsewhere)
	  std::cerr << "sigma0 undefined: ignoring\n";
	  return; 
	}
	if (this->isspecialcase(scale,mzeig,blk)) {
	  if (verbose_>1)
	    std::cout << "identity sigma0/detect special: adding " << scale << '\n';
	  objlist_(blkc).add_FID(FID,scale);
	  return;
	}
	if (verbose_>1) {
	  std::cout << "sigma0\n" << sigma0;
	  std::cout << "detect\n" << detect_(mzeig,blk) << '\n';
	}
	if (isherm)
	  objlist_(blkc).add_FID_hermitian(FID,scale,sigma0,detect_(mzeig,blk));
	else
	  objlist_(blkc).add_FID(FID,scale,sigma0,detect_(mzeig,blk));
      }
      
    };

    class MetaFID {
    public:
      template<class PropBase,class Source,class TypeSD> MetaFID(const PropBase& objv, const Source& genv, int nobsv, double startv, double endv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(),int verbosev =0)
	: metaobjp_(new RawMetaFID<Source,PropBase,TypeSD>(objv,genv,nobsv,startv,endv,sigma0v,detectv,verbosev)) {}

	template<class PropBase,class Source,class TypeSD> MetaFID(const PropBase& objv, const Source& Hv, double periodv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(),int verbosev =0)
	  : metaobjp_(new RawMetaFID<Source,PropBase,TypeSD>(objv,Hv,periodv,sigma0v,detectv,verbosev)) {}
	  
	  template<class PropBase,class Source,class TypeSD> MetaFID(const PropBase& objv, const Source& Hv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(),int verbosev =0)
	    : metaobjp_(new RawMetaFID<Source,PropBase,TypeSD>(objv,Hv,sigma0v,detectv,verbosev)) {}

      void add_FID(BaseList<complex>& FID, complex scale) { metaobjp_->add_FID(FID,scale); }
      void larmor(double larmorv) { metaobjp_->larmor(larmorv); }
    private:
      smartptr<BaseFID_,false> metaobjp_;
    };

  template<class PropUse, class Source, class TypeSD> void add_MetaFID(const PropUse& objv, BaseList<complex>& FID, complex scale, const Source& genv, int nobsv, double startv, double endv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(), int verbosev =0)
  {
    RawMetaFID<Source,PropUse,TypeSD> metaobj(objv,genv,nobsv,startv,endv,sigma0v,detectv,verbosev);
    metaobj.add_FID(FID,scale);
  }

  template<class PropUse, class TypeH, class TypeSD> void add_MetaFID(const PropUse& objv, BaseList<complex>& FID, complex scale, const TypeH& Hv, double periodv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(), int verbosev =0)
  {
    RawMetaFID<TypeH,PropUse,TypeSD> metaobj(objv,Hv,periodv,sigma0v,detectv,verbosev);
    metaobj.add_FID(FID,scale);
  }

  template<class PropUse, class TypeH, class TypeSD> void add_MetaFID(const PropUse& objv, BaseList<complex>& FID, complex scale, const TypeH& Hv, const TypeSD& sigma0v, const TypeSD& detectv =TypeSD(), int verbosev =0)
  {
    RawMetaFID<TypeH,PropUse,TypeSD> metaobj(objv,Hv,0.0,sigma0v,detectv,verbosev);
    metaobj.add_FID(FID,scale);
  }
 
  //Implementation details

template<class Source, class PropUse, class TypeSD> 
bool BaseMetaObj<Source,PropUse,TypeSD>::inciter(bool& makeherm)
{
  size_t rind,cind;
  bool ismiddle;
  
  if (!(mziterp_->next(rind,cind,mzeigSD_,ismiddle))) 
    return false;

  makeherm=blkspec_.isherm || (ismiddle && (dirn_!='M'));

  if (rind==cind)
    set_HUs(rind,'M');
  else {
    if (rind==lastcind_) {
      shuffle('C');
      set_HUs(cind,'C');
    }
    else {
      if (cind==lastrind_) {
	shuffle('R');
	set_HUs(rind,'R');
      }
      else {
	set_HUs(rind,'R');
	set_HUs(cind,'C');
      }
    }
  }
  if (this->verbose_>1)
    std::cout << "Row: " << lastrind_ << "  Col: " << lastcind_ << "  Hermitian: " << (makeherm ? "Yes\n" : "No\n");
  return true;
}

  template<class Source, class PropUse, class TypeSD>
    template<class PropBase> BaseMetaObj<Source,PropUse,TypeSD>::BaseMetaObj(const PropBase& objv, const Source& genv, int nobsv, double startv, double endv, const TypeSD& sigma0v, const TypeSD& detectv, int verbosev)
      : 	Hpgenp_(&genv),
	opgen_base_(genv.generator()),
	sigma0_(sigma0v), detect_(detectv),
    verbose_(verbosev),
	start_(startv), end_(endv),
		objlist_(opgen_base_.eigblocks(),objv,mxflag::normal) {
      
    common_init(nobsv,endv-startv,true);
  }

  template<class Source, class PropUse, class TypeSD> void
    RawMetaSpectrum<Source,PropUse,TypeSD>::set_sigma0dets(size_t mzeigSD)
  {
    if (this->verbose_>1)
      std::cout << "Setting sigma0/detect for mz index: " << mzeigSD << '\n';
    for (size_t blkc=useblocks_;blkc--;) {
      set_sigma0det(blkc,mzeigSD,Bool2Type<SignalGenerator_traits<PropUse>::allowED>(),Bool2Type<SignalGenerator_traits<PropUse>::allownonED>());
      if (verbose_>1)
	std::cout << "Object " << (blkc+eigoffset_) << '\n' << objlist_(blkc+eigoffset_);
    }
  }

  template<class Source, class PropUse, class TypeSD>
  void RawMetaSpectrum<Source,PropUse,TypeSD>::reset()
    {
      if (!blkspec_) { //No overlap between sigm0 and detect!
	finished_=true;
	return;
      }
    finished_=false;
    havenext_=false;
    mziterp_->reset();
    if (!next())
      throw Failed("RawMetaSpectrum: Failed to initialise");
    reset_blocks(useblocks_-1);
  }

  namespace {
    complex conj_(const complex& a) { return conj(a); }
    double conj_(double a) { return a; }
  };

  template<class Source, class PropUse, class TypeSD>
    bool RawMetaSpectrum<Source,PropUse,TypeSD>::operator() (typename PropUse::amplitude_type& amp, double& freq) {
    if (finished_)
      return false;

    if (havenext_) {
      amp=nextamp_;
      freq=nextfreq_;
      havenext_=false;
      return true;
    }

    while (!(*curobjp_)(amp,freq)) {
      if (blockpos_==0) {
	if (!next()) {
	  finished_=true;
	  return false;
	}
	reset_blocks(useblocks_-1);
      }
      else
	reset_blocks(blockpos_-1);
    }
    amp*=totweight_; //any weighting from mz symmetry
    if (makeherm_) { //if auto-hermitian, line up flipped amp,freq pair
      havenext_=true;
      nextamp_=conj_(amp);
      nextfreq_=-freq;
    }
    return true;
  }    

  template<class Source, class PropUse, class TypeSD>
  void RawMetaSpectrum<Source,PropUse,TypeSD>::set_sigma0det(size_t blkc, size_t mzeigSD, Bool2Type<true>, Bool2Type<false>)
  {
    if (!!detect_)
      throw InternalError("set_sigma0det<true,false> called with detect operator");
    const size_t blk(blkc+eigoffset_);    
    const cmatrix& subsigma0(sigma0_(mzeigSD,blk));
    if (subsigma0.size())
      objlist_(blkc).observe(subsigma0);
    else
      objlist_(blkc).reset();
  }

  template<class Source, class PropUse, class TypeSD>
  void RawMetaSpectrum<Source,PropUse,TypeSD>::set_sigma0det(size_t blkc, size_t mzeigSD, Bool2Type<true>, Bool2Type<true>)
  {
    const size_t blk(blkc+eigoffset_);
    complex dummy;
    if (verbose_ && this->isspecialcase(dummy,mzeigSD,blk))
      diagonalsigma0detectmissed_warning.raise();
    if (!detect_)
      set_sigma0det(blkc,mzeigSD,Bool2Type<true>(),Bool2Type<false>());
    else
      set_sigma0det(blkc,mzeigSD,Bool2Type<false>(),Bool2Type<true>());
  }

  template<class Source, class PropUse, class TypeSD>
  void RawMetaSpectrum<Source,PropUse,TypeSD>::set_sigma0det(size_t blkc, size_t mzeigSD, Bool2Type<false>, Bool2Type<true>)
  {
    const size_t blk(blkc+eigoffset_);
    const cmatrix& subsigma0(sigma0_(mzeigSD,blk));
    if (subsigma0.size())
      objlist_(blkc).observe(subsigma0,detect_(mzeigSD,blk));
    else
      objlist_(blkc).reset();
  }

  template<class Source, class PropUse, class TypeSD> void
  BaseMetaObj<Source,PropUse,TypeSD>::set_U(PropUse& obj, Int2Type<2>, Bool2Type<false>, Bool2Type<true>, char which) {
    if (which=='M')
      obj.set_U(passlist_.front(),period_);
    else
      obj.set_U(which,passlist_.front(),period_);
  }

  template<class Source, class PropUse, class TypeSD> void
  BaseMetaObj<Source,PropUse,TypeSD>::set_U(PropUse& obj, Int2Type<3>, Bool2Type<false>, Bool2Type<true>, char which) {
    if (which=='M')
      obj.set_Us(passlist_,period_);
    else
      obj.set_Us(which,passlist_,period_);
  }      

  template<class Source, class PropUse, class TypeSD> void
  BaseMetaObj<Source,PropUse,TypeSD>::set_U(PropUse& obj, Int2Type<2>, Bool2Type<true>, Bool2Type<true>, char which) {
    cmatrix& source(passlist_.front());
//     cmatrix Utmp;
//     cmatrix* usep=&source;
//     if (reduction_factor_>1) {
//       pow(Utmp,source,reduction_factor_);
//       usep=&Utmp;
//     }
    if (which=='M')
      obj.set_U(source);
    else
      obj.set_U(which,source);
  }

  template<class Source, class PropUse, class TypeSD> void
  BaseMetaObj<Source,PropUse,TypeSD>::set_U(PropUse& obj, Int2Type<3>, Bool2Type<true>, Bool2Type<true>, char which) {
    if (which=='M')
      obj.set_Us(passlist_);
    else
      obj.set_Us(which,passlist_);
  }      

  template<class Source,class PropUse,class TypeSD>
  template<bool isTD> void BaseMetaObj<Source,PropUse,TypeSD>::set_H_2nd_(PropUse& obj, Bool2Type<isTD> istd, size_t mzeig, size_t blk,char which,Bool2Type<false>) const
  {
    const BaseList<double> Hzeeman(opgen_base_.Hzeeman(mzeig,blk));
    if (opgen_base_.mzblocks()==1) { //if only 1 mz block then can ignore "twisting" of eigenbasis
      hermitian_eigenvalues_second(seigs_,(*Hpgenp_)(mzeig,blk),Hzeeman);
      set_H_(obj,istd,seigs_,which);
    }
    else {
      typename Source::suboutput_type V;
      hermitian_eigensystem_second(V,seigs_,(*Hpgenp_)(mzeig,blk),opgen_base_.Hzeeman_structure(mzeig,blk),Hzeeman);
      set_H_(obj,istd,V,seigs_,which);
    }
  }      

  template<class Source, class PropUse, class TypeSD>
    template<bool isTD,bool allowU> bool BaseMetaObj<Source,PropUse,TypeSD>::set_HU(PropUse& obj, Bool2Type<true>, Bool2Type<allowU>, Bool2Type<isTD> istd, size_t mzeig, size_t blk, char which) {
    switch (opgen_base_.effective_quadrupole_order()) {
    case 1:
      set_H_(obj,istd,(*Hpgenp_)(mzeig,blk),which);
      break;
    case 0: {
      typename Source::suboutput_type Huse((*Hpgenp_)(mzeig,blk));
      Huse+=opgen_base_.Hzeeman(mzeig,blk);
      set_H_(obj,istd,Huse,which);
      if (obslarmor_)
	obj.larmor(obslarmor_);
    }
      break;
    case 2:
      set_H_2nd_(obj,istd,mzeig,blk,which,Bool2Type<Ham_traits<Source>::isdiagonal>());
      break;
    default:
      throw Failed("BaseMetaObj: unhandled quadrupole order");
    }
    return true;
  }

  template<class Source, class PropUse, class TypeSD> void
  BaseMetaObj<Source,PropUse,TypeSD>::set_HUs(size_t mzeig, char which)
  {
    bool collapse=false;
    if (mzeig==mzblocks_) {
      switch (which) {
      case 'R':
	if (lastcind_==mzeig-1)
	  collapse=true;
	break;
      case 'C':
	if (lastrind_==mzeig-1)
	  collapse=true;
      }
      if (!collapse)
	throw InternalError("set_HUs");
    }
    switch (which) {
    case 'M':
      lastrind_=mzeig;
      lastcind_=mzeig;
      break;
    case 'R':
      lastrind_=mzeig;
      break;
    case 'C':
      lastcind_=mzeig;
      break;
    default:
      throw InvalidParameter("BaseMetaObj::set_HUs");
    }
    if (verbose_>1)
      std::cout << "mzeigH(" << which << "): " << mzeig << '\n';
    for (size_t blkc=useblocks_;blkc--;) {
      const size_t blk=eigoffset_+blkc;
      PropUse& obj_=objlist_(blkc);
      if (collapse) {
	if (verbose_>1)
	  std::cout << "k=" << blk << ": setting row=col\n";
	obj_.copy(which);
      }
      else {
	const size_t s(opgen_base_.size(mzeig,blk));
	if (verbose_>1)
	  std::cout << "k=" << blk << ": " << s << 'x' << s << '\n';
	if (s) { //non-zero block?
	  if (!set_HU(obj_ ,Bool2Type<SignalGenerator_traits<PropUse>::allowH && Propagator_traits<Source>::allowsetH >(), Bool2Type<SignalGenerator_traits<PropUse>::allowU>(), Bool2Type<SignalGenerator_traits<PropUse>::timedomain>(),mzeig,blk,which)) {
	    set_U(obj_, Int2Type<SignalGenerator_traits<PropUse>::input_dimensionality>(), Bool2Type<SignalGenerator_traits<PropUse>::timedomain>(), Bool2Type<SignalGenerator_traits<PropUse>::allowU>(),which);
	  }
	}
	else {
	  if (which=='M')
	    obj_.clear();
	  else
	    obj_.clear(which);
	}
      }
    }
  }

  template<class Source, class PropUse, class TypeSD>
  void
    BaseMetaObj<Source,PropUse,TypeSD>::initSD()
  {
    if (!sigma0_)
      throw Undefined("BaseMetaObj: initial matrix");
    if (opgen_base_.eigblocks()!=sigma0_.eigblocks())
      throw Mismatch("MetaObj: Hamiltonian and sigma0 have different eigenvalue structure");

    mzblocksSD_=sigma0_.actual_mzblocks();
    blkspec_=sigma0_.blockstructure();

    if (!!detect_) {
      if (!arematching(sigma0_,detect_))
	throw Mismatch("MetaObj: sigma0 / detect don't match");	  

      blkspec_&= detect_.blockstructure();
      if (!blkspec_)
	std::cerr << "Warning: sigma0 and detection operator do not overlap:\n" << sigma0_.blockstructure() << detect_.blockstructure();
    }

    if (verbose_>1)
      ::std::cout << "Block structure: " << blkspec_;

    if (blkspec_.isdiagonal())
      dirn_='M';
    else
      dirn_=blkspec_.isupper() ? 'R' : 'C';

    mziterp_.reset(new block_pattern::iterator(blkspec_));
  }

  template<class Source, class PropUse, class TypeSD> void
  BaseMetaObj<Source,PropUse,TypeSD>::common_init(int nobsv, double periodv, bool havevalidv)
  {
    lastrind_=lastcind_=-1;
    useblocks_=opgen_base_.eigblocks(); //by default use all eigenvalues
    eigoffset_=0;
    mzblocks_=opgen_base_.actual_mzblocks();
    usemzsym_=opgen_base_.usemzsym();
    //    reduction_factor_=1;
    mzeigSD_=0; //!< should be redundant, set by inciter
    obslarmor_=0.0;

    if (nobsv<=0 || periodv<0.0)
      throw InvalidParameter("MetaObj initialise");
    nobs_=nobsv;
    static const bool isinhomo= !SignalGenerator_traits<PropUse>::allowH && !SignalGenerator_traits<PropUse>::allowU;
    static const bool need_period=!isinhomo && ((!SignalGenerator_traits<PropUse>::allowH) ^ (SignalGenerator_traits<PropUse>::timedomain));
    if ((periodv==0.0) && need_period)
      throw InvalidParameter("MetaObj: period must be specified (and >0)");
  
    havevalidpgen_=havevalidv;
    period_=periodv;

    initSD();
   }

#endif
