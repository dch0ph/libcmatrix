//internal file only

template<bool> struct lcm_cache_type_ { typedef Type2Type< BlockedMatrix<double> > type; };
template<> struct lcm_cache_type_<true> { typedef Type2Type< ListList<double> > type; };

class SpinOpGenerator : public SpinOpGeneratorBase {
  public:
  SpinOpGenerator(const basespin_system&, int flagsv =0, int verbosev =0);
  SpinOpGenerator(const basespin_system&, const CrystalStructure&, int flagsv =0, int verbosev =0);
  SpinOpGenerator(const basespin_system&, const BaseList<nuclei_spec>&, int flagsv =0, int verbosev =0);
  SpinOpGenerator(const basespin_system&, const CrystalStructure&, const BaseList<nuclei_spec>&, int flagsv =0, int verbosev =0);
  SpinOpGenerator(const basespin_system&, const CrystalStructure&, const char* label, int flagsv =0, int verbosev =0);
  SpinOpGenerator(const HamiltonianStructure&, const CrystalStructure&, int verbosev =0);

  SpinOpGeneratorBase* clone() const { return new SpinOpGenerator(*this); }

  typedef BlockedMatrix<complex> output_type;
  typedef Matrix<complex> suboutput_type;
  typedef double base_type; //operators *can* be real
  static const bool allowproductop=true;

  void proton_frequency(double);

    const spin& operator()(size_t n) const {
      return (*sysp_)(n);
    }

    size_t size(size_t mzblk, size_t eigblk) const;
    size_t eigweight(size_t) const { return 1; }
  const BaseList<size_t>& diagonal_structure() const { return sizes_; }

    bool ishomonuclear(size_t j,size_t k) const
    { return ((*sysp_)(j)==(*sysp_)(k)); }
    
  void add(BlockedMatrix<complex>&, block_pattern&, const productoperator_spec&) const;
  void mla(BlockedMatrix<complex>&, block_pattern&, double, const operator_spec&) const;
  void mla_Iz(ListList<double>&, double scale, size_t sel) const;
    void mla_Fz(ListList<double>&,double,nuclei_spec) const;
  const ListList<double>& diag_Fz(nuclei_spec whichn) const { return fzcache_(*this,whichn()); }
  void rawdiag_Fz(ListList<double>&,size_t) const;

  template<class T> void add_A2_homo(T& dest, const typename Ham_traits<T>::coupling_type& coup, size_t j, size_t k,Bool2Type<false>) const {
    add_A2(dest,&spin_dipolar,coup,j,k,Type2Type< BlockedMatrix<double> >());
  }
  
  template<class T> void add_A2_homo(T&, const typename Ham_traits<T>::coupling_type&, size_t, size_t,Bool2Type<true>) const {
    throw Failed("add_A2: can't add homonuclear to diagonal Hamiltonian");
  }

  template<class T> void add_A2(T& dest, const typename Ham_traits<T>::coupling_type& coup, size_t j, size_t k, bool forceweak =false) const {
    if (forceweak || !ishomonuclear(j,k))
      add_A2(dest,&diag_spin_truncated_dipolar,coup,j,k,Type2Type< ListList<double> >());
    else
      add_A2_homo(dest,coup,j,k,Bool2Type<Ham_traits<T>::isdiagonal>());
  }

  template<class T> void add_A0_strong(T& dest, double coup, size_t j, size_t k,Bool2Type<false>) const {
      add_A0(dest,&spin_J,coup,j,k,Type2Type< BlockedMatrix<double> >());
  }

  template<class T> void add_A0_strong(T&, double, size_t, size_t,Bool2Type<true>) const {
    throw Failed("add_A0: can't add strong coupling to diagonal Hamiltonian");
  }

  template<class T> void add_A0(T& dest, double coup, size_t j, size_t k ,bool forceweak =false) const {
    if (forceweak || !ishomonuclear(j,k))
      add_A0(dest,&diag_spin_weakJ,coup,j,k,Type2Type< ListList<double> >());
    else
      add_A0_strong(dest,coup,j,k,Bool2Type<Ham_traits<T>::isdiagonal>());
  }

  template<class T> void add_ns_(T&, const space_T&, size_t, size_t, ns_flag, Bool2Type<true>) const { throw Failed("Diagonal Hamiltonian incompatible with non-secular interactions"); }

  template<class T> void add_ns_(T&, const space_T&, size_t, size_t, ns_flag, Bool2Type<false>) const;

  template<class T> void add_ns(T& dest, const space_T& coup, size_t i, size_t j, ns_flag nstype) const
  { add_ns_(dest,coup,i,j,nstype,Bool2Type<Ham_traits<T>::isdiagonal>()); }

  template<class T> void add_ns(T&, double, size_t, size_t, ns_flag) const { throw Failed("Non-secular Hamiltonian cannot be initialised from scalar interactions"); }

  template<class T> void add_Hquadrupolar2(T&, const space_T&, size_t) const;
  template<class T> void add_Hquadrupolar2(T&, double, size_t) const { throw Failed("Can't add second order quadrupole to this Hamiltonian"); }
  void add_Hquadrupolar2(List<DiagonalSpinningHamiltonian>&, const space_T&, size_t) const;

  template<class T> void add_Hquadrupolar(T& dest, const typename Ham_traits<T>::coupling_type& coup, size_t j) const;
  template<class T> void add_Hcs(T& dest, const typename Ham_traits<T>::coupling_type&, size_t j) const;
        
    void print(std::ostream&) const;

    template<class T> void create(BlockedMatrix<T>& dest, const block_pattern&) const;

  template<class M> void create(M& dest) const
  { SpinOpGeneratorBase::create(dest); }

    void create_ns(List< MultiMatrix<double,3> >& dest) const {
      dest.create(mzblocks()); //later functions will create contents
    }

    const basespin_system& spinsystem() const { return *sysp_; }

  size_t index(size_t mzblk,size_t eigblk) const {
    if (eigblk)
      throw Failed("SpinOpGenerator: no eigenvalue structure");
    return mzblk;
  }

  //void coherence_filter(BlockedMatrix<bool>&, block_pattern&, const spinorder_spec&) const;
  //void coherence_filter_ip(BlockedMatrix<bool>&, block_pattern&, const spinorder_spec&, state_t) const;
  //  void mla_Fz_compressed(ListList<double>&, double, state_t) const;

  private:
  const smartptr<const basespin_system> sysp_;

  mutable ScratchList<complex,5> tmp_;
  mutable ScratchList<complex,5> coeffs_;

  void addQ2_(List<DiagonalSpinningHamiltonian>&, size_t, double, const space_T&, const ListList<double>&) const;
  
  typedef UnionHolder< 3, ListList<double>, BlockedMatrix<double>, List< MultiMatrix<double,3> > > cache_type;
  enum { BLOCKEDLIST=1, BLOCKEDMATRIX=2, BLOCKEDTENSOR=3 };

  mutable FzCache_ fzcache_;

    //These caches are protected by mutexes
    mutable Mutex<ThreadingActive> lock;
    mutable Matrix<cache_type> A0store,A2store;
    mutable List< ListList<double> > CSstore;

  mutable Mutex<ThreadingActive> iterlock;
  mutable smartptr<CrystalStructure_iterator> cstructiterp;

  void common_init(const BaseList<nuclei_spec>& =BaseList<nuclei_spec>());
  //  void makepartitioning();
    block_pattern getdiagtype(size_t nuc, char op) const;

  void mla_(cmatrix&,double,const operator_spec&, const operator_spec&) const;
  void mla_(cmatrix&,double,const operator_spec&) const;

  template<class T1,class T2> void chop_mla(BlockedMatrix<T1>&, double, const Matrix<T2>&, const block_pattern&) const;
  template<class T1,class T2> void chop_mla(BlockedMatrix<T1>& dest, double scale, const Matrix<T2>& source) const { chop_mla(dest,scale,source,block_pattern(active_blocking_info())); }
  template<class T1,class T2> void chop_mla(Matrix<T1>&, double, const Matrix<T2>&) const;// const block_pattern&) const;
    template<class T> void chop_mla(BlockedMatrix<T>&, double, const List<double>&) const;
    template<class T> void chop_mla(Matrix<T>&, double, const List<double>&) const;
    void chop_mla(ListList<double>&, double, const List<double>&) const;
    template<class T,class C> void chop_mla(List<T>&, const C&, const List<double>&) const;
    template<class C> void chop_mla(List<double>&, const C&, const List<double>&) const;
    template<class C,class M> void chop_mla(RealSpinningHamiltonian&, const C&, const M&) const;
    template<class C,class M> void chop_mla(SpinningHamiltonian&, const C&, const M&) const;
    template<class C,class M> void chop_mla(DiagonalSpinningHamiltonian&, const C&, const M&) const;
  template<class M,class C> void chop_mla(List<M>&, const C&, const rmatrix&) const;

    void chop_ns(List< MultiMatrix<double,3> >&, const MultiMatrix<double,3>&) const;

  template<class T> void add_ns_(T&, const space_T&, const BaseList< MultiMatrix<double,3> >& ) const { throw Failed("Can't add non-secular interaction to this Hamiltonian"); }

  void add_ns_(List<SpinningHamiltonian>& dest, const space_T& coup, const BaseList< MultiMatrix<double,3> >& source) const;
  void add_ns_(SpinningHamiltonian& dest, const space_T& coup, const BaseList< MultiMatrix<double,3> >& source) const;
  void add_ns_(BlockedMatrix<complex>& dest, const space_T& coup, const BaseList< MultiMatrix<double,3> >& source) const;
  void add_ns_(Matrix<complex>& dest, const space_T& coup, const BaseList< MultiMatrix<double,3> >& source) const;

    template<class T,class C,class W> static void add_(List<T>& dest, const C& coup, const W& a) {
      for (size_t k=a.size();k--;)
	dest(k).add(coup,a(k));
    }

    static void add_(List<double>& dest, double coup, const ListList<double>& a) {
      ::libcmatrix::mla(dest,coup,a.row());
    }
    template<class T,class W> static void add_(BlockedMatrix<T>& dest, double coup, const W& a)
    { ::libcmatrix::mla(dest,coup,a); }

    static void add_(ListList<double>& dest, double coup, const ListList<double>& a)
    { ::libcmatrix::mla(dest,coup,a); }

    template<class T,class M> void add_(Matrix<T>&, double, const M&) const;
    
  template<class T,class F,class C,class ResultType> void add_A2(T&,F func, const C&, size_t,size_t,Type2Type<ResultType>) const;
  template<class T,class F,class ResultType> void add_A0(T&,F func, double, size_t,size_t,Type2Type<ResultType>) const;

  template<class T,class F,class C,class ResultType> void addraw(T&, F func, const C&, size_t,size_t,Type2Type<ResultType>) const;
    template<class T,class C> void addCS(T&, const C&, size_t) const;
    template<class T,class C> void addrawCS(T&, const C&, size_t) const;
    template<class T,class C> void addrawQ(T&, const C&, size_t) const;
 
  void makeraw_ns(List< MultiMatrix<double,3> >& dest, size_t, size_t, ns_flag) const;
  
    static void spydefined(std::ostream&, const char*, const Matrix<cache_type>&);
    static void spydefined(std::ostream&, const char*, const List< ListList<double> >&);

    void dump() const;    
  void update();
  
  void makemzinds(Matrix<size_t>&, List<size_t>&, const BaseList<size_t>&);

  };

  template<> struct spinsystem_glue_<SpinOpGenerator> {
    static const basespin_system& spinsystem(const SpinOpGenerator& a)
    { return a.spinsystem(); }

    static size_t ncells(const SpinOpGenerator& a)
    { return a.ncells(); }

    static const bool ispurereal=true;
    
    template<class M> static bool arematching(const SpinOpGenerator& a, const HamiltonianStore<M>& store) {
      return (store.ncells()==a.ncells()) && (a.nspins()==(store.ncells()*store.nspins()));
    }
  };
  
  inline std::ostream& operator<< (std::ostream& ostr, const SpinOpGenerator& a) {
    a.print(ostr);
    return ostr;
  }
     
