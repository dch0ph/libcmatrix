/*! \file 
  \brief Chebyshev and Pade approximant propagation routines */

#undef LCM_SUPPRESS_VIEWS
#include "PartitionedMatrix.h"
#include "cmatrix.h"
#include "NMR.h"
#include "BlockedMatrix.h"

namespace libcmatrix {

  void hermitian_multiply(PartitionedMatrix<double>& d, const PartitionedMatrix<double>& a, PartitionedMatrix<double>& b) { 
    raw_multiply(d,a,b,false,true);
  }

  void hermitian_multiply(PartitionedMatrix<complex>& d, const PartitionedMatrix<complex>& a, PartitionedMatrix<complex>& b) { 
    raw_multiply(d,a,b,false,true);
  }

  void matrix_partition::swappresent(matrix_partition& a) 
  {
    if (!arematching(*this,a))
      throw Failed("swappresent");
    present_.swap(a.present_);
    std::swap(density_,a.density_);
  }

  const size_t eigensystem_controller::max_chebyshev_iterations=11;

  std::ostream& operator<< (std::ostream& ostr, const matrix_partition& a)
  {
    if (a.isdiagonal_) {
      ostr << "Offsets: " << a.roffs_ << '\n';
      //      if (a.totalblocks_)
      //	ostr << "Antisymmetric partition of " << a.totalblocks_ << " block matrix\n";
    }
    else {
      ostr << "Row offsets: " << a.roffs_ << '\n';
      ostr << "Column offsets: " << a.coffs_ << '\n';
    }
    spy(ostr,a.present_);
    ostr << "Density (whole block): " << 100.0*a.density() << "%\n";
    return ostr;
  }

  std::ostream& operator<< (std::ostream& ostr, const diagonal_partition& a)
  {
    return ostr << "Offsets: " << a.offs_  << "  Density (whole block): " << 100.0*a.density() << "%\n";
  }
  
  std::ostream& operator<< (std::ostream& ostr, const matrix_partition_set& a)
  {
    for (size_t i=0;i<a.size();i++) {
      if (a.size()>1)
	ostr << "Block " << i << '\n';
      ostr << a(i);
    }
    return ostr;
  }

  double diagonal_partition::density() const { return 1.0/rows(); }

  double diagonal_partition_set::density() const 
  {
    size_t totr=0;
    for (size_t i=size();i--;)
      totr+=(*this)(i).rows();
    return 1.0/totr;
  }
  
  double matrix_partition_set::density() const
  {
    double totused=0.0;
    size_t totblocks=0;
    for (size_t i=size();i--;) {
      const matrix_partition& part((*this)(i));      
      const size_t blocks=part.rows()*part.cols();
      totused+=size_t(0.5+part.density()*blocks);
      totblocks+=blocks;
    }
    return totused/totblocks;      
  }

  double matrix_partition::density() const
  {
    if (density_<0.0) {
      size_t nonzero=0;
      for (size_t r=present_.rows();r--;) {
	const size_t nr=rows(r);
	for (size_t c=present_.cols();c--;) {
	  if (present_(r,c))
	    nonzero+=nr*cols(c);
	}
      }
      density_=double(nonzero)/(rows()*cols());
    }
    return density_;
  }

  void diagonal_partition::makepartition(ScratchList<size_t>& offs, const BaseList<size_t>& sizes)
  {
    const size_t blks=sizes.size();
    offs.create(blks+1);
    size_t sofar=0;
    for (size_t i=0;i<blks;i++) {
      offs(i)=sofar;
      sofar+=sizes(i);
    }
    offs.back()=sofar;
  }
 
  diagonal_partition::diagonal_partition(const BaseList<size_t>& sizes)
   {
     makepartition(offs_,sizes);
   }

  void matrix_partition_base::create(const BaseList<size_t>& sizes)//, size_t totalblocksv)
  {
    diagonal_partition::makepartition(roffs_,sizes);
    coffs_=roffs_;
    isdiagonal_=true;
  }
    
  void matrix_partition_base::create(const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes)
  {
    isdiagonal_=false;
    diagonal_partition::makepartition(roffs_,rsizes);
    diagonal_partition::makepartition(coffs_,csizes);
  }

  void matrix_partition::create(const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes)
  {
    density_=-1.0;
    matrix_partition_base::create(rsizes,csizes);
    present_.create(rsizes.size(),csizes.size(),false);
  }

  void matrix_partition::create(const BaseList<size_t>& sizes)//, size_t totalblocksv)
  {
    density_=-1.0;
    matrix_partition_base::create(sizes);//,totalblocksv);
    present_.create(sizes.size(),sizes.size(),false);
  }

  void matrix_partition::adddiagonal()
  {
    for (size_t r=present_.rows();r--;) {
      if (rows(r))
	present_(r,r)=1;
    }
  }

  void matrix_partition_set::adddiagonal()
  {
    for (size_t i=size();i--;)
      (*this)(i).adddiagonal();
  }

  void matrix_partition::cleanempty()
  {
    for (size_t r=present_.rows();r--;) {
      if (rows(r)==0) {
	for (size_t c=present_.cols();c--;)
	  present_(r,c)=false;
      }
    }
    for (size_t c=present_.cols();c--;) {
      if (cols(c)==0) {
	for (size_t r=present_.rows();r--;)
	  present_(r,c)=false;
      }
    }
  }    
      
  matrix_partition::matrix_partition(const Matrix<bool>& present, const BaseList<size_t>& rsizes, const BaseList<size_t>& csizes)
  {
    create(rsizes,csizes);
    if ((present.rows()!=present_.rows()) || (present.cols()!=present_.cols()))
      throw Mismatch("matrix_partition");
    present_=present;
    cleanempty();
  }

  matrix_partition::matrix_partition(const Matrix<bool>& present, const BaseList<size_t>& sizes)
  {
    if (!issquare(present))
      throw NotSquare("matrix_partition");
    create(sizes);
    if (present.rows()!=present_.rows())
      throw Mismatch("matrix_partition");
    present_=present;
    cleanempty();
  }

  namespace {
    inline std::pair<double,double> eigenvalue_range_(const Matrix<double>& a) { return eigenvalue_range(a); }
    inline std::pair<double,double> eigenvalue_range_(const cmatrix& a) { return hermitian_eigenvalue_range(a); }
  }

  matrix_partition& matrix_partition::operator= (const matrix_partition_base& a) 
  { 
    matrix_partition_base::operator= (a);
    density_=0.0;
    present_.create(a.row_size(),a.col_size(),false);
    return *this;
  }
  
  matrix_partition& matrix_partition::operator|= (const matrix_partition& a)
  {
    if (!*this)
      return (*this=a);
    apply_ip2(doesor_ip<bool,bool>(),present_,a.present_);
    return *this;
  }

  namespace {
    template<class T> inline const Matrix<T>& asmatrix(const Matrix<T>& a) { return a; }
    template<class T> inline Matrix<T>& asmatrix(Matrix<T>& a) { return a; }
    template<class T> inline const Matrix<T>& asmatrix(const PartitionedMatrix<T>& a) { return a.asmatrix(); }
    template<class T> inline Matrix<T>& asmatrix(PartitionedMatrix<T>& a) { return a.asmatrix(); }

    template<class T> void add_diag_ip(Matrix<T>& U, double v)
    {
      for (size_t i=U.rows();i--;)
	U(i,i)+=v;
    }

    template<class T> void add_diag_ip(PartitionedMatrix<T>& U, double v)
    {
      add_diag_ip(U.asmatrix(),v);
      for (size_t i=U.partitioning().row_size();i--;)
	U.partitioning().present(i,i)=true;
    }

    template<class T,class PartType> void validate_partitioning(const Matrix<T>& H, const PartType& part) 
    {
      if (cmatrix_eigensystem_controller.verbose) {
	if (!part.validate(H,cmatrix_eigensystem_controller.tolerance)) {
	  std::cerr << part;
	  spy(std::cerr,H,cmatrix_eigensystem_controller.tolerance);
	  std::cerr.flush();
	  throw Failed("propagator: partitioning failure\n");
	}
	std::cout << "propagator: matrix partitioning validated\n";
      }
    }

    void expand_partitioning(cmatrix& d, const cmatrix& a, const matrix_partition& part)
    {
      d=complex(0.0); 
      const range sel(0,part.rows()-1);
      d(sel,sel)=a;
      size_t start=d.rows();
      size_t curind=0;
      for (size_t i=part.rows();i--;) {
	if (start<a.rows()) {
	  if (start==a.rows()-1)
	    return;
	  throw Mismatch("expand_partitioning");
	}
	const size_t nels=part.rows(curind);
	const range sel(part.row_range(curind++));
	const range destsel(start-nels,start-1);
	d(destsel,destsel)=a(sel,sel);
	start-=nels;
      }
    }

    const double cheby_s=0.5;
    const double cheby_a0=0.938469807240831;
    const ScratchList<complex,11> cheby_as(complex(0,0.484536915349748),
					    complex(-0.0612080469173653,0),
					    complex(0,-0.00512745998917449),
					    complex(0.000321472952728575,0),
					    complex(0,1.6107254482715e-5),
					    complex(-6.7213692572377e-007,0),
					    complex(0,-2.40317346555261e-8),
					    complex(7.51644630959522e-10,0),
					    complex(0,2.0893535178658e-11),
					    complex(-5.2263547216456e-13,0),
					    complex(0,2.47676511895986e-16));
  }

  template<typename T> struct lessfunc : public ::std::binary_function<T,T,bool> {
    inline bool operator()(const T& a, const T& b) const { return (a<b); }
  };
  template<> struct lessfunc<complex> : public ::std::binary_function<complex,complex,bool> {
    inline bool operator()(const complex& a, const complex& b) const { return (norm(a)<norm(b)); }
  };
    
  template<class M> void raw_chebyshev_propagator(cmatrix& U, M& Hnorm, double dt)
  {
    M Hprev,Htmp;
    M Hcur(Hnorm);
    Hnorm*=2.0; //!< avoids factor of 2 multiply later
    const M& Hnorm2(Hnorm);

    const size_t maxn=cmatrix_eigensystem_controller.chebyshev_iterations;     
    for (size_t k=1;;k++) {
      const complex an(cheby_as(k-1));
      if (cmatrix_eigensystem_controller.verbose) {
	std::cout << "Iteration " << k << ": adding " << an << " of: max=" << max(asmatrix(Hcur).row(),lessfunc<LCM_VAL(M) >()) << std::endl;
	if (cmatrix_eigensystem_controller.verbose>1)
	  std::cout << Hcur << '\n';
      }
      if (k>1) {
	if (imag(an)==0.0)
	  mla(U,real(an),asmatrix(Hcur));
 	else
 	  mla(U,an,asmatrix(Hcur));
      }
      else {
	U=asmatrix(Hcur);
	U*=an;
	add_diag_ip(U,cheby_a0);
      }
      if (k==maxn)
	break;
      hermitian_multiply(Htmp,Hnorm2,Hcur);
      if (!Hprev)
	add_diag_ip(Htmp,-1.0);
      else
	Htmp-=Hprev;
      Hprev.swap(Hcur);
      Hcur.swap(Htmp);
    }
  }

  Warning<> chebyshev_convergence_warning("chebyshev_propagator: convergence criterion not met.  Try decreasing integration timestep by a factor of ",&lcm_base_warning);
  Warning<> chebyshev_nonunitary_warning("non-unitary propagator / eigenvectors detected",&lcm_base_warning);

  void docheck_(Matrix<complex>& dest, const Matrix<complex>& U)
  {
    conj_transpose_multiply(dest,U,U);
  }

  void docheck_(Matrix<double>& dest, const Matrix<double>& U)
  {
    transpose_multiply(dest,U,U);
  }

  template<typename T> bool isunitary_(const Matrix<T>& U, const char* name)
  {
    if (cmatrix_eigensystem_controller.tolerance==0)
      return true;
    Matrix<T> checkm;
    docheck_(checkm,U);
    const bool failed=hasoffdiagonal(checkm,cmatrix_eigensystem_controller.tolerance);
    if (failed) {
      if (!(cmatrix_eigensystem_controller.throwexception))
	chebyshev_nonunitary_warning.raise();
      if (cmatrix_eigensystem_controller.verbose>1)
	std::cerr << name << "'*" << name << ":\n" << checkm;
      if (cmatrix_eigensystem_controller.throwexception)
	chebyshev_nonunitary_warning.raiseas(BaseWarning::RaiseException);
    }
    else {
      switch (cmatrix_eigensystem_controller.verbose) {
      case 0: break;
      case 1:
	std::cout << "Passed orthogonality check\n";
	break;
      default:
	std::cout << name << "'*" << name << ":\n" << checkm;
      }
    }
    return !failed;
  }      

  bool isunitary(const cmatrix& U, const char* name)
  {
    return isunitary_(U,name);
  }

  bool isunitary(const rmatrix& U, const char* name)
  {
    return isunitary_(U,name);
  }

  char stabfacchar[10];

  template<class T> void chebyshev_propagator_(cmatrix& U, const Matrix<T>& H, double dt, const matrix_partition* partp)//, bool acc =true)
  {
    dt*=2*M_PI;
    const size_t maxn=cmatrix_eigensystem_controller.chebyshev_iterations;
    if ((maxn>cheby_as.size()) || (maxn<1))
      throw InvalidParameter("chebyshev_propagator: iterations");
    if (dt<=0.0)
      throw InvalidParameter("chebyshev_propagator: dt");
    if (cmatrix_eigensystem_controller.tolerance) {
      const std::pair<double,double> erange(eigenvalue_range_(H));
      const double diff=std::max(fabs(erange.second),fabs(erange.first));
      const double stabfac=diff*dt;
      if (cmatrix_eigensystem_controller.verbose) {
	std::cout << "Estimated maximum eigenvalue range: " << erange.first << " to " << erange.second << '\n';
	std::cout << "Stability factor: " << stabfac << '\n';
      }
      if (stabfac>1.0) {
	snprintf(stabfacchar,sizeof(stabfacchar)-1,"%.2g\n",stabfac);
	if (cmatrix_eigensystem_controller.throwexception)
	  chebyshev_convergence_warning.raiseas(BaseWarning::RaiseException,stabfacchar);
	else
	  chebyshev_convergence_warning.raise(stabfacchar);
      }
    }
    const double scale=dt/cheby_s;
    if (H.empty())
      throw Undefined("chebyshev_propagator");

    if (partp) {
      validate_partitioning(H,*partp);
      PartitionedMatrix<T> Hnorm(H,*partp,mxflag::normal);

      Hnorm*=-scale;
      raw_chebyshev_propagator(U,Hnorm,dt);
      if (H.rows()>U.rows()) {
	cmatrix tmp(H.rows(),H.cols());
	expand_partitioning(tmp,U,*partp);
	U.swap(tmp);
      }
    }
    else {
      Matrix<T> Hnorm(H,mxflag::normal);
      Hnorm*=-scale;
      raw_chebyshev_propagator(U,Hnorm,dt);
    }      
    if (!isunitary(U) && (cmatrix_eigensystem_controller.verbose>1) )
      std::cerr << "Input matrix\n" << H;
  }
  
  void chebyshev_propagator(cmatrix& U, const Matrix<double>& H, double dt, const matrix_partition* partp) { chebyshev_propagator_(U,H,dt,partp); }

  void chebyshev_propagator(cmatrix& U, const Matrix<complex>& H, double dt, const matrix_partition* partp) { chebyshev_propagator_(U,H,dt,partp); }

  namespace {
    inline double lreal(double x) { return x; }
    inline double lreal(const complex& z) { return z.real(); }
  }

  template<class T> void propagator_(cmatrix& U, const Matrix<T>& H, double dt, const diagonal_partition* partp)
  {
    if (!partp) {
      propagator(U,H,dt);
      return;
    }
    if (!issquare(H))
      throw NotSquare("propagator (partitioned)");

    validate_partitioning(H,*partp);
    const size_t trows=H.rows();
    U.create(trows,trows,complex(0.0));
    Matrix<T> Htmp;
    cmatrix Utmp;
    for (size_t r=partp->rows();r--;) {
      const range sel(partp->index_range(r));
      if (sel.size()==1) {
	const size_t ind=sel.start();
	U(ind,ind)=propagator(lreal(H(ind,ind)),dt);
      }
      else {
	Htmp=RawMatrix<T>(H,sel);
	propagator(Utmp,Htmp,dt);
	RawMatrix<complex> Usel(U,sel);
	Usel=Utmp;
      }
    }
  }

  void propagator(cmatrix& U, const Matrix<double>& H, double dt, const diagonal_partition* partp) { propagator_(U,H,dt,partp); }

  void propagator(cmatrix& U, const Matrix<complex>& H, double dt, const diagonal_partition* partp) { propagator_(U,H,dt,partp); }

  template<class T> void propagator_(BlockedMatrix<complex>& U, const BlockedMatrix<T>& H, double dt, const diagonal_partition_set& part)
  {
    if (part.empty()) {
      if (U.empty())
	propagator(U,H,dt);
      else {
	BlockedMatrix<complex> Utmp;
	propagator(Utmp,H,dt);
	U&=Utmp;
      }
      return;
    }
    if (part.size()!=H.size())
      throw Mismatch("propagator: Hamiltonian and partitioning",H.size(),part.size());
    U.duplicate_structure(H);
    for (size_t i=H.size();i--;)
      propagator(U(i),H(i),dt,&(part(i)));
  }

  void propagator(BlockedMatrix<complex>& U, const BlockedMatrix<double>& H, double dt, const diagonal_partition_set& part) { propagator_(U,H,dt,part); }

  void propagator(BlockedMatrix<complex>& U, const BlockedMatrix<complex>& H, double dt, const diagonal_partition_set& part) { propagator_(U,H,dt,part); }

  //Pade

  template<class T> inline double abs_(const T& x) { return abs(x); }
  inline double abs_(double x) { return std::fabs(x); }

  template<class T> double infinity_norm(const Matrix<T>& A)
  {
    double maxval=0.0;  
    for (size_t i=A.rows();i--;) {
      double rowsum=0.0;
    for (size_t j=A.cols();j--;)
      rowsum+=abs_(A(i,j));
    if (rowsum>maxval)
      maxval=rowsum;
    }
    return maxval;
  }
  
  inline double ln2(double x)
  {
    static const double iln2(1.0/std::log(2.0));
    return std::log(x)*iln2;
  }
  
  template<class T1,class T2> inline void multiply_ip_tmp(T1& A, const T2& B, T1& tmp)
  {
    multiply(tmp,A,B);
    A.swap(tmp);
  }
  
  template<class Mout, class Min> void exppade_ip_(Mout& E, Min& A, size_t order)
  {
    if (!issquare(A))
      throw NotSquare("exppade_ip");
    
    static List<double> padecoeffs;
    if (order<3)
      throw InvalidParameter("exppade_ip: order must be >2");
    if (padecoeffs.size()!=order+1) {    
      padecoeffs.create(order+1);
      padecoeffs(0U)=1.0;
      for (size_t k=1;k<=order;k++) {
	const int ik(k);
	padecoeffs(k)=padecoeffs(k-1)*(order+1-ik)/(ik*(2*order+1-ik));
      }
    }
    //  std::cout << "Pade coefficients: " << padecoeffs << '\n';
    const double s=infinity_norm(asmatrix(A));
    if (s==0.0)
      throw Failed("exppade_ip: null input matrix");
    int is((int)s);
    if (s>0.5) {
      is=int(ln2(s)+2);
      if (is>0)
	A*=(1.0/(1<<is)); 
      else
	is=0;
    }
    //  std::cout << "Infinity norm: " << s << "  Squares: " << is << '\n';
    
    Min A2,tmp,Q,P;
    multiply(A2,A,A);
    
    bool odd=true;
    for (size_t k=order-1;k--;) {    
      Min& use(odd ? Q : P);
      if (!use) {
	const double coeff=padecoeffs( (&use==&Q) ? order : order-1);
	use=A2;
	use*=coeff;
      }
      else 
	multiply_ip_tmp(use,A2,tmp);
      add_diag_ip(use,padecoeffs(k));
      odd=!odd;
    }
    
    Min& use(odd ? Q : P);
    multiply_ip_tmp(use,A,tmp);
    Q-=P;

    E=asmatrix(P);
    Mout& Qout(asmatrix(Q));
    multiply_inv_ip2(Qout,E); //!< N.B. *both* input matrices are overwritten
    E*=2.0;
    add_diag_ip(E,1.0);
    if (odd)
      negate_ip(E);
    for (size_t k=is;k--;)
      multiply_ip_tmp(E,E,asmatrix(tmp));
  }
    
  void exppade_ip(Matrix<complex>& A, size_t iter) { exppade_ip_(A,A,iter); }
  void exppade_ip(Matrix<double>& A, size_t iter) { exppade_ip_(A,A,iter); }

  size_t lcm_pade_iterations=LCM_DEFAULT_PADE_ITERATIONS;

  void pade_propagatorL(cmatrix& U, const cmatrix& L, double dt, const matrix_partition* partp)
  {
    if (dt<=0.0)
      throw InvalidParameter("pade_propagator: dt");
    //    dt*=2*M_PI;
    if (L.empty())
      throw Undefined("pade_propagator");

    if (partp) {
      validate_partitioning(L,*partp);
      PartitionedMatrix<complex> Lscaled(L,*partp,mxflag::normal);
      Lscaled*=dt;
      exppade_ip_(U,Lscaled,lcm_pade_iterations);
//       if (L.rows()>Lscaled.rows()) {
// 	U.create(L.rows(),L.cols());
// 	expand_partitioning(U,Lscaled.asmatrix(),*partp);
//       }
//       else
// 	U.swap(Lscaled.asmatrix()); //!< Lscaled will be an inconsistent state (but who cares?)
    }
    else {
      multiply(U,dt,L);
      exppade_ip(U,lcm_pade_iterations);
    }      
    isunitary(U);
  }

}


//   void matrix_partition::reduce(size_t& blk, size_t& offset, size_t ind, const BaseList<size_t>& offs)
//   {
//     bool bad=true;
//     for (blk=offs.size();blk--;) {
//       if (ind>=offs(blk)) {
// 	if (bad)
// 	  throw BadIndex("matrix_partition::reduce",ind,offs.back());	  
// 	offset=ind-offs(blk);
// 	return;
//       }
//       bad=false;
//     }	
//     throw InternalError("matrix_partition::reduce");
//   }

//   class row_iterator_base {
//   public:
//     row_iterator_base(const PartitionedMatrix<T>&, size_t, size_t);
//     bool operator==(const row_iterator_base& x) const { return (curcol==x.curcol); }
//     bool operator!=(const row_iterator_base& x) const { return (curcol!=x.curcol); }
//   protected:
//     const BlockedMatrix<T>& data;
//     const BaseList<size_t> coloffsets;
//     const size_t offsetr;
//     size_t curcol;
//     size_t curblkc;
//     void up();
//   };
  
//   struct const_row_iterator : public row_iterator_base, ::std::iterator< ::std::forward_iterator_tag,T,ptrdiff_t,const T*,const T&> {
//     const_row_iterator(const PartitionedMatrix<T>& a, size_t r, size_t c)
//       : row_iterator_base(a,r,c) {}
//     const T& operator*() const { return crow[cindex_(c_)]; }
//     const_iterator& operator++() { this->up(); return *this; }
//     const_iterator operator++(int) { const_iterator tmp(*this); this->up(); return tmp; }
//   };
  
