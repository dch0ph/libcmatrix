//internal file BlockedOperator/Filter

class BlockedFilter;

   class BlockedOperator : public BaseStructure {
    public:
     BlockedOperator() : opgenp_(NULL) {}

     template<class OpGen> BlockedOperator(const OpGen& opgen, const operator_spec& spec, double scale =1.0)
     //       : Hindexer_(opgen.structure().indexer())
       : opgenp_(&opgen)
     {
       	opgen.mla(store_,blkspec_,scale,spec);
	static_cast< BaseStructure& >(*this)= BaseStructure(blkspec_.blocks,opgen.eigblocks());
      }

     template<class OpGen> BlockedOperator(const OpGen& opgen, const productoperator_spec& spec)
     //       : Hindexer_(opgen.structure().indexer())
       : opgenp_(&opgen)
     {
       opgen.add(store_,blkspec_,spec);
       static_cast< BaseStructure& >(*this)= BaseStructure(blkspec_.blocks,opgen.eigblocks());
     }
     
      bool operator!() const { return store_.empty(); }

      const cmatrix& operator()(size_t mzeig, size_t blk) const {
	return store_(index(mzeig,blk)); 
      } 

      cmatrix& operator()(size_t mzeig, size_t blk) {
	return store_(index(mzeig,blk)); 
      } 

     const complex& identityscale(size_t mzeig, size_t blk) const;

     const SpinOpGeneratorBase& generator() const { 
       if (!opgenp_)
	 throw Undefined("BlockedOperator: generator");
       return *opgenp_;
     }

     template<class OpGen> void mla(const OpGen& opgen,double scale, const operator_spec& spec) {
       if (!*this)
	 throw Undefined("BlockedOperator::mla");
       opgen.mla(store_,blkspec_,scale,spec);
     }

     template<class OpGen> void mla(const OpGen& opgen,double scale, const operator_spec& spec1, const operator_spec& spec2) {
       if (!*this)
	 throw Undefined("BlockedOperator::mla");
       opgen.mla(store_,blkspec_,scale,spec1,spec2);
     }
       
      const cmatrix& operator()() const {
	if (store_.size()>1)
	  throw Failed("BlockedOperator(): invalid with block structure");
	return store_.front();
      }

      cmatrix& operator()() {
	if (store_.size()>1)
	  throw Failed("BlockedOperator(): invalid with block structure");
	return store_.front();
      }
      
      const block_pattern& blockstructure() const {
	return blkspec_;
      }
     const BaseStructure& structure() const { return static_cast<const BaseStructure&>(*this); }

     BlockedOperator& operator*= (double scale) { store_*=scale; return *this; }
     BlockedOperator& operator*= (complex scale) { store_*=scale; return *this; }
     void print(std::ostream&, bool includeblocking =false) const;

      void unitary_simtrans(const BlockedMatrix<complex>&);
      void unitary_isimtrans(const BlockedMatrix<complex>&);

     void invalidatecache() { identityscale_.clear(); }
     void makeidentitycache() const;
      bool ishermitian() const { return blkspec_.isherm; }

      //Functions to be used carefully as separate use can lead to inconsistent objects
     void copy_props(const BlockedOperator&, bool dotranspose =false);
    BlockedMatrix<complex>& row() { return store_; }
    const BlockedMatrix<complex>& row() const { return store_; }

    template<class T2> BlockedOperator& operator+=(const T2& v)
      { 
	store_+=v;
	invalidatecache();
	return *this;
      }
    template<class T2> BlockedOperator& operator-=(const T2& v)
      { 
	store_-=v;
	invalidatecache();
	return *this;
      }

     //     BlockedOperator& operator= (const BlockedOperator&);
     BlockedOperator& operator+=(const BlockedOperator&);
     BlockedOperator& operator-=(const BlockedOperator&);

    template<class T2> void emultiply(const T2& v) { store_.emultiply(v); }

    void apply(const BlockedFilter&);
     void scale(double, const BlockedFilter&);
     void scale(const complex& , const BlockedFilter&);
     void clear() { store_.clear(); invalidatecache(); }

     void swap(BlockedOperator&);
     //     size_t Hblock(size_t mzn, size_t eign) const { return Hindexer_(mzn,eign); }

     void full(cmatrix&) const; //!< full, unblocked matrix
     cmatrix full() const { cmatrix tmp(mxflag::temporary); full(tmp); return tmp; }

    private:
     block_pattern blkspec_; //!< block structure for *operator*
      BlockedMatrix<complex> store_;
     const SpinOpGeneratorBase* opgenp_;
     //     Indexer<2> Hindexer_; //!< Indexer for *Hamiltonian*
     mutable List<complex> identityscale_;
    
    void dump() const;
     static complex findscale_(const Matrix<complex>&);
     template<typename T> void scale_(T, const BlockedFilter&);
  };

  inline BlockedOperator operator+ (const BlockedOperator& a, const BlockedOperator& b)
    {
      BlockedOperator c(a);
      c+=b;
      return c;
    }

  inline BlockedOperator operator- (const BlockedOperator& a, const BlockedOperator& b)
    {
      BlockedOperator c(a);
      c-=b;
      return c;
    }

 //rather crude check
inline bool arematching(const BlockedOperator& a, const BlockedOperator& b)
{ return a.blockstructure().ismatching(b.blockstructure()); }

void mla(BlockedOperator&, double, const BlockedOperator&);
void mla(BlockedOperator&, complex, const BlockedOperator&);

inline void multiply(BlockedOperator& d, const BlockedOperator& a, const BlockedOperator& b)
    {
      if (!arematching(a,b))
	throw Mismatch("BlockedOperator::multiply");
      d.copy_props(a);
      multiply(d.row(),a.row(),b.row());
    }

inline BlockedOperator operator* (const BlockedOperator& a, const BlockedOperator& b)
    {
      BlockedOperator c;
      multiply(c,a,b);
      return c;
    }

// inline void spy(std::ostream& ostr, const BlockedOperator& a)
// {
//   spy(ostr,a.row());
// }

//   inline complex trace(const BlockedOperator& a)
//     {
//       return trace(a.row());
//     }

//   inline double hermitian_trace(const BlockedOperator& a)
//     {
//       return hermitian_trace(a.row());
//     }

complex trace_multiply(const BlockedOperator&, const BlockedOperator&);
//double hermitian_trace_multiply(const BlockedOperator&, const BlockedOperator&);

  template<class T2> BlockedOperator operator+ (const BlockedOperator& a,const T2& b)
    { 
      BlockedOperator d(a);
      d+=b;
      return d;
    }

  template<class T2> BlockedOperator operator- (const BlockedOperator& a,const T2& b)
    { 
      BlockedOperator d(a);
      d-=b;
      return d;
    }

  inline std::ostream& operator<< (std::ostream& ostr, const BlockedOperator& a)
    {
      a.print(ostr);
      return ostr;
    }

  inline void unitary_simtrans(BlockedOperator& d, const BlockedOperator& a, const BlockedMatrix<complex>& U)
  {
    d=a;
    d.unitary_simtrans(U);
  }

  inline void unitary_isimtrans(BlockedOperator& d, const BlockedOperator& a, const BlockedMatrix<complex>& U)
  {
    d=a;
    d.unitary_isimtrans(U);
  }

  void transpose(BlockedOperator&, const BlockedOperator&);
  void conj_transpose(BlockedOperator&, const BlockedOperator&);

  inline BlockedOperator conj_transpose(const BlockedOperator& a) {
    BlockedOperator b;
    conj_transpose(b,a);
    return b;
  }

  inline BlockedOperator transpose(const BlockedOperator& a) {
    BlockedOperator b;
    transpose(b,a);
    return b;
  }

double norm(const BlockedOperator&);

//define how detect and sigma are combined
 inline complex NMR_trace(const BlockedOperator& detect, const BlockedOperator& sigma) {
   return (!detect) ? complex(norm(sigma)) : trace_multiply(detect,sigma);
 }

inline complex NMR_trace(const BlockedOperator& sigma) { return complex(norm(sigma)); }

// inline complex NMR_hermitian_trace(const BlockedOperator& detect, const BlockedOperator& sigma) { return hermitian_trace_multiply(detect,sigma); }

//   template<> struct Operator_traits<BlockedOperator> {
//     typedef BlockedMatrix<complex> output_type;
//     typedef BlockedMatrix<complex> genoutput_type;
//     typedef cmatrix suboutput_type;
//     typedef cmatrix gensuboutput_type;
//     typedef cmatrix& selectreturn_type;
//     static const bool isnotPOD=true;
//     static const bool trivialclone=true;
//   };

class SpinOpGenerator;

class BlockedFilter : public BaseStructure {
public:
  //template<class OpGen> BlockedFilter(const OpGen&, const BaseList<int>&);
  //template<class OpGen,class F> BlockedFilter(const OpGen&, const BaseList<size_t>&, const F&);
  BlockedFilter(const SpinOpGenerator&, const block_pattern&, const spinorder_spec&);

  const BlockedMatrix<bool>& row() const { return store_; }
  bool operator!() const { return !store_; }
  size_t size() const { return store_.size(); }
  const BaseStructure& structure() const { return static_cast<const BaseStructure&>(*this); }
  bool ismatching(const BlockedOperator& a) const { return (store_.size()==a.row().size()) && (structure()==a.structure()); }
  bool isnonzero() const;
  const block_pattern& blockstructure() const {	return blkspec_; }
  void full(Matrix<bool>&) const;
  Matrix<bool> full() const { Matrix<bool> tmp(mxflag::temporary); full(tmp); return tmp; }

  const Matrix<bool>& operator()(size_t mzeig, size_t blk) const {
    return store_(index(mzeig,blk)); 
  } 
  
  Matrix<bool>& operator()(size_t mzeig, size_t blk) {
    return store_(index(mzeig,blk)); 
  } 

private:
  BlockedMatrix<bool> store_;
  block_pattern blkspec_; //!< block structure for *operator*
  const SpinOpGeneratorBase& opgen_;

  void coherence_filter_ip(const ListList<double>& Fz, const spinorder_spec& spec, state_t cohers);
};

std::ostream& operator<< (std::ostream&, const BlockedFilter&);

