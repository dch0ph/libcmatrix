#undef LCM_SUPPRESS_VIEWS
#include <cstdlib>
#include "CrystalSystem.h"
#include "lcm_MetaPropagation.h"
#include "cmatrix.h"
#include "rmatrix.h"
#include "space_T.h"
#include "ScratchList.h"
#include "matlabio.h"
#include "ListList.h"
#include "NMR.h"

namespace libcmatrix {

  using ::std::sqrt;

  const size_t invalid_value=-1;

/****************************************************************/

  CrystalSymmetriseBase::CrystalSymmetriseBase(size_t N, size_t M, size_t neigs,int flags)
    : neigs_(neigs), ncells_(N), total_(M*N),
    useeigsym_(flags & MetaFlags::UseEigSymmetry),
    scalefacs(neigs*neigs+1,DEADBEEF)
{
  const size_t max(neigs*neigs);
  if (neigs % N)
    throw InvalidParameter("CrystalSymmetriseBase");
  for (size_t j=1;j<=max;j++) {
    if ( (max % j)==0)
      scalefacs(j)=N/sqrt(double(j));
  }
  maxstates_=long(1)<<total_;
}

//Basic 1D symmetriser
class SimpleSymmetrise : public CrystalSymmetriseBase {
public:
  SimpleSymmetrise(int N, int M, int flags_);

  complex symmetrise(Matrix<double>& store, magicsize_t,magicsize_t,int k) const { return symmetrise_(store,k); }
  complex symmetrise(Matrix<complex>& store, magicsize_t,magicsize_t,int k) const { return symmetrise_(store,k); }

  CrystalSymmetriseBase* clone() const { return new SimpleSymmetrise(*this); }
 
  void addtoevals(BaseList<complex>& evals, size_t blki, size_t blkj, size_t r, size_t c) const;
  bool haspermutation() const { return false; }

private:
  List<complex> eigfacs;
  template<typename T> complex symmetrise_(Matrix<T>&, int k) const;
};

struct simple_creator {
  explicit simple_creator(size_t N_,size_t M_) : M(M_), shift(M_*N_-M_), mask((1<<M)-1) {}
  const size_t M;
  const size_t shift;
  const state_t mask;
  magicsize_t push_states(List<state_t>& state_list, BaseList<bool> used, state_t current) const
  {
    size_t count=0;
    const state_t start=current;
    do {
      used(current)=true; //flag state used
      count++;
      state_list.push_back(current);
      current=(current>>M) | ((current & mask)<<shift); //apply translation
    } while (current!=start); //stop when we return to beginning
    return count;
  }
};

  void CrystalSymmetriseBase::makeindices()
{
  if (linkedstates_.empty())
    throw Failed("CrystalSymmetriseBase::makeindices called before linked states created");

  const long totdim=maxstates(); //total number of states

  //indexes take a state number and return the set(block) number and position (index) within set
  state_to_block_.create(totdim,invalid_value);
  state_to_index_.create(totdim,invalid_value);

  for (size_t blk=permutation_blocks();blk--;) {
    const BaseList<state_t> clist(linkedstates_(blk));
    size_t j=clist.size();
    for (;j--;) {
      state_to_block_(clist(j))=blk;
      state_to_index_(clist(j))=j;
    }
  }
}

void CrystalSymmetriseBase::make_eigfacs(BaseList<complex> eigfacs)
{
  const size_t N=eigfacs.size();
  for (size_t j=N;j--;)
    eigfacs(j)=expi(2*M_PI*double(j)/N); //store e^2pi i n/N
}

SimpleSymmetrise::SimpleSymmetrise(int N, int M_, int flags_)
  : CrystalSymmetriseBase(N,M_,N,flags_), eigfacs(N)
{
  useeigs = useeigsym_ ? 1+N/2 : N; //number of active eigenvalues
  isreal_.create(useeigs,false);
  isreal_.front()=true;
  if (!useeigsym_ || ((N & 1)==0))
    isreal_(N/2)=true;

  size_t j,k;

  make_eigfacs(eigfacs);
 
  haseig.create(N+1,useeigs,false);
  //haseig(j,k) is true if a set of size j contains eigenvalue k, where j is a factor of N
  for (j=1;j<=N;j++) {
    if ( (N % j)==0) {
      const size_t step=N/j;
      for (k=0;k<useeigs;k+=step)
	haseig(j,k)=true;
    }
  }
  simple_creator creator(N,M_);
  getsymmetry(creator); //get sets of symmetry-linked states
}
  
/***************************************************/

CrystalOpGenerator::CrystalOpGenerator(const spinhalf_system& sys_, const CrystalStructure& cstructv, int flagsv, int verbosev)
  : SpinOpGeneratorBase(cstructv,sys_.nspins(),flagsv,verbosev),
    sysp_(N_==1 ? &sys_ : sys_.clone(N_),N_==1 ? mxflag::nondynamic : mxflag::normal) 
{
  create();
}
      
CrystalOpGenerator::CrystalOpGenerator(const spinhalf_system& sys_, const CrystalStructure& cstructv, const char *nuclabel, int flagsv, int verbosev)
  :  SpinOpGeneratorBase(cstructv,sys_.nspins(),flagsv,verbosev),
     sysp_(N_==1 ? &sys_ : sys_.clone(N_),N_==1 ? mxflag::nondynamic : mxflag::normal)
{
  create(ExplicitList<1,nuclei_spec>(nuclabel));
}

  CrystalOpGenerator::CrystalOpGenerator(const spinhalf_system& sys_, const CrystalStructure& cstructv, const BaseList<nuclei_spec>& blocknucsv, int flagsv, int verbosev)
    : SpinOpGeneratorBase(cstructv,sys_.nspins(),flagsv,verbosev),
      sysp_(N_==1 ? &sys_ : sys_.clone(N_),N_==1 ? mxflag::nondynamic : mxflag::normal)
{
  create(blocknucsv);
}

CrystalOpGenerator::CrystalOpGenerator(const HamiltonianStructure& structurev, const CrystalStructure& cstructv, int verbosev)
  : SpinOpGeneratorBase(cstructv,structurev,verbosev),
    sysp_(new spinhalf_system(structurev.spinsystem(),N_),mxflag::normal)
{
  if (quadrupole_order()!=1)
    throw Failed("CrystalOpGenerator doesn't support 2nd order quadrupoles");

  create(structurev.blockingnuclei());
 }

template<typename T> void CrystalOpGenerator::create(BlockedMatrix<T>& dest, const block_pattern& blkspec) const
{
  block_pattern::iterator iter(blkspec);
  //  const size_t limit=actual_mzblocks();
  const size_t eigblks(eigstr_.eigblocks());
  const size_t mzblks=blkspec.blocks;
  const size_t totblks(eigblks*mzblks);
//   size_t mzlevels=blkspec.mzlevels;
//   if (blkspec.hasmiddle)
//     mzlevels=(mzlevels+1)/2; //actual cutoff
  
  ScratchList<size_t> rstr(totblks);
  ScratchList<size_t> cstr(totblks);

  size_t r,c,mz;
  bool ismiddle;

  while (iter.next(r,c,mz,ismiddle)) {
//     if (r>=mzlevels) {
//       r=c;
//       assert(blkspec.hasmiddle);
//     }
//     else {
//       if (c>=mzlevels) {
// 	c=r;
// 	assert(blkspec.hasmiddle);
//       }
//     }
 //    if (ismiddle) {
//       if (c>r)
// 	c=r;
//       else
// 	r=c;
//     }
    for (size_t k=eigblks;k--;) {
      const size_t ind(eigstr_.index(mz,k));
      rstr(ind)=eigsizes(r,k);
      cstr(ind)=eigsizes(c,k);
    }
  }
  if (verbose()>1)
    std::cout << "Creating empty BlockedMatrix:  rows: " << rstr << "  cols: " << cstr << '\n';
  dest.create(rstr,cstr);
}

  //iterates over a set of blocks
class CrystalSystem_iterator {
public:  
  CrystalSystem_iterator(const CrystalOpGenerator&); //zero coherence blocks
  CrystalSystem_iterator(const CrystalOpGenerator&, const block_pattern&, size_t nuc, int coher, bool isupper =false);
  bool next(size_t&, size_t&);

  //are we blocked by coherence?
  //  bool isblocked() const { return sym.isblocked(); }

  template<typename T> void addsym(MatrixTensor<complex>&, size_t blki, size_t blkj, int l, int m, int k, bool, Type2Type<T>);
  template<typename T> void addsym(BlockedMatrixTensor<complex>&, size_t blki, size_t blkj, int l, int m, bool, Type2Type<T>);
  template<class T> void addsym(BlockedMatrix<complex>&, size_t, size_t, Type2Type<T>, bool);
  template<class T> void addsym(cmatrix&, size_t, size_t, int k, const T&, bool);
  template<class T> void addsym(complex&, size_t, size_t, int k, const T&, bool);

  void add(BlockedMatrix<complex>&, const BaseList<complex>&, size_t, size_t) const;
  void add_hermitian(BlockedMatrix<complex>&, const BaseList<complex>&, size_t, size_t) const;
  
  void add(cmatrix&, complex, size_t, size_t, int) const;
  void add_hermitian(cmatrix&, complex, size_t, size_t, int) const;

  void mla(BlockedMatrix<complex>&, double,const operator_spec&);
  bool Ipevalues(DynamicList<complex>&, size_t blki,size_t blkj,int);
  bool Fpevalues(DynamicList<complex>&, size_t blki,size_t blkj,size_t nuc =NULL_NUCLEUS);

  void docreate();

  bool isdiagonal() const {
    return (rmzind==cmzind);
  }
  void ensure_null(BlockedMatrixTensor<complex>&, int l) const;
  void ensure_null(BlockedMatrix<complex>&) const;
  
  size_t rows(size_t k) const {
    return rblkstr(k); }
  size_t cols(size_t k) const {
    return cblkstr(k); }

  template<class T> bool Hcoupling_matrix(const Matrix<T>& couplings,int,int,size_t blki,size_t blkj);  

  bool A2_matrix(const Matrix<space_T>& tensors,size_t blki,size_t blkj,int m);
  bool A2_matrix(const rmatrix& couplings,size_t blki,size_t blkj) { return Hcoupling_matrix(couplings,2,-1,blki,blkj); }
  bool A2_matrix(const Matrix<space_T>& couplings,size_t blki,size_t blkj) { return Hcoupling_matrix(couplings,2,-1,blki,blkj); }
  bool A0_matrix(const rmatrix& couplings,size_t blki,size_t blkj) { return Hcoupling_matrix(couplings,1,1,blki,blkj); }

  size_t get_block(int eval) const {
    return sym.eigstr_.index(mzcount,eval);
  }

  rmatrix& get_store(size_t rs, size_t cs, const Type2Type<double>&)
  { rmatrix& res(rstores(rs,cs));
    if (!res)
      res.create(rs,cs);
    return res;
  }

  cmatrix& get_store(size_t rs, size_t cs, const Type2Type<complex>&)
  { cmatrix& res(cstores(rs,cs));
    if (!res)
      res.create(rs,cs);
    return res;
  }

  double symmetrise0(const Type2Type<double>&) const {
    return ncells*H00;
  }
  
  complex symmetrise0(const Type2Type<complex>&) const {
    return float_t(ncells)*H00c;
  }

  double symmetrise0(size_t rs,size_t cs, const Type2Type<double>&) const {
    return scalefacs(rs*cs)*sum(rstores(rs,cs)); }
  
  complex symmetrise0(size_t rs,size_t cs, const Type2Type<complex>&) const {
    return scalefacs(rs*cs)*sum(cstores(rs,cs)); }
     
private:
  const CrystalOpGenerator& sym;
  const CrystalSymmetriseBase& symmetriser;
  size_t useeigs;
  const List<double>& scalefacs;
  const List<size_t>& state_to_block;
  const List<size_t>& state_to_index;
  int coher_; //coherence being considered
  BaseList<size_t> nucmzinds;
  size_t ncells;
  size_t neigs;
  bool finished; //true if there are no blocks left
  bool isupper; //true if problem is symmetric and we are only iterating over upper diagonal
  size_t rmzind,cmzind; //current row/col block
  size_t mzcount;
  bool bumpmz;
  BaseList<size_t> cwhich; //current set (row,column) of states
  BaseList<size_t> rwhich;
  size_t brablks,ketblks; //number of bra/ket states
  size_t rows_,cols_; //size of mz2 block
  size_t crow,ccol; //current bra/ket index
  Matrix<rmatrix> rstores;
  Matrix<cmatrix> cstores;
  double H00;
  complex H00c;
  DynamicList<complex> evalues;
  BaseList<size_t> rblkstr;
  BaseList<size_t> cblkstr;
  int verbose;
  //  bool ismiddle;
  block_pattern::iterator mziter;

  void resetmz2();
  bool nextmz2() {
    mzcount++;
    size_t mzeigSD; //!< ignored - could be merged with mzcount, but don't want to touch!
    if (mziter.next(rmzind,cmzind,mzeigSD)) {
      resetmz2();
      return true;
    }
    return false;
  }

  void init(size_t selnuc);
};

class CrystalSystem_diagiterator {
 public:
  CrystalSystem_diagiterator(const CrystalOpGenerator&);
  bool next(size_t&);

  size_t get_block(int eval) const {
    return sym.eigstr_.index(mzind,eval);
  }

private:
  const CrystalOpGenerator& sym;
  BaseList<size_t> rwhich; //current set of states
  int blki; //current index
  size_t mzind;
  bool bumpmz;

  void resetmz2();
};

CrystalSystem_diagiterator::CrystalSystem_diagiterator(const CrystalOpGenerator& Sym)
  : sym(Sym)
{
  mzind=Sym.actual_mzblocks()-1;
  resetmz2();
}

void
CrystalSystem_diagiterator::resetmz2()
{
  rwhich.create(sym.blockindices(mzind));
  blki=rwhich.length()-1;
  bumpmz=false;
}

  bool
  CrystalSystem_diagiterator::next(size_t& blk) {
    if (bumpmz) {
      mzind--;
      resetmz2();
    }
    if (blki<0)
      return false; //finished
    blk = rwhich(blki--);
    if ((blki<0) && mzind)
      bumpmz=true;
    return true;
  }
  
CrystalSystem_iterator::CrystalSystem_iterator(const CrystalOpGenerator& Sym)
  : sym(Sym), symmetriser(*(Sym.symmetriserp_)), 
    scalefacs(symmetriser.scale_factors()),
    state_to_block(symmetriser.state_to_block()),
    state_to_index(symmetriser.state_to_index()),
    coher_(0),
    isupper(false),
    mziter(Sym.diag_blkspec_)
{
  init(Sym.defaultnucleus_);
}

CrystalSystem_iterator::CrystalSystem_iterator(const CrystalOpGenerator& Sym, const block_pattern& blkspec, size_t selnuc, int coherv, bool isupper_)
  : sym(Sym), symmetriser(*(Sym.symmetriserp_)), 
    scalefacs(symmetriser.scale_factors()),
    state_to_block(symmetriser.state_to_block()),
    state_to_index(symmetriser.state_to_index()),
    coher_(coherv),
    isupper(isupper_),
    mziter(blkspec)
{
  if (isupper && (coherv!=0))
    throw InternalError("CrystalSystem_iterator: isupper incompatible with non-zero coherences");
  init(selnuc);
}

void CrystalSystem_iterator::init(size_t selnuc)
{
  useeigs=symmetriser.actual_eigenvalues();
  ncells=symmetriser.ncells();
  neigs=symmetriser.eigenvalues();
  finished=false;
  rstores.create(neigs+1,neigs+1);
  cstores.create(neigs+1,neigs+1);
  //  evalues.create(symmetriser.actual_eigenvalues(),mxflag::normal);
  evalues.create(symmetriser.actual_eigenvalues()); //!< removed normal flag (May 10) Why there?
  verbose=sym.verbose();
  if (selnuc)
    nucmzinds.create(sym.allmzinds_.row(sym.nuctoindex(selnuc)));
  bumpmz=false;
  nextmz2();
  mzcount=0;
}

  //reset after changing rmzind/cmzind
void
CrystalSystem_iterator::resetmz2()
{
  size_t usermz=rmzind;
  size_t usecmz=cmzind;
  const size_t maxind=sym.mzsizes.size();
  if (rmzind>=maxind)
    usermz=cmzind;
  else {
    if (cmzind>=maxind)
      usecmz=rmzind;
  }
//     if (rmzind>cmzind)
//       usermz=cmzind;
//     else
//       usecmz=rmzind;
//   }

  rblkstr.create(sym.eigsizes.row(usermz)); //get block structure for given mz pair
  cblkstr.create(sym.eigsizes.row(usecmz));
  
  rows_=sym.mzsizes(usermz); //size of block
  cols_=sym.mzsizes(usecmz);
  
  //if single state set, store set number for quick reference
  //store junk value otherwise to catch attempts to use firstrblk in other cases
//   if (rows_==1) 
//     firstrblk=int(sym.mzblocking_(usermz).front());
//   else
//     firstrblk=-1;

  //must use actual indices here!
  rwhich.create(sym.blockindices(rmzind)); //update rwhich etc. to point to current set of bra states
  cwhich.create(sym.blockindices(cmzind));
  brablks=rwhich.length();
  ketblks=cwhich.length();
  crow=ccol=0;
}

bool
CrystalSystem_iterator::next(size_t& rblk, size_t& cblk)
{
  int blkcoher;
  do {
    if (bumpmz) {
      if (!nextmz2())
	finished=true;
      else {
	resetmz2();
	bumpmz=false;
      }
    }

    if (finished)
      return false; //iterator is exhausted?
    
    rblk=rwhich(crow); //read current position
    cblk=cwhich(ccol);

    ccol++; //increment col index
    if (ccol==ketblks) {
      crow++; //increment row index
      ccol= isupper ? crow : 0; //upper diagonal only?
      if (crow==brablks) //finished state block?
	bumpmz=true;
    }
    if (nucmzinds.empty())
      return true;
    blkcoher=(int)nucmzinds(cblk)-(int)nucmzinds(rblk);
  } while (blkcoher!=coher_);
  return true;
}

  //construct dipolar matrix for given set of bra/ket states
  //ignore couplings to spins of type nuc (unless NULL_NUCLEUS)

template<class T> bool
CrystalSystem_iterator::Hcoupling_matrix(const Matrix<T>& couplings, int scalez, int scalexy, size_t blki, size_t blkj)
{
  if (sym.nspins()!=couplings.cols())
    throw Mismatch("Hcoupling_matrix: coupling matrix doesn't match CrystalOpGenerator");

  const ListList<state_t>& linkedstates(symmetriser.linkedstates());
  const BaseList<state_t> cstates(linkedstates(blkj));
  const size_t blkisize=linkedstates.size(blki);
  const size_t blkjsize=linkedstates.size(blkj);

  rmatrix& H(rstores(blkisize,blkjsize)); //where we'll store the temporary matrix
  bool issingle=(blkisize==1) && (blkjsize==1); //single element?

  reducer_<double,T> reduce;

  //clear temporary matrix
  if (issingle)
    H00=0.0;
  else {
    if (!H) //if temp. matrix doesn't exist, make it
      H.create(blkisize,blkjsize);
    H=0;
  }

  const basespin_system& sys=*(sym.sysp_);
  const double scalez2=scalez*0.5;
  const size_t total(sym.nspins());
  const size_t M(sym.nspins_cell());

  for (size_t j=0;j<M;j++) {
    const state_t maskj=maskelement(total,j);
    
    for (size_t k=0;k<M;k++) { //loop over spin pairs

      const bool ishomo=(sys(j)==sys(k)); //homonuclear coupling
      
      for (size_t cell=0;cell<ncells;cell++) { //loop over cells
	if ((j>k) || cell) { //don't count both (0,j)-(0,k) and (0,k)-(0,j)
	  const size_t sk=sym.cell_to_spin(cell,k); //index of spin (cell,k)
	  const state_t masksk=maskelement(total,sk);

	  const T& curcoup((j<sk) ? couplings(j,sk) : couplings(sk,j));
	  if (!curcoup)
	    continue;
	  const double coup=(cell ? 0.25 : 0.5)*reduce(curcoup);
	  if (coup==0.0)
	    continue;
	  
	  const state_t mask=maskj | masksk; //has bits corresponding to spins j & sk set
	  
	  for (size_t jj=blkjsize;jj--;) {
	    const state_t cstate=cstates(jj); //current ket state
	    const bool isflip=!(cstate & maskj) ^ !(cstate & masksk);
	    if (blki==blkj) { //zz components
	      const double del= isflip ? -coup : coup;
	      if (issingle) 
		::libcmatrix::mla(H00,scalez2,del);
	      else
		::libcmatrix::mla(H(jj,jj),scalez2,del);
	    }
	    if (ishomo && isflip) {
	      const state_t nstate=cstate ^ mask; //flip spin states
	      if (state_to_block(nstate)==blki ) { //is new spin state within bra states?
		if (issingle)
		  H00+=scalexy*coup;
		else
		  H(state_to_index(nstate),jj)+=scalexy*coup;
	      }
	    }
	  }
	}
      }
    }
  }
  if (verbose>1) {
    std::cout << "H" << blki << "," << blkj;
    if (issingle)
      std::cout  << ": " << H00 << "\n";
    else
      std::cout << ":\n" << H << "\n";
  }
  return issingle;
}
  //return true if result is single element

  //corresponding function for spinning Hamiltonian (Fourier component m)
bool 
CrystalSystem_iterator::A2_matrix(const Matrix<space_T>& tensors,size_t blki, size_t blkj, int m)
{
  if (sym.nspins_cell()!=tensors.rows())
    throw Mismatch("Hcoupling_matrix");
  
  const ListList<state_t>& linkedstates(symmetriser.linkedstates());
  const BaseList<state_t> cstates(linkedstates(blkj));

  const size_t blkisize=linkedstates.size(blki);
  const size_t blkjsize=linkedstates.size(blkj);

  cmatrix& H=cstores(blkisize,blkjsize);
  const bool issingle=(blkisize==1) && (blkjsize==1);

  if (issingle)
    H00c=0.0;
  else {
    if (!H) H.create(blkisize,blkjsize);
    H=0;
  }

  const basespin_system& sys=*(sym.sysp_);
  const size_t total(sym.nspins());
  const size_t M(sym.nspins_cell());

  for (size_t j=0;j<M;j++) {
    const state_t maskj=maskelement(total,j);

    for (size_t k=0;k<M;k++) {
      const bool ishomo=(sys(j)==sys(k));

      for (size_t cell=0;cell<ncells;cell++) {
	if ((j>k) || cell) {
	  const size_t sk=sym.cell_to_spin(cell,k);
	  const state_t masksk=maskelement(total,sk);

	  const space_T& curtens((j<sk) ? tensors(j,sk) : tensors(sk,j));

	  if (!curtens)
	    continue;

	  const state_t mask=maskj | masksk;
	  const complex coup=(cell ? 0.25 : 0.5)*curtens(2,m);

	  for (size_t jj=blkjsize;jj--;) {
	    const state_t cstate=cstates(jj);
	    const bool isflip = !(cstate & maskj) ^ !(cstate & masksk);
	    if (blki==blkj) {
	      const complex del= isflip ? -coup : coup;
	      if (issingle) 
		H00c+=del;
	      else
		H(jj,jj)+=del;
	    }
	    if (ishomo && isflip) {
	      const state_t nstate=cstate ^ mask;
	      if (state_to_block(nstate)==blki ) {
		if (issingle)
		  H00c-=coup;
		else
		  H(state_to_index(nstate),jj)-=coup;
	      }
	    }
	  }
	}
      }
    }
  }
  if (verbose>1) {
    std::cout << "H(" << m << ")" << blki << "," << blkj;
    if (issingle)
      std::cout  << ": " << H00c << "\n";
    else
      std::cout << ":\n" << H << "\n";
  }
  
  return issingle;
}

  //return value corresponding to eigenvalue m of H transformed into symmetrised basis
template<class T> complex 
SimpleSymmetrise::symmetrise_(Matrix<T>& H, int m) const
{
  if (m==0)
    throw Failed("symmetrise: don't use for k=0!");

  const size_t blkisize=H.rows();
  const size_t blkjsize=H.cols();

  if (!haseig(blkisize,m) || !haseig(blkjsize,m))  //don't have that eigenvalue for this no. of bra/ket states
    throw Failed("symmetrise: eigenvalue doesn't exist in given block");

  complex s(0.0);
      
  for (size_t ii=blkisize;ii--;) {
    const BaseList<T> Hii=H.row(ii);
    size_t wh=(m*(ncells_-ii)) % ncells_;
    
    for (size_t jj=0;jj<blkjsize;jj++) {
      //      const size_t wh=(m*(jj-ii+N)) % N;
      mla(s,Hii(jj),eigfacs(wh));
      wh+=m;
      if (wh>=ncells_) wh-=ncells_;
    }
  }
  return s*scalefacs(blkisize*blkjsize);
}


  //return 2 mz for given nucleus type and state
int
CrystalOpGenerator::mz2val(size_t nuc,state_t state) const
{
  int retval=0;
  const size_t total(nspins());

  for (size_t j=0;j<nspins_cell();j++) { //loop over spins in unit cell
    if ((*sysp_)(j).nucleus()!=nuc)
      continue;

    for (size_t i=0;i<N_;i++) { //loop over cells
      if (state & maskelement(total,cell_to_spin_(i,j))) 
	//count -1/2 for bit set, +1/2 for not set
	retval--;
      else
	retval++;
    }
  }
  return retval;
}

  //return 2 Fz for given state (val) and total number of spins
  // = (bits unset)-(bits set) = total bits - 2 (bits set)
int
CrystalOpGenerator::mz2val(state_t val) const
{
  int retval(nspins());
  
  while (val) { //while still some bits set
    if (val & 1)
      retval-=2;
    
    val>>=1;
  }
  return retval;
}

  //compute symmetrised Iz for given spin
ListList<double>
diag_Iz(const CrystalOpGenerator& sys, size_t spinn)
{
  return sys.diag_Iz(spinn);
}

void
CrystalOpGenerator::mla_Fz(ListList<double>& dest, double scale, nuclei_spec whichn) const 
{
  const size_t nuc=whichn();
  if (nuc!=NULL_NUCLEUS)
    ::libcmatrix::mla(dest,scale,diag_Fz(nuc));
  else {
    if (tzops_.empty())
      throw Failed("Can't compute unspecific Fz here");
    for (size_t j=nspins_cell();j--;)
      ::libcmatrix::mla(dest,scale,tzops_(j));
  }
}

void
CrystalOpGenerator::rawdiag_Fz(ListList<double>& dest, size_t nuc) const
{
  const basespin_system& sys(spinsystem());
  List<size_t> whichspins(nspins_cell());
  whichspins.create(0);
  for (size_t i=nspins_cell();i--;) {
    if (sys(i).nucleus()==nuc)
      whichspins.push_back(i);
  }
  make_tzop(dest,whichspins);
}

  // construct symmetrised z operator
void
CrystalOpGenerator::mla_Iz(ListList<double>& tzop, double scale, size_t spinn) const
{
  if (tzops_.empty())
    throw Failed("Can't compute Iz when permutation active");
  ::libcmatrix::mla(tzop,scale,tzops_(spinn));
}

void
CrystalOpGenerator::make_tzop(ListList<double>& tzop, const BaseList<size_t>& spins) const
{
  create(tzop);
#ifndef NDEBUG
  tzop=DEADBEEF; //easier to spot problems
#endif

  const size_t n(spins.size());
  ScratchList<state_t> masks(n);
  const size_t total(nspins());

  size_t j,k,ablk;
  for (j=n;j--;)
    masks(j)=maskelement(total,spins(j));

  CrystalSystem_diagiterator iter(*this);
  const ListList<state_t>& linkedstates(symmetriserp_->linkedstates());
  const size_t useeigs=symmetriserp_->actual_eigenvalues();

  while (iter.next(ablk)) { //loop over state sets
    const BaseList<state_t> states(linkedstates(ablk));
    const size_t blksize=states.length();  
    const BaseList<int> ptrs(eigptrs.row(ablk)); //stash value of eigptrs for set

    double sum=0.0;
    for (j=n;j--;) {
      const state_t& mask(masks(j));
      for (k=blksize;k--;) //loop over states
	sum+=Izelement_raw(states(k) & mask);
    }

    sum*=N_/blksize;
    for (k=useeigs;k--;) {//store result (sum) in correct position (ptrs(k)) of each eigenvalue block (k)
      const size_t ind(iter.get_block(k));
      const int ptr(ptrs(k));
      if (ptr>=0) //easy way to check that eigenvalue is present
	(tzop(ind))(size_t(ptr))=sum;
    }
  }
  if (verbose()>1)
    std::cout << "Computed Fz for spins " << spins << ": " << tzop << '\n';
}

/* given system spanning reduced Hilbert space construct
   "reverse index" taking state to index within sub-space */
void
makerindex(List<size_t>& rind,const spinhalf_system& sys)
{
  const BaseList<state_t> states=sys.ketstates();

  rind.create(sys.size());
  rind=-1; //fill with dummy
  
  for (size_t i=states.length();i--;) //loop over sub-space
    rind(states(i))=i;
}


  //return sum over columns in source in dest
template<class T> void
colsum(List<T>& dest,const Matrix<T>& source)
{
  if (!source)
    throw Undefined("colsum");

  size_t r=source.rows();
  dest=source.row(--r);
  for (;r--;)
    dest+=source.row(r);
}
    
void CrystalOpGenerator::makemzinds(Matrix<size_t>& allmzinds, List<size_t>& maxmzinds, const BaseList<size_t>& nucs)
{
  const size_t ntypes=nucs.size();
  const ListList<state_t>& linkedstates(symmetriserp_->linkedstates());
  const size_t permblocks=linkedstates.size();
  
  allmzinds.create(ntypes,permblocks);
  maxmzinds.create(ntypes);
    
  size_t nuci,j;
  for (nuci=0;nuci<ntypes;nuci++) { //calculate 2mz for all blocks and all nuclei
    const size_t tmpnuc=indextonuc(nuci);
    BaseList<size_t> curmz2(allmzinds.row(nuci));
    const size_t maxI2(sysp_->nspins(tmpnuc));
    for (j=permblocks;j--;)
      curmz2(j)=(maxI2-mz2val(tmpnuc,linkedstates(j,0)))/2;
    maxmzinds(nuci)=size_t(maxI2+1.5);
  }
}

//Common part of CrystalOpGenerator constructor (private)
void
CrystalOpGenerator::create(const BaseList<nuclei_spec>& nucs)
{
  if (M_<1 || N_<1)
    throw InvalidParameter("CrystalOpGenerator"); 

  if (cstruct.dimensions()<2 && !cstruct.haspermutation())
    symmetriserp_.reset(new SimpleSymmetrise(N_,M_,flags()));
  else {
#ifdef LCM_ENABLE_GENERICPERIODIC
    symmetriserp_.reset(new GeneralSymmetrise(cstruct,M_,flags(),verbose()));
#else
    throw Failed("library compiled without support for periodicity in >1 dimension");
#endif
  }
  
  const size_t useeigs=symmetriserp_->actual_eigenvalues();
  const ListList<state_t>& linkedstates(symmetriserp_->linkedstates());
  const size_t permblocks=linkedstates.size();
  if (verbose()) {
    std::cout << "Found " << permblocks << " sets of symmetry-linked states\n";
    if (verbose()>1)
      std::cout << linkedstates << '\n';
  }

  init_blockstr(useeigs,nucs);

  defaultnucleus_= (indextonuc.size()==1) ? indextonuc.front() : NULL_NUCLEUS;

  //find sizes of eigenvalue blocks
  const size_t nblocks(actual_mzblocks());
  eigsizes.create(nblocks,useeigs,size_t(0));
  eigptrs.create(permblocks,useeigs,-1);

  size_t j,k;

  List<int> eigcount(useeigs);
  const CrystalSymmetriseBase& symmetriser(*symmetriserp_);

  List< List<size_t> > subsizes(useeigs*mzblocking_.size());

  for (size_t indexi=mzblocking_.size();indexi--;) { //loop over mz blocks

    const bool included(indexi<nblocks);
    const ListList<size_t> curmzblocking(mzblocking_(indexi));
    eigcount=0;

    for (size_t minorj=curmzblocking.size();minorj--;) {
      const BaseList<size_t> cwhich(curmzblocking(minorj));
      for (j=cwhich.length();j--;) { //loop over each set
	size_t blk=cwhich(j); 
	const magicsize_t shape=symmetriser.blockshape(blk);
	const BaseList<bool> whichvalues(symmetriser.which_eigenvalues(shape));
	
	//eigptrs(blk,k) gives the index within the symmetrised block for eigenvalue k generated from set blk
	for (k=useeigs;k--;) {
	  if (whichvalues(k)) {
	    eigptrs(blk,k)=eigcount(k)++;
	    if (included) {
	      List<size_t>& cursubsize(subsizes(index(indexi,k)));
	      if (cursubsize.empty())
		cursubsize.create(curmzblocking.size(),size_t(0U));
	      cursubsize(minorj)++;
	    }
	  }
	}
      }
    }
    if (included) {
      BaseList<size_t> csizes(eigsizes.row(indexi));
      csizes+=eigcount;
    }
  }
  if (verbose()) {
    std::cout << "Eigenvalue distribution:\n" << eigsizes;
    double flopc=0;
    const BaseList<size_t> asrow(eigsizes.row());
    for (j=asrow.size();j--;) {
      const size_t csize(asrow(j));
      flopc+=csize*csize*csize;
    }
    std::cout << "Flop count: " << flopc << '\n';
    if (verbose()>1)
      std::cout << "Eigenvalue pointers:\n" << eigptrs;
  }

  mzsizes.create(nblocks);
  for (j=nblocks;j--;) {
    const BaseList<size_t> blks(blockindices(j));
    size_t mzsize=0;
    for (k=blks.length();k--;) {
      const size_t blk(blks(k));
      mzsize+=linkedstates.size(blk);
    }
    mzsizes(j)=mzsize;
  }
  if (verbose())
    std::cout << "Sizes of mz blocks: " << mzsizes << '\n';

  diagstr_.create(nblocks*useeigs); //diagonal structure
  for (size_t mz=nblocks;mz--;) { 
    for (k=useeigs;k--;) {
      const size_t ind(eigstr_.index(mz,k));
      diagstr_(ind)=eigsizes(mz,k);
    }
  } 
  if (verbose()>1)
    std::cout << "Overall diagonal structure: " << diagstr_ << '\n';

  if (flags() & MetaFlags::UsePartitioning) {
    partitioned_diagstr_=subsizes;
    init_partitioning();
  }

  //don't built z operators if permutation symmetry active
  if (!symmetriserp_->haspermutation()) {
    tzops_.create(M_);
    for (size_t m=0;m<M_;m++) {
      make_tzop(tzops_(m),m);
      if (verbose())
	std::cout << "Iz_" << m << ": " << tzops_(m) << '\n';
    }
  }
  diag_blkspec_=block_pattern(*this);
}

void CrystalSymmetriseBase::print(std::ostream& ostr) const
{
  ostr << "Linked states: " << linkedstates_ << "\n";
  ostr << "Eigenvalue pattern:\n" << haseig;
}

std::ostream& operator<< (std::ostream& ostr, const CrystalOpGenerator& a)
{
  if (a.N_==0) return ostr << "<undefined>\n";

  ostr << "Total spin system: " << *(a.sysp_) << "\n";
  ostr << "Number of cells: " << a.N_ << "\n";

  const SpinOpGeneratorBase& asbase(a);
  asbase.print(ostr);
  if (a.verbose()) {
    a.symmetriserp_->print(ostr);
    ostr << "Eigenvalue distribution:\n" << a.eigsizes << "\n";
    if (a.verbose()>1)
      ostr << "Ptrs:\n" << a.eigptrs << "\n";
  }
  return ostr;
}

void
CrystalSystem_iterator::ensure_null(BlockedMatrix<complex>& Hmats) const
{
  if (!Hmats) {
    sym.create(Hmats);
    Hmats=complex(0.0);
  }
}

void
CrystalSystem_iterator::ensure_null(BlockedMatrixTensor<complex>& Htens, int l) const
{
  if (!Htens) {
    sym.create(Htens);
    Htens=complex(0.0);
  }
  if (!(Htens.front().have_rank(l)))
    Htens.ensure_rank(l,complex(0.0));
}
  
  //add symmetrised block (blki,blkj) to spinning Hamiltonian component m
template<typename T> void
CrystalSystem_iterator::addsym(MatrixTensor<complex>& Htens, size_t blki,size_t blkj,int l, int m, int eval,bool issingle, Type2Type<T> realcomplex)
{
  if (issingle && (eval!=0))
    return;

  BaseList<int> reigptrs(sym.eigptrs.row(blki));
  BaseList<int> ceigptrs(sym.eigptrs.row(blkj));

  //k=0 is special case
  if (eval==0) {
    const size_t blkisize(symmetriser.blocksize(blki));
    const size_t blkjsize(symmetriser.blocksize(blkj));
    const complex val= complex(issingle ? symmetrise0(realcomplex) : symmetrise0(blkisize,blkjsize,realcomplex));
    if (Htens.ismatrix()) {
      cmatrix& Hsym=Htens.matrix(l,m);
      Hsym(reigptrs.front(),ceigptrs.front())+=val;
    }
    else
      Htens.element(l,m)+=val;
  }
  else {
    const magicsize_t blkishape(symmetriser.blockshape(blki));
    const magicsize_t blkjshape(symmetriser.blockshape(blkj));
    const int dr=reigptrs(eval);
    const int dc=ceigptrs(eval);

    if ((dr>=0) && (dc>=0)) { //if eigenvalue exists
      const size_t lrows(symmetriser.blocksize(blki));
      const size_t lcols(symmetriser.blocksize(blkj));
      const complex val=symmetriser.symmetrise(get_store(lrows,lcols,realcomplex),blkishape,blkjshape,eval);
      if (Htens.ismatrix()) {
	cmatrix& Hsym=Htens.matrix(l,m);
	Hsym(dr,dc)+=val;
      }
      else {
	assert((dr==0) && (dc==0));
	Htens.element(l,m)+=val;
      }
    }
  }
}

template<typename T> void
CrystalSystem_iterator::addsym(BlockedMatrixTensor<complex>& Htens, size_t blki,size_t blkj,int l,int m, bool issingle, Type2Type<T> realcomplex)
{
  if (issingle) {
    const size_t ind(get_block(0));
    addsym(Htens(ind),blki,blkj,l,m,0,true,realcomplex);
    return;
  }
  for (size_t eval=useeigs;eval--;) {
    const size_t ind(get_block(eval));
    addsym(Htens(ind),blki,blkj,l,m,eval,false,realcomplex);
  }
}

template<class T> void
CrystalSystem_iterator::addsym(complex& Hk,size_t blki,size_t blkj,int eval,const T& type_, bool issingle)
{
  assert((sym.eigptrs(blki,eval)==0) && (sym.eigptrs(blkj,eval)==0));

  if (issingle) {
    if (eval==0)
      Hk+=symmetrise0(type_);
    return;
  }

  const size_t lrows(symmetriser.blocksize(blki));
  const size_t lcols(symmetriser.blocksize(blkj));
  if (eval==0)
    Hk+=symmetrise0(lrows,lcols,type_);
  else {
    const magicsize_t rs=symmetriser.blockshape(blki);
    const magicsize_t cs=symmetriser.blockshape(blkj);
    Matrix<T>& store(get_store(lrows,lcols,type_));
    Hk+=symmetriser.symmetrise(store,rs,cs,eval);
  }
}

template<class T> void
CrystalSystem_iterator::addsym(cmatrix& Hk,size_t blki,size_t blkj,int eval,const T& type_, bool issingle)
{
  const int dr=sym.eigptrs(blki,eval);
  const int dc=sym.eigptrs(blkj,eval);

  if (issingle) {
    if (eval==0)
      Hk(dr,dc)+=symmetrise0(type_);
    return;
  }

  const size_t lrows=symmetriser.blocksize(blki);
  const size_t lcols=symmetriser.blocksize(blkj);
  if (eval==0)
    Hk(dr,dc)+=symmetrise0(lrows,lcols,type_);
  else {
    if ((dc>=0) && (dr>=0)) {
      const magicsize_t rs=symmetriser.blockshape(blki);
      const magicsize_t cs=symmetriser.blockshape(blkj);
      Hk(dr,dc)+=symmetriser.symmetrise(get_store(lrows,lcols,type_),rs,cs,eval);
    }
  }
}

template<class T> void
CrystalSystem_iterator::addsym(BlockedMatrix<complex>& Hmats,size_t blki,size_t blkj, Type2Type<T> type_, bool issingle)
{
  if (issingle) {
    const size_t ind(get_block(0));
    addsym(Hmats(ind),blki,blkj,0,type_,true);
    return;
  }
  BaseList<int> reigptrs(sym.eigptrs.row(blki));
  BaseList<int> ceigptrs(sym.eigptrs.row(blkj));
  
  cmatrix& Hsym0(Hmats(get_block(0)));
  Hsym0(reigptrs.front(),ceigptrs.front())+=symmetrise0(symmetriser.blocksize(blki),symmetriser.blocksize(blkj),type_);

  const size_t lrows=symmetriser.blocksize(blki);
  const size_t lcols=symmetriser.blocksize(blkj);
  const magicsize_t rs=symmetriser.blockshape(blki);
  const magicsize_t cs=symmetriser.blockshape(blkj);
  Matrix<T>& store(get_store(lrows,lcols,type_));

  for (size_t eval=1;eval<useeigs;eval++) {
    const int dr=reigptrs(eval);
    const int dc=ceigptrs(eval);

    if ((dc>=0) && (dr>=0)) { 
      const size_t ind(get_block(eval));
      cmatrix& Hsym(Hmats(ind));
      Hsym(size_t(dr),size_t(dc))+=symmetriser.symmetrise(store,rs,cs,eval);
    }
  }
}

void
CrystalOpGenerator::add_A2(BlockedMatrixTensor<complex>& Htens, const Matrix<space_T>& tensors) const
{
  size_t blki,blkj;
  CrystalSystem_iterator iter(*this);

  iter.ensure_null(Htens,2);
  
  while (iter.next(blki,blkj)) {
    for (int m=-2;m<=2;m++) {
      const bool issingle=iter.A2_matrix(tensors,blki,blkj,m);
      iter.addsym(Htens,blki,blkj,2,m,issingle,Type2Type<complex>());
    }
  }
}

void
CrystalOpGenerator::add_A0(BlockedMatrixTensor<complex>& Htens, const Matrix<double>& A0s) const
{
  size_t blki,blkj;
  CrystalSystem_iterator iter(*this);

  iter.ensure_null(Htens,0);
  
  while (iter.next(blki,blkj)) {
    const bool issingle=iter.A0_matrix(A0s,blki,blkj);
    iter.addsym(Htens,blki,blkj,0,0,issingle,Type2Type<double>());
  }
}

template<class T> void
CrystalOpGenerator::add_A2(BlockedMatrix<complex>& Hmats, const Matrix<T>& couplings) const
{
  CrystalSystem_iterator iter(*this);
  iter.ensure_null(Hmats);
  size_t blki,blkj;

  while (iter.next(blki,blkj)) {
    const bool issingle=iter.A2_matrix(couplings,blki,blkj);
    iter.addsym(Hmats,blki,blkj,Type2Type<double>(),issingle);
  }  
}

void
CrystalOpGenerator::add_A0(BlockedMatrix<complex>& Hmats, const rmatrix& couplings) const
{
  CrystalSystem_iterator iter(*this);
  iter.ensure_null(Hmats);
  size_t blki,blkj;

  while (iter.next(blki,blkj)) {
    const bool issingle=iter.A0_matrix(couplings,blki,blkj);
    iter.addsym(Hmats,blki,blkj,Type2Type<double>(),issingle);
  }  
}

  //create symmetrised shift Hamiltonian, given shifts and complete set of symmetrised z operators
  ListList<double> diag_Hcs(const CrystalOpGenerator& sys, const BaseList<double>& shifts)
  {
    ListList<double> Hcs(mxflag::temporary);
    const size_t M=sys.nspins_cell();
    if (M!=shifts.size())
      throw Mismatch("diag_Hcs: size of shift vector");
        
    for (size_t m=0;m<M;m++)
      mla(Hcs,shifts(m),sys.diag_Iz(m));
    
    return Hcs;
  }

static void multiply_full(cmatrix& dest, const complex& v, const BaseList<double>& H)
{
  size_t n=H.size();
  dest.create(n,n);
  dest=complex(0.0);
  for (;n--;)
    mla(dest(n,n),H(n),v);
}

static void mla(Tensor<cmatrix>& d, const Tensor<complex>& A, const BaseList<double>& H)
{
  for (int l=A.rank();l>=0;l--) {
    if (A.have_rank(l)) { //rank present?
      if (!d.have_rank(l)) {
	d.ensure_rank(l);
	for (int m=-l;m<=l;m++) 
	  multiply_full(d(l,m),A(l,m),H);
      }
      else {
	for (int m=-l;m<=l;m++) 
	  mla(d(l,m),A(l,m),H);
      }
    }
  }
}

  void
  CrystalOpGenerator::add_Hcs(BlockedMatrixTensor<complex>& Htens, const space_T& CSA, size_t m) const
{
  if (Htens.size()!=diagstr_.size())
    throw Mismatch("add_Hcs");

  for (size_t k=Htens.size();k--;) //loop over eigenvalues
    switch (diagstr_(k)) {
    case 0: break;
    case 1: {
      space_T& curHisol=Htens(k).element();
      ::libcmatrix::mla(curHisol,diag_Iz(m)(k).front(),CSA);
    }
      break;
    default: {
      Tensor<cmatrix>& curHtens=Htens(k).matrix();
      ::libcmatrix::mla(curHtens,CSA,diag_Iz(m)(k));
    }
    }
}

//special case for single state block (k=0 only)
// void add_Hcs(space_T& Hisol, const CrystalOpGenerator& sys, const BaseList<space_T>& CSAs)
// {
//   (void)sys.state(); //ensure single state
//   const size_t M=sys.nspins();
//   if (M!=CSAs.length())
//     throw Mismatch("add_Hcs");
  
//   for (size_t m=M;m--;)
//     mla(Hisol,sys.diag_Iz(m)(0U,0U),CSAs(m));
// }

void CrystalOpGenerator::add_Hcs(cmatrix& dest, double shift, size_t m, int k) const
{
  ::libcmatrix::mla(dest,shift,diag_Iz(m)(k));
}

 template<class T1,class T2,class T3> void multiply(BlockedMatrix<T1>& d, const T2& v, const ListList<T3>& a)
   {
     LCM_STATIC_CHECK( LCM_DIM(T2)==0, BlockedMatrix_multiply);
     const size_t n=a.size();
     ScratchList<size_t> sizes(n);
     size_t i;
     for (i=n;i--;)
       sizes(i)=a.size(i);
     d.create(sizes);
     for (i=n;i--;)
       multiply(d(i),v,a(i));
   }


void CrystalOpGenerator::add_Hcs(BlockedMatrix<complex>& dest, double shift, size_t m) const
{

  ::libcmatrix::mla(dest,shift,diag_Iz(m));
}

void CrystalOpGenerator::add_Hcs(ListList<double>& dest, double shift, size_t m) const
{
  ::libcmatrix::mla(dest,shift,diag_Iz(m));
}

void SimpleSymmetrise::addtoevals(BaseList<complex>& evals, size_t blki, size_t blkj, size_t r, size_t c) const
{
  const size_t blkisize(linkedstates_.size(blki));
  const size_t blkjsize(linkedstates_.size(blkj));
  const double scale=scalefacs(blkisize*blkjsize); //normalisation factor
  
  evals(0)+=scale; //eigenvalue 0
  
  const int mstep=ncells_+c-r;

  for (size_t m=1;m<useeigs;m++) {
    if (haseig(blkisize,m) && haseig(blkjsize,m)) { //if eigenvalue present
      const size_t wh=(m*mstep) % ncells_;
      mla(evals(m),scale,eigfacs(wh));
    }
  }
}

  //Symmetrised Ip
bool CrystalSystem_iterator::Ipevalues(DynamicList<complex>& evals, size_t blki, size_t blkj, int spini)
{
  const BaseList<state_t> cstates(symmetriser.linkedstates()(blkj));
  const size_t total(sym.nspins());

  if (verbose>2) {
    const BaseList<state_t> rstates(symmetriser.linkedstates()(blki));
    std::cout << "Col states: "; printbin(std::cout,cstates,total); std::cout << '\n';
    std::cout << "Row states: "; printbin(std::cout,rstates,total); std::cout << '\n';
  }

  const state_t flipmask=~maskelement(total,spini);
  evals=0.0; //clear result

  bool haveany=false;
  for (size_t c=cstates.size();c--;) {
    const state_t nstate=cstates(c) & flipmask; // if bit is already clear, nstate cannot be in blki
    if (state_to_block(nstate)==blki) {
      symmetriser.addtoevals(evals,blki,blkj,state_to_index(nstate),c);
      haveany=true;
    }
  }
  return haveany;
}

bool CrystalSystem_iterator::Fpevalues(DynamicList<complex>& evals,size_t blki,size_t blkj,size_t nuc)
{
  const BaseList<state_t> cstates(symmetriser.linkedstates()(blkj));
  const size_t total(sym.nspins());

  if (verbose>2) {
    const BaseList<state_t> rstates(symmetriser.linkedstates()(blki));
    std::cout << "Col states: "; printbin(std::cout,cstates,total); std::cout << '\n'; 
    std::cout << "Row states: "; printbin(std::cout,rstates,total); std::cout << '\n'; 
  }

  bool haveany=false;
  evals=0.0;

  const size_t M(sym.nspins_cell());

  for (size_t c=cstates.size();c--;) { //loop over ket states
    const state_t ostate=cstates(c);
    for (size_t spini=0;spini<M;spini++) { //loop over unit cell
      if (nuc==NULL_NUCLEUS || sym(spini).nucleus()==nuc) { //considering this spin?
	const state_t flipmask=~maskelement(total,spini);
	const state_t nstate=ostate & flipmask; // if bit is clear, nstate cannot be in blki
	if (state_to_block(nstate)==blki) {
	  symmetriser.addtoevals(evals,blki,blkj,state_to_index(nstate),c);
	  haveany=true;
	}
      }
    }
  }
  return haveany;
}

  //accumulate symmetrised elements into final matrix set
void
CrystalSystem_iterator::add(BlockedMatrix<complex>& dest, const BaseList<complex>& evals, size_t blki, size_t blkj) const
{
  BaseList<int> reigptrs(sym.eigptrs.row(blki));
  BaseList<int> ceigptrs(sym.eigptrs.row(blkj));

  for (size_t eval=useeigs;eval--;) {
    int dr=reigptrs(eval);
    int dc=ceigptrs(eval);
    if ((dr>=0) && (dc>=0)) {//eigenvalue exists?
      const size_t ind(get_block(eval));
      dest(ind)(size_t(dr),size_t(dc))+=evals(eval);
    }
  }
}

  //ditto but for Hermitian operator where iterated over upper diagonal only
void
CrystalSystem_iterator::add_hermitian(BlockedMatrix<complex>& dest, const BaseList<complex>& evals, size_t blki, size_t blkj) const
{
  BaseList<int> reigptrs(sym.eigptrs.row(blki));
  BaseList<int> ceigptrs(sym.eigptrs.row(blkj));

  for (size_t eval=useeigs;eval--;) {
    int dr=reigptrs(eval);
    int dc=ceigptrs(eval);
    if ((dr>=0) && (dc>=0)) {
      assert(dr!=dc);
      const size_t ind(get_block(eval));
      const complex& v=evals(eval);
      dest(ind)(size_t(dr),size_t(dc))+=v;
      dest(ind)(size_t(dc),size_t(dr))+=conj(v);
    }
  }
}

  //as above, but for single eigenvalue
void
CrystalSystem_iterator::add(cmatrix& dest, complex v, size_t blki, size_t blkj, int eval) const
{
  const int dr=sym.eigptrs(blki,eval);
  const int dc=sym.eigptrs(blkj,eval);
  if ((dr>=0) && (dc>=0))
    dest(size_t(dr),size_t(dc))+=v;
}

void
CrystalSystem_iterator::add_hermitian(cmatrix& dest, complex v, size_t blki, size_t blkj, int eval) const
{
  const int dr=sym.eigptrs(blki,eval);
  const int dc=sym.eigptrs(blkj,eval);
  if ((dr>=0) && (dc>=0)) {
    dest(size_t(dr),size_t(dc))+=v;
    dest(size_t(dc),size_t(dr))+=conj(v);
  }
}

void CrystalOpGenerator::mla_Iz(ListList<double>& dest, double scale, const operator_spec& spec) const
{
  if (spec.op!='z')
    throw InvalidParameter("mla_Iz: only valid for z operators");
  if (spec.nuc==NULL_NUCLEUS)
    mla_Iz(dest,scale,spec.number);
  else
    mla_Fz(dest,scale,spec.nuc);
}

  //ditto, all eigenvalues
void CrystalOpGenerator::mla(BlockedMatrix<complex>& dest, block_pattern& blkspec, double scale, const operator_spec& opspec) const
{
  if (opspec.op=='z') { //ignore 'c' as not valid
    ListList<double> destz;
    mla_Iz(destz,scale,opspec);
    dest+=destz;
    blkspec=block_pattern(*this);
    return;
  }

  blkspec=block_pattern(*this,opspec);

  if (!dest) {
    create(dest,blkspec);
    dest=complex(0.0);
  }
  else {
    if (dest.length()!=totalblocks())
      throw Mismatch("rawmla");
  }

  size_t nuc(opspec.nuc);
  if (nuc==NULL_NUCLEUS)
    nuc=(*sysp_)(opspec.number).nucleus();
  CrystalSystem_iterator iter(*this,blkspec,nuc,1);
  iter.mla(dest,scale,opspec);
}

namespace {
  const char PUREERR[]="productoperator_spec: scale factor must be pure real or pure imaginary";
}

void CrystalOpGenerator::add(BlockedMatrix<complex>& dest, block_pattern& blkspec, const productoperator_spec& spec) const 
{ 
  static const complex iconst(0.0,1.0);
  for (size_t i=spec.size();i--;) {
    const BaseList<operator_spec> opspec(spec.specs(i));
    if (opspec.size()!=1)
      throw Failed("Can't use multi-spin productoperators with CrystalSystem");

    const complex& cscale(spec.scales(i));
    if (imag(cscale)) {
      if (real(cscale))
	throw Failed(PUREERR);
      BlockedMatrix<complex> tmp;
      mla(tmp,blkspec,imag(cscale),opspec.front());
      ::libcmatrix::mla(dest,iconst,tmp);
    }
    else
      mla(dest,blkspec,real(cscale),opspec.front());
  }
}

void CrystalSystem_iterator::mla(BlockedMatrix<complex>& dest, double scale, const operator_spec& opspec)
{
  size_t blki,blkj;
  const size_t nuc=opspec.nucleus(sym.spinsystem());
  const bool notherm(sym.isblocked(nuc)); //if blocked by coherence, don't construct both sides of hermitian operators

  const bool isI=(opspec.nuc==NULL_NUCLEUS);

  while (next(blki,blkj)) {

    const bool haveany = isI 
      ? Ipevalues(evalues,blki,blkj,opspec.number)
      : Fpevalues(evalues,blki,blkj,opspec.nuc);
    
    if (verbose>1) {
      std::cout << "I/Fp" << blki << "," << blkj << ": ";
      if (haveany)
	std::cout << evalues << '\n';
      else
	std::cout << "<all zero>\n";
    }
    if (!haveany)
      continue;

    switch (opspec.op) {
    case 'x':
      evalues*=0.5*scale;
      if (notherm)
	add(dest,evalues,blki,blkj);
      else
	add_hermitian(dest,evalues,blki,blkj);
      break;

    case 'y':
      evalues*=complex(0,-0.5*scale);      
      if (notherm)
	add(dest,evalues,blki,blkj);
      else
	add_hermitian(dest,evalues,blki,blkj);
      break;

    case '+':
      if (scale!=1.0) evalues*=scale;
      add(dest,evalues,blki,blkj);
      break;

    case '-':
      if (scale!=1.0) evalues*=scale;
      conj_ip(evalues);
      add(dest,evalues,blkj,blki);
      break;

    default:
      throw InvalidParameter("Operator must be x, y, + , or -");
    }
  }
}

void sum_permuted(rmatrix& dest,const rmatrix& source,const BaseList<size_t>& next, int N)
{
  dest=source;
  rmatrix tmp;

  for (size_t i=1;i<N;i++) {
    tmp=dest(next,next);
    dest=tmp;
    dest+=source;
  }
}

void sum_permuted(List<double>& dest,const BaseList<double>& source,const BaseList<size_t>& next, int N)
{
  dest=source;
  List<double> tmp;

  for (int i=1;i<N;i++) {
    tmp=dest(next);
    dest=tmp;
    tmp+=source;
  }
}

//Specialise default add_Hamiltonian
template<class OutType, class StoreType> void add_coupling_Hamiltonians(OutType& dest, const CrystalOpGenerator& opgen, const HamiltonianStore<StoreType>& Hstore) 
{
  Matrix<double> A0s;
  Matrix<StoreType> A2s;

  Hstore.split_couplings(A0s,A2s);
  if (!!A0s) {
    if (opgen.verbose()>1)
      std::cout << "A0 couplings\n" << A0s;
    opgen.add_A0(dest,A0s);
    if (opgen.verbose()>1)
      std::cout << "After A0 add\n" << dest;
  }
  if (opgen.verbose()>1)
    std::cout << "A2 couplings\n" << A2s;
  opgen.add_A2(dest,A2s);
}

template<> void BlockedSpinningHamiltonianImpl<complex,CrystalOpGenerator>::interactions(const HamiltonianStore<space_T>& store)
{
  const size_t verbose(opgen_.verbose());
  BlockedMatrixTensor<complex> tmp;

  add_Hamiltonian(tmp,opgen_,store);
  size_t n(tmp.size());
  List<SpinningHamiltonian>& H(*this);
  if (H.size()!=n) {
    if (H.empty())
      H.create(n,SpinningHamiltonian(rotor_speed_,rotor_phase_,rinfo_));
    else
      throw Mismatch("BlockedSpinningHamiltonian::interactions");
  }
  for (;n--;) {
    MatrixTensor<complex>& source(tmp(n));
    if (verbose>1) 
      std::cout << "MatrixTensor(" << n << ")\n" << source  << '\n';
    SpinningHamiltonian& Hspin(H(n));
    switch (source.rows()) {
    case 0: break;
    case 1:
      Hspin.tensor(source.element());
      if (!!Hbase_)
	Hspin+=Hbase_(n);
      break;
    default:
      Hspin.tensor(source.matrix());
      source.clear();
      if (!!Hbase_)
	Hspin+=Hbase_(n);
      break;
    }
  }    
}

template<> struct Ham_traits< BlockedMatrixTensor<complex> > {
  typedef space_T coupling_type;
};

LCM_INSTANTIATE_TYPED_HAMILTONIANS(CrystalOpGenerator,complex)
LCM_INSTANTIATE_SPINOPGEN_CREATE(CrystalOpGenerator)

template void add_Hamiltonian<BlockedMatrix<complex>, CrystalOpGenerator, Tensor<complex> >(BlockedMatrix<complex>&, CrystalOpGenerator const&, HamiltonianStore< Tensor<libcmatrix::complex> > const&);

//template<> BlockedFilter::BlockedFilter(const CrystalOpGenerator&, const BaseList<int>&);

} //namespace libcmatrix

// static void translate_raw(rmatrix &dest,const rmatrix &source,const BaseList<size_t> &TT)
// {
//   const int n=TT.length();
//   if (!issquare(source) || n!=source.rows()) throw Mismatch("translate_raw");

//   dest.create(n,n);
//   for (int i=n;i--;)
//     for (int j=n;j--;)
//       dest(TT(i),TT(j))=source(i,j);
// }

// void Tmatrix(cmatrix& TV, const ListList<state_t>& TT, int total, const BaseList<size_t>& which)
// {
//   int j;

//   const bool havesel=(which.length()!=0);
//   const size_t blks=havesel ? which.length() : TT.blocks();

//   const long fulldim=long(1)<<total;

//   long totdim=0; //size of output matrix
//   if (havesel) { 
//     for (j=blks;j--;) //sum over selected blocks
//       totdim+=TT.size(which(j));
//   }
//   else
//     totdim=fulldim;

//   TV.create(totdim,totdim,0.0); //construct matrix to be filled

//   int point=0;
//   for (j=blks;j--;) { //loop over blocks
//     const BaseList<state_t> clist=TT(havesel ? which(j) : j); //get state list
//     const size_t size=clist.length();
//     const double scale=1.0/sqrt((double)size); //normalisation factor
    
//     for (size_t n=size;n--;) {
//       for (size_t m=size;m--;)
// 	TV(clist(m),point)=expi( (2*M_PI*m*n)/size)*scale;
//       point++;
//     }
//   }
// }

  //create tensor of correct size (spinning) for eigenvalue, k, default value, val
// void
// CrystalOpGenerator::clear(Tensor<cmatrix>& ctens, int k, double val) const
// {
//   if (ismatrix(k)) {
//     const size_t rows=rblkstr(k);
//     const size_t cols=cblkstr(k);
    
//     ctens.create(2,mxflag::maximum);
    
//     for (int m=-2;m<=2;m++) // create (2,0) matrix, even if unused
//       ctens(2,m).create(rows,cols,val);
//   }
//   else //if eigenvalue doesn't exist, kill off 
//     ctens.clear();
// }

//   //space_T used to store "single element" H's in spinning case
// void
// CrystalOpGenerator::clear(space_T& Hisol, int k, double val) const
// {
//   if (!ismatrix(k)) {
//     Hisol.create(2,mxflag::maximum);
//     Hisol=complex(val);
//   }
//   else
//     Hisol.clear();
// }
    
//   //create tensors of correct size (spinning)  
// void
// CrystalOpGenerator::clear(List< Tensor<cmatrix> >& Htens, BaseList<space_T> Hisol, double val, double valisol) const
// {
//   const bool use_isol=(Hisol.length()!=0);

//   Htens.create(useeigs);
//   for (size_t j=useeigs;j--;) {
//     if (ismatrix(j))
//       clear(Htens(j),j,val);
//     else {
//       if (use_isol)
// 	clear(Hisol(j),j,valisol);
//     }
//   }
// }

// void
// CrystalOpGenerator::clear(BaseList<cmatrix> Hmats, BaseList<double> Hisol, double val, double valisol) const
// {
//   if (Hmats.size()!=useeigs)
//     throw Mismatch("CrystalOpGenerator::clear");

//   const size_t n=Hisol.size();
//   if (n) {
//     if (n!=useeigs) 
//       throw Mismatch("CrystalOpGenerator::clear");
//     Hisol=valisol;
//   }

//   for (size_t j=useeigs;j--;) {
//     if (ismatrix(j))
//       clear(Hmats(j),j,val);
//     else
//       Hmats(j).clear();
//   }
// }

//void
// CrystalOpGenerator::clear(BaseList<cmatrix> Hmats, double val) const
// {
//   if (Hmats.size()!=useeigs)
//     throw Mismatch("CrystalOpGenerator::clear");

//   for (size_t j=useeigs;j--;)
//     clear(Hmats(j),j,val);
// }

// void
// CrystalOpGenerator::clear(List<cmatrix>& Hmats, List<double>& Hisol, double val, double valisol) const
// {
//   Hmats.create(useeigs);
//   Hisol.create(useeigs);
//   clear(static_cast<BaseList<cmatrix>& >(Hmats),Hisol,val,valisol);
// }

// void
// CrystalOpGenerator::clear(List<cmatrix>& Hmats, BaseList<double> Hisol, double val, double valisol) const
// {
//   Hmats.create(useeigs);
//   clear(static_cast<BaseList<cmatrix>& >(Hmats),Hisol,val,valisol);
// }

// void
// CrystalOpGenerator::clear(BlockedMatrix<complex>& Hmats, double val) const
// {
//   if (isdiagonal())
//     Hmats.create(rblkstr,val);
//   else
//     Hmats.create(rblkstr,cblkstr,val);
// }

// void
// CrystalOpGenerator::clear(cmatrix& d, int k, double val) const
// {
//   if (exists(k))
//     d.create(rows(k),cols(k),val);
//   else
//     d.clear();
// }
  
  //add symmetrised block (blki,blkj) eigenvalue eval to (static) Hamiltonian
// void
// CrystalOpGenerator::place(cmatrix& Hk, complex& Hisol, size_t blki, size_t blkj, int eval) const
// {
//   const int dr=sym.eigptrs(blki,eval);
//   const int dc=sym.eigptrs(blkj,eval);
//   const size_t rs=sym.T.size(blki);
//   const size_t cs=sym.T.size(blkj);

//   if (eval==0) {
//     if (ismatrix(0U))
//       Hk(dr,dc)=csymmetrise0(rs,cs);
//     else
//       Hisol=csymmetrise0(rs,cs);
//   }
//   else {
//     if ((dr>=0) && (dc>=0)) {
//       if (ismatrix(eval))
// 	Hk(dr,dc)=csymmetrise(rs,cs,eval);
//       else {
// 	Hisol=csymmetrise(rs,cs,eval);
//       }
//     }
//   }
// }
 
//calculate matrix that symmetrises w.r.t translation (not useful for large systems!)
// void 
// CrystalOpGenerator::get_diag(cmatrix& D) const 
// {
//   List<size_t> rind;
//   makerindex(rind,*sysp_); //create state->index reverse index

//   D.create(sysp_->size(),sysp_->size(),0.0); //empty matrix

//   CrystalSystem_diagiterator iter(*this);
//   const ListList<state_t>& linkedstates(symmetriserp_->linkedstates());

//   size_t blk;
//   size_t point=0;
//   while (iter.next(blk)) { //loop over state sets
//     const BaseList<state_t> states(linkedstates(blk)); //state list
//     const size_t blksize=states.size();
//     const double scalef=1.0/sqrt(double(blksize)); //normalisation factor
//     for (size_t j=0;j<blksize;j++) {
//       for (size_t k=0;k<blksize;k++)
// 	D(rind(states(k)),point)=eigfacs((j*k) % N_)*scalef;
//       point++;
//     }
//   }
//   if (point!=D.cols()) //sanity check - have we filled D?
//     throw Failed("get_diag");
// }

//As above but restricted to single eigenvalue (usek)
// void
// CrystalOpGenerator::get_diag(cmatrix& D,int usek) const 
// {
//   List<size_t> rind;
//   makerindex(rind,sys);

//   D.create(rows(),diag_str_(usek),0.0);

//   CrystalSystem_diagiterator iter(*this);

//   size_t blk;
//   size_t point=0;
//   while (iter.next(blk)) {
//     const BaseList<state_t> states(linkedstates_(blk));
//     const size_t blksize=states.size();
//     if (sym.haseig(blksize,usek)) {
//       const double scalef=1.0/sqrt(double(blksize));
//       for (size_t k=0;k<blksize;k++)
// 	D(rind(states(k)),point)=sym.eigfacs( (usek*k) % N)*scalef;
//       point++;
//     }
//   }
//   if (point!=D.cols())
//     throw Failed();
// }

// void
// CrystalOpGenerator::create(cmatrix& a, int k) const 
// {
//   const size_t rs=rblkstr(k);
//   if (isdiagonal())
//     a.create(rs,rs);
//   else
//     a.create(rs,cblkstr(k));
// }

// void
// CrystalOpGenerator::create(MatrixTensor<complex>& a, int k) const 
// {
//   if (!isdiagonal())
//     throw Failed("Can only create MatrixTensor objects for diagonal CrystalOpGenerator");
//   a.create(rblkstr(k),2,mxflag::maximum);
// }

// void
// CrystalOpGenerator::create(BlockedMatrix<complex>& Hmats) const
// {
//   if (isdiagonal())
//     Hmats.create(rblkstr);
//   else
//     Hmats.create(rblkstr,cblkstr);
// }

// void
// CrystalOpGenerator::create(BlockedMatrixTensor<complex>& Htens) const
// {
//   if (!isdiagonal())
//     throw Failed("Can only create MatrixTensor objects for diagonal CrystalOpGenerator");
//   Htens.create(rblkstr,2,mxflag::maximum);
// }

// void
// CrystalOpGenerator::add_Hcoupling(complex& H, P_GENSPINNINGMATRIX func, const Matrix<space_T>& tensors, int m) const
// {
//   if (firstrblk<0 || !(this->*func)(tensors,firstrblk,firstrblk,m))
//     throw Failed("sym_Hdipolar_element: not single state block!");
//   H+=symmetrise0(Type2Type<complex>());
// }

  //as above but for single eigenvalue
// template<class Type> inline void CrystalOpGenerator::addtoeval(Type& s, size_t blkisize, size_t blkjsize, size_t r, size_t c, int m)
// {
//   if (haseig(blkisize,m) && haseig(blkjsize,m)) {
//     const int wh=(m*(N+c-r)) % N;
//     mla(s,scalefacs(blkisize,blkjsize),eigfacs(wh));
//   }
// }
