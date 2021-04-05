/* Code for creating subblocks of (product) spin operators for spin-1/2 only systems */

#undef LCM_SUPPRESS_VIEWS
#include <ctype.h>
#include "spinhalf_system.h"
#include "rmatrix.h"
#include "ScratchList.h"

namespace libcmatrix {

  void spinhalf_system::nucleus(size_t n, const spin& nuc)
  {
    if (nuc.deg()>2)
      throw Failed("spinhalf_system::nucleus: spin-1/2 only allowed");
    else
      List<spin>::operator()(n)=nuc;
  }

  void spinhalf_system::docreate(size_t n)
  {
    iscomplete_=true;
    _nspins=n;
    cstates.create(1<<n);
    for (size_t i=cols();i--;)
      cstates(i)=i;
    _rows=cols();
    ensurereverse();
  }

  spinhalf_system::spinhalf_system(int n,const char* label,float_t mz)
    : basespin_system(n,spin(label))
{
  _nspins=n;
  cstates=::libcmatrix::mzstates(n,mz);
  _rows=cols();
  ensurereverse();
}

spinhalf_system::spinhalf_system(int n,const char* label,float_t rmz,float_t cmz)
  : basespin_system(n,spin(label))
{
  _nspins=n;
  cstates=::libcmatrix::mzstates(n,cmz);
  rstates=::libcmatrix::mzstates(n,rmz);
  ensurereverse();
}

spinhalf_system::spinhalf_system(const basespin_system& sys, size_t rep)
  : basespin_system(sys,rep)
{
  if (!sys.isspinhalfonly())
    throw Failed("spinhalf_system: can only be constructed from pure spin-1/2 systems");
  docreate(nspins());  
}

void spinhalf_system::docreate(size_t n,const BaseList<state_t> &inds)
{
  iscomplete_=false;
  _nspins=n;
  if (max(inds)>=(state_t(1)<<n))
    throw BadIndex("spinhalf_system");
  cstates=inds;
  _rows=cols();
  ensurereverse();
}

void spinhalf_system::docreate(size_t n,const BaseList<state_t> &rinds,const BaseList<state_t> &cinds)
{
  iscomplete_=false;
  _nspins=n;
  const state_t maxdeg=state_t(1)<<n;
  if ( (max(rinds)>=maxdeg) || (max(cinds)>=maxdeg) )
    throw BadIndex("spinhalf_system::create");
  cstates=cinds;
  rstates=rinds;
  _rows=rstates.length();
  ensurereverse();
}

void spinhalf_system::ensurematrix(cmatrix &dest) const
{
  if (isdefined(dest)) {
    if (dest.rows()!=rows() || dest.cols()!=cols())
      throw Mismatch("spinhalf_system::create");
  }
  else {
    dest.create(rows(),cols());
    dest=0.0;
  }
}

void spinhalf_system::ensurereverse()
{
  reverse.create(1 << _nspins);
  reverse=-1;
  const List<state_t> &states = isdiagonal() ? cstates : rstates;
  for (size_t i=states.length();i--;)
    reverse(states(i))=i;
}

state_t spinhalf_system::apply_permutation(state_t orig, const Permutation& permvec)
{
  state_t state(0);
  const size_t nspin=permvec.size();
  for (size_t i=nspin;i--;) {
    state|= tostate(orig & 1,nspin,permvec(i));
    orig>>=1;
  }
  return state;
}

void spinhalf_system::permutation_vectorH(BaseList<state_t>& ninds, const Permutation& permvec) const
{
  validate_diagonal();
  const size_t n=rows();
  if (ninds.length()!=n)
    throw Mismatch("permutation_vectorH");
  for (size_t i=n;i--;) {
    const int nstate=reverse(rawpermute(cstates(i),permvec));
    if (nstate<0)
      throw Failed("permutation_matrixH: permuted state outside sub-space");
    else
      ninds(i)=nstate;
  }
}

/* Iz etc. only work for diagonal blocks */

void spinhalf_system::rawmla_Iz(BaseList<double> &dest,double scale,size_t n) const
{
  const state_t umask=mask(n);
  if (!isdiagonal())
    throw Failed("mla_Iz: diagonal blocks only");
  if (dest.length()!=cols())
    throw Mismatch("mla_Iz");
  scale*=0.5;
  for (size_t i=cols();i--;) {
    if (cstates(i) & umask)
      dest(i)-=scale;
    else
      dest(i)+=scale;
  }
}

void spinhalf_system::rawmla_Iz(cmatrix &dest,double scale,state_t umask) const
{
  scale*=0.5;
  if (isdiagonal()) {
    for (size_t i=cols();i--;) {
      if ( cstates(i) & umask)
	dest(i,i)-=scale;
      else
	dest(i,i)+=scale;
    }
  }
  else {
    for (size_t i=cols();i--;) {
      const state_t cstate=cstates(i);
      const int nindex=reverse(cstate);
      if (nindex>=0) {
	if ( cstate & umask)
	  dest(nindex,i)-=scale;
	else
	  dest(nindex,i)+=scale;
      }
    }
  }
}

void spinhalf_system::mla_Ip(cmatrix &dest,double scale,state_t umask) const
{
  for (size_t i=cols();i--;) {
    state_t cstate=cstates(i);
    if (cstate & umask) {
      const int nwhich=reverse( cstate - umask);
      if (nwhich>=0)
	dest(nwhich,i)+=scale;
    }
  }
}

void spinhalf_system::mla_Im(cmatrix &dest,double scale,state_t umask) const
{
  for (size_t i=cols();i--;) {
    state_t cstate=cstates(i);
    if ((cstate & umask)==0) {
      const int nwhich=reverse( cstate | umask);
      if (nwhich>=0) dest(nwhich,i)+=scale;
    }
  }
}

void spinhalf_system::print(std::ostream& ostr) const
{
  ostr << *this;
}

void spinhalf_system::mla_I(cmatrix &dest,double scale,size_t n,char op) const
{
  ensurematrix(dest);
  const state_t umask=mask(n);

  int i;

  switch (op) {
  case '+':
    mla_Ip(dest,scale,umask);
    return;
  case '-':
    mla_Im(dest,scale,umask);
    return;
  case 'x':
    mla_Ix(dest,scale,umask);
    return;
  case 'y':
    mla_Iy(dest,scale,umask);
    return;
  case 'z':
    rawmla_Iz(dest,scale,umask);
    return;
  case 'I':
    if (isdiagonal()) {
      for (i=cols();i--;)
	dest(i,i)+=scale;
    }
    else {
      for (i=cols();i--;) {
	const int nindex=reverse(cstates(i));
	if (nindex>=0)
	  dest(nindex,i)+=scale;
      }
    }
    return;
  case 'a': case 'b':
    const state_t which=(op=='a') ? umask : 0;

    if (isdiagonal()) {
      for (i=cols();i--;) {
	if ( (cstates(i) & umask)==which)
	  dest(i,i)+=scale;
      }
    }
    else {
      for (i=cols();i--;) {
	if ( (cstates(i) & umask)==which) {
	  const int nindex=reverse(cstates(i));
	  if (nindex>=0)
	    dest(nindex,i)+=scale;
	}
      }
    }
    return;
  }
  throw InvalidParameter("mla_I");
}

void spinhalf_system::mla_Ix(cmatrix &dest,double scale,state_t umask) const
{
  scale*=0.5;

  for (size_t i=cols();i--;) {
    const state_t nstate=cstates(i) ^ umask;
    const int nwhich=reverse(nstate);
    if (nwhich>=0)
      dest(nwhich,i)+=scale;
  }
}

void spinhalf_system::mla_Iy(cmatrix &dest,double scale,state_t umask) const
{
  const complex val(0.0,0.5*scale);

  for (size_t i=cols();i--;) {
    const state_t nstate=cstates(i) ^ umask;
    const int nwhich=reverse(nstate);
    if (nwhich>=0) {
      if (nstate & umask)
	dest(nwhich,i)+=val;
      else
	dest(nwhich,i)-=val;
    }
  }
}

void spinhalf_system::rawmla_Iz(BaseList<double> &dest,double scale,size_t j,size_t k) const
{
  if (!isdiagonal()) 
    throw Failed("mla_Iz: diagonal blocks only");
  if (dest.length()!=cols())
    throw Mismatch("mla_Iz");
  scale*=0.25;
  const state_t maskj=mask(j);
  const state_t maskk=mask(k);
  for (size_t i=cols();i--;) {
    const bool issetj=(( cstates(i) & maskj)!=0);
    const bool issetk=(( cstates(i) & maskk)!=0);
    if (issetj ^ issetk)
      dest(i)-=scale;
    else
      dest(i)+=scale;
  }
}

static const complex ihalf(0,0.5);
static const complex mihalf(0,-0.5);

static bool advance(state_t &cstate,complex &val,state_t umask,char cop)
{
  switch (cop) {
  case 'x':
    cstate ^=umask;
    val*=0.5;
    return true;
  case 'y':
    val*=(cstate & umask) ? mihalf : ihalf;
    cstate ^=umask;
    return true;
  case 'z':
    val*=(cstate & umask) ? -0.5 : 0.5;
    return true;
  case '-':
    if (cstate & umask) 
      return false;
    cstate|=umask;
    return true;
  case '+':
    if (cstate & umask) {
      cstate-=umask;
      return true;
    }
    return false;
  case 'a':
    return (cstate & umask) ? false : true;
  case 'b':
    return (cstate & umask) ? true : false;
  }
  throw InvalidParameter("advance");
}

void spinhalf_system::mla_I(cmatrix &dest,double scale,size_t n,char opn,size_t m,char opm) const
{
  if (n>=_nspins || m>=_nspins)
    throw BadIndex("mla_I");

  ensurematrix(dest);

  const state_t maskn=mask(n);
  const state_t maskm=mask(m);
  
  for (size_t i=cols();i--;) {
    state_t cstate=cstates(i);
    complex val(scale);
    if (advance(cstate,val,maskn,opn) && advance(cstate,val,maskm,opm)) {
      const int nwhich=reverse(cstate);
      if (nwhich>=0)
	dest(nwhich,i)+=val;
    }
  }
}

void spinhalf_system::mla_I(cmatrix &dest,double scale,const BaseList<char> &ops) const
{
  if (_nspins!=ops.length())
    throw Mismatch("mla_I");

  ensurematrix(dest);
  
  for (size_t i=cols();i--;) {
    state_t cstate=cstates(i);
    complex val(scale);
    for (size_t j=_nspins;j--;) {
      if (!advance(cstate,val,mask(j),ops(j))) {
	val=0.0;
	break;
      }
    }
    const int nwhich=reverse(cstate);
    if (nwhich>=0)
      dest(nwhich,i)+=val;
  }
}

void printbin(std::ostream& ostr,state_t state,int bits)
{
  state_t mask=state_t(1)<< (bits-1);
  for (;mask;mask>>=1)
    ostr << ((state & mask) ? '1' : '0');
}

void printbin(std::ostream& ostr,const BaseList<state_t> &list,int bits)
{
  ostr << "[ ";
  const size_t n=list.length();
  for (size_t i=0;i<n;i++) {
    printbin(ostr,list(i),bits);
    ostr << ((i==n-1) ? " ]" : ", ");
  }
}

std::ostream& operator << (std::ostream& ostr,const spinhalf_system& sys)
{
  if (!(ostr.flags() & std::ios::hex))
    //ostr << "Spins: " << sys.nspins() << "\n";
    ostr << static_cast<const basespin_system&>(sys) << '\n';
  const bool isdiag=sys.isdiagonal();
  ostr << (isdiag ? "Diagonal states: " : "Row states: ");
  printbin(ostr,sys.brastates(),sys.nspins());
  if (!isdiag) {
    ostr << "\nColumn states: ";
    printbin(ostr,sys.ketstates(),sys.nspins());
  }
  return ostr << std::endl;
}

bool isdiagop(char op) /* identity doesn't count! */
{
  switch (tolower(op)) {
  case 'z': case 'a': case 'b':
    return true;
  }
  return false;
}

double setval(char op)
{
  switch (op) {
  case 'z':
    return -0.5;
  case 'a':
    return 0.0;
  case 'b':
    return 1.0;
  }
  throw InvalidParameter("setval");
}

double notsetval(char op)
{
  switch (op) {
  case 'z':
    return 0.5;
  case 'a':
    return 1.0;
  case 'b':
    return 0.0;
  }
  throw InvalidParameter("notsetval");
}

complex Ielement(int nspins,state_t bra,state_t ket,size_t n,char op)
{
  complex val(1.0);
  if (!advance(ket,val,maskelement(nspins,n),op))
    return complex(0.0);
  return (bra==ket) ? val : complex(0.0);
}
  
complex Ielement(int nspins,state_t bra,state_t ket,size_t n,char nop,size_t m,char mop)
{
  complex val(1.0);
  if (!advance(ket,val,maskelement(nspins,m),mop))
    return complex(0.0);
  if (!advance(ket,val,maskelement(nspins,n),nop))
    return complex(0.0);
  return (bra==ket) ? val : complex(0.0);
}

complex Ielement(int nspins,state_t bra,state_t ket,const BaseList<char> &ops)
{
  if (ops.length()!=nspins)
    throw Mismatch("Ielement");
  complex val(1.0);
  for (size_t i=nspins;i--;) {
    if (!advance(ket,val,maskelement(nspins,i),ops(i)))
      return complex(0.0);
  }
  return (bra==ket) ? val : complex(0.0);
}

int getmz2(float_t mz)
{
  mz*=2;
  int mz2=int(mz);
  if (fabs((mz-mz2))>1e-4)
    throw InvalidParameter("getmz2: mz must be an integer or half-integer");
  return mz2;
}

static void mz2states(BaseList<state_t> &states,int &scount,state_t sofar,int left,int current)
{
  if (left==0) {
    states(scount++)=sofar;
    return;
  }
  if (left==current) {
    const state_t toadd= sofar | ((1<<current)-1);
    states(scount++)=toadd;
    return;
  }
  current--;
  mz2states(states,scount,sofar,left,current);
  mz2states(states,scount,sofar | (1<<current),left-1,current);
}

List<state_t> mzstates(size_t nspins,float_t mz)
{
  const int mz2=getmz2(mz);
  const int bits=(nspins-mz2)/2;

  if (bits<0 || nspins<bits)
    throw InvalidParameter("mzstates: mz out of range");

  ScratchList<state_t> states_d(1<<nspins);
  int nstates=0;
  mz2states(states_d,nstates,0,bits,nspins);
  return List<state_t>(states_d.truncate(nstates));
}

List<state_t> spinhalf_system::Lstates() const
{
  const BaseList<state_t> bras=brastates();
  const BaseList<state_t> kets=ketstates();

  const size_t nbras=bras.length();
  const size_t nkets=kets.length();
  
  List<state_t> fstates(mxflag::temporary);
  fstates.create(nbras*nkets);

  size_t dptr=0;
  for (size_t nb=0;nb<nbras;nb++) {
    const state_t braL= bras(nb) << _nspins;
    for (size_t nk=0;nk<nkets;nk++)
      fstates(dptr++)= kets(nk) | braL;
  }
  return fstates;
}

state_t insert(state_t ostate,size_t n,size_t nbits,size_t bitvals)
{
  const state_t mask= (1<<n)-1;
  return ( (ostate & ~mask)<<nbits) | (bitvals << n) | (ostate & mask);
}
    
template<typename T> void _expandsuperop(Matrix<T>& d, const spinhalf_system& ax,size_t index,const Matrix<T>& a)
{
  if (a.rows()!=4 || a.cols()!=4)
    throw Mismatch("expandsuperop");

  const size_t ns=ax.nspins();
  const size_t dim=1 << (2*ns);
  
  d.create(dim,dim,T(0.0));

  size_t rawins[4];
  BaseList<size_t> ins(4,rawins);

  for (state_t ostate=dim >> 2;ostate--;) {
    ins(0)=insert(insert(ostate,index,1,0),index+ns,1,0);
    ins(1)=insert(insert(ostate,index,1,0),index+ns,1,1);
    ins(2)=insert(insert(ostate,index,1,1),index+ns,1,0);
    ins(3)=insert(insert(ostate,index,1,1),index+ns,1,1);
    d(ins,ins)=a;
  }  
}

void spinhalf_system::expandsuperop(rmatrix& d,size_t index,const rmatrix &a) const
{
  _expandsuperop(d,*this,index,a);
}

void spinhalf_system::expandsuperop(cmatrix& d,size_t index,const cmatrix &a) const
{
  _expandsuperop(d,*this,index,a);
}

List<state_t> spinhalf_system::mzstates(const char *label,float_t mz) const
{
  if (!iscomplete())
    throw Failed("spinhalf_system::mzstates only valid for full Hilbert space");

  size_t i;

  const size_t nuc=labeltonuc(label);
  const size_t ntot=nspins();
  ScratchList<bool> issel(ntot,false);
  size_t ptr=0;
  
  const spinhalf_system& ax=*this;

  for (i=0;i<ntot;i++) {
    if (ax(i).nucleus()==nuc) {
      issel(i)=true;
      ptr++;
    }
  }
  if (ptr==0)
    throw Failed("No nuclei of this type!");

  const List<state_t> basestates=::libcmatrix::mzstates(ptr,mz);
  const size_t nbase=basestates.length();
  const size_t ostates=1<< (ntot-ptr);

  List<state_t> fstates(mxflag::temporary);
  fstates.create(nbase*ostates);

  size_t dptr=0;
  state_t sbit;
  for (i=0;i<nbase;i++) {
    for (state_t oth=0;oth<ostates;oth++) {
      state_t cstate=0;
      state_t ostate=oth;
      state_t bstate=basestates(i);
      for (size_t j=ntot;j--;) {
	if (issel(j)) {
	  sbit=bstate & 1;
	  bstate >>=1;
	}
	else {
	  sbit = ostate & 1;
	  ostate >>=1;
	}
	cstate = ((sbit << ntot) |cstate) >> 1;
      }
      fstates(dptr++)=cstate;
    }
  }  
  return fstates;
}

void spinhalf_system::dump() const 
{
  std::cout << *this;
}

}//namespace libcmatrix
