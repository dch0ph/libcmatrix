#undef LCM_SUPPRESS_VIEWS
#include <cstring>
#include <ctype.h>
#include "spin_system.h"
#include "ScratchList.h"

namespace libcmatrix {

const cmatrix& spinhalfop(char c)
{
  static const complex Si[]={complex(1.0),complex(0.0),complex(0.0),complex(1.0)};
  static const complex Sa[]={complex(1.0),complex(0.0),complex(0.0),complex(0.0)};
  static const complex Sb[]={complex(0.0),complex(0.0),complex(0.0),complex(1.0)};
  static const complex Sx[]={complex(0.0),complex(0.5),complex(0.5),complex(0.0)};
  static const complex Sy[]={complex(0.0),complex(0,-0.5),complex(0,0.5),complex(0.0)};
  static const complex Sz[]={complex(0.5),complex(0.0),complex(0.0),complex(-0.5)};
  static const complex Sp[]={complex(0.0),complex(1.0),complex(0.0),complex(0.0)};
  static const complex Sm[]={complex(0.0),complex(0.0),complex(1.0),complex(0.0)};
  
  static const cmatrix Mi(2,2,Si);
  static const cmatrix Ma(2,2,Sa);
  static const cmatrix Mb(2,2,Sb);
  static const cmatrix Mx(2,2,Sx);
  static const cmatrix My(2,2,Sy);
  static const cmatrix Mz(2,2,Sz);
  static const cmatrix Mp(2,2,Sp);
  static const cmatrix Mm(2,2,Sm);
  
  switch (tolower(c)) {
  case 'i': return Mi;
  case 'a': return Ma;
  case 'b': return Mb;
  case 'x': return Mx;
  case 'y': return My;
  case 'z': return Mz;
  case '+': case 'c': return Mp;
  case '-': return Mm;
  default:
    throw InvalidParameter("spinhalfop: invalid operator type");
  }
}

  void spinopST(cmatrix& d, size_t m, size_t r, size_t c, double scale)
  {
    if (r>=m || c>=m)
      throw BadIndex("create single transition operator",r,m,c,m);
    d.create(m,m,complex(0.0));
    d(r,c)=scale;
  }

  void spinopST(List<double>& d, size_t m, size_t r, double scale)
  {
    if (r>=m)
      throw BadIndex("create single transition operator",r,m);
    d.create(m,0.0);
    d(r)=scale;
  }
    
void spinop(cmatrix &d,int m,char c,double scale, bool restrict)
{
  if (restrict && (m & 1))
    throw InvalidParameter("spinop: central transition restriction only valid for half-integer nuclei");
  if (restrict || (m==2)) {
    d=spinhalfop(c);
//     if (m>2) {
//       switch (c) {
//       case '+': case '-': case 'c': case 'x': case 'y':
// 	scale*=0.5*m;
//       }
//     }
    if (scale!=1.0)
      d*=scale;
    return;
  }
  if (c=='a' || c=='b' || m<1)
    throw InvalidParameter("spinop");

  d.create(m,m); d=0.0;

  if (c=='i') {
    for (;m--;)
      d(m,m)=scale;
    return;
  }
  if (m==1)
    throw Failed("Only operator for spin=0 is 'i'");

  double qnum=(m-1)/2.0;
  double ii=qnum*(qnum+1);
  double q=qnum-1;
  double x;

  int i;

  switch (c) {
  case 'z':
    q=-qnum;
    for (;m--;q+=1.0)
      d(m,m)=scale*q;
    return;
  case 'c':
    if (m % 2)
      throw Failed("'c' only valid for half-integer nuclei");
    m>>=1;
    d(m-1,m)=scale*std::sqrt(ii+0.25); //!< equivalent to +1 on central transition
    return;
  case 'x':
    for (i=0;i<m-1;i++) {
      d(i,i+1)=(d(i+1,i)=scale*std::sqrt(ii-q*(q+1))/2);
      q-=1.0;
    }
    return;
  case 'y':
    for (i=0;i<m-1;i++) {
      x=scale*std::sqrt(ii-q*(q+1))/2;
      d(i+1,i)=complex(0,x);
      d(i,i+1)=complex(0,-x);
      q-=1.0;
    }
    return;
  case '+':
    for (i=0;i<m-1;i++) {
      d(i,i+1)=scale*std::sqrt(ii-q*(q+1));
      q-=1.0;
    }
    return;
  case '-':
    for (i=0;i<m-1;i++) {
      d(i+1,i)=scale*std::sqrt(ii-q*(q+1));
      q-=1.0;
    }
    return;
  }
  throw InvalidParameter("spinop: invalid operator");
}
 
void spin_system::I(cmatrix &dest,size_t index,char op,double scale) const
{
  if (op=='z') {
    ScratchList<double> tmp(size());
    tmp=0.0;
    mla_Iz(tmp,scale,index);
    full(dest,tmp);
  }
  else {
    cmatrix tmp;
    const spin& cspin((*this)(index));
    spinop(tmp,cspin.deg(),op,scale,cspin.isrestricted());
    expandop_tmp(dest,index,tmp);
  }
}

void spin_system::ST(cmatrix &dest,size_t index,size_t r,size_t c,double scale) const
{
  cmatrix tmp;
  spinopST(tmp,(*this)(index).deg(),r,c,scale);
  expandop_tmp(dest,index,tmp);
}

spin_system::spin_system(const basespin_system& sys, size_t rep)
  : basespin_system(sys,rep)
{}

bool spin_system::isspinhalfonly() const
{
  for (size_t n=nspins();n--;) {
    if ((*this)(n).deg()!=2)
      return false;
  }
  return true;
}

void spin_system::mla_I(cmatrix& final,double scale,size_t n,const BaseList<char>& ops) const
{
  cmatrix tmp2,tmp,dest;
  const spin& cspin((*this)(n));
  const int degv=cspin.deg();

  const size_t nops=ops.length();
  for (size_t i=0;i<nops;i++) {
    const char op=tolower(ops(i));
    if (op=='i')
      continue;
    if (isdefined(tmp)) {
      spinop(tmp2,degv,op,1.0,cspin.isrestricted());
      multiply(dest,tmp,tmp2);
      dest.swap(tmp);
    }
    else
      spinop(tmp,degv,op,1.0,cspin.isrestricted());
  }
  if (isdefined(tmp))
    expandop_tmp(dest,n,tmp);
  else
    dest.identity(size());
  mla(final,scale,dest);
}
    
template<typename T> void inline kronecker_ip(cmatrix& dest, const T& a, cmatrix& tmp)
{
  kronecker(tmp,dest,a);
  dest.swap(tmp);
}

void spin_system::mla_I(cmatrix &final,double scale,const BaseList<char> &ops) const
{
  const size_t spins=nspins();
  if (spins!=ops.length())
    throw Mismatch("mla_I");

  cmatrix tmp,dest,singlespinop;
  const spin_system& sys=*this;

  int acc=1;
  for (size_t i=0;i<spins;i++) {
    const spin& cspin(sys(i));
    if (tolower(ops(i))=='i')
      acc*=cspin.effective_deg();
    else {
      if (acc>1)
	kronecker_ip(dest,acc,tmp);
      spinop(singlespinop,cspin.deg(),ops(i),scale,cspin.isrestricted());
      scale=1.0;
      kronecker_ip(dest,singlespinop,tmp);
      acc=1;
    }
  }

  if (acc>1)
    kronecker_ip(dest,acc,tmp);
  final+=dest;
}

void spin_system::mla_I(cmatrix &dest,double scale,size_t i,char opi,size_t j,char opj) const
{
  const size_t spins=nspins();
  if (i==j) {
    if (i>=spins)
      throw BadIndex("mla_I");
    char ops[2];
    ops[0]=opi;
    ops[1]=opj;
    mla_I(dest,scale,i,BaseList<char>(2,ops));
  }
  else {
    if (i>=spins || j>=spins)
      throw BadIndex("mla_I: bad spin");

    ScratchList<char> x(spins);
    x='i';
    x(i)=opi;
    x(j)=opj;
    
    mla_I(dest,scale,x);
  }
}

void zspinop(List<double>& z,int m,double scale, bool restrict)
{
  if (m<0)
    throw InvalidParameter("zspinop");
  if (restrict && (m & 1))
    throw InvalidParameter("zspinop: central transition restriction only valid for half-integer nuclei");
  if (restrict || (m==2)) {
    scale*=0.5;
    z.create(2U);
    z(0U)=scale;
    z(1U)=-scale;
  }
  else {
    z.create(m);
    double q=(1-m)/2.0;
    for (;m--;q+=1.0)
      z(m)=scale*q;
  }
}

void spin_system::mla_I(cmatrix& dest,double scale,size_t n,char op) const
{
  if (isdefined(dest)) {
    cmatrix tmp;
    I(tmp,n,op,scale);
    dest+=tmp;
  }
  else
    I(dest,n,op,scale);
}

void spin_system::mla_ST(cmatrix& dest,double scale,size_t n,size_t r, size_t c) const
{
  if (isdefined(dest)) {
    cmatrix tmp;
    ST(tmp,n,r,c,scale);
    dest+=tmp;
  }
  else
    ST(dest,n,r,c,scale);
}

  void spin_system::mla_ST(List<double>& dest,double scale,size_t n,size_t r) const
{
  List<double> tmp;
  spinopST(tmp,(*this)(n).effective_deg(),r,scale);
  mla_expandop(dest,n,tmp);
}

void spin_system::rawmla_Iz(BaseList<double>& dest,double scale,size_t n) const
{
  List<double> tmp;
  const spin& cspin((*this)(n));
  zspinop(tmp,cspin.deg(),scale,cspin.isrestricted());
  mla_expandop(dest,n,tmp);
}

void spin_system::rawmla_Iz(BaseList<double> &dest,double scale,size_t i,size_t j) const
{
  List<double> tmps;
  List<double> tmp1(size(),0.0);
  const spin& cspini((*this)(i));
  const int degi=cspini.deg();
  zspinop(tmps,degi,1.0,cspini.isrestricted());
  mla_expandop(tmp1,i,tmps);
  
  if (i==j)
    tmp1*=tmp1;
  else {
    List<double> tmp2(size(),0.0);
    const spin& cspinj((*this)(j));
    const int degj=cspinj.deg();
    if (degi!=degj)
      zspinop(tmps,degj,1.0,cspinj.isrestricted());
    mla_expandop(tmp2,j,tmps);
    tmp1*=tmp2;
  }
  mla(dest,scale,tmp1);
}

void spin_system::get_states(BaseList<size_t>& states, state_t ostate) const
{
  size_t n=nspins();
  if (n!=states.length())
    throw Mismatch("get_states");
  for (;n--;) {
    const size_t deg=(*this)(n).effective_deg();
    const size_t state=ostate % deg;
    states(n)=state;
    ostate=(ostate-state)/deg;
  }
}

state_t spin_system::from_states(const BaseList<size_t>& states) const
{
  state_t state(0);
  const size_t spins=nspins();
  if (spins!=states.length())
    throw Mismatch("from_states");
  const spin_system& sys=*this;
  for (size_t n=0;n<spins;n++)
    state= state*(sys(n).effective_deg())+states(n);
  return state;
}

state_t spin_system::rawpermute(state_t orig, const Permutation& permvec) const
{
  const size_t n=nspins();
  ScratchList<size_t,2*SCRATCH_SIZE> basescr(2*n);  
  BaseList<size_t> origstates(n,basescr.vector());
  BaseList<size_t> finalstates(n,basescr.vector()+n);
  
  get_states(origstates,orig);
  permvec.apply(finalstates,origstates);
  return from_states(finalstates);
}

void spin_system::dump() const 
{
  std::cout << (*this);
}

}//namespace libcmatrix
