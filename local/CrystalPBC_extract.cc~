#undef LCM_SUPPRESS_VIEWS
#include "CrystalSystem.h"
#include "space_T.h"

#ifdef ENABLE_GENERICPERIODIC
#include "cmatrix_hash_map.h"
#endif

namespace libcmatrix {

  const size_t invalid_value=-1;
  const size_t MAXID_LEN=128;

  double CrystalGeometry::dipolar_coupling(const basespin_system& sys,size_t j,size_t sk) const
  {
    double d;
    Euler PAS;
    dipolar_parameters(d,PAS,sys,j,sk);
    return real(rotate(spatial_tensor(d),2,0,PAS));
  }

  space_T CrystalGeometry::dipolar_tensor(const basespin_system& sys,size_t j,size_t sk) const
  {
    double d;
    Euler PAS;
    dipolar_parameters(d,PAS,sys,j,sk);
    return rotate(spatial_tensor(d),PAS);
  }

  double CrystalGeometry::dipolar_coupling_unscaled(const basespin_system& sys,size_t j,size_t sk) const
  {
    const spherical sphco((*this)(j,sk));
    if (sphco.r==0.0)
      throw Failed("dipolar_coupling_unscaled: spins are co-incident");
    return ::libcmatrix::dipolar_coupling(sys(j).gamma(),sys(spin_to_spin(sk)).gamma(),sphco.r*1e-10);
  }

  void CrystalGeometry::dipolar_parameters(double& d, Euler& PAS,const basespin_system& sys,size_t j,size_t sk) const
  {
    const spherical sphco((*this)(j,sk));
    if (sphco.r==0.0)
      throw Failed("dipolar_parameters: spins are co-incident");
    d=::libcmatrix::dipolar_coupling(sys(j).gamma(),sys(spin_to_spin(sk)).gamma(),sphco.r*1e-10);
#ifdef LCM_USE_VECTOR_TO_EULER
    PAS=vector_to_Euler(sphco);
#else
    PAS=Euler(0,sphco.theta,sphco.phi);
#endif
  }
  
  template<int N> void CrystalGeometry::mla(vector3& d,int n,Int2Type<N> ndim) const
  {
    if (n==0)
      return;
    const int nmax(dimension(ndim));
    const vector3& a(cell_vectors_(ndims_+N-3));
    if (2*n>nmax)
      n-=nmax;
    d.x+=n*a.x;
    d.y+=n*a.y;
    d.z+=n*a.z;
  }

  void CrystalGeometry::print(std::ostream& ostr) const
  {
    static const char labels[]="abc";
    for (size_t m=0;m<unit_cell_.size();m++)
      ostr << "Spin " << m << ": " << unit_cell_(m) << " A\n";
    ostr << '\n';
    for (size_t m=0;m<cell_vectors_.size();m++)
      ostr << "Vector " << labels[m] << ": " << cell_vectors_(m) << " A x " << dimension(m) << '\n';
  }
  
  void CrystalGeometry::reverse(BaseList<size_t> d,size_t& k,size_t sk) const {
    const size_t n=d.size();
    if (n!=dimensions())
      throw Mismatch("CrystalGeometry::reverse");
    size_t ns[3];
    indexer_.reverse(ns[0],ns[1],ns[2],k,sk);
    d=ns+3-n;
//     const size_t offset=3-n;
//     for (size_t i=n;i--;)
//       d(i)=ns[offset+i];
  }

  vector3 CrystalGeometry::operator()(size_t j) const
  {
    if (j>=nspins())
      throw BadIndex("CrystalGeometry()");
    size_t nk,na,nb,nc;
    indexer_.reverse(na,nb,nc,nk,j);
    vector3 pos(unit_cell_(nk));
    this->mla(pos,na,Int2Type<0>());
    this->mla(pos,nb,Int2Type<1>());
    this->mla(pos,nc,Int2Type<2>());
    return pos;
  }

  //MUST MATCH corresponding function in CrystalStructure

vector3 CrystalGeometry::operator()(size_t j,size_t k) const
{
  if (j>=nspins_cell())
    throw BadIndex("CrystalGeometry()");

  size_t nk,na,nb,nc;
  indexer_.reverse(na,nb,nc,nk,k);
  const bool higher=(dimensions()>1);
  
  const size_t cell=cell_to_cell_(na,nb,nc);
  vector3 diff(unit_cell_(nk)-unit_cell_(j));

  if (cell && (nk>j)) {
    const size_t swapcell0(inverse(na,Int2Type<0>()));
    const size_t swapcell1(inverse(nb,Int2Type<1>()));      
    const size_t swapcell2(inverse(nc,Int2Type<2>()));
    if (cell_to_cell_(swapcell0,swapcell1,swapcell2)==cell) {
      //equidistant forward and back - swap if nk>j
      na=swapcell0;
      nb=swapcell1;
      nc=swapcell2;
      diff=-diff;
    }      
//     else {
//       size_t newj,newk;
//       permute(newj,newk,j,nk);
//       if (newk<newj)
// 	diff=unit_cell_(newk)-unit_cell_(newj);
//     }
//     }
//     else {
//       if (swapcell0<na) //find shortest
//  	na=swapcell0;
//       if (swapcell1<nb)
//  	nb=swapcell1;
//       if (swapcell2<nc)
//  	nc=swapcell2;
//     }
  }

//   if (j!=nk) {
//     //bool isswap=(2*nc==dimension(Int2Type<2>()));
//     //if (higher)
//     //  isswap &= ((2*nb==dimension(Int2Type<1>())) & (2*na==dimension(Int2Type<0>())));
//     if ((nk<j) && isswap)
//       diff=unit_cell_(j)-unit_cell_(nk);
//     else
//       diff=unit_cell_(nk)-unit_cell_(j);
//   }
  if (higher) {
    this->mla(diff,na,Int2Type<0>());
    this->mla(diff,nb,Int2Type<1>());
  }
  this->mla(diff,nc,Int2Type<2>());
#ifndef NDEBUG
  std::cout << "Coupling " << j << "-" << k << " resolved to vector " << j << '-' << na << ',' << nb << ',' << nc << ',' << nk << ": " << diff << '\n';
#endif
  return diff;
}

//   void CrystalGeometry::permute(size_t& newj, size_t& newk, size_t j, size_t k) const
//   {
//     if (perm_.empty()) {
//       newj=j;
//       newk=k;
//       return;
//     }
//     newj=perm_(j);
//     newk=perm_(k);
//   }
    
  CrystalGeometry::CrystalGeometry(size_t n, const vector3& vectora, const BaseList<vector3>& unit_cellv)//, const Permutation& permv)
    : ndims_(1), 
      indexer_(1,1,n,unit_cellv.size()),
      cell_to_cell_(1,1,n),
      unit_cell_(unit_cellv), cell_vectors_(1,vectora) {
      //    perm_(permv) {
  verify();
}
  
  CrystalGeometry::CrystalGeometry(const BaseList<size_t>& dims, const BaseList<vector3>& vectors, const BaseList<vector3>& unit_cellv)//, const Permutation& permv)
    : ndims_(dims.size()), unit_cell_(unit_cellv), cell_vectors_(vectors) {
    //      perm_(permv)        {
    if (ndims_>3)
      throw InvalidParameter("CrystalGeometry: >3 dimensions!");
    if (vectors.size()!=ndims_)
      throw Mismatch("CrystalGeometry");
    const size_t na = (ndims_>2) ? dims.front() : 1;
    const size_t nb = (ndims_>1) ? dims(ndims_-2) : 1;
    const size_t nc = ndims_ ? dims.back() : 1;
    indexer_.setdims(na,nb,nc,nspins_cell());
    cell_to_cell_.setdims(na,nb,nc);
    verify();
  }

void CrystalGeometry::verify()
{
//   if (!perm_.empty() && (perm_.size()!=nspins_cell()))
//     throw Mismatch("CrystalGeometry");
  for (size_t n=dimensions();n--;) {
    if ((dimension(n)>1) && (cell_vectors_(n).norm()==0.0))
      throw InvalidParameter("cell vector of zero length");
  }
}

/****************************************************************/

 //    size_t log2(state_t s)
//     {
//       size_t res=0;
//       while ((s & 1)==0) {
// 	res++;
// 	s>>=1;
//       }
//       return res;
//     }

#ifdef LCM_ENABLE_GENERICPERIODIC

class bitroller {
public:
  bitroller() : donothing(true) {}

  bitroller(size_t bitsv, size_t totbitsv, size_t totalv)
    : permp(NULL) {
    init(bitsv,totbitsv,totalv);
  }

  bitroller(const Permutation& perm, size_t totalv)
    : permp(&perm) {
    init(perm.size(),perm.size(),totalv);
  }

  state_t operator()(state_t in) const;

private:  
  size_t bits,totbits,total,shift;
  const Permutation* permp;
  state_t chunkmask,mask;
  bool donothing;

  void init(size_t bitsv, size_t totbitsv, size_t totalv);       
};

void bitroller::init(size_t bitsv, size_t totbitsv, size_t totalv)
{
  bits=bitsv;
  totbits=totbitsv;
  total=totalv;
  
  if ((bits>totbits) || (totbits>total) || (total % totbits) || (totbits % bits))
    throw InvalidParameter("bitroller");
  shift=totbits-bits;
  chunkmask=(state_t(1)<<totbits)-1;
  mask=(state_t(1)<<bits)-1;    
  donothing=((shift==0) && (!permp || (permp->order()==1)));
}

state_t bitroller::operator()(state_t in) const 
{
  if (donothing)
    return in;
  state_t out=0;
  for (size_t chunkbase=0;chunkbase<total;chunkbase+=totbits) {
    state_t chunk= in & chunkmask;
    in >>= totbits;
    if (permp)
      chunk=spinhalf_system::apply_permutation(chunk,*permp);
    else
      chunk = (chunk >> bits) | ((chunk & mask) << shift);
    out |= chunk << chunkbase;
  }
  return out;
}

struct special_shape {
  std::string id;
  magicsize_t shape;

  special_shape() : shape(0) {}
  special_shape(const char* id_, magicsize_t shape_)
    : id(id_), shape(shape_) {}
};

struct cache_arrays {
  mutable List< List<complex> > eigfacs;
  mutable List<cmatrix> smalltransforms;
  List< List<size_t> > allowed;
  //  List< List<bool> > lisreal;
  
  cache_arrays(size_t n)
    : eigfacs(n), smalltransforms(n), allowed(n) {}
  void create(size_t);
  void make_transform(cmatrix&, size_t[]) const;
};

struct general_creator {
  general_creator(const CrystalStructure&, int M, GeneralSymmetrise&, const cache_arrays&);

  magicsize_t push_states(List<state_t>& state_list, BaseList<bool> used, state_t current);
  magicsize_t push_entangled(List<state_t>& state_list, const BaseList<state_t>&);
  magicsize_t push_entangled3(List<state_t>& state_list, const BaseList<state_t>&);

  static bool areentangled(int&, int&, const ListList<size_t>&, size_t, size_t);
  state_t roll(state_t current, size_t) const;
  static size_t findstep(const BaseList<size_t>& longperm, const BaseList<size_t>& perm);
  size_t findk(const complex&, size_t) const;
  void create_perm(List<size_t>&, size_t& innercount, size_t& outercount, const BaseList<state_t>& states, const bitroller& innerroll, const bitroller& outerroll) const;

  magicsize_t maxshape() const { return indexer.size(); }

  magicsize_t makespecial();
  static void makeidstr(char idstr[MAXID_LEN],const ListList<size_t>&);

  GeneralSymmetrise& symmetriser;
  const Indexer<3>& indexer;
  size_t total;
  bitroller roll0,roll1,roll2;
  const cache_arrays& arrays;

  stringhash hasher;
  typedef std::map<size_t,special_shape> specialmap_t;
  specialmap_t specialmap;

  size_t specialcount;
  List<state_t> entangled;
  List<int> reverse;
  List<size_t> cycle;
  List<int> dummyeigrow;
};

general_creator::general_creator(const CrystalStructure& cstruct, int M, GeneralSymmetrise& symmetriserv, const cache_arrays& arrays_)
  : symmetriser(symmetriserv), indexer(symmetriser.indexer), total(cstruct.ncells()*M), arrays(arrays_)
 {
   if (symmetriser.perm.empty()) {
    roll0=bitroller(M*cstruct.indexer().multiplier(0),total,total);
    roll1=bitroller(M*cstruct.indexer().multiplier(1),M*cstruct.indexer().multiplier(0),total);
    roll2=bitroller(M,M*cstruct.indexer().multiplier(1),total);
   }
  else {
    if (M!=symmetriser.perm.size())
      throw Mismatch("GeneralSymmetriser: permutation vs. spins-per-cell");
    roll0=bitroller(M*cstruct.indexer().multiplier(1),total,total);
    roll1=bitroller(M,M*cstruct.indexer().multiplier(1),total);
    roll2=bitroller(symmetriser.perm,total);
  }
   specialcount=0;
   dummyeigrow.create(symmetriser.eigindex.cols(),-1);
}

size_t general_creator::findstep(const BaseList<size_t>& longperm, const BaseList<size_t>& perm)
{
  const BaseList<size_t>::const_iterator where(::std::find(longperm.begin(),longperm.end(),perm.front()));
  if (where==longperm.end()) {
    std::cerr << "Failed to find overlap in entangled states: No " << perm.front() << " in " << longperm << std::endl;
    throw InternalError("findstep:1");
  }
  const size_t step=where-longperm.begin()+1; //k step
#ifndef NDEBUG
  const size_t longn=longperm.size()+1;
  const size_t permn=perm.size()+1;
  for (size_t m=2;m<=permn;m++) {
    const size_t lpos=(m*step) % longn;
    const size_t lval= lpos ? longperm(lpos-1) : 0;
    const size_t pval= (m==permn) ? 0 : perm(m-1);
    if (pval!=lval) {
      std::cerr << "Failed to find cyclic relationship between " << perm << " and " << longperm << std::endl;
      throw InternalError("findstep:2");	
    }
  }
#endif
  return step;
}

bool general_creator::areentangled(int& eshorti, int& elongi, const ListList<size_t>& perms, size_t longeri, size_t shorteri)
{
  if (perms.size(shorteri)>perms.size(longeri))
    ::std::swap(shorteri,longeri);
  const BaseList<size_t> shorter(perms(shorteri));
  if (shorter.empty())
    return false;
  const BaseList<size_t> longer(perms(longeri));
  const BaseList<size_t>::const_iterator end(longer.end());
  const BaseList<size_t>::const_iterator begin(longer.begin());
  bool allpresent=true;
  bool allmissing=true;
  for (size_t j=shorter.size();j--;) {
    if (std::find(begin,end,shorter(j))==end)
      allpresent=false;
    else
      allmissing=false;
  }
  if (allpresent ^ allmissing) {
    if (allmissing)
      return false;
    if ((eshorti<0) || (perms.size(eshorti)>shorter.size()))
      eshorti=shorteri;
   if ((elongi<0) || (perms.size(elongi)<longer.size()))
      elongi=longeri;
    return true;
  }
  std::cerr << "Lists partially intersect: " << longer << " & " << shorter << '\n';
  throw InternalError("areentangled");
}
    
state_t general_creator::roll(state_t current, size_t which) const
{
  switch (which) {
  case 0:
    return roll0(current);
  case 1:
    return roll1(current);
  case 2:
    return roll2(current);
  }
  throw InvalidParameter("roll");
}

void general_creator::create_perm(List<size_t>& perm, size_t& innercount, size_t& outercount, const BaseList<state_t>& states, const bitroller& innerroll, const bitroller& outerroll) const
{
  const size_t nstates(states.size());
  perm.reserve(nstates);
  perm.create(0);

  state_t current=states.front();
  const state_t start1=current;
  outercount=0;
  do {
    const state_t start2=current;
    innercount=0;
    do {
      const int whichst(reverse(current));
      assert(whichst>=0);
      perm.push_back(whichst);
      current=innerroll(current);
      innercount++;
    } while (current!=start2);
    current=outerroll(current);
    outercount++;
  } while (current!=start1);
  if (perm.size()!=nstates) {
    std::cerr << "create_perm: permutation on 2 dimensions did not generate " << nstates << " states. Found: " << perm << std::endl;
    throw InternalError("create_perm");
  }  
}

size_t general_creator::findk(const complex& z, size_t n) const
{ 
  List<complex>& ceigfacs(arrays.eigfacs(n));
  for (size_t j=n;j--;) {
    if (norm(z-ceigfacs(j))<1e-6)
      return j;
  }
  std::cerr << z << " doesn't match any of the " << n << " th roots of 1\n";
  throw Failed("findk");
}

magicsize_t general_creator::makespecial()
{
  const magicsize_t shape=indexer.size()+specialcount;

  if (symmetriser.verbose)
    std::cout << "Creating special transform " << shape << '\n';
  
  symmetriser.transforms.push_back(cmatrix()); //create space
  symmetriser.extraeigindex.push_back(dummyeigrow);
  specialcount++;
  return shape;
}

magicsize_t general_creator::push_entangled3(List<state_t>& state_list, const BaseList<state_t>& states)
{
  const int verbose=symmetriser.verbose;
  const magicsize_t shape=makespecial();

  List<size_t> xperm;
  size_t count1,count2;
  create_perm(xperm,count2,count1,states,roll2,roll1);

  const magicsize_t subshape=symmetriser.getshape(1,count1,count2);
  cmatrix& ctransform(symmetriser.transforms(shape));
  ctransform=symmetriser.transforms(subshape);//transformation matrix

  const size_t k0max(indexer.dimension(Int2Type<0>()));
  const size_t k1max(indexer.dimension(Int2Type<1>()));
  const size_t k1step=k1max/count1;
  const size_t k2max(indexer.dimension(Int2Type<2>()));
  const size_t k2step=k2max/count2;

  //create z permutation and store states
  const size_t nstates(states.size());
  List<int> zperm(nstates,-1);
  for (size_t j=0;j<nstates;j++) {
    state_list.push_back(states(j));
    const int which(reverse(roll0(states(j))));
    assert(which>=0);
    if (zperm(which)>=0) {
      std::cerr << "z permutation has wrong structure" << std::endl;
      throw Failed("create_perm");
    }
    else
      zperm(which)=j;
  }
  if (verbose>1)
    std::cout << "z permutation: " << zperm << '\n';

  size_t index=0;
  const double normfac=real(ctransform(0U,0U));
  for (size_t k1=0;k1<k1max;k1+=k1step) {
    for (size_t k2=0;k2<k2max;k2+=k2step) {
      const BaseList<complex>& Vrow(ctransform.row(index)); //current mixing 
      size_t k0=findk(Vrow(zperm.front())/normfac,k0max);
#ifndef NDEBUG
      for (size_t j=1;j<nstates;j++) {
	const size_t lk0=findk(Vrow(zperm(j))/Vrow(j),k0max);
	if (k0!=lk0) {
	  std::cerr << "Mismatch between z eigenvalue at state " << j << ": " << k0 << " vs. " << lk0 << std::endl;
	  throw Failed("push_entangled3");
	}
      }
#endif
      symmetriser.addeigenvalue(index,shape,indexer(k0,k1,k2));
    }
  }
  return shape;
}

void general_creator::makeidstr(char idstr[MAXID_LEN],const ListList<size_t>& perms)
{
  for (size_t i=0;i<perms.size();i++) {
    if (i)
      *idstr++=' '; //space between runs
    const BaseList<size_t> cperm(perms(i));
    for (size_t j=0;j<cperm.size();j++)
      *idstr++='A'+cperm(j);
  }
  *idstr='\0';
}

magicsize_t general_creator::push_entangled(List<state_t>& state_list, const BaseList<state_t>& states)
{
  const int verbose=symmetriser.verbose;
  if (verbose) {
    std::cout << "Entangled states: ";
    printbin(std::cout,states,total); std::cout << '\n';
  }
  const size_t nstates(states.size());
  reverse.create(symmetriser.maxstates(),-1);
  //construct indices
  size_t j;
  for (j=states.size();j--;)
    reverse(states(j))=j;
  ListList<size_t> permstructure;
  size_t maxdim=0;
  const char names[]="zyx";
  const size_t initstate=states.front();
  cycle.create(nstates);
  for (j=0;j<3;j++) {
    cycle.create(0);
    size_t current=initstate;
    for (;;) {
      current=roll(current,j);
      if (current==initstate)
	break;
      const int ind=reverse(current);
      if (ind<0)
	throw InternalError("push_entangled:2");
      cycle.push_back(ind);
    }
    if (verbose>1)
      std::cout << "T" << names[j] << ": " << cycle << '\n';
    const size_t dim=cycle.length()+1;
    if (dim>maxdim)
      maxdim=dim;
    permstructure.push_back(cycle);
  }
  if (nstates % maxdim) {
    std::cerr << "Doesn't divide number of states!\n";
    throw InternalError("push_entangled:1");
  }
  size_t entanglenum=0;
  enum { EXY=1, EXZ=2, EYZ=4 };
  int shorti=-1,longi=-1;
  int sparedim=-1;
  if (areentangled(shorti,longi,permstructure,0,1)) {
    entanglenum|=EYZ;
    sparedim=2;
  }
  if (areentangled(shorti,longi,permstructure,0,2)) {
    entanglenum|=EXZ;
    sparedim=1;
  }
  if (areentangled(shorti,longi,permstructure,1,2)) {
    entanglenum|=EXY;
    sparedim=0;
  }
  switch (entanglenum) {
  case 0:
    if (maxdim<nstates)
      return push_entangled3(state_list,states);
    throw InternalError("push_entangled:3");
  case EXY: case EXZ: case EYZ:
    break;
  case 7:
    //return push_entangled3(state_list,states);
    sparedim=-1;
    break;
  default:
    std::cerr << "Unexpected entanglenum: " << entanglenum << std::endl;
    throw InternalError("push_entangled:4");
  }

  char idstr[MAXID_LEN];
  makeidstr(idstr,permstructure);

  if (verbose)
    std::cout << "Maximum run: " << maxdim << "  Entangle type: " << entanglenum << "  Structure: " << idstr << '\n';

  const BaseList<size_t> longinds(permstructure(longi));
  const size_t n=longinds.size()+1;
  ScratchList<state_t> lstates(n);
  lstates.front()=states.front();
  for (j=n-1;j--;)
    lstates(j+1)=states(longinds(j));
  if (verbose>1)
    std::cout << "Base states: " << lstates << '\n';
  const size_t cursize(state_list.size());
  for (;;) {
    state_list.push_back(lstates);
    if (sparedim>=0) {
      for (j=n;j--;)
	lstates(j)=roll(lstates(j),sparedim);
    }
    if (lstates.front()==initstate)
      break;
  }
  if (state_list.size()-cursize!=nstates) {
    std::cerr << "Mismatch between initial and final number of states\n";
    throw InternalError("push_entangled:5");
  }
  
  const size_t hashedid(hasher(idstr));
  specialmap_t::iterator curpos(specialmap.find(hashedid));
  magicsize_t shape=0;
  if (curpos==specialmap.end()) {
    shape=makespecial();
    specialmap[hashedid]=special_shape(idstr,shape); //store index in map
    cmatrix& ctransform(symmetriser.transforms(shape));
    //    symmetriser.cache(n);
    ScratchList<size_t> ks(3); //ks(3,0) creates two element array!
    ks=size_t(0);
    size_t index=0;

    if (sparedim>=0) {
      const size_t sparen(1+permstructure.size(sparedim));
      const size_t sparemax(indexer.dimension(sparedim));
      const size_t sparestep=sparemax/sparen;
      ScratchList<size_t> ns(size_t(1),sparen,n); //DODGY
      arrays.make_transform(ctransform,ns.vector());
      const size_t longmax(indexer.dimension(longi));
      const size_t longstep=longmax/n;
      const size_t shortmax(indexer.dimension(shorti));
      const size_t shortstep=longstep*findstep(longinds,permstructure(shorti));

      size_t& shortk(ks(shorti));
      size_t& longk(ks(longi));
      size_t& sparek(ks(size_t(sparedim)));
      for (sparek=0;sparek<sparemax;sparek+=sparestep) {
	shortk=0;
	for (longk=0;longk<longmax;longk+=longstep) {
	  symmetriser.addeigenvalue(index,shape,indexer(ks(0U),ks(1U),ks(2U)));
	  shortk+=shortstep;
	  if (shortk>=shortmax)
	    shortk-=shortmax;
	}	  
      }
    }
    else {
      ctransform=arrays.smalltransforms(n);
      
      ScratchList<size_t> ksteps(3);
      const size_t longmax=indexer.dimension(longi);
      const size_t longstep= longmax/(1+longinds.size());
      size_t m;
      for (m=3;m--;) {
	if (m==longi)
	  ksteps(m)=longstep;
	else {
	  const size_t genkstep=longstep*findstep(longinds,permstructure(m));
	  if (longmax>indexer.dimension(m)) {
	    const size_t divfac=longmax/indexer.dimension(m);
	    if (genkstep % divfac)
	      throw InternalError("push_entangled:7");
	    ksteps(m)=genkstep/divfac;
	  }
	  else {
	    const size_t multfac=indexer.dimension(m)/longmax;
	    ksteps(m)=genkstep*multfac;
	  }
	}
      }
      if (verbose>1)
	std::cout << "k steps: " << ksteps << '\n';
      
      for (size_t k=n;k--;) {
 	symmetriser.addeigenvalue(index,shape,indexer(ks(0U),ks(1U),ks(2U)));
 	ks+=ksteps;
 	for (m=3;m--;) {
 	  if (ks(m)>=indexer.dimension(m))
 	    ks(m)-=indexer.dimension(m);
 	}
      }
    }
    if (verbose>1)
      std::cout << "Transformation matrix:\n" << ctransform;
  }
  else {
    shape=curpos->second.shape;
    if (verbose)
      std::cout << "Matches shape: " << shape << '\n';
  }
  return shape;
}

magicsize_t general_creator::push_states(List<state_t>& state_list, BaseList<bool> used, state_t initstate)
{
  size_t count0=0;
  size_t count1,count2;
  size_t current=initstate;
  bool isentangled=false;
  const state_t start0=current;
  const size_t initsize=state_list.size();
  for (size_t n0=indexer.dimension(0);n0--;) {
    count0++;
    count1=0;
    const state_t start1=current;
    for (size_t n1=indexer.dimension(1);n1--;) {
      count1++;
      count2=0;
      const state_t start2=current;
      for (size_t n2=indexer.dimension(2);n2--;) {
	count2++;
	if (used(current))
	  isentangled=true;
	else {
	  used(current)=true;
	  state_list.push_back(current);
	}
	current=roll2(current);
	if (current==start2)
	  break;
      }
      current=roll1(current);
      if (current==start1)
	break;
    }
    current=roll0(current);
    if (current==start0)
      break;
  }
  if (isentangled) {
    //copy out entangled states
    entangled=BaseList<state_t>(state_list.size()-initsize,state_list.vector()+initsize);
    state_list.resize(initsize);    
    return push_entangled(state_list,entangled);
  }
  return symmetriser.getshape(count0,count1,count2);
}

  void GeneralSymmetrise::addeigenvalue(size_t& index, magicsize_t shape,size_t inteig)
{
  const int exteig(eigreverse(inteig));
  if (exteig>=0) {
    if (shape>=eigindex.rows())
      extraeigindex(shape-eigindex.rows(),exteig)=index;
    else
      eigindex(shape,exteig)=index;
  }
  index++;
}

template<class T> complex 
GeneralSymmetrise::symmetrise_(Matrix<T>& H, magicsize_t blkishape, magicsize_t blkjshape, int m) const
{
  if (m==0)
    throw Failed("symmetrise: don't use for k=0!");

  const int indexi=eigindex(blkishape,m);
  const int indexj=eigindex(blkjshape,m);
  if ((indexi<0) || (indexj<0))
    throw Failed("symmetrise: eigenvalue doesn't exist in given block");

  const size_t blkisize=H.rows();
  const size_t blkjsize=H.cols(); 
  const BaseList<complex> Vrj(transforms(blkishape).row(size_t(indexi)));
  const BaseList<complex> Vcj(transforms(blkjshape).row(size_t(indexj)));
  
  complex s(0.0);
  for (size_t ii=blkisize;ii--;) {
    const BaseList<T> Hii(H.row(ii));
    complex tmp(0.0);
    for (size_t jj=blkjsize;jj--;)
      mla(tmp,Hii(jj),Vcj(jj));
    conj_mla(s,Vrj(ii),tmp);
  }
  return s*scalefacs(blkisize*blkjsize);
}

//Construct overall transform matrix for set of eigenvalue sizes
void cache_arrays::make_transform(cmatrix& dest, size_t ns[]) const
{
  for (size_t j=3;j--;) {
    if (ns[j]!=1)
      dest=kronecker(smalltransforms(ns[j]),dest);
  }
  if (!dest)
    dest.create(1,1,complex(1.0));
}

void cache_arrays::create(size_t n)
{
    List<complex>& ceigfacs(eigfacs(n));
    if (!(ceigfacs.empty()))
      return; //already created?
    ceigfacs.create(n);
    CrystalSymmetriseBase::make_eigfacs(ceigfacs);

    cmatrix& ctransform(smalltransforms(n));
    ctransform.create(n,n,complex(1.0));
    for (size_t j=1;j<n;j++)
      for (size_t k=1;k<n;k++)
	ctransform(j,k)=ceigfacs((j*k) % n);
    
    List<size_t>& callowed(allowed(n));
    callowed.reserve(n);
    for (size_t j=1;j<=n;j++) {
      if ((n % j)==0)
	callowed.push_back(j);
    }
    
    //    List<bool>& cisreal(lisreal(n));
//     cisreal.create(n,false);
//     cisreal.front()=true;
//     if ((n & 1)==0)
//       cisreal(n/2)=true;
  }

  namespace {
    struct ExistsIfPositive : public ::std::unary_function<int,bool> {
      bool operator()(int v) const { return (v>=0); }
    };
  }

GeneralSymmetrise::GeneralSymmetrise(const CrystalStructure& cstruct, int M_, int flags_, int verbose_)
    : CrystalSymmetriseBase(cstruct.ncells(),M_,cstruct.ncells()*cstruct.permutation().order(),flags_),
      perm(cstruct.permutation()),
      verbose(verbose_)
  {
    if (perm.empty())
      indexer.setdims(cstruct.dimension(Int2Type<0>()),cstruct.dimension(Int2Type<1>()),cstruct.dimension(Int2Type<2>()));
    else {
      if (cstruct.dimension(Int2Type<0>())>1)
	throw Failed("GeneralSymmetrise: can't specify 3 translation axes and permutation");
      indexer.setdims(cstruct.dimension(Int2Type<1>()),cstruct.dimension(Int2Type<2>()),perm.order());
    }

    size_t j;
    size_t maxdim=1;
    for (j=0;j<3;j++) {
      if (indexer.dimension(j)>maxdim)
	maxdim=indexer.dimension(j);
    }
    if (verbose>1)
      std::cout << "Maximum dimension: " << maxdim << '\n';
 
    const size_t N0=indexer.dimension(Int2Type<0>());
    const size_t N1=indexer.dimension(Int2Type<1>());
    const size_t N2=indexer.dimension(Int2Type<2>());
    
    if (useeigsym_) { //compute eigenvalues to use
      eigreverse.create(neigs_,-1);
      useeigs=0;
 //      const size_t useN0=1+N0/2;
//       const size_t useN1=1+N1/2;
//       const size_t useN2=1+N2/2;
      for (size_t n0=0;n0<N0;n0++) {
	const size_t rn0=n0 ? N0-n0 : 0;
	for (size_t n1=0;n1<N1;n1++) {
	  const size_t rn1=n1 ? N1-n1 : 0;
	  for (size_t n2=0;n2<N2;n2++) {
	    const size_t rn2=n2 ? N2-n2 : 0;
	    const size_t eig=indexer(n0,n1,n2);
	    const size_t reig=indexer(rn0,rn1,rn2);
	    if (eig<=reig) {
	      isreal_.push_back(eig==reig);
	      eigreverse(eig)=useeigs++;
	      eigtranslate.push_back(eig);
	    }
	  }
	}
      }
      if (verbose>1) {
	std::cout << "ext->int translation: " << eigtranslate << '\n';
	std::cout << "int->ext translation: " << eigreverse << '\n';
      }
    }
    else {
      eigtranslate.create(neigs_);
      for (j=neigs_;j--;)
	eigtranslate(j)=j;
      eigreverse=eigtranslate;
      useeigs=neigs_;
    }
    if (verbose)
      std::cout << "Eigenvalues in use: " << useeigs << '\n';
    
    cache_arrays arrays(maxdim+1);
    for (j=3;j--;)
      arrays.create(indexer.dimension(j));

    const size_t lmaxshape=indexer.size();
    transforms.create(lmaxshape);
    eigindex.create(lmaxshape,useeigs,-1);

    size_t ns[3];
    size_t& n0=ns[0];
    size_t& n1=ns[1];
    size_t& n2=ns[2];
    
    for (size_t ni0=arrays.allowed(N0).size();ni0--;) {
      n0=arrays.allowed(N0)(ni0);
      const size_t step0=N0/n0;
      for (size_t ni1=arrays.allowed(N1).size();ni1--;) {
	n1=arrays.allowed(N1)(ni1);
	const size_t step1=N1/n1;
	for (size_t ni2=arrays.allowed(N2).size();ni2--;) {
	  n2=arrays.allowed(N2)(ni2);
	  const size_t step2=N2/n2;
	  const magicsize_t shape=getshape(n0,n1,n2);
	  if (verbose>1)
	    std::cout << "Making transforms for shape " << shape << ": " << n0 << 'x' << n1 << 'x' << n2 << '\n';
	  arrays.make_transform(transforms(shape),ns);
	  if (verbose>1)
	    std::cout << transforms(shape);
	  size_t index=0;
	  for (size_t k0=0;k0<N0;k0+=step0)
	    for (size_t k1=0;k1<N1;k1+=step1)
	      for (size_t k2=0;k2<N2;k2+=step2)
		addeigenvalue(index,shape,indexer(k0,k1,k2)); //another place to screw up
	}
      }
    }
    
    general_creator creator(cstruct,M_,*this,arrays);
    getsymmetry(creator); //get sets of symmetry-linked states

    //combine eigindex and extraeigindex
    const size_t baser=eigindex.rows();
    Matrix<int> tmp(baser+extraeigindex.size(),eigindex.cols());
    const BaseList<int> asrow(eigindex.row());
    tmp.row().truncate(asrow.size())=asrow;
    for (j=extraeigindex.size();j--;)
      tmp.row(j+baser)=extraeigindex(j);
    eigindex.swap(tmp);

    apply(haseig,ExistsIfPositive(),eigindex);

    if (verbose) {
      std::cout << "Is real: " << isreal_ << '\n';
      if (verbose>1) {
	std::cout << "Shapes: " << blockshapes_ << '\n';
	std::cout << "eigindex (shape x exteig)\n" << eigindex << '\n';
      }
    }
  }
  
void GeneralSymmetrise::addtoevals(BaseList<complex>& evals, size_t blki, size_t blkj, size_t r, size_t c) const
{
  const magicsize_t blkishape(blockshapes_(blki));
  const magicsize_t blkjshape(blockshapes_(blkj));
  const size_t blkisize(linkedstates_.size(blki));
  const size_t blkjsize(linkedstates_.size(blkj));
  
  const double scale=scalefacs(blkisize*blkjsize); //normalisation factor
  
  evals(0)+=scale; //eigenvalue 0

  const BaseList<complex> transformi(transforms(blkishape).row(r));
  const BaseList<complex> transformj(transforms(blkjshape).row(c)); 
  const BaseList<int> indexis(eigindex.row(blkishape));
  const BaseList<int> indexjs(eigindex.row(blkjshape));

  for (size_t m=1;m<useeigs;m++) {
    //    const size_t inteig(eigreverse(m));
    const int indexi(indexis(m));
    const int indexj(indexjs(m));
    if ((indexi>=0) && (indexj>=0)) { //if eigenvalue present
      const complex Vtot(conj_multiply(transformi(indexi),transformj(indexj)));
      mla(evals(m),scale,Vtot);
    }
  }
}

#endif

} //namespace libcmatrix
