#include "basespin_system.h"
#include "spinhalf_system.h"
#include "rmatrix.h"
#include "ScratchList.h"
#include <string>
#include <map>
//#include "cmatrix_hash_map.h"

namespace libcmatrix {

  static const size_t gamma_internal_nucs=114;
  static const size_t harris_internal_nucs=128;
  static const size_t internal_nucs=harris_internal_nucs; //!< needs to be highest of available (checked on first use)
  const size_t MAX_NUCLEUS=internal_nucs+LCM_MAX_USER_NUC;
  const size_t DEFAULT_NUCLEUS=1; // 1H is default 
  const size_t H_NUCLEUS=1;
    
  namespace {

    char error_buf[256];

    void throw_notpresent(size_t nuc) 
    {
      LCM_SNPRINTF(error_buf,sizeof(error_buf),"mla_F: nucleus %s not present in spin system",nuctolabel(nuc));
      throw Failed(error_buf);
    }

    //! summarises nucleus properties
    struct nucleus_properties {
      size_t deg; //!< degeneracy
      const char* label; //!< pointer to nucleus name
      double gamma; //!< magnetogyric ratio
    };
    
    //! user defined nucleus type (internal type)
    struct user_nuc {
      std::string label; //!< nucleus name
      nucleus_properties props; //!< nucleus properties
      
      user_nuc(const char*, size_t, double);
    };
    
  static const nucleus_properties* used_properties_table=NULL; //!< last table used

    //! construct user defined nucleus
    /** ::Failed exception is thrown if nucleus name is already used */
    user_nuc::user_nuc(const char* labelv, size_t deg, double gamma)
      : label(labelv)
    {
      try { 
	(void)labeltonuc(labelv);
      }
      catch (Failed&) {
	props.deg=deg;
	props.label=label.c_str();
	props.gamma=gamma;	
	return;
      }
      throw Failed("user_nuc: name is already used");      
    }
    
    List<user_nuc> usernucs; //!< list of user defined nuclei
    typedef std::map<std::string,size_t> usernucs_map_t;

    usernucs_map_t& getusernucs_map() 
    {
      static usernucs_map_t usernucs_map; //!< index of user defined nuclei
      return usernucs_map;
    }
    
    const nucleus_properties* get_nucleus_properties()
    {
      if (!used_properties_table)
	set_nucleus_properties(LCM_DEFAULT_NUCLEUS_PROPERTIES);
      return used_properties_table;
    }

    //! Table of ::nucleus_properties for known nuclei
    /** gamma values from Harris and Mann, CRC handbook, GAMMA Isotopes.list
	N.B. many values are calculated and are not accurate to all s.f. given! */
    const nucleus_properties nucleus_props_GAMMA[gamma_internal_nucs]={
      {0,NULL,0},
      {2,"1H", 267515255},
      {2,"13C", 6.7263e7},
      {3, "2H", 4.1064e7},
      {3,"14N", 1.9324e7},
      {2,"15N",-2.7107e7},
      {3,"6Li", 3.9366e7},
      {4,"7Li",10.396e7},
      {4,"9Be",-3.75911e+07},
      {7,"10B", 2.8748e7},
      {4,"11B", 8.5827e7},
      {6,"17O",-3.6266e7},
      {2,"19F",25.1666e7},
      {4,"21Ne",-2.11175e+07},
      {4,"23Na", 7.076e7},
      {6,"25Mg",-1.63773e+07},
      {6,"27Al", 6.9706e7},
      {2,"29Si",-5.314e7},
      {2,"31P", 10.829e7},
      {4,"33S",2.05345e+07},
      {4,"35Cl", 2.6212e7},
      {4,"37Cl", 2.182e7},
      {4,"39K",1.24822e+07},
      {4,"41K",6.85085e+06},
      {8,"43Ca",-1.80013e+07},
      {8,"45Sc",6.49851e+07},
      {6,"47Ti",-1.5079e7},
      {8,"49Ti",-1.5083e7},
      {13,"50V",2.66713e+07},
      {8,"51V",  7.032e7},
      {4,"53Cr",-1.51197e+07},
      {6,"55Mn",6.61965e+07},
      {2,"57Fe", 0.8644e7},
      {8,"59Co",6.34734e+07},
      {4,"61Ni",-2.39054e+07},
      {4,"63Cu", 7.0906e7},
      {4,"65Cu", 7.5958e7},
      {6,"67Zn",1.67383e+07},
      {4,"69Ga",6.42062e+07},
      {4,"71Ga",8.15844e+07},
      {10,"73Ge",-9.33124e+06},
      {4,"75As",4.58065e+07},
      {2,"77Se", 5.101e7},
      {4,"79Br", 6.7021e7},
      {4,"81Br", 7.2245e7},
      {6,"85Rb",2.58289e+07},
      {4,"87Rb",8.7534e+07},
      {10,"87Sr",-1.15944e+07},
      {2,"89Y", -1.3106e7},
      {6,"91Zr",-2.48682e+07},
      {10,"93Nb",6.54772e+07},
      {6,"95Mo",1.74337e+07},
      {6,"97Mo",-1.78007e+07},
      {10,"99Tc",6.02122e+07},
      {6,"99Ru",-1.23432e+07},
      {6,"101Ru",-1.38334e+07},
      {2,"103Rh",-0.842e7},
      {2,"107Ag",-1.0825e7},
      {2,"109Ag",-1.2445e7},
      {2,"111Cd",-5.6727e7},
      {2,"113Cd",-5.934e7},
      {10,"113In",5.84946e+07},
      {10,"115In",5.86203e+07},
      {2,"117Sn",-9.5301e7},
      {2,"119Sn",-9.9708e7},
      {6,"121Sb",6.4019e+07},
      {8,"123Sb",3.46674e+07},
      {2,"123Te",-7.00087e+07},
      {2,"125Te",-8.44011e+07},
      {6,"127I",  5.3521e7},
      {2,"129Xe",-7.3995e7},
      {4,"131Xe", 2.1935e7},
      {8,"133Cs",3.50873e+07},
      {4,"135Ba",2.6575e+07},
      {4,"137Ba",2.97287e+07},
      {11,"138La",3.52959e+07},
      {8,"139La",3.7789e+07},
      {6,"141Pr",7.83578e+07},
      {8,"143Nd",-1.45448e+07},
      {8,"145Nd",-8.94815e+06},
      {8,"147Sm",-1.10428e+07},
      {8,"149Sm",-8.79839e+06},
      {6,"151Eu",6.63463e+07},
      {6,"153Eu",2.92954e+07},
      {4,"155Gd",-1.02164e+07},
      {4,"157Gd",-1.2771e+07},
      {4,"159Tb",6.06668e+07},
      {6,"161Dy",-8.81176e+06},
      {6,"163Dy",1.22603e+07},
      {8,"165Ho",5.48756e+07},
      {8,"167Er",-7.73135e+06},
      {2,"169Tm",-2.213e7},
      {2,"171Yb", 4.7117e7},
      {6,"173Yb",-1.29796e+07},
      {8,"175Lu",3.05156e+07},
      {8,"177Hf",8.34644e+06},
      {10,"179Hf",-4.99957e+06},
      {8,"181Ta",3.20754e+07},
      {2,"183W",  1.1131e7},
      {6,"185Re",6.02576e+07},
      {6,"187Re",6.08654e+07},
      {2,"187Os", 0.616e7},
      {4,"189Os",2.07725e+07},
      {4,"191Ir",4.59576e+06},
      {4,"193Ir",5.00492e+06},
      {2,"195Pt", 5.7505e7},
      {4,"197Au",4.62517e+06},
      {2,"199Hg", 4.769e7},
      {4,"201Hg",-1.76884e+07},
      {2,"203Tl",1.53083e+08},
      {2,"205Tl",15.438e7},
      {2,"207Pb", 5.5968e7},
      {10,"209Bi",4.29871e+07},
      {8,"235U",4.7883e+06}
    };

    // Nuclear gyromagnetic ratios, source IUPAC Recommendations 2001, Robin K. Harris et al via Simpson isotopes.c

 const nucleus_properties nucleus_props_Harris[harris_internal_nucs]={
  {0,NULL,0},
{2,"1H", 26.7522128e7},
{2,"13C", 6.7282840e7},
{3,"2H", 4.1066279e7},
{3,"14N", 1.9337792e7},
{2,"15N", -2.7126180e7},
{2,"3H", 28.5349779e7},
{2,"3He", -20.3801580e7},
{3,"6Li", 3.9371709e7},
{4,"7Li", 10.3977013e7},
{4,"9Be", 3.7596660e7},
{7,"10B", 2.8746786e7},
{4,"11B", 8.5847044e7},
{6,"17O", -3.6280800e7},
{2,"19F", 25.1814800e7},
{4,"21Ne", -2.1130800e7},
{4,"23Na", 7.0808493e7},
{6,"25Mg", -1.6388700e7},
{6,"27Al", 6.9762715e7},
{2,"29Si", -5.3190000e7},
{2,"31P", 10.8394000e7},
{4,"33S", 2.0556850e7},
{4,"35Cl", 2.6241980e7},
{4,"37Cl", 2.1843680e7},
{4,"39K", 1.2500608e7},
{9,"40K", -1.5542854e7},
{4,"41K", 0.6860681e7},
{8,"43Ca", -1.8030690e7},
{8,"45Sc", 6.5087973e7},
{6,"47Ti", -1.5105000e7},
{8,"49Ti", -1.5109500e7},
{13,"50V", 2.6706490e7},
{8,"51V", 7.0455117e7},
{4,"53Cr", -1.5152000e7},
{6,"55Mn", 6.6452546e7},
{2,"57Fe", 0.8680620e7},
{8,"59Co", 6.3320000e7},
{4,"61Ni", -2.3948000e7},
{4,"63Cu", 7.1117890e7},
{4,"65Cu", 7.6043500e7},
{6,"67Zn", 1.6766880e7},
{4,"69Ga", 6.4388550e7},
{4,"71Ga", 8.1811710e7},
{10,"73Ge", -0.9360303e7},
{4,"75As", 4.5961630e7},
{2,"77Se", 5.1253857e7},
{4,"79Br", 6.7256160e7},
{4,"81Br", 7.2497760e7},
{10,"83Kr", -1.0331000e7},
{6,"85Rb", 2.5927050e7},
{4,"87Rb", 8.7864000e7},
{10,"87Sr", -1.1639376e7},
{2,"89Y", -1.3162791e7},
{6,"91Zr", -2.4974300e7},
{10,"93Nb", 6.5674000e7},
{6,"95Mo", -1.7510000e7},
{6,"97Mo", -1.7880000e7},
{10,"99Tc", 6.0460000e7},
{6,"99Ru", -1.2290000e7},
{6,"101Ru", -1.3770000e7},
{2,"103Rh", -0.8468000e7},
{6,"105Pd", -1.2300000e7},
{2,"107Ag", -1.0889181e7},
{2,"109Ag", -1.2518634e7},
{2,"111Cd", -5.6983131e7},
{2,"113Cd", -5.9609153e7},
{10,"113In", 5.8845000e7},
{10,"115In", 5.8972000e7},
{2,"115Sn", -8.8013000e7},
{2,"117Sn", -9.5887900e7},
{2,"119Sn", -10.0317000e7},
{6,"121Sb", 6.4435000e7},
{8,"123Sb", 3.4892000e7},
{2,"123Te", -7.0590980e7},
{2,"125Te", -8.5108404e7},
{6,"127I", 5.3895730e7},
{2,"129Xe", -7.4521030e7},
{4,"131Xe", 2.2090760e7},
{8,"133Cs", 3.5332539e7},
{4,"135Ba", 2.6575000e7},
{4,"137Ba", 2.9929500e7},
{11,"138La", 3.5572390e7},
{8,"139La", 3.8083318e7},
{6,"141Pr", 8.1907000e7},
{8,"143Nd", -1.4570000e7},
{8,"145Nd", -0.8980000e7},
{8,"147Pm", 3.6130000e7},
{8,"147Sm", -1.1150000e7},
{8,"149Sm", -0.9192000e7},
{6,"151Eu", 6.6510000e7},
{6,"153Eu", 2.9369000e7},
{4,"155Gd", -0.8213200e7},
{4,"157Gd", -1.0769000e7},
{4,"159Tb", 6.4310000e7},
{6,"161Dy", -0.9201000e7},
{6,"163Dy", 1.2890000e7},
{8,"165Ho", 5.7100000e7},
{8,"167Er", -0.7715700e7},
{2,"169Tm", -2.2180000e7},
{2,"171Yb", 4.7288000e7},
{6,"173Yb", -1.3025000e7},
{8,"175Lu", 3.0552000e7},
{15,"176Lu", 2.1684000e7},
{8,"177Hf", 1.0860000e7},
{10,"179Hf", -0.6821000e7},
{8,"181Ta", 3.2438000e7},
{2,"183W", 1.1282403e7},
{6,"185Re", 6.1057000e7},
{6,"187Re", 6.1682000e7},
{2,"187Os", 0.6192895e7},
{4,"189Os", 2.1071300e7},
{4,"191Ir", 0.4812000e7},
{4,"193Ir", 0.5227000e7},
{2,"195Pt", 5.8385000e7},
{4,"197Au", 0.4730600e7},
{2,"199Hg", 4.8457916e7},
{4,"201Hg", -1.7887690e7},
{2,"203Tl", 15.5393338e7},
{2,"205Tl", 15.6921808e7},
{2,"207Pb", 5.5804600e7},
{10,"209Bi", 4.5444000e7},
{2,"209Po", 7.4000000e7},
{4,"227Ac", 3.5000000e7},
{6,"229Th", 0.4000000e7},
{4,"231Pa", 3.2100000e7},
{8,"235U", -0.5200000e7},
{6,"237Np", 3.1000000e7},
{2,"239Pu", 0.9720000e7},
};

    //! return ::nucleus_properties for nucleus type \a nuc
    const nucleus_properties& get_properties(size_t nuc)
    {
      if (nuc==NULL_NUCLEUS)
	throw InvalidParameter("get_properties: can't request properties of generic (NULL) nucleus");
      if (nuc<internal_nucs) {
	const nucleus_properties* propstable(get_nucleus_properties());
	return *(propstable+nuc);
      }
      nuc-=internal_nucs;
      if (nuc>=usernucs.size())
	throw Undefined("get_properties: nucleus number is undefined");
      return usernucs(nuc).props;
    }    
    
    //! throw ::Undefined exception if \a nuc is not defined
    inline void validate_nuc(size_t nuc) { (void)get_properties(nuc); }
  
  } //anon namespace

    void set_nucleus_properties(const char* table)
    {
      if (used_properties_table)
	throw Failed("set_nucleus_properties: cannot switch nucleus properties table after first use");
      size_t inucs=0;
      if (strcmp(table,"HarrisIUPAC")==0) {
	used_properties_table=nucleus_props_Harris;	
	inucs=harris_internal_nucs;
      }
      else if (strcmp(table,"HarrisMann")==0) {
	used_properties_table=nucleus_props_GAMMA;
	inucs=gamma_internal_nucs;
      }
      else
	throw Failed("set_nucleus_properties: unknown properties table (current values are HarrisIUPAC and HarrisMann");
      if (inucs>internal_nucs)
	throw InternalError("set_nucleus_properties: internal_nucs is invalid");
    }

  
    //! checked store for nucleus id
    nuclei_spec::nuclei_spec(size_t value_)
      : value(value_) 
    {
      if (value_!=NULL_NUCLEUS)
	validate_nuc(value_);
    }

    /** \note The search is very crude, but table is ordered with
	common nuclei at the start and creating a hashed/sorted
	table would add to start up time */

    size_t labeltonuc(const char* lab)
    {
      const usernucs_map_t& usernucs_map(getusernucs_map());

      if (!usernucs_map.empty()) {
	const usernucs_map_t::const_iterator userp(usernucs_map.find(lab));
	if (userp!=usernucs_map.end())
	  return userp->second;
      }
      const nucleus_properties* propstable(get_nucleus_properties());
      for (size_t i=1;i<internal_nucs;i++) {
	if (!strcmp(lab,(propstable[i].label)))
	  return i;
      }      
      throw Failed("labeltonuc: unknown nucleus name");
    }

    const char *nuctolabel(size_t nuc)
    {
      return get_properties(nuc).label;
    }

    /**
       ::Failed exception is thrown if definition failed (name already used) */
    void define_nucleus(const char* label, size_t deg, double gamma)
    {
      usernucs_map_t& usernucs_map(getusernucs_map());
      usernucs.push_back(user_nuc(label,deg,gamma));
      usernucs_map[label]=internal_nucs+usernucs.size()-1;
    }

  size_t maximumnuc() { return internal_nucs+usernucs.size(); }

double spin::gamma() const
{
  return get_properties(nucleus_).gamma;
}

size_t spin::deg() const
{
  return get_properties(nucleus_).deg;
}

float spin::qn() const
{
  return 0.5*(get_properties(nucleus_).deg-1);
}

void spin::nucleus(size_t nuc, bool restrictv)
{
  const size_t m=get_properties(nuc).deg;
  if (restrictv) {
    if (m & 1)
      throw InvalidParameter("spin: central transition restriction is only valid for half-integer nuclei");
    if (m==2)
      restrictv=false; //!< silently remove restrict flag to ensure only one "type" of spin-1/2
  }
  restrict_=restrictv;
  //  validate_nuc(nuc);
  //const size_t m=nucleus_props[nuc].deg;
  //if (halfonly_ && (nucleus_props[nuc].deg!=2)) throw Failed("spin::nucleus: restricted to spin-1/2");
  nucleus_=nuc;
  //  _deg=m;
}

// int logtwo(long n)
// {
//   if (n<=0) throw Failed("logtwo input must be >0");
//   int ret=0;
//   while (!(n & 1)) {
//     ret++;
//     n>>=1;
//   }
//   if (n>1) throw Failed("logtwo input not power of 2");
//   return ret;
// }

size_t degeneracy(float v)
{
  v=2*v+1;
  size_t res=int(v+0.5);
  if (res<1 || fabs(v-res)>1e-4)
    throw InvalidParameter("degeneracy: Spin must be integer or half integer");
  return res;
}

void basespin_system::nucleus(size_t n, const spin& nuc)
{
  List<spin>::operator()(n)=nuc;
}

size_t basespin_system::nspins(nuclei_spec which) const
{
  const size_t nuc=which();

  size_t res=0;
  const List<spin>& spinlist(*this);
  for (size_t i=nspins();i--;) {
    if (spinlist(i).nucleus()==nuc)
      res++;
  }
  return res;
}

size_t basespin_system::size() const
{
  size_t m=1;
  const List<spin>& spinlist(*this);
  for (size_t i=nspins();i--;) {
    m*=spinlist(i).effective_deg();
    if (m>NMRSIM_MAX_STATES) { //!< trap here - size() will always be called, and prevent more obscure errors later
      char errmess[256];
      snprintf(errmess,sizeof(errmess),"Spin system is ridiculously large - contains more than %i states!",NMRSIM_MAX_STATES);
      throw Failed(errmess);
    }
  }
  return m;
}

basespin_system::basespin_system(int nspin, const spin& def)
  : List<spin>(nspin,def) {}

void basespin_system::isotope(size_t from,size_t to, const spin& nuc)
{
  if (from>to || to>nspins())
    throw BadIndex("basespin_system::isotope");
  //  const size_t nuc=nucspec();
  for (size_t i=from;i<=to;i++)
    nucleus(i,nuc);
}
  
double gamma(const nuclei_spec& nspec)
{
  return get_properties(nspec()).gamma;
}

std::ostream& operator << (std::ostream& ostr,const spin& a)
{
  ostr << a.isotope();
  if (a.restrict_)
    ostr << "(c)"; //!< central transition indicator
  return ostr;
}

std::ostream& operator << (std::ostream& ostr,const basespin_system& a)
{
  for (size_t i=0;i<a.nspins();i++)
    ostr << a(i) << ' ';
  return ostr;
}

void basespin_system::permutation_vectorH(BaseList<state_t>& ninds, const Permutation& permvec) const
{
  validate_diagonal();
  size_t n=rows();
  if (n!=ninds.length())
    throw Mismatch("permutation_vectorH");
  for (;n--;) 
    ninds(n)=rawpermute(n,permvec);
}

bool Permutation::isvalid(const BaseList<size_t>& permvec) 
{
  const size_t n(permvec.size());
  ScratchList<bool> used(n,false);
  for (size_t i=n;i--;) {
    const size_t which(permvec(i));
    if ((which>=n) || (used(which)))
      return false;
    
    used(which)=true;
  }
  return true;
}

size_t countnonzero(state_t n)
{
  size_t nz=0;
  while (n) {
    if (n & 1)
      nz++;
    n>>=1;
  }
  return nz;
}

template<class T> void spinorder_spec::add_(size_t n, state_t which, const BaseList<T>& ordersv)
{
  registerorders(which,compressorders(n,ordersv));
}
	
void spinorder_spec::create(size_t n, const BaseList<size_t>& ordersv)
{
  init(n,false);
  const state_t whichmask=(1<<n)-1;
  add_(n,whichmask,ordersv);
}

void spinorder_spec::create(size_t n, const BaseList<int>& ordersv)
{
  init(n,true);
  const state_t whichmask=(1<<n)-1;
  add_(n,whichmask,ordersv);
}

void spinorder_spec::init(size_t nspinsv, bool iscoherv)
{
  nspins=nspinsv;
  iscoher=iscoherv;
  totalspinmask=0;
  //  totalordermask=0; 
  //  posmask=(1<<nspins)-1;
}

state_t spinorder_spec::coherencemask(int coher) const
{
  if (coher<0) {
    if (coher<-int(nspins))
      throw InvalidParameter("coherencemask: coherence is out of range");
    return state_t(1) << (nspins-1-coher);
  }
  else {
    if (coher>nspins)
      throw InvalidParameter("coherencemask: coherence is out of range");
    return state_t(1) << coher;
  }
}

void spinorder_spec::dumpcoherences(std::ostream& ostr, state_t n) const
{
  bool done=false;
  const int inspins(nspins);
  for (int co=-inspins;co<=inspins;co++) {
    if (n & coherencemask(co)) {
      if (done)
	ostr << ',';
      else
	done=true;
      ostr << co;
    }
  }
}

void spinorder_spec::dumpbits(std::ostream& ostr, state_t n, size_t start, bool rev) const
{  
  bool done=false;
  if (rev) {
    state_t startmask=maskelement(nspins,0);
    while (startmask) {
      if (n & startmask) {
	if (done)
	  ostr << ',';
	else
	  done=true;
	ostr << start;
      }
      start++;
      startmask >>=1;
    }
  }
  else {
    while (n) {
      if (n & 1) {
	if (done)
	  ostr << ',';
	else
	  done=true;
	ostr << start;
      }    
      n >>= 1;
      start++;
    }
  }
}

void spinorder_spec::print(std::ostream& ostr) const
{
  if (specs.empty())
    ostr << "<empty>";
  else {
    for (size_t i=0;i<specs.size();i++) {
      if (i)
	ostr << "  ";
      const subspec_t& curspec(specs(i));
      ostr << "Spins: ";
      dumpbits(ostr,curspec.first,cmatrix_ostream_controller(ostr).indexbase,true);
      ostr << "  " << (iscoher ? "Coherences" : "Orders") << ": ";
      if (curspec.second) {
	if (iscoher)
	  dumpcoherences(ostr,curspec.second);
	else
	  dumpbits(ostr,curspec.second);
      }
      else
	ostr << "none";
    }
  }
}

bool spinorder_spec::overlappingranges() const
{
  state_t sofar=0;
  for (size_t i=specs.size();i--;) {
    if (sofar & specs(i).first)
      return true;
    sofar |= specs(i).first;
  }
  return false;
}

  // if (max==0)
  //   throw InvalidParameter("spinorder_spec: no spins involved!");

  //  if (whichmask & totalspinmask)
  //  throw Failed("Invalid specification of spin order - repeated spin numbers");


state_t spinorder_spec::compressorders(size_t max, const BaseList<size_t>& ordersv) const
{
  if (iscoher)
    throw Failed("spinorder_spec: trying to add spin orders to coherence order specification");

  state_t orders=0;
  char buf[100];
  for (size_t i=ordersv.size();i--;) {
    const size_t n=ordersv(i);
    if (n>max) {
      snprintf(buf,sizeof(buf),"Invalid specification of spin order - order (%" LCM_PRI_SIZE_T_MODIFIER "u) exceeds number of spins selected (%" LCM_PRI_SIZE_T_MODIFIER "u)",n,max);
      throw Failed(buf);
    }
    const state_t mask=1<<n;
    //    if (orders & (mask | totalordermask)) {
    if (orders & mask) {
      snprintf(buf,sizeof(buf),"Invalid specification of spin order - repeated order (%" LCM_PRI_SIZE_T_MODIFIER "u)",n);
      throw Failed(buf);
    } 
    orders|=mask;
  }
  return orders;
}

state_t spinorder_spec::compressorders(size_t max, const BaseList<int>& ordersv) const
{
  if (!iscoher)
    throw Failed("spinorder_spec: trying to add coherences to spin order specification");

  state_t orders=0;
  char buf[100];
  const int imax(max);
  for (size_t i=ordersv.size();i--;) {
    const int n=ordersv(i);
    if ((n>imax) || (n<-imax)) {
      snprintf(buf,sizeof(buf),"Invalid specification of coherence order - order (%i) exceeds number of spins selected (%" LCM_PRI_SIZE_T_MODIFIER "u)",n,max);
      throw Failed(buf);
    }
    const state_t mask=coherencemask(n);
    //    if (orders & (mask | totalordermask)) {
    if (orders & mask) {
      snprintf(buf,sizeof(buf),"Invalid specification of coherence order - repeated order (%i)",n);
      throw Failed(buf);
    } 
    orders|=mask;
  }
  return orders;
}

void spinorder_spec::registerorders(state_t whichmask, state_t orders)
{
  for (size_t i=specs.size();i--;) {
    if (specs(i).first==whichmask) {
      if (specs(i).second==orders)
	return; //!< already specified (ignore)      
      throw Failed("Invalid specification of spin / coherence order - same range of spins specified twice with inconsistent coherence / order set");
    }
  }

  //  totalordermask|=orders;  
  totalspinmask|=whichmask;
  specs.push_back(subspec_t(whichmask,orders));
}

state_t spinorder_spec::compress_spinindices(size_t nspinsv,const BaseList<size_t>& whichv, state_t ltotalspinmask)
{
  if (whichv.empty())
    throw InvalidParameter("spinorder_spec: no spins involved!");
  char buf[100];

  state_t which=0;
  for (size_t i=whichv.size();i--;) {
    const size_t n=whichv(i);
    if (n>=nspinsv) {
      snprintf(buf,sizeof(buf),"Invalid specification of spin order - spin index (%" LCM_PRI_SIZE_T_MODIFIER "u) out of range (%" LCM_PRI_SIZE_T_MODIFIER "u)",n+cmatrix_ostream_controller(std::cerr).indexbase,nspinsv);
      throw Failed(buf);
    }
    const state_t mask=maskelement(nspinsv,whichv(i));
    if (which & (mask | ltotalspinmask)) {
      snprintf(buf,sizeof(buf),"Invalid specification of spin order - repeated spin index (%" LCM_PRI_SIZE_T_MODIFIER "u)",n+cmatrix_ostream_controller(std::cerr).indexbase);
      throw Failed(buf);
    } 
    which|=mask;
  }
  return which;
}

void spinorder_spec::add(const BaseList<size_t>& whichv, const BaseList<size_t>& ordersv)
{
  const state_t which=compress_spinindices(nspins,whichv,totalspinmask);
  add_(whichv.size(),which,ordersv);
}

void spinorder_spec::add(const BaseList<size_t>& whichv, const BaseList<int>& ordersv)
{
  const state_t which=compress_spinindices(nspins,whichv,totalspinmask);
  add_(whichv.size(),which,ordersv);
}

bool spinorder_spec::isactive(size_t n) const
{
  const state_t mask=maskelement(nspins,n);
  if (mask==0) {
    char buf[100];
    snprintf(buf,sizeof(buf),"spinorder_spec: spin index wildly out of range (%" LCM_PRI_SIZE_T_MODIFIER "u)",n);
    throw InvalidParameter(buf);
  }
  return totalspinmask & mask;
}

bool spinorder_spec::operator()(state_t bra, state_t ket) const
{
  if (specs.empty())
    return false;
  
  const state_t diff=(bra ^ ket);
  if (diff & ~totalspinmask)
    return false; //!< other spins not same state
    
  for (size_t i=specs.size();i--;) {
    const subspec_t curspec(specs(i));
    state_t selmask=0;
    size_t count=0;
    if (iscoher) {
      const state_t lbra = bra & curspec.first;
      const state_t lket = ket & curspec.first;
      const int coherket =int(countnonzero(lket));
      const int coherbra =int(countnonzero(lbra));
      const int coher = coherket - coherbra;
#ifndef NDEBUG
      std::cout << ket << " (" << coherket << ") - " << bra << " (" << coherbra << ") against spec " << i;
#endif
      selmask = coherencemask(coher);
    }
    else {
      state_t ldiff=diff & curspec.first;
      while (ldiff) {
	count+=(ldiff & 1);
	ldiff >>=1;
      }
      selmask= 1 << count;
#ifndef NDEBUG
      std::cout << bra << ' ' << ket << " against spec " << i << " (" << count << ")";
#endif
    }
    const bool issel=selmask & curspec.second;
#ifndef NDEBUG
    std::cout << ": " << issel << '\n';
#endif
    if (!issel)
      return false;
  }
  return true;
}

bool spinorder_spec::iscoherence(double coher, state_t sel) const
{  
  int icoher=(coher>0) ? int(coher+0.5) : int(coher-0.5); 
  if (fabs(icoher-coher)>0.2)
    throw InvalidParameter("iscoherence: argument is not close to integer");
  return sel & coherencemask(icoher);
}

Permutation::Permutation(const BaseList<size_t>& permvec)
  : DynamicList<size_t>(permvec,mxflag::normal), order_(0) 
{
  if (!isvalid(*this))
    throw Failed("Invalid permutation vector");
}

Permutation::Permutation(const char* permvec)
  : DynamicList<size_t>(strlen(permvec),mxflag::normal), order_(0) 
{
  for (size_t m=size();m--;) {
    if (!isdigit(permvec[m]))
      throw InvalidParameter("Permutation()");
    (*this)(m)=permvec[m]-'0';
  }
  if (!isvalid(*this))
    throw Failed("Invalid permutation vector");
}

bool Permutation::isidentity(const BaseList<size_t>& a)
{
  for (size_t i=a.size();i--;) {
    if (a(i)!=i)
      return false;
  }
  return true;
}

size_t Permutation::order() const
{
  if (order_>0)
    return order_;
  const size_t n(this->size());
  List<size_t> lorder(n);
  List<size_t> ordertmp;
  size_t i;
  for (i=n;i--;)
    lorder(i)=i;
  for (i=1;i<=n;i++) {
    apply(ordertmp,lorder);
    lorder.swap(ordertmp);
    if (isidentity(lorder))
      return i;
  }
  throw Failed("Permutation::order()");
}

std::ostream& operator<< (std::ostream& ostr, const Permutation& a)
{
  const size_t n=a.size();
  if (n==0)
    return ostr << "(none)";
  const char* sep = (n>9) ? " " : "";
  for (size_t i=0;i<n;i++)
    ostr << a(i) << sep;
  return ostr;
}
    
rmatrix basespin_system::permutation_matrixH(const Permutation& permvec) const
{
  const size_t n=rows();
  rmatrix d(mxflag::temporary);
  d.create(n,n,0.0);
  for (size_t i=n;i--;)
    d(rawpermute(i,permvec),i)=1.0;
  return d;
}

rmatrix basespin_system::permutation_matrixL(const Permutation& permvec) const
{
  const size_t n=rows();
  ScratchList<state_t> nstates(n);
  permutation_vectorH(nstates,permvec);
  const state_t Ldim=n*n;
  rmatrix d(mxflag::temporary);
  d.create(Ldim,Ldim,0.0);
  size_t ptr=0;
  for (size_t j=0;j<n;j++) {
    const state_t nLbase=n*nstates(j);
    for (size_t k=0;k<n;k++)
      d(nLbase+nstates(k),ptr++)=1.0;
  }
  return d;
}

void basespin_system::permutation_vectorL(BaseList<state_t>& nLstates,const Permutation& permvec) const
{
  const size_t n=rows();
  if (nLstates.length()!=n*n)
    throw Mismatch("permutation_vectorL");
  ScratchList<state_t> nstates(n);
  permutation_vectorH(nstates,permvec);
  size_t ptr=0;
  for (size_t j=0;j<n;j++) {
    const state_t nLbase=n*nstates(j);
    for (size_t k=0;k<n;k++)
      nLstates(ptr++)=nLbase+nstates(k);
  }
}

size_t basespin_system::dim_range(size_t start,size_t endv) const
{
  if (endv<start)
    throw BadIndex("basespin_system::dim_range");

  size_t dim=1;
  const basespin_system& sys(*this);
  for (size_t i=start;i<endv;i++)
    dim*=sys(i).effective_deg();
  return dim;
}

void basespin_system::mla_expandop(BaseList<double>& Z, size_t index, const BaseList<double>& z) const
{
  if (nspins()==1) {
    Z+=z;
    return;
  }
  const size_t zl=z.length();
  if (zl!=(*this)(index).effective_deg())
    throw Mismatch("basespin_system::mla_expandop");

  const size_t dimbef=dim_range(0,index);
  const size_t dimaft=dim_range(index+1,nspins());
  
  if (Z.length()!=dimaft*dimbef*zl)
    throw Mismatch("basespin_system::mla_expandop");

  size_t count=0;
  for (size_t i=dimbef;i--;) {
    for (size_t j=0;j<zl;j++) {
      const double el=z(j);
      for (size_t k=dimaft;k--;)
	Z(count++)+=el;
    }
  }
}

template<typename T> void expandop_(const basespin_system& sys, Matrix<T>& dest,size_t index,const Matrix<T>& a)
{
  if (!issquare(a))
    throw NotSquare("basespin_system::expandop");
  if (a.rows()!=sys(index).effective_deg())
    throw Mismatch("basespin_system::expandop");

  const size_t dimaft=sys.dim_range(index+1,sys.nspins());
  const size_t dimbef=sys.dim_range(0,index);

  if (dimbef>1) {
    kronecker(dest,dimbef,a);
    if (dimaft>1) {
      Matrix<T> tmp;
      kronecker(tmp,dest,dimaft);
      dest.swap(tmp);
    }
  }
  else
    kronecker(dest,a,dimaft);
}

void basespin_system::expandop(cmatrix& dest, size_t index, const cmatrix& a) const
{
  expandop_(*this,dest,index,a);
}

void basespin_system::expandop_tmp(cmatrix& dest, size_t index, cmatrix& a) const
{
  if (nspins()==1)
    dest.swap(a);
  else
    expandop_(*this,dest,index,a);
}

void basespin_system::expandop(rmatrix& dest, size_t index, const rmatrix& a) const
{
  expandop_(*this,dest,index,a);
}

void basespin_system::mla_F(cmatrix& dest,double scale,size_t nuc,char op) const
{
  const basespin_system& sys=*this;
  bool found=false;
  for (size_t i=nspins();i--;) {
    if (nuc==NULL_NUCLEUS || sys(i).nucleus()==nuc) {
      found=true;
      mla_I(dest,scale,i,op);
    }
  }
  if (!found) 
    throw_notpresent(nuc);
}

char basespin_system::squeezeST(size_t r, size_t c)
{
  switch (r) {
  case 0:
    switch (c) {
    case 0:
      return 'a';
    case 1:
      return '+';
    }
    break;
  case 1:
    switch (c) {
    case 0:
      return '-';
    case 1:
      return 'b';
    }
  }
  throw InvalidParameter("squeezeST: invalid spin-1/2 single transition operator");
}

void basespin_system::mla_Iz(List<double>& dest,double scale,size_t which) const 
{
  if (dest.length()==0) {
    dest.create(rows()); 
    dest=0.0;
  } 
  mla_Iz( static_cast< BaseList<double>& >(dest),scale,which);
}

void basespin_system::mla_Iz(List<double>& dest,double scale,size_t which1, size_t which2) const 
{
  if (dest.length()==0) {
    dest.create(rows());
    dest=0.0;
  }
  mla_Iz( static_cast< BaseList<double>& >(dest),scale,which1,which2);
}

void basespin_system::mla_Fz(List<double>& dest,double scale,size_t which) const 
{
  if (dest.length()==0) {
    dest.create(rows());
    dest=0.0;
  }
  mla_Fz( static_cast< BaseList<double>& >(dest),scale,which);
}

void basespin_system::mla_I(cmatrix& dest, double scale, size_t n, size_t r, size_t c) const
{
  if ((*this)(n).effective_deg()==2)
    mla_I(dest,scale,n,squeezeST(r,c));
  else
    mla_ST(dest,scale,n,r,c);
}

void basespin_system::mla_F(cmatrix& dest,double scale,size_t nuc,size_t r, size_t c) const
{
  const basespin_system& sys=*this;
  if (spin(nuc).effective_deg()==2) { //if spin-1/2, convert ST operator to standard
    mla_F(dest,scale,nuc,squeezeST(r,c));
    return;
  }
  bool found=false;
  if (r==c) {
    List<double> tmpd;
    for (size_t i=nspins();i--;) {
      if ((nuc==NULL_NUCLEUS) || (sys(i).nucleus()==nuc)) {
	found=true;
	mla_ST(tmpd,scale,i,r);
      }
    }
    if (!dest)
      full(dest,tmpd);
    else {
      cmatrix tmp;
      full(tmp,tmpd);
      dest+=tmp;
    }
  }
  else {
    for (size_t i=nspins();i--;) {
      if (nuc==NULL_NUCLEUS || sys(i).nucleus()==nuc) {
	mla_ST(dest,scale,i,r,c);
	found=true;
      }
    }
  }
  if (!found)
    throw_notpresent(nuc);
}

void basespin_system::mla_Fz(BaseList<double> &dest,double scale,size_t nuc) const
{
  const basespin_system& sys=*this;
  bool found=false;
  for (size_t i=nspins();i--;) {
    if (nuc==NULL_NUCLEUS || sys(i).nucleus()==nuc) {
      found=true;
      mla_Iz(dest,scale,i);
    }
  }
  if (!found)
    throw_notpresent(nuc);
}

cmatrix I(const basespin_system &sys,size_t n,char op)
{
  cmatrix A(mxflag::temporary);
  sys.mla_I(A,1.0,n,op);
  return A;
}

cmatrix I(const basespin_system &sys,size_t n,size_t r,size_t c)
{
  cmatrix A(mxflag::temporary);
  sys.mla_ST(A,1.0,n,r,c);
  return A;
}

cmatrix I(const basespin_system &sys,const BaseList<char> &ops)
{
  cmatrix A(mxflag::temporary);
  sys.mla_I(A,1.0,ops);
  return A;
}

List<double> diag_Iz(const basespin_system &sys,size_t n)
{
  List<double> Zop(mxflag::temporary);
  sys.mla_Iz(Zop,1.0,n);
  return Zop;
}

cmatrix I(const basespin_system &sys,size_t n,const BaseList<char>& ops)
{
  cmatrix A(mxflag::temporary);
  sys.mla_I(A,1.0,n,ops);
  return A;
}

cmatrix I(const basespin_system& sys,size_t i,char opi,size_t j,char opj)
{
  cmatrix A(mxflag::temporary);
  sys.mla_I(A,1.0,i,opi,j,opj);
  return A;
}

List<double> diag_Iz(const basespin_system& sys,size_t i,size_t j)
{
  if (!sys.isdiagonal())
    throw Failed("Only valid for diagonal blocks");
  List<double> tmp(mxflag::temporary);
  sys.mla_Iz(tmp,1.0,i,j);
  return tmp;
}

List<double> diag_Fz(const basespin_system& sys,nuclei_spec nuc)
{
  if (!sys.isdiagonal())
    throw Failed("Only valid for diagonal blocks");
  List<double> tmp(mxflag::temporary);
  sys.mla_Fz(tmp,1.0,nuc());
  return tmp;
}

cmatrix F(const basespin_system& a,nuclei_spec nucspec,char op)
{ 
  cmatrix d(mxflag::temporary);
  const size_t nuc=nucspec();
  bool found=false;
  for (size_t i=a.nspins();i--;) {
    if ((nuc==NULL_NUCLEUS) || (a(i).nucleus()==nuc)) {
      found=true;
      a.mla_I(d,1.0,i,op);
    }
  }
  if (!found)
    throw_notpresent(nuc);
  return d;
}

cmatrix F(const basespin_system& a,nuclei_spec nucspec,size_t r,size_t c)
{ 
  cmatrix d(mxflag::temporary);
  const size_t nuc=nucspec();
  bool found=false;
  for (size_t i=a.nspins();i--;) {
    if ((nuc==NULL_NUCLEUS) || (a(i).nucleus()==nuc)) {
      found=true;
      a.mla_ST(d,1.0,i,r,c);
    }
  }
  if (!found)
    throw_notpresent(nuc);
  return d;
}

}//namespace libcmatrix
