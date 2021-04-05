#include "powder.h"
#include <cmath>

namespace libcmatrix {

  using ::std::sin;

  inline void validate_scale(const Euler& dest, double scale, size_t which)
  {
    if (scale<0.0) {
      std::cerr << "-ve scaling factor on powder orientation " << which << ": " << dest << ' ' << scale << std::endl;
      throw Failed("PowderMethod: negative powder scaling factor");
    }
  }

  void PowderMethod::orientation(Euler& dest, double& scale, size_t which) const
  { if (which>=max) 
      throw BadIndex("orientation");
    else {
      index_to_orientation(dest,scale,which); 
      validate_scale(dest,scale,which);
    }
  }
  
  void PowderMethod::create(size_t maxv)
  { 
    if (maxv<1)
      throw InvalidParameter("PowderMethod::create");
    max=maxv;
    current=0;
  }

bool PowderMethod::next(Euler& dest,double& weight)
{
  current++;
  if (current<=max) {
    index_to_orientation(dest,weight,current-1);
    validate_scale(dest,weight,current-1);
    return true;
  }
  if (current>max+1)
    throw Failed("PowderMethod::next");
  return false;
}

 //  range_t PowderMethod::range_type(const char* range)
//   {
//     if (strcmp(range,"sphere")==0)
//       return sphere;

//     if (strcmp(range,"hemisphere")==0)
//       return hemisphere;

//     if (strcmp(range,"octant")==0)
//       return octant;

//     throw Failed("PowderMethod::range_type: unknown range qualifier");
//   }

double PowderMethod::alphamax(range_t range)
{
  switch (range) {
  case sphere:
    return 2*M_PI;
  case hemisphere:
    return 2*M_PI;
  case octant:
    return M_PI/2;
  }
  throw InternalError("PowderMethod::alphamax");
}

double PowderMethod::betamax(range_t range)
{
  switch (range) {
  case sphere:
    return M_PI;
  case hemisphere:
    return M_PI/2;
  case octant:
    return M_PI/2;
  }
  throw InternalError("betamax");
}

  void PowderMethod::include_parameters(double& step, double& offset, double range, size_t steps, include_t type)
  {
    if (steps<2)
      throw InvalidParameter("PowderMethod: too few powder steps (<2)");
    step=range/((type==both) ? (steps-1) : steps);
    switch (type) {
    case start: case both:
      offset=0.0;
      break;
    case middle:
      offset=0.5;
      break;
    case end:
      offset=1.0;
      break;
    }
  }

void PowderMethod::alpha_parameters(double& step, double& offset, size_t steps, range_t range, include_t inctype)
{
  if (steps==1)
    offset=step=0.0;
  else
    include_parameters(step,offset,alphamax(range),steps,(range==octant) ? inctype : start);
}

  double PlanarGrid::getbeta(size_t whichb) const
  { 
    return (whichb+boffset)*bstep;
  }

  PlanarGrid::PlanarGrid(size_t astepsv, size_t bstepsv, range_t range, include_t inctype)
    : PowderMethod(astepsv*bstepsv,2), asteps(astepsv), bsteps(bstepsv)
{
  if (asteps<1 || bsteps<1)
    throw InvalidParameter("PlanarGrid()");

  alpha_parameters(astep,aoffset,asteps,range,inctype);
  issingle=(bsteps==1);

  if (issingle) {
    angles_=(asteps==1) ? 0 : 1; //correct no. of averaging angles
    normal=1.0/orientations();
  }
  else {
    include_parameters(bstep,boffset,betamax(range),bsteps,inctype);
    if (asteps==1)
      angles_=1;
    double total=0;
    for (size_t i=bsteps;i--;)
      total+=sin(getbeta(i));
    normal=1.0/(asteps*total);
  }
}

void PlanarGrid::index_to_orientation(Euler& dest,double &w,size_t which) const
{
  if (issingle) {
    if (asteps>1)
      dest.alpha=astep*(which+aoffset);
    if (bsteps>1)
      dest.beta=bstep*(which+boffset);
    w=normal;
  }
  else {
    dest.alpha=astep*(aoffset+(which % asteps));
    dest.beta=getbeta(which / asteps);
    w=normal*sin(dest.beta);
  }
}

  SphericalGrid::SphericalGrid(size_t astepsv,size_t bstepsv, range_t range, include_t inctype)
    : PowderMethod(astepsv*bstepsv,2)
  {
    asteps=astepsv;
    if (asteps!=1)
      alpha_parameters(astep,aoffset,asteps,range,inctype);
    angles_= ((astepsv==1) ? 0 : 1)+((bstepsv==1) ? 0 : 1);
    
    bsteps=bstepsv;
    if (bsteps!=1)
      include_parameters(bstep,boffset,(2/M_PI)*betamax(range),bsteps,inctype);
    
    normal=1.0/orientations();
  }

void SphericalGrid::index_to_orientation(Euler& dest,double& w,size_t which) const
{
  if (bsteps>1) {
    const size_t whichb = which / asteps;
    dest.beta=acos(1.0-(whichb+boffset)*bstep);
  }
  if (asteps>1) {
    const size_t whicha = which % asteps;
    dest.alpha=astep*(aoffset+whicha);
  }
  w=normal;
}

ExplicitSampling::ExplicitSampling(const Matrix<double>& weightings)
  : PowderMethod(weightings.rows(),(weightings.cols()==3) ? 2 : 3), //avoid exception triggered by bad input matrix
    angweights(weightings)
{
  switch (angweights.cols()) {
  case 3: case 4:
    break;
  default:
    throw Failed("ExplicitSampling: angle set must have 3 or 4 columns");
  }
  //  create(angweights.rows());
  
  double total=0;
  for (size_t i=orientations();i--;)
    total+=angweights(i,angles_);
  normal=1/total;
}

void ExplicitSampling::index_to_orientation(Euler& dest,double& w,size_t which) const
{
  dest.alpha=angweights(which,0);
  dest.beta=angweights(which,1);
  if (angles_>2)
    dest.gamma=angweights(which,2);
  w=normal*angweights(which,angles_);
}

// This is rather stupid...
int fibonacci(int n)
{
  if (n>2)
    return fibonacci(n-1)+fibonacci(n-2);
  switch (n) {
  case 1: case 2:
    return 1;
  case 0:
    return 0;
  }
  throw InvalidParameter("fibonacci");
}

size_t PlanarZCW::N_to_orientations(size_t N)
{
  return fibonacci(N+2)-1;
}

int PlanarZCW::orientations_to_N(size_t orients)
{
  if (orients==1)
    return 0;
  size_t N=1;
  size_t prevorients=fibonacci(N+1);
  size_t torients=fibonacci(N+2);
  for (;;N++) {
    if (torients-1==orients)
      return N;
    if (torients-1>orients)
      return -1;
    size_t tmp=torients;
    torients=prevorients+torients;
    prevorients=tmp;
  }
}

  PlanarZCW::PlanarZCW(size_t n,range_t range, include_t inctype)
{
  if (n==0) {
    g2=1;
    create(1);
    return;
  }
  if (n<1)
    throw InvalidParameter("PlanarZCW");
  g1=fibonacci(n);
  g2=fibonacci(n+2);
  create(g2);
  alpha_parameters(ascale,aoffset,g2,range,inctype);
  include_parameters(bscale,boffset,betamax(range),g2,inctype);

  double total=0;
  for (size_t i=orientations();i--;)
    total+=sin(getbeta(i));
  normal=1/total;
}

double PlanarZCW::getbeta(size_t which) const
{
  const LCM_LONGLONG m=which*(LCM_LONGLONG)(g1);
  if (m<0)
    throw InternalError("PlanarZCW::getbeta");
  return bscale*(boffset+(m % g2));
}
 
void PlanarZCW::index_to_orientation(Euler& dest, double& w, size_t which) const
{
  if (g2==1) {
    w=1.0;
    return;
  }
  dest.alpha=ascale*(aoffset+which);
  dest.beta=getbeta(which);
  w=normal*sin(dest.beta);
}

Warning<> WithGamma::singlecombine_warning("WithGamma: adding gamma angle integration to single angle averaging method makes little physical sense",&lcm_base_warning);

WithGamma::WithGamma(const PowderMethod& meth_, size_t gsteps_) 
  : PowderMethod(meth_.orientations()*gsteps_,meth_.angles()+1), 
    methp(meth_.clone()),
    gsteps(gsteps_), scale(2*M_PI/gsteps_)
{
  switch (angles_) {
  case 1: case 3: 
    break;
  case 2:
    singlecombine_warning.raise();
    break;
  default:
    throw Failed("WithGamma: can't add gamma angle to powder method which already includes gamma angle!");
  }
}

void WithGamma::index_to_orientation(Euler& dest, double& w, size_t which) const
{
  methp->index_to_orientation(dest,w,which/gsteps);
  dest.gamma= scale*(which % gsteps);
  w/=gsteps;
}

  SphericalZCW::SphericalZCW(size_t n,range_t range, include_t inctype)
{
  if (n==0) {
    g2=1;
    create(1);
    return;
  }
  if (n<1)
    throw InvalidParameter("SphericalZCW");
  g1=fibonacci(n);
  g2=fibonacci(n+2);
  create(g2);
  alpha_parameters(ascale,aoffset,g2,range,inctype);
  include_parameters(bscale,boffset,(2/M_PI)*betamax(range),g2,inctype);
  normal=1.0/orientations();
}
 
void SphericalZCW::index_to_orientation(Euler &dest,double &w,size_t which) const
{
  if (g2==1) {
    w=1.0;
    return;
  }
  const LCM_LONGLONG m=which*(LCM_LONGLONG)(g1);
  if (m<0)
    throw InternalError("SphericalZCW::getbeta");
  dest.alpha=ascale*(aoffset+(m % g2));
  dest.beta=acos(1.0-(which+boffset)*bscale);
  w=normal;
}

namespace {
  //! values courtesy of Matthias Ernst (ETH, Zurich)
  const size_t ZCW3_SETS=41;
  const size_t ZCW3_mod1[ZCW3_SETS] = {  10,   20,   30,   40,   50,  100,  150,  200,
                          300,  400,  500,  600,  700,  800,  900, 1000,
                         1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
                         1900, 2000, 2100, 2200, 2300, 2400, 2500, 2750,
                         3000, 3500, 4000, 4500, 5000, 6000, 7000,
                        10000, 50000};

  const size_t ZCW3_mod2[ZCW3_SETS] = {   3,    3,    7,    3,    5,   15,   35,   55,
                           89,  127,   97,  103,  145,  189,  233,  313,
                          523,  573,  447,  159,  139,  205,  321,  291,
                          671,  297,  395,  697,  527,  549,  363,  739,
                          637,  647,  403,  437, 1197,  531, 1889,
                         3189, 9027};

  const size_t ZCW3_mod3[ZCW3_SETS] = {   5,    7,   11,   15,   13,   47,   63,   81,
                          137,  187,  229,  265,  223,  257,  355,  477,
                          391,  181,  191,  553,  621,  551,  789,  481,
                          829,  479,  993,  887,  827,  841,  917, 1131,
                          933, 1069,  945, 1555, 1715, 1891, 2747,
                         4713, 14857};
  double ZCW3_normal[ZCW3_SETS] = {0.0};
}

size_t ZCW3::N_to_orientations(size_t N)
{
  if (N==0)
    return 1;
  return (N>ZCW3_SETS) ? 0 : ZCW3_mod1[N-1];
}

int ZCW3::orientations_to_N(size_t orients)
{
  static const BaseList<size_t> aslist(ZCW3_SETS,const_cast<size_t*>(ZCW3_mod1));
  if (orients==1)
    return 0;
  const size_t mod1value=orients;
  const size_t* where=std::lower_bound(aslist.begin(),aslist.end(),mod1value);
  const size_t ind=where-ZCW3_mod1+1;
  return ((where==aslist.end()) || (*where!=mod1value)) ? -1 : ind;
}

  ZCW3::ZCW3(size_t n,range_t range, include_t inctype, double gammamax)
  : PowderMethod(1,3) //! initial orientations is dummy
{
  if (gammamax<=0.0)
    throw InvalidParameter("ZCW3 (gammamax)");
  if (n==0) {
    mod1=1;
    create(1);
    return;
  }
  if ((n<1) || (n>ZCW3_SETS))
    throw InvalidParameter("ZCW3 (n)");
  n--;
  mod1=ZCW3_mod1[n];
  mod2=ZCW3_mod2[n];
  mod3=ZCW3_mod3[n];
  create(mod1);
  alpha_parameters(ascale,aoffset,mod1,range,inctype);
  include_parameters(bscale,boffset,betamax(range),mod1,inctype);
  gscale=gammamax/mod1;

  if (ZCW3_normal[n]==0.0) {
    double total=0;
    for (size_t i=orientations();i--;)
      total+=sin(getbeta(i));
    ZCW3_normal[n]=1/total;
  }
  normal=ZCW3_normal[n];
}

double ZCW3::getbeta(size_t which) const
{
  return bscale*(which+boffset);
}
 
void ZCW3::index_to_orientation(Euler& dest, double& w, size_t i) const
{
  if (mod1==1) {
    w=1.0;
    return;
  }
  const LCM_LONGLONG mod2i=i*(LCM_LONGLONG)(mod2);
  const LCM_LONGLONG mod3i=i*(LCM_LONGLONG)(mod3);
  if ((mod2i<0) || (mod3i<0))
    throw InternalError("ZCW3::index_to_orientation");
  dest.alpha=ascale*(aoffset+(mod2i % mod1));
  dest.beta=getbeta(i);
  dest.gamma=gscale*(mod3i % mod1);
  w=normal*sin(dest.beta);
}


}//namespace libcmatrix
