#undef LCM_SUPPRESS_VIEWS
#include <ctype.h>
#include <cstring>
#include "NMR.h"
#include "ScratchList.h"
#include "geometry.h"

namespace libcmatrix {

  int verbose_NMR=1;

static const double tol=1e-8;

  static double char_to_phase(const char* ax)
  {
    switch (ax[0]) {
    case 'x':
      if (ax[1]=='\0')
	return 0.0;
      break;
    case 'y':
      if (ax[1]=='\0')
	return M_PI/2;
      break;
    }
    throw InvalidParameter("phase_spec: bad phase specification");
  }
  
  phase_spec::phase_spec(char ax)
  {
    switch (ax) {
    case 'x': value=0.0; return;
    case 'y': value=M_PI/2; return;
    }
    throw InvalidParameter("phase_spec: bad phase specification");
  }
    
  phase_spec::phase_spec(const char* ax)
  {
    try {
      value = (ax[0]=='-') ? char_to_phase(ax+1)+M_PI : char_to_phase(ax);
    } catch (InvalidParameter&) {
      char* tail;
      value=strtod(ax,&tail);
      if (*tail && !isspace(*tail))
	throw InvalidParameter("phase_spec: failed to parse numeric phase");
      value*=M_PI/180.0;
    }
  }

  char operator_intersection(char a, char b)
  {
    if (a==b)
      return a;
    char test;
    if (a=='x')
      test=b;
    else {
      if (b=='x')
	test=a;
      else
	return '\0';
    }
    return ((test=='-') || (test=='+')) ? test : '\0';
  }

cmatrix spin_J(const basespin_system& sys, size_t i, size_t j)
{
  cmatrix H(mxflag::temporary);		    
  sys.mla_I(H,1.0,i,'z',j,'z');
  if (sys(i)==sys(j)) {
    sys.mla_I(H,0.5,i,'+',j,'-');
    sys.mla_I(H,0.5,i,'-',j,'+');
  }
  return H;
}

static const double dip_factor=-1e-7*hbar/(2.0*M_PI);

double dipolar_coupling(double gam1,double gam2,double r)
{
  return dip_factor*gam1*gam2/(r*r*r);
}

double dipolar_coupling_to_r(double d, double gam1, double gam2)
{
  return std::pow(std::fabs(dip_factor*gam1*gam2/d), 1.0/3);
}

static void forcerighthanded_bycolumn(Matrix<double>& V)
{
  const vector3 eigx(V(0,0),V(1,0),V(2,0));
  const vector3 eigy(V(0,1),V(1,1),V(2,1));
  const vector3 eigz(cross(eigx,eigy));
  V(0,2)=eigz.x;
  V(1,2)=eigz.y;
  V(2,2)=eigz.z;
}

static void reverseeigenvector_bycolumn(Matrix<double>& V, size_t col)
{
  if (col>=V.cols())
    throw BadIndex("reverseeigenvector_bycolumn",col,V.cols());
  for (size_t r=3;r--;) 
    V(r,col)=-V(r,col);
}

void cartesian_to_eigensystem_symmetric(double& XX, double& YY, double& ZZ, Matrix<double>& Vsorted, const Matrix<double>& A, ordering_convention_t convention, const Euler_controller& ctrl)
{
   Matrix<double> V;
   ScratchList<double,3> eigs(3);
   hermitian_eigensystem(V,eigs,A);

   if (ctrl.eigenvalue_flip_pattern) {
     const size_t revpattern=ctrl.eigenvalue_flip_pattern;
     if (revpattern & 1)
       reverseeigenvector_bycolumn(V,0);
     if (revpattern & 2)
       reverseeigenvector_bycolumn(V,1);
     if (revpattern & 4)
       reverseeigenvector_bycolumn(V,2);
   }
     
   if (ctrl.verbose>1) {
     std::cout << "Unsorted eigenvalues: " << eigs << '\n';
     std::cout << "Unsorted eigenvectors\n" << V << '\n';  
   }

   ScratchList<size_t,3> order(0U,0U,0U);
   const double isotropic=(eigs(0U)+eigs(1U)+eigs(2U))/3.0;

   switch (convention) {
   case convention_Haeberlen: case convention_NQR: {
     ScratchList<double,3> tmpeigs(3);
       for (size_t i=3;i--;)
	 tmpeigs(i)=fabs(eigs(i)-isotropic);
     indexed_sort(order,tmpeigs);     
     if (convention==convention_Haeberlen)
       std::swap(order(0U),order(1U)); //!< order is y < x < z !
   }
     break;

   case convention_IUPACIncreasing: case convention_IUPACDecreasing: 
     indexed_sort(order,eigs);     
     if (convention==convention_IUPACDecreasing)
       std::swap(order(size_t(0)),order(size_t(2))); //!< swap from increasing to decreasing order
     break;

   default:
     throw InvalidParameter("cartesian_to_eigensystem_symmetric: unrecognised ordering convention");
   }

   if (ctrl.verbose>1)
     std::cout << "Ordering of principal components: " << order << '\n';

   Vsorted=V(range(),order); //!< swapped rows for columns Feb 2013!

   if (ctrl.verbose>1)
     std::cout << "Sorted eigenvectors (rotation matrix) before RHS forcing\n" << Vsorted << '\n';

   forcerighthanded_bycolumn(Vsorted);   

   XX=eigs(order(size_t(0)));
   YY=eigs(order(size_t(1)));
   ZZ=eigs(order(size_t(2)));
}

  void cartesian_to_PAS_symmetric(double& iso11, double& aniso22, double& asym33, Euler& PAS, const Matrix<double>& A, ordering_convention_t convention, const Euler_controller& ctrl, bool doscale)
 {
   double XX,YY,ZZ;
   Matrix<double> Vsorted;

   cartesian_to_eigensystem_symmetric(XX,YY,ZZ,Vsorted,A,convention,ctrl);
   
   if (ctrl.verbose>1)
     std::cout << "Sorted eigenvalues: " << XX << ", " << YY << ", " << ZZ << '\n';
   if (ctrl.verbose)
     std::cout << "Sorted eigenvectors (rotation matrix)\n" << Vsorted << '\n';     

   PAS=DCM_to_Euler(Vsorted,ctrl);

   const double asympar=fabs((XX-YY)/ZZ);
   if (asympar<ctrl.asymmetry_tolerance) {
     const double av=0.5*(XX+YY);
     XX=YY=av;
     PAS.alpha=0.0;
     if (ctrl.verbose>1)
       std::cout << "Forcing XX=YY, alpha =0\n";
   }

   if (isorderingIUPAC(convention)) {
     iso11=XX;
     aniso22=YY;
     asym33=ZZ;
   }
   else {
     cartesian_to_anisotropy(iso11,aniso22,asym33,XX,YY,ZZ,convention);

     if (doscale) {
       static const double sqrt3halves=std::sqrt(3.0/2.0);
       static const double minussqrt3=-std::sqrt(3.0);
       iso11*=minussqrt3;
       aniso22*=sqrt3halves;
     }
   }
 }

static const double sqrt3rd=std::sqrt(1.0/3);
const double sqrt2thrd=std::sqrt(2.0/3);
static const double sqrt6th=std::sqrt(1.0/6);

Warning<InvalidParameter> NMR_asymmetry_warning("anisotropy_to_cartesian: asymmetry out of range",&lcm_serious_warning);

bool anisotropy_to_cartesian(double& xx,double& yy, double& zz,double iso,double aniso,double asym, ordering_convention_t convention)
{
  const bool ok=(asym<1+tol) && (asym>(-1.0-tol));
  if (!ok && NMR_asymmetry_warning.enabled()) {
    char buf[64];
    snprintf(buf,sizeof(buf)," (%g)",asym);
    NMR_asymmetry_warning.raise(buf);
  }
  xx=iso-0.5*aniso*(1+asym);
  yy=iso-0.5*aniso*(1-asym);
  zz=iso+aniso;
  switch (convention) {
  case convention_Haeberlen:
    break;
  case convention_NQR:
    std::swap(xx,yy);
    break;
  default:
    throw InvalidParameter("anisotropy_to_cartesian: only 'Haerberlen' and 'NQR' conventions supported");
  }
  //  if (!ok && verbose_NMR)
  //  ::std::cerr << "Warning: anisotropy_to_cartesian: asymmetry out of range: " << asym << '\n';
  return ok;
}

Warning<InvalidParameter> NMR_axis_ordering_warning("cartesian_to_anisotropy: ordering of x, y, z inconsistent with convention used",&lcm_serious_warning);

bool cartesian_to_anisotropy(double& iso, double& aniso, double& asym, double xx, double yy, double zz, ordering_convention_t convention)
{
  iso=(xx+yy+zz)/3.0;
  xx-=iso;
  yy-=iso;
  zz-=iso;
  //std::cout << "x,y,z: " << xx << ',' << yy << ',' << zz << '\n';
  const double fxx=std::fabs(xx);
  bool ok=(fxx<=std::fabs(zz));
  if (ok) {
    switch (convention) {
    case convention_Haeberlen:      
      ok=(std::fabs(yy)<=fxx);
      break;
    case convention_NQR:
      ok=(fxx<=std::fabs(yy));
      break;
    default:
      throw InvalidParameter("cartesian_to_anisotropy: only 'Haerberlen' and 'NQR' conventions supported");
    }
  }
  if (!ok)
    NMR_axis_ordering_warning.raise();
  aniso=zz;
  asym=std::fabs((xx-yy)/zz);
  return ok;
}

Warning<> tensor_ignore_warning("default contents passed to Tensor<T> are being ignored",&lcm_base_warning);

space_T spatial_tensor(double iso11,double d22,double asym33, ordering_convention_t convention)
{
  if (isorderingIUPAC(convention))
    return A2(iso11,d22,asym33);
  double xx,yy,zz;
  anisotropy_to_cartesian(xx,yy,zz,-sqrt3rd*iso11,d22*sqrt2thrd,asym33,convention);
  return A2(xx,yy,zz);
}

space_T spatial_tensor(double d,double asym, ordering_convention_t convention)
{
  if (isorderingIUPAC(convention))
    throw Failed("IUPAC ordering not appropriate for traceless tensors");
  double xx,yy,zz;
  anisotropy_to_cartesian(xx,yy,zz,0.0,d*sqrt2thrd,asym,convention);
  return A2(xx,yy,zz,2);
}

List<double> diag_spin_truncated_dipolar(const basespin_system &a,size_t i,size_t j)
  {
    List<double> Hdip;
    a.mla_Iz(Hdip,2.0,i,j);
    return Hdip;
  }

  List<double> diag_spin_dipolar(const basespin_system &a,size_t i,size_t j)
  {
    if (a(i)==a(j))
      throw Failed("diag_spin_dipolar: invalid for homonuclear spin pair");
    List<double> Hdip;
    a.mla_Iz(Hdip,2.0,i,j);
    return Hdip;
  }

cmatrix spin_dipolar(const basespin_system& sys, size_t i, size_t j)
{
  cmatrix H(mxflag::temporary);
  sys.mla_I(H,2.0,i,'z',j,'z');
  if (sys(i)==sys(j)) {
    sys.mla_I(H,-0.5,i,'+',j,'-');
    sys.mla_I(H,-0.5,i,'-',j,'+');
  }
  return H;
}

  List<double> diag_spin_weakJ(const basespin_system& a,size_t i,size_t j)
  {
    return diag_Iz(a,i,j);
  }

List<double> diag_spin_quadrupolar(const basespin_system& a, size_t i)
{
  List<double> res(diag_spin_quadrupolar_Cq(a,i));
  res*=(1.0/3.0);
  return res;
}

List<double> diag_spin_quadrupolar_Cq(const basespin_system& a,size_t i, size_t order)
{
  const double Iv=a(i).qn();
  if (Iv<0.99999)
    throw Failed("diag_spin_quadrupolar: spin not quadrupolar");

  List<double> tmp(mxflag::temporary);
  a.mla_Iz(tmp,1.0,i);    
  const double I2=Iv*(Iv+1);

  switch (order) {
  case 0: {
    //    const double tr_val=I2/3.0;
    for (size_t j=tmp.size();j--;) {
      double& el=tmp(j);
      el=3.0*el*el-I2;
    }
  }
    break;
  case 1: {
    const double tr_val=4*I2-1;
    for (size_t j=tmp.size();j--;) {
      double& el=tmp(j);
      el*=tr_val-8.0*el*el;
    }    
  }
    break;
  case 2: {
    const double tr_val=2*I2-1;
    for (size_t j=tmp.size();j--;) {
      double& el=tmp(j);
      el*=tr_val-2.0*el*el;
    }    
  }
    break;
  default:
    throw InvalidParameter("diag_spin_quadupolar::order");
  }
  return tmp;
} 

cmatrix spin_quadrupolar(const basespin_system& a,size_t i)
{
  cmatrix tmp(mxflag::temporary);
  full(tmp,diag_spin_quadrupolar(a,i));
  return tmp;
}

cmatrix spin_CS(const basespin_system &sys,size_t i)
{
  cmatrix total(mxflag::temporary);
  sys.mla_I(total,1.0,i,'z');
  return total;
}

void rotatez_ip(cmatrix& U,const BaseList<double>& Fz,double angle)
{
  if (!issquare(U))
    throw NotSquare("rotatez_ip");
  const size_t n=U.rows();
  if (n!=Fz.size()) 
    throw Mismatch("rotatez_ip");

  ScratchList<complex> zrot(n);
  rotatezfacs(zrot,Fz,angle);

  U.unitary_simtrans(zrot);
}

void rotatezfacs(BaseList<complex> zrot,const BaseList<double>& Fz,double angle)
{
  const size_t n=Fz.size();
  if (zrot.size()!=n)
    throw Mismatch("rotatezfacs");
  if (angle==0)
    zrot=1.0;
  else
    for (size_t i=n;i--;)
      zrot(i)=expi(-angle*Fz(i));  
}

void rotatezfacs(cmatrix& d, const BaseList<double>& Fz, const BaseList<double>& angles)
{
  size_t n=angles.size();
  d.create(n,Fz.size());
  for (;n--;) {
    BaseList<complex> dr(d.row(n));
    rotatezfacs(dr,Fz,angles(n));
  }
}

cmatrix rotatez(const cmatrix& U, const BaseList<double>& Fz, double angle)
{
  cmatrix dU(U,mxflag::temporary);
  if (angle)
    rotatez_ip(dU,Fz,angle);
  return dU;
}

cmatrix spin_rf(const basespin_system& a,nuclei_spec whichn,phase_spec phase, double scale)
{
  const size_t nuc=whichn();
  cmatrix d(mxflag::temporary);
  double sphase,cphase;
  cmatrix_sincos(phase(),sphase,cphase);
  if (fabs(cphase)>tol)
    a.mla_F(d,cphase*scale,nuc,'x');
  if (fabs(sphase)>tol)
    a.mla_F(d,sphase*scale,nuc,'y');
  return d;
}

  
cmatrix Upulse(const basespin_system& a,nuclei_spec whichn,double angle,phase_spec phase)
{
  const size_t nuc=whichn();
  cmatrix U(mxflag::temporary);
  Upulse(U,real(F(a,nuc,'x')),diag_Fz(a,nuc),angle,phase);
  return U;
}

void Upulse(cmatrix& total,const rmatrix& Fx,const BaseList<double>& Fz,double angle,phase_spec phase)
{
  hermitian_expi(total,Fx,-angle);
  if (phase()) 
    rotatez_ip(total,Fz,phase());
}

cmatrix Upulse(const rmatrix& Fx,const BaseList<double>& Fz,double angle,phase_spec phase)
{
  cmatrix d(mxflag::temporary);
  Upulse(d,Fx,Fz,angle,phase); 
  return d;
}

// cmatrix Usoftpulse(const basespin_system &a,const cmatrix &H,nuclei_spec whichn,double time,double freq,double angle,phase_spec pspec)
// {
//   PulseGenerator pgen(a,whichn);
//   return pgen(H,time,freq,angle,pspec);
// }

// void Usoftpulse(cmatrix& U,const cmatrix& Fx,const BaseList<double>& Fz,const cmatrix& H,double time,double freq,double angle,phase_spec phase)
// {
//   PulseGenerator pgen(Fx,Fz);
//   pgen(U,H,time,freq,angle,phase);
// }

void propagator(cmatrix& d,const cmatrix& a,double t)
{
  if (a.rows()==1)
    d.create(1,1,propagator(real(a(0U,0U)),t));
  else
    hermitian_exp(d,a,complex(0,MINUS_TWO_PI*t));
}

void propagator(cmatrix& d,const rmatrix& a,double t)
{
  if (a.rows()==1)
    d.create(1,1,propagator(a(0U,0U),t));
  else
    hermitian_expi(d,a,MINUS_TWO_PI*t);
}

cmatrix propagator(const cmatrix& a,double t)
{
  cmatrix d(mxflag::temporary);
  propagator(d,a,t);
  return d;
}

cmatrix propagator(const rmatrix &a,double t)
{
  cmatrix d(mxflag::temporary);
  propagator(d,a,t);
  return d;
}

void propagator(List<complex>& d,const BaseList<double>& a,double t)
{
  d.create(a.size());
  propagator(static_cast< BaseList<complex> >(d),a,t);
}

List<complex> propagator(const BaseList<double> a,double t)
{
  List<complex> d(a.size(),mxflag::temporary);
  propagator(static_cast< BaseList<complex> >(d),a,t);
  return d;
}

void propagator(BaseList<complex> d,const BaseList<double>& a,double t)
{
  size_t n=a.size();
  if (d.size()!=n
      ) throw Mismatch("propagator");
  t*=MINUS_TWO_PI;
  for (;n--;)
    d(n)=expi(t*a(n));
}

//find eigensystem of (Hilbert space) propagator
void diag_propagator(cmatrix& D, List<double>& eff_freq, const cmatrix& U,double dt, double ptol)
{
  eff_freq.create(U.rows());
  diag_propagator(D,static_cast< BaseList<double> >(eff_freq),U,dt,ptol);
}

void diag_propagator(cmatrix& D,BaseList<double> eff_freq, const cmatrix& U,double dt, double ptol)
{
  ScratchList<complex> cycle_eigs(eff_freq.size());
  if (ptol)
    eigensystem(D,cycle_eigs,U,ptol);
  else
    eigensystem(D,cycle_eigs,U);
  diag_propagator(eff_freq,cycle_eigs,dt);
}

void diag_propagator(BaseList<double> eff_freq,const BaseList<complex>& U,double dt)
{
  if (dt==0.0)
    throw Failed("diag_propagator: zero time interval");
  const size_t dim=U.size();
  if (dim!=eff_freq.size()) 
    throw Mismatch("diag_propagator");
  const double scale=1.0/(dt*MINUS_TWO_PI);
  for (size_t j=dim;j--;)
    eff_freq(j)=scale*arg(U(j));
}

void diag_propagator(List<double>& eff_freq,const BaseList<complex>& U,double dt)
{
  eff_freq.create(U.size());
  diag_propagator( static_cast< BaseList<double> >(eff_freq),U,dt);
}

double diag_propagator(const complex& U,double dt)
{
  if (dt==0.0)
    throw Failed("diag_propagator: zero time interval");
  return arg(U)/(MINUS_TWO_PI*dt);
}

// The complex factors are stored along the diagonal to avoid a temporary!
static void _evolution_matrix(cmatrix &to)
{
  size_t i;
  for (i=to.rows();i--;) {
    const complex &wi=to(i,i);
    for (size_t j=0;j<i;j++)
      to(j,i)=conj(to(i,j)=conj_multiply(to(j,j),wi));
  }
  for (i=to.rows();i--;) to(i,i)=complex(1.0);
}

void evolution_matrix(cmatrix &to,const BaseList<double> &eigs,double t)
{
  const size_t n=eigs.size();
  to.create(n,n);

  if (t==0) {
    to=1.0;
    return;
  }
  t*=MINUS_TWO_PI;

  for (size_t i=n;i--;) to(i,i)=expi(t*eigs(i));
  _evolution_matrix(to);
}

void evolution_matrix(cmatrix &to,const BaseList<complex> &w)
{
  const size_t n=w.size();
  to.create(n,n);

  for (size_t i=n;i--;) {
    to(i,i)=1.0;
    for (size_t j=0;j<i;j++)
      to(j,i)=conj(to(i,j)=conj_multiply(w(j),w(i)));
  }
}

void evolution_matrix(cmatrix &to,const BaseList<complex> &wL,const BaseList<complex> &wR)
{
  const size_t nr=wL.size();
  const size_t nc=wR.size();
  to.create(nr,nc);

  for (size_t r=nr;r--;) {
    const complex &rfac=wL(r);
    BaseList<complex> tor=to.row(r);
    for (size_t c=nc;c--;)
      tor(c)=conj_multiply(wR(c),rfac);
  }
}

void evolution_matrix(BaseList<complex> &to,const complex &wL,const BaseList<complex> &wR)
{
  size_t c=wR.size();
  if (to.size()!=c)
    throw Mismatch("evolution_matrix");
  for (;c--;)
    to(c)=conj_multiply(wL,wR(c));
}

void evolution_matrix(BaseList<complex> &to,const BaseList<complex> &wL,const complex &wR)
{
  size_t r=wL.size();
  if (to.size()!=r)
    throw Mismatch("evolution_matrix");
  for (;r--;)
    to(r)=conj_multiply(wL(r),wR);
}

// void propagate(cmatrix& CU,const BinaryFunction<cmatrix,double,double>& propmeth,double start,double end,double intdt)
// {
//   if (intdt<=0.0) throw InvalidParameter("propagate: invalid time integration step");

//   const double diff=end-start; 
//   int steps=int(0.5+diff/intdt);
//   if (steps<1) steps=1;
  
//   const double dt=diff/steps;

//   cmatrix tmp,tmp2;

//   for (int n=0;n<steps;n++) {
//     const double t=start+n*dt;
// #ifndef NDEBUG
//     // cout << "Propagator from " << t*1e3 << " to " << (t+dt)*1e3 << " ms\n";
// #endif
//     if (isdefined(CU)) {
//       propmeth(tmp,t,t+dt);
// #ifndef NDEBUG
//       // cout << tmp << std::endl;
// #endif
//       multiply(tmp2,tmp,CU);
//       CU=tmp2;
//     }
//     else {
//       propmeth(CU,t,t+dt);
// #ifndef NDEBUG
//       // cout << CU << std::endl;
// #endif
//     }
//   }
// }

void propagators(BaseList<cmatrix> CUs,const PropGen_t& propmeth,double start,double end)
{
  cmatrix tmp;
  const size_t nprops=CUs.size();
  if (nprops==0)
    throw Failed("propagators: empty destination vector");
  
  const double tstep=(end-start)/nprops;
  double t=start;

  for (size_t j=0;j<nprops;j++) {
    const double nextt=(j==nprops-1) ? end : start+(j+1)*tstep;
    propmeth(j ? tmp : CUs.front(),t,nextt);
    if (j)
      multiply(CUs(j),tmp,CUs(j-1));
    t=nextt;
  }
}

void propagators(cmatrix& CUs,size_t nprops,const DiagPropGen_t& propmeth,double start,double end)
{
  if (nprops==0)
    throw Failed("propagators: no propagators to calculate");

  const double tstep=(end-start)/nprops;
  List<complex> Utmp;

  for (size_t j=0;j<nprops;j++) {
    propmeth(Utmp,start,start+(j+1)*tstep);
    if (!CUs)
      CUs.create(nprops,Utmp.size());    
    CUs.row(j)=Utmp;
  }
}

void coherencematrix_ip(Matrix<bool>& d,const BaseList<double>& zops, const BaseList<int>& cohers,int maxcoher)
{
  const size_t n=zops.size();
  if (!d)
    d.create(n,n,true);
  if (!maxcoher)
    maxcoher=coherencematrix_range(zops,cohers);

  ScratchList<bool> maskin(2*maxcoher+1,false);
  for (size_t j=cohers.size();j--;)
    maskin(cohers(j)+maxcoher)=true;
  
  for (size_t r=n;r--;) {
    const double offr=0.5+zops(r)+maxcoher;
    for (size_t c=n;c--;) {
      const bool isin=maskin(size_t(offr-zops(c)));      
      d(r,c)&=isin;
    }
  }
}

  int coherencematrix_range(const BaseList<double>& zops, const BaseList<int>& cohers)
  {
    int maxcoher=int(0.5+max(zops)-min(zops));
    for (size_t i=cohers.size();i--;) {
      const int coher(cohers(i));       
      if (coher>maxcoher || coher<-maxcoher)
	throw InvalidParameter("coherencematrix: invalid coherence");
    }
    return maxcoher;
  }
 
Matrix<bool> coherencematrix(const basespin_system& sys, const BaseList<int>& cohers)
{
  Matrix<bool> d(mxflag::temporary);
  coherencematrix_ip(d,diag_Fz(sys),cohers);
  return d;
}

Matrix<bool> coherencematrix(const basespin_system& sys, const BaseList<size_t>& which, const ListList<int>& cohers)
{
  Matrix<bool> d(mxflag::temporary);
  size_t n=which.size();
  if (n!=cohers.size())
    throw Mismatch("coherencematrix");
  for (;n--;)
    coherencematrix_ip(d,diag_Fz(sys,which(n)),cohers(n));
  return d;
}

Matrix<bool> coherencematrix(const basespin_system& sys, const BaseList<size_t>& which, const BaseList<int>& cohers)
{
  Matrix<bool> d(mxflag::temporary);
  size_t n=which.size();
  if (n!=cohers.size())
    throw Mismatch("coherencematrix");
  for (;n--;)
    coherencematrix_ip(d,diag_Fz(sys,which(n)),BaseList<int>(1,const_cast<int*>(&(cohers(n)))));
  return d;
}

Matrix<bool> coherencematrix(const basespin_system& sys, nuclei_spec which, const BaseList<int>& cohers)
{
  Matrix<bool> d(mxflag::temporary);
  coherencematrix_ip(d,diag_Fz(sys,which),cohers);
  return d;
}

int getL(size_t len)
{
  switch (len) {
  case 5: return 2;
  case 3: return 1;
  case 1: return 0;
  }
  throw InvalidParameter("Hamiltonian doesn't appear to be rank 2");
}

const double sqrtfac=-std::sqrt(0.5);

int validateL(const space_T& A,size_t len)
{
  const int useL=getL(len);
  int r=A.max_rank();
  if (r==0)
    return 0;
  for (;r>0;r--) {
    if ((!A.have_rank(r)) ^ (r!=useL))
      throw InvalidParameter("add: interaction has wrong rank");
  }
  return useL;
}

Warning<> propagation_emptyhamiltonian_warning("propagator object: Hamiltonian was empty!",&lcm_base_warning);

// void add_Hamiltonian_ns(cmatrix& Htot,const space_T& A, const BaseList<cmatrix>& H)
// {
//   const int useL=validateL(A,H.size());

//   for (int m=-useL;m<=useL;m++) {
//     complex coeff=A(2,-m);
//     if (m)
//       mla(Htot,A(2,-m),H(m+useL));
//     else {
//       double rcoeff=real(coeff);
//       if (A.have_rank(0))
// 	rcoeff+=sqrtfac*real(A(0,0));
//       mla(Htot,rcoeff,H(useL));
//     }
//   }
// }

}//namespace libcmatrix
