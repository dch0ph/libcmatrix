#include "config.h" 
#include <cmath>
// needed to ensure M_PI defined
#include "geometry.h"
// needed for cmatrix_sincos
#include "cmatrix_complex.h"

namespace libcmatrix {

static const double deg_to_rad=M_PI/180.0;

std::ostream& operator << (std::ostream &ostr,const vector3 &a)
{
  return ostr << '(' << a.x << ", " << a.y << ", " << a.z << ')';
}

std::ostream& operator << (std::ostream &ostr,const spherical &a)
{
  return ostr << '(' << a.r << ", " << a.theta/deg_to_rad << ", " << a.phi/deg_to_rad << ')';
}

double dot(const vector3 &A,const vector3 &B)
{
  return A.x*B.x+A.y*B.y+A.z*B.z;
}

double vector_angle(const vector3 &A,const vector3 &B)
{
  const double costheta = dot(A,B)/(A.length()*B.length());
  errno=0;
  const double retval = acos(costheta);
  if (errno==0)
    return retval;
  if (costheta >= 1.0)
    return 0.0;
  if (costheta <= -1.0)
    return M_PI;
  throw Failed("vector_angle: unknown error returned by acos");
}

vector3 cross(const vector3 &A,const vector3 &B)
{
  return vector3(A.y*B.z -A.z*B.y, -A.x*B.z +A.z*B.x, A.x*B.y -A.y*B.x);
}

// Klaus SR p. 447 - designated "pseudo-active" : used as active rotation A' = R * A, but Euler angles defined in Rose convention (passive rotations of axis system)
void rotation_matrix(rmatrix3& R,const Euler& ang)
{
  double ca,sa,cb,sb,cg,sg;

  cmatrix_sincos(ang.alpha,sa,ca);
  cmatrix_sincos(ang.beta,sb,cb);
  cmatrix_sincos(ang.gamma,sg,cg);

  R[0][0]=ca*cb*cg-sa*sg;
  R[0][1]=sa*cb*cg+ca*sg;
  R[0][2]=-sb*cg;

  R[1][0]=-ca*cb*sg-sa*cg;
  R[1][1]=-sa*cb*sg+ca*cg;
  R[1][2]=sb*sg;

  R[2][0]=ca*sb;
  R[2][1]=sa*sb;
  R[2][2]=cb;
}

void rotation_matrix(Matrix<double>& R,const Euler& ang)
{
  double ca,sa,cb,sb,cg,sg;
  R.create(3,3);

  cmatrix_sincos(ang.alpha,sa,ca);
  cmatrix_sincos(ang.beta,sb,cb);
  cmatrix_sincos(ang.gamma,sg,cg);

  R(0,0)=ca*cb*cg-sa*sg;
  R(0,1)=sa*cb*cg+ca*sg;
  R(0,2)=-sb*cg;

  R(1,0)=-ca*cb*sg-sa*cg;
  R(1,1)=-sa*cb*sg+ca*cg;
  R(1,2)=sb*sg;

  R(2,0)=ca*sb;
  R(2,1)=sa*sb;
  R(2,2)=cb;
}

void vector3::rotate(const rmatrix3& R)
{
  const double ox=x;
  const double oy=y;

  x=R[0][0]*ox+R[0][1]*oy+R[0][2]*z;
  y=R[1][0]*ox+R[1][1]*oy+R[1][2]*z;
  z=R[2][0]*ox+R[2][1]*oy+R[2][2]*z;
}

//!< pre-multiplication A' = R*A i.e. assuming R is an active rotation
vector3 rotate(const vector3& A,const rmatrix3& R)
{
  return vector3 (R[0][0]*A.x+R[0][1]*A.y+R[0][2]*A.z,
		 R[1][0]*A.x+R[1][1]*A.y+R[1][2]*A.z,
		 R[2][0]*A.x+R[2][1]*A.y+R[2][2]*A.z);
}

void vector3::rotate(const Euler& rot)
{
  rmatrix3 R;
  rotation_matrix(R,rot);
  rotate(R);
}

vector3 rotate(const vector3& A,const Euler& rot)
{
  rmatrix3 R;
  rotation_matrix(R,rot);
  return rotate(A,R);
}

vector3::vector3(const spherical& a)
{
  double cphi,sphi,ctheta,stheta;
  cmatrix_sincos(a.theta,stheta,ctheta);
  cmatrix_sincos(a.phi,sphi,cphi);
  const double projr=a.r*stheta;
  x=projr*cphi;
  y=projr*sphi;
  z=a.r*ctheta;
}

spherical::spherical(const vector3 &a)
{
  r=std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
  if (r) {
    theta=acos(a.z/r);
    phi=atan2(a.y,a.x);
  }
  else
    theta=phi=0.0;
}

inline const BaseList<double> aslist(const rmatrix3& a) { return BaseList<double>(9,const_cast<double*>(a[0])); }
inline BaseList<double> aslist(rmatrix3& a) { return BaseList<double>(9,a[0]); }

void rmatrix3_to_rmatrix(Matrix<double>& d, const rmatrix3& a)
{
  d.create(3,3);
  d.row()=aslist(a);
}

Euler Euler::operator- () const
{
  if (beta)
    return Euler(M_PI-gamma, beta, M_PI-alpha);
  else {
    static const double TWOPI=2*M_PI;

    double total=-(alpha+gamma);
    if (total>=TWOPI)
      total-=TWOPI;
    else {
      if (total<=-TWOPI)
	total+=TWOPI;
    }
    return Euler(total,0.0,0.0);
  }
}

Euler vector_to_Euler(const spherical& sphco)
{
  return Euler(0.0,sphco.theta,M_PI-sphco.phi);
}

Euler_controller::Euler_controller() :
  verbose(0),
  trig_tolerance(1e-8),
  gimbal_tolerance(1e-6), //!< this needs to be relatively high
  verify_tolerance(0.0), //!< off by default
  asymmetry_tolerance(1e-6), //!< tolerance for forcing asymmetry to be zero
  eigenvalue_flip_pattern(0) //!< eigenvector flip pattern (debug only)
{}

Euler_controller cmatrix_euler_controller; //!< default values object

double Euler_controller::filter_trigvalue(double v) const
{
  if (v>1.0) {
    if (v-1.0<trig_tolerance)
      return 1.0;
  }
  else {
    if (v<-1.0) {
      if (-v-1.0<trig_tolerance)
	return -1.0;
    }
    else
      return v; //!< value in range
  }
  throw Failed("Euler_controller: invalid direction cosine (<-1.0 or >1.0)");
}

//! ZYZ convention - R must do an active rotation and columns must form an orthonormal right-handed set, though Euler angles are defined in terms of a passive rotation of axis system

Euler DCM_to_Euler(const Matrix<double>& R, const Euler_controller& ctrl)
{
  if (!issquare(R) || (R.rows()!=3))
    throw InvalidParameter("DCM_to_Euler");
  if (ctrl.verify_tolerance<0.0)
    throw InvalidParameter("DCM_to_Euler: verify_tolerance cannot be <0");
 
  double cosb=ctrl.filter_trigvalue(R(size_t(2),size_t(2)));
  double beta=acos(cosb);  //!< acos is defined to return 0...pi
  const double sb=std::sin(beta);
  double alpha,gamma;

  if (sb<ctrl.gimbal_tolerance) {
    if (cosb<0) { //!< clean up values
      cosb=-1.0;
      beta=M_PI;
    }
    else {
      cosb=1.0;
      beta=0.0;
    }
    const double cosalpha=ctrl.filter_trigvalue(R(size_t(0),size_t(0))/cosb); //!< gimbal lock - only sensitive to alpha + gamma. Return gamma = 0 (works for beta=0 and pi)
    const double sinalpha=-R(size_t(1),size_t(0));
    alpha=atan2(sinalpha,cosalpha);    
    return Euler( beta ? -alpha : alpha, 0.0,0.0); //!< "normalise" to beta =0
  }

//   if (transpose) {
//     alpha=atan2(R(1U,2U),R(0U,2U));
//     gamma=atan2(R(2U,1U),-R(2U,0U));
//   }
//   else {
  alpha=atan2(R(2U,1U),R(2U,0U));
  gamma=atan2(R(1U,2U),-R(0U,2U));

  const Euler euler(alpha,beta,gamma);

  if (ctrl.verify_tolerance) {
    Matrix<double> newR;
    rotation_matrix(newR,euler);
    const double normdiff=norm(newR-R);
    if (normdiff>ctrl.verify_tolerance) {
      if (ctrl.verbose) {
	std::cerr << "Initial direction cosine matrix:\n" << R;
	std::cerr << "Euler angles: " << euler << '\n';
	std::cerr << "Recreated direction cosine matrix:\n" << newR << std::endl;
      }
      throw Failed("DCM_to_Euler: recreated direction cosine matrix differs from input by more than tolerance");
    }
  }
  
  return euler;
}

}//namespace libcmatrix
