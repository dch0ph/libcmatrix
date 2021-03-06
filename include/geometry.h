#ifndef _geometry_h_
#define _geometry_h_

#include <cmath>
#include <iostream>
#include "Euler.h"
#include "Matrix.h"

namespace libcmatrix {

  struct vector3;

struct spherical {
  double r;
  double theta;
  double phi;

  spherical() {}
  spherical(double _r,double _theta,double _phi =0) : r(_r), theta(_theta), phi(_phi)  {}
  explicit spherical(const vector3&);
};

std::ostream& operator << (std::ostream &,const spherical &);

typedef double rmatrix3 [3][3];
 void rmatrix3_to_rmatrix(Matrix<double>&, const rmatrix3&);

struct vector3 {
  double x,y,z;

  vector3() {}
  vector3(double _x,double _y,double _z =0) : x(_x), y(_y), z(_z) {}
  explicit vector3(const spherical&);

  double norm() const { return x*x+y*y+z*z; }
  double length() const { return std::sqrt(x*x+y*y+z*z); }

  void operator+= (const vector3 &a) { x+=a.x; y+=a.y; z+=a.z; }
  void operator-= (const vector3 &a) { x-=a.x; y-=a.y; z-=a.z; }
  void operator*= (double a) { x*=a; y*=a; z*=a; }
  void operator/= (double a) { x/=a; y/=a; z/=a; }

  void rotate(const Euler&);
  void rotate(const rmatrix3&);
};

inline vector3 operator- (const vector3 &a) { return vector3(-a.x,-a.y,-a.z); }
inline vector3 operator+ (const vector3 &a, const vector3 &b) { return vector3(a.x+b.x,a.y+b.y,a.z+b.z); }
inline vector3 operator- (const vector3 &a, const vector3 &b) { return vector3(a.x-b.x,a.y-b.y,a.z-b.z); }

inline vector3 operator* (const vector3 &a,double z) { return vector3(a.x*z,a.y*z,a.z*z); }
inline vector3 operator* (double z,const vector3 &a) { return vector3(a.x*z,a.y*z,a.z*z); }

inline vector3 operator/ (const vector3 &a,double z) { return vector3(a.x/z,a.y/z,a.z/z); }

double dot(const vector3&, const vector3&);
vector3 cross(const vector3&, const vector3&);
std::ostream& operator << (std::ostream&, const vector3&);

void rotation_matrix(rmatrix3 &,const Euler &);
void rotation_matrix(Matrix<double> &,const Euler &);
vector3 rotate(const vector3 &,const Euler &);
vector3 rotate(const vector3 &,const rmatrix3 &);

double vector_angle(const vector3&,const vector3&);

Euler vector_to_Euler(const spherical&); //!< return Euler angles associated with (symmetric) dipolar interaction oriented along vector (passive rotation of axis into PAS, active rotation of tensor originally along z)


} //namespace libcmatrix

#endif



