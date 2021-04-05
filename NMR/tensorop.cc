/* Return specific tensor operators */

#include "tensorop.h"
#include "basespin_system.h"
#include "space_T.h"

static const double sqrthalf=sqrt(0.5);
static const double sqrtthird=sqrt(1.0/3);
static const double sqrtsixth=sqrt(1.0/6);

namespace libcmatrix {

cmatrix T1(const basespin_system &a,size_t spini,int order)
{
  cmatrix tmp(mxflag::temporary);
  switch (order) {
  case -1:
    a.mla_I(tmp,sqrthalf,spini,'-');
    break;
  case 0:
    a.mla_I(tmp,1.0,spini,'z');
    break;
  case 1:
    a.mla_I(tmp,-sqrthalf,spini,'+');
    break;
  default:
    throw BadIndex("T1");
  }
  return tmp;
}


cmatrix T2(const basespin_system& a,size_t i,size_t j,int rank,int order)
{
  cmatrix result(mxflag::temporary);

  switch (rank) {
  case 0:
    if (order!=0)
      throw BadIndex("T2");
    a.mla_I(result,-sqrtthird,i,'z',j,'z');
    a.mla_I(result,-sqrtthird,i,'x',j,'x');
    a.mla_I(result,-sqrtthird,i,'y',j,'y');
    break;
  case 1:
    switch (order) {
    case -1:
      a.mla_I(result,-0.5,i,'-',j,'z');
      a.mla_I(result,+0.5,i,'z',j,'-');
      break;
    case 0:
      a.mla_I(result,-0.5*sqrthalf,i,'+',j,'-');
      a.mla_I(result,+0.5*sqrthalf,i,'-',j,'+');
      break;
    case 1:
      a.mla_I(result,-0.5,i,'+',j,'z');
      a.mla_I(result,+0.5,i,'z',j,'+'); 
      break;
    default:
      throw BadIndex("T2");
    }
    break;

  case 2:
    switch (order) {
    case -2:
      a.mla_I(result,0.5,i,'-',j,'-');
      break;
    case -1:
      a.mla_I(result,0.5,i,'-',j,'z');
      a.mla_I(result,0.5,i,'z',j,'-');
      break;
    case 0:
      a.mla_I(result,2*sqrtsixth,i,'z',j,'z');
      a.mla_I(result,-sqrtsixth,i,'x',j,'x');
      a.mla_I(result,-sqrtsixth,i,'y',j,'y');
      break;
    case 1:
      a.mla_I(result,-0.5,i,'+',j,'z');
      a.mla_I(result,-0.5,i,'z',j,'+');
      break;
    case 2:
      a.mla_I(result,0.5,i,'+',j,'+');
      break;
    default:
      throw BadIndex("T2");
    }
    break;
  default:    
    throw BadRank("T2");
  }
  return result;
}

}//namespace libcmatrix
