/* Single transition operators */

#include <ctype.h>
#include <cstring>
#include "singletransop.h"

namespace libcmatrix {

cmatrix singletransop(int m,size_t r,size_t s,char c)
{
  if (r<1 || r>m || s<1 || s>m || s==r) throw InvalidParameter("singletransop");

  if (m==2) return spinhalfop(c);

  r--; s--;

  cmatrix d(m,m,0.0,mxflag::temporary);

  switch (tolower(c)) {
  case 'x':
    d(r,s)=0.5;
    d(s,r)=0.5;
    break;
  case 'y':
    d(r,s)=complex(0,-0.5);
    d(s,r)=complex(0, 0.5);
    break;
  case 'z':
    d(r,s)=0.5;
    d(s,r)=-0.5;
    break;
  case '+':
    d(r,s)=1.0;
    break;
  case '-':
    d(s,r)=1.0;
    break;
  default:
    throw InvalidParameter("Invalid operator");
  }
  return d;
}

rmatrix singletransop(int m,size_t r)
{
  if (r<1 || r>m) throw InvalidParameter("singletransop");

  rmatrix d(m,m,0.0,mxflag::temporary);
  d(r-1,r-1)=1.0;
  return d;
}

cmatrix Isingle(const basespin_system& Q,size_t index,size_t r,size_t s,char c)
{
  cmatrix d(mxflag::temporary);
  Q.expandop(d,index,singletransop(Q(index).deg(),r,s,c));
  return d;
}

rmatrix Isingle(const basespin_system& Q,size_t index,size_t r)
{
  rmatrix d(mxflag::temporary);
  Q.expandop(d,index,singletransop(Q(index).deg(),r));
  return d;
}

cmatrix Fsingle(const basespin_system& a,nuclei_spec nucspec,size_t r,size_t s,char op)
{ 
  cmatrix d(mxflag::temporary);
  cmatrix tmp;
  const size_t nuc=nucspec();
  for (size_t i=a.nspins();i--;) {
    if (nuc==NULL_NUCLEUS || a(i).nucleus()==nuc) {
      a.expandop(tmp,i,singletransop(a(i).deg(),r,s,op));
      d+=tmp;
    }
  }
  if (!d) return cmatrix(a.rows(),a.cols(),complex(0.0));
  return d;
}

rmatrix Fsingle(const basespin_system& a,nuclei_spec nucspec,size_t r)
{ 
  rmatrix d(mxflag::temporary);
  rmatrix tmp;
  const size_t nuc=nucspec();
  for (size_t i=a.nspins();i--;) {
    if (nuc==NULL_NUCLEUS || a(i).nucleus()==nuc) {
      a.expandop(tmp,i,singletransop(a(i).deg(),r));
      d+=tmp;
    }
  }
  if (!d) return rmatrix(a.rows(),a.cols(),0.0);
  return d;
}

}//namespace libcmatrix
