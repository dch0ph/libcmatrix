#include "cmatrix.h"
#include "rmatrix.h"
#include "basespin_system.h"

namespace libcmatrix {
  cmatrix singletransop(int,size_t,size_t,char);
  rmatrix singletransop(int,size_t);
  cmatrix Isingle(const basespin_system&, size_t,size_t,size_t,char);
  rmatrix Isingle(const basespin_system&, size_t,size_t);
  cmatrix Fsingle(const basespin_system&, nuclei_spec, size_t,size_t,char);
  rmatrix Fsingle(const basespin_system&, nuclei_spec, size_t);
  inline cmatrix Fsingle(const basespin_system& sys,size_t r, size_t s, char c) { return Fsingle(sys,NULL_NUCLEUS,r,s,c); }
  inline rmatrix Fsingle(const basespin_system& sys, size_t r) { return Fsingle(sys,NULL_NUCLEUS,r); }
}
