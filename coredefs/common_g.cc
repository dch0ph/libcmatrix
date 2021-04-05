/* components that are always compiled with symbols */

#include "basedefs.h"

namespace libcmatrix {

  void MatrixException::hook() const
  {
    (void)isspace(*what()); //!< harmless function to prevent elimination
  }

}
