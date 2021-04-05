#include "basedefs.h"

using namespace libcmatrix;
using namespace std;

int main()
{
  try {
    throw BadIndex("test",3,2);
  }
//   catch (const MatrixException& exc) {
//     cerr << exc;
//     return 1;
//   }
  catch (const std::exception& exc) {
    cerr << exc.what() << '\n';
    return 1;
  }
  return 0;
}
