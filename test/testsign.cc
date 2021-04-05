
#include <iostream.h>

const unsigned ui2=2;
const int i2=2;
const int im2=-2;

inline void testsign(unsigned n)
{
  cout << "greater than (unsigned int)2: " << ((n>ui2) ? "y\n" : "n\n");
  cout << "greater than 2: " << ((n>i2) ? "y\n" : "n\n");
  cout << "greater than -2: " << ((n>im2) ? "y\n" : "n\n");
}

int main()
{
  cout << "1:\n";
  testsign(1);
  cout << "\n";

  cout << "-1:\n";
  testsign(-1);
  cout << "\n";

  cout << "(unsigned int)1:\n";
  testsign( (unsigned int)3);
  cout << "\n";

  cout << "(unsigned int)-1:\n";
  testsign( (unsigned int)(-3) );
  cout << "\n";

  return 0;
}
