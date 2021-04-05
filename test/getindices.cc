#include "ttyio.h"

using namespace libcmatrix;
using namespace std;

void getcheckindices(const char* str,size_t min,size_t max)
{
  cout << "Parsing '" << str << "': " << getuniqueindices(str,min,max) << '\n';
}

int main()
{
  const int MAXINDS=20;

  getcheckindices("- ",1,MAXINDS);
  getcheckindices("1-21",1,MAXINDS);
  getcheckindices("1-3,5,7-8",1,MAXINDS);
  getcheckindices("1,2,5,7-8 5",1,MAXINDS);
  getcheckindices("8-7",1,MAXINDS);
  getcheckindices(" 2, -1",1,MAXINDS);
  return 0;
}
