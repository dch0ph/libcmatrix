//#include "List.h"
#include <vector>
#include <memory>
#include <iostream>

//using namespace libcmatrix;

struct intpair {
  intpair(int x, int y) : x_(x), y_(y) { std::cout << "Constructed: "; this->print(); };   
  intpair(const intpair &) = delete;
  void operator= (const intpair&) = delete;
  void operator= (const intpair&& a) { x_ = a.x_; y_ = a.y_; }
  void print(std::ostream& ostr =std::cout) const
  { ostr << "x=" << x_ << ", y=" << y_ << '\n'; };
  int x_, y_;
};

class expression : public std::vector< std::unique_ptr<const intpair> > {
public:
  expression() {}
};

//template class std::vector< std::unique_ptr<const intpair> >;

std::ostream& operator<< (std::ostream& ostr, const intpair& a) { a.print(); return ostr; }

int main()
{
  intpair fred(2,3);
  fred.print();  
  //  std::vector< std::unique_ptr<intpair> > pairlist;
  expression pairlist;
  //  pairlist.push_back(new intpair(5,6));
  pairlist.emplace_back(new intpair(4,5));
  std::cout << *(pairlist[0]) << '\n';
  return 0;
}
