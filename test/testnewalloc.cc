#include "cmatrix.h"
#include <iostream>

struct dodgy;
std::ostream& operator<< (std::ostream&, const dodgy&);

static int initcount=0;
static int destroycount=0;
static int copycount=0;
//static int assigncount=0;

struct dodgy {
  int val;

  explicit dodgy(int valv =1234) : val(valv) { 
    std::cerr << "Make: " << valv << '\n';
    initcount++;
  }
  dodgy(const dodgy& a) { val=a.val;
    std::cerr << "Copy: " << *this << '\n';
    copycount++;
  }
//   dodgy& operator= (const dodgy& a) {
//     std::cerr << "Assigned: " << *this << '\n';
//     assigncount++;
//     val=a.val;
//     return *this;
//   }
  ~dodgy() {
    //if (val!=1234)
      std::cerr << "Uninitialised dodgy: " << *this << '\n';
    destroycount++;
  }

  
};

std::ostream& operator<< (std::ostream& ostr, const dodgy& a)
{
  return ostr << a.val;
}

using namespace libcmatrix;

template<class T> bool stress(T& obj, int left)
{
  typedef typename T::value_type Type;
  char junk[]="DEADBEEFDEADBEEFDE"; //mess up stack a bit
  try {
    T newlist(3,mxflag::normal);
    std::cout << newlist << '\n';
    //newlist=dodgy(5);
    //std::cout << newlist << '\n';
    newlist.create(4,Type(99));
    newlist.push_back(Type(100));
    std::cout << newlist << '\n';

    //T newlist2(5,Type(3),mxflag::normal);
    //std::cout << newlist2 << '\n';

  } catch (MatrixException& exc) {
    std::cerr << exc << '\n';
    return false;
  }
  if (--left>0)
    stress(obj,--left);
  return true;
}

int main()
{
  List<dodgy> dlist1;
  stress(dlist1,3);

  std::cout << "Created: " << initcount << '\n';
  std::cout << "Copied: " << copycount << '\n';
  std::cout << "Destroyed: " << destroycount << '\n';
  //  std::cout << "Assigned: " << copycount << '\n';

  cmatrix b(4,4);
  Matrix<int> x(5,5); //this won't show up

  for (size_t j=4;j--;) {
    cmatrix a(3,3,4.0);
    std::cout << a << '\n';
    matrix_traits<complex>::allocator::print(std::cout);
  }

  return 0;
}
      
