#include "basedefs.h"
#include "List.h"
#include "MultiMatrix.h"
#include "Matrix.h"
#include "PermutedMatrix.h"
#include "PermutedMultiMatrix.h"

using namespace libcmatrix;
using namespace std;

 template<class T> std::ostream& operator<< (std::ostream& ostr, const PermutedMultiMatrix<T,2>& a) {
    Print_< PermutedMultiMatrix<T,2> ,2>::print(a,ostr);
    return ostr;
 }

int main()
{
  try {
    ListList<size_t> inds(ExplicitList<3,size_t>(2,3,2),ExplicitList<7,size_t>(0,1,0,2,1,1,0));
    const BaseList<size_t> sel1=inds(0);
    const BaseList<size_t> sel2=inds(1);
    const BaseList<size_t> sel3=inds(2);

    cout << "Index selection: " << inds << '\n';
    DirectSum_iterator<size_t> cur(inds);
    DirectSum_iterator<size_t> end(inds,true);
    while (cur!=end) cout << *cur++ << ' ';
    cout << '\n';

    MultiMatrix<int,1> Alist(2,3);
    cout << Alist << endl;
    MultiMatrix<int,1> Blist(3);
    MultiMatrix<int,2> Amatrix(2,3,4);
    cout << Amatrix << endl;
    const MultiMatrix<int,3> Bmatrix(2,3,4);
    cout << Bmatrix << endl;

    Matrix<int> Bsub(Bmatrix(0));
    cout << "Bsub\n" << Bsub << '\n';
    Bsub*=2;
    cout << "Bsub\n" << Bsub << '\n';

    MultiMatrix<float,1> mymat1(4,float(5.0));
    cout << mymat1 << endl;
    mymat1-=1.0;
    cout << "After subtract 1.0: " << mymat1 << endl;

    List<double> mylist(mymat1);
    cout << "Converted into List<double>: " << mylist << endl;

    MultiMatrix<float,2> A(4,3);
    for (int i=0;i<A.dimension(0);i++)
      for (int j=0;j<A.dimension(1);j++)
	  A(i,j)=i+10*j;
    cout << A << endl;

    MultiMatrix<double,2> B(4,3,3.0);
    B+=A;
    cout << "+B (3.0)\n" << B << endl;
    MultiMatrix<double,2> C=A+B;
    cout << "A+B\n" << C << endl;

    const MultiMatrix<float,2> cA(A);

    cout << "A(1,1): " << A(1,1) << endl;
    cout << "A(1,:)(1): " << A(1,range())(1) << endl;
    cout << "A(1): " << A(1) << endl;
    cout << "A():\n" << A() << endl;
    cout << "A(sel1,0): " << A(sel1,0) << endl;
    //cout << "A(sel1): " << A(sel1) << endl;
    cout << "(const) A(2,sel2): " << cA(2,sel2) << endl;
    //    cout << "A(0,0,sel3): " << A(0,0,sel3) << endl;
    cout << "A(sel1,sel2): " << A(sel1,sel2) << endl;
    cout << "A(sel1,:): " << A(sel1,range()) << endl;
    A(sel1,sel2)+=float(10.0);
    cout << "A(sel1,sel2)+=10.0\n" << A << endl;
    //cA(sel1,sel2)=float(0.0);
    
    MultiMatrix<bool,2> mymat2(3,4);
    mymat2=true;
    mymat2(1,2)=false;
    cout << mymat2 << endl;

    const PermutedMultiMatrix<bool,2> pmm(mymat2);
    //cout << "MultiMatrix dimensionality: " << type_traits<MultiMatrix<bool,2> >::dimensionality << '\n';
    //cout << "PermutedMultiMatrix dimensionality: " << type_traits< PermutedMultiMatrix<bool,2> >::dimensionality << '\n';
    cout << "Permuted version (multimatrix):\n" << pmm;
    cout << "(2,1): " << pmm(2,1) << "  (1,2): " << pmm(1,2) << "\n\n";

    PermutedMatrix<bool> pm(mymat2());
    cout << "Permuted version (matrix):\n" <<  pm;
    cout << "(2,1): " << pm(2,1) << "  (1,2): " << pm(1,2) << "\n\n";
    return 0;

  } catch (MatrixException& err) {
    cerr << err;
    return 1;
  }
}

