#include "spin_system.h"
#include "spinhalf_system.h"
#include "ListList.h"
#include "NMR.h"

using namespace std;
using namespace libcmatrix;

int main()
{
  try {
  const int nspins=3;

  cout << "Test of spinhalf_system using " << nspins << " spins\n";

  spin_system a(nspins,"1H");

  const Permutation perm(ScratchList<size_t>(0,2,1));

  cout.precision(2);

  List<state_t> pvec1(a.rows()); a.permutation_vectorH(pvec1,perm);
  cout << "Permutation vector: " << pvec1 << endl;
  cout << "Permutation matrix (Hilbert)\n" << a.permutation_matrixH(perm);
  cout << "Permutation matrix (Liouville)\n"; spy(cout,a.permutation_matrixL(perm));

  cout << "Select rows/columns 0-3 for spinhalf_system\n";
  state_t sel[]={0,1,2,3};
  spinhalf_system ab(nspins,"1H",BaseList<state_t>(4,sel));
  cout << ab << endl;

  List<state_t> pvec2(ab.rows()); ab.permutation_vectorH(pvec2,perm);
  cout << "Permutation vector: " << pvec2 << endl;
  cout << "Permutation matrix (Hilbert)\n"; spy(cout,ab.permutation_matrixH(perm));
  cout << "Permutation matrix (Liouville)\n"; spy(cout,ab.permutation_matrixL(perm));

  List<double> Fz=diag_Fz(a);

  cout << "Each `block' operator should be the corresponding submatrix of the `normal' operator\n";

  if (ab.isdiagonal()) {
    cout << "Normal Fz:\n" << Fz << endl;
    cout << "Block Fz:\n" << diag_Fz(ab) << endl;
  }
  
    ListList<size_t> blist=find_blocks(Fz,0.0);
    cout << "\nBlock structure from find_blocks on Fz: " << blist << endl;  

  cout << "Normal Iz0:\n" << I(a,0,'z') << endl;
  cout << "Block Iz0:\n" << I(ab,0,'z') << endl;

  cout << "Normal Iy1:\n" << I(a,1,'y') << endl;
  cout << "Block Iy1:\n" << I(ab,1,'y') << endl;

  cout << "Normal Ix2:\n" << I(a,2,'x') << endl;
  cout << "Block Ix2:\n" << I(ab,2,'x') << endl;

  //  cout << "Normal Fx:\n" << HRF(a,1000.0,0) << endl;
  //cout << "Block Fx:\n" << HRF(ab,1000.0,0) << endl;

  cout << "Normal I1zI2z:\n" << I(a,1,'z',2,'z') << endl;
  cout << "Block I1zI2z:\n" << I(ab,1,'z',2,'z') << endl;

  if (ab.isdiagonal()) {
    cout << "Normal I1zI2z:\n" << diag_Iz(a,1,2) << endl;
    cout << "Block I1zI2z:\n" << diag_Iz(ab,1,2) << endl;
  }

  cout << "Normal I0xI1+:\n" << I(a,0,'x',1,'+') << endl;
  cout << "Block I0xI1+:\n" << I(ab,0,'x',1,'+') << endl;

  cout << "Normal I1zI2x:\n" << I(a,1,'z',2,'x') << endl;
  cout << "Block I1zI2x:\n" << I(ab,1,'z',2,'x') << endl;

  cout << "Normal I1+I2-:\n" << I(a,1,'+',2,'-') << endl;
  cout << "Block I1+I2-:\n" << I(ab,1,'+',2,'-') << endl;
  }
  catch (MatrixException &exc) {
    cerr << exc << endl;
  }

  return(0);
}



  
