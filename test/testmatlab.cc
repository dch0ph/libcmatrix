#include "matlabio.h"
#include "cmatrix.h"
#include "rmatrix.h"
#include "MultiMatrix.h"
#include "BlockedMatrix.h"

using namespace std;
using namespace libcmatrix;

int main()
{
  try {

  rmatrix A(2,3,0.0);

  A(0,0)=1.0; A(0,1)=2.0;
  A(1,0)=3.0; A(1,1)=4.0;

  WriteMATLAB("Aorig",A,"A");
  WriteMATLAB("Aorig5",A,"A",5);
  WriteMATLAB("Along5",A,"Alongname",5);
  
  rmatrix tmpr;
  cmatrix tmpc;

  ReadMATLAB(tmpr,"Aorig");
  cout << "Read V4 into real\n" << tmpr << endl;

  ReadMATLAB(tmpc,"Aorig");
  cout << "Read V4 into complex\n" << tmpc << endl;

  tmpr.kill(); tmpc.kill();
  ReadMATLAB(tmpr,"Aorig5");
  cout << "Read V5 into real\n" << tmpr << endl;

  ReadMATLAB(tmpc,"Aorig5");
  cout << "Read V5 into complex\n" << tmpc << endl;

  cmatrix B=A*expi(M_PI/4);
  {
    matlab_controller ctrl("Borig",4);
    ctrl.write(B,"B");
    ctrl.write(B,"B2");
  }

  WriteMATLAB("B5new",B,"B",5);

  cout << "Original:\n" << B << endl;

  tmpc.kill(); ReadMATLAB(tmpc,"Borig");
  cout << "Read\n" << tmpc << endl;

  tmpc.kill(); ReadMATLAB(tmpc,"a5");
  cout << "Read A5PC\n" << tmpc << endl;

  tmpc.kill(); ReadMATLAB(tmpc,"a4pc");
  cout << "Read A4PC\n" << tmpc << endl;

  tmpc.kill(); ReadMATLAB(tmpc,"b5");
  cout << "Read B5PC\n" << tmpc << endl;

  tmpc.kill(); ReadMATLAB(tmpc,"b4pc");
  cout << "Read B4PC\n" << tmpc << endl;

  MultiMatrix<double,3> c(2,3,4);
  for (int i=c.dimension(0);i--;)
    for (int j=c.dimension(1);j--;)
      for (int k=c.dimension(2);k--;)
	c(i,j,k)=100*i+10*j+k;
  cout << "c(0): " << c(0) << "c(1): " << c(1) <<'\n';

  WriteMATLAB("C5",c);
  MultiMatrix<double,3> d;
  ReadMATLAB(d,"C5");
  cout << "d(0): " << d(0) << "d(1): " << d(1) <<'\n';

  BlockedMatrix<double> blkd(ScratchList<size_t>(2,3,4));
  blkd(0)=2.0;
  blkd(1)=3.0;
  blkd(2)=4.0;
  cout << blkd;
  
  {
    matlab_controller ctrl("fred",5,matlab_controller::singleprecision);
    ctrl.write(blkd,"density");
    {
      matlab_controller::composite comp(ctrl,"density2",matlab_controller::STRUCT);
      comp.write(blkd(0U),"one");
      comp.write<rmatrix>(blkd(1U),"two");
      comp.write<const rmatrix& >(blkd(2U),"three"); //!< this is the right way - avoids a copy
    }
    ctrl.write("a string","test");
    ctrl.write(ExplicitList<1,double>(6.5),"time");
  }
  {
    matlab_controller ctrl("fred");
    List< Matrix<double> > blke;
    ctrl.read(blke);
    BlockedMatrix<double> bme(blke);
    cout << bme;
    List<double> times;
    ctrl.read(times);
    cout << times.front() << '\n';
    List<char> str;
    ctrl.read(str);
    cout << str.vector() << '\n';
  }

  } catch (const MatrixException& exc) {
    cerr << exc;
    return 1;
  }
  return 0;
}


 
