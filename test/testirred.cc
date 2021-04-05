/* Checking that we've got the definitions of spherical tensors correct! */

#include "NMR.h"
#include "powder.h"

using namespace std;
using namespace libcmatrix;

int main()
{
  try {
  const double sw=15000;

  const double isoshift=0;
  const double aniso=5000;
  const double asym=0.0;

  const space_T ACSA_PAS=spatial_tensor(isoshift,aniso,asym);
  std::cout << ACSA_PAS << endl;

  double xx,yy,zz;
  anisotropy_to_cartesian(xx,yy,zz,isoshift,aniso,asym);
  std::cout << xx << "  " << yy << "  " << zz << endl;
  std::cout << A2(xx,yy,zz) << endl;

  const space_T ADIP_PAS=spatial_tensor(aniso);

  const Euler ranang(M_PI/3,M_PI/5,M_PI/7);

  const int L=2;
  
  const cmatrix wignerels=D(L,ranang);

  cout << "Rank " << L << " rotation matrix\n" << wignerels << endl;

  space_T ACSA_R=rotate(ACSA_PAS,ranang);
  cout << "Normal rotate:\n" << ACSA_R << endl;
  cout << "Matrix rotate:\n" << rotate(ACSA_PAS,wignerels) << endl;

  List<cmatrix> DS(L+1);
  DS(L)=wignerels;
  cout << "All matrix rotate:\n" << rotate(ACSA_PAS,DS) << endl;

  cout << "orig rotate(2,1): " << ACSA_R(2,1) << endl;
  cout << "normal rotate(2,1): " << rotate(ACSA_PAS,2,1,ranang) << endl;
  cout << "matrix rotate(2,1): " << rotate(ACSA_PAS,1,wignerels) << endl;

  PlanarZCW twoangle(20);
  PlanarGrid betaonly(1,40);

  Euler powder;
  double weight;

  const int npoints=1000;
  
  List<double> spectrum(npoints);
  const double binres=sw/npoints;

  /* CS */
  spectrum=0.0;
  while (betaonly.next(powder,weight)) {
    double f=real(rotate(ACSA_PAS,2,0,powder)+ACSA_PAS(0,0));
    const cmatrix rotM(D(2,powder));
    const space_T ACSA(rotate(ACSA_PAS,rotM));
    double f2=sumzero(ACSA);
    cout << powder << " " << f << " " << f2 << '\n';
    if (f<0) f+=sw;
    const int bin=int(f/binres);
    spectrum(bin)+=weight;
  }

  write_vector("CSA",spectrum);
  }
  catch (MatrixException &exc) {
    cerr << exc << endl;
  }

  return 0;
}
