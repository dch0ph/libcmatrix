#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "rmatrix.h"
#include "powder.h"

using namespace std;
using namespace libcmatrix;

void writefile(const char *fname,PowderMethod &meth)
{
  FILE *fp=fopen(fname,"w");
  double total=0;
  double weight;
  Euler dest(M_PI/4,M_PI/4,M_PI/4);
  while (meth.next(dest,weight)) {
    fprintf(fp,"%g %g %g\n",dest.alpha,dest.beta,weight);
    total+=weight;
  }
  fclose(fp);
  printf("%s weight sum: %g\n",fname,total);
}

int main()
{
  PlanarGrid pg(30,40,sphere);
  writefile("PG3040",pg);

  PlanarGrid pgh(30,20,hemisphere);
  writefile("PGH3020",pgh);

  PlanarGrid pgo(20,25,octant);
  writefile("PGO2025",pgo);

  PlanarGrid pga(30,1);
  writefile("PGA",pga);

  PlanarGrid pgb(1,40);
  writefile("PGB",pgb);

  SphericalGrid sgo(20,25,octant);
  writefile("SGO2025",sgo);

  SphericalGrid sga(30,1);
  writefile("SGA",sga);

  SphericalGrid sgb(1,40);
  writefile("SGB",sgb);

  int n;

  for (n=8;n<15;n++) {
    char buf[64];
    
    PlanarZCW pzcw(n,octant);
    sprintf(buf,"PZCW%i",n);
    writefile(buf,pzcw);
  }

  for (n=8;n<15;n++) {
    char buf[64];
    
    SphericalZCW szcw(n,hemisphere);
    sprintf(buf,"SZCW%i",n);
    writefile(buf,szcw);
  }

  rmatrix weights;
  read_matrix(weights,"PZCW12");

  ExplicitSampling exps(weights);
  writefile("ES",exps);
  return 0;
}
  
