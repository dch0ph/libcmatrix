#include "cmatrix.h"
#include "ttyio.h"
#include "Propagation.h"

using namespace std;
using namespace libcmatrix;

const double BADNORM=1e-2;

FILE *open_file(const char *base,const char *suf)
{
  char tmp[64];
  sprintf(tmp,"%s%s",base,suf);
  FILE *fp=fopen(tmp,"r");
  if (fp==NULL) {
    cerr << "Failed to open " << tmp << endl;
    exit(1);
  }
  return fp;
}

void read_file(List<complex> &FID,const char *file,double sw,int npoints)
{
  double amp,freq;

  FILE *fp=fopen(file,"r");
  //if (fpraw) {
   // read_vector(FID,fpraw);
   // fclose(fpraw);
   // if (FID.length()!=npoints) throw Mismatch("read_file");
   // return;
 // }

  FID.create(npoints);

  //FILE *fpamp=open_file(file,"Amp");
  //FILE *fpfreq=open_file(file,"Freq");

  FID=complex(0.0);
  int count=0;
  while (fscanf(fp,"%lg %lg",&amp,&freq)==2) {
    add_FID(FID,amp,freq/sw);
    count++;
  }
  cout << " (" << count << " items  sum " << real(FID(0)) << ")";
  fclose(fp);

  char tmp[128];
  sprintf(tmp,"%sFID",file);
  write_vector(tmp,FID);
}

int main(int argc,const char *argv[])
{
  int count=1;

  const int npoints=getint(argc,argv,count,"Number of FID points? ",256);
  const double sw=getint(argc,argv,count,"Spectral width? ",5000);

  const int nfiles=argc-count;
  if (nfiles<1) {
    cerr << "Need more than 1 file!\n";
    return 1;
  }

  List<complex> FIDbase;
  List<complex> FIDtest;
  double normref;
  int failcount=0;

  for (int i=0;i<nfiles;i++) {
    cout << (i==0 ? "Reading " : "Comparing ") << argv[i+count];
    cout.flush();
    if (i==0) {
      read_file(FIDbase,argv[i+count],sw,npoints);
      normref=sqrt(norm(FIDbase));
      cout << " (" << normref << ")\n";
      if (normref==0.0) {
	cerr << "No intensity in spectrum (spectral width too small?)\n";
	return 1;
      }
    }
    else {
      read_file(FIDtest,argv[i+count],sw,npoints);
      FIDtest-=FIDbase;
      const double checkval=sqrt(norm(FIDtest))/normref;
      if (checkval>BADNORM)
      	failcount++;
      cout << ": " << checkval << endl;
    }
  }
  cout << "Failures: " << failcount << '\n';
  return (failcount ? 1 : 0);
}
