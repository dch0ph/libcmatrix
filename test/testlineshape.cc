// Direct lineshape addition

#include "Lineshapes.h"
#include "simpsonio.h"
#include "ScratchList.h"

using namespace std;
using namespace libcmatrix;

void makehist(const char* fname, const LineshapeSpec& paras, const HistogramSpec& histspec,double lw)
{
  char totname[256];
  const double sw=histspec.range();
  static const ExplicitList<4,double> freqs(0.0,-0.4,0.8,1.2);
  LorentzianGaussian shapegen(lw,0.0,paras);
  List<double> tmp;
  for (size_t shape=0;shape<3;shape++) {
    strcpy(totname,fname);
    double frac;
    switch (shape) {
    case 0:
      frac=LorentzianGaussian::PureLorentzian;
      strcat(totname,"_L");
      break;
    case 1:
      frac=LorentzianGaussian::PureGaussian;
      strcat(totname,"_G");
      break;
    case 2:
      frac=0.5;
      strcat(totname,"_M");
      break;
    }
    strcat(totname,".spe");
    shapegen.fraction(frac);
    List<complex> res(histspec.nbins,complex(0.0));
    LineshapeHistogram<complex> hist(res,shapegen,sw);
    for (size_t line=freqs.size();line--;) {
      const double freq=sw*freqs(line);
      hist.add(complex(line+1.0),freq);
    }
    std::cout << "Writing output file: " << totname << "\n\n";
    write_simpson(totname,res,sw);
  }
  shapegen.print();
  std::cout << std::endl;
}
  
int main()
{
  const size_t npts=512;
  const double sw=1e3; 
  const double lw=40.0;
  const HistogramSpec histspec(npts,sw);

  char fname[256];
  LineshapeSpec paras;
  int docache=1;
  for (size_t dovar=2;dovar--;) {
    for (size_t docut=2;docut--;) {
      paras.cutoff= docut ? 1e-1 : 0.0;
      for (size_t coarse=0;coarse<2;coarse++) {
	paras.resolution_steps=coarse ? 1 : 10;
	paras.cache_maximum=docache ? 10 : 0;
	for (size_t dofold=2;dofold--;) {
	  paras.flags=LineshapeSpec::verbose | (dofold ? 0 : LineshapeSpec::nofold) | (dovar ? LineshapeSpec::variable : 0);
	  sprintf(fname,"lineshape_%s_%i_%s_%s_%s",
		  (docut ? "cut" : "nocut"),
		  paras.resolution_steps,
		  (docache ? "cache" : "nocache"),
		  (dovar ? "variable" : "fixed"),
		  (dofold ? "fold" : "nofold"));
	  makehist(fname,paras,histspec,lw);
	}
      }
    }
  }
  return 0;
}
