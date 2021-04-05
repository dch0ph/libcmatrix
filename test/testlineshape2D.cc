// Direct lineshape addition in two dimensions

#include "Lineshapes.h"
#include "matlabio.h"

using namespace std;
using namespace libcmatrix;

template<typename T> class BaseHistogram2D {
public:
  BaseHistogram2D(Matrix<T>& datav, double rsw, double csw)
    : data(datav), rowspec_(datav.rows(),rsw), colspec_(datav.cols(),csw) {
    create();
  }
  BaseHistogram2D(Matrix<T>& datav, double rsw, double rminv, double csw, double cminv)
    : data(datav), rowspec_(datav.rows(),rsw,rminv), colspec_(datav.cols(),csw,cminv) {
    create();
  }

  const HistogramSpec& rowspecification() const { return rowspec_; }
  const HistogramSpec& colspecification() const { return colspec_; }
  size_t lost() const { return lostv; }
  T lostsum() const { return lostsumv; }

  virtual ~BaseHistogram2D() {}
  virtual void add(const T&, double, double) =0;

protected:
  Matrix<T>& data;
  HistogramSpec rowspec_,colspec_;
  double rscalef,rdeltaf;
  double cscalef,cdeltaf;
  size_t lostv;
  T lostsumv;

private:
  BaseHistogram2D<T>(const BaseHistogram2D<T>&);
  BaseHistogram2D<T>& operator=(const BaseHistogram2D<T>&); //!< make private to prevent copies

  void create() {
    create_(rscalef,rdeltaf,rowspec_);
    create_(cscalef,cdeltaf,colspec_);
    lostv=0U;
    lostsumv=T(0);
  }

  void create_(double& scalef, double& deltaf, const HistogramSpec& spec) {
    scalef=spec.nbins/spec.rangef;
    deltaf=1.0/scalef;
  }
};

template<typename T> class LineshapeHistogram2D : public BaseHistogram2D<T>
{
  public:

  using BaseHistogram2D<T>::data;
  using BaseHistogram2D<T>::lostv;
  using BaseHistogram2D<T>::lostsumv;

  LineshapeHistogram2D(Matrix<T>& datav, BaseLineshape& rlshapev, double rmaxv, BaseLineshape& clshapev, double cmaxv)
   : BaseHistogram2D<T>(datav,rmaxv,cmaxv), rlineshape_(rlshapev), clineshape_(clshapev) {}

  LineshapeHistogram2D(Matrix<T>& datav, BaseLineshape& rlshapev, double rmaxv, double rminv, BaseLineshape& clshapev, double cmaxv, double cminv)
    : BaseHistogram2D<T>(datav,rmaxv,rminv,cmaxv,cminv), rlineshape_(rlshapev), clineshape_(clshapev) {}
    
  void add(const T& a, double rf, double cf);

private:
  BaseLineshape& rlineshape_;
  BaseLineshape& clineshape_;  //!< non-const ref not essential but prevents inappropriate copy and reduces risk of referencing a temporary object
  List<double> rtmp_,ctmp_;
};

template<typename T> void LineshapeHistogram2D::add(const T& a, double rf, double cf)
{
  const int roffset=rlineshape_(rtmp_,this->rowspecification(),rf);
  if (roffset>=0) {
    assert(roffset+rtmp_.size()<=data.rows());
    const int coffset=clineshape_(ctmp_,this->colspecification(),cf);
    if (coffset>=0) {
      assert(coffset+ctmp_.size()<=data.cols());
      const range colsel(coffset,coffset+ctmp_.size()-1);
      for (size_t r=rtmp_.size();r--;) {
	BaseList<T> dest(data(size_t(r+roffset),colsel));
	mla(dest,a*rtmp_(r),ctmp_);
      }
      return;
    }
  }
  lostv++;
  lostsumv+=a;
}


int main()
{
  const size_t cols=256;
  const size_t rows=48;
  const double csw=1e3;
  const double rsw=1e3;
  const double rlw=120.0;
  const double clw=180.0;
  const double rfrac=LorentzianGaussian::PureGaussian;
  const double cfrac=LorentzianGaussian::PureLorentzian;

  LineshapeSpec paras;
  paras.cutoff=0.0;
  paras.resolution_steps=5;
  //paras.flags=LineshapeSpec::verbose; //!< fold

  LorentzianGaussian rshapegen(rlw,rfrac,paras);
  LorentzianGaussian cshapegen(clw,cfrac,paras); //!< use same lineshape paras for both dimensions

  rmatrix data(rows,cols,0.0);
  LineshapeHistogram2D<double> hist(data,rshapegen,rsw,cshapegen,csw);
  hist.add(1.0,rsw*0.2,csw*0.2);
  hist.add(2.0,rsw*0.2,csw*0.6);
  hist.add(3.0,rsw*0.6,csw*0.6);
  hist.add(4.0,rsw*0.6,csw*0.2);
  
  matlab_controller matlabctrl("lineshape2D",4);
  matlabctrl.write(data,"lineshape");

  return 0;
}

  
