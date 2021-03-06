#ifndef simpsonio_h_
#define simpsonio_h_

#include <cstdio>
#include "List.h"
#include "Matrix.h"
#include "cmatrix_complex.h"
#include "Warnings.h"

namespace libcmatrix {

enum { FD_TYPE_FID=0, FD_TYPE_SPE=1 };
enum { FD_FORMAT_TEXT, FD_FORMAT_BINARY };
enum { FD_PREC_SINGLE, FD_PREC_DOUBLE };

struct simpsonFD {
  const complex* data;
  int np,ni,type,format,prec;
  double ref,ref1,sw,sw1;
  double sfrq,sfrq1;

  simpsonFD() : data(NULL), np(0), ni(0), type(FD_TYPE_FID), format(FD_FORMAT_TEXT), prec(FD_PREC_SINGLE), ref(0.0), ref1(0.0), sw(0.0), sw1(0.0), sfrq(0.0), sfrq1(0.0) {}

  explicit simpsonFD(const BaseList<complex>& a, double sw_ =0.0, bool isspectrum =true) : data(a.vector()), np(a.length()), ni(0), type(isspectrum ? FD_TYPE_SPE : FD_TYPE_FID), format(FD_FORMAT_TEXT), prec(FD_PREC_SINGLE), ref(0.0), ref1(0.0), sw(sw_), sw1(0.0), sfrq(0.0), sfrq1(0.0) {}

  explicit simpsonFD(const Matrix<complex>& a, double sw_ =0.0, double sw1_ =0.0, bool isspectrum =true) : data(a.vector()), np(a.cols()), ni(a.rows()), type(isspectrum ? FD_TYPE_SPE : FD_TYPE_FID), format(FD_FORMAT_TEXT), prec(FD_PREC_SINGLE), ref(0.0), ref1(0.0), sw(sw_), sw1(sw1_), sfrq(0.0), sfrq1(0.0) {}

  bool istimedomain() const { return (type==FD_TYPE_FID); }
};

 void read_simpson(List<complex>&, simpsonFD&, FILE*);
  void read_simpson(List<complex>&, const char*);
 void read_simpson(List<complex>&, simpsonFD&, const char*);
 void read_simpson(Matrix<complex>&, simpsonFD&, const char*);
 void read_simpson(Matrix<complex>&, const char*);
 void write_simpson(FILE*, const simpsonFD&);
 void write_simpson(const char*, const simpsonFD&);
 void write_simpson(const char*, const BaseList<complex>&, double sw, bool isspectrum =true);
 void write_simpson(const char*, const Matrix<complex>& , double sw, double sw1, bool isspectrum =true);
 void write_simpson_rows(const char*, const Matrix<complex>& , double sw, bool isspectrum =true);
 void write_simpson_rows(const char*, const Matrix<complex>& , double sw, bool isspectrum , double sfrq);
 void write_simpson_rows(const char*, const Matrix<complex>& , const BaseList<const char*>&, double sw, bool isspectrum =true);
 void write_simpson_rows(const char*, const Matrix<complex>& , const BaseList<const char*>&, double sw, bool isspectrum , double sfrq);
  List<char> insert_simpson(const char* fname,const char* ins, char sep ='\0');

  class simpson_controller {
  public:
    simpson_controller(const char*, size_t =100);
    simpson_controller(const char*, const BaseList<const char*>&);
    void write(size_t,const BaseList<complex>&, double, bool =true) const;
    void write(size_t,const simpsonFD&) const;
    static Warning<> ignoring_binary_flag_warning;

    void write_rows(const Matrix<complex>& a, double sw, bool isspectrum, double sfrq =0.0);
  private:
    int nsplit;
    std::string split;
    mutable List<char> tmp_;
    bool insert;
    int prec_;
    size_t max_;
    List<std::string> labels_;

    void findsplit(const char*, size_t);  
  };
}

#endif //!INCLUDED
