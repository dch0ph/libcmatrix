/* Declarations for functions in utils directory */
#ifndef _cmatrix_utils_h_
#define _cmatrix_utils_h_

#include "List.h"
#include "cmatrix.h"
#include "rmatrix.h"
#include <sys/types.h>

namespace libcmatrix {

bool ispowerof2(int);

const int FT_FORWARD=-1;
const int FT_BACKWARD=1;

void ft_create_table(BaseList<complex> &,int =FT_FORWARD);

void fft_ip(BaseList<complex> &,int =FT_FORWARD,double =1.0);
List<complex> fft(const BaseList<complex>&, int =FT_FORWARD,double =1.0);
List<complex> fft(const BaseList<double>&, int =FT_FORWARD,double =1.0);

void ft(BaseList<complex> &,const BaseList<complex> &,int =FT_FORWARD,double =1.0);
inline void ft(List<complex> &d,const BaseList<complex> &a,int dir =FT_FORWARD,double scale =1.0)
{ d.create(a.length()); ft( static_cast< BaseList<complex>& >(d),a,dir,scale); }

void real_ft(BaseList<double> &,const BaseList<complex> &,int =FT_FORWARD,double =1.0);
inline void real_ft(List<double> &d,const BaseList<complex> &a,int dir =FT_FORWARD,double scale =1.0)
{ d.create(a.length()); real_ft( static_cast< BaseList<double>& >(d),a,dir,scale); }

void ft(BaseList<complex> &,const BaseList<complex> &,const BaseList<complex> &,double =1.0);
void ft(List<complex> &,const BaseList<complex> &,const BaseList<complex> &,double =1.0);

List<complex> ft(const BaseList<complex> &a,int dir =FT_FORWARD,double scale =1.0);

void real_ft(BaseList<double> &,const BaseList<complex> &,const BaseList<complex> &,double =1.0);
inline void real_ft(List<double> &d,const BaseList<complex> &a,const BaseList<complex> &facs,double scale =1.0)
{ d.create(a.length()); real_ft( static_cast< BaseList<double>& >(d),a,facs,scale); }

void fft_ip(cmatrix &,int =FT_FORWARD,double =1.0);
cmatrix fft(const cmatrix &,int =FT_FORWARD,double =1.0);

void ft(cmatrix &,const cmatrix &,int =FT_FORWARD,double =1.0);
cmatrix ft(const cmatrix &a,int dir =FT_FORWARD,double scale =1.0);

cmatrix phasefft(const cmatrix &,int =FT_FORWARD,double =1.0,double =1.0);
void phasefft_ip(cmatrix &,int =FT_FORWARD,double =1.0,double =1.0);

void ampfft_ip(cmatrix &,cmatrix &,int =FT_FORWARD,double =1.0,double =1.0);
void ampfft_ip(cmatrix &,int =FT_FORWARD,double =1.0,double =1.0);
cmatrix ampfft(const cmatrix &,int =FT_FORWARD,double =1.0,double =1.0);

typedef int randomseed_t;

class RandomGenerator
{
 public:
  RandomGenerator();
  RandomGenerator(randomseed_t);
  double operator() ();
  static randomseed_t get_seed_random();
  static randomseed_t get_seed_default();
  void set_seed(randomseed_t);
 private:
  long ix1,ix2,ix3;
  double r[98];
  randomseed_t seedval;
};

template<class Gen =RandomGenerator> class GaussianGenerator
  {
    public:
    GaussianGenerator() : rgen_() { set_correlationtime(0.0); }
    GaussianGenerator(double cortimv, double chartimv) : rgen_() { set_correlationtime(cortimv, chartimv); }
    GaussianGenerator(randomseed_t seedv) : rgen_(seedv) { set_correlationtime(0.0); }
    GaussianGenerator(double cortimv, double chartimv, randomseed_t seedv) : rgen_(seedv) { set_correlationtime(cortimv,chartimv); }
    void set_seed(randomseed_t seedv) { rgen_.set_seed(seedv); reset(); }
    double operator() ();
    double operator() (double);
    void reset();
    void set_correlationtime(double, double =0.0);

    private:
    Gen rgen_;
    bool iset;
    double gset;
    double cortim_;
    double chartim_;
    double weightfac_;
    double lme2_;
    double prev_;

    void set_chartime(double);
  };

double random(double, bool =false); // must have arg to distinguish from std::random
void set_seed(); //!< set seed using current time
void set_seed(randomseed_t);
double gauss(double =1.0, bool =false);
//double nr_random(randomseed_t*, bool =false);

void write_hypercomplex(FILE *,const cmatrix &,const cmatrix &,const char * =NULL,unsigned =(mxflag::doublep | mxflag::binary));
int write_hypercomplex(const char *,const cmatrix &,const cmatrix &,const char * =NULL,unsigned =(mxflag::doublep | mxflag::binary));

void exponential_multiply_ip(cmatrix&, double,double, size_t =1);
cmatrix exponential_multiply(const cmatrix& , double rtc,double ctc, size_t =1);

void exponential_multiply_ip(rmatrix&, double,double, size_t =1);
rmatrix exponential_multiply(const rmatrix&, double rtc,double ctc, size_t =1);

void exponential_multiply_ip(BaseList<double>,double);
List<double> exponential_multiply(const BaseList<double>&, double);

void exponential_multiply_ip(BaseList<complex>, double);
List<complex> exponential_multiply(const BaseList<complex>&, double);

void gaussian_multiply_ip(cmatrix&, double, double, double, double, size_t =1);
void gaussian_multiply_ip(BaseList<complex>, double lbg, double dt);

inline double phase_correction(double f, double zero, double first, double pivot =0.5) { return zero+first*(f-pivot); }

void phase_correct_ip(BaseList<complex>, double zero, double first, double pivot =0.5);
void phase_correct_ip(cmatrix&, double zero, double first, double pivot =0.5);
void phase_correct_ip(cmatrix&, double zero2, double first2, double pivot2, double zero1, double first1, double pivot1 =0.0);
void ampphase_correct_ip(cmatrix&, cmatrix&, double zero2, double first2, double pivot2, double zero1, double first1, double pivot1 =0.0);

off_t file_length(FILE*);
off_t file_length(const char*);
FILE* file_open(const char* s, const char* mode); //!< like fpopen(\a s,\a mode) but "a" does 'expected thing'

bool isfilebinary(FILE *);

class file_position_lock {
 public:
  file_position_lock(FILE* fpv) 
    : fp(fpv), pos(std::ftell(fpv)) {}
  ~file_position_lock() {
    if (fp)
      std::fseek(fp,pos,SEEK_SET); //!< restore to locked position
  }
  void unlock() { fp=NULL; } //!< flag OK
 private:
  FILE* fp;
  long pos;
};


} //namespace libcmatrix

#endif
