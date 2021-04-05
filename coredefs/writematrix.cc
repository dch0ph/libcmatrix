#include "cmatrix.h"
#include "rmatrix.h"

namespace libcmatrix {

extern const char RMATRIX_MAGIC[];
extern const char RMATRIX_MAGICB[];

/* Allow sloppy use of comments, skip blank lines */
void write_comment(FILE *fp,const char *c,const char *intro)
{
  if (!c)
    return;

  do {
    for (;*c=='\n' && *c;c++);
    if (!*c)
      break;
    if (intro)
      fputs(intro,fp);
    if (*c=='#')
      c++;
    for (;*c!='\n' && *c;c++)
      putc(*c,fp);
    putc('\n',fp);
  } while (*c);
}


// Write vector or line from data set

void write_vector(FILE *fp,const BaseList<complex> &vec,int flags)
{
  int i;
  const int points=vec.length();
  const complex *data=vec.vector();

  if (flags & mxflag::binary) {
    int fail=0;
    if (flags & mxflag::doublep) {
#ifdef LCM_COMPLEX_CHEAT
      fail=fwrite( (void *)(data),sizeof(complex),points,fp);
#else
      double tmp[2];

      for (i=points;i--;data++) {
	tmp[0]=real(*data);
	tmp[1]=imag(*data);
	fail+=fwrite( static_cast<void *>(tmp),2*sizeof(double),1,fp);
      }
#endif
    }
    else {
      float tmp[2];

      for (i=points;i--;data++) {
	tmp[0]=real(*data); tmp[1]=imag(*data);
	fail+=fwrite( static_cast<void *>(tmp),2*sizeof(float),1,fp);
      }
    }
    if (fail!=points)
      throw Failed("Truncated write");
  }
  else {
    int precision=(flags & mxflag::doublep) ? LCM_DOUBLE_PRES : LCM_SINGLE_PRES;
    char sepchar,termchar;

    if (flags & mxflag::block) {
      sepchar=' '; termchar=' ';
    }
    else {
      sepchar=(flags & mxflag::rmat) ? ' ' : '\n';
      termchar='\n';
    }
    for (i=points;i--;data++)
	fprintf(fp,"%.*g%c%.*g%c",precision,real(*data),sepchar,precision,imag(*data),termchar);
  }
}

template<class T> void write_matrix_(const char *fname,const Matrix<T>& a,const char *comment,int flags)
{
  FileHandle fp(fname,flags & mxflag::binary ? "wb" : "w");
  write_matrix(fp(),a,comment,flags);
}


template<class T> inline void _write_vector(const char *fname,const BaseList<T> &vec,int flags)
{
  FileHandle fp(fname,flags & mxflag::binary ? "wb" : "w");
  write_vector(fp(),vec,flags);
}

void write_vector(const char *fname,const BaseList<double>& vec,int flags)
{
  _write_vector(fname,vec,flags);
}

void write_vector(const char *fname,const BaseList<complex> &vec,int flags)
{
  _write_vector(fname,vec,flags);
}

// Write cmatrix structure

void write_matrix(FILE *fp,const cmatrix &a,const char *comment,int flags)
{
  if (!(flags & mxflag::rmat))
    write_comment(fp,comment,"%");
  else {
    fprintf(fp,"%s RC\n",
	    flags & mxflag::binary ? RMATRIX_MAGICB : RMATRIX_MAGIC);
    write_comment(fp,comment,"#");
    fprintf(fp,"%" LCM_PRI_SIZE_T_MODIFIER "u %" LCM_PRI_SIZE_T_MODIFIER "u\n",a.rows(),2*a.cols());

// Note that RMAT format doesn't support double precision binary
    if (flags & mxflag::binary) flags&=~mxflag::doublep;
  }
  
  if ((flags & mxflag::binary) || ((flags & mxflag::norowsep) & !(flags & mxflag::block))) {
    write_vector(fp,a.row(),flags);
    return;
  }

  const size_t rows=a.rows();
  for (size_t i=0;i<rows;i++) {
    write_vector(fp,a.row(i),flags);
    if (i<rows-1)
      putc('\n',fp);
  }
  if (flags & mxflag::block)
    putc('\n',fp);
}

// Write vector or line from data set

void write_vector(FILE* fp,const BaseList<double>& vec,int flags)
{
  int i;
  const double *data=vec.vector();
  const size_t points=vec.length();

  if (flags & mxflag::binary) {
    int fail=0;

    if (flags & mxflag::doublep)
      fail=fwrite(data,sizeof(double),points,fp);
    else {
      float tmp;

      for (i=points;i--;) {
	tmp=*data++;
	fail+=fwrite(&tmp,sizeof(float),1,fp);
      }
    }
    if (fail!=points)
      throw Failed("Truncated write");
  }
  else {
    int precision=(flags & mxflag::doublep) ? LCM_DOUBLE_PRES : LCM_SINGLE_PRES;
    char termchar=(flags & mxflag::block) ? ' ' : '\n';

    for (i=points;i--;)
	fprintf(fp,"%.*g%c",precision,*data++,termchar);
  }
}

void write_matrix(FILE *fp,const rmatrix& a, const char *comment,int flags)
{
  if (!(flags & mxflag::rmat))
    write_comment(fp,comment,"%");
  else {
    fprintf(fp,"%s RR\n",
	    flags & mxflag::binary ? RMATRIX_MAGICB : RMATRIX_MAGIC);
    write_comment(fp,comment,"#");
    fprintf(fp,"%" LCM_PRI_SIZE_T_MODIFIER "u %" LCM_PRI_SIZE_T_MODIFIER "u\n",a.rows(),a.cols());

// Note that RMAT format doesn't support double precision binary
    if (flags & mxflag::binary)
      flags&=~mxflag::doublep;
  }
  
  if ((flags & mxflag::binary) || ((flags & mxflag::norowsep) & !(flags & mxflag::block))) {
    write_vector(fp,a.row(),flags);
    return;
  }

  const size_t rows=a.rows();
  for (size_t i=0;i<rows;i++) {
    write_vector(fp,a.row(i),flags);
    if (i<rows-1)
      putc('\n',fp);
  }
  if (flags & mxflag::block) 
    putc('\n',fp);
}

void write_matrix(const char* fname, const rmatrix& a, const char* comment, int flags)
{
  write_matrix_(fname,a,comment,flags);
}

void write_matrix(const char* fname, const cmatrix& a, const char* comment, int flags)
{
  write_matrix_(fname,a,comment,flags);
}

void write_hypercomplex(FILE *fp,const cmatrix &co,const cmatrix &si,const char *comment,unsigned flags)
{
  const size_t rows=co.rows();
  const size_t cols=co.cols();

  if (rows!=si.rows() || cols!=si.cols())
    throw Mismatch("write_hypercomplex");

  // Note that only the RMAT format is output
  fprintf(fp,"%s CC\n",
	  flags & mxflag::binary ? RMATRIX_MAGICB : RMATRIX_MAGIC);
  write_comment(fp,comment,"#");
  fprintf(fp,"%" LCM_PRI_SIZE_T_MODIFIER "u %" LCM_PRI_SIZE_T_MODIFIER "u\n",2*rows,2*cols);
    
  if (flags & mxflag::binary) flags&=~mxflag::doublep;

  int writes=!(flags & (mxflag::norowsep | mxflag::binary)); 

  for (size_t i=0;i<2*rows;i++) {
    write_vector(fp,(i % 2) ? si.row(i/2) : co.row(i/2),flags);
    if (writes && i<rows-1)
      putc('\n',fp);
  }
}

void write_hypercomplex(const char *fname,const cmatrix &co,const cmatrix &si,const char *comment,unsigned flags)
{
  FileHandle fp(fname,flags & mxflag::binary ? "wb" : "w");
  write_hypercomplex(fp(),co,si,comment,flags);
}

}//namespace libcmatrix
