#include <cstring>
#include <cstdio>
#include "cmatrix.h"
#include "rmatrix.h"
#include <cstdlib>
#include <sys/stat.h>
#include "cmatrix_utils.h"

namespace libcmatrix {

extern const char RMATRIX_MAGIC[]="RMAT";
extern const char RMATRIX_MAGICB[]="RMATB";

static void read_comment(FILE*);

/* returns -1 in rows for block format matrices (i.e. cols determined)
   routines may refuse such matrices */

  namespace {

    inline bool iseol(char c) { return (c==EOF) || (c=='\n') || (c=='\r'); }

    void read_desc_file_(FILE *fp,int& rs,int& cs,bool& isbinary,bool& iscomplex,bool& integer)
    {
      char magic_str[8],extra_str[3];
      
      integer=false;
      iscomplex=false;
      
      const size_t curpos=ftell(fp);
      
      int items=fscanf(fp,"%7s %2s ",magic_str,extra_str);
      
      bool fail=false;

      if (items<1)
	fail=true;
      else {      
	if (strcmp(magic_str,RMATRIX_MAGIC)==0)
	  isbinary=false;
	else {
	  if (strcmp(magic_str,RMATRIX_MAGICB)==0)
	    isbinary=true; 
	  else {
	    if (strcmp(magic_str,"P2")==0) {
	      integer=true; isbinary=false;
	    }
	    else {
	      if (strcmp(magic_str,"P5")==0) {
		integer=true; isbinary=false;
	      }
	      else
		fail=true;
	    }
	  }
	}
      }

      if (!fail && !integer) {
	if (items<2)
	  throw FileCorrupt();
	
	switch (extra_str[1]) {
	case 'c':case 'C':
	  iscomplex=true;
	  break;
	case 'r':case 'R':
	  iscomplex=false;
	  break;
	default:
	  throw FileCorrupt();
	}
      }
      
      if (fail) {
	fseek(fp,curpos,SEEK_SET);
	
	/* Refuse to read raw binary files */
	if ((isbinary=isfilebinary(fp)))
	  throw Failed("Can't read raw binary file");
	
	read_comment(fp);
	
	size_t offset=ftell(fp);
	char c;
	
	for (cs=0;;) {
	  do {
	    c=getc(fp);
	  } while (c==' ' || c=='\t');
	  if (iseol(c))
	    break;
	  cs++;
	  do {
	    c=getc(fp);
	  } while (!isspace(c) && (c!=EOF));
	  if (iseol(c))
	    break;
	}
	fseek(fp,offset,SEEK_SET);
	
	if (cs<1)//!< fail if all white space
	  throw Failed("No data found");
	rs=-1;
	return;
      }
      
      /* Only if not from pipe */
      read_comment(fp);
      
      if (integer) {
	int max; 
	if (fscanf(fp,"%i %i %i \n",&cs,&rs,&max)<3)
	  throw FileCorrupt();
      }
      else {
	if (fscanf(fp,"%i %i \n",&rs,&cs)<2)
	  throw FileCorrupt();
      }
      
      if (iscomplex) {
	if (cs & 1)
	  throw FileCorrupt("Odd number of data entries impossible for complex data");
	cs/=2;
      }
      
      if (rs<1 || cs<1)
	throw FileCorrupt();
    }
  }

const int BUFFER_SIZE=512;
/* Check for binary file - pinched from gnuplot */

bool isfilebinary(FILE *fp)
{
  size_t where;
  unsigned char buffer[BUFFER_SIZE];
  int i,odd,len;
  unsigned char *c;
  
  if ((where=ftell(fp))==-1)
    return false;
  
  len=fread((char *)buffer,sizeof(char),BUFFER_SIZE,fp);
  if (len<=0)
    return false;
  
  fseek(fp,where,0);

  for (i=0,odd=0,c=buffer;i<len;i++,c++) {
    if (!*c)
      return true;
    if ((*c & 128) ||/* Meta-characters--we hope it's not formatting */
	(*c == 127)|| /* DEL */
	(*c < 32 && 
	 *c != '\n' && *c != '\r' && *c != '\b' &&
	 *c != '\t' && *c != '\f' && *c != 27 /*ESC*/))
      odd++;
  }
  return(odd*10>len);
}
  
static void read_comment(FILE *fp)
{
  char wr_linebuf[LINE_MAXLEN];
  
  for (;;) {
    const char c=getc(fp);
    if (c=='#' || c=='%')
      fgets(wr_linebuf,LINE_MAXLEN,fp);
    else
      break;
  }
  /* ungetc is dangerous here */
  fseek(fp,-1L,SEEK_CUR);
}

// Assume that integer file is purely real
static size_t read_binint_c(FILE* fp, complex* data, size_t n)
{
  size_t count;
  int c;

  for (count=0;count<n && (c=getc(fp))!=EOF;)
    data[count++]=complex(c,0);
  return count;
}

static bool scanfloat(double& dest, FILE* fp)
{
  char buf[30];
  if (fscanf(fp,"%30s",buf)!=1)
    return false;
  char* tail;
  dest=strtod(buf,&tail);
  if (tail && *tail) {
    char message[75];
    snprintf(message,sizeof(message)-1,"Failed to parse %s as floating point number", buf);
    throw Failed(message);
  }
  return true;
}

static size_t read_ascii_c(FILE* fp, complex* data, size_t n)
{
  double store;
  size_t count;

  for (count=0;count<n;count++) {
    if (scanfloat(store,fp))
      data[count]=complex(store,0);
    else
      break;
  }
  return count;
}

static size_t read_cascii(FILE* fp, complex* data, size_t n)
{
  double storer,storei;
  size_t count;

  for (count=0;count<n;count++) {
    if (scanfloat(storer,fp)) {
      if (scanfloat(storei,fp))
	data[count]=complex(storer,storei);
      else
	throw Failed("Read as complex failed: odd number of entries");
    }
    else
      break;
  }
  return count;
}

static size_t read_binfloat_c(FILE *fp,complex *data,size_t n)
{
  size_t count;
  float store;

  for (count=0;count<n && (fread(&store,sizeof(float),1,fp)==1);)
    data[count++]=complex(store,0);
  return count;
}

static size_t read_cbinfloat(FILE *fp,complex *data,size_t n)
{
  size_t count;
  float store[2];

  for (count=0;count<n && (fread(store,sizeof(float),2,fp)==2);)
    data[count++]=complex(store[0],store[1]);
  return count;
}

template<class T> void read_matrix_(Matrix<T>& a, const char* fname)
{
  FileHandle fp(fname,"r");
  read_matrix(a,fp());
//   FILE* fp=fopen(fname,"r");
//   if (fp==NULL)
//     throw OpenFailed(fname);
//   try {
//     read_matrix(a,fp);
//   } catch (...) {
//     fclose(fp);
//     throw;
//   }
//  fclose(fp);
}

void read_matrix(cmatrix& a,const char* fname)
{
  read_matrix_(a,fname);
}

void read_matrix(rmatrix& a,const char* fname)
{
  read_matrix_(a,fname);
}

void read_matrix(cmatrix& a,FILE *fp)
{
  int rs,cs;
  bool integer,isbinary,iscomplex;

  read_desc_file_(fp,rs,cs,isbinary,iscomplex,integer);

  // refuse block matrix format unless single ASCII column
  if (rs<0) {
//     if (cs!=1)
//       throw Failed("Can't read general ASCII block matrix into complex");
    rmatrix tmp;
    read_matrix(tmp,fp);
    a=tmp;
    return;
  }

  size_t (*read_func)(FILE *,complex *,size_t);

  if (integer)
    read_func= isbinary ? read_binint_c : read_ascii_c;
  else {
    if (iscomplex)
      read_func= isbinary ? read_cbinfloat : read_cascii;
    else
      read_func= isbinary ? read_binfloat_c : read_ascii_c;
  }

  a.create(rs,cs);

  const size_t toread=rs*cs;
  if (read_func(fp,a.vector(),toread)!=toread) {
    a.kill();
    throw Failed("File truncated?");
  }
}
 
typedef size_t (*PRF)(FILE *,double *,size_t);

// Assume that integer file is purely real
static size_t read_binint_r(FILE *fp,double *data,size_t n)
{
  size_t count;

  for (count=0;count<n && (data[count]=getc(fp))!=EOF;count++);
  return count;
}

static size_t read_ascii_r(FILE* fp, double* data, size_t n)
{
  size_t count;

  for (count=0;count<n && scanfloat(data[count],fp);count++);
  return count;
}

static size_t read_binfloat_r(FILE *fp,double *data,size_t n)
{
  size_t count;
  float store;

  for (count=0;count<n && (fread(&store,sizeof(float),1,fp)==1);)
    data[count++]=store;
  return count;
}

const size_t BLK_INIT=1024;

static const double* read_raw_vector(FILE* fp, PRF read_func, int& total)
{
  double *block, *point;

  int blk_size=BLK_INIT;
  if (!(block=(double*)malloc((size_t)(blk_size*sizeof(double)))))
    throw Failed("read_raw_vector: malloc failure");
  
  total=0;
  for (point=block;;) {
    const int read=read_func(fp,point,blk_size);
    total+=read;
    if (read<blk_size)
      break;
    blk_size=total;
    //double block size 
    if (!(block=(double*)realloc((char*)block,(size_t)(2*total*sizeof(double)))))
      throw Failed("read_raw_vector: realloc failure");
    point=block+total;
  }
  return block;
}

static int getnfloats(FILE* fp)
{
  struct stat buf;
  fstat(fileno(fp),&buf);
  const int lenleft=buf.st_size-ftell(fp);
  if (lenleft % sizeof(float))
    throw Failed("Truncated file?");
  return lenleft/sizeof(float);
}

void read_vector(List<double>& a,FILE* fp)
{
  if (isfilebinary(fp)) {
    const int n=getnfloats(fp);
    a.create(n);
    read_binfloat_r(fp,a.vector(),n);
  }
  else
    read_vector_ascii(a,fp);
}

void read_vector_ascii(List<double>& a, FILE* fp)
{
  int n;
  const double* data=read_raw_vector(fp,read_ascii_r,n);
  a.create(n);
  a=const_cast<double*>(data);
  free( (void*)(data));
}

void read_vector(List<complex>& a, FILE* fp)
{
  if (isfilebinary(fp)) {
    int n=getnfloats(fp);
    if (n % 2)
      throw Failed("Odd number of data items - not complex data?");
    n/=2;
    a.create(n);
    read_cbinfloat(fp,a.vector(),n);
  }
  else
    read_vector_ascii(a,fp);
}

void read_vector_ascii(List<complex>& a, FILE* fp)
{
  int n;
  const double* data=read_raw_vector(fp,read_ascii_r,n);
  if (n % 2)
    throw Failed("Odd number of data items - not complex data?");
  copyascomplex(a,BaseList<double>(n,const_cast<double*>(data)));
  free( (void *)data);
}

void read_matrix(rmatrix& a,FILE *fp)
{
  int rs,cs;
  bool integer,isbinary,iscomplex;

  read_desc_file_(fp,rs,cs,isbinary,iscomplex,integer);

  if (iscomplex) 
    throw Failed("Can't read complex matrix into real matrix");

  PRF read_func;

  if (isbinary)
    read_func= integer ? read_binint_r : read_binfloat_r;
  else
    read_func= read_ascii_r;

  if (rs<1) {
    int total;

    const double* where=read_raw_vector(fp,read_func,total);

    if (total % cs) {
      free( (void*)where);
      throw Failed("Truncated file?");
    }
    rs=total/cs;

    a.create(rs,cs,where);
    free( (void *)where);
    return;
  }

  a.create(rs,cs);

  const size_t toread=rs*cs;
  if (read_func(fp,a.vector(),toread)!=toread) {
    a.kill();
    throw Failed("File truncated?");
  }
}

template<class T> inline void read_vector_(List<T> &a,const char *fname)
{
  FileHandle fp(fname,"r");
  read_vector(a,fp());
//   FILE* fp=fopen(fname,"r");
//   if (fp==NULL)
//     throw OpenFailed(fname);
//   try {
//     read_vector(a,fp);
//   } catch (MatrixException& exc) {
//     fclose(fp);
//     throw;
//   }
//   fclose(fp);
}

void read_vector(List<double> &a,const char *fname)
{
  read_vector_(a,fname);
}

void read_vector(List<complex> &a,const char *fname)
{
  read_vector_(a,fname);
}

}//namespace libcmatrix
