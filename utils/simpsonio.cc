/* read and write (ASCII) simpson files */

#include <cstdio>
#include <cctype>
#include <cstdlib>
#include "List.h"
#include "Matrix.h"
#include "cmatrix_complex.h"
#include "simpsonio.h"
#include "ScratchList.h"

namespace libcmatrix {

static char errmesg[LINE_MAXLEN];

inline void parse(double& v,const char *val,size_t lineno)
{
  if (sscanf(val,"%lg",&v)!=1) {
    LCM_SNPRINTF(errmesg,sizeof(errmesg),
	     "read_simpson: unable to parse '%s' as double at line %i",val,(int)lineno);
    throw Failed(errmesg);
  }
}

  namespace {
    void striptrail(char* buf) {      
      char* last=buf+strlen(buf)-1;
      while (isspace(*last) && (last>=buf)) 
	last--;
      last[1]='\0';
    }
  }

void read_simpson(List<complex>& a,simpsonFD& fd,FILE *fp)
{
  char buf[LINE_MAXLEN];
  size_t lineno=1;

  if (!fgets(buf,sizeof(buf),fp)) {
    LCM_SNPRINTF(errmesg,sizeof(errmesg),"read_simpson: unable to read line %i",(int)lineno);
    throw Failed(errmesg);
  }
  striptrail(buf);
  if (strcmp(buf,"SIMP")!=0) { //don't check for newline (may be DOS)
    LCM_SNPRINTF(errmesg,sizeof(errmesg),"read_simpson: invalid file type: %s",buf);
    throw Failed(errmesg);
  }

  bool got_sw=false;
  bool got_np=false;
  bool got_data=false;

  while (!feof(fp) && !got_data) {
    lineno++;
    if (!fgets(buf,sizeof(buf),fp)) {
      if (got_sw && got_data && got_np) break;
      strcpy(errmesg,"read_simpson: missing ");
      buf[0]=0;
      if (!got_sw) strcat(buf," SW=<spectral width>,");
      if (!got_data) strcat(buf," <data section>,");
      if (!got_np) strcat(buf," NP=<number of points>,");
      throw Failed(errmesg);
    }

    striptrail(buf);

    char* com=buf;
    while (isspace(*com)) com++;
    if (!(*com)) continue; //ignore blank lines

    char *val=strchr(com,'=');
    if (val) {
      *val++='\0';
      
      if (strcmp(com,"NP")==0) {
	fd.np=atoi(val);
	if (fd.np<=0)
	  throw Failed("read_simpson: number of points must be larger than zero");
	got_np=true;
      }
      else if (strcmp(com,"NI")==0) {
	fd.ni=atoi(val);
	if (fd.ni<=0)
	  throw Failed("read_simpson: number of increments must be larger than zero");
      }
      else if (strcmp(com,"SW")==0) {
	parse(fd.sw,val,lineno);
	got_sw=true;
      }
      else if (strcmp(com,"SW1")==0)
	parse(fd.sw1,val,lineno);
      else if (strcmp(com,"TYPE")==0) {
	if (strcmp(val,"FID")==0)
	  fd.type=FD_TYPE_FID;
	else if (strcmp(val,"SPE")==0) 
	  fd.type=FD_TYPE_SPE;
	else {
	  LCM_SNPRINTF(errmesg,sizeof(errmesg),"read_simpson: unknown file type: %s",val);
	  throw Failed(errmesg);
	}
      }
      else if (strcmp(com,"FORMAT")==0) {
	if (strcmp(val,"TEXT")==0)
	  fd.format=FD_FORMAT_TEXT;
	else if (strcmp(val,"BINARY")==0)
	  fd.format=FD_FORMAT_BINARY;
	else {
	  LCM_SNPRINTF(errmesg,sizeof(errmesg),"read_simpson: unknown format: %s",val);
	  throw Failed(errmesg);
	}
      }
      else if (strcmp(com,"PREC")==0) {
	if (strcmp(com,"SINGLE")==0)
	  fd.prec=FD_PREC_SINGLE;
	else if (strcmp(com,"DOUBLE")==0)
	  fd.prec=FD_PREC_DOUBLE;
	else {
	  LCM_SNPRINTF(errmesg,sizeof(errmesg),"read_simpson: unknown precision: %s",val);
	  throw Failed(errmesg);
	}
      }
      else if (strcmp(com,"REF")==0)
	parse(fd.ref,val,lineno);
      else if (strcmp(com,"REF1")==0)
	parse(fd.ref1,val,lineno);
      else if (strcmp(com,"SFRQ")==0)
	parse(fd.sfrq,val,lineno);
      else if (strcmp(com,"SFRQ1")==0)
	parse(fd.sfrq1,val,lineno);
      //ignore unrecognised tags
    }
    else { //no '='
      if (strcmp(com,"DATA"))
	throw Failed("read_simpson: data section does not start with DATA");
      if (!got_np)
	throw Failed("read_simpson: number of data points not set!");
      if (fd.format!=FD_FORMAT_TEXT)
	throw Failed("read_simpson: binary data read unsupported");
      const size_t ntot=fd.np * (fd.ni ? fd.ni : 1);

      a.create(ntot);
      double re,im;
      for (size_t i=0;i<ntot;i++) {
	lineno++;
	if (fscanf(fp,"%lg%lg",&re,&im)!=2) {
	  LCM_SNPRINTF(errmesg,sizeof(errmesg),"read_simpson: failed to read point %i at line %i",(int)i+1,(int)lineno);
	  throw Failed(errmesg);
	}
	a(i)=complex(re,im);
      }
      if (!fscanf(fp,"%s",buf) || strcmp(buf,"END")) {
	LCM_SNPRINTF(errmesg,sizeof(errmesg),"read_simpson: data section finishes with '%s' not END at line %i",buf,(int)lineno);
	throw Failed(errmesg);
      }
      got_data=true;
    }
  }
}

void read_simpson(List<complex>& a,simpsonFD& fd,const char* fname)
{
   FILE *fp=fopen(fname,"rb");
   if (!fp) 
     throw OpenFailed("read_simpson");
   try {
     read_simpson(a,fd,fp);
   } catch (...) {
     fclose(fp);
     throw;
   }
}

void read_simpson(List<complex>& a,const char* fname)
{
  simpsonFD fd;
  read_simpson(a,fd,fname);
}

void read_simpson(Matrix<complex>& a,simpsonFD& fd,const char* fname)
{
  FILE *fp(fopen(fname,"rb"));
  if (!fp)
    throw OpenFailed("read_simpson");
  try {
    List<complex> tmp(mxflag::temporary);
    read_simpson(tmp,fd,fp);
    a.create(fd.ni ? fd.ni : 1,fd.np,tmp.begin());
  }
  catch (MatrixException&) {
    fclose(fp);
    throw;
  }
  fclose(fp);
}

void read_simpson(Matrix<complex>& a,const char* fname)
{
  simpsonFD fd;
  read_simpson(a,fd,fname);
}

Warning<> simpson_controller::ignoring_binary_flag_warning("write_simpson: ignoring binary format flag",&lcm_io_warning);

void write_simpson(FILE* fp,const simpsonFD& fd)
{
  fprintf(fp,"SIMP\n");
  fprintf(fp,"NP=%i\n",fd.np);
  fprintf(fp,"SW=%.12g\n",fd.sw);
#if LCM_ENABLE_SFRQ
  if (fd.sfrq)
    fprintf(fp,"SFRQ=%.12g\n",fd.sfrq);
#endif
  if (fd.ref) fprintf(fp,"REF=%.12g\n",fd.ref);
  if (fd.ni) fprintf(fp,"NI=%i\n",fd.ni);
  if (fd.sw1) fprintf(fp,"SW1=%.12g\n",fd.sw1);
#if LCM_ENABLE_SFRQ
  if (fd.sfrq1)
    fprintf(fp,"SFRQ1=%.12g\n",fd.sfrq1);
#endif
  if (fd.ref1)
    fprintf(fp,"REF1=%.12g\n",fd.ref1);
  fprintf(fp,"TYPE=%s\n",(fd.type==FD_TYPE_FID) ? "FID" : "SPE");
  if (fd.format==FD_FORMAT_BINARY)
    simpson_controller::ignoring_binary_flag_warning.raise();
  fprintf(fp,"FORMAT=TEXT\n");
  fprintf(fp,"DATA\n");
  
  const size_t ntot=fd.np * (fd.ni ? fd.ni : 1);
  for (size_t i=0;i<ntot;i++) {
    const complex& ai=fd.data[i];
    fprintf(fp,"%.9g %.9g\n",real(ai),imag(ai));
  }
  fprintf(fp,"END\n");
}  

void write_simpson(const char* fname, const simpsonFD& fd)
{
  FILE* fp=fopen(fname,"wb"); // open in binary mode otherwise Windows-based programs will insert CR which confuses simplot
  if (fp==NULL) {
    LCM_SNPRINTF(errmesg,sizeof(errmesg),"write_open: failed to open '%s' for writing",fname);
    throw Failed(errmesg);
  }
  write_simpson(fp,fd);
  fclose(fp);
}

void write_simpson(const char* fname, const BaseList<complex>& a, double sw, bool isspectrum)
{
  write_simpson(fname,simpsonFD(a,sw,isspectrum));
}

void write_simpson(const char* fname, const Matrix<complex>& a, double sw, double sw1, bool isspectrum)
{
  write_simpson(fname,simpsonFD(a,sw,sw1,isspectrum));
}

List<char> insert_simpson(const char* fname,const char* ins, char sep)
{
  int nsplit=strlen(fname)-4;
  List<char> tmp_(nsplit+6+strlen(ins));
  char* tmp(tmp_.vector());
  const char* split=fname+nsplit;
  if ((nsplit>0) && (*split=='.')) {
    memcpy(tmp,fname,nsplit);
    if (sep)
      tmp[nsplit++]=sep;
    sprintf(tmp+nsplit,"%s%s",ins,split);
  }
  else {
    if (sep)
      sprintf(tmp,"%s%c%s",fname,sep,ins);
    else
      sprintf(tmp,"%s%s",fname,ins);
  }
  return tmp_;
}

simpson_controller::simpson_controller(const char* fname, size_t maxv)
  : max_(maxv)
{
  if (maxv==0)
    throw InvalidParameter("simpson_controller: rows cannot be zero!");
  findsplit(fname,4);
  if (maxv>=1000)
    throw InvalidParameter("simpson_controller: row output cannot be used for >999 rows");
  if (maxv>100)
    prec_=3;
  else if (maxv>10)
    prec_=2;
  else
    prec_=1;
}

void simpson_controller::findsplit(const char* fname, size_t len)
{
  nsplit=strlen(fname)-4;
  tmp_.create(nsplit+4+len);
  const char* lsplit=fname+nsplit;
  insert = (nsplit>0) && (*lsplit=='.');
  if (insert)
    split=lsplit;
  else
    nsplit+=4; //point to end of string
  strcpy(tmp_.vector(),fname);
}

simpson_controller::simpson_controller(const char* fname, const BaseList<const char*>& labelsv)
  : prec_(0)
{
  if (labelsv.empty())
    throw InvalidParameter("simpson_controller: empty labels list");
  max_=labelsv.size();
  labels_.create(max_);
  size_t maxlen=0;
  for (size_t i=max_;i--;) {
    labels_(i)=labelsv(i);
    const size_t len=strlen(labelsv(i));
    if (len==0) 
      throw InvalidParameter("simpson_controller: labels cannot be empty");
    if (len>maxlen)
      maxlen=len;
  }
  findsplit(fname,maxlen+1);
}

void simpson_controller::write(size_t n, const BaseList<complex>& a, double sw, bool isspectrum) const
{
  const simpsonFD spec(a,sw,isspectrum);
  write(n,spec);
}

void simpson_controller::write(size_t n, const simpsonFD& FD) const
{
  if (n>=max_)
    throw InvalidParameter("simpson_controller: row exceeds previously declared maximum");
  char* tmp=tmp_.vector()+nsplit;
  if (labels_.empty()) {
    if (insert)
      sprintf(tmp,"%0*i%s",prec_,int(n),split.c_str());
    else
      sprintf(tmp,"%0*i",prec_,int(n));
  }
  else {
    if (insert)
      sprintf(tmp,"%s%s",labels_(n).c_str(),split.c_str());
    else
      strcpy(tmp,labels_(n).c_str());
  }
  write_simpson(tmp_.vector(),FD);
}

void simpson_controller::write_rows(const Matrix<complex>& a, double sw, bool isspectrum, double sfrq)
{
  if (a.rows()!=max_)
    throw InvalidParameter("simpson_controller: data matrix of unexpected size");
  for (size_t i=a.rows();i--;) {
    if (sfrq) {
      simpsonFD fd(a.row(i),sw,isspectrum);
      fd.sfrq=sfrq;
      write(i,fd);
    }
    else 
      write(i,a.row(i),sw,isspectrum);
  }
}

void write_simpson_rows(const char* fname, const Matrix<complex>& a, double sw, bool isspectrum)
{
  simpson_controller ctrl(fname,a.rows());
  ctrl.write_rows(a,sw,isspectrum);
}

void write_simpson_rows(const char* fname, const Matrix<complex>& a, double sw, bool isspectrum, double sfrq)
{
  simpson_controller ctrl(fname,a.rows());
  ctrl.write_rows(a,sw,isspectrum,sfrq);
}

void write_simpson_rows(const char* fname, const Matrix<complex>& a, const BaseList<const char*>& labelsv, double sw, bool isspectrum)
{
  simpson_controller ctrl(fname,labelsv);
  ctrl.write_rows(a,sw,isspectrum);
}

void write_simpson_rows(const char* fname, const Matrix<complex>& a, const BaseList<const char*>& labelsv, double sw, bool isspectrum, double sfrq)
{
  simpson_controller ctrl(fname,labelsv);
  ctrl.write_rows(a,sw,isspectrum,sfrq);
}

}//namespace libcmatrix
