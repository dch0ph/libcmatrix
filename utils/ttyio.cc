#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <sys/types.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef _MSC_VER
#include <io.h>
#endif
#include <fcntl.h>
#include <sys/stat.h>
#include "ttyio.h"
#include "List.h"
#include "ScratchList.h"
#ifdef LCM_ENABLE_ARGSITER
#include "args_iter.h"
#endif

namespace libcmatrix {

  namespace {
    const char OUTRANGE[]="input out of range";
    
    inline void chewwhite(char*& str) {
      while (isspace(*str)) { str++; } 
    }

    long getunsigned(char** strp,int min,int max)
    {
      char* str=*strp;
      char* nstr;
      const long readval=strtol(str,&nstr,10);
      if (nstr==str) {
	std::cerr << "getindices: failed to parse number: " << str << std::endl;
	return -1;
      }
      if (readval<min || readval>max) {
	std::cerr <<  "getindices: index outside of range" << std::endl;
	return -1;
      }
      str=nstr;
      chewwhite(str);
      *strp=str;
      return readval;
    }
  }

size_t pushindices(List<size_t>& inds,const char* cstr,int min,int max)
{
  if (!cstr)
    throw InvalidParameter("pushindices: NULL input string");
  if (min>max || min<0 || max<0)
    throw InvalidParameter("pushindices: invalid index range");

  char* str=const_cast<char *>(cstr); //arguments to strtoX are not qualified const
  const size_t startcount=inds.length();
  bool firstitem=true;
  for (;;firstitem=false) {
    long endval; 

    chewwhite(str);
    if (*str=='\0')
      return inds.length()-startcount;
    if (*str=='-') { /* item starts with - */
      if (!firstitem) 
	break;
      chewwhite(++str);
      if (*str) {
	std::cerr << "pushindices: invalid character after initial -\n";
	break;
      }
      for (int i=min;i<=max;i++)
	inds.push_back(i); /* special case of "-" */
      return inds.length()-startcount;
    }
    long readval=getunsigned(&str,min,max);
    if (readval<0)
      break;
    if (*str=='-') {
      str++;
      endval=getunsigned(&str,min,max);
      if (endval<0)
	break;
      for (int i=readval;i<=endval;i++) 
	inds.push_back(i);
    }
    else
      inds.push_back(readval);
    if (*str==',')
      str++;
  }
  inds.resize(startcount); //read failed; reset list
  return 0;
}

List<size_t> getuniqueindices(const char* str,int min,int max)
{
  const int maxvals=max-min+1;
  List<size_t> inds(maxvals,mxflag::temporary);
  inds.resize(0);
  pushindices(inds,str,min,max);
  ScratchList<bool> hasval(maxvals,false);
  for (size_t i=inds.length();i--;) {
    bool& has=hasval(inds(i)-min);
    if (has)
      return List<size_t>();
    has=true;
  }
  return inds;
}

void getline(char *line,size_t max,::std::istream& istr)
{
  istr.getline(line,max,'\n');
  if (istr.fail()) {
    std::cerr << "End of input\n";
    throw Failed("getline: end of input");
  }
}

FILE *getwritefile(const char *message,int flags)
{
  char linebuf[LINE_MAXLEN];
  return getwritefile(message,linebuf,LINE_MAXLEN,flags);
}

FILE *getwritefile(const char *message,char *buffer,size_t max,int flags)
{
  if (flags & mxflag::defaultstdout)
    flags|=mxflag::allownull;
  if (!buffer)
    return getwritefile(message,flags);

  FILE *fp;

  for (;;) {
    getstring(message,buffer,max);
    if (!buffer[0]) {
      if (flags & mxflag::allownull)
	return (flags & mxflag::defaultstdout) ? stdout : NULL;
    }
    else {
      if ( !(flags & mxflag::confirmoverwrite) || !isreadable(buffer) || getlogical("File exists, overwrite? ")) {
	if ((fp=fopen(buffer,flags & mxflag::binary ? "wb" : "w")))
	  return fp;
	std::cout << "Cannot open " << buffer << " for writing.\n";
      }
    }
  }
}

FILE *getreadfile(const char *message,char *buffer,size_t max,int flags)
{
  FILE *fp;
  
  for (;;) {
    getstring(message,buffer,max);
    if (!buffer[0]) {
      if (flags & mxflag::allownull) return NULL;
    }
    else {
      if ((fp=fopen(buffer,flags & mxflag::openbinary ? "rb" : "r")))
	return fp;
      std::cout << "Can't find " << buffer << ".\n";
    }
  }
}

FILE *getreadfile(const char *message,int flags)
{
  char linebuf[LINE_MAXLEN];
  return getreadfile(message,linebuf,LINE_MAXLEN,flags);
}

bool getlogical(const char *s)
{
  char linebuf[LINE_MAXLEN];
  for (;;) {
    std::cout << s;
    getline(linebuf,LINE_MAXLEN);
    switch (linebuf[0]) {
    case 'y':case 'Y': return true;
    case 'n':case 'N': return false;
    default:
      std::cout << "Please answer the question.\n";
    }
  }
}	

bool getlogical(const char *s,bool def)
{
  char linebuf[LINE_MAXLEN];
  for (;;) {
    std::cout << s << "[" << (def ? 'y' : 'n') << "] ";
    getline(linebuf,LINE_MAXLEN);
    switch (linebuf[0]) {
    case 'y':case 'Y':case '1': return true;
    case 'n':case 'N':case '0': return false;
    case '\0': return def;
    default:
      std::cout << "Please answer the question.\n";
    }
  }
}	

int getint(const char *s)
{
  char linebuf[LINE_MAXLEN];
  int value;
  
  for (;;) {
    std::cout << s;
    getline(linebuf,LINE_MAXLEN);
    if (sscanf(linebuf,"%i",&value)==1)
      return value;
    printf("Please enter an integer.\n");
  }
}

int getint(const char *s,int imin,int imax)
{
  if (imin>imax)
    throw InvalidParameter("getint: min must be less than max");
  
  for (;;) {
    const int i=getint(s);
    if (i<=imax && i>=imin) 
      return i;
    printf("Between %d and %d please.\n",imin,imax);
  }
}

int getint(const char *s,int defval)
{
  char linebuf[LINE_MAXLEN];
  int value;

  std::cout << s << "[" << defval << "] ";
  getline(linebuf,LINE_MAXLEN);
  return (sscanf(linebuf,"%d",&value)!=1) ? defval : value;
}

int getint(const char *s,int imin,int imax,int idef)
{
  if (imin>imax || idef<imin || idef>imax)
    throw InvalidParameter("getint: min must be less than max");

  for (;;) {
    const int i=getint(s,idef);
    if (i<=imax && i>=imin)
      return i;
    std::cout << "Between " << imin << " and " << imax << " please.\n";
  }
}

float getfloat(const char *s)
{
  char linebuf[LINE_MAXLEN];
  float value;
  
  for (;;) {
    std::cout << s;
    getline(linebuf,LINE_MAXLEN);
    if (sscanf(linebuf,"%g",&value)==1)
      return value;
    std::cout << "Please enter a float.\n";
  }
}

float getfloat(const char *s,float fmin,float fmax)
{
  if (fmin>fmax)
    throw InvalidParameter("getfloat: min must be less than max");
  
  for (;;) {
    const float f=getfloat(s);
    if (f<=fmax && f>=fmin)
      return f;
    std::cout << "Between " << fmin << " and " << fmax << " please.\n";
  }
}

float getfloat(const char *s,float defval)
{
  char linebuf[LINE_MAXLEN];
  float value;
  
  std::cout << s << "[" << defval << "] ";
  getline(linebuf,LINE_MAXLEN);
  return (sscanf(linebuf,"%g",&value)!=1) ? defval : value;
}

float getfloat(const char* s,float fmin,float fmax,float fdef)
{
  if (fmin>fmax || fdef<fmin || fdef>fmax)
    throw InvalidParameter("getfloat: min must be less than max");

  for (;;) {
    const float f=getfloat(s,fdef);
    if (f<=fmax && f>=fmin)
      return f;
    std::cout << "Between " << fmin << " and " << fmax << " please.\n";
  }
}

char *getstring(const char *s,char *t,size_t max,const char *def)
{
  std::cout << s;
  if (def)
    std::cout << "[" << def << "] ";

  getline(t,max);
  if (def && !t[0])
    strncpy(t,def,max);
  return t;
}

size_t getoption(const char *prompt,const char *options)
{
  char linebuf[LINE_MAXLEN];
  for (;;) {
    std::cout << prompt;
    getline(linebuf,LINE_MAXLEN);
    const char charf=linebuf[0];

    if (charf) {
      const char *found=strchr(options,charf);
      if (found)
	return found-options;
    }
    std::cout << "Reply with one character from: " << options << std::endl;
  }
}

size_t getoption(const char *prompt,const char *options,size_t def)
{
  char linebuf[LINE_MAXLEN];
  if (def>=strlen(options)) 
    throw InvalidParameter("getoption");
  const char* found;

  for (;;) {
    std::cout << prompt << "[" << options[def] << "] ";
    getline(linebuf,LINE_MAXLEN);
    const char charf=linebuf[0];

    if (!charf)
      return def;
    if ((found=strchr(options,charf)))
      return found-options;

    std::cout << "Reply with one character from: " << options << std::endl;
  }
}

bool isreadable(const char* fname)
{
  int fd=open(fname,O_RDONLY);
  if (fd>=0) {
    close(fd);
    return true;
  }
  return false;
}

off_t file_length(const char* fname)
{
  struct stat buf;
  return stat(fname,&buf) ? -1 : buf.st_size;
}

off_t file_length(FILE* fp)
{
  struct stat buf;
  if (fstat(fileno(fp),&buf))
    return -1;
  return buf.st_size;
}

const char* getbasename(const char* s)
{
  const char *c;
  if (!s)
    return NULL;
  return (c=strrchr(s,'/')) ? ++c : s;
}

FILE* file_open(const char* s, const char* mode)  
{
  if (!mode || !s)
    throw InvalidParameter("file_open");
  FILE* fp=NULL;
  switch (*mode) {
  case 'r': case 'w':
    fp=fopen(s,mode);
    break;
  case 'a': {
    const int fildes=open(s,O_CREAT | O_WRONLY,0666);
    if (fildes>=0) {
      fp=fdopen(fildes,"wb");
      if (fp)
	fseek(fp,0,SEEK_END);
    }
  }
    break;
  }
  if (!fp) {
    char buf[1024];
    snprintf(buf,sizeof(buf)-20,"file_open: %s (%s)",s,mode);
    throw Failed(buf);
  }
  return fp;
}

int parseint(const char *str)
{
  int value;
  if (sscanf(str,"%i",&value)!=1) {
    std::cerr << "Unparsable as integer: " << str << std::endl;
    throw Failed("Unparsable integer");
  }
  return value;
}

int getint(int argc,const char *argv[],int& count,const char *prompt)
{
  if (count<argc) {
    int value=parseint(argv[count++]);
    std::cout << prompt << value << std::endl;
    return value;
  }
  return getint(prompt);
}

float parsefloat(const char *str)
{
  float value;
  if (sscanf(str,"%g",&value)!=1) {
    std::cerr << "Unparsable as floating point number: " << str << std::endl;
    throw Failed("Unparseable float");
  }
  return value;
}
 
float getfloat(int argc,const char *argv[],int& count,const char *prompt)
{
  if (count<argc) {
    float value=parsefloat(argv[count++]);
    std::cout << prompt << value << std::endl;
    return value;
  }
  return getfloat(prompt);
}


int getint(int argc,const char *argv[],int& count,const char *prompt,int def)
{
  if (count<argc) {
    if (!strcmp(argv[count],"-")) {
      count++;
      std::cout << prompt << def << std::endl;
      return def;
    }
    else
      return getint(argc,argv,count,prompt);
  }
  return getint(prompt,def);
}

float getfloat(int argc,const char *argv[],int& count,const char *prompt,float def)
{
  if (count<argc) {
    if (!strcmp(argv[count],"-")) {
      count++;
      std::cout << prompt << def << std::endl;
      return def;
    }
    else
      return getfloat(argc,argv,count,prompt);
  }
  return getfloat(prompt,def);
}


char* getstring(int argc,const char *argv[],int& count,const char *prompt,char *dest,size_t max,const char *def)
{
  if (count<argc) {
    const char *usep=argv[count++];

    if (!strcmp(usep,"-")) {
      if (def)
	usep=def;
      else {
        dest[0]='\0';
        std::cout << prompt << "<none>" << std::endl;
        return dest;
      }
    }
    std::cout << prompt << usep << std::endl;
    return strncpy(dest,usep,max);
  }
  return getstring(prompt,dest,max,def);
}

int getint(int argc,const char *argv[],int& count,const char *prompt,int min,int max)
{
  if (count<argc) {
    if (max<min)
      throw InvalidParameter("getint: min must be less than max");

    int value=getint(argc,argv,count,prompt);
    if (value<min || value>max) {
      std::cerr << "Value (" << value << ") outside permitted range (" << min << "-" << max << ")" << std::endl;
      throw Failed(OUTRANGE);
    }
    return value;
  }
  return getint(prompt,min,max);
}

float getfloat(int argc,const char *argv[],int& count,const char *prompt,float min,float max)
{
  if (count<argc) {
    if (max<min)
      throw InvalidParameter("getfloat: min must be less than max");
    float value=getfloat(argc,argv,count,prompt);
    if (value<min || value>max) {
      std::cerr << "Value (" << value << ") outside permitted range (" << min << "-" << max << ")" << std::endl;
      throw Failed(OUTRANGE);
    }
    return value;
  }
  return getfloat(prompt,min,max);
}

int getint(int argc,const char *argv[],int& count,const char *prompt,int min,int max,int def)
{
  if (count<argc) {
    if (max<min || def<min || def>max)
      throw InvalidParameter("getint: min must be less than max");
    const int value=getint(argc,argv,count,prompt,def);
    if (value<min || value>max) {
      std::cerr << "Value (" << value << ") outside permitted range (" << min << "-" << max << ")" << std::endl;
      throw Failed(OUTRANGE);
    }
    return value;
  }
  return getint(prompt,min,max,def);
}

float getfloat(int argc,const char *argv[],int& count,const char *prompt,float min,float max,float def)
{
  if (count<argc) {
    if (max<min || def<min || def>max)
      throw InvalidParameter("getfloat: min must be less than max");
    const float value=getfloat(argc,argv,count,prompt,def);
    if (value<min || value>max) {
      std::cerr << "Value (" << value << ") outside permitted range (" << min << "-" << max << ")" << std::endl;
      throw Failed(OUTRANGE);
    }
    return value;
  }
  return getfloat(prompt,min,max,def);
}

FILE *getreadfile(int argc,const char *argv[],int &count,const char *prompt,char *dest,size_t max,int flags)
{
  if (count<argc) {
    if ((flags & mxflag::allownull) && !strcmp(argv[count],"-")) {
      count++;
      return NULL;
    }
    getstring(argc,argv,count,prompt,dest,max);
    FILE *fp=fopen(dest,flags & mxflag::openbinary ? "rb" : "r");
    if (!fp) {
      std::cerr << "Failed to find " << dest << "\n";
      throw Failed("Open failed");
    }
    return fp;
  }
  return getreadfile(prompt,dest,max,flags);
}

FILE *getwritefile(int argc,const char *argv[],int &count,const char *prompt,char *dest,size_t max,int flags)
{
  if (flags & mxflag::defaultstdout) 
    flags|=mxflag::allownull;
  if (count<argc) {
    if (!strcmp(argv[count],"-")) {
      if (!(flags & mxflag::allownull)) {
	std::cerr << "Filename must be supplied.\n";
	throw Failed("Invalid filename");
      }
      count++;
      return flags & mxflag::defaultstdout ? stdout : NULL;
    }
    const char *usep=argv[count++];
    if (dest)
      strncpy(dest,usep,max);
    FILE *fp=fopen(usep,flags & mxflag::openbinary ? "wb" : "w");
    if (!fp) {
      std::cerr << "Failed to open " << usep << "\n";
      throw Failed("Open failed");
    }
    std::cout << prompt << usep << std::endl;
    return fp;
  }
  return getwritefile(prompt,dest,max,flags);
}

FILE *getwritefile(int argc,const char *argv[],int &count,const char *prompt,int flags)
{
  return getwritefile(argc,argv,count,prompt,NULL,0,flags);
}

  namespace {
    bool parse_logical(const char *arg,const char *prompt,bool def,bool have_def)
    {
      bool value;
      
      switch (arg[0]) {
      case 'y': case 'Y': value=true; break;
      case 'n': case 'N': value=false; break;
      case '-':
	if (have_def)
	  value=def;
	else {
	  std::cerr << "(" << prompt << "): No default response." << std::endl;
	  throw Failed("No default response");
	}
	break;
      default:
	std::cerr << "(" << prompt << ") unparseable as logical: " << arg << std::endl;
	throw Failed("Unparseable logical");
      }
      std::cout << prompt << (value ? 'y' : 'n') << std::endl;
      return value;
    }
  }

bool getlogical(int argc,const char *argv[],int &count,const char *prompt)
{
  return  count<argc ? parse_logical(argv[count++],prompt,false,false) : getlogical(prompt);
}

bool getlogical(int argc,const char *argv[],int &count,const char *prompt,bool def)
{
  return  count<argc ? parse_logical(argv[count++],prompt,def,true) : getlogical(prompt,def);
}

  namespace {
    size_t parse_option(char arg,const char *prompt,const char *options,size_t def,bool have_def,::std::ostream& ostr,bool haveostr)
    {
      if (have_def && (def>=strlen(options)))
	throw InvalidParameter("parse_option: default out of range!");
      
      if (arg=='-') {
	if (have_def)
	  return def;
	std::cerr << "(" << prompt << "): No default response." << std::endl;
	throw Failed("No default response");
      }
      
      const char *found=strchr(options,arg);
      
      if (!found) {
	std::cerr << "(" << prompt << "): input must be one of : " << options << std::endl;
	throw Failed("Invalid option");
      }
      
      if (haveostr)
	ostr << prompt << arg << std::endl;
      return found-options;
    }    
  }

size_t getoption(int argc,const char *argv[],int &count,const char *prompt,const char *options,size_t def)
{
  return count<argc ? parse_option(argv[count++][0],prompt,options,def,true,::std::cout,true) : getoption(prompt,options,def);
}

size_t getoption(int argc,const char *argv[],int &count,const char *prompt,const char *options)
{
  return count<argc ? parse_option(argv[count++][0],prompt,options,0,false,::std::cout,true) : getoption(prompt,options);
}

void
write_arg(FILE* fp,int argc,const char* argv[])
{
  fputc('%',fp);
  for (int i=0;i<argc;i++)
    fprintf(fp,"%s ",argv[i]);
  fputc('\n',fp);
}

#ifdef LCM_ENABLE_ARGSITER

 template<> bool readstream(bool& v,std::istream& in) {
   char tmp[10];
   in.width(10);
   if (in >> tmp) {
     switch (tmp[0]) {
     case 'Y':case 'y': v=true; return true;
     case 'N':case 'n': v=false; return true;
     }
   }
   return false;
 }

args_iter::args_iter(size_t start,size_t end,const char *argv_[]) :  cpos(start), epos(end), argv(argv_), istr(::std::cin), ostr(::std::cout), haveostr(true) { init(); }

args_iter::args_iter(size_t start,size_t end,const char *argv_[],::std::istream& istr_) : cpos(start), epos(end), argv(argv_), istr(istr_), ostr(::std::cout), haveostr(false) { init(); }

args_iter::args_iter(size_t start,size_t end,const char *argv_[],::std::istream& istr_, ::std::ostream& ostr_) : cpos(start), epos(end), argv(argv_), istr(istr_), ostr(ostr_), haveostr(true) { init(); }

args_iter::args_iter() : cpos(0), epos(0), istr(::std::cin), ostr(::std::cout), haveostr(true) { init(); }

args_iter::args_iter(::std::istream& istr_) : cpos(0), epos(0), istr(istr_), ostr(::std::cout), haveostr(false) { init(); }

args_iter::args_iter(::std::istream& istr_, ::std::ostream& ostr_) : cpos(0), epos(0), istr(istr_), ostr(ostr_), haveostr(true) { init(); }

std::istream& args_iter::operator()(char *dest,size_t min,size_t max,const char *prompt,const char *def)
{
  do {
    if (doprompt()) {
      ostr << prompt;
      if (def)
	ostr << " [" << def << "] ";
    }
    std::istream& cistr=update();

    cistr.width(max);
    
    const bool ok=(cistr >> dest);
    next();
    if (!ok || dest[0]=='\0' || !strcmp(dest,"-")) {
      if (def)
	strncpy(dest,def,max);
      else
	dest[0]='\0';
    }
    if (strlen(dest)<min)
      std::cerr << "Input must have at least " << min << "character" << ((min>1) ? "s" : "") << std::endl;
    else {
      if (!doprompt()) 
	ostr << prompt << (dest[0] ? dest : "<none>") << ::std::endl;
      return cistr;
    }
  } while (doprompt());
  throw Failed("getstring: bad input");
}

std::istream& args_iter::operator()(size_t& which, const char *prompt,const char *options)
{
  char c;
  if (!doprompt()) {    
    std::istream& cistr=update();
    if (cistr >> c) {
      next();
      which=parse_option(c,prompt,options,0,false,haveostr ? ostr : ::std::cout,haveostr);
      return cistr;
    }
    throw Failed("args_iter::option");
  }

  for (;;) {
    ostr << prompt;
    std::istream& cistr=update();
    if (cistr >> c) {
      const char *found=strchr(options,c);
      if (found) {
	which=found-options;
	return cistr;
      }
    }
    std::cerr << "Respond with one character of " << options << std::endl;
  }
}
    
std::istream& args_iter::operator()(size_t& which, const char *prompt,const char *options, size_t def)
{
  if (def>=strlen(options)) 
    throw InvalidParameter("getoption");

  char c;

  if (!doprompt()) {    
    std::istream& cistr=update();
    if (cistr >> c) {
      next();
      which=parse_option(c,prompt,options,def,true,haveostr ? ostr : ::std::cout,haveostr);
      return cistr;
    }
    throw Failed("args_iter::option");
  }

  for (;;) {
    ostr << prompt << "[" << options[def] << "] ";
    std::istream& cistr=update();
    if (cistr >> c) {
      const char *found=strchr(options,c);
      if (found) {
	which=found-options;
	return cistr;
      }
      std::cerr << "Respond with one character of " << options << std::endl;
      continue;
    }
    else {
      which=def;
      return cistr;
    }
  }
}

  namespace {
    bool iswrite(const char *mode)
    {
      switch (mode[0]) {
      case 'r': return false;
      case 'w': return true;
      }
      throw InvalidParameter("File open mode must start with r or w");
    }
  }

std::istream& args_iter::operator()(FILE* & fp,const char *prompt,const char *mode,char *dest,size_t max,int flags)
{
  const bool isw=iswrite(mode);
  const bool allown=(flags & mxflag::allownull);
  const bool confirm=isw ? (flags & mxflag::confirmoverwrite) : false;

  do {
    std::istream& cistr=operator()(dest,allown ? 0 : 1,max,prompt,allown ? "" : NULL);
    if (dest[0]=='\0') {
      fp=NULL;
      return cistr;
    }
    if (confirm && isreadable(dest) && (!doprompt() || !getlogical("File exists, overwrite? ")))
      std::cerr << "Will not overwrite " << dest << std::endl;
    else {
      if ((fp=fopen(dest,mode)))
	return cistr;
      std::cerr << "Failed to open " << dest << std::endl;
    }
  } while (doprompt());
  throw Failed("File open");
}

std::istream& args_iter::update() {
  if (cpos<epos)
    attach(argv[cpos]);
  else {
    getline(tmp,LINE_MAXLEN,istr);
    attach(tmp);
  }
  return operator()();
}

//would use str with istringstream, but this seems broken...       
void args_iter::attach(const char *inp)
{
  if (cstrp)
    delete cstrp;
  cstrp=new LCM_STRTYPE(inp);
}

void args_iter::init() {
  cstrp=NULL;
}

//LCM_ENABLE_ARGSITER
#endif

}//namespace libcmatrix
