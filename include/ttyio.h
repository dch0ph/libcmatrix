#ifndef _ttyio_h_
#define _ttyio_h_

#include <cstdio>
#include "List.h"

namespace libcmatrix {

  void getline(char *line,size_t max,::std::istream& =::std::cin);
  int parseint(const char *);
  float parsefloat(const char *);

  bool isreadable(const char *);
  const char *getbasename(const char *);
  void write_arg(FILE *,int,const char *argv[]);

  List<size_t> getuniqueindices(const char* str,int min,int max);
  size_t pushindices(List<size_t>& inds,const char* cstr,int min,int max);

FILE *getreadfile(const char *,char *,size_t,int =0);
FILE *getreadfile(const char *,int =0);
FILE *getwritefile(const char *,char *,size_t,int =0);
FILE *getwritefile(const char *,int =0);

bool getlogical(const char *);
bool getlogical(const char *,bool);

int getint(const char *);
int getint(const char *,int def);
int getint(const char *,int min,int max);
int getint(const char *,int min,int max,int def);

float getfloat(const char *);
float getfloat(const char *,float min,float max);
float getfloat(const char *,float min,float max,float def);
float getfloat(const char *,float def);

size_t getoption(const char *,const char *);
size_t getoption(const char *,const char *,size_t);

char *getstring(const char *,char *,size_t,const char * =NULL);

bool getlogical(int,const char *argv[],int&,const char *);
bool getlogical(int,const char *argv[],int&,const char *,bool);

int getint(int,const char *argv[],int &,const char *);
int getint(int,const char *argv[],int &,const char *,int);
int getint(int,const char *argv[],int &,const char *,int,int);
int getint(int,const char *argv[],int &,const char *,int,int,int);

float getfloat(int,const char *argv[],int &,const char *);
float getfloat(int,const char *argv[],int &,const char *,float);
float getfloat(int,const char *argv[],int &,const char *,float,float);
float getfloat(int,const char *argv[],int &,const char *,float,float,float);

char *getstring(int,const char *argv[],int &,const char *,char *,size_t,const char * =NULL);

size_t getoption(int,const char *argv[],int &,const char *,const char *);
size_t getoption(int,const char *argv[],int &,const char *,const char *,size_t);

FILE *getreadfile(int,const char *argv[],int &,const char *,char *,size_t,int =0);
FILE *getwritefile(int,const char *argv[],int &,const char *,char *,size_t,int =0);
FILE *getwritefile(int,const char *argv[],int &,const char *,int =0);

} //namespace libcmatrix

#endif
