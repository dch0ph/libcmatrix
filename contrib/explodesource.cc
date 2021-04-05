/* Split source files into sections with common header
   C++ hardly necessary, but may as well since compilation environment must be sane */

#include <iostream>
#include <fstream>

using namespace std;
const int MAXLINE=8192;
const int MAXFNAME=2048;

char* chewwhite(char* ptr)
{
  while (isspace(*ptr))
    ptr++;
  return ptr;
}

void fnamesplit(char*& head, char*& tail, char* in)
{
  head=chewwhite(in);
  char* sep=strrchr(head,'/');
  if (sep) {
    *sep='\0';
    tail=sep+1;
  }
  else {
    tail=head;
    head=NULL;
  }
}


int main(int argc, char* argv[])
{
  if (argc!=2) {
    cerr << argv[0] << ": <input file>\n";
    return 1;
  }

  ifstream istr(argv[1]);
  if (!istr) {
    cerr << argv[0] << ": failed to open: " << argv[1] << '\n';
    return 2;
  }

  char buffer[MAXLINE];
  long endheader;
  
  while (istr.getline(buffer,sizeof(buffer))) {
    if (strncmp(buffer,"//EXPLODE",9)==0)
      break;
    endheader=istr.tellg();
  }

  if (istr.eof()) {
    cerr << argv[0] << ": failed to find an //EXPLODE directive\n";
    return 3;
  }

  char* head;
  char* ftail;
  fnamesplit(head,ftail,argv[1]);

  char newfname[MAXFNAME];
  for (;;) {
    //start with //EXPLODE line in buffer
    char* source=chewwhite(buffer+9);
    if (*source=='\0') {
      cerr << argv[0] << ": missing destination file in //EXPLODE\n";
      return 4;
    }

    if (head) {
      if (snprintf(newfname,MAXFNAME,"%s/%s",head,source)>MAXFNAME) {
	cerr << argv[0] << ": buffer overflow\n";
	return 5;
      }
      source=newfname;
    }

    ofstream ostr(source);
    if (!ostr) {
      cerr << argv[0] << ": failed to open output file: " << source << '\n';
      return 6;
    }
    cout << source << ' ';
    const long curpos=istr.tellg();
    istr.seekg(0); //rewind to header
    while (istr.tellg()<endheader) {
      istr.getline(buffer,MAXLINE);
      ostr << buffer << '\n';
    }
    istr.seekg(curpos); //return to current segment
    while (istr.getline(buffer,MAXLINE)) {
      if (strncmp(buffer,"//EXPLODE",9)==0)
	break;
      ostr << buffer << '\n';
    }
    ostr.close(); //would be done automatically, but no harm in being explicit
    if (istr.eof())
      break;
  }
  cout << '\n';
  return 0;
}
