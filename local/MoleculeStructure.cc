/* This contains completely random stuff used in some programs */

#include <cstdio>
#include <ctype.h>
#include <cmath>
#include "MoleculeStructure.h"

namespace libcmatrix {

  using ::std::sqrt;
  using ::std::cos;

  namespace {
    class file_reader {
    public:
      explicit file_reader(const char* fname) 
	: fp(fopen(fname,"r"))
      {
	if (fp==NULL)
	  throw Failed("file_reader: open");
      }
      FILE* operator()() { return fp; }
      ~file_reader() { fclose(fp); }
    private:
      FILE* fp;
    };
    
    inline void checknewline(FILE* fp)
    {
      const char ch=getc(fp);
      if (ch!='\n')
	ungetc(ch,fp);
    }

    inline void push_ifnonnull(List<size_t>& atomlist, char* start)
    {
      const size_t atom=atoi(start); //atom.serial number	
      if (atom)
	atomlist.push_back(atom);
    }
    
   inline char* chewwhite(char* str) {
     while (isspace(*str)) { str++; } 
     return str;
   }

   inline const char* chewwhite(const char* str) {
     while (isspace(*str)) { str++; } 
     return str;
   }
    
    void getfixedline(char* ptr, size_t max, FILE* fp)
    {
      char* const end=ptr+max-1;
      char ch;
      while (ptr<end) {
	ch =getc(fp);
	switch (ch) {
	case '\n': case EOF:
	  *ptr='\0';
	  return;
	case '\r':
	  checknewline(fp);
	  *ptr='\0';
	  return;
	  // 	case ' ':
 	  //ch='\0'; // replace spaces with NULL
 	  //break;
	default:
	  ch=toupper(ch); //enforce upper case
	  break;
	}
	*ptr++=ch;
      }
      //record is full, swallow rest of line
      do {
	ch=getc(fp);
      } while( (ch!='\n') && (ch!='\r') && (ch!=EOF) );
      
      if (ch=='\r')
	checknewline(fp);

      *ptr ='\0';
    }
  }

  void MoleculeStructure::read_molecule(const char* fname, char selatom)
  {
    file_reader freader(fname);
    read_molecule(freader(),selatom);
  }

  void MoleculeStructure::read_PDB(const char* fname, char selatom, bool findunitcell)
  {
    file_reader freader(fname);
    read_PDB(freader(),selatom,findunitcell);
  }

  void MoleculeStructure::read_PDB_single(const char* fname, int slice, char selatom, bool findunitcell)
  {
    file_reader freader(fname);
    read_PDB_single(freader(),slice,selatom,findunitcell);
  }

  void MoleculeStructure::read_PDB_single(FILE* fp, int slice, char selatom, bool findunitcell)
  {
    if (slice < 1)
      throw InvalidParameter("read_PDB_single: model number must be >0");
    read_PDB_(fp, slice, selatom, findunitcell);
  }

  bool PDBatom::istype(const char* ref) const //compare atom types from user's input and PDB-file (e.g. C and C2 are the same, F2 and F2A are the same etc.)
  {
    if (!ref || (*ref=='\0'))
      throw InvalidParameter("PDB::istype");

    const int n=strlen(ref)-1;
    
    if (ref[n]!='*')
      return (strcmp(ref,type)==0);

    if (n<1)
      throw InvalidParameter("PDB::istype");

    if (strncmp(ref,type,n)) //lead doesn't match?
      return false;

    if (type[n]=='\0')
      return true; //exact match

    const bool trail_is_number=isdigit(type[n]);
    const bool matching_number=isdigit(ref[n-1]);

    return matching_number ? !trail_is_number : trail_is_number;
  }

  List<size_t> MoleculeStructure::find(const char* type, bool incell0) const
  {
    List<size_t> result;
    const size_t start = incell0 ? cell0_start() : 0;
    const size_t lend = incell0 ? cell0_end() : size();
    for (size_t i=start;i<lend;i++) {
      if ((*this)(i).istype(type))
	result.push_back(i);
    }
    return result;
  }

  Warning<> MoleculeStructure::connectivity_warning("MoleculeStructure: connectivity information refers to (pseudo)atoms that were not imported",&lcm_io_warning);

  void MoleculeStructure::tidy_connections()
  {
    if (ncells_==0)
      throw Failed("MoleculeStructure::tidy_connections: unit cell size undetermined");

    size_t maxs=0;
    for (size_t i=size();i--;) {
      const size_t snumber=(*this)(i).serial;
      if (snumber>maxs)
	maxs=snumber;
    }
    {
      List<int>& mutable_index=const_cast<MoleculeStructure*>(this)->serial_to_internal_;
      mutable_index.create(maxs+1,-1);
      vector3 total(0,0,0);
      List<vector3> totals(ncells_,vector3(0,0,0));
      for (size_t i=size();i--;) {
	const PDBatom& curatom((*this)(i));
	mutable_index(curatom.serial)=i;
	total+=curatom.coords;
	totals(i / cellsize_)+=curatom.coords;
      }
      centre_=total*(1.0/size());
      
      if (ncells_>1) {
	const double scalef=1.0/cellsize_;
	double mindist=1e30;
	for (size_t i=ncells_;i--;) {
	  vector3 lcentre(totals(i)*scalef);
	  const double dist=(lcentre-centre_).length();
	  //	if (verbose>1)
	  // std::cout << "Centre of cell " << i << " at " << lcentre << " (distance=" << dist << " A)\n";
	  if (dist<mindist) {
	    mindist=dist;
	    cell0_=i;
	  }
	}
	if (verbose_)
	  std::cout << "Cell " << cell0_ << " is closest to origin at " << centre_ << " (distance of " << mindist << " A)\n";
      }
      else
	cell0_=0;
    }
    
    if (!has_connections())
      return;
    
    {
      List<int>& mutable_index=const_cast<MoleculeStructure*>(this)->whichmolecule_;
      mutable_index.create(size(),-1);
      bool warn=false;
      for (size_t m=connections.size();m--;) {
	const List<size_t>& curmol(connections(m));
	for (size_t n=curmol.size();n--;) {
	  const int index(serial_to_internal_(curmol(n)));
	  if (index>=0) //not all atoms in CONECT may be present in structure
	    mutable_index(index)=m;
	  else
	    warn=true;
	}
      }
      if (warn)
	connectivity_warning.raise();
    }
  }

  Warning<> MoleculeStructure::missing_axis_warning("read_molecule: axis definition missing.",&lcm_io_warning);

  void MoleculeStructure::read_molecule(FILE* fp, char selatom)
  {
    char scr[128];
    
    if (fgets(scr,sizeof(scr),fp)==NULL)
      throw Failed("MoleculeStructure::read_molecule: failed to read 1st line - file corrupt?");
    clear();
    
    int axis_count=0;
    vector3 pos;
    size_t serial=1;
    
    for (int lico=1;fgets(scr,sizeof(scr),fp);lico++) {
      const int items=sscanf(scr+10,"%lg %lg %lg",&(pos.x),&(pos.y),&(pos.z));
      if (items==0 || items==EOF)
	break;
      if (items<3) {
	sprintf(scr,"read_molecule: incomplete data line/corrupted file on line %i:\n",lico);
	throw Failed(scr);
      }
            
      switch (scr[0]) {
      case 'o':
	axis_count++;
	axisdef[0]=pos;
	break;
      case 'a':case 'b': case 'c':
	axis_count++;
	axisdef[1+scr[0]-'a']=pos;
	break;
      case ' ':
	break;
      default:
	if (!selatom || (scr[0]==selatom))
	  push_back(PDBatom(pos,scr[0],serial));
	serial++; //counts atom in *input* file
      }
    }      

    switch (axis_count) {
    case 0:
      missing_axis_warning.raise();
      break;
    case 4:
      haveaxes_=true;
      break;
    default:
      throw Failed("read_molecule: axis definition corrupt");
    }
  }

// double cos_length(double r1,double r2,double ang12)
// {
//   return sqrt(r1*r1+r2*r2-2*r1*r2*cos(ang12));
// }

// double cos_angle(double r12,double r13,double r23)
// {
//   return acos((r12*r12+r13*r13-r23*r23)/(2*r12*r13));
// }

 void MoleculeStructure::addConnections(const BaseList<size_t>& newcon)
  {
    for (size_t i=newcon.size();i--;) {  //looking for what molecule to add
      const size_t curatom=newcon(i);
      for (size_t m=connections.size();m--;) {
	List<size_t>& curmol(connections(m));
	if (ismember(curmol,curatom)) {
	  for (size_t j=newcon.size();j--;) { //add all new connections
	    if (!ismember(curmol,newcon(j)))
	      curmol.push_back(newcon(j));
	  }
	  return;
	}
      }
    }
    //molecule is not found i.e. all connections are new
    connections.push_back(newcon);
  }

  Warning<> MoleculeStructure::division_warning("atoms incommensurate with cell size",&lcm_io_warning);

  void MoleculeStructure::read_PDB(FILE* fp, char selatom, bool findunitcell)
  {
    read_PDB_(fp, 0, selatom, findunitcell);
  }

  static void copystring(char* dest, const char* buffer, size_t start, size_t end)
  {
    const char* labelstart=chewwhite(buffer+start);
    const char* labelend=buffer+end;
    while ((labelstart<labelend) && (*labelstart!=' '))
      *dest++=*labelstart++;
    *dest='\0';
  }

//rough and ready parser for PDB-files (from standard guide)
//see www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html


  void MoleculeStructure::read_PDB_(FILE* fp, int slice, char selatom, bool findunitcell)
  {
    char buffer[82];
    const double cutoff=9000.0;
    
    PDBatom single_entry;
    List<size_t> atomlist;
    bool foundslice = (slice < 1); //!< always true if not selected 1 model
    
    cellsize_=0;
    while (!feof(fp)) { //analyse each string
      getfixedline(buffer,sizeof(buffer),fp);
      if (strncmp(buffer,"END   ",6)==0) {//!< fixed 2017-11-07 from ENDM, which seemed to be neither END nor ENDMDL
	if (slice > 0) {
	  if (!foundslice)
	    throw Failed("Failed to find requested model in PDB file");
	  else
	    throw Failed("END encountered before ENDMDL!");
	}
	break;
      }
      if ((slice > 0) && (strncmp(buffer, "ENDMDL",6)==0) && foundslice) //!< stop if at end of selected model
	break;
      
      if (!strncmp(buffer,"CRYST1",6)) { //!< CRYST1 occurs before MODEL
	lattice_info_.a = atof(buffer+6);
	lattice_info_.b = atof(buffer+15);
	lattice_info_.c = atof(buffer+24);
	if (lattice_info_.a * lattice_info_.b * lattice_info_.c <= 0.0) {
	  std::cerr << "Warning: invalid lattice parameters find in CRYST1 line - ignoring uit cell information\n";
	  lattice_info_.a = 0.0; // indicate no unit cell info
	}
	else {
	  lattice_info_.alpha = atof(buffer+33);
	  lattice_info_.beta = atof(buffer+40);
	  lattice_info_.gamma = atof(buffer+47);
	}
	char sbuffer[12];
	copystring(sbuffer, buffer, 55, 66);
	space_group_info_.sgroup = sbuffer;
	space_group_info_.Z = atoi(buffer+66);
	continue;
      }
      if (slice > 0) {
	if (strncmp(buffer, "MODEL ",6)==0) {
	  int mdlnum = atoi(buffer+10);
#ifndef NDEBUG
	  std::cout << "Found model: " << mdlnum << '\n';
#endif
	  foundslice = (mdlnum == slice);
	}
	if (!foundslice)
	  continue;
      }
      if (!strncmp(buffer,"HETATM",6) || !strncmp(buffer,"ATOM  ",6)) {
	if (buffer[13]=='Q')
	  continue; //ignore pseudo-atoms
	//atoi etc. don't check if record has spilled outside allocated columns
	
	copystring(single_entry.type,buffer,12,20);
	if (selatom && (single_entry.type[0]!=selatom)) //not selected?
	  continue;
	single_entry.serial=atoi(buffer+6);  //enter (index) (PDB serial number);
	single_entry.residue=atoi(buffer+22);
	single_entry.coords.x=atof(buffer+30);
	single_entry.coords.y=atof(buffer+38);
	single_entry.coords.z=atof(buffer+46);
	if ((single_entry.coords.x>cutoff) && (single_entry.coords.y>cutoff) && (single_entry.coords.z>cutoff))
	  continue; //ignore XPLOR pseudo atom
	
	const size_t curcount=size();
	if (!empty() && (cellsize_==0) && findunitcell && (strcmp(front().type,single_entry.type)==0)) //starting to repeat entries?
	  cellsize_=curcount;
	
	if (verbose_>1) {
	  if (cellsize_ && ((curcount % cellsize_)==0))
	    std::cout << "--- cell " << (curcount / cellsize_) << " ---\n";
	  std::cout << single_entry.serial << ": " <<  single_entry.type;
	  if (single_entry.residue)
	    std::cout << '(' << single_entry.residue << ')';
	  std::cout << " at " << single_entry.coords << '\n';
	}
	
	push_back(single_entry);
	continue;
      }
      if (!strcmp(buffer,"CONECT")) {
	atomlist.create(0); //reset
	push_ifnonnull(atomlist,buffer+6);
	push_ifnonnull(atomlist,buffer+11);
	push_ifnonnull(atomlist,buffer+16);
	push_ifnonnull(atomlist,buffer+21);
	push_ifnonnull(atomlist,buffer+26);
	addConnections(atomlist);
	if (verbose_>1)
	  std::cout<< "Connected atoms: " << atomlist << '\n';
	continue;
      }
    }
    if (cellsize_==0) {
      cellsize_=size();
      ncells_=1;
    }
    else {
      if (size() % cellsize_) {
	char buf[256];
	snprintf(buf,sizeof(buf)," (%lu atoms do not divide into cell of %lu atoms)",(long unsigned)(size()),(long unsigned)cellsize_);
	division_warning.raise(buf);
      }
      else
	ncells_ = size() / cellsize_;
    }
    tidy_connections();
    if (has_connections() && verbose_) {
      std::cout << connections.size() << " molecules found\n";
      if (verbose_>1)
	std::cout<< connections<< '\n';
    }
  }

  const vector3& MoleculeStructure::axis(char which) const
  {
    if (!has_axes())
      throw Failed("MoleculeStructure::axis: axis information missing");
    switch (which) {
    case 'o':
      return axisdef[0];
    case 'a':
      return axisdef[1];
    case 'b':
      return axisdef[2];
    case 'c':
      return axisdef[3];
    }
    throw InvalidParameter("MoleculeStructure::axis: axis must be one of o,a,b,c");
  }
      

}//namespace libcmatrix

//old PDB read routines (nicked from somewhere)


// static long ReadPDBCoord(int offset)
// {
//     int len,neg;
//     long result;
//     char *ptr;
//     char ch;

//     result = 0;
//     neg = 0;
//     len = 8;

//     ptr = Record+offset;
//     while( len-- )
//     {   ch = *ptr++;
// 	if( (ch>='0') && (ch<='9') )
// 	{   result = (10*result)+(ch-'0');
// 	} else if( ch=='-' )
// 	    neg = 1;
//     }

//     /* Handle Chem3D PDB Files! */
//     if( Record[offset+3]=='.' )
// 	result /= 10;
//     return( neg? -result : result );
// }

//   void MoleculeStructure::read_PDB(FILE* fp, char selatom)
//   {
//     int count=0;

//     long dx,dy,dz;
//     size_t resid;
    
//     do {
//       FetchRecord(fp);
      
//       if (!strncmp("ATOM",Record+1,4)) {
// 	/* Ignore Pseudo Atoms!! */
// 	if ( (Record[13]==' ') && (Record[14]=='Q') )
// 	  continue; 
	
// 	dx = ReadPDBCoord(31);
// 	dy = ReadPDBCoord(39);
// 	dz = ReadPDBCoord(47);
	
// 	/* Ignore XPLOR Pseudo Atoms!! */
// 	if ( (dx==9999000L) && (dy==9999000L) && (dz==9999000L) )
// 	  continue;
	
// 	char type=toupper(Record[14]);

// 	count++;
// 	if (selatom && type!=selatom)
// 	  continue;
	
// 	sscanf(Record+22,"%i",&resid);
// 	push_back(PDBatom(vector3(dx/1000.0,dy/1000.0,dz/1000.0),type,count,resid));
//       }
//       else
// 	if(!strncmp("ENDM",Record+1,4))
// 	  break;
      
//     } while (!feof(fp));
//   }
