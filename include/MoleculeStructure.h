#ifndef lcm_moleculestructure_h_
#define lcm_moleculestructure_h_

#include "geometry.h"
#include "List.h"
#include "Warnings.h"
#include <string>

namespace libcmatrix {

  struct PDBatom {
    vector3 coords;
    char type[8];
    size_t serial;
    size_t residue;

    PDBatom() { *type='\0'; }
    PDBatom(const vector3& coords_, char type_, size_t serial_ =0, size_t residue_ =0) 
      : coords(coords_), serial(serial_), residue(residue_) {
      type[0]=type_;
      type[1]='\0';
    }

    PDBatom(const vector3& coords_, const char* type_, size_t serial_ =0, size_t residue_ =0) 
      : coords(coords_), serial(serial_), residue(residue_) {
      strncpy(type,type_,8);
    }

    bool istype(const char*) const;
  };
  
  struct space_group_info_t {
    space_group_info_t()
      : Z(0) {}

    space_group_info_t(const char* str, size_t Zv)
      : sgroup(str), Z(Zv) {}

    std::string sgroup;
    size_t Z;
  }; 

  struct lattice_info_t {
    lattice_info_t()
      : a(0.0), b(0.0), c(0.0), alpha(0.0), beta(0.0), gamma(0.0) {}

    lattice_info_t(double av, double bv, double cv, double alphav, double betav, double gammav)
      : a(av), b(bv), c(cv), alpha(alphav), beta(betav), gamma(gammav) {
      if (a*b*c<=0.0)
	throw InvalidParameter("lattice_info: invalid lattice parameters");
    }

    double a, b, c;
    double alpha;
    double beta;
    double gamma;
  };
    
  class MoleculeStructure : public List<PDBatom> {
  public:
    
    
    MoleculeStructure(int verbosev =0) : verbose_(verbosev), haveaxes_(false),
					 cellsize_(0), ncells_(0), cell0_(-1) {}

    void read_molecule(const char*, char ='\0');
    void read_molecule(FILE*, char ='\0');
    void read_PDB(FILE*, char ='\0', bool findunitcell =true);
    void read_PDB(const char*, char ='\0', bool findunitcell =true);
    void read_PDB_single(FILE*, int, char ='\0', bool findunitcell =true); //!< added 2017-11-07 to allow reading single model
    void read_PDB_single(const char*, int, char ='\0', bool findunitcell =true);

    int snumber_to_index(size_t n) const { return serial_to_internal_(n); }//return internal index of atom with supplied serial number
    List<size_t> find(const char*, bool incell0 =false) const; //return internal indexes of atoms of given type (optionally restricted to "middle" cell)

    size_t cell0_start() const { return cell0_*cellsize_; }
    size_t cell0_end() const { return (cell0_+1)*cellsize_; }
    const vector3& centre() const { return centre_; }

    int molecule(size_t m) const {//return internal index of molecule, which contain atom with given internal index (-1 if not included in a molecule)
      return has_connections() ? whichmolecule_(m) : -1;
    }
    
    bool has_connections() const { return !(connections.empty()); }
    bool has_axes() const { return haveaxes_; }
    const vector3& axis(char) const;
    bool has_lattice() const { return (lattice_info_.a != 0.0); }
    const lattice_info_t& lattice_parameters() const {
      if (has_lattice())
	return lattice_info_;
      else
	throw Failed("MoleculeStructure: no lattice information found");
    }
    const space_group_info_t& space_group() const {
      if (has_lattice())
	return space_group_info_;
      else
	throw Failed("MoleculeStructure: no lattice information found");
    }
    
    void clear() {
      List<PDBatom>::clear();
      haveaxes_=false;
    }

    static Warning<> connectivity_warning;
    static Warning<> missing_axis_warning;
    static Warning<> division_warning;

  private:
    int verbose_;
    bool haveaxes_;
    vector3 axisdef[4];
    lattice_info_t lattice_info_;
    space_group_info_t space_group_info_;

    List< List<size_t> > connections;//contain serial numbers, grouped by molecules
    List<int> whichmolecule_;
    List<int> serial_to_internal_;
    vector3 centre_;
    size_t cellsize_,ncells_,cell0_;

    void addConnections(const BaseList<size_t>&);
    void tidy_connections();
    void read_PDB_(FILE*, int slice, char ='\0', bool findunitcell =true);
  };

}//namespace libcmatrix

#endif
