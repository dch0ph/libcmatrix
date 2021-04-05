/*! \file
  \brief write (read?) of "composite" Matlab V5 formats (cell and struct) */

#include "matlabio.h"
#include "cmatrix_utils.h"
#include <cstring>

namespace libcmatrix {

  Matlab5WriteStoreBase::Matlab5WriteStoreBase(const char* namev) { 
    ::std::strncpy(name,namev,matlab_controller_maxstructlength); 
  }
  
   void matlab_controller::composite::ensureknown() const {
     if ((max_==0) && (class_==mxSTRUCT_CLASS))
       throw Failed("matlab_controller::composite: cannot create composite object inside struct of unknown size");
    }

    void matlab_controller::composite::ensurewrite() const {
      ensureactive();
      ctrl_.ensurewrite();
    }

    void matlab_controller::composite::ensureread() const {
      ensureactive();
      ctrl_.ensureread();
    }

    void matlab_controller::composite::ensureactive() const {
      if (!active_)
	throw Failed("matlab_controller::composite: composite object has been closed");
    }

  Warning<> matlab_controller::composite::ignoring_item_name_warning("matlab_controller::composite: item name is ignored for cell save",&lcm_io_warning);
  Warning<> matlab_controller::composite::write_overrun_warning("matlab_controller::composite: more writes than expected",&lcm_io_warning);

    void matlab_controller::composite::checkwritestart(const char* name) const
    {
      ensurewrite();
      const bool havename=name && name[0];
      if (havename) {
	if (!disablecheck_ && (class_==mxCELL_CLASS))
	  ignoring_item_name_warning.raise();
      }
      else {
	if (class_==mxSTRUCT_CLASS)
	  throw InvalidParameter("matlab_controller::composite: item name must be given for struct save");
	if (strlen(name)>=matlab_controller_maxstructlength)
	  throw InvalidParameter("matlab_controller::composite: struct names must be less 32 characters");
      }
      if (max_) {
	if (class_==mxSTRUCT_CLASS) {
	  file_position_lock lock(ctrl_.fp);
	  ctrl_.fseek(namespos_+count_*matlab_controller_maxstructlength);
	  fputs(name,ctrl_.file_pointer()); //!< write null-terminated string
	  fputc('\0',ctrl_.file_pointer());
	}
	if (count_>=max_) {
	  if (disablecheck_)
	    write_overrun_warning.raise();
	  else
	    write_overrun_warning.raiseas(BaseWarning::RaiseException);
	}
      }
      checkposition();
    }

  Warning<> matlab_controller::composite::intruding_write_warning("matlab_controller::composite: intruding write?  File is likely to be corrupted",&lcm_io_warning);

void matlab_controller::composite::checkposition() const
{
  if (!disablecheck_ && (ctrl_.ftell()!=checkpos_))
    intruding_write_warning.raise();
}

void matlab_controller::composite::validate(matlab_t type)
{
  switch (type) {
  case STRUCT:
    class_=mxSTRUCT_CLASS;
    break;
  case CELL:
    class_=mxCELL_CLASS;
    break;
  default:
    throw InvalidParameter("matlab_controller::composite: type must be CELL or STRUCT");
  }  
  count_=0;
}

matlab_controller::composite::composite(matlab_controller& ctrl, const char* name, matlab_t type, bool nocheck)
  : ctrl_(ctrl), disablecheck_(nocheck), max_(0)
{
  write_header(type,ScratchList<size_t>(0U,0U),name);
}

matlab_controller::composite::composite(composite& comp, const char* name, matlab_t type, bool nocheck)
  : ctrl_(comp.ctrl_), disablecheck_(nocheck), max_(0)
{
  comp.ensureknown();
  comp.checkwritestart(name); //!< write name
  validate(type);
  write_header(type,ScratchList<size_t>(0U,0U),(class_==mxSTRUCT_CLASS) ? "" : name);
  comp.increment();//!< increment object count in parent
}

void matlab_controller::composite::write_header(matlab_t type, const BaseList<size_t>& dims, const char* name)
{
  ctrl_.ensurewrite();
  validate(type);
  FILE* fp=ctrl_.file_pointer();
  if (class_==mxSTRUCT_CLASS) {
    initpos_=WriteMATLAB5MatrixHeader(fp,mxSTRUCT_CLASS,false,ScratchList<size_t>(size_t(1),size_t(1)),name);
    static const int32_t maxstructlength=matlab_controller_maxstructlength;
    WriteMATLAB5Tag(fp,miINT32,4,(void*)(&maxstructlength),true);
    namespos_=ctrl_.ftell()+8;
    if (max_) {
      const UINT32_t tagsize=WriteMATLAB5Tag_raw(fp,miINT8,maxstructlength*max_,false);
      assert(tagsize==8U);
      ctrl_.fseek(namespos_+max_*maxstructlength); //!< add tag size + space required by names
    }
  }
  else
    initpos_=WriteMATLAB5MatrixHeader(fp,class_,false,dims,name);
  checkupdate();
  active_=true;
}

matlab_controller::composite::composite(matlab_controller& ctrl, const char* name, matlab_t type, int entries, bool nocheck)
  : ctrl_(ctrl), disablecheck_(nocheck), max_(entries)
{
  if (entries<1)
    throw InvalidParameter("matlab_controller::composite");  
  write_header(type,ScratchList<size_t>(size_t(1),size_t(entries)),name); 
}

matlab_controller::composite::composite(composite& comp, const char* name, matlab_t type, int entries, bool nocheck)
  : ctrl_(comp.ctrl_), disablecheck_(nocheck), max_(entries)
{
  comp.ensureknown();
  if (entries<1)
    throw InvalidParameter("matlab_controller::composite");  
  comp.checkwritestart(name); //!< write name
  validate(type);
  write_header(type,ScratchList<size_t>(size_t(1),size_t(entries)),(class_==mxSTRUCT_CLASS) ? "" : name); 
  comp.increment();
}

  namespace {
    template<class T> T product_(const BaseList<T>& a) {
      T prod(1);
      for (size_t i=a.size();i--;)
	prod*=a(i);
      return prod;
    }
  }

matlab_controller::composite::composite(matlab_controller& ctrl, const char* name, matlab_t type, const BaseList<size_t>& dims, bool nocheck)
  : ctrl_(ctrl), disablecheck_(nocheck)
{
  max_=product_(dims);
  if (max_==0)
    throw InvalidParameter("matlab_controller::composite");
  write_header(type,dims,name);
}

matlab_controller::composite::composite(composite& comp, const char* name, matlab_t type, const BaseList<size_t>& dims, bool nocheck)
  : ctrl_(comp.ctrl_), disablecheck_(nocheck)
{
  comp.ensureknown();
  max_=product_(dims);
  if (max_==0)
    throw InvalidParameter("matlab_controller::composite");
  write_header(type,dims,name);
  comp.increment();
}

  void matlab_controller::composite::increment()
  {
    count_++;
    checkupdate();
  }

matlab_controller::composite::composite(matlab_controller& ctrl)
  : ctrl_(ctrl)
{
  read_init();
}

  Warning<> matlab_controller::composite::name_length_warning("matlab_controller::composite: non-standard name length detected",&lcm_io_warning);

  void matlab_controller::composite::read_init()
  {
    count_=0;
    namelen_=0;
    header_info info;
    FILE* fp(ctrl_.fp);
    initpos_=::std::ftell(fp);
    //std::cout << "Current file pos: " << initpos_ << "\n";
    const long endpos=ctrl_.read_size(); //!< not sure what to do with this!
    if (endpos==0)
      throw InternalError("matlab_controller::composite::read_init");
    ctrl_.read_matrix_header(info);
    class_=static_cast<matlab_class_t>(info.dataclass);
    if (!ctrl_.next()) //!< skip array name tag
      throw Failed("matlab_controller::composite: name missing");
    //    Matlab5Tag arrayname;
    //arrayname.skip(fp,ctrl_.isbigend);
    switch (class_) {
    case mxCELL_CLASS:
      max_=product_(info.dims);
      break;
    case mxSTRUCT_CLASS: {
      class_=mxSTRUCT_CLASS;
      Matlab5Tag namelen(fp,ctrl_.isbigend);
      const size_t len=namelen.row<int32_t>().front();
      if (len!=matlab_controller_maxstructlength)
	name_length_warning.raise();
      Matlab5Tag names;
      names.read_header(fp,ctrl_.isbigend);    
      assert(!(names.iscompressed));
      namespos_=ctrl_.ftell();
      assert((names.nbytes % len)==0);
      max_=names.nbytes / len;
      ctrl_.fseek(namespos_+max_*len);
    }
      break;
    default:
      throw Failed("matlab_controller::composite: item is not a composite object");
    }
    active_=true;
  }

  matlab_controller::composite::composite(composite& comp)
    : ctrl_(comp.ctrl_)
  {
    read_init();
    comp.count_++; //!< increment count in parent
  }

 bool matlab_controller::composite::peek(header_info& info)
{
  if (count_==max_)
    return false;    
  if (!ctrl_.peek(info))
      throw InternalError("matlab_controller::composite::peek: corrupt file?");
  if (class_==mxSTRUCT_CLASS) {
    file_position_lock lock(ctrl_.fp);
    ctrl_.fseek(namespos_+count_*matlab_controller_maxstructlength);
    fread(namebuf,matlab_controller_maxstructlength,1,ctrl_.fp);
    info.name=namebuf;
  }
  return true;
}

bool matlab_controller::composite::next() 
{
  if (count_==max_)
    return false;
  if (!ctrl_.next())
    throw InternalError("matlab_controller::composite::next");
  count_++;
  return true;
}

bool matlab_controller::composite::find(header_info& info, const char* name)
{
  file_position_lock lock(ctrl_.fp); //!< restore original file pos unless match found
  while (peek(info)) {
    if (strcmp(name,info.name)==0) {
      lock.unlock();
      return true;
    }
    if (!next())
      return false;
  }
  return false;
}

  Warning<> matlab_controller::composite::write_mismatch_warning("matlab_controller::composite: items written doesn't match original declaration (file is likely to be corrupt)",&lcm_io_warning);

bool matlab_controller::composite::ok_to_close() const
{
  if (!active_ || ctrl_.isreading())
    return true;

  if (max_ || (class_!=mxSTRUCT_CLASS))
    return true;

  return (writestack_.size()>0);
}

void matlab_controller::composite::close()
{
  if (!active_)
    return;

  FILE* fp(ctrl_.file_pointer());
  if (ctrl_.isreading()) {
    if (count_!=max_) {
      ctrl_.fseek(initpos_);
      ctrl_.next(); //!< advance to end of composite
    }
  }
  else {
    if (max_ || (class_!=mxSTRUCT_CLASS)) {
      if (!disablecheck_ && max_ && (max_!=count_))
	write_mismatch_warning.raise();
      TidyMATLAB5Header(fp,initpos_);
    }
    else {
      const size_t nitems=writestack_.size();
      if (nitems==0)
	throw Failed("matlab_controller::composite: nothing to write!"); //!< don't have much choice here - can't carry on since file will be junk
      checkposition();
      ctrl_.fseek(namespos_-8); //!< -8 to move to header start
      const UINT32_t tagsize=WriteMATLAB5Tag_raw(fp,miINT8,matlab_controller_maxstructlength*nitems,false);
      assert(tagsize==8U);
      for (size_t i=0;i<nitems;i++)
	fwrite(writestack_(i)->name,matlab_controller_maxstructlength,1,fp);
      for (size_t i=0;i<nitems;i++)
	(*writestack_(i))(ctrl_);
      TidyMATLAB5Header(fp,initpos_);
    }
  }
  active_=false;
}

} // namespace libcmatrix
