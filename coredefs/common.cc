/* Shared bits and pieces */

#undef LCM_SUPPRESS_VIEWS
// flag that Level 1 (list) data types should be instantiated
#define LCM_LEVEL1_INSTANT
#include "List.h"
#include "IndirectList.h"
#include "CachedAllocator.h"
#include "MultiMatrix.h"
#include <sstream>
#include "Warnings.h"
#include "simple_counter.h"

namespace libcmatrix {
  
  const char cmatrix_abi_version[]="3.15.0";

  FixedAllocator::FixedAllocator(size_t sizev, size_t totalv)
    : size_(sizev), stack_(totalv/sizev)
  {
    stack_.resize(0);
    failedcount=successcount=0;
  }

#ifdef LCM_NEED_128BIT_ALIGN
#define LCM_ALIGNMENT 16
#else
#define LCM_ALIGNMENT 0
#endif

  static DefaultAllocator<LCM_ALIGNMENT> rawfixedallocator;

  FixedAllocator::~FixedAllocator() {
    while (!(stack_.empty())) {
      //      operator delete(stack_.back());
      rawfixedallocator.deallocate(stack_.back(),size_);
      stack_.pop_back();
    }
  }
  
  void* FixedAllocator::allocate()
  {
    if (!(stack_.empty())) {
      void* newp=stack_.back();
      stack_.pop_back();
      successcount++;
      return newp;
    }
    failedcount++;
    //    return operator new(size_);
    return rawfixedallocator.allocate(size_);
  }

  void FixedAllocator::deallocate(void* p)
  {
    if (stack_.size()!=stack_.capacity())
      stack_.push_back(p);
    else
      rawfixedallocator.deallocate(p,size_);
      //operator delete(p);
  }
  
  void FixedAllocator::print(std::ostream& ostr) const 
  {
    ostr << "Allocator for blocks of size " << size_ << " bytes: current no.=" << stack_.size();
    ostr << "  success=" << successcount << "  failed=" << failedcount << "\n";
  }
  
  CachedAllocator_::CachedAllocator_(int maxkv, int totalkv)
    : pLastAlloc_(NULL), pLastDealloc_(NULL), 
      destroyed_(false) 
  {
      setcachelimits(maxkv,totalkv);
  }

  CachedAllocator_::~CachedAllocator_()
  { //bit dirty
    const pool_t::iterator end=pool_.end();
    pool_t::iterator cur=pool_.begin();
    destroyed_=true;
    while (cur!=end) {
      FixedAllocator*& curv=*cur;
      delete curv;
      curv=NULL;
      cur++;
    }
  }

  void CachedAllocator_::setcachelimits(int maxkv, int totalkv)
  {
    if ((maxkv>totalkv) || (totalkv<0))
      throw InvalidParameter("CachedAllocator: bad cache limits");
    max_=1024*maxkv;
    total_=1024*totalkv;
  }
  
  void* CachedAllocator_::allocate(size_t length_)
  {
    if (length_>max_) {
      oversize_claimed_words_+=(length_>>2);
      //      return operator new(length_);
      return rawfixedallocator.allocate(length_); // 18/2/2016 - replaced operator new
    }
    
    if (!pLastAlloc_ || (pLastAlloc_->BlockSize()!=length_)) {
      pool_t::iterator iter = std::lower_bound(pool_.begin(), pool_.end(), length_, 
					       CompareFixedAllocatorSize());
      if (iter==pool_.end() || (*iter)->BlockSize()!=length_) {
	iter=pool_.insert(iter, new FixedAllocator(length_,total_));
	pLastDealloc_ = *(pool_.begin()); //pLastDealloc may be invalid - reset to something valid
      }
      pLastAlloc_ = *iter;
    }
    return pLastAlloc_->allocate();
  }

  Warning<> CachedAllocator_::post_destruction_warning("deallocate called after CachedAllocator destruction",&lcm_base_warning,LCM_DEBUG_WARNING);
  
  void CachedAllocator_::deallocate(void* p, size_t length_) 
  {
    if (destroyed_) {
      //     post_destruction_warning.raise();
      return; //Don't bother to release empty memory
    }
    
    if (length_>max_) { //Not critical if max_ has been changed
      oversize_released_words_+=(length_>>2);
      rawfixedallocator.deallocate(p,length_); // 18/2/2016 - replaced operator delete
      //      operator delete(p);
      return;
    }
    if (!pLastDealloc_ || (pLastDealloc_->BlockSize()!=length_)) {
      pool_t::iterator iter = std::lower_bound(pool_.begin(), pool_.end(), length_, 
					       CompareFixedAllocatorSize());
      assert(iter != pool_.end());
      assert((*iter)->BlockSize() == length_);
      pLastDealloc_=*iter;
    }
    return pLastDealloc_->deallocate(p);
  }
  
  void CachedAllocator_::print(std::ostream& ostr) const
  {
    const pool_t::const_iterator end=pool_.end();
    pool_t::const_iterator cur=pool_.begin();
    while (cur!=end) {
      (*cur)->print(ostr);
      ++cur;
    }
    ostr << "Oversize words claimed: " << oversize_claimed_words_ << "  released: " << oversize_released_words_ << "  diff: " << (oversize_claimed_words_-oversize_released_words_) << '\n';
  }

  size_t CachedAllocator_::allocated_words() const
  {
    size_t alloc=0;
    const pool_t::const_iterator end=pool_.end();
    pool_t::const_iterator cur=pool_.begin();
    while (cur!=end) {
      const FixedAllocator& curalloc(*(*cur));
      alloc+=((curalloc.BlockSize()) >> 2)*curalloc.BlockCount();
      ++cur;
    }
    return alloc+oversize_claimed_words_-oversize_released_words_;
  }

  std::ostream& operator << (std::ostream& ostr, const MatrixException& obj)
  {
    ostr << obj.type;
    if (!obj.errm.empty())
      ostr << ": " << obj.errm;
    return ostr << std::endl;
  }

  namespace {
    struct MakeError {
      MakeError(const char* errm_)
	: ostr(std::ostringstream::out), haveerr(errm_ && errm_[0]) {
	if (haveerr)
	  ostr << errm_ << " (";
      }
      
      void flush(std::string& dest) {
	if (haveerr)
	  ostr << ')';
	dest=ostr.str();
      }
      mutable std::ostringstream ostr;
      bool haveerr;
    };
    
    template<typename T> MakeError& operator<< (MakeError& ostr, const T& a) {
      ostr.ostr << a;
      return ostr;
    }
  }
	
  BadIndex::BadIndex(const char* errm_, int arg_, size_t max_)
    : MatrixException("Bad index"),
      std::out_of_range(errm_) 
  {
    MakeError ostr(errm_);
    const int indexbase=cmatrix_ostream_controller(std::cerr).indexbase;
    ostr << "requested element " << (arg_+indexbase) << " when object ";
    if (max_==0)
      ostr << "is empty";
    else {
      ostr << "has only " << max_ << " element";
      if (max_>1)
	ostr << 's';
    }
    ostr.flush(errm);
  }

  BadIndex::BadIndex(const char* errm_, int arg1_, size_t max1_, int arg2_, size_t max2_)
    : MatrixException("Bad index"),
      std::out_of_range(errm_) 
  {
    MakeError ostr(errm_);
    const int indexbase=cmatrix_ostream_controller(std::cerr).indexbase;
    ostr << "requested element " << (indexbase+arg1_) << ',' << (indexbase+arg2_) << " when ";
    if (max1_ && max2_)
      ostr << "matrix is " << max1_ << 'x' << max2_;
    else
      ostr << "matrix is undefined";
    ostr.flush(errm);
  }

  BadIndex::BadIndex(const char* errm_, int arg1_, size_t max1_, int arg2_, size_t max2_, int arg3_, size_t max3_)
    : MatrixException("Bad index"),
      std::out_of_range(errm_) 
  {
    MakeError ostr(errm_);
    const int indexbase=cmatrix_ostream_controller(std::cerr).indexbase;
    ostr << "requested element " << (arg1_+indexbase) << ',' << (indexbase+arg2_) << ',' << (indexbase+arg3_) << " when ";
    if (max1_ && max2_ && max3_)
      ostr << "matrix is " << arg1_ << 'x' << arg2_ << 'x' << arg3_;
    else
      ostr << "matrix is undefined";
    ostr.flush(errm);
  }

  Mismatch::Mismatch(const char* errm_, size_t v1, size_t v2)
    : MatrixException("Mismatch"),
      std::domain_error(errm_)
  {
    MakeError ostr(errm_);
    ostr << v1 << " vs. " << v2;
    ostr.flush(errm);
  }
      
  Mismatch::Mismatch(const char* errm_, size_t r1, size_t c1, size_t r2, size_t c2)
    : MatrixException("Mismatch"),
      std::domain_error(errm_)
  {
    MakeError ostr(errm_);
    ostr << r1 << 'x' << c1 << " matrix vs. " << r2 << 'x' << c2;
    ostr.flush(errm);
  }
      
  //assumes size_t is unsigned
bool isvalid_indexlist(const BaseList<size_t>& a,size_t maxv)
{
  for (size_t n=a.length();n--;) {
    if (a(n)>=maxv)
      return false;
  }
  return true;
}

//   void lcm_verify_init_(int& iac, int& ibc, size_t ac, size_t bc)
//   {
//     if (iac<0) {
//       iac=a.cols();
//       ibc=b.cols();
//     }
//     else {
//       if ((iac<a.cols()) || (ibc<b.cols()))
// 	throw InvalidParameter("multiply: storage steps smaller than number of columns");
//     }
//   }

  void slice::create(int startv, int sizev, int stridev)
  { 
    start_=startv;
    size_=sizev;
    stride_=stridev;
    if (stride_==0)
      throw InvalidParameter("slice(): stride cannot be zero");
    if (size_) {
      const int endv=startv+(sizev-1)*stridev;
      max_=(start_>endv) ? start_ : endv;
      if (startv<0 || sizev<0 || endv<0)
	throw InvalidParameter("slice(): indices must be non-negative");
    }
    else
      max_=start_;      
  }

slice& slice::operator*= (int scale)
{
  if (scale<=0)
    throw InvalidParameter("slice: scaling factor must be +ve");
  start_*=scale; stride_*=scale; max_*=scale;
  return *this;
}

slice& slice::operator+= (int offset) 
{
  if (offset<0)
    throw InvalidParameter("slice: offset cannot be -ve");
  start_+=offset; max_+=offset;
  return *this;
}
      
// static const int MAX_ERROR=8;

// static const char* error_names[MAX_ERROR]={"<No error>",
// 				  "Operation failed.",
// 				  "Write failed.",
// 				  "File open failed.",
// 				  "Corrupt file?",
// 				  "Unknown file type.",
// 				  "Invalid parameter for operation.",
// 				  "Fatal internal error."
// };

// const char* error_name(int err)
// {
//   const char* result= (err>=0 && err<MAX_ERROR) ? error_names[err] : NULL;
//   if (!result)
//     throw InvalidParameter("error_name");
//   return result;
// }

MultiIterator::MultiIterator(const BaseList<size_t>& dims)
  : dims_(dims), mults_(dims.size()), cur_(dims.size())
{
  //  if (std::find(dims.begin(),dims.end(),0)!=dims.end())
  total_=1;
  if (dims_.size()) {
    //    throw InvalidParameter("MultiIterator: zero dimensions");
    for (size_t i=dims.size();i--;) {
      mults_(i)=total_;
      total_*=dims_(i);
    }
    if (total_==0)
      throw InvalidParameter("MultiIterator: dimensions must be non-zero");
  }
  reset();
}

void MultiIterator::reset()
{
  if (!dims_.empty())
    cur_=size_t(0);
  finished_=false;
}

  void MultiIterator::reverse(BaseList<size_t> dest, size_t ind) const
  {
    if (ind>=total_)
      throw Mismatch("MultiIterator: index exceeds total number of steps",ind,total_);
    const size_t n=dims_.size();
//     if (n<1)
//       throw Undefined("MultiIterator::reverse");
    if (dest.size()!=n)
      throw Mismatch("MultiIterator: destination vector vs. actual",dest.size(),n);
    if (n) {
      size_t i=0;
      for (;i<n-1;i++) {
	const size_t cm=mults_(i);
	const size_t ci=ind / cm;
	dest(i)=ci;
	ind-=ci*cm;
      }
      dest(i)=ind;
    }
  }

bool MultiIterator::next(BaseList<size_t> dest)
{
  if (finished_)
    return false;
  dest=cur_;
  if (dest.empty()) 
    finished_=true;
  else {
    size_t cdim=0;
    while (++cur_(cdim)==dims_(cdim)) {
      cur_[cdim]=0;
      if (++cdim==size()) {
	finished_=true;
	break;
      }
    }
  }
  return true;
}    

BaseWarning lcm_base_warning(BaseWarning::Always); //!< default is to show all warnings
BaseWarning lcm_io_warning(lcm_base_warning);
BaseWarning lcm_serious_warning(BaseWarning::RaiseException); //!< default to raise exception

Warning<> Mutex<true>::mutexnotfunctional_warning("no mutexes as no functional libpthread",&lcm_serious_warning);

void BaseWarning::print() const
{
  ostr_ << "'" << name_ << "' warning ";
  if (count_) {
    ostr_ << "raised " << count_ << " time";
    if (count_>1)
	ostr_ << 's';
    ostr_ << '\n';
  }
  else
    ostr_ << "never raised\n";
}

BaseWarning::BaseWarning(warning_t typev)
  : parentp_(NULL), type_(typev),
    ostr_(std::cerr) //!< needed so reference is initialised but should never be used
{
  type(typev);
  reset();
}

BaseWarning::BaseWarning(BaseWarning& parentv, warning_t typev) 
  : parentp_(&parentv), type_(typev), //!< must have parent so no need to check
    ostr_(std::cerr) //!< needed so reference is initialised but should never be used
{
  reset();
}

BaseWarning::BaseWarning(const char* namev, BaseWarning* parentpv, warning_t typev, std::ostream& ostrv) 
  : parentp_(parentpv), name_(namev),
    ostr_(ostrv)
{
  type(typev); 
  reset();
}

void BaseWarning::reset()
{
  count_=0;
  stringsmap_.clear();
}

BaseWarning::warning_t BaseWarning::evaluated_type() const
{
  if (type_==Inherit) {
    assert(parentp_!=NULL);
    return parentp_->evaluated_type();
  }
  return type_;
}

void BaseWarning::type(warning_t typev)
{
  if ((typev==Inherit) && !parentp_)
    throw InvalidParameter("BaseWarning::type - no parent to inherit from");
  type_=typev;
}
  
void BaseWarning::recursive_increment()
{
  count_++;
  if (parentp_)
    parentp_->recursive_increment();
}

std::ostream& BaseWarning::print_message(const char* extra) const
{  
  ostr_ << "Warning: " << name_;
  if (extra)
    ostr_ << extra;
  return ostr_;
}

void BaseWarning::raise_generic(warning_t typev, const char* extra, bool uniquev)
{
  recursive_increment();
  switch (typev) {
  case Silent: return;
  case Always:
    print_message(extra) << '\n';
    return;
  case FirstOnly:
    if (uniquev) {
      const std::string sextra(extra);
      if (stringsmap_.count(sextra)==0) {
	print_message(extra) << " (future identical warnings suppressed)\n";
	stringsmap_.insert(sextra);
      }
    }
    else {
      if (count_==1)
	print_message(extra) << " (future warnings of this type suppressed)\n";
    }
    return;
  default:
    throw InternalError("BaseWarning::raise_generic");
  }      
}

  Warning<InternalError> simple_counter::negative_warning("counter has gone negative!",&lcm_base_warning);

  simple_counter::simple_counter(long limitv, long usedv)
    : used_(usedv)
  {
    limit(limitv);
    if (usedv<0)
      throw InvalidParameter("simple_counter");
  }

  void simple_counter::limit(long valv)
  {
    if (valv<0)
      throw InvalidParameter("simple_counter: limit");
    limit_=valv;
  }

  simple_counter& simple_counter::operator+=(long valv)
  {
    used_+=valv;
    if (used_<0)
      negative_warning.raise();
    return *this;
  }

  Warning<> dynamiclist_alias_warning("aliasing DynamicList",&lcm_base_warning,LCM_DEBUG_WARNING);

void base_parallel_controller::writelog(const char* mess) const
{
  if (verbose_ || (log!=NULL)) {
    FILE* fp= log ? log : stdout;
    fprintf(fp,"Process %i: %s (at time +%g minutes)\n",id_,mess,timer_()/60.0);
    fflush(fp);
  }
}

  Warning<> base_parallel_controller::logopenfailed_warning("Failed to open log file for parallel controller. Logging output will be discarded or written to standard out depending on verbose level",&lcm_base_warning);

  void base_parallel_controller::closelog()
  {
    if (log)
      fclose(log);
    log=NULL;
  }

  void base_parallel_controller::openlog(const char* name)
  {
    log=fopen(name,"w");
    if (!log)
      logopenfailed_warning.raise(name);
  }

  base_parallel_controller::~base_parallel_controller()
  {
    if (log)
      fclose(log);
  }

}//namespace libcmatrix
