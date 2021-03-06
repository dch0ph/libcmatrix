#ifndef LCM_Warnings_h_
#define LCM_Warnings_h_

/*! \file
 \brief  Objects for handling/tracking warnings
*/
#include <set>
#include <string>
#include <cstdio>
#include "basedefs.h"

namespace libcmatrix {

  class BaseWarning {
  public:
    enum warning_t { Ignore, Silent, Always, Inherit, FirstOnly, RaiseException };
    explicit BaseWarning(warning_t =Always);
    explicit BaseWarning(BaseWarning&, warning_t =Inherit);
    explicit BaseWarning(const char* namev, BaseWarning* parentpv =NULL, warning_t typev =Inherit, std::ostream& ostrv =std::cerr); 
    virtual ~BaseWarning() {}
    size_t count() const { return count_; }
    void reset();
    bool enabled() const { return (evaluated_type()!=Ignore); }
    warning_t type() const { return type_; }
    warning_t evaluated_type() const;
    void type(warning_t);
    virtual std::ostream& print_message(const char* =NULL) const; //!< print raw message
    void print() const; //!< print status
    void flush() { //!< only output if called
      if (count_)
	print();
    }
  protected:
    BaseWarning* parentp_;
    std::string name_;
    warning_t type_;
    std::ostream& ostr_;
    size_t count_;

    typedef std::set<std::string> record_t;
    record_t stringsmap_;

    void raise_generic(warning_t, const char* =NULL, bool uniquev =false); //!< cases that don't depend on exception type
    void recursive_increment(); //!< recursively increment count
  };

  extern BaseWarning lcm_base_warning; //!< master control for libcmatrix warnings
  extern BaseWarning lcm_io_warning; //!< subset IO warnings
  extern BaseWarning lcm_serious_warning; //!< master control for "nasty" warnings

  template<class ExceptClass =Failed> class Warning : public BaseWarning {
  public:
    Warning(const char* namev, BaseWarning* parentpv, warning_t typev =Inherit, std::ostream& ostrv =std::cerr) 
      : BaseWarning(namev,parentpv,typev,ostrv) {}
    void raise(const char* extramv =NULL, bool uniquev =false) { raise_(evaluated_type(),extramv, uniquev); } //!< raise warning
  void raisesafe(const char* extramv =NULL); //!< raise warning, but without throwing exception
    void raiseas(warning_t typev, const char* extramv =NULL, bool uniquev =false) { raise_(typev==Inherit ? evaluated_type() : typev,extramv,uniquev); }
  
  private:
  void raise_(warning_t, const char*, bool);
  };

  template<class ExceptClass> void Warning<ExceptClass>::raisesafe(const char* extram)
    {
      warning_t typev=evaluated_type();
      if (typev==RaiseException) {
	std::fputs("[suppressing exception] ",stderr);
	typev=Always;
      }
      raise_(typev,extram,false);
    }

  template<class ExceptClass> void Warning<ExceptClass>::raise_(warning_t typev, const char* extram, bool uniquev)
  {
    switch (typev) {
    case Ignore:
      return;

    case RaiseException:
      throw ExceptClass(name_.c_str());

    default:
      raise_generic(typev,extram,uniquev);
    }
  }

#ifdef NDEBUG
#define LCM_DEBUG_WARNING BaseWarning::Ignore
#else
#define LCM_DEBUG_WARNING BaseWarning::Inherit
#endif

}//namespace libcmatrix

#endif
