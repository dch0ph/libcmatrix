#ifndef cmatrix_hash_map_h_
#define cmatrix_hash_map_h_

//Quick and dirty string hash map
//This should not be used in new code

#ifdef _MSC_VER
#include <hash_map>
#define LCM_HASHER stdext::hash_value
#else
#include <ext/hash_map>
#define LCM_HASHER hasher_
#endif
#include <map>

namespace libcmatrix {
  struct stringhash : public std::unary_function<const char*,size_t> {
#ifndef _MSC_VER
    __gnu_cxx::hash<const char*> hasher_;
#endif
    size_t operator()(const std::string& a) const { return LCM_HASHER(a.c_str()); }
    size_t operator()(const char* a) const { return LCM_HASHER(a); }
  };

//NB this isn't a fully encapsulated version of map

  template<class Hasher =stringhash, class StoreKeyType =std::string> struct hasher_key_t {
    typedef typename Hasher::argument_type arg_t;    
    hasher_key_t(const Hasher& hasher, const arg_t& v) : hashed_key(hasher(v)), name(v) {}
    size_t hashed_key;
    const StoreKeyType name;

    bool operator== (const hasher_key_t& a) const { return (hashed_key==a.hashed_key); }
    bool operator!= (const hasher_key_t& a) const { return (hashed_key!=a.hashed_key); }
    bool operator< (const hasher_key_t& a) const { return (hashed_key<a.hashed_key); }
    bool operator> (const hasher_key_t& a) const { return (hashed_key>a.hashed_key); }
    bool operator<= (const hasher_key_t& a) const { return (hashed_key<=a.hashed_key); }
    bool operator>= (const hasher_key_t& a) const { return (hashed_key>=a.hashed_key); }
  };

  template<class StoreType, class Hasher =stringhash, class StoreKeyType =std::string> class hashed_map_key : public std::map< hasher_key_t<Hasher>, StoreType> {
    Hasher hasher_;
    typedef typename Hasher::argument_type arg_t;
    typedef hasher_key_t<Hasher,StoreKeyType> key_t;
    typedef std::map<key_t,StoreType> map_t;

  public:
    StoreType& operator[](const arg_t& name) { return map_t::operator[](key_t(hasher_,name)); }
    typename map_t::const_iterator find(const arg_t& name) const { return map_t::find(key_t(hasher_,name)); }
    typename map_t::iterator find(const arg_t& name) { return map_t::find(key_t(hasher_,name)); }
    size_t count(const arg_t& name) const { return map_t::count(key_t(hasher_,name)); }
  };

  template<class StoreType,class Hasher =stringhash> class hashed_map : public std::map<size_t, StoreType> {
    Hasher hasher_;
    typedef typename Hasher::argument_type arg_t;
  public:
    typedef std::map<size_t,StoreType> map_t;
    //    typedef StoreType datatype;
    StoreType& operator[](const arg_t& name) { return map_t::operator[](hasher_(name)); }
    typename map_t::const_iterator find(const arg_t& name) const { return map_t::find(hasher_(name)); }
    typename map_t::iterator find(const arg_t& name) { return map_t::find(hasher_(name)); }
    size_t count(const arg_t& name) const { return map_t::count(hasher_(name)); }
    size_t hash(const arg_t& name) const { return hasher_(name); }
  };
  
}

#endif
