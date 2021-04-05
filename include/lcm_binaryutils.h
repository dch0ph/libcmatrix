/* Routines for dealing with binary data */

#ifndef _binaryutils_h_
#define _binaryutils_h_

#include <sys/types.h>
#include "config.h"

#ifdef HAVE_STDINT_H
#include <stdint.h>
#else
    typedef signed __int64 int64_t;
    typedef signed __int32 int32_t;
    typedef signed __int16 int16_t;
    typedef signed __int8 int8_t;
    typedef unsigned __int64 uint64_t;
    typedef unsigned __int32 uint32_t;
    typedef unsigned __int16 uint16_t;
    typedef unsigned __int8 uint8_t;
#undef NEED_UNDERSCORE
#endif
#ifdef NEED_UNDERSCORE
#define UINT8_t  u_int8_t
#define UINT16_t u_int16_t
#define UINT32_t u_int32_t
#define UINT64_t u_int64_t
#else
#define UINT8_t  uint8_t
#define UINT16_t uint16_t
#define UINT32_t uint32_t
#define UINT64_t uint64_t
#endif

namespace libcmatrix {

inline bool ambigendian()
{
  static const int x=1;
  return ( (*(char *)&x) != 1);
}

 template<class Sel> void endianswap_(char*, Sel);
 template<> inline void endianswap_(char* in,Int2Type<2>) {
   ::std::swap(in[0],in[1]);
 }
 template<> inline void endianswap_(char* in,Int2Type<4>) {
   ::std::swap(in[0],in[3]);
   ::std::swap(in[1],in[2]);
 }
 template<> inline void endianswap_(char* in,Int2Type<8>) {
   ::std::swap(in[0],in[7]);
   ::std::swap(in[1],in[6]);
   ::std::swap(in[2],in[5]);
   ::std::swap(in[3],in[4]);
 }

 template<class T> inline void endianswap(T& in) {
   endianswap_( reinterpret_cast<char*>(&in),Int2Type<sizeof(T)>());
 }

} //namespace libcmatrix

#endif
