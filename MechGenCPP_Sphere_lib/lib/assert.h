/*
 * This file is intended to supersede the system assert.h files
 * Unlike other header files, "assert.h" may usefully be included
 * multiple times, with and without NDEBUG defined.
 */

#include <stdio.h> // for printf
#include <stdlib.h> // for abort
//#include <sys/cdefs.h>

#ifdef WIN32
#define __func__ __FUNCTION__
#endif

// asserts higher than MAX_ASSERT_LEVEL are compiled out
// (to avoid slowing execution)
// uncomment next line to make this variable global
// #undef MAX_ASSERT_LEVEL
#ifndef MAX_ASSERT_LEVEL
  #define MAX_ASSERT_LEVEL 2
#endif

// for consistency with functionality of assert.h
#ifdef NDEBUG // completely turns off assert statements
// if debug is turned off create vacuous versions
// of everything that user might use
// (a list of what we currently actually support;
// other stuff below is experimental)
#define assert(e) ((void)0)
#define assert1(e) ((void)0)
#define assert2(e) ((void)0)
#define assert3(e) ((void)0)
#ifdef WIN32
#define assert_printf(args,...) ((void)0)
#define assert_printf1(args,...) ((void)0)
#define assert_printf2(args,...) ((void)0)
#define assert_printf3(args,...) ((void)0)
#else
#define assert_printf(args...) ((void)0)
#define assert_printf1(args...) ((void)0)
#define assert_printf2(args...) ((void)0)
#define assert_printf3(args...) ((void)0)
#endif
#define assert_almost_eq(a) ((void)0)
#define assert_eq(a,b) ((void)0)
#define assert_ne(a,b) ((void)0)
#define assert_gt(a,b) ((void)0)
#define assert_lt(a,b) ((void)0)
#define assert_ge(a,b) ((void)0)
#define assert_le(a,b) ((void)0)
#define assert_isnum(a) ((void)0)

#else // ifndef NDEBUG

// override system assert.h
//#define assert_fileLine(e, file, line) \
//    ((void)printf ("%s:%u: failed assertion `%s'\n", file, line, e), abort())
//void eprintf_fileLine(const char *func, const char *file, int line_number,
//  const char *format, ...);
#define dassert_fileLine(e, file, line, func) \
    (void)(printf("ERROR: file %s, line %d, function %s:\n\tfailed assertion: (%s)\n", file, line, func,e),abort())
#ifdef WIN32
#define dassert_printf_fileLine(e, file, line, func, args,...) \
    (void)(printf("ERROR: file %s, line %d, function %s:\n\tfailed assertion: (%s)\n\t", file, line, func,e), printf(args), printf("\n"), abort())
#else
#define dassert_printf_fileLine(e, file, line, func, args...) \
    (void)(printf("ERROR: file %s, line %d, function %s:\n\tfailed assertion: (%s)\n\t", file, line, func,e), printf(args), printf("\n"), abort())
#endif

// comment out the next line if __builtin_expect causes problems
#define USE_GCC_OPTIMIZATION
#ifndef USE_GCC_OPTIMIZATION
// override system assert.h
//#define assert(e)  \
//    ((void) ((e) ? (void)0 : assert_fileLine (#e, __FILE__, __LINE__)))
#define dassert_(e)  \
    ((void) ((e) ? (void)0 : dassert_fileLine(#e, __FILE__, __LINE__, __func__)))
#ifdef WIN32
#define dassert_printf_(e, args,...)  \
    ((void) ((e) ? (void)0 : dassert_printf_fileLine(#e, __FILE__, __LINE__, __func__,##args)))
#else
#define dassert_printf_(e, args...)  \
    ((void) ((e) ? (void)0 : dassert_printf_fileLine(#e, __FILE__, __LINE__, __func__,##args)))
#endif
#else // ifdef USE_GCC_OPTIMIZATION
// optimized version of preceding
//#define assert(e)  \
//    (__builtin_expect(!(e), 0) ? assert_fileLine (#e, __FILE__, __LINE__) : (void)0)
#define dassert_(e)  \
    (__builtin_expect(!(e), 0) ? dassert_fileLine (#e, __FILE__, __LINE__, __func__) : (void)0)
#ifdef WIN32
#define dassert_printf(e, args,...)  \
    (__builtin_expect(!(e), 0) ? dassert_printf_fileLine (#e, __FILE__, __LINE__, __func__,##args) : (void)0)
#else
#define dassert_printf(e, args...)  \
    (__builtin_expect(!(e), 0) ? dassert_printf_fileLine (#e, __FILE__, __LINE__, __func__,##args) : (void)0)
#endif
#endif // USE_GCC_OPTIMIZATION

#if(MAX_ASSERT_LEVEL>=1)
  #define assert1 dassert_
  #define assert_printf1 dassert_printf
#else
  #define assert1(e)
#ifdef WIN32
  #define assert_printf1(args,...)
#else
  #define assert_printf1(args...)
#endif
#endif
#if(MAX_ASSERT_LEVEL>=2)
  #define assert2 dassert_
  #define assert  dassert_
  #define assert_printf2 dassert_printf
  #define assert_printf  dassert_printf
#else
  #define assert2(e)
  #define assert(e)
#ifdef WIN32
  #define assert_printf2(args,...)
  #define assert_printf(args,...)
#else
  #define assert_printf2(args...)
  #define assert_printf(args...)
#endif
#endif
#if(MAX_ASSERT_LEVEL>=3)
  #define assert3 dassert_
  #define assert_printf3 dassert_printf
#else
  #define assert3(e)
#ifdef WIN32
  #define assert_printf3(args,...)
#else
  #define assert_printf3(args...)
#endif
#endif

// asserting specific relationships

void assert_error_double(const char* file, int line, const char* func,
  const char* op, const char* lhs_str, const char* rhs_str,
  double lhs, double rhs);
#define declare_assert_errmsg(t1,t2) \
  void assert_error(const char* file, int line, const char* func, \
    const char* op, const char* lhs_str, const char* rhs_str, \
    t1 lhs, t2 rhs); 
declare_assert_errmsg(double,double); // this seems enough for all numbers
declare_assert_errmsg(int,int); // but maybe this is more efficient
declare_assert_errmsg(const char*,const char*);
// put in assert_string.h:
//#include "assert.h"
//#include<string>
//declare_assert_errmsg(const string&,const string&);

extern "C" {
int fcmp(double x1, double x2, double epsilon);
}
#ifndef USE_GCC_OPTIMIZATION
#define builtin_expect(a,b) (a)
#else
#ifdef WIN32
#define builtin_expect(a,b) (a)
#define __builtin_expect(a,b) builtin_expect(a,b)
#else
#define builtin_expect(a,b) __builtin_expect(a,b)
#endif
#endif
#define assert_not_almost_eq(lhs,rhs) \
  (fcmp(lhs, rhs, 1e-14) \
    ? (void)0 \
    : assert_error_double(__FILE__, __LINE__, __func__, " !=~= ", #lhs, #rhs, lhs, rhs))
#define assert_almost_eq(lhs,rhs) \
  (builtin_expect(fcmp(lhs, rhs, 1e-14),0) \
    ? assert_error_double(__FILE__, __LINE__, __func__, " =~= ", #lhs, #rhs, lhs, rhs) \
    : (void)0)
#define assert_divides(lhs,rhs) \
  (builtin_expect(rhs%lhs,0) \
    ? assert_error(__FILE__, __LINE__, __func__, "(divides)", #lhs, #rhs, lhs, rhs) \
    : (void)0)
#define assert_streq(lhs,rhs) \
  (builtin_expect(strcmp(lhs,rhs),0) \
    ? assert_error(__FILE__, __LINE__, __func__, "==", #lhs, #rhs, lhs, rhs) \
    : (void)0)
#define assert_op(op,lhs,rhs) \
  (builtin_expect(!(lhs op rhs),0) \
    ? assert_error(__FILE__, __LINE__, __func__, #op, #lhs, #rhs, lhs, rhs) \
    : (void)0)

#define assert_eq(a,b) assert_op(==,a,b);
#define assert_ne(a,b) assert_op(!=,a,b);
#define assert_gt(a,b) assert_op(>,a,b);
#define assert_lt(a,b) assert_op(<,a,b);
#define assert_ge(a,b) assert_op(>=,a,b);
#define assert_le(a,b) assert_op(<=,a,b);
// should have custom implementation
#define assert_isnum(a) assert_eq(a,a);

#endif // NDEBUG
