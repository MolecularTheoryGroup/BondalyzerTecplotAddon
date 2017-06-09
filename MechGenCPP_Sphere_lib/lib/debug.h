#ifndef DEBUG_H
#define DEBUG_H

#ifdef WIN32
#define snprintf _snprintf
#define __func__ __FUNCTION__
#endif

// on GNU systems dprintf is defined in stdio.h;
// the following line insures that this file
// is included *before* we override this symbol.
#include <stdio.h>
//////////////////////////////////////////
// for user documentation see debug.cpp //
//////////////////////////////////////////

// stuff that applies for both C-style and C++-style debug

// debug higher than this level is compiled out
// (to avoid slowing execution)
// uncomment next line to make this variable global
// #undef MAX_DEBUG_LEVEL
#ifndef MAX_DEBUG_LEVEL
  #define MAX_DEBUG_LEVEL 2
#endif
// a value of 0 turns debug completely off
#if(MAX_DEBUG_LEVEL>=0)
  #define DEBUG_ON_
#else
  #undef DEBUG_ON_
#endif

// singleton
class DebugLevel
{
private:
  static int level;
private:
  //DebugLevel(){}; // private constructor
  DebugLevel(const DebugLevel&); // prevent copy-construction
  DebugLevel& operator=(const DebugLevel&); // prevent assignment
public:
  DebugLevel(int level_){ level=level_; } // "constructor" sets (default) level.
  static void set(int level_) { level=level_; }
  static int get() { return level; }
};

// stuff specific to C-style debug

void errmsg_printf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
void eprintf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
void Wprintf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
void dprintf_fileLine(int dlevel,
                      const char *func, const char *file, int line_number, const char *format, ...);

// void eprintf_fileLine(const char* file, int line, ...);
#ifdef WIN32
// #define eprintf(args,...) error_fileLine(__func__,__FILE__,__LINE__ , ## args)
#define errmsg_printf(args,...) \
  errmsg_printf_fileLine(__func__, __FILE__, __LINE__, ## args);
#define eprintf(args,...) \
  errmsg_printf_fileLine(__func__, __FILE__, __LINE__, ## args);
#define Wprintf(args,...) \
  Wprintf_fileLine(__func__, __FILE__, __LINE__, ## args);
#else
// #define eprintf(args...) error_fileLine(__func__,__FILE__,__LINE__ , ## args)
#define errmsg_printf(args...) \
  errmsg_printf_fileLine(__func__, __FILE__, __LINE__, ## args);
#define eprintf(args...) \
  errmsg_printf_fileLine(__func__, __FILE__, __LINE__, ## args);
#define Wprintf(args...) \
  Wprintf_fileLine(__func__, __FILE__, __LINE__, ## args);
#endif
#define declare_invalid_value_error(t1) \
  void invalid_value_error_fileLine(const char* file, int line, const char* func, \
    const char* type, const char* expr, t1 val);
declare_invalid_value_error(double);
declare_invalid_value_error(int);
declare_invalid_value_error(const char*);
//#include<string>
//declare_invalid_value_error(const string&);
#define unsupported_value_error(val) invalid_value_error_fileLine( \
  __FILE__, __LINE__, __func__, "unsupported", #val, val);
#define invalid_value_error(val) invalid_value_error_fileLine( \
  __FILE__, __LINE__, __func__, "invalid", #val, val);

const char *dog_basename (const char *name);

inline bool debug_printFileLine(int dlevel, const char* func,
  const char *file, int line_number)
{
  if(DebugLevel::get() >= dlevel)
  {
    fprintf(stderr, "DEBUG(%d) %s(), %s:%d: ",
      dlevel, func, dog_basename(file), line_number);
    return true;
  }
  return false;
}

// version of dprintfn that returns true
//#define dbprintfn(dlevel, str, args...) debugn(dlevel) && \
//  printf(str, ## args), printf("\n"), true

// compile out unused debug
#ifndef DEBUG_ON
#define DEBUG_ON 1
#endif
#ifdef WIN32
#if(DEBUG_ON!=0)
  #define dprintfn(dlevel, args,...) if(DebugLevel::get() >= dlevel) \
    dprintf_fileLine(dlevel,__func__,  __FILE__, __LINE__, ## args);
  #define debugn(dlevel) debug_printFileLine(dlevel, __func__, __FILE__, __LINE__)
  #define dbprintfn(dlevel, str, args,...) debugn(dlevel) && \
    printf(str, ## args), printf("\n"), true
#else
  #define dprintfn(args,...)
  #define dbprintfn(args,...) false
  #define debugn(dlevel) false
#endif
#if(MAX_DEBUG_LEVEL>=1)
  #define dprintf1(args,...) if(DebugLevel::get() >= 1) \
    dprintf_fileLine(1, __func__, __FILE__, __LINE__, ## args);
  #define dprint1(var) if(DebugLevel::get() >= 1) \
    dprintvar_fileLine(1,__func__,__FILE__,__LINE__,#var,var);
  #define debug1 debug_printFileLine(1, __func__, __FILE__, __LINE__)
#else
  #define dprintf1(args,...)
  #define dprint1(args,...)
  #define debug1 false
#endif
#if(MAX_DEBUG_LEVEL>=2)
  #define dprintf2(args,...) if(DebugLevel::get() >= 2) \
    dprintf_fileLine(2, __func__, __FILE__, __LINE__, ## args);
  #define dprint2(var) if(DebugLevel::get() >= 2) \
    dprintvar_fileLine(2,__func__,__FILE__,__LINE__,#var,var);
  #define debug2 debug_printFileLine(2, __func__, __FILE__, __LINE__)
#else
  #define dprintf2(args,...)
  #define dprint2(args,...)
  #define debug2 false
#endif
#if(MAX_DEBUG_LEVEL>=3)
  #define dprintf3(args,...) if(DebugLevel::get() >= 3) \
    dprintf_fileLine(3, __func__, __FILE__, __LINE__, ## args);
  #define dprint3(var) if(DebugLevel::get() >= 3) \
    dprintvar_fileLine(3,__func__,__FILE__,__LINE__,#var,var);
  #define debug3 debug_printFileLine(3, __func__, __FILE__, __LINE__)
#else
  #define dprintf3(args,...)
  #define dprint3(args,...)
  #define debug3 false
#endif
#else
#if(DEBUG_ON!=0)
#define dprintfn(dlevel, args...) if(DebugLevel::get() >= dlevel) \
    dprintf_fileLine(dlevel,__func__,  __FILE__, __LINE__, ## args);
#define debugn(dlevel) debug_printFileLine(dlevel, __func__, __FILE__, __LINE__)
#define dbprintfn(dlevel, str, args...) debugn(dlevel) && \
    printf(str, ## args), printf("\n"), true
#else
#define dprintfn(args...)
#define dbprintfn(args...) false
#define debugn(dlevel) false
#endif
#if(MAX_DEBUG_LEVEL>=1)
#define dprintf1(args...) if(DebugLevel::get() >= 1) \
    dprintf_fileLine(1, __func__, __FILE__, __LINE__, ## args);
#define dprint1(var) if(DebugLevel::get() >= 1) \
    dprintvar_fileLine(1,__func__,__FILE__,__LINE__,#var,var);
#define debug1 debug_printFileLine(1, __func__, __FILE__, __LINE__)
#else
#define dprintf1(args...)
#define dprint1(args...)
#define debug1 false
#endif
#if(MAX_DEBUG_LEVEL>=2)
#define dprintf2(args...) if(DebugLevel::get() >= 2) \
    dprintf_fileLine(2, __func__, __FILE__, __LINE__, ## args);
#define dprint2(var) if(DebugLevel::get() >= 2) \
    dprintvar_fileLine(2,__func__,__FILE__,__LINE__,#var,var);
#define debug2 debug_printFileLine(2, __func__, __FILE__, __LINE__)
#else
#define dprintf2(args...)
#define dprint2(args...)
#define debug2 false
#endif
#if(MAX_DEBUG_LEVEL>=3)
#define dprintf3(args...) if(DebugLevel::get() >= 3) \
    dprintf_fileLine(3, __func__, __FILE__, __LINE__, ## args);
#define dprint3(var) if(DebugLevel::get() >= 3) \
    dprintvar_fileLine(3,__func__,__FILE__,__LINE__,#var,var);
#define debug3 debug_printFileLine(3, __func__, __FILE__, __LINE__)
#else
#define dprintf3(args...)
#define dprint3(args...)
#define debug3 false
#endif
#endif

// debug level of dprintf and dprint is intended to be contextually defined;
// define a default debug level for dprintf
// (user is expected to override this definition)
// 
#define dprintf dprintf2
#define dprint dprint2

#define declare_dprintvar_fileLine(type) \
  void dprintvar_fileLine(int,const char*,const char*,int,const char*,type);
declare_dprintvar_fileLine(int);
declare_dprintvar_fileLine(double);
declare_dprintvar_fileLine(const char*);

#ifdef WIN32
#define println(arg, args,...) \
  printf(arg, ## args); printf("\n");
#else
#define println(arg, args...) \
  printf(arg, ## args); printf("\n");
#endif
void printvar_name_val(const char*name,int val);
void printvar_name_val(const char*name,double val);
void printvar_name_val_nl(const char*name,int val);
void printvar_name_val_nl(const char*name,double val);
#define printvn(var) printvar_name_val_nl(#var,var);
#define printv(var) printvar_name_val(#var,var);
#define printvc(var) printvar_name_val(#var,var); printf(", ");

#endif
