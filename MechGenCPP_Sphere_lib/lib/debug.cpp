#include <stdio.h>
#include <stdlib.h> /* for exit and abort */
#include <stdarg.h>
#ifdef WIN32
#include <io.h>
#include <process.h>
#else
#include <unistd.h> /* for isatty */
#endif
#include "debug.h"

// implemented here because the new
// standard basename() in libgen.h eliminates the const
// (because standard basename returns non-const? -- seems
// like an attempted correction in the wrong direction...).
const char *dog_basename (const char *name)
{
  const char *base;
  for (base = name; *name; name++)
    if (*name == '/') base = name + 1;
  return base;
}

// === documentation of debug code ===
//
// For an example, see the test code:
//   dogpack/lib/tests/test_debug.cpp
//
// eprintf passes its arguments to eprintf_fileLine;
// derr causes eprintf_fileLine to be called when it goes
//   out of scope.
// 
// I have errmsg_printf_fileLine defined to call a script
//   eprintf_script which (by uncommenting a line)
//   can be configured to email the error message to the user.
// 
// To use debug, you should
//   #include "debug.h"
// which defines eprintf and dprintf.
// 
// If you insist on using streams (which I would discourage,
// since including the stream libraries increases the time
// necessary to compile an object module by over a factor of 6),
// you should
//   #include "DebugStream.h"
// to get dout (which works like cout) and
//   #include "ErrorStream.h"
// to get derr (which works like cerr, except that
// when derr goes out of scope eprintf_fileLine gets called).
// As a convenience, the line
//   #include "DebugStreams.h"
// will include both DebugStream.h and ErrorStream.h.
// 
// I have implemented three levels of debug.
// 
//   dprintf1 (or dout1) can be used for important debug,
//   dprintf2 (or dout2) can be used for standard debug, and
//   dprintf3 (or dout3) is for verbose debug
//     (such as from the inside of a loop).
//
// By default dprintf is aliased to dprintf2 and dout is
// aliased to dout2. The expectation is that these symbols are
// contextually defined per file with these default values.
// Adding the line
//   #define dout dout3
// *after* including debug.h would override this default
// 
// In debug.h notice the line
//   #define MAX_DEBUG_LEVEL 2
// If you wish to see level-3 debug you should change this to
//   #define MAX_DEBUG_LEVEL 3
// If you make this change in this file it will cause level-3
// debug to be shown universally.
// If you include the line
//   #define MAX_DEBUG_LEVEL 3
// in a particular file (*before* including debug.h)
// it will redefine the maximum debug level for that particular file.
// 
// The purpose of this constant is to provide a mechanism that
// (1) allows for file-level debug level settings and
// (2) allows debug to be removed at compile time.
// 
// The setting
//   #define MAX_DEBUG_LEVEL 2
// will cause any appearances of dprintf3 to be removed
// at compile time.  Appearances of dout3 will not be
// removed, but will be replaced with "if(0) ...",
// which (if the optimizer has any intelligence  at all)
// will incur no run-time cost.
// 
// To change the debug level universally
// at *run* time, you should call, e.g.,
// 
//   DebugLevel::set(1)
// 
// if you wish to change the debug level from its default
// value of 2.
//
// === end of documentation of debug code ===

/* To do:
   - modify to append to a log file
   - create personalized assert which email the user.
 */

int DebugLevel::level=2;  // set default debug level.

// int debug_level=2; /* default debug level */
void dprintf_fileLine(int dlevel,
  const char *func, const char *file, int line_number,
  const char *format, ...)
{
  // This test is unnecessary since the macro 
  // already takes care of it.
  //if(DebugLevel::get() < dlevel) return;
  va_list args;
  va_start(args, format);
  fprintf(stderr, "DEBUG(%d) %s(), %s:%d: ",
    dlevel, func, dog_basename(file), line_number);
  /* print out remainder of message */
  vfprintf(stderr, format, args);
  va_end(args);
  printf("\n");
}

void Wprintf_fileLine(const char *func, const char *file, int line_number,
  const char *format, ...)
{
  //if(DebugLevel::get() < dlevel) return;
  va_list args;
  va_start(args, format);
  fprintf(stderr, "WARNING, %s(), %s:%d:\n\t",
    func, dog_basename(file), line_number);
  /* print out remainder of message */
  vfprintf(stderr, format, args);
  va_end(args);
  printf("\n");
}

void eprintf_fileLine(const char *func, const char *file, int line_number,
  const char *format, ...)
{
  va_list args;
  va_start(args, format);
  fprintf(stderr, "ERROR, %s(), %s:%d:\n\t",
    func, file, line_number);
  /* print out remainder of message */
  vfprintf(stderr, format, args);
  va_end(args);
  printf("\n");
  abort();
}

/*  An example of a typical escaped string:
 *  echo "\\\"\`'()[]{}<>~@#$|\&*"
 */
static void escape_special_shell_characters(
  const char* message, char* escaped_msg)
{
  char c = '\0';
  const char* message_ptr;
  char* escaped_msg_ptr;

  message_ptr = message;
  escaped_msg_ptr = escaped_msg;

  *escaped_msg_ptr++ = '"';
  while(c==*message_ptr++)
  {
    switch(c)
    {
      /* escape characters that are special even inside double quotes */
      case '$':
      case '`':
      case '"':
      case '\\':
        *escaped_msg_ptr++ = '\\';
    }
    *escaped_msg_ptr++ = c;
  }
  *escaped_msg_ptr++ = '"';
  *escaped_msg_ptr++ = '\0';
}


/*
 *  quit and optionally send email message
 */
void errmsg_printf_fileLine(const char *func, const char *file, int line_number,
  const char *format, ...)
{
  va_list args;
  va_start(args, format);
  fprintf(stderr, "ERROR in function %s, file %s, line %d: \n\t",
    func, file, line_number);
  /* print out remainder of message */
  vfprintf(stderr, format, args);
  va_end(args);
  // append terminating newline so user does not have to do it
  fprintf(stderr, "\n");

  // create and send message
  const int BUFFER_SIZE = 1024;
  // if the output is detached from a terminal email
  // the error message to the user.
  int exit_status = 0;
  if(!isatty(fileno(stdout)))
  {
    char message[BUFFER_SIZE];
    char escaped_msg[2*BUFFER_SIZE+2];
    char shell_command[128+2*BUFFER_SIZE];
    char *str_pos;
    int pos;
    pos = snprintf(message, BUFFER_SIZE,
      "ERROR in function %s, file %s, line %d: \n\t",func, file, line_number);
    va_start(args, format);
    pos+=vsnprintf(message+pos,BUFFER_SIZE-pos, format, args);
    va_end(args);
    // append terminating newline to message
    snprintf(message+pos, BUFFER_SIZE-pos, "\n");

    escape_special_shell_characters(message,escaped_msg);
    str_pos = shell_command;
    str_pos+= sprintf(str_pos, 
      "${DOGPACK}/scripts/errmsg_printf_script ");
    sprintf(str_pos, "%s", escaped_msg);
    
    fprintf(stderr, "executing shell command: %s", shell_command);
    exit_status = system(shell_command);
  }

  abort();
}

#include "VoidStream.h"
// is this "instantiation" of this vacuous class necessary?
// I imagine not, since all of its (nonexistent) data is static
VoidStream voidStream;

#include <iostream>
using namespace std;
#define implement_invalid_value_error(t1) \
  void invalid_value_error_fileLine(const char* file, int line, const char* func, \
    const char* type, const char* expr, t1 val) \
  { \
    std::cerr<< "ERROR in file " << file << ", line " << line  \
      << ", function " << func  \
      <<"\n\t" << type << " value: " << expr << " = " << val << endl; \
      abort(); \
  }

implement_invalid_value_error(double);
implement_invalid_value_error(int);
implement_invalid_value_error(const char*);
#include<string>
implement_invalid_value_error(const string&);

#define implement_dprintvar_fileLine(code,type) \
  void dprintvar_fileLine(int dlevel, const char* func, const char* file, int line, \
    const char* name, type val) \
  { \
    dprintf_fileLine(dlevel,func,file,line, \
      "%s == " code,name,val); \
  }
implement_dprintvar_fileLine("%s",const char*);
implement_dprintvar_fileLine("%d",int);
implement_dprintvar_fileLine("%24.16e",double);

void printvar_name_val_nl(const char*name,int val)
{ printf("%24d =: %s\n",val,name); }
void printvar_name_val_nl(const char*name,double val)
{ printf("%24.16e =: %s\n",val,name); }

void printvar_name_val(const char*name,int val)
{ printf("%s = %d",name,val); }
void printvar_name_val(const char*name,double val)
{ printf("%s = %24.16e",name,val); }

// void printvar_name_val_nl(const char*name,int val)
// { printf("%24s := %d\n",name,val); }
// void printvar_name_val_nl(const char*name,double val)
// { printf("%24s := %24.16e\n",name,val); }
