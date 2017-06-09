#include <string>
#include <iostream>
#include <stdlib.h>
#include "debug.h"
using namespace std;
//#include "stdio.h"

// so that we can print doubles to desired precision
//
void assert_error_double(const char* file, int line, const char* func,
  const char* op, const char* lhs_str, const char* rhs_str,
  double lhs, double rhs)
{
  fprintf(stderr,"ERROR in file %s, line %d, function %s"
      "\n\tassertion failed: %s %s %s, i.e., %24.16e %s %24.16e\n",
    file, line, func, lhs_str, op, rhs_str, lhs, op, rhs);
  abort();
}

#define implement_assert_errmsg(t1,t2) \
  void assert_error(const char* file, int line, const char* func, \
    const char* op, const char* lhs_str, const char* rhs_str, \
    t1 lhs, t2 rhs) \
  { \
    std::cerr<< "ERROR in file " << file << ", line " << line  \
      << ", function " << func  \
      <<"\n\tassertion failed: " << lhs_str << op << rhs_str \
      << ", i.e., " << lhs << op << rhs << endl; \
      abort(); \
  }

implement_assert_errmsg(int,int);
implement_assert_errmsg(const char*,const char*);
implement_assert_errmsg(const string&,const string&);
//implement_assert_errmsg(string,string);

implement_assert_errmsg(double,double);
//void assert_error(const char* file, int line, const char* func,
//  const char* op, const char* lhs_str, const char* rhs_str,
//  double lhs, double rhs)
//{
//  fprintf(stderr, "ERROR, %s(), %s:%d:\n\t"
//    "assertion failed: %s %s %s, i.e., %.16e %s %.16e\n",
//    func, file, line, 
//    lhs_str, op, rhs_str, lhs, op, rhs);
//    abort();
//}

// void assert_error(const char* file, int line, const char* func,
//   const char* op, const char* lhs_str, const char* rhs_str,
//   double lhs, double rhs)
// {
//   std::cerr<< "ERROR in file " << file << ", line " << line 
//     << ", function " << func 
//     <<"\n\tassertion failed: " << lhs_str << op << rhs_str
//     << ", i.e., " << lhs << op << rhs << endl;
//   abort();
// }

//  std::cerr<< "ERROR in file " << __FILE__ << ", line " << __LINE__
//           << ", function " << __func__
//           <<"\n\tassertion failed: " #lhs "==" #rhs ", lhs="
//           <<lhs<<", rhs="<<rhs<<std::endl;
