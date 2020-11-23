#pragma once

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

//#include "ThirdPartyHeadersBegin.h"
#  include <sstream>
#  include <stdio.h>
#  include <string>
//#include "ThirdPartyHeadersEnd.h"


#ifdef _MSC_VER

/*  Wrapper for snprintf until Microsoft supports the C++ ISO 11 Standard (post-VS2013).
 *
*/

#define snprintf std_snprintf

/* Format string according to format specification up to size-1 characters
 * and terminate buffer with null.
 * Return number of characters that would have been written if count had been sufficiently large,
 * not counting the terminating null.  If encoding error occurs, return negative number.
 * Function succeeds if return value is positive and less than size.
 * @param str 
 * Buffer to hold the formatted string.
 * @param size
 * Size of buffer.
 * @format 
 * String format specification.
 * @ap
 * List of variable arguments.
 *
*/
namespace{

inline int std_vsnprintf(char* str, size_t size, const char* format, va_list ap)
{
    int count = -1;

    if (size != 0)
        count = _vsnprintf_s(str, size, _TRUNCATE, format, ap);
    if (count == -1)
        count = _vscprintf(format, ap);

    return count;
}

inline int std_snprintf(char* str, size_t size, const char* format, ...)
{
    int count;
    va_list ap;

    va_start(ap, format);
    count = std_vsnprintf(str, size, format, ap);
    va_end(ap);

    return count;
}
} // namespace

#endif // _MSC_VER


//
// to_string
//  Convert from numeric type to string.
//
// TODO (JCB, 11/27/2013): combine with loadensight\UtilitiesString.h when revisit consolidation of string handling functions.
// Use ISO 11 Standard when upgrade to compiler that supports it.

template <typename T>
std::string to_string(T value)
{
    std::ostringstream o;
    o << value;
    return o.str();
}

