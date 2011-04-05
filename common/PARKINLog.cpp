// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-03-23 td
// last changed:
//
#include <cstdio>           // vsprintf
#include <cstdarg>          // va_list, va_start, etc.
#include <sstream>          // stringstream

#include "PARKINLog.h"

// using namespace PARKIN;

//---------------------------------------------------------------------------
void
PARKIN::printl(
    dlib::logger&           lu,
    dlib::log_level const&  lvl,
    std::string const       fmt,
                            ...
    )
{
    char            buffer[8196];
    std::string     line;
    std::va_list    args;

    va_start(args, fmt);
    std::vsprintf(buffer, fmt.c_str(), args);
    va_end(args);

    std::stringstream stream(buffer);
    while ( std::getline(stream,line) ) lu << lvl << line;

    return;
}
//---------------------------------------------------------------------------
void
PARKIN::print_parkin_logger_header(
    std::ostream&           o,
    std::string const&      n,
    dlib::log_level const&  l,
    dlib::uint64            j
    )
{
    o << std::setw(8) << n << " | ";
    return;
}
//---------------------------------------------------------------------------
