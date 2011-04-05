// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-03-23 td
// last changed:
//
#ifndef __PARKIN_LOG_H
#define __PARKIN_LOG_H

#include <iostream>
#include <string>

#include "addpkg/dlib/logger.h"


namespace PARKIN
{
    void printl(
                    dlib::logger&           lu,
                    dlib::log_level const&  lvl,
                    std::string const       fmt,
                                            ...
               );

     void print_parkin_logger_header(
                    std::ostream&           o,
                    std::string const&      n,
                    dlib::log_level const&  l,
                    dlib::uint64            j
                    );

}
#endif // __PARKIN_LOG_H
