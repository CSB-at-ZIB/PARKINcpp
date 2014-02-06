// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 06.04.2010 td
// last changed:
//

#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <cmath>
#include <limits>

#include <addpkg/dlib/logger.h>


const long PARKINCPP_VERSION = 2000900; // <=> 02 00 09 00
const char PARKINCPP_DOTTED_VERSION[] = "2.0.9 ("  __DATE__ " " __TIME__ ")\0";


// pi ~ 3.1415926535 8979323846 2643383279 5028841971 7

/*
const double constPI      = 3.1415926535897932384626433832795;
const double constEPS     = 3.0e-16;
const double constSqrtEPS = sqrt(constEPS);
*/

const double EPMACH     = std::numeric_limits<double>::epsilon();
const double sqrtEPMACH = std::sqrt(EPMACH);
const double SMALL      = std::sqrt( std::numeric_limits<double>::min() / EPMACH );
const double GREAT      = std::sqrt( std::numeric_limits<double>::max() / 10.0 );

namespace dlib
{
    const dlib::log_level LGABBY( 40, "GABBY");
    const dlib::log_level LTALK ( 60, "TALK ");
    const dlib::log_level LVERB ( 80, "VERB ");
}

#endif
