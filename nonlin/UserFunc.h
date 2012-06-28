// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-03-23 td
// last changed:
//
#ifndef __USER_FUNC_H
#define __USER_FUNC_H

#include <cstdlib> // for rand(), RAND_MAX, ... see below!

#include "linalg/Matrix.h"
#include "linalg/Vector.h"


namespace PARKIN
{

    class UserFunc
    {
        public:
            virtual ~UserFunc() { }
            virtual Vector fcn(Vector const& x, int& ifail) { Vector f; ifail=987; return f; }
            virtual Matrix jac(Vector const& x, int& ifail) { Matrix J; ifail=987; return J; }
    };

    ///

    struct IOpt
    {
        // Real    ;
        int       iscal, mode, iterm, jacgen, boundeddamp;
        int       nonlin, rscal, itmax;
        int       mprerr, mprmon, mprsol, mprtim;
        int       transf;
        Vector    itrans;
        bool      lpos, norowscal, qrank1, qstat;
        bool      zerodat;

        //          _mprmon =   0      1      2      3      4       5       6
        //  dlib::log_level =  LNONE  LINFO  LVERB  LTALK  LGABBY  LDEBUG  LTRACE
        // logical units for logging; log levels by _mprerr, _mprmon, _mprsol, _mptim

        IOpt() :
            iscal(0), mode(0), iterm(0), jacgen(0), boundeddamp(0),
            nonlin(0), rscal(0), itmax(50),
            mprerr(1), mprmon(0), mprsol(0), mprtim(0),
            transf(0), itrans(),
            lpos(false), norowscal(false), qrank1(false), qstat(false),
            zerodat(false)
        { }
    };


    ///
    ///
    ///

    double randu();  // return a random number ( uniformly distributed in [0,1[ )
    double randn();

}
#endif // __USER_FUNC_H
