// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-10 td
// last changed:
//
%module linalg

%{
#include <common/Constants.h>
#include <common/Types.h>
#include <linalg/Matrix.h>
#include <linalg/Vector.h>

#include <linalg/QRDecomp.h>
#include <linalg/QRconDecomp.h>

#include <sstream>
%}

//
%constant double EPMACH     = std::numeric_limits<double>::epsilon();
%constant double sqrtEPMACH = std::sqrt(EPMACH);
%include <std_string.i>

//
%ignore PARKIN::Matrix::operator();
%ignore PARKIN::Matrix::operator=;
%ignore PARKIN::Mref::operator=;

%ignore PARKIN::Matrix::operator();
%ignore PARKIN::Vector::operator=;
%ignore PARKIN::Vref::operator=;

//
%include "common/Types.h"
%include "linalg/GenericMatrixImpl.h"
%include "linalg/Matrix.h"
%include "linalg/Vector.h"


//
// %ignore PARKIN::dcmpErr;
%ignore PARKIN::GenericDecomp::operator=;
//%ignore PARKIN::QRDecomp::operator=;
//%ignore PARKIN::QRconDecomp::operator=;

%include "linalg/GenericDecomp.h"
%include "linalg/GenericQRPseudoInv.h"
%include "linalg/QRDecomp.h"
%include "linalg/QRconDecomp.h"


//
%extend PARKIN::Matrix {

    Real __getitem__(PyObject* args)
    {
        int j, k;
        int success = PyArg_ParseTuple(args, "ii", &j, &k);
        if (!success) return 0.0;
        if ( j < 0 ) j += (*self).nr();
        if ( k < 0 ) k += (*self).nc();
        return (*self)(1+j,1+k);
    }

    Real __getitem__(long l)
    {
        // long m = (*self).nr();
        long n = (*self).nc();
        int j,k;
        j = l / n;
        k = l % n;
        return (*self)(1+j,1+k);
    }

    void __setitem__(PyObject* args, Real val)
    {
        int j, k;
        int success = PyArg_ParseTuple(args, "ii", &j, &k);
        if (!success) return;
        if ( j < 0 ) j += (*self).nr();
        if ( k < 0 ) k += (*self).nc();
        (*self)(1+j,1+k) = val;
    }

    void __setitem__(long l, Real val)
    {
        // long m = (*self).nr();
        long n = (*self).nc();
        int j,k;
        j = l / n;
        k = l % n;
        (*self)(1+j,1+k) = val;
    }

    std::string __str__()
    {
        std::ostringstream s;
        s << (*self);
        return s.str();
    }

};

//
%extend PARKIN::Vector {

    PARKIN::Vector __getslice__(long j1, long j2)
    {
        long m = (*self).nr();
        if ( (j1==0) && (j2>=m) ) return (*self).row(1,m);
        if ( j1 < 0 ) j1 += m;
        if ( j2 >= m ) j2 = m-1;
        if ( j2 < 0 )  j2 += m-1;
        return (*self).row(1+j1,1+j2);
    }

    Real __getitem__(long j)
    {
        if ( j < 0 ) j += (*self).nr();
        return (*self)(1+j);
    }

    void __setitem__(long j, Real val)
    {
        if ( j < 0 ) j += (*self).nr();
        (*self)(j+1) = val;
    }

    long __len__()
    {
        return (*self).nr();
    }

    std::string __str__()
    {
        std::ostringstream s;
        s << (*self);
        return s.str();
    }

};
