// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-03-13 td
// last changed:
//
#include "linalg/Matrix.h"
#include "BioSystem.h"
#include "BioSystemVAR.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
BioSystemVAR::BioSystemVAR() :
        _dim(0), _qq(0), _nz(0), _obj(0), _rhs()
{
}
//----------------------------------------------------------------------------
BioSystemVAR::~BioSystemVAR()
{
}
//----------------------------------------------------------------------------
void
BioSystemVAR::setObj(BioSystem& obj)
{
    _obj = &obj;
    _rhs = _obj -> getODE();
    _dim = _rhs.getSpecies().size() + 1;
    _qq  = _obj -> getOptPar().size();
    _nz  = _dim * (_qq + 1);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
BioSystemVAR::computeDerivatives(Real const t, Real* y, Real* dy, int* info)
{
    Expression::Param&          opt = _obj->getOptPar();
    // int                         qq = opt.size();   // ((long)(*nq) - n) / n;
    double*                     yy = y + _dim;

    //_nz = _dim*(qq + 1);
    y[0]  =
    yy[0] = t;

    _rhs.df( opt, yy, y, _dim, _qq, dy + _dim );
    _rhs.f( y, _dim, dy );

    *info = 0;
}
//----------------------------------------------------------------------------
void
BioSystemVAR::computeMassMatrix(Real const t, Real* y,
                                Real* B, int* ir, int* ic)
{
    const Real one = 1.0;

    for (int j = 1; j <= _nz; ++j)
    {
        *B++ = one;
        *ir++ = *ic++ = j;
    }
}
//----------------------------------------------------------------------------
void
BioSystemVAR::computeJacobian(Real const t, Real* y, Real* dy,
                            Real* J, int* ldJ, int* full_or_band,
                            int* info)
{
    Expression::Param&          sys  = _obj->getSysPar();
    BioRHS::Species const&      spec = _rhs.getSpecies();
    StrIterConst                sBeg = spec.begin();
    StrIterConst                sEnd = spec.end();

    sys["odeTime"] = t;
    y++; // skip first component as it is reserved as time variable
    for (StrIterConst it = sBeg; it != sEnd; ++it) sys[*it] = *y++;

    Matrix Fz = _rhs.Jf( sys );
    long   n  = Fz.nr();
    // long   gap = (long)(*ldJ) - n;

    if ( *full_or_band == 0 )
    {

        long  q1 = (long)(_nz) / n;  // q1 := q + 1; loop below counts til m < q1 !

        for (long k = 1; k <= n; ++k)
        {
            for (long j = 1; j <= n; ++j)
            {
                double tmp = Fz(j,k);

                for (long m = 0; m < q1; ++m)
                {
                    J[ (*ldJ)*(k-1 + m*n) + (j-1 + m*n) ] = tmp;
                                //  J(j + q*n, k + q*n) = tmp
                }
            }
        }

    }
    else
    {
        // long  nn = 2*n + 1L;
        long  q1 = (long)(_nz) / n;  // q1 := q + 1; loop below counts til m < q1 !

        for (long k = 1; k <= n; ++k)
        {
            long off = n + 1L - k;

            for (long j = 1; j <= n; ++j)
            {
                double tmp = Fz(j,k);

                for (long m = 0; m < q1; ++m)
                {
                    J[ (*ldJ)*(k-1 + m*n) + (j-1 + off) ] = tmp;
                                //  J(j + q*n, k + q*n) = tmp
                }
            }
        }
    }

    *info = 0;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

