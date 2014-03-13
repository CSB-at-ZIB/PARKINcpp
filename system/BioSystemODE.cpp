// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-03-13 td
// last changed:
//
#include "BioSystem.h"
#include "BioSystemODE.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
BioSystemODE::BioSystemODE() :
        _dim(0), _nz(0), _obj(0), _rhs() // , _adm()
{
}
//----------------------------------------------------------------------------
BioSystemODE::~BioSystemODE()
{
}
//----------------------------------------------------------------------------
void
BioSystemODE::setObj(BioSystem& obj)
{
    _obj = &obj;
    _rhs = _obj -> getODE();
    // _adm = _obj -> getMedication();
    _dim = _rhs.getSpecies().size()+1;
    _nz  = _dim;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void
BioSystemODE::computeDerivatives(Real const t, Real* y, Real* dy, int* info)
{
    y[0] = t;

    _rhs.f( y, _dim, dy );

    // _adm.computeMedication( y, _dim, dy );

    *info = 0;
}
//----------------------------------------------------------------------------
void
BioSystemODE::computeMassMatrix(Real const t, Real* y,
                                Real* B, int* ir, int* ic)
{
    // ??? _adm.computeMedication( y, _dim, dy );

    _rhs.b( y, &_nz, B, ir, ic );
}
//----------------------------------------------------------------------------
void
BioSystemODE::computeJacobian(Real const t, Real* y, Real* dy,
                            Real* J, int* ldJ, int* full_or_band,
                            int* info)
{
    Expression::Param&          sys  = _obj->getSysPar();
    BioRHS::Species const&      spec = _rhs.getSpecies();
    StrIterConst                sBeg = spec.begin();
    StrIterConst                sEnd = spec.end();

    // ??? _adm.computeMedication( y, _dim, dy );

    sys["odeTime"] = t;
    y++; // skip first component as it is reserved as time variable
    for (StrIterConst it = sBeg; it != sEnd; ++it) sys[*it] = *y++;

    _rhs.Jf( sys, _dim, J, *ldJ );

    *info = 0;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

