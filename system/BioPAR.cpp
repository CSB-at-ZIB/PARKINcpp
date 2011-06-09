// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-17 td
// last changed:
//

#include "BioPAR.h"

using namespace PARKIN;

//----------------------------------------------------------------------------
BioPAR::BioPAR(BioSystem* biosys) :
    _bioSystem(biosys),
    _parameter(biosys->getParameters()),
    _optPar()
    // , _invPTrafo(id)
{
    // _bioSystem -> setParameters(_parameter);
}
//----------------------------------------------------------------------------
BioPAR::BioPAR(BioSystem* biosys, BioSystem::Parameter const& param) :
    _bioSystem(biosys),
    _parameter(param),
    _optPar()
    // , _invPTrafo(id)
{
    // _bioSystem -> setParameters(_parameter);
}
//----------------------------------------------------------------------------
BioPAR::~BioPAR()
{
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
Vector
BioPAR::fcn(Vector const& x, int& ifail)
{
    // Expression::Param       optPar;
    BioSystem::StrIterConst pBeg = _parameter.begin();
    BioSystem::StrIterConst pEnd = _parameter.end();
    long                    k = 0;

    // if ( ifail == 998 ) mode = "init";

    for (BioSystem::StrIterConst it = pBeg; it != pEnd; ++it)
    {
        _optPar[*it] = x(++k);
    }

    Vector v;
    v = _bioSystem -> computeModel( _optPar );
    ifail = 0;

    return v;
}
//----------------------------------------------------------------------------
Matrix
BioPAR::jac(Vector const& x, int& ifail)
{
    // Expression::Param       optPar;
    BioSystem::StrIterConst pBeg = _parameter.begin();
    BioSystem::StrIterConst pEnd = _parameter.end();
    long                    k = 0;

    for (BioSystem::StrIterConst it = pBeg; it != pEnd; ++it)
    {
        _optPar[*it] = x(++k);
    }

    Matrix J;

    J = _bioSystem -> computeJacobian( _optPar );
    ifail = 0;

    return J;
}
//----------------------------------------------------------------------------
