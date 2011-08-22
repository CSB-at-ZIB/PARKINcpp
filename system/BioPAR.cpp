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
    _parameter(),
    _optPar(),
    _optIdx()
    // , _invPTrafo(id)
{
    setCurrentParameter( biosys->getParameters() );
}
//----------------------------------------------------------------------------
BioPAR::BioPAR(BioSystem* biosys, BioSystem::Parameter const& param) :
    _bioSystem(biosys),
    _parameter(),
    _optPar(),
    _optIdx()
    // , _invPTrafo(id)
{
    setCurrentParameter( param );
}
//----------------------------------------------------------------------------
BioPAR::~BioPAR()
{
}
//----------------------------------------------------------------------------
void
BioPAR::setCurrentParameter(BioSystem::Parameter const& param)
{
    BioSystem::StrIterConst pBeg = param.begin();
    BioSystem::StrIterConst pEnd = param.end();
    long                    k = 0;

    _parameter = param;

    _optPar.clear();
    _optIdx.clear();
    for (BioSystem::StrIterConst it = pBeg; it != pEnd; ++it)
    {
        _optPar[*it] = 0.0;
        _optIdx[*it] = ++k;
    }
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

    // if ( ifail == 998 ) mode = "adaptive";

    for (BioSystem::StrIterConst it = pBeg; it != pEnd; ++it)
    {
        _optPar[*it] = x(++k);
    }

    Vector v;
    v = _bioSystem -> computeModel( _optPar );
    ifail = _bioSystem -> getComputeErrorFlag();

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

    J = _bioSystem -> computeJacobian( _optPar, _optIdx );
    ifail = _bioSystem -> getComputeErrorFlag();

    return J;
}
//----------------------------------------------------------------------------
