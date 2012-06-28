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
    , _zerodat(false)
    , _absRes()
    , _relRes()
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
    , _zerodat(false)
    , _absRes()
    , _relRes()
{
    setCurrentParameter( param );
}
//----------------------------------------------------------------------------
BioPAR::~BioPAR()
{
}
//----------------------------------------------------------------------------
void
BioPAR::setEvalMode(bool zerodat)
{
    _zerodat = zerodat;
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

    //

    if ( _zerodat == true )
    {
        Vector d = _bioSystem -> getMeasurements();
        Vector w = _bioSystem -> getMeasurementWeights();

        // for (long j = 1; j <= w.nr(); ++j)
        // {
        //     w(j) = std::max( std::fabs(d(j)), w(j) );
        // }

        // _absRes = v - d;
        // _relRes = _absRes.dotdiv(w);
        v = v - d;
        v = v.dotdiv(w);

//std::cerr << std::endl;
//std::cerr << "*** BioPAR::fcn():" << std::endl;
//std::cerr << "  v = " << std::endl;
//std::cerr << v << std::endl;

        // return (v - d).dotdiv(w);
    }

    //


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

    if ( _zerodat == true )
    {
        const Real  relPert = _bioSystem->getSolverRTol(); // 1.0e-05;
        const Real  absPert = _bioSystem->getSolverATol(); // 1.0e-05;

        Vector      v;
        Vector      v0 = _bioSystem -> computeModel( _optPar );
        Vector      w  = _bioSystem -> getMeasurementWeights();
        // Vector      d  = _bioSystem -> getMeasurements();

        // v0 = v0 - d;
        J.zeros(v0.nr(), x.nr());
        k = 0;

        for (BioSystem::StrIterConst it = pBeg; it != pEnd; ++it)
        {
            Expression::Param   hPar = _optPar;
            // Real                sign = (_optPar[*it] < 0.0) ? -1.0 : 1.0;
            // Real                dp   = sign * _optPar[*it] * relPert + absPert;
            Real                dp   = relPert * std::fabs(_optPar[*it]) + absPert;

            hPar[*it] += dp;

            v = _bioSystem -> computeModel( hPar );
            ifail = _bioSystem -> getComputeErrorFlag();

            ++k;

            if ( ifail == 0 )
            {
                // v = v - d;
                v = (1.0/dp)*(v - v0);

                J.set_colm(k) = v.dotdiv(w);
            }

            hPar[*it] -= dp;
        }

//std::cerr << std::endl;
//std::cerr << "*** BioPAR::jac() :" << std::endl;
//std::cerr << "  J = " << std::endl;
//std::cerr << J << std::endl;

        return J;
    }

    J = _bioSystem -> computeJacobian( _optPar, _optIdx );
    ifail = _bioSystem -> getComputeErrorFlag();

    return J;
}
//----------------------------------------------------------------------------
