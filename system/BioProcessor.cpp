// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-09-13 td
// last changed:
//

#include "BioProcessor.h"

using namespace PARKIN;

//---------------------------------------------------------------------------
BioProcessor::BioProcessor(BioSystem* biosys, std::string const& method) :
    _biosys(biosys), _biopar(biosys),
    _method(method), _iopt(),
    _curSpecies( biosys->getSpecies() ),
    _optPar(), // _optIdx(),
    _pw(),
    _nlscon(), _nlsconWk(),
    _parkin(), _parkinWk()
{
}
//---------------------------------------------------------------------------
BioProcessor::~BioProcessor()
{
}
//---------------------------------------------------------------------------
BioProcessor::BioProcessor(BioProcessor const& other) :
    _biosys(other._biosys), _biopar(other._biopar),
    _method(other._method), _iopt(other._iopt),
    _curSpecies(other._curSpecies),
    _optPar(other._optPar), // _optIdx(other._optIdx),
    _pw(other._pw),
    _nlscon(), _nlsconWk(),
    _parkin(), _parkinWk()
{
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void
BioProcessor::setProcessingMethod(std::string const& method)
{
    _method = method;
}
//---------------------------------------------------------------------------
void
BioProcessor::setIOpt(IOpt const& iopt)
{
    _iopt = iopt;
}
//---------------------------------------------------------------------------
void
BioProcessor::setCurrentParamValues(Expression::Param const& par)
{
    Expression::ParamIterConst pBeg = par.begin();
    Expression::ParamIterConst pEnd = par.end();

    _optPar = par;

    _pw.clear();

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
    {
        _pw[it->first] = 1.0;
    }
}
//---------------------------------------------------------------------------
Expression::Param const&
BioProcessor::getCurrentParamValues()
{
    return _optPar;
}
//---------------------------------------------------------------------------
void
BioProcessor::setCurrentParamWeights(Expression::Param& par)
{
    Expression::ParamIterConst pBeg = _optPar.begin();
    Expression::ParamIterConst pEnd = _optPar.end();

    _pw.clear();

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
    {
        _pw[it->first] = par[it->first];
    }
}
//---------------------------------------------------------------------------
Expression::Param const&
BioProcessor::getCurrentParamWeights()
{
    return _pw;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
BioProcessor::TrajectoryMap
BioProcessor::computeModel()
{
    BioSystem::StrIterConst sBeg = _curSpecies.begin();
    BioSystem::StrIterConst sEnd = _curSpecies.end();
    Vector                  model;
    TrajectoryMap           trajMap;

    trajMap.clear();

    model = _biosys->computeModel( _optPar, "adaptive" );

    if ( _biosys->getComputeErrorFlag() != 0 )
    {
        return trajMap;
    }

    long k = 0;
    long K = model.nr();

    while ( k < K )
    {
        for (BioSystem::StrIterConst it = sBeg; it != sEnd; ++it)
        {
            trajMap[*it].push_back( model(++k) );
        }
    }

    return trajMap;
}
//---------------------------------------------------------------------------
BioProcessor::TrajectoryMap
BioProcessor::computeSensitivityTrajectories()
{
    int                        jacgen = _iopt.jacgen;   // 1: variational eqn
                                                        // 2: num.diff.
                                                        // 3: num.diff.(with feedback)
    BioSystem::StrIterConst    sBeg = _curSpecies.begin();
    BioSystem::StrIterConst    sEnd = _curSpecies.end();
    Expression::ParamIterConst pBeg = _optPar.begin();
    Expression::ParamIterConst pEnd = _optPar.end();
    Matrix                     mat;
    TrajectoryMap              trajMap;
    int                        ifail = 0;

    trajMap.clear();

    if ( jacgen == 1 )
    {
        mat = _biosys->computeJacobian( _optPar, "adaptive" );

        ifail = _biosys->getComputeErrorFlag();
    }
    else if ( jacgen == 2 )
    {
        mat = computeJac( ifail );
    }
    else if ( jacgen == 3 )
    {
        mat = computeJcf( ifail );
    }
    else
    {
        ifail = -1;
    }


    if ( ifail != 0 )
    {
        return trajMap;
    }

    long j = 0;
    long J = mat.nr();
    long k = 0;

    while ( j < J )
    {
        for (BioSystem::StrIterConst it = sBeg; it != sEnd; ++it)
        {
            ++j; k = 0;
            for (Expression::ParamIterConst itPar = pBeg;
                                            itPar != pEnd; ++itPar)
            {
                std::string s = *it + " / " + itPar->first;

                trajMap[s].push_back( mat(j, ++k) );
            }
        }
    }

    return trajMap;
}
//---------------------------------------------------------------------------
Vector
BioProcessor::getAdaptiveTimepoints()
{
    return _biosys->getOdeTrajectoryTimePoints();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int
BioProcessor::prepareDetailedSensitivities(Vector const& tp)
{
    return 0;
}
//---------------------------------------------------------------------------
std::vector<Matrix>
BioProcessor::getSensitivityMatrices()
{
    std::vector<Matrix> sensMat;

    return sensMat;
}
//---------------------------------------------------------------------------
std::vector<QRconDecomp>
BioProcessor::getSensitivityDecomps()
{
    std::vector<QRconDecomp> sensDcmp;

    return sensDcmp;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int
BioProcessor::identifyParameters()
{
    return 0;
}
//---------------------------------------------------------------------------
Expression::Param
BioProcessor::getIdentificationResults()
{
    Expression::Param par;

    return par;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Matrix
BioProcessor::computeJac(int& ifail)
{
    Expression::ParamIterConst pBeg = _optPar.begin();
    Expression::ParamIterConst pEnd = _optPar.end();

    Real   ajmin = SMALL;
    Real   ajdelta = std::sqrt(10.0*EPMACH);

    Matrix mat;
    Vector fh;
    Vector f = _biosys->computeModel( _optPar, "adaptive" );

    ifail = _biosys->getComputeErrorFlag();

    if ( ifail != 0 )
    {
        return mat;
    }

    long k = 0;
    mat.zeros( f.nr(), _optPar.size() );

    for (Expression::ParamIterConst itPar = pBeg;
                                    itPar != pEnd; ++itPar)
    {
        ++k;

        Real w  = itPar->second;

        int  su = ( w < 0.0 ) ? -1 : 1;
        Real u  = std::max( std::max(std::fabs(w), ajmin), _pw[itPar->first] );

        u *= (ajdelta * su);
        _optPar[itPar->first] = w + u;

        fh = _biosys->computeModel( _optPar );

        _optPar[itPar->first] = w;
        mat.set_colm(k) = (1.0/u) * Matrix(fh - f);
    }

    return mat;
}
//---------------------------------------------------------------------------
Matrix
BioProcessor::computeJcf(int& ifail)
{
    Expression::ParamIterConst pBeg = _optPar.begin();
    Expression::ParamIterConst pEnd = _optPar.end();

    Real   etadif = 1.0e-6;
    Real   etaini = 1.0e-6;
    Real   epdiff = std::sqrt(10.0 * EPMACH);
    Real   etamax = std::sqrt(epdiff);
    Real   etamin = epdiff*etamax;

    Matrix mat;
    Vector fu;
    Vector f = _biosys->computeModel( _optPar, "adaptive" );

    ifail = _biosys->getComputeErrorFlag();
    if ( ifail != 0 )
    {
        return mat;
    }

    Vector v;  v.ones( f.nr() );
    Vector eta = etaini * v;

    long m = f.nr();
    long n = _optPar.size();

    mat.zeros(m,n);
    long k = 0;

    for (Expression::ParamIterConst it = pBeg;
                                    it != pEnd; ++it)
    {
        bool is = false;
        bool qexit = false;
        bool qfine = false;

        ++k;

        while ( !qfine )
        {
            Real w  = _optPar[it->first];
            int  su = ( w < 0.0 ) ? -1 : 1;
            Real u  = eta(k) * _pw[it->first] * su;

            _optPar[it->first] = w + u;

            fu = _biosys->computeModel( _optPar );

            _optPar[it->first] = w;

            ifail = _biosys->getComputeErrorFlag();
            if ( ifail != 0 )
            {
                qexit = true; break;
            }

            Real sumd = 0.0;
            for (long j = 1; j <= m; ++j)
            {
                Real hg = std::max( std::fabs( f(j) ), std::fabs( fu(j) ) );
                Real fhj = fu(j) - f(j);
                if ( hg != 0.0 ) sumd += (fhj/hg)*(fhj/hg);

                mat(j,k) = fhj / u;
            }

            sumd = std::sqrt( sumd / m );
            qfine = true;

            if ( (sumd != 0.0) && (is == false) )
            {
                eta(k) = std::min(
                            etamax,
                            std::max( etamin, eta(k)*std::sqrt(etadif/sumd) )
                         );
                is = true;
                qfine = (sumd >= etamin);
            }
        }

        if ( qexit ) break;
    }

    return mat;
}
//---------------------------------------------------------------------------
