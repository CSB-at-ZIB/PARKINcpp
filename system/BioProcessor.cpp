// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-09-13 td
// last changed:
//

#include "BioProcessor.h"

using namespace PARKIN;

//---------------------------------------------------------------------------
BioProcessor::BioProcessor(BioSystem* biosys, int method) :
    _biosys(biosys), _biopar(biosys),
    _method(method), _iopt(),
    _curSpecies( biosys->getSpecies() ),
    _optPar(), // _optIdx(),
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
    _nlscon(), _nlsconWk(),
    _parkin(), _parkinWk()
{
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void
BioProcessor::setProcessingMethod(int method)
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
    _optPar = par;
}
//---------------------------------------------------------------------------
Expression::Param const&
BioProcessor::getCurrentParamValues()
{
    return _optPar;
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
    BioSystem::StrIterConst    sBeg = _curSpecies.begin();
    BioSystem::StrIterConst    sEnd = _curSpecies.end();
    Expression::ParamIterConst pBeg = _optPar.begin();
    Expression::ParamIterConst pEnd = _optPar.end();
    Matrix                     mat;
    TrajectoryMap              trajMap;

    trajMap.clear();

    mat = _biosys->computeJacobian( _optPar, "adaptive" );

    if ( _biosys->getComputeErrorFlag() != 0 )
    {
        return trajMap;
    }

    long j = 0;
    long k = 0;
    long J = mat.nr();

    while ( j < J )
    {
        for (BioSystem::StrIterConst it = sBeg; it != sEnd; ++it)
        {
            ++j; k = 0;
            for (Expression::ParamIterConst itPar = pBeg; itPar != pEnd; ++itPar)
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
