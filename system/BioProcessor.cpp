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
    _speThres(), _parThres(),
    _trajMap(), _sensTraj(0),
    _sensiMat(), _sensiDcmp(),
    _nlscon(), _nlsconWk(),
    _parkin(), _parkinWk()
{
}
//---------------------------------------------------------------------------
BioProcessor::~BioProcessor()
{
    delete _sensTraj;
}
//---------------------------------------------------------------------------
BioProcessor::BioProcessor(BioProcessor const& other) :
    _biosys(other._biosys), _biopar(other._biopar),
    _method(other._method), _iopt(other._iopt),
    _curSpecies(other._curSpecies),
    _optPar(other._optPar), // _optIdx(other._optIdx),
    _speThres(other._speThres), _parThres(other._parThres),
    _trajMap(), _sensTraj(0),
    _sensiMat(other._sensiMat), _sensiDcmp(other._sensiDcmp),
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
    // bool                        lpos = _iopt.lpos;
    Expression::ParamIterConst  pBeg = par.begin();
    Expression::ParamIterConst  pEnd = par.end();

    /*
    if ( lpos == true )
    {
        _optPar.clear();

        for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
        {
            Real tmp = it->second;

            _optPar[it->first] = (tmp > 0.0) ? std::log(tmp) : -1.0e-38;
        }
    }
    else
    */
    {
        _optPar = par;
    }


    _parThres.clear();

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
    {
        _parThres[it->first] = 0.0;
    }
}
//---------------------------------------------------------------------------
Expression::Param const&
BioProcessor::getCurrentParamValues()
{
    /*
    bool                        lpos = _iopt.lpos;
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    Expression::Param           par;

    par.clear();

    if ( lpos == true )
    {
        for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
        {
            Real tmp = it->second;

            par[it->first] = std::exp( tmp );
        }
    }
    else
    {
        par = _optPar;
    }

    return par;
    */

    return _optPar;
}
//---------------------------------------------------------------------------
void
BioProcessor::setCurrentParamThres(Expression::Param const& par)
{
    Expression::ParamIterConst pBeg = _optPar.begin();
    Expression::ParamIterConst pEnd = _optPar.end();

    _parThres.clear();

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
    {
        _parThres[it->first] = it->second;
    }
}
//---------------------------------------------------------------------------
Expression::Param const&
BioProcessor::getCurrentParamThres()
{
    return _parThres;
}
//---------------------------------------------------------------------------
void
BioProcessor::setCurrentSpeciesThres(Expression::Param& par)
{
    BioSystem::StrIterConst sBeg = _curSpecies.begin();
    BioSystem::StrIterConst sEnd = _curSpecies.end();

    _speThres.clear();

    for (BioSystem::StrIterConst it = sBeg; it != sEnd; ++it)
    {
        _speThres[*it] = par[*it];
    }
}
//---------------------------------------------------------------------------
Expression::Param const&
BioProcessor::getCurrentSpeciesThres()
{
    return _speThres;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
BioProcessor::TrajectoryMap
BioProcessor::computeModel()
{
    BioSystem::StrIterConst sBeg = _curSpecies.begin();
    BioSystem::StrIterConst sEnd = _curSpecies.end();
    Vector                  model;
    TrajectoryMap           trajectory;

    trajectory.clear();

    model = _biosys->computeModel( _optPar, "adaptive" );

    if ( _biosys->getComputeErrorFlag() != 0 )
    {
        return trajectory;
    }

    long k = 0;
    long K = model.nr();

    while ( k < K )
    {
        for (BioSystem::StrIterConst it = sBeg; it != sEnd; ++it)
        {
            trajectory[*it].push_back( model(++k) );
        }
    }

    return trajectory;
}
//---------------------------------------------------------------------------
BioProcessor::TrajectoryMap
BioProcessor::computeSensitivityTrajectories()
{
    bool                        lpos = _iopt.lpos;
    int                         jacgen = _iopt.jacgen;   // 1: variational eqn
                                                         // 2: num.diff.
                                                         // 3: num.diff.(with feedback)
    BioSystem::StrIterConst     sBeg = _curSpecies.begin();
    BioSystem::StrIterConst     sEnd = _curSpecies.end();
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    // TrajectoryMap               trajMap;
    Matrix                      mat;
    int                         ifail = 0;

    _trajMap.clear();

    if ( jacgen == 1 )
    {
        mat = _biosys->computeJacobian( _optPar, "adaptive" );
        ifail = _biosys->getComputeErrorFlag();

        delete _sensTraj;
        _sensTraj = _biosys->getEvaluationTrajectories();

        _linftyModel = _biosys->getLinftyModel();
    }
    else if ( jacgen == 2 )
    {
        mat = computeJac( "adaptive", ifail );
    }
    else if ( jacgen == 3 )
    {
        mat = computeJcf( "adaptive", ifail );
    }
    else
    {
        ifail = -1;
    }


    if ( ifail != 0 )
    {
        return _trajMap;
    }


    if ( lpos == true )
    {
        Vector y;
        long   j = mat.nc();

        y.zeros( j );
        j = 0;

        for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
        {
            y(++j) = it->second;
        }

        mat = mat * y.diag();
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

                _trajMap[s].push_back( mat(j, ++k) );
            }
        }
    }

    return _trajMap;
}
//---------------------------------------------------------------------------
Vector
BioProcessor::getAdaptiveTimepoints()
{
    return _biosys->getOdeTrajectoryTimePoints();
}
//---------------------------------------------------------------------------
BioProcessor::TrajectoryMap
BioProcessor::getScaledSensitivityTrajectories()
{
    BioSystem::StrIterConst     sBeg = _curSpecies.begin();
    BioSystem::StrIterConst     sEnd = _curSpecies.end();
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    Expression::Param           parScale;
    Expression::Param           speScale;
    TrajectoryMap               scaledTrajMap;

    scaledTrajMap.clear();

    parScale = computeParameterScales();
    speScale = computeSpeciesScales();

    for (BioSystem::StrIterConst it = sBeg; it != sEnd; ++it)
    {
        std::string s = *it + " / ";
        Real        yScale = speScale[*it];

        if ( yScale > 0.0 )
        {
            for (Expression::ParamIterConst itPar = pBeg;
                                            itPar != pEnd; ++itPar)
            {
                s += itPar->first;

                long T = _trajMap[s].size();
                Real pScale = parScale[itPar->first];

                for (long tp = 0; tp < T; ++tp)
                {
                    scaledTrajMap[s][tp] =
                        _trajMap[s][tp] * ( pScale / yScale );
                }
            }
        }
    }

    return scaledTrajMap;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int
BioProcessor::prepareDetailedSensitivities(Vector const& tp)
{
    bool                lpos = _iopt.lpos;
    int                 jacgen = _iopt.jacgen;
    Expression::Param   parScale = computeParameterScales();
    Expression::Param   speScale = computeSpeciesScales();
    unsigned            mcon = 0;
    unsigned            irank = 0;
    Real                cond = 1.0/( _biosys->getSolverRTol() );
    long                T = tp.nr();
    long                m = _curSpecies.size();
    long                q = _optPar.size();

    _sensiMat.clear();
    _sensiDcmp.clear();

    Expression::ParamIterConst itPar = _optPar.begin();
    BioSystem::StrIterConst    itSpe = _curSpecies.begin();
    Vector                     one_fw, pw;

    one_fw.zeros(m);
    pw.zeros(q);

    for (long k = 1; k <= m; ++k)
    {
        Real tmp = speScale[*itSpe++];
        one_fw(k) = (tmp > 0.0) ? 1.0/tmp : 0.0;
    }
    if ( lpos == true )
    {
        for (long l = 1; l <= q; ++l)
        {
            pw(l) = parScale[itPar->first] * (itPar->second);
            ++itPar;
        }
    }
    else
    {
        for (long l = 1; l <= q; ++l)
        {
            pw(l) = parScale[itPar->first];
            ++itPar;
        }
    }

    if ( jacgen == 1 )
    {
        if ( _sensTraj == 0 )
        {
            return -1;
        }

        for (long j = 1; j <= T; ++j)
        {
            Vector v = _sensTraj->eval( tp(j) );
            Matrix mat;
            long   n = m;

            mat.zeros(m,q);

            if ( v.nr() == m*(q+1) )
            {
                for (long k = 1; k <= m; ++k)
                {
                    for (long l = 1; l <= q; ++l)
                    {
                        mat(k,l) = v(++n);
                    }
                }
            }

            mat = one_fw.diag() * mat * pw.diag();

            _sensiMat.push_back( mat );
            _sensiDcmp.push_back( mat.factorQRcon(mcon, irank, cond) );
        }
    }
    else if ( jacgen == 2 )
    {
        // _biosys->setEmptyMeasurementList();
        _biosys->setMeasurementTimePoints( tp );

        int ifail = 0;
        Matrix Jac = computeJac( "", ifail );

        if ( (ifail != 0) || (Jac.nr() != m*T) )
        {
            return -2;
        }

        for (long j = 1; j <= T; ++j)
        {
            Matrix mat = Jac.rowm( (j-1)*m + 1, j*m );

            mat = one_fw.diag() * mat * pw.diag();

            _sensiMat.push_back( mat );
            _sensiDcmp.push_back( mat.factorQRcon(mcon, irank, cond) );
        }
    }
    else if ( jacgen == 3 )
    {
        // _biosys->setEmptyMeasurementList();
        _biosys->setMeasurementTimePoints( tp );

        int ifail = 0;
        Matrix Jac = computeJcf( "", ifail );

        if ( (ifail != 0) || ( Jac.nr() != m*T) )
        {
            return -3;
        }

        for (long j = 1; j <= T; ++j)
        {
            Matrix mat = Jac.rowm( (j-1)*m + 1, j*m );

            mat = one_fw.diag() * mat * pw.diag();

            _sensiMat.push_back( mat );
            _sensiDcmp.push_back( mat.factorQRcon(mcon, irank, cond) );
        }
    }
    else
    {
        return -4;
    }

    return 0;
}
//---------------------------------------------------------------------------
std::vector<Matrix>
BioProcessor::getSensitivityMatrices()
{
    return _sensiMat;
}
//---------------------------------------------------------------------------
std::vector<QRconDecomp>
BioProcessor::getSensitivityDecomps()
{
    return _sensiDcmp;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
int
BioProcessor::identifyParameters(Real xtol)
{
    int                         rc = 0;
    unsigned                    j, m;
    Expression::Param           pScal;
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    BioSystem::Parameter        pname;
    Vector                      x, xscal;
    Vector                      fobs, fscal;

    fobs = _biosys->getMeasurements();
    fscal = _biosys->getMeasurementWeights();

    j = 0;
    pname.clear();
    x.zeros( _optPar.size() );

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
    {
        pname.push_back( it->first );
        x(++j) = it->second;
    }


    j = 0;
    pScal = computeParameterScales();
    xscal.zeros( pScal.size() );

    for (Expression::ParamIterConst it = pScal.begin();
                                    it != pScal.end();
                                    ++it)
    {
        xscal(++j) = it->second;
    }


    m = fobs.nr();

    _biopar = BioPAR( _biosys, pname );



    if ( _method == "parkin" )
    {
        _parkin.setProblem( &_biopar );
        _parkin.initialise( m,
                            x, xscal,
                            fobs, fscal,
                            xtol, _iopt,
                            _parkinWk
                          );

        rc = _parkin.run();

        // _parkinWk = _parkin.getWk();
        _idResult = _parkin.getSolution();
    }
    else if ( _method == "nlscon" )
    {
        _nlscon.setProblem( &_biopar );
        _nlscon.initialise( m,
                            x, xscal,
                            fobs, fscal,
                            xtol, _iopt,
                            _nlsconWk
                          );

        rc = _nlscon.run();

        // _nlsconWk = _nlscon.getWk();
        _idResult = _nlscon.getSolution();
    }
    else
    {
        rc = -999;
    }

    return rc;
}
//---------------------------------------------------------------------------
Expression::Param
BioProcessor::getIdentificationResults()
{
    long                    J = _idResult.nr();
    BioSystem::Parameter    pname = _biopar.getCurrentParameter();
    Expression::Param       par;

    par.clear();

    if ( pname.size() != (unsigned)J )
    {
        return par;
    }

    for (long j = 1; j <= J; ++j)
    {
        std::string s = pname[j-1];
        par[s] = _idResult(j);
    }

    return par;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Expression::Param
BioProcessor::computeParameterScales()
{
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    Expression::Param           parScale;

    parScale.clear();

    for (Expression::ParamIterConst itPar = pBeg;
                                    itPar != pEnd; ++itPar)
    {
        std::string s = itPar->first;

        parScale[s] =
            std::max(   std::fabs( _optPar[s] ) ,
                        _parThres[s]
                    );
    }

    return parScale;
}
//---------------------------------------------------------------------------
Expression::Param
BioProcessor::computeSpeciesScales()
{
    BioSystem::StrIterConst sBeg = _curSpecies.begin();
    BioSystem::StrIterConst sEnd = _curSpecies.end();
    Expression::Param       speScale;

    speScale.clear();

    for (BioSystem::StrIterConst it = sBeg; it != sEnd; ++it)
    {
        speScale[*it] =
            std::max(   _linftyModel[*it],
                        _speThres[*it]
                    );
    }

    return speScale;
}
//---------------------------------------------------------------------------
Matrix
BioProcessor::computeJac(std::string mode, int& ifail)
{
    Expression::ParamIterConst pBeg = _optPar.begin();
    Expression::ParamIterConst pEnd = _optPar.end();
    Expression::Param          parScale = computeParameterScales();

    Real   ajmin = SMALL;
    Real   ajdelta = std::sqrt(10.0*EPMACH);

    _linftyModel.clear();

    Matrix mat;
    Vector fh;
    Vector f = _biosys->computeModel( _optPar, mode );

    ifail = _biosys->getComputeErrorFlag();
    if ( ifail != 0 )
    {
        return mat;
    }

    _linftyModel = _biosys->getLinftyModel();

    long k = 0;
    mat.zeros( f.nr(), _optPar.size() );

    for (Expression::ParamIterConst itPar = pBeg;
                                    itPar != pEnd; ++itPar)
    {
        ++k;

        Real w  = itPar->second;

        int  su = ( w < 0.0 ) ? -1 : 1;
        Real u  = std::max( std::max(std::fabs(w), ajmin), parScale[itPar->first] );

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
BioProcessor::computeJcf(std::string mode, int& ifail)
{
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    Expression::Param           parScale = computeParameterScales();

    Real   etadif = 1.0e-6;
    Real   etaini = 1.0e-6;
    Real   epdiff = std::sqrt(10.0 * EPMACH);
    Real   etamax = std::sqrt(epdiff);
    Real   etamin = epdiff*etamax;

    _linftyModel.clear();

    Matrix mat;
    Vector fu;
    Vector f = _biosys->computeModel( _optPar, mode );

    ifail = _biosys->getComputeErrorFlag();
    if ( ifail != 0 )
    {
        return mat;
    }

    _linftyModel = _biosys->getLinftyModel();

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
            Real u  = eta(k) * parScale[it->first] * su;

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
