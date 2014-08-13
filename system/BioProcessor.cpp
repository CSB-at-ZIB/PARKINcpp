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
    _method(method), _logstream(0), _iopt(),
    _curSpecies( _biosys->getSpecies() ),
    _optPar(), // _optIdx(),
    _speThres(), _parThres(),
    _trajMap(), _sensTraj(0),
    _sensiMat(), _sensiDcmp(),
    _nlscon(), _nlsconWk(),
    _parkin(), _parkinWk(),
    _idResult(),
    _niter(-1), _kappa(0.0)
{
    setLogStream( /* std::clog */ );

    BioSystem::StrIterConst sBeg = _curSpecies.begin();
    BioSystem::StrIterConst sEnd = _curSpecies.end();

    _speThres.clear();

    for (BioSystem::StrIterConst it = sBeg; it != sEnd; ++it)
    {
        _speThres[*it] = 0.0;
    }
}
//---------------------------------------------------------------------------
BioProcessor::~BioProcessor()
{
    delete _sensTraj;
}
//---------------------------------------------------------------------------
BioProcessor::BioProcessor(BioProcessor const& other) :
    _biosys(other._biosys), _biopar(other._biosys),
    _method(other._method), _logstream(0), _iopt(other._iopt),
    _curSpecies(other._curSpecies),
    _optPar(other._optPar), // _optIdx(other._optIdx),
    _speThres(other._speThres), _parThres(other._parThres),
    _trajMap(), _sensTraj(0),
    _sensiMat(other._sensiMat), _sensiDcmp(other._sensiDcmp),
    _nlscon(), _nlsconWk(),
    _parkin(), _parkinWk(),
    _idResult(other._idResult),
    _niter(other._niter), _kappa(other._kappa)
{
    setLogStream( /* std::clog */ );
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
BioProcessor::setLogStream(std::ostream& logstream)
{
    _logstream = &logstream;

    /*
    if ( _method == "parkin" )
    {
        _parkin.setLogStream(logstream);
    }
    else if ( _method == "nlscon" )
    {
        _nlscon.setLogStream(logstream);
    }
    */
}
//---------------------------------------------------------------------------
void
BioProcessor::setIOpt(IOpt const& iopt)
{
    _iopt = iopt;
}
//---------------------------------------------------------------------------
IOpt
BioProcessor::getIOpt()
{
    return _iopt;
}
//---------------------------------------------------------------------------
void
BioProcessor::setParameterConstraints(Vector const& itrans,
                                      Vector const& xlb, Vector const& xub)
{
    if ( itrans.nr() > 0 )
    {
        _iopt.transf = 1;
        _iopt.itrans = itrans;
        _nlsconWk.xlb = xlb;
        _nlsconWk.xub = xub;
    }
    else
    {
        _iopt.transf = 0;
    }
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
    Expression::ParamIterConst pBeg = par.begin();
    Expression::ParamIterConst pEnd = par.end();

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
BioProcessor::setCurrentSpeciesThres(Expression::Param const& par)
{
    Expression::ParamIterConst sBeg = par.begin();
    Expression::ParamIterConst sEnd = par.end();

    _speThres.clear();

    for (Expression::ParamIterConst it = sBeg; it != sEnd; ++it)
    {
        _speThres[it->first] = it->second;
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
    // bool                        lpos = _iopt.lpos;
    int                         transf = _iopt.transf;
    int                         jacgen = _iopt.jacgen;   // 1: variational eqn
                                                         // 2: num.diff.
                                                         // 3: num.diff.(with feedback)
    BioSystem::StrIterConst     sBeg = _curSpecies.begin();
    BioSystem::StrIterConst     sEnd = _curSpecies.end();
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    // TrajectoryMap               trajMap;
    // Vector                      itrans = _iopt.itrans;
    // Vector                      xlb = _nlsconWk.xlb;
    // Vector                      xub = _nlsconWk.xub;
    Matrix                      mat;
    int                         ifail = 0;


//    std::string tstr;
//
//    switch ( jacgen )
//    {
//        case 1 : tstr = "Variational equation";     break;
//        case 2 : tstr = "Num. diff.";               break;
//        case 3 : tstr = "Num. diff. w/ feedback";   break;
//        default: tstr = " N/A ";                    break;
//    }
//    std::cerr << std::endl;
//    std::cerr << "*** BioProcessor::computeSensitivityTrajectories() ***" << std::endl;
//    std::cerr << " jacgen = " << jacgen << "(" << tstr << ")" << std::endl;
//    std::cerr << "***" << std::endl;


    _trajMap.clear();
    _linftyModel.clear();

    _biosys -> computeModel( _optPar, "adaptive" );
    ifail = _biosys->getComputeErrorFlag();
    if ( ifail != 0 )
    {
        return _trajMap;
    }
    _linftyModel = _biosys->getLinftyModel();


    if ( jacgen == 1 )
    {
        mat = _biosys->computeJacobian( _optPar, "adaptive" );
        ifail = _biosys->getComputeErrorFlag();

        // delete _sensTraj;
        // _sensTraj = _biosys->getEvaluationTrajectories();

        // _linftyModel = _biosys->getLinftyModel();
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


    // if ( lpos == true )
    if ( (transf > 0) && (jacgen == 1) )
    {
        Vector dy;
        long   j = mat.nc();

        dy.zeros( j );
        j = 1;

        for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
        {
            Real tmp = _optPar[it->first];

            dy(j) = dertransf_p( tmp, j );

            ++j;
        }

        mat = mat * dy.diag();
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
    /// return _biosys->getOdeTrajectoryTimePoints(); // not quite working?!?
    return _biosys->getMeasurementTimePoints();
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
        Real yScale = speScale[*it];

        if ( yScale > 0.0 )
        {
            for (Expression::ParamIterConst itPar = pBeg;
                                            itPar != pEnd; ++itPar)
            {
                std::string s = *it + " / " + itPar->first;

                long T = _trajMap[s].size();
                Real pScale = parScale[itPar->first];

                scaledTrajMap[s].resize(T);

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
    // bool                lpos = _iopt.lpos;
    int                 transf = _iopt.transf;
    int                 jacgen = _iopt.jacgen;
    int                 ifail;
    Expression::Param   parScale = computeParameterScales();
    Expression::Param   speScale;  //  = computeSpeciesScales();
    unsigned            mcon = 0;
    unsigned            irankMax = 0;
 // Real                condMax = 1.0/EPMACH; // ( _biosys->getSolverRTol() );
    Real                condMax = 1.0/( _biosys->getSystemTol() );
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

//std::cerr << std::endl;
//std::cerr << "entering: ### BioProcessor::prepareDetailedSensitivities() ###\n";
//std::cerr << "        m = " << m << std::endl;
//std::cerr << "        q = " << q << std::endl;
//std::cerr << "        T = " << T << std::endl;
//std::cerr << "   jacgen = " << jacgen << std::endl;
////std::cerr << "     lpos = " << ((lpos == true) ? "true" : "false") << std::endl;
//std::cerr << "   transf = " << transf << std::endl;

    for (long l = 1; l <= q; ++l)
    {
        pw(l) = parScale[itPar->first];
        ++itPar;
    }
    // if ( lpos == true )
    if ( (transf > 0) && (jacgen == 1) )
    {
        itPar = _optPar.begin();

        for (long l = 1; l <= q; ++l)
        {
            Real tmp = _optPar[itPar->first];
            Real dy  = 1.0;

            dy = dertransf_p( tmp, l);

            pw(l) = parScale[itPar->first] * dy;
            ++itPar;
        }
    }

    ///

    _biosys->computeModel( _optPar, "adaptive" );

    ifail = _biosys->getComputeErrorFlag();

    if ( ifail != 0 )
    {
        return /* -1000 */ ifail;
    }

    _linftyModel.clear();
    _linftyModel = _biosys->getLinftyModel();

    speScale = computeSpeciesScales();
    itSpe = _curSpecies.begin();

    for (long k = 1; k <= m; ++k)
    {
        Real tmp = speScale[*itSpe++];
        one_fw(k) = (tmp > 0.0) ? (1.0/tmp) : 0.0;
    }

    ///

    if ( jacgen == 1 )
    {
        // if ( _sensTraj == 0 )
        // {
            _biosys->computeJacobian( _optPar, "adaptive" );

            ifail = _biosys->getComputeErrorFlag();

            if ( ifail != 0 )
            {
                return /* -1001 */ ifail;
            }

            delete _sensTraj;
            _sensTraj = _biosys->getEvaluationTrajectories();
        // }

        // speScale = computeSpeciesScales();
        // itSpe = _curSpecies.begin();
        // for (long k = 1; k <= m; ++k)
        // {
        //     Real tmp = speScale[*itSpe++];
        //     one_fw(k) = (tmp > 0.0) ? (1.0/tmp) : 0.0;
        // }

        for (long j = 1; j <= T; ++j)
        {
            Vector v = _sensTraj->eval( tp(j) );
            Matrix mat;
            long   n = m+1;     // small(!!) remainder: first entry is reserved
                                // (for the time variable...)

            mat.zeros(m,q);

/*
std::cerr << "\n";
std::cerr << "*******************************************" << std::endl;
std::cerr << "* BioProc::prepareDetailedSensitivities() *" << std::endl;
std::cerr << "*******************************************" << std::endl;
std::cerr << " T       = " << T << std::endl;
std::cerr << " v.nr()  = " << v.nr() << std::endl;
std::cerr << " n*(q+1) = " << n*(q+1) << std::endl;
std::cerr << " n       = " << n << std::endl;
std::cerr << " m       = " << m << std::endl;
std::cerr << " q       = " << q << std::endl;
std::cerr << " v.t()   = " << v.t() << std::endl;
std::cerr << "*******************************************" << std::endl;
*/

            if ( v.nr() == n*(q+1) )
            {
                for (long k = 1; k <= m; ++k)
                {
                    Real r = one_fw(k);

                    for (long l = 1; l <= q; ++l)
                    {
                        Real s = pw(l);

                        // mat(k,l) = r * v(++n) * s;
                        mat(k,l) = r * v( n + (n*(l-1) + k+1) ) * s;
                    }
                }
            }

            // mat = one_fw.diag() * mat * pw.diag();

            unsigned irank = irankMax;
            Real     cond = condMax;

            _sensiMat.push_back( mat );
            _sensiDcmp.push_back( mat.factorQRcon(mcon, irank, cond) );
        }
    }
    else if ( jacgen == 2 )
    {
        // _biosys->setEmptyMeasurementList();
        _biosys->setMeasurementTimePoints( tp );

        int ifail = 0;
        Matrix Jac = computeJac( "init", ifail );

        if ( (ifail != 0) || (Jac.nr() != m*T) )
        {
            return (ifail == 0) ? -1002 : ifail;
        }

        // speScale = computeSpeciesScales();
        // itSpe = _curSpecies.begin();
        // for (long k = 1; k <= m; ++k)
        // {
        //     Real tmp = speScale[*itSpe++];
        //     one_fw(k) = (tmp > 0.0) ? (1.0/tmp) : 0.0;
        // }

        for (long j = 1; j <= T; ++j)
        {
            Matrix mat;

            // mat.zeros(m,q);
            mat = Jac.rowm( (j-1)*m + 1, j*m );

            for (long k = 1; k <= m; ++k)
            {
                Real r = one_fw(k);

                for (long l = 1; l <= q; ++l)
                {
                    Real s = pw(l);

                    mat(k,l) *= (r * s);
                }
            }

            // mat = one_fw.diag() * mat * pw.diag();

            unsigned irank = irankMax;
            Real     cond = condMax;

            _sensiMat.push_back( mat );
            _sensiDcmp.push_back( mat.factorQRcon(mcon, irank, cond) );
        }
    }
    else if ( jacgen == 3 )
    {
        // _biosys->setEmptyMeasurementList();
        _biosys->setMeasurementTimePoints( tp );

        int ifail = 0;
        Matrix Jac;
//        if (transf > 0)  /// q'n'd hack for the time being: if transf, then feedbck off
//        {
//            Jac = computeJac( "init", ifail );
//        }
//        else
        {
            Jac = computeJcf( "init", ifail );
        }

        if ( (ifail != 0) || ( Jac.nr() != m*T) )
        {
            return (ifail == 0) ? -1003 : ifail;
        }

        // speScale = computeSpeciesScales();
        // itSpe = _curSpecies.begin();
        // for (long k = 1; k <= m; ++k)
        // {
        //     Real tmp = speScale[*itSpe++];
        //     one_fw(k) = (tmp > 0.0) ? (1.0/tmp) : 0.0;
        // }

        for (long j = 1; j <= T; ++j)
        {
            Matrix mat;

            // mat.zeros(m,q);
            mat = Jac.rowm( (j-1)*m + 1, j*m );

            for (long k = 1; k <= m; ++k)
            {
                Real r = one_fw(k);

                for (long l = 1; l <= q; ++l)
                {
                    Real s = pw(l);

                    mat(k,l) *= (r * s);
                }
            }

            // mat = one_fw.diag() * mat * pw.diag();

            unsigned irank = irankMax;
            Real     cond = condMax;

            _sensiMat.push_back( mat );
            _sensiDcmp.push_back( mat.factorQRcon(mcon, irank, cond) );
        }
    }
    else
    {
        return -1004;
    }

//std::cerr << std::endl;
//std::cerr << " pw = " << std::endl;
//std::cerr << pw.t() << std::endl;
//std::cerr << " one_fw = " << std::endl;
//std::cerr << one_fw.t() << std::endl;
//std::cerr << "### BioProcessor::prepareDetailedSensitivities() ###\n";

    return 0;
}
//---------------------------------------------------------------------------
BioProcessor::MatrixList
BioProcessor::getSensitivityMatrices()
{
    return _sensiMat;
}
//---------------------------------------------------------------------------
BioProcessor::QRconDecompList
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
    // Expression::Param           parScal;
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    BioSystem::Parameter        pname;
    Vector                      x, xscal;
    Vector                      fobs, fscal;

    fobs = _biosys->getMeasurements();
    fscal = _biosys->getMeasurementWeights();

    m = fobs.nr();

    if ( _iopt.zerodat == true )
    {
        fobs.zeros(m);
        fscal.ones(m);
    }

/*
std::cerr << std::endl;
std::cerr << "***" << std::endl;
std::cerr << "*** BioProcessor::identifyParameters():" << std::endl;
std::cerr << "***" << std::endl;
std::cerr << "  fobs.t() = " << std::endl;
std::cerr << fobs.t() << std::endl;
std::cerr << "  fscal.t() = " << std::endl;
std::cerr << fscal.t() << std::endl;
*/

    j = 1;
    pname.clear();
    x.zeros( _optPar.size() );
    xscal.zeros( _optPar.size() );

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
    {
        std::string s = it->first;

        pname.push_back( s );
            x(j) = _optPar[s];      // it->second;
        xscal(j) = _parThres[s];

        ++j;
    }

    /*
    j = 0;
    // parScal = computeParameterScales();
    xscal.zeros( _optPar.size() );

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
    {
        std::string s = it->first;

        xscal(++j) = _parThres[s]; // it->second;
    }
    */

    // _biopar = BioPAR( _biosys, pname );
    _biopar.setCurrentParameter( pname );
    _biopar.setEvalMode( _iopt.zerodat );


    if ( _method == "parkin" )
    {
        _parkinWk.itmax = _iopt.itmax;
        // _parkinWk.cond = 1.0 / (xtol * 1.0e+1); // ( _biosys->getSolverRTol() );
        _parkinWk.cond = 1.0 / ( _biosys->getSystemTol() );

        _parkin.setLogStream( *_logstream );

        _parkin.setProblem( &_biopar );
        rc = _parkin.initialise( m,
                                 x, xscal,
                                 fobs, fscal,
                                 xtol, _iopt,
                                 _parkinWk
                               );
        if ( rc != 0 ) return rc;

        rc = _parkin.run();

        // _parkinWk = _parkin.getWk();
        _niter = _parkin.getWk().iter;
        _kappa = _parkin.getWk().kappa;
        _idResult = _parkin.getSolution();
    }
    else if ( _method == "nlscon" )
    {

        *_logstream << std::endl;
        *_logstream << "          --------------------------------" << std::endl;
        *_logstream << "          Internal scaling of measurements" << std::endl;
        *_logstream << "          --------------------------------" << std::endl;
        *_logstream << std::endl;

        char line[256]; *line = '\0';

        for (unsigned j = 1; j <= m; ++j)
        {
            std::sprintf(line, "%s %7d %13.6e", line, j, fscal(j) );

            if ( j%2 == 0 )
            {
                *_logstream << line << std::endl;
                *line = '\0';
            }
        }

        if ( *line != '\0' )
        {
            *_logstream << line << std::endl;
        }

        *_logstream << std::endl;


        // _nlsconWk.xlb = _xlb;
        // _nlsconWk.xub = _xub;
        _nlsconWk.nitmax = _iopt.itmax;
        // _nlsconWk.cond    = 1.0 / (xtol * 1.0e+1);
        //_nlsconWk.fcstart = 1.0e-1;
        //_nlsconWk.fcmin   = 1.0e-2;

        _nlsconWk.cond = 1.0 / ( _biosys->getSystemTol() );
        // _nlsconWk.cond   = 1.0 / ( _biosys->getSolverRTol() );
        // _nlsconWk.fcmin  = 1.0e-4;


        _nlscon.setLogStream( *_logstream );

        _nlscon.setProblem( &_biopar );
        rc = _nlscon.initialise( m,
                                 x, xscal,
                                 fobs, fscal,
                                 xtol, _iopt,
                                 _nlsconWk
                               );
        if ( rc != 0 ) return rc;


        rc = _nlscon.run();
        // _nlsconWk = _nlscon.getWk();
        _idResult = _nlscon.getSolution();


        if ( (rc == 0) || (rc == 1) || (rc == 2) || (rc == 3) )
        {
            _nlscon.analyse();
        }

        ///

        *_logstream << std::endl;
        *_logstream << std::endl;
        *_logstream << "         ------------------------------------------" << std::endl;
        *_logstream << "         NLSCON: Parameter Map ( no. ---> ID-name )" << std::endl;
        *_logstream << "         ------------------------------------------" << std::endl;
        for (unsigned j = 0; j < pname.size(); ++j)
        {
            *_logstream << std::setw(8) << std::right << j+1 << " : "
                                        << std::left << pname[j]
                                        << std::endl;
        }
        *_logstream << std::endl;


        std::vector<Vector> piter = _nlscon.getSolutionIter();
        _niter = _nlscon.getWk().niter;
        _kappa = _nlscon.getWk().skap;

        *_logstream << std::endl;
        *_logstream << "         ------------------------------------------" << std::endl;
        *_logstream << "         NLSCON: Inverse Problem Solution Iteration" << std::endl;
        *_logstream << "         ------------------------------------------" << std::endl;
        for (long j = 0; j <= _niter; ++j)
        {
            *_logstream << "it = " << j << "\n" << piter[j].t();
        }
        *_logstream << std::endl;

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
int
BioProcessor::getNIter()
{
    return _niter;
}
//---------------------------------------------------------------------------
Real
BioProcessor::getKappa()
{
    return _kappa;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Expression::Param
BioProcessor::computeParameterScales()
{
    // int                         transf = _iopt.transf;
    Vector                      itrans = _iopt.itrans;
    Expression::ParamIterConst  pBeg   = _optPar.begin();
    Expression::ParamIterConst  pEnd   = _optPar.end();
    Expression::Param           parScale;
    long                        k;

    parScale.clear();
    k = 1;
    if ( itrans.nr() < 1 )
    {
        itrans.zeros( _optPar.size() );
    }

    for (Expression::ParamIterConst itPar = pBeg;
                                    itPar != pEnd; ++itPar)
    {
        std::string s = itPar->first;
        Real optPar   = _optPar[s];
        Real parThres = _parThres[s];

        /*
        if ( transf > 0 )
        {
            optPar = transform_p( optPar, k );
            parThres = transform_p( parThres, k );

            ++k;
        }
        */

        parScale[s] =
            std::max(   std::fabs( optPar ) , parThres );

        if ( itrans(k) > 0.0 )  // switch off scaling if any transformation is used
        {
            parScale[s] = 1.0;
        }

        ++k;
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
BioProcessor::computeJac(std::string const& mode, int& ifail)
{
    int                        transf = _iopt.transf;
    Expression::ParamIterConst pBeg = _optPar.begin();
    Expression::ParamIterConst pEnd = _optPar.end();
    Expression::Param          parScale = computeParameterScales();

    Real   ajmin = SMALL;
    Real   ajdelta = std::sqrt(10.0*EPMACH);

    Matrix mat;
    Vector fh;
    Vector f = _biosys->computeModel( _optPar, mode );

    ifail = _biosys->getComputeErrorFlag();
    if ( ifail != 0 )
    {
        return mat;
    }

    // _linftyModel.clear();
    // _linftyModel = _biosys->getLinftyModel();


    long k = 1;
    mat.zeros( f.nr(), _optPar.size() );

//std::cerr << std::endl;
//std::cerr << "### BioProcessor::computeJac(...) ###" << std::endl;
//std::cerr << std::endl;
//std::cerr << "  Parameter:" << std::endl;

    for (Expression::ParamIterConst itPar = pBeg;
                                    itPar != pEnd; ++itPar, ++k)
    {
        std::string s = itPar->first;
        Real w        = _optPar[s];
        Real wSave    = w;
        if ( transf > 0 )
        {
            w = transform_p(w, k);
        }

        int  su = ( w < 0.0 ) ? -1 : 1;
        Real u  = std::max( std::max(std::fabs(w), ajmin), parScale[s] );

        u *= (ajdelta * su);

        _optPar[s] = w + u;

        if ( transf > 0 )
        {
            _optPar[s] = backtrans_p(w + u, k);
        }

        fh = _biosys->computeModel( _optPar );

        _optPar[s] = wSave;
        mat.set_colm(k) = (1.0/u) * Matrix(fh - f);

//std::cerr << std::right << std::setw(4) << k << ": ";
//std::cerr << std::left << std::setw(15) << itPar->first << " ";
//std::cerr << mat.colm(k).t();
//std::cerr << std::setw(22) << " " << "w=" << w << ",  u=" << u << ",  w+u=" << w+u << std::endl;
//std::cerr << std::setw(22) << " " << fh.t();
//std::cerr << std::setw(22) << " " << f.t();
//std::cerr << std::setw(22) << " " << Matrix(fh - f).t();
    }

//std::cerr << std::endl;
//std::cerr << "###" << std::endl;

    return mat;
}
//---------------------------------------------------------------------------
Matrix
BioProcessor::computeJcf(std::string const& mode, int& ifail)
{
    int                         transf = _iopt.transf;
    Expression::ParamIterConst  pBeg = _optPar.begin();
    Expression::ParamIterConst  pEnd = _optPar.end();
    Expression::Param           parScale = computeParameterScales();

    Real   etadif = 1.0e-6;
    Real   etaini = 1.0e-6;
    Real   epdiff = std::sqrt(10.0 * EPMACH);
    Real   etamax = std::sqrt(epdiff);
    Real   etamin = epdiff*etamax;

    Matrix mat;
    Vector fu;
    Vector f = _biosys->computeModel( _optPar, mode );

    ifail = _biosys->getComputeErrorFlag();
    if ( ifail != 0 )
    {
        return mat;
    }

    // _linftyModel.clear();
    // _linftyModel = _biosys->getLinftyModel();

    long m = f.nr();
    long n = _optPar.size();
    // Vector v;  v.ones( n );
    Vector eta;

    eta.zeros( n );
    for (long k = 1; k <= n; ++k)
    {
        eta(k) = etaini;
    }

    mat.zeros(m,n);
    long k = 1;

//std::cerr << std::endl;
//std::cerr << "### BioProcessor::computeJcf(...) ###" << std::endl;
//std::cerr << std::endl;
//std::cerr << "  Parameter:" << std::endl;

    for (Expression::ParamIterConst it = pBeg;
                                    it != pEnd; ++it, ++k)
    {
        std::string s = it->first;
        Real sumd, hg, fhj;
        Real w = 0.0, u = 0.0;
        Real wSave = 0.0;
        int  su = 0;
        bool is = false;
        bool qexit = false;
        bool qfine = false;

        while ( !qfine )
        {
            w  = _optPar[s];
            wSave = w;
            if ( transf > 0 )
            {
                w = transform_p(w, k);
            }
            su = ( w < 0.0 ) ? -1 : 1;
            u  = eta(k) * parScale[s] * su;

            _optPar[s] = w + u;
            if ( transf > 0 )
            {
                _optPar[s] = backtrans_p(w + u, k);
            }

            fu = _biosys->computeModel( _optPar );

            _optPar[s] = wSave;

            ifail = _biosys->getComputeErrorFlag();

            if ( ifail != 0 )
            {
                qexit = true; break;
            }

            sumd = 0.0;
            for (long j = 1; j <= m; ++j)
            {
                hg = std::max( std::fabs( f(j) ), std::fabs( fu(j) ) );
                fhj = fu(j) - f(j);
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

//std::cerr << std::right << std::setw(4) << k << ": ";
//std::cerr << std::left << std::setw(15) << it->first << " ";
//std::cerr << mat.colm(k).t();
//std::cerr << std::setw(22) << " " << "eta(k)=" << eta(k) << std::endl;
//std::cerr << std::setw(22) << " " << "w=" << w << ",  u=" << u << ",  w+u=" << w+u << std::endl;
//std::cerr << std::setw(22) << " " << fu.t();
//std::cerr << std::setw(22) << " " << f.t();
//std::cerr << std::setw(22) << " " << Matrix(fu - f).t();
    }

//std::cerr << std::endl;
//std::cerr << "###" << std::endl;

    return mat;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Real
BioProcessor::dertransf_p(Real p, long k)  // dy = phi'(u) if u = phi^(-1)(xparam)
{
    Real utmp;
    Real dy     =  1.0;
    Real itrans = _iopt.itrans(k);
    Real xlb    = _nlsconWk.xlb(k);
    Real xub    = _nlsconWk.xub(k);

    // dy(j) = 1.0;

    if ( itrans == 1.0 )
    {
        // utmp = std::log( p );
        // dy   = std::exp( utmp );

        dy = std::max( p, 0.0 );
    }
    else if ( itrans == 2.0 )
    {
        utmp = 1.0 - xlb + p;
        utmp  = std::sqrt( -1.0 + utmp*utmp );

        dy = utmp / std::sqrt( 1.0 + utmp*utmp );
    }
    else if ( itrans == 3.0 )
    {
        utmp = 1.0 + xub - p;
        utmp = std::sqrt( -1.0 + utmp*utmp );

        dy = - utmp / std::sqrt( 1.0 + utmp*utmp );
    }
    else if ( itrans == 4.0 )
    {
        utmp = (p - xlb) / (xub - xlb);
        utmp = std::asin( -1.0 + 2.0*utmp );

        dy = 0.5 * (xub - xlb) * std::cos( utmp );
    }

    return dy;
}
//---------------------------------------------------------------------------
Real
BioProcessor::transform_p(Real p, long k)   // u = phi^(-1)(xparm)
{
    Real ptmp;
    Real u      = p;
    Real itrans = _iopt.itrans(k);
    Real xlb    = _nlsconWk.xlb(k);
    Real xub    = _nlsconWk.xub(k);

    if ( itrans == 1.0 )
    {
        u = -1.0e38;
        if ( p > 0.0 )
        {
            u = std::log( p );
        }
    }
    else if ( itrans == 2.0 )
    {
        ptmp = 1.0 - xlb + p;
        u = std::sqrt( -1.0 + ptmp*ptmp );
    }
    else if ( itrans == 3.0 )
    {
        ptmp = 1.0 + xub - p;
        u = std::sqrt( -1.0 + ptmp*ptmp );
    }
    else if ( itrans == 4.0 )
    {
        ptmp = (p - xlb) / (xub - xlb);
        u = std::asin( -1.0 + 2.0 * ptmp );
    }

    return u;
}
//---------------------------------------------------------------------------
Real
BioProcessor::backtrans_p(Real u, long k)   // xparam = phi(u)
{
    Real x      = u;
    Real itrans = _iopt.itrans(k);
    Real xlb    = _nlsconWk.xlb(k);
    Real xub    = _nlsconWk.xub(k);

    if ( itrans == 1.0 )
    {
        x = std::exp( u );                                     //   0 < x
    }
    else if ( itrans == 2.0 )
    {
        x = -1.0 + xlb + std::sqrt( 1.0 + u*u );               // xlb <= x
    }
    else if ( itrans == 3.0 )
    {
        x =  1.0 + xub - std::sqrt( 1.0 + u*u );               //        x <= xub
    }
    else if ( itrans == 4.0 )
    {
        x = xlb + 0.5*(xub - xlb) * ( 1.0 + std::sin( u ) );   // xlb <= x <= xub
    }

    return x;
}
//---------------------------------------------------------------------------
