// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-08 td
// last changed:
//

#include <addpkg/dlib/time_this.h>
#include "BioSystem.h"
#include "BioSystemODE.h"
#include "BioSystemVAR.h"
#include "linalg/QRconDecomp.h"
#include "odelib/LIMEX_A.h"
#include "odelib/DOP853.h"
#include "odelib/METAN_A.h"

using namespace PARKIN;
// / *
// BioSystem*  BioSystemWrapper::_obj = 0;
// BioRHS      BioSystemWrapper::_ode = BioRHS();
// * /
//---------------------------------------------------------------------------
/// c'tor
BioSystem::BioSystem( Real tStart, Real tEnd ) :
    _adm(),
    _ode(), _iniCond(),
    _iniPar(), _sysPar(), _optPar(),
    _parValue(0),
    _linftyModel(),
    _tpMeas(), _measData(),
    _synData(),
    _systemOde( new BioSystemODE() ),
    _variationalOde( new BioSystemVAR() ),
    //$$$ _odeSolver( new DOP853() ),
    _odeSolver( new LIMEX_A() ),
    _odeErrorFlag(0),
    _totmeasData(0),
    _tInterval(),
    _systemTol( 10.0*EPMACH )
{
    _linftyModel.clear();
    _tpMeas.clear();
    _measData.clear();
    _tInterval.clear();
    _tInterval.push_back( tStart );
    _tInterval.push_back( tEnd );
    _iniCond.clear();
    _iniCond.push_back( BioRHS() );
    _iniCond.push_back( BioRHS() );
    dynamic_cast<BioSystemODE*>(_systemOde) -> setObj(*this);
    dynamic_cast<BioSystemVAR*>(_variationalOde) -> setObj(*this);
}
//---------------------------------------------------------------------------
/// c'tor
BioSystem::BioSystem( Vector const&  tInterval ) :
    _adm(),
    _ode(), _iniCond(),
    _iniPar(), _sysPar(), _optPar(),
    _parValue(0),
    _linftyModel(),
    _tpMeas(), _measData(),
    _synData(),
    _systemOde( new BioSystemODE() ),
    _variationalOde( new BioSystemVAR() ),
    //$$$ _odeSolver( new DOP853() ),
    _odeSolver( new LIMEX_A() ),
    _odeErrorFlag(0),
    _totmeasData(0),
    _tInterval(),
    _systemTol( 10.0*EPMACH )
{
    _linftyModel.clear();
    _tpMeas.clear();
    _measData.clear();
    _tInterval.clear();
    _iniCond.clear();
    for (long j = 1; j <= tInterval.nr(); ++j)
    {
        _tInterval.push_back( tInterval(j) );
        _iniCond.push_back( BioRHS() );
    }
    dynamic_cast<BioSystemODE*>(_systemOde) -> setObj(*this);
    dynamic_cast<BioSystemVAR*>(_variationalOde) -> setObj(*this);
}
//---------------------------------------------------------------------------
/// c'tor
BioSystem::BioSystem( ExpressionMap const&      eMap,
                      ODESolver::Grid const&    tPoints,
                      MeasurementList const&    meas,
                      ODESolver::Grid const&    tInterval
                    ) :
    _adm(),
    _ode(eMap), _iniCond(),
    _iniPar(), _sysPar(), _optPar(),
    _parValue(0),
    _linftyModel(),
    _tpMeas(tPoints), _measData(meas),
    _synData(),
    _systemOde( new BioSystemODE() ),
    _variationalOde( new BioSystemVAR() ),
    //$$$ _odeSolver( new DOP853() ),
    _odeSolver( new LIMEX_A() ),
    _odeErrorFlag(0),
    _tInterval(tInterval),
    _systemTol( 10.0*EPMACH )
{
    long n = 0;

    _linftyModel.clear();

    for (unsigned j = 0; j < _measData.size(); ++j)
        n += _measData[j].size();

    _totmeasData = n;

    setIdentityEvents();

    dynamic_cast<BioSystemODE*>(_systemOde) -> setObj(*this);
    dynamic_cast<BioSystemVAR*>(_variationalOde) -> setObj(*this);
}
//---------------------------------------------------------------------------
/// d'tor
BioSystem::~BioSystem()
{
    delete[] _parValue;
    delete   _systemOde;
    delete   _variationalOde;
    delete   _odeSolver;
}
//---------------------------------------------------------------------------
/// copy c'tor
 BioSystem::BioSystem(BioSystem const& s)
 {
     if (this != &s)
     {
         _adm     = s._adm;
         _ode     = s._ode;
         _iniCond = s._iniCond;

         _iniPar = s._iniPar;
         _sysPar = s._sysPar;
         _optPar = s._optPar;

         long q = s._ode.getParameters().size();

         _parValue = new double[q];

         _linftyModel = s._linftyModel;
         _tpMeas      = s._tpMeas;
         _measData    = s._measData;
         _synData.clear();      // = s._synData;
         _jacobian.clear();     // = s._jacobian;
         _jac = Matrix();

         _systemOde = new BioSystemODE();
         _variationalOde = new BioSystemVAR();

         //$$$ _odeSolver = new DOP853();
         _odeSolver = new LIMEX_A();
         //_odeSolver = s._odeSolver->clone();
         _odeSolver->setDebugFlag( s._odeSolver->getDebugFlag() );
         _odeSolver->setInterpolationFlag( s._odeSolver->getInterpolationFlag() );
         _odeSolver->setRTol( s._odeSolver->getRTol() );
         _odeSolver->setATol( s._odeSolver->getATol() );
         _odeErrorFlag = 0;

         _totmeasData = s._totmeasData;
         _tInterval   = s._tInterval;
         _systemTol   = s._systemTol;

         dynamic_cast<BioSystemODE*>(_systemOde) -> setObj(*this);
         dynamic_cast<BioSystemVAR*>(_variationalOde) -> setObj(*this);
     }
 }
//---------------------------------------------------------------------------
/// assignment
// BioSystem const&
// BioSystem::operator= (BioSystem const& s)
// {
//     if (this != &s)
//     {
//         _odeExpr = s._odeExpr;
//         _sysPar = s._sysPar;
//         _iniPar = s._iniPar;
//     }
//     return *this;
// }
//---------------------------------------------------------------------------
////---------------------------------------------------------------------------
//Vector
//BioSystem::fcn(Vector const& x, int& ifail)
//{
//    ifail = 0;
//    return computeModel(x);
//}
////---------------------------------------------------------------------------
//Matrix
//BioSystem::jac(Vector const& x, int& ifail)
//{
//    ifail = 0;
//    return computeJacobian(x);
//}
////---------------------------------------------------------------------------
void
BioSystem::setSolverDebugFlag(int flag)
{
    _odeSolver -> setDebugFlag(flag);
}
//---------------------------------------------------------------------------
void
BioSystem::setSolverInterpolationFlag(int cubint)
{
    _odeSolver -> setInterpolationFlag(cubint);
}
//---------------------------------------------------------------------------
ODESolverId
BioSystem::getSolverId() const
{
    return _odeSolver -> getId();
}
//---------------------------------------------------------------------------
void
BioSystem::setSolver(ODESolverId solverid)
{
    ODESolver* newSolver = 0;

    switch( solverid )
    {
        case ODE_SOLVER_METAN_A:
                // newSolver = new DOP853( *_odeSolver );
                newSolver = new METAN_A();
                break;

        case ODE_SOLVER_DOP853:
                // newSolver = new DOP853( *_odeSolver );
                newSolver = new DOP853();
                break;

        case ODE_SOLVER_DEFAULT:
        case ODE_SOLVER_LIMEX_A:
                // newSolver = new LIMEX_A( *_odeSolver );
                newSolver = new LIMEX_A();
                break;

        default:
                break;
    }

    if ( newSolver != 0 )
    {
        newSolver -> setRTol( _odeSolver -> getRTol() );
        newSolver -> setATol( _odeSolver -> getATol() );
        newSolver -> setIniStep( _odeSolver -> getIniStep() );
        newSolver -> setMaxStep( _odeSolver -> getMaxStep() );
        newSolver -> setDebugFlag( _odeSolver -> getDebugFlag() );
        newSolver -> setInterpolationFlag( _odeSolver -> getInterpolationFlag() );

        delete _odeSolver;

        _odeSolver = newSolver;
        _odeErrorFlag = 0;
    }
}
//---------------------------------------------------------------------------
Real
BioSystem::getSolverRTol() const
{
    return _odeSolver -> getRTol();
}
//---------------------------------------------------------------------------
void
BioSystem::setSolverRTol(Real tol)
{
    _odeSolver -> setRTol(tol);
}
//---------------------------------------------------------------------------
Real
BioSystem::getSolverATol() const
{
    return _odeSolver -> getATol();
}
//---------------------------------------------------------------------------
void
BioSystem::setSolverATol(Real tol)
{
    _odeSolver -> setATol(tol);
}
//---------------------------------------------------------------------------
Real
BioSystem::getSolverIniStep() const
{
    return _odeSolver -> getIniStep();
}
//---------------------------------------------------------------------------
void
BioSystem::setSolverIniStep(Real inistep)
{
    _odeSolver -> setIniStep(inistep);
}
//---------------------------------------------------------------------------
Real
BioSystem::getSolverMaxStep() const
{
    return _odeSolver -> getMaxStep();
}
//---------------------------------------------------------------------------
void
BioSystem::setSolverMaxStep(Real maxstep)
{
    _odeSolver -> setMaxStep(maxstep);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
std::string
BioSystem::getSolverMessage(int rc) const
{
    return _odeSolver -> getErrorMessage(rc);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Real
BioSystem::getSystemTol() const
{
    return _systemTol;
}
//---------------------------------------------------------------------------
void
BioSystem::setSystemTol(Real tol)
{
    _systemTol = tol;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Expression::Param&
BioSystem::getSysPar()
{
    return _sysPar;
}
//---------------------------------------------------------------------------
Expression::Param&
BioSystem::getOptPar()
{
    return _optPar;
}
//---------------------------------------------------------------------------
Expression::Param&
BioSystem::getLinftyModel()
{
    return _linftyModel;
}
//---------------------------------------------------------------------------
BioSystem::Species
BioSystem::getSpecies() const
{
    BioSystem::Species ret;
    BioSystem::Species species = _ode.getSpecies();
    StrIterConst       sBeg = species.begin();
    StrIterConst       sEnd = species.end();

    ret.clear();

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        ret.push_back( *it );
    }

    return ret;
}
//---------------------------------------------------------------------------
void
BioSystem::resetSpecies(Species const& species)
{
    long Tb = _tInterval.size();

    _ode.resetSpecies(species);

    for (long j = 0; j < Tb; ++j) _iniCond[j].resetSpecies(species);

    // _measData.clear(); _totmeasData = 0;
}
//---------------------------------------------------------------------------
BioSystem::Parameter
BioSystem::getParameters() const
{
    BioSystem::Parameter ret;
    BioSystem::Parameter param = _ode.getParameters();
    StrIterConst         pBeg = param.begin();
    StrIterConst         pEnd = param.end();

    ret.clear();

    for (StrIterConst it = pBeg; it != pEnd; ++it)
    {
        ret.push_back( *it );
    }

    return ret;
}
//---------------------------------------------------------------------------
void
BioSystem::setParameters(Parameter const& parameter)
{
    _ode.setParameters(parameter);

    long Tb = _tInterval.size();
    long q = _ode.getParameters().size();

    for (long j = 0; j < Tb; ++j) _iniCond[j].setParameters(parameter);

    delete[] _parValue;
    _parValue = new double[q];
}
//---------------------------------------------------------------------------
Vector
BioSystem::getBreakpoints() const
{
    return Vector( _tInterval );
}
//---------------------------------------------------------------------------
void
BioSystem::setBreakpoints(Vector const& tInterval)
{
    if ( 1 < tInterval.nr() )
    {
        _tInterval.clear();

        for (long j = 1; j <= tInterval.nr(); ++j)
        {
            _tInterval.push_back( tInterval(j) );
        }

        setIdentityEvents();
    }
}
//---------------------------------------------------------------------------
BioSystem::ExpressionMap
BioSystem::getEvent(long j) const
{
    long Tb = _tInterval.size();

    if ( (0 <= j) && (j < Tb) ) return _iniCond[j].getRHS();

    return ExpressionMap();
}
//---------------------------------------------------------------------------
void
BioSystem::setEvent(long j, ExpressionMap const& emap)
{
    long Tb = _tInterval.size();

    if ( (0 <= j) && (j < Tb) )
    {
        _iniCond[j].setRHS(emap);
    }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void
BioSystem::setMedicationList(MedicationList medList)
{
    int             nAdm    = medList.size();
    Species const&  species = _ode.getSpecies();
    StrIterConst    sBeg    = species.begin();
    StrIterConst    sEnd    = species.end();

    for (int j = 0; j < nAdm; ++j)
    {
        StrIterConst it = std::find(sBeg, sEnd, medList[j].adm_comp);

        medList[j].index = 0;

        if ( it != sEnd )
        {
            medList[j].index = (it - sBeg);  //trick to get index into vec!
            medList[j].index += 1;           //index 0 reserved for time t
        }
    }

    _adm.setMedicationList(medList);
}
//---------------------------------------------------------------------------
void
BioSystem::addMedication(Medication med)
{
    Species const&  species = _ode.getSpecies();
    StrIterConst    sBeg    = species.begin();
    StrIterConst    sEnd    = species.end();

    StrIterConst it = std::find(sBeg, sEnd, med.adm_comp);

    med.index = 0;

    if ( it != sEnd )
    {
        med.index = (it - sBeg);  // neat trick to get index into vector!!
        med.index += 1;           // index 0 reserved for time t
    }

    _adm.addMedication(med);
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void
BioSystem::setODESystem(ExpressionMap const& eMap)
{
    _ode.setRHS(eMap);

    setEmptyMeasurementList();

    setIdentityEvents();

    // BioSystemWrapper::setObj(*this);
}
//---------------------------------------------------------------------------
void
BioSystem::setODESystem(ExpressionMap const& eMap, ExprTypeMap const& tMap)
{
    _ode.setRHS(eMap);
    _ode.setRHSType(tMap);

    setEmptyMeasurementList();

    setIdentityEvents();

    // BioSystemWrapper::setObj(*this);
}
//---------------------------------------------------------------------------
void
BioSystem::setODETypes(ExprTypeMap const& tMap)
{
    _ode.setRHSType(tMap);
}
//---------------------------------------------------------------------------
void
BioSystem::setIdentityEvents()
{
    Species const&      species = _ode.getSpecies();
    Parameter const&    parameter = _ode.getParameters();
    BioRHS              rhs = BioRHS(species);

    rhs.setParameters(parameter);
    _iniCond.clear();

    for (long j = 0; j < (long)_tInterval.size(); ++j)
    {
        // BioRHS rhs = BioRHS(species);
        // rhs.setParameters(parameter);

        _iniCond.push_back( rhs );
    }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
BioSystem::MeasurementList const&
BioSystem::getMeasurementList()
{
    return _measData;
}
//---------------------------------------------------------------------------
BioSystem::MeasurementList const&
BioSystem::getSimTrajectoryList()
{
    return _synData;
}
//---------------------------------------------------------------------------
void
BioSystem::setMeasurementList(
                                Vector const& tp,
                                MeasurementList const& meas
                             )
{
    long T = meas.size();
    long n = 0;

    if ( tp.nr() == T )
    {
        for (long j = 0; j < T; ++j) n += meas[j].size();
        _totmeasData = n;
        _measData = meas;
        _tpMeas.clear();
        for (long j = 1; j <= T; ++j) _tpMeas.push_back( tp(j) );
    }
}
//---------------------------------------------------------------------------
void
BioSystem::setMeasurementList(MeasurementList const& meas)
{
    long T = meas.size();
    long n = 0;
    for (long j = 0; j < T; ++j) n += meas[j].size();
    _totmeasData = n;
    _measData = meas;
    _tpMeas.clear();
    for (long j = 0; j < T; ++j) _tpMeas.push_back(j);
}
//---------------------------------------------------------------------------
void
BioSystem::setEmptyMeasurementList()
{
    long           T = _tpMeas.size();
    Species const& species = _ode.getSpecies();
    StrIterConst   sBeg = species.begin();
    StrIterConst   sEnd = species.end();

    if ( T <= 0 )
    {
        T = 100;
        _tpMeas.clear();
        for (long j = 0; j < T; ++j) _tpMeas.push_back(j+1);
    }

    _measData.clear();

    for (long tp = 0; tp < T; ++tp)
    {
        _measData.push_back( MeasurementPoint() );

        for (StrIterConst it = sBeg; it != sEnd; ++it)
        {
            _measData[tp][*it] = std::make_pair<Real,Real>(0.0, GREAT);
        }
    }

    _totmeasData = T*species.size();
}
//---------------------------------------------------------------------------
Vector
BioSystem::getMeasurements() const
{
    long                T = _tpMeas.size();
    std::vector<Real>   meas;

    meas.clear();

    for (long tp = 0; tp < T; ++tp)
    {
        MeasurementPoint    mp = _measData[tp];

        MeasIterConst       mBeg = mp.begin();
        MeasIterConst       mEnd = mp.end();

        for (MeasIterConst it = mBeg; it != mEnd; ++it)
        {
            std::string spec = it->first;

            meas.push_back( ( mp[spec] ).first );
        }
    }

    return Vector( meas );
}
//---------------------------------------------------------------------------
Vector
BioSystem::getMeasurementWeights() const
{
    long                T = _tpMeas.size();
    std::vector<Real>   weights;

    weights.clear();

    for (long tp = 0; tp < T; ++tp)
    {
        MeasurementPoint    mp = _measData[tp];

        MeasIterConst       mBeg = mp.begin();
        MeasIterConst       mEnd = mp.end();

        for (MeasIterConst it = mBeg; it != mEnd; ++it)
        {
            std::string spec = it->first;

            weights.push_back( ( mp[spec] ).second );
        }
    }

    return Vector( weights );
}
//---------------------------------------------------------------------------
//void
//BioSystem::setMeasurementTimePoints(ODESolver::Grid const& tp)
//{
//    _tpMeas = tp;
//}
//---------------------------------------------------------------------------
void
BioSystem::setMeasurementTimePoints(Vector const& tp)
{
    _tpMeas.clear();
    for (long j = 1; j <= tp.nr(); ++j) _tpMeas.push_back( tp(j) );

    setEmptyMeasurementList();
}
//---------------------------------------------------------------------------
Vector
BioSystem::getMeasurementTimePoints() const
{
    return Vector( _tpMeas );
}
//---------------------------------------------------------------------------
Vector //ODESolver::Grid&
BioSystem::getOdeTrajectoryTimePoints() const
{
    //$$$ return Vector( dynamic_cast<DOP853*>(_odeSolver)->getSolutionGridPoints() );
    return Vector( _odeSolver -> getSolutionGridPoints() );
}
//---------------------------------------------------------------------------
Real
BioSystem::getParamValue(std::string const& name) // const
{
    return _sysPar[name];
}
//---------------------------------------------------------------------------
void
BioSystem::setParamValue(std::string const& name, Real value)
{
    _sysPar[name] = value;
}
//---------------------------------------------------------------------------
void
BioSystem::setParamValues(Expression::Param const& par)
{
    Expression::ParamIterConst pBeg = par.begin();
    Expression::ParamIterConst pEnd = par.end();

    for (Expression::ParamIterConst it = pBeg; it != pEnd; ++it)
    {
        _sysPar[it->first] = it->second;
    }
}
//---------------------------------------------------------------------------
Real
BioSystem::getInitialValue(std::string const& name) // const
{
    return _iniPar[name];
}
//---------------------------------------------------------------------------
void
BioSystem::setInitialValue(std::string const& name, Real value)
{
    _iniPar[name] = value;
}
//---------------------------------------------------------------------------
void
BioSystem::setInitialValues(Expression::Param const& inipar)
{
    Expression::ParamIterConst iBeg = inipar.begin();
    Expression::ParamIterConst iEnd = inipar.end();

    for (Expression::ParamIterConst it = iBeg; it != iEnd; ++it)
    {
        _iniPar[it->first] = it->second;
    }
}
//---------------------------------------------------------------------------
void
BioSystem::setInitialValues(Vector const& val)
{
    Species const&  spec = _ode.getSpecies();
    StrIterConst    sBeg = spec.begin();
    StrIterConst    sEnd = spec.end();

    if ( val.nr() >= (long)spec.size() )
    {
        long j = 1;
        for (StrIterConst it = sBeg; it != sEnd; ++it)
            _iniPar[*it] = val(j++);
    }
}
//---------------------------------------------------------------------------
ODETrajectory*
BioSystem::getEvaluationTrajectories()
{
    //$$$ return dynamic_cast<DOP853*>(_odeSolver)->getRawTrajectory();
    //if ( _odeSolver->getId() == ODE_SOLVER_LIMEX_A )
    //{
    //    return dynamic_cast<LIMEX_A*>(_odeSolver)->getRawTrajectory();
    //}
    return _odeSolver -> getRawTrajectory();
    //return 0;
}
//---------------------------------------------------------------------------
ODESolver::Trajectory
BioSystem::getOdeTrajectory() const
{
    ODESolver::Trajectory& tra = _odeSolver -> getSolutionTrajectory();
    ODESolver::Trajectory ret;
    //$$$ return dynamic_cast<DOP853*>(_odeSolver)->getSolutionTrajectory();
    // return _odeSolver -> getSolutionTrajectory();
    for (long j = 1; j < (long)tra.size(); ++j)
    {
        ret[j-1] = tra[j];
    }
    return ret;
}
//---------------------------------------------------------------------------
Vector
BioSystem::getOdeTrajectory(long j) const
{
    //$$$ ODESolver::Trajectory& tra = dynamic_cast<DOP853*>(_odeSolver) ->
    //$$$                                                getSolutionTrajectory();
    ODESolver::Trajectory& tra = (_odeSolver -> getSolutionTrajectory());

    if ( (0 <= j) && (j < (long)tra.size()-1) ) return Vector( tra[j+1] );

    return Vector();
}
//---------------------------------------------------------------------------
Vector
BioSystem::getSimTrajectoryPoints(long j)
{
    Species const& species = _ode.getSpecies();
    long           n       = species.size();
    long           T       = _tpMeas.size();
    Vector         v;

    if ( (0 <= j) && (j < n) )
    {
        v.zeros(T);
        for (long tp = 1; tp <= T; ++tp)
            v(tp) = ( _synData[tp-1][species[j]] ).first;
    }

    return v;
}
//---------------------------------------------------------------------------
Vector
BioSystem::getSimTrajectoryPoints(std::string spec)
{
    Species const&  species = _ode.getSpecies();
    StrIterConst    sBeg    = species.begin();
    StrIterConst    sEnd    = species.end();
    long            T       = _tpMeas.size();
    Vector          v;

    StrIterConst it = std::find(sBeg, sEnd, spec);

    if ( it != sEnd )
    {
        v.zeros(T);
        for (long tp = 1; tp <= T; ++tp)
            v(tp) = ( _synData[tp-1][*it] ).first;
    }

    return v;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void
BioSystem::iniODE()
{
    Species const&  spec = _ode.getSpecies();
    StrIterConst    sBeg = spec.begin();
    StrIterConst    sEnd = spec.end();

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        _sysPar[*it] = _iniPar[*it];
    }

    _sysPar["odeTime"] = _tInterval[0];
}
//---------------------------------------------------------------------------
void
BioSystem::iniODE(long n, double* y)
{
    long            j = 0;
    Species const&  spec = _ode.getSpecies();
    StrIterConst    sBeg = spec.begin();
    StrIterConst    sEnd = spec.end();

    y[j++] = _tInterval[0];

    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        if (j < n) y[j++] = _iniPar[*it];
    }
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Vector
BioSystem::computeModel(Expression::Param const& var, std::string mode)
{
    Species const&      species   = _ode.getSpecies();
    Parameter const&    parameter = _ode.getParameters();
    long                k = 0;
    long                n = species.size() + 1;
    long                q = parameter.size();
    long                T; //  = _tpMeas.size();
    long                Tb = _tInterval.size();
    StrIterConst        sBeg = species.begin();
    StrIterConst        sEnd = species.end();
    StrIterConst        pBeg = parameter.begin();
    StrIterConst        pEnd = parameter.end();
    ODESolver::Grid     iniValues(n);
    double*             y = new double[n];
    Vector              model; // (n*T);
    // long                tmin = std::min( _tpMeas.size(), _measData.size() );

    // iniODE();
    iniODE(n,y);

    //    Expression::ParamIterConst vBeg = var.begin();
    //    Expression::ParamIterConst vEnd = var.end();
    //    for (Expression::ParamIterConst it = vBeg; it != vEnd; ++it)
    //        _sysPar[it->first] = it->second;
    setParamValues(var);

    _optPar = var;


    k = 0;
    for (StrIterConst it = pBeg; it != pEnd; ++it)
    {
        if (k < q) _parValue[k++] = _sysPar[*it];
    }

    _ode.setParBase(_parValue);
    for (long j = 1; j < Tb; ++j)
    {
        _iniCond[j-1].setParBase(_parValue);
    }

    if ( mode == "adaptive" )
    {
        _tpMeas.clear();
        T = -1;
    }

    iniValues[0] = _tInterval[0];

    //$$$ dynamic_cast<DOP853*>(_odeSolver) ->
//    dynamic_cast<LIMEX_A*>(_odeSolver) ->
//                            setODESystem(
//                                            BioSystemWrapper::fcnODE,
//                                            BioSystemWrapper::jacODE,
//                                            _tInterval[0],
//                                            iniValues, _tpMeas,
//                                            _tInterval[Tb-1]
//                                            // BioSystemWrapper::outODE
//                                        );
//
//    BioSystemWrapper::setObj(*this);
    _odeSolver -> setODESystem(
                                *_systemOde,
                                _tInterval[0],
                                iniValues, _tpMeas,
                                _tInterval[Tb-1]
                               );

    dynamic_cast<BioSystemODE*>(_systemOde) -> setObj(*this);


// TIME_THIS_TO( std::cerr << "*** BioSystem::computeModel() : ";

    for (long j = 1; j < Tb; ++j)
    {
        // _iniCond[j-1].f( _sysPar, n, y );
        _iniCond[j-1].f( y, n, y );

        if ( _odeSolver->getId() == ODE_SOLVER_LIMEX_A )
        {
            dynamic_cast<LIMEX_A*>(_odeSolver) -> resetSuccessiveCallFlag();
        }

        //$$$ std::cerr << "\n*** dop853: rc = ";
        // std::cerr << "\n*** BioSystem::computeModel: limex_a: rc = ";
        _odeErrorFlag = _odeSolver -> integrate(n, y, _tInterval[j-1], _tInterval[j]);
        //$$$ _odeErrorFlag = dynamic_cast<LIMEX_A*>(_odeSolver) -> integrateWithoutInterpolation(n, y, _tInterval[j-1], _tInterval[j]);
        // std::cerr << _odeErrorFlag << std::endl;

        /*
        double* ytmp = y;
        for (StrIterConst it = sBeg; it != sEnd; ++it)
        {
            _sysPar[*it] = *ytmp++;
        }
        */
    }

// , std::cerr )

    ODESolver::Trajectory sim;

    if ( _tpMeas.size() <= 0 )
    {
        ODESolver::Grid grid = _odeSolver -> getAdaptiveGridPoints();
        T = grid.size();
        for (long tp = 0; tp < T; ++tp)
        {
            _tpMeas.push_back( grid[tp] );
        }
        sim = _odeSolver -> getAdaptiveSolution();
    }
    else
    {
        T = _tpMeas.size();
        sim = _odeSolver -> getDataTrajectory();
    }


    _linftyModel.clear();
    _synData.clear();

    k = 1;  // starting from k = 1, since k = 0 is the (dummy) time variable!
    for (StrIterConst it = sBeg; it != sEnd; ++it)
    {
        _linftyModel[*it] = std::fabs( sim[k++][0] );
    }

    for (long tp = 0; tp < T; ++tp)
    {
        if ( tp > 0 )
        {
            k = 1; // starting from k = 1, see above comment!
            for (StrIterConst it = sBeg; it != sEnd; ++it)
            {
                Real tmp = std::fabs( sim[k++][tp] );

                if ( tmp > _linftyModel[*it] )
                {
                    _linftyModel[*it] = tmp;
                }
            }
        }

        k = 1;
        _synData.push_back( MeasurementPoint() );

        for (StrIterConst it = sBeg; it != sEnd; ++it)
        {
            if ( k < (long)sim.size() )
            {
                _synData[tp][*it] = std::make_pair<Real,Real>( sim[k++][tp], 1.0 );
            }
            else
            {
                _synData[tp][*it] = std::make_pair<Real,Real>( 0.0, GREAT );
            }
        }
    }

    if ( mode == "adaptive" ) { _measData = _synData; _totmeasData = (n-1)*T; }
    if ( mode == "init" ) { _measData = _synData; _totmeasData = (n-1)*T; }

    k = 0;
    model.zeros(_totmeasData);

    for (long tp = 0; tp < T; ++tp)
    {
        MeasIterConst mBeg = _measData[tp].begin();
        MeasIterConst mEnd = _measData[tp].end();

        for (MeasIterConst it = mBeg; it != mEnd; ++it)
        {
            std::string spec = it->first;

            model(++k) = ( _synData[tp][spec] ).first;
        }
    }

    delete[] y;

    return model;
}
//---------------------------------------------------------------------------
Matrix
BioSystem::computeJacobian(Expression::Param const& var, std::string mode)
{
    Expression::ParamIterConst vBeg = var.begin();
    Expression::ParamIterConst vEnd = var.end();
    Expression::Param          idx;
    long                       k = 0;

    idx.clear();
    for (Expression::ParamIterConst it = vBeg; it != vEnd; ++it)
    {
        idx[it->first] = ++k;
    }

    return computeJacobian(var, idx, mode);
}
//---------------------------------------------------------------------------
Matrix
BioSystem::computeJacobian(Expression::Param const& var,
                           Expression::Param&       idx,
                           std::string              mode)
{
    Species const&      species   = _ode.getSpecies();
    Parameter const&    parameter = _ode.getParameters();
    long                j = 0;
    long                k = 0;
    long                n = species.size() + 1L;  // one more for virtual time variable
    long                q = parameter.size();
    long                qq = var.size();
    long                T; //  = _tpMeas.size();
    long                Tb = _tInterval.size();
    long                N = n + n*qq;  // n + n*q;
    std::vector<Real>   iniValues(N, 0.0);
    double*             yy = new double[N];
    double*             Zp = yy + n;
    // Matrix              Z; // (n,qq);
    // Matrix              jacobian(n*T, qq);
    Matrix              jacobian; // (_totmeasData, qq);
    StrIterConst        sBeg = species.begin();
    StrIterConst        sEnd = species.end();
    StrIterConst        pBeg = parameter.begin();
    StrIterConst        pEnd = parameter.end();
/*
    Real                rTolSave = _odeSolver -> getRTol();
    Real                aTolSave = _odeSolver -> getATol();
*/
    // Z.zeros(n,qq);
    // jacobian.zeros(_totmeasData,qq);

    for (long jj = 0; jj < N; ++jj) yy[jj] = 0.0;
    // for (long jj = 0; jj < Tb; ++jj) _iniCond[jj].setParameters(var);
    // iniODE();
    iniODE(n,yy);

    Expression::ParamIterConst vBeg = var.begin();
    Expression::ParamIterConst vEnd = var.end();
    //    for (Expression::ParamIterConst it = vBeg; it != vEnd; ++it)
    //    {
    //        _sysPar[it->first] = it->second;
    //    }
    setParamValues(var);

    _optPar = var;

//std::cerr << "*** BioSystem::computeJacobian(...) : current parameter values ***" << std::endl;
//    k = 0;
//    for (StrIterConst it = pBeg; it != pEnd; ++it)
//    {
//std::cerr << std::setw(5) << k++ << ": '" << *it << "' = " << _sysPar[*it] << std::endl;
//    }
//std::cerr << "***" << std::endl;

    k = 0;
    for (StrIterConst it = pBeg; it != pEnd; ++it)
    {
        if (k < q) _parValue[k++] = _sysPar[*it];
    }

    _ode.setParBase(_parValue);
    for (long jj = 1; jj < Tb; ++jj)
    {
        _iniCond[jj-1].setParBase(_parValue);
    }

    /// ODESolver::Grid tp;

    if ( mode == "adaptive" )
    {
        _tpMeas.clear();
        T = -1;
    }

    iniValues[0] = _tInterval[0];

/*
    _odeSolver -> setRTol(1.0e-2);
    _odeSolver -> setATol(1.0e-3);
*/

//    dynamic_cast<LIMEX_A*>(_odeSolver) ->
//                            setODESystem(
//                                            BioSystemWrapper::fcnVar,
//                                            BioSystemWrapper::jacVar,
//                                            _tInterval[0],
//                                            iniValues, _tpMeas,
//                                            _tInterval[Tb-1]
//                                            , (int) n
//                                            // BioSystemWrapper::outVar
//                                        );
//    // dynamic_cast<LIMEX_A*>(_odeSolver) -> resetSuccessiveCallFlag();
//    BioSystemWrapper::setObj(*this);

    _odeSolver -> setODESystem(
                                    *_variationalOde,
                                    _tInterval[0],
                                    iniValues, _tpMeas,
                                    _tInterval[Tb-1]
                                    , (int) n
                               );

    dynamic_cast<BioSystemVAR*>(_variationalOde) -> setObj(*this);


// TIME_THIS_TO( std::cerr << "*** BioSystem::computeJacobian() : ";

    for (long jj = 1; jj < Tb; ++jj)
    {
        // double* yytmp;

        // Z = _iniCond[jj-1].df( Z, _sysPar );
       // _iniCond[jj-1].f( _sysPar, n, yy );
       _iniCond[jj-1].df( var, Zp, yy, n, qq, Zp );
       _iniCond[jj-1].f( yy, n, yy );
       /*
        yytmp = yy + n;
        for (long r = 1; r <= n; ++r)
            for (long s = 1; s <= qq; ++s)
            {
                *yytmp++ = Z(r,s);
            }
        */

        if ( _odeSolver->getId() == ODE_SOLVER_LIMEX_A )
        {
            dynamic_cast<LIMEX_A*>(_odeSolver) -> resetSuccessiveCallFlag();
        }

        //$$$ std::cerr << "*** dop853: rc = " <<
        // std::cerr << "\n*** BioSystem::computeJacobian() : limex_a: rc = " <<
/*
        _odeErrorFlag = _odeSolver ->
                            integrate(N, yy, _tInterval[jj-1], _tInterval[jj]);
*/
        _odeErrorFlag = _odeSolver ->
                            integrateSensitivitySystem(n, N, yy, _tInterval[jj-1], _tInterval[jj]);
        // std::cerr << std::endl;

        /*
        yytmp = yy + n;
        for (StrIterConst it = sBeg; it != sEnd; ++it)
        {
            _sysPar[*it] = *yytmp++;
        }
        for (long r = 1; r <= n; ++r)
            for (long s = 1; s <= qq; ++s)
            {
                Z(r,s) = *yytmp++;
            }
        */
    }

// , std::cerr )


    ODESolver::Trajectory sim;

    if ( _tpMeas.size() <= 0 )
    {
        ODESolver::Grid grid = _odeSolver -> getAdaptiveGridPoints();
        T = grid.size();
        for (long tp = 0; tp < T; ++tp)
        {
            _tpMeas.push_back( grid[tp] );
        }
        sim = _odeSolver -> getAdaptiveSolution();
    }
    else
    {
        T = _tpMeas.size();
        sim = _odeSolver -> getDataTrajectory();
    }

/*
    _odeSolver -> setRTol( rTolSave );
    _odeSolver -> setATol( aTolSave );
*/

    _linftyModel.clear();
    _synData.clear();
    _jacobian.clear();

    /// n--; // from here on, the time variable is gone since "sim" contains only species
    /// 15.03.13 td: *** NO!  ODESolver is now dull again: ***
    ///              _No_ knowledge about the system(s) it shall integrate!!!
    k = 1;
    for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
    {
        _linftyModel[*itSpe] = std::fabs( sim[k++][0] );
    }

    for (long tp = 0; tp < T; ++tp)
    {
        if ( tp > 0 )
        {
            k = 1;
            for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
            {
                Real tmp = std::fabs( sim[k++][tp] );

                if ( tmp > _linftyModel[*itSpe] )
                {
                    _linftyModel[*itSpe] = tmp;
                }
            }
        }

        k = 1;
        _synData.push_back( MeasurementPoint() );

        for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
        {
            if ( k < n )
            {
                _synData[tp][*itSpe] = std::make_pair<Real,Real>( sim[k++][tp], 1.0 );
            }
            else
            {
                _synData[tp][*itSpe] = std::make_pair<Real,Real>( 0.0, GREAT );
            }
        }

        // now skip to the results of the variational equations (k = n+qq is correct here!)
        // 13.09.12 td: !! see comment below !!  (Now k = n is correct.  Really!)
        // 15.03.13 td: The comments from Sep 13th are still in place,
        //              and correct as it seems ...
        k = n; // + qq;
        _jacobian.push_back( MeasurementPoint() );

        // 13.09.12 td: Note, here  we first iterate over the parameters because there had
        //              been a reordering in yy (for LIMEX), and therefore in the sim array.
        //              And, most importantly (!!!), the k++ offsets are changed: k = n from
        //              above is now correct.  Indeed!

        // for (StrIterConst itPar = pBeg; itPar != pEnd; ++itPar)
        for (Expression::ParamIterConst itPar = vBeg; itPar != vEnd; ++itPar)
        {
            k++;

            for (StrIterConst itSpe = sBeg; itSpe != sEnd; ++itSpe)
            {
                std::string s = *itSpe + " / " + (itPar->first);

                if ( k < (long)sim.size() )
                {
                    _jacobian[tp][s] = std::make_pair<Real,Real>( sim[k++][tp], 1.0 );
                }
                else
                {
                    _jacobian[tp][s] = std::make_pair<Real,Real>( 0.0, GREAT );
                }
            }
        }
    }

    if ( mode == "adaptive" ) { _measData = _synData; _totmeasData = (n-1)*T; }

    j = 0;
    jacobian.zeros(_totmeasData, qq);

    /*
    if ( mode == "adaptive" )
    {
        for (long tp = 0; tp < T; ++tp)
        {
            for (StrIterConst it = sBeg; it != sEnd; ++it)
            {
                ++j; // k = 0;
                for (Expression::ParamIterConst itPar = vBeg;
                                                itPar != vEnd; ++itPar)
                {
                    std::string s;
                    s = (*it) + " / " + (itPar->first);
                    k = idx[(itPar->first)];

                    jacobian( j, k ) = (_jacobian[tp][s]).first;
                }
            }
        }
    }
    else
    */
    {
        for (long tp = 0; tp < T; ++tp)
        {
            MeasIterConst mBeg = _measData[tp].begin();
            MeasIterConst mEnd = _measData[tp].end();

            for (MeasIterConst it = mBeg; it != mEnd; ++it)
            {
                ++j; // k = 0;
                for (Expression::ParamIterConst itPar = vBeg;
                                                itPar != vEnd; ++itPar)
                {
                    std::string s;
                    s = (it->first) + " / " + (itPar->first);
                    k = idx[(itPar->first)];

                    jacobian( j, k ) = (_jacobian[tp][s]).first;
                }
            }
        }
    }

    delete[] yy;

    return jacobian;
}
//---------------------------------------------------------------------------
int
BioSystem::getComputeErrorFlag()
{
    return _odeErrorFlag;
}
//---------------------------------------------------------------------------
QRconDecomp
BioSystem::computeSensitivity(Expression::Param&    var,
                              Expression::Param&    vscal,
                              std::string           mode)
{
    Expression::ParamIterConst  vBeg = var.begin();
    Expression::ParamIterConst  vEnd = var.end();
    Real                        rtol = _odeSolver->getRTol();
    long                        qq = var.size();
    long                        k = 0;
    Vector                      pw;
    // Matrix                      jac;

    pw.ones(qq);

    if ( qq == (long)vscal.size() )
    {
        k = 0;
        for (Expression::ParamIterConst itVar = vBeg;
                                        itVar != vEnd; ++itVar)
        {
            pw(++k) = std::max( std::fabs(vscal[itVar->first]), rtol );
        }
    }

    if ( mode == "inner" )
    {
        _jac = computeJacobian(var);
    }
    else
    {
        Vector fh, f;
        Real   ajmin = SMALL;
        Real   ajdelta = sqrtEPMACH;
        Real   u, w;
        int    su;

        f = computeModel(var, mode);

        k = 0;
        _jac.zeros( f.nr(), qq );

        for (Expression::ParamIterConst itVar = vBeg;
                                        itVar != vEnd; ++itVar)
        {
            ++k;

            w  = itVar->second;

            su = ( w < 0.0 ) ? -1 : 1;
            u  = std::max( std::max(std::fabs(w), ajmin), pw(k) );
            u *= ajdelta * su;
            var[itVar->first] = w + u;

            fh = computeModel(var);

            var[itVar->first] = w;
            _jac.set_colm(k) = (1.0/u) * Matrix(fh - f);
        }
    }

    _jac = _jac * pw.diag();
    // if ( qScal ) computeRowScaling();

    return _jac.factorQRcon( 0, 0, 1.0/std::sqrt(rtol) );
}
//---------------------------------------------------------------------------
Matrix
BioSystem::getSensitivityMatrix()
{
    return _jac;
}
////---------------------------------------------------------------------------
//BioSystem::MeasurementList const&
//BioSystem::getSensitivityList()
//{
//    return _jacobian;
//}
////---------------------------------------------------------------------------

///
///
///
/*
//---------------------------------------------------------------------------
//$$$ extern "C"
//$$$ void
//$$$ BioSystemWrapper::fcnODE(
//$$$                         unsigned   n,
//$$$                         double     t,
//$$$                         double*    y,
//$$$                         double*    dy,
//$$$                         double*    cd  // 'client data'
//$$$                        )
extern "C"
void
BioSystemWrapper::fcnODE(
                              int*      n,
                              int*      nz,
                              double*   t,
                              double*   y,
                              double*   dy,
                              double*   B,
                              int*      ir,
                              int*      ic,
                              int*      info
                        )
{
//    const double                one = 1.0;
    // Expression::Param&          sys  = _obj->getSysPar();
    // BioRHS                      ode  = _obj->getODE();
    // / *
    // BioRHS::Species const&      spec = _ode.getSpecies();
    // StrIterConst                sBeg = spec.begin();
    // StrIterConst                sEnd = spec.end();
    // * /

    // *nz = *n;
    y[0] = *t;
    // sys["odeTime"] = *t;
    // for (StrIterConst it = sBeg; it != sEnd; ++it) sys[*it] = *y++;

    // Vector dz =
    // _ode.f( sys, *n, dy );
    _ode.f( y, *n, dy );

    // for (long j = 1; j <= dz.nr(); ++j) *(dy++) = dz(j);

    // if ( *B != one )
//    {
//        for (int j = 1; j <= *nz; ++j)
//        {
//            *B++ = one;
//            *ir++ = *ic++ = j;
//        }
//    }
    _ode.b( y, nz, B, ir, ic);

    //

    *info = 0;
}
//---------------------------------------------------------------------------
extern "C"
void
BioSystemWrapper::jacODE(
                            int*    n,
                            double* t,
                            double* y,
                            double* dy,
                            double* J,
                            int*    ldJ,
                            int*    ml,
                            int*    mu,
                            int*    full_or_band,
                            int*    info
                        )
{
    Expression::Param&      sys  = _obj->getSysPar();
    // BioRHS                  ode  = _obj->getODE();
    BioRHS::Species const&  spec = _ode.getSpecies();
    StrIterConst            sBeg = spec.begin();
    StrIterConst            sEnd = spec.end();

    sys["odeTime"] = *t;
    y++; // skip first component as it is reserved as time variable
    for (StrIterConst it = sBeg; it != sEnd; ++it) sys[*it] = *y++;

    _ode.Jf( sys, *n, J, *ldJ );
    // Matrix Fz = _obj -> rhsJac( sys );
    // long   gap = (long)(*ldJ) - Fz.nr();
    //
    // for (long k = 1; k <= Fz.nc(); ++k)
    // {
    //     for (long j = 1; j <= Fz.nr(); ++j)
    //     {
    //         *(J++) = Fz(j,k);
    //     }
    //
    //     J += gap;
    // }

    *info = 0;
}
//---------------------------------------------------------------------------
//extern "C"
//void
//BioSystemWrapper::outODE(
//                         long       nr,
//                         double     told,
//                         double     t,
//                         double*    y,
//                         unsigned   n,
//                         int*       irtrn
//                        )
//{
////    double* z = y;
//    BioSystem::TimePoints& tp = _obj->getOdeTrajectoryTimePoints();
//
//    if ( tp.empty() || (tp.back() == told) )
//    {
//        ODESolver::Trajectory&  tra = _obj->getOdeTrajectories();
//        BioSystem::Species      spe = _obj->getSpecies();
//        StrIterConst            beg = spe.begin();
//        StrIterConst            end = spe.end();
//
//        tp.push_back(t);
//
//        for (StrIterConst it = beg; it != end; ++it)
//            tra[*it].push_back( *(y++) );
//    }
//
//    //_obj->interpolate();
//
////    std::cerr << std::scientific;
////    std::cerr << "#nr=" << nr << "  t=" << t;
////    std::cerr << "\t z[0]=" << z[0] << "\t z[1]=" << z[1] << std::endl;
//
//    *irtrn = 0;
//}
//---------------------------------------------------------------------------
//$$$ extern "C"
//$$$ void
//$$$ BioSystemWrapper::fcnVar(
//$$$                          unsigned   nq,
//$$$                          double     t,
//$$$                          double*    yyu,
//$$$                          double*    dydyu,
//$$$                          double*    cd
//$$$                         )
extern "C"
void
BioSystemWrapper::fcnVar(
                         int*       nq,
                         int*       nqz,
                         double*    t,
                         double*    yyu,
                         double*    dydyu,
                         double*    B,
                         int*       ir,
                         int*       ic,
                         int*       info
                        )
{
    const double                one = 1.0;
    Expression::Param&          opt  = _obj->getOptPar();
//Expression::Param&          sys  = _obj->getSysPar();
//BioRHS::Species const&      spe  = _ode.getSpecies();
//StrIterConst                sBeg = spe.begin();
//StrIterConst                sEnd = spe.end();
    long                        n    = _ode.getSpecies().size() + 1L;
    long                        qq   = opt.size();          // ((long)(*nq) - n) / n;
    double*                     yy   = yyu + n;
//Matrix                      Z(n,qq);
//Vector                      y(n);

    *nqz = *nq;
    yy[0] =
    yyu[0] = *t;
//sys["odeTime"] = *t;
//yyu++;
//for (StrIterConst it = sBeg; it != sEnd; ++it) sys[*it] = *yyu++;
//
//for (long k = 1; k <= qq; ++k)
//{
//    for (long j = 1; j <= n; ++j)
//    {
//        Z(j,k) = *yy++;
//    }
//}
//
//Z = _ode.df( opt, Z, sys );
//y = _ode.f( sys );

    _ode.df( opt, yy, yyu, n, qq, dydyu+n );
    _ode.f( yyu, n, dydyu );

//for (long j = 1; j <= n; ++j)
//{
//    *dydyu++ = y(j);
//}
//
// / *
//std::cerr << std::endl;
//std::cerr << "### BioSystemWrapper::fcnVar() ###" << std::endl;
//std::cerr << "   Z = " << std::endl;
//std::cerr << Z.t() << std::endl;
// * /
//
//for (long k = 1; k <= qq; ++k)
//{
//    for (long j = 1; j <= n; ++j)
//    {
//        *dydyu++ = Z(j,k);
//    }
//}


    for (int j = 1; j <= *nqz; ++j)
    {
        *B++ = one;
        *ir++ = *ic++ = j;
    }

    *info = 0;
}
//---------------------------------------------------------------------------
extern "C"
void
BioSystemWrapper::jacVar(
                            int*    nq,
                            double* t,
                            double* yu,
                            double* dyu,
                            double* Ju,
                            int*    ldJu,
                            int*    ml,
                            int*    mu,
                            int*    full_or_band,
                            int*    info
                        )
{
    Expression::Param       sys  = _obj->getSysPar();
    BioRHS::Species const&  spec = _ode.getSpecies();
    StrIterConst            sBeg = spec.begin();
    StrIterConst            sEnd = spec.end();

    sys["odeTime"] = *t;
    yu++; // skip first component since it is reserved as time variable
    for (StrIterConst it = sBeg; it != sEnd; ++it) sys[*it] = *yu++;

    Matrix Fz = _ode.Jf( sys );
    long   n  = Fz.nr();
    // long   gap = (long)(*ldJu) - n;

    if ( *full_or_band == 0 )
    {

        long  q1 = (long)(*nq) / n;    // q1 := q + 1; loop below counts til m < q1 !

        for (long k = 1; k <= n; ++k)
        {
            for (long j = 1; j <= n; ++j)
            {
                double tmp = Fz(j,k);

                for (long m = 0; m < q1; ++m)
                {
                    Ju[ (*ldJu)*(k-1 + m*n) + (j-1 + m*n) ] = tmp;
                                            //  Ju(j + q*n, k + q*n) = tmp
                }
            }
        }

    }
    else // if ( ((long)(*mu) == n) && ((long)(*ml) == n) )
    {
// / *
//std::cerr << std::endl;
//std::cerr << "##### BioSystemWrapper::jacVar() #####" << std::endl;
//std::cerr << "  *ldJu = " << *ldJu << std::endl;
//std::cerr << "    *mu = " << *mu << std::endl;
//std::cerr << "    *ml = " << *ml << std::endl;
//std::cerr << "      n = " << n << std::endl;
//std::cerr << std::endl;
// * /
        // long  nn = 2*n + 1L;
        long  q1 = (long)(*nq) / n;    // q1 := q + 1; loop below counts til m < q1 !

        for (long k = 1; k <= n; ++k)
        {
            long off = n + 1L - k;

            for (long j = 1; j <= n; ++j)
            {
                double tmp = Fz(j,k);

                for (long m = 0; m < q1; ++m)
                {
                    Ju[ (*ldJu)*(k-1 + m*n) + (j-1 + off) ] = tmp;
                                            //  Ju(j + q*n, k + q*n) = tmp
                }
            }
        }
    }
// / *
//    else
//    {
//        // long q1 = (long)(*nq) / n;
//
//        for (long k = 1; k <= *nq; ++k)
//        {
//            long upp = std::max(          1L, k - (long)(*mu) );
//            long low = std::min( (long)(*nq), k + (long)(*ml) );
//            long off = (long)(*mu) + 1 - k;
//
//            long  q = 1 + (k-1) % n;
//            long  m = 1 + (k-1) / n;
//            long mm = m + n;
//            long jj = 0;
//
//            // for (long jj = 1, j = upp; j <= low; ++j)
//            for (long j = upp; j <= low; ++j)
//            {
//                double tmp = 0.0;
//                // long m = (k-1) / n;
//
//                // if ( (m*n < j) && (j <= (m+1)*n) ) tmp = Fz( jj++, 1 + (k-1)%n );
//                if ( (m <= j) && (j <= mm) ) tmp = Fz( ++jj, q );
//
//                Ju[ (*ldJu)*(k-1) + (j-1+off) ] = tmp;
//                                                    //  Ju(j+off, k) = tmp
//            }
//        }
//    }
// * /

    *info = 0;
}
//---------------------------------------------------------------------------
//extern "C"
//void
//BioSystemWrapper::outVar(
//                         long       nr,
//                         double     told,
//                         double     t,
//                         double*    yu,
//                         unsigned   nq,
//                         int*       irtrn
//                        )
//{
////    double* z = y;
//    BioSystem::TimePoints& tp = _obj->getVarTrajectoryTimePoints();
//
//    if ( tp.empty() || (tp.back() == told) )
//    {
//        ODESolver::Trajectory&  vtra = _obj->getVarTrajectories();
//        BioSystem::Species      spec = _obj->getSpecies();
//        BioSystem::Parameter    parm = _obj->getParameters();
//        StrIterConst            sBeg = spec.begin();
//        StrIterConst            sEnd = spec.end();
//        StrIterConst            pBeg = parm.begin();
//        StrIterConst            pEnd = parm.end();
//
//        tp.push_back(t);
//
//        for (StrIterConst itSpec = sBeg; itSpec != sEnd; ++itSpec)
//            for (StrIterConst itParm = pBeg; itParm != pEnd; ++itParm)
//            {
//                std::string s = *itSpec + " / " + *itParm;
//
//                vtra[s].push_back( *(yu++) );
//            }
//    }
////    std::cerr << std::scientific;
////    std::cerr << "#nr=" << nr << "  t=" << t;
////    std::cerr << "\t z[0]=" << z[0] << "\t z[1]=" << z[1] << std::endl;
//
//    *irtrn = 0;
//}
//---------------------------------------------------------------------------
*/
