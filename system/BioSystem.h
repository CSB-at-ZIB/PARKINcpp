// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-03 td
// last changed:
//
#ifndef __BIO_SYSTEM_H
#define __BIO_SYSTEM_H

#include <utility>  // for `std::pair'
#include <string>
#include <vector>
#include <map>

#include "linalg/Matrix.h"
#include "linalg/Vector.h"

#include "odelib/ODESolver.h"
#include "odelib/ODETrajectory.h"
#include "odelib/FirstOrderODESystem.h"
#include "Expression.h"
#include "BioRHS.h"

namespace PARKIN
{
    typedef std::map< std::string, std::pair<Real,Real> >    MeasurementPoint;

    class QRconDecomp;

    class BioSystem
    {

        public:

            typedef BioRHS::Species                         Species;
            typedef BioRHS::Parameter                       Parameter;
            typedef BioRHS::ExprTypeMap                     ExprTypeMap;
            typedef BioRHS::ExpressionMap                   ExpressionMap;
            typedef std::vector< MeasurementPoint >         MeasurementList;

            /// c'tor
            BioSystem(  Real tStart = 0.0, Real tEnd = 1.0 );
            BioSystem(  Vector const& tInterval );

            BioSystem(  ExpressionMap const&            eList,
                        ODESolver::Grid const&          tPoints,
                        MeasurementList const&          meas,
                        ODESolver::Grid const&          tInterval
                     );
            /// d'tor
            virtual ~BioSystem();
            /// copy c'tor
            BioSystem(BioSystem const& biosys);
            /// assignment
//            BioSystem const& operator= (BioSystem const& biosys);

            ///

//            virtual Vector
//            fcn(Vector const& x, int& ifail);
//
//            virtual Matrix
//            jac(Vector const& x, int& ifail);

            void
            setSolverDebugFlag(int flag);

            void
            setSolverInterpolationFlag(int cubint);

            ///

            Real
            getSolverRTol() const;

            void
            setSolverRTol(Real tol);

            Real
            getSolverATol() const;

            void
            setSolverATol(Real tol);

            Real
            getSystemTol() const;

            void
            setSystemTol(Real tol);

            //

            Expression::Param&
            getSysPar();

            Expression::Param&
            getOptPar();

            //

            Expression::Param&
            getLinftyModel();

            //

            Species
            getSpecies() const;

            void
            resetSpecies(Species const& species);   // Species serve also as index
                                                    // to BioRHS, therefore we will
                                                    // have a resetSpecies method here!
            //

            Parameter
            getParameters() const;

            void
            setParameters(Parameter const& parameter);

            //

            Vector
            getBreakpoints() const;

            void
            setBreakpoints(Vector const& tInterval);

            //

            ExpressionMap
            getEvent(long j) const;

            void
            setEvent(long j, ExpressionMap const& emap);

            //

            BioRHS const&
            getODE() const
            { return _ode; }

            //

            ExpressionMap const&
            getVarExpr() const
            { return _ode.getdRHS(); }
            //
            ExpressionMap const&
            getODEExpr() const
            { return _ode.getRHS(); }

            ExprTypeMap const&
            getODETypes() const
            { return _ode.getRHSType(); }

            //

            void
            setODESystem(ExpressionMap const& elist);

            void
            setODESystem(ExpressionMap const& elist, ExprTypeMap const& tlist);

            void
            setODETypes(ExprTypeMap const& tlist);

            //

            MeasurementList const&
            getMeasurementList();

            MeasurementList const&
            getSimTrajectoryList();

            // void
            // setMeasurementList(ODESolver::Grid const& tp, MeasurementList const& meas);

            void
            setMeasurementList(Vector const& tp, MeasurementList const& meas);

            void
            setMeasurementList(MeasurementList const& meas);

            void
            setEmptyMeasurementList();

            Vector
            getMeasurements() const;

            Vector
            getMeasurementWeights() const;

            //

//            void
//            setMeasurementTimePoints(ODESolver::Grid const& tp);

            void
            setMeasurementTimePoints(Vector const& tp);

            Vector
            getMeasurementTimePoints() const;

            Vector  // ODESolver::Grid&
            getOdeTrajectoryTimePoints() const;

            //

            Real
            getParamValue(std::string const& name); // const;

            void
            setParamValue(std::string const& name, Real value);

            void
            setParamValues(Expression::Param const& par);

            //

            Real
            getInitialValue(std::string const& name); // const;

            void
            setInitialValue(std::string const& name, Real value);

            void
            setInitialValues(Expression::Param const& inipar);

            void
            setInitialValues(Vector const& val);

            ///

            ODETrajectory*
            getEvaluationTrajectories();

            ODESolver::Trajectory const&
            getOdeTrajectory() const;

            Vector
            getOdeTrajectory(long j) const;

            //

            Vector
            getSimTrajectoryPoints(long j);

            Vector
            getSimTrajectoryPoints(std::string spec);

            ///

            void
            iniODE();

            void
            iniODE(long n, double* y);

            ///

            Vector
            computeModel(Expression::Param const&   var,
                         std::string                mode = "");
                         // mode    -   "adaptive"   : generate artificial data in _measData  **and**  _tpMeas,
                         //             "init"       : compute _sysData at _tpMeas; and initialise _measData with these values,
                         //             "any-string" : only _synData is changed

            Matrix
            computeJacobian(Expression::Param const& var,
                            std::string              mode = "");
                            // mode -   "adaptive"   : generate Jacobian adaptively in _jacobian; **and** simultaneously change _tpMeas, _measData,
                            //          "any-string" : take _tpMeas as computation base, i.e. compute Jacobian at _tpMeas
            Matrix
            computeJacobian(Expression::Param const& var,
                            Expression::Param&       idx,
                            std::string              mode = "");
                            // mode -   "adaptive"   : generate Jacobian adaptively in _jacobian; **and** simultaneously change _tpMeas, _measData,
                            //          "any-string" : take _tpMeas as computation base, i.e. compute Jacobian at _tpMeas


            int
            getComputeErrorFlag();

            ///

            typedef MeasurementPoint::const_iterator                    MeasIterConst;
            typedef std::vector< std::string >::const_iterator          StrIterConst;


        private:

            //
            //      Note:   These two methods have been taken out for good.
            //      -----   See class 'BioProcessor' instead!
            //

            QRconDecomp
            computeSensitivity(Expression::Param&   var,
                               Expression::Param&   vscal,
                               std::string          mode = "initbla");
                               // mode  -   "inner"  : use variational equation to compute Jacobian
                               //           "outer"  : use difference quotient
                               //           "initbla": use difference quotient with newly artificial data

            ///

            Matrix
            getSensitivityMatrix();

//        public:
//            BioSystem::MeasurementList const&
//            getSensitivityList();

            ///

        private:

            void setIdentityEvents();

            BioRHS                          _ode;
            std::vector<BioRHS>             _iniCond;
            Expression::Param               _iniPar;
            Expression::Param               _sysPar;
            Expression::Param               _optPar;
            double*                         _parValue;

            Expression::Param               _linftyModel;

            ODESolver::Grid                 _tpMeas;
            MeasurementList                 _measData;      // e.g. _measData[3][species] = pair< value, std.dev >
            MeasurementList                 _synData;
            //std::vector<MeasurementPoint>   _jacobian;
            MeasurementList                 _jacobian;
            Matrix                          _jac;           // Jacobian (sensitivity matrix) as stagged blocks of
                                                            // solutions to the variational equation

            FirstOrderODESystem*            _systemOde;
            FirstOrderODESystem*            _variationalOde;
            ODESolver*                      _odeSolver;
            int                             _odeErrorFlag;

            long                            _totmeasData;
            ODESolver::Grid                 _tInterval;
            // double                          _tStart, _tEnd;
            Real                            _systemTol;
    };

    ///
    ///
    ///
    ///

    class BioSystemWrapper
    {
            typedef BioSystem::StrIterConst StrIterConst;

        public:
            //$$$ static void fcnVar(unsigned nq, double t, double* yu, double* dyu, double* cd);
            static void fcnVar(int* nq, int* nqz,
                               double* t, double* yu, double* dyu,
                               double* B, int *ir, int* ic,
                               int* info);
            static void jacVar(int* nq,
                               double* t, double* yu, double* dyu,
                               double* Ju, int* ldJu,
                               int* ml, int* mu, int* full_or_band,
                               int* info);

            //$$$ static void fcnODE(unsigned n, double t, double* y, double* dy, double* cd);
            static void fcnODE(int* n, int *nz,
                               double* t, double* y, double* dy,
                               double* B, int* ir, int* ic,
                               int* info);
            static void jacODE(int *n,
                               double* t, double* y, double* dy,
                               double* J, int* ldJ,
                               int* ml, int* mu, int* full_or_band,
                               int* info);

            static void setObj(BioSystem& obj)
            { _obj = &obj; _ode = _obj->getODE(); }

        private:
            static BioSystem*   _obj;
            static BioRHS       _ode;
    };

}
#endif // __BIO_SYSTEM_H
