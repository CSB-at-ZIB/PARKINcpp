// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-07 td
// last changed:
//
#ifndef __ODE_SOLVER_H
#define __ODE_SOLVER_H

#include <vector>
#include <map>
// #include <set>
#include "common/Constants.h"
#include "common/Types.h"
#include "FirstOrderODESystem.h"
#include "ODETrajectory.h"

namespace PARKIN
{
    ///

    enum ODESolverId
    {
        ODE_SOLVER_DEFAULT =  0 ,
        ODE_SOLVER_LIMEX_A =  1 ,
        ODE_SOLVER_DOP853  =  2 ,
        ODE_SOLVER_METAN_A =  3
    };

    ///

    class ODESolver
    {
        public:
            typedef std::vector< Real >                         Grid;
            typedef std::map< unsigned, std::vector<Real> >     Trajectory;
            typedef Grid::const_iterator                        GridIterConst;


            virtual void setODESystem(
                                            FirstOrderODESystem& ode,
                                            double              t0,
                                            Grid const&         y0,
                                            double              tEnd,
                                            int                  bandwidth = 0
                                       ) = 0;

            virtual void setODESystem(
                                            FirstOrderODESystem& ode,
                                            double              t0,
                                            Grid const&         y0,
                                            Grid const&         refGrid,
                                            double              tEnd,
                                            int                  bandwidth = 0
                                       ) = 0;

            virtual ~ODESolver() { }

            virtual ODESolver*      clone()                          = 0;

            virtual int             integrate()                      = 0;
            virtual int             integrate(unsigned n, double* yIni,
                                                double xLeft, double xRight) = 0;
            virtual int integrateSensitivitySystem(unsigned nDAE)   = 0;
            virtual int integrateSensitivitySystem(unsigned nDAE,
                                                    unsigned n, double* yIni,
                                                    double tLeft, double tRight) = 0;
            virtual Grid&           getAdaptiveGridPoints()           = 0;
            virtual Trajectory&     getAdaptiveSolution()             = 0;
            virtual Grid&           getSolutionGridPoints()           = 0;
            virtual Trajectory&     getSolutionTrajectory()           = 0;
            virtual Grid&           getDataGridPoints()               = 0;
            virtual Trajectory&     getDataTrajectory()               = 0;
            virtual ODETrajectory*  getRawTrajectory()                = 0;

            virtual std::string     getErrorMessage(int rc)           = 0;

            ///

            void    setDebugFlag(int flag) { _debugflag = flag; }
            int     getDebugFlag() const { return _debugflag; }

            void    setInterpolationFlag(int cubint) { _cubintflag = cubint; }
            int     getInterpolationFlag() const { return _cubintflag; }

            ///
            ODESolverId getId() const { return _solverId; }
            ///

            void setRTol(Real rtol) { _rtol = rtol; }
            void setATol(Real atol) { _atol = atol; }
            void setIniStep(Real inistep) { _inistep = inistep; }
            void setMaxStep(Real maxstep) { _maxstep = maxstep; }

            Real getRTol() const { return _rtol; }
            Real getATol() const { return _atol; }
            Real getIniStep() const { return _inistep; }
            Real getMaxStep() const { return _maxstep; }

        protected:
            explicit ODESolver(ODESolverId solverId) :
                                        _solverId(solverId),
                                        _debugflag(0),
                                        _cubintflag(-1),
                                        _atol(10*EPMACH),
                                        _rtol(1.0e-12),
                                        _inistep(1.0e-4),
                                        _maxstep(0.0)
            { }

            // The following precision setting seems to work only semi-stable:
            // ODESolver() : _atol(EPMACH), _rtol(1.0e-9) { }

            ODESolverId _solverId;
            int         _debugflag;
            int         _cubintflag;
            Real        _atol, _rtol;
            Real        _inistep, _maxstep;
    };

}
#endif // __ODE_SOLVER_H
