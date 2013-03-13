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

namespace PARKIN
{

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
            // virtual ODETrajectory*  getRawTrajectory()                = 0;

            ///

            void    setDebugFlag(int flag) { _debugflag = flag; }
            int     getDebugFlag() { return _debugflag; }

            void    setInterpolationFlag(int cubint) { _cubintflag = cubint; }
            int     getInterpolationFlag() { return _cubintflag; }

            ///

            void setRTol(Real rtol) { _rtol = rtol; }
            void setATol(Real atol) { _atol = atol; }

            Real getRTol() { return _rtol; }
            Real getATol() { return _atol; }

        protected:
            ODESolver() : _debugflag(0), _cubintflag(-1), _atol(10*EPMACH), _rtol(1.0e-12) { }

            // The following precision setting seems to work only semi-stable:
            // ODESolver() : _atol(EPMACH), _rtol(1.0e-9) { }

            int     _debugflag;
            int     _cubintflag;
            Real    _atol, _rtol;
    };

}
#endif // __ODE_SOLVER_H
