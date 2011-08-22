// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-07 td
// last changed:
//
#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <vector>
#include <map>
// #include <set>
#include "common/Constants.h"
#include "common/Types.h"

namespace PARKIN
{

    class ODESolver
    {
        public:
            typedef std::vector< Real >                         Grid;
            typedef std::map< unsigned, std::vector<Real> >     Trajectory;
            typedef Grid::const_iterator                        GridIterConst;


            virtual ~ODESolver() { }


            virtual int         integrate()                            = 0;
            virtual int         integrate(unsigned n, double* yIni,
                                          double xLeft, double xRight) = 0;
            virtual Grid&       getAdaptiveGridPoints()                = 0;
            virtual Trajectory& getAdaptiveSolution()                  = 0;


            void setRTol(Real rtol) { _rtol = rtol; }
            void setATol(Real atol) { _atol = atol; }

            Real getRTol() { return _rtol; }
            Real getATol() { return _atol; }


        protected:
            ODESolver() : _atol(10*EPMACH), _rtol(1.0e-12) { }

            // The following precision setting seems to work only semi-stable:
            // ODESolver() : _atol(EPMACH), _rtol(1.0e-9) { }

            Real _atol, _rtol;
    };

}
#endif // ODE_SOLVER_H
