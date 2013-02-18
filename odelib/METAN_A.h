// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-01-22 td
// last changed:
//
#ifndef __METAN_A_H
#define __METAN_A_H

#include "FirstOrderODESystem.h"
#include "ODESolver.h"


namespace PARKIN
{

    class METAN_A : public ODESolver
    {
        public:
            METAN_A();
            virtual ~METAN_A();

            virtual void setODESystem(
                                        FirstOrderODESystem& ode,
                                        double              t0,
                                        Grid const&         y0,
                                        double              tEnd
                                       );


            virtual int integrate();
            virtual int integrate(unsigned n, double* yIni,
                                   double tLeft, double tRight);

            virtual Grid&       getAdaptiveGridPoints();
            virtual Trajectory& getAdaptiveSolution();

            Grid&           getSolutionGridPoints()   { return _solPoints; }
            Trajectory&     getSolutionTrajectory()   { return _solution; }
            Grid&           getDataGridPoints()       { return _datPoints; }
            Trajectory&     getDataTrajectory()       { return _data; }


        private:
            Grid        _solPoints;
            Trajectory  _solution;
            Grid        _datPoints;
            Trajectory  _data;
    };

}

#endif // __METAN_A_H
