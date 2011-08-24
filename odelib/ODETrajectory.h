// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-08-23 td
// last changed:
//
#ifndef __ODE_TRAJECTORY_H
#define __ODE_TRAJECTORY_H

#include <vector>
#include "common/Types.h"


namespace PARKIN
{

    struct TrajData
    {
        Real                t;
        std::vector<Real>   y;

        TrajData() :
            t(0.0), y()
        { }
    };

    ///

    class ODETrajectory
    {
        public:

            virtual ~ODETrajectory() { }

            virtual void               clear()                             = 0;
            virtual void               insert(double t, int n, double* y)  = 0;
            virtual std::vector<Real>  eval(double t)                      = 0;

            void setDim(int n) { _n = n; }

        protected:
            ODETrajectory() : _n(0), _trajectory() { }

            int                      _n;
            std::vector< TrajData >  _trajectory;
    };

}
#endif // __ODE_TRAJECTORY_H
