// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-04-04 td
// last changed:
//
#ifndef __LINEAR_TRAJECTORY_H
#define __LINEAR_TRAJECTORY_H

#include "ODETrajectory.h"

namespace PARKIN
{

    /*
    struct LinearData
    {
        double               tLeft, tRight;
        unsigned             nDim;

        LinearData() :
            tLeft(0.0), tRight(0.0), nDim(0)
        { }
    };
    */

    ///

    class LinearTrajectory : public ODETrajectory
    {
        public:
            explicit LinearTrajectory(int n);
            virtual  ~LinearTrajectory();

            virtual LinearTrajectory*   clone() { return new LinearTrajectory(*this); }

            virtual void               clear();
            virtual void               insert(double t, int n, double* y);
            virtual std::vector<Real>  eval(double t);

            //

            void insertLin(double t, int n, double* y);
            void appendLin(double t, int n, double* y);

        private:
            std::vector<Real> evalLinear(double t,
                                            TrajData const& trajLeft,
                                            TrajData const& trajRight);

            unsigned                       _iLow;
            unsigned                       _iHigh;
    };

}
#endif // __LINEAR_TRAJECTORY_H
