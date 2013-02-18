// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-02-18 td
// last changed:
//
#ifndef __CUBIC_HERMITE_TRAJECTORY_H
#define __CUBIC_HERMITE_TRAJECTORY_H

#include "ODETrajectory.h"

namespace PARKIN
{

    /*
    struct CubicHermiteData
    {
        double               tLeft, tRight;
        unsigned             nDim;

        CubicHermiteData() :
            tLeft(0.0), tRight(0.0), nDim(0)
        { }
    };
    */

    ///

    class CubicHermiteTrajectory : public ODETrajectory
    {
        public:
            explicit CubicHermiteTrajectory(int n);
            virtual  ~CubicHermiteTrajectory();

            virtual CubicHermiteTrajectory*   clone() { return new CubicHermiteTrajectory(*this); }

            virtual void               clear();
            virtual void               insert(double t, int n, double* y);
            virtual std::vector<Real>  eval(double t);

            //

            void insertHerm(double t, int n, double* y, double* dy);
            void appendHerm(double t, int n, double* y, double* dy);

        private:
            std::vector<Real> evalHermite(double t,
                                            TrajData const& trajLeft,
                                            TrajData const& trajRight);

            unsigned                         _iLow;
            unsigned                         _iHigh;
    };

}
#endif // __CUBIC_HERMITE_TRAJECTORY_H
