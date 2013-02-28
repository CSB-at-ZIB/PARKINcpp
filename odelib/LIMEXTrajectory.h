// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-08-23 td
// last changed:
//
#ifndef __LIMEX_TRAJECTORY_H
#define __LIMEX_TRAJECTORY_H

#include "ODETrajectory.h"

namespace PARKIN
{

    struct LIMEXHermiteData
    {
        double               t1, t2;
        unsigned             kOrder, nDim;
        std::vector< Real >  coeff;

        LIMEXHermiteData() :
            t1(0.0), t2(0.0), kOrder(0), nDim(0), coeff()
        { }
    };

    ///

    class LIMEXTrajectory : public ODETrajectory
    {
        public:
            explicit LIMEXTrajectory(int n);
            virtual  ~LIMEXTrajectory();

            virtual LIMEXTrajectory*   clone() { return new LIMEXTrajectory(*this); }

            virtual void               clear();
            virtual void               insert(double t, int n, double* y);
            virtual std::vector<Real>  eval(double t);

            //

            void insertHerm(double t, int n, double* y);
            void appendHerm(double t, int n, double* y,
                            int k, int N, double coeff[],
                            double t1, double t2);

        private:
            std::vector<Real> evalHermite(double t, LIMEXHermiteData const& herm);

            unsigned                         _iLow;
            unsigned                         _iHigh;
            std::vector< LIMEXHermiteData >  _hermite;

    };

}
#endif // __LIMEX_TRAJECTORY_H
