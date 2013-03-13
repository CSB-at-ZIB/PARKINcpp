// Copyright (C) 2010 - 2012
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2012-11-20 td
// last changed:
//
#ifndef __FIRST_ORDER_ODE_SYSTEM_H
#define __FIRST_ORDER_ODE_SYSTEM_H

#include "common/Types.h"  // for definition of "Real"

namespace PARKIN
{
    class FirstOrderODESystem
    {
        public:
            virtual ~FirstOrderODESystem() { }
            virtual void computeDerivatives( Real const t,
                                               Real* y, Real* dy,
                                               int* info
                                             ) { *info = -987; }
            virtual void computeJacobian( Real const t,
                                            Real* y, Real* dy,
                                            Real* J, int* ldJ,
                                            int* full_or_band, int* info
                                          ) { *info = -987; }
            virtual void computeMassMatrix( Real const t, Real* y,
                                              Real* B, int* ir, int* ic
                                            ) { }

            virtual int getSystemDimension() { return 0; }
            virtual int getMassMatrixNz() { return 0; }
    };
}

#endif // __FIRST_ORDER_ODE_SYSTEM_H
