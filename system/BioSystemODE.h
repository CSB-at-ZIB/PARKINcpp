// Copyright (C) 2010 - 2013
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2013-03-12 td
// last changed:
//
#ifndef __BIO_SYSTEM_ODE_H
#define __BIO_SYSTEM_ODE_H


#include "odelib/FirstOrderODESystem.h"
#include "BioRHS.h"


namespace PARKIN
{

    class BioSystem;

    class BioSystemODE : public FirstOrderODESystem
    {

        typedef std::vector< std::string >::const_iterator  StrIterConst;

        public:
            BioSystemODE();
            virtual ~BioSystemODE();

            virtual void computeDerivatives( Real const t,
                                               Real* y, Real* dy,
                                               int* info
                                             );
            virtual void computeJacobian( Real const t,
                                            Real* y, Real* dy,
                                            Real* J,int* ldJ,
                                            int* full_or_band, int* info
                                          );
            virtual void computeMassMatrix( Real const t, Real* y,
                                              Real* B, int* ir, int* ic
                                            );

            virtual int getSystemDimension() { return _dim; }
            virtual int getMassMatrixNz() { return _nz; }

            ///

            void setObj(BioSystem& obj);

        private:
            int         _dim;
            int         _nz;
            BioSystem*  _obj;
            BioRHS      _rhs;
    };

}

#endif // __BIO_SYSTEM_ODE_H
