// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-12-17 td
// last changed:
//
#ifndef __BIO_PAR_H
#define __BIO_PAR_H

#include "common/Types.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"

#include "nonlin/GaussNewton.h"

#include "BioSystem.h"


namespace PARKIN
{

    class BioPAR : public UserFunc
    {
        public:

            typedef double (*PTrafo)(double);
            typedef double (*InvPTrafo)(double);

            ///

            BioPAR(BioSystem* biosys);
            BioPAR(BioSystem* biosys, BioSystem::Parameter const& param);

            virtual ~BioPAR();

            //

            virtual Vector fcn(Vector const& x, int& ifail);
            virtual Matrix jac(Vector const& x, int& ifail);

            //

            void
            setParameter(BioSystem::Parameter const& param)
            {
                // _bioSystem->setParameters(param);
                _parameter = param;
                _optPar.clear();
            }

            BioSystem::Parameter const&
            getParameter()
            { return _bioSystem->getParameters(); }

            ///

            Vector
            getMeasurements()
            { return _bioSystem->getMeasurements(); }

            Vector
            getMeasurementWeights()
            { return _bioSystem->getMeasurementWeights(); }

//            BioSystem::MeasurementList
//            getSensitivityList()
//            { return _bioSystem->getSensitivityList(); }

        private:


            BioSystem*              _bioSystem;
            BioSystem::Parameter    _parameter;
            Expression::Param       _optPar;
            // PTrafo                  _pTrafo;
            // PTrafo                  _dpTrafo;
            // InvPTrafo               _invPTrafo;
    };

    ///

//    double id(double x);
//    double did(double /*x*/);
//    double invid(double x);
//
//    double posi(double x);
//    double dposi(double x);
//    double invposi(double x);

}

#endif // __BIO_PAR_H
