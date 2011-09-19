// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-09-13 td
// last changed:
//
#ifndef __BIO_PROCESSOR_H
#define __BIO_PROCESSOR_H


#include <common/Types.h>

#include <linalg/Matrix.h>
#include <linalg/QRconDecomp.h>

#include <nonlin/GaussNewton.h>
#include <nonlin/YeOldeParkinCore.h>

#include "BioSystem.h"
#include "BioPAR.h"


namespace PARKIN
{

    class BioProcessor
    {
        public:

            typedef std::map< std::string, std::vector<Real> >  TrajectoryMap;

            ///

            // c'tor
            BioProcessor(BioSystem* biosys, std::string const& method = "parkin");
            // d'tor
            ~BioProcessor();

            // copy c'tor
            BioProcessor(BioProcessor const& other);

            ///

            void
            setProcessingMethod(std::string const& method);

            void
            setIOpt(IOpt const& iopt);

            //

            void
            setCurrentParamValues(Expression::Param const& par);

            Expression::Param const&
            getCurrentParamValues();

            void
            setCurrentParamThres(Expression::Param& par);

            Expression::Param const&
            getCurrentParamThres();

            //

            void
            setCurrentSpeciesThres(Expression::Param& spe);

            Expression::Param const&
            getCurrentSpeciesThres();

            //

            TrajectoryMap
            computeModel();

            TrajectoryMap
            computeSensitivityTrajectories();

            Vector
            getAdaptiveTimepoints();

            TrajectoryMap
            getScaledSensitivityTrajectories();

            //

            int
            prepareDetailedSensitivities(Vector const& tp);

            std::vector<Matrix>
            getSensitivityMatrices();

            std::vector<QRconDecomp>
            getSensitivityDecomps();

            //

            int
            identifyParameters();

            Expression::Param
            getIdentificationResults();

            //

        protected:

        private:
            // assignment
            BioProcessor const& operator= (BioProcessor const& other); // { return *this; }

            Expression::Param computeParameterScales();
            Expression::Param computeSpeciesScales();
            Matrix computeJac(std::string mode, int& ifail);
            Matrix computeJcf(std::string mode, int& ifail);

            BioSystem*                  _biosys;
            BioPAR                      _biopar;
            std::string                 _method;
            IOpt                        _iopt;

            BioSystem::Species          _curSpecies;
            Expression::Param           _optPar;
            // Expression::Param           _optIdx;
            Expression::Param           _speThres;
            Expression::Param           _parThres;

            Expression::Param           _linftyModel;

            TrajectoryMap               _trajMap;
            ODETrajectory*              _sensTraj;

            std::vector<Matrix>         _sensiMat;
            std::vector<QRconDecomp>    _sensiDcmp;

            GaussNewton                 _nlscon;
            GaussNewtonWk               _nlsconWk;

            YeOldeParkinCore            _parkin;
            YeOldeParkinWk              _parkinWk;
    };

}
#endif // __BIO_PROCESSOR_H
