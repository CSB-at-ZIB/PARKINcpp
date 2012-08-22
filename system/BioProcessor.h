// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-09-13 td
// last changed:
//
#ifndef __BIO_PROCESSOR_H
#define __BIO_PROCESSOR_H


#include "common/Types.h"

#include "linalg/Matrix.h"
#include "linalg/QRconDecomp.h"

#include "nonlin/UserFunc.h"
#include "nonlin/GaussNewton.h"
#include "nonlin/YeOldeParkinCore.h"

#include "BioSystem.h"
#include "BioPAR.h"


namespace PARKIN
{
    // class IOpt;

    class BioProcessor
    {
        public:

            typedef std::map< std::string, std::vector<Real> >  TrajectoryMap;
            typedef std::vector< Matrix >                       MatrixList;
            typedef std::vector< QRconDecomp >                  QRconDecompList;

            ///

            // c'tor
            BioProcessor(BioSystem*         biosys,
                         std::string const& method = "parkin");
            // d'tor
            ~BioProcessor();

            // copy c'tor
            BioProcessor(BioProcessor const& other);

            ///

            void
            setProcessingMethod(std::string const& method);

            void
            setLogStream(std::ostream& logstream = std::clog);

            void
            setIOpt(IOpt const& iopt);

            IOpt
            getIOpt();

            void
            setParameterConstraints(Vector const& itrans, Vector const& xlb, Vector const& xub);

            //

            void
            setCurrentParamValues(Expression::Param const& par);

            Expression::Param const&
            getCurrentParamValues();

            void
            setCurrentParamThres(Expression::Param const& par);

            Expression::Param const&
            getCurrentParamThres();

            //

            void
            setCurrentSpeciesThres(Expression::Param const& spe);

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

            MatrixList
            getSensitivityMatrices();

            QRconDecompList
            getSensitivityDecomps();

            //

            int
            identifyParameters(Real xtol = 1.0e-4);

            Expression::Param
            getIdentificationResults();

            //

        protected:

        private:
            // assignment
            BioProcessor const& operator= (BioProcessor const& other); // { return *this; }

            Expression::Param computeParameterScales();
            Expression::Param computeSpeciesScales();
            Matrix computeJac(std::string const& mode, int& ifail);
            Matrix computeJcf(std::string const& mode, int& ifail);
            Real dertransf_p(Real p, long k);
            Real transform_p(Real p, long k);
            Real backtrans_p(Real u, long k);

            BioSystem*                  _biosys;
            BioPAR                      _biopar;
            std::string                 _method;
            std::ostream*                _logstream;
            IOpt                        _iopt;

            BioSystem::Species          _curSpecies;
            Expression::Param           _optPar;
            // Expression::Param           _optIdx;
            Expression::Param           _speThres;
            Expression::Param           _parThres;

            Expression::Param           _linftyModel;

            TrajectoryMap               _trajMap;
            ODETrajectory*              _sensTraj;

            MatrixList                  _sensiMat;
            QRconDecompList             _sensiDcmp;

            GaussNewton                 _nlscon;
            GaussNewtonWk               _nlsconWk;

            YeOldeParkinCore            _parkin;
            YeOldeParkinWk              _parkinWk;

            Vector                      _idResult;
    };

}
#endif // __BIO_PROCESSOR_H
