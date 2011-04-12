// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 12.03.2010 td
// last changed:
//

#ifndef __QR_DECOMP_H
#define __QR_DECOMP_H

#include "common/Constants.h"

#include "GenericDecomp.h"
#include "QRCholesky.h"
#include "QRMoorePenrose.h"

namespace PARKIN
{

    class QRDecomp : public GenericDecomp
    {
            // enum PInv
            // { QRCholesky, QRWilkinson, QRPenrose };

        public:
            // C'tors
            explicit QRDecomp(Real condMax = 1.0/sqrtEPMACH, int method = 0) :
                _qrA(), /* _qrAH(),*/ _diag(), _pivot(),
                _rank(0), _rankMax(0), _cond(0.0), _condMax(condMax),
                _pinv(0)
            { if (method == 0) _pinv = new QRCholesky(); else _pinv = new QRMoorePenrose(); }

            explicit QRDecomp(long rankMax, Real condMax = 1.0/sqrtEPMACH, int method = 0) :
                _qrA(), /* _qrAH(),*/ _diag(), _pivot(),
                _rank(rankMax), _rankMax(rankMax), _cond(0.0), _condMax(condMax),
                _pinv(0)
            { if (method == 0) _pinv = new QRCholesky(); else _pinv = new QRMoorePenrose(); }

            // Destructor
            virtual ~QRDecomp() { delete _pinv; _pinv = 0; }

            // decompose with Householder reflections
            virtual ErrType const& decompose(Matrix const& A);

            // Solve for right-hand side vector b
            virtual ErrType const& solve(Matrix const& A, Vector const& b, Vector& x);
            virtual ErrType const& solve(Vector& b, Vector& x) const;

            // Extra methods replacing the kred < 0 calls of deccon/solcon
            ErrType const& solveR(Vector const& b, Vector& x) const;
            ErrType const& setNewRank(long rank);

            //
            Real projectIntoSubspace(Vector const& v) const;

            //
            void getFirstFactor(Matrix& Q) const;
            void getSecondFactor(Matrix& R) const;
            void getThirdFactor(Matrix& P) const;
            Matrix getFirstFactor() const;
            Matrix getSecondFactor() const;
            Matrix getThirdFactor() const;

            //
            void setPseudoInverseMethod(GenericQRPseudoInv* pinv)
            { delete _pinv; _pinv = pinv; }

            //
            void setMaxRank(long rankMax)
            { _rankMax = rankMax; }

            //
            void setMaxCond(Real condMax)
            { _condMax = condMax; }

            // Return pseudo-rank of matrix
            long getRank() const
            { return _rank; }

            // Return sub-condition of least squares part
            virtual Real getSubCond() const
            { return _cond; }

            // Return diagonal part Rjj of decomposition
            Vector const& getDiag() const
            { return _diag; }

            // Return index vector storing permutation of columns
            Vector const& getPivot() const
            { return _pivot; }

            Matrix const& getMatH() const
            { return _pinv -> getMat(); }

            Matrix const& getMat() const
            { return _qrA; }

            // Copy constructor
            QRDecomp(QRDecomp const& qrdcmp) :
                GenericDecomp(qrdcmp)
            {
                if (this != &qrdcmp)
                {
                    _qrA   = qrdcmp._qrA;
                    // _qrAH  = qrdcmp._qrAH;
                    _diag  = qrdcmp._diag;
                    _pivot = qrdcmp._pivot;

                    _rank    = qrdcmp._rank;
                    _rankMax = qrdcmp._rankMax;
                    _cond    = qrdcmp._cond;
                    _condMax = qrdcmp._condMax;

                    _pinv = qrdcmp._pinv -> clone();
                    // NOTE: No (!!!) delete here since we are
                    // just c'ting a new object; hence _pinv
                    // could point to anything ...
                }
            }

            // Assignment operator
            QRDecomp const& operator= (QRDecomp const& qrdcmp)
            {
                // if (_pinv == 0) { _pinv = new QRCholesky(); }
                if (this != &qrdcmp)
                {
                    _valid = qrdcmp._valid;

                    _qrA   = qrdcmp._qrA;
                    // _qrAH  = qrdcmp._qrAH;
                    _diag  = qrdcmp._diag;
                    _pivot = qrdcmp._pivot;

                    _rank    = qrdcmp._rank;
                    _rankMax = qrdcmp._rankMax;
                    _cond    = qrdcmp._cond;
                    _condMax = qrdcmp._condMax;

                    if (_pinv != 0) delete _pinv;

                    _pinv = qrdcmp._pinv -> clone();
                }
                return *this;
            }

        private:
//            // Copy constructor
//            QRDecomp(QRDecomp const&) { }
//            // Assignment operator
//            QRDecomp const& operator= (QRDecomp const&) { return *this; }

            //
            bool computeHouseholderReflections(long m, long n);
            bool computePseudoInverse();

        protected:
            Matrix  _qrA;
            // Matrix  _qrAH;
            Vector  _diag;
            Vector  _pivot;

            long    _rank;
            long    _rankMax;
            Real    _cond;
            Real    _condMax;

            GenericQRPseudoInv* _pinv;
    };

}
#endif
