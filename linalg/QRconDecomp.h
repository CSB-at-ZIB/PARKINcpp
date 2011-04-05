// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-05 td
// last changed:
//
#ifndef __QRCON_DECOMP_H
#define __QRCON_DECOMP_H

#include "common/Constants.h"
#include "QRDecomp.h"

namespace PARKIN
{

    class QRconDecomp : public QRDecomp
    {
            // enum PInv
            // { QRCholesky, QRWilkinson, QRPenrose };

        public:
            // C'tors
            explicit QRconDecomp(long mcon = 0, Real condMax = 1.0/sqrtEPMACH, int method = 0) :
                QRDecomp(0, condMax, method),
                _mcon(mcon), _rankc(mcon), _rankcMax(mcon), _condc(0.0)
            { }
            explicit QRconDecomp(long mcon, long rankMax, Real condMax = 1.0/sqrtEPMACH, int method = 0) :
                QRDecomp(rankMax, condMax, method),
                _mcon(mcon), _rankc(mcon), _rankcMax(mcon), _condc(0.0)
            { }

            // Destructor
            virtual ~QRconDecomp() { }

            // decompose with Householder reflections
            virtual ErrType const& decompose(Matrix const& A);

            // Solve for right-hand side vector b
            virtual ErrType const& solve(Matrix const& A, Vector const& b, Vector& x);
            virtual ErrType const& solve(Vector& b, Vector& x) const;

            // Extra methods replacing the kred < 0 calls of deccon/solcon
            ErrType const& solveR(Vector const& b, Vector& x) const;
            ErrType const& setNewRank(long rank);

            //
            void getFirstFactor(Matrix& Q) const;
            void getSecondFactor(Matrix& R) const;

            //
            void setMaxRankc(long rankcMax)
            { _rankcMax = rankcMax; }

            // Return pseudo-rank of matrix
            long getRankc() const
            { return _rankc; }

            // Return sub-condition of least squares part
            virtual Real getSubCondc() const
            { return _condc; }


        private:
//            // Copy constructor
//            QRconDecomp(QRconDecomp const&) { }
//            // Assignment operator
//            QRconDecomp const& operator= (QRconDecomp const&) { return *this; }

            //
            bool computeHouseholderReflections(long m, long n);
            bool computePseudoInverse();

        protected:
            long    _mcon;
            long    _rankc;
            long    _rankcMax;
            Real    _condc;
    };

}
#endif // __QRCON_DECOMP_H
