// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-01-19 td
// last changed:
//
#ifndef __QR_CHOLESKY_H
#define __QR_CHOLESKY_H

#include "GenericQRPseudoInv.h"

namespace PARKIN
{
    /// Pseudo-inverse computation by Cholesky decomposition

    class QRCholesky : public GenericQRPseudoInv
    {
        public:
            // default c'tor
            // QRCholesky() { }

            // d'tor
            virtual ~QRCholesky() { }

            // help for copy c'tor / assignment
            virtual QRCholesky* clone() { return new QRCholesky(*this); }

            //

            virtual bool prepare(Matrix const& qrA, Vector const& diag, long rank);
            virtual bool solvePInv(long rank, Vector& x);
            virtual bool solveR(Matrix const& qrA, Vector const& diag,
                                long rank, Vector const& b,
                                Vector& v);

            virtual Real projectIntoSubspace(long rank, Vector& v);
    };
}

#endif // __QR_CHOLESKY_H
