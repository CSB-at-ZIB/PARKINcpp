// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-01-19 td
// last changed:
//
#ifndef __GENERIC_QR_PSEUDO_INV_H
#define __GENERIC_QR_PSEUDO_INV_H

#include "common/Types.h"
#include "Matrix.h"
#include "Vector.h"

namespace PARKIN
{
    /// Interface class for computing rank-deficient pseudo-inverse

    class GenericQRPseudoInv
    {
        public:
            // c'tor: take the default one!
            // GenericQRPseudoInv() : _qrAH() { }

            // d'tor
            virtual ~GenericQRPseudoInv() { }

            // help for copy c'tor / assignment
            virtual GenericQRPseudoInv* clone() = 0;

            //

            virtual bool prepare(Matrix const& qrA, Vector const& diag, long rank) = 0;
            virtual bool solvePInv(long rank, Vector& x) = 0;
            virtual bool solveR(Matrix const& qrA, Vector const& diag,
                                long rank, Vector const& b,
                                Vector& v) = 0;

            virtual Real projectIntoSubspace(long rank, Vector& v) = 0;

            //

            Matrix const& getMat() const { return _qrAH; }

        protected:
            Matrix  _qrAH;
    };
}

#endif // __GENERIC_QR_PSEUDO_INV_H
