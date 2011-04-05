// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-01-19 td
// last changed:
//
#ifndef __QR_MOORE_PENROSE_H
#define __QR_MOORE_PENROSE_H

#include "GenericQRPseudoInv.h"

namespace PARKIN
{
    /// Pseudo-inverse computation by (classic) Moore-Penrose inversion

    class QRMoorePenrose : public GenericQRPseudoInv
    {
        public:
            // default c'tor
            // QRMoorePenrose() { }

            // d'tor
            virtual ~QRMoorePenrose() { }

            // help for copy c'tor / assignment
            virtual QRMoorePenrose* clone() { return new QRMoorePenrose(*this); }

            //

            virtual bool prepare(Matrix const& qrA, Vector const& diag, long rank);
            virtual bool solvePInv(long rank, Vector& x);
            virtual bool solveR(Matrix const& qrA, Vector const& diag,
                                long rank, Vector const& b,
                                Vector& v);

            virtual Real projectIntoSubspace(long rank, Vector& v);
    };
}

#endif // __QR_MOORE_PENROSE_H
