// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-03 td
// last changed:
//
#ifndef __MATRIX_H
#define __MATRIX_H

#include "common/Constants.h"
#include "common/Types.h"

#include "GenericMatrixImpl.h"
// #include "Vector.h"


namespace PARKIN
{
    //
    class Mref;
    class Vector;
    class QRDecomp;
    class QRconDecomp;


    class Matrix
    {
        public:
            Matrix();
            Matrix(long m, long n);
            Matrix(Matrix const& A);
            Matrix(Vector const& v);
            virtual ~Matrix();
            Matrix const& operator= (Matrix const& A);
            //
            long nr() const;
            long nc() const;
            //
            Real&       operator() (long j, long k);
            Real const& operator() (long j, long k) const;
            //
            Matrix t() const;
            //
            Matrix rowm(long j) const;
            Matrix rowm(long j1, long j2) const;
            Vector colm(long k) const;
            Matrix colm(long k1, long k2) const;
            Matrix subm(long j1, long j2, long k1, long k2) const;
            //
            Mref set_rowm(long j);
            Mref set_rowm(long j1, long j2);
            Mref set_colm(long k);
            Mref set_colm(long k1, long k2);
            Mref set_subm(long j1, long j2, long k1, long k2);
            // set/create random matrix
            void randm(long m, long n);
            void zero();
            void zeros(long m, long n);
            void ones(long m, long n);
            void scale_columns(Vector const& v);
            // possible nrmtypes so far: "l1", "l2", "inf", "fro"
            Real norm(std::string nrmtype) const;
            Real trace() const;
            //
            QRconDecomp factorQRcon(
                long mcon = 0,
                long rank = 0,
                Real cond = 1.0/EPMACH,
                int  meth = 0
                ) const;
            QRDecomp factorQR(
                long rank = 0,
                Real cond = 1.0/EPMACH,
                int  meth = 0
                ) const;
            //
            GenericMatrixImpl* impl() const;

        private:
            GenericMatrixImpl* _impl;

            // friend Matrix operator* (Matrix const& A, Matrix const& B);
            // friend std::ostream& operator<< (std::ostream& os, Matrix const& A);
    };

    // ------------------------------------------------------------------------

    class Mref
    {
        public:
            ~Mref() { }
            void operator= (Matrix const& M)
            {
                (*_pM).impl()->write_subm( *M.impl(), _j1,_j2,_k1,_k2);
            }

        private:
            friend class Matrix;

            Mref(Matrix const* M, long j1, long j2, long k1, long k2) :
                _pM(M), _j1(j1), _j2(j2), _k1(k1), _k2(k2)
            { }

            Matrix const*   _pM;
            long            _j1, _j2, _k1, _k2;
    };


    /// =======================================================================

    Matrix operator+ (Matrix const& A, Matrix const& B);
    Matrix operator* (Real const c, Matrix const& A);
    Matrix operator* (Matrix const& A, Real const c);
    Matrix operator* (Matrix const& A, Matrix const& B);
    std::ostream& operator<< (std::ostream& os, Matrix const& A);

    Matrix add(Matrix const& A, Matrix const& B);
    Matrix mult(Matrix const& A, Matrix const& B);
}
#endif // __MATRIX_H
