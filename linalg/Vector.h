// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-10-27 td
// last changed:
//
#ifndef __VECTOR_H
#define __VECTOR_H

#include <vector>

#include "common/Types.h"

#include "GenericMatrixImpl.h"
// #include "Matrix.h"


namespace PARKIN
{

    class Vref;
    class Matrix;

    class Vector
    {
        public:
            Vector();
            explicit Vector(long m);
            Vector(std::vector<Real> const& v);
            Vector(Vector const& v);
            virtual ~Vector();
            Vector const& operator= (Vector const& v);
            //
            long nr() const;
            long nc() const;
            //
            Real&       operator() (long j);
            Real const& operator() (long j) const;
            //
            Matrix t() const;
            Matrix diag() const;
            //
            Vector row(long j1, long j2) const;
            //
            Vref set_row(long j1, long j2);
            // set/create random vector
            void rand(long m);
            void zero();
            void zeros(long m);
            void ones(long m);
            // possible nrmtypes so far: "l1", "l2", "inf", "fro"
            Real norm(std::string nrmtype) const;
            //Real length_squared() const
            //
            Vector dotmult(Vector const& v) const;
            Vector dotdiv(Vector const& v) const;
            //
            GenericMatrixImpl* impl() const;

        private:
            GenericMatrixImpl* _impl;

            // friend Vector operator+ (Vector const& v, Vector const& w);
            // friend Vector operator- (Vector const& v, Vector const& w);
            // friend Vector operator* (Matrix const& A, Vector const& x);
            // friend Real operator* (Vector const& v, Vector const& w);
            // friend std::ostream& operator<< (std::ostream& os, Vector const& v);

            // friend Real norm(Vector const& v);
            // friend Real length_squared(Vector const& v);
    };

    // ------------------------------------------------------------------------

    class Vref
    {
        public:
            ~Vref() { }
            void operator= (Vector const& v)
            {
                (*_pv).impl()->write_rowm( *v.impl(), _j1,_j2);
            }

        private:
            friend class Vector;

            Vref(Vector const* v, long j1, long j2) :
                _pv(v), _j1(j1), _j2(j2)
            { }

            Vector const*   _pv;
            long            _j1, _j2;
    };

    /// =======================================================================

    Vector operator+ (Vector const& v, Vector const& w);
    Vector operator- (Vector const& v, Vector const& w);
    Vector operator* (Real const c, Vector const& v);
    Vector operator* (Vector const& v, Real const c);
    Vector operator* (Matrix const& A, Vector const& x);
    Matrix operator* (Vector const& x, Matrix const& A);
    Real operator* (Vector const& v, Vector const& w);
    std::ostream& operator<< (std::ostream& os, Vector const& v);

    Real norm(Vector const& v);
    Real length_squared(Vector const& v);

    //
    Vector mult(Matrix const& A, Vector const& x);

}
#endif // __VECTOR_H
