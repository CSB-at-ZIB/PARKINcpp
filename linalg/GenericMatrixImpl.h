// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-10-26 td
// last changed:
//
#ifndef __GENERIC_MATRIX_IMPL_H
#define __GENERIC_MATRIX_IMPL_H

#include "common/Types.h"

namespace PARKIN
{

    class GenericMatrixImpl
    {
        public:
            // d'tor
            virtual ~GenericMatrixImpl() {}

            // assign other matrix to this
            virtual void assign(GenericMatrixImpl const* a) = 0;

            // return number of rows and columns, resp.
            virtual long nr() const = 0;
            virtual long nc() const = 0;

            // access to element (j,k) of matrix (const and non-const version)
            virtual Real const& at(long j, long k) const = 0;
            virtual Real&       at(long j, long k) = 0;

            // swap rows and cols of this
            virtual void transpose() = 0;
            // ... of a matrix and store it to this
            virtual void transpose(GenericMatrixImpl const& a) = 0;

            // take only the diagonal elements of this (and scap the rest!)
            virtual void diag() = 0;
            // ... of a vector (!!) and store it to this
            virtual void diagm(GenericMatrixImpl const& a) = 0;

            // copy specific rows / cols to this (destructive!)
            virtual void rowm(GenericMatrixImpl const& a, long j          ) = 0;
            virtual void rowm(GenericMatrixImpl const& a, long j1, long j2) = 0;
            virtual void colm(GenericMatrixImpl const& a, long k          ) = 0;
            virtual void colm(GenericMatrixImpl const& a, long k1, long k2) = 0;
            virtual void subm(GenericMatrixImpl const& a, long j1, long j2, long k1, long k2) = 0;

            // set all entries to zero
            virtual void zero() = 0;

            // set (new) (m x n)-matrix to zero and one, resp.
            virtual void zeros(long m, long n) = 0;
            virtual void ones(long m, long n) = 0;

            // set (new) (m x n)-matrix to random values
            virtual void randm(long m, long n) = 0;

            // scale the columns of matrix by a vector
            virtual void scale_columns(GenericMatrixImpl const& a) = 0;

            // compute the matrix norm(nrmtype) with nrmtype: "l1", "l2", "inf", "fro"
            virtual Real norm(std::string nrmtype) const = 0;
            virtual Real length_squared() const = 0;

            // compute the trace, i.e. the sum of the diagonal elements of the matrix
            virtual Real trace() const = 0;

            // add two matrices and write the result to this
            virtual void add(GenericMatrixImpl const& a, GenericMatrixImpl const& b) = 0;
            // add two matrices and write the result to this
            virtual void subtract(GenericMatrixImpl const& a, GenericMatrixImpl const& b) = 0;
            // multiply two matrices and write the result to this
            virtual void mult(Real const c, GenericMatrixImpl const& a) = 0;
            virtual void mult(GenericMatrixImpl const& a, Real const c) = 0;
            virtual void mult(GenericMatrixImpl const& a, GenericMatrixImpl const& b) = 0;

            // multiply two matrices (of same size) element- aka point-wise, and write the result to this
            virtual void pointwise_multiply(GenericMatrixImpl const& a, GenericMatrixImpl const& b) = 0;
            // divide two matrices (of same size) element by element (aka point-wise), and write the result to this
            virtual void pointwise_divide(GenericMatrixImpl const& a, GenericMatrixImpl const& b) = 0;

            // ...

            // cast a GenericMatrixImpl to its derived classes (const and non-const)
            template <typename T>
            T const& down_cast() const
            {
                return dynamic_cast< T const& >( *instance() );
            }
            template <typename T>
            T& down_cast()
            {
                return dynamic_cast< T& >( *instance() );
            }
            virtual GenericMatrixImpl const* instance() const { return this; }
            virtual GenericMatrixImpl*       instance()       { return this; }

            // put this to output stream
            virtual std::ostream& put(std::ostream& os) const = 0;

            // write given rows / cols of some matrix to this
            virtual void write_rowm(GenericMatrixImpl const& a, long j          ) = 0;
            virtual void write_rowm(GenericMatrixImpl const& a, long j1, long j2) = 0;
            virtual void write_colm(GenericMatrixImpl const& a, long k          ) = 0;
            virtual void write_colm(GenericMatrixImpl const& a, long k1, long k2) = 0;
            virtual void write_subm(GenericMatrixImpl const& a, long j1, long j2, long k1, long k2) = 0;

        protected:
            //  c'tor (no private members, thus no data to initialise!)
            GenericMatrixImpl() {}

        private:
            // no copy c'tor
            GenericMatrixImpl(GenericMatrixImpl const&) {}
            // no assignment operator
            GenericMatrixImpl const& operator= (GenericMatrixImpl const&) { return *this; }

    };

}
#endif // __GENERIC_MATRIX_IMPL_H
