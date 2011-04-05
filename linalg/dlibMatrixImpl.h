// Copyright (C) 2010
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2010-11-03 td
// last changed:
//
#ifndef __DLIB_MATRIX_IMPL_H
#define __DLIB_MATRIX_IMPL_H

#define DLIB_USE_BLAS
//#if defined(_WIN32) || defined(_WIN64)
//    #undef DLIB_USE_BLAS
//#endif
#include <addpkg/dlib/matrix.h>

#include "GenericMatrixImpl.h"

namespace PARKIN
{

    typedef dlib::matrix<Real>      dMatrix;
    typedef dlib::matrix<Real,0,1>  dVector;

    //

    class dlibMatrixImpl : public GenericMatrixImpl
    {
        public:
            //
            dlibMatrixImpl() : _A(dMatrix()) {}
            dlibMatrixImpl(long m, long n) : _A(dMatrix(m,n)) {}
            virtual ~dlibMatrixImpl() {}

            //
            virtual long nr() const { return _A.nr(); }
            virtual long nc() const { return _A.nc(); }

            //
            virtual void assign(GenericMatrixImpl const* a)
            {
                _A = (*a).down_cast<dlibMatrixImpl>().ref();
            }

            //
            virtual void transpose()
            {
                dlib::trans(_A);
            }
            virtual void transpose(GenericMatrixImpl const& a)
            {
                _A = dlib::trans( a.down_cast<dlibMatrixImpl>().ref() );
            }

            //
            virtual void diag()
            {
                dlib::diag(_A);
            }
            virtual void diagm(GenericMatrixImpl const& a)
            {
                _A = dlib::diagm( a.down_cast<dlibMatrixImpl>().ref() );
            }

            //
            virtual Real&       at(long j, long k)       { return _A(j-1,k-1); }
            virtual Real const& at(long j, long k) const { return _A(j-1,k-1); }

            //
            virtual void zero()                { _A = 0.0; }
            virtual void randm(long m, long n) { _A = dlib::randm(m,n); }
            virtual void zeros(long m, long n) { _A = dlib::zeros_matrix<Real>(m,n); }
            virtual void ones(long m, long n)  { _A = dlib::ones_matrix<Real>(m,n); }

            virtual void scale_columns(GenericMatrixImpl const& a)
            {
                dlib::scale_columns( _A, a.down_cast<dlibMatrixImpl>().ref() );
            }

            virtual void colm(GenericMatrixImpl const& a, long k)
            {
                _A = dlib::colm(
                                    a.down_cast<dlibMatrixImpl>().ref(),
                                    k-1
                               );
            }
            virtual void colm(GenericMatrixImpl const& a, long k1, long k2)
            {
                _A = dlib::colm(
                                    a.down_cast<dlibMatrixImpl>().ref(),
                                    dlib::range(k1-1,k2-1)
                               );
            }
            virtual void rowm(GenericMatrixImpl const& a, long j)
            {
                _A = dlib::rowm(
                                    a.down_cast<dlibMatrixImpl>().ref(),
                                    j-1
                               );
            }
            virtual void rowm(GenericMatrixImpl const& a, long j1, long j2)
            {
                _A = dlib::rowm(
                                    a.down_cast<dlibMatrixImpl>().ref(),
                                    dlib::range(j1-1,j2-1)
                               );
            }
            virtual void subm(GenericMatrixImpl const& a, long j1, long j2, long k1, long k2)
            {
                _A = dlib::subm(
                                    a.down_cast<dlibMatrixImpl>().ref(),
                                    dlib::range(j1-1,j2-1), dlib::range(k1-1,k2-1)
                               );
            }

            //
            virtual Real norm(std::string nrmtype) const
            {
                if ( nrmtype == "l1" ) return dlib::max( dlib::sum_rows(dlib::abs(_A)) );
                else
                if ( nrmtype == "l2" ) { dMatrix U,S,V; dlib::svd2(false,false,_A,U,S,V); return dlib::max(S); }
                else
                if ( nrmtype == "inf" ) return dlib::max( dlib::sum_cols(dlib::abs(_A)) );
                else
                if ( nrmtype == "fro" ) return sqrt( dlib::sum(dlib::squared(dlib::abs(_A))) );
                else
                    return -1.0;
            }
            virtual Real length_squared() const
            {
                return dlib::sum( dlib::squared(_A) );
            }

            //
            virtual Real trace() const
            {
                return dlib::trace(_A);
            }

            //
            virtual void add(GenericMatrixImpl const& a, GenericMatrixImpl const& b)
            {
                _A = a.down_cast<dlibMatrixImpl>().ref() + b.down_cast<dlibMatrixImpl>().ref();
            }
            //
            virtual void subtract(GenericMatrixImpl const& a, GenericMatrixImpl const& b)
            {
                _A = a.down_cast<dlibMatrixImpl>().ref() - b.down_cast<dlibMatrixImpl>().ref();
            }
            //
            virtual void mult(Real const c, GenericMatrixImpl const& a)
            {
                _A = c * a.down_cast<dlibMatrixImpl>().ref();
            }
            virtual void mult(GenericMatrixImpl const& a, Real const c)
            {
                _A = a.down_cast<dlibMatrixImpl>().ref() * c;
            }
            virtual void mult(GenericMatrixImpl const& a, GenericMatrixImpl const& b)
            {
                _A = a.down_cast<dlibMatrixImpl>().ref() * b.down_cast<dlibMatrixImpl>().ref();
            }

            //
            virtual void pointwise_multiply(GenericMatrixImpl const& a, GenericMatrixImpl const& b)
            {
                _A = dlib::pointwise_multiply(
                                                a.down_cast<dlibMatrixImpl>().ref(),
                                                b.down_cast<dlibMatrixImpl>().ref()
                                             );
            }
            virtual void pointwise_divide(GenericMatrixImpl const& a, GenericMatrixImpl const& b)
            {
                _A = dlib::pointwise_multiply(
                                                a.down_cast<dlibMatrixImpl>().ref(),
                                                dlib::reciprocal( b.down_cast<dlibMatrixImpl>().ref() )
                                             );
            }


            //
            virtual std::ostream& put(std::ostream& os) const { return dlib::operator<<(os,_A); }

            ///

            virtual void write_rowm(GenericMatrixImpl const& a, long j)
            {
                dlib::set_rowm(
                                    _A,
                                    j-1
                              )
                    = a.down_cast<dlibMatrixImpl>().ref();
            }
            virtual void write_rowm(GenericMatrixImpl const& a, long j1, long j2)
            {
                dlib::set_rowm(
                                    _A,
                                    dlib::range(j1-1,j2-1)
                              )
                    = a.down_cast<dlibMatrixImpl>().ref();
            }

            virtual void write_colm(GenericMatrixImpl const& a, long k)
            {
                dlib::set_colm(
                                    _A,
                                    k-1
                              )
                    = a.down_cast<dlibMatrixImpl>().ref();
            }
            virtual void write_colm(GenericMatrixImpl const& a, long k1, long k2)
            {
                dlib::set_colm(
                                    _A,
                                    dlib::range(k1-1,k2-1)
                              )
                    = a.down_cast<dlibMatrixImpl>().ref();
            }

            virtual void write_subm(GenericMatrixImpl const& a, long j1, long j2, long k1, long k2)
            {
                dlib::set_subm(
                                    _A,
                                    dlib::range(j1-1,j2-1),
                                    dlib::range(k1-1,k2-1)
                              )
                    = a.down_cast<dlibMatrixImpl>().ref();
            }

        private:
            dMatrix const& ref() const { return _A; }

            dlibMatrixImpl(dlibMatrixImpl const&) { }
            dlibMatrixImpl const& operator= (dlibMatrixImpl const&) { return *this; }

            dMatrix _A;
    };

}

#endif // __DLIB_MATRIX_IMPL_H
