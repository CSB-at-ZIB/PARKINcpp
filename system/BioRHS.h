// Copyright (C) 2010 - 2011
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2011-02-21 td
// last changed:
//
#ifndef __BIO_RHS_H
#define __BIO_RHS_H

#include <map>
#include <vector>
#include <string>

#include "common/Types.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"

#include "Expression.h"


namespace PARKIN
{

    class BioRHS
    {
        public:

            typedef std::map< std::string, int >            ExprTypeMap;
            typedef std::map< std::string, Expression >     ExpressionMap;
            typedef std::vector< std::string >              Species;
            typedef std::vector< std::string >              Parameter;

            //

            BioRHS();
            // explicit BioRHS(long n);
            BioRHS(ExpressionMap const& emap);
            BioRHS(Species const& species);
            BioRHS(ExpressionMap const& emap, Parameter const& param);

            ~BioRHS();

            BioRHS(BioRHS const& rhs);
            BioRHS const& operator= (BioRHS const& rhs);

            //

            void
            setRHSType(ExprTypeMap const& tmap);
            //
            ExprTypeMap const&
            getRHSType() const { return _rhs_type; }

            //

            void
            setRHS(ExpressionMap const& emap);
            //
            ExpressionMap const&
            getRHS() const { return _rhs; }
            //
            ExpressionMap const&
            getdRHS() const { return _drhs; }

            void
            resetSpecies(Species const& species);
            //
            Species const&
            getSpecies() const;

            void
            setParameters(Expression::Param const& param);
            void
            setParameters(Parameter const& param);
            //
            Parameter const&
            getParameters() const;
            //
            void
            setParBase(double* p) { if (_data_table != 0) _data_table[1] = p; }

            //

            //

            void
            b(double* y, int* nz, double* B, int* ir, int* ic);

            //

            Vector
            f(Expression::Param& par, long n = 0, double* dy = 0);
            Vector
            f(double* y, long n, double* dy);

            Matrix
            Jf(Expression::Param& par, long n = 0, double* J = 0, long ldJ = 0);

            Matrix
            df(Expression::Param const& var, Matrix const& Zp, Expression::Param& par);
            Matrix
            df(Expression::Param const& var, double* Zp, double* y,
               long n, long q, double* dZ);

            //

            typedef ExpressionMap::const_iterator               EMapIterConst;
            typedef ExprTypeMap::const_iterator                 TMapIterConst;
            typedef std::vector< std::string >::const_iterator  StrIterConst;

        private:

            void computeDerivativeExpression();

            ExprTypeMap     _rhs_type;
            ExpressionMap   _rhs;
            ExpressionMap   _drhs;
            Species         _species;
            Parameter       _parameter;
            long            _num_table_entries;
            double**        _data_table;
    };

}
#endif // __BIO_RHS_H
